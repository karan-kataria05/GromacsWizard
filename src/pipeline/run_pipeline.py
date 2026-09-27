#!/usr/bin/env python3
"""
GROMACS Protein-Ligand MD Pipeline Automation (CHARMM36-CGenFF/SwissParam, GPU)
=================================================================================
Matches your Chimera + SwissParam + GROMACS SOP. Automates everything from
"copy the SwissParam zip contents + REC.pdb into a working folder" onward --
the Chimera/SwissParam ligand prep itself (steps 1-5 of your SOP) is manual/GUI
and can't be scripted.

Assumes you already have, in your working directory:
  - REC.pdb          : DockPrep'd protein structure (chain-cleaned)
  - LIG.pdb          : ligand coordinates, same frame as REC.pdb
  - LIG.itp          : ligand topology (from the SwissParam zip)
  - LIG.prm          : optional extra parameters, if SwissParam generated one
  - <name>.ff/        : the SwissParam-downloaded force field folder (any name,
                         auto-detected -- e.g. 'charmm36-feb2026_cgenff-5.0.ff')

Usage:
  python3 run_pipeline.py --prod-ns 100
  python3 run_pipeline.py --prod-ns 100 --step md          # run/resume one step
  python3 run_pipeline.py --prod-ns 100 --force            # rerun everything

Full step order: prep -> box -> solvate -> ions -> em -> index -> nvt -> npt -> md -> analysis
Every step is checkpointed -- rerun the same command after a crash/Ctrl+C and it
picks up where it left off.
"""

import argparse
import json
import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).resolve().parent
MDP_DIR = SCRIPT_DIR / "mdp"

STEPS = ["prep", "box", "solvate", "ions", "em", "index", "nvt", "npt", "md", "analysis"]

GROUP_LINE_RE = re.compile(r'^\s*(\d+)\s+(\S+)\s*:\s*\d+\s+atoms?', re.MULTILINE | re.IGNORECASE)

# Matches GROMACS mdrun -v progress lines, e.g.:
#   "step 1234000, remaining wall clock time:   842 s"
STEP_RE = re.compile(r'step\s+(\d+)(?:,\s*remaining wall clock time:\s*([\d.]+)\s*s)?', re.IGNORECASE)


def read_nsteps_from_mdp(mdp_path: Path):
    """Pull the 'nsteps = N' value out of an .mdp file, for driving a progress bar."""
    for line in mdp_path.read_text().splitlines():
        line = line.split(";")[0].strip()
        if line.lower().startswith("nsteps"):
            return int(line.split("=")[1].strip())
    return None


# --------------------------------------------------------------------------
# Infra: logging, running commands, checkpointing
# --------------------------------------------------------------------------

def setup_logging(workdir: Path):
    log_path = workdir / "pipeline.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    return logging.getLogger("pipeline")


def run(cmd, log, cwd=None, input_text=None):
    log.info(f"RUNNING: {' '.join(cmd)}" + (f"  (stdin: {input_text!r})" if input_text else ""))
    result = subprocess.run(cmd, cwd=cwd, input=input_text, text=True, capture_output=True)
    if result.stdout:
        log.info(result.stdout[-4000:])
    if result.returncode != 0:
        log.error(result.stderr[-5000:])
        raise RuntimeError(
            f"Command failed (exit {result.returncode}): {' '.join(cmd)}\nSee pipeline.log for full output."
        )
    if result.stderr:
        log.debug(result.stderr[-1000:])
    return result


def run_streaming(cmd, log, cwd=None, total_steps=None):
    """Like run(), but streams output live -- use for mdrun, where the job can
    take hours. If total_steps is known (from the .mdp's nsteps), shows a tqdm
    progress bar driven off GROMACS's own 'step N, remaining wall clock time'
    lines instead of dumping raw text. If total_steps is None (e.g. EM, which
    has no fixed step target), falls back to raw live echo."""
    log.info(f"RUNNING (live): {' '.join(cmd)}")
    process = subprocess.Popen(
        cmd, cwd=cwd, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1,
    )

    pbar = None
    if total_steps:
        pbar = tqdm(total=total_steps, unit="step",
                    desc=Path(cwd).name if cwd else "mdrun")
    last_step = 0

    for line in process.stdout:
        logging.getLogger("pipeline").debug(line.rstrip())
        if pbar is not None:
            m = STEP_RE.search(line)
            if m:
                step = int(m.group(1))
                pbar.update(max(0, step - last_step))
                last_step = step
                if m.group(2):
                    pbar.set_postfix(remaining_s=m.group(2))
        else:
            print(line, end="")  # raw live echo when we don't know the step target

    if pbar is not None:
        pbar.close()

    process.wait()
    if process.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {process.returncode}): {' '.join(cmd)}\n"
            f"Scroll up in the terminal, or check the step's own .log file (e.g. EM.log) for details."
        )
    return process.returncode


def skip_if_exists(path: Path, log, force=False):
    if path.exists() and not force:
        log.info(f"SKIP (already exists): {path}")
        return True
    return False


# --------------------------------------------------------------------------
# Group-index parsing helpers (for make_ndx / genrestr / analysis tools)
# --------------------------------------------------------------------------

def parse_groups(stdout_text: str) -> dict:
    """Parse a gmx make_ndx group table. Later occurrences win (final state)."""
    groups = {}
    for m in GROUP_LINE_RE.finditer(stdout_text):
        groups[m.group(2)] = int(m.group(1))
    return groups


def max_group_index(stdout_text: str):
    idxs = [int(m.group(1)) for m in GROUP_LINE_RE.finditer(stdout_text)]
    return max(idxs) if idxs else None


def find_group(groups: dict, *candidates):
    """First exact (case-insensitive) match, else first substring match."""
    lower_map = {k.lower(): v for k, v in groups.items()}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    for cand in candidates:
        for name, idx in groups.items():
            if cand.lower() in name.lower():
                return idx
    return None


# --------------------------------------------------------------------------
# Step: prep (pdb2gmx, ligand merge, topology patch)
# --------------------------------------------------------------------------

def find_ff_dir(search_dir: Path) -> Path:
    candidates = [c for c in sorted(search_dir.glob("*.ff")) if c.is_dir()]
    if not candidates:
        raise FileNotFoundError(
            f"No '*.ff' force-field folder found in {search_dir}. "
            "Copy the SwissParam-downloaded force field folder there first."
        )
    if len(candidates) > 1:
        names = ", ".join(c.name for c in candidates)
        raise RuntimeError(f"Multiple *.ff folders found ({names}) -- pass --ff-name explicitly.")
    return candidates[0]


def fix_moleculetype_name(itp_path: Path, log, target_name="LIG"):
    lines = itp_path.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.strip().startswith("[") and "moleculetype" in line.lower():
            j = i + 1
            while j < len(lines) and (not lines[j].strip() or lines[j].strip().startswith(";")):
                j += 1
            if j < len(lines):
                parts = lines[j].split()
                if parts and parts[0] != target_name:
                    old_name = parts[0]
                    lines[j] = lines[j].replace(old_name, target_name, 1)
                    log.info(f"Renamed moleculetype '{old_name}' -> '{target_name}' in {itp_path.name}")
            break
    itp_path.write_text("\n".join(lines) + "\n")


def merge_gro(rec_gro: Path, lig_gro: Path, out_gro: Path, log):
    rec_lines = rec_gro.read_text().splitlines()
    lig_lines = lig_gro.read_text().splitlines()
    rec_n = int(rec_lines[1].strip())
    lig_n = int(lig_lines[1].strip())
    total = rec_n + lig_n
    merged = (
        [rec_lines[0], f"{total:5d}"]
        + rec_lines[2:2 + rec_n]
        + lig_lines[2:2 + lig_n]
        + [rec_lines[-1]]
    )
    out_gro.write_text("\n".join(merged) + "\n")
    log.info(f"Merged {rec_n} protein atoms + {lig_n} ligand atoms = {total} -> {out_gro}")


def patch_topology(top_path: Path, ff_dir_name: str, log):
    lines = top_path.read_text().splitlines()

    ff_idx = next((i for i, l in enumerate(lines) if "forcefield.itp" in l and ff_dir_name in l), None)
    if ff_idx is None:
        ff_idx = next(i for i, l in enumerate(lines) if "forcefield.itp" in l)

    block = ["", "; Include ligand topology", '#include "LIG.itp"']
    for j, b in enumerate(block):
        lines.insert(ff_idx + 1 + j, b)

    posres_idx = next((i for i, l in enumerate(lines) if 'include "posre.itp"' in l), None)
    if posres_idx is not None:
        endif_idx = posres_idx + 1  # the #endif right after #include "posre.itp"
        block = ["", "; Ligand position restraints", "#ifdef POSRES", '#include "posre_LIG.itp"', "#endif"]
        for j, b in enumerate(block):
            lines.insert(endif_idx + 1 + j, b)
    else:
        log.warning("Could not find protein's POSRES block -- add the ligand posre_LIG.itp block manually.")

    lines.append("LIG                 1")
    top_path.write_text("\n".join(lines) + "\n")
    log.info(f"Patched topology: {top_path}")


def step_prep(args, dirs, log):
    out_gro = dirs["prep"] / "complex.gro"
    if skip_if_exists(out_gro, log, args.force):
        return

    rec_pdb = Path(args.rec).resolve()
    lig_pdb = Path(args.lig_pdb).resolve()
    lig_itp = Path(args.lig_itp).resolve()
    for f in (rec_pdb, lig_pdb, lig_itp):
        if not f.exists():
            raise FileNotFoundError(f"Required input not found: {f}")

    ff_src = find_ff_dir(Path(args.ff_source).resolve())
    ff_dest = dirs["prep"] / ff_src.name
    if not ff_dest.exists():
        shutil.copytree(ff_src, ff_dest)
    ff_name = ff_src.name[:-3]  # strip '.ff'
    log.info(f"Using force field: {ff_name} (from {ff_src})")

    run(
        ["gmx", "pdb2gmx", "-f", str(rec_pdb), "-o", "conf_protein.gro",
         "-p", "topol.top", "-i", "posre.itp", "-water", "tip3p", "-ff", ff_name, "-ignh"],
        log, cwd=dirs["prep"],
    )

    run(["gmx", "editconf", "-f", str(lig_pdb), "-o", "LIG.gro"], log, cwd=dirs["prep"])

    shutil.copy(lig_itp, dirs["prep"] / "LIG.itp")
    fix_moleculetype_name(dirs["prep"] / "LIG.itp", log)
    if args.lig_prm:
        lig_prm = Path(args.lig_prm)
        if lig_prm.exists():
            shutil.copy(lig_prm, dirs["prep"] / "LIG.prm")

    merge_gro(dirs["prep"] / "conf_protein.gro", dirs["prep"] / "LIG.gro", out_gro, log)
    patch_topology(dirs["prep"] / "topol.top", ff_name, log)

    log.info(
        "\n*** MANUAL CHECKPOINT ***\n"
        "Open prep/topol.top and prep/complex.gro and verify:\n"
        "  1. '#include \"LIG.itp\"' sits right after the forcefield.itp include\n"
        "  2. 'LIG   1' is the last line of [ molecules ]\n"
        "  3. complex.gro's atom count = protein atoms + ligand atoms (see log above)\n"
    )


# --------------------------------------------------------------------------
# Step: box / solvate / ions
# --------------------------------------------------------------------------

def step_box(args, dirs, log):
    out = dirs["prep"] / "box.gro"
    if skip_if_exists(out, log, args.force):
        return
    run(["gmx", "editconf", "-f", "complex.gro", "-o", "box.gro", "-d", "1.0", "-bt", "triclinic"],
        log, cwd=dirs["prep"])


def step_solvate(args, dirs, log):
    out = dirs["prep"] / "box_sol.gro"
    if skip_if_exists(out, log, args.force):
        return
    run(["gmx", "solvate", "-cp", "box.gro", "-cs", "spc216.gro", "-p", "topol.top", "-o", "box_sol.gro"],
        log, cwd=dirs["prep"])


def step_ions(args, dirs, log):
    out = dirs["prep"] / "box_sol_ion.gro"
    if skip_if_exists(out, log, args.force):
        return
    shutil.copy(MDP_DIR / "ions.mdp", dirs["prep"] / "ions.mdp")
    run(["gmx", "grompp", "-f", "ions.mdp", "-c", "box_sol.gro", "-p", "topol.top",
         "-o", "ION.tpr", "-maxwarn", "2"], log, cwd=dirs["prep"])
    run(["gmx", "genion", "-s", "ION.tpr", "-p", "topol.top", "-conc", "0.1",
         "-neutral", "-o", "box_sol_ion.gro"], log, cwd=dirs["prep"], input_text="SOL\n")


# --------------------------------------------------------------------------
# Step: EM
# --------------------------------------------------------------------------

def step_em(args, dirs, log):
    out = dirs["em"] / "EM.gro"
    if skip_if_exists(out, log, args.force):
        return
    shutil.copy(MDP_DIR / "EM.mdp", dirs["em"] / "EM.mdp")
    run(["gmx", "grompp", "-f", "EM.mdp", "-c", "../prep/box_sol_ion.gro", "-p", "../prep/topol.top",
         "-o", "EM.tpr", "-maxwarn", "2"], log, cwd=dirs["em"])
    # EM has no fixed step target (stops on convergence) -- raw live echo.
    run_streaming(["gmx", "mdrun", "-v", "-deffnm", "EM"], log, cwd=dirs["em"])


# --------------------------------------------------------------------------
# Step: index (dynamic group lookup -> posre_LIG.itp + index.ndx)
# --------------------------------------------------------------------------

def step_index(args, dirs, log):
    out = dirs["prep"] / "index.ndx"
    if skip_if_exists(out, log, args.force):
        return

    # 1. Ligand heavy-atom group (for ligand position restraints)
    result = run(["gmx", "make_ndx", "-f", "LIG.gro", "-o", "index_LIG.ndx"],
                 log, cwd=dirs["prep"], input_text="0 & ! a H*\nq\n")
    heavy_idx = max_group_index(result.stdout)
    if heavy_idx is None:
        raise RuntimeError("Could not determine the new heavy-atom group index from make_ndx output.")
    log.info(f"Ligand heavy-atom group created at index {heavy_idx}")

    run(["gmx", "genrestr", "-f", "LIG.gro", "-n", "index_LIG.ndx", "-o", "posre_LIG.itp",
         "-fc", "1000", "1000", "1000"], log, cwd=dirs["prep"], input_text=f"{heavy_idx}\n")

    # 2. Combined Protein|LIG group, built on EM.gro (has all default groups present)
    em_gro = "../em/EM.gro"
    query = run(["gmx", "make_ndx", "-f", em_gro, "-o", "tmp_query.ndx"],
                log, cwd=dirs["prep"], input_text="q\n")
    groups = parse_groups(query.stdout)
    protein_idx = find_group(groups, "Protein")
    lig_idx = find_group(groups, "LIG")
    if protein_idx is None or lig_idx is None:
        raise RuntimeError(f"Could not find Protein/LIG groups in default listing. Found: {groups}")
    log.info(f"Protein group={protein_idx}, LIG group={lig_idx}")

    final = run(["gmx", "make_ndx", "-f", em_gro, "-o", "index.ndx"],
                log, cwd=dirs["prep"], input_text=f"{protein_idx} | {lig_idx}\nq\n")
    final_groups = parse_groups(final.stdout)
    (dirs["prep"] / "groups.json").write_text(json.dumps(final_groups, indent=2))

    log.info(
        "\n*** MANUAL CHECKPOINT ***\n"
        f"Final index.ndx groups: {final_groups}\n"
        "Confirm your NVT.mdp / NPT.mdp / MD.mdp 'tc-grps' line uses names that actually\n"
        "exist above (commonly 'Protein_LIG' and 'Water_and_ions' -- but verify against\n"
        "this exact list before running grompp for NVT).\n"
    )


# --------------------------------------------------------------------------
# Step: NVT / NPT / MD (production)
# --------------------------------------------------------------------------

def gpu_mdrun_cmd(deffnm, args, extra=None):
    cmd = ["gmx", "mdrun", "-v", "-deffnm", deffnm,
           "-nb", "gpu", "-pme", "gpu", "-bonded", "gpu",
           "-gpu_id", str(args.gpu_id), "-nt", str(args.nt)]
    if args.pin:
        cmd += ["-pin", "on"]
    if extra:
        cmd += extra
    return cmd


def step_nvt(args, dirs, log):
    out = dirs["nvt"] / "NVT.gro"
    if skip_if_exists(out, log, args.force):
        return
    shutil.copy(MDP_DIR / "NVT.mdp", dirs["nvt"] / "NVT.mdp")
    run(["gmx", "grompp", "-f", "NVT.mdp", "-c", "../em/EM.gro", "-r", "../em/EM.gro",
         "-p", "../prep/topol.top", "-n", "../prep/index.ndx", "-maxwarn", "2", "-o", "NVT.tpr"],
        log, cwd=dirs["nvt"])
    nvt_nsteps = read_nsteps_from_mdp(dirs["nvt"] / "NVT.mdp")
    run_streaming(gpu_mdrun_cmd("NVT", args), log, cwd=dirs["nvt"], total_steps=nvt_nsteps)


def step_npt(args, dirs, log):
    out = dirs["npt"] / "NPT.gro"
    if skip_if_exists(out, log, args.force):
        return
    shutil.copy(MDP_DIR / "NPT.mdp", dirs["npt"] / "NPT.mdp")
    run(["gmx", "grompp", "-f", "NPT.mdp", "-c", "../nvt/NVT.gro", "-r", "../nvt/NVT.gro",
         "-t", "../nvt/NVT.cpt", "-p", "../prep/topol.top", "-n", "../prep/index.ndx",
         "-maxwarn", "2", "-o", "NPT.tpr"], log, cwd=dirs["npt"])
    npt_nsteps = read_nsteps_from_mdp(dirs["npt"] / "NPT.mdp")
    run_streaming(gpu_mdrun_cmd("NPT", args), log, cwd=dirs["npt"], total_steps=npt_nsteps)


def step_md(args, dirs, log):
    out = dirs["md"] / "MD.gro"
    if skip_if_exists(out, log, args.force):
        cpt = dirs["md"] / "MD.cpt"
        if cpt.exists():
            log.info("MD.gro missing but checkpoint found -- resuming production run.")
            nsteps = read_nsteps_from_mdp(dirs["md"] / "MD.mdp")
            extra = ["-cpi", "MD.cpt"]
            if args.maxh:
                extra += ["-maxh", str(args.maxh)]
            run_streaming(gpu_mdrun_cmd("MD", args, extra=extra), log, cwd=dirs["md"], total_steps=nsteps)
        return

    nsteps = int(args.prod_ns * 500000)  # dt = 0.002 ps
    template = (MDP_DIR / "MD_template.mdp").read_text()
    (dirs["md"] / "MD.mdp").write_text(template.replace("{NSTEPS}", str(nsteps)))

    run(["gmx", "grompp", "-f", "MD.mdp", "-c", "../npt/NPT.gro", "-t", "../npt/NPT.cpt",
         "-p", "../prep/topol.top", "-n", "../prep/index.ndx", "-maxwarn", "2", "-o", "MD.tpr"],
        log, cwd=dirs["md"])
    extra = ["-maxh", str(args.maxh)] if args.maxh else None
    run_streaming(gpu_mdrun_cmd("MD", args, extra=extra), log, cwd=dirs["md"], total_steps=nsteps)


# --------------------------------------------------------------------------
# Publication plotting helpers
# --------------------------------------------------------------------------

PLOT_LEGEND_LABEL = "Abrusogenin-GSK3b"  # Edit this for your system.


def read_xvg(path: Path):
    """Read the first two numeric columns from a GROMACS XVG file."""
    x, y = [], []
    for line in path.read_text(errors="replace").splitlines():
        if not line or line.startswith(("@", "#")):
            continue
        fields = line.split()
        if len(fields) >= 2:
            try:
                x.append(float(fields[0]))
                y.append(float(fields[1]))
            except ValueError:
                continue
    if not x:
        raise RuntimeError(f"No numeric data found in {path}")
    return x, y


def make_publication_plots(analysis_dir: Path, log):
    """Create a five-panel PNG from the completed analysis XVG files."""
    try:
        import matplotlib
        matplotlib.use("Agg")  # headless WSL-safe image generation
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Publication plots need matplotlib. In WSL run: pip3 install matplotlib") from exc

    plt.style.use("ggplot")

    plot_dir = analysis_dir / "Publication plots"
    plot_dir.mkdir(exist_ok=True)

    rmsd_x, rmsd_y = read_xvg(analysis_dir / "rmsd.xvg")
    rmsf_x, rmsf_y = read_xvg(analysis_dir / "rmsf.xvg")
    rg_x, rg_y = read_xvg(analysis_dir / "gyrate1.xvg")   # your script called this rg.xvg
    hb_x, hb_y = read_xvg(analysis_dir / "hb.xvg")        # your script called this hbnum.xvg
    sasa_x, sasa_y = read_xvg(analysis_dir / "sasa.xvg")

    # Convert ps to ns (rg only -- rmsd, hb, sasa already ns via -tu ns)
    rg_x = [v / 1000 for v in rg_x]

    fig, axes = plt.subplots(5, 1, figsize=(12, 22))

    axes[0].plot(rmsd_x, rmsd_y, color="cornflowerblue", linewidth=1.8)
    axes[0].set_title("RMSD", fontsize=16, fontweight="bold")
    axes[0].set_xlabel("Time (ns)"); axes[0].set_ylabel("RMSD (nm)")
    axes[0].legend([PLOT_LEGEND_LABEL]); axes[0].grid(True)

    axes[1].plot(rmsf_x, rmsf_y, color="seagreen", linewidth=1.5)
    axes[1].set_title("RMSF", fontsize=16, fontweight="bold")
    axes[1].set_xlabel("Residues"); axes[1].set_ylabel("RMSF (nm)")
    axes[1].legend([PLOT_LEGEND_LABEL]); axes[1].grid(True)

    axes[2].plot(rg_x, rg_y, color="peru", linewidth=1.5)
    axes[2].set_title("RG", fontsize=16, fontweight="bold")
    axes[2].set_xlabel("Time (ns)"); axes[2].set_ylabel("RG (nm)")
    axes[2].legend([PLOT_LEGEND_LABEL]); axes[2].grid(True)

    axes[3].plot(hb_x, hb_y, color="mediumpurple", linewidth=1.0, alpha=0.8)
    axes[3].set_title("HB", fontsize=16, fontweight="bold")
    axes[3].set_xlabel("Time (ns)"); axes[3].set_ylabel("HB (count)")
    axes[3].legend([PLOT_LEGEND_LABEL]); axes[3].grid(True)

    axes[4].plot(sasa_x, sasa_y, color="goldenrod", linewidth=1.6)
    axes[4].set_title("SASA", fontsize=16, fontweight="bold")
    axes[4].set_xlabel("Time (ns)"); axes[4].set_ylabel("SASA (nm²)")
    axes[4].legend([PLOT_LEGEND_LABEL]); axes[4].grid(True)

    plt.tight_layout()
    output = plot_dir / "MD_Analysis_5Panel.png"
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    log.info(f"Publication plot created: {output}")


def ensure_sasa(args, dirs, log, protein_idx):
    """Generate the one analysis input the previous pipeline did not create."""
    out = dirs["analysis"] / "sasa.xvg"
    if skip_if_exists(out, log, args.force):
        return
    run(["gmx", "sasa", "-s", "../md/MD.tpr", "-f", "MD_center.xtc", "-n", "../prep/index.ndx",
         "-o", "sasa.xvg", "-tu", "ns"], log, cwd=dirs["analysis"], input_text=f"{protein_idx}\n")


# --------------------------------------------------------------------------
# Step: analysis (recenter, RMSD, RMSF, H-bonds, Rg, SASA, energy, plots)
# --------------------------------------------------------------------------

def step_analysis(args, dirs, log):
    out = dirs["analysis"] / "rmsd.xvg"
    groups_path = dirs["prep"] / "groups.json"
    if not groups_path.exists():
        raise RuntimeError("prep/groups.json not found -- run the 'index' step first.")
    groups = json.loads(groups_path.read_text())
    protein_idx = find_group(groups, "Protein")
    lig_idx = find_group(groups, "LIG")
    backbone_idx = find_group(groups, "Backbone")
    system_idx = find_group(groups, "System")

    # Existing analyses can still be upgraded with SASA and the new plot.
    if skip_if_exists(out, log, args.force):
        ensure_sasa(args, dirs, log, protein_idx)
        make_publication_plots(dirs["analysis"], log)
        return

    md_tpr = "../md/MD.tpr"
    md_xtc = "../md/MD.xtc"
    ndx = "../prep/index.ndx"

    # 1. Recenter/rewrap
    run(["gmx", "trjconv", "-s", md_tpr, "-f", md_xtc, "-n", ndx,
         "-o", "MD_center.xtc", "-center", "-pbc", "mol", "-ur", "compact"],
        log, cwd=dirs["analysis"], input_text=f"{protein_idx}\n{system_idx}\n")

    # 2. First frame as a PDB
    run(["gmx", "trjconv", "-s", md_tpr, "-f", "MD_center.xtc", "-n", ndx,
         "-o", "start.pdb", "-dump", "0"], log, cwd=dirs["analysis"], input_text=f"{system_idx}\n")

    # 3. RMSD (fit on Backbone, compute for LIG)
    run(["gmx", "rms", "-s", md_tpr, "-f", "MD_center.xtc", "-n", ndx, "-o", "rmsd.xvg", "-tu", "ns"],
        log, cwd=dirs["analysis"], input_text=f"{backbone_idx}\n{lig_idx}\n")

    # 4. RMSF (Backbone)
    run(["gmx", "rmsf", "-s", md_tpr, "-f", "MD_center.xtc", "-n", ndx, "-o", "rmsf.xvg"],
        log, cwd=dirs["analysis"], input_text=f"{backbone_idx}\n")

    # 5. H-bonds (Protein vs LIG)
    run(["gmx", "hbond", "-s", md_tpr, "-f", "MD_center.xtc", "-n", ndx, "-num", "hb.xvg", "-tu", "ns"],
        log, cwd=dirs["analysis"], input_text=f"{protein_idx}\n{lig_idx}\n")

    # 6. Radius of gyration (Protein)
    run(["gmx", "gyrate", "-s", md_tpr, "-f", "MD_center.xtc", "-n", ndx, "-o", "gyrate1.xvg"],
        log, cwd=dirs["analysis"], input_text=f"{protein_idx}\n")

    # 7. Solvent-accessible surface area for the protein.
    ensure_sasa(args, dirs, log, protein_idx)

    # 8. Energy terms -- name-based selection (GROMACS accepts term names directly
    #    in recent versions). If this errors, open pipeline.log for the printed
    #    term table and switch to numeric indices instead.
    run(["gmx", "energy", "-f", "../md/MD.edr", "-o", "energy1.xvg"],
        log, cwd=dirs["analysis"], input_text="Potential\nTemperature\nPressure\nDensity\n\n")

    make_publication_plots(dirs["analysis"], log)

    log.info(
        "\n*** ANALYSIS COMPLETE ***\n"
        "Files in analysis/: MD_center.xtc, start.pdb, rmsd.xvg, rmsf.xvg, hb.xvg,\n"
        "gyrate1.xvg, sasa.xvg, energy1.xvg, and Publication plots/MD_Analysis_5Panel.png.\n"
    )


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

STEP_FUNCS = {
    "prep": step_prep, "box": step_box, "solvate": step_solvate, "ions": step_ions,
    "em": step_em, "index": step_index, "nvt": step_nvt, "npt": step_npt,
    "md": step_md, "analysis": step_analysis,
}


def main():
    parser = argparse.ArgumentParser(description="GROMACS CHARMM36/CGenFF protein-ligand pipeline")
    parser.add_argument("--rec", default="REC.pdb")
    parser.add_argument("--lig-pdb", default="LIG.pdb")
    parser.add_argument("--lig-itp", default="LIG.itp")
    parser.add_argument("--lig-prm", default=None)
    parser.add_argument("--ff-source", default=".", help="Directory containing the SwissParam *.ff folder")
    parser.add_argument("--prod-ns", type=float, default=100, help="Production run length in ns")
    parser.add_argument("--gpu-id", default="0")
    parser.add_argument("--nt", type=int, default=16, help="Total thread count for mdrun")
    parser.add_argument("--pin", action="store_true", default=True)
    parser.add_argument("--maxh", type=float, default=None, help="Max wall-hours for the MD step (optional)")
    parser.add_argument("--workdir", default=".")
    parser.add_argument("--step", choices=STEPS, default=None)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    workdir = Path(args.workdir).resolve()
    dirs = {name: workdir / name for name in ["prep", "em", "nvt", "npt", "md", "analysis"]}
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    log = setup_logging(workdir)
    log.info(f"Pipeline starting. Production length target: {args.prod_ns} ns")

    for step_name in ([args.step] if args.step else STEPS):
        log.info(f"\n{'='*60}\nSTEP: {step_name}\n{'='*60}")
        try:
            STEP_FUNCS[step_name](args, dirs, log)
        except Exception as e:
            log.error(f"Pipeline stopped at step '{step_name}': {e}")
            # GUI-specific, line-buffer-safe signal.  Keep this separate from
            # normal logs so the GUI can finish even if WSL keeps a pipe open.
            print("PIPELINE_GUI_STATUS: FAILED", flush=True)
            sys.exit(1)

    log.info("Pipeline completed successfully.")
    # A flushed sentinel lets the Windows GUI finish immediately after the
    # scientific work is complete, without waiting for WSL pipe cleanup.
    print("PIPELINE_GUI_STATUS: COMPLETE", flush=True)


if __name__ == "__main__":
    main()
