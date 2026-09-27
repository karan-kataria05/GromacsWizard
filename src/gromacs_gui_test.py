#!/usr/bin/env python3
"""
GROMACS Pipeline GUI Wizard (Windows, CustomTkinter)
=====================================================
A step-by-step graphical front-end for run_pipeline.py, built for a
non-coding user. It never runs `gmx` itself -- it checks WSL/GROMACS are
ready, walks the user through picking input files and run options, then
launches `python3 run_pipeline.py ...` inside WSL and streams the live
output into an on-screen log console with per-step progress bars.

INSTALLED LAYOUT (what the Setup.exe puts down, all under one folder):
  GromacsWizard.exe
  pipeline/
      run_pipeline.py
      mdp/                (ions.mdp, EM.mdp, NVT.mdp, NPT.mdp, MD_template.mdp)
      forcefield/
          charmm36-feb2026_cgenff-5.0.ff/   (bundled, shared by every project)

PER-PROJECT FOLDER (whatever folder the user picks in the wizard) only needs:
  REC.pdb, LIG.pdb, LIG.itp, optionally LIG.prm
The big *.ff force field is no longer copied per-project -- it's read once
from the install location.

Run directly with:            python gromacs_gui_test.py
Install the one dependency:   pip install customtkinter

Build into a standalone .exe (run ON WINDOWS, not in WSL) with Nuitka:
  pip install nuitka customtkinter
  python -m nuitka --onefile --standalone --windows-console-mode=disable ^
      --enable-plugin=tk-inter --windows-icon-from-ico=GromacsWizard.ico ^
      --output-filename=GromacsWizard.exe gromacs_gui_test.py
  -> GromacsWizard.exe (in the current folder)
"""

import os
import json
import re
import shutil
import subprocess
import sys
import threading
import queue
import webbrowser
from pathlib import Path

import customtkinter as ctk
from tkinter import filedialog, messagebox

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

STEPS = ["prep", "box", "solvate", "ions", "em", "index", "nvt", "npt", "md", "analysis"]
STEP_LABELS = {
    "prep": "Prep (pdb2gmx + ligand merge)",
    "box": "Build box",
    "solvate": "Solvate",
    "ions": "Add ions",
    "em": "Energy minimization",
    "index": "Build index groups",
    "nvt": "NVT equilibration",
    "npt": "NPT equilibration",
    "md": "Production MD",
    "analysis": "Analysis (RMSD/RMSF/H-bonds/etc.)",
}

PERCENT_RE = re.compile(r"(\d{1,3})%\|")
STEP_HEADER_RE = re.compile(r"^STEP:\s*(\w+)\s*$")
SIM_STEP_RE = re.compile(r"step\s+(\d+)(?:,\s*remaining wall clock time:\s*([\d.]+)\s*s)?", re.IGNORECASE)

FONT_TITLE = ("Segoe UI", 22, "bold")
FONT_SUBTITLE = ("Segoe UI", 13)
FONT_BODY = ("Segoe UI", 12)
FONT_MONO = ("Consolas", 11)

# Change this one value to the exact author/credit wording you want shown.
CREATOR_CREDIT = "Created by Karan Kataria"

# Store lightweight GUI preferences separately from simulation data.
SETTINGS_PATH = Path(os.environ.get("APPDATA", str(Path.home()))) / "GromacsPipelineWizard" / "settings.json"


def load_last_workdir() -> Path | None:
    """Return the last valid run folder, if one was saved."""
    try:
        saved = json.loads(SETTINGS_PATH.read_text(encoding="utf-8"))
        path = Path(saved.get("last_workdir", ""))
        return path if path.is_dir() else None
    except (OSError, ValueError, TypeError):
        return None


def save_last_workdir(workdir: Path) -> None:
    """Persist the folder selection without making a write failure fatal."""
    try:
        SETTINGS_PATH.parent.mkdir(parents=True, exist_ok=True)
        SETTINGS_PATH.write_text(json.dumps({"last_workdir": str(workdir.resolve())}, indent=2), encoding="utf-8")
    except OSError:
        pass


# --------------------------------------------------------------------------
# WSL helpers (unchanged logic from the console wizard, adapted to be
# non-blocking / callable from a worker thread)
# --------------------------------------------------------------------------

def wsl_run(bash_line, timeout=None):
    return subprocess.run(["wsl", "bash", "-lc", bash_line], text=True,
                           capture_output=True, timeout=timeout)


def windows_to_wsl_path(win_path: Path) -> str:
    win_path = win_path.resolve()
    drive = win_path.drive.rstrip(":").lower()
    rest = str(win_path)[len(win_path.drive):].replace("\\", "/").lstrip("/")
    return f"/mnt/{drive}/{rest}"


# --------------------------------------------------------------------------
# App/bundled-resource paths
#
# IMPORTANT: this app is built with Nuitka's --onefile mode. In that mode,
# __file__ points at a TEMPORARY extraction folder that changes every launch
# and disappears on exit -- it is NOT where GromacsWizard.exe actually lives
# on disk. Nuitka's own docs are explicit about this: sys.argv[0] is the
# real/original executable path; __file__ is the temp bootstrap path. So we
# must use sys.argv[0] (via the "__compiled__" check Nuitka injects) to find
# resources installed next to the real .exe, and only fall back to __file__
# when running as a plain, uncompiled .py script.
# --------------------------------------------------------------------------

def get_app_dir() -> Path:
    if "__compiled__" in globals():
        return Path(sys.argv[0]).resolve().parent
    return Path(__file__).resolve().parent


APP_DIR = get_app_dir()
PIPELINE_SCRIPT = APP_DIR / "pipeline" / "run_pipeline.py"
FORCEFIELD_PARENT_DIR = APP_DIR / "pipeline" / "forcefield"


def bundled_ff_name() -> str:
    """Name of the shared force field folder shipped with this install, for
    display purposes only. Falls back gracefully if running unbuilt/uninstalled."""
    try:
        candidates = [p for p in FORCEFIELD_PARENT_DIR.glob("*.ff") if p.is_dir()]
        if candidates:
            return candidates[0].name
    except OSError:
        pass
    return "(not found -- reinstall the app)"


# --------------------------------------------------------------------------
# Shared wizard state
# --------------------------------------------------------------------------

class WizardState:
    def __init__(self):
        self.workdir: Path | None = load_last_workdir()
        self.rec = None
        self.lig_pdb = None
        self.lig_itp = None
        self.lig_prm = None
        self.run_mode = "full"       # "full" or "single"
        self.single_step = STEPS[0]
        self.force = False
        self.prod_ns = "100"
        self.gpu_id = "0"
        self.nt = "16"
        self.maxh = ""
        self.pin = True


# --------------------------------------------------------------------------
# Base page
# --------------------------------------------------------------------------

class WizardPage(ctk.CTkFrame):
    """Base class for a wizard screen. Subclasses set self.can_go_next()
    and build their widgets in build()."""

    title_text = ""
    subtitle_text = ""

    def __init__(self, master, app):
        super().__init__(master, fg_color="transparent")
        self.app = app
        self.state: WizardState = app.wiz_state

        header = ctk.CTkFrame(self, fg_color="transparent")
        header.pack(fill="x", pady=(0, 14))
        ctk.CTkLabel(header, text=self.title_text, font=FONT_TITLE, anchor="w").pack(fill="x")
        if self.subtitle_text:
            ctk.CTkLabel(header, text=self.subtitle_text, font=FONT_SUBTITLE,
                         anchor="w", text_color="gray70", wraplength=740, justify="left").pack(fill="x", pady=(2, 0))

        self.body = ctk.CTkFrame(self, fg_color="transparent")
        self.body.pack(fill="both", expand=True)
        self.build()

    def build(self):
        pass

    def on_show(self):
        """Called every time the page becomes visible."""
        pass

    def can_go_next(self) -> bool:
        return True


# --------------------------------------------------------------------------
# Page 1: Welcome
# --------------------------------------------------------------------------

class WelcomePage(WizardPage):
    title_text = "GROMACS Pipeline Wizard"
    subtitle_text = ("This wizard runs your protein-ligand MD pipeline (prep through analysis) "
                      "inside WSL, using the run_pipeline.py + force field installed alongside "
                      "this application. Click through each step -- nothing runs until you "
                      "confirm on the review screen.")

    def build(self):
        card = ctk.CTkFrame(self.body, corner_radius=14)
        card.pack(fill="x", pady=10)
        lines = [
            "1.  Check your system (WSL + GROMACS)",
            "2.  Pick the folder with your run's input files",
            "3.  Choose which steps to run",
            "4.  Set simulation options (length, GPU, threads)",
            "5.  Review everything and confirm",
            "6.  Watch it run, with live progress and a log",
        ]
        for line in lines:
            ctk.CTkLabel(card, text=line, font=FONT_BODY, anchor="w").pack(fill="x", padx=20, pady=6)
        ctk.CTkLabel(self.body, text=CREATOR_CREDIT, font=FONT_SUBTITLE,
                     text_color="gray60").pack(anchor="e", pady=(8, 0))


# --------------------------------------------------------------------------
# Page 2: System check
# --------------------------------------------------------------------------

class CheckPage(WizardPage):
    title_text = "Step 1 -- System Check"
    subtitle_text = "Making sure WSL, GROMACS and Python are ready inside your Linux environment."

    # Each failing check gets a one-click "Install" button that opens a real
    # terminal window and runs the actual install command. We can't make
    # these fully silent -- `sudo` needs a password and installing WSL
    # itself needs a Windows admin prompt -- but the user never has to type
    # or copy a command themselves.
    INSTALL_SPECS = {
        "wsl": "Install WSL",
        "distro": "Install Ubuntu",
        "gmx": "Install GROMACS",
        "py": "Install Python deps",
    }

    def build(self):
        self.rows = {}
        self.install_buttons = {}
        for key, label in [("wsl", "WSL installed"), ("distro", "A Linux distro is set up"),
                            ("gmx", "GROMACS (gmx) available in WSL"),
                            ("py", "python3 + tqdm available in WSL")]:
            row = ctk.CTkFrame(self.body, corner_radius=10)
            row.pack(fill="x", pady=5)
            dot = ctk.CTkLabel(row, text="\u25CF", text_color="gray50", font=("Segoe UI", 16))
            dot.pack(side="left", padx=(14, 8), pady=10)
            ctk.CTkLabel(row, text=label, font=FONT_BODY, anchor="w").pack(side="left", pady=10)
            msg = ctk.CTkLabel(row, text="", font=("Segoe UI", 11), text_color="gray60", anchor="e")
            msg.pack(side="right", padx=14, pady=10)
            btn = ctk.CTkButton(row, text=self.INSTALL_SPECS[key], width=140, fg_color="#8e44ad",
                                hover_color="#732d91", command=lambda k=key: self.run_install(k))
            # not packed yet -- only shown once a check comes back failed
            self.install_buttons[key] = btn
            self.rows[key] = (dot, msg)

        self.recheck_btn = ctk.CTkButton(self.body, text="Run Checks Again", command=self.start_checks)
        self.recheck_btn.pack(pady=16)

        self.help_box = ctk.CTkTextbox(self.body, height=140, font=FONT_MONO, wrap="word")
        self.help_box.pack(fill="both", expand=True, pady=(4, 0))
        self.help_box.insert("1.0", "Install / troubleshooting tips will appear here if a check fails.")
        self.help_box.configure(state="disabled")

        self._passed = False

    # -- one-click installers ------------------------------------------------
    # Each opens a visible terminal window running the real command, so any
    # sudo password prompt or Windows admin (UAC) prompt shows up there.

    def _open_console(self, cmd_list, note):
        try:
            flags = getattr(subprocess, "CREATE_NEW_CONSOLE", 0)
            subprocess.Popen(cmd_list, creationflags=flags)
            self.append_help(note)
        except Exception as e:
            self.append_help(f"Could not open an install window automatically ({e}). "
                              "You'll need to run the command shown above manually instead.")

    def run_install(self, key):
        if key == "wsl":
            if not messagebox.askyesno(
                "Install WSL",
                "This opens an elevated PowerShell and runs 'wsl --install'. Windows will ask you to "
                "approve admin access, and your PC will likely need a restart afterward. Continue?"):
                return
            self._open_console(
                ["powershell", "-Command",
                 "Start-Process powershell -ArgumentList '-NoExit -Command wsl --install' -Verb RunAs"],
                "Opened an elevated PowerShell to install WSL. Approve the admin prompt, let it finish, "
                "restart your PC if it asks you to, then come back and click 'Run Checks Again'.")
        elif key == "distro":
            self._open_console(
                ["wsl", "--install", "-d", "Ubuntu"],
                "Opened a window installing Ubuntu under WSL. It'll ask you to create a username/password "
                "the first time -- follow the prompts, then click 'Run Checks Again'.")
        elif key == "gmx":
            bash_cmd = ("sudo apt update && sudo apt install -y gromacs; echo; "
                        "echo Done -- close this window and click Run Checks Again in the wizard.; "
                        "read -n 1 -s")
            self._open_console(
                ["wsl", "bash", "-lc", bash_cmd],
                "Opened a WSL terminal installing GROMACS. Type your WSL password if it's asked for "
                "(it won't show on screen -- that's normal), let it finish, then click 'Run Checks Again'.")
        elif key == "py":
            bash_cmd = ("sudo apt install -y python3 python3-pip && "
                        "(pip3 install tqdm --break-system-packages || pip3 install tqdm); "
                        "echo; echo Done -- close this window and click Run Checks Again in the wizard.; "
                        "read -n 1 -s")
            self._open_console(
                ["wsl", "bash", "-lc", bash_cmd],
                "Opened a WSL terminal installing python3/pip/tqdm. Type your WSL password if asked, "
                "let it finish, then click 'Run Checks Again'.")

    def on_show(self):
        if not getattr(self, "_ran_once", False):
            self._ran_once = True
            self.start_checks()

    def set_row(self, key, ok, msg=""):
        dot, msg_label = self.rows[key]
        dot.configure(text_color="#2ecc71" if ok else "#e74c3c")
        msg_label.configure(text=msg)
        btn = self.install_buttons.get(key)
        if btn is not None:
            if ok is False:
                btn.pack(side="right", padx=(4, 0))
            else:
                btn.pack_forget()

    def append_help(self, text):
        self.help_box.configure(state="normal")
        self.help_box.insert("end", text + "\n\n")
        self.help_box.configure(state="disabled")

    def start_checks(self):
        self.recheck_btn.configure(state="disabled", text="Checking...")
        self.help_box.configure(state="normal")
        self.help_box.delete("1.0", "end")
        self.help_box.configure(state="disabled")
        for key in self.rows:
            self.set_row(key, None, "checking...")
            self.rows[key][0].configure(text_color="gray50")
        threading.Thread(target=self._run_checks, daemon=True).start()

    def _run_checks(self):
        all_ok = True

        wsl_ok = shutil.which("wsl") is not None
        self.after_set(self.set_row, "wsl", wsl_ok, "found" if wsl_ok else "not found")
        if not wsl_ok:
            self.after_append(
                "WSL isn't installed. In an Administrator PowerShell run:\n"
                "  wsl --install\n"
                "then restart your PC and open 'Ubuntu' once from the Start menu to finish setup."
            )
            self.after_set(self.set_row, "distro", False, "skipped")
            self.after_set(self.set_row, "gmx", False, "skipped")
            self.after_set(self.set_row, "py", False, "skipped")
            self.after_finish(False)
            return

        try:
            result = subprocess.run(["wsl", "-l", "-q"], text=True, capture_output=True, timeout=15)
            distros = [d.strip() for d in result.stdout.splitlines() if d.strip()]
        except Exception:
            distros = []
        distro_ok = len(distros) > 0
        self.after_set(self.set_row, "distro", distro_ok, distros[0] if distro_ok else "none found")
        if not distro_ok:
            all_ok = False
            self.after_append("No Linux distro is installed under WSL. In PowerShell run:\n  wsl --install -d Ubuntu\n"
                               "then open Ubuntu once from the Start menu to finish setup.")

        # Try the plain command first, then fall back to sourcing GMXRC from
        # common install locations -- catches source-built (e.g. CUDA)
        # GROMACS installs whose GMXRC is only sourced from ~/.bashrc, which
        # non-interactive login shells (like this one) skip.
        gmx_probe = (
            "gmx --version >/dev/null 2>&1 && gmx --version && exit 0; "
            "for f in /usr/local/gromacs*/bin/GMXRC /opt/gromacs*/bin/GMXRC "
            "$HOME/gromacs*/bin/GMXRC /usr/local/bin/GMXRC; do "
            "[ -f \"$f\" ] && . \"$f\" && gmx --version >/dev/null 2>&1 && gmx --version && exit 0; "
            "done; exit 1"
        )
        try:
            gmx = wsl_run(gmx_probe, timeout=20)
            gmx_ok = gmx.returncode == 0
        except Exception:
            gmx_ok = False
        version_line = ""
        if gmx_ok and gmx.stdout.strip():
            version_line = gmx.stdout.strip().splitlines()[0]
        self.after_set(self.set_row, "gmx", gmx_ok, version_line if gmx_ok else "not found")
        if not gmx_ok:
            all_ok = False
            self.after_append(
                "GROMACS isn't available inside WSL (checked PATH and common source-build locations).\n\n"
                "Already built GROMACS yourself (e.g. a CUDA build) and it works in your normal terminal? "
                "Don't use the Install button below -- that installs Ubuntu's CPU-only apt package, which "
                "won't have GPU support and may shadow your own build. Instead this is almost always a PATH "
                "issue: your GMXRC source line is probably only in ~/.bashrc, which non-interactive shells "
                "(like this check, and the pipeline runs) skip. Fix it by adding the same 'source "
                ".../GMXRC' line to ~/.profile instead (or in addition).\n\n"
                "Don't have GROMACS at all yet? Inside your WSL Ubuntu terminal run:\n"
                "  sudo apt update && sudo apt install -y gromacs\n"
                "For GPU (CUDA) support you'll need to build GROMACS from source with CUDA enabled -- "
                "see https://manual.gromacs.org/current/install-guide/index.html"
            )

        try:
            py = wsl_run(
                'python3 -c "import tqdm" 2>/dev/null && echo OK || '
                '(pip3 install --user tqdm --break-system-packages >/dev/null 2>&1 || '
                'pip3 install --user tqdm >/dev/null 2>&1 && '
                'python3 -c "import tqdm" && echo OK)',
                timeout=60)
            py_ok = py.returncode == 0 and "OK" in py.stdout
        except Exception:
            py_ok = False
        self.after_set(self.set_row, "py", py_ok, "ready" if py_ok else "missing tqdm")
        if not py_ok:
            all_ok = False
            self.after_append(
                "python3 or the 'tqdm' package isn't ready inside WSL. Run inside WSL:\n"
                "  sudo apt install -y python3 python3-pip\n"
                "  pip3 install tqdm --break-system-packages"
            )

        self.after_finish(all_ok and distro_ok and gmx_ok and py_ok)

    # -- thread-safe UI helpers --------------------------------------------
    def after_set(self, fn, *args):
        self.after(0, lambda: fn(*args))

    def after_append(self, text):
        self.after(0, lambda: self.append_help(text))

    def after_finish(self, ok):
        def _finish():
            self._passed = ok
            self.recheck_btn.configure(state="normal", text="Run Checks Again")
            self.app.refresh_nav()
        self.after(0, _finish)

    def can_go_next(self):
        if not self._passed:
            return messagebox.askyesno(
                "Checks not all green",
                "Not everything passed. You can still continue, but the pipeline will likely fail. Continue anyway?"
            )
        return True


# --------------------------------------------------------------------------
# Page 3: Folder + input files
#
# NOTE: the *.ff force-field folder is no longer scanned for here -- it now
# ships bundled with the application install (see FORCEFIELD_PARENT_DIR
# above) and is shared by every project. This page only needs the actual
# per-run inputs.
# --------------------------------------------------------------------------

class FolderPage(WizardPage):
    title_text = "Step 2 -- Choose Your Run Folder"
    subtitle_text = ("Pick the folder that has REC.pdb, LIG.pdb, and LIG.itp (the SwissParam-generated "
                      "ligand topology). The wizard will scan it automatically.")

    def build(self):
        row = ctk.CTkFrame(self.body, fg_color="transparent")
        row.pack(fill="x", pady=(0, 10))
        self.folder_entry = ctk.CTkEntry(row, placeholder_text="No folder selected", font=FONT_BODY)
        self.folder_entry.pack(side="left", fill="x", expand=True, padx=(0, 8))
        ctk.CTkButton(row, text="Browse...", width=110, command=self.browse).pack(side="left")
        ctk.CTkButton(row, text="Scan Folder", width=110, command=self.scan).pack(side="left", padx=(8, 0))

        self.status_frame = ctk.CTkFrame(self.body, corner_radius=14)
        self.status_frame.pack(fill="both", expand=True, pady=10)
        self.status_rows = {}
        for key, label in [("rec", "REC.pdb (protein)"), ("lig_pdb", "LIG.pdb (ligand coords)"),
                            ("lig_itp", "LIG.itp (ligand topology)"), ("lig_prm", "LIG.prm (optional)")]:
            r = ctk.CTkFrame(self.status_frame, fg_color="transparent")
            r.pack(fill="x", padx=16, pady=6)
            dot = ctk.CTkLabel(r, text="\u25CF", text_color="gray50", font=("Segoe UI", 16))
            dot.pack(side="left", padx=(0, 8))
            ctk.CTkLabel(r, text=label, font=FONT_BODY, width=220, anchor="w").pack(side="left")
            val = ctk.CTkLabel(r, text="not scanned yet", font=FONT_BODY, text_color="gray60", anchor="w")
            val.pack(side="left", fill="x", expand=True)
            self.status_rows[key] = (dot, val)

        ff_row = ctk.CTkFrame(self.status_frame, fg_color="transparent")
        ff_row.pack(fill="x", padx=16, pady=(10, 6))
        ctk.CTkLabel(ff_row, text="Force field (bundled with the app):", font=FONT_BODY,
                     width=260, anchor="w", text_color="gray60").pack(side="left")
        ctk.CTkLabel(ff_row, text=bundled_ff_name(), font=FONT_BODY, text_color="gray60").pack(side="left")

    def on_show(self):
        # Reuse the saved folder after restarting the app or finishing a step.
        if self.state.workdir and self.state.workdir.is_dir():
            current = self.folder_entry.get().strip()
            if current != str(self.state.workdir):
                self.folder_entry.delete(0, "end")
                self.folder_entry.insert(0, str(self.state.workdir))
            self.scan()

    def browse(self):
        path = filedialog.askdirectory(title="Select your run folder")
        if path:
            self.folder_entry.delete(0, "end")
            self.folder_entry.insert(0, path)
            self.scan()

    def scan(self):
        path_str = self.folder_entry.get().strip()
        if not path_str:
            messagebox.showwarning("No folder", "Choose a folder first.")
            return
        workdir = Path(path_str)
        if not workdir.is_dir():
            messagebox.showerror("Not found", f"'{workdir}' is not a folder.")
            return
        self.state.workdir = workdir
        save_last_workdir(workdir)

        found = {
            "rec": workdir / "REC.pdb" if (workdir / "REC.pdb").exists() else None,
            "lig_pdb": workdir / "LIG.pdb" if (workdir / "LIG.pdb").exists() else None,
            "lig_itp": workdir / "LIG.itp" if (workdir / "LIG.itp").exists() else None,
            "lig_prm": workdir / "LIG.prm" if (workdir / "LIG.prm").exists() else None,
        }

        for key in ("rec", "lig_pdb", "lig_itp"):
            ok = found[key] is not None
            dot, val = self.status_rows[key]
            dot.configure(text_color="#2ecc71" if ok else "#e74c3c")
            val.configure(text=found[key].name if ok else "MISSING -- required")
            setattr(self.state, key, found[key])

        dot, val = self.status_rows["lig_prm"]
        if found["lig_prm"]:
            dot.configure(text_color="#2ecc71")
            val.configure(text=found["lig_prm"].name)
        else:
            dot.configure(text_color="gray50")
            val.configure(text="not present (optional)")
        self.state.lig_prm = found["lig_prm"]

        self.app.refresh_nav()

    def can_go_next(self):
        # Also support a folder path typed or pasted into the input field.
        if self.folder_entry.get().strip() != str(self.state.workdir or ""):
            self.scan()
        missing = [k for k in ("rec", "lig_pdb", "lig_itp") if getattr(self.state, k) is None]
        if missing:
            messagebox.showerror("Missing files", f"Still missing: {', '.join(missing)}. "
                                                    "Add them to the folder and click 'Scan Folder' again.")
            return False
        return True


# --------------------------------------------------------------------------
# Page 4: Step selection
# --------------------------------------------------------------------------

class StepSelectPage(WizardPage):
    title_text = "Step 3 -- What Do You Want to Run?"
    subtitle_text = "Run the whole pipeline start to finish, or just resume/re-run a single step."

    def build(self):
        self.mode_var = ctk.StringVar(value="full")
        self.segment = ctk.CTkSegmentedButton(
            self.body, values=["Full pipeline", "Single step"],
            command=self._on_mode_change)
        self.segment.set("Full pipeline")
        self.segment.pack(fill="x", pady=(0, 16))

        self.step_menu_frame = ctk.CTkFrame(self.body, fg_color="transparent")
        ctk.CTkLabel(self.step_menu_frame, text="Which step?", font=FONT_BODY).pack(side="left", padx=(0, 10))
        self.step_menu = ctk.CTkOptionMenu(
            self.step_menu_frame, values=[f"{s} -- {STEP_LABELS[s]}" for s in STEPS],
            command=lambda v: setattr(self.state, "single_step", v.split(" -- ")[0]))
        self.step_menu.pack(side="left")

        self.force_var = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(self.body, text="Force rerun (ignore existing checkpoints, redo from scratch)",
                        variable=self.force_var,
                        command=lambda: setattr(self.state, "force", self.force_var.get())).pack(anchor="w", pady=20)

        info = ctk.CTkLabel(self.body, font=FONT_SUBTITLE, text_color="gray60", justify="left", wraplength=740,
                            text=("Every step is checkpointed -- if a run crashes or you close the wizard, "
                                  "running 'Full pipeline' again picks up where it left off, without redoing "
                                  "finished steps (unless Force rerun is checked)."))
        info.pack(fill="x")

    def on_show(self):
        """Synchronize reusable controls with state after a finished run."""
        if self.state.run_mode == "single":
            self.segment.set("Single step")
            self.step_menu.set(f"{self.state.single_step} -- {STEP_LABELS[self.state.single_step]}")
            self.step_menu_frame.pack(fill="x", pady=(0, 10))
        else:
            self.segment.set("Full pipeline")
            self.step_menu_frame.pack_forget()

    def _on_mode_change(self, choice):
        if choice == "Full pipeline":
            self.state.run_mode = "full"
            self.step_menu_frame.pack_forget()
        else:
            self.state.run_mode = "single"
            self.step_menu_frame.pack(fill="x", pady=(0, 10))


# --------------------------------------------------------------------------
# Page 5: Run options
# --------------------------------------------------------------------------

class RunOptionsPage(WizardPage):
    title_text = "Step 4 -- Simulation Options"
    subtitle_text = "These control the production MD run and GPU usage."

    def build(self):
        grid = ctk.CTkFrame(self.body, corner_radius=14)
        grid.pack(fill="x", pady=6)

        self.entries = {}
        fields = [
            ("prod_ns", "Production run length (ns)", "100"),
            ("gpu_id", "GPU id", "0"),
            ("nt", "Total thread count for mdrun", "16"),
            ("maxh", "Max wall-hours for MD (blank = no limit)", ""),
        ]
        for i, (key, label, placeholder) in enumerate(fields):
            ctk.CTkLabel(grid, text=label, font=FONT_BODY, width=300, anchor="w").grid(
                row=i, column=0, padx=16, pady=10, sticky="w")
            entry = ctk.CTkEntry(grid, placeholder_text=placeholder, width=180)
            entry.insert(0, getattr(self.state, key))
            entry.grid(row=i, column=1, padx=16, pady=10, sticky="w")
            self.entries[key] = entry

        self.pin_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(self.body, text="Pin threads to cores (-pin on) -- recommended",
                        variable=self.pin_var).pack(anchor="w", pady=16)

    def can_go_next(self):
        for key in ("prod_ns", "gpu_id", "nt"):
            val = self.entries[key].get().strip()
            if not val:
                messagebox.showerror("Missing value", f"Please fill in '{key}'.")
                return False
        try:
            float(self.entries["prod_ns"].get().strip())
            int(self.entries["nt"].get().strip())
            if self.entries["maxh"].get().strip():
                float(self.entries["maxh"].get().strip())
        except ValueError:
            messagebox.showerror("Invalid value", "Production ns, thread count, and max hours must be numbers.")
            return False
        self.state.prod_ns = self.entries["prod_ns"].get().strip()
        self.state.gpu_id = self.entries["gpu_id"].get().strip()
        self.state.nt = self.entries["nt"].get().strip()
        self.state.maxh = self.entries["maxh"].get().strip()
        self.state.pin = self.pin_var.get()
        return True


# --------------------------------------------------------------------------
# Page 6: Review
# --------------------------------------------------------------------------

class ReviewPage(WizardPage):
    title_text = "Step 5 -- Review & Confirm"
    subtitle_text = "Check everything below, then click Next to start the run."

    def build(self):
        self.box = ctk.CTkTextbox(self.body, font=FONT_MONO, wrap="word")
        self.box.pack(fill="both", expand=True)
        self.box.configure(state="disabled")

    def on_show(self):
        s = self.state
        run_desc = "Full pipeline" if s.run_mode == "full" else f"Single step: {s.single_step}"
        lines = [
            f"Working folder : {s.workdir}",
            f"REC.pdb        : {s.rec.name if s.rec else '-'}",
            f"LIG.pdb        : {s.lig_pdb.name if s.lig_pdb else '-'}",
            f"LIG.itp        : {s.lig_itp.name if s.lig_itp else '-'}",
            f"LIG.prm        : {s.lig_prm.name if s.lig_prm else '(none)'}",
            f"Force field    : {bundled_ff_name()} (bundled)",
            "",
            f"Run mode       : {run_desc}",
            f"Force rerun    : {'yes' if s.force else 'no'}",
            "",
            f"Production ns  : {s.prod_ns}",
            f"GPU id         : {s.gpu_id}",
            f"Threads (-nt)  : {s.nt}",
            f"Max wall-hours : {s.maxh or '(no limit)'}",
            f"Pin threads    : {'yes' if s.pin else 'no'}",
        ]
        self.box.configure(state="normal")
        self.box.delete("1.0", "end")
        self.box.insert("1.0", "\n".join(lines))
        self.box.configure(state="disabled")


# --------------------------------------------------------------------------
# Page 7: Run
# --------------------------------------------------------------------------

class RunPage(WizardPage):
    title_text = "Running Pipeline"
    subtitle_text = "Live output is streamed below. You can leave this running in the background."

    def build(self):
        self.step_rows = {}
        step_list = ctk.CTkFrame(self.body, corner_radius=14)
        step_list.pack(fill="x", pady=(0, 10))
        for step in STEPS:
            row = ctk.CTkFrame(step_list, fg_color="transparent")
            row.pack(fill="x", padx=14, pady=4)
            dot = ctk.CTkLabel(row, text="\u25CB", width=20, text_color="gray50", font=("Segoe UI", 14))
            dot.pack(side="left")
            ctk.CTkLabel(row, text=STEP_LABELS[step], font=FONT_BODY, width=260, anchor="w").pack(side="left")
            bar = ctk.CTkProgressBar(row, width=300)
            bar.set(0)
            bar.pack(side="left", padx=10)
            self.step_rows[step] = (dot, bar)

        self.log_box = ctk.CTkTextbox(self.body, font=FONT_MONO, wrap="none")
        self.log_box.pack(fill="both", expand=True, pady=(6, 6))
        self.log_box.configure(state="disabled")

        btn_row = ctk.CTkFrame(self.body, fg_color="transparent")
        btn_row.pack(fill="x")
        self.status_label = ctk.CTkLabel(btn_row, text="Starting...", font=FONT_BODY)
        self.status_label.pack(side="left")
        self.cancel_btn = ctk.CTkButton(btn_row, text="Cancel", fg_color="#c0392b",
                                        hover_color="#922b21", command=self.cancel, width=100)
        self.cancel_btn.pack(side="right")

        self._process = None
        self._log_queue = queue.Queue()
        self._current_step = None
        self._started = False
        self._finished = False

    def on_show(self):
        if not self._started:
            self._started = True
            self.start()

    def reset_for_new_run(self):
        """Make this reusable page ready for another pipeline step."""
        self._started = False
        self._finished = False
        self._process = None
        self._current_step = None
        self._log_queue = queue.Queue()
        self._last_was_progress = False
        for dot, bar in self.step_rows.values():
            dot.configure(text="\u25CB", text_color="gray50")
            bar.set(0)
        self.log_box.configure(state="normal")
        self.log_box.delete("1.0", "end")
        self.log_box.configure(state="disabled")
        self.status_label.configure(text="Starting...")
        self.cancel_btn.configure(text="Cancel", fg_color="#c0392b", hover_color="#922b21", command=self.cancel)

    def append_log(self, text):
        self.log_box.configure(state="normal")
        self.log_box.insert("end", text)
        self.log_box.see("end")
        self.log_box.configure(state="disabled")

    def append_log_smart(self, text):
        is_progress = bool(PERCENT_RE.search(text))
        self.log_box.configure(state="normal")
        if is_progress and getattr(self, "_last_was_progress", False):
            self.log_box.delete("end-2l", "end-1l")
        self.log_box.insert("end", text)
        self.log_box.see("end")
        self.log_box.configure(state="disabled")
        self._last_was_progress = is_progress

    def mark_step(self, step, status):
        """status: 'running' | 'done' | 'failed'"""
        if step not in self.step_rows:
            return
        dot, bar = self.step_rows[step]
        if status == "running":
            dot.configure(text="\u25CF", text_color="#3498db")
        elif status == "done":
            dot.configure(text="\u25CF", text_color="#2ecc71")
            bar.set(1)
        elif status == "failed":
            dot.configure(text="\u25CF", text_color="#e74c3c")

    def build_command(self):
        s = self.state
        wsl_workdir = windows_to_wsl_path(s.workdir)
        wsl_pipeline_script = windows_to_wsl_path(PIPELINE_SCRIPT)
        wsl_ff_parent_dir = windows_to_wsl_path(FORCEFIELD_PARENT_DIR)

        if not PIPELINE_SCRIPT.is_file():
            raise FileNotFoundError(
                f"run_pipeline.py not found at {PIPELINE_SCRIPT}. The installation looks broken -- "
                "try reinstalling GROMACS Wizard."
            )

        # Make sure gmx is on PATH for this non-interactive shell even if the
        # user's GMXRC source line only lives in ~/.bashrc (common for
        # source-built/CUDA installs, which non-interactive shells skip).
        gmx_bootstrap = (
            "command -v gmx >/dev/null 2>&1 || "
            "for f in /usr/local/gromacs*/bin/GMXRC /opt/gromacs*/bin/GMXRC "
            "$HOME/gromacs*/bin/GMXRC /usr/local/bin/GMXRC; do "
            "[ -f \"$f\" ] && . \"$f\" && break; done"
        )
        parts = [
            f"cd '{wsl_workdir}' &&",
            f"{gmx_bootstrap};",
            f"python3 -u '{wsl_pipeline_script}'",
            f"--rec '{s.rec.name}'",
            f"--lig-pdb '{s.lig_pdb.name}'",
            f"--lig-itp '{s.lig_itp.name}'",
            f"--ff-source '{wsl_ff_parent_dir}'",
            f"--prod-ns {s.prod_ns}",
            f"--gpu-id {s.gpu_id}",
            f"--nt {s.nt}",
        ]
        if s.lig_prm:
            parts.append(f"--lig-prm '{windows_to_wsl_path(s.lig_prm)}'")
        if s.maxh:
            parts.append(f"--maxh {s.maxh}")
        if s.run_mode == "single":
            parts.append(f"--step {s.single_step}")
        if s.force:
            parts.append("--force")
        return " ".join(parts)

    def start(self):
        try:
            bash_line = self.build_command()
        except Exception as e:
            self.append_log(f"[error] Could not build command: {e}\n")
            self.status_label.configure(text="Failed to start")
            return
        self.append_log(f"$ {bash_line}\n\n")
        threading.Thread(target=self._run, args=(bash_line,), daemon=True).start()
        self.after(80, self._poll_queue)

    def _run(self, bash_line):
        try:
            self._process = subprocess.Popen(
                ["wsl", "bash", "-lc", bash_line], text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1,
                # WSL/GROMACS may emit UTF-8 progress glyphs that Windows'
                # default cp1252 decoder cannot read. Never let one byte kill
                # the worker thread that reports pipeline completion.
                encoding="utf-8", errors="replace")
        except Exception as e:
            self._log_queue.put(("error", f"Could not launch WSL: {e}\n"))
            self._log_queue.put(("done", 1))
            return

        for line in self._process.stdout:
            self._log_queue.put(("line", line))
        self._process.wait()
        self._log_queue.put(("done", self._process.returncode))

    def _poll_queue(self):
        try:
            while True:
                kind, payload = self._log_queue.get_nowait()
                if kind == "line":
                    self.append_log_smart(payload if payload.endswith("\n") else payload + "\n")
                    # run_pipeline.py prints this flushed marker after every
                    # step command has returned.  Do not wait for WSL's pipe
                    # cleanup before taking the user to the Done screen.
                    if "PIPELINE_GUI_STATUS: COMPLETE" in payload:
                        self._finish(0)
                        return
                    if "PIPELINE_GUI_STATUS: FAILED" in payload:
                        self._finish(1)
                        return
                    m = STEP_HEADER_RE.match(payload.strip())
                    if m:
                        if self._current_step:
                            self.mark_step(self._current_step, "done")
                        self._current_step = m.group(1)
                        self.mark_step(self._current_step, "running")
                        self.status_label.configure(text=f"Running: {STEP_LABELS.get(self._current_step, self._current_step)}")
                    pm = PERCENT_RE.search(payload)
                    if pm and self._current_step:
                        _, bar = self.step_rows.get(self._current_step, (None, None))
                        if bar is not None:
                            bar.set(min(100, int(pm.group(1))) / 100)
                elif kind == "error":
                    self.append_log(f"[error] {payload}")
                elif kind == "done":
                    self._finish(payload)
                    return
        except queue.Empty:
            pass
        self.after(80, self._poll_queue)

    def _finish(self, returncode):
        # A status marker can arrive just before the WSL subprocess-reaped
        # event.  Handle the first one only.
        if self._finished:
            return
        self._finished = True
        if returncode == 0:
            if self._current_step:
                self.mark_step(self._current_step, "done")
            self.status_label.configure(text="Finished successfully")
        else:
            if self._current_step:
                self.mark_step(self._current_step, "failed")
            self.status_label.configure(text=f"Failed (exit code {returncode})")
        self.cancel_btn.configure(text="Close", fg_color="#555555", hover_color="#444444", command=self.close_only)
        self.app.on_run_finished(returncode == 0)

    def cancel(self):
        if self._process and self._process.poll() is None:
            if not messagebox.askyesno("Cancel run", "Stop the pipeline now? Progress already made is checkpointed."):
                return
            try:
                subprocess.run(["wsl", "bash", "-lc", "pkill -f run_pipeline.py; pkill -f 'gmx mdrun'"])
            except Exception:
                pass
            self.status_label.configure(text="Cancelling...")

    def close_only(self):
        idx = self.app.pages.index(next(p for p in self.app.pages if isinstance(p, FinishedPage)))
        self.app.show_page(idx)

# --------------------------------------------------------------------------
# Page 8: Finished
# --------------------------------------------------------------------------

class FinishedPage(WizardPage):
    title_text = "Done"
    subtitle_text = ""

    def build(self):
        self.msg = ctk.CTkLabel(self.body, text="", font=("Segoe UI", 16))
        self.msg.pack(pady=30)
        btn_row = ctk.CTkFrame(self.body, fg_color="transparent")
        btn_row.pack()
        self.next_step_btn = ctk.CTkButton(btn_row, text="Run Next Step", command=self.app.run_next_step)
        self.next_step_btn.pack(side="left", padx=8)
        ctk.CTkButton(btn_row, text="Change Step / Options", command=self.app.resume_at_step_selection).pack(side="left", padx=8)
        ctk.CTkButton(btn_row, text="Open Run Folder", command=self.open_folder).pack(side="left", padx=8)
        self.plots_btn = ctk.CTkButton(btn_row, text="Open Publication Plots",
                                       command=self.open_publication_plots)
        self.plots_btn.pack(side="left", padx=8)
        ctk.CTkButton(btn_row, text="Start New Run", command=self.app.restart).pack(side="left", padx=8)
        ctk.CTkButton(btn_row, text="Exit", fg_color="#555555", hover_color="#444444",
                     command=self.app.destroy).pack(side="left", padx=8)

    def on_show(self):
        ok = self.app.last_run_ok
        self.msg.configure(
            text=("\u2705  Pipeline completed successfully." if ok else
                  "\u274C  Pipeline stopped with an error -- check the log on the previous screen or pipeline.log."))
        can_continue = ok and self.state.run_mode == "single" and self.state.single_step != STEPS[-1]
        self.next_step_btn.configure(state="normal" if can_continue else "disabled")
        plot_dir = self.state.workdir / "analysis" / "Publication plots" if self.state.workdir else None
        self.plots_btn.configure(state="normal" if plot_dir and plot_dir.is_dir() else "disabled")

    def open_folder(self):
        if self.state.workdir:
            os.startfile(str(self.state.workdir))  # Windows only, matches target platform

    def open_publication_plots(self):
        if not self.state.workdir:
            return
        plot_dir = self.state.workdir / "analysis" / "Publication plots"
        if plot_dir.is_dir():
            os.startfile(str(plot_dir))
        else:
            messagebox.showinfo("Publication plots", "Run the Analysis step first to create publication plots.")


# --------------------------------------------------------------------------
# Main app: page container + nav bar
# --------------------------------------------------------------------------

class WizardApp(ctk.CTk):
    PAGES = [WelcomePage, CheckPage, FolderPage, StepSelectPage, RunOptionsPage, ReviewPage, RunPage, FinishedPage]

    def __init__(self):
        super().__init__()
        self.title(f"GROMACS Pipeline Wizard — {CREATOR_CREDIT}")
        self.geometry("880x680")
        self.minsize(760, 600)

        self.wiz_state = WizardState()
        self.last_run_ok = False
        self.index = 0

        self.container = ctk.CTkFrame(self, fg_color="transparent")
        self.container.pack(fill="both", expand=True, padx=24, pady=(20, 0))

        nav = ctk.CTkFrame(self, fg_color="transparent")
        nav.pack(fill="x", padx=24, pady=16)
        self.back_btn = ctk.CTkButton(nav, text="< Back", width=100, command=self.go_back)
        self.back_btn.pack(side="left")
        self.next_btn = ctk.CTkButton(nav, text="Next >", width=100, command=self.go_next)
        self.next_btn.pack(side="right")
        self.progress_label = ctk.CTkLabel(nav, text="", text_color="gray60")
        self.progress_label.pack(side="right", padx=16)

        self.pages = []
        for cls in self.PAGES:
            page = cls(self.container, self)
            page.place(relx=0, rely=0, relwidth=1, relheight=1)
            self.pages.append(page)

        self.show_page(0)

    def show_page(self, index):
        self.index = index
        page = self.pages[index]
        page.tkraise()
        page.on_show()
        is_run_page = isinstance(page, RunPage)
        is_finished_page = isinstance(page, FinishedPage)
        self.back_btn.configure(state="disabled" if (index == 0 or is_run_page or is_finished_page) else "normal")
        self.next_btn.configure(state="disabled" if (is_run_page or is_finished_page) else "normal")
        if is_finished_page:
            self.next_btn.pack_forget()
            self.back_btn.pack_forget()
        else:
            self.next_btn.pack(side="right")
            if index != 0:
                self.back_btn.pack(side="left")
        self.progress_label.configure(text=f"Step {index + 1} of {len(self.pages)}")

    def refresh_nav(self):
        pass  # hook for pages to request a nav re-evaluation; buttons re-check on click

    def go_next(self):
        page = self.pages[self.index]
        if not page.can_go_next():
            return
        if self.index + 1 < len(self.pages):
            self.show_page(self.index + 1)

    def go_back(self):
        if self.index > 0:
            self.show_page(self.index - 1)

    def on_run_finished(self, ok):
        self.last_run_ok = ok
        if ok:
            self.after(600, lambda: self.show_page(self.pages.index(next(p for p in self.pages if isinstance(p, FinishedPage)))))
        # On failure, stay on RunPage so the user can read/copy the actual
        # WSL/Python error. RunPage's Close button still opens Done manually.

    def restart(self):
        self.wiz_state = WizardState()
        for page in self.pages:
            page.state = self.wiz_state
        self.show_page(0)

    def _page_index(self, page_type):
        return self.pages.index(next(page for page in self.pages if isinstance(page, page_type)))

    def resume_at_step_selection(self):
        """Keep the selected folder and options, but change the requested step."""
        self.show_page(self._page_index(StepSelectPage))

    def run_next_step(self):
        """After a successful single step, launch its following step directly."""
        if self.wiz_state.run_mode != "single":
            return
        try:
            next_index = STEPS.index(self.wiz_state.single_step) + 1
            self.wiz_state.single_step = STEPS[next_index]
        except (ValueError, IndexError):
            return
        run_page = next(page for page in self.pages if isinstance(page, RunPage))
        run_page.reset_for_new_run()
        self.show_page(self._page_index(RunPage))


def main():
    if os.name != "nt":
        print("Note: this GUI is designed to run on Windows (it calls WSL). "
              "You can still explore the interface, but pipeline launches will fail here.")
    app = WizardApp()
    app.mainloop()


if __name__ == "__main__":
    main()
