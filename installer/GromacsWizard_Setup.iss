; ============================================================================
; GROMACS Wizard -- Inno Setup installer script
; ============================================================================
; BEFORE COMPILING THIS SCRIPT:
;   1. Build GromacsWizard.exe with Nuitka (see BUILD_STEPS.txt) and place it,
;      together with GromacsWizard.ico, in the ..\App folder next to this
;      script.
;   2. Make sure ..\App\pipeline\run_pipeline.py, ..\App\pipeline\mdp\*.mdp,
;      and ..\App\pipeline\forcefield\charmm36-feb2026_cgenff-5.0.ff\ (the
;      full 17.6 MB folder) are all in place.
;   3. Open this .iss file in Inno Setup (Tools > Compile, or press F9/Ctrl+F9).
;      Output goes to .\Output\GromacsWizard_Setup_1.0.0.exe
; ============================================================================

#define MyAppName "GROMACS Wizard"
#define MyAppVersion "1.0.0"
#define MyAppPublisher "Karan Kataria"
#define MyAppExeName "GromacsWizard.exe"
#define MyAppSourceDir "..\App"

[Setup]
AppId={{B7E2B6B0-5C2A-4C1E-9B8B-7B6D9B8F3E11}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={autopf}\GromacsWizard
DefaultGroupName=GROMACS Wizard
DisableProgramGroupPage=yes
OutputBaseFilename=GromacsWizard_Setup_{#MyAppVersion}
OutputDir=Output
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
; The bundled force field folder is ~17-18 MB; onefile GromacsWizard.exe is
; ~40-50 MB. Nothing here needs admin-only system changes, but installing to
; Program Files by default still asks for elevation -- that's expected.
PrivilegesRequired=admin
ArchitecturesInstallIn64BitMode=x64compatible
SetupIconFile={#MyAppSourceDir}\GromacsWizard.ico
UninstallDisplayIcon={app}\{#MyAppExeName}

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Additional shortcuts:"; Flags: unchecked
Name: "startmenuicon"; Description: "Create a &Start Menu shortcut"; GroupDescription: "Additional shortcuts:"; Flags: checkedonce

[Files]
; The application itself
Source: "{#MyAppSourceDir}\{#MyAppExeName}"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#MyAppSourceDir}\GromacsWizard.ico"; DestDir: "{app}"; Flags: ignoreversion

; The pipeline script (shared, single copy for every project)
Source: "{#MyAppSourceDir}\pipeline\run_pipeline.py"; DestDir: "{app}\pipeline"; Flags: ignoreversion

; The mdp templates
Source: "{#MyAppSourceDir}\pipeline\mdp\*"; DestDir: "{app}\pipeline\mdp"; Flags: ignoreversion recursesubdirs createallsubdirs

; The bundled base force field (CHARMM36/CGenFF) -- shared by every project,
; the user never has to copy this into a project folder again.
Source: "{#MyAppSourceDir}\pipeline\forcefield\*"; DestDir: "{app}\pipeline\forcefield"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\GROMACS Wizard"; Filename: "{app}\{#MyAppExeName}"; Tasks: startmenuicon
Name: "{group}\Uninstall GROMACS Wizard"; Filename: "{uninstallexe}"; Tasks: startmenuicon
Name: "{autodesktop}\GROMACS Wizard"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "Launch GROMACS Wizard now"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
; Clean up anything the app writes at runtime, if present.
Type: filesandordirs; Name: "{app}\pipeline\__pycache__"
