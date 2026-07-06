import tkinter as tk
from tkinter import ttk, messagebox
import subprocess, sys, os, threading


if getattr(sys, 'frozen', False):
    BASE = os.path.dirname(sys.executable)
else:
    BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
DEFAULTS = {
    'test_case' : 'dambreak',
    'finaltime' : '4.0',
    'save'      : '100',
    'mode_rk'   : '1',
    'mesh_mode' : 'lagrangian',
    'dx'        : '0.05',
}

def ecrire_config(params):
    dx      = float(params['dx'])
    h       = round(1.2 * dx, 10)
    mode_ale = 'wass'
    vtk_dir  = f"{params['test_case']}_{params['mesh_mode']}_{mode_ale}"

    # ── system/control_dict.py ───────────────────────────────
    chemin = os.path.join(BASE, 'system', 'control_dict.py')
    with open(chemin, 'w', encoding='utf-8') as f:
        f.write("import numpy as np\n\n")
        f.write(f"test_case = '{params['test_case']}'\n")
        f.write(f"finaltime = {params['finaltime']}\n")
        f.write(f"save      = {params['save']}\n")
        f.write(f"mode_rk   = {params['mode_rk']}\n")
        f.write(f"vtk_dir   = '{vtk_dir}'\n")

    # ── constant/ale_properties.py ───────────────────────────
    chemin = os.path.join(BASE, 'constant', 'ale_properties.py')
    with open(chemin, 'w', encoding='utf-8') as f:
        f.write("import numpy as np\n\n")
        f.write(f"mesh_mode = '{params['mesh_mode']}'\n")
        f.write(f"mode_ale  = '{mode_ale}'\n")
        f.write(f"beta_geom = 0.001\n")
        f.write(f"alpha     = 0.01\n")

    # ── constant/sph_properties.py ───────────────────────────
    chemin = os.path.join(BASE, 'constant', 'sph_properties.py')
    with open(chemin, 'w', encoding='utf-8') as f:
        f.write("import numpy as np\n\n")
        f.write(f"n_layers    = 3\n")
        f.write(f"dx          = {dx}\n")
        f.write(f"h           = {h}\n")
        f.write(f"CFL         = 0.1\n")
        f.write(f"FLUID, WALL = 0, 1\n")

    return vtk_dir

def lancer_simulation(params, log_widget, btn):

    # Fonction d'écriture dans le widget de logs
    # Doit être appelée depuis n'importe quel thread via root.after
    def log(msg):
        log_widget.config(state='normal')
        log_widget.insert('end', msg)
        log_widget.see('end')          # scroll automatique vers le bas
        log_widget.config(state='disabled')

    def run():
        # Désactiver le bouton pendant la simulation
        btn.config(state='disabled', text='Simulation en cours…')
        try:
            vtk_dir = ecrire_config(params)
            log(f'✔ Config écrite — dossier VTK : {vtk_dir}\n')
            log('─' * 50 + '\n')

            # Chemin vers run_solver.py
            solver_script = os.path.join(BASE, 'src', 'run_solver.py')

            # Lancement en sous-processus
            # -u = unbuffered : les print() arrivent immédiatement
            proc = subprocess.Popen(
                [sys.executable, '-u', solver_script],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,  # stderr redirigé vers stdout
                text=True,
                cwd=BASE                   # répertoire de travail = www/
            )

            # Lecture ligne par ligne en temps réel
            for line in proc.stdout:
                log(line)

            proc.wait()  # attendre la fin du processus

            if proc.returncode == 0:
                log('─' * 50 + '\n')
                log('✔ Simulation terminée avec succès.\n')
            else:
                log('─' * 50 + '\n')
                log(f'✘ Erreur (code {proc.returncode})\n')

        except Exception as e:
            log(f'✘ Erreur inattendue : {e}\n')
        finally:
            # Toujours réactiver le bouton, même en cas d'erreur
            btn.config(state='normal', text='▶  Lancer la simulation')

    # Lancement dans un thread séparé
    threading.Thread(target=run, daemon=True).start()
    
BG     = '#1e1e2e'   # fond principal
FG     = '#cdd6f4'   # texte
ACCENT = '#89b4fa'   # titres et bordures
ENTRY  = '#313244'   # fond des champs de saisie
VERT   = '#a6e3a1'   # couleur des logs

def build_gui():
    root = tk.Tk()
    root.title("SPHIRIT — Lanceur de simulation")
    root.resizable(False, False)
    root.configure(bg=BG)

    # ── Icône (implémentation future) ────────────────────────
    # Décommenter quand assets/sphirit.ico est disponible :
    # root.iconbitmap(os.path.join(BASE, 'assets', 'sphirit.ico'))

    # ── Style ttk ────────────────────────────────────────────
    style = ttk.Style()
    style.theme_use('clam')
    style.configure('TLabel',
                    background=BG, foreground=FG,
                    font=('Consolas', 10))
    style.configure('TEntry',
                    fieldbackground=ENTRY, foreground=FG,
                    font=('Consolas', 10))
    style.configure('TCombobox',
                    fieldbackground=ENTRY, foreground=FG,
                    font=('Consolas', 10))
    style.configure('TLabelframe',
                    background=BG, foreground=ACCENT,
                    font=('Consolas', 10, 'bold'))
    style.configure('TLabelframe.Label',
                    background=BG, foreground=ACCENT)

    pad = {'padx': 10, 'pady': 5}

    # ── Titre ────────────────────────────────────────────────
    tk.Label(root, text="SPHIRIT",
             bg=BG, fg=ACCENT,
             font=('Consolas', 20, 'bold')).grid(
                 row=0, column=0, columnspan=2, pady=(16, 2))

    tk.Label(root,
             text="Smoothed Particle Hydrodynamics Integrated Resource for Instruction and Teaching— Launcher",
             bg=BG, fg=FG,
             font=('Consolas', 9)).grid(
                 row=1, column=0, columnspan=2, pady=(0, 14))

    # ── Cadre paramètres ─────────────────────────────────────
    frame = ttk.LabelFrame(root, text="  Cellule 1 — Paramètres  ", padding=14)
    frame.grid(row=2, column=0, columnspan=2, padx=20, pady=6, sticky='ew')

    # Dictionnaire qui stocke les variables tkinter de chaque champ
    vars_ = {}

    # Helpers pour créer les lignes du formulaire proprement
    def ajouter_label(parent, ligne, texte):
        ttk.Label(parent, text=texte).grid(
            row=ligne, column=0, sticky='w', **pad)

    def ajouter_combo(parent, ligne, cle, valeurs, defaut):
        var = tk.StringVar(value=defaut)
        ttk.Combobox(parent, textvariable=var,
                     values=valeurs, state='readonly',
                     width=22).grid(row=ligne, column=1, sticky='w', **pad)
        vars_[cle] = var

    def ajouter_entry(parent, ligne, cle, defaut):
        var = tk.StringVar(value=defaut)
        ttk.Entry(parent, textvariable=var,
                  width=24).grid(row=ligne, column=1, sticky='w', **pad)
        vars_[cle] = var

    # Ligne 0 — test_case
    ajouter_label(frame, 0, 'test_case')
    ajouter_combo(frame, 0, 'test_case',
                  ['dambreak', 'TGV', 'funnel'],
                  DEFAULTS['test_case'])

    # Ligne 1 — finaltime
    ajouter_label(frame, 1, 'finaltime  (s)')
    ajouter_entry(frame, 1, 'finaltime', DEFAULTS['finaltime'])

    # Ligne 2 — save
    ajouter_label(frame, 2, 'save  (pas VTK)')
    ajouter_entry(frame, 2, 'save', DEFAULTS['save'])

    # Ligne 3 — mode_rk
    ajouter_label(frame, 3, 'mode_rk')
    ajouter_combo(frame, 3, 'mode_rk',
                  ['1 — Euler', '2 — RK2'],
                  '1 — Euler')

    # Ligne 4 — mesh_mode
    ajouter_label(frame, 4, 'mesh_mode')
    ajouter_combo(frame, 4, 'mesh_mode',
                  ['lagrangian', 'eulerian', 'wass'],
                  DEFAULTS['mesh_mode'])

    # Ligne 5 — dx
    ajouter_label(frame, 5, 'dx  (m)')
    ajouter_entry(frame, 5, 'dx', DEFAULTS['dx'])

    # ── Bouton Lancer ─────────────────────────────────────────
    btn = tk.Button(root,
                    text='▶  Lancer la simulation',
                    bg=ACCENT, fg='#1e1e2e',
                    font=('Consolas', 11, 'bold'),
                    relief='flat', cursor='hand2',
                    padx=20, pady=10)
    btn.grid(row=3, column=0, columnspan=2, pady=14)

    # ── Zone de logs ──────────────────────────────────────────
    log_frame = ttk.LabelFrame(root, text="  Logs  ", padding=6)
    log_frame.grid(row=4, column=0, columnspan=2,
                   padx=20, pady=(0, 18), sticky='ew')

    log_text = tk.Text(log_frame,
                       height=16, width=70,
                       bg='#11111b', fg=VERT,
                       font=('Consolas', 9),
                       state='disabled', relief='flat')
    scroll = ttk.Scrollbar(log_frame, command=log_text.yview)
    log_text.configure(yscrollcommand=scroll.set)
    log_text.grid(row=0, column=0)
    scroll.grid(row=0, column=1, sticky='ns')

    # ── Callback du bouton ────────────────────────────────────
    def on_lancer():
        # Validation et extraction des valeurs
        try:
            params = {
                'test_case' : vars_['test_case'].get(),
                'finaltime' : float(vars_['finaltime'].get()),
                'save'      : int(vars_['save'].get()),
                # "1 — Euler" → on prend uniquement le premier caractère
                'mode_rk'   : int(vars_['mode_rk'].get()[0]),
                'mesh_mode' : vars_['mesh_mode'].get(),
                'dx'        : float(vars_['dx'].get()),
            }
        except ValueError as e:
            messagebox.showerror("Paramètre invalide",
                                 f"Valeur incorrecte : {e}")
            return

        lancer_simulation(params, log_text, btn)

    btn.configure(command=on_lancer)
    root.mainloop()


# ── Point d'entrée ────────────────────────────────────────────
if __name__ == '__main__':
    build_gui()