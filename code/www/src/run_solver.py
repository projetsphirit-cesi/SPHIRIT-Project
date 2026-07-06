import os, sys, importlib, time

# ── Résolution du chemin racine www/ ────────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE)
#─────────────────────────────────────────────────────────────────────────────────────────────────

# ── Import des modules de configuration ─────────────────────────────────────────────────────────
import system.control_dict           as ctrl
import constant.physical_properties  as phys
import constant.sph_properties       as sph 
import constant.ale_properties       as ale 
import constant.numerical_properties as num 
from src.Solver import run                  
#─────────────────────────────────────────────────────────────────────────────────────────────────

# ── Rechargement forcé (les fichiers viennent d'être réécrits) ──────────────────────────────────
for mod in [ctrl, phys, sph, ale, num]:     
    importlib.reload(mod)
    
sys.stdout.reconfigure(line_buffering=True) 
#─────────────────────────────────────────────────────────────────────────────────────────────────

#── Initialisation du cas test ───────────────────────────────────────────────────────────────────
case_module = importlib.import_module(f'cases.{ctrl.test_case}.setup')
state = case_module.init(
    vars(ctrl), vars(phys), vars(sph), vars(ale), vars(num)
)
print(f'Init OK : {len(state["pos"])} particules — dt={state["dt"]:.2e}s', flush=True)
#─────────────────────────────────────────────────────────────────────────────────────────────────

#── Lancement de la simulation ───────────────────────────────────────────────────────────────────
t0 = time.perf_counter()
run(state, vars(ctrl), vars(phys), vars(sph))
print(f'Simulation terminée en {time.perf_counter()-t0:.1f} s', flush=True)
#─────────────────────────────────────────────────────────────────────────────────────────────────