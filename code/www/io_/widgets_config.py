import ipywidgets as widgets
from IPython.display import display

import system.control_dict           as ctrl
import system.sph_schemes            as schemes
import constant.ale_properties       as ale
import constant.sph_properties       as sph
import constant.numerical_properties as num

def _make_control_dict_widgets():
    # ── Widgets ───────────────────────────────────────────────────────────
    test_case_dropdown = widgets.Dropdown(
        options        =['dambreak','TGV','funnel'],
        value          ='dambreak',
        description    ='Test Case:',
        style          ={'description_width': '150px'},
    )
    
    finaltime_text            = widgets.FloatText(
        value          = 4.0,
        description    ='Final time (s):',
        style          ={'description_width': '150px'},
    )
    
    save_text          =widgets.IntText(
        value          =50,
        description    ='Pas de sauvegarde:',
        style          ={'description_width': '150px'},
    )
    
    mode_rk_dropdown   = widgets.Dropdown(
        options        =['Euler','RK2'],
        description    ='Time integration:',
        style          ={'description_width': '150px'},
    )
    
    # ── Assemblage ────────────────────────────────────────────────────────
    box = widgets.VBox([
        test_case_dropdown,
        finaltime_text,
        save_text,
        mode_rk_dropdown,
    ])

    # ── Dict de retour ────────────────────────────────────────────────────
    w = {
        'test_case_dropdown' : test_case_dropdown,
        'finaltime_text'     : finaltime_text,
        'save_text'          : save_text,
        'mode_rk_dropdown'   : mode_rk_dropdown,
    }

    return box, w

def _make_ale_properties_widgets():
    # ── Widgets ───────────────────────────────────────────────────────────
    mesh_mode_dropdown  = widgets.Dropdown(
        options         =['lagrangian','eulerian','wass'],
        value           ='lagrangian',
        description     ='Mesh mode:',
        style           ={'description_width': '150px'},
    )
    
    mode_ale_dropdown   = widgets.Dropdown(
        options         =['wass','loyd','shepard','coupled'],
        value           ='wass',
        description     ='ALE mode:',
        style           ={'description_width': '150px'},
    )
    
    beta_geom_slider    = widgets.FloatSlider(
    value = 0.001, min = 0.0001, max = 0.1, step = 0.001,
    description         ='ALE Method β:',
    style               ={'description_width': '150px'}
    )
    
    alpha_slider        = widgets.FloatSlider(
    value = 0.01, min = 0.0, max = 10.0, step = 0.1,
    description         ='ALE Method α:',
    style               ={'description_width': '150px'}
    )
    
    # ── Assemblage ────────────────────────────────────────────────────────
    box = widgets.VBox([
        mesh_mode_dropdown,
        mode_ale_dropdown,
        beta_geom_slider,
        alpha_slider,
    ])
    
    # ── Dict de retour ────────────────────────────────────────────────────
    w = {
        'mesh_mode_dropdown' : mesh_mode_dropdown,
        'mode_ale_dropdown'  : mode_ale_dropdown,
        'beta_geom_slider'   : beta_geom_slider,
        'alpha_slider'       : alpha_slider,
    }
    
    return box, w

def _make_sph_properties_widgets():
    # ── Widgets ───────────────────────────────────────────────────────────
    dx_slider           = widgets.FloatSlider(
        value = 0.05, min = 0.01, max = 0.15, step = 0.005,
        description     ='Discrétisation (dx):',
        style           ={'description_width': '150px'},
        readout_format  ='.3f',
    )
    
    h_ratio_slider      = widgets.FloatSlider(
        value =  1.2, min = 1.2, max = 2.5, step = 0.1,
        description     = 'Ratio h/dx:',
        style           ={'description_width': '150px'},
        readout_format  ='.1f',
    ) 
    
    cfl_slider          = widgets.FloatSlider(
        value = 0.2, min = 0.01, max = 0.5, step = 0.01,
        description     ='CFL:',
        style           ={'description_width': '150px'},
        readout_format  ='.2f',
    )
    
    # ── Assemblage ────────────────────────────────────────────────────────
    box = widgets.VBox([
        dx_slider,
        h_ratio_slider,
        cfl_slider,
    ])
    
    # ── Dict de retour ────────────────────────────────────────────────────
    w= {
        'dx_slider'         : dx_slider,
        'h_ratio_slider'    : h_ratio_slider,
        'cfl_slider'        : cfl_slider,
    }
     
    return box, w

def _make_numerical_properties_schemes_widgets():
    # ── Widgets ───────────────────────────────────────────────────────────
    mode_p_refinement_dropdown = widgets.Dropdown(
        options                = {'00 - defaut' : 0,'01 - MLS gradient' : 1},
        value                  = 0,
        description            ='Odre spatial :',
        style                  ={'description_width': '150px'},
    )
    
    visco_on_toggle            = widgets.ToggleButton(
        value                  = False,
        description            ='Viscosité : ',
        indent                 = False,
        style                  ={'description_width': '150px'},
    )
    
    shepard_corr_toggle        = widgets.ToggleButton(
        value                  = False,
        description            ='Shepard corr : ',
        indent                 = False,
        style                  ={'description_width': '150px'},
    )
    
    shift_part_toggle          = widgets.ToggleButton(
        value                  = False,
        description            ='Shift particles : ',
        indent                 = False,
        style                  ={'description_width': '150px'},
    )
    
    random_pos_toggle          = widgets.ToggleButton(
        value                  = True,
        description            ='Random pos : ',
        indent                 = False,
        style                  ={'description_width': '150px'},
    )
    
    # ── Assemblage ────────────────────────────────────────────────────────
    box = widgets.VBox([
        mode_p_refinement_dropdown,
        visco_on_toggle,
        shepard_corr_toggle,
        shift_part_toggle,
        random_pos_toggle,
    ])
    
    # ── Dict de retour ────────────────────────────────────────────────────
    w = {
        'mode_p_refinement_dropdown' : mode_p_refinement_dropdown,
        'visco_on_toggle'            : visco_on_toggle,
        'shepard_corr_toggle'        : shepard_corr_toggle,
        'shift_part_toggle'          : shift_part_toggle,
        'random_pos_toggle'          : random_pos_toggle,
    }
    
    return box, w

def _apply_config (w):
    # ── Control dict ──────────────────────────────────────────────────────
    ctrl.test_case         = w['test_case_dropdown'].value
    ctrl.finaltime         = w['finaltime_text'].value
    ctrl.save              = w['save_text'].value
    ctrl.mode_rk           = w['mode_rk_dropdown'].value
    
    # ── ALE properties ────────────────────────────────────────────────────
    ale.mesh_mode          = w['mesh_mode_dropdown'].value
    ale.mode_ale           = w['mode_ale_dropdown'].value
    ale.beta_geom          = w['beta_geom_slider'].value
    ale.alpha              = w['alpha_slider'].value
    
    # ── SPH properties ────────────────────────────────────────────────────
    sph.dx                 = w['dx_slider'].value
    sph.h_ratio            = w['h_ratio_slider'].value
    sph.cfl                = w['cfl_slider'].value
    
    # ── Numerical properties ──────────────────────────────────────────────
    num.mode_p_refinement  = w['mode_p_refinement_dropdown'].value
    num.visco_on            = w['visco_on_toggle'].value
    num.shepard_corr        = w['shepard_corr_toggle'].value
    num.shift_part          = w['shift_part_toggle'].value
    num.random_pos          = w['random_pos_toggle'].value
    
    ctrl.vtk_dir = f"{ctrl.test_case}_{ale.mesh_mode}_{ale.mode_ale}"
    print(f"Cas test     : {ctrl.test_case}")
    print(f"Mode maillage: {ale.mesh_mode}")
    print(f"Durée        : {ctrl.finaltime} s")
    print(f"Dossier VTK  : {ctrl.vtk_dir}")    

def _display_widgets():
    sim_box, w_sim = _make_control_dict_widgets()
    ale_box, w_ale = _make_ale_properties_widgets()
    sph_box, w_sph = _make_sph_properties_widgets()
    num_box, w_num = _make_numerical_properties_schemes_widgets()
    
    w = {**w_sim, **w_ale, **w_sph, **w_num}
    
    output_config = widgets.Output()
    
    apply_btn           = widgets.Button(
        description     ="Appliquer la configuration",
        button_style    ='primary',
        icon            ='check',
    )
    
    def _on_apply(b):
        output_config.clear_output()
        with output_config:
            _apply_config(w)
            
    apply_btn.on_click(_on_apply)
    
    root = widgets.VBox([
        widgets.HTML("<b>── Simulation ──</b>"),  sim_box,
        widgets.HTML("<b>── Maillage / ALE ──</b>"), ale_box,
        widgets.HTML("<b>── SPH ──</b>"),         sph_box,
        widgets.HTML("<b>── Options numériques ──</b>"), num_box,
        apply_btn,
        output_config
    ])
    
    display(root)
    with output_config:
        _apply_config(w)