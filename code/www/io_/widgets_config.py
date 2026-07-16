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
    
    finaltime_text            = widgets.Text(
        value          = '4.0',
        description    ='Final time (s):',
        style          ={'description_width': '150px'},
    )
    
    save_text          =widgets.Text(
        value          ='50',
        description    ='Pas de sauvegarde:',
        style          ={'description_width': '150px'},
    )
    
    mode_rk_dropdown   = widgets.Dropdown(
        options        =['Euler','RK2'],
        description    ='Time integration:',
        style          ={'description_width': '150px'},
    )
    
    # ── Assemblage ────────────────────────────────────────────────────────
    box = widgets.GridBox([
        test_case_dropdown,
        finaltime_text,
        save_text,
        mode_rk_dropdown,
        ],
        layout = widgets.Layout(grid_template_columns = "1fr 1fr",width = "100%"),
    )

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
    box = widgets.GridBox([
        mesh_mode_dropdown,
        mode_ale_dropdown,
        beta_geom_slider,
        alpha_slider,
    ],
        layout = widgets.Layout(grid_template_columns = "1fr 1fr",width = "100%"),
    )
    
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
    box = widgets.GridBox([
        dx_slider,
        h_ratio_slider,
        cfl_slider,
    ],
        layout = widgets.Layout(grid_template_columns = "1fr 1fr",width = "100%"),
    )
    
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
        description            ='Viscosité',
        indent                 = False,
        layout                 = widgets.Layout(width = 'auto', min_width = '0'),
    )
    
    shepard_corr_toggle        = widgets.ToggleButton(
        value                  = False,
        description            ='Shepard corr',
        indent                 = False,
        layout                 = widgets.Layout(width = 'auto', min_width = '0'),
    )
    
    shift_part_toggle          = widgets.ToggleButton(
        value                  = False,
        description            ='Shift particles',
        indent                 = False,
        layout                 = widgets.Layout(width = 'auto', min_width = '0'),
    )
    
    random_pos_toggle          = widgets.ToggleButton(
        value                  = True,
        description            ='Random pos',
        indent                 = False,
        layout                 = widgets.Layout(width = 'auto', min_width = '0'),
    )
    
    # ── Assemblage ────────────────────────────────────────────────────────
    box = widgets.GridBox([
        mode_p_refinement_dropdown,
        widgets.Label(''),
        visco_on_toggle,
        shepard_corr_toggle,
        shift_part_toggle,
        random_pos_toggle,
    ],
        layout = widgets.Layout(grid_template_columns = "1fr 1fr",width = "100%", overflow = 'hidden'),
    )
    
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
    try:
        ctrl.finaltime = float(w['finaltime_text'].value)
    except ValueError:
        print("⚠  finaltime invalide — valeur conservée")
    try:
        ctrl.save = int(w['save_text'].value)
    except ValueError:
        print("⚠  save invalide — valeur conservée")
    ctrl.mode_rk           = 1 if w['mode_rk_dropdown'].value == 'Euler' else 2
    
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
    
    ctrl.config_validated = True

def _make_section_title(label):
    return widgets.HTML(
        f'<div style="'
        f'display:flex; align-items:center; '
        f'margin:12px 0 4px 0; gap:8px;">'
        f'<hr style="flex:1; border:none; border-top:1px solid #3c3c3c; margin:0;">'
        f'<span style="'
        f'font-weight:bold; font-family:monospace; '
        f'font-size:13px; color:#969696; white-space:nowrap;">'
        f'{label}'
        f'</span>'
        f'<hr style="flex:1; border:none; border-top:1px solid #3c3c3c; margin:0;">'
        f'</div>',
        layout=widgets.Layout(width='100%'),
    )

def _display_widgets():
    
    display(widgets.HTML("""<style>
    .jp-OutputArea-output,
    .jp-Cell-outputArea        { background-color : #1e1e1e !important;
                                 overflow-x        : hidden  !important; }
    .widget-vbox,
    .widget-gridbox            { background-color : #1e1e1e !important; }
    .widget-label,
    .widget-readout            { color : #cccccc !important; }
    .jp-OutputArea-output pre  { color : #cccccc !important; }

    .widget-dropdown select    { background-color : #3c3c3c !important;
                                 color            : #cccccc !important;
                                 border           : 1px solid #555 !important;
                                 border-radius    : 6px    !important; }
    .widget-text input,
    .widget-floattext input,
    .widget-inttext input      { background-color : #3c3c3c !important;
                                 color            : #cccccc !important;
                                 border           : 1px solid #555 !important;
                                 border-radius    : 6px    !important; }
    .widget-toggle-button      { background-color : #3c3c3c !important;
                                 color            : #cccccc !important;
                                 border           : 1px solid #555 !important;
                                 border-radius    : 6px    !important; }
    .widget-toggle-button.mod-active {
                                 background-color : #0078d4 !important;
                                 color            : #ffffff !important;
                                 border-color     : #005fa3 !important; }

    .widget-button             { border-radius : 6px !important; }

    .jp-OutputArea-output,
    .jp-Cell-outputArea,
    .widget-vbox,
    .widget-gridbox            { overflow-x : hidden !important; }
    .widget-floattext input::-webkit-inner-spin-button,
    .widget-floattext input::-webkit-outer-spin-button,
    .widget-inttext   input::-webkit-inner-spin-button,
    .widget-inttext   input::-webkit-outer-spin-button {
        -webkit-appearance : none !important;
        margin             : 0    !important;
    }
    </style>"""))
    
    sim_box, w_sim = _make_control_dict_widgets()
    ale_box, w_ale = _make_ale_properties_widgets()
    sph_box, w_sph = _make_sph_properties_widgets()
    num_box, w_num = _make_numerical_properties_schemes_widgets()
    
    w = {**w_sim, **w_ale, **w_sph, **w_num}
    
    output_config = widgets.Output(
        layout=widgets.Layout(
            border        = '1px solid #3c3c3c',
            padding       = '8px 12px',
            border_radius = '6px',
            )
    )

    apply_btn = widgets.Button(
        description  = 'Appliquer la configuration',
        button_style = 'primary',
        icon         = 'check',
        layout       = widgets.Layout(
            width         = '100%',
            height        = '100%',
            border_radius = '6px',
        ),
    )
    apply_btn.style.font_size    = '14px'
    apply_btn.style.font_weight  = 'bold'
    apply_btn.style.button_color = '#0078d4'

    def _on_apply(b):
        output_config.clear_output(wait=True)
        with output_config:
            _apply_config(w)

    apply_btn.on_click(_on_apply)

    btn_wrapper = widgets.VBox(
        [apply_btn],
        layout       = widgets.Layout(
            width    = '100%',
            height   = '100%',
            overflow = 'hidden',
        ),
    )

    footer = widgets.GridBox(
        [btn_wrapper, output_config],
        layout                    = widgets.Layout(
            grid_template_columns = '1fr 2fr',
            grid_template_rows    = '1fr',
            width                 = '100%',
            align_items           = 'stretch',
            gap                   = '8px',
            overflow              = 'hidden',
        ),
    )

    root = widgets.VBox([
        _make_section_title('Simulation'),          sim_box,
        _make_section_title('Maillage / ALE'),      ale_box,
        _make_section_title('SPH'),                 sph_box,
        _make_section_title('Options numériques'),  num_box,
        footer,
        ],
        layout       = widgets.Layout(
            width    = '100%',
            overflow = 'hidden',
            padding  = '12px',
        ),
    )

    display(root)
    with output_config:
        _apply_config(w)
    ctrl.config_validated = False