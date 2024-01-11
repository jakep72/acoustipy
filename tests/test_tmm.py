from acoustipy import AcousticTMM
import matplotlib.pyplot as plt 
import pytest

def lengths(normal_layer, diffuse_layer):
    normal_array = normal_layer[0].shape
    diffuse_array = diffuse_layer[0].shape

    return(normal_array, diffuse_array)

def test_tmm_init():
    AcousticTMM()

def test_air_layer():
    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    normal_layer = structure.Add_Air_Layer()
    structure_diffuse = AcousticTMM(incidence='Diffuse',fmin=4, fmax=3200, fs=4, angles=[0,79,1], air_temperature=20)
    diffuse_layer = structure_diffuse.Add_Air_Layer()
    assert lengths(normal_layer, diffuse_layer) == ((2,2,800), (2,2,800,79))

def test_db_layer():
    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    normal_layer = structure.Add_DB_Layer(thickness=10,
                                          flow_resistivity=100000)
    structure_diffuse = AcousticTMM(incidence='Diffuse',fmin=4, fmax=3200, fs=4, angles=[0,79,1], air_temperature=20)
    diffuse_layer = structure_diffuse.Add_DB_Layer(thickness=10,
                                          flow_resistivity=100000)
    assert lengths(normal_layer, diffuse_layer) == ((2,2,800), (2,2,800,79))

def test_dbm_layer():
    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    normal_layer = structure.Add_DBM_Layer(thickness=10,
                                           flow_resistivity=100000)
    structure_diffuse = AcousticTMM(incidence='Diffuse',fmin=4, fmax=3200, fs=4, angles=[0,79,1], air_temperature=20)
    diffuse_layer = structure_diffuse.Add_DBM_Layer(thickness=10,
                                           flow_resistivity=100000)
    assert lengths(normal_layer, diffuse_layer) == ((2,2,800), (2,2,800,79))

def test_jca_layer():
    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    normal_layer = structure.Add_JCA_Layer(thickness=10,
                                     flow_resistivity=100000,
                                     porosity=.9,
                                     tortuosity=2.5,
                                     viscous_characteristic_length=10,
                                     thermal_characteristic_length=10)
    structure_diffuse = AcousticTMM(incidence='Diffuse',fmin=4, fmax=3200, fs=4, angles=[0,79,1], air_temperature=20)
    diffuse_layer = structure_diffuse.Add_JCA_Layer(thickness=10,
                                     flow_resistivity=100000,
                                     porosity=.9,
                                     tortuosity=2.5,
                                     viscous_characteristic_length=10,
                                     thermal_characteristic_length=10)
    assert lengths(normal_layer, diffuse_layer) == ((2,2,800), (2,2,800,79))

def test_jcal_layer():
    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    normal_layer = structure.Add_JCAL_Layer(thickness=10,
                                     flow_resistivity=100000,
                                     porosity=.9,
                                     tortuosity=2.5,
                                     viscous_characteristic_length=10,
                                     thermal_characteristic_length=10,
                                     thermal_permeability=60)
    structure_diffuse = AcousticTMM(incidence='Diffuse',fmin=4, fmax=3200, fs=4, angles=[0,79,1], air_temperature=20)
    diffuse_layer = structure_diffuse.Add_JCAL_Layer(thickness=10,
                                     flow_resistivity=100000,
                                     porosity=.9,
                                     tortuosity=2.5,
                                     viscous_characteristic_length=10,
                                     thermal_characteristic_length=10,
                                     thermal_permeability=60)
    assert lengths(normal_layer, diffuse_layer) == ((2,2,800), (2,2,800,79))

def test_jcapl_layer():
    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    normal_layer = structure.Add_JCAPL_Layer(thickness=10,
                                     flow_resistivity=100000,
                                     porosity=.9,
                                     tortuosity=2.5,
                                     viscous_characteristic_length=10,
                                     thermal_characteristic_length=10,
                                     thermal_permeability=60,
                                     thermal_tortuosity=10,
                                     viscous_tortuosity=10)
    structure_diffuse = AcousticTMM(incidence='Diffuse',fmin=4, fmax=3200, fs=4, angles=[0,79,1], air_temperature=20)
    diffuse_layer = structure_diffuse.Add_JCAPL_Layer(thickness=10,
                                     flow_resistivity=100000,
                                     porosity=.9,
                                     tortuosity=2.5,
                                     viscous_characteristic_length=10,
                                     thermal_characteristic_length=10,
                                     thermal_permeability=60,
                                     thermal_tortuosity=10,
                                     viscous_tortuosity=10)
    assert lengths(normal_layer, diffuse_layer) == ((2,2,800), (2,2,800,79))

def test_horoshenkov_layer():
    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    normal_layer = structure.Add_Horoshenkov_Layer(thickness=10,
                                            porosity=0.9,
                                            median_pore_size=147,
                                            pore_size_distribution=.325)
    structure_diffuse = AcousticTMM(incidence='Diffuse',fmin=4, fmax=3200, fs=4, angles=[0,79,1], air_temperature=20)
    diffuse_layer = structure_diffuse.Add_Horoshenkov_Layer(thickness=10,
                                            porosity=0.9,
                                            median_pore_size=147,
                                            pore_size_distribution=.325)
    assert lengths(normal_layer, diffuse_layer) == ((2,2,800), (2,2,800,79))

def test_biot_limp_layer():
    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    normal_layer = structure.Add_Biot_Limp_Layer('JCA',
                                          thickness=25.4,
                                          flow_resistivity=50000,
                                          mass_density=200,
                                          porosity=.90,
                                          tortuosity=1.2,
                                          viscous_characteristic_length=80,
                                          thermal_characteristic_length=110)
    structure_diffuse = AcousticTMM(incidence='Diffuse',fmin=4, fmax=3200, fs=4, angles=[0,79,1], air_temperature=20)
    diffuse_layer = structure_diffuse.Add_Biot_Limp_Layer('JCA',
                                          thickness=25.4,
                                          flow_resistivity=50000,
                                          mass_density=200,
                                          porosity=.90,
                                          tortuosity=1.2,
                                          viscous_characteristic_length=80,
                                          thermal_characteristic_length=110)
    assert lengths(normal_layer, diffuse_layer) == ((2,2,800), (2,2,800,79))

def test_biot_rigid_layer():
    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    normal_layer = structure.Add_Biot_Rigid_Layer('JCA',
                                          thickness=25.4,
                                          flow_resistivity=50000,
                                          mass_density=200,
                                          porosity=.90,
                                          tortuosity=1.2,
                                          viscous_characteristic_length=80,
                                          thermal_characteristic_length=110)
    structure_diffuse = AcousticTMM(incidence='Diffuse',fmin=4, fmax=3200, fs=4, angles=[0,79,1], air_temperature=20)
    diffuse_layer = structure_diffuse.Add_Biot_Rigid_Layer('JCA',
                                          thickness=25.4,
                                          flow_resistivity=50000,
                                          mass_density=200,
                                          porosity=.90,
                                          tortuosity=1.2,
                                          viscous_characteristic_length=80,
                                          thermal_characteristic_length=110)
    assert lengths(normal_layer, diffuse_layer) == ((2,2,800), (2,2,800,79))

def test_screen_layer():
    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    normal_layer = structure.Add_Resistive_Screen(thickness=2.54,
                                           flow_resistivity=50000,
                                           porosity=.95)
    structure_diffuse = AcousticTMM(incidence='Diffuse',fmin=4, fmax=3200, fs=4, angles=[0,79,1], air_temperature=20)
    diffuse_layer = structure_diffuse.Add_Resistive_Screen(thickness=2.54,
                                           flow_resistivity=50000,
                                           porosity=.95)
    assert lengths(normal_layer, diffuse_layer) == ((2,2,800), (2,2,800,79))

def test_maa_mpp_layer():
    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    normal_layer = structure.Add_MAA_MPP_Layer(thickness=2.54,
                                        pore_diameter=0.5,
                                        c_to_c_dist=1)
    structure_diffuse = AcousticTMM(incidence='Diffuse',fmin=4, fmax=3200, fs=4, angles=[0,79,1], air_temperature=20)
    diffuse_layer = structure_diffuse.Add_MAA_MPP_Layer(thickness=2.54,
                                        pore_diameter=0.5,
                                        c_to_c_dist=1)
    assert lengths(normal_layer, diffuse_layer) == ((2,2,800), (2,2,800,79))

def test_ef_mpp_layer():
    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    normal_layer = structure.Add_MPP_EF_Layer(thickness=2.54,
                                        pore_diameter=0.5,
                                        c_to_c_dist=1)
    structure_diffuse = AcousticTMM(incidence='Diffuse',fmin=4, fmax=3200, fs=4, angles=[0,79,1], air_temperature=20)
    diffuse_layer = structure_diffuse.Add_MPP_EF_Layer(thickness=2.54,
                                        pore_diameter=0.5,
                                        c_to_c_dist=1)
    assert lengths(normal_layer, diffuse_layer) == ((2,2,800), (2,2,800,79))

# def test_layer_to_database():
#     structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
#     structure.Add_DB_Layer(25.4, 50000,save_layer=True,layer_name='test_DB')

def test_air_props():
    structure = AcousticTMM(air_temperature=25)
    density = structure.density_temp
    soundspeed = structure.soundspeed_temp
    gamma = structure.gamma_temp
    viscosity = structure.viscosity_temp
    pr = structure.Pr_temp
    Z0 = structure.Z0
    assert (density, soundspeed, gamma, viscosity, pr, Z0) == (1.1868857215,
                                                           346.16,
                                                           1.4015208336019964,
                                                           1.8440511979580485e-05,
                                                           0.714729563521373,
                                                           410.85236135444)

def test_frequencies():
    with pytest.raises(ValueError):
        AcousticTMM(fmin=0)
    
    with pytest.raises(ValueError):
        AcousticTMM(fs=0)
    
    with pytest.raises(ValueError):
        AcousticTMM(fmin=4, fmax=4)

    structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs=4, air_temperature=20)
    lin_freq_len = len(structure.frequency)
    ang_freq_len = len(structure.ang_freq)
    k0_len = len(structure.k0)

    assert (lin_freq_len, ang_freq_len, k0_len) == (800,800,800)

def test_angles():
    with pytest.raises(ValueError):
        AcousticTMM(angles=[2,1,1])

    with pytest.raises(ValueError):
        AcousticTMM(angles=[0,79,3])

def test_transfer_matrix():
    normal_structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs = 4, air_temperature=20)

    layer_normal = normal_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    normal_transfer_matrix = normal_structure.assemble_structure(layer_normal)

    diffuse_structure = AcousticTMM(incidence='Diffuse', fmin=4, fmax=3200, fs = 4, angles=[0,79,1])

    layer_diffuse = diffuse_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    diffuse_transfer_matrix = diffuse_structure.assemble_structure(layer_diffuse)

    assert lengths(normal_transfer_matrix, diffuse_transfer_matrix) == ((2,2,800),(2,2,800,79))

def test_reflection():
    normal_structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs = 4, air_temperature=20)

    layer_normal = normal_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    normal_transfer_matrix = normal_structure.assemble_structure(layer_normal)

    diffuse_structure = AcousticTMM(incidence='Diffuse', fmin=4, fmax=3200, fs = 4, angles=[0,79,1])

    layer_diffuse = diffuse_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    diffuse_transfer_matrix = diffuse_structure.assemble_structure(layer_diffuse)

    normal_r = normal_structure.reflection(normal_transfer_matrix)

    diffuse_r = diffuse_structure.reflection(diffuse_transfer_matrix)

    assert (normal_r.shape, diffuse_r.shape) == ((800,2),(800,2))

def test_absorption():
    normal_structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs = 4, air_temperature=20)

    layer_normal = normal_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    normal_transfer_matrix = normal_structure.assemble_structure(layer_normal)

    diffuse_structure = AcousticTMM(incidence='Diffuse', fmin=4, fmax=3200, fs = 4, angles=[0,79,1])

    layer_diffuse = diffuse_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    diffuse_transfer_matrix = diffuse_structure.assemble_structure(layer_diffuse)

    normal_abs = normal_structure.absorption(normal_transfer_matrix)

    diffuse_abs = diffuse_structure.absorption(diffuse_transfer_matrix)

    assert (normal_abs.shape, diffuse_abs.shape) == ((800,2),(800,2))

def test_tl():
    normal_structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs = 4, air_temperature=20)

    layer_normal = normal_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    normal_transfer_matrix = normal_structure.assemble_structure(layer_normal)

    diffuse_structure = AcousticTMM(incidence='Diffuse', fmin=4, fmax=3200, fs = 4, angles=[0,79,1])

    layer_diffuse = diffuse_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    diffuse_transfer_matrix = diffuse_structure.assemble_structure(layer_diffuse)

    normal_tl = normal_structure.transmission_loss(normal_transfer_matrix)

    diffuse_tl = diffuse_structure.transmission_loss(diffuse_transfer_matrix)

    assert (normal_tl.shape, diffuse_tl.shape) == ((800,2),(800,2))

def test_octaves():
    normal_structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs = 4, air_temperature=20)

    layer_normal = normal_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    normal_transfer_matrix = normal_structure.assemble_structure(layer_normal)

    diffuse_structure = AcousticTMM(incidence='Diffuse', fmin=4, fmax=3200, fs = 4, angles=[0,79,1])

    layer_diffuse = diffuse_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    diffuse_transfer_matrix = diffuse_structure.assemble_structure(layer_diffuse)

    normal_abs = normal_structure.absorption(normal_transfer_matrix)

    diffuse_abs = diffuse_structure.absorption(diffuse_transfer_matrix)

    normal_bands = normal_structure.octave_bands(normal_abs, kind='THIRD_OCTAVE')

    diffuse_bands = diffuse_structure.octave_bands(diffuse_abs, kind='THIRD_OCTAVE')

    assert (normal_bands.shape, diffuse_bands.shape) == ((24,2),(24,2))

def test_ffa():
    normal_structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs = 4, air_temperature=20)

    layer_normal = normal_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    normal_transfer_matrix = normal_structure.assemble_structure(layer_normal)

    diffuse_structure = AcousticTMM(incidence='Diffuse', fmin=4, fmax=3200, fs = 4, angles=[0,79,1])

    layer_diffuse = diffuse_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    diffuse_transfer_matrix = diffuse_structure.assemble_structure(layer_diffuse)

    normal_abs = normal_structure.absorption(normal_transfer_matrix)

    diffuse_abs = diffuse_structure.absorption(diffuse_transfer_matrix)

    normal_bands = normal_structure.octave_bands(normal_abs, kind='THIRD_OCTAVE')

    diffuse_bands = diffuse_structure.octave_bands(diffuse_abs, kind='THIRD_OCTAVE')

    normal_ffa = normal_structure.FFA(normal_bands)

    diffuse_ffa = diffuse_structure.FFA(diffuse_bands)

    assert (normal_ffa, diffuse_ffa) == (0.157, 0.206)

def test_saa():
    normal_structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs = 4, air_temperature=20)

    layer_normal = normal_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    normal_transfer_matrix = normal_structure.assemble_structure(layer_normal)

    diffuse_structure = AcousticTMM(incidence='Diffuse', fmin=4, fmax=3200, fs = 4, angles=[0,79,1])

    layer_diffuse = diffuse_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    diffuse_transfer_matrix = diffuse_structure.assemble_structure(layer_diffuse)

    normal_abs = normal_structure.absorption(normal_transfer_matrix)

    diffuse_abs = diffuse_structure.absorption(diffuse_transfer_matrix)

    normal_bands = normal_structure.octave_bands(normal_abs, kind='THIRD_OCTAVE')

    diffuse_bands = diffuse_structure.octave_bands(diffuse_abs, kind='THIRD_OCTAVE')

    normal_saa = normal_structure.SAA(normal_bands)

    diffuse_saa = diffuse_structure.SAA(diffuse_bands)

    assert (normal_saa, diffuse_saa) == (0.161, 0.209)

def test_plots(monkeypatch):
    monkeypatch.setattr(plt, 'show', lambda: None)
    normal_structure = AcousticTMM(incidence='Normal', fmin=4, fmax=3200, fs = 4, air_temperature=20)

    layer_normal = normal_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    normal_transfer_matrix = normal_structure.assemble_structure(layer_normal)

    diffuse_structure = AcousticTMM(incidence='Diffuse', fmin=4, fmax=3200, fs = 4, angles=[0,79,1])

    layer_diffuse = diffuse_structure.Add_DBM_Layer(thickness=10,
                                     flow_resistivity=100000)
    
    diffuse_transfer_matrix = diffuse_structure.assemble_structure(layer_diffuse)

    normal_abs = normal_structure.absorption(normal_transfer_matrix)

    diffuse_abs = diffuse_structure.absorption(diffuse_transfer_matrix)

    normal_bands = normal_structure.octave_bands(normal_abs, kind='THIRD_OCTAVE')

    diffuse_bands = diffuse_structure.octave_bands(diffuse_abs, kind='THIRD_OCTAVE')

    normal_structure.plot_curve([normal_abs,normal_bands],["absorption","third octave"])

    diffuse_structure.plot_curve([diffuse_abs,diffuse_bands],["absorption","third octave"])
