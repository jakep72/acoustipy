from acoustipy import AcousticTMM
from acoustipy import AcousticID
import os
import pytest

def test_param_init():
    with pytest.raises(ValueError):
        AcousticID()

def test_inverse(tmpdir):
    gap_file = tmpdir.mkdir("sub")
    no_gap_file = tmpdir.mkdir("sub2")

    structure = AcousticTMM(incidence='Normal',air_temperature = 20)

    layer1 = structure.Add_JCA_Layer(30, 46182,.917,2.1,83,128)

    air = structure.Add_Air_Layer(thickness = 100)

    s1 = structure.assemble_structure(layer1)

    A1 = structure.reflection(s1)

    no_gap = os.path.join(no_gap_file,'no_gap.csv')
    structure.to_csv(no_gap,A1)

    s2 = structure.assemble_structure(layer1,air)

    A2 = structure.reflection(s2)

    gap = os.path.join(gap_file,'gap.csv')
    structure.to_csv(gap,A2)

    inv = AcousticID(mount_type='Dual',no_gap_file=no_gap_file.join("no_gap.csv"),gap_file = gap_file.join("gap.csv"),input_type='reflection',air_temperature=20)

    res = inv.Inverse(thickness=30,
                      flow_resistivity=46182,
                      porosity=.917,
                      air_gap=100,
                      uncertainty = .10)
    stats = inv.stats(res)

    assert stats['r_value'] > .99

def test_indirect(tmpdir):
    gap_file = tmpdir.mkdir("sub")
    no_gap_file = tmpdir.mkdir("sub2")

    structure = AcousticTMM(incidence='Normal',air_temperature = 20)

    layer1 = structure.Add_JCA_Layer(30, 46182,.917,2.1,83,128)

    air = structure.Add_Air_Layer(thickness = 100)

    s1 = structure.assemble_structure(layer1)

    A1 = structure.reflection(s1)

    no_gap = os.path.join(no_gap_file,'no_gap.csv')
    structure.to_csv(no_gap,A1)

    s2 = structure.assemble_structure(layer1,air)

    A2 = structure.reflection(s2)

    gap = os.path.join(gap_file,'gap.csv')
    structure.to_csv(gap,A2)

    inv = AcousticID(mount_type='Dual',no_gap_file=no_gap_file.join("no_gap.csv"),gap_file = gap_file.join("gap.csv"),input_type='reflection',air_temperature=20)

    res = inv.Indirect(thickness=30,
                       porosity=.917,
                       flow_resistivity=46182,
                       air_gap=100)
    stats = inv.stats(res)

    assert stats['r_value'] > .99

def test_hybrid(tmpdir):
    gap_file = tmpdir.mkdir("sub")
    no_gap_file = tmpdir.mkdir("sub2")

    structure = AcousticTMM(incidence='Normal',air_temperature = 20)

    layer1 = structure.Add_JCA_Layer(30, 46182,.917,2.1,83,128)

    air = structure.Add_Air_Layer(thickness = 100)

    s1 = structure.assemble_structure(layer1)

    A1 = structure.reflection(s1)

    no_gap = os.path.join(no_gap_file,'no_gap.csv')
    structure.to_csv(no_gap,A1)

    s2 = structure.assemble_structure(layer1,air)

    A2 = structure.reflection(s2)

    gap = os.path.join(gap_file,'gap.csv')
    structure.to_csv(gap,A2)

    inv = AcousticID(mount_type='Dual',no_gap_file=no_gap_file.join("no_gap.csv"),gap_file = gap_file.join("gap.csv"),input_type='reflection',air_temperature=20)

    res = inv.Hybrid(30,.917,air_gap=100,uncertainty = .10)
    stats = inv.stats(res)

    assert stats['r_value'] > .99