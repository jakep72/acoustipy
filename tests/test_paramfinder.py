from acoustipy import AcousticTMM
from acoustipy import AcousticID
import os
import pytest
import torch


class TestAcousticIDInit:
    """Tests for AcousticID initialization and validation."""
    
    def test_param_init_no_file(self):
        """Test that initializing without required files raises ValueError."""
        with pytest.raises(ValueError):
            AcousticID()
    
    def test_param_init_invalid_mount_type(self, tmpdir):
        """Test that invalid mount_type is handled gracefully."""
        # Create a dummy file
        no_gap_file = tmpdir.join("no_gap.csv")
        no_gap_file.write("100,0.5\n200,0.6\n")
        
        # Invalid mount type should still work but won't load data
        inv = AcousticID(mount_type='InvalidType', no_gap_file=str(no_gap_file))
        assert inv.opt_type == 'InvalidType'
    
    def test_param_init_no_gap_absorption(self, tmpdir):
        """Test initialization with No Gap mount type and absorption data."""
        no_gap_file = tmpdir.join("no_gap.csv")
        no_gap_file.write("100,0.5\n200,0.6\n")
        
        inv = AcousticID(mount_type='No Gap', no_gap_file=str(no_gap_file), input_type='absorption')
        assert inv.opt_type == 'No Gap'
        assert inv.input_type == 'absorption'
        assert inv.no_gap_data is not None

    def test_param_init_gap_absorption(self, tmpdir):
        """Test initialization with Gap mount type and absorption data."""
        gap_file = tmpdir.join("gap.csv")
        gap_file.write("100,0.5\n200,0.6\n")
        
        inv = AcousticID(mount_type='Gap', gap_file=str(gap_file), input_type='absorption')
        assert inv.opt_type == 'Gap'
        assert inv.gap_data is not None


# Keep backward compatibility with old test name
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

def test_ml(tmpdir):
    no_gap_file = tmpdir.mkdir("sub2")

    structure = AcousticTMM(incidence='Normal',air_temperature = 20)

    layer1 = structure.Add_JCA_Layer(30, 46182,.917,2.1,83,128)

    s1 = structure.assemble_structure(layer1)

    A1 = structure.absorption(s1)

    no_gap = os.path.join(no_gap_file,'no_gap.csv')
    structure.to_csv(no_gap,A1)


    inv = AcousticID(mount_type='No Gap',no_gap_file=no_gap_file.join("no_gap.csv"), air_temperature=20)

    # Run ML with verbose=False to reduce test output noise
    res = inv.ML(thickness=30, verbose=False)
    stats = inv.stats(res)

    assert stats['r_value'] > .99


class TestAcousticIDValidation:
    """Tests for parameter validation and error handling."""
    
    def test_inverse_uncertainty_bounds(self, tmpdir):
        """Test that uncertainty values outside [0, 1] are handled."""
        gap_file = tmpdir.mkdir("sub")
        no_gap_file = tmpdir.mkdir("sub2")

        structure = AcousticTMM(incidence='Normal', air_temperature=20)
        layer1 = structure.Add_JCA_Layer(30, 46182, .917, 2.1, 83, 128)
        air = structure.Add_Air_Layer(thickness=100)

        s1 = structure.assemble_structure(layer1)
        A1 = structure.reflection(s1)
        no_gap = os.path.join(no_gap_file, 'no_gap.csv')
        structure.to_csv(no_gap, A1)

        s2 = structure.assemble_structure(layer1, air)
        A2 = structure.reflection(s2)
        gap = os.path.join(gap_file, 'gap.csv')
        structure.to_csv(gap, A2)

        inv = AcousticID(
            mount_type='Dual',
            no_gap_file=no_gap_file.join("no_gap.csv"),
            gap_file=gap_file.join("gap.csv"),
            input_type='reflection',
            air_temperature=20
        )

        # Test with uncertainty > 1 (should revert to default)
        res = inv.Inverse(
            thickness=30,
            flow_resistivity=46182,
            porosity=.917,
            air_gap=100,
            uncertainty=1.5  # Invalid - should revert to 0.01
        )
        # Should still produce valid results
        assert 'tortuosity' in res
        assert res['tortuosity'] > 0
    
    def test_indirect_requires_dual_mount(self, tmpdir):
        """Test that Indirect method requires Dual mount type."""
        no_gap_file = tmpdir.mkdir("sub")
        structure = AcousticTMM(incidence='Normal', air_temperature=20)
        layer1 = structure.Add_JCA_Layer(30, 46182, .917, 2.1, 83, 128)
        
        s1 = structure.assemble_structure(layer1)
        A1 = structure.reflection(s1)
        no_gap = os.path.join(no_gap_file, 'no_gap.csv')
        structure.to_csv(no_gap, A1)

        inv = AcousticID(
            mount_type='No Gap',
            no_gap_file=no_gap_file.join("no_gap.csv"),
            input_type='reflection',
            air_temperature=20
        )

        with pytest.raises(ValueError, match="Dual mount types are required"):
            inv.Indirect(thickness=30, porosity=.917, air_gap=100)
    
    def test_indirect_requires_reflection_or_surface(self, tmpdir):
        """Test that Indirect method rejects absorption data."""
        gap_file = tmpdir.mkdir("sub")
        no_gap_file = tmpdir.mkdir("sub2")

        structure = AcousticTMM(incidence='Normal', air_temperature=20)
        layer1 = structure.Add_JCA_Layer(30, 46182, .917, 2.1, 83, 128)
        air = structure.Add_Air_Layer(thickness=100)

        s1 = structure.assemble_structure(layer1)
        A1 = structure.absorption(s1)  # Using absorption instead of reflection
        no_gap = os.path.join(no_gap_file, 'no_gap.csv')
        structure.to_csv(no_gap, A1)

        s2 = structure.assemble_structure(layer1, air)
        A2 = structure.absorption(s2)
        gap = os.path.join(gap_file, 'gap.csv')
        structure.to_csv(gap, A2)

        inv = AcousticID(
            mount_type='Dual',
            no_gap_file=no_gap_file.join("no_gap.csv"),
            gap_file=gap_file.join("gap.csv"),
            input_type='absorption',  # This should fail
            air_temperature=20
        )

        with pytest.raises(ValueError, match="Absorption data cannot be used"):
            inv.Indirect(thickness=30, porosity=.917, air_gap=100)


class TestAcousticIDResults:
    """Tests for result dictionary structure and values."""
    
    def test_inverse_result_keys(self, tmpdir):
        """Test that Inverse returns all expected keys."""
        gap_file = tmpdir.mkdir("sub")
        no_gap_file = tmpdir.mkdir("sub2")

        structure = AcousticTMM(incidence='Normal', air_temperature=20)
        layer1 = structure.Add_JCA_Layer(30, 46182, .917, 2.1, 83, 128)
        air = structure.Add_Air_Layer(thickness=100)

        s1 = structure.assemble_structure(layer1)
        A1 = structure.reflection(s1)
        no_gap = os.path.join(no_gap_file, 'no_gap.csv')
        structure.to_csv(no_gap, A1)

        s2 = structure.assemble_structure(layer1, air)
        A2 = structure.reflection(s2)
        gap = os.path.join(gap_file, 'gap.csv')
        structure.to_csv(gap, A2)

        inv = AcousticID(
            mount_type='Dual',
            no_gap_file=no_gap_file.join("no_gap.csv"),
            gap_file=gap_file.join("gap.csv"),
            input_type='reflection',
            air_temperature=20
        )

        res = inv.Inverse(
            thickness=30,
            flow_resistivity=46182,
            porosity=.917,
            air_gap=100,
            uncertainty=.10
        )

        expected_keys = [
            'thickness', 'flow resistivity', 'porosity', 'tortuosity',
            'viscous characteristic length', 'thermal characteristic length',
            'air gap', 'error'
        ]
        for key in expected_keys:
            assert key in res, f"Missing key: {key}"
    
    def test_physical_constraints(self, tmpdir):
        """Test that identified parameters are physically reasonable."""
        gap_file = tmpdir.mkdir("sub")
        no_gap_file = tmpdir.mkdir("sub2")

        structure = AcousticTMM(incidence='Normal', air_temperature=20)
        layer1 = structure.Add_JCA_Layer(30, 46182, .917, 2.1, 83, 128)
        air = structure.Add_Air_Layer(thickness=100)

        s1 = structure.assemble_structure(layer1)
        A1 = structure.reflection(s1)
        no_gap = os.path.join(no_gap_file, 'no_gap.csv')
        structure.to_csv(no_gap, A1)

        s2 = structure.assemble_structure(layer1, air)
        A2 = structure.reflection(s2)
        gap = os.path.join(gap_file, 'gap.csv')
        structure.to_csv(gap, A2)

        inv = AcousticID(
            mount_type='Dual',
            no_gap_file=no_gap_file.join("no_gap.csv"),
            gap_file=gap_file.join("gap.csv"),
            input_type='reflection',
            air_temperature=20
        )

        res = inv.Inverse(
            thickness=30,
            flow_resistivity=46182,
            porosity=.917,
            air_gap=100,
            uncertainty=.10
        )

        # Physical constraints
        assert res['tortuosity'] >= 1.0, "Tortuosity must be >= 1"
        assert res['porosity'] > 0 and res['porosity'] <= 1, "Porosity must be in (0, 1]"
        assert res['viscous characteristic length'] <= res['thermal characteristic length'], \
            "VCL should be <= TCL for most materials"