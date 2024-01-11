# Inverse JCA Characterization
Use an inverse characterization method to identify parameters of the Johnson-Champaux-Allard model given data obtained from an impedance tube.

```python
from acoustipy import acousticTMM, AcousticID

# Create an AcousticTMM object to generate toy impedance tube data
structure = acousticTMM(incidence='Normal',air_temperature = 20)

# Define the JCA and air gap material parameters for the toy data
layer1 = structure.Add_JCA_Layer(thickness = 30, flow_resistivity = 46879, porosity = .93, tortuosity = 1.7, viscous_characteristic_length = 80, thermal_characteristic_length = 105)
air = structure.Add_Air_Layer(thickness = 375)

# Generate rigid backed absorption data and save to a csv file
s1 = structure.assemble_structure(layer1)
A1 = structure.absorption(s1)
structure.to_csv('no_gap',A1)

# Generate air backed absorption data and save to a csv file
s2 = structure.assemble_structure(layer1,air)
A2 = structure.absorption(s2)
structure.to_csv('gap',A2)

# Create an AcousticID object, specifying to mount types, data files, and data types
inv = AcousticID(mount_type='Dual',no_gap_file="no_gap.csv", gap_file = 'gap.csv',air_temperature=20,input_type='absorption')

# Call the Inverse method to find the tortuosity, viscous, and thermal characteristic lengths of the material
res = inv.Inverse(30, 47000,.926,air_gap=375,uncertainty=.2,verbose=True)

# Display summary statistics about the optimization
stats = inv.stats(res)
print(stats)

# Plot the results of the found parameters compared to the toy input data
inv.plot_comparison(res)

# Save the optimization results to a csv
inv.to_csv("params.csv",res)

stats = {'slope': 1.000037058594857, 'intercept': 9.276088883464206e-05, 'r_value': 0.9999999674493408, 'p_value': 0.0, 'std_err': 8.732362148426126e-06}
```

<br></br>

![](../assets/ex_material_identification_inverse.png)