from src.acoustipy.TMM import AcousticTMM
from src.acoustipy.Params import AcousticID

# Create an AcousticTMM object to generate toy impedance tube data
# structure = AcousticTMM(incidence='Normal',air_temperature = 20)

# # Define the JCA and air gap material parameters for the toy data
# layer1 = structure.Add_JCA_Layer(thickness = 30, flow_resistivity = 25000, porosity = .93, tortuosity = 1.7, viscous_characteristic_length = 80, thermal_characteristic_length = 105)
# air = structure.Add_Air_Layer(thickness = 375)

# # Generate rigid backed absorption data and save to a csv file
# s1 = structure.assemble_structure(layer1)
# A1 = structure.absorption(s1)
# structure.to_csv('no_gap',A1)

# # Generate air backed absorption data and save to a csv file
# s2 = structure.assemble_structure(layer1,air)
# A2 = structure.absorption(s2)
# structure.to_csv('gap',A2)

# # Create an AcousticID object, specifying to mount types, data files, and data types
# inv = AcousticID(mount_type='Dual',no_gap_file="no_gap.csv", gap_file = 'gap.csv',air_temperature=20,input_type='absorption')

# # Call the Inverse method to find the tortuosity, viscous, and thermal characteristic lengths of the material
# res = inv.Inverse(30, 26000,.93,air_gap=375,uncertainty=.05,verbose=True)
# print(res)
# # Display summary statistics about the optimization
# stats = inv.stats(res)
# print(stats)

# Plot the results of the found parameters compared to the toy input data
# inv.plot_comparison(res)

# # Save the optimization results to a csv
# inv.to_csv("params.csv",res)

structure = AcousticTMM(incidence='Normal',air_temperature = 20)

layer1 = structure.Add_JCA_Layer(30, 46182,.917,2.1,83,128)

air = structure.Add_Air_Layer(thickness = 100)

s1 = structure.assemble_structure(layer1)

A1 = structure.reflection(s1)

# no_gap = os.path.join(no_gap_file,'no_gap.csv')
structure.to_csv('no_gap',A1)

s2 = structure.assemble_structure(layer1,air)

A2 = structure.reflection(s2)

# gap = os.path.join(gap_file,'gap.csv')
structure.to_csv('gap',A2)

inv = AcousticID(mount_type='Dual',no_gap_file="no_gap.csv",gap_file ="gap.csv",input_type='reflection',air_temperature=20)

res = inv.Hybrid(thickness=30,porosity=.917,air_gap=100,uncertainty = .10)
print(res)
stats = inv.stats(res)
print(stats)
inv.plot_comparison(res)