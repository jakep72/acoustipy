# Using impedance tube data
Estimate the diffuse sound field absorption coefficients of a structure
using the normal incidence reflection coefficients obtained via impedance
tube. 

```python
from acoustipy import AcousticTMM

# Generate synthetic "impedance tube" data
structure = AcousticTMM(incidence='Normal', air_temperature = 20)

# Define the JCA and air gap material parameters for the synthetic data
layer1 = structure.Add_JCA_Layer(25.4, 60000, .91, 1.3, 85, 110)
air = structure.Add_Air_Layer(thickness = 100)

# Generate rigid backed reflection data and save to a csv file
s1 = structure.assemble_structure(layer1)
r1 = structure.reflection(s1)
structure.to_csv('no_gap_tube',r1)

# Generate air backed absorption data and save to a csv file
s2 = structure.assemble_structure(layer1,air)
r2 = structure.reflection(s2)
structure.to_csv('gap_tube',r2)

# Create a new structure to estimate the diffuse field absorption coefficients
estimated_structure = AcousticTMM(incidence='Diffuse', air_temperature = 20)

# Load the synthetic impedance tube data and define a new air layer
estimated_layer = estimated_structure.Add_Layer_From_Tube(no_gap_file = 'no_gap_tube.csv', gap_file='gap_tube.csv', sample_thickness=25.4, air_gap_thickness=100)
air2 = estimated_structure.Add_Air_Layer(thickness = 400)

# Generate narrow band absorption data for the new structure
s3 = estimated_structure.assemble_structure(estimated_layer,air2)
absorption = structure.absorption(s3)

# Calculate the 3rd octave bands absorption coefficients
bands = structure.octave_bands(absorption)

# Calculate the four frequency average absorption
ffa_estimate = structure.FFA(bands)
print(ffa_estimate)

>>> 0.717

# Compare the estimated FFA to that calculated from the original layer data
air3 = structure.Add_Air_Layer(thickness = 400)
s4 = structure.assemble_structure(layer1,air3)
abs1 = structure.absorption(s4)
bands1 = structure.octave_bands(abs1)
ffa_original = structure.FFA(bands1)
print(ffa_original)

>>> 0.717
```














