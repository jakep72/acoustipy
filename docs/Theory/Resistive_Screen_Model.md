The Resistive Screen model implemented in acoustipy is from the paper below, neglecting the frame mechanical properties:

Mathieu Gaborit, Olivier Dazel, Peter Göransson; A simplified model for thin acoustic screens. J. Acoust. Soc. Am. 1 July 2018; 144 (1): EL76–EL81. https://doi.org/10.1121/1.5047929

The equations for the dynamic mass density $(\tilde{\rho})$ and dynamic bulk modulus $(\widetilde{K})$ can be found below:

\[
\tilde{\rho} = \frac{\rho_{0}}{\phi}+\frac{\sigma}{j{\omega}}
\]

\[
\widetilde{K} = \frac{P_{0}}{\phi}
\]

The [Add_Resistive_Screen](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM.Add_Resistive_Screen) method then converts The dynamic mass density and bulk modulus to the characteristic impedence $(Z_{c})$ and wavenumber $(k_{c})$ for use in the layer transfer matrix.

\[
Z_{c} = \sqrt{\tilde{\rho}\widetilde{K}}
\]

\[
k_{c} = {\omega}\sqrt{\frac{\tilde{\rho}}{\widetilde{K}}}
\]

## Model Parameters:

##### Using the following nomenclature --- Symbol = [Units] (name)

\[
\sigma = \Bigg[\frac{Pa*s}{m^2}\Bigg]\tag{static airflow resistivity}
\]

\[
\phi = \Bigg[unitless\Bigg]\tag{porosity}
\]
## Defining Other Symbols:

##### Using the following nomenclature --- Symbol = [Units] (name) or Symbol = equation = [Units] (name)

\[
  \rho_{0} = \Bigg[\frac{kg}{m^3}\Bigg]\tag{air density}
\]

\[
\omega = 2{\pi}f = \Bigg[\frac{radians}{s}\Bigg]\tag{angular frequency}
\]