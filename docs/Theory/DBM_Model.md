The Delaney-Bazley-Mike model is an empirical, one parameter  model consisting of the static airflow resistivity $(\sigma)$ of the porous material.  It is generally considered a better fit of the data from Delaney-Bazley.

The equations for the complex characteristic impedance $(Z_{c})$ and wavenumber $(k_{c})$ can be found below:

\[
Z_{c} = \rho_{0} c_{0} 
  \Bigg[ 1+0.0699  \left(\frac{f}{\sigma}\right)^{-0.632} 
        - j 0.107  \left(\frac{f}{\sigma}\right)^{-0.632} \Bigg]
\]

\[
k_{c} = \displaystyle\frac{\omega}{c_{0}} 
  \Bigg[ 1 + 0.109  \left(\frac{f}{\sigma}\right)^{-0.618} 
         - j 0.160  \left(\frac{f}{\sigma}\right)^{-0.618} \Bigg]
\]

As with the Delaney-Bazley Model, the acoustipy implementation converts the characteristic impedence and wavenumber to the dynamic mass density $(\tilde{\rho})$ and dynamic bulk modulus $(\widetilde{K})$ via the equations below, as the [Add_DBM_Layer](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM.Add_DBM_Layer) method calls the internal [_calc_dynamics](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM._calc_dynamics) method for consistency with the other equivalent fluid models.

\[
\tilde{\rho} = \frac{Z_{c}k_{c}}{\omega}
\]

\[
\widetilde{K} = \frac{{\omega}Z_{c}}{k_{c}}
\]

The dynamic mass density and bulk modulus are then converted back to the characteristic impedence and wavenumber for use in the layer transfer matrix via:

\[
Z_{c} = \sqrt{\tilde{\rho}\widetilde{K}}
\]

\[
k_{c} = {\omega}\sqrt{\frac{\tilde{\rho}}{\widetilde{K}}}
\]

The model is valid in the frequency range defined below:

\[
0.01 < \displaystyle{\frac{f}{\sigma}} < 1.00
\]

## Model Parameters:

##### Using the following nomenclature --- Symbol = [Units] (name)

\[
\sigma = \Bigg[\frac{Pa*s}{m^2}\Bigg]\tag{static airflow resistivity}
\]

## Defining Other Symbols:

##### Using the following nomenclature --- Symbol = [Units] (name) or Symbol = equation = [Units] (name)

\[
\rho_{0} = \Bigg[\frac{kg}{m^3}\Bigg]\tag{air density}
\]

\[
f = \Bigg[Hz\Bigg]\tag{linear frequency}
\]

\[
\omega = 2{\pi}f = \Bigg[\frac{radians}{s}\Bigg]\tag{angular frequency}
\]