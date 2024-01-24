The Delaney-Bazley model is an empirical, one parameter model consisting of the static airflow resistivity $(\sigma)$ of the porous material.

The equations for the complex characteristic impedance $(Z_{c})$ and wavenumber $(k_{c})$ can be found below:

\[
Z_{c} = \rho_{0} c_{0} 
  \Bigg[ 1+0.0571  \left(\frac{\rho_{0}f}{\sigma}\right)^{-0.754} 
        - j 0.0870  \left(\frac{\rho_{0}f}{\sigma}\right)^{-0.732} \Bigg]
\]

\[
k_{c} = \displaystyle\frac{\omega}{c_{0}} 
  \Bigg[ 1 + 0.0978  \left(\frac{\rho_{0}f}{\sigma}\right)^{-0.700} 
         - j 0.1890  \left(\frac{\rho_{0}f}{\sigma}\right)^{-0.595} \Bigg]
\]

In the acoustipy implementation, the characteristic impedence and wavenumber are converted to the dynamic mass density $(\tilde{\rho})$ and dynamic bulk modulus $(\widetilde{K})$ via the equations below, as the [Add_DB_Layer](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM.Add_DB_Layer) method calls the internal [_calc_dynamics](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM._calc_dynamics) method for consistency with the other equivalent fluid models.

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