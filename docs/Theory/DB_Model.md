The Delaney-Bazley model is an empirical, one parameter parameter model reliant on the static airflow resistivity $(\sigma)$ of the porous material.

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

The model is valid in the frequency range defined below:

\[
0.01 < \displaystyle{\frac{f}{\sigma}} < 1.00
\]

## Parameter Units:

\[
    \sigma = \Bigg[\frac{Pa*s}{m^2}\Bigg]
\]