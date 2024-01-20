The Delaney-Bazley-Mike model is an empirical, one parameter parameter model reliant on the static airflow resistivity $(\sigma)$ of the porous material.  It is generally considered a better model than Delaney-Bazley.

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

The model is valid in the frequency range defined below:

\[
0.01 < \displaystyle{\frac{f}{\sigma}} < 1.00
\]

## Parameter Units:

\[
    \sigma = \Bigg[\frac{Pa*s}{m^2}\Bigg]
\]