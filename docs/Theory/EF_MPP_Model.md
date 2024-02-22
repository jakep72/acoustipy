# Equivalent Fluid MPP Model

The equivalent fluid microperforate model estimates the input parameters to the [JCA](https://jakep72.github.io/acoustipy/Theory/JCA_Model/) model based on the perforate diameter $(d)$, center-to-center distance $(b)$, and thickness $(t)$.

The internal [_calc_dynamics](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM._calc_dynamics) method is then used to determine the dynamic mass density $(\tilde{\rho})$ and dynamic bulk modulus $(\widetilde{K})$ using the [JCA](https://jakep72.github.io/acoustipy/Theory/JCA_Model/) model.

# Estimation of JCA parameters

\[
\phi = \frac{\pi}{4}\Bigg(\frac{d}{b}\Bigg)^2
\]

\[
\sigma = \frac{32\eta}{\phi d^2}
\]

\[
\Lambda = \frac{d}{2}
\]

\[
\Lambda^{\prime} = \frac{d}{2}
\]

\[
\tau = 1+\frac{2*fok}{t}
\]

where $fok$ is:

\[
fok = \frac{4d}{3\pi}(1-1.13eps-0.09eps^2+0.27eps^3)
\]

and $eps$ is:

\[
eps = 2\sqrt{\frac{\phi}{\pi}}
\]

The [Add_MPP_EF_Layer](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM.Add_MPP_EF_Layer) method then converts The modified dynamic mass density and bulk modulus to the characteristic impedence $(Z_{c})$ and wavenumber $(k_{c})$ for use in the layer transfer matrix.

\[
Z_{c} = \sqrt{\tilde{\rho}\widetilde{K}}
\]

\[
k_{c} = {\omega}\sqrt{\frac{\tilde{\rho}}{\widetilde{K}}}
\]

## Model Parameters:

##### Using the following nomenclature --- Symbol = [Units] (name)

###### Model Specific
\[
d = \Bigg[m\Bigg]\tag{perforate diameter}
\]

\[
b = \Bigg[m\Bigg]\tag{center-to-center distance}
\]

\[
t = \Bigg[m\Bigg]\tag{layer thickness}
\]

###### JCA Parameters

\[
\sigma = \Bigg[\frac{Pa*s}{m^2}\Bigg]\tag{static airflow resistivity}
\]

\[
\phi = \Bigg[unitless\Bigg]\tag{porosity}
\]

\[
\tau = \Bigg[unitless\Bigg]\tag{tortuosity}
\]

\[
\Lambda = \Bigg[{\mu}m\Bigg]\tag{viscous characteristic length}
\]

\[
\Lambda^{\prime} = \Bigg[{\mu}m\Bigg]\tag{thermal characteristic length}
\]



## Defining Other Symbols:

##### Using the following nomenclature --- Symbol = [Units] (name)


\[
\eta = \Bigg[Pa*s\Bigg]\tag{viscosity of air}
\]