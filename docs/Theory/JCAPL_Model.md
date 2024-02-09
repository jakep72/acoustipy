The Johnson-Champoux-Allard-Pride-Lafarge model is an eight parameter model consisting of the static airflow resistivity $(\sigma)$, porosity $(\phi)$, tortuosity $(\tau)$, viscous characteristic length $(\Lambda)$, thermal characteristic length $(\Lambda^\prime)$, thermal permeability $(k_{0}^\prime)$, thermal tortuosity $(\alpha_{0}^\prime)$, and viscous tortuosity $(\alpha_{0})$.

The acoustipy implementation for the JCA, JCAL, and JCAPL models are all based on the implementation from [APMR](https://apmr.matelys.com/PropagationModels/MotionlessSkeleton/JohnsonChampouxAllardPrideLafargeModel.html).  The equations described below can be found in the [_calc_dynamics](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM._calc_dynamics) method.

# Dynamic Mass Density
\[
\tilde{\rho} = \frac{\rho_{0}\tilde{\alpha}(\omega)}{\phi}
\]

\[
\tilde{\alpha}(\omega) = \tau\Bigg[1+\frac{\tilde{F}(\omega)}{j\bar{\omega}}\Bigg]
\]

\[
\tilde{F}(\omega) = 1-P+P\sqrt{1+\frac{M}{2P^2}j\bar{\omega}}
\]

\[
\bar{\omega} = \frac{{\omega}{\rho_{0}}{\tau}}{{\sigma}{\phi}}
\]

\[
M = \frac{8{\eta}{\tau}}{{\sigma}{\phi}{\Lambda}^2}
\]

\[
P = \frac{M}{4\Bigg(\frac{\alpha_{0}}{\tau}-1\Bigg)}
\]

# Dynamic Bulk Modulus
\[
\widetilde{K} = \frac{{\gamma}P_{0}}{{\phi}\tilde{\beta}(\omega)}
\]

\[
\tilde{\beta}(\omega) = \gamma-(\gamma - 1)\Bigg[1+\frac{\tilde{F^{\prime}}(\omega)}{j\bar{\omega^\prime}}\Bigg]^{-1}
\]

\[
\tilde{F^{\prime}}(\omega) = 1-P^{\prime}+P^{\prime}\sqrt{1+\frac{M^{\prime}}{2P^{\prime{2}}}j\bar{\omega^{\prime}}}
\]

\[
\bar{\omega^{\prime}} = \frac{{\omega}{\rho_{0}}{Pr}{k_{0}^\prime}}{{\eta}{\phi}}
\]

\[
M^{\prime} = \frac{8{k_{0}^\prime}}{{\phi}{\Lambda^{\prime{2}}}}
\]

\[
P^{\prime} = \frac{M^\prime}{4\Bigg(\alpha_{0}^\prime-1\Bigg)}
\]

The [Add_JCAPL_Layer](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM.Add_JCAPL_Layer) method then converts The dynamic mass density and bulk modulus to the characteristic impedence $(Z_{c})$ and wavenumber $(k_{c})$ for use in the layer transfer matrix.

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

\[
\tau = \Bigg[unitless\Bigg]\tag{tortuosity}
\]

\[
\Lambda = \Bigg[{\mu}m\Bigg]\tag{viscous characteristic length}
\]

\[
\Lambda^{\prime} = \Bigg[{\mu}m\Bigg]\tag{thermal characteristic length}
\]

\[
k_{0}^{\prime} = \Bigg[m^2\Bigg]\tag{thermal permeability}
\]

\[
\alpha_{0}^{\prime} = \Bigg[unitless\Bigg]\tag{thermal tortuosity}
\]

\[
\alpha_{0} = \Bigg[unitless\Bigg]\tag{viscous tortuosity}
\]



## Defining Other Symbols:

##### Using the following nomenclature --- Symbol = [Units] (name)

\[
\rho_{0} = \Bigg[\frac{kg}{m^3}\Bigg]\tag{air density}
\]

\[
\omega = \Bigg[\frac{radians}{s}\Bigg]\tag{angular frequency}
\]

\[
\eta = \Bigg[Pa*s\Bigg]\tag{viscosity of air}
\]

\[
\Pr = \Bigg[unitless\Bigg]\tag{Prandtl Number}
\]