The Horoshenkov model consists of three parameters -- the porosity $(\phi)$, median pore size $(\bar s)$, and pore size distribution $(\sigma_{s})$.

From these three parameters, the rest of the JCAL parameters -- the static airflow resistivity $(\sigma)$, porosity $(\phi)$, tortuosity $(\tau)$, viscous characteristic length $(\Lambda)$, thermal characteristic length $(\Lambda^\prime)$, and thermal permeability $(k_{0}^\prime)$ -- can be calculated as seen below:

\[
\tau = \exp\Bigg(4[\sigma_{s}ln(2)]^2\Bigg)
\]

\[
\sigma = \frac{8\eta\tau}{\phi{\bar s^2}}\exp\Bigg(6[\sigma_{s}ln(2)]^2\Bigg)
\]

\[
\Lambda = \bar s \exp\Bigg(-\frac{5}{2}[\sigma_{s}ln(2)]^2\Bigg)
\]

\[
\Lambda^\prime = \bar s \exp\Bigg(\frac{3}{2}[\sigma_{s}ln(2)]^2\Bigg)
\]

\[
k_{0}^\prime = \frac{\phi{\bar s^2}}{8\tau}\exp\Bigg(-6[\sigma_{s}ln(2)]^2\Bigg)
\]

The dynamic mass density and bulk modulus are then determined within the [_calc_dynamics](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM._calc_dynamics) method using the calculated JCAL parameters, following the implemention found here: [JCAL](https://jakep72.github.io/acoustipy/Theory/JCAL_Model/).

The [Add_Horoshenkov_Layer](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM.Add_Horoshenkov_Layer) method then converts The dynamic mass density and bulk modulus to the characteristic impedence $(Z_{c})$ and wavenumber $(k_{c})$ for use in the layer transfer matrix.

\[
Z_{c} = \sqrt{\tilde{\rho}\widetilde{K}}
\]

\[
k_{c} = {\omega}\sqrt{\frac{\tilde{\rho}}{\widetilde{K}}}
\]

## Model Parameters:

##### Using the following nomenclature --- Symbol = [Units] (name)


\[
\phi = \Bigg[unitless\Bigg]\tag{porosity}
\]

\[
\bar s = \Bigg[{\mu}m\Bigg]\tag{median pore size}
\]

\[
\sigma_{s} = \Bigg[unitless\Bigg]\tag{pore size distribution}
\]




## Defining Other Symbols:

##### Using the following nomenclature --- Symbol = [Units] (name)

\[
\eta = \Bigg[Pa*s\Bigg]\tag{viscosity of air}
\]