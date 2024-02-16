The Biot-Rigid model modifies the dynamic mass density $(\tilde{\rho}_{eq})$ determined via an equivalent fluid model such as: [DB](https://jakep72.github.io/acoustipy/Theory/DB_Model/), [DBM](https://jakep72.github.io/acoustipy/Theory/DBM_Model/), [JCA](https://jakep72.github.io/acoustipy/Theory/JCA_Model/), [JCAL](https://jakep72.github.io/acoustipy/Theory/JCAL_Model/), and [JCAPL](https://jakep72.github.io/acoustipy/Theory/JCAPL_Model/).

The acoustipy implementation follows from eq. 23 in [Bécot and Jaouen](https://doi.org/10.1121/1.4826175).

\[
\frac{1}{\tilde{\rho}_{rigid}} = \frac{1}{\phi\tilde{\rho}_{eq}}+\frac{\gamma^2}{\phi\tilde{\rho}}+\frac{(1-\phi)}{\phi}\frac{\gamma}{\tilde{\rho}}
\]

where $\gamma$ and $\tilde{\rho}$ are defined as:

\[
\gamma = \frac{\rho_{0}}{\tilde{\rho}_{eq}}-1
\]

\[
\tilde{\rho} = \rho_{1}+\phi\rho_{0}-\frac{\rho_{0}^2}{\tilde{\rho}_{eq}}
\]

The [Add_Biot_Rigid_Layer](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM.Add_Biot_Rigid_Layer) method then converts The modified dynamic mass density and bulk modulus to the characteristic impedence $(Z_{c})$ and wavenumber $(k_{c})$ for use in the layer transfer matrix.

\[
Z_{c} = \sqrt{\tilde{\rho}_{rigid}\widetilde{K}}
\]

\[
k_{c} = {\omega}\sqrt{\frac{\tilde{\rho}_{rigid}}{\widetilde{K}}}
\]


## Model Parameters:

##### Using the following nomenclature --- Symbol = [Units] (name)

\[
\rho_{1} = \Bigg[\frac{kg}{m^3}\Bigg]\tag{volumetric density}
\]

The following parameters are used to find the dynamic mass density $(\tilde{\rho}_{eq})$ from the specified equivalent fluid model.

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

