# MAA MPP Model

The Maa microperforate model determines the characteristic impedence $(Z_{c})$ of the layer based on the perforate diameter $(d)$, center-to-center distance $(b)$, and thickness $(t)$.

Unlike other material models, the Maa MPP model does not calculate a characteristic wavenumber $(k_{c})$.  Instead, a separate transfer matrix is used.

# Determination of Characteristic Impedence

\[
Z_{c} = r+j\omega m
\]

where $r$ is:

\[
r = \frac{32\eta t}{\phi d^2} r_{1}
\]

\[
r_{1} = \sqrt{1+\frac{x^2}{32}}+\frac{\sqrt{2}}{32} x \frac{d}{t}
\]

$m$ is:

\[
m = \frac{\rho_{0}t}{\phi} m_{1}
\]

\[
m_{1} = 1+\frac{1}{\sqrt{1+\frac{x^2}{2}}}+\frac{0.85d}{t}
\]

and 

\[
x = \frac{d}{2}\sqrt{\frac{\omega\rho_{0}}{\eta}}
\]

\[
\phi = \frac{\pi}{4}\Bigg(\frac{d}{b}\Bigg)^2
\]

The characteristic impedence $(Z_{c})$ is then used directly in the layer transfer matrix via the [Add_MAA_MPP_Layer](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM.Add_MAA_MPP_Layer) method.


## Model Parameters:

##### Using the following nomenclature --- Symbol = [Units] (name)

\[
d = \Bigg[m\Bigg]\tag{perforate diameter}
\]

\[
b = \Bigg[m\Bigg]\tag{center-to-center distance}
\]

\[
t = \Bigg[m\Bigg]\tag{layer thickness}
\]


## Defining Other Symbols:

##### Using the following nomenclature --- Symbol = [Units] (name)

\[
\rho_{0} = \Bigg[\frac{kg}{m^3}\Bigg]\tag{air density}
\]

\[
\eta = \Bigg[Pa*s\Bigg]\tag{viscosity of air}
\]

\[
\omega = \Bigg[\frac{radians}{s}\Bigg]\tag{angular frequency}
\]