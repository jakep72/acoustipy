The frequency-dependent reflection coefficients $(R)$ are calculated directly from the total transfer matrix $(T_{t})$ of a multilayered structure.  These coefficients are a measure of how much sound is reflected off the surface of a structure.

The coefficients can be calculated under both normal and diffuse sound field conditions.  Under a normal incidence sound field, the sound impinges on the surface from a single, perpendicular angle.  In the diffuse field case, the incident sound theoretically strikes the surface of the material from all possible angles -- though the acoustipy implementation defaults to angles between 0 and 79, as seen in literature on the topic.

The acoustipy implementation for both cases can be found [here](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM.reflection).

### Normal Incidence



Starting from the total transfer matrix:

\[
T_{t} = 
\begin{bmatrix}
T_{11} & T_{12}\\
T_{21} & T_{22} \\
\end{bmatrix}
\]

First, the surface impedence $(Z_{s})$ is calculated:

\[
Z_{s} = \frac{T_{11}}{T_{21}}  
\]

Then the reflection coeffients are:

\[
R = \frac{Z_{s}-Z_{0}}{Z_{s}-Z_{0}}  
\]

where $Z_{0}$ is the characteristic impedence of air:

\[
Z_{0} = \rho_{0} c_{0}
\]

and $\rho_{0}$ is the density of air and $c_{0}$ is the speed of sound in air.

### Diffuse Incidence

Under the diffuse sound field condition, the calculation of surface impedence $(Z_{s})$ is the same as the normal incidence condition.

The reflection coefficients at each angle are then:

\[
r = \frac{Z_{s}\cos(\theta)-Z_{0}}{Z_{s}\cos(\theta)+Z_{0}}
\]

which yields a vector of shape $[f, \theta]$.  To collapse this vector to shape $[f,1]$, Paris' formula is used as shown below.

\[
R = \frac{\sum r\cos(\theta)\sin(\theta)}{\sum \cos(\theta)\sin(\theta)}
\]