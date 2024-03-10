The frequency-dependent transmission coefficients $(\tau)$ are calculated directly from the total transfer matrix $(T_{t})$ of a multilayered structure.  These coefficients are a measure of how much sound passes through a structure.

The coefficients can be calculated under both normal and diffuse sound field conditions.  Under a normal incidence sound field, the sound impinges on the surface from a single, perpendicular angle.  In the diffuse field case, the incident sound theoretically strikes the surface of the material from all possible angles -- though the acoustipy implementation defaults to angles between 0 and 79, as seen in literature on the topic.

The transmission coefficients can then be used to calculate transmission loss $(TL)$, which is also a frequency dependent metric.

The acoustipy implementation for both cases can be found [here](https://jakep72.github.io/acoustipy/AcousticTMM/#src.acoustipy.TMM.AcousticTMM.transmission_loss).

### Normal Incidence

Starting from the total transfer matrix:

\[
T_{t} = 
\begin{bmatrix}
T_{11} & T_{12}\\
T_{21} & T_{22} \\
\end{bmatrix}
\]

Then, the transmission coefficients $(\tau)$ are:

\[
\tau = \frac{2e^{jk_{0}t}}{T_{11}+\frac{T_{12}}{Z_{0}}+Z_{0}T_{21}+T_{22}}  
\]

where $t$ is the total thickness of the structure and $k_{0}$ is the wavenumber:

\[
k_{0} = \frac{\omega}{c_{0}}
\]

where $\omega$ is the angular frequency and $c_{0}$ is the speed of sound in air.

and where $Z_{0}$ is the characteristic impedence of air:

\[
Z_{0} = \rho_{0} c_{0}
\]

and $\rho_{0}$ is the density of air.

Then the transmission loss is:

\[
TL = 10\log_{10}\frac{1}{|\tau|^2}
\]


### Diffuse Incidence

Under the diffuse sound field condition, the calculation of the transmission coefficients is:

\[
\tau = \frac{2e^{jk_{0}t}}{T_{11}+\frac{T_{12}\cos\theta}{Z_{0}}+\frac{Z_{0}T_{21}}{\cos\theta}+T_{22}}  
\]

which yields a vector of shape $[f, \theta]$.  To collapse this vector to shape $[f,1]$, Paris' formula is used as shown below.

\[
T = \frac{\sum \tau\cos(\theta)\sin(\theta)}{\sum \cos(\theta)\sin(\theta)}
\]

the transmission loss is then calculated in the same manner as the normal incidence condition:

\[
TL = 10\log_{10}\frac{1}{|T|^2}
\]