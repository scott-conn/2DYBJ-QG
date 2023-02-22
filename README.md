# 2DYBJ-QG
A pseudospectral solver for the coupled YBJ-QG system in 2D

## The Coupled YBJ-QG System
The coupled YBJ-QG system describes the evolution of a NIW field in the presence of a QG eddy field [[1]](#1) as well as the back-reaction of the waves on the evolution of the QG field [[2]](#2) [[3]](#3). The 2D version of this system, where we make a vertical plane wave assumption for the NIW field, was formulated in [[4]](#4). The model equations are:

$$\begin{align}
\partial_t\phi + J(\psi,\phi)+\frac{i\zeta}{2}\phi-\frac{i\eta}{2}\nabla^2\phi &= -\nu_\phi\nabla^4\phi + F,\\
\partial_tq +J(\psi,q) &= -\nu_q\nabla^4q,\\
q &= \nabla^2\psi + \frac{1}{f_0}\left[\frac{1}{4}\nabla^2|\phi|^2 + \frac{i}{2}J(\phi^*,\phi)\right],
\end{align}$$

where $\phi$ is the back-rotated, complexified NIW velocity, $\psi$ is the QG streamfunction, $\zeta=\nabla^2\psi$ is the QG vorticity, $\eta=N_0^2/f_0m^2$ is the dispersivity (with $N_0$ being a reference stratification and $m$ being the vertical wavenumber of the NIWs), $f_0$ is the reference Coriolis parameter, $\nu_\phi$ is the wave hyperviscosity, $F$ is a forcing term, $q$ is the QGPV and $\nu_q$ is the QG hyperviscosity.

We use dedalus to solve these equations on a periodic domain.

The test file contains code to run the two examples from [[4]](#4).

## References
<a id="1">[1]</a> 
Young, W. R., Ben Jelloul, M. (1997).
Propagation of near-inertial oscillations through a geostrophic flow
Journal of Marine Research, 55(4), 735-766.

<a id="2">[2]</a> 
Xie, J-H., and Vanneste, J. (2015).
A generalised-Lagrangian-mean model of the interactions between near-inertial waves and mean flow
Journal of Fluid Mechanics, 774, 143-169.

<a id="3">[3]</a> 
Wagner, G. L., and W. R. Young. (2016).
A three-component model for the coupled evolution of near-inertial waves, quasi-geostrophic flow and the near-inertial second harmonic
Journal of Fluid Mechanics, 802, 806-837.

<a id="4">[4]</a> 
Rocha, C., Wagner, G., & Young, W. (2018).
Stimulated generation: Extraction of energy from balanced flow by near-inertial waves
Journal of Fluid Mechanics, 847, 417-451.
