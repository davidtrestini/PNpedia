# Angular momentum for nonspinning compact binaries on elliptic orbits

As pointed out in [arXiv:2504.13245v1](https://arxiv.org/abs/2504.13245v1), there are two distinct notions of angular momentum. The conservative angular momentum $`\mathbf{J}_\mathrm{cons}`$ is a constant of motion under the conservative equations of motion. The binding angular momentum $`\mathbf{J}`$ enters the flux balance law $`d \mathbf{J}/d t = - \mathcal{F}_{\mathbf{J}}`$. The difference between the two is called a Schott term, $`\mathbf{J}_\mathrm{Schott}`$, which is nonvanishing for circular orbits starting at 4PN. The expression of the binding angular momentum in terms of the orbital (radial and azimuthal) frequencies depends on an arbitrary scale $b_0$, which is related to the choice of slicing in relating near-zone and far zone quantities. This arbitrary constant also appears in the flux, and drops out of the balance law. For nonprecessing binaries, only the direction of the angular momentum is constant and orthogonal to the orbital plane: only the norm $`J = |\mathbf{J}|`$ is needed.

For nonspinning compact binaries on circular orbits:
* the file ``angular_momentum_conservative.txt`` contains the conservative angular momentum in terms of the Blanchet orbital parameters $(x,\iota)$, which are directly related to the radial and azimuthal frequencies (see below)
* the file  ``lambda_0_expansion.txt`` contains the small $e$ expansion  of the function $\lambda_0(e)$ up to $\mathcal{O}(e^8)$

The conservative angular momentum is given at 4PN accuracy. The Schott term, and thus the binding angular momentum, are not known at 4PN. At 3PN, the Schott term is vanishing and the binding angular momentum is thus identical to the conservative angular momentum. 

## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless Blanchet parameter $x = G m \omega /c^3$, where $\omega$ is the dimensionful azimuthal frequency
* ``\[Iota]`` is the dimensionless Blanchet parameter $\iota = \frac{3 x}{\omega/n-1}$, where $n$ is the dimensionful radial frequency
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{(m_1 + m_2)^2}$
* ``b_0`` is an arbitary constant linked to the choice of foliation

We also introduce the enhancement function

$$\Lambda(e) = \frac{1}{16}\sum_{p=1}^{\infty} p^6 \ln\left(\frac{p}{2}\right) ({}_{p} \hat{\mathrm{I}}_{ij}) ({}_{-p} \hat{\mathrm{I}}_{ij}) ,$$


where $e\in [0,1)$ and we use the Einstein summation on the spatial indices $i,j \in \{x,y,z\}$. The Fourier coefficients are explicitly expressed in terms of the Bessel functions $J_p$ and their derivatives $J_p'$ as follows:

$$\begin{aligned}
{}_{p}\hat{\mathrm{I}}_{xx} &= - \frac{2}{3} \frac{3-e^2}{e}\frac{J_p(pe)}{p^2} + \frac{2(1-e^2)}{e} \frac{J_p'(e)}{p}\\
{}_{p}\hat{\mathrm{I}}_{xy} &= {}_{p} \hat{\mathrm{I}}_{yx}  = 2 \mathrm{i} \sqrt{1-e^2} \left[- \frac{1-e^2}{e^2} \frac{J_{p}(p e)}{p} + \frac{1}{e} \frac{J_p'(p e)}{p ^2}\right] \\
{}_{p}\hat{\mathrm{I}}_{yy} &=  \frac{2}{3} \frac{3-2e^2}{e^2}\frac{J_p(e)}{p^2} -  \frac{2(1-e^2)}{e} \frac{J_p'(p e)}{p} \\
{}_{p}\hat{\mathrm{I}}_{zz} &= \frac{2}{3} \frac{J_p(p e)}{p^2}\\
{}_{p}\hat{\mathrm{I}}_{xz} &=  {}_{p} \hat{\mathrm{I}}_{yz} = {}_{p} \hat{\mathrm{I}}_{zx} = {}_{p} \hat{\mathrm{I}}_{zy} = 0
\end{aligned}$$


We then define function $\lambda_0(e)$ such that 

$$\Lambda_0(e) = - \frac{3}{2(1-e^2)^{7/2}}\left[\ln(1-e^2)\left(1 + \frac{73}{24}e^2 + \frac{37}{96}e^4\right) + e^2 \lambda_0(e)\right] .$$

The function $\lambda_0(e)$ is not known in closed form, but replacing it by its small $e$ expansion leads to very accurate results; see Sec. IV.D of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1). The small $e$ expansion of $\lambda_0(e)$ is provided up to $\mathcal{O}(e^8)$ in ``lambda_0_expansion.txt``.

In the files, these functions are evaluated for $e=\sqrt{1-\iota}$. We denote $\lambda_0(\sqrt{1-\iota})$ as ``\[Lambda]0[Sqrt[1 - \[Iota]]``, and its derivative $\lambda_0'(\sqrt{1-\iota})$ is denoted ``Derivative[1][\[Lambda]0][Sqrt[1 - \[Iota]]``.

The conservative angular momentum was obtained:
* at 4PN in
    * (5.11)-(5.12) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1)

The small $e$ expansion of $\lambda_0(e)$ was obtained
* up to $\mathcal{O}(e^8)$ in
    * (4.41) and (4.43) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1)

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
