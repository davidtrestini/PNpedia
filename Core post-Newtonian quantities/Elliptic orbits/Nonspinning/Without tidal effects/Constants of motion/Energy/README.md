# Energy for nonspinning compact binaries on elliptic orbits

As pointed out in [arXiv:2504.13245v1](https://arxiv.org/abs/2504.13245v1), there are two distinct notions of energy. The conservative energy $E_\mathrm{cons}$ is a constant of motion under the conservative equations of motion. The binding energy $E$ enters the flux balance law $d E/d t = - F_E$. The difference between the two is called a Schott term, $E_\mathrm{Schott}$, which is nonvanishing for circular orbits starting at 4PN. The expression of the binding energy in terms of the orbital frequency depends on an arbitrary scale $b_0$, which is related to the choice of slicing in relating near-zone and far zone quantities. This arbitrary constant also appears in the flux, and drops out of the balance law.

For nonspinning compact binaries on circular orbits:
* the file ``energy_conservative.txt`` contains the conservative energy in terms of the Blanchet orbital parameters $(x,\iota)$, which are directly related to the radial and azimuthal frequencies (see below)

The conservative energy is given at 4PN accuracy. The Schott term, and thus the binding energy, are not known at 4PN. At 3PN, the Schott term is vanishing and the binding energy is thus identical to the conservative energy. 

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

We also introduce the enhancement function $\Lambda(e) = \frac{1}{16}\sum_{p=1}^{\infty} p^2 \ln\left(\frac{p}{2}\right)  \phantom{\mathrm{I}}_p \mathrm{I}_{ij}$, where $e\in [0,1)$ and we use the Einstein summation on the spatial indices $i,j \in \{x,y,z\}$. The Fourier coefficients are explicitly expressed in terms of the Bessel functions $J_p(x)$ as follows:

${\,}_{p} \hat{\mathrm{I}}_{xx} = - \frac{2}{3}\,\frac{3-e^2}{e}\,\frac{J_p(pe)}{p^2} + \frac{2(1-e^2)}{e} \frac{J_p'(e)}{p}$
${\,}_{p} \hat{\mathrm{I}}_{xy} = {\,}_{p} \hat{\mathrm{I}}_{yx} = 2 \mathrm{i} \sqrt{1-e^2} \Big[- \frac{1-e^2}{e^2} \frac{J_{p}(p e)}{p} + \frac{1}{e} \frac{J_p'(p e)}{p ^2}\Big] $
${\,}_{p} \hat{\mathrm{I}}_{yy} =  \frac{2}{3} \frac{3-2e^2}{e^2}\frac{J_p(e)}{p^2} -  \frac{2(1-e^2)}{e} \frac{J_p'(p e)}{p}$
${\,}_{p} \hat{\mathrm{I}}_{xz} =  {\,}_{p} \hat{\mathrm{I}}_{yz} = {\,}_{p} \hat{\mathrm{I}}_{zx} = {\,}_{p} \hat{\mathrm{I}}_{zy = 0}$
${\,}_{p} \hat{\mathrm{I}}_{zz} = \frac{2}{3} \frac{J_p(p e)}{p^2}$

We then define function $\lambda_0(e)$ such that 
$\Lambda_0(e) = - \frac{3}{2(1-e^2)^{7/2}}\Big[\ln(1-e^2)\left(1 + \frac{73}{24}e^2 + \frac{37}{96}e^4\right) + e^2 \lambda_0(e)\Big]$
The function $\lambda_0(e)$ is not known in closed form, but replacing it by its small-$e$ expansion leads to very accurate results; see Sec. IV.D of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1). The small-$e$ expansion of $\lambda_0(e)$ is provided up to $e^8$ in ``lambda_0_expansion.txt``.

In the files, these functions are evaluated for $e=\sqrt{1-\iota}$. We denote $\lambda_0(\sqrt{1-\iota})$ as ``\[Lambda]0[Sqrt[1 - \[Iota]]``, and its derivative $\lambda_0'(\sqrt{1-\iota})$ is denoted ``Derivative[1][\[Lambda]0][Sqrt[1 - \[Iota]]``.

The conservative energy was obtained:
* at 4PN in
    * (5.11)--(5.12) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1)

The small-$e$ expansion of $\lambda_0(e)$ was obtained
* up to $e^8$ in
    * (4.41) and (4.43) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1)

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
