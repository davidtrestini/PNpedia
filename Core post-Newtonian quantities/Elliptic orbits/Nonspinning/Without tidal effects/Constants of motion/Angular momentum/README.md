# Angular momentum for nonspinning compact binaries on elliptic orbits

As pointed out in [arXiv:2504.13245v1](https://arxiv.org/abs/2504.13245v1), there are two distinct notions of angular momentum. The conservative angular momentum $`\mathbf{J}_\mathrm{cons}`$ is a constant of motion under the conservative equations of motion. The binding angular momentum $`\mathbf{J}`$ enters the flux balance law $`\frac{\mathrm{d} \mathbf{J}}{\mathrm{d}t}  = - \mathcal{F}_{\mathbf{J}}`$. The difference between the two is called a Schott term, $`\mathbf{J}_\mathrm{Schott}`$, which is nonvanishing for elliptic orbits (after orbit averaging) starting at 4PN. The expression of the binding angular momentum in terms of the orbital (radial and azimuthal) frequencies depends on an arbitrary scale $b_0$, which is related to the choice of slicing in relating near-zone and far zone quantities. This arbitrary constant also appears in the flux, and drops out of the balance law. For nonprecessing binaries, only the direction of the angular momentum is constant and orthogonal to the orbital plane: only the norm $`J = |\mathbf{J}|`$ is needed.

For nonspinning compact binaries on elliptic orbits:
* the file ``angular_momentum_conservative.txt`` contains the conservative angular momentum in terms of the Blanchet orbital parameters $(x,\iota)$, which are directly related to the radial and azimuthal frequencies (see below)
* the file  ``lambda0_expansion.txt`` contains the small $e$ expansion  of the function $\lambda_0(e)$ up to $\mathcal{O}(e^8)$

The conservative angular momentum is given at 4PN accuracy. The Schott term, and thus the binding angular momentum, are not known at 4PN. At 3PN, the Schott term is vanishing and the binding angular momentum is thus identical to the conservative angular momentum. 

The small $e$ expansion of $\lambda_0(e)$ is provided up to $\mathcal{O}(e^8)$ in ``lambda0_expansion.txt``.

## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless Blanchet parameter $x = \left(\frac{G m \omega}{c^3}\right)^{\frac{2}{3}}$, where $\omega$ is the dimensionful azimuthal frequency
* ``\[Iota]`` is the dimensionless Blanchet parameter $\iota = \frac{3 x}{\omega/n-1}$, where $n$ is the dimensionful radial frequency
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{(m_1 + m_2)^2}$

We also introduce the following enhancement function defined in (4.15), (4.28), and (4.40) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1):
* ``\[Lambda]0[e]`` corresponds to $\lambda_0(e)$

where $e$ is a dummy (eccentricity-like) variable. The derivative of the enhancement functions with respect to $e$ is denoted with an apostrophe: ``\[Lambda]0'[e]``.


The function $\lambda_0(e)$ is not known in closed form, but replacing it by its small $e$ expansion leads to very accurate results; see Sec. IV.D of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1). 

## Sources

The conservative angular momentum was obtained:
* at 4PN in
    * (5.11)-(5.12) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1)

The small $e$ expansion of $\lambda_0(e)$ was obtained
* up to $\mathcal{O}(e^8)$ in
    * (4.41) and (4.43) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1)

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
