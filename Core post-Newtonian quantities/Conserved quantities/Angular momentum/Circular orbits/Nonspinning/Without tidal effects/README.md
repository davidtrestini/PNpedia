# Angular momentum for nonspinning compact binaries on circular orbits

The file ``angular_momentum_conservative.txt`` contains the conservative angular momentum for nonspinning compact binaries on circular orbits in terms of the orbital frequency $x$. It is a constant of motion when using the conservative equations of motion.

The file ``angular_momentum.txt`` contains the binding angular momentum for nonspinning compact binaries on circular orbits in terms of the orbital frequency $x$. It differs from the conservative angular momentum by a Schott term due to dissipative effects. It depends on the arbitrary scale $b_0$, which is related to the choice of slicing in relating near-zone and far zone quantities. The waveform frequency of the $(2,2)$ mode can be related to the orbital frequency, and this relation also depends on $b_0$. Thus, the binding energy in terms of the waveform frequency does not feature the $b_0$ constant, and has the same functional form as the conservative energy in terms of the orbital frequency. See [arXiv:2504.13245v1](https://arxiv.org/abs/2504.13245v1) for details.

## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless orbital frequency $x = G m \omega /c^3$, where $\omega$ is the dimensionful orbital frequency
* ``x22`` is the dimensionless orbital frequency $x_{22} = G m \omega_{22} /c^3$, where $\omega_{22}$ is the dimensionful (half-)frequency of the $(2,2)$ mode
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{(m_1 + m_2)^2}$
* ``b_0`` is an arbitary constant linked to the choice of foliation

The result is given at 4.5PN accuracy.

## Sources

The conservative energy was obtained:
* at 4PN in
    * (5.8) of [arXiv:1711.00283v2](https://arxiv.org/abs/1711.00283v2)
The binding angular momentum was obtained:
* at 4PN in
    * (5.7b), (6.1b) and (6.3b) of [arXiv:2504.13245v1](https://arxiv.org/abs/2504.13245v1)

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]