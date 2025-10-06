# Energy for nonspinning compact binaries on circular orbits

As pointed out in [arXiv:2504.13245v1](https://arxiv.org/abs/2504.13245v1), there are two distinct notions of energy. The conservative energy $E_\mathrm{cons}$ is a constant of motion under the conservative equations of motion. The binding energy $E$ enters the flux balance law $d E/d t = - F_E$. The difference between the two is called a Schott term, $E_\mathrm{Schott}$, which is nonvanishing for circular orbits starting at 4PN. The expression of the binding energy in terms of the orbital frequency depends on an arbitrary scale $b_0$, which is related to the choice of slicing in relating near-zone and far zone quantities. This arbitrary constant also appears in the flux, and drops out of the balance law.

For nonspinning compact binaries on circular orbits:
* the file ``energy_conservative.txt`` contains the conservative energy in terms of the orbital frequency $x$
* the file ``energy_binding.txt`` contains the binding energy in terms of the orbital frequency $x$
* the file ``energy_Schott.txt`` contains the Schott term in terms of the orbital frequency $x$

The result is given at 4.5PN accuracy.

## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless orbital frequency $x = G m \omega /c^3$, where $\omega$ is the dimensionful orbital frequency
* ``x22`` is the dimensionless waveform frequency $x_{22} = G m \omega_{22} /c^3$, where $\omega_{22}$ is the dimensionful (half-)frequency of the $(2,2)$ mode
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{(m_1 + m_2)^2}$
* ``b_0`` is an arbitary constant linked to the choice of foliation

## Sources

The conservative energy was obtained:
* at 4PN in
    * (5.5) of [arXiv:1401.4548v2](https://arxiv.org/abs/1401.4548v2)
    * (5.6) of [arXiv:1711.00283v2](https://arxiv.org/abs/1711.00283v2)
The Schott term was obtained:
* at 4PN in
    * (5.7a) of [arXiv:2504.13245v1](https://arxiv.org/abs/2504.13245v1)

Note that all 4.5PN terms vanish for circular orbits, both in the conservative and binding energies.

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
