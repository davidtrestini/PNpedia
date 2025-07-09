# Energy for spinning, non-precessing, compact binaries without tidal deformation on circular orbits

As pointed out in [arXiv:2504.13245v1](https://arxiv.org/abs/2504.13245v1), there are two distinct notions of energy. The conservative energy $E_\mathrm{cons}$ is a constant of motion under the conservative equations of motion. The binding energy $E$ enters the flux balance law $\mathrm{d}E/\mathrm{d}t = - F_E$. The difference between the two is called a Schott term, $E_\mathrm{Schott}$, which is nonvanishing for circular orbits starting at 4PN. The expression of the binding energy in terms of the orbital frequency depends on an arbitrary scale $b_0$, which is related to the choice of slicing in relating near-zone and far zone quantities. This arbitrary constant also appears in the flux, and drops out of the balance law.

For spinning, non-precessing (i.e. aligned or antialigned spin) compact binaries without tidal deformation on circular orbits:
* the file ``energy_conservative.txt`` contains the conservative energy in terms of the orbital frequency $x$
* the file ``energy_binding.txt`` contains the binding energy in terms of the orbital frequency $x$
* the file ``energy_Schott.txt`` contains the Schott term in terms of the orbital frequency $x$

The result is given at 4PN accuracy.

## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless orbital frequency $x = G m \omega /c^3$, where $\omega$ is the dimensionful orbital frequency
* ``x22`` is the dimensionless waveform frequency $x_{22} = G m \omega_{22} /c^3$, where $\omega_{22}$ is the dimensionful (half-)frequency of the $(2,2)$ mode
* ``m`` us the total mass of the binary, $m = m_1+m_2$
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{m^2}$
* ``\[Delta]`` is the relative mass difference, $\delta = \frac{m_1-m_2}{m}$ such that $\delta^2=1-4\nu$
* ``b_0`` is an arbitary constant linked to the choice of foliation
* we adopt PN-counting convention which assumes that the black holes are maximally spinning, where thesspins $S_1$ and $S_2$ having dimension $[ML^3T^{-2}]$
* we define the dimensionless spin parameters $\chi_1 = \frac{S_1}{G m_1^2}$ and $\chi_2 = \frac{S_2}{G m_2^2}$, such that maximally spinning black holes satisfy $\chi_{1,2} = 1$.
* we define $S = S_1 + S_2 = G (m_1^2 \chi_1+m_2^2 \chi_2)$ and $\Sigma = m\left(\frac{S_2}{m_2}-\frac{S_1}{m_1}\right) = G m\left(m_2 \chi_2 - m_1 \chi_1\right)$
* ``s`` is one of the reduced spin parameters, $s = \frac{S}{G m^2} = \frac{m_1^2 \chi_1 + m_2^2 \chi_2}{m^2}$
* ``\[Sigma]`` is the other reduced spin parameter, $\sigma = \frac{\Sigma}{G m^2} = \frac{m_2 \chi_2 - m_1 \chi_1}{m}$
* we introduce the spin-induced quadrupolar deformability of the bodies, $\kappa_1$ and $\kappa_2$, where $\kappa_{1,2}=1$ for black holes

## Sources

This result was obtained:
* at 4PN
    * in the nonspinning sector (conservative piece) in
        * (5.5) of [arXiv:1401.4548v2](https://arxiv.org/abs/1401.4548v2)
        * (5.6) of [arXiv:1711.00283v2](https://arxiv.org/abs/1711.00283v2)
    * in the nonspinning sector (Schott term) in
        * (5.7a), (6.1a) and (6.3a) of [arXiv:2504.13245v1](https://arxiv.org/abs/2504.13245v1)
    * in the spin-orbit sector in
        * (12) of [arXiv:2201.05138v1](https://arxiv.org/abs/2201.05138v1)
    * in the spin-spin sector in
        * (12) of [arXiv:2201.05138v1](https://arxiv.org/abs/2201.05138v1)
    * in the cubic-in-spin sector in
        * (5.18) of [arXiv:1411.4118v2](https://arxiv.org/abs/1411.4118v2)
    * in the quartic-in-spin sector in
        * (62) of [arXiv:1712.08603v2](https://arxiv.org/pdf/1712.08603v2)

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
