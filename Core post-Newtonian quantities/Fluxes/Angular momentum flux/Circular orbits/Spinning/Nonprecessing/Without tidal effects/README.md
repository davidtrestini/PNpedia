# Flux of energy for spinning, non-precessing, structureless compact binaries on elliptic orbits

The file ``flux.txt`` contains the flux of energy for spinning, non-precessing (i.e. aligned or antialigned spin), structureless (no tidal effects) compact binaries on elliptic orbits

## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless frequency of the (2,2) mode of the GW
* ``m`` is the total mass of the binary, $m = m_1+m_2$
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{m^2}$
* ``\[Delta]`` is the relative mass difference, $\delta = \frac{m_1-m_2}{m}$ such that $\delta^2=1-4\nu$
* we adopt PN-counting convention which assumes that the black holes are maximally spinning, where the spins $S_1$ and $S_2$ having dimension $[ML^3T^{-2}]$
* we define the dimensionless spin parameters $\chi_1 = \frac{S_1}{G m_1^2}$ and $\chi_2 = \frac{S_2}{G m_2^2}$, such that maximally spinning black holes satisfy $\chi_{1,2} = 1$.
* we define $S = S_1 + S_2 = G (m_1^2 \chi_1+m_2^2 \chi_2)$ and $\Sigma = m\left(\frac{S_2}{m_2}-\frac{S_1}{m_1}\right) = G m\left(m_2 \chi_2 - m_1 \chi_1\right)$
* ``s`` is one of the reduced spin parameters, $s = \frac{S}{G m^2} = \frac{m_1^2 \chi_1 + m_2^2 \chi_2}{m^2}$
* ``\[Sigma]`` is the other reduced spin parameter, $\sigma = \frac{\Sigma}{G m^2} = \frac{m_2 \chi_2 - m_1 \chi_1}{m}$

The result is given at 3.5PN accuracy, including all non-spinning, spin-orbit, spin-spin, and cubic-in-spin terms.

## Sources

The angular momentum flux for circular orbits is obtained by dividing the energy flux by the frequency of the (2,2) mode, $\omega_{22} = \frac{c^3 x^{3/2}}{G m}$. The energy flux was obtained:
* at 4PN
    * in the nonspinning sector in
        * (6.11) of [arXiv:2304.11186v4](https://arxiv.org/abs/2304.11186v4)
        * (4) of [arXiv:2304.11185v4](https://arxiv.org/abs/2304.11185v4)
    * in the spin-orbit sector in
        * (13) of [arXiv:2201.05138v1](https://arxiv.org/abs/2201.05138v1)
    * in the spin-spin sector in
        * (13) of [arXiv:2201.05138v1](https://arxiv.org/abs/2201.05138v1)
    * in the cubic-in-spin sector in
        * (6.18)-(6.19) of [arXiv:1411.4118v2](https://arxiv.org/abs/1411.4118v2)

However, at 4PN, the quartic-in-spin term is unknown, so only the 3.5PN result is presented.
## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
