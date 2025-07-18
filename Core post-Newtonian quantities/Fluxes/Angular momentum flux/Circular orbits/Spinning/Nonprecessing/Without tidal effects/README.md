# Flux of angular momentum for spinning, non-precessing, structureless compact binaries on circular orbits

For spinning, non-precessing (i.e. aligned or antialigned spin), structureless (no tidal effects) compact binaries on circular orbits, the flux of angular momentum through future null infinity is given:
* in terms of the orbital frequency $x$ in the file ``flux_infinity_x.txt``
* in terms of the waveform frequency $x_{22}$ in the file ``flux_infinity_x22.txt``

Note that when expressed in terms of the orbital frequency, the flux depends on the arbitrary scale $b_0$, which is related to the choice of slicing in relating near-zone and far zone quantities. This arbitrary constant also appears in the binding angular momentum, and drops out of the balance law.

The result is given at 4PN accuracy, including all non-spinning, spin-orbit, spin-spin, cubic-in-spin and quartic-in-spin terms.

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

## Sources

The full flux of angular momentum for circular orbits is obtained by dividing the flux of energy by the frequency of the (2,2) mode, $\omega_{22} = \frac{c^3 x^{3/2}}{G m}$, see [arXiv:2504.13245v1](https://arxiv.org/abs/2504.13245v1). For references, see the [readme](https://github.com/davidtrestini/PNpedia/tree/main/Core%20post-Newtonian%20quantities/Fluxes/Energy%20flux/Circular%20orbits/Spinning/Nonprecessing/Without%20tidal%20effects) for the flux of energy for nonspinning compact binaires of circular orbits. 

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
