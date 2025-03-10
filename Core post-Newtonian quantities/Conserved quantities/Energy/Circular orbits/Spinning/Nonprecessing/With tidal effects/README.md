# Energy for spinning, non-precessing, compact binaries with tidal deformation on circular orbits

The file ``energy.txt`` contains the conserved energy for spinning, non-precessing (i.e. aligned or antialigned spin) compact binaries with tidal deformation on circular orbits.

## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless frequency of the (2,2) mode of the GW
* ``m`` us the total mass of the binary, $m = m_1+m_2$
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{m^2}$
* ``\[Delta]`` is the relative mass difference, $\delta = \frac{m_1-m_2}{m}$ such that $\delta^2=1-4\nu$
* we adopt PN-counting convention which assumes that the black holes are maximally spinning, where the spins $S_1$ and $S_2$ having dimension $[ML^3T^{-2}]$
* we define the dimensionless spin parameters $\chi_1 = \frac{S_1}{G m_1^2}$ and $\chi_2 = \frac{S_2}{G m_2^2}$, such that maximally spinning black holes satisfy $\chi_{1,2} = 1$.
* we define $S = S_1 + S_2 = G (m_1^2 \chi_1+m_2^2 \chi_2)$ and $\Sigma = m\left(\frac{S_2}{m_2}-\frac{S_1}{m_1}\right) = G m\left(m_2 \chi_2 - m_1 \chi_1\right)$
* ``s`` is one of the reduced spin parameters, $s = \frac{S}{G m^2} = \frac{m_1^2 \chi_1 + m_2^2 \chi_2}{m^2}$
* ``\[Sigma]`` is the other reduced spin parameter, $\sigma = \frac{\Sigma}{G m^2} = \frac{m_2 \chi_2 - m_1 \chi_1}{m}$
* we introduce the spin-induced quadrupolar deformability of the bodies, $\kappa_1$ and $\kappa_2$, where $\kappa_{1,2}=1$ for black holes
* ``\[Kappa]p`` is defined as $\kappa_+ = \kappa_1 + \kappa_2$
* ``\[Kappa]m`` is defined as $\kappa_- = \kappa_1 - \kappa_2$
* we introduce the spin-induced octupolar deformability of the bodies, $\lambda_1$ and $\lambda_2$, where $\kappa_{1,2}=1$ for black holes
* ``\[Lambda]p`` is defined as $\lambda_+ = \lambda_1 + \lambda_2$
* ``\[Lambda]m`` is defined as $\lambda_- = \lambda_1 - \lambda_2$

The result is given at 3.5PN accuracy.

## Sources

This result was obtained:
* at 4PN
    * in the nonspinning sector in
        * (5.5) of [arXiv:1401.4548v2](https://arxiv.org/abs/1401.4548v2)
        * (5.6) of [arXiv:1711.00283v2](https://arxiv.org/abs/1711.00283v2)
    * in the spin-orbit sector in
        * (12) of [arXiv:2201.05138v1](https://arxiv.org/abs/2201.05138v1)
    * in the spin-spin sector in
        * (12) of [arXiv:2201.05138v1](https://arxiv.org/abs/2201.05138v1)
    * in the cubic-in-spin sector in
        * (5.18) of [arXiv:1411.4118v2](https://arxiv.org/abs/1411.4118v2)
However, at 4PN, the quartic-in-spin term is unknown, so only the 3.5PN result is presented. 

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
