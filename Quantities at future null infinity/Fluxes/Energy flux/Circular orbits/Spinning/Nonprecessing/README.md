# Flux of energy for nonspinning compact binaries on elliptic orbits

The file ``flux.txt`` contains the *orbit averaged* flux of energy for nonspinning compact binaries on elliptic orbits.

## Notations

We use the following notations:
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless frequency of the (2,2) mode of the GW
* ``m`` us the total mass of the binary, $m = m_1+m_2$
* ``\[Nu]`` is the symmetric mass ratio, $\nu = m_1 m_2 / m^2$
* we adopt PN-counting convention which assumes that the black holes are maximally spinning, with spins $S_1$ and $S_2$ having dimension $[ML^3T^{-2}]$
* we define the dimensionless spin parameters $s_1 = S_1/(G m_1^2)$ and $s_2 = S_2/(G m_2^2)$, such that maximally spinning black holes satisfy $s_{1,2} = 1$.
* we define $S = S_1 + S_2 = G (m_1^2 s_1+m_2^2 s_2)$ and $\Sigma = m\left(\frac{S_2}{m_2}-\frac{S_1}{m_1}) = G m\left(m_2 s_2 - m_1 s_1\right)$
* ``s`` is one of the reduced spin parameters, $s = \frac{S}{G m^2} = \frac{m_1^2 s_1 + m_2^2 s_2}{m^2}$
* ``\[sigma]`` is the other reduced spin parameter, $\sigma = \frac{\Sigma}{G m^2} = \frac{m_2 s_2 - m_1 s_1}{m}$

The result is given at 4PN accuracy, including all non-spining, spin-orbit, spin-spin, and cubic-in-spin terms.

## Sources

This result was obtained:
* at 3PN in
    * (8.8) and (8.11) of [arXiv:0711.0302v2](https://arxiv.org/abs/0711.0302v2)
## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
