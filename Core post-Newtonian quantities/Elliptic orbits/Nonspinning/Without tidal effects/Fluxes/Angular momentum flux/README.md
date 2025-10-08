# Flux of angular momentum for nonspinning compact binaries on elliptic orbits

The file ``flux.txt`` contains the *orbit averaged* flux of angular momentum for nonspinning compact binaries on elliptic orbits.

## Notations

We use the following notations:
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``m`` is the total mass of the binary, $m = m_1 + m_2$
* ``x`` is the dimensionless frequency of the (2,2) mode of the GW
* ``\[Nu]`` is the symmetric mass ratio, $\nu = m_1 m_2 / m^2$
* ``et`` is the time-eccentricity $e_t$ from the quasi-Keplerian parametrization **in modified harmonic coordinates**

We also introduce the following *enhancement functions* defined in (5.28) of [arXiv:0908.3854v2](https://arxiv.org/abs/0908.3854v2):
* ``\[CurlyPhi]$tilde[et]`` corresponds to $\tilde\varphi(e_t)$
* ``\[Psi]$tilde[et]`` is defined as $\tilde\psi(e_t) = \frac{13696}{8191}\tilde\alpha(e_t) - \frac{16403}{24573}\tilde\beta(e_t) - \frac{112}{24573}\tilde\gamma(e_t)$
*  ``\[Zeta]$tilde[et]`` is defined as $\tilde\zeta(e_t) = -\frac{1424}{4081}\tilde\theta(e_t) + \frac{16403}{12243}\tilde\beta(e_t) + \frac{16}{1749}\tilde\gamma(e_t)$
*  ``\[Kappa]$tilde[et]`` is defined as $\tilde\kappa(e_t) = \frac{1 + \frac{229}{32}e_t^2 + \frac{327}{64} e_t^4 + \frac{69}{256}e_t^6}{(1-e_t^2)^{5}} + \frac{59920}{116761} \tilde\chi(e_t)$

The intermediate enhancement functions $\tilde\varphi(e_t)$, $\tilde\alpha(e_t)$, $\tilde\beta(e_t)$, $\tilde\gamma(e_t)$, $\tilde\theta(e_t)$, $\tilde\chi(e_t)$ are defined in (92-102) of arXiv:1607.05409v3 in terms of the Fourier decomposition of the multipolar moments, which are given in (65-70) and (A1-A8) of arXiv:1607.05409v3 [note that there is a typo in (A8), which should acquire a global minus sign]. Another equivalent formulation of the Fourier decomposition of the multipolar moments at Newtonian order is given by (A3-A5) [note however the typo in (A5a)].

It is common to perform the small eccentricity ($e_t \ll 1$) expansion  of the enhancement functions. These are given, up to neglected $\mathcal{O}(e_t^8)$ terms, by the following files:
* ``varphi_tilde_expanded.txt`` corresponds to $\tilde\varphi(e_t)$, see (7.2a) of [arXiv:0908.3854v2](https://arxiv.org/abs/0908.3854v2) or (B7b) of [arXiv:1906.06263v2](https://arxiv.org/abs/1906.06263v2)
* ``psi_tilde_expanded.txt`` corresponds to $\tilde\psi(e_t)$, see (7.2b) of [arXiv:0908.3854v2](https://arxiv.org/abs/0908.3854v2) or (B7d) of [arXiv:1906.06263v2](https://arxiv.org/abs/1906.06263v2)
* ``zeta_tilde_expanded.txt`` corresponds to $\tilde\zeta(e_t)$, see (7.2c) of [arXiv:0908.3854v2](https://arxiv.org/abs/0908.3854v2) or (B7h) of [arXiv:1906.06263v2](https://arxiv.org/abs/1906.06263v2)
* ``kappa_tilde_expanded.txt`` corresponds to $\tilde\kappa(e_t)$, see (7.1d) of [arXiv:0908.3854v2](https://arxiv.org/abs/0908.3854v2) or (B7f) of [arXiv:1906.06263v2](https://arxiv.org/abs/1906.06263v2)

The result is given at 3PN accuracy.

## Sources

This result was obtained:
* at 3PN in
    * (4.10) and (5.29) of [arXiv:0711.0302v2](https://arxiv.org/abs/0711.0302v2)
    
Note that the treatment of the memory contribution of [arXiv:0711.0302v2](https://arxiv.org/abs/0711.0302v2) is incomplete, and is completed in Appendix A of [arXiv:2410.12898v2](https://arxiv.org/abs/2410.12898v2)

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]