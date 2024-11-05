# Flux of energy for nonspinning compact binaries on elliptic orbits

The file ``flux.txt`` contains the flux of energy for nonspinning compact binaries on circular orbits.

## Notations

We use the following notations:
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless frequency of the (2,2) mode of the GW
* ``\[Nu]`` is the symmetric mass ratio
* ``et`` is the time-eccentricity $e_t$ from the quasi-Keplerian parametrization

We also introduce the following *enhancement functions* defined in (6.1) of arXiv:0711.0250v2:
* ``\[CurlyPhi][et]`` is defined as $\varphi(e_t)$
* ``\[Psi][et]`` is defined as $\psi(e_t) = \frac{13696}{8191}\alpha(e_t) - \frac{16403}{24573}\beta(e_t) - \frac{112}{24573}\gamma(e_t)$
*  ``\[Zeta][et]`` is defined as $\zeta(e_t) = -\frac{1424}{4081}\theta(e_t) + \frac{16403}{12243}\beta(e_t) + \frac{16}{1749}\gamma(e_t)$
*  ``\[Kappa][et]`` is defined as $\kappa(e_t) = \frac{1 + \frac{85}{6}e_t^2 + \frac{5171}{192} e_t^4 + \frac{1751}{192}e_t^6 + \frac{297}{1024} e_t^8}{(1-e_t^2)^{13/2}} + \frac{59920}{116761} \chi(e_t)$

The intermediate enhancement functions $\varphi(e_t)$, $\alpha(e_t)$, $\beta(e_t)$, $\gamma(e_t)$, $\theta(e_t)$, $\chi(et_)$ are defined in (92-102) of arXiv:1607.05409v3 in terms of the Fourier decomposition of the multipolar moments, which are given in (65-70) and (A1-A8) of arXiv:1607.05409v3 [note that there is a typo in (A8), which should acquire a global minus sign]. Another equivalent formulation of the Fourier decomposition of the multipolar moments at Newtonian order is given by (A3-A5) [note however the typo in (A5a)].

It is common to perform the small eccentricity ($e_t \ll 1$) expansion  of the enhancement functions. These are given, up to neglected $\mathcal(O)(e_t^8)$ terms, by the following files:
* ``varphi_expanded.txt`` corresponds to $\varphi(e_t)$, see (7.1a) of arXiv:0908.3854v2 or (B7a) of arXiv:1906.06263v2
* ``psi_expanded.txt`` corresponds to $\psi(e_t)$, see (7.1 b) of arXiv:0908.3854v2 or (B7c) of arXiv:1906.06263v2
* ``zeta_expanded.txt`` corresponds to $\zeta(e_t)$, see (7.1 c) of arXiv:0908.3854v2 or (B7g) of arXiv:1906.06263v2
* ``kappa_expanded.txt`` corresponds to $\kappa_(e_t)$, see (7.1 d) of arXiv:0908.3854v2 or (B7e) of arXiv:1906.06263v2

The result is given at 4.5PN accuracy.

## Sources

This result was obtained:
* at 4.5PN in
    * (6.11) of arXiv:2304.11186v4
    * (4) of arXiv:2304.11185v4
* at 3.5PN in
    * (12.9) of arXiv:gr-qc/0105098v3
* at 3PN in 
    * (4.11) of arXiv:2406.03457v2

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
