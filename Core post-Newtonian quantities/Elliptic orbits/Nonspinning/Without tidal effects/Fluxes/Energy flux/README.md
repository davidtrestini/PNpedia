# Flux of energy for nonspinning compact binaries on elliptic orbits

The file ``flux_infinity_x_et.txt`` contains the *orbit averaged* flux of energy at $\mathcal{I}^+$ for nonspinning compact binaries on elliptic orbits, expressed in terms of frequency parameter $x$ and the harmonic time eccentricity $e_t$.

The file ``flux_infinity_x_iota.txt`` contains the *orbit averaged* flux of energy at $\mathcal{I}^+$ for nonspinning compact binaries on elliptic orbits, expressed in terms of Blanchet frequency parameters $(x,\iota)$.

The file ``flux_infinity_epsilon_j.txt`` contains the *orbit averaged* flux of energy at $\mathcal{I}^+$ for nonspinning compact binaries on elliptic orbits, expressed in terms of dimensionless energy and angular momentum $(\varepsilon, j)$.

The results are given at 3PN accuracy.

## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless Blanchet parameter $x = \left(\frac{G m \omega}{c^3}\right)^{\frac{2}{3}}$, where $\omega$ is the dimensionful azimuthal frequency
* ``\[Iota]`` is the dimensionless Blanchet parameter $\iota = \frac{3 x}{\omega/n-1}$, where $n$ is the dimensionful radial frequency* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{(m_1+m_2)^2}$
* ``\[CurlyEpsilon]`` is the dimensionless energy parameter $\varepsilon = -\frac{2E}{m \nu c^2}$
* ``j`` is the dimensionless angular momentum parameter $j = -\frac{2 J^2 E}{G^2 m^5 \nu^3}$
* ``et`` is the time-eccentricity $e_t$ from the quasi-Keplerian parametrization **in modified harmonic coordinates**

We also introduce the following *enhancement functions* defined in (6.1) of [arXiv:0711.0250v2](https://arxiv.org/abs/0711.0250v2):
* ``\[CurlyPhi][e]`` corresponds to $\varphi(e)$
* ``\[Psi][e]`` is defined as $\psi(e) = \frac{13696}{8191}\alpha(e) - \frac{16403}{24573}\beta(e) - \frac{112}{24573}\gamma(e)$
*  ``\[Zeta][e]`` is defined as $\zeta(e) = -\frac{1424}{4081}\theta(e) + \frac{16403}{12243}\beta(e) + \frac{16}{1749}\gamma(e)$
*  ``\[Kappa][e]`` is defined as $\kappa(e) = \frac{1 + \frac{85}{6}e^2 + \frac{5171}{192} e^4 + \frac{1751}{192}e^6 + \frac{297}{1024} e^8}{(1-e^2)^{13/2}} + \frac{59920}{116761} \chi(e)$
where $e$ is a dummy (eccentricity-like) variable. The intermediate enhancement functions $\varphi(e_t)$, $\alpha(e_t)$, $\beta(e_t)$, $\gamma(e_t)$, $\theta(e_t)$, $\chi(e_t)$ are defined in (92-102) of arXiv:1607.05409v3 in terms of the Fourier decomposition of the multipolar moments, which are given in (65-70) and (A1-A8) of arXiv:1607.05409v3 [note that there is a typo in (A8), which should acquire a global minus sign]. Another equivalent formulation of the Fourier decomposition of the multipolar moments at Newtonian order is given by (A3-A5) [note however the typo in (A5a)].

We introduce yet another enhancement defined in (4.15), (4.28), and (4.40) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1):
* ``\[Lambda]0[e]`` corresponds to $\lambda_0(e)$

The derivatives of the various enhancement functions with respect to $e$ are denoted with an apostrophe: ``\[CurlyPhi]'[e]``, ``\[Psi]'[e]``, ``\[Zeta]'[e]``, ``\[Kappa]'[e]`` and ``\[Lambda]0'[e]``.

The *(modified) harmonic* time-eccentricity is expressed in terms of energy and angular momentum at 3PN in (25d); this translates to $e_t(\varepsilon,j)$ at 3PN as given in ``et.txt``.

It is common to perform the small eccentricity ($e \ll 1$) expansion  of the enhancement functions. These are given, up to neglected $\mathcal{O}(e^8)$ terms, by the following files:
* ``varphi_expanded.txt`` corresponds to $\varphi(e)$, see (7.1a) of [arXiv:0908.3854v2](https://arxiv.org/abs/0908.3854v2) or (B7a) of [arXiv:1906.06263v2](https://arxiv.org/abs/1906.06263v2)
* ``psi_expanded.txt`` corresponds to $\psi(e)$, see (7.1b) of [arXiv:0908.3854v2](https://arxiv.org/abs/0908.3854v2) or (B7c) of [arXiv:1906.06263v2](https://arxiv.org/abs/1906.06263v2)
* ``zeta_expanded.txt`` corresponds to $\zeta(e)$, see (7.1c) of [arXiv:0908.3854v2](https://arxiv.org/abs/0908.3854v2) or (B7g) of [arXiv:1906.06263v2](https://arxiv.org/abs/1906.06263v2)
* ``kappa_expanded.txt`` corresponds to $\kappa(e)$, see (7.1d) of [arXiv:0908.3854v2](https://arxiv.org/abs/0908.3854v2) or (B7e) of [arXiv:1906.06263v2](https://arxiv.org/abs/1906.06263v2)
* ``lambda_0.txt`` corresponds to $\lambda_0(e)$, see (4.41) and (4.43) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1)

## Sources

This result was obtained:
* at 3PN
    * in terms of $(x,e_t)$ in
        * (8.8) and (8.11) of [arXiv:0711.0302v2](https://arxiv.org/abs/0711.0302v2)
    * in terms of $(x,\iota)$ in
        * (8.2a) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1)
    * in terms of $(\varepsilon,j)$ in
        * (8.1a) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1)

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]