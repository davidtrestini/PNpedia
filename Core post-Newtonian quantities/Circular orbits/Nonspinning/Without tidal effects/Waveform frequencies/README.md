# Waveform frequencies associated to each $(\ell,m)$ mode in terms of the orbital frequency

For nonspinning compact binaries on circular orbits, we report here the relation between the waveform frequencies and orbital frequency.

The orbital frequency $\omega$ is represented by the dimensionless parameter $x= \left(\frac{G m \omega}{c^3}\right)^{\frac{2}{3}}$. 
The metric can be decomposed as

$$ h = h_+ - \mathrm{i} h_\times = \frac{1}{R}\sum_{\ell=2}^\infty \sum_{m=-\ell}^\ell |h_{\ell m}| \mathrm{e}^{\mathrm{i} m \phi_{\ell m}}$$

where $|h_{\ell m}|$ are real amplitudes and $\phi_{\ell m}$ are real phases. We define the waveform frequency $\omega_{\ell m} = \frac{\mathrm{d} \phi_{\ell m}}{ \mathrm{d} t}$ and the associated dimensionless parameter $x_{\ell m} =  \left(\frac{G m \omega_{\ell m}}{c^3}\right)^{\frac{2}{3}}$.

We express $x_{\ell m}$ is terms of $x$ for the following values of $(\ell, m)$
* $x_{22}$ is expressed at 4.5PN accuracy in ``x22.txt``

## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless orbital frequency $x = G m \omega /c^3$, where $\omega$ is the dimensionful orbital frequency
* ``x22`` is the dimensionless waveform frequency $x_{22} = G m \omega_{22} /c^3$, where $\omega_{22}$ is the dimensionful (half-)frequency of the $(2,2)$ mode
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{(m_1 + m_2)^2}$
* ``b0`` is an arbitary constant lengthscale $b_0$, which is linked to the choice of foliation. 

## Sources

This result was obtained:
* at 4.5PN in
    * (6.10) of [arXiv:2304.11186v4](https://arxiv.org/abs/2304.11186v4)
    * (19) of [arXiv:2407.00366v2](https://arxiv.org/abs/2407.00366v2)
    * (6.5) of [arXiv:2504.13245v2](https://arxiv.org/abs/2504.13245v2) has a typo: $-32/15$ should be $576/5$

## Endorsers

[Loïc Honet](https://github.com/honetloic) [[0009-0007-2863-6085](https://orcid.org/0009-0007-2863-6085)]
[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]