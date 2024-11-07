# Waveform for nonspinning compact binaries on circular orbits

We decompose the waveform into $(\ell,m)$ modes as follows:

$$h = h_+ - i h_- = \sum_{\ell=0}^{\infty} \sum_{m=-\ell}^{\ell} h_{\ell m} e^{i m \psi}$$

The phase $\psi$ (given in ``phase.txt``) and the amplitudes $h_{\ell m}$ (given in ``h22.txt``,``h21.txt``, etc.) are expressed in terms of the dimensionless waveform frequency $x$. The time evolution of the waveform frequency is itself given ``chirp.txt``, in which $x$ is expressed in terms of the dimensionless time variable $\tau$.

The phase $\phi(x)$ and the chirp $\x(\tau)$ are given at 4.5PN accuracy. The $h_{22}$ mode is given with 4PN accuracy, and all the others are given at 3.5PN accuracy.

## Notations

We use the following notations:
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless frequency of the (2,2) mode of the GW
* ``m`` us the total mass of the binary, $m = m_1+m_2$
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{m^2}$
* ``\[Delta]`` is the relative mass difference, $\delta = \frac{m_1-m_2}{m}$ such that $\delta^2=1-4\nu$
* ``[Tau]`` is the dimensionless time variable $\tau = \frac{\nu c^3(t-t_0)}{5 G m}$

## Sources

The phase was obtained:
* at 4.5PN in
    * (8) of [arXiv:2304.11185v4](https://arxiv.org/abs/2304.11185v4)
    
The chirp was obtained:
* at 4.5PN in 
    * (6) of [arXiv:2304.11185v4](https://arxiv.org/abs/2304.11185v4)
    
The $h_22$ mode was obtained:
* at 4PN in:
    * (6.17) of [arXiv:2304.11186v4](https://arxiv.org/abs/2304.11186v4)
    * (11) of [arXiv:2304.11185v4](https://arxiv.org/abs/2304.11185v4)

The $h_{l0}$ modes were obtained:
* at 3.5PN in:
    * (39) and (26) of [arXiv:2410.23950v2](https://arxiv.org/abs/2410.23950v2)
* at 3PN in:
    * (4.3) and (4.1) of [arXiv:0812.0069v2](https://arxiv.org/abs/0812.0069v2)

The other $h_{lm}$ modes were obtained
* at 3.5PN in:
    * (3.4) and (3.2) of [arXiv:2210.15602v2](https://arxiv.org/abs/2210.15602v2)

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
