# Waveform for spinning compact binaries on circular orbits, including tidal effects

We decompose the waveform into $(\ell,m)$ modes as follows:

$$h = h_+ - i h_- = \frac{1}{R}\sum_{\ell=0}^{\infty} \sum_{m=-\ell}^{\ell} h_{\ell m} e^{- i m \psi}\ \_{-2}Y_{\ell m}(\theta,\phi)$$

The phase $\psi$ and the amplitudes $h_{\ell m}$ (given in ``h_2_2.txt``,``h_2_1.txt``, etc.) are expressed in terms of the dimensionless waveform frequency $x$. The time evolution of the waveform frequency is itself given ``chirp.txt``, in which $x$ is expressed in terms of the dimensionless time variable $\tau =  \frac{\nu c^3(t-t_0)}{5 G m}$. The $h_{\ell m}$ modes for $m<0$ are not presented, because they are trivially related to the $m>0$ modes through the relation $h_{\ell (-m)} = (-1)^\ell (h_{\ell m})^*$, where the star denotes a complex conjugate.

Expressions for the phase $\psi(x)$ and chirp $x(\tau)$ are known with 3.5PN when ignoring the contribution of the horizon fluxes or tidal heating, and are given, respectively, by ``phase_ignoring_horizon_flux.txt`` and ``chirp_ignoring_horizon_flux.txt``. Since the horizon flux contribution start at 2.5PN, the full phase and chirp are only known at 2PN, and are given in ``phase.txt`` and ``chirp.txt``. The $h_{lm}$ modes for $m>0$ are given with 3.5PN accuracy, and the $h_{l0}$ modes are unknown.


## Notations

We use the following notations:
* ``I`` is the imaginary unit $i$
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the dimensionless frequency of the (2,2) mode of the GW
* ``m`` us the total mass of the binary, $m = m_1+m_2$
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{m^2}$
* ``\[Delta]`` is the relative mass difference, $\delta = \frac{m_1-m_2}{m}$ such that $\delta^2=1-4\nu$
* ``[Tau]`` is the dimensionless time variable $\tau = \frac{\nu c^3(t-t_0)}{5 G m}$
* we adopt PN-counting convention which assumes that the black holes are maximally spinning, where the spins $S_1$ and $S_2$ having dimension $[ML^3T^{-2}]$
* we define the dimensionless spin parameters $\chi_1 = \frac{S_1}{G m_1^2}$ and $\chi_2 = \frac{S_2}{G m_2^2}$, such that maximally spinning black holes satisfy $\chi_{1,2} = 1$.
* we define $S = S_1 + S_2 = G (m_1^2 \chi_1+m_2^2 \chi_2)$ and $\Sigma = m\left(\frac{S_2}{m_2}-\frac{S_1}{m_1}\right) = G m\left(m_2 \chi_2 - m_1 \chi_1\right)$
* ``s`` is one of the reduced spin parameters, $s = \frac{S}{G m^2} = \frac{m_1^2 \chi_1 + m_2^2 \chi_2}{m^2}$
* ``\[Sigma]`` is the other reduced spin parameter, $\sigma = \frac{\Sigma}{G m^2} = \frac{m_2 \chi_2 - m_1 \chi_1}{m}$

## Sources

The phase (ignoring horizon fluxes) was obtained:
* in the nonspinning sector:
    * at 4.5PN in
        * (8) of [arXiv:2304.11185v4](https://arxiv.org/abs/2304.11185v4)
* in the spin-orbit sector:
    * at 4PN in
        * (14) of [arXiv:2201.05138v1](https://arxiv.org/abs/2201.05138v1)
* in the spin-spin sector:
    * at 4PN in
        * (14) of [arXiv:2201.05138v1](https://arxiv.org/abs/2201.05138v1)
* in the cubic-in-spin sector
    * at 4PN in
        * (6.20)-(6.21) of [arXiv:1411.4118v2](https://arxiv.org/abs/1411.4118v2)
    
The chirp was obtained:
* in the nonspinning sector:
    * at 4.5PN in 
        * (6) of [arXiv:2304.11185v4](https://arxiv.org/abs/2304.11185v4)
    * at 3.5PN in the spinning sector by computing it directly from the fluxes and energy at 3.5PN, see the relevant sections of PNpedia
    
The $h_{\ell m}$ modes (for $m>0$) were obtained:
* in the nonspinning sector:
    * at 3.5PN in:
        * (3.4) and (3.2) of [arXiv:2210.15602v2](https://arxiv.org/abs/2210.15602v2)
* in the spinning sector:
    * at 3.5PN in:
        * the ancillary file ``modes_PNexp_full_35PN.dat.m`` of [arXiv:2210.15602v2](https://arxiv.org/abs/2210.15602v2)  

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
