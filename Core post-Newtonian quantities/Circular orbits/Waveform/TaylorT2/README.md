# Waveform for circular orbits: TaylorT2 approximant

We decompose the waveform into $(\ell,m)$ modes as follows:

$$h = h_+ - i h_- = \frac{G m}{c^2 R}\sum_{\ell=0}^{\infty} \sum_{m=-\ell}^{\ell} h_{\ell m} e^{- i m \psi}\ \_{-2}Y_{\ell m}(\theta,\phi)$$

The waveform phase $\psi = \mathrm{d} \omega_{22}/\mathrm{d} t$ and the amplitudes $h_{\ell m}$ are expressed in terms of the dimensionless waveform frequency $x_{22} = (G m \omega_{22}/c^3)^{2/3}$. We also introduce the ``chirp'' as the right-hand side of the evolution equation for the waveform frequency: $\dot{x}_{22} = \cdots$. 

The waveform phase $\psi(x_{22})$ is given in the ``Phase'' folder and is  decomposed as follows:
* ``phase_NS.txt'' is the phase of two nonspinning black holes
* ``phase_S1.txt'' is the linear-in-spin contribution to the phase of two spinning black holes
* ``phase_S2.txt'' is the quadratic-in-spin contribution to the phase of two spinning black holes
* ``phase_S3.txt'' is the cubic-in-spin contribution to the phase of two spinning black holes
* ``phase_tidal.txt'' is the nonspinning tidal contribution to the phase of two compact objects

The chirp $\dot{x}\_{22}(x\_{22})$ is given in the ``Chirp'' folder and is  decomposed as follows:
* ``chirp_NS.txt'' is the phase of two nonspinning black holes
* ``chirp_S1.txt'' is the linear-in-spin contribution to the chirp of two spinning black holes
* ``chirp_S2.txt'' is the quadratic-in-spin contribution to the chirp of two spinning black holes
* ``chirp_S3.txt'' is the cubic-in-spin contribution to the chirp of two spinning black holes

The waveform amplitudes $h_{\ell m}$ are given in the ``Amplitudes'' folder. 

Each mode is given in a subfolder ``h_l_m''. Within each of these folders:
* ``h_l_m_NS.txt'' are the amplitudes for two nonspinning black holes
* ``h_l_m_tidal.txt'' are the nonspinning tidal contributions to the amplitude for two nonspinning compact objects
The $h_{\ell m}$ modes for $m<0$ are not presented, because they are trivially related to the $m>0$ modes through the relation $h_{\ell (-m)} = (-1)^\ell (h_{\ell m})^*$, where the star denotes a complex conjugate.

The time evolution of the waveform frequency is itself given ``chirp.txt'', in which $x$ is expressed in terms of the dimensionless time variable $\tau =  \frac{\nu c^3(t-t_0)}{5 G m}$. 

The expression for $h_{22}$ is given with 4PN accuracy. Expressions for the phase $\psi(x_{22})$, chirp $\dot{x}\_{22}(x\_{22})$ and the other $h_{lm}$ modes are presented with 3.5PN accuracy.

## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``m`` is the total mass of the binary, $m = m_1+m_2$
* ``x22`` is the dimensionless waveform frequency $x_{22} = G m \omega_{22} /c^3$, where $\omega_{22}$ is the dimensionful (half-)frequency of the $(2,2)$ mode
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{(m_1 + m_2)^2}$
* ``\[Delta]`` is the relative mass difference, $\delta = \frac{m_1-m_2}{m}$ such that $\delta^2=1-4\nu$
* ``b0`` is the arbitary constant $b_0$ linked to the choice of foliation
* we adopt PN-counting convention which assumes that the black holes are maximally spinning, where the spins $S_1$ and $S_2$ having dimension $[ML^3T^{-2}]$
* we define the dimensionless spin parameters $\chi_1 = \frac{S_1}{G m_1^2}$ and $\chi_2 = \frac{S_2}{G m_2^2}$, such that maximally spinning black holes satisfy $\chi_{1,2} = 1$.
* for the spinning sector:
    * we define $S = S_1 + S_2 = G (m_1^2 \chi_1+m_2^2 \chi_2)$ and $\Sigma = m\left(\frac{S_2}{m_2}-\frac{S_1}{m_1}\right) = G m\left(m_2 \chi_2 - m_1 \chi_1\right)$
    * ``s`` is one of the reduced spin parameters, $s = \frac{S}{G m^2} = \frac{m_1^2 \chi_1 + m_2^2 \chi_2}{m^2}$
    * ``\[Sigma]`` is the other reduced spin parameter, $\sigma = \frac{\Sigma}{G m^2} = \frac{m_2 \chi_2 - m_1 \chi_1}{m}$

## Sources

The phase $\psi(x_{22})$ was obtained:
* for nonspinning black holes:
    * at 4.5PN
        * in (8) of [arXiv:2304.11185v4](https://arxiv.org/abs/2304.11185v4) [ignoring horizon flux contributions]
* in the spinning sector for black holes
    * at 4PN
        * in (14) of [arXiv:2201.05138v1](https://arxiv.org/abs/2201.05138v1) for the linear-in-spin and quartic-in-spin sectors only
        * in (6.20)-(6.21) of [arXiv:1411.4118v2](https://arxiv.org/abs/1411.4118v2) for the cubic-in-spin sector only
    * at 3.5PN
        * in (49) of [arXiv:2410.23950v2](https://arxiv.org/abs/2410.23950v2) 
* in the tidal sector
    * at the $\text{N}^{2.5}\text{LO}$ in:
        * in (4.8) of [arXiv:2412.14249v2](https://arxiv.org/abs/2412.14249v2)

The chirp $\dot{x}\_{22}(x\_{22})$ was obtained:
* for nonspinning black holes:
    * at 4.5PN
        * by deduction from (8) of [arXiv:2304.11185v4](https://arxiv.org/abs/2304.11185v4) [ignoring horizon flux contributions]
    * at 3.5PN
        * in (3.19) of [arXiv:0812.0069v2](https://arxiv.org/abs/0812.0069v2)
* in the spinning sector for black holes
    * at 3.5PN
        * in (45) of [arXiv:2410.23950v2](https://arxiv.org/abs/2410.23950v2) 
        
The $h_{\ell m}(x_{22})$ modes were obtained:
* for nonspinning black holes:
    * for the $(2,2)$ mode
        * at 4PN
            * in (11) of [arXiv:2304.11185v4](https://arxiv.org/abs/2304.11185v4) 
    * for the other oscillatory ($m>0$) modes
        * at 3.5PN
            * in (3.4) and (3.2) of [arXiv:2210.15602v2](https://arxiv.org/abs/2210.15602v2)
    * for the nonoscillatory ($m=0$) modes
        * at 3.5PN
            * in (38) of [arXiv:2410.23950v2](https://arxiv.org/abs/2410.23950v2)
        * at 3PN
            * in (4.3) of [arXiv:0812.0069v2](https://arxiv.org/abs/0812.0069v2)
* in the spinning sector for black holes
    * for the oscillatory ($m>0$) modes      
        * at 3.5PN in:
            * the ancillary file ``modes_PNexp_full_35PN.dat.m`` of [arXiv:2210.15602v2](https://arxiv.org/abs/2210.15602v2)  
    * for the nonoscillatory ($m=0$) modes
        * in (50) and the Supplemental Material of [arXiv:2410.23950v2](https://arxiv.org/abs/2410.23950v2)
* in the tidal sector
    * at the $\text{N}^{2.5}\text{LO}$ in:
	    * (4.13) of [arXiv:2412.14249v2](https://arxiv.org/abs/2412.14249v2) 

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)] --- except tidal contributions
[Eve Dones](https://github.com/evedones) [[0009-0003-0239-4584](https://orcid.org/0009-0003-0239-4584)] --- tidal contributions only
