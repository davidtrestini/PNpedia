# Waveform for non-spinning, non-precessing, compact binaries with adiabatic tidal deformation on circular orbits

We decompose the waveform into $(\ell,m)$ modes as follows:

$$h = h_+ - i h_- = \frac{G m}{c^2 R} \sum_{\ell=0}^{\infty} \sum_{m=-\ell}^{\ell} h_{\ell m} e^{- i m \psi}\ \_{-2}Y\_{\ell m}(\theta,\phi) \.$$

We decompose the amplitudes as

$$h\_{\ell m} = h^{\text{non-tidal}}\_{\ell m} + h^{\text{tidal}}\_{\ell m} \.$$ 

We also decompose the phase as 

$$\psi = \psi^{\text{non-tidal}} + \psi^{\text{tidal}} = \phi^{\text{non-tidal}} + \phi^{\text{tidal}} - \frac{2 G \mathcal{M} \omega}{c^3}\log{\left(\frac{\omega}{\omega_0}\right)}  $$

where $\mathcal{M} = m + \frac{E}{c^2}$ is the ADM mass, and $E=E^{\text{non-tidal}} + E^{\text{tidal}}$ the binding energy.

The orbital phase $\phi^{\text{tidal}}$ (given in ``orbital_phase_tidal.txt``) and non-vanishing amplitudes $h^{\text{tidal}}\_{\ell m}$ (given in ``h_tidal_2_2.txt``,``h_tidal_2_1.txt``, etc.) are expressed in terms of the dimensionless waveform frequency $x$. The $h^{\text{tidal}}\_{\ell m}$ modes for $m\<0$ are not presented, because they are trivially related to the $m>0$ modes through the relation $h_{\ell (-m)} = (-1)^\ell (h_{\ell m})^*$, where the star denotes a complex conjugate.

Note that we provide the result in the tidal sector for the orbital phase $\phi^{\text{tidal}}$, and not the phase $\psi^{\text{tidal}}$, to match the result published in [arXiv:2412.14249v2](https://arxiv.org/abs/2412.14249v2). 

Results are given at the $\text{N}^{2.5}\text{LO}$ accuracy, that is at relative 2.5PN order beyond the leading tidal term. 

For state-of-the-art results in the non-tidal sector ($h^{\text{non-tidal}}\_{\ell m}$ and $\psi^{non-tidal}$), see other branches of the repository.

Currently, presented results do not include effects due to dynamical or dissipative tides. 

## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``x`` is the orbital frequency $x$
* ``m`` us the total mass of the binary, $m = m_1+m_2$
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{m^2}$
* ``\[Delta]`` is the relative mass difference, $\delta = \frac{m_1-m_2}{m}$ such that $\delta^2=1-4\nu$
* we introduce the mass-type quadrupolar tidal deformabilities of the bodies, $\mu^{(2)}\_1$ and $\mu^{(2)}\_2$, where $\mu^{(2)}\_{1,2} = 0$ for stationary black holes.
* ``\[Mu]2p`` is defined as $\mu^{(2)}_+ = 1/2\left(\frac{m_2}{m_1} \mu^{(2)}_1 + \frac{m_1}{m_2} \mu^{(2)}_2\right)$
* ``\[Mu]2m`` is defined as $\mu^{(2)}_- = 1/2\left(\frac{m_2}{m_1} \mu^{(2)}_1 - \frac{m_1}{m_2} \mu^{(2)}_2\right)$
* ``\[Mu]t2p`` is defined as $\tilde{\mu}^{(2)}\_+ = \left(\frac{c^2}{G m}\right)^5 G \mu^{(2)}\_+$
* ``\[Mu]t2m`` is defined as $\tilde{\mu}^{(2)}\_- = \left(\frac{c^2}{G m}\right)^5 G \mu^{(2)}\_-$
* we introduce the current-type quadrupolar tidal deformabilities of the bodies, $\sigma^{(2)}\_1$ and $\sigma^{(2)}\_2$, where $\sigma^{(2)}\_{1,2} = 0$ for stationary black holes.
* ``\[Sigma]2p`` is defined as $\sigma^{(2)}_+ = 1/2\left(\frac{m_2}{m_1} \sigma^{(2)}_1 + \frac{m_1}{m_2} \sigma^{(2)}_2\right)$
* ``\[Sigma]2m`` is defined as $\sigma^{(2)}_- = 1/2\left(\frac{m_2}{m_1} \sigma^{(2)}_1 - \frac{m_1}{m_2} \sigma^{(2)}_2\right)$
* ``\[Sigma]t2p`` is defined as $\tilde{\sigma}^{(2)}\_+ = \left(\frac{c^2}{G m}\right)^5 G \sigma^{(2)}\_+$
* ``\[Sigma]t2m`` is defined as $\tilde{\sigma}^{(2)}\_- = \left(\frac{c^2}{G m}\right)^5 G \sigma^{(2)}\_-$
* we introduce the mass-type octupolar tidal deformabilities of the bodies, $\mu^{(3)}\_1$ and $\mu^{(3)}\_2$, where $\mu^{(3)}\_{1,2} = 0$ for stationary black holes.
* ``\[Mu]3p`` is defined as $\mu^{(3)}_+ = 1/2\left(\frac{m_2}{m_1} \mu^{(3)}_1 + \frac{m_1}{m_2} \mu^{(3)}_2\right)$
* ``\[Mu]3m`` is defined as $\mu^{(3)}_- = 1/2\left(\frac{m_2}{m_1} \mu^{(3)}_1 - \frac{m_1}{m_2} \mu^{(3)}_2\right)$
* ``\[Mu]t3p`` is defined as $\tilde{\mu}^{(3)}\_+ = \left(\frac{c^2}{G m}\right)^7 G \mu^{(3)}\_+$
* ``\[Mu]t3m`` is defined as $\tilde{\mu}^{(3)}\_- = \left(\frac{c^2}{G m}\right)^7 G \mu^{(3)}\_-$

## Sources

The non vanishing $h^{\text{tidal}}\_{\ell m}$ modes (for $m>0$) were obtained:
* at the $\text{N}^{2.5}\text{LO}$ in:
	* (4.13) of [arXiv:2412.14249v2](https://arxiv.org/abs/2412.14249v2)
  * the ancillary file ``Amplitude_modes_adiabatic_tides_2.5PN.m`` of [arXiv:2412.14249v2](https://arxiv.org/abs/2412.14249v2)  

The orbital phase $\phi^{tidal}$ was obtained:
* at the $\text{N}^{2.5}\text{LO}$ in:
	* (4.8) of [arXiv:2412.14249v2](https://arxiv.org/abs/2412.14249v2) 
  * the ancillary file ``Amplitude_modes_adiabatic_tides_2.5PN.m`` of [arXiv:2412.14249v2](https://arxiv.org/abs/2412.14249v2)  

## Endorsers

[Eve Dones](https://github.com/evedones) [[0009-0003-0239-4584](https://orcid.org/0009-0003-0239-4584)]
