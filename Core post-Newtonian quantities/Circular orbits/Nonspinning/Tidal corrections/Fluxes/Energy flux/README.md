# Flux of energy for non-spinning, non-precessing, compact binaries with adiabatic tidal deformation on circular orbits

We decompose flux of energy at infinity as follows:

$$\mathcal{F} = \mathcal{F}^{\text{non-tidal}} + \mathcal{F}^{\text{tidal}} \.$$

The flux  $\mathcal{F}^{\text{tidal}}$ (given in ``flux_tidal_infinity.txt``) is expressed in terms of the dimensionless waveform frequency $x$.

The result is given at the $\text{N}^{2.5}\text{LO}$ accuracy, that is at relative 2.5PN order beyond the leading tidal term.

For state-of-the-art results for the flux in the non-tidal sector $\mathcal{F}^{\text{non-tidal}}$, see other branches of the repository.

Currently, the presented result do not include effects due to dynamical or dissipative tides. 

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

This result was obtained:
* at the $\text{N}^{2.5}\text{LO}$ in
    * (6.2) of [arXiv:2005.13367v4](https://arxiv.org/abs/2005.13367v4)
    * (4.3) of [arXiv:2412.14249v2](https://arxiv.org/abs/2412.14249v2)
    * the ancillary file ``Amplitude_modes_adiabatic_tides_2.5PN.m`` of [arXiv:2412.14249v2](https://arxiv.org/abs/2412.14249v2)


## Endorsers

[Eve Dones](https://github.com/evedones) [[0009-0003-0239-4584](https://orcid.org/0009-0003-0239-4584)]
