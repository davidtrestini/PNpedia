# Energy for compact binaries on circular orbits expressed in terms of the waveform frequency


The energy associated to compact binaries on circular orbits is expressed in terms of the waveform frequency $\omega_{22}$, defined here as the half-frequency of the (2,2) mode, through the dimensionless variable $x_{22}=(G m \omega_{22}/c^3)^{2/3}$.

The binding energy is decomposed as follows:
* ``energy_NS_binding.txt'' is the binding energy for two nonspinning black holes
* ``energy_S1.txt'' is the linear-in-spin contribution to the energy for spinning black holes
* ``energy_S2.txt'' is the quadratic-in-spin contribution to the energy for two spinning black holes
* ``energy_S3.txt'' is the cubic-in-spin contribution  to the energy for two spinning black holes
* ``energy_S4.txt'' is the quartic-in-spin contribution  to the energy for two spinning black holes
* ``energy_tidal.txt'' is the extra contribution due to tidal effects for nonspinning finite-size compact objects
* ``energy_SID.txt'' is the extra contribution due spin-induced deformability effects for spinning finite-size compact objects

The result is given:
* at 4.5PN for the nonspinning sector (in the case of black holes) --- notice that the 4.5PN coefficient is in fact vanishing
* at 4PN the the linear, quadratic, cubic, and quartic in spin sectors (in the case of black holes) --- notice that the 4PN coefficient of the quartic-in-spin contribution is in fact vanishing
* at $\text{N}^{2.5}\text{LO}$ for the nonspinning tidal sector (i.e, 2.5PN orders beyond the leading tidal contributions)
* at 4PN for the spin-induced deformability effects

## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``m`` is the total mass of the binary, $m = m_1+m_2$
* ``x`` is the dimensionless orbital frequency $x = G m \omega /c^3$, where $\omega$ is the dimensionful orbital frequency
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{(m_1 + m_2)^2}$
* ``\[Delta]`` is the relative mass difference, $\delta = \frac{m_1-m_2}{m}$ such that $\delta^2=1-4\nu$
* ``b0`` is the arbitary constant $b_0$ linked to the choice of foliation
* we adopt PN-counting convention which assumes that the black holes are maximally spinning, where the spins $S_1$ and $S_2$ having dimension $[ML^3T^{-2}]$
* we define the dimensionless spin parameters $\chi_1 = \frac{S_1}{G m_1^2}$ and $\chi_2 = \frac{S_2}{G m_2^2}$, such that maximally spinning black holes satisfy $\chi_{1,2} = 1$.
* for the spinning sector:
    * we define $S = S_1 + S_2 = G (m_1^2 \chi_1+m_2^2 \chi_2)$ and $\Sigma = m\left(\frac{S_2}{m_2}-\frac{S_1}{m_1}\right) = G m\left(m_2 \chi_2 - m_1 \chi_1\right)$
    * ``s`` is one of the reduced spin parameters, $s = \frac{S}{G m^2} = \frac{m_1^2 \chi_1 + m_2^2 \chi_2}{m^2}$
    * ``\[Sigma]`` is the other reduced spin parameter, $\sigma = \frac{\Sigma}{G m^2} = \frac{m_2 \chi_2 - m_1 \chi_1}{m}$
* for the spin-induced deformability sector:
    * we introduce the following spin-induced deformability parameters: (i) quadrupolar $\kappa_1$ and $\kappa_2$; (ii) and octupolar $\lambda_1$ and $\lambda_2$. If body $A\in\{1,2\}$ is a black hole, then $\kappa_A=\lambda_A=1$.
    * ``\[CapitalDelta]\[Kappa]p`` is defined as $\Delta\kappa_+ = \frac{(\kappa_1-1) + (\kappa_2-1)}{2}$ (for black holes, $\Delta\kappa_+=0$)
    * ``\[CapitalDelta]\[Kappa]m`` is defined as $\Delta\kappa_- = \frac{(\kappa_1-1) - (\kappa_2-1)}{2}$ (for black holes, $\Delta\kappa_-=0$)
    * ``\[CapitalDelta]\[Lambda]p`` is defined as $\Delta\lambda_+ = \frac{\lambda_1 + \lambda_2}{2}$ (for black holes, $\Delta\lambda_+=0$)
    * ``\[CapitalDelta]\[Lambda]m`` is defined as $\Delta\lambda_- = \frac{\lambda_1 - \lambda_2}{2}$ (for black holes, $\Delta\lambda_-=0$)
* for the tidal sector:
    * we introduce the mass-type quadrupolar tidal deformabilities of the bodies, $\mu^{(2)}\_1$ and $\mu^{(2)}\_2$, where $\mu^{(2)}\_{1,2} = 0$ for stationary black holes.
    * we define $\mu^{(2)}_+ = \frac{1}{2}\left(\frac{m_2}{m_1} \mu^{(2)}_1 + \frac{m_1}{m_2} \mu^{(2)}_2\right)$
    * we define $\mu^{(2)}_- = \frac{1}{2}\left(\frac{m_2}{m_1} \mu^{(2)}_1 - \frac{m_1}{m_2} \mu^{(2)}_2\right)$
    * ``\[Mu]t2p`` is defined as $\tilde{\mu}^{(2)}\_+ = \left(\frac{c^2}{G m}\right)^5 G \mu^{(2)}\_+$
    * ``\[Mu]t2m`` is defined as $\tilde{\mu}^{(2)}\_- = \left(\frac{c^2}{G m}\right)^5 G \mu^{(2)}\_-$
    * we introduce the current-type quadrupolar tidal deformabilities of the bodies, $\sigma^{(2)}\_1$ and $\sigma^{(2)}\_2$, where $\sigma^{(2)}\_{1,2} = 0$ for stationary black holes.
    * we define $\sigma^{(2)}_+ = \frac{1}{2}\left(\frac{m_2}{m_1} \sigma^{(2)}_1 + \frac{m_1}{m_2} \sigma^{(2)}_2\right)$
    * we define $\sigma^{(2)}_- = \frac{1}{2}\left(\frac{m_2}{m_1} \sigma^{(2)}_1 - \frac{m_1}{m_2} \sigma^{(2)}_2\right)$
    * ``\[Sigma]t2p`` is defined as $\tilde{\sigma}^{(2)}\_+ = \left(\frac{c^2}{G m}\right)^5 G \sigma^{(2)}\_+$
    * ``\[Sigma]t2m`` is defined as $\tilde{\sigma}^{(2)}\_- = \left(\frac{c^2}{G m}\right)^5 G \sigma^{(2)}\_-$
    * we introduce the mass-type octupolar tidal deformabilities of the bodies, $\mu^{(3)}\_1$ and $\mu^{(3)}\_2$, where $\mu^{(3)}\_{1,2} = 0$ for stationary black holes.
    * we define $\mu^{(3)}_+= \frac{1}{2}\left(\frac{m_2}{m_1} \mu^{(3)}_1 + \frac{m_1}{m_2} \mu^{(3)}_2\right)$
    * we define $\mu^{(3)}_- = \frac{1}{2}\left(\frac{m_2}{m_1} \mu^{(3)}_1 - \frac{m_1}{m_2} \mu^{(3)}_2\right)$
    * ``\[Mu]t3p`` is defined as $\tilde{\mu}^{(3)}\_+ = \left(\frac{c^2}{G m}\right)^7 G \mu^{(3)}\_+$
    * ``\[Mu]t3m`` is defined as $\tilde{\mu}^{(3)}\_- = \left(\frac{c^2}{G m}\right)^7 G \mu^{(3)}\_-$

## Sources

This result was obtained:
* in the nonspinning sector
    * at 4.5PN 
         * for the conservative piece in
            * (5.5) of [arXiv:1401.4548v2](https://arxiv.org/abs/1401.4548v2)
            * (5.6) of [arXiv:1711.00283v2](https://arxiv.org/abs/1711.00283v2)
        * for the Schott term in
            * (5.7a) of [arXiv:2504.13245v1](https://arxiv.org/abs/2504.13245v1)
    * at 3.5PN in
        * (41) of [arXiv:2410.23950v2](https://arxiv.org/abs/arXiv:2410.23950v2)
* in the spinning sector
   * in the spin-orbit sector in
        * (12) of [arXiv:2201.05138v1](https://arxiv.org/abs/2201.05138v1)
    * in the spin-spin sector in
        * (12) of [arXiv:2201.05138v1](https://arxiv.org/abs/2201.05138v1)
    * in the cubic-in-spin sector in
        * (5.18) of [arXiv:1411.4118v2](https://arxiv.org/abs/1411.4118v2)
    * in the quartic-in-spin sector (for black holes only) in
        * (62) of [arXiv:1712.08603v2](https://arxiv.org/pdf/1712.08603v2) [the 4PN contribution is vanishing]

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
[Loïc Honet](https://github.com/honetloic) [[0009-0007-2863-6085](https://orcid.org/0009-0007-2863-6085)]