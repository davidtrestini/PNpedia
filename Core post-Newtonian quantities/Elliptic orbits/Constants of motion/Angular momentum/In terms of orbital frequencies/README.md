# Angular momentum for compact binaries on elliptic orbits

As pointed out in [arXiv:2504.13245v2](https://arxiv.org/abs/2504.13245v2), there are two distinct notions of angular momentum. The conservative angular momentum $`\mathbf{J}_\mathrm{cons}`$ is a constant of motion under the conservative equations of motion. The binding angular momentum $`\mathbf{J}`$ enters the flux balance law $`\frac{\mathrm{d} \mathbf{J}}{\mathrm{d}t}  = - \mathcal{F}_{\mathbf{J}}`$. The difference between the two is called a Schott term, $`\mathbf{J}_\mathrm{Schott}`$, which is nonvanishing for elliptic orbits (after orbit averaging) starting at 4PN. The expression of the binding angular momentum in terms of the orbital (radial and azimuthal) frequencies depends on an arbitrary scale $b_0$, which is related to the choice of slicing in relating near-zone and far zone quantities. This arbitrary constant also appears in the flux, and drops out of the balance law. For nonprecessing binaries, the direction of the angular momentum is constant and orthogonal to the orbital plane: only the norm $`J = |\mathbf{J}|`$ is needed.

For nonprecessing compact binaries on elliptic orbits:
* the file ``angular_momentum_NS_conservative.txt`` contains the conservative, nonspinning contributions to the 4PN angular momentum in terms of the Blanchet frequency parameters $(x,\iota)$, which are directly related to the radial and azimuthal frequencies (see below)
* the file ``angular_momentum_S1.txt`` contains the linear-in-spin contribution to the 4PN angular momentum in terms of the Blanchet frequency parameters $(x,\iota)$ 
* the file ``angular_momentum_S2.txt`` contains the quadratic-in-spin contribution to the 4PN angular momentum in terms of the Blanchet frequency parameters $(x,\iota)$ 
* the file ``angular_momentum_S3.txt`` contains the cubic-in-spin contribution to the 4PN angular momentum in terms of the Blanchet frequency parameters $(x,\iota)$ 
* the file ``angular_momentum_S4.txt`` contains the quartic-in-spin contribution to the 4PN angular momentum in terms of the Blanchet frequency parameters $(x,\iota)$ 
* the file ``angular_momentum_S5.txt`` contains the $S^5$ contribution to the 4PN angular momentum in terms of the Blanchet frequency parameters $(x,\iota)$ 
* the file ``angular_momentum_S6.txt`` contains the $S^6$ contribution to the 4PN angular momentum in terms of the Blanchet frequency parameters $(x,\iota)$ 
* the file ``angular_momentum_SID.txt`` contains the spin-induced deformability contributions to the 4PN angular momentum in terms of the Blanchet frequency parameters $(x,\iota)$ 
* the file  ``lambda0_expanded.txt`` contains the small $e$ expansion  of the function $\lambda_0(e)$ up to $\mathcal{O}(e^8)$

The Schott term, and thus the binding angular momentum, are not known at 4PN. At 3PN, the Schott term is vanishing and the binding angular momentum is thus identical to the conservative angular momentum.


## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``EulerGamma`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``m`` is the total mass of the binary, $m = m_1+m_2$
* ``x`` is the dimensionless Blanchet parameter $x = \left(\frac{G m \omega}{c^3}\right)^{\frac{2}{3}}$, where $\omega$ is the dimensionful azimuthal frequency
* ``\[Iota]`` is the dimensionless Blanchet parameter $\iota = \frac{3 x}{\omega/n-1}$, where $n$ is the dimensionful radial frequency
* ``\[Nu]`` is the symmetric mass ratio, $\nu = \frac{m_1 m_2}{(m_1 + m_2)^2}$
* ``\[Delta]`` is the relative mass difference, $\delta = \frac{m_1-m_2}{m}$ such that $\delta^2=1-4\nu$
* we define the dimensionless spin parameters $\chi_1 = \frac{S_1}{G m_1^2}$ and $\chi_2 = \frac{S_2}{G m_2^2}$, such that maximally spinning black holes satisfy $\chi_{1,2} = 1$.
* for the spinning sector:
    * we define $S = S_1 + S_2 = G (m_1^2 \chi_1+m_2^2 \chi_2)$ and $\Sigma = m\left(\frac{S_2}{m_2}-\frac{S_1}{m_1}\right) = G m\left(m_2 \chi_2 - m_1 \chi_1\right)$
    * ``s`` is one of the reduced spin parameters, $s = \frac{S}{G m^2} = \frac{m_1^2 \chi_1 + m_2^2 \chi_2}{m^2}$
    * ``\[Sigma]`` is the other reduced spin parameter, $\sigma = \frac{\Sigma}{G m^2} = \frac{m_2 \chi_2 - m_1 \chi_1}{m}$
* for the spin-induced deformability sector:
    * we introduce the following spin-induced deformability parameters: (i) quadrupolar $\kappa_1$ and $\kappa_2$; (ii) octupolar $\lambda_1$ and $\lambda_2$; and (iii) hexadecapolar $\mu_1$ and $\mu_2$. If body $A\in\{1,2\}$ is a black hole, then $\kappa_A=\lambda_A=\mu_A=1$.
    * ``\[CapitalDelta]\[Kappa]p`` is defined as $\Delta\kappa_+ = \frac{(\kappa_1-1) + (\kappa_2-1)}{2}$ (for black holes, $\Delta\kappa_+=0$)
    * ``\[CapitalDelta]\[Kappa]m`` is defined as $\Delta\kappa_- = \frac{(\kappa_1-1) - (\kappa_2-1)}{2}$ (for black holes, $\Delta\kappa_-=0$)
    * ``\[CapitalDelta]\[Lambda]p`` is defined as $\Delta\lambda_+ = \frac{(\lambda_1-1) + (\lambda_2-1)}{2}$ (for black holes, $\Delta\lambda_+=0$)
    * ``\[CapitalDelta]\[Lambda]m`` is defined as $\Delta\lambda_- = \frac{(\lambda_1-1) - (\lambda_2-1)}{2}$ (for black holes, $\Delta\lambda_-=0$)
    * ``\[CapitalDelta]\[Mu]p`` is defined as $\Delta\mu_+ = \frac{(\mu_1-1) + (\mu_2-1)}{2}$ (for black holes, $\Delta\mu_+=0$)
    * ``\[CapitalDelta]\[Mu]m`` is defined as $\Delta\mu_- = \frac{(\mu_1-1) - (\mu_2-1)}{2}$ (for black holes, $\Delta\mu_-=0$)
    
We also introduce the following enhancement function defined in (4.15), (4.28), and (4.40) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1):
* ``\[Lambda]0[e]`` corresponds to $\lambda_0(e)$

where $e$ is a dummy (eccentricity-like) variable. The derivative of the enhancement functions with respect to $e$ is denoted with an apostrophe: ``\[Lambda]0'[e]``.


The function $\lambda_0(e)$ is not known in closed form, but replacing it by its small $e$ expansion leads to very accurate results; see Sec. IV.D of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1).

## Sources

The conservative angular momentum was obtained:
* in the nonspinning sector
    * at 4PN in
        * (5.11)-(5.12) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1)
    * at 3PN in
        * (4.40b) of [arXiv:1503.01374v2](https://arxiv.org/abs/1503.01374v2)
* in the spinning sector
    * in the spin-orbit sector
        * at 4PN in
            * (4.7a) of [arXiv:2608.31072v2](https://arxiv.org/abs/2608.31072v2)
    * in the spin-spin sector in
        * at 4PN in
            * (4.7b) of [arXiv:2608.31072v2](https://arxiv.org/abs/2608.31072v2)
    * in the cubic-in-spin sector in
        * at 4PN in
            * (4.7c) of [arXiv:2608.31072v2](https://arxiv.org/abs/2608.31072v2)
    * in the quartic-in-spin sector (for black holes only) in
         * at 4PN in
            * (4.7d) of [arXiv:2608.31072v2](https://arxiv.org/abs/2608.31072v2)
    * in the $S^5$ sector (for black holes only) in
         * at 4PN in
            * (4.7e) of [arXiv:2608.31072v2](https://arxiv.org/abs/2608.31072v2)
    * in the $S^6$ sector (for black holes only) in
         * at 4PN in
            * (4.7f) of [arXiv:2608.31072v2](https://arxiv.org/abs/2608.31072v2)
    * in the spin induced deformability sector in
         * at 4PN in
            * the Supplemental Material of [arXiv:2608.31072v2](https://arxiv.org/abs/2608.31072v2)

 

The small $e$ expansion of $\lambda_0(e)$ was obtained
* up to $\mathcal{O}(e^8)$ in
    * (4.41) and (4.43) of [arXiv:2511.10735v1](https://arxiv.org/abs/2511.10735v1)

The Schott contributions (which enter at 4PN and need to be added to obtain the binding angular momentum) are unknown.

## Endorsers

[David Trestini](https://github.com/davidtrestini) [[0000-0002-4140-0591](https://orcid.org/0000-0002-4140-0591)]
