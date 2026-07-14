# Conventions

## Orbital frequency vs waveform frequency

There are also two notions of frequency: (i) the orbital frequency $\omega$, and (ii) the waveform frequency $\omega_{22}$, defined here as the half-frequency of the (2,2) mode. These frequencies typically enter expressions through the dimensionless variables $x=(G m \omega/c^3)^{2/3}$ and $x_{22}=(G m \omega_{22}/c^3)^{2/3}$, which also count PN order. There differ by small 4PN terms : $x_{22} = x (1+\mathcal{O}(x^4))$

## Conservative energy versus binding energy

As pointed out in [arXiv:2504.13245v1](https://arxiv.org/abs/2504.13245v1), there are also two distinct notions of energy. The conservative energy $E_\mathrm{cons}$ is a constant of motion under the conservative equations of motion. The binding energy $E$ enters the flux balance law $d E/d t = - F_E$. The difference between the two is called a Schott term, $E_\mathrm{Schott}$, which is nonvanishing for circular orbits starting at 4PN.

The expression of the binding energy $E$ in terms of the orbital frequency $x$ depends on an arbitrary scale $b_0$, which is related to the choice of slicing in relating near-zone and far zone quantities. Conversely, the expression of the binding energy $E$ in terms of the waveform frequency $x_{22}$ does not depend on the scale $b_0$.

## What is presented

In the section ``In terms of orbital frequency'', we provide:
* the *conservative* energy in terms of orbital frequency, $E_\text{cons}(x)$ [does not depend on $b_0$]
* the Schott term in terms of orbital frequency, $E_\text{Schott}(x)$ [depends on $b_0$]
* the binding energy in terms of orbital frequency is simply obtained by summing these two contributions: $E(x) = E_\text{cons}(x)+E_\text{Schott}(x)$
In the section ``In terms of waveform frequency'', we provide:
* the *binding* energy in terms of waveform frequency, $E(x_{22})$ [does not depend on $b_0$]