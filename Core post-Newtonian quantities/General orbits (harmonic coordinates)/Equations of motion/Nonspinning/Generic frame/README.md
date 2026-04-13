# Equations of motion of the particle 1 for nonspinning compact binaries on general orbits in a general frame in harmonic coordinates

For the dissipative part of nonspinning compact binaries on general orbits in a general frame in harmonic coordinates
* the file ``2p5PN_Eoms_nonspinning.txt`` contains the 2.5PN contribution to the equations of motion
* the file ``3p5PN_Eoms_nonspinning.txt`` contains the 3.5PN contribution to the equations of motion
* the file ``I_4PN_nonlocal_nonspinning.txt`` contains the function $I_{\text{4PN}}[i,j]$ defining the integrand of the hereditary tail contribution in the equations of motion at 4PN (See the definition below)
* the file ``4p5PN_local_Eoms_nonspinning.txt`` contains the local 4.5PN contribution to the equations of motion
* the file ``I_4p5PN_nonlocal_nonspinning.txt`` contains the function $I_{\text{4.5PN}}[i]$ defining the integrand of the hereditary 4.5PN contribution to the equations of motion (See the definition below)


## Notations

We use the following notations:
* ``Pi`` is $\pi \approx 3.14$
* ``\[Gamma]E`` is the Euler's constant $\gamma_\text{E} \approx 0.58$
* ``qbar`` is  $4 \pi e^{\gamma_\text{E}} $
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``m1`` and ``m2`` are the masses of body 1 and body 2
* ``\[Delta][i,j]`` is the flat metric defined on the 3 dimensional euclidian space
* ``y1[i]`` and ``y2[i]`` are the positions of body 1 and body 2
* ``n12[i]`` is the unit separation vector between the two bodies, $n_{12}^i = (y_1^i-y_2^i)/r_{12}$
* ``v1[i]`` and ``v2[i]`` are the velocities of body 1 and body 2
* ``v12[i]`` is the relative Velocity, $v_{12}^i = v_{1}^i - v_{2}^i$
* Scalar products are denoted ``U[i]V[-i] = UV`` :  n12y1, n12y2, n12v1, n12v2, y1y1, y1y2, y2y2, y1v1, y1v2, y2v1, y2v2, v1v1, v1v2, v2v2
* ``b0`` is an arbitrary scale used for the tail contribution
* ``\[Epsilon]`` is d - 3 the dimensional regularization small parameter
* ``l0`` is an arbitrary length associated with dimensional regularization

We have also introduced the functions defining the hereditary contributions :
* $I_{\text{4PN}}$, such that the 4PN contribution to the equations of motion reads $ a^i_{\text{4PN}} = \frac{4 G (m1+m2)}{5 c^8} y_1^j \int_0^{+\infty} d\tau \, \ln\left(\frac{c \tau}{2 b0}\right) \Big[I_{\text{4PN}}[i,j](t-\tau)+I_{\text{4PN}}[i,j](t+\tau)\Big]$, where $I_{\text{4PN}}[i,j]=\frac{d^7}{dt^7}(m_1 y_1^i y_1^j + m_2 y_2^i y_2^j -\frac{1}{3} \delta^{ij}(m_2 y_{2}y_{2} + m_{1} y_{1} y_{1}))$ is the time derivative of the quadrupole at Newtonian order with $\delta^{ij}$ the Kronecker delta.
* $I_{\text{4.5PN}}$, such that the hereditary 4.5PN contribution to the equations of motion reads $ a^i_{\text{4.5PN}}= \frac{2 G }{c^9} \int_0^{+\infty} d\tau \, \ln\left(\frac{c \tau \sqrt{\pi}}{l_0}\right) \Big[I_{\text{4.5PN}}[i](t-\tau)-I_{\text{4.5PN}}[i](t+\tau)\Big]$, where $I_{\text{4.5PN}}[i]=\frac{d^4}{dt^4}(\frac{G^3 n_{12}^{i}}{r_{12}}(m_2^3 m_1 - m_1^3 m_2))$

## Sources

The equations of motion were obtained in harmonic coordinates :
* at 3.5PN in 
(5.4) of [https://arxiv.org/abs/gr-qc/0201001](https://arxiv.org/pdf/gr-qc/0412018)
* at 4.5PN in
(5.5)-(5.6) of [arXiv:2601.06743](https://arxiv.org/abs/2601.06743)