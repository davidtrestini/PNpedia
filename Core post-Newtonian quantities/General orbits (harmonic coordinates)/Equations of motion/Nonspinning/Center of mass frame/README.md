# Relative acceleration between particle 1 and particle 2 ($a_{1}^i-a_{2}^i$) for nonspinning compact binaries on general orbits in the center of mass frame in harmonic coordinates

For the dissipative part of nonspinning compact binaries on general orbits in the center of mass frame in harmonic coordinates
* the file ``2p5PN_CM_Eoms_nonspinning.txt`` contains the 2.5PN contribution to the relative acceleration in the center of mass frame
* the file ``3p5PN_CM_Eoms_nonspinning.txt`` contains the 3.5PN contribution to the relative acceleration in the center of mass frame
* the file ``I_4PN_CM_nonlocal_nonspinning.txt`` contains the function $I_{\text{CM}\,,\text{4PN}}[i,j]$ defining the integrand of the hereditary tail contribution in the relative acceleration at 4PN (See the definition below)
* the file ``4p5PN_CM_mat_Eoms_nonspinning.txt`` contains the matter part of the 4.5PN contribution to the relative acceleration in the center of mass frame
* the file ``4p5PN_CM_rad_Eoms_nonspinning.txt`` contains the radiation part of the 4.5PN contribution to the relative acceleration in the center of mass frame
* the file ``FluxP3p5PN.txt`` contains the impulsion flux $F_{P}$ entering the 4.5PN radiation contribution to the relative acceleration in the center of mass frame


## Notations

We use the following notations:
* ``G`` is Newton's constant of gravitation
* ``c`` is the speed of light
* ``\[Delta][i,j]`` is the flat metric defined on the 3 dimensional euclidian space
* ``m1`` and ``m2`` are the masses of body 1 and body 2
* ``m`` is the total mass of the binary, $m=m_{1}+m_{2}$
* ``\[Nu]`` is the symmetric mass ratio,  $\nu = m_1 m_2 / m^2$
* ``\[CapitalDelta]`` is the relative mass difference, $\Delta= (m_1-m_2)/m$ 
* ``y1[i]`` and ``y2[i]`` are the positions of body 1 and body 2
* ``r`` is the distance between the two bodies
* ``n[i]`` is the unit separation vector  between the two bodies, $n^i = (y_1^i-y_2^i)/r$
* ``v1[i]`` and ``v2[i]`` are the velocities of body 1 and body 2
* ``v[i]`` is the relative Velocity, $v^i = v_1^i - v_2^i$
* ``v`` is the norm of the relative Velocity
* ``b0`` is an arbitrary scale used for the tail contribution
* ``\[CapitalPi][i]`` is the integral of the impulsion flux,  $\Pi^i = \int_{-\infty}^{t}  dt' \,F_{P}^i(t')$ 
* Scalar products are denoted ``U[i]V[-i] = UV`` :  $nv$, $\Pi v$, $\Pi n$


We have also introduced the functions defining the hereditary contributions :
* $I_{\text{CM}\,\text{4PN}}$, such that the 4PN contribution to the relative acceleration in the center of mass frame reads $ a^i_{\text{CM}\,\text{4PN}} = -\frac{4 G m}{5 c^8} y_1^j \int_0^{+\infty} d\tau \, \ln\left(\frac{c \tau}{2 b0}\right) \Big[I_{\text{CM}\,\text{4PN}}[i,j](t-\tau)+I_{\text{CM}\,\text{4PN}}(t+\tau)\Big]$, where $I_{\text{CM}\,\text{4PN}}[i,j]=\frac{d^7}{dt^7}(m \nu r^2 ( x^i x^j -\frac{1}{3} \delta^{ij} x^2))$ is the time derivative of the quadrupole at Newtonian order in the center of mass frame with $\delta^{ij}$ being the Kronecker delta.
* The total 4.5PN contribution to the relative acceleration in the center of mass frame is the sum of the matter and radiative parts.


## Sources

The equations of motion in the center of mass frame were obtained in harmonic coordinates :
* at 3.5PN in 
(5.9) of [https://arxiv.org/abs/gr-qc/0201001](https://arxiv.org/pdf/gr-qc/0412018)
* at 4.5PN in
(6.10) and (6.12) of [arXiv:2601.06743](https://arxiv.org/abs/2601.06743)