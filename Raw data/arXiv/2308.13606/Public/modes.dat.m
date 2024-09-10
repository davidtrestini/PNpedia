(* ::Package:: *)

(* ::Chapter:: *)
(*Supplementary file to "Spin effects in gravitational waveforms and fluxes for binaries on eccentric orbits to the third post-Newtonian order"*)


(* ::Text:: *)
(*by Quentin Henry and Mohammed Khalil*)


(* ::Text:: *)
(*This file contains the waveform modes for eccentric orbits and aligned spins, in harmonic coordinates using the covariant spin-supplementary condition.*)
(**)
(*The spin contributions are included to 3PN, while the nonspinning contributions are mostly included to 2PN, except in the (2,1) mode, which includes the 2.5PN nonspinning part, and the oscillatory memory contributions, which include the 3PN nonspinning part.*)
(**)
(*The full 3PN nonspinning modes for eccentric orbits are provided in the Supplemental Material of the following references: [Mishra, Arun and Iyer arXiv:1501.07096], [Boetzel et al. arXiv: 1904.11814], [Ebersold et al. DOI: 10.1103/PhysRevD.100.084043]*)


(*
Notation:
We use geometric units with c = G = 1
\[Nu] = m1 m2 / M^2 is the symmetric mass ratio 
\[Delta] = (m1 - m2) / M is the antisymmetric mass ratio 
M = m1 + m2 is the total mass

\[Epsilon] = 1 is a PN counting parameter
SO = 1 is a spin order counting parameter

x = (M \[CapitalOmega])^(2/3), where \[CapitalOmega] is the azimuthal orbital frequency
e = e_t is the time eccentricity
\[ScriptL] is the mean anomaly
b is an arbitrary gauge parameter in the tail part
\[Gamma]E = 0.577... is the EulerGamma constant
R is the distance from the observer to the source
rd and \[Phi]d denote dr/dt and d\[Phi]/dt, repsectively

\[Chi]1 = S1 / m1^2, and \[Chi]2 = S2 / m2^2 are the dimensionless spins
\[Chi]S = 1/2 (\[Chi]1 + \[Chi]2)
\[Chi]A = 1/2 (\[Chi]1 - \[Chi]2)
\[Kappa]1 and \[Kappa]2 are the spin quadrupole constants, which equal 1 for black holes
\[Kappa]S = 1/2 ((\[Kappa]1 - 1) \[Chi]1^2 + (\[Kappa]2 - 1) \[Chi]2^2)
\[Kappa]A = 1/2 ((\[Kappa]1 - 1) \[Chi]1^2 - (\[Kappa]2 - 1) \[Chi]2^2)

The modes are factored as
h[l, m] = 8 M \[Nu] / R Sqrt[\[Pi]/5] E^(-i m \[Phi]) H[l, m]
with H[l, m] split into instantaneous, tail, oscillatory memory, and DC memory contributions
H[l, m] = Hinst[l, m] + Htail[l, m] + HoscMem[l, m] + HDCMem[l, m]
*)


(* ::Section::Closed:: *)
(*Instantaneous contributions*)


Hinst[2, 0] = (-M + r*rd^2 + r^3*\[Phi]d^2)/(Sqrt[6]*r) + 
    (\[Epsilon]^2*(-7*M^2*(-10 + \[Nu]) + M*r*rd^2*(15 + 32*\[Nu]) + 
       M*r^3*(-37 + 20*\[Nu])*\[Phi]d^2 - 9*r^2*(-1 + 3*\[Nu])*
        (rd^2 + r^2*\[Phi]d^2)^2))/(14*Sqrt[6]*r^2) - 
    (\[Epsilon]^4*(M^3*(6056 + 2534*\[Nu] + 316*\[Nu]^2) - 
       3*r^3*(83 - 589*\[Nu] + 1111*\[Nu]^2)*(rd^2 + r^2*\[Phi]d^2)^3 + 
       M^2*(-2*r*rd^2*(-619 + 2789*\[Nu] + 934*\[Nu]^2) + 
         r^3*(283 + 4075*\[Nu] - 1165*\[Nu]^2)*\[Phi]d^2) + 
       3*M*(r^2*rd^4*(-557 + 664*\[Nu] + 1712*\[Nu]^2) + 
         6*r^4*rd^2*(-113 - 42*\[Nu] + 488*\[Nu]^2)*\[Phi]d^2 + 
         2*r^6*(-169 - 37*\[Nu] + 577*\[Nu]^2)*\[Phi]d^4)))/
     (504*Sqrt[6]*r^3) + SO*((Sqrt[2/3]*M^2*\[Epsilon]^3*\[Phi]d*
        (\[Delta]*\[Chi]A + \[Chi]S + \[Nu]*\[Chi]S))/r - 
      (M^2*\[Epsilon]^5*\[Phi]d*(M*(\[Delta]*(86 - 31*\[Nu])*\[Chi]A + 
           (86 + 369*\[Nu] + 6*\[Nu]^2)*\[Chi]S) + r^3*\[Phi]d^2*
          (3*\[Delta]*(12 + 5*\[Nu])*\[Chi]A + (36 - 263*\[Nu] + 90*\[Nu]^2)*
            \[Chi]S) + r*rd^2*(-3*\[Delta]*(14 + 5*\[Nu])*\[Chi]A + 
           (-42 + 295*\[Nu] + 180*\[Nu]^2)*\[Chi]S)))/(14*Sqrt[6]*r^2)) + 
    SO^2*(-(Sqrt[3/2]*M^3*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
          \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
          \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/(2*r^3) + 
      (M^3*\[Epsilon]^6*(2*M*(104*\[Kappa]S - 143*\[Kappa]S*\[Nu] + 
           36*\[Kappa]S*\[Nu]^2 + 13*\[Delta]*\[Kappa]A*(8 + 5*\[Nu]) + 
           104*\[Chi]A^2 - 526*\[Nu]*\[Chi]A^2 + 72*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(104 + 51*\[Nu])*\[Chi]A*\[Chi]S + 104*\[Chi]S^2 + 
           212*\[Nu]*\[Chi]S^2 + 56*\[Nu]^2*\[Chi]S^2) - 
         r^3*\[Phi]d^2*(\[Delta]*\[Kappa]A*(11 - 14*\[Nu]) + 
           \[Kappa]S*(11 - 36*\[Nu] + 52*\[Nu]^2) + 11*\[Chi]A^2 - 
           268*\[Nu]*\[Chi]A^2 + 104*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(11 + 42*\[Nu])*\[Chi]A*\[Chi]S + 11*\[Chi]S^2 + 
           308*\[Nu]*\[Chi]S^2 + 112*\[Nu]^2*\[Chi]S^2) + 
         r*rd^2*(\[Delta]*\[Kappa]A*(-71 + 20*\[Nu]) + 
           \[Kappa]S*(-71 + 162*\[Nu] + 8*\[Nu]^2) - 71*\[Chi]A^2 - 
           116*\[Nu]*\[Chi]A^2 + 16*\[Nu]^2*\[Chi]A^2 - 
           2*\[Delta]*(71 + 36*\[Nu])*\[Chi]A*\[Chi]S - 71*\[Chi]S^2 + 
           328*\[Nu]*\[Chi]S^2 + 224*\[Nu]^2*\[Chi]S^2)))/(28*Sqrt[6]*r^4))
 
Hinst[2, 1] = (I/3)*M*\[Delta]*\[Epsilon]*\[Phi]d + 
    ((I/84)*M*\[Delta]*\[Epsilon]^3*\[Phi]d*(4*M*(-9 + 11*\[Nu]) + 
       r*(rd^2*(-66 + 20*\[Nu]) + (2*I)*r*rd*(83 + 2*\[Nu])*\[Phi]d + 
         r^2*(19 - 24*\[Nu])*\[Phi]d^2)))/r + 
    ((I/1512)*M*\[Delta]*\[Epsilon]^5*\[Phi]d*
      (10*M^2*(31 - 205*\[Nu] + 111*\[Nu]^2) - 
       2*M*r*(-8*rd^2*(-202 - 587*\[Nu] + 177*\[Nu]^2) + 
         I*r*rd*(-3167 - 5278*\[Nu] + 201*\[Nu]^2)*\[Phi]d + 
         r^2*(-197 + 5*\[Nu] + 660*\[Nu]^2)*\[Phi]d^2) + 
       3*r^2*(rd^4*(-241 + 550*\[Nu] - 264*\[Nu]^2) - 
         (2*I)*r*rd^3*(-265 + 526*\[Nu] + 18*\[Nu]^2)*\[Phi]d - 
         3*r^2*rd^2*(75 - 560*\[Nu] + 68*\[Nu]^2)*\[Phi]d^2 + 
         (2*I)*r^3*rd*(308 - 1607*\[Nu] + 111*\[Nu]^2)*\[Phi]d^3 + 
         r^4*(152 - 692*\[Nu] + 333*\[Nu]^2)*\[Phi]d^4)))/r^2 - 
    (((3*I)/4)*M^4*SO^3*\[Epsilon]^6*(\[Delta]*\[Kappa]A*\[Chi]A + 
       \[Chi]A^3 - 4*\[Nu]*\[Chi]A^3 + \[Kappa]S*(\[Chi]A - 
         2*\[Nu]*\[Chi]A) + \[Kappa]A*\[Chi]S - 4*\[Kappa]A*\[Nu]*\[Chi]S + 
       3*\[Chi]A*\[Chi]S^2 - 8*\[Nu]*\[Chi]A*\[Chi]S^2 + 
       \[Delta]*\[Chi]S*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
         (3 - 4*\[Nu])*\[Chi]A^2 + \[Chi]S^2)))/r^4 + 
    (M^3*SO^2*\[Epsilon]^5*((-I)*r*\[Phi]d*(\[Kappa]A*(-1 + 8*\[Nu]) + 
         4*(-5 + 13*\[Nu])*\[Chi]A*\[Chi]S + \[Delta]*
          (\[Kappa]S*(-1 + 6*\[Nu]) + 2*(-5 + 6*\[Nu])*\[Chi]A^2 - 
           10*\[Chi]S^2)) - 4*rd*(\[Kappa]A - 2*\[Kappa]A*\[Nu] + 
         2*(1 - 2*\[Nu])*\[Chi]A*\[Chi]S + \[Delta]*(\[Kappa]S + \[Chi]A^2 + 
           \[Chi]S^2))))/(6*r^3) + 
    SO*(((-I/2)*M^2*\[Epsilon]^2*(\[Chi]A + \[Delta]*\[Chi]S))/r^2 + 
      ((I/84)*M^2*\[Epsilon]^4*(2*M*(77 + 59*\[Nu])*\[Chi]A + 
         22*M*\[Delta]*(7 + \[Nu])*\[Chi]S + 
         r*rd^2*((105 - 52*\[Nu])*\[Chi]A + 15*\[Delta]*(7 - 4*\[Nu])*
            \[Chi]S) + 4*r^3*\[Phi]d^2*(3*(-7 + 22*\[Nu])*\[Chi]A + 
           \[Delta]*(-21 + 4*\[Nu])*\[Chi]S) + (2*I)*r^2*rd*\[Phi]d*
          ((-147 + 83*\[Nu])*\[Chi]A - \[Delta]*(147 + 13*\[Nu])*\[Chi]S)))/
       r^3 + ((I/3024)*M^2*\[Epsilon]^6*
        (2*M^2*((-5661 - 17156*\[Nu] + 231*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(-5661 - 9172*\[Nu] + 775*\[Nu]^2)*\[Chi]S) + 
         2*M*r*(r^2*\[Phi]d^2*((765 + 667*\[Nu] + 7773*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(765 + 3731*\[Nu] - 1715*\[Nu]^2)*\[Chi]S) + 
           rd^2*((936 + 8960*\[Nu] - 1932*\[Nu]^2)*\[Chi]A + 
             4*\[Delta]*(234 + 3604*\[Nu] - 1111*\[Nu]^2)*\[Chi]S) + 
           I*r*rd*\[Phi]d*((2043 + 37*\[Nu] + 10635*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(2043 + 2597*\[Nu] + 139*\[Nu]^2)*\[Chi]S)) - 
         3*r^2*(4*r^4*\[Phi]d^4*((252 - 1315*\[Nu] + 1140*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(252 - 857*\[Nu] - 172*\[Nu]^2)*\[Chi]S) - 
           (4*I)*r*rd^3*\[Phi]d*((-315 + 502*\[Nu] + 60*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(-315 + 926*\[Nu] + 4*\[Nu]^2)*\[Chi]S) - 
           12*r^2*rd^2*\[Phi]d^2*((189 - 1042*\[Nu] + 385*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(189 - 586*\[Nu] + 161*\[Nu]^2)*\[Chi]S) + 
           (4*I)*r^3*rd*\[Phi]d^3*((936 - 4895*\[Nu] + 2706*\[Nu]^2)*
              \[Chi]A + \[Delta]*(936 - 1075*\[Nu] + 586*\[Nu]^2)*\[Chi]S) - 
           rd^4*((567 - 1232*\[Nu] + 1416*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(567 - 1024*\[Nu] + 1384*\[Nu]^2)*\[Chi]S))))/r^4)
 
Hinst[2, 2] = (M + r*(I*rd + r*\[Phi]d)^2)/(2*r) + 
    (\[Epsilon]^2*(21*M^2*(-10 + \[Nu]) - 27*r^2*(-1 + 3*\[Nu])*
        ((-I)*rd + r*\[Phi]d)*(I*rd + r*\[Phi]d)^3 + 
       M*r*(-3*rd^2*(15 + 32*\[Nu]) + (10*I)*r*rd*(5 + 27*\[Nu])*\[Phi]d + 
         r^2*(11 + 156*\[Nu])*\[Phi]d^2)))/(84*r^2) + 
    (\[Epsilon]^4*(6*M^3*(3028 + 1267*\[Nu] + 158*\[Nu]^2) + 
       9*r^3*(83 - 589*\[Nu] + 1111*\[Nu]^2)*(rd - I*r*\[Phi]d)^4*
        ((-I)*rd + r*\[Phi]d)^2 + 
       M^2*r*(-6*rd^2*(-619 + 2789*\[Nu] + 934*\[Nu]^2) + 
         (8*I)*r*rd*(-773 - 3767*\[Nu] + 2852*\[Nu]^2)*\[Phi]d + 
         r^2*(-11891 - 36575*\[Nu] + 13133*\[Nu]^2)*\[Phi]d^2) - 
       3*M*r^2*(-3*rd^4*(-557 + 664*\[Nu] + 1712*\[Nu]^2) + 
         (4*I)*r*rd^3*(-863 + 1462*\[Nu] + 2954*\[Nu]^2)*\[Phi]d + 
         6*r^2*rd^2*(-33 + 1014*\[Nu] + 232*\[Nu]^2)*\[Phi]d^2 + 
         (6*I)*r^3*rd*(-433 - 721*\[Nu] + 1703*\[Nu]^2)*\[Phi]d^3 + 
         2*r^4*(-835 - 19*\[Nu] + 2995*\[Nu]^2)*\[Phi]d^4)))/(3024*r^3) + 
    SO*(((-I/6)*M^3*\[Epsilon]^6*(M + r*(2*rd^2 - (10*I)*r*rd*\[Phi]d + 
           7*r^2*\[Phi]d^2))*(\[Delta]*\[Chi]A + \[Chi]S - 2*\[Nu]*\[Chi]S))/
       r^4 + (M^2*\[Epsilon]^3*((-I)*rd*(3*\[Delta]*\[Chi]A + 
           (3 - 8*\[Nu])*\[Chi]S) + r*\[Phi]d*(-3*\[Delta]*\[Chi]A + 
           (-3 + 5*\[Nu])*\[Chi]S)))/(3*r^2) + 
      (M^2*\[Epsilon]^5*((-8*I)*M*rd*(\[Delta]*(-55 + 19*\[Nu])*\[Chi]A + 
           (-55 + 100*\[Nu] - 86*\[Nu]^2)*\[Chi]S) + (2*I)*r*rd^3*
          (3*\[Delta]*(-9 + 14*\[Nu])*\[Chi]A + (-27 + 30*\[Nu] - 4*\[Nu]^2)*
            \[Chi]S) - r^4*\[Phi]d^3*(\[Delta]*(120 + 83*\[Nu])*\[Chi]A + 
           3*(40 - 161*\[Nu] + 78*\[Nu]^2)*\[Chi]S) + r^2*rd^2*\[Phi]d*
          (\[Delta]*(18 + 275*\[Nu])*\[Chi]A + (18 - 315*\[Nu] + 188*\[Nu]^2)*
            \[Chi]S) - (2*I)*r^3*rd*\[Phi]d^2*(\[Delta]*(51 + 62*\[Nu])*
            \[Chi]A + (51 - 264*\[Nu] + 440*\[Nu]^2)*\[Chi]S) + 
         M*r*\[Phi]d*(\[Delta]*(238 - 141*\[Nu])*\[Chi]A + 
           (238 - 181*\[Nu] + 474*\[Nu]^2)*\[Chi]S)))/(84*r^3)) + 
    SO^2*((3*M^3*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
         4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
           2*\[Chi]A*\[Chi]S)))/(4*r^3) - 
      (M^3*\[Epsilon]^6*(2*M*(\[Delta]*\[Kappa]A*(438 + 103*\[Nu]) + 
           \[Kappa]S*(438 - 773*\[Nu] + 108*\[Nu]^2) + 
           2*(219 - 863*\[Nu] + 108*\[Nu]^2)*\[Chi]A^2 + 
           2*\[Delta]*(438 + 145*\[Nu])*\[Chi]A*\[Chi]S + 
           6*(73 + 44*\[Nu] - 28*\[Nu]^2)*\[Chi]S^2) + 
         r*(-(r^2*\[Phi]d^2*(\[Delta]*\[Kappa]A*(153 + 430*\[Nu]) + 
              \[Kappa]S*(153 + 124*\[Nu] - 804*\[Nu]^2) + 153*\[Chi]A^2 - 
              28*\[Nu]*\[Chi]A^2 - 1608*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(153 + 
                38*\[Nu])*\[Chi]A*\[Chi]S + 153*\[Chi]S^2 - 508*\[Nu]*\[Chi]S^
                2 - 112*\[Nu]^2*\[Chi]S^2)) - (2*I)*r*rd*\[Phi]d*
            (\[Delta]*\[Kappa]A*(117 + 211*\[Nu]) + \[Kappa]S*
              (117 - 23*\[Nu] - 390*\[Nu]^2) + 117*\[Chi]A^2 - 
             439*\[Nu]*\[Chi]A^2 - 780*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(117 - 223*\[Nu])*\[Chi]A*\[Chi]S + 117*\[Chi]S^2 - 
             475*\[Nu]*\[Chi]S^2 + 56*\[Nu]^2*\[Chi]S^2) + 
           rd^2*(\[Delta]*\[Kappa]A*(291 - 308*\[Nu]) + \[Kappa]S*
              (291 - 890*\[Nu] + 24*\[Nu]^2) + 291*\[Chi]A^2 - 
             940*\[Nu]*\[Chi]A^2 + 48*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(291 - 364*\[Nu])*\[Chi]A*\[Chi]S + 291*\[Chi]S^2 - 
             952*\[Nu]*\[Chi]S^2 + 224*\[Nu]^2*\[Chi]S^2))))/(168*r^4))
 
Hinst[3, 0] = ((I/2)*M*rd*\[Epsilon]^2*(-1 + 3*\[Nu])*\[Phi]d)/Sqrt[42] + 
    ((I/36)*M*rd*\[Epsilon]^4*\[Phi]d*(2*M*(-82 + 239*\[Nu] + 55*\[Nu]^2) - 
       3*(r*rd^2*(-13 + 38*\[Nu] + 10*\[Nu]^2) + 
         4*r^3*(-23 + 70*\[Nu] + 2*\[Nu]^2)*\[Phi]d^2)))/(Sqrt[42]*r) + 
    ((I/4)*Sqrt[3/14]*M^3*rd*SO^2*\[Epsilon]^6*\[Phi]d*
      (\[Kappa]S - 5*\[Kappa]S*\[Nu] + 6*\[Kappa]S*\[Nu]^2 + \[Chi]A^2 - 
       7*\[Nu]*\[Chi]A^2 + 12*\[Nu]^2*\[Chi]A^2 + \[Chi]S^2 - 
       19*\[Nu]*\[Chi]S^2 + 16*\[Nu]^2*\[Chi]S^2 + 
       \[Delta]*(\[Kappa]A - 3*\[Kappa]A*\[Nu] + 2*(1 - 11*\[Nu])*\[Chi]A*
          \[Chi]S)))/r^2 + 
    SO*(((-I)*Sqrt[2/21]*M^2*rd*\[Epsilon]^3*\[Nu]*\[Chi]S)/r^2 - 
      (I*M^3*\[Epsilon]^6*(M + 2*r*rd^2 - r^3*\[Phi]d^2)*
        (\[Delta]*\[Chi]A + \[Chi]S - 2*\[Nu]*\[Chi]S))/(Sqrt[42]*r^4) - 
      ((I/12)*M^2*rd*\[Epsilon]^5*(12*r*rd^2*\[Nu]*(\[Delta]*\[Chi]A - 
           4*\[Chi]S) + 2*M*(\[Delta]*(-12 + 5*\[Nu])*\[Chi]A + 
           (-12 + 97*\[Nu] + 4*\[Nu]^2)*\[Chi]S) + 3*r^3*\[Phi]d^2*
          (\[Delta]*(18 - 31*\[Nu])*\[Chi]A + (18 - 179*\[Nu] + 32*\[Nu]^2)*
            \[Chi]S)))/(Sqrt[42]*r^3))
 
Hinst[3, 1] = (\[Delta]*\[Epsilon]*(6*r*(rd - I*r*\[Phi]d)^2*
        (rd + I*r*\[Phi]d) + M*(-12*rd + (7*I)*r*\[Phi]d)))/(12*Sqrt[14]*r) + 
    (\[Delta]*\[Epsilon]^3*((6*I)*r^2*(-5 + 19*\[Nu])*((-I)*rd + r*\[Phi]d)^2*
        (I*rd + r*\[Phi]d)^3 + M^2*(rd*(218 - 172*\[Nu]) + 
         (2*I)*r*(-101 + 43*\[Nu])*\[Phi]d) + 
       3*M*r*(4*rd^3*(8 + 17*\[Nu]) - I*r*rd^2*(33 + 62*\[Nu])*\[Phi]d + 
         6*r^2*rd*(2 + 9*\[Nu])*\[Phi]d^2 - (4*I)*r^3*(-9 + 14*\[Nu])*
          \[Phi]d^3)))/(72*Sqrt[14]*r^2) + 
    (M^3*SO^2*\[Epsilon]^5*(2*rd*(\[Kappa]A*(-17 + 52*\[Nu]) + 
         2*(-17 + 52*\[Nu])*\[Chi]A*\[Chi]S + \[Delta]*
          (\[Kappa]S*(-17 + 18*\[Nu]) + (-17 + 36*\[Nu])*\[Chi]A^2 - 
           17*\[Chi]S^2)) - I*r*\[Phi]d*(\[Kappa]A*(-13 + 68*\[Nu]) + 
         2*(-13 + 68*\[Nu])*\[Chi]A*\[Chi]S + \[Delta]*
          (\[Kappa]S*(-13 + 42*\[Nu]) + (-13 + 84*\[Nu])*\[Chi]A^2 - 
           13*\[Chi]S^2))))/(24*Sqrt[14]*r^3) + 
    SO*(((I/24)*M^2*\[Epsilon]^4*(4*M*(-1 + 5*\[Nu])*(\[Chi]A + 
           \[Delta]*\[Chi]S) + r*(-2*rd^2*((-6 + 25*\[Nu])*\[Chi]A + 
             3*\[Delta]*(-2 + 5*\[Nu])*\[Chi]S) + r^2*\[Phi]d^2*
            (3*(-8 + 29*\[Nu])*\[Chi]A - \[Delta]*(24 + 31*\[Nu])*\[Chi]S) + 
           (2*I)*r*rd*\[Phi]d*((-6 + 31*\[Nu])*\[Chi]A - 
             \[Delta]*(6 + 35*\[Nu])*\[Chi]S))))/(Sqrt[14]*r^3) + 
      ((I/432)*M^2*\[Epsilon]^6*((-6*I)*r^3*rd^3*\[Phi]d*
          ((-186 + 745*\[Nu] + 354*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(-186 + 191*\[Nu] - 554*\[Nu]^2)*\[Chi]S) - 
         6*r^6*\[Phi]d^4*((-90 + 28*\[Nu] + 579*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(-90 + 800*\[Nu] - 551*\[Nu]^2)*\[Chi]S) + 
         6*r^2*rd^4*((30 - 187*\[Nu] + 318*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(30 - 101*\[Nu] + 122*\[Nu]^2)*\[Chi]S) + 
         2*M^2*((252 - 1277*\[Nu] + 96*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(252 - 1279*\[Nu] + 376*\[Nu]^2)*\[Chi]S) + 
         9*r^4*rd^2*\[Phi]d^2*((32 - 451*\[Nu] + 230*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(32 + 691*\[Nu] + 626*\[Nu]^2)*\[Chi]S) + 
         (6*I)*r^5*rd*\[Phi]d^3*((-12 - 341*\[Nu] + 315*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(-12 - 91*\[Nu] + 1213*\[Nu]^2)*\[Chi]S) + 
         2*M*r*(I*r*rd*\[Phi]d*((-1104 + 3073*\[Nu] + 4179*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(-1104 + 3071*\[Nu] - 3083*\[Nu]^2)*\[Chi]S) + 
           2*r^2*\[Phi]d^2*((312 - 1654*\[Nu] + 1005*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(312 + 1846*\[Nu] - 655*\[Nu]^2)*\[Chi]S) - 
           2*rd^2*((105 - 541*\[Nu] + 924*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(105 - 77*\[Nu] + 494*\[Nu]^2)*\[Chi]S))))/
       (Sqrt[14]*r^4))
 
Hinst[3, 2] = -(Sqrt[5/7]*M*\[Epsilon]^2*(-1 + 3*\[Nu])*\[Phi]d*
       (I*rd + 4*r*\[Phi]d))/12 - (M*\[Epsilon]^4*\[Phi]d*
      (2*M*((5*I)*rd*(-82 + 239*\[Nu] + 55*\[Nu]^2) + 
         r*(167 - 925*\[Nu] + 1615*\[Nu]^2)*\[Phi]d) - 
       3*r*((5*I)*rd^3*(-13 + 38*\[Nu] + 10*\[Nu]^2) + 
         12*r*rd^2*(-23 + 70*\[Nu] + 10*\[Nu]^2)*\[Phi]d - 
         (60*I)*r^2*rd*(-8 + 25*\[Nu] + \[Nu]^2)*\[Phi]d^2 + 
         2*r^3*(-13 - 25*\[Nu] + 355*\[Nu]^2)*\[Phi]d^3)))/(216*Sqrt[35]*r) + 
    SO*((Sqrt[5/7]*M^2*\[Epsilon]^3*\[Nu]*(I*rd + 4*r*\[Phi]d)*\[Chi]S)/
       (3*r^2) + ((I/6)*Sqrt[5/7]*M^3*\[Epsilon]^6*
        (M + r*(2*rd^2 - (10*I)*r*rd*\[Phi]d + 7*r^2*\[Phi]d^2))*
        (\[Delta]*\[Chi]A + \[Chi]S - 2*\[Nu]*\[Chi]S))/r^4 + 
      ((I/72)*Sqrt[5/7]*M^2*\[Epsilon]^5*
        (2*M*(rd*(\[Delta]*(-12 + 5*\[Nu])*\[Chi]A + (-12 + 97*\[Nu] + 4*
                \[Nu]^2)*\[Chi]S) - (2*I)*r*\[Phi]d*
            (\[Delta]*(-12 + 23*\[Nu])*\[Chi]A - (12 + 53*\[Nu] + 8*\[Nu]^2)*
              \[Chi]S)) + 3*r*(4*rd^3*\[Nu]*(\[Delta]*\[Chi]A - 4*\[Chi]S) + 
           (4*I)*r*rd^2*\[Nu]*\[Phi]d*(-5*\[Delta]*\[Chi]A + 
             (15 + 4*\[Nu])*\[Chi]S) + (4*I)*r^3*\[Nu]*\[Phi]d^3*
            (-17*\[Delta]*\[Chi]A + (-5 + 16*\[Nu])*\[Chi]S) + 
           r^2*rd*\[Phi]d^2*(\[Delta]*(-30 + 17*\[Nu])*\[Chi]A + 
             (-30 + 189*\[Nu] + 32*\[Nu]^2)*\[Chi]S))))/r^3) + 
    (Sqrt[5/7]*M^3*SO^2*\[Epsilon]^6*
      (2*M*\[Nu]*(\[Kappa]S - \[Chi]A^2 - 9*\[Chi]S^2 + 24*\[Nu]*\[Chi]S^2 + 
         \[Delta]*(\[Kappa]A - 10*\[Chi]A*\[Chi]S)) - 
       r*(2*r^2*\[Phi]d^2*(\[Delta]*\[Kappa]A*(-6 + 11*\[Nu]) + 
           \[Kappa]S*(-6 + 23*\[Nu] - 36*\[Nu]^2) - 6*\[Chi]A^2 + 
           49*\[Nu]*\[Chi]A^2 - 72*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(-6 + 37*\[Nu])*\[Chi]A*\[Chi]S - 6*\[Chi]S^2 + 
           49*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2) + 
         I*r*rd*\[Phi]d*(\[Delta]*\[Kappa]A*(3 + 11*\[Nu]) + 
           \[Kappa]S*(3 + 5*\[Nu] + 18*\[Nu]^2) + 3*\[Chi]A^2 - 
           41*\[Nu]*\[Chi]A^2 + 36*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(3 - 5*\[Nu])*\[Chi]A*\[Chi]S + 3*\[Chi]S^2 + 
           19*\[Nu]*\[Chi]S^2 + 16*\[Nu]^2*\[Chi]S^2) - 
         4*rd^2*\[Nu]*(\[Kappa]S - \[Chi]A^2 - \[Chi]S^2 + 
           8*\[Nu]*\[Chi]S^2 + \[Delta]*(\[Kappa]A - 2*\[Chi]A*\[Chi]S)))))/
     (24*r^4)
 
Hinst[3, 3] = (Sqrt[5/42]*\[Delta]*\[Epsilon]*(2*r*(-rd + I*r*\[Phi]d)^3 + 
       M*(4*rd - (7*I)*r*\[Phi]d)))/(4*r) + 
    (Sqrt[5/42]*\[Delta]*\[Epsilon]^3*(6*r^2*(-5 + 19*\[Nu])*
        (rd - I*r*\[Phi]d)^4*(rd + I*r*\[Phi]d) + 
       2*M^2*(rd*(-109 + 86*\[Nu]) + (3*I)*r*(101 - 43*\[Nu])*\[Phi]d) + 
       3*M*r*(-4*rd^3*(8 + 17*\[Nu]) + (3*I)*r*rd^2*(33 + 62*\[Nu])*\[Phi]d + 
         6*r^2*rd*(14 + 31*\[Nu])*\[Phi]d^2 - (12*I)*r^3*(1 + 4*\[Nu])*
          \[Phi]d^3)))/(72*r^2) - (Sqrt[15/14]*M^3*SO^2*\[Epsilon]^5*
      (2*rd - (7*I)*r*\[Phi]d)*(\[Kappa]A*(-1 + 4*\[Nu]) + 
       2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S + 
       \[Delta]*(\[Kappa]S*(-1 + 2*\[Nu]) + (-1 + 4*\[Nu])*\[Chi]A^2 - 
         \[Chi]S^2)))/(8*r^3) + 
    SO*(((-I/8)*Sqrt[5/42]*M^2*\[Epsilon]^4*
        (4*M*(-1 + 5*\[Nu])*(\[Chi]A + \[Delta]*\[Chi]S) + 
         r*(-2*rd^2*((-6 + 25*\[Nu])*\[Chi]A + 3*\[Delta]*(-2 + 5*\[Nu])*
              \[Chi]S) + r^2*\[Phi]d^2*((-24 + 119*\[Nu])*\[Chi]A + 
             3*\[Delta]*(-8 + 11*\[Nu])*\[Chi]S) + (2*I)*r*rd*\[Phi]d*
            ((-18 + 77*\[Nu])*\[Chi]A + 3*\[Delta]*(-6 + 13*\[Nu])*
              \[Chi]S))))/r^3 - ((I/144)*M^2*\[Epsilon]^6*
        (10*M^2*((252 - 1277*\[Nu] + 96*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(252 - 1279*\[Nu] + 376*\[Nu]^2)*\[Chi]S) - 
         3*r^2*(3*r^2*rd^2*\[Phi]d^2*((-480 + 1711*\[Nu] + 2322*\[Nu]^2)*
              \[Chi]A + \[Delta]*(-480 + 1889*\[Nu] - 1514*\[Nu]^2)*
              \[Chi]S) - 10*rd^4*((30 - 187*\[Nu] + 318*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(30 - 101*\[Nu] + 122*\[Nu]^2)*\[Chi]S) + 
           (2*I)*r*rd^3*\[Phi]d*((90 - 1321*\[Nu] + 5118*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(90 - 319*\[Nu] + 714*\[Nu]^2)*\[Chi]S) + 
           2*r^4*\[Phi]d^4*((350 - 1616*\[Nu] + 883*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(350 - 1844*\[Nu] + 769*\[Nu]^2)*\[Chi]S) + 
           (2*I)*r^3*rd*\[Phi]d^3*((660 - 4061*\[Nu] + 2643*\[Nu]^2)*
              \[Chi]A + \[Delta]*(660 - 2899*\[Nu] + 4789*\[Nu]^2)*
              \[Chi]S)) + 2*M*r*(-10*rd^2*((105 - 541*\[Nu] + 924*\[Nu]^2)*
              \[Chi]A + \[Delta]*(105 - 77*\[Nu] + 494*\[Nu]^2)*\[Chi]S) + 
           2*r^2*\[Phi]d^2*((1320 - 8938*\[Nu] + 8709*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(1320 - 422*\[Nu] + 2777*\[Nu]^2)*\[Chi]S) + 
           (3*I)*r*rd*\[Phi]d*((2000 - 9147*\[Nu] + 8911*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(2000 - 3173*\[Nu] + 5273*\[Nu]^2)*\[Chi]S))))/
       (Sqrt[210]*r^4))
 
Hinst[4, 0] = (\[Epsilon]^2*(-1 + 3*\[Nu])*
      (7*M^2 + 6*r^2*(rd^2 + r^2*\[Phi]d^2)^2 - 
       M*(18*r*rd^2 + 13*r^3*\[Phi]d^2)))/(84*Sqrt[2]*r^2) + 
    (\[Epsilon]^4*(8*M^3*(314 - 987*\[Nu] + 195*\[Nu]^2) - 
       12*r^3*(23 - 159*\[Nu] + 291*\[Nu]^2)*(rd^2 + r^2*\[Phi]d^2)^3 - 
       M^2*(r*rd^2*(580 - 4066*\[Nu] + 8730*\[Nu]^2) + 
         r^3*(3013 - 10324*\[Nu] + 5044*\[Nu]^2)*\[Phi]d^2) + 
       M*(12*r^2*rd^4*(-83 + 30*\[Nu] + 762*\[Nu]^2) + 
         2*r^4*rd^2*(-789 + 84*\[Nu] + 8012*\[Nu]^2)*\[Phi]d^2 + 
         r^6*(909 - 4776*\[Nu] + 7108*\[Nu]^2)*\[Phi]d^4)))/
     (3696*Sqrt[2]*r^3) - (M^2*SO*\[Epsilon]^5*\[Phi]d*
      (M*\[Delta]*(-10 + 29*\[Nu])*\[Chi]A + M*(-10 + 17*\[Nu] + 40*\[Nu]^2)*
        \[Chi]S - 2*r*rd^2*\[Nu]*(\[Delta]*\[Chi]A + (-15 + 44*\[Nu])*
          \[Chi]S) - r^3*\[Phi]d^2*(\[Delta]*(-9 + 26*\[Nu])*\[Chi]A + 
         (-9 + 15*\[Nu] + 37*\[Nu]^2)*\[Chi]S)))/(42*Sqrt[2]*r^2) + 
    (M^3*SO^2*\[Epsilon]^6*
      (-8*r*rd^2*(\[Kappa]S*(-3 + 17*\[Nu] - 18*\[Nu]^2) - 3*\[Chi]A^2 + 
         23*\[Nu]*\[Chi]A^2 - 36*\[Nu]^2*\[Chi]A^2 - 3*\[Chi]S^2 + 
         11*\[Nu]*\[Chi]S^2 + \[Delta]*(-3 + 11*\[Nu])*(\[Kappa]A + 
           2*\[Chi]A*\[Chi]S)) + 
       4*M*(\[Kappa]S*(-9 + 43*\[Nu] - 54*\[Nu]^2) - 9*\[Chi]A^2 + 
         61*\[Nu]*\[Chi]A^2 - 108*\[Nu]^2*\[Chi]A^2 - 9*\[Chi]S^2 + 
         25*\[Nu]*\[Chi]S^2 + \[Delta]*(-9 + 25*\[Nu])*(\[Kappa]A + 
           2*\[Chi]A*\[Chi]S)) + r^3*\[Phi]d^2*
        (\[Kappa]S*(33 - 157*\[Nu] + 198*\[Nu]^2) + 33*\[Chi]A^2 - 
         223*\[Nu]*\[Chi]A^2 + 396*\[Nu]^2*\[Chi]A^2 + 33*\[Chi]S^2 - 
         91*\[Nu]*\[Chi]S^2 - \[Delta]*(-33 + 91*\[Nu])*
          (\[Kappa]A + 2*\[Chi]A*\[Chi]S))))/(168*Sqrt[2]*r^4)
 
Hinst[4, 1] = (M*\[Delta]*\[Epsilon]^3*(-1 + 2*\[Nu])*\[Phi]d*
      ((-12*I)*M + r*((6*I)*rd^2 + 10*r*rd*\[Phi]d + (11*I)*r^2*\[Phi]d^2)))/
     (84*Sqrt[10]*r) + SO*((Sqrt[5/2]*M^2*\[Epsilon]^4*\[Nu]*
        ((-12*I)*M + r*((6*I)*rd^2 + 10*r*rd*\[Phi]d + (11*I)*r^2*\[Phi]d^2))*
        (\[Chi]A - \[Delta]*\[Chi]S))/(168*r^3) - 
      ((I/3696)*M^2*\[Epsilon]^6*(2*M^2*((440 - 6801*\[Nu] + 1428*\[Nu]^2)*
            \[Chi]A + \[Delta]*(440 + 3193*\[Nu] - 300*\[Nu]^2)*\[Chi]S) + 
         2*M*r*(2*r^2*\[Phi]d^2*((66 - 158*\[Nu] + 3389*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(66 - 238*\[Nu] - 301*\[Nu]^2)*\[Chi]S) + 
           I*r*rd*\[Phi]d*(5*(264 - 311*\[Nu] + 823*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(1320 - 9093*\[Nu] - 59*\[Nu]^2)*\[Chi]S) + 
           rd^2*((-440 - 3318*\[Nu] + 3024*\[Nu]^2)*\[Chi]A - 
             2*\[Delta]*(220 - 3067*\[Nu] + 240*\[Nu]^2)*\[Chi]S)) - 
         r^2*(6*rd^4*\[Nu]*((-301 + 298*\[Nu])*\[Chi]A + 3*\[Delta]*
              (71 + 50*\[Nu])*\[Chi]S) + 2*r^4*\[Nu]*\[Phi]d^4*
            ((-1396 + 2593*\[Nu])*\[Chi]A + \[Delta]*(1220 + 551*\[Nu])*
              \[Chi]S) - (2*I)*r*rd^3*\[Nu]*\[Phi]d*
            (25*(-65 + 74*\[Nu])*\[Chi]A + \[Delta]*(1009 + 982*\[Nu])*
              \[Chi]S) + r^2*rd^2*\[Nu]*\[Phi]d^2*((-4973 + 4634*\[Nu])*
              \[Chi]A + \[Delta]*(2773 + 4102*\[Nu])*\[Chi]S) + 
           (2*I)*r^3*rd*\[Phi]d^3*(5*(396 - 353*\[Nu] + 1763*\[Nu]^2)*
              \[Chi]A + \[Delta]*(1980 - 14779*\[Nu] + 2393*\[Nu]^2)*
              \[Chi]S))))/(Sqrt[10]*r^4))
 
Hinst[4, 2] = -(Sqrt[5]*\[Epsilon]^2*(-1 + 3*\[Nu])*
       (7*M^2 - 6*r^2*((-I)*rd + r*\[Phi]d)*(I*rd + r*\[Phi]d)^3 + 
        3*M*r*(-6*rd^2 + (9*I)*r*rd*\[Phi]d + r^2*\[Phi]d^2)))/(252*r^2) - 
    (\[Epsilon]^4*(40*M^3*(314 - 987*\[Nu] + 195*\[Nu]^2) + 
       60*r^3*(23 - 159*\[Nu] + 291*\[Nu]^2)*(rd - I*r*\[Phi]d)^4*
        ((-I)*rd + r*\[Phi]d)^2 + 
       M^2*r*(-10*rd^2*(290 - 2033*\[Nu] + 4365*\[Nu]^2) + 
         (12*I)*r*rd*(967 - 4615*\[Nu] + 5935*\[Nu]^2)*\[Phi]d + 
         r^2*(1987 - 11200*\[Nu] + 12960*\[Nu]^2)*\[Phi]d^2) - 
       3*M*r^2*(-20*rd^4*(-83 + 30*\[Nu] + 762*\[Nu]^2) + 
         (2*I)*r*rd^3*(-1853 + 1730*\[Nu] + 13230*\[Nu]^2)*\[Phi]d - 
         2*r^2*rd^2*(549 - 2140*\[Nu] + 2140*\[Nu]^2)*\[Phi]d^2 + 
         (4*I)*r^3*rd*(-454 - 315*\[Nu] + 5980*\[Nu]^2)*\[Phi]d^3 + 
         r^4*(1577 - 7940*\[Nu] + 9920*\[Nu]^2)*\[Phi]d^4)))/
     (11088*Sqrt[5]*r^3) + (M^2*SO*\[Epsilon]^5*
      ((2*I)*M*rd*(\[Delta]*(-100 + 333*\[Nu])*\[Chi]A + 
         (-100 + 577*\[Nu] - 864*\[Nu]^2)*\[Chi]S) + 
       4*M*r*\[Phi]d*(\[Delta]*(-70 + 237*\[Nu])*\[Chi]A + 
         (-70 + 253*\[Nu] - 156*\[Nu]^2)*\[Chi]S) - 
       3*r*((2*I)*rd^3*(\[Delta]*(-20 + 63*\[Nu])*\[Chi]A + 
           (-20 + 107*\[Nu] - 144*\[Nu]^2)*\[Chi]S) + 4*r^3*\[Phi]d^3*
          (\[Delta]*(5 + 2*\[Nu])*\[Chi]A + (5 - 7*\[Nu] - 41*\[Nu]^2)*
            \[Chi]S) + 4*r*rd^2*\[Phi]d*(\[Delta]*(-20 + 67*\[Nu])*\[Chi]A + 
           (-20 + 63*\[Nu] - 16*\[Nu]^2)*\[Chi]S) - I*r^2*rd*\[Phi]d^2*
          (\[Delta]*(-10 + 89*\[Nu])*\[Chi]A + (-10 - 339*\[Nu] + 
             1048*\[Nu]^2)*\[Chi]S))))/(504*Sqrt[5]*r^3) + 
    (Sqrt[5]*M^3*SO^2*\[Epsilon]^6*
      (-2*M*(\[Kappa]S*(-6 + 29*\[Nu] - 36*\[Nu]^2) - 6*\[Chi]A^2 + 
         41*\[Nu]*\[Chi]A^2 - 72*\[Nu]^2*\[Chi]A^2 - 6*\[Chi]S^2 + 
         17*\[Nu]*\[Chi]S^2 + \[Delta]*(-6 + 17*\[Nu])*(\[Kappa]A + 
           2*\[Chi]A*\[Chi]S)) + 
       r*(r^2*\[Phi]d^2*(\[Kappa]S*(5 - 11*\[Nu] + 30*\[Nu]^2) + 
           5*\[Chi]A^2 - 21*\[Nu]*\[Chi]A^2 + 60*\[Nu]^2*\[Chi]A^2 + 
           5*\[Chi]S^2 - \[Nu]*\[Chi]S^2 - \[Delta]*(-5 + \[Nu])*
            (\[Kappa]A + 2*\[Chi]A*\[Chi]S)) + 
         4*rd^2*(\[Kappa]S*(-2 + 11*\[Nu] - 12*\[Nu]^2) - 2*\[Chi]A^2 + 
           15*\[Nu]*\[Chi]A^2 - 24*\[Nu]^2*\[Chi]A^2 - 2*\[Chi]S^2 + 
           7*\[Nu]*\[Chi]S^2 + \[Delta]*(-2 + 7*\[Nu])*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S)) - I*r*rd*\[Phi]d*
          (\[Kappa]S*(-13 + 85*\[Nu] - 78*\[Nu]^2) - 13*\[Chi]A^2 + 
           111*\[Nu]*\[Chi]A^2 - 156*\[Nu]^2*\[Chi]A^2 - 13*\[Chi]S^2 + 
           59*\[Nu]*\[Chi]S^2 + \[Delta]*(-13 + 59*\[Nu])*
            (\[Kappa]A + 2*\[Chi]A*\[Chi]S)))))/(168*r^4)
 
Hinst[4, 3] = ((I/12)*M*\[Delta]*\[Epsilon]^3*(-1 + 2*\[Nu])*\[Phi]d*
      (4*M + r*(-2*rd^2 + (10*I)*r*rd*\[Phi]d + 23*r^2*\[Phi]d^2)))/
     (Sqrt[70]*r) + SO*(((I/24)*Sqrt[5/14]*M^2*\[Epsilon]^4*\[Nu]*
        (4*M + r*(-2*rd^2 + (10*I)*r*rd*\[Phi]d + 23*r^2*\[Phi]d^2))*
        (\[Chi]A - \[Delta]*\[Chi]S))/r^3 + 
      ((I/1584)*M^2*\[Epsilon]^6*(2*M^2*((440 - 6801*\[Nu] + 1428*\[Nu]^2)*
            \[Chi]A + \[Delta]*(440 + 3193*\[Nu] - 300*\[Nu]^2)*\[Chi]S) + 
         3*r^2*(-2*rd^4*\[Nu]*((-301 + 298*\[Nu])*\[Chi]A + 
             3*\[Delta]*(71 + 50*\[Nu])*\[Chi]S) + (2*I)*r*rd^3*\[Nu]*\[Phi]d*
            (25*(-65 + 74*\[Nu])*\[Chi]A + \[Delta]*(1009 + 982*\[Nu])*
              \[Chi]S) + 3*r^2*rd^2*\[Nu]*\[Phi]d^2*
            ((-2587 + 3526*\[Nu])*\[Chi]A + \[Delta]*(1267 + 1978*\[Nu])*
              \[Chi]S) + 2*r^4*\[Nu]*\[Phi]d^4*((-3304 + 11617*\[Nu])*
              \[Chi]A + \[Delta]*(1016 + 4847*\[Nu])*\[Chi]S) + 
           (2*I)*r^3*rd*\[Phi]d^3*(5*(308 - 727*\[Nu] + 421*\[Nu]^2)*
              \[Chi]A + \[Delta]*(1540 - 7981*\[Nu] - 753*\[Nu]^2)*
              \[Chi]S)) + 2*M*r*((3*I)*r*rd*\[Phi]d*
            (5*(264 - 311*\[Nu] + 823*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(1320 - 9093*\[Nu] - 59*\[Nu]^2)*\[Chi]S) + 
           rd^2*((-440 - 3318*\[Nu] + 3024*\[Nu]^2)*\[Chi]A - 
             2*\[Delta]*(220 - 3067*\[Nu] + 240*\[Nu]^2)*\[Chi]S) + 
           2*r^2*\[Phi]d^2*((1826 - 19530*\[Nu] + 20145*\[Nu]^2)*\[Chi]A + 
             \[Delta]*(1826 + 1534*\[Nu] + 567*\[Nu]^2)*\[Chi]S))))/
       (Sqrt[70]*r^4))
 
Hinst[4, 4] = (Sqrt[5/7]*\[Epsilon]^2*(-1 + 3*\[Nu])*
      (7*M^2 + 6*r^2*(rd - I*r*\[Phi]d)^4 + 
       3*M*r*(-6*rd^2 + (18*I)*r*rd*\[Phi]d + 17*r^2*\[Phi]d^2)))/(72*r^2) + 
    (\[Epsilon]^4*(40*M^3*(314 - 987*\[Nu] + 195*\[Nu]^2) - 
       60*r^3*(23 - 159*\[Nu] + 291*\[Nu]^2)*(rd - I*r*\[Phi]d)^5*
        (rd + I*r*\[Phi]d) + M^2*r*(-10*rd^2*(290 - 2033*\[Nu] + 
           4365*\[Nu]^2) + (24*I)*r*rd*(967 - 4615*\[Nu] + 5935*\[Nu]^2)*
          \[Phi]d + r^2*(53143 - 199660*\[Nu] + 127500*\[Nu]^2)*\[Phi]d^2) - 
       3*M*r^2*(-20*rd^4*(-83 + 30*\[Nu] + 762*\[Nu]^2) + 
         (4*I)*r*rd^3*(-1853 + 1730*\[Nu] + 13230*\[Nu]^2)*\[Phi]d + 
         2*r^2*rd^2*(-6141 + 8980*\[Nu] + 31500*\[Nu]^2)*\[Phi]d^2 - 
         (8*I)*r^3*rd*(-976 + 1745*\[Nu] + 3150*\[Nu]^2)*\[Phi]d^3 + 
         r^4*(613 - 920*\[Nu] + 6420*\[Nu]^2)*\[Phi]d^4)))/
     (3168*Sqrt[35]*r^3) + (M^2*SO*\[Epsilon]^5*
      ((-3*I)*r^3*rd*\[Phi]d^2*(\[Delta]*(-250 + 849*\[Nu])*\[Chi]A + 
         (-250 + 1221*\[Nu] - 1512*\[Nu]^2)*\[Chi]S) - 
       (2*I)*M*rd*(\[Delta]*(-100 + 333*\[Nu])*\[Chi]A + 
         (-100 + 577*\[Nu] - 864*\[Nu]^2)*\[Chi]S) - 
       6*r^4*\[Phi]d^3*(\[Delta]*(-65 + 282*\[Nu])*\[Chi]A + 
         (-65 + 263*\[Nu] - 291*\[Nu]^2)*\[Chi]S) + 12*r^2*rd^2*\[Phi]d*
        (\[Delta]*(-40 + 129*\[Nu])*\[Chi]A + (-40 + 201*\[Nu] - 252*\[Nu]^2)*
          \[Chi]S) + (6*I)*r*rd^3*(\[Delta]*(-20 + 63*\[Nu])*\[Chi]A + 
         (-20 + 107*\[Nu] - 144*\[Nu]^2)*\[Chi]S) + 
       2*M*r*\[Phi]d*(\[Delta]*(130 - 513*\[Nu])*\[Chi]A + 
         (130 - 757*\[Nu] + 1224*\[Nu]^2)*\[Chi]S)))/(72*Sqrt[35]*r^3) - 
    (Sqrt[5/7]*M^3*SO^2*\[Epsilon]^6*(-1 + 3*\[Nu])*
      (12*M + r*(-8*rd^2 + (26*I)*r*rd*\[Phi]d + 53*r^2*\[Phi]d^2))*
      (\[Kappa]S*(-1 + 2*\[Nu]) - \[Chi]A^2 + 4*\[Nu]*\[Chi]A^2 - \[Chi]S^2 - 
       \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/(48*r^4)
 
Hinst[5, 0] = ((-I/36)*M*rd*\[Epsilon]^4*(1 - 5*\[Nu] + 5*\[Nu]^2)*\[Phi]d*
      (22*M - 6*r*rd^2 - 21*r^3*\[Phi]d^2))/(Sqrt[462]*r) + 
    ((I/12)*M^2*rd*SO*\[Epsilon]^5*\[Nu]*(22*M - 6*r*rd^2 - 21*r^3*\[Phi]d^2)*
      (\[Delta]*\[Chi]A + (-1 + 2*\[Nu])*\[Chi]S))/(Sqrt[462]*r^3)
 
Hinst[5, 1] = (\[Delta]*\[Epsilon]^3*(-1 + 2*\[Nu])*
      (2*M^2*(205*rd - (86*I)*r*\[Phi]d) + 120*r^2*(-rd + I*r*\[Phi]d)^3*
        ((-I)*rd + r*\[Phi]d)^2 + (3*I)*M*r*((160*I)*rd^3 + 
         132*r*rd^2*\[Phi]d + (160*I)*r^2*rd*\[Phi]d^2 + 97*r^3*\[Phi]d^3)))/
     (288*Sqrt[385]*r^2) - ((I/432)*M^2*SO*\[Epsilon]^6*
      (M^2*((90 - 622*\[Nu] + 966*\[Nu]^2)*\[Chi]A + 
         2*\[Delta]*(45 - 262*\[Nu] + 385*\[Nu]^2)*\[Chi]S) - 
       3*r^2*(I*r^3*rd*\[Phi]d^3*(35*(3 - 22*\[Nu] + 36*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(105 + 494*\[Nu] - 1268*\[Nu]^2)*\[Chi]S) + 
         r^4*\[Phi]d^4*((180 - 1009*\[Nu] + 1227*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(180 - 95*\[Nu] - 601*\[Nu]^2)*\[Chi]S) + 
         (4*I)*r*rd^3*\[Phi]d*(5*(3 - 20*\[Nu] + 30*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(15 + 106*\[Nu] - 262*\[Nu]^2)*\[Chi]S) - 
         4*rd^4*((15 - 92*\[Nu] + 126*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(15 - 64*\[Nu] + 70*\[Nu]^2)*\[Chi]S) - 
         3*r^2*rd^2*\[Phi]d^2*((35 - 198*\[Nu] + 244*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(35 - 382*\[Nu] + 612*\[Nu]^2)*\[Chi]S)) + 
       M*r*(r^2*\[Phi]d^2*((609 - 3352*\[Nu] + 3966*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(609 + 58*\[Nu] - 2854*\[Nu]^2)*\[Chi]S) + 
         (2*I)*r*rd*\[Phi]d*(5*(78 - 499*\[Nu] + 717*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(390 + 677*\[Nu] - 2759*\[Nu]^2)*\[Chi]S) - 
         2*rd^2*((255 - 1592*\[Nu] + 2226*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(255 - 1144*\[Nu] + 1330*\[Nu]^2)*\[Chi]S))))/
     (Sqrt[385]*r^4)
 
Hinst[5, 2] = (M*\[Epsilon]^4*(1 - 5*\[Nu] + 5*\[Nu]^2)*\[Phi]d*
      (M*((22*I)*rd + 41*r*\[Phi]d) - 3*r*((2*I)*rd^3 + 6*r*rd^2*\[Phi]d - 
         (3*I)*r^2*rd*\[Phi]d^2 + 11*r^3*\[Phi]d^3)))/(108*Sqrt[55]*r) + 
    (M^2*SO*\[Epsilon]^5*\[Nu]*((-22*I)*M*rd + (6*I)*r*rd^3 - 
       41*M*r*\[Phi]d + 18*r^2*rd^2*\[Phi]d - (9*I)*r^3*rd*\[Phi]d^2 + 
       33*r^4*\[Phi]d^3)*(\[Delta]*\[Chi]A + (-1 + 2*\[Nu])*\[Chi]S))/
     (36*Sqrt[55]*r^3)
 
Hinst[5, 3] = (\[Delta]*\[Epsilon]^3*(-1 + 2*\[Nu])*
      (M^2*(-410*rd + (516*I)*r*\[Phi]d) - (120*I)*r^2*(rd - I*r*\[Phi]d)^4*
        ((-I)*rd + r*\[Phi]d) + 3*M*r*(160*rd^3 - (396*I)*r*rd^2*\[Phi]d - 
         240*r^2*rd*\[Phi]d^2 - (51*I)*r^3*\[Phi]d^3)))/(288*Sqrt[330]*r^2) + 
    ((I/144)*M^2*SO*\[Epsilon]^6*
      (M^2*((90 - 622*\[Nu] + 966*\[Nu]^2)*\[Chi]A + 
         2*\[Delta]*(45 - 262*\[Nu] + 385*\[Nu]^2)*\[Chi]S) + 
       M*r*(r^2*\[Phi]d^2*((849 - 5360*\[Nu] + 7590*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(849 - 1774*\[Nu] + 418*\[Nu]^2)*\[Chi]S) + 
         (6*I)*r*rd*\[Phi]d*(5*(38 - 243*\[Nu] + 349*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(190 - 643*\[Nu] + 601*\[Nu]^2)*\[Chi]S) - 
         2*rd^2*((255 - 1592*\[Nu] + 2226*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(255 - 1144*\[Nu] + 1330*\[Nu]^2)*\[Chi]S)) - 
       3*r^2*(3*r^2*rd^2*\[Phi]d^2*((45 - 298*\[Nu] + 444*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(45 - 2*\[Nu] - 148*\[Nu]^2)*\[Chi]S) - 
         4*rd^4*((15 - 92*\[Nu] + 126*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(15 - 64*\[Nu] + 70*\[Nu]^2)*\[Chi]S) + 
         r^4*\[Phi]d^4*((20 - 317*\[Nu] + 751*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(20 + 13*\[Nu] + 91*\[Nu]^2)*\[Chi]S) + 
         (4*I)*r*rd^3*\[Phi]d*(5*(9 - 56*\[Nu] + 78*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(45 - 142*\[Nu] + 114*\[Nu]^2)*\[Chi]S) + 
         I*r^3*rd*\[Phi]d^3*(25*(3 - 14*\[Nu] + 12*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(75 - 958*\[Nu] + 1516*\[Nu]^2)*\[Chi]S))))/
     (Sqrt[330]*r^4)
 
Hinst[5, 4] = -(M*\[Epsilon]^4*(1 - 5*\[Nu] + 5*\[Nu]^2)*\[Phi]d*
       ((22*I)*M*rd - (6*I)*r*rd^3 + 82*M*r*\[Phi]d - 36*r^2*rd^2*\[Phi]d + 
        (99*I)*r^3*rd*\[Phi]d^2 + 174*r^4*\[Phi]d^3))/(72*Sqrt[165]*r) + 
    (M^2*SO*\[Epsilon]^5*\[Nu]*((22*I)*M*rd - (6*I)*r*rd^3 + 82*M*r*\[Phi]d - 
       36*r^2*rd^2*\[Phi]d + (99*I)*r^3*rd*\[Phi]d^2 + 174*r^4*\[Phi]d^3)*
      (\[Delta]*\[Chi]A + (-1 + 2*\[Nu])*\[Chi]S))/(24*Sqrt[165]*r^3)
 
Hinst[5, 5] = (\[Delta]*\[Epsilon]^3*(-1 + 2*\[Nu])*
      (24*r^2*(rd - I*r*\[Phi]d)^5 + 2*M^2*(41*rd - (86*I)*r*\[Phi]d) + 
       3*M*r*(-32*rd^3 + (132*I)*r*rd^2*\[Phi]d + 208*r^2*rd*\[Phi]d^2 - 
         (143*I)*r^3*\[Phi]d^3)))/(96*Sqrt[66]*r^2) - 
    ((I/144)*M^2*SO*\[Epsilon]^6*
      (M^2*((90 - 622*\[Nu] + 966*\[Nu]^2)*\[Chi]A + 
         2*\[Delta]*(45 - 262*\[Nu] + 385*\[Nu]^2)*\[Chi]S) + 
       3*r^2*(4*rd^4*((15 - 92*\[Nu] + 126*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(15 - 64*\[Nu] + 70*\[Nu]^2)*\[Chi]S) - 
         15*r^2*rd^2*\[Phi]d^2*((41 - 258*\[Nu] + 364*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(41 - 154*\[Nu] + 156*\[Nu]^2)*\[Chi]S) - 
         (4*I)*r*rd^3*\[Phi]d*((75 - 464*\[Nu] + 642*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(75 - 298*\[Nu] + 310*\[Nu]^2)*\[Chi]S) + 
         r^4*\[Phi]d^4*((300 - 2347*\[Nu] + 4041*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(300 - 869*\[Nu] + 1085*\[Nu]^2)*\[Chi]S) + 
         I*r^3*rd*\[Phi]d^3*((675 - 4414*\[Nu] + 6492*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(675 - 2558*\[Nu] + 2780*\[Nu]^2)*\[Chi]S)) + 
       M*r*(-2*rd^2*((255 - 1592*\[Nu] + 2226*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(255 - 1144*\[Nu] + 1330*\[Nu]^2)*\[Chi]S) + 
         (2*I)*r*rd*\[Phi]d*((870 - 5563*\[Nu] + 7989*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(870 - 3743*\[Nu] + 4349*\[Nu]^2)*\[Chi]S) + 
         r^2*\[Phi]d^2*((1329 - 9376*\[Nu] + 14838*\[Nu]^2)*\[Chi]A + 
           \[Delta]*(1329 - 5438*\[Nu] + 6962*\[Nu]^2)*\[Chi]S))))/
     (Sqrt[66]*r^4)
 
Hinst[6, 0] = (\[Epsilon]^4*(1 - 5*\[Nu] + 5*\[Nu]^2)*
     (-172*M^3 + 120*r^3*(rd^2 + r^2*\[Phi]d^2)^3 + 
      M^2*(806*r*rd^2 + 463*r^3*\[Phi]d^2) - 
      3*M*(200*r^2*rd^4 + 372*r^4*rd^2*\[Phi]d^2 + 137*r^6*\[Phi]d^4)))/
    (792*Sqrt[273]*r^3)
 
Hinst[6, 1] = ((I/2376)*M^2*SO*\[Epsilon]^6*\[Nu]*
     (410*M^2 - 2*M*r*(310*rd^2 - (343*I)*r*rd*\[Phi]d + 359*r^2*\[Phi]d^2) + 
      3*r^2*(40*rd^4 - (56*I)*r*rd^3*\[Phi]d + 108*r^2*rd^2*\[Phi]d^2 - 
        (196*I)*r^3*rd*\[Phi]d^3 + 103*r^4*\[Phi]d^4))*
     ((-1 + 3*\[Nu])*\[Chi]A - \[Delta]*(-1 + \[Nu])*\[Chi]S))/(Sqrt[26]*r^4)
 
Hinst[6, 2] = (\[Epsilon]^4*(1 - 5*\[Nu] + 5*\[Nu]^2)*
     (516*M^3 + 360*r^3*(rd - I*r*\[Phi]d)^4*((-I)*rd + r*\[Phi]d)^2 + 
      M^2*r*(-2418*rd^2 + (2920*I)*r*rd*\[Phi]d - 145*r^2*\[Phi]d^2) - 
      3*M*r^2*(-600*rd^4 + (1040*I)*r*rd^3*\[Phi]d - 252*r^2*rd^2*\[Phi]d^2 + 
        (1050*I)*r^3*rd*\[Phi]d^3 + 233*r^4*\[Phi]d^4)))/(4752*Sqrt[65]*r^3)
 
Hinst[6, 3] = ((-I/1584)*M^2*SO*\[Epsilon]^6*\[Nu]*
     (410*M^2 + 2*M*r*(-310*rd^2 + (1029*I)*r*rd*\[Phi]d + 
        929*r^2*\[Phi]d^2) - 3*r^2*(-40*rd^4 + (168*I)*r*rd^3*\[Phi]d + 
        228*r^2*rd^2*\[Phi]d^2 + (28*I)*r^3*rd*\[Phi]d^3 + 
        513*r^4*\[Phi]d^4))*((-1 + 3*\[Nu])*\[Chi]A - 
      \[Delta]*(-1 + \[Nu])*\[Chi]S))/(Sqrt[65]*r^4)
 
Hinst[6, 4] = (\[Epsilon]^4*(1 - 5*\[Nu] + 5*\[Nu]^2)*
     (-516*M^3 + 360*r^3*(rd - I*r*\[Phi]d)^5*(rd + I*r*\[Phi]d) + 
      M^2*r*(2418*rd^2 - (5840*I)*r*rd*\[Phi]d - 3587*r^2*\[Phi]d^2) + 
      15*M*r^2*(-120*rd^4 + (416*I)*r*rd^3*\[Phi]d + 468*r^2*rd^2*\[Phi]d^2 - 
        (108*I)*r^3*rd*\[Phi]d^3 + 113*r^4*\[Phi]d^4)))/(3960*Sqrt[78]*r^3)
 
Hinst[6, 5] = ((I/144)*M^2*SO*\[Epsilon]^6*\[Nu]*
     (82*M^2 + 2*M*r*(-62*rd^2 + (343*I)*r*rd*\[Phi]d + 701*r^2*\[Phi]d^2) + 
      3*r^2*(8*rd^4 - (56*I)*r*rd^3*\[Phi]d - 180*r^2*rd^2*\[Phi]d^2 + 
        (364*I)*r^3*rd*\[Phi]d^3 + 547*r^4*\[Phi]d^4))*
     ((-1 + 3*\[Nu])*\[Chi]A - \[Delta]*(-1 + \[Nu])*\[Chi]S))/(Sqrt[429]*r^4)
 
Hinst[6, 6] = (\[Epsilon]^4*(1 - 5*\[Nu] + 5*\[Nu]^2)*
     (172*M^3 + 120*r^3*(I*rd + r*\[Phi]d)^6 + 
      M^2*r*(-806*rd^2 + (2920*I)*r*rd*\[Phi]d + 3269*r^2*\[Phi]d^2) + 
      15*M*r^2*(40*rd^4 - (208*I)*r*rd^3*\[Phi]d - 444*r^2*rd^2*\[Phi]d^2 + 
        (494*I)*r^3*rd*\[Phi]d^3 + 281*r^4*\[Phi]d^4)))/(720*Sqrt[143]*r^3)


(* ::Section::Closed:: *)
(*Tail contributions*)


(* ::Text:: *)
(*The tail part of the waveform modes is expanded in eccentricity to O(e^6), and expressed in terms of (x,e,\[ScriptL]).*)


Htail[2, 0] = x^(5/2)*\[Epsilon]^3*
     (e*(E^(I*\[ScriptL])*(((11*I)/12)/Sqrt[6] + Pi/(2*Sqrt[6]) - 
          (I*\[Gamma]E)/Sqrt[6] - (I*Log[2*b])/Sqrt[6] - 
          (I/2)*Sqrt[3/2]*Log[x]) + (((-11*I)/12)/Sqrt[6] + Pi/(2*Sqrt[6]) + 
          (I*\[Gamma]E)/Sqrt[6] + (I*Log[2*b])/Sqrt[6] + 
          (I/2)*Sqrt[3/2]*Log[x])/E^(I*\[ScriptL])) + 
      e^2*(E^((2*I)*\[ScriptL])*(((11*I)/6)/Sqrt[6] + Pi/Sqrt[6] - 
          I*Sqrt[2/3]*\[Gamma]E - (2*I)*Sqrt[2/3]*Log[2] - 
          I*Sqrt[2/3]*Log[b] - I*Sqrt[3/2]*Log[x]) + 
        (((-11*I)/6)/Sqrt[6] + Pi/Sqrt[6] + I*Sqrt[2/3]*\[Gamma]E + 
          (2*I)*Sqrt[2/3]*Log[2] + I*Sqrt[2/3]*Log[b] + I*Sqrt[3/2]*Log[x])/
         E^((2*I)*\[ScriptL])) + 
      e^3*((((11*I)/96)/Sqrt[6] - Pi/(16*Sqrt[6]) - ((I/8)*\[Gamma]E)/
           Sqrt[6] - ((I/8)*Log[2*b])/Sqrt[6] - (I/16)*Sqrt[3/2]*Log[x])/
         E^(I*\[ScriptL]) + E^(I*\[ScriptL])*(((-11*I)/96)/Sqrt[6] - 
          Pi/(16*Sqrt[6]) + ((I/8)*\[Gamma]E)/Sqrt[6] + ((I/8)*Log[2*b])/
           Sqrt[6] + (I/16)*Sqrt[3/2]*Log[x]) + E^((3*I)*\[ScriptL])*
         (((33*I)/32)*Sqrt[3/2] + (9*Sqrt[3/2]*Pi)/16 - ((9*I)/8)*Sqrt[3/2]*
           \[Gamma]E - ((9*I)/8)*Sqrt[3/2]*Log[6*b] - ((27*I)/16)*Sqrt[3/2]*
           Log[x]) + (((-33*I)/32)*Sqrt[3/2] + (9*Sqrt[3/2]*Pi)/16 + 
          ((9*I)/8)*Sqrt[3/2]*\[Gamma]E + ((9*I)/8)*Sqrt[3/2]*Log[6*b] + 
          ((27*I)/16)*Sqrt[3/2]*Log[x])/E^((3*I)*\[ScriptL])) + 
      e^6*((((88*I)/45)*Sqrt[2/3] - (16*Sqrt[2/3]*Pi)/15 - 
          ((32*I)/15)*Sqrt[2/3]*\[Gamma]E - ((32*I)/5)*Sqrt[2/3]*Log[2] - 
          ((32*I)/15)*Sqrt[2/3]*Log[b] - ((16*I)/5)*Sqrt[2/3]*Log[x])/
         E^((4*I)*\[ScriptL]) + E^((4*I)*\[ScriptL])*
         (((-88*I)/45)*Sqrt[2/3] - (16*Sqrt[2/3]*Pi)/15 + 
          ((32*I)/15)*Sqrt[2/3]*\[Gamma]E + ((32*I)/5)*Sqrt[2/3]*Log[2] + 
          ((32*I)/15)*Sqrt[2/3]*Log[b] + ((16*I)/5)*Sqrt[2/3]*Log[x]) + 
        E^((6*I)*\[ScriptL])*(((297*I)/80)*Sqrt[3/2] + (81*Sqrt[3/2]*Pi)/40 - 
          ((81*I)/20)*Sqrt[3/2]*\[Gamma]E - ((81*I)/10)*Sqrt[3/2]*Log[2] - 
          ((81*I)/20)*Sqrt[3/2]*Log[3*b] - ((243*I)/40)*Sqrt[3/2]*Log[x]) + 
        (((-297*I)/80)*Sqrt[3/2] + (81*Sqrt[3/2]*Pi)/40 + 
          ((81*I)/20)*Sqrt[3/2]*\[Gamma]E + ((81*I)/10)*Sqrt[3/2]*Log[2] + 
          ((81*I)/20)*Sqrt[3/2]*Log[3*b] + ((243*I)/40)*Sqrt[3/2]*Log[x])/
         E^((6*I)*\[ScriptL]) + E^((2*I)*\[ScriptL])*(((11*I)/144)/Sqrt[6] + 
          Pi/(24*Sqrt[6]) - ((I/12)*\[Gamma]E)/Sqrt[6] - 
          ((I/6)*Log[2])/Sqrt[6] - ((I/12)*Log[b])/Sqrt[6] - 
          ((I/8)*Log[x])/Sqrt[6]) + (((-11*I)/144)/Sqrt[6] + 
          Pi/(24*Sqrt[6]) + ((I/12)*\[Gamma]E)/Sqrt[6] + 
          ((I/6)*Log[2])/Sqrt[6] + ((I/12)*Log[b])/Sqrt[6] + 
          ((I/8)*Log[x])/Sqrt[6])/E^((2*I)*\[ScriptL])) + 
      e^4*(E^((4*I)*\[ScriptL])*(((22*I)/9)*Sqrt[2/3] + (4*Sqrt[2/3]*Pi)/3 - 
          ((8*I)/3)*Sqrt[2/3]*\[Gamma]E - (8*I)*Sqrt[2/3]*Log[2] - 
          ((8*I)/3)*Sqrt[2/3]*Log[b] - (4*I)*Sqrt[2/3]*Log[x]) + 
        (((-22*I)/9)*Sqrt[2/3] + (4*Sqrt[2/3]*Pi)/3 + ((8*I)/3)*Sqrt[2/3]*
           \[Gamma]E + (8*I)*Sqrt[2/3]*Log[2] + ((8*I)/3)*Sqrt[2/3]*Log[b] + 
          (4*I)*Sqrt[2/3]*Log[x])/E^((4*I)*\[ScriptL]) + 
        (((11*I)/18)/Sqrt[6] - Pi/(3*Sqrt[6]) - (I/3)*Sqrt[2/3]*\[Gamma]E - 
          ((2*I)/3)*Sqrt[2/3]*Log[2] - (I/3)*Sqrt[2/3]*Log[b] - 
          (I*Log[x])/Sqrt[6])/E^((2*I)*\[ScriptL]) + E^((2*I)*\[ScriptL])*
         (((-11*I)/18)/Sqrt[6] - Pi/(3*Sqrt[6]) + (I/3)*Sqrt[2/3]*\[Gamma]E + 
          ((2*I)/3)*Sqrt[2/3]*Log[2] + (I/3)*Sqrt[2/3]*Log[b] + 
          (I*Log[x])/Sqrt[6])) + 
      e^5*((((297*I)/512)*Sqrt[3/2] - (81*Sqrt[3/2]*Pi)/256 - 
          ((81*I)/128)*Sqrt[3/2]*\[Gamma]E - ((81*I)/128)*Sqrt[3/2]*
           Log[6*b] - ((243*I)/256)*Sqrt[3/2]*Log[x])/E^((3*I)*\[ScriptL]) + 
        E^((3*I)*\[ScriptL])*(((-297*I)/512)*Sqrt[3/2] - (81*Sqrt[3/2]*Pi)/
           256 + ((81*I)/128)*Sqrt[3/2]*\[Gamma]E + ((81*I)/128)*Sqrt[3/2]*
           Log[6*b] + ((243*I)/256)*Sqrt[3/2]*Log[x]) + 
        E^(I*\[ScriptL])*(((11*I)/2304)/Sqrt[6] + Pi/(384*Sqrt[6]) - 
          ((I/192)*\[Gamma]E)/Sqrt[6] - ((I/192)*Log[2*b])/Sqrt[6] - 
          ((I/128)*Log[x])/Sqrt[6]) + (((-11*I)/2304)/Sqrt[6] + 
          Pi/(384*Sqrt[6]) + ((I/192)*\[Gamma]E)/Sqrt[6] + 
          ((I/192)*Log[2*b])/Sqrt[6] + ((I/128)*Log[x])/Sqrt[6])/
         E^(I*\[ScriptL]) + E^((5*I)*\[ScriptL])*(((34375*I)/4608)/Sqrt[6] + 
          (3125*Pi)/(768*Sqrt[6]) - (((3125*I)/384)*\[Gamma]E)/Sqrt[6] - 
          (((3125*I)/384)*Log[10*b])/Sqrt[6] - (((3125*I)/256)*Log[x])/
           Sqrt[6]) + (((-34375*I)/4608)/Sqrt[6] + (3125*Pi)/(768*Sqrt[6]) + 
          (((3125*I)/384)*\[Gamma]E)/Sqrt[6] + (((3125*I)/384)*Log[10*b])/
           Sqrt[6] + (((3125*I)/256)*Log[x])/Sqrt[6])/
         E^((5*I)*\[ScriptL]))) + SO*x^4*\[Epsilon]^6*
     (\[Delta]*\[Chi]A*(e^2*(E^((2*I)*\[ScriptL])*(((52*I)/9)*Sqrt[2/3] + 
            (16*Sqrt[2/3]*Pi)/3 - ((32*I)/3)*Sqrt[2/3]*\[Gamma]E - 
            ((64*I)/3)*Sqrt[2/3]*Log[2] - ((32*I)/3)*Sqrt[2/3]*Log[b] - 
            (16*I)*Sqrt[2/3]*Log[x]) + (((-52*I)/9)*Sqrt[2/3] + 
            (16*Sqrt[2/3]*Pi)/3 + ((32*I)/3)*Sqrt[2/3]*\[Gamma]E + 
            ((64*I)/3)*Sqrt[2/3]*Log[2] + ((32*I)/3)*Sqrt[2/3]*Log[b] + 
            (16*I)*Sqrt[2/3]*Log[x])/E^((2*I)*\[ScriptL])) + 
        e^4*(E^((2*I)*\[ScriptL])*(((50*I)/27)*Sqrt[2/3] + (32*Sqrt[2/3]*Pi)/
             9 - ((64*I)/9)*Sqrt[2/3]*\[Gamma]E - ((128*I)/9)*Sqrt[2/3]*
             Log[2] - ((64*I)/9)*Sqrt[2/3]*Log[b] - ((32*I)/3)*Sqrt[2/3]*
             Log[x]) + (((-50*I)/27)*Sqrt[2/3] + (32*Sqrt[2/3]*Pi)/9 + 
            ((64*I)/9)*Sqrt[2/3]*\[Gamma]E + ((128*I)/9)*Sqrt[2/3]*Log[2] + 
            ((64*I)/9)*Sqrt[2/3]*Log[b] + ((32*I)/3)*Sqrt[2/3]*Log[x])/
           E^((2*I)*\[ScriptL]) + E^((4*I)*\[ScriptL])*
           (((680*I)/27)*Sqrt[2/3] + (176*Sqrt[2/3]*Pi)/9 - 
            ((352*I)/9)*Sqrt[2/3]*\[Gamma]E - ((352*I)/3)*Sqrt[2/3]*Log[2] - 
            ((352*I)/9)*Sqrt[2/3]*Log[b] - ((176*I)/3)*Sqrt[2/3]*Log[x]) + 
          (((-680*I)/27)*Sqrt[2/3] + (176*Sqrt[2/3]*Pi)/9 + 
            ((352*I)/9)*Sqrt[2/3]*\[Gamma]E + ((352*I)/3)*Sqrt[2/3]*Log[2] + 
            ((352*I)/9)*Sqrt[2/3]*Log[b] + ((176*I)/3)*Sqrt[2/3]*Log[x])/
           E^((4*I)*\[ScriptL])) + e*(E^(I*\[ScriptL])*(((71*I)/18)/Sqrt[6] + 
            (13*Pi)/(3*Sqrt[6]) - ((13*I)/3)*Sqrt[2/3]*\[Gamma]E - 
            ((13*I)/3)*Sqrt[2/3]*Log[2*b] - ((13*I)*Log[x])/Sqrt[6]) + 
          (((-71*I)/18)/Sqrt[6] + (13*Pi)/(3*Sqrt[6]) + ((13*I)/3)*Sqrt[2/3]*
             \[Gamma]E + ((13*I)/3)*Sqrt[2/3]*Log[2*b] + ((13*I)*Log[x])/
             Sqrt[6])/E^(I*\[ScriptL])) + 
        e^3*(E^((3*I)*\[ScriptL])*(((137*I)/16)*Sqrt[3/2] + 
            (57*Sqrt[3/2]*Pi)/8 - ((57*I)/4)*Sqrt[3/2]*\[Gamma]E - 
            ((57*I)/4)*Sqrt[3/2]*Log[6*b] - ((171*I)/8)*Sqrt[3/2]*Log[x]) + 
          (((-137*I)/16)*Sqrt[3/2] + (57*Sqrt[3/2]*Pi)/8 + 
            ((57*I)/4)*Sqrt[3/2]*\[Gamma]E + ((57*I)/4)*Sqrt[3/2]*Log[6*b] + 
            ((171*I)/8)*Sqrt[3/2]*Log[x])/E^((3*I)*\[ScriptL]) + 
          E^(I*\[ScriptL])*(((451*I)/144)/Sqrt[6] + (113*Pi)/(24*Sqrt[6]) - 
            (((113*I)/12)*\[Gamma]E)/Sqrt[6] - (((113*I)/12)*Log[2*b])/
             Sqrt[6] - (((113*I)/8)*Log[x])/Sqrt[6]) + 
          (((-451*I)/144)/Sqrt[6] + (113*Pi)/(24*Sqrt[6]) + 
            (((113*I)/12)*\[Gamma]E)/Sqrt[6] + (((113*I)/12)*Log[2*b])/
             Sqrt[6] + (((113*I)/8)*Log[x])/Sqrt[6])/E^(I*\[ScriptL])) + 
        e^6*((((1316*I)/135)*Sqrt[2/3] - (56*Sqrt[2/3]*Pi)/45 - 
            ((112*I)/45)*Sqrt[2/3]*\[Gamma]E - ((112*I)/15)*Sqrt[2/3]*
             Log[2] - ((112*I)/45)*Sqrt[2/3]*Log[b] - ((56*I)/15)*Sqrt[2/3]*
             Log[x])/E^((4*I)*\[ScriptL]) + E^((4*I)*\[ScriptL])*
           (((-1316*I)/135)*Sqrt[2/3] - (56*Sqrt[2/3]*Pi)/45 + 
            ((112*I)/45)*Sqrt[2/3]*\[Gamma]E + ((112*I)/15)*Sqrt[2/3]*
             Log[2] + ((112*I)/45)*Sqrt[2/3]*Log[b] + ((56*I)/15)*Sqrt[2/3]*
             Log[x]) + E^((6*I)*\[ScriptL])*(((531*I)/10)*Sqrt[3/2] + 
            (189*Sqrt[3/2]*Pi)/5 - ((189*I)/5)*Sqrt[6]*\[Gamma]E - 
            ((378*I)/5)*Sqrt[6]*Log[2] - ((189*I)/5)*Sqrt[6]*Log[3*b] - 
            ((567*I)/5)*Sqrt[3/2]*Log[x]) + (((-531*I)/10)*Sqrt[3/2] + 
            (189*Sqrt[3/2]*Pi)/5 + ((189*I)/5)*Sqrt[6]*\[Gamma]E + 
            ((378*I)/5)*Sqrt[6]*Log[2] + ((189*I)/5)*Sqrt[6]*Log[3*b] + 
            ((567*I)/5)*Sqrt[3/2]*Log[x])/E^((6*I)*\[ScriptL]) + 
          E^((2*I)*\[ScriptL])*(((389*I)/54)/Sqrt[6] + (91*Pi)/(9*Sqrt[6]) - 
            ((91*I)/9)*Sqrt[2/3]*\[Gamma]E - ((182*I)/9)*Sqrt[2/3]*Log[2] - 
            ((91*I)/9)*Sqrt[2/3]*Log[b] - (((91*I)/3)*Log[x])/Sqrt[6]) + 
          (((-389*I)/54)/Sqrt[6] + (91*Pi)/(9*Sqrt[6]) + ((91*I)/9)*Sqrt[2/3]*
             \[Gamma]E + ((182*I)/9)*Sqrt[2/3]*Log[2] + ((91*I)/9)*Sqrt[2/3]*
             Log[b] + (((91*I)/3)*Log[x])/Sqrt[6])/E^((2*I)*\[ScriptL])) + 
        e^5*(E^((3*I)*\[ScriptL])*(((-123*I)/256)*Sqrt[3/2] + 
            (261*Sqrt[3/2]*Pi)/128 - ((261*I)/64)*Sqrt[3/2]*\[Gamma]E - 
            ((261*I)/64)*Sqrt[3/2]*Log[6*b] - ((783*I)/128)*Sqrt[3/2]*
             Log[x]) + (((123*I)/256)*Sqrt[3/2] + (261*Sqrt[3/2]*Pi)/128 + 
            ((261*I)/64)*Sqrt[3/2]*\[Gamma]E + ((261*I)/64)*Sqrt[3/2]*
             Log[6*b] + ((783*I)/128)*Sqrt[3/2]*Log[x])/
           E^((3*I)*\[ScriptL]) + E^(I*\[ScriptL])*
           (((13703*I)/3456)/Sqrt[6] + (3373*Pi)/(576*Sqrt[6]) - 
            (((3373*I)/288)*\[Gamma]E)/Sqrt[6] - (((3373*I)/288)*Log[2*b])/
             Sqrt[6] - (((3373*I)/192)*Log[x])/Sqrt[6]) + 
          (((-13703*I)/3456)/Sqrt[6] + (3373*Pi)/(576*Sqrt[6]) + 
            (((3373*I)/288)*\[Gamma]E)/Sqrt[6] + (((3373*I)/288)*Log[2*b])/
             Sqrt[6] + (((3373*I)/192)*Log[x])/Sqrt[6])/E^(I*\[ScriptL]) + 
          E^((5*I)*\[ScriptL])*(((634375*I)/6912)/Sqrt[6] + 
            (78125*Pi)/(1152*Sqrt[6]) - (((78125*I)/576)*\[Gamma]E)/Sqrt[6] - 
            (((78125*I)/576)*Log[10*b])/Sqrt[6] - (((78125*I)/384)*Log[x])/
             Sqrt[6]) + (((-634375*I)/6912)/Sqrt[6] + (78125*Pi)/
             (1152*Sqrt[6]) + (((78125*I)/576)*\[Gamma]E)/Sqrt[6] + 
            (((78125*I)/576)*Log[10*b])/Sqrt[6] + (((78125*I)/384)*Log[x])/
             Sqrt[6])/E^((5*I)*\[ScriptL]))) + 
      \[Chi]S*(e*(((Pi*(13 - 8*\[Nu]))/(3*Sqrt[6]) - (I/3)*Sqrt[2/3]*
             \[Gamma]E*(-13 + 8*\[Nu]) + ((I/18)*(-71 + 52*\[Nu]))/Sqrt[6] - 
            (I/3)*Sqrt[2/3]*(-13 + 8*\[Nu])*Log[2*b] + 
            (I*(13 - 8*\[Nu])*Log[x])/Sqrt[6])/E^(I*\[ScriptL]) + 
          E^(I*\[ScriptL])*((Pi*(13 - 8*\[Nu]))/(3*Sqrt[6]) + 
            (I/3)*Sqrt[2/3]*\[Gamma]E*(-13 + 8*\[Nu]) - 
            ((I/18)*(-71 + 52*\[Nu]))/Sqrt[6] + (I/3)*Sqrt[2/3]*
             (-13 + 8*\[Nu])*Log[2*b] + (I*(-13 + 8*\[Nu])*Log[x])/
             Sqrt[6])) + e^3*(((-3*Sqrt[3/2]*Pi*(-19 + 14*\[Nu]))/8 - 
            ((3*I)/4)*Sqrt[3/2]*\[Gamma]E*(-19 + 14*\[Nu]) + 
            (I/16)*Sqrt[3/2]*(-137 + 118*\[Nu]) - ((3*I)/4)*Sqrt[3/2]*
             (-19 + 14*\[Nu])*Log[6*b] - ((9*I)/8)*Sqrt[3/2]*(-19 + 14*\[Nu])*
             Log[x])/E^((3*I)*\[ScriptL]) + E^((3*I)*\[ScriptL])*
           ((-3*Sqrt[3/2]*Pi*(-19 + 14*\[Nu]))/8 + ((3*I)/4)*Sqrt[3/2]*
             \[Gamma]E*(-19 + 14*\[Nu]) - (I/16)*Sqrt[3/2]*
             (-137 + 118*\[Nu]) + ((3*I)/4)*Sqrt[3/2]*(-19 + 14*\[Nu])*
             Log[6*b] + ((9*I)/8)*Sqrt[3/2]*(-19 + 14*\[Nu])*Log[x]) + 
          ((Pi*(113 - 58*\[Nu]))/(24*Sqrt[6]) + 
            (((11*I)/144)*(-41 + 22*\[Nu]))/Sqrt[6] - 
            ((I/12)*\[Gamma]E*(-113 + 58*\[Nu]))/Sqrt[6] - 
            ((I/12)*(-113 + 58*\[Nu])*Log[2*b])/Sqrt[6] - 
            ((I/8)*(-113 + 58*\[Nu])*Log[x])/Sqrt[6])/E^(I*\[ScriptL]) + 
          E^(I*\[ScriptL])*((Pi*(113 - 58*\[Nu]))/(24*Sqrt[6]) - 
            (((11*I)/144)*(-41 + 22*\[Nu]))/Sqrt[6] + 
            ((I/12)*\[Gamma]E*(-113 + 58*\[Nu]))/Sqrt[6] + 
            ((I/12)*(-113 + 58*\[Nu])*Log[2*b])/Sqrt[6] + 
            ((I/8)*(-113 + 58*\[Nu])*Log[x])/Sqrt[6])) + 
        e^5*(((-9*Sqrt[3/2]*Pi*(-29 + 4*\[Nu]))/128 - ((9*I)/64)*Sqrt[3/2]*
             \[Gamma]E*(-29 + 4*\[Nu]) - ((3*I)/256)*Sqrt[3/2]*
             (-41 + 136*\[Nu]) - ((9*I)/64)*Sqrt[3/2]*(-29 + 4*\[Nu])*
             Log[6*b] - ((27*I)/128)*Sqrt[3/2]*(-29 + 4*\[Nu])*Log[x])/
           E^((3*I)*\[ScriptL]) + E^((3*I)*\[ScriptL])*
           ((-9*Sqrt[3/2]*Pi*(-29 + 4*\[Nu]))/128 + ((9*I)/64)*Sqrt[3/2]*
             \[Gamma]E*(-29 + 4*\[Nu]) + ((3*I)/256)*Sqrt[3/2]*
             (-41 + 136*\[Nu]) + ((9*I)/64)*Sqrt[3/2]*(-29 + 4*\[Nu])*
             Log[6*b] + ((27*I)/128)*Sqrt[3/2]*(-29 + 4*\[Nu])*Log[x]) + 
          ((-15625*Pi*(-5 + 4*\[Nu]))/(1152*Sqrt[6]) - 
            (((15625*I)/576)*\[Gamma]E*(-5 + 4*\[Nu]))/Sqrt[6] + 
            (((3125*I)/6912)*(-203 + 184*\[Nu]))/Sqrt[6] - 
            (((15625*I)/576)*(-5 + 4*\[Nu])*Log[10*b])/Sqrt[6] - 
            (((15625*I)/384)*(-5 + 4*\[Nu])*Log[x])/Sqrt[6])/
           E^((5*I)*\[ScriptL]) + E^((5*I)*\[ScriptL])*
           ((-15625*Pi*(-5 + 4*\[Nu]))/(1152*Sqrt[6]) + 
            (((15625*I)/576)*\[Gamma]E*(-5 + 4*\[Nu]))/Sqrt[6] - 
            (((3125*I)/6912)*(-203 + 184*\[Nu]))/Sqrt[6] + 
            (((15625*I)/576)*(-5 + 4*\[Nu])*Log[10*b])/Sqrt[6] + 
            (((15625*I)/384)*(-5 + 4*\[Nu])*Log[x])/Sqrt[6]) + 
          ((Pi*(3373 - 1748*\[Nu]))/(576*Sqrt[6]) - 
            ((I/288)*\[Gamma]E*(-3373 + 1748*\[Nu]))/Sqrt[6] + 
            ((I/3456)*(-13703 + 7528*\[Nu]))/Sqrt[6] - 
            ((I/288)*(-3373 + 1748*\[Nu])*Log[2*b])/Sqrt[6] - 
            ((I/192)*(-3373 + 1748*\[Nu])*Log[x])/Sqrt[6])/E^(I*\[ScriptL]) + 
          E^(I*\[ScriptL])*((Pi*(3373 - 1748*\[Nu]))/(576*Sqrt[6]) + 
            ((I/288)*\[Gamma]E*(-3373 + 1748*\[Nu]))/Sqrt[6] - 
            ((I/3456)*(-13703 + 7528*\[Nu]))/Sqrt[6] + 
            ((I/288)*(-3373 + 1748*\[Nu])*Log[2*b])/Sqrt[6] + 
            ((I/192)*(-3373 + 1748*\[Nu])*Log[x])/Sqrt[6])) + 
        e^2*(E^((2*I)*\[ScriptL])*((Sqrt[2/3]*Pi*(16 - 11*\[Nu]))/3 + 
            ((2*I)/3)*Sqrt[2/3]*\[Gamma]E*(-16 + 11*\[Nu]) - 
            ((I/9)*(-104 + 85*\[Nu]))/Sqrt[6] - ((64*I)/3)*Sqrt[2/3]*Log[2] + 
            ((22*I)/3)*Sqrt[2/3]*\[Nu]*Log[4] + ((2*I)/3)*Sqrt[2/3]*
             (-16 + 11*\[Nu])*Log[b] + I*Sqrt[2/3]*(-16 + 11*\[Nu])*Log[x]) + 
          ((Sqrt[2/3]*Pi*(16 - 11*\[Nu]))/3 - ((2*I)/3)*Sqrt[2/3]*\[Gamma]E*
             (-16 + 11*\[Nu]) + ((I/9)*(-104 + 85*\[Nu]))/Sqrt[6] + 
            ((64*I)/3)*Sqrt[2/3]*Log[2] + ((32*I)/3)*Sqrt[2/3]*Log[b] + 
            (16*I)*Sqrt[2/3]*Log[x] - ((11*I)/3)*Sqrt[2/3]*\[Nu]*
             Log[16*b^2*x^3])/E^((2*I)*\[ScriptL])) + 
        e^4*(((-8*Sqrt[2/3]*Pi*(-22 + 17*\[Nu]))/9 - ((16*I)/9)*Sqrt[2/3]*
             \[Gamma]E*(-22 + 17*\[Nu]) + ((4*I)/27)*Sqrt[2/3]*
             (-170 + 151*\[Nu]) - ((16*I)/3)*Sqrt[2/3]*(-22 + 17*\[Nu])*
             Log[2] - ((16*I)/9)*Sqrt[2/3]*(-22 + 17*\[Nu])*Log[b] - 
            ((8*I)/3)*Sqrt[2/3]*(-22 + 17*\[Nu])*Log[x])/
           E^((4*I)*\[ScriptL]) + E^((4*I)*\[ScriptL])*
           ((-8*Sqrt[2/3]*Pi*(-22 + 17*\[Nu]))/9 + ((16*I)/9)*Sqrt[2/3]*
             \[Gamma]E*(-22 + 17*\[Nu]) - ((4*I)/27)*Sqrt[2/3]*
             (-170 + 151*\[Nu]) + ((16*I)/3)*Sqrt[2/3]*(-22 + 17*\[Nu])*
             Log[2] + ((16*I)/9)*Sqrt[2/3]*(-22 + 17*\[Nu])*Log[b] + 
            ((8*I)/3)*Sqrt[2/3]*(-22 + 17*\[Nu])*Log[x]) + 
          E^((2*I)*\[ScriptL])*((Pi*(64 - 29*\[Nu]))/(9*Sqrt[6]) + 
            (I/9)*Sqrt[2/3]*\[Gamma]E*(-64 + 29*\[Nu]) - 
            ((I/54)*(-200 + 67*\[Nu]))/Sqrt[6] + ((2*I)/9)*Sqrt[2/3]*
             (-64 + 29*\[Nu])*Log[2] + (I/9)*Sqrt[2/3]*(-64 + 29*\[Nu])*
             Log[b] + ((I/3)*(-64 + 29*\[Nu])*Log[x])/Sqrt[6]) + 
          ((Pi*(64 - 29*\[Nu]))/(9*Sqrt[6]) - (I/9)*Sqrt[2/3]*\[Gamma]E*
             (-64 + 29*\[Nu]) + ((I/54)*(-200 + 67*\[Nu]))/Sqrt[6] + 
            ((128*I)/9)*Sqrt[2/3]*Log[2] + ((64*I)/9)*Sqrt[2/3]*Log[b] + 
            ((32*I)/3)*Sqrt[2/3]*Log[x] - (((29*I)/9)*\[Nu]*Log[16*b^2*x^3])/
             Sqrt[6])/E^((2*I)*\[ScriptL])) + 
        e^6*(E^((4*I)*\[ScriptL])*((28*Sqrt[2/3]*Pi*(-2 + 7*\[Nu]))/45 - 
            ((56*I)/45)*Sqrt[2/3]*\[Gamma]E*(-2 + 7*\[Nu]) + 
            ((14*I)/135)*Sqrt[2/3]*(-94 + 113*\[Nu]) - ((56*I)/15)*Sqrt[2/3]*
             (-2 + 7*\[Nu])*Log[2] - ((56*I)/45)*Sqrt[2/3]*(-2 + 7*\[Nu])*
             Log[b] - ((28*I)/15)*Sqrt[2/3]*(-2 + 7*\[Nu])*Log[x]) + 
          ((28*Sqrt[2/3]*Pi*(-2 + 7*\[Nu]))/45 + ((56*I)/45)*Sqrt[2/3]*
             \[Gamma]E*(-2 + 7*\[Nu]) - ((14*I)/135)*Sqrt[2/3]*
             (-94 + 113*\[Nu]) + ((56*I)/15)*Sqrt[2/3]*(-2 + 7*\[Nu])*
             Log[2] + ((56*I)/45)*Sqrt[2/3]*(-2 + 7*\[Nu])*Log[b] + 
            ((28*I)/15)*Sqrt[2/3]*(-2 + 7*\[Nu])*Log[x])/
           E^((4*I)*\[ScriptL]) + ((-27*Sqrt[3/2]*Pi*(-28 + 23*\[Nu]))/20 - 
            ((27*I)/10)*Sqrt[3/2]*\[Gamma]E*(-28 + 23*\[Nu]) + 
            ((9*I)/40)*Sqrt[3/2]*(-236 + 217*\[Nu]) - ((27*I)/5)*Sqrt[3/2]*
             (-28 + 23*\[Nu])*Log[2] - ((27*I)/10)*Sqrt[3/2]*(-28 + 23*\[Nu])*
             Log[3*b] - ((81*I)/20)*Sqrt[3/2]*(-28 + 23*\[Nu])*Log[x])/
           E^((6*I)*\[ScriptL]) + E^((6*I)*\[ScriptL])*
           ((-27*Sqrt[3/2]*Pi*(-28 + 23*\[Nu]))/20 + ((27*I)/10)*Sqrt[3/2]*
             \[Gamma]E*(-28 + 23*\[Nu]) - ((9*I)/40)*Sqrt[3/2]*
             (-236 + 217*\[Nu]) + ((27*I)/5)*Sqrt[3/2]*(-28 + 23*\[Nu])*
             Log[2] + ((27*I)/10)*Sqrt[3/2]*(-28 + 23*\[Nu])*Log[3*b] + 
            ((81*I)/20)*Sqrt[3/2]*(-28 + 23*\[Nu])*Log[x]) + 
          E^((2*I)*\[ScriptL])*((Pi*(182 - 97*\[Nu]))/(18*Sqrt[6]) + 
            ((I/9)*\[Gamma]E*(-182 + 97*\[Nu]))/Sqrt[6] - 
            ((I/108)*(-778 + 455*\[Nu]))/Sqrt[6] - ((182*I)/9)*Sqrt[2/3]*
             Log[2] + (((97*I)/9)*\[Nu]*Log[4])/Sqrt[6] + 
            ((I/9)*(-182 + 97*\[Nu])*Log[b])/Sqrt[6] + 
            ((I/6)*(-182 + 97*\[Nu])*Log[x])/Sqrt[6]) + 
          ((Pi*(182 - 97*\[Nu]))/(18*Sqrt[6]) - ((I/9)*\[Gamma]E*
              (-182 + 97*\[Nu]))/Sqrt[6] + ((I/108)*(-778 + 455*\[Nu]))/
             Sqrt[6] + ((182*I)/9)*Sqrt[2/3]*Log[2] + ((91*I)/9)*Sqrt[2/3]*
             Log[b] + (((91*I)/3)*Log[x])/Sqrt[6] - 
            (((97*I)/18)*\[Nu]*Log[16*b^2*x^3])/Sqrt[6])/
           E^((2*I)*\[ScriptL]))))
 
Htail[2, 1] = x^3*\[Delta]*\[Epsilon]^4*(7/9 + (I/3)*Pi - (2*\[Gamma]E)/3 + 
      Log[1/(64*b^6*x^9)]/9 + e*(E^(I*\[ScriptL])*(7/9 + (I/3)*Pi - 
          (2*\[Gamma]E)/3 + Log[1/(64*b^6*x^9)]/9) + 
        (7/3 + I*Pi - 2*\[Gamma]E - (14*Log[2])/3 - 2*Log[b] - 3*Log[x])/
         E^(I*\[ScriptL])) + 
      e^3*((595/72 + ((85*I)/24)*Pi - (85*\[Gamma]E)/12 - (1267*Log[2])/36 + 
          (27*Log[3])/4 - (85*Log[b])/12 - (85*Log[x])/8)/
         E^((3*I)*\[ScriptL]) + E^(I*\[ScriptL])*(133/72 + ((17*I)/24)*Pi - 
          (19*\[Gamma]E)/12 - (55*Log[2])/12 - (19*Log[b])/12 - 
          (19*Log[x])/8) + (119/72 + ((17*I)/24)*Pi - (17*\[Gamma]E)/12 + 
          (55*Log[2])/12 - (27*Log[3])/4 - (17*Log[b])/12 - (17*Log[x])/8)/
         E^(I*\[ScriptL]) + E^((3*I)*\[ScriptL])*(49/72 + ((43*I)/72)*Pi - 
          (7*\[Gamma]E)/12 - (13*Log[2])/36 - (7*Log[b])/12 - 
          (7*Log[x])/8)) + e^4*(49/24 + ((83*I)/96)*Pi - (7*\[Gamma]E)/4 + 
        (59*Log[2])/12 - (243*Log[3])/32 - (7*Log[b])/4 + 
        (1981/144 + ((283*I)/48)*Pi - (283*\[Gamma]E)/24 + (135*Log[2])/8 + 
          (27*Log[3])/32 - (15625*Log[5])/576 - (283*Log[b])/24 - 
          (283*Log[x])/16)/E^((4*I)*\[ScriptL]) + E^((2*I)*\[ScriptL])*
         (7/3 + ((49*I)/72)*Pi - 2*\[Gamma]E - (52*Log[2])/9 - 2*Log[b] - 
          3*Log[x]) + (49/54 + ((7*I)/18)*Pi - (7*\[Gamma]E)/9 - 
          (293*Log[2])/9 + (135*Log[3])/8 - (7*Log[b])/9 - (7*Log[x])/6)/
         E^((2*I)*\[ScriptL]) + E^((4*I)*\[ScriptL])*
         (175/432 + ((175*I)/192)*Pi - (25*\[Gamma]E)/72 - Log[2]/8 + 
          (27*Log[3])/64 - (25*Log[b])/72 - (25*Log[x])/48) - 
        (21*Log[x])/8) + e^2*(14/9 + ((2*I)/3)*Pi - (4*\[Gamma]E)/3 - 
        4*Log[2] - (4*Log[b])/3 + (14/3 + (2*I)*Pi - 4*\[Gamma]E - 
          (4*Log[2])/3 - (27*Log[3])/4 - 4*Log[b] - 6*Log[x])/
         E^((2*I)*\[ScriptL]) + E^((2*I)*\[ScriptL])*(7/9 + ((5*I)/12)*Pi - 
          (2*\[Gamma]E)/3 - (2*Log[2*b])/3 - Log[x]) - 2*Log[x]) + 
      e^5*((25403/1152 + ((3629*I)/384)*Pi - (3629*\[Gamma]E)/192 - 
          (183623*Log[2])/2880 - (3843*Log[3])/80 + (15625*Log[5])/576 - 
          (3629*Log[b])/192 - (3629*Log[x])/128)/E^((5*I)*\[ScriptL]) + 
        E^((3*I)*\[ScriptL])*(3647/1152 + ((69*I)/128)*Pi - 
          (521*\[Gamma]E)/192 - (4223*Log[2])/576 - (27*Log[3])/64 - 
          (521*Log[b])/192 - (521*Log[x])/128) + 
        (3913/1728 + ((185*I)/192)*Pi - (559*\[Gamma]E)/288 - 
          (3697*Log[2])/96 + (297*Log[3])/16 - (559*Log[b])/288 - 
          (559*Log[x])/192)/E^(I*\[ScriptL]) + E^(I*\[ScriptL])*
         (427/192 + ((541*I)/576)*Pi - (61*\[Gamma]E)/32 + 
          (1819*Log[2])/288 - 9*Log[3] - (61*Log[b])/32 - (183*Log[x])/64) + 
        E^((5*I)*\[ScriptL])*(-581/3456 + ((8191*I)/5760)*Pi + 
          (83*\[Gamma]E)/576 + (5231*Log[2])/2880 + (27*Log[3])/64 + 
          (83*Log[b])/576 + (83*Log[x])/384) + 
        (-1771/1152 - ((253*I)/384)*Pi + (253*\[Gamma]E)/192 + 
          (15181*Log[2])/192 - (81*Log[3])/8 - (15625*Log[5])/576 + 
          (253*Log[b])/192 + (253*Log[x])/128)/E^((3*I)*\[ScriptL])) + 
      e^6*(175/72 + ((887*I)/864)*Pi - (25*\[Gamma]E)/12 - (553*Log[2])/12 + 
        (729*Log[3])/32 - (25*Log[b])/12 + (49609/1440 + ((7087*I)/480)*Pi - 
          (7087*\[Gamma]E)/240 + (46697*Log[2])/2160 + (125631*Log[3])/2560 + 
          (15625*Log[5])/4608 - (5764801*Log[7])/69120 - (7087*Log[b])/240 - 
          (7087*Log[x])/160)/E^((6*I)*\[ScriptL]) + E^((4*I)*\[ScriptL])*
         (9919/2160 + ((509*I)/2880)*Pi - (1417*\[Gamma]E)/360 - 
          (1331*Log[2])/120 - (189*Log[3])/320 - (1417*Log[b])/360 - 
          (1417*Log[x])/240) + (2261/864 + ((5141*I)/4608)*Pi - 
          (323*\[Gamma]E)/144 + (12101*Log[2])/144 - (7857*Log[3])/512 - 
          (15625*Log[5])/512 - (323*Log[b])/144 - (323*Log[x])/96)/
         E^((2*I)*\[ScriptL]) + E^((2*I)*\[ScriptL])*
         (2023/864 + ((2353*I)/2304)*Pi - (289*\[Gamma]E)/144 + 
          (1255*Log[2])/144 - (1413*Log[3])/128 - (289*Log[b])/144 - 
          (289*Log[x])/96) - (25*Log[x])/8 + E^((6*I)*\[ScriptL])*
         (-1757/1440 + ((153793*I)/69120)*Pi + (251*\[Gamma]E)/240 + 
          (5971*Log[2])/2160 + (243*Log[3])/512 + (15625*Log[5])/13824 + 
          (251*Log[b])/240 + (251*Log[x])/160) + 
        (-15337/2160 - ((2191*I)/720)*Pi + (2191*\[Gamma]E)/360 - 
          (741*Log[2])/8 - (8001*Log[3])/160 + (15625*Log[5])/192 + 
          (2191*Log[b])/360 + (2191*Log[x])/240)/E^((4*I)*\[ScriptL]))) + 
    SO*x^(7/2)*\[Epsilon]^5*(\[Chi]A*(-7/6 - (I/2)*Pi + \[Gamma]E + 
        e*(E^(I*\[ScriptL])*(-7/6 - (I/2)*Pi + \[Gamma]E + Log[2] + Log[b] + 
            (3*Log[x])/2) + (-7/2 - ((3*I)/2)*Pi + 3*\[Gamma]E + 7*Log[2] + 
            3*Log[b] + (9*Log[x])/2)/E^(I*\[ScriptL])) + 
        e^3*(E^((3*I)*\[ScriptL])*(-49/48 - ((43*I)/48)*Pi + 
            (7*\[Gamma]E)/8 + (13*Log[2])/24 + (7*Log[b])/8 + 
            (21*Log[x])/16) + E^(I*\[ScriptL])*(-161/48 - ((21*I)/16)*Pi + 
            (23*\[Gamma]E)/8 + (59*Log[2])/8 + (23*Log[b])/8 + 
            (69*Log[x])/16) + (-203/48 - ((29*I)/16)*Pi + (29*\[Gamma]E)/8 - 
            (27*Log[2])/8 + (81*Log[3])/8 + (29*Log[b])/8 + (87*Log[x])/16)/
           E^(I*\[ScriptL]) + (-595/48 - ((85*I)/16)*Pi + (85*\[Gamma]E)/8 + 
            (1267*Log[2])/24 - (81*Log[3])/8 + (85*Log[b])/8 + 
            (255*Log[x])/16)/E^((3*I)*\[ScriptL])) + 
        e^5*(E^((5*I)*\[ScriptL])*(581/2304 - ((8191*I)/3840)*Pi - 
            (83*\[Gamma]E)/384 - (5231*Log[2])/1920 - (81*Log[3])/128 - 
            (83*Log[b])/384 - (83*Log[x])/256) + 
          (-2989/768 - ((427*I)/256)*Pi + (427*\[Gamma]E)/128 - 
            (35407*Log[2])/384 + (81*Log[3])/8 + (15625*Log[5])/384 + 
            (427*Log[b])/128 + (1281*Log[x])/256)/E^((3*I)*\[ScriptL]) + 
          E^(I*\[ScriptL])*(-1981/384 - ((817*I)/384)*Pi + 
            (283*\[Gamma]E)/64 - (1087*Log[2])/192 + (27*Log[3])/2 + 
            (283*Log[b])/64 + (849*Log[x])/128) + E^((3*I)*\[ScriptL])*
           (-4039/768 - ((965*I)/768)*Pi + (577*\[Gamma]E)/128 + 
            (4327*Log[2])/384 + (81*Log[3])/128 + (577*Log[b])/128 + 
            (1731*Log[x])/256) + (-6853/1152 - ((325*I)/128)*Pi + 
            (979*\[Gamma]E)/192 + (3645*Log[2])/64 - (729*Log[3])/32 + 
            (979*Log[b])/192 + (979*Log[x])/128)/E^(I*\[ScriptL]) + 
          (-25403/768 - ((3629*I)/256)*Pi + (3629*\[Gamma]E)/128 + 
            (183623*Log[2])/1920 + (11529*Log[3])/160 - (15625*Log[5])/384 + 
            (3629*Log[b])/128 + (10887*Log[x])/256)/E^((5*I)*\[ScriptL])) + 
        e^6*(-77/12 - ((3133*I)/1152)*Pi + (11*\[Gamma]E)/2 + 68*Log[2] - 
          (3645*Log[3])/128 + (11*Log[b])/2 + E^((6*I)*\[ScriptL])*
           (1757/960 - ((153793*I)/46080)*Pi - (251*\[Gamma]E)/160 - 
            (5971*Log[2])/1440 - (729*Log[3])/1024 - (15625*Log[5])/9216 - 
            (251*Log[b])/160 - (753*Log[x])/320) + 
          (959/2880 + ((137*I)/960)*Pi - (137*\[Gamma]E)/480 + 
            (4041*Log[2])/32 + (47601*Log[3])/640 - (78125*Log[5])/768 - 
            (137*Log[b])/480 - (137*Log[x])/320)/E^((4*I)*\[ScriptL]) + 
          (33*Log[x])/4 + E^((2*I)*\[ScriptL])*(-3283/576 - 
            ((3497*I)/1536)*Pi + (469*\[Gamma]E)/96 - (803*Log[2])/96 + 
            (4239*Log[3])/256 + (469*Log[b])/96 + (469*Log[x])/64) + 
          E^((4*I)*\[ScriptL])*(-20713/2880 - ((3643*I)/3840)*Pi + 
            (2959*\[Gamma]E)/480 + (2677*Log[2])/160 + (729*Log[3])/1280 + 
            (2959*Log[b])/480 + (2959*Log[x])/320) + 
          (-4165/576 - ((9493*I)/3072)*Pi + (595*\[Gamma]E)/96 - 
            (9685*Log[2])/96 + (14499*Log[3])/1024 + (46875*Log[5])/1024 + 
            (595*Log[b])/96 + (595*Log[x])/64)/E^((2*I)*\[ScriptL]) + 
          (-49609/960 - ((7087*I)/320)*Pi + (7087*\[Gamma]E)/160 - 
            (46697*Log[2])/1440 - (376893*Log[3])/5120 - (15625*Log[5])/
             3072 + (5764801*Log[7])/46080 + (7087*Log[b])/160 + 
            (21261*Log[x])/320)/E^((6*I)*\[ScriptL])) + Log[2*b*x^(3/2)] + 
        e^4*(-14/3 - ((127*I)/64)*Pi + 4*\[Gamma]E + (729*Log[3])/64 + 
          E^((4*I)*\[ScriptL])*(-175/288 - ((175*I)/128)*Pi + 
            (25*\[Gamma]E)/48 + (3*Log[2])/16 - (81*Log[3])/128 + 
            (25*Log[b])/48 + (25*Log[x])/32) + E^((2*I)*\[ScriptL])*
           (-49/12 - ((4*I)/3)*Pi + (7*\[Gamma]E)/2 + (55*Log[2])/6 + 
            (7*Log[b])/2 + (21*Log[x])/4) + (-175/36 - ((25*I)/12)*Pi + 
            (25*\[Gamma]E)/6 + (299*Log[2])/6 - (81*Log[3])/4 + 
            (25*Log[b])/6 + (25*Log[x])/4)/E^((2*I)*\[ScriptL]) + 
          (-1981/96 - ((283*I)/32)*Pi + (283*\[Gamma]E)/16 - 
            (405*Log[2])/16 - (81*Log[3])/64 + (15625*Log[5])/384 + 
            (283*Log[b])/16 + (849*Log[x])/32)/E^((4*I)*\[ScriptL]) + 
          Log[(b^4*x^6)/16]) + e^2*(-35/12 - ((5*I)/4)*Pi + (5*\[Gamma]E)/2 + 
          (13*Log[2])/2 + (5*Log[b])/2 + (15*Log[x])/4 + E^((2*I)*\[ScriptL])*
           (-7/6 - ((5*I)/8)*Pi + \[Gamma]E + Log[2] + Log[b] + 
            (3*Log[x])/2) + (-7 - (3*I)*Pi + 6*\[Gamma]E + (81*Log[3])/8 + 
            Log[4*b^6*x^9])/E^((2*I)*\[ScriptL]))) + 
      \[Delta]*\[Chi]S*(-7/6 - (I/2)*Pi + \[Gamma]E + 
        e*(E^(I*\[ScriptL])*(-7/6 - (I/2)*Pi + \[Gamma]E + Log[2] + Log[b] + 
            (3*Log[x])/2) + (-7/2 - ((3*I)/2)*Pi + 3*\[Gamma]E + 7*Log[2] + 
            3*Log[b] + (9*Log[x])/2)/E^(I*\[ScriptL])) + 
        e^3*(E^((3*I)*\[ScriptL])*(-49/48 - ((43*I)/48)*Pi + 
            (7*\[Gamma]E)/8 + (13*Log[2])/24 + (7*Log[b])/8 + 
            (21*Log[x])/16) + E^(I*\[ScriptL])*(-161/48 - ((21*I)/16)*Pi + 
            (23*\[Gamma]E)/8 + (59*Log[2])/8 + (23*Log[b])/8 + 
            (69*Log[x])/16) + (-203/48 - ((29*I)/16)*Pi + (29*\[Gamma]E)/8 - 
            (27*Log[2])/8 + (81*Log[3])/8 + (29*Log[b])/8 + (87*Log[x])/16)/
           E^(I*\[ScriptL]) + (-595/48 - ((85*I)/16)*Pi + (85*\[Gamma]E)/8 + 
            (1267*Log[2])/24 - (81*Log[3])/8 + (85*Log[b])/8 + 
            (255*Log[x])/16)/E^((3*I)*\[ScriptL])) + 
        e^5*(E^((5*I)*\[ScriptL])*(581/2304 - ((8191*I)/3840)*Pi - 
            (83*\[Gamma]E)/384 - (5231*Log[2])/1920 - (81*Log[3])/128 - 
            (83*Log[b])/384 - (83*Log[x])/256) + 
          (-2989/768 - ((427*I)/256)*Pi + (427*\[Gamma]E)/128 - 
            (35407*Log[2])/384 + (81*Log[3])/8 + (15625*Log[5])/384 + 
            (427*Log[b])/128 + (1281*Log[x])/256)/E^((3*I)*\[ScriptL]) + 
          E^(I*\[ScriptL])*(-1981/384 - ((817*I)/384)*Pi + 
            (283*\[Gamma]E)/64 - (1087*Log[2])/192 + (27*Log[3])/2 + 
            (283*Log[b])/64 + (849*Log[x])/128) + E^((3*I)*\[ScriptL])*
           (-4039/768 - ((965*I)/768)*Pi + (577*\[Gamma]E)/128 + 
            (4327*Log[2])/384 + (81*Log[3])/128 + (577*Log[b])/128 + 
            (1731*Log[x])/256) + (-6853/1152 - ((325*I)/128)*Pi + 
            (979*\[Gamma]E)/192 + (3645*Log[2])/64 - (729*Log[3])/32 + 
            (979*Log[b])/192 + (979*Log[x])/128)/E^(I*\[ScriptL]) + 
          (-25403/768 - ((3629*I)/256)*Pi + (3629*\[Gamma]E)/128 + 
            (183623*Log[2])/1920 + (11529*Log[3])/160 - (15625*Log[5])/384 + 
            (3629*Log[b])/128 + (10887*Log[x])/256)/E^((5*I)*\[ScriptL])) + 
        e^6*(-77/12 - ((3133*I)/1152)*Pi + (11*\[Gamma]E)/2 + 68*Log[2] - 
          (3645*Log[3])/128 + (11*Log[b])/2 + E^((6*I)*\[ScriptL])*
           (1757/960 - ((153793*I)/46080)*Pi - (251*\[Gamma]E)/160 - 
            (5971*Log[2])/1440 - (729*Log[3])/1024 - (15625*Log[5])/9216 - 
            (251*Log[b])/160 - (753*Log[x])/320) + 
          (959/2880 + ((137*I)/960)*Pi - (137*\[Gamma]E)/480 + 
            (4041*Log[2])/32 + (47601*Log[3])/640 - (78125*Log[5])/768 - 
            (137*Log[b])/480 - (137*Log[x])/320)/E^((4*I)*\[ScriptL]) + 
          (33*Log[x])/4 + E^((2*I)*\[ScriptL])*(-3283/576 - 
            ((3497*I)/1536)*Pi + (469*\[Gamma]E)/96 - (803*Log[2])/96 + 
            (4239*Log[3])/256 + (469*Log[b])/96 + (469*Log[x])/64) + 
          E^((4*I)*\[ScriptL])*(-20713/2880 - ((3643*I)/3840)*Pi + 
            (2959*\[Gamma]E)/480 + (2677*Log[2])/160 + (729*Log[3])/1280 + 
            (2959*Log[b])/480 + (2959*Log[x])/320) + 
          (-4165/576 - ((9493*I)/3072)*Pi + (595*\[Gamma]E)/96 - 
            (9685*Log[2])/96 + (14499*Log[3])/1024 + (46875*Log[5])/1024 + 
            (595*Log[b])/96 + (595*Log[x])/64)/E^((2*I)*\[ScriptL]) + 
          (-49609/960 - ((7087*I)/320)*Pi + (7087*\[Gamma]E)/160 - 
            (46697*Log[2])/1440 - (376893*Log[3])/5120 - (15625*Log[5])/
             3072 + (5764801*Log[7])/46080 + (7087*Log[b])/160 + 
            (21261*Log[x])/320)/E^((6*I)*\[ScriptL])) + Log[2*b*x^(3/2)] + 
        e^4*(-14/3 - ((127*I)/64)*Pi + 4*\[Gamma]E + (729*Log[3])/64 + 
          E^((4*I)*\[ScriptL])*(-175/288 - ((175*I)/128)*Pi + 
            (25*\[Gamma]E)/48 + (3*Log[2])/16 - (81*Log[3])/128 + 
            (25*Log[b])/48 + (25*Log[x])/32) + E^((2*I)*\[ScriptL])*
           (-49/12 - ((4*I)/3)*Pi + (7*\[Gamma]E)/2 + (55*Log[2])/6 + 
            (7*Log[b])/2 + (21*Log[x])/4) + (-175/36 - ((25*I)/12)*Pi + 
            (25*\[Gamma]E)/6 + (299*Log[2])/6 - (81*Log[3])/4 + 
            (25*Log[b])/6 + (25*Log[x])/4)/E^((2*I)*\[ScriptL]) + 
          (-1981/96 - ((283*I)/32)*Pi + (283*\[Gamma]E)/16 - 
            (405*Log[2])/16 - (81*Log[3])/64 + (15625*Log[5])/384 + 
            (283*Log[b])/16 + (849*Log[x])/32)/E^((4*I)*\[ScriptL]) + 
          Log[(b^4*x^6)/16]) + e^2*(-35/12 - ((5*I)/4)*Pi + (5*\[Gamma]E)/2 + 
          (13*Log[2])/2 + (5*Log[b])/2 + (15*Log[x])/4 + E^((2*I)*\[ScriptL])*
           (-7/6 - ((5*I)/8)*Pi + \[Gamma]E + Log[2] + Log[b] + 
            (3*Log[x])/2) + (-7 - (3*I)*Pi + 6*\[Gamma]E + (81*Log[3])/8 + 
            Log[4*b^6*x^9])/E^((2*I)*\[ScriptL]))))
 
Htail[2, 2] = x^(5/2)*\[Epsilon]^3*((-11*I)/3 + 2*Pi + (4*I)*\[Gamma]E + 
      (8*I)*Log[2] + (4*I)*Log[b] + 
      e*(((-121*I)/24 + (11*Pi)/4 + ((11*I)/2)*\[Gamma]E - ((5*I)/2)*Log[2] + 
          ((27*I)/2)*Log[3] + ((11*I)/2)*Log[b] + ((33*I)/4)*Log[x])/
         E^(I*\[ScriptL]) + E^(I*\[ScriptL])*((-143*I)/24 + (13*Pi)/4 + 
          ((13*I)/2)*\[Gamma]E + ((29*I)/2)*Log[2] + ((13*I)/2)*Log[b] + 
          ((39*I)/4)*Log[x])) + 
      e^3*(E^(I*\[ScriptL])*((-649*I)/192 + (59*Pi)/32 + 
          ((59*I)/16)*\[Gamma]E - ((693*I)/16)*Log[2] + ((351*I)/8)*Log[3] + 
          ((59*I)/16)*Log[b] + ((177*I)/32)*Log[x]) + 
        ((-671*I)/192 + (61*Pi)/32 + ((61*I)/16)*\[Gamma]E + 
          ((2541*I)/16)*Log[2] - ((1377*I)/16)*Log[3] + ((61*I)/16)*Log[b] + 
          ((183*I)/32)*Log[x])/E^(I*\[ScriptL]) + 
        ((-2035*I)/192 + (185*Pi)/32 + ((185*I)/16)*\[Gamma]E - 
          ((5573*I)/48)*Log[2] + ((81*I)/8)*Log[3] + ((3125*I)/48)*Log[5] + 
          ((185*I)/16)*Log[b] + ((555*I)/32)*Log[x])/E^((3*I)*\[ScriptL]) + 
        E^((3*I)*\[ScriptL])*((-2629*I)/192 + (703*Pi)/96 + 
          ((239*I)/16)*\[Gamma]E + ((1661*I)/48)*Log[2] + 
          ((239*I)/16)*Log[b] + ((717*I)/32)*Log[x])) + 
      e^5*(((-1309*I)/3072 + (119*Pi)/512 + ((119*I)/256)*\[Gamma]E + 
          ((518197*I)/768)*Log[2] + ((14121*I)/64)*Log[3] - 
          ((334375*I)/768)*Log[5] + ((119*I)/256)*Log[b] + 
          ((357*I)/512)*Log[x])/E^((3*I)*\[ScriptL]) + E^((3*I)*\[ScriptL])*
         ((4037*I)/3072 + (469*Pi)/1536 - ((367*I)/256)*\[Gamma]E - 
          ((90781*I)/768)*Log[2] + ((3105*I)/32)*Log[3] - 
          ((367*I)/256)*Log[b] - ((1101*I)/512)*Log[x]) + 
        ((-17017*I)/4608 + (1547*Pi)/768 + ((1547*I)/384)*\[Gamma]E - 
          ((99479*I)/128)*Log[2] + ((25245*I)/128)*Log[3] + 
          ((40625*I)/192)*Log[5] + ((1547*I)/384)*Log[b] + 
          ((1547*I)/256)*Log[x])/E^(I*\[ScriptL]) + E^(I*\[ScriptL])*
         ((-17039*I)/4608 + (1465*Pi)/768 + ((1549*I)/384)*\[Gamma]E + 
          ((162653*I)/384)*Log[2] - ((15993*I)/64)*Log[3] + 
          ((1549*I)/384)*Log[b] + ((1549*I)/256)*Log[x]) + 
        ((-199661*I)/9216 + (18151*Pi)/1536 + ((18151*I)/768)*\[Gamma]E - 
          ((821453*I)/3840)*Log[2] - ((3879*I)/16)*Log[3] + 
          ((3125*I)/64)*Log[5] + ((823543*I)/3840)*Log[7] + 
          ((18151*I)/768)*Log[b] + ((18151*I)/512)*Log[x])/
         E^((5*I)*\[ScriptL]) + E^((5*I)*\[ScriptL])*((-268411*I)/9216 + 
          (37257*Pi)/2560 + ((24401*I)/768)*\[Gamma]E + ((94071*I)/1280)*
           Log[2] + ((459*I)/1280)*Log[3] + ((24401*I)/768)*Log[b] + 
          ((24401*I)/512)*Log[x])) + e^2*((-11*I)/3 + 2*Pi + 
        (4*I)*\[Gamma]E - (22*I)*Log[2] + (27*I)*Log[3] + (4*I)*Log[b] + 
        ((-22*I)/3 + 4*Pi + (8*I)*\[Gamma]E + (75*I)*Log[2] - (27*I)*Log[3] + 
          (8*I)*Log[b] + (12*I)*Log[x])/E^((2*I)*\[ScriptL]) + 
        E^((2*I)*\[ScriptL])*((-55*I)/6 + 5*Pi + (10*I)*\[Gamma]E + 
          (23*I)*Log[2] + (10*I)*Log[b] + (15*I)*Log[x]) + (6*I)*Log[x]) + 
      e^4*((-11*I)/3 + 2*Pi + (4*I)*\[Gamma]E + ((543*I)/2)*Log[2] - 
        ((621*I)/4)*Log[3] + (4*I)*Log[b] + E^((2*I)*\[ScriptL])*
         ((-143*I)/72 + (11*Pi)/8 + ((13*I)/6)*\[Gamma]E - 
          ((147*I)/2)*Log[2] + ((531*I)/8)*Log[3] + ((13*I)/6)*Log[b] + 
          ((13*I)/4)*Log[x]) + ((-187*I)/72 + (17*Pi)/12 + 
          ((17*I)/6)*\[Gamma]E - ((2521*I)/6)*Log[2] + ((351*I)/4)*Log[3] + 
          ((3125*I)/24)*Log[5] + ((17*I)/6)*Log[b] + ((17*I)/4)*Log[x])/
         E^((2*I)*\[ScriptL]) + ((-2189*I)/144 + (199*Pi)/24 + 
          ((199*I)/12)*\[Gamma]E + ((745*I)/4)*Log[2] + ((981*I)/8)*Log[3] - 
          ((3125*I)/24)*Log[5] + ((199*I)/12)*Log[b] + ((199*I)/8)*Log[x])/
         E^((4*I)*\[ScriptL]) + E^((4*I)*\[ScriptL])*((-2893*I)/144 + 
          (125*Pi)/12 + ((263*I)/12)*\[Gamma]E + ((611*I)/12)*Log[2] + 
          ((263*I)/12)*Log[b] + ((263*I)/8)*Log[x]) + (6*I)*Log[x]) + 
      e^6*((-11*I)/3 + (1145*Pi)/576 + (4*I)*\[Gamma]E - 
        ((45641*I)/36)*Log[2] + ((22329*I)/64)*Log[3] + 
        ((184375*I)/576)*Log[5] + (4*I)*Log[b] + 
        ((-2189*I)/576 + (199*Pi)/96 + ((199*I)/48)*\[Gamma]E + 
          ((22351*I)/16)*Log[2] + ((18747*I)/64)*Log[3] - 
          ((303125*I)/384)*Log[5] + ((199*I)/48)*Log[b] + 
          ((199*I)/32)*Log[x])/E^((2*I)*\[ScriptL]) + E^((2*I)*\[ScriptL])*
         ((-2233*I)/576 + (689*Pi)/384 + ((203*I)/48)*\[Gamma]E + 
          ((30605*I)/48)*Log[2] - ((49167*I)/128)*Log[3] + 
          ((203*I)/48)*Log[b] + ((203*I)/32)*Log[x]) + 
        ((176*I)/45 - (32*Pi)/15 - ((64*I)/15)*\[Gamma]E - 
          ((28649*I)/30)*Log[2] - ((538029*I)/640)*Log[3] + 
          ((59375*I)/128)*Log[5] + ((823543*I)/1920)*Log[7] - 
          ((64*I)/15)*Log[b] - ((32*I)/5)*Log[x])/E^((4*I)*\[ScriptL]) + 
        E^((4*I)*\[ScriptL])*((352*I)/45 - (427*Pi)/240 - 
          ((128*I)/15)*\[Gamma]E - ((5513*I)/30)*Log[2] + 
          ((44289*I)/320)*Log[3] - ((128*I)/15)*Log[b] - ((64*I)/5)*Log[x]) + 
        ((-29491*I)/960 + (2681*Pi)/160 + ((2681*I)/80)*\[Gamma]E + 
          ((175631*I)/144)*Log[2] + ((58563*I)/640)*Log[3] + 
          ((3125*I)/576)*Log[5] - ((823543*I)/1920)*Log[7] + 
          ((2681*I)/80)*Log[b] + ((8043*I)/160)*Log[x])/
         E^((6*I)*\[ScriptL]) + E^((6*I)*\[ScriptL])*((-40183*I)/960 + 
          (115751*Pi)/5760 + ((3653*I)/80)*\[Gamma]E + ((8467*I)/80)*Log[2] + 
          ((459*I)/640)*Log[3] + ((3653*I)/80)*Log[b] + ((10959*I)/160)*
           Log[x]) + (6*I)*Log[x]) + (6*I)*Log[x]) + 
    SO*x^4*\[Epsilon]^6*(\[Delta]*\[Chi]A*((44*I)/9 - (8*Pi)/3 - 
        ((16*I)/3)*\[Gamma]E - ((32*I)/3)*Log[2] - ((16*I)/3)*Log[b] + 
        e*(((338*I)/9 - (32*Pi)/3 - ((64*I)/3)*\[Gamma]E + 
            ((88*I)/3)*Log[2] - (72*I)*Log[3] - ((64*I)/3)*Log[b] - 
            (32*I)*Log[x])/E^(I*\[ScriptL]) + E^(I*\[ScriptL])*
           ((515*I)/18 - (37*Pi)/3 - ((74*I)/3)*\[Gamma]E - 
            ((226*I)/3)*Log[2] - ((74*I)/3)*Log[b] - (37*I)*Log[x])) + 
        e^3*(E^(I*\[ScriptL])*((1577*I)/18 - (155*Pi)/6 - 
            ((155*I)/3)*\[Gamma]E + (648*I)*Log[2] - ((1197*I)/2)*Log[3] - 
            ((155*I)/3)*Log[b] - ((155*I)/2)*Log[x]) + 
          ((13517*I)/144 - (643*Pi)/24 - ((643*I)/12)*\[Gamma]E - 
            ((8885*I)/4)*Log[2] + ((4887*I)/4)*Log[3] - ((643*I)/12)*Log[b] - 
            ((643*I)/8)*Log[x])/E^(I*\[ScriptL]) + 
          ((2033*I)/18 - (233*Pi)/6 - ((233*I)/3)*\[Gamma]E + 
            ((14660*I)/9)*Log[2] - ((459*I)/2)*Log[3] - ((12625*I)/18)*
             Log[5] - ((233*I)/3)*Log[b] - ((233*I)/2)*Log[x])/
           E^((3*I)*\[ScriptL]) + E^((3*I)*\[ScriptL])*((18163*I)/144 - 
            (4057*Pi)/72 - ((1373*I)/12)*\[Gamma]E - ((12707*I)/36)*Log[2] - 
            ((1373*I)/12)*Log[b] - ((1373*I)/8)*Log[x])) + 
        e^5*(E^((3*I)*\[ScriptL])*((86047*I)/576 - (8249*Pi)/144 - 
            ((965*I)/12)*\[Gamma]E + ((294775*I)/144)*Log[2] - 
            ((28161*I)/16)*Log[3] - ((965*I)/12)*Log[b] - 
            ((965*I)/8)*Log[x]) + E^(I*\[ScriptL])*((549695*I)/3456 - 
            (25105*Pi)/576 - ((26569*I)/288)*\[Gamma]E - ((2287877*I)/288)*
             Log[2] + ((37095*I)/8)*Log[3] - ((26569*I)/288)*Log[b] - 
            ((26569*I)/192)*Log[x]) + ((358813*I)/2304 - (18347*Pi)/384 - 
            ((18347*I)/192)*\[Gamma]E - ((7061101*I)/576)*Log[2] - 
            ((30843*I)/8)*Log[3] + ((4414375*I)/576)*Log[5] - 
            ((18347*I)/192)*Log[b] - ((18347*I)/128)*Log[x])/
           E^((3*I)*\[ScriptL]) + ((17861*I)/108 - (6905*Pi)/144 - 
            ((6905*I)/72)*\[Gamma]E + ((169853*I)/12)*Log[2] - 
            ((27081*I)/8)*Log[3] - ((581375*I)/144)*Log[5] - 
            ((6905*I)/72)*Log[b] - ((6905*I)/48)*Log[x])/E^(I*\[ScriptL]) + 
          ((479843*I)/1728 - (14897*Pi)/144 - ((14897*I)/72)*\[Gamma]E + 
            ((2945681*I)/720)*Log[2] + ((69333*I)/16)*Log[3] - 
            ((65875*I)/48)*Log[5] - ((2266201*I)/720)*Log[7] - 
            ((14897*I)/72)*Log[b] - ((14897*I)/48)*Log[x])/
           E^((5*I)*\[ScriptL]) + E^((5*I)*\[ScriptL])*((2588395*I)/6912 - 
            (294901*Pi)/1920 - ((200621*I)/576)*\[Gamma]E - 
            ((968783*I)/960)*Log[2] - ((1809*I)/320)*Log[3] - 
            ((200621*I)/576)*Log[b] - ((200621*I)/384)*Log[x])) + 
        e^2*((949*I)/18 - (47*Pi)/3 - ((94*I)/3)*\[Gamma]E + 
          ((850*I)/3)*Log[2] - (279*I)*Log[3] - ((94*I)/3)*Log[b] + 
          ((2449*I)/36 - (131*Pi)/6 - ((131*I)/3)*\[Gamma]E - 
            (633*I)*Log[2] + (279*I)*Log[3] - ((131*I)/3)*Log[b] - 
            ((131*I)/2)*Log[x])/E^((2*I)*\[ScriptL]) + E^((2*I)*\[ScriptL])*
           ((2357*I)/36 - (175*Pi)/6 - ((175*I)/3)*\[Gamma]E - 
            ((551*I)/3)*Log[2] - ((175*I)/3)*Log[b] - ((175*I)/2)*Log[x]) - 
          (47*I)*Log[x]) + e^4*((4199*I)/36 - (205*Pi)/6 - 
          ((205*I)/3)*\[Gamma]E - (4527*I)*Log[2] + ((5193*I)/2)*Log[3] - 
          ((205*I)/3)*Log[b] + E^((2*I)*\[ScriptL])*((13099*I)/108 - 
            (931*Pi)/24 - ((641*I)/9)*\[Gamma]E + ((3593*I)/3)*Log[2] - 
            ((8565*I)/8)*Log[3] - ((641*I)/9)*Log[b] - ((641*I)/6)*Log[x]) + 
          ((13517*I)/108 - (679*Pi)/18 - ((679*I)/9)*\[Gamma]E + 
            ((60797*I)/9)*Log[2] - (1521*I)*Log[3] - ((147875*I)/72)*Log[5] - 
            ((679*I)/9)*Log[b] - ((679*I)/6)*Log[x])/E^((2*I)*\[ScriptL]) + 
          ((38773*I)/216 - (2327*Pi)/36 - ((2327*I)/18)*\[Gamma]E - 
            ((5847*I)/2)*Log[2] - ((12531*I)/8)*Log[3] + ((147875*I)/72)*
             Log[5] - ((2327*I)/18)*Log[b] - ((2327*I)/12)*Log[x])/
           E^((4*I)*\[ScriptL]) + E^((4*I)*\[ScriptL])*((48161*I)/216 - 
            (6949*Pi)/72 - ((3691*I)/18)*\[Gamma]E - ((11059*I)/18)*Log[2] - 
            ((3691*I)/18)*Log[b] - ((3691*I)/12)*Log[x]) - 
          ((205*I)/2)*Log[x]) + e^6*((9277*I)/48 - (98177*Pi)/1728 - 
          ((455*I)/4)*\[Gamma]E + ((2730371*I)/108)*Log[2] - 
          ((403947*I)/64)*Log[3] - ((11930875*I)/1728)*Log[5] - 
          ((455*I)/4)*Log[b] + E^((4*I)*\[ScriptL])*((670999*I)/4320 - 
            (5777*Pi)/72 - ((22133*I)/360)*\[Gamma]E + ((1204859*I)/360)*
             Log[2] - ((175641*I)/64)*Log[3] - ((22133*I)/360)*Log[b] - 
            ((22133*I)/240)*Log[x]) + ((154219*I)/864 - (38773*Pi)/720 - 
            ((38773*I)/360)*\[Gamma]E + ((1315615*I)/72)*Log[2] + 
            ((10219689*I)/640)*Log[3] - ((3405625*I)/384)*Log[5] - 
            ((48612361*I)/5760)*Log[7] - ((38773*I)/360)*Log[b] - 
            ((38773*I)/240)*Log[x])/E^((4*I)*\[ScriptL]) + 
          E^((2*I)*\[ScriptL])*((346537*I)/1728 - (58835*Pi)/1152 - 
            ((17195*I)/144)*\[Gamma]E - ((1875131*I)/144)*Log[2] + 
            ((984627*I)/128)*Log[3] - ((17195*I)/144)*Log[b] - 
            ((17195*I)/96)*Log[x]) + ((351173*I)/1728 - (17695*Pi)/288 - 
            ((17695*I)/144)*\[Gamma]E - ((1267109*I)/48)*Log[2] - 
            ((449049*I)/64)*Log[3] + ((18231875*I)/1152)*Log[5] - 
            ((17695*I)/144)*Log[b] - ((17695*I)/96)*Log[x])/
           E^((2*I)*\[ScriptL]) + ((1214389*I)/2880 - (77423*Pi)/480 - 
            ((77423*I)/240)*\[Gamma]E - ((45187303*I)/2160)*Log[2] - 
            ((1763019*I)/640)*Log[3] - ((91625*I)/1728)*Log[5] + 
            ((48612361*I)/5760)*Log[7] - ((77423*I)/240)*Log[b] - 
            ((77423*I)/160)*Log[x])/E^((6*I)*\[ScriptL]) + 
          E^((6*I)*\[ScriptL])*((350165*I)/576 - (4024913*Pi)/17280 - 
            ((137003*I)/240)*\[Gamma]E - ((387979*I)/240)*Log[2] - 
            ((9531*I)/640)*Log[3] - ((137003*I)/240)*Log[b] - 
            ((137003*I)/160)*Log[x]) - ((1365*I)/8)*Log[x]) - (8*I)*Log[x]) + 
      \[Chi]S*(((-44*I)/9)*(-1 + \[Nu]) + (8*Pi*(-1 + \[Nu]))/3 + 
        ((16*I)/3)*\[Gamma]E*(-1 + \[Nu]) + ((32*I)/3)*(-1 + \[Nu])*Log[2] + 
        ((16*I)/3)*(-1 + \[Nu])*Log[b] + (8*I)*(-1 + \[Nu])*Log[x] + 
        e*(((2*Pi*(-16 + 7*\[Nu]))/3 + ((4*I)/3)*\[Gamma]E*(-16 + 7*\[Nu]) - 
            ((2*I)/9)*(-169 + 79*\[Nu]) - ((4*I)/3)*(-22 + 13*\[Nu])*Log[2] + 
            (36*I)*(-2 + \[Nu])*Log[3] + ((4*I)/3)*(-16 + 7*\[Nu])*Log[b] + 
            (2*I)*(-16 + 7*\[Nu])*Log[x])/E^(I*\[ScriptL]) + 
          E^(I*\[ScriptL])*((Pi*(-37 + 22*\[Nu]))/3 + ((2*I)/3)*\[Gamma]E*
             (-37 + 22*\[Nu]) - (I/18)*(-515 + 296*\[Nu]) + 
            ((2*I)/3)*(-113 + 62*\[Nu])*Log[2] + ((2*I)/3)*(-37 + 22*\[Nu])*
             Log[b] + I*(-37 + 22*\[Nu])*Log[x])) + 
        e^2*((Pi*(-47 + 26*\[Nu]))/3 + ((2*I)/3)*\[Gamma]E*(-47 + 26*\[Nu]) - 
          (I/18)*(-949 + 502*\[Nu]) - ((2*I)/3)*(-425 + 182*\[Nu])*Log[2] + 
          (9*I)*(-31 + 14*\[Nu])*Log[3] + ((2*I)/3)*(-47 + 26*\[Nu])*Log[b] + 
          I*(-47 + 26*\[Nu])*Log[x] + ((2449*I)/36 - (25*I)*\[Nu] + 
            Pi*(-131/6 + 6*\[Nu]) + (I/3)*\[Gamma]E*(-131 + 36*\[Nu]) + 
            I*(-633 + 262*\[Nu])*Log[2] - (9*I)*(-31 + 14*\[Nu])*Log[3] + 
            (I/3)*(-131 + 36*\[Nu])*Log[b] + (I/2)*(-131 + 36*\[Nu])*Log[x])/
           E^((2*I)*\[ScriptL]) + E^((2*I)*\[ScriptL])*
           ((I/36)*(2357 - 1184*\[Nu]) + (Pi*(-175 + 88*\[Nu]))/6 + 
            (I/3)*\[Gamma]E*(-175 + 88*\[Nu]) + ((19*I)/3)*(-29 + 14*\[Nu])*
             Log[2] + (I/3)*(-175 + 88*\[Nu])*Log[b] + 
            (I/2)*(-175 + 88*\[Nu])*Log[x])) + 
        e^3*(E^(I*\[ScriptL])*((5*Pi*(-124 + 65*\[Nu]))/24 + 
            ((5*I)/12)*\[Gamma]E*(-124 + 65*\[Nu]) - 
            (I/144)*(-12616 + 6473*\[Nu]) + (648*I - ((3235*I)/12)*\[Nu])*
             Log[2] + ((9*I)/4)*(-266 + 115*\[Nu])*Log[3] + 
            ((5*I)/12)*(-124 + 65*\[Nu])*Log[b] + ((5*I)/8)*(-124 + 65*\[Nu])*
             Log[x]) + ((Pi*(-932 + 155*\[Nu]))/24 + (I/12)*\[Gamma]E*
             (-932 + 155*\[Nu]) - (I/144)*(-16264 + 4711*\[Nu]) - 
            (I/36)*(-58640 + 24887*\[Nu])*Log[2] + ((27*I)/4)*
             (-34 + 15*\[Nu])*Log[3] + ((125*I)/18)*(-101 + 38*\[Nu])*
             Log[5] + (I/12)*(-932 + 155*\[Nu])*Log[b] + 
            (I/8)*(-932 + 155*\[Nu])*Log[x])/E^((3*I)*\[ScriptL]) + 
          ((Pi*(-643 + 303*\[Nu]))/24 + (I/12)*\[Gamma]E*(-643 + 303*\[Nu]) - 
            (I/144)*(-13517 + 6555*\[Nu]) + (I/4)*(-8885 + 3677*\[Nu])*
             Log[2] - ((27*I)/4)*(-181 + 74*\[Nu])*Log[3] + 
            (I/12)*(-643 + 303*\[Nu])*Log[b] + (I/8)*(-643 + 303*\[Nu])*
             Log[x])/E^(I*\[ScriptL]) + E^((3*I)*\[ScriptL])*
           (((-41*I)/144)*(-443 + 205*\[Nu]) + (I/12)*\[Gamma]E*
             (-1373 + 625*\[Nu]) + (Pi*(-4057 + 1879*\[Nu]))/72 + 
            ((97*I)/36)*(-131 + 59*\[Nu])*Log[2] + (I/12)*(-1373 + 625*\[Nu])*
             Log[b] + (I/8)*(-1373 + 625*\[Nu])*Log[x])) + 
        e^4*((I/36)*(4199 - 2072*\[Nu]) + (5*Pi*(-41 + 20*\[Nu]))/6 + 
          ((5*I)/3)*\[Gamma]E*(-41 + 20*\[Nu]) + 
          (-4527*I + ((5474*I)/3)*\[Nu])*Log[2] - ((9*I)/2)*
           (-577 + 230*\[Nu])*Log[3] + ((5*I)/3)*(-41 + 20*\[Nu])*Log[b] + 
          ((5*I)/2)*(-41 + 20*\[Nu])*Log[x] + (Pi*(-2327/36 + (47*\[Nu])/9) + 
            (I/18)*\[Gamma]E*(-2327 + 188*\[Nu]) - (I/216)*
             (-38773 + 8656*\[Nu]) + (I/6)*(-17541 + 6410*\[Nu])*Log[2] + 
            ((3*I)/8)*(-4177 + 1436*\[Nu])*Log[3] - ((125*I)/72)*
             (-1183 + 454*\[Nu])*Log[5] + (I/18)*(-2327 + 188*\[Nu])*Log[b] + 
            (I/12)*(-2327 + 188*\[Nu])*Log[x])/E^((4*I)*\[ScriptL]) + 
          ((Pi*(-679 + 313*\[Nu]))/18 + (I/9)*\[Gamma]E*(-679 + 313*\[Nu]) - 
            (I/108)*(-13517 + 6467*\[Nu]) - (I/9)*(-60797 + 23174*\[Nu])*
             Log[2] + (117*I)*(-13 + 5*\[Nu])*Log[3] + ((125*I)/72)*
             (-1183 + 454*\[Nu])*Log[5] + (I/9)*(-679 + 313*\[Nu])*Log[b] + 
            (I/6)*(-679 + 313*\[Nu])*Log[x])/E^((2*I)*\[ScriptL]) + 
          E^((2*I)*\[ScriptL])*((I/9)*\[Gamma]E*(-641 + 335*\[Nu]) + 
            (Pi*(-931 + 458*\[Nu]))/24 - (I/108)*(-13099 + 6709*\[Nu]) - 
            (I/3)*(-3593 + 1462*\[Nu])*Log[2] + ((3*I)/8)*
             (-2855 + 1198*\[Nu])*Log[3] + (I/9)*(-641 + 335*\[Nu])*Log[b] + 
            (I/6)*(-641 + 335*\[Nu])*Log[x]) + E^((4*I)*\[ScriptL])*
           ((I/18)*\[Gamma]E*(-3691 + 1564*\[Nu]) + (Pi*(-6949 + 3076*\[Nu]))/
             72 - (I/216)*(-48161 + 20984*\[Nu]) + 
            (I/18)*(-11059 + 4726*\[Nu])*Log[2] + (I/18)*(-3691 + 1564*\[Nu])*
             Log[b] + (I/12)*(-3691 + 1564*\[Nu])*Log[x])) + 
        e^6*((I/48)*(9277 - 4446*\[Nu]) + ((35*I)/4)*\[Gamma]E*
           (-13 + 6*\[Nu]) + (Pi*(-98177 + 45410*\[Nu]))/1728 - 
          (I/108)*(-2730371 + 1009202*\[Nu])*Log[2] + ((81*I)/64)*
           (-4987 + 1798*\[Nu])*Log[3] + ((125*I)/1728)*
           (-95447 + 36086*\[Nu])*Log[5] + ((35*I)/4)*(-13 + 6*\[Nu])*
           Log[b] + ((105*I)/8)*(-13 + 6*\[Nu])*Log[x] + 
          ((191*Pi*(-203 + 138*\[Nu]))/720 + ((191*I)/360)*\[Gamma]E*
             (-203 + 138*\[Nu]) - (I/4320)*(-771095 + 462234*\[Nu]) + 
            (I/360)*(6578075 - 2206338*\[Nu])*Log[2] - ((27*I)/640)*
             (-378507 + 126878*\[Nu])*Log[3] + ((625*I)/384)*
             (-5449 + 1842*\[Nu])*Log[5] + ((343*I)/5760)*
             (-141727 + 48054*\[Nu])*Log[7] + ((191*I)/360)*
             (-203 + 138*\[Nu])*Log[b] + ((191*I)/240)*(-203 + 138*\[Nu])*
             Log[x])/E^((4*I)*\[ScriptL]) + 
          ((5*Pi*(-3539 + 1510*\[Nu]))/288 + ((5*I)/144)*\[Gamma]E*
             (-3539 + 1510*\[Nu]) - (I/1728)*(-351173 + 161314*\[Nu]) + 
            (I/144)*(-3801327 + 1326734*\[Nu])*Log[2] + ((3*I)/64)*
             (-149683 + 55874*\[Nu])*Log[3] - ((625*I)/1152)*
             (-29171 + 10358*\[Nu])*Log[5] + ((5*I)/144)*(-3539 + 1510*\[Nu])*
             Log[b] + ((5*I)/96)*(-3539 + 1510*\[Nu])*Log[x])/
           E^((2*I)*\[ScriptL]) + ((Pi*(-77423 - 4294*\[Nu]))/480 - 
            (I/240)*\[Gamma]E*(77423 + 4294*\[Nu]) - ((11*I)/2880)*
             (-110399 + 12194*\[Nu]) + (I/2160)*(-45187303 + 13670374*\[Nu])*
             Log[2] + ((27*I)/640)*(-65297 + 25426*\[Nu])*Log[3] + 
            ((125*I)/1728)*(-733 + 154*\[Nu])*Log[5] - ((343*I)/5760)*
             (-141727 + 48054*\[Nu])*Log[7] - (I/240)*(77423 + 4294*\[Nu])*
             Log[b] - (I/160)*(77423 + 4294*\[Nu])*Log[x])/
           E^((6*I)*\[ScriptL]) + E^((2*I)*\[ScriptL])*
           ((I/144)*\[Gamma]E*(-17195 + 8142*\[Nu]) + 
            (Pi*(-58835 + 28902*\[Nu]))/1152 - (I/1728)*(-346537 + 
              168258*\[Nu]) + (I/144)*(-1875131 + 729054*\[Nu])*Log[2] - 
            ((9*I)/128)*(-109403 + 42190*\[Nu])*Log[3] + 
            (I/144)*(-17195 + 8142*\[Nu])*Log[b] + 
            (I/96)*(-17195 + 8142*\[Nu])*Log[x]) + E^((4*I)*\[ScriptL])*
           ((I/360)*\[Gamma]E*(-22133 + 16246*\[Nu]) + 
            (Pi*(-57770 + 24989*\[Nu]))/720 - (I/4320)*(-670999 + 
              392474*\[Nu]) - (I/360)*(-1204859 + 471362*\[Nu])*Log[2] + 
            ((3*I)/320)*(-292735 + 117986*\[Nu])*Log[3] + 
            (I/360)*(-22133 + 16246*\[Nu])*Log[b] + 
            (I/240)*(-22133 + 16246*\[Nu])*Log[x]) + E^((6*I)*\[ScriptL])*
           ((I/240)*\[Gamma]E*(-137003 + 52338*\[Nu]) - 
            (I/2880)*(-1750825 + 697614*\[Nu]) + 
            (Pi*(-4024913 + 1708922*\[Nu]))/17280 + 
            (I/720)*(-1163937 + 459074*\[Nu])*Log[2] + ((27*I)/640)*
             (-353 + 66*\[Nu])*Log[3] + (I/240)*(-137003 + 52338*\[Nu])*
             Log[b] + (I/160)*(-137003 + 52338*\[Nu])*Log[x])) + 
        e^5*(((Pi*(-29794 + 249*\[Nu]))/288 + (I/144)*\[Gamma]E*
             (-29794 + 249*\[Nu]) - (I/3456)*(-959686 + 157587*\[Nu]) - 
            (I/720)*(-2945681 + 1117149*\[Nu])*Log[2] - ((33*I)/32)*
             (-4202 + 1513*\[Nu])*Log[3] + ((125*I)/96)*(-1054 + 427*\[Nu])*
             Log[5] + ((343*I)/720)*(-6607 + 2103*\[Nu])*Log[7] + 
            (I/144)*(-29794 + 249*\[Nu])*Log[b] + (I/96)*(-29794 + 249*\[Nu])*
             Log[x])/E^((5*I)*\[ScriptL]) + E^((3*I)*\[ScriptL])*
           (((5*I)/48)*\[Gamma]E*(-772 + 431*\[Nu]) + 
            (Pi*(-16498 + 7609*\[Nu]))/288 - (I/1152)*(-172094 + 
              90997*\[Nu]) - ((5*I)/144)*(-58955 + 23477*\[Nu])*Log[2] + 
            ((27*I)/32)*(-2086 + 855*\[Nu])*Log[3] + ((5*I)/48)*
             (-772 + 431*\[Nu])*Log[b] + ((5*I)/32)*(-772 + 431*\[Nu])*
             Log[x]) + ((Pi*(-6905 + 3107*\[Nu]))/144 + (I/72)*\[Gamma]E*
             (-6905 + 3107*\[Nu]) - (I/1728)*(-285776 + 135287*\[Nu]) - 
            (I/24)*(-339706 + 127061*\[Nu])*Log[2] + ((27*I)/16)*
             (-2006 + 735*\[Nu])*Log[3] + ((125*I)/288)*(-9302 + 3551*\[Nu])*
             Log[5] + (I/72)*(-6905 + 3107*\[Nu])*Log[b] + 
            (I/48)*(-6905 + 3107*\[Nu])*Log[x])/E^(I*\[ScriptL]) + 
          ((Pi*(-18347 + 9440*\[Nu]))/384 + (I/192)*\[Gamma]E*
             (-18347 + 9440*\[Nu]) - ((13*I)/2304)*(-27601 + 14026*\[Nu]) + 
            (I/576)*(-7061101 + 2551240*\[Nu])*Log[2] + ((9*I)/32)*
             (-13708 + 4913*\[Nu])*Log[3] - ((625*I)/576)*
             (-7063 + 2524*\[Nu])*Log[5] + (I/192)*(-18347 + 9440*\[Nu])*
             Log[b] + (I/128)*(-18347 + 9440*\[Nu])*Log[x])/
           E^((3*I)*\[ScriptL]) + E^(I*\[ScriptL])*
           ((Pi*(-25105 + 12306*\[Nu]))/576 + (I/288)*\[Gamma]E*
             (-26569 + 12744*\[Nu]) - (I/3456)*(-549695 + 268902*\[Nu]) + 
            (I/288)*(-2287877 + 904272*\[Nu])*Log[2] - ((3*I)/32)*
             (-49460 + 19363*\[Nu])*Log[3] + (I/288)*(-26569 + 12744*\[Nu])*
             Log[b] + (I/192)*(-26569 + 12744*\[Nu])*Log[x]) + 
          E^((5*I)*\[ScriptL])*((I/576)*\[Gamma]E*(-200621 + 80336*\[Nu]) + 
            (Pi*(-294901 + 127124*\[Nu]))/1920 - (I/6912)*(-2588395 + 
              1074478*\[Nu]) + ((23*I)/960)*(-42121 + 17304*\[Nu])*Log[2] + 
            ((27*I)/320)*(-67 + 8*\[Nu])*Log[3] + (I/576)*
             (-200621 + 80336*\[Nu])*Log[b] + (I/384)*(-200621 + 80336*\[Nu])*
             Log[x]))))
 
Htail[3, 0] = SO*x^4*\[Epsilon]^6*\[Chi]S*
    (e*(E^(I*\[ScriptL])*(((-5*I)/3)*Sqrt[2/21]*\[Nu] - (Pi*\[Nu])/Sqrt[42] + 
         I*Sqrt[2/21]*\[Gamma]E*\[Nu] + I*Sqrt[2/21]*\[Nu]*Log[2*b] + 
         I*Sqrt[3/14]*\[Nu]*Log[x]) + (((-5*I)/3)*Sqrt[2/21]*\[Nu] + 
         (Pi*\[Nu])/Sqrt[42] + I*Sqrt[2/21]*\[Gamma]E*\[Nu] + 
         I*Sqrt[2/21]*\[Nu]*Log[2*b] + I*Sqrt[3/14]*\[Nu]*Log[x])/
        E^(I*\[ScriptL])) + 
     e^3*(((((5*I)/12)*\[Nu])/Sqrt[42] - (Pi*\[Nu])/(8*Sqrt[42]) - 
         ((I/4)*\[Gamma]E*\[Nu])/Sqrt[42] - ((I/4)*\[Nu]*Log[2*b])/Sqrt[42] - 
         (I/8)*Sqrt[3/14]*\[Nu]*Log[x])/E^(I*\[ScriptL]) + 
       E^(I*\[ScriptL])*((((5*I)/12)*\[Nu])/Sqrt[42] + 
         (Pi*\[Nu])/(8*Sqrt[42]) - ((I/4)*\[Gamma]E*\[Nu])/Sqrt[42] - 
         ((I/4)*\[Nu]*Log[2*b])/Sqrt[42] - (I/8)*Sqrt[3/14]*\[Nu]*Log[x]) + 
       E^((3*I)*\[ScriptL])*(((-45*I)/4)*Sqrt[3/14]*\[Nu] - 
         (27*Sqrt[3/14]*Pi*\[Nu])/8 + ((27*I)/4)*Sqrt[3/14]*\[Gamma]E*\[Nu] + 
         ((27*I)/4)*Sqrt[3/14]*\[Nu]*Log[6*b] + ((81*I)/8)*Sqrt[3/14]*\[Nu]*
          Log[x]) + (((-45*I)/4)*Sqrt[3/14]*\[Nu] + (27*Sqrt[3/14]*Pi*\[Nu])/
          8 + ((27*I)/4)*Sqrt[3/14]*\[Gamma]E*\[Nu] + ((27*I)/4)*Sqrt[3/14]*
          \[Nu]*Log[6*b] + ((81*I)/8)*Sqrt[3/14]*\[Nu]*Log[x])/
        E^((3*I)*\[ScriptL])) + 
     e^6*(E^((4*I)*\[ScriptL])*(((256*I)/9)*Sqrt[2/21]*\[Nu] + 
         (128*Sqrt[2/21]*Pi*\[Nu])/15 - ((256*I)/15)*Sqrt[2/21]*\[Gamma]E*
          \[Nu] - ((256*I)/5)*Sqrt[2/21]*\[Nu]*Log[2] - 
         ((256*I)/15)*Sqrt[2/21]*\[Nu]*Log[b] - ((128*I)/5)*Sqrt[2/21]*\[Nu]*
          Log[x]) + (((256*I)/9)*Sqrt[2/21]*\[Nu] - (128*Sqrt[2/21]*Pi*\[Nu])/
          15 - ((256*I)/15)*Sqrt[2/21]*\[Gamma]E*\[Nu] - 
         ((256*I)/45)*Sqrt[2/21]*\[Nu]*Log[512] - ((256*I)/15)*Sqrt[2/21]*
          \[Nu]*Log[b] - ((128*I)/5)*Sqrt[2/21]*\[Nu]*Log[x])/
        E^((4*I)*\[ScriptL]) + ((-81*I)*Sqrt[3/14]*\[Nu] + 
         (243*Sqrt[3/14]*Pi*\[Nu])/10 + ((243*I)/5)*Sqrt[3/14]*\[Gamma]E*
          \[Nu] + ((81*I)/5)*Sqrt[3/14]*\[Nu]*Log[64] + 
         ((243*I)/5)*Sqrt[3/14]*\[Nu]*Log[3*b] + ((729*I)/10)*Sqrt[3/14]*
          \[Nu]*Log[x])/E^((6*I)*\[ScriptL]) + E^((6*I)*\[ScriptL])*
        ((-81*I)*Sqrt[3/14]*\[Nu] - (243*Sqrt[3/14]*Pi*\[Nu])/10 + 
         ((243*I)/5)*Sqrt[3/14]*\[Gamma]E*\[Nu] + ((81*I)/10)*Sqrt[3/14]*
          \[Nu]*Log[4096] + ((243*I)/5)*Sqrt[3/14]*\[Nu]*Log[3*b] + 
         ((729*I)/10)*Sqrt[3/14]*\[Nu]*Log[x]) + E^((2*I)*\[ScriptL])*
        ((((-5*I)/9)*\[Nu])/Sqrt[42] - (Pi*\[Nu])/(6*Sqrt[42]) + 
         ((I/3)*\[Gamma]E*\[Nu])/Sqrt[42] + ((I/9)*\[Nu]*Log[64])/Sqrt[42] + 
         ((I/3)*\[Nu]*Log[b])/Sqrt[42] + ((I/2)*\[Nu]*Log[x])/Sqrt[42]) + 
       ((((-5*I)/9)*\[Nu])/Sqrt[42] + (Pi*\[Nu])/(6*Sqrt[42]) + 
         ((I/3)*\[Gamma]E*\[Nu])/Sqrt[42] + ((I/9)*\[Nu]*Log[64])/Sqrt[42] + 
         ((I/3)*\[Nu]*Log[b])/Sqrt[42] + ((I/2)*\[Nu]*Log[x])/Sqrt[42])/
        E^((2*I)*\[ScriptL])) + 
     e^5*((((405*I)/64)*Sqrt[3/14]*\[Nu] - (243*Sqrt[3/14]*Pi*\[Nu])/128 - 
         ((243*I)/64)*Sqrt[3/14]*\[Gamma]E*\[Nu] - ((243*I)/64)*Sqrt[3/14]*
          \[Nu]*Log[6*b] - ((729*I)/128)*Sqrt[3/14]*\[Nu]*Log[x])/
        E^((3*I)*\[ScriptL]) + E^((3*I)*\[ScriptL])*
        (((405*I)/64)*Sqrt[3/14]*\[Nu] + (243*Sqrt[3/14]*Pi*\[Nu])/128 - 
         ((243*I)/64)*Sqrt[3/14]*\[Gamma]E*\[Nu] - ((243*I)/64)*Sqrt[3/14]*
          \[Nu]*Log[6*b] - ((729*I)/128)*Sqrt[3/14]*\[Nu]*Log[x]) + 
       E^(I*\[ScriptL])*((((-5*I)/288)*\[Nu])/Sqrt[42] - 
         (Pi*\[Nu])/(192*Sqrt[42]) + ((I/96)*\[Gamma]E*\[Nu])/Sqrt[42] + 
         ((I/96)*\[Nu]*Log[2*b])/Sqrt[42] + ((I/64)*\[Nu]*Log[x])/Sqrt[42]) + 
       ((((-5*I)/288)*\[Nu])/Sqrt[42] + (Pi*\[Nu])/(192*Sqrt[42]) + 
         ((I/96)*\[Gamma]E*\[Nu])/Sqrt[42] + ((I/96)*\[Nu]*Log[2*b])/
          Sqrt[42] + ((I/64)*\[Nu]*Log[x])/Sqrt[42])/E^(I*\[ScriptL]) + 
       E^((5*I)*\[ScriptL])*((((-78125*I)/576)*\[Nu])/Sqrt[42] - 
         (15625*Pi*\[Nu])/(384*Sqrt[42]) + (((15625*I)/192)*\[Gamma]E*\[Nu])/
          Sqrt[42] + (((15625*I)/192)*\[Nu]*Log[10*b])/Sqrt[42] + 
         (((15625*I)/128)*\[Nu]*Log[x])/Sqrt[42]) + 
       ((((-78125*I)/576)*\[Nu])/Sqrt[42] + (15625*Pi*\[Nu])/(384*Sqrt[42]) + 
         (((15625*I)/192)*\[Gamma]E*\[Nu])/Sqrt[42] + 
         (((15625*I)/192)*\[Nu]*Log[10*b])/Sqrt[42] + 
         (((15625*I)/128)*\[Nu]*Log[x])/Sqrt[42])/E^((5*I)*\[ScriptL])) + 
     e^4*((((20*I)/9)*Sqrt[2/21]*\[Nu] - (2*Sqrt[2/21]*Pi*\[Nu])/3 - 
         ((4*I)/3)*Sqrt[2/21]*\[Gamma]E*\[Nu] - ((4*I)/9)*Sqrt[2/21]*\[Nu]*
          Log[64] - ((4*I)/3)*Sqrt[2/21]*\[Nu]*Log[b] - 
         (2*I)*Sqrt[2/21]*\[Nu]*Log[x])/E^((2*I)*\[ScriptL]) + 
       E^((4*I)*\[ScriptL])*(((-320*I)/9)*Sqrt[2/21]*\[Nu] - 
         (32*Sqrt[2/21]*Pi*\[Nu])/3 + ((64*I)/3)*Sqrt[2/21]*\[Gamma]E*\[Nu] + 
         (64*I)*Sqrt[2/21]*\[Nu]*Log[2] + ((64*I)/3)*Sqrt[2/21]*\[Nu]*
          Log[b] + (32*I)*Sqrt[2/21]*\[Nu]*Log[x]) + 
       (((-320*I)/9)*Sqrt[2/21]*\[Nu] + (32*Sqrt[2/21]*Pi*\[Nu])/3 + 
         ((64*I)/3)*Sqrt[2/21]*\[Gamma]E*\[Nu] + ((64*I)/9)*Sqrt[2/21]*\[Nu]*
          Log[512] + ((64*I)/3)*Sqrt[2/21]*\[Nu]*Log[b] + 
         (32*I)*Sqrt[2/21]*\[Nu]*Log[x])/E^((4*I)*\[ScriptL]) + 
       E^((2*I)*\[ScriptL])*(((20*I)/9)*Sqrt[2/21]*\[Nu] + 
         (2*Sqrt[2/21]*Pi*\[Nu])/3 - ((4*I)/3)*Sqrt[2/21]*\[Gamma]E*\[Nu] - 
         ((2*I)/9)*Sqrt[2/21]*\[Nu]*Log[4096*b^6*x^9])) + 
     e^2*((((-20*I)/3)*Sqrt[2/21]*\[Nu] + 2*Sqrt[2/21]*Pi*\[Nu] + 
         (4*I)*Sqrt[2/21]*\[Gamma]E*\[Nu] + ((4*I)/3)*Sqrt[2/21]*\[Nu]*
          Log[64] + (4*I)*Sqrt[2/21]*\[Nu]*Log[b] + (2*I)*Sqrt[6/7]*\[Nu]*
          Log[x])/E^((2*I)*\[ScriptL]) + E^((2*I)*\[ScriptL])*
        (((-20*I)/3)*Sqrt[2/21]*\[Nu] - 2*Sqrt[2/21]*Pi*\[Nu] + 
         (4*I)*Sqrt[2/21]*\[Gamma]E*\[Nu] + ((2*I)/3)*Sqrt[2/21]*\[Nu]*
          Log[4096*b^6*x^9])))
 
Htail[3, 1] = x^3*\[Delta]*\[Epsilon]^4*(97/(360*Sqrt[14]) + 
     ((I/12)*Pi)/Sqrt[14] - \[Gamma]E/(6*Sqrt[14]) - Log[2*b]/(6*Sqrt[14]) - 
     Log[x]/(4*Sqrt[14]) + e*(E^(I*\[ScriptL])*(97/(360*Sqrt[14]) + 
         ((I/12)*Pi)/Sqrt[14] - \[Gamma]E/(6*Sqrt[14]) - 
         Log[2*b]/(6*Sqrt[14]) - Log[x]/(4*Sqrt[14])) + 
       (-97/(40*Sqrt[14]) - (((3*I)/4)*Pi)/Sqrt[14] + (3*\[Gamma]E)/
          (2*Sqrt[14]) + (17*Log[2])/(6*Sqrt[14]) + (3*Log[b])/(2*Sqrt[14]) + 
         (9*Log[x])/(4*Sqrt[14]))/E^(I*\[ScriptL])) + 
     e^6*(-1067/(288*Sqrt[14]) - (((1919*I)/1728)*Pi)/Sqrt[14] + 
       (55*\[Gamma]E)/(24*Sqrt[14]) + (967*Log[2])/(24*Sqrt[14]) - 
       (1215*Log[3])/(64*Sqrt[14]) + (55*Log[b])/(24*Sqrt[14]) + 
       (55*Log[x])/(16*Sqrt[14]) + E^((4*I)*\[ScriptL])*
        ((-57133*Sqrt[7/2])/43200 + (((47*I)/288)*Pi)/Sqrt[14] + 
         (589*Sqrt[7/2]*\[Gamma]E)/720 + (2717*Log[2])/(240*Sqrt[14]) + 
         (297*Log[3])/(320*Sqrt[14]) + (589*Sqrt[7/2]*Log[b])/720 + 
         (589*Sqrt[7/2]*Log[x])/480) + ((-220093*Sqrt[7/2])/28800 - 
         ((2269*I)/960)*Sqrt[7/2]*Pi + (2269*Sqrt[7/2]*\[Gamma]E)/480 - 
         (10907*Sqrt[7/2]*Log[2])/4320 - (501309*Log[3])/(10240*Sqrt[14]) - 
         (59375*Log[5])/(18432*Sqrt[14]) + (3411821*Sqrt[7/2]*Log[7])/
          276480 + (2269*Sqrt[7/2]*Log[b])/480 + (2269*Sqrt[7/2]*Log[x])/320)/
        E^((6*I)*\[ScriptL]) + (475009/(43200*Sqrt[14]) + 
         (((4897*I)/1440)*Pi)/Sqrt[14] - (4897*\[Gamma]E)/(720*Sqrt[14]) + 
         (1455*Log[2])/(16*Sqrt[14]) + (15957*Log[3])/(320*Sqrt[14]) - 
         (15625*Log[5])/(192*Sqrt[14]) - (4897*Log[b])/(720*Sqrt[14]) - 
         (4897*Log[x])/(480*Sqrt[14]))/E^((4*I)*\[ScriptL]) + 
       E^((6*I)*\[ScriptL])*(119407/(14400*Sqrt[14]) - 
         (((591407*I)/276480)*Pi)/Sqrt[14] - (1231*\[Gamma]E)/
          (240*Sqrt[14]) - (17351*Log[2])/(2160*Sqrt[14]) - 
         (243*Sqrt[7/2]*Log[3])/2048 - (96875*Log[5])/(55296*Sqrt[14]) - 
         (1231*Log[b])/(240*Sqrt[14]) - (1231*Log[x])/(160*Sqrt[14])) + 
       E^((2*I)*\[ScriptL])*(-15229/(4320*Sqrt[14]) - (((9593*I)/9216)*Pi)/
          Sqrt[14] + (157*\[Gamma]E)/(72*Sqrt[14]) - (325*Log[2])/
          (72*Sqrt[14]) + (2133*Log[3])/(256*Sqrt[14]) + 
         (157*Log[b])/(72*Sqrt[14]) + (157*Log[x])/(48*Sqrt[14])) + 
       (-69161/(17280*Sqrt[14]) - ((3217*I)/18432)*Sqrt[7/2]*Pi + 
         (713*\[Gamma]E)/(288*Sqrt[14]) - (22607*Log[2])/(288*Sqrt[14]) + 
         (29619*Log[3])/(2048*Sqrt[14]) + (59375*Log[5])/(2048*Sqrt[14]) + 
         (713*Log[b])/(288*Sqrt[14]) + (713*Log[x])/(192*Sqrt[14]))/
        E^((2*I)*\[ScriptL])) + e^2*((-97*Sqrt[7/2])/360 - 
       (I/12)*Sqrt[7/2]*Pi + (Sqrt[7/2]*\[Gamma]E)/6 + 
       (5*Log[2])/(2*Sqrt[14]) + (Sqrt[7/2]*Log[b])/6 + 
       (Sqrt[7/2]*Log[x])/4 + E^((2*I)*\[ScriptL])*(97/(144*Sqrt[14]) - 
         ((I/48)*Pi)/Sqrt[14] - (5*\[Gamma]E)/(12*Sqrt[14]) - 
         (5*Log[2*b])/(12*Sqrt[14]) - (5*Log[x])/(8*Sqrt[14])) + 
       (-97/(16*Sqrt[14]) - (((15*I)/8)*Pi)/Sqrt[14] + 
         (15*\[Gamma]E)/(4*Sqrt[14]) + (29*Log[2])/(12*Sqrt[14]) + 
         (81*Log[3])/(16*Sqrt[14]) + (15*Log[b])/(4*Sqrt[14]) + 
         (45*Log[x])/(8*Sqrt[14]))/E^((2*I)*\[ScriptL])) + 
     e^3*(E^(I*\[ScriptL])*((-1067*Sqrt[7/2])/2880 - (((55*I)/96)*Pi)/
          Sqrt[14] + (11*Sqrt[7/2]*\[Gamma]E)/48 + (149*Log[2])/
          (48*Sqrt[14]) + (11*Sqrt[7/2]*Log[b])/48 + (11*Sqrt[7/2]*Log[x])/
          32) + E^((3*I)*\[ScriptL])*(4171/(2880*Sqrt[14]) - 
         (((65*I)/288)*Pi)/Sqrt[14] - (43*\[Gamma]E)/(48*Sqrt[14]) - 
         (193*Log[2])/(144*Sqrt[14]) - (43*Log[b])/(48*Sqrt[14]) - 
         (43*Log[x])/(32*Sqrt[14])) + (-6499/(2880*Sqrt[14]) - 
         (((67*I)/96)*Pi)/Sqrt[14] + (67*\[Gamma]E)/(48*Sqrt[14]) - 
         (125*Log[2])/(48*Sqrt[14]) + (81*Log[3])/(16*Sqrt[14]) + 
         (67*Log[b])/(48*Sqrt[14]) + (67*Log[x])/(32*Sqrt[14]))/
        E^(I*\[ScriptL]) + (-33659/(2880*Sqrt[14]) - (((347*I)/96)*Pi)/
          Sqrt[14] + (347*\[Gamma]E)/(48*Sqrt[14]) + (4601*Log[2])/
          (144*Sqrt[14]) - (81*Log[3])/(16*Sqrt[14]) + 
         (347*Log[b])/(48*Sqrt[14]) + (347*Log[x])/(32*Sqrt[14]))/
        E^((3*I)*\[ScriptL])) + e^4*(-2813/(960*Sqrt[14]) - 
       (((337*I)/384)*Pi)/Sqrt[14] + (29*\[Gamma]E)/(16*Sqrt[14]) - 
       (121*Log[2])/(48*Sqrt[14]) + (729*Log[3])/(128*Sqrt[14]) + 
       (29*Log[b])/(16*Sqrt[14]) + (87*Log[x])/(32*Sqrt[14]) + 
       E^((2*I)*\[ScriptL])*((-97*Sqrt[7/2])/180 - (((149*I)/288)*Pi)/
          Sqrt[14] + (Sqrt[7/2]*\[Gamma]E)/3 + (41*Log[2])/(9*Sqrt[14]) + 
         (Sqrt[7/2]*Log[b])/3 + (Sqrt[7/2]*Log[x])/2) + 
       E^((4*I)*\[ScriptL])*(47821/(17280*Sqrt[14]) - (((449*I)/768)*Pi)/
          Sqrt[14] - (493*\[Gamma]E)/(288*Sqrt[14]) - 
         (69*Log[2])/(32*Sqrt[14]) - (27*Sqrt[7/2]*Log[3])/256 - 
         (493*Log[b])/(288*Sqrt[14]) - (493*Log[x])/(192*Sqrt[14])) + 
       (-2813/(2160*Sqrt[14]) - (((29*I)/72)*Pi)/Sqrt[14] + 
         (29*\[Gamma]E)/(36*Sqrt[14]) + (1757*Sqrt[7/2]*Log[2])/288 - 
         (459*Log[6])/(32*Sqrt[14]) + (29*Log[b])/(36*Sqrt[14]) + 
         (29*Log[x])/(24*Sqrt[14]))/E^((2*I)*\[ScriptL]) + 
       (-39091/(1920*Sqrt[14]) - (((403*I)/64)*Pi)/Sqrt[14] + 
         (403*\[Gamma]E)/(32*Sqrt[14]) - (397*Log[2])/(32*Sqrt[14]) - 
         (81*Log[3])/(128*Sqrt[14]) + (59375*Log[5])/(2304*Sqrt[14]) + 
         (403*Log[b])/(32*Sqrt[14]) + (1209*Log[x])/(64*Sqrt[14]))/
        E^((4*I)*\[ScriptL])) + 
     e^5*(E^((5*I)*\[ScriptL])*(678709/(138240*Sqrt[14]) - 
         ((3887*I)/23040)*Sqrt[7/2]*Pi - (6997*\[Gamma]E)/(2304*Sqrt[14]) - 
         (67369*Log[2])/(11520*Sqrt[14]) - (27*Sqrt[7/2]*Log[3])/256 - 
         (6997*Log[b])/(2304*Sqrt[14]) - (6997*Log[x])/(1536*Sqrt[14])) + 
       (21631/(9216*Sqrt[14]) + (((1115*I)/1536)*Pi)/Sqrt[14] - 
         (1115*\[Gamma]E)/(768*Sqrt[14]) - (57499*Log[2])/(768*Sqrt[14]) + 
         (297*Log[3])/(32*Sqrt[14]) + (59375*Log[5])/(2304*Sqrt[14]) - 
         (1115*Log[b])/(768*Sqrt[14]) - (1115*Log[x])/(512*Sqrt[14]))/
        E^((3*I)*\[ScriptL]) + E^(I*\[ScriptL])*(-15229/(4608*Sqrt[14]) - 
         (((2171*I)/2304)*Pi)/Sqrt[14] + (785*\[Gamma]E)/(384*Sqrt[14]) - 
         (11405*Log[2])/(1152*Sqrt[14]) + (27*Log[6])/(4*Sqrt[14]) + 
         (785*Log[b])/(384*Sqrt[14]) + (785*Log[x])/(256*Sqrt[14])) + 
       (-46657/(13824*Sqrt[14]) - (((787*I)/768)*Pi)/Sqrt[14] + 
         (2405*\[Gamma]E)/(1152*Sqrt[14]) + (19001*Log[2])/(384*Sqrt[14]) - 
         (999*Log[6])/(64*Sqrt[14]) + (2405*Log[b])/(1152*Sqrt[14]) + 
         (2405*Log[x])/(768*Sqrt[14]))/E^(I*\[ScriptL]) + 
       E^((3*I)*\[ScriptL])*(-53447/(9216*Sqrt[14]) - (((493*I)/1536)*Pi)/
          Sqrt[14] + (2755*\[Gamma]E)/(768*Sqrt[14]) + 
         (14353*Log[2])/(2304*Sqrt[14]) + (27*Sqrt[7/2]*Log[3])/256 + 
         (2755*Log[b])/(768*Sqrt[14]) + (2755*Log[x])/(512*Sqrt[14])) + 
       (-515749/(15360*Sqrt[14]) - (((5317*I)/512)*Pi)/Sqrt[14] + 
         (5317*\[Gamma]E)/(256*Sqrt[14]) + (762217*Log[2])/(11520*Sqrt[14]) + 
         (15417*Log[3])/(320*Sqrt[14]) - (59375*Log[5])/(2304*Sqrt[14]) + 
         (5317*Log[b])/(256*Sqrt[14]) + (15951*Log[x])/(512*Sqrt[14]))/
        E^((5*I)*\[ScriptL])))
 
Htail[3, 2] = SO*x^4*\[Epsilon]^6*\[Chi]S*(((-80*I)/9)*Sqrt[5/7]*\[Nu] + 
     (8*Sqrt[5/7]*Pi*\[Nu])/3 + ((16*I)/3)*Sqrt[5/7]*\[Gamma]E*\[Nu] + 
     ((16*I)/9)*Sqrt[5/7]*\[Nu]*Log[64] + ((16*I)/3)*Sqrt[5/7]*\[Nu]*Log[b] + 
     (8*I)*Sqrt[5/7]*\[Nu]*Log[x] + e^2*(((-320*I)/9)*Sqrt[5/7]*\[Nu] + 
       (32*Sqrt[5/7]*Pi*\[Nu])/3 + ((64*I)/3)*Sqrt[5/7]*\[Gamma]E*\[Nu] - 
       ((40*I)/3)*Sqrt[5/7]*\[Nu]*Log[2] + (54*I)*Sqrt[5/7]*\[Nu]*Log[3] + 
       ((64*I)/3)*Sqrt[5/7]*\[Nu]*Log[b] + (32*I)*Sqrt[5/7]*\[Nu]*Log[x] + 
       E^((2*I)*\[ScriptL])*(((-230*I)/9)*Sqrt[5/7]*\[Nu] + 
         (23*Sqrt[5/7]*Pi*\[Nu])/3 + ((46*I)/3)*Sqrt[5/7]*\[Gamma]E*\[Nu] + 
         ((14*I)/3)*Sqrt[35]*\[Nu]*Log[2] + ((46*I)/3)*Sqrt[5/7]*\[Nu]*
          Log[b] + (23*I)*Sqrt[5/7]*\[Nu]*Log[x]) + 
       (((-530*I)/9)*Sqrt[5/7]*\[Nu] + (53*Sqrt[5/7]*Pi*\[Nu])/3 + 
         ((106*I)/3)*Sqrt[5/7]*\[Gamma]E*\[Nu] + (30*I)*Sqrt[35]*\[Nu]*
          Log[2] - (54*I)*Sqrt[5/7]*\[Nu]*Log[3] + ((106*I)/3)*Sqrt[5/7]*
          \[Nu]*Log[b] + (53*I)*Sqrt[5/7]*\[Nu]*Log[x])/
        E^((2*I)*\[ScriptL])) + 
     e^3*(E^((3*I)*\[ScriptL])*(((-305*I)/8)*Sqrt[5/7]*\[Nu] + 
         (1661*Sqrt[5/7]*Pi*\[Nu])/144 + ((183*I)/8)*Sqrt[5/7]*\[Gamma]E*
          \[Nu] + ((505*I)/72)*Sqrt[35]*\[Nu]*Log[2] + ((183*I)/8)*Sqrt[5/7]*
          \[Nu]*Log[b] + ((549*I)/16)*Sqrt[5/7]*\[Nu]*Log[x]) + 
       E^(I*\[ScriptL])*(((-395*I)/8)*Sqrt[5/7]*\[Nu] + 
         (237*Sqrt[5/7]*Pi*\[Nu])/16 + ((237*I)/8)*Sqrt[5/7]*\[Gamma]E*
          \[Nu] - ((793*I)/24)*Sqrt[5/7]*\[Nu]*Log[2] + ((351*I)/4)*Sqrt[5/7]*
          \[Nu]*Log[3] + ((237*I)/8)*Sqrt[5/7]*\[Nu]*Log[b] + 
         ((711*I)/16)*Sqrt[5/7]*\[Nu]*Log[x]) + 
       (((-4055*I)/72)*Sqrt[5/7]*\[Nu] + (811*Sqrt[5/7]*Pi*\[Nu])/48 + 
         ((811*I)/24)*Sqrt[5/7]*\[Gamma]E*\[Nu] + ((3289*I)/8)*Sqrt[5/7]*
          \[Nu]*Log[2] - ((1377*I)/8)*Sqrt[5/7]*\[Nu]*Log[3] + 
         ((811*I)/24)*Sqrt[5/7]*\[Nu]*Log[b] + ((811*I)/16)*Sqrt[5/7]*\[Nu]*
          Log[x])/E^(I*\[ScriptL]) + (((-8045*I)/72)*Sqrt[5/7]*\[Nu] + 
         (1609*Sqrt[5/7]*Pi*\[Nu])/48 + ((1609*I)/24)*Sqrt[5/7]*\[Gamma]E*
          \[Nu] - ((19717*I)/72)*Sqrt[5/7]*\[Nu]*Log[2] + 
         ((81*I)/4)*Sqrt[5/7]*\[Nu]*Log[3] + ((15625*I)/72)*Sqrt[5/7]*\[Nu]*
          Log[5] + ((1609*I)/24)*Sqrt[5/7]*\[Nu]*Log[b] + 
         ((1609*I)/16)*Sqrt[5/7]*\[Nu]*Log[x])/E^((3*I)*\[ScriptL])) + 
     e^4*(((-650*I)/9)*Sqrt[5/7]*\[Nu] + (65*Sqrt[5/7]*Pi*\[Nu])/3 + 
       ((130*I)/3)*Sqrt[5/7]*\[Gamma]E*\[Nu] + ((2032*I)/3)*Sqrt[5/7]*\[Nu]*
        Log[2] - ((621*I)/2)*Sqrt[5/7]*\[Nu]*Log[3] + ((130*I)/3)*Sqrt[5/7]*
        \[Nu]*Log[b] + (65*I)*Sqrt[5/7]*\[Nu]*Log[x] + 
       E^((4*I)*\[ScriptL])*(((-1480*I)/27)*Sqrt[5/7]*\[Nu] + 
         (611*Sqrt[5/7]*Pi*\[Nu])/36 + ((296*I)/9)*Sqrt[5/7]*\[Gamma]E*
          \[Nu] + ((638*I)/9)*Sqrt[5/7]*\[Nu]*Log[2] + ((296*I)/9)*Sqrt[5/7]*
          \[Nu]*Log[b] + ((148*I)/3)*Sqrt[5/7]*\[Nu]*Log[x]) + 
       E^((2*I)*\[ScriptL])*(((-1795*I)/27)*Sqrt[5/7]*\[Nu] + 
         (79*Sqrt[5/7]*Pi*\[Nu])/4 + ((359*I)/9)*Sqrt[5/7]*\[Gamma]E*\[Nu] - 
         (61*I)*Sqrt[5/7]*\[Nu]*Log[2] + ((531*I)/4)*Sqrt[5/7]*\[Nu]*Log[3] + 
         ((359*I)/9)*Sqrt[5/7]*\[Nu]*Log[b] + ((359*I)/6)*Sqrt[5/7]*\[Nu]*
          Log[x]) + (((-2045*I)/27)*Sqrt[5/7]*\[Nu] + 
         (409*Sqrt[5/7]*Pi*\[Nu])/18 + ((409*I)/9)*Sqrt[5/7]*\[Gamma]E*
          \[Nu] - ((9659*I)/9)*Sqrt[5/7]*\[Nu]*Log[2] + ((351*I)/2)*Sqrt[5/7]*
          \[Nu]*Log[3] + ((15625*I)/36)*Sqrt[5/7]*\[Nu]*Log[5] + 
         ((409*I)/9)*Sqrt[5/7]*\[Nu]*Log[b] + ((409*I)/6)*Sqrt[5/7]*\[Nu]*
          Log[x])/E^((2*I)*\[ScriptL]) + (((-5330*I)/27)*Sqrt[5/7]*\[Nu] + 
         (533*Sqrt[5/7]*Pi*\[Nu])/9 + ((1066*I)/9)*Sqrt[5/7]*\[Gamma]E*
          \[Nu] + ((314*I)/3)*Sqrt[35]*\[Nu]*Log[2] + ((279*I)/4)*Sqrt[35]*
          \[Nu]*Log[3] - ((15625*I)/36)*Sqrt[5/7]*\[Nu]*Log[5] + 
         ((1066*I)/9)*Sqrt[5/7]*\[Nu]*Log[b] + ((533*I)/3)*Sqrt[5/7]*\[Nu]*
          Log[x])/E^((4*I)*\[ScriptL])) + 
     e^6*(((-50*I)/3)*Sqrt[35]*\[Nu] + (4321*Sqrt[35]*Pi*\[Nu])/864 + 
       (10*I)*Sqrt[35]*\[Gamma]E*\[Nu] - ((12497*I)/27)*Sqrt[35]*\[Nu]*
        Log[2] + ((22329*I)/32)*Sqrt[5/7]*\[Nu]*Log[3] + 
       ((921875*I)/864)*Sqrt[5/7]*\[Nu]*Log[5] + (10*I)*Sqrt[35]*\[Nu]*
        Log[b] + (15*I)*Sqrt[35]*\[Nu]*Log[x] + E^((2*I)*\[ScriptL])*
        (((-23605*I)/216)*Sqrt[5/7]*\[Nu] + (19079*Sqrt[5/7]*Pi*\[Nu])/576 + 
         ((4721*I)/72)*Sqrt[5/7]*\[Gamma]E*\[Nu] + ((109649*I)/72)*Sqrt[5/7]*
          \[Nu]*Log[2] - ((49167*I)/64)*Sqrt[5/7]*\[Nu]*Log[3] + 
         ((4721*I)/72)*Sqrt[5/7]*\[Nu]*Log[b] + ((4721*I)/48)*Sqrt[5/7]*\[Nu]*
          Log[x]) + (((-2995*I)/24)*Sqrt[5/7]*\[Nu] + 
         (599*Sqrt[5/7]*Pi*\[Nu])/16 + ((599*I)/8)*Sqrt[5/7]*\[Gamma]E*
          \[Nu] + ((310403*I)/72)*Sqrt[5/7]*\[Nu]*Log[2] + 
         ((44019*I)/32)*Sqrt[5/7]*\[Nu]*Log[3] - ((1515625*I)/576)*Sqrt[5/7]*
          \[Nu]*Log[5] + ((599*I)/8)*Sqrt[5/7]*\[Nu]*Log[b] + 
         ((1797*I)/16)*Sqrt[5/7]*\[Nu]*Log[x])/E^((2*I)*\[ScriptL]) + 
       (((-1801*I)/27)*Sqrt[5/7]*\[Nu] + (1801*Pi*\[Nu])/(18*Sqrt[35]) + 
         (((1801*I)/9)*\[Gamma]E*\[Nu])/Sqrt[35] - 
         (((162563*I)/9)*\[Nu]*Log[2])/Sqrt[35] - 
         (((1074573*I)/64)*\[Nu]*Log[3])/Sqrt[35] + ((296875*I)/192)*
          Sqrt[5/7]*\[Nu]*Log[5] + ((823543*I)/576)*Sqrt[7/5]*\[Nu]*Log[7] + 
         (((1801*I)/9)*\[Nu]*Log[b])/Sqrt[35] + (((1801*I)/6)*\[Nu]*Log[x])/
          Sqrt[35])/E^((4*I)*\[ScriptL]) + E^((6*I)*\[ScriptL])*
        (((-7589*I)/72)*Sqrt[5/7]*\[Nu] + (310351*Pi*\[Nu])/(1728*Sqrt[35]) + 
         (((7589*I)/24)*\[Gamma]E*\[Nu])/Sqrt[35] + 
         (((49463*I)/72)*\[Nu]*Log[2])/Sqrt[35] - (((459*I)/64)*\[Nu]*Log[3])/
          Sqrt[35] + (((7589*I)/24)*\[Nu]*Log[b])/Sqrt[35] + 
         (((7589*I)/16)*\[Nu]*Log[x])/Sqrt[35]) + E^((4*I)*\[ScriptL])*
        (((-1027*I)/9)*Sqrt[5/7]*\[Nu] + (44659*Pi*\[Nu])/(288*Sqrt[35]) + 
         (((1027*I)/3)*\[Gamma]E*\[Nu])/Sqrt[35] - 
         (((7324*I)/9)*\[Nu]*Log[2])/Sqrt[35] + (((11187*I)/8)*\[Nu]*Log[3])/
          Sqrt[35] + (((1027*I)/3)*\[Nu]*Log[b])/Sqrt[35] + 
         (((1027*I)/2)*\[Nu]*Log[x])/Sqrt[35]) + 
       (((-39259*I)/72)*Sqrt[5/7]*\[Nu] + (39259*Pi*\[Nu])/(48*Sqrt[35]) + 
         (((39259*I)/24)*\[Gamma]E*\[Nu])/Sqrt[35] + ((1409251*I)/216)*
          Sqrt[5/7]*\[Nu]*Log[2] + (((116883*I)/64)*\[Nu]*Log[3])/Sqrt[35] + 
         ((15625*I)/864)*Sqrt[5/7]*\[Nu]*Log[5] - ((823543*I)/576)*Sqrt[7/5]*
          \[Nu]*Log[7] + (((39259*I)/24)*\[Nu]*Log[b])/Sqrt[35] + 
         (((39259*I)/16)*\[Nu]*Log[x])/Sqrt[35])/E^((6*I)*\[ScriptL])) + 
     e*(E^(I*\[ScriptL])*(((-145*I)/9)*Sqrt[5/7]*\[Nu] + 
         (29*Sqrt[5/7]*Pi*\[Nu])/6 + ((29*I)/3)*Sqrt[5/7]*\[Gamma]E*\[Nu] + 
         ((61*I)/3)*Sqrt[5/7]*\[Nu]*Log[2] + ((29*I)/3)*Sqrt[5/7]*\[Nu]*
          Log[b] + ((29*I)/2)*Sqrt[5/7]*\[Nu]*Log[x]) + 
       (((-35*I)/9)*Sqrt[35]*\[Nu] + (7*Sqrt[35]*Pi*\[Nu])/6 + 
         ((7*I)/3)*Sqrt[35]*\[Gamma]E*\[Nu] + ((17*I)/3)*Sqrt[5/7]*\[Nu]*
          Log[2] + (27*I)*Sqrt[5/7]*\[Nu]*Log[3] + ((7*I)/3)*Sqrt[35]*\[Nu]*
          Log[b] + ((7*I)/2)*Sqrt[35]*\[Nu]*Log[x])/E^(I*\[ScriptL])) + 
     e^5*(E^((5*I)*\[ScriptL])*(((-265385*I)/3456)*Sqrt[5/7]*\[Nu] + 
         (31663*Pi*\[Nu])/(256*Sqrt[35]) + ((53077*I)/1152)*Sqrt[5/7]*
          \[Gamma]E*\[Nu] + (((192707*I)/384)*\[Nu]*Log[2])/Sqrt[35] - 
         (((459*I)/128)*\[Nu]*Log[3])/Sqrt[35] + ((53077*I)/1152)*Sqrt[5/7]*
          \[Nu]*Log[b] + ((53077*I)/768)*Sqrt[5/7]*\[Nu]*Log[x]) + 
       E^((3*I)*\[ScriptL])*(((-33665*I)/384)*Sqrt[5/7]*\[Nu] + 
         (58259*Sqrt[5/7]*Pi*\[Nu])/2304 + ((6733*I)/128)*Sqrt[5/7]*\[Gamma]E*
          \[Nu] - ((117227*I)/1152)*Sqrt[5/7]*\[Nu]*Log[2] + 
         ((3105*I)/16)*Sqrt[5/7]*\[Nu]*Log[3] + ((6733*I)/128)*Sqrt[5/7]*
          \[Nu]*Log[b] + ((20199*I)/256)*Sqrt[5/7]*\[Nu]*Log[x]) + 
       E^(I*\[ScriptL])*(((-155345*I)/1728)*Sqrt[5/7]*\[Nu] + 
         (31153*Sqrt[5/7]*Pi*\[Nu])/1152 + ((31069*I)/576)*Sqrt[5/7]*
          \[Gamma]E*\[Nu] + ((594941*I)/576)*Sqrt[5/7]*\[Nu]*Log[2] - 
         ((15993*I)/32)*Sqrt[5/7]*\[Nu]*Log[3] + ((31069*I)/576)*Sqrt[5/7]*
          \[Nu]*Log[b] + ((31069*I)/384)*Sqrt[5/7]*\[Nu]*Log[x]) + 
       (((-13985*I)/1152)*Sqrt[35]*\[Nu] + (2797*Sqrt[35]*Pi*\[Nu])/768 + 
         ((2797*I)/384)*Sqrt[35]*\[Gamma]E*\[Nu] + ((2504465*I)/1152)*
          Sqrt[5/7]*\[Nu]*Log[2] + ((4239*I)/32)*Sqrt[35]*\[Nu]*Log[3] - 
         ((1671875*I)/1152)*Sqrt[5/7]*\[Nu]*Log[5] + ((2797*I)/384)*Sqrt[35]*
          \[Nu]*Log[b] + ((2797*I)/256)*Sqrt[35]*\[Nu]*Log[x])/
        E^((3*I)*\[ScriptL]) + (((-24035*I)/1728)*Sqrt[35]*\[Nu] + 
         (4807*Sqrt[35]*Pi*\[Nu])/1152 + ((4807*I)/576)*Sqrt[35]*\[Gamma]E*
          \[Nu] - ((127095*I)/64)*Sqrt[5/7]*\[Nu]*Log[2] + 
         ((25245*I)/64)*Sqrt[5/7]*\[Nu]*Log[3] + ((203125*I)/288)*Sqrt[5/7]*
          \[Nu]*Log[5] + ((4807*I)/576)*Sqrt[35]*\[Nu]*Log[b] + 
         ((4807*I)/384)*Sqrt[35]*\[Nu]*Log[x])/E^(I*\[ScriptL]) + 
       (((-164555*I)/3456)*Sqrt[35]*\[Nu] + (32911*Sqrt[35]*Pi*\[Nu])/2304 + 
         ((32911*I)/1152)*Sqrt[35]*\[Gamma]E*\[Nu] - 
         (((4364051*I)/1152)*\[Nu]*Log[2])/Sqrt[35] - ((7767*I)/8)*Sqrt[5/7]*
          \[Nu]*Log[3] + ((15625*I)/96)*Sqrt[5/7]*\[Nu]*Log[5] + 
         ((823543*I)/1152)*Sqrt[7/5]*\[Nu]*Log[7] + ((32911*I)/1152)*Sqrt[35]*
          \[Nu]*Log[b] + ((32911*I)/768)*Sqrt[35]*\[Nu]*Log[x])/
        E^((5*I)*\[ScriptL])))
 
Htail[3, 3] = x^3*\[Delta]*\[Epsilon]^4*((-291*Sqrt[3/70])/8 - 
     ((9*I)/4)*Sqrt[15/14]*Pi + (9*Sqrt[15/14]*\[Gamma]E)/2 + 
     (9*Sqrt[15/14]*Log[6*b])/2 + (27*Sqrt[15/14]*Log[x])/4 + 
     e*((-4559/(24*Sqrt[210]) - ((47*I)/4)*Sqrt[5/42]*Pi + 
         (47*Sqrt[5/42]*\[Gamma]E)/2 - (27*Sqrt[15/14]*Log[6])/2 + 
         16*Sqrt[10/21]*Log[64] + (47*Sqrt[5/42]*Log[b])/2 + 
         (47*Sqrt[15/14]*Log[x])/4)/E^(I*\[ScriptL]) + 
       E^(I*\[ScriptL])*(-1843/(8*Sqrt[210]) - ((19*I)/4)*Sqrt[15/14]*Pi + 
         (19*Sqrt[15/14]*\[Gamma]E)/2 + (11*Sqrt[15/14]*Log[2])/2 + 
         (27*Sqrt[15/14]*Log[3])/2 + (19*Sqrt[15/14]*Log[b])/2 + 
         (57*Sqrt[15/14]*Log[x])/4)) + e^2*(-1649/(8*Sqrt[210]) - 
       ((17*I)/4)*Sqrt[15/14]*Pi + (17*Sqrt[15/14]*\[Gamma]E)/2 + 
       108*Sqrt[30/7]*Log[2] - (135*Sqrt[15/14]*Log[6])/2 + 
       (17*Sqrt[15/14]*Log[b])/2 + (51*Sqrt[15/14]*Log[x])/4 + 
       ((-3007*Sqrt[5/42])/48 - ((155*I)/8)*Sqrt[5/42]*Pi + 
         (155*Sqrt[5/42]*\[Gamma]E)/4 - (6091*Sqrt[5/42]*Log[2])/16 + 
         (3125*Sqrt[5/42]*Log[5])/16 + (27*Sqrt[105/2]*Log[6])/16 + 
         (155*Sqrt[5/42]*Log[b])/4 + (155*Sqrt[15/14]*Log[x])/8)/
        E^((2*I)*\[ScriptL]) + E^((2*I)*\[ScriptL])*((-3977*Sqrt[5/42])/48 - 
         ((205*I)/8)*Sqrt[5/42]*Pi + (205*Sqrt[5/42]*\[Gamma]E)/4 + 
         (61*Sqrt[5/42]*Log[2])/4 + (459*Sqrt[15/14]*Log[3])/16 + 
         (205*Sqrt[5/42]*Log[b])/4 + (205*Sqrt[15/14]*Log[x])/8)) + 
     e^3*(E^(I*\[ScriptL])*(-43747/(192*Sqrt[210]) - ((451*I)/32)*Sqrt[5/42]*
          Pi + (451*Sqrt[5/42]*\[Gamma]E)/16 + (15587*Sqrt[5/42]*Log[2])/16 - 
         (675*Sqrt[15/14]*Log[3])/4 + (451*Sqrt[5/42]*Log[b])/16 + 
         (451*Sqrt[15/14]*Log[x])/32) + (-44717/(192*Sqrt[210]) - 
         ((461*I)/32)*Sqrt[5/42]*Pi + (461*Sqrt[5/42]*\[Gamma]E)/16 - 
         (1343*Sqrt[105/2]*Log[2])/16 + (135*Sqrt[105/2]*Log[3])/8 + 
         (3125*Sqrt[15/14]*Log[5])/16 + (461*Sqrt[5/42]*Log[b])/16 + 
         (461*Sqrt[15/14]*Log[x])/32)/E^(I*\[ScriptL]) + 
       ((-10573*Sqrt[3/70])/64 - ((327*I)/32)*Sqrt[15/14]*Pi + 
         (327*Sqrt[15/14]*\[Gamma]E)/16 + (673*Sqrt[105/2]*Log[2])/16 + 
         (639*Sqrt[15/14]*Log[3])/4 - (3125*Sqrt[15/14]*Log[5])/16 + 
         (327*Sqrt[15/14]*Log[b])/16 + (981*Sqrt[15/14]*Log[x])/32)/
        E^((3*I)*\[ScriptL]) + E^((3*I)*\[ScriptL])*(-44329/(64*Sqrt[210]) - 
         ((457*I)/32)*Sqrt[15/14]*Pi + (457*Sqrt[15/14]*\[Gamma]E)/16 + 
         (7*Sqrt[105/2]*Log[2])/16 + (423*Sqrt[15/14]*Log[3])/8 + 
         (457*Sqrt[15/14]*Log[b])/16 + (1371*Sqrt[15/14]*Log[x])/32)) + 
     e^4*((-5723*Sqrt[3/70])/64 - ((177*I)/32)*Sqrt[15/14]*Pi + 
       (177*Sqrt[15/14]*\[Gamma]E)/16 - (23743*Sqrt[15/14]*Log[2])/16 + 
       (23679*Sqrt[15/14]*Log[3])/64 + (53125*Sqrt[15/14]*Log[5])/128 + 
       (177*Sqrt[15/14]*Log[b])/16 + (531*Sqrt[15/14]*Log[x])/32 + 
       (-879887/(1152*Sqrt[210]) - ((9071*I)/192)*Sqrt[5/42]*Pi + 
         (9071*Sqrt[5/42]*\[Gamma]E)/96 - (137041*Sqrt[5/42]*Log[2])/96 - 
         (124461*Sqrt[15/14]*Log[3])/256 + (3125*Sqrt[105/2]*Log[5])/128 + 
         (117649*Sqrt[35/6]*Log[7])/768 + (9071*Sqrt[5/42]*Log[b])/96 + 
         (9071*Sqrt[5/42]*Log[x])/64)/E^((4*I)*\[ScriptL]) + 
       (-1843/(8*Sqrt[210]) - ((19*I)/4)*Sqrt[15/14]*Pi + 
         (19*Sqrt[15/14]*\[Gamma]E)/2 + (2915*Sqrt[15/14]*Log[2])/2 + 
         (13005*Sqrt[15/14]*Log[3])/32 - (28125*Sqrt[15/14]*Log[5])/32 + 
         (19*Sqrt[15/14]*Log[b])/2 + (57*Sqrt[15/14]*Log[x])/4)/
        E^((2*I)*\[ScriptL]) + E^((4*I)*\[ScriptL])*
        (-426509/(384*Sqrt[210]) - ((2509*I)/256)*Sqrt[35/6]*Pi + 
         (4397*Sqrt[5/42]*\[Gamma]E)/32 - (115*Sqrt[5/42]*Log[2])/32 + 
         (23103*Sqrt[15/14]*Log[3])/256 + (4397*Sqrt[5/42]*Log[b])/32 + 
         (4397*Sqrt[15/14]*Log[x])/64) + E^((2*I)*\[ScriptL])*
        ((-3977*Sqrt[7/30])/144 - ((41*I)/24)*Sqrt[35/6]*Pi + 
         (41*Sqrt[35/6]*\[Gamma]E)/12 + (21935*Sqrt[5/42]*Log[2])/12 - 
         (10845*Sqrt[15/14]*Log[3])/32 + (41*Sqrt[35/6]*Log[b])/12 + 
         (41*Sqrt[35/6]*Log[x])/8)) + 
     e^5*(((-270727*Sqrt[5/42])/4608 - ((13955*I)/768)*Sqrt[5/42]*Pi + 
         (13955*Sqrt[5/42]*\[Gamma]E)/384 + (1543201*Sqrt[5/42]*Log[2])/128 + 
         (84519*Sqrt[15/14]*Log[3])/128 - (34375*Sqrt[15/14]*Log[5])/16 + 
         (13955*Sqrt[5/42]*Log[b])/384 + (13955*Sqrt[5/42]*Log[x])/256)/
        E^(I*\[ScriptL]) + E^(I*\[ScriptL])*((-271697*Sqrt[5/42])/4608 - 
         ((14005*I)/768)*Sqrt[5/42]*Pi + (14005*Sqrt[5/42]*\[Gamma]E)/384 - 
         (3463843*Sqrt[5/42]*Log[2])/384 + (105669*Sqrt[15/14]*Log[3])/128 + 
         (146875*Sqrt[5/42]*Log[5])/64 + (14005*Sqrt[5/42]*Log[b])/384 + 
         (14005*Sqrt[5/42]*Log[x])/256) + (-10652443/(9216*Sqrt[210]) - 
         ((109819*I)/1536)*Sqrt[5/42]*Pi + (109819*Sqrt[5/42]*\[Gamma]E)/
          768 + (4370657*Sqrt[7/30]*Log[2])/768 + (68031*Sqrt[3/70]*Log[3])/
          32 - (3125*Sqrt[5/42]*Log[5])/32 - (117649*Sqrt[35/6]*Log[7])/256 + 
         (109819*Sqrt[5/42]*Log[b])/768 + (109819*Sqrt[5/42]*Log[x])/512)/
        E^((5*I)*\[ScriptL]) + E^((5*I)*\[ScriptL])*
        (-15972893/(9216*Sqrt[210]) - (((819431*I)/1536)*Pi)/Sqrt[210] + 
         (164669*Sqrt[5/42]*\[Gamma]E)/768 - (99943*Log[2])/(768*Sqrt[210]) + 
         (3357*Sqrt[21/10]*Log[3])/32 + (164669*Sqrt[5/42]*Log[b])/768 + 
         (164669*Sqrt[5/42]*Log[x])/512) + ((-33853*Sqrt[5/42])/1024 - 
         ((1745*I)/512)*Sqrt[15/14]*Pi + (1745*Sqrt[15/14]*\[Gamma]E)/256 - 
         (2146325*Sqrt[5/42]*Log[2])/256 - (69579*Sqrt[15/14]*Log[3])/32 + 
         (90625*Sqrt[15/14]*Log[5])/64 + (117649*Sqrt[35/6]*Log[7])/256 + 
         (1745*Sqrt[15/14]*Log[b])/256 + (5235*Sqrt[15/14]*Log[x])/512)/
        E^((3*I)*\[ScriptL]) + E^((3*I)*\[ScriptL])*
        ((-1261*Sqrt[35/6])/1024 - ((505*I)/512)*Sqrt[15/14]*Pi + 
         (65*Sqrt[105/2]*\[Gamma]E)/256 + (808277*Sqrt[5/42]*Log[2])/256 - 
         (9801*Sqrt[15/14]*Log[3])/16 + (65*Sqrt[105/2]*Log[b])/256 + 
         (195*Sqrt[105/2]*Log[x])/512)) + 
     e^6*((-97*Sqrt[105/2])/32 - ((15*I)/16)*Sqrt[105/2]*Pi + 
       (15*Sqrt[105/2]*\[Gamma]E)/8 + (205307*Sqrt[5/42]*Log[2])/8 + 
       (531*Sqrt[105/2]*Log[3])/4 - (821875*Sqrt[5/42]*Log[5])/64 + 
       (15*Sqrt[105/2]*Log[b])/8 + (45*Sqrt[105/2]*Log[x])/16 + 
       ((-92053*Sqrt[3/70])/160 - ((2847*I)/16)*Sqrt[3/70]*Pi + 
         (2847*Sqrt[3/70]*\[Gamma]E)/8 - (22505*Sqrt[35/6]*Log[2])/8 + 
         (27868131*Sqrt[3/70]*Log[3])/2048 - (15625*Sqrt[5/42]*Log[5])/2048 + 
         (823543*Sqrt[35/6]*Log[7])/2048 + (2847*Sqrt[3/70]*Log[b])/8 + 
         (8541*Sqrt[3/70]*Log[x])/16)/E^((6*I)*\[ScriptL]) + 
       E^((6*I)*\[ScriptL])*(-1697791/(640*Sqrt[210]) - 
         (((831779*I)/1024)*Pi)/Sqrt[210] + (17503*Sqrt[3/70]*\[Gamma]E)/32 - 
         (10051*Log[2])/(32*Sqrt[210]) + (37053*Sqrt[3/70]*Log[3])/32 + 
         (17503*Sqrt[3/70]*Log[b])/32 + (52509*Sqrt[3/70]*Log[x])/64) + 
       (-184979/(576*Sqrt[210]) - ((1907*I)/96)*Sqrt[5/42]*Pi + 
         (1907*Sqrt[5/42]*\[Gamma]E)/48 - (1178449*Sqrt[5/42]*Log[2])/48 - 
         (10619001*Sqrt[15/14]*Log[3])/2048 + (80115625*Sqrt[5/42]*Log[5])/
          6144 + (2000033*Sqrt[35/6]*Log[7])/2048 + (1907*Sqrt[5/42]*Log[b])/
          48 + (1907*Sqrt[5/42]*Log[x])/32)/E^((2*I)*\[ScriptL]) + 
       E^((4*I)*\[ScriptL])*((132211*Sqrt[7/30])/2880 + (((17659*I)/192)*Pi)/
          Sqrt[210] - (1363*Sqrt[7/30]*\[Gamma]E)/48 + (1243967*Log[2])/
          (48*Sqrt[210]) - (166653*Sqrt[3/70]*Log[3])/32 - 
         (1363*Sqrt[7/30]*Log[b])/48 - (1363*Sqrt[7/30]*Log[x])/32) + 
       E^((2*I)*\[ScriptL])*((-54029*Sqrt[7/30])/1152 - 
         ((17599*I)/6144)*Sqrt[35/6]*Pi + (557*Sqrt[35/6]*\[Gamma]E)/96 - 
         (1577341*Sqrt[5/42]*Log[2])/96 + (3272733*Sqrt[15/14]*Log[3])/2048 + 
         (8021875*Sqrt[5/42]*Log[5])/2048 + (557*Sqrt[35/6]*Log[b])/96 + 
         (557*Sqrt[35/6]*Log[x])/64) + (50537/(2880*Sqrt[210]) + 
         (((521*I)/96)*Pi)/Sqrt[210] - (521*\[Gamma]E)/(48*Sqrt[210]) + 
         (7322227*Log[2])/(48*Sqrt[210]) + 8811*Sqrt[6/35]*Log[3] - 
         (171875*Sqrt[5/42]*Log[5])/64 - (2000033*Sqrt[7/30]*Log[7])/192 - 
         (521*Log[b])/(48*Sqrt[210]) - (521*Log[x])/(32*Sqrt[210]))/
        E^((4*I)*\[ScriptL])))


(* ::Section::Closed:: *)
(*Oscillatory memory contributions*)


(* ::Text:: *)
(*The oscillatory memory is expanded in eccentricity to O(e^6), and expressed in terms of (x,e,\[ScriptL]).*)
(*In addition to the 3PN spin contributions, we also include the full 3PN nonspinning part.*)


HoscMem[2, 0] = x^(7/2)*\[Epsilon]^5*
    ((((16*I)/7)*Sqrt[6]*e*(-1 + E^((2*I)*\[ScriptL]))*\[Nu])/
      E^(I*\[ScriptL]) + (((647*I)/42)*e^2*(-1 + E^((4*I)*\[ScriptL]))*\[Nu])/
      (Sqrt[6]*E^((2*I)*\[ScriptL])) + (((80*I)/63)*Sqrt[2/3]*e^3*
       (-8 - 21*E^((2*I)*\[ScriptL]) + 21*E^((4*I)*\[ScriptL]) + 
        8*E^((6*I)*\[ScriptL]))*\[Nu])/E^((3*I)*\[ScriptL]) + 
     ((I/336)*e^4*(-9413 - 14264*E^((2*I)*\[ScriptL]) + 
        14264*E^((6*I)*\[ScriptL]) + 9413*E^((8*I)*\[ScriptL]))*\[Nu])/
      (Sqrt[6]*E^((4*I)*\[ScriptL])) + 
     ((I/1260)*e^5*(-49471 - 51925*E^((2*I)*\[ScriptL]) - 
        161350*E^((4*I)*\[ScriptL]) + 161350*E^((6*I)*\[ScriptL]) + 
        51925*E^((8*I)*\[ScriptL]) + 49471*E^((10*I)*\[ScriptL]))*\[Nu])/
      (Sqrt[6]*E^((5*I)*\[ScriptL])) + 
     ((I/5040)*e^6*(-279128 - 207983*E^((2*I)*\[ScriptL]) - 
        454160*E^((4*I)*\[ScriptL]) + 454160*E^((8*I)*\[ScriptL]) + 
        207983*E^((10*I)*\[ScriptL]) + 279128*E^((12*I)*\[ScriptL]))*\[Nu])/
      (Sqrt[6]*E^((6*I)*\[ScriptL])))
 
HoscMem[2, 2] = x^(5/2)*\[Epsilon]^3*
     (((-13*I)/252)*e^2*E^((2*I)*\[ScriptL])*\[Nu] - 
      ((13*I)/126)*e^3*E^(I*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Nu] - 
      (I/1008)*e^4*(39 - 70*E^((2*I)*\[ScriptL]) + 169*E^((4*I)*\[ScriptL]))*
       \[Nu] - ((I/3024)*e^5*(13 - 555*E^((2*I)*\[ScriptL]) - 
         225*E^((4*I)*\[ScriptL]) + 767*E^((6*I)*\[ScriptL]))*\[Nu])/
       E^(I*\[ScriptL]) - ((I/12096)*e^6*(26 + 1320*E^((2*I)*\[ScriptL]) - 
         1509*E^((4*I)*\[ScriptL]) - 1352*E^((6*I)*\[ScriptL]) + 
         4485*E^((8*I)*\[ScriptL]))*\[Nu])/E^((2*I)*\[ScriptL])) + 
    SO*x^3*\[Epsilon]^4*(((13*I)/378)*e^2*E^((2*I)*\[ScriptL])*\[Nu]*
       (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S) + 
      ((13*I)/189)*e^3*E^(I*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
       (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S) + 
      (I/1512)*e^4*(39 - 44*E^((2*I)*\[ScriptL]) + 169*E^((4*I)*\[ScriptL]))*
       \[Nu]*(-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S) + 
      ((I/4536)*e^5*(13 - 711*E^((2*I)*\[ScriptL]) - 
         69*E^((4*I)*\[ScriptL]) + 767*E^((6*I)*\[ScriptL]))*\[Nu]*
        (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/E^(I*\[ScriptL]) + 
      ((I/18144)*e^6*(26 + 1554*E^((2*I)*\[ScriptL]) - 
         1695*E^((4*I)*\[ScriptL]) - 338*E^((6*I)*\[ScriptL]) + 
         4485*E^((8*I)*\[ScriptL]))*\[Nu]*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu])*\[Chi]S))/E^((2*I)*\[ScriptL])) + 
    \[Epsilon]^5*(x^(7/2)*((((-4*I)/63)*e*(-1 + 3*E^((2*I)*\[ScriptL]))*
          \[Nu])/E^(I*\[ScriptL]) - ((I/1512)*e^3*\[Nu]*
          (-342 - 162*E^((2*I)*\[ScriptL]) + E^((4*I)*\[ScriptL])*
            (7449 - 12032*\[Nu]) + E^((6*I)*\[ScriptL])*
            (-6465 + 12032*\[Nu])))/E^((3*I)*\[ScriptL]) - 
        ((I/3024)*e^2*\[Nu]*(-390 - 456*E^((2*I)*\[ScriptL]) + 
           E^((4*I)*\[ScriptL])*(-6273 + 12110*\[Nu])))/
         E^((2*I)*\[ScriptL]) - ((I/12096)*e^4*\[Nu]*
          (-4498 - 1476*E^((2*I)*\[ScriptL]) + 3*E^((4*I)*\[ScriptL])*
            (-8623 + 11902*\[Nu]) + E^((6*I)*\[ScriptL])*
            (-9670 + 21724*\[Nu]) + E^((8*I)*\[ScriptL])*
            (-84705 + 155558*\[Nu])))/E^((4*I)*\[ScriptL]) - 
        ((I/36288)*e^5*\[Nu]*(-21372 - 4212*E^((2*I)*\[ScriptL]) + 
           E^((6*I)*\[ScriptL])*(614925 - 1023564*\[Nu]) + 
           93*E^((8*I)*\[ScriptL])*(-1825 + 3324*\[Nu]) + 
           E^((4*I)*\[ScriptL])*(-13635 + 11876*\[Nu]) + 
           E^((10*I)*\[ScriptL])*(-385821 + 702556*\[Nu])))/
         E^((5*I)*\[ScriptL]) - ((I/725760)*e^6*\[Nu]*
          (-660582 - 43776*E^((2*I)*\[ScriptL]) + 20*E^((4*I)*\[ScriptL])*
            (-10497 + 5899*\[Nu]) - 60*E^((8*I)*\[ScriptL])*
            (-87781 + 146170*\[Nu]) + 60*E^((6*I)*\[ScriptL])*
            (-105495 + 163828*\[Nu]) + 4*E^((10*I)*\[ScriptL])*
            (-1466547 + 2634340*\[Nu]) + 3*E^((12*I)*\[ScriptL])*
            (-3769863 + 6816610*\[Nu])))/E^((6*I)*\[ScriptL])) + 
      SO^2*x^(7/2)*(((13*I)/4536)*e^2*E^((2*I)*\[ScriptL])*\[Nu]*
         (9*\[Delta]*\[Kappa]A + \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 
          92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 
          23*\[Chi]S^2 + 32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2) + 
        (I/18144)*e^4*(39 - 18*E^((2*I)*\[ScriptL]) + 
          169*E^((4*I)*\[ScriptL]))*\[Nu]*(9*\[Delta]*\[Kappa]A + 
          \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 92*\[Nu]*\[Chi]A^2 + 
          2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 23*\[Chi]S^2 + 
          32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2) + 
        ((I/54432)*e^5*(13 - 867*E^((2*I)*\[ScriptL]) + 
           87*E^((4*I)*\[ScriptL]) + 767*E^((6*I)*\[ScriptL]))*\[Nu]*
          (9*\[Delta]*\[Kappa]A + \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 
           92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 
           23*\[Chi]S^2 + 32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2))/
         E^(I*\[ScriptL]) + ((I/217728)*e^6*(26 + 1788*E^((2*I)*\[ScriptL]) - 
           1725*E^((4*I)*\[ScriptL]) + 676*E^((6*I)*\[ScriptL]) + 
           4485*E^((8*I)*\[ScriptL]))*\[Nu]*(9*\[Delta]*\[Kappa]A + 
           \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 92*\[Nu]*\[Chi]A^2 + 
           2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 23*\[Chi]S^2 + 
           32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2))/E^((2*I)*\[ScriptL]) - 
        ((13*I)/2268)*e^3*E^(I*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
         (-9*\[Delta]*\[Kappa]A + 9*\[Kappa]S*(-1 + 2*\[Nu]) + 23*\[Chi]A^2 - 
          92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(23 - 16*\[Nu])*\[Chi]A*\[Chi]S + 
          23*\[Chi]S^2 - 32*\[Nu]*\[Chi]S^2 + 8*\[Nu]^2*\[Chi]S^2))) + 
    \[Epsilon]^6*(x^4*(((-29*I)/126)*e^2*E^((2*I)*\[ScriptL])*Pi*\[Nu] - 
        ((29*I)/63)*e^3*E^(I*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*Pi*
         \[Nu] - (I/1512)*e^4*(261 + 473*E^((2*I)*\[ScriptL]) + 
          1131*E^((4*I)*\[ScriptL]))*Pi*\[Nu] - 
        ((I/1512)*e^5*(29 - 3121*E^((2*I)*\[ScriptL]) + 
           1381*E^((4*I)*\[ScriptL]) + 1711*E^((6*I)*\[ScriptL]))*Pi*\[Nu])/
         E^(I*\[ScriptL]) - ((I/24192)*e^6*(232 + 
           23076*E^((2*I)*\[ScriptL]) - 13897*E^((4*I)*\[ScriptL]) + 
           36892*E^((6*I)*\[ScriptL]) + 40020*E^((8*I)*\[ScriptL]))*Pi*\[Nu])/
         E^((2*I)*\[ScriptL])) + SO*x^4*((I/3024)*e^3*E^(I*\[ScriptL])*
         (-1 + E^((2*I)*\[ScriptL]))*\[Nu]*(\[Delta]*(20017 - 32172*\[Nu])*
           \[Chi]A + (20017 - 58604*\[Nu] + 16320*\[Nu]^2)*\[Chi]S) + 
        (I/6048)*e^2*E^((2*I)*\[ScriptL])*\[Nu]*
         (\[Delta]*(20537 - 32380*\[Nu])*\[Chi]A + 
          (20537 - 59228*\[Nu] + 16424*\[Nu]^2)*\[Chi]S) + 
        (I/72576)*e^4*\[Nu]*(-(\[Delta]*(39*E^((4*I)*\[ScriptL])*
              (-19577 + 31996*\[Nu]) + 4*E^((2*I)*\[ScriptL])*
              (-71809 + 90264*\[Nu]) + 3*(-57451 + 95476*\[Nu]))*\[Chi]A) + 
          (39*E^((4*I)*\[ScriptL])*(19577 - 58076*\[Nu] + 16232*\[Nu]^2) + 
            4*E^((2*I)*\[ScriptL])*(71809 - 173864*\[Nu] + 45240*\[Nu]^2) + 
            3*(57451 - 172692*\[Nu] + 48440*\[Nu]^2))*\[Chi]S) + 
        ((I/72576)*e^5*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
          (-(\[Delta]*(18977 - 31756*\[Nu] + 10*E^((2*I)*\[ScriptL])*(
                -198933 + 307628*\[Nu]) + E^((4*I)*\[ScriptL])*(-1132123 + 
                1878596*\[Nu]))*\[Chi]A) + (-18977 + 57356*\[Nu] - 
             16112*\[Nu]^2 + 2*E^((2*I)*\[ScriptL])*(994665 - 2823964*\[Nu] + 
               778160*\[Nu]^2) + E^((4*I)*\[ScriptL])*(1132123 - 3398980*
                \[Nu] + 953104*\[Nu]^2))*\[Chi]S))/E^(I*\[ScriptL]) - 
        ((I/290304)*e^6*(\[Delta]*\[Nu]*(-36914 + E^((4*I)*\[ScriptL])*
              (2109756 - 3763332*\[Nu]) + 63096*\[Nu] + 
             18*E^((2*I)*\[ScriptL])*(-200561 + 322628*\[Nu]) + 
             3*E^((8*I)*\[ScriptL])*(-2166235 + 3645492*\[Nu]) + 
             E^((6*I)*\[ScriptL])*(-5679950 + 8033592*\[Nu]))*\[Chi]A - 
           \[Nu]*(36914 - 113464*\[Nu] + 32016*\[Nu]^2 - 
             12*E^((4*I)*\[ScriptL])*(175813 - 559523*\[Nu] + 159204*
                \[Nu]^2) + 18*E^((2*I)*\[ScriptL])*(200561 - 587332*\[Nu] + 
               163336*\[Nu]^2) + 3*E^((8*I)*\[ScriptL])*(2166235 - 6576596*
                \[Nu] + 1849656*\[Nu]^2) + 2*E^((6*I)*\[ScriptL])*
              (2839975 - 7520924*\[Nu] + 2023608*\[Nu]^2))*\[Chi]S))/
         E^((2*I)*\[ScriptL])) + SO^3*x^4*
       (((-13*I)/3402)*e^2*E^((2*I)*\[ScriptL])*\[Nu]*
         (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
             \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
            (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
             \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
            3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
            (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)) - 
        ((13*I)/1701)*e^3*E^(I*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
         (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
             \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
            (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
             \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
            3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
            (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)) - 
        (I/13608)*e^4*(39 + 8*E^((2*I)*\[ScriptL]) + 
          169*E^((4*I)*\[ScriptL]))*\[Nu]*
         (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
             \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
            (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
             \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
            3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
            (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)) - 
        ((I/40824)*e^5*(13 - 1023*E^((2*I)*\[ScriptL]) + 
           243*E^((4*I)*\[ScriptL]) + 767*E^((6*I)*\[ScriptL]))*\[Nu]*
          (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
              \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
              \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/
         E^(I*\[ScriptL]) - ((I/163296)*e^6*(26 + 2022*E^((2*I)*\[ScriptL]) - 
           1599*E^((4*I)*\[ScriptL]) + 1690*E^((6*I)*\[ScriptL]) + 
           4485*E^((8*I)*\[ScriptL]))*\[Nu]*
          (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
              \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
              \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/
         E^((2*I)*\[ScriptL])))
 
HoscMem[3, 1] = x^3*\[Epsilon]^4*
     ((-22*Sqrt[14]*e*E^(I*\[ScriptL])*\[Delta]*\[Nu])/135 - 
      (22*Sqrt[14]*e^2*(-1 + E^((2*I)*\[ScriptL]))*\[Delta]*\[Nu])/135 - 
      (e^3*(-308 + 8735*E^((2*I)*\[ScriptL]) + 2772*E^((4*I)*\[ScriptL]))*
        \[Delta]*\[Nu])/(1080*Sqrt[14]*E^(I*\[ScriptL])) - 
      (e^4*(-616 - 33597*E^((2*I)*\[ScriptL]) + 24357*E^((4*I)*\[ScriptL]) + 
         9856*E^((6*I)*\[ScriptL]))*\[Delta]*\[Nu])/
       (3240*Sqrt[14]*E^((2*I)*\[ScriptL])) - 
      (e^5*(-4158 - 32365*E^((2*I)*\[ScriptL]) + 
         374448*E^((4*I)*\[ScriptL]) + 202581*E^((6*I)*\[ScriptL]) + 
         96250*E^((8*I)*\[ScriptL]))*\[Delta]*\[Nu])/
       (25920*Sqrt[14]*E^((3*I)*\[ScriptL])) - 
      (e^6*(-9856 - 52915*E^((2*I)*\[ScriptL]) - 
         1591890*E^((4*I)*\[ScriptL]) + 804325*E^((6*I)*\[ScriptL]) + 
         550960*E^((8*I)*\[ScriptL]) + 299376*E^((10*I)*\[ScriptL]))*\[Delta]*
        \[Nu])/(64800*Sqrt[14]*E^((4*I)*\[ScriptL]))) + 
    SO*x^(7/2)*\[Epsilon]^5*((e*E^(I*\[ScriptL])*\[Nu]*
        ((-5729 + 19712*\[Nu])*\[Chi]A + \[Delta]*(-5729 + 2464*\[Nu])*
          \[Chi]S))/(1620*Sqrt[14]) + (e^2*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
        ((-5729 + 19712*\[Nu])*\[Chi]A + \[Delta]*(-5729 + 2464*\[Nu])*
          \[Chi]S))/(1620*Sqrt[14]) + 
      (e^3*\[Nu]*((5729 - 19712*\[Nu] + 9*E^((4*I)*\[ScriptL])*
            (-5729 + 19712*\[Nu]) + 4*E^((2*I)*\[ScriptL])*
            (-46249 + 159472*\[Nu]))*\[Chi]A + \[Delta]*(5729 - 2464*\[Nu] + 
           9*E^((4*I)*\[ScriptL])*(-5729 + 2464*\[Nu]) + 
           4*E^((2*I)*\[ScriptL])*(-46249 + 19934*\[Nu]))*\[Chi]S))/
       (12960*Sqrt[14]*E^(I*\[ScriptL])) + 
      (e^5*\[Nu]*((154683 - 532224*\[Nu] + 8640*E^((6*I)*\[ScriptL])*
            (-1013 + 3494*\[Nu]) + 625*E^((8*I)*\[ScriptL])*
            (-5729 + 19712*\[Nu]) - 128*E^((2*I)*\[ScriptL])*
            (-10462 + 36061*\[Nu]) + 54*E^((4*I)*\[ScriptL])*
            (-344511 + 1188608*\[Nu]))*\[Chi]A + 
         \[Delta]*(154683 - 66528*\[Nu] + 2160*E^((6*I)*\[ScriptL])*
            (-4052 + 1747*\[Nu]) + 625*E^((8*I)*\[ScriptL])*
            (-5729 + 2464*\[Nu]) - 16*E^((2*I)*\[ScriptL])*
            (-83696 + 36061*\[Nu]) + 54*E^((4*I)*\[ScriptL])*
            (-344511 + 148576*\[Nu]))*\[Chi]S))/(622080*Sqrt[14]*
        E^((3*I)*\[ScriptL])) + (e^4*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
        ((-5729 + 19712*\[Nu] + 16*E^((4*I)*\[ScriptL])*
            (-5729 + 19712*\[Nu]) + E^((2*I)*\[ScriptL])*
            (-351971 + 1213088*\[Nu]))*\[Chi]A + 
         \[Delta]*(-5729 + 2464*\[Nu] + 16*E^((4*I)*\[ScriptL])*
            (-5729 + 2464*\[Nu]) + E^((2*I)*\[ScriptL])*(-351971 + 
             151636*\[Nu]))*\[Chi]S))/(19440*Sqrt[14]*E^((2*I)*\[ScriptL])) + 
      (e^6*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
        ((-91664 + 315392*\[Nu] + 486*E^((8*I)*\[ScriptL])*
            (-5729 + 19712*\[Nu]) + 27*E^((2*I)*\[ScriptL])*
            (-23707 + 81696*\[Nu]) + 3*E^((4*I)*\[ScriptL])*
            (-6347263 + 21890464*\[Nu]) + E^((6*I)*\[ScriptL])*
            (-8809174 + 30364672*\[Nu]))*\[Chi]A + 
         \[Delta]*(-91664 + 39424*\[Nu] + 486*E^((8*I)*\[ScriptL])*
            (-5729 + 2464*\[Nu]) + 27*E^((2*I)*\[ScriptL])*
            (-23707 + 10212*\[Nu]) + 3*E^((4*I)*\[ScriptL])*
            (-6347263 + 2736308*\[Nu]) + E^((6*I)*\[ScriptL])*
            (-8809174 + 3795584*\[Nu]))*\[Chi]S))/(388800*Sqrt[14]*
        E^((4*I)*\[ScriptL]))) + \[Epsilon]^6*
     (x^4*((-121*\[Delta]*\[Nu])/(45*Sqrt[14]) - 
        (e^2*\[Delta]*\[Nu]*(65714 + E^((2*I)*\[ScriptL])*
            (335463 - 6616*\[Nu]) + E^((4*I)*\[ScriptL])*
            (-167097 + 6616*\[Nu])))/(11880*Sqrt[14]*E^((2*I)*\[ScriptL])) - 
        (e*\[Delta]*\[Nu]*(39732 + E^((2*I)*\[ScriptL])*(-138607 + 
             20168*\[Nu])))/(11880*Sqrt[14]*E^(I*\[ScriptL])) + 
        (e^3*\[Delta]*\[Nu]*(-848870 + E^((2*I)*\[ScriptL])*
            (-1158395 + 20168*\[Nu]) + E^((6*I)*\[ScriptL])*
            (1870811 + 35320*\[Nu]) + E^((4*I)*\[ScriptL])*
            (2343930 + 563612*\[Nu])))/(95040*Sqrt[14]*
          E^((3*I)*\[ScriptL])) + (e^4*\[Delta]*\[Nu]*(-2010855 + 
           4*E^((2*I)*\[ScriptL])*(-497155 + 5042*\[Nu]) + 
           4*E^((8*I)*\[ScriptL])*(1031179 + 51460*\[Nu]) + 
           6*E^((6*I)*\[ScriptL])*(420121 + 247072*\[Nu]) - 
           3*E^((4*I)*\[ScriptL])*(3983461 + 569480*\[Nu])))/
         (142560*Sqrt[14]*E^((4*I)*\[ScriptL])) + 
        (e^6*\[Delta]*\[Nu]*(-95514177 - 10*E^((4*I)*\[ScriptL])*
            (7780026 + 151307*\[Nu]) + E^((2*I)*\[ScriptL])*
            (-52452347 + 322688*\[Nu]) + 3*E^((12*I)*\[ScriptL])*
            (61570609 + 4525184*\[Nu]) - 30*E^((6*I)*\[ScriptL])*
            (16558956 + 5511281*\[Nu]) + 10*E^((8*I)*\[ScriptL])*
            (8495720 + 9338059*\[Nu]) + E^((10*I)*\[ScriptL])*
            (14922117 + 59572670*\[Nu])))/(2851200*Sqrt[14]*
          E^((6*I)*\[ScriptL])) + (e^5*\[Delta]*\[Nu]*(-99854216 - 
           8*E^((4*I)*\[ScriptL])*(15366264 + 430285*\[Nu]) + 
           E^((2*I)*\[ScriptL])*(-74658597 + 544536*\[Nu]) + 
           24*E^((8*I)*\[ScriptL])*(2416139 + 2878583*\[Nu]) + 
           E^((10*I)*\[ScriptL])*(197166351 + 12981176*\[Nu]) + 
           2*E^((6*I)*\[ScriptL])*(106584163 + 61386984*\[Nu])))/
         (4561920*Sqrt[14]*E^((5*I)*\[ScriptL]))) + 
      SO^2*x^4*(-(e*E^(I*\[ScriptL])*\[Nu]*(2772*\[Kappa]A*(-1 + 4*\[Nu]) + 
            (17372 - 73737*\[Nu] + 39424*\[Nu]^2)*\[Chi]A*\[Chi]S + 
            \[Delta]*(2772*\[Kappa]S*(-1 + 2*\[Nu]) + (8686 - 28336*\[Nu])*
               \[Chi]A^2 + (8686 - 10657*\[Nu] + 2464*\[Nu]^2)*\[Chi]S^2)))/
         (2430*Sqrt[14]) - (e^2*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
          (2772*\[Kappa]A*(-1 + 4*\[Nu]) + (17372 - 73737*\[Nu] + 
             39424*\[Nu]^2)*\[Chi]A*\[Chi]S + \[Delta]*
            (2772*\[Kappa]S*(-1 + 2*\[Nu]) + (8686 - 28336*\[Nu])*\[Chi]A^2 + 
             (8686 - 10657*\[Nu] + 2464*\[Nu]^2)*\[Chi]S^2)))/
         (2430*Sqrt[14]) - (e^4*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
          (9*(616 + 41605*E^((2*I)*\[ScriptL]) + 9856*E^((4*I)*\[ScriptL]))*
            \[Kappa]A*(-1 + 4*\[Nu]) + 2*(17372 - 73737*\[Nu] + 
             39424*\[Nu]^2 + 16*E^((4*I)*\[ScriptL])*(17372 - 73737*\[Nu] + 
               39424*\[Nu]^2) + 5*E^((2*I)*\[ScriptL])*(234187 - 994977*
                \[Nu] + 532544*\[Nu]^2))*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*(616 + 41605*E^((2*I)*\[ScriptL]) + 9856*
                E^((4*I)*\[ScriptL]))*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (17372 - 56672*\[Nu] - 64*E^((4*I)*\[ScriptL])*(-4343 + 
                 14168*\[Nu]) - 5*E^((2*I)*\[ScriptL])*(-234187 + 
                 765532*\[Nu]))*\[Chi]A^2 + (17372 - 21314*\[Nu] + 4928*
                \[Nu]^2 + 32*E^((4*I)*\[ScriptL])*(8686 - 10657*\[Nu] + 
                 2464*\[Nu]^2) + 5*E^((2*I)*\[ScriptL])*(234187 - 
                 287674*\[Nu] + 66568*\[Nu]^2))*\[Chi]S^2)))/
         (58320*Sqrt[14]*E^((2*I)*\[ScriptL])) - 
        (e^3*\[Nu]*(9*(-308 + 11199*E^((2*I)*\[ScriptL]) + 
             2772*E^((4*I)*\[ScriptL]))*\[Kappa]A*(-1 + 4*\[Nu]) + 
           (-17372 + 73737*\[Nu] - 39424*\[Nu]^2 + 9*E^((4*I)*\[ScriptL])*
              (17372 - 73737*\[Nu] + 39424*\[Nu]^2) + 6*E^((2*I)*\[ScriptL])*
              (105011 - 446256*\[Nu] + 238912*\[Nu]^2))*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*(-308 + 11199*E^((2*I)*\[ScriptL]) + 2772*
                E^((4*I)*\[ScriptL]))*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (-8686 + E^((2*I)*\[ScriptL])*(315033 - 1030308*\[Nu]) + 
               E^((4*I)*\[ScriptL])*(78174 - 255024*\[Nu]) + 28336*\[Nu])*
              \[Chi]A^2 + (-8686 + 10657*\[Nu] - 2464*\[Nu]^2 + 9*
                E^((4*I)*\[ScriptL])*(8686 - 10657*\[Nu] + 2464*\[Nu]^2) + 
               E^((2*I)*\[ScriptL])*(315033 - 387096*\[Nu] + 89592*\[Nu]^2))*
              \[Chi]S^2)))/(19440*Sqrt[14]*E^(I*\[ScriptL])) - 
        (e^5*\[Nu]*(18*(-4158 - 39757*E^((2*I)*\[ScriptL]) + 
             643224*E^((4*I)*\[ScriptL]) + 269109*E^((6*I)*\[ScriptL]) + 
             96250*E^((8*I)*\[ScriptL]))*\[Kappa]A*(-1 + 4*\[Nu]) + 
           (-27*(17372 - 73737*\[Nu] + 39424*\[Nu]^2) + 
             625*E^((8*I)*\[ScriptL])*(17372 - 73737*\[Nu] + 39424*\[Nu]^2) + 
             108*E^((6*I)*\[ScriptL])*(280289 - 1191294*\[Nu] + 637888*
                \[Nu]^2) - 4*E^((2*I)*\[ScriptL])*(1118819 - 4753674*\[Nu] + 
               2544448*\[Nu]^2) + 6*E^((4*I)*\[ScriptL])*(12052972 - 51240087*
                \[Nu] + 27444224*\[Nu]^2))*\[Chi]A*\[Chi]S + 
           \[Delta]*(18*(-4158 - 39757*E^((2*I)*\[ScriptL]) + 643224*
                E^((4*I)*\[ScriptL]) + 269109*E^((6*I)*\[ScriptL]) + 96250*
                E^((8*I)*\[ScriptL]))*\[Kappa]S*(-1 + 2*\[Nu]) - 
             2*(117261 + E^((2*I)*\[ScriptL])*(1118819 - 3657644*\[Nu]) - 
               382536*\[Nu] + 625*E^((8*I)*\[ScriptL])*(-4343 + 
                 14168*\[Nu]) + 27*E^((6*I)*\[ScriptL])*(-280289 + 
                 916964*\[Nu]) + 6*E^((4*I)*\[ScriptL])*(-3013243 + 
                 9862768*\[Nu]))*\[Chi]A^2 + (E^((2*I)*\[ScriptL])*
                (-2237638 + 2748856*\[Nu] - 636112*\[Nu]^2) - 27*
                (8686 - 10657*\[Nu] + 2464*\[Nu]^2) + 625*
                E^((8*I)*\[ScriptL])*(8686 - 10657*\[Nu] + 2464*\[Nu]^2) + 54*
                E^((6*I)*\[ScriptL])*(280289 - 344468*\[Nu] + 79736*
                  \[Nu]^2) + 6*E^((4*I)*\[ScriptL])*(6026486 - 7408607*
                  \[Nu] + 1715264*\[Nu]^2))*\[Chi]S^2)))/
         (933120*Sqrt[14]*E^((3*I)*\[ScriptL])) - 
        (e^6*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
          (9*(9856 + 75091*E^((2*I)*\[ScriptL]) + 2486761*E^((4*I)*
                \[ScriptL]) + 1047456*E^((6*I)*\[ScriptL]) + 
             299376*E^((8*I)*\[ScriptL]))*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(16*(17372 - 73737*\[Nu] + 39424*\[Nu]^2) + 
             486*E^((8*I)*\[ScriptL])*(17372 - 73737*\[Nu] + 39424*\[Nu]^2) + 
             18*E^((6*I)*\[ScriptL])*(1637564 - 6957819*\[Nu] + 3724288*
                \[Nu]^2) + E^((2*I)*\[ScriptL])*(2113697 - 8979687*\[Nu] + 
               4805824*\[Nu]^2) + E^((4*I)*\[ScriptL])*(69927587 - 297217377*
                \[Nu] + 159152704*\[Nu]^2))*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*(9856 + 75091*E^((2*I)*\[ScriptL]) + 2486761*
                E^((4*I)*\[ScriptL]) + 1047456*E^((6*I)*\[ScriptL]) + 299376*
                E^((8*I)*\[ScriptL]))*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (277952 + E^((4*I)*\[ScriptL])*(69927587 - 228782012*\[Nu]) + 
               E^((2*I)*\[ScriptL])*(2113697 - 6908372*\[Nu]) - 906752*
                \[Nu] - 1944*E^((8*I)*\[ScriptL])*(-4343 + 14168*\[Nu]) - 72*
                E^((6*I)*\[ScriptL])*(-409391 + 1338416*\[Nu]))*\[Chi]A^2 + 
             (32*(8686 - 10657*\[Nu] + 2464*\[Nu]^2) + 972*
                E^((8*I)*\[ScriptL])*(8686 - 10657*\[Nu] + 2464*\[Nu]^2) + 36*
                E^((6*I)*\[ScriptL])*(818782 - 1005859*\[Nu] + 232768*
                  \[Nu]^2) + E^((2*I)*\[ScriptL])*(2113697 - 2596214*\[Nu] + 
                 600728*\[Nu]^2) + E^((4*I)*\[ScriptL])*(69927587 - 
                 85942394*\[Nu] + 19894088*\[Nu]^2))*\[Chi]S^2)))/
         (1166400*Sqrt[14]*E^((4*I)*\[ScriptL]))))
 
HoscMem[3, 3] = x^3*\[Epsilon]^4*
     ((-269*e^3*E^((3*I)*\[ScriptL])*\[Delta]*\[Nu])/(648*Sqrt[210]) - 
      (269*e^4*E^((2*I)*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Delta]*
        \[Nu])/(216*Sqrt[210]) - (e^5*E^(I*\[ScriptL])*
        (5649 - 12754*E^((2*I)*\[ScriptL]) + 13719*E^((4*I)*\[ScriptL]))*
        \[Delta]*\[Nu])/(5184*Sqrt[210]) - 
      (e^6*(-538 - 1044*E^((2*I)*\[ScriptL]) - 11061*E^((4*I)*\[ScriptL]) + 
         12643*E^((6*I)*\[ScriptL]))*\[Delta]*\[Nu])/(2592*Sqrt[210])) + 
    SO*x^(7/2)*\[Epsilon]^5*((e^3*E^((3*I)*\[ScriptL])*\[Nu]*
        (4*(-121 + 538*\[Nu])*\[Chi]A + \[Delta]*(-484 + 269*\[Nu])*\[Chi]S))/
       (972*Sqrt[210]) + (e^4*E^((2*I)*\[ScriptL])*
        (-1 + E^((2*I)*\[ScriptL]))*\[Nu]*(4*(-121 + 538*\[Nu])*\[Chi]A + 
         \[Delta]*(-484 + 269*\[Nu])*\[Chi]S))/(324*Sqrt[210]) + 
      (e^5*E^(I*\[ScriptL])*\[Nu]*
        ((E^((2*I)*\[ScriptL])*(42041 - 186848*\[Nu]) + 
           168*(-121 + 538*\[Nu]) + 408*E^((4*I)*\[ScriptL])*
            (-121 + 538*\[Nu]))*\[Chi]A + \[Delta]*
          (E^((2*I)*\[ScriptL])*(42041 - 23356*\[Nu]) + 
           42*(-484 + 269*\[Nu]) + 102*E^((4*I)*\[ScriptL])*
            (-484 + 269*\[Nu]))*\[Chi]S))/(15552*Sqrt[210]) + 
      (e^6*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
        ((32*(-121 + 538*\[Nu]) + 752*E^((4*I)*\[ScriptL])*
            (-121 + 538*\[Nu]) + E^((2*I)*\[ScriptL])*(-22949 + 
             102272*\[Nu]))*\[Chi]A + \[Delta]*(8*(-484 + 269*\[Nu]) + 
           188*E^((4*I)*\[ScriptL])*(-484 + 269*\[Nu]) + E^((2*I)*\[ScriptL])*
            (-22949 + 12784*\[Nu]))*\[Chi]S))/(15552*Sqrt[210])) + 
    \[Epsilon]^6*(x^4*((11*\[Delta]*\[Nu])/(27*Sqrt[210]) + 
        (e*(33 + E^((2*I)*\[ScriptL]))*\[Delta]*\[Nu])/
         (18*Sqrt[210]*E^(I*\[ScriptL])) - 
        (e^2*(-713 - 570*E^((2*I)*\[ScriptL]) + 595*E^((4*I)*\[ScriptL]))*
          \[Delta]*\[Nu])/(180*Sqrt[210]*E^((2*I)*\[ScriptL])) - 
        (e^4*\[Delta]*\[Nu]*(-612744 - 303083*E^((2*I)*\[ScriptL]) - 
           275649*E^((4*I)*\[ScriptL]) - 5*E^((6*I)*\[ScriptL])*
            (-419447 + 349080*\[Nu]) + 5*E^((8*I)*\[ScriptL])*
            (-346099 + 349080*\[Nu])))/(47520*Sqrt[210]*
          E^((4*I)*\[ScriptL])) - (e^3*\[Delta]*\[Nu]*(-1056132 - 
           662013*E^((2*I)*\[ScriptL]) - 1219680*E^((4*I)*\[ScriptL]) + 
           5*E^((6*I)*\[ScriptL])*(-268351 + 354998*\[Nu])))/
         (142560*Sqrt[210]*E^((3*I)*\[ScriptL])) - 
        (e^5*\[Delta]*\[Nu]*(-24493326 - 8776944*E^((2*I)*\[ScriptL]) - 
           9541686*E^((4*I)*\[ScriptL]) + E^((8*I)*\[ScriptL])*
            (71473275 - 58127180*\[Nu]) + 30*E^((10*I)*\[ScriptL])*
            (-3130156 + 2922795*\[Nu]) + 6*E^((6*I)*\[ScriptL])*
            (-12171484 + 5975745*\[Nu])))/(1140480*Sqrt[210]*
          E^((5*I)*\[ScriptL])) - (e^6*\[Delta]*\[Nu]*(-39587427 - 
           8702001*E^((2*I)*\[ScriptL]) - 12173502*E^((4*I)*\[ScriptL]) + 
           E^((8*I)*\[ScriptL])*(98457378 - 82937820*\[Nu]) + 
           E^((10*I)*\[ScriptL])*(83203062 - 69887220*\[Nu]) + 
           E^((6*I)*\[ScriptL])*(1148696 - 6567340*\[Nu]) + 
           2*E^((12*I)*\[ScriptL])*(-88612207 + 79696190*\[Nu])))/
         (1140480*Sqrt[210]*E^((6*I)*\[ScriptL]))) + 
      SO^2*x^4*(-(e^3*E^((3*I)*\[ScriptL])*\[Nu]*
           (2421*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(5323 - 27108*\[Nu] + 
              17216*\[Nu]^2)*\[Chi]A*\[Chi]S + \[Delta]*
             (2421*\[Kappa]S*(-1 + 2*\[Nu]) + (5323 - 24748*\[Nu])*\[Chi]A^
                2 + (5323 - 8176*\[Nu] + 2152*\[Nu]^2)*\[Chi]S^2)))/
         (11664*Sqrt[210]) - (e^4*E^((2*I)*\[ScriptL])*
          (-1 + E^((2*I)*\[ScriptL]))*\[Nu]*(2421*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(5323 - 27108*\[Nu] + 17216*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(2421*\[Kappa]S*(-1 + 2*\[Nu]) + (5323 - 24748*\[Nu])*
              \[Chi]A^2 + (5323 - 8176*\[Nu] + 2152*\[Nu]^2)*\[Chi]S^2)))/
         (3888*Sqrt[210]) - (e^6*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
          (9*(538 + 4810*E^((2*I)*\[ScriptL]) + 12643*E^((4*I)*\[ScriptL]))*
            \[Kappa]A*(-1 + 4*\[Nu]) + 2*(10646 - 54216*\[Nu] + 
             34432*\[Nu]^2 + 47*E^((4*I)*\[ScriptL])*(5323 - 27108*\[Nu] + 
               17216*\[Nu]^2) + 5*E^((2*I)*\[ScriptL])*(18994 - 96849*\[Nu] + 
               61568*\[Nu]^2))*\[Chi]A*\[Chi]S + \[Delta]*
            (9*(538 + 4810*E^((2*I)*\[ScriptL]) + 12643*E^((4*I)*\[ScriptL]))*
              \[Kappa]S*(-1 + 2*\[Nu]) + (10646 + E^((2*I)*\[ScriptL])*
                (94970 - 442520*\[Nu]) - 49496*\[Nu] - 47*
                E^((4*I)*\[ScriptL])*(-5323 + 24748*\[Nu]))*\[Chi]A^2 + 
             (2*(5323 - 8176*\[Nu] + 2152*\[Nu]^2) + 47*E^((4*I)*\[ScriptL])*
                (5323 - 8176*\[Nu] + 2152*\[Nu]^2) + 10*E^((2*I)*\[ScriptL])*
                (9497 - 14609*\[Nu] + 3848*\[Nu]^2))*\[Chi]S^2)))/
         (46656*Sqrt[210]) - (e^5*E^(I*\[ScriptL])*\[Nu]*
          (9*(1883 - 3534*E^((2*I)*\[ScriptL]) + 4573*E^((4*I)*\[ScriptL]))*
            \[Kappa]A*(-1 + 4*\[Nu]) + 2*(7*(5323 - 27108*\[Nu] + 17216*
                \[Nu]^2) + 17*E^((4*I)*\[ScriptL])*(5323 - 27108*\[Nu] + 
               17216*\[Nu]^2) - 6*E^((2*I)*\[ScriptL])*(11663 - 59373*\[Nu] + 
               37696*\[Nu]^2))*\[Chi]A*\[Chi]S + \[Delta]*
            (9*(1883 - 3534*E^((2*I)*\[ScriptL]) + 4573*E^((4*I)*\[ScriptL]))*
              \[Kappa]S*(-1 + 2*\[Nu]) + (37261 + E^((4*I)*\[ScriptL])*
                (90491 - 420716*\[Nu]) - 173236*\[Nu] + 6*
                E^((2*I)*\[ScriptL])*(-11663 + 54188*\[Nu]))*\[Chi]A^2 + 
             (7*(5323 - 8176*\[Nu] + 2152*\[Nu]^2) + 17*E^((4*I)*\[ScriptL])*
                (5323 - 8176*\[Nu] + 2152*\[Nu]^2) - 6*E^((2*I)*\[ScriptL])*
                (11663 - 17906*\[Nu] + 4712*\[Nu]^2))*\[Chi]S^2)))/
         (31104*Sqrt[210])))
 
HoscMem[4, 0] = x^(7/2)*\[Epsilon]^5*
    ((((4*I)/105)*Sqrt[2]*e*(-1 + E^((2*I)*\[ScriptL]))*\[Nu])/
      E^(I*\[ScriptL]) + (((143*I)/1680)*e^2*(-1 + E^((4*I)*\[ScriptL]))*
       \[Nu])/(Sqrt[2]*E^((2*I)*\[ScriptL])) + 
     ((I/1890)*e^3*(-211 - 567*E^((2*I)*\[ScriptL]) + 
        567*E^((4*I)*\[ScriptL]) + 211*E^((6*I)*\[ScriptL]))*\[Nu])/
      (Sqrt[2]*E^((3*I)*\[ScriptL])) + 
     ((I/8064)*e^4*(-1235 - 1928*E^((2*I)*\[ScriptL]) + 
        1928*E^((6*I)*\[ScriptL]) + 1235*E^((8*I)*\[ScriptL]))*\[Nu])/
      (Sqrt[2]*E^((4*I)*\[ScriptL])) + 
     ((I/9450)*e^5*(-2019 - 2200*E^((2*I)*\[ScriptL]) - 
        6825*E^((4*I)*\[ScriptL]) + 6825*E^((6*I)*\[ScriptL]) + 
        2200*E^((8*I)*\[ScriptL]) + 2019*E^((10*I)*\[ScriptL]))*\[Nu])/
      (Sqrt[2]*E^((5*I)*\[ScriptL])) + 
     ((I/67200)*e^6*(-20168 - 15773*E^((2*I)*\[ScriptL]) - 
        34160*E^((4*I)*\[ScriptL]) + 34160*E^((8*I)*\[ScriptL]) + 
        15773*E^((10*I)*\[ScriptL]) + 20168*E^((12*I)*\[ScriptL]))*\[Nu])/
      (Sqrt[2]*E^((6*I)*\[ScriptL])))
 
HoscMem[4, 2] = x^(5/2)*\[Epsilon]^3*
     ((((-13*I)/3024)*e^2*E^((2*I)*\[ScriptL])*\[Nu])/Sqrt[5] - 
      (((13*I)/1512)*e^3*E^(I*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Nu])/
       Sqrt[5] - ((I/12096)*e^4*(39 - 70*E^((2*I)*\[ScriptL]) + 
         169*E^((4*I)*\[ScriptL]))*\[Nu])/Sqrt[5] - 
      ((I/36288)*e^5*(13 - 555*E^((2*I)*\[ScriptL]) - 
         225*E^((4*I)*\[ScriptL]) + 767*E^((6*I)*\[ScriptL]))*\[Nu])/
       (Sqrt[5]*E^(I*\[ScriptL])) - 
      ((I/145152)*e^6*(26 + 1320*E^((2*I)*\[ScriptL]) - 
         1509*E^((4*I)*\[ScriptL]) - 1352*E^((6*I)*\[ScriptL]) + 
         4485*E^((8*I)*\[ScriptL]))*\[Nu])/(Sqrt[5]*E^((2*I)*\[ScriptL]))) + 
    SO*x^3*\[Epsilon]^4*((((13*I)/4536)*e^2*E^((2*I)*\[ScriptL])*\[Nu]*
        (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/Sqrt[5] + 
      (((13*I)/2268)*e^3*E^(I*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
        (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/Sqrt[5] + 
      ((I/18144)*e^4*(39 - 44*E^((2*I)*\[ScriptL]) + 
         169*E^((4*I)*\[ScriptL]))*\[Nu]*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu])*\[Chi]S))/Sqrt[5] + 
      ((I/54432)*e^5*(13 - 711*E^((2*I)*\[ScriptL]) - 
         69*E^((4*I)*\[ScriptL]) + 767*E^((6*I)*\[ScriptL]))*\[Nu]*
        (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/
       (Sqrt[5]*E^(I*\[ScriptL])) + 
      ((I/217728)*e^6*(26 + 1554*E^((2*I)*\[ScriptL]) - 
         1695*E^((4*I)*\[ScriptL]) - 338*E^((6*I)*\[ScriptL]) + 
         4485*E^((8*I)*\[ScriptL]))*\[Nu]*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu])*\[Chi]S))/(Sqrt[5]*E^((2*I)*\[ScriptL]))) + 
    \[Epsilon]^5*(x^(7/2)*(((-I/189)*e*(-1 + 3*E^((2*I)*\[ScriptL]))*\[Nu])/
         (Sqrt[5]*E^(I*\[ScriptL])) - ((I/199584)*e^3*\[Nu]*
          (-3762 - 1782*E^((2*I)*\[ScriptL]) + E^((4*I)*\[ScriptL])*
            (275055 - 710146*\[Nu]) + E^((6*I)*\[ScriptL])*
            (-264231 + 710146*\[Nu])))/(Sqrt[5]*E^((3*I)*\[ScriptL])) - 
        ((I/399168)*e^2*\[Nu]*(-4290 - 5016*E^((2*I)*\[ScriptL]) + 
           E^((4*I)*\[ScriptL])*(-262119 + 711004*\[Nu])))/
         (Sqrt[5]*E^((2*I)*\[ScriptL])) - ((I/1596672)*e^4*\[Nu]*
          (-49478 - 16236*E^((2*I)*\[ScriptL]) + 3*E^((4*I)*\[ScriptL])*
            (-287969 + 708716*\[Nu]) + 2*E^((6*I)*\[ScriptL])*
            (-243221 + 726760*\[Nu]) + E^((8*I)*\[ScriptL])*
            (-3442263 + 9222460*\[Nu])))/(Sqrt[5]*E^((4*I)*\[ScriptL])) - 
        ((I/4790016)*e^5*\[Nu]*(-235092 - 46332*E^((2*I)*\[ScriptL]) + 
           E^((6*I)*\[ScriptL])*(23528307 - 61881090*\[Nu]) + 
           99*E^((8*I)*\[ScriptL])*(-71153 + 195502*\[Nu]) + 
           E^((4*I)*\[ScriptL])*(-343101 + 708430*\[Nu]) + 
           E^((10*I)*\[ScriptL])*(-15637875 + 41817962*\[Nu])))/
         (Sqrt[5]*E^((5*I)*\[ScriptL])) - ((I/95800320)*e^6*\[Nu]*
          (-7266402 - 481536*E^((2*I)*\[ScriptL]) + 20*E^((4*I)*\[ScriptL])*
            (-212025 + 353786*\[Nu]) - 240*E^((8*I)*\[ScriptL])*
            (-844467 + 2213690*\[Nu]) + 60*E^((6*I)*\[ScriptL])*
            (-3859449 + 9935450*\[Nu]) + 4*E^((10*I)*\[ScriptL])*
            (-60248397 + 164745050*\[Nu]) + 3*E^((12*I)*\[ScriptL])*
            (-152510193 + 407214260*\[Nu])))/(Sqrt[5]*
          E^((6*I)*\[ScriptL]))) + SO^2*x^(7/2)*
       ((((13*I)/54432)*e^2*E^((2*I)*\[ScriptL])*\[Nu]*
          (9*\[Delta]*\[Kappa]A + \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 
           92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 
           23*\[Chi]S^2 + 32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2))/
         Sqrt[5] + ((I/217728)*e^4*(39 - 18*E^((2*I)*\[ScriptL]) + 
           169*E^((4*I)*\[ScriptL]))*\[Nu]*(9*\[Delta]*\[Kappa]A + 
           \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 92*\[Nu]*\[Chi]A^2 + 
           2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 23*\[Chi]S^2 + 
           32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2))/Sqrt[5] + 
        ((I/653184)*e^5*(13 - 867*E^((2*I)*\[ScriptL]) + 
           87*E^((4*I)*\[ScriptL]) + 767*E^((6*I)*\[ScriptL]))*\[Nu]*
          (9*\[Delta]*\[Kappa]A + \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 
           92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 
           23*\[Chi]S^2 + 32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2))/
         (Sqrt[5]*E^(I*\[ScriptL])) + ((I/2612736)*e^6*
          (26 + 1788*E^((2*I)*\[ScriptL]) - 1725*E^((4*I)*\[ScriptL]) + 
           676*E^((6*I)*\[ScriptL]) + 4485*E^((8*I)*\[ScriptL]))*\[Nu]*
          (9*\[Delta]*\[Kappa]A + \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 
           92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 
           23*\[Chi]S^2 + 32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2))/
         (Sqrt[5]*E^((2*I)*\[ScriptL])) - (((13*I)/27216)*e^3*
          E^(I*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
          (-9*\[Delta]*\[Kappa]A + 9*\[Kappa]S*(-1 + 2*\[Nu]) + 
           23*\[Chi]A^2 - 92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(23 - 16*\[Nu])*
            \[Chi]A*\[Chi]S + 23*\[Chi]S^2 - 32*\[Nu]*\[Chi]S^2 + 
           8*\[Nu]^2*\[Chi]S^2))/Sqrt[5])) + 
    \[Epsilon]^6*(x^4*((((-29*I)/1512)*e^2*E^((2*I)*\[ScriptL])*Pi*\[Nu])/
         Sqrt[5] - (((29*I)/756)*e^3*E^(I*\[ScriptL])*
          (-1 + E^((2*I)*\[ScriptL]))*Pi*\[Nu])/Sqrt[5] - 
        ((I/18144)*e^4*(261 + 473*E^((2*I)*\[ScriptL]) + 
           1131*E^((4*I)*\[ScriptL]))*Pi*\[Nu])/Sqrt[5] - 
        ((I/18144)*e^5*(29 - 3121*E^((2*I)*\[ScriptL]) + 
           1381*E^((4*I)*\[ScriptL]) + 1711*E^((6*I)*\[ScriptL]))*Pi*\[Nu])/
         (Sqrt[5]*E^(I*\[ScriptL])) - ((I/290304)*e^6*
          (232 + 23076*E^((2*I)*\[ScriptL]) - 13897*E^((4*I)*\[ScriptL]) + 
           36892*E^((6*I)*\[ScriptL]) + 40020*E^((8*I)*\[ScriptL]))*Pi*\[Nu])/
         (Sqrt[5]*E^((2*I)*\[ScriptL]))) + 
      SO*x^4*(((I/199584)*e^2*E^((2*I)*\[ScriptL])*\[Nu]*
          (\[Delta]*(206184 - 474241*\[Nu])*\[Chi]A + 
           4*(51546 - 177424*\[Nu] + 59441*\[Nu]^2)*\[Chi]S))/Sqrt[5] + 
        ((I/99792)*e^3*E^(I*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
          (\[Delta]*(204754 - 473669*\[Nu])*\[Chi]A + 
           2*(102377 - 353990*\[Nu] + 118739*\[Nu]^2)*\[Chi]S))/Sqrt[5] + 
        ((I/798336)*e^4*\[Nu]*(-(\[Delta]*(-607112 + 1418147*\[Nu] + 
              8*E^((2*I)*\[ScriptL])*(-106733 + 238883*\[Nu]) + 
              13*E^((4*I)*\[ScriptL])*(-203544 + 473185*\[Nu]))*\[Chi]A) + 
           4*(151778 - 528840*\[Nu] + 177751*\[Nu]^2 + 
             13*E^((4*I)*\[ScriptL])*(50886 - 176632*\[Nu] + 59309*\[Nu]^2) + 
             E^((2*I)*\[ScriptL])*(213466 - 716433*\[Nu] + 238982*\[Nu]^2))*
            \[Chi]S))/Sqrt[5] + ((I/2395008)*e^5*(-1 + E^((2*I)*\[ScriptL]))*
          \[Nu]*(-(\[Delta]*(201894 - 472525*\[Nu] + 10*E^((2*I)*\[ScriptL])*(
                -2009238 + 4644485*\[Nu]) + E^((4*I)*\[ScriptL])*(-11946066 + 
                27892703*\[Nu]))*\[Chi]A) + 2*(-100947 + 352274*\[Nu] - 
             118453*\[Nu]^2 + 2*E^((2*I)*\[ScriptL])*(5023095 - 17341448*
                \[Nu] + 5818105*\[Nu]^2) + E^((4*I)*\[ScriptL])*
              (5973033 - 20804758*\[Nu] + 6992159*\[Nu]^2))*\[Chi]S))/
         (Sqrt[5]*E^(I*\[ScriptL])) + ((I/9580032)*e^6*\[Nu]*
          (-(\[Delta]*(-400928 + E^((4*I)*\[ScriptL])*(23784192 - 
                56532531*\[Nu]) + 943906*\[Nu] + 6*E^((2*I)*\[ScriptL])*(
                -6283156 + 14661725*\[Nu]) + 3*E^((8*I)*\[ScriptL])*(
                -23173480 + 54322643*\[Nu]) + 2*E^((6*I)*\[ScriptL])*(
                -27260420 + 62273449*\[Nu]))*\[Chi]A) + 
           (8*(50116 - 175708*\[Nu] + 59155*\[Nu]^2) + 
             12*E^((2*I)*\[ScriptL])*(3141578 - 10929155*\[Nu] + 3673772*
                \[Nu]^2) + 12*E^((8*I)*\[ScriptL])*(5793370 - 20242456*
                \[Nu] + 6808831*\[Nu]^2) - 3*E^((4*I)*\[ScriptL])*
              (7928064 - 28014157*\[Nu] + 9448472*\[Nu]^2) + 
             4*E^((6*I)*\[ScriptL])*(13630210 - 46557223*\[Nu] + 15589276*
                \[Nu]^2))*\[Chi]S))/(Sqrt[5]*E^((2*I)*\[ScriptL]))) + 
      SO^3*x^4*((((-13*I)/40824)*e^2*E^((2*I)*\[ScriptL])*\[Nu]*
          (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
              \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
              \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/Sqrt[5] - 
        (((13*I)/20412)*e^3*E^(I*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*
          \[Nu]*(9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
              \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
              \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/Sqrt[5] - 
        ((I/163296)*e^4*(39 + 8*E^((2*I)*\[ScriptL]) + 
           169*E^((4*I)*\[ScriptL]))*\[Nu]*
          (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
              \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
              \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/Sqrt[5] - 
        ((I/489888)*e^5*(13 - 1023*E^((2*I)*\[ScriptL]) + 
           243*E^((4*I)*\[ScriptL]) + 767*E^((6*I)*\[ScriptL]))*\[Nu]*
          (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
              \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
              \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/
         (Sqrt[5]*E^(I*\[ScriptL])) - ((I/1959552)*e^6*
          (26 + 2022*E^((2*I)*\[ScriptL]) - 1599*E^((4*I)*\[ScriptL]) + 
           1690*E^((6*I)*\[ScriptL]) + 4485*E^((8*I)*\[ScriptL]))*\[Nu]*
          (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
              \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
              \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/
         (Sqrt[5]*E^((2*I)*\[ScriptL]))))
 
HoscMem[4, 4] = x^(5/2)*\[Epsilon]^3*
     (((5*I)/6912)*Sqrt[5/7]*e^4*E^((4*I)*\[ScriptL])*\[Nu] + 
      ((5*I)/1728)*Sqrt[5/7]*e^5*E^((3*I)*\[ScriptL])*
       (-1 + E^((2*I)*\[ScriptL]))*\[Nu] + ((5*I)/13824)*Sqrt[5/7]*e^6*
       E^((2*I)*\[ScriptL])*(11 - 27*E^((2*I)*\[ScriptL]) + 
        21*E^((4*I)*\[ScriptL]))*\[Nu]) + SO*x^3*\[Epsilon]^4*
     (((5*I)/10368)*Sqrt[5/7]*e^4*E^((4*I)*\[ScriptL])*\[Nu]*
       (2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S) + 
      ((5*I)/2592)*Sqrt[5/7]*e^5*E^((3*I)*\[ScriptL])*
       (-1 + E^((2*I)*\[ScriptL]))*\[Nu]*(2*\[Delta]*\[Chi]A - 
        (-2 + \[Nu])*\[Chi]S) - ((5*I)/20736)*Sqrt[5/7]*e^6*
       E^((2*I)*\[ScriptL])*(11 - 26*E^((2*I)*\[ScriptL]) + 
        21*E^((4*I)*\[ScriptL]))*\[Nu]*(-2*\[Delta]*\[Chi]A + 
        (-2 + \[Nu])*\[Chi]S)) + \[Epsilon]^5*
     (x^(7/2)*(((I/9)*\[Nu])/Sqrt[35] + 
        ((I/45)*e*(7 + 15*E^((2*I)*\[ScriptL]))*\[Nu])/
         (Sqrt[35]*E^(I*\[ScriptL])) + 
        ((I/4320)*e^2*(1037 + 993*E^((2*I)*\[ScriptL]) + 
           3255*E^((4*I)*\[ScriptL]))*\[Nu])/(Sqrt[35]*
          E^((2*I)*\[ScriptL])) + ((I/1080)*e^3*
          (394 + 317*E^((2*I)*\[ScriptL]) + 239*E^((4*I)*\[ScriptL]) + 
           1650*E^((6*I)*\[ScriptL]))*\[Nu])/(Sqrt[35]*
          E^((3*I)*\[ScriptL])) - ((I/4561920)*e^4*\[Nu]*
          (-2492886 - 1581360*E^((2*I)*\[ScriptL]) - 
           1709400*E^((4*I)*\[ScriptL]) + 426096*E^((6*I)*\[ScriptL]) + 
           5*E^((8*I)*\[ScriptL])*(-2754285 + 288332*\[Nu])))/
         (Sqrt[35]*E^((4*I)*\[ScriptL])) - ((I/1140480)*e^5*\[Nu]*
          (-922900 - 436524*E^((2*I)*\[ScriptL]) - 
           495264*E^((4*I)*\[ScriptL]) - 590216*E^((6*I)*\[ScriptL]) + 
           E^((8*I)*\[ScriptL])*(1792359 - 1449910*\[Nu]) + 
           5*E^((10*I)*\[ScriptL])*(-1303539 + 289982*\[Nu])))/
         (Sqrt[35]*E^((5*I)*\[ScriptL])) - ((I/45619200)*e^6*\[Nu]*
          (-54176298 - 17196410*E^((2*I)*\[ScriptL]) - 
           22597630*E^((4*I)*\[ScriptL]) - 23166000*E^((6*I)*\[ScriptL]) + 
           55*E^((8*I)*\[ScriptL])*(-1249333 + 1459660*\[Nu]) + 
           15*E^((12*I)*\[ScriptL])*(-31510139 + 10201620*\[Nu]) - 
           4*E^((10*I)*\[ScriptL])*(-55303447 + 45587775*\[Nu])))/
         (Sqrt[35]*E^((6*I)*\[ScriptL]))) + SO^2*x^(7/2)*
       (((5*I)/124416)*Sqrt[5/7]*e^4*E^((4*I)*\[ScriptL])*\[Nu]*
         (-9*\[Delta]*\[Kappa]A + 9*\[Kappa]S*(-1 + 2*\[Nu]) + 23*\[Chi]A^2 - 
          92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(23 - 16*\[Nu])*\[Chi]A*\[Chi]S + 
          23*\[Chi]S^2 - 32*\[Nu]*\[Chi]S^2 + 8*\[Nu]^2*\[Chi]S^2) + 
        ((5*I)/31104)*Sqrt[5/7]*e^5*E^((3*I)*\[ScriptL])*
         (-1 + E^((2*I)*\[ScriptL]))*\[Nu]*(-9*\[Delta]*\[Kappa]A + 
          9*\[Kappa]S*(-1 + 2*\[Nu]) + 23*\[Chi]A^2 - 92*\[Nu]*\[Chi]A^2 + 
          2*\[Delta]*(23 - 16*\[Nu])*\[Chi]A*\[Chi]S + 23*\[Chi]S^2 - 
          32*\[Nu]*\[Chi]S^2 + 8*\[Nu]^2*\[Chi]S^2) + 
        ((5*I)/248832)*Sqrt[5/7]*e^6*E^((2*I)*\[ScriptL])*
         (11 - 25*E^((2*I)*\[ScriptL]) + 21*E^((4*I)*\[ScriptL]))*\[Nu]*
         (-9*\[Delta]*\[Kappa]A + 9*\[Kappa]S*(-1 + 2*\[Nu]) + 23*\[Chi]A^2 - 
          92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(23 - 16*\[Nu])*\[Chi]A*\[Chi]S + 
          23*\[Chi]S^2 - 32*\[Nu]*\[Chi]S^2 + 8*\[Nu]^2*\[Chi]S^2))) + 
    \[Epsilon]^6*(x^4*((((19*I)/1152)*e^4*E^((4*I)*\[ScriptL])*Pi*\[Nu])/
         Sqrt[35] + (((19*I)/288)*e^5*E^((3*I)*\[ScriptL])*
          (-1 + E^((2*I)*\[ScriptL]))*Pi*\[Nu])/Sqrt[35] + 
        ((I/414720)*e^6*E^((2*I)*\[ScriptL])*(37620 - 
           76793*E^((2*I)*\[ScriptL]) + 71820*E^((4*I)*\[ScriptL]))*Pi*\[Nu])/
         Sqrt[35]) + SO*x^4*(((I/152064)*e^4*E^((4*I)*\[ScriptL])*\[Nu]*
          (\[Delta]*(20960 - 63921*\[Nu])*\[Chi]A + 
           4*(5240 - 17008*\[Nu] + 7887*\[Nu]^2)*\[Chi]S))/Sqrt[35] + 
        ((I/114048)*e^5*E^((3*I)*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*
          \[Nu]*(\[Delta]*(65630 - 192863*\[Nu])*\[Chi]A + 
           2*(32815 - 103698*\[Nu] + 47597*\[Nu]^2)*\[Chi]S))/Sqrt[35] + 
        ((I/912384)*e^6*E^((2*I)*\[ScriptL])*\[Nu]*
          (-(\[Delta]*(33*(-22960 + 64721*\[Nu]) - 8*E^((2*I)*\[ScriptL])*(
                -208460 + 581701*\[Nu]) + E^((4*I)*\[ScriptL])*(-1430480 + 
                4071023*\[Nu]))*\[Chi]A) + 
           (E^((2*I)*\[ScriptL])*(-1667680 + 5076631*\[Nu] - 2297104*
                \[Nu]^2) + 132*(5740 - 17608*\[Nu] + 7987*\[Nu]^2) + 
             4*E^((4*I)*\[ScriptL])*(357620 - 1104504*\[Nu] + 502381*
                \[Nu]^2))*\[Chi]S))/Sqrt[35]) + 
      SO^3*x^4*(((5*I)/93312)*Sqrt[5/7]*e^4*E^((4*I)*\[ScriptL])*\[Nu]*
         (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
             \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
            (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
             \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
            3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
            (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)) + 
        ((5*I)/23328)*Sqrt[5/7]*e^5*E^((3*I)*\[ScriptL])*
         (-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
         (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
             \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
            (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
             \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
            3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
            (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)) + 
        ((5*I)/186624)*Sqrt[5/7]*e^6*E^((2*I)*\[ScriptL])*
         (11 - 24*E^((2*I)*\[ScriptL]) + 21*E^((4*I)*\[ScriptL]))*\[Nu]*
         (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
             \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
            (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
             \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
            3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
            (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2))))
 
HoscMem[5, 1] = x^3*\[Epsilon]^4*
     ((-43*e*E^(I*\[ScriptL])*\[Delta]*\[Nu])/(216*Sqrt[385]) - 
      (43*e^2*(-1 + E^((2*I)*\[ScriptL]))*\[Delta]*\[Nu])/(216*Sqrt[385]) - 
      (e^3*(-43 + 1300*E^((2*I)*\[ScriptL]) + 387*E^((4*I)*\[ScriptL]))*
        \[Delta]*\[Nu])/(1728*Sqrt[385]*E^(I*\[ScriptL])) - 
      (e^4*(-43 - 2466*E^((2*I)*\[ScriptL]) + 1821*E^((4*I)*\[ScriptL]) + 
         688*E^((6*I)*\[ScriptL]))*\[Delta]*\[Nu])/(2592*Sqrt[385]*
        E^((2*I)*\[ScriptL])) - (e^5*(-1161 - 9520*E^((2*I)*\[ScriptL]) + 
         112926*E^((4*I)*\[ScriptL]) + 60912*E^((6*I)*\[ScriptL]) + 
         26875*E^((8*I)*\[ScriptL]))*\[Delta]*\[Nu])/
       (82944*Sqrt[385]*E^((3*I)*\[ScriptL])) - 
      (e^6*(-688 - 3895*E^((2*I)*\[ScriptL]) - 118770*E^((4*I)*\[ScriptL]) + 
         60775*E^((6*I)*\[ScriptL]) + 41680*E^((8*I)*\[ScriptL]) + 
         20898*E^((10*I)*\[ScriptL]))*\[Delta]*\[Nu])/
       (51840*Sqrt[385]*E^((4*I)*\[ScriptL]))) + SO*x^(7/2)*\[Epsilon]^5*
     ((43*e*E^(I*\[ScriptL])*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(324*Sqrt[385]) + 
      (43*e^2*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(324*Sqrt[385]) + 
      (e^3*(-43 + 1472*E^((2*I)*\[ScriptL]) + 387*E^((4*I)*\[ScriptL]))*\[Nu]*
        ((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*\[Chi]S))/
       (2592*Sqrt[385]*E^(I*\[ScriptL])) + 
      (e^4*(-43 - 2724*E^((2*I)*\[ScriptL]) + 2079*E^((4*I)*\[ScriptL]) + 
         688*E^((6*I)*\[ScriptL]))*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(3888*Sqrt[385]*
        E^((2*I)*\[ScriptL])) + (e^5*(-1161 - 10552*E^((2*I)*\[ScriptL]) + 
         150318*E^((4*I)*\[ScriptL]) + 70200*E^((6*I)*\[ScriptL]) + 
         26875*E^((8*I)*\[ScriptL]))*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(124416*Sqrt[385]*
        E^((3*I)*\[ScriptL])) + (e^6*(-688 - 4325*E^((2*I)*\[ScriptL]) - 
         147300*E^((4*I)*\[ScriptL]) + 82855*E^((6*I)*\[ScriptL]) + 
         48560*E^((8*I)*\[ScriptL]) + 20898*E^((10*I)*\[ScriptL]))*\[Nu]*
        ((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*\[Chi]S))/
       (77760*Sqrt[385]*E^((4*I)*\[ScriptL]))) + 
    \[Epsilon]^6*(x^4*((-13*\[Delta]*\[Nu])/(63*Sqrt[385]) - 
        (e^2*\[Delta]*\[Nu]*(15587 + 56*E^((4*I)*\[ScriptL])*
            (-2355 + 5156*\[Nu]) - 7*E^((2*I)*\[ScriptL])*
            (-24963 + 41248*\[Nu])))/(39312*Sqrt[385]*E^((2*I)*\[ScriptL])) - 
        (e*\[Delta]*\[Nu]*(9789 + 7*E^((2*I)*\[ScriptL])*
            (-17254 + 41807*\[Nu])))/(39312*Sqrt[385]*E^(I*\[ScriptL])) - 
        (e^4*\[Delta]*\[Nu]*(449319 + E^((2*I)*\[ScriptL])*
            (611242 - 292649*\[Nu]) + E^((8*I)*\[ScriptL])*
            (-2648878 + 4529777*\[Nu]) + 3*E^((6*I)*\[ScriptL])*
            (-2701778 + 7717255*\[Nu]) - 3*E^((4*I)*\[ScriptL])*
            (-3994993 + 9129631*\[Nu])))/(471744*Sqrt[385]*
          E^((4*I)*\[ScriptL])) - (e^3*\[Delta]*\[Nu]*
          (194987 + E^((2*I)*\[ScriptL])*(392582 - 292649*\[Nu]) + 
           E^((6*I)*\[ScriptL])*(-1324070 + 2571233*\[Nu]) + 
           E^((4*I)*\[ScriptL])*(-5742651 + 16138108*\[Nu])))/
         (314496*Sqrt[385]*E^((3*I)*\[ScriptL])) - 
        (e^6*\[Delta]*\[Nu]*(40778088 + E^((4*I)*\[ScriptL])*
            (70260975 - 89235230*\[Nu]) + E^((2*I)*\[ScriptL])*
            (31691461 - 9364768*\[Nu]) - 840*E^((6*I)*\[ScriptL])*
            (-1839552 + 4585619*\[Nu]) + 3*E^((12*I)*\[ScriptL])*
            (-68188307 + 90318326*\[Nu]) + 10*E^((8*I)*\[ScriptL])*
            (-89386741 + 254900422*\[Nu]) + E^((10*I)*\[ScriptL])*
            (-375664593 + 1130560760*\[Nu])))/(18869760*Sqrt[385]*
          E^((6*I)*\[ScriptL])) - (e^5*\[Delta]*\[Nu]*(21778913 + 
           E^((4*I)*\[ScriptL])*(67666026 - 108284176*\[Nu]) + 
           E^((2*I)*\[ScriptL])*(22257906 - 7901523*\[Nu]) + 
           E^((10*I)*\[ScriptL])*(-116720742 + 175517881*\[Nu]) + 
           3*E^((8*I)*\[ScriptL])*(-91270681 + 265880048*\[Nu]) + 
           E^((6*I)*\[ScriptL])*(-798034940 + 2243839626*\[Nu])))/
         (15095808*Sqrt[385]*E^((5*I)*\[ScriptL]))) + 
      SO^2*x^4*((-43*e*E^(I*\[ScriptL])*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(3888*Sqrt[385]) - 
        (43*e^2*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
          (9*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(23 - 108*\[Nu] + 64*\[Nu]^2)*
            \[Chi]A*\[Chi]S + \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (23 - 92*\[Nu])*\[Chi]A^2 + (23 - 32*\[Nu] + 8*\[Nu]^2)*
              \[Chi]S^2)))/(3888*Sqrt[385]) - 
        (e^3*(-43 + 1644*E^((2*I)*\[ScriptL]) + 387*E^((4*I)*\[ScriptL]))*
          \[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(23 - 108*\[Nu] + 64*\[Nu]^2)*
            \[Chi]A*\[Chi]S + \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (23 - 92*\[Nu])*\[Chi]A^2 + (23 - 32*\[Nu] + 8*\[Nu]^2)*
              \[Chi]S^2)))/(31104*Sqrt[385]*E^(I*\[ScriptL])) - 
        (e^4*(-43 - 2982*E^((2*I)*\[ScriptL]) + 2337*E^((4*I)*\[ScriptL]) + 
           688*E^((6*I)*\[ScriptL]))*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(46656*Sqrt[385]*
          E^((2*I)*\[ScriptL])) - (e^5*(-1161 - 11584*E^((2*I)*\[ScriptL]) + 
           191838*E^((4*I)*\[ScriptL]) + 79488*E^((6*I)*\[ScriptL]) + 
           26875*E^((8*I)*\[ScriptL]))*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(1492992*Sqrt[385]*
          E^((3*I)*\[ScriptL])) - (e^6*(-688 - 4755*E^((2*I)*\[ScriptL]) - 
           178410*E^((4*I)*\[ScriptL]) + 107515*E^((6*I)*\[ScriptL]) + 
           55440*E^((8*I)*\[ScriptL]) + 20898*E^((10*I)*\[ScriptL]))*\[Nu]*
          (9*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(23 - 108*\[Nu] + 64*\[Nu]^2)*
            \[Chi]A*\[Chi]S + \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (23 - 92*\[Nu])*\[Chi]A^2 + (23 - 32*\[Nu] + 8*\[Nu]^2)*
              \[Chi]S^2)))/(933120*Sqrt[385]*E^((4*I)*\[ScriptL]))))
 
HoscMem[5, 3] = x^3*\[Epsilon]^4*
     ((-31*e^3*E^((3*I)*\[ScriptL])*\[Delta]*\[Nu])/(1296*Sqrt[330]) - 
      (31*e^4*E^((2*I)*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Delta]*
        \[Nu])/(432*Sqrt[330]) - (e^5*E^(I*\[ScriptL])*
        (651 - 1481*E^((2*I)*\[ScriptL]) + 1581*E^((4*I)*\[ScriptL]))*
        \[Delta]*\[Nu])/(10368*Sqrt[330]) - 
      (e^6*(-124 - 207*E^((2*I)*\[ScriptL]) - 2583*E^((4*I)*\[ScriptL]) + 
         2914*E^((6*I)*\[ScriptL]))*\[Delta]*\[Nu])/(10368*Sqrt[330])) + 
    SO*x^(7/2)*\[Epsilon]^5*((31*e^3*E^((3*I)*\[ScriptL])*\[Nu]*
        ((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*\[Chi]S))/
       (1944*Sqrt[330]) + (31*e^4*E^((2*I)*\[ScriptL])*
        (-1 + E^((2*I)*\[ScriptL]))*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(648*Sqrt[330]) + 
      (e^5*E^(I*\[ScriptL])*(651 - 1357*E^((2*I)*\[ScriptL]) + 
         1581*E^((4*I)*\[ScriptL]))*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(15552*Sqrt[330]) + 
      (e^6*(-124 - 579*E^((2*I)*\[ScriptL]) - 2211*E^((4*I)*\[ScriptL]) + 
         2914*E^((6*I)*\[ScriptL]))*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(15552*Sqrt[330])) + 
    \[Epsilon]^6*(x^4*(-(\[Delta]*\[Nu])/(189*Sqrt[330]) - 
        (e*(-123 + 134*E^((2*I)*\[ScriptL]))*\[Delta]*\[Nu])/
         (2016*Sqrt[330]*E^(I*\[ScriptL])) - 
        (e^2*(-1531 - 1215*E^((2*I)*\[ScriptL]) + 3290*E^((4*I)*\[ScriptL]))*
          \[Delta]*\[Nu])/(10080*Sqrt[330]*E^((2*I)*\[ScriptL])) - 
        (e^4*\[Delta]*\[Nu]*(-848172 - 406939*E^((2*I)*\[ScriptL]) - 
           348582*E^((4*I)*\[ScriptL]) + 50*E^((8*I)*\[ScriptL])*
            (-276347 + 485436*\[Nu]) - 5*E^((6*I)*\[ScriptL])*
            (-3050263 + 4854360*\[Nu])))/(1572480*Sqrt[330]*
          E^((4*I)*\[ScriptL])) - (e^3*\[Delta]*\[Nu]*(-2835027 - 
           1712178*E^((2*I)*\[ScriptL]) - 3948165*E^((4*I)*\[ScriptL]) + 
           10*E^((6*I)*\[ScriptL])*(-2515391 + 4865644*\[Nu])))/
         (9434880*Sqrt[330]*E^((3*I)*\[ScriptL])) - 
        (e^6*\[Delta]*\[Nu]*(-56661579 - 12008997*E^((2*I)*\[ScriptL]) - 
           16806426*E^((4*I)*\[ScriptL]) + E^((8*I)*\[ScriptL])*
            (728924394 - 1160358360*\[Nu]) + E^((6*I)*\[ScriptL])*
            (53182588 - 96297320*\[Nu]) - 12*E^((10*I)*\[ScriptL])*
            (-51887827 + 84664930*\[Nu]) + 8*E^((12*I)*\[ScriptL])*
            (-168374936 + 284079355*\[Nu])))/(37739520*Sqrt[330]*
          E^((6*I)*\[ScriptL])) - (e^5*\[Delta]*\[Nu]*(-138314397 - 
           48121398*E^((2*I)*\[ScriptL]) - 51649962*E^((4*I)*\[ScriptL]) + 
           30*E^((10*I)*\[ScriptL])*(-96411197 + 164709720*\[Nu]) + 
           12*E^((6*I)*\[ScriptL])*(-122612419 + 169394820*\[Nu]) - 
           5*E^((8*I)*\[ScriptL])*(-412828047 + 659176112*\[Nu])))/
         (150958080*Sqrt[330]*E^((5*I)*\[ScriptL]))) + 
      SO^2*x^4*((-31*e^3*E^((3*I)*\[ScriptL])*\[Nu]*
          (9*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(23 - 108*\[Nu] + 64*\[Nu]^2)*
            \[Chi]A*\[Chi]S + \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (23 - 92*\[Nu])*\[Chi]A^2 + (23 - 32*\[Nu] + 8*\[Nu]^2)*
              \[Chi]S^2)))/(23328*Sqrt[330]) - 
        (31*e^4*E^((2*I)*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
          (9*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(23 - 108*\[Nu] + 64*\[Nu]^2)*
            \[Chi]A*\[Chi]S + \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (23 - 92*\[Nu])*\[Chi]A^2 + (23 - 32*\[Nu] + 8*\[Nu]^2)*
              \[Chi]S^2)))/(7776*Sqrt[330]) - 
        (e^5*E^(I*\[ScriptL])*(217 - 411*E^((2*I)*\[ScriptL]) + 
           527*E^((4*I)*\[ScriptL]))*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(62208*Sqrt[330]) - 
        (e^6*(-124 - 951*E^((2*I)*\[ScriptL]) - 1839*E^((4*I)*\[ScriptL]) + 
           2914*E^((6*I)*\[ScriptL]))*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(186624*Sqrt[330])))
 
HoscMem[5, 5] = x^3*\[Epsilon]^4*
     ((5*e^5*E^((5*I)*\[ScriptL])*\[Delta]*\[Nu])/(1152*Sqrt[66]) + 
      (25*e^6*E^((4*I)*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Delta]*
        \[Nu])/(1152*Sqrt[66])) + SO*x^(7/2)*\[Epsilon]^5*
     ((-5*e^5*E^((5*I)*\[ScriptL])*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(1728*Sqrt[66]) - 
      (25*e^6*E^((4*I)*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
        ((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*\[Chi]S))/
       (1728*Sqrt[66])) + \[Epsilon]^6*
     (x^4*((3*Sqrt[3/22]*\[Delta]*\[Nu])/35 + 
        (e*(914 + 1773*E^((2*I)*\[ScriptL]))*\[Delta]*\[Nu])/
         (2016*Sqrt[66]*E^(I*\[ScriptL])) + 
        (e^2*(7870 + 8909*E^((2*I)*\[ScriptL]) + 21845*E^((4*I)*\[ScriptL]))*
          \[Delta]*\[Nu])/(10080*Sqrt[66]*E^((2*I)*\[ScriptL])) + 
        (e^3*(20724 + 19659*E^((2*I)*\[ScriptL]) + 
           19794*E^((4*I)*\[ScriptL]) + 74653*E^((6*I)*\[ScriptL]))*\[Delta]*
          \[Nu])/(16128*Sqrt[66]*E^((3*I)*\[ScriptL])) + 
        (e^4*(248125 + 191760*E^((2*I)*\[ScriptL]) + 
           213882*E^((4*I)*\[ScriptL]) + 121805*E^((6*I)*\[ScriptL]) + 
           1105260*E^((8*I)*\[ScriptL]))*\[Delta]*\[Nu])/
         (120960*Sqrt[66]*E^((4*I)*\[ScriptL])) - 
        (e^6*\[Delta]*\[Nu]*(-20655258 - 9236019*E^((2*I)*\[ScriptL]) - 
           11150165*E^((4*I)*\[ScriptL]) - 11854128*E^((6*I)*\[ScriptL]) - 
           14895270*E^((8*I)*\[ScriptL]) - 68*E^((10*I)*\[ScriptL])*
            (-451411 + 58240*\[Nu]) + 4*E^((12*I)*\[ScriptL])*
            (-32905349 + 990080*\[Nu])))/(4193280*Sqrt[66]*
          E^((6*I)*\[ScriptL])) - (e^5*\[Delta]*\[Nu]*(-161250310 - 
           97871085*E^((2*I)*\[ScriptL]) - 110328400*E^((4*I)*\[ScriptL]) - 
           121370210*E^((6*I)*\[ScriptL]) + 46703670*E^((8*I)*\[ScriptL]) + 
           E^((10*I)*\[ScriptL])*(-866689133 + 9395568*\[Nu])))/
         (50319360*Sqrt[66]*E^((5*I)*\[ScriptL]))) + 
      SO^2*x^4*((5*e^5*E^((5*I)*\[ScriptL])*\[Nu]*
          (9*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(23 - 108*\[Nu] + 64*\[Nu]^2)*
            \[Chi]A*\[Chi]S + \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (23 - 92*\[Nu])*\[Chi]A^2 + (23 - 32*\[Nu] + 8*\[Nu]^2)*
              \[Chi]S^2)))/(20736*Sqrt[66]) + (25*e^6*E^((4*I)*\[ScriptL])*
          (-1 + E^((2*I)*\[ScriptL]))*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(20736*Sqrt[66])))
 
HoscMem[6, 0] = 0
 
HoscMem[6, 2] = x^(7/2)*\[Epsilon]^5*
     (((-I/59136)*e^2*E^((2*I)*\[ScriptL])*\[Nu]*(-2783 + 8904*\[Nu]))/
       Sqrt[65] - ((I/29568)*e^3*E^(I*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*
        \[Nu]*(-2783 + 8904*\[Nu]))/Sqrt[65] - 
      ((I/2128896)*e^5*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
        (8349 - 26712*\[Nu] + 177*E^((4*I)*\[ScriptL])*(-2783 + 8904*\[Nu]) + 
         98*E^((2*I)*\[ScriptL])*(-7711 + 24888*\[Nu])))/
       (Sqrt[65]*E^(I*\[ScriptL])) - ((I/2128896)*e^4*\[Nu]*
        (27*(-2783 + 8904*\[Nu]) + 117*E^((4*I)*\[ScriptL])*
          (-2783 + 8904*\[Nu]) + E^((2*I)*\[ScriptL])*
          (-68926 + 231168*\[Nu])))/Sqrt[65] - 
      ((I/8515584)*e^6*\[Nu]*(E^((4*I)*\[ScriptL])*(1282418 - 
           4119864*\[Nu]) + 6*(-2783 + 8904*\[Nu]) + 
         1035*E^((8*I)*\[ScriptL])*(-2783 + 8904*\[Nu]) + 
         24*E^((2*I)*\[ScriptL])*(-60797 + 195846*\[Nu]) + 
         8*E^((6*I)*\[ScriptL])*(-222629 + 729582*\[Nu])))/
       (Sqrt[65]*E^((2*I)*\[ScriptL]))) + SO*x^4*\[Epsilon]^6*
     (((I/88704)*e^2*E^((2*I)*\[ScriptL])*\[Nu]*(-2783 + 8904*\[Nu])*
        (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/Sqrt[65] + 
      ((I/44352)*e^3*E^(I*\[ScriptL])*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
        (-2783 + 8904*\[Nu])*(-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/
       Sqrt[65] + ((I/3193344)*e^4*\[Nu]*(27*(-2783 + 8904*\[Nu]) + 
         117*E^((4*I)*\[ScriptL])*(-2783 + 8904*\[Nu]) + 
         20*E^((2*I)*\[ScriptL])*(-5951 + 19572*\[Nu]))*
        (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/Sqrt[65] + 
      ((I/3193344)*e^5*(-1 + E^((2*I)*\[ScriptL]))*\[Nu]*
        (8349 - 26712*\[Nu] + 177*E^((4*I)*\[ScriptL])*(-2783 + 8904*\[Nu]) + 
         E^((2*I)*\[ScriptL])*(-855866 + 2759568*\[Nu]))*
        (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/
       (Sqrt[65]*E^(I*\[ScriptL])) + ((I/12773376)*e^6*\[Nu]*
        (6*(-2783 + 8904*\[Nu]) + 1035*E^((8*I)*\[ScriptL])*
          (-2783 + 8904*\[Nu]) - 36*E^((4*I)*\[ScriptL])*
          (-27619 + 88242*\[Nu]) + 30*E^((2*I)*\[ScriptL])*
          (-53647 + 172704*\[Nu]) + E^((6*I)*\[ScriptL])*
          (-2432254 + 7920192*\[Nu]))*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu])*\[Chi]S))/(Sqrt[65]*E^((2*I)*\[ScriptL])))
 
HoscMem[6, 4] = x^(7/2)*\[Epsilon]^5*
     (((-I/101376)*e^4*E^((4*I)*\[ScriptL])*\[Nu]*(-703 + 2244*\[Nu]))/
       Sqrt[78] - ((I/25344)*e^5*E^((3*I)*\[ScriptL])*
        (-1 + E^((2*I)*\[ScriptL]))*\[Nu]*(-703 + 2244*\[Nu]))/Sqrt[78] - 
      ((I/202752)*e^6*E^((2*I)*\[ScriptL])*\[Nu]*
        (-7733 + E^((2*I)*\[ScriptL])*(17656 - 56388*\[Nu]) + 24684*\[Nu] + 
         21*E^((4*I)*\[ScriptL])*(-703 + 2244*\[Nu])))/Sqrt[78]) + 
    SO*x^4*\[Epsilon]^6*(((I/152064)*e^4*E^((4*I)*\[ScriptL])*\[Nu]*
        (-703 + 2244*\[Nu])*(-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/
       Sqrt[78] + ((I/38016)*e^5*E^((3*I)*\[ScriptL])*
        (-1 + E^((2*I)*\[ScriptL]))*\[Nu]*(-703 + 2244*\[Nu])*
        (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/Sqrt[78] + 
      ((I/304128)*e^6*E^((2*I)*\[ScriptL])*\[Nu]*
        (-7733 + E^((2*I)*\[ScriptL])*(16953 - 54144*\[Nu]) + 24684*\[Nu] + 
         21*E^((4*I)*\[ScriptL])*(-703 + 2244*\[Nu]))*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu])*\[Chi]S))/Sqrt[78])
 
HoscMem[6, 6] = (((25*I)/110592)*e^6*E^((6*I)*\[ScriptL])*x^(7/2)*
      \[Epsilon]^5*\[Nu]*(-19 + 64*\[Nu]))/Sqrt[143] - 
    (((25*I)/165888)*e^6*E^((6*I)*\[ScriptL])*SO*x^4*\[Epsilon]^6*\[Nu]*
      (-19 + 64*\[Nu])*(-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/Sqrt[143]
 
HoscMem[7, 1] = x^4*\[Epsilon]^6*
    ((-5*e*E^(I*\[ScriptL])*\[Delta]*\[Nu]*(-5023 + 16296*\[Nu]))/
      (5189184*Sqrt[2]) - (5*e^2*(-1 + E^((2*I)*\[ScriptL]))*\[Delta]*\[Nu]*
       (-5023 + 16296*\[Nu]))/(5189184*Sqrt[2]) - 
     (5*e^4*(-1 + E^((2*I)*\[ScriptL]))*\[Delta]*\[Nu]*
       (-5023 + 16296*\[Nu] + 16*E^((4*I)*\[ScriptL])*(-5023 + 16296*\[Nu]) + 
        70*E^((2*I)*\[ScriptL])*(-8143 + 24972*\[Nu])))/
      (62270208*Sqrt[2]*E^((2*I)*\[ScriptL])) - 
     (5*e^3*\[Delta]*\[Nu]*(5023 - 16296*\[Nu] + 9*E^((4*I)*\[ScriptL])*
         (-5023 + 16296*\[Nu]) + 18*E^((2*I)*\[ScriptL])*
         (-18693 + 56896*\[Nu])))/(41513472*Sqrt[2]*E^(I*\[ScriptL])) - 
     (e^6*(-1 + E^((2*I)*\[ScriptL]))*\[Delta]*\[Nu]*
       (8*(-5023 + 16296*\[Nu]) + 243*E^((8*I)*\[ScriptL])*
         (-5023 + 16296*\[Nu]) + 7*E^((2*I)*\[ScriptL])*
         (-71207 + 218964*\[Nu]) + 28*E^((4*I)*\[ScriptL])*
         (-779033 + 2338491*\[Nu]) + 3*E^((6*I)*\[ScriptL])*
         (-2449103 + 7495656*\[Nu])))/(124540416*Sqrt[2]*
       E^((4*I)*\[ScriptL])) - (5*e^5*\[Delta]*\[Nu]*
       (135621 + E^((2*I)*\[ScriptL])*(2219764 - 6796608*\[Nu]) - 
        439992*\[Nu] + 625*E^((8*I)*\[ScriptL])*(-5023 + 16296*\[Nu]) + 
        108*E^((6*I)*\[ScriptL])*(-158191 + 479472*\[Nu]) + 
        18*E^((4*I)*\[ScriptL])*(-2796539 + 8327368*\[Nu])))/
      (1992646656*Sqrt[2]*E^((3*I)*\[ScriptL])))
 
HoscMem[7, 3] = x^4*\[Epsilon]^6*
    ((-5*e^3*E^((3*I)*\[ScriptL])*\[Delta]*\[Nu]*(-12539 + 28896*\[Nu]))/
      (6918912*Sqrt[6]) - (5*e^4*E^((2*I)*\[ScriptL])*
       (-1 + E^((2*I)*\[ScriptL]))*\[Delta]*\[Nu]*(-12539 + 28896*\[Nu]))/
      (2306304*Sqrt[6]) - (5*e^6*(-1 + E^((2*I)*\[ScriptL]))*\[Delta]*\[Nu]*
       (-25078 + 57792*\[Nu] + 70*E^((2*I)*\[ScriptL])*
         (-4915 + 11532*\[Nu]) + 47*E^((4*I)*\[ScriptL])*
         (-12539 + 28896*\[Nu])))/(27675648*Sqrt[6]) - 
     (5*e^5*E^(I*\[ScriptL])*\[Delta]*\[Nu]*(-263319 + E^((2*I)*\[ScriptL])*
         (414302 - 945168*\[Nu]) + 606816*\[Nu] + 51*E^((4*I)*\[ScriptL])*
         (-12539 + 28896*\[Nu])))/(55351296*Sqrt[6]))
 
HoscMem[7, 5] = x^4*\[Epsilon]^6*
    (-(e^5*E^((5*I)*\[ScriptL])*\[Delta]*\[Nu]*(-6547 + 14664*\[Nu]))/
      (1797120*Sqrt[66]) - (e^6*E^((4*I)*\[ScriptL])*
       (-1 + E^((2*I)*\[ScriptL]))*\[Delta]*\[Nu]*(-6547 + 14664*\[Nu]))/
      (359424*Sqrt[66]))
 
HoscMem[7, 7] = 0


(* ::Section::Closed:: *)
(*DC memory contributions*)


(* ::Text:: *)
(*The DC memory contributions are expanded to O(e^6), and expressed in terms of x, e and an initial eccentricity ei.*)
(**)
(*We first include the integrand \dot{h}_{lm}^{DCmem}, which is given by Eq. (2.14) in the paper, then the provide the integration result for h_{lm}^{DCmem}, computed using*)
(*h_{lm}^{DCmem} = \int_{e_0}^{e(T_R)} de \dot{h}_{lm}^{DCmem} / \dot{e}*)
(*where \dot{e} is computed from the fluxes at infinity. Therefore, our result for the DC memory does not account for the horizon fluxes, but the integrand is valid without needing to account for the horizon absorption.*)
(**)
(*For convenience, we also provide separately the circular-orbit limit of the eccentric-orbit results, since the 3PN DC memory was the only contribution not previously known at that order for circular orbits and aligned spins (it was derived at 2PN in [Mitman et. al. arXiv:2208.04356]).*)


(* ::Subsection::Closed:: *)
(*DC memory integrand*)


hdotDCMem[2, 0] = x^5*((-128*Sqrt[(2*Pi)/15]*\[Nu]^2)/7 - 
      (2504*e^2*Sqrt[(2*Pi)/15]*\[Nu]^2)/21 - 
      (482*e^4*Sqrt[(10*Pi)/3]*\[Nu]^2)/7 - (434*e^6*Sqrt[(10*Pi)/3]*\[Nu]^2)/
       3) + x^6*\[Epsilon]^2*((-4*Sqrt[(2*Pi)/15]*\[Nu]^2*(-1219 + 12*\[Nu]))/
       63 + (4*e^2*Sqrt[(2*Pi)/15]*\[Nu]^2*(2052 + 2365*\[Nu]))/63 + 
      (e^4*Sqrt[(2*Pi)/15]*\[Nu]^2*(-9539 + 18050*\[Nu]))/21 + 
      (e^6*Sqrt[(2*Pi)/15]*\[Nu]^2*(-162681 + 169505*\[Nu]))/63) + 
    SO*x^(15/2)*\[Epsilon]^5*((2*Sqrt[(2*Pi)/15]*\[Nu]^2*
        (\[Delta]*(1028 - 3145*\[Nu])*\[Chi]A + 
         (1028 + 2855*\[Nu] + 4004*\[Nu]^2)*\[Chi]S))/63 + 
      (e^2*Sqrt[(2*Pi)/15]*\[Nu]^2*(3*\[Delta]*(113109 - 70606*\[Nu])*
          \[Chi]A + (339327 - 609162*\[Nu] + 253288*\[Nu]^2)*\[Chi]S))/63 + 
      (e^4*Sqrt[Pi/30]*\[Nu]^2*(\[Delta]*(11037783 - 6245000*\[Nu])*\[Chi]A + 
         (11037783 - 19642280*\[Nu] + 7354248*\[Nu]^2)*\[Chi]S))/126 + 
      (e^6*Sqrt[Pi/30]*\[Nu]^2*(\[Delta]*(93493801 - 51506762*\[Nu])*
          \[Chi]A + (93493801 - 163320138*\[Nu] + 60123112*\[Nu]^2)*\[Chi]S))/
       252) + \[Epsilon]^3*(x^(13/2)*((-512*Sqrt[2/15]*Pi^(3/2)*\[Nu]^2)/7 - 
        (18664*Sqrt[2/15]*e^2*Pi^(3/2)*\[Nu]^2)/21 - 
        (85786*Sqrt[2/15]*e^4*Pi^(3/2)*\[Nu]^2)/21 - 
        (2064521*e^6*Pi^(3/2)*\[Nu]^2)/(84*Sqrt[30])) + 
      SO*x^(13/2)*((e^6*Sqrt[Pi/30]*\[Nu]^2*(1098781*\[Delta]*\[Chi]A + 
           (1098781 - 1152476*\[Nu])*\[Chi]S))/36 + 
        (e^4*Sqrt[(5*Pi)/6]*\[Nu]^2*(81833*\[Delta]*\[Chi]A + 
           (81833 - 91476*\[Nu])*\[Chi]S))/42 + 
        (e^2*Sqrt[(2*Pi)/15]*\[Nu]^2*(60433*\[Delta]*\[Chi]A + 
           (60433 - 77468*\[Nu])*\[Chi]S))/63 + 
        (2*Sqrt[(2*Pi)/15]*\[Nu]^2*(505*\[Delta]*\[Chi]A + 
           (505 - 1124*\[Nu])*\[Chi]S))/21)) + 
    \[Epsilon]^4*(x^7*((2*Sqrt[(2*Pi)/15]*\[Nu]^2*(25376 + 126207*\[Nu] + 
           151236*\[Nu]^2))/6237 + (e^2*Sqrt[(2*Pi)/15]*\[Nu]^2*
          (2888258 + 949671*\[Nu] + 1752368*\[Nu]^2))/2079 + 
        (e^4*Sqrt[Pi/30]*\[Nu]^2*(22823810 + 9087473*\[Nu] + 
           12944500*\[Nu]^2))/1386 + (e^6*Sqrt[Pi/30]*\[Nu]^2*
          (467575078 + 223997529*\[Nu] + 266898380*\[Nu]^2))/8316) + 
      SO^2*x^7*((-4*Sqrt[(2*Pi)/15]*\[Nu]^2*(64*\[Delta]*\[Kappa]A + 
           64*\[Kappa]S - 128*\[Kappa]S*\[Nu] + 63*\[Chi]A^2 - 
           256*\[Nu]*\[Chi]A^2 + 126*\[Delta]*\[Chi]A*\[Chi]S + 
           63*\[Chi]S^2 + 4*\[Nu]*\[Chi]S^2))/7 - 
        (4*e^2*Sqrt[(6*Pi)/5]*\[Nu]^2*(388*\[Delta]*\[Kappa]A + 
           388*\[Kappa]S - 776*\[Kappa]S*\[Nu] + 383*\[Chi]A^2 - 
           1552*\[Nu]*\[Chi]A^2 + 766*\[Delta]*\[Chi]A*\[Chi]S + 
           383*\[Chi]S^2 + 20*\[Nu]*\[Chi]S^2))/7 - 
        (e^4*Sqrt[(10*Pi)/3]*\[Nu]^2*(-15070*\[Kappa]S*(-1 + 2*\[Nu]) + 
           14881*\[Chi]A^2 - 60280*\[Nu]*\[Chi]A^2 + 14881*\[Chi]S^2 + 
           756*\[Nu]*\[Chi]S^2 + 2*\[Delta]*(7535*\[Kappa]A + 
             14881*\[Chi]A*\[Chi]S)))/21 - (2*e^6*Sqrt[(2*Pi)/15]*\[Nu]^2*
          (-126508*\[Kappa]S*(-1 + 2*\[Nu]) + 124933*\[Chi]A^2 - 
           506032*\[Nu]*\[Chi]A^2 + 124933*\[Chi]S^2 + 6300*\[Nu]*\[Chi]S^2 + 
           2*\[Delta]*(63254*\[Kappa]A + 124933*\[Chi]A*\[Chi]S)))/21)) + 
    \[Epsilon]^6*(SO^2*x^8*((2*Sqrt[(2*Pi)/15]*\[Nu]^2*
          (4*\[Delta]*\[Kappa]A*(65 + 251*\[Nu]) + \[Kappa]S*
            (260 + 484*\[Nu] - 720*\[Nu]^2) + 103*\[Chi]A^2 + 
           100*\[Nu]*\[Chi]A^2 - 1440*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(103 + 3474*\[Nu])*\[Chi]A*\[Chi]S + 103*\[Chi]S^2 + 
           6436*\[Nu]*\[Chi]S^2 - 3840*\[Nu]^2*\[Chi]S^2))/21 - 
        (e^2*Sqrt[(2*Pi)/15]*\[Nu]^2*(\[Delta]*\[Kappa]A*
            (746268 - 742260*\[Nu]) + 12*\[Kappa]S*(62189 - 186233*\[Nu] + 
             54082*\[Nu]^2) + 1566985*\[Chi]A^2 - 6665938*\[Nu]*\[Chi]A^2 + 
           1297968*\[Nu]^2*\[Chi]A^2 + 70*\[Delta]*(44771 - 67942*\[Nu])*
            \[Chi]A*\[Chi]S + 1566985*\[Chi]S^2 - 4357942*\[Nu]*\[Chi]S^2 + 
           1530424*\[Nu]^2*\[Chi]S^2))/189 - (e^4*Sqrt[Pi/30]*\[Nu]^2*
          (\[Delta]*\[Kappa]A*(9159276 - 7584816*\[Nu]) + 
           12*\[Kappa]S*(763273 - 2158614*\[Nu] + 634260*\[Nu]^2) + 
           18118659*\[Chi]A^2 - 76831650*\[Nu]*\[Chi]A^2 + 
           15222240*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(18118659 - 
             22912582*\[Nu])*\[Chi]A*\[Chi]S + 18118659*\[Chi]S^2 - 
           41468150*\[Nu]*\[Chi]S^2 + 12867160*\[Nu]^2*\[Chi]S^2))/126 - 
        (e^6*Sqrt[Pi/30]*\[Nu]^2*(12*\[Kappa]S*(21207839 - 58794087*\[Nu] + 
             17739582*\[Nu]^2) + 488810695*\[Chi]A^2 - 2072916694*\[Nu]*
            \[Chi]A^2 + 425749968*\[Nu]^2*\[Chi]A^2 + 488810695*\[Chi]S^2 - 
           1024283842*\[Nu]*\[Chi]S^2 + 297419368*\[Nu]^2*\[Chi]S^2 - 
           2*\[Delta]*(6*\[Kappa]A*(-21207839 + 16378409*\[Nu]) + 
             (-488810695 + 570978878*\[Nu])*\[Chi]A*\[Chi]S)))/756) + 
      SO*x^8*((4*Sqrt[2/15]*Pi^(3/2)*\[Nu]^2*(1017*\[Delta]*\[Chi]A + 
           (1017 - 2276*\[Nu])*\[Chi]S))/21 - (e^2*Sqrt[(2*Pi)/15]*\[Nu]^2*
          (Pi*(-406274*\[Delta]*\[Chi]A - 406274*\[Chi]S + 522640*\[Nu]*
              \[Chi]S) - (3969*I)*(32*\[Delta]*\[Chi]A - 15*\[Nu]*\[Chi]S)*
            (2*ArcCoth[5] - Log[3/2])))/63 + (e^4*Sqrt[Pi/30]*\[Nu]^2*
          (24778459*Pi*\[Delta]*\[Chi]A + Pi*(24778459 - 27458780*\[Nu])*
            \[Chi]S + (65691*I)*(32*\[Delta]*\[Chi]A - 15*\[Nu]*\[Chi]S)*
            (2*ArcCoth[5] - Log[3/2])))/252 + 
        (e^6*Sqrt[Pi/30]*\[Nu]^2*(4*Pi*(642138509*\[Delta]*\[Chi]A + 
             (642138509 - 663392766*\[Nu])*\[Chi]S) + (2725569*I)*
            (32*\[Delta]*\[Chi]A - 15*\[Nu]*\[Chi]S)*(2*ArcCoth[5] - 
             Log[3/2])))/6048))
 
hdotDCMem[4, 0] = x^5*((-163*e^4*Sqrt[Pi/10]*\[Nu]^2)/42 - 
      (49*e^6*Sqrt[Pi/10]*\[Nu]^2)/6 - (32*Sqrt[(2*Pi)/5]*\[Nu]^2)/315 - 
      (211*e^2*Sqrt[(2*Pi)/5]*\[Nu]^2)/315) + x^6*\[Epsilon]^2*
     ((e^6*Sqrt[Pi/10]*(4745247 - 21440090*\[Nu])*\[Nu]^2)/13860 + 
      (e^4*Sqrt[Pi/10]*(1716039 - 7562150*\[Nu])*\[Nu]^2)/13860 + 
      (Sqrt[Pi/10]*(30399 - 103100*\[Nu])*\[Nu]^2)/10395 - 
      (e^2*Sqrt[(2*Pi)/5]*\[Nu]^2*(-155049 + 645119*\[Nu]))/10395) + 
    SO*x^(15/2)*\[Epsilon]^5*
     ((e^6*Sqrt[Pi/10]*\[Nu]^2*(3*\[Delta]*(-61221299 + 347563626*\[Nu])*
          \[Chi]A + (-183663897 + 1291942544*\[Nu] - 961013102*\[Nu]^2)*
          \[Chi]S))/41580 + (e^4*Sqrt[Pi/10]*\[Nu]^2*
        (\[Delta]*(-45047937 + 257928856*\[Nu])*\[Chi]A + 
         (-45047937 + 339328444*\[Nu] - 267068588*\[Nu]^2)*\[Chi]S))/41580 + 
      (e^2*Sqrt[(2*Pi)/5]*\[Nu]^2*(13*\[Delta]*(-62478 + 354127*\[Nu])*
          \[Chi]A + (-812214 + 6894589*\[Nu] - 5876556*\[Nu]^2)*\[Chi]S))/
       10395 + (Sqrt[Pi/10]*\[Nu]^2*(\[Delta]*(-80247 + 368188*\[Nu])*
          \[Chi]A - 3*(26749 - 275028*\[Nu] + 253568*\[Nu]^2)*\[Chi]S))/
       10395) + \[Epsilon]^3*
     (x^(13/2)*((-128*Sqrt[2/5]*Pi^(3/2)*\[Nu]^2)/315 - 
        (104*Sqrt[2/5]*e^2*Pi^(3/2)*\[Nu]^2)/21 - 
        (479*Sqrt[2/5]*e^4*Pi^(3/2)*\[Nu]^2)/21 - 
        (6232357*e^6*Pi^(3/2)*\[Nu]^2)/(45360*Sqrt[10])) + 
      SO*x^(13/2)*((e^6*Sqrt[Pi/10]*\[Nu]^2*(91345*\[Delta]*\[Chi]A + 
           (91345 - 534632*\[Nu])*\[Chi]S))/540 + 
        (e^4*Sqrt[Pi/10]*\[Nu]^2*(33589*\[Delta]*\[Chi]A + 
           (33589 - 206192*\[Nu])*\[Chi]S))/630 + 
        (e^2*Sqrt[(2*Pi)/5]*\[Nu]^2*(4829*\[Delta]*\[Chi]A + 
           (4829 - 33268*\[Nu])*\[Chi]S))/945 + 
        (2*Sqrt[(2*Pi)/5]*\[Nu]^2*(107*\[Delta]*\[Chi]A + (107 - 1324*\[Nu])*
            \[Chi]S))/945)) + \[Epsilon]^4*
     (x^7*((Sqrt[(2*Pi)/5]*\[Nu]^2*(-2902797 + 12988674*\[Nu] + 
           1902920*\[Nu]^2))/405405 + (e^2*Sqrt[Pi/10]*\[Nu]^2*
          (-95981238 + 368585403*\[Nu] + 319367500*\[Nu]^2))/810810 + 
        (e^4*Sqrt[Pi/10]*\[Nu]^2*(-371976636 + 977010747*\[Nu] + 
           3205659440*\[Nu]^2))/1081080 + (e^6*Sqrt[Pi/10]*\[Nu]^2*
          (-606923607 - 132850998*\[Nu] + 12887095970*\[Nu]^2))/1081080) + 
      SO^2*x^7*((-8*Sqrt[(2*Pi)/5]*\[Nu]^2*(8*\[Delta]*\[Kappa]A + 
           8*\[Kappa]S - 16*\[Kappa]S*\[Nu] + 7*\[Chi]A^2 - 
           32*\[Nu]*\[Chi]A^2 + 14*\[Delta]*\[Chi]A*\[Chi]S + 7*\[Chi]S^2 + 
           4*\[Nu]*\[Chi]S^2))/315 - (2*e^2*Sqrt[(2*Pi)/5]*\[Nu]^2*
          (115*\[Delta]*\[Kappa]A + 115*\[Kappa]S - 230*\[Kappa]S*\[Nu] + 
           103*\[Chi]A^2 - 460*\[Nu]*\[Chi]A^2 + 206*\[Delta]*\[Chi]A*
            \[Chi]S + 103*\[Chi]S^2 + 48*\[Nu]*\[Chi]S^2))/63 - 
        (e^6*Sqrt[(2*Pi)/5]*\[Nu]^2*(-20779*\[Kappa]S*(-1 + 2*\[Nu]) + 
           18679*\[Chi]A^2 - 83116*\[Nu]*\[Chi]A^2 + 18679*\[Chi]S^2 + 
           8400*\[Nu]*\[Chi]S^2 + \[Delta]*(20779*\[Kappa]A + 
             37358*\[Chi]A*\[Chi]S)))/315 + (e^4*Sqrt[Pi/10]*\[Nu]^2*
          (24769*\[Kappa]S*(-1 + 2*\[Nu]) - 22249*\[Chi]A^2 + 
           99076*\[Nu]*\[Chi]A^2 - 22249*\[Chi]S^2 - 10080*\[Nu]*\[Chi]S^2 - 
           \[Delta]*(24769*\[Kappa]A + 44498*\[Chi]A*\[Chi]S)))/630)) + 
    \[Epsilon]^6*(SO^2*x^8*((e^6*Sqrt[Pi/10]*\[Nu]^2*
          (-3*\[Delta]*\[Kappa]A*(-516496079 + 4102554728*\[Nu]) + 
           3*\[Kappa]S*(516496079 - 5135546886*\[Nu] + 7424880960*\[Nu]^2) + 
           803511137*\[Chi]A^2 - 14529852416*\[Nu]*\[Chi]A^2 + 
           44549285760*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(803511137 - 
             6469176232*\[Nu])*\[Chi]A*\[Chi]S + 803511137*\[Chi]S^2 - 
           1622544596*\[Nu]*\[Chi]S^2 - 7277631808*\[Nu]^2*\[Chi]S^2))/
         498960 + (e^4*Sqrt[Pi/10]*\[Nu]^2*(\[Delta]*\[Kappa]A*
            (61720029 - 480398094*\[Nu]) + \[Kappa]S*(61720029 - 
             603838152*\[Nu] + 868999620*\[Nu]^2) + 33710641*\[Chi]A^2 - 
           576769060*\[Nu]*\[Chi]A^2 + 1737999240*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(33710641 - 251488160*\[Nu])*\[Chi]A*\[Chi]S + 
           33710641*\[Chi]S^2 - 61049824*\[Nu]*\[Chi]S^2 - 
           304572368*\[Nu]^2*\[Chi]S^2))/83160 + 
        (e^2*Sqrt[Pi/10]*\[Nu]^2*(\[Delta]*\[Kappa]A*(6793239 - 
             49516860*\[Nu]) + \[Kappa]S*(6793239 - 63103338*\[Nu] + 
             89237736*\[Nu]^2) + 4324355*\[Chi]A^2 - 62806232*\[Nu]*
            \[Chi]A^2 + 178475472*\[Nu]^2*\[Chi]A^2 + 14*\[Delta]*
            (617765 - 3745432*\[Nu])*\[Chi]A*\[Chi]S + 4324355*\[Chi]S^2 - 
           6927236*\[Nu]*\[Chi]S^2 - 35254912*\[Nu]^2*\[Chi]S^2))/62370 + 
        (Sqrt[Pi/10]*\[Nu]^2*(-3*\[Delta]*\[Kappa]A*(-66049 + 349372*\[Nu]) + 
           3*\[Kappa]S*(66049 - 481470*\[Nu] + 610152*\[Nu]^2) + 
           199511*\[Chi]A^2 - 1743200*\[Nu]*\[Chi]A^2 + 3660912*\[Nu]^2*
            \[Chi]A^2 + 2*\[Delta]*(199511 - 651676*\[Nu])*\[Chi]A*\[Chi]S + 
           199511*\[Chi]S^2 - 358196*\[Nu]*\[Chi]S^2 - 919456*\[Nu]^2*
            \[Chi]S^2))/31185) + 
      SO*x^8*((4*Sqrt[2/5]*Pi^(3/2)*\[Nu]^2*(235*\[Delta]*\[Chi]A + 
           (235 - 2732*\[Nu])*\[Chi]S))/945 - (2*e^2*Sqrt[(2*Pi)/5]*\[Nu]^2*
          (-5527*Pi*\[Delta]*\[Chi]A + Pi*(-5527 + 37824*\[Nu])*\[Chi]S - 
           (441*I)*(4*\[Delta]*\[Chi]A - 15*\[Nu]*\[Chi]S)*
            (2*ArcCoth[5] - Log[3/2])))/315 + 
        (e^4*Sqrt[Pi/10]*\[Nu]^2*(685369*Pi*\[Delta]*\[Chi]A + 
           Pi*(685369 - 4204276*\[Nu])*\[Chi]S + (14598*I)*
            (4*\[Delta]*\[Chi]A - 15*\[Nu]*\[Chi]S)*(2*ArcCoth[5] - 
             Log[3/2])))/1260 + (e^6*Sqrt[Pi/10]*\[Nu]^2*
          (322536290*Pi*\[Delta]*\[Chi]A + 2*Pi*(161268145 - 948224699*\[Nu])*
            \[Chi]S + (2725569*I)*(4*\[Delta]*\[Chi]A - 15*\[Nu]*\[Chi]S)*
            (2*ArcCoth[5] - Log[3/2])))/136080))
 
hdotDCMem[6, 0] = x^6*\[Epsilon]^2*
     ((e^4*Sqrt[(5*Pi)/273]*(16109 - 62412*\[Nu])*\[Nu]^2)/1848 + 
      (e^6*Sqrt[(5*Pi)/273]*(6757 - 25886*\[Nu])*\[Nu]^2)/264 + 
      (e^2*Sqrt[Pi/1365]*(2131 - 8463*\[Nu])*\[Nu]^2)/231 + 
      (Sqrt[Pi/1365]*(839 - 3612*\[Nu])*\[Nu]^2)/1386) + 
    x^7*\[Epsilon]^4*(-(Sqrt[Pi/1365]*\[Nu]^2*(982361 - 5074830*\[Nu] + 
          5601960*\[Nu]^2))/124740 - (e^2*Sqrt[Pi/1365]*\[Nu]^2*
        (2162947 - 11302362*\[Nu] + 14504616*\[Nu]^2))/16632 - 
      (e^4*Sqrt[Pi/1365]*\[Nu]^2*(16548751 - 88343430*\[Nu] + 
         120961008*\[Nu]^2))/22176 - (e^6*Sqrt[Pi/1365]*\[Nu]^2*
        (353520905 - 1913159178*\[Nu] + 2701311984*\[Nu]^2))/133056) + 
    SO*x^(15/2)*\[Epsilon]^5*((e^6*Sqrt[Pi/1365]*\[Nu]^2*
        (\[Delta]*(-4735519 + 2058105*\[Nu])*\[Chi]A + 
         (-4735519 + 41828155*\[Nu] - 89904192*\[Nu]^2)*\[Chi]S))/2772 + 
      (e^4*Sqrt[Pi/1365]*\[Nu]^2*(3*\[Delta]*(-253495 + 70868*\[Nu])*
          \[Chi]A + (-760485 + 6907612*\[Nu] - 15267000*\[Nu]^2)*\[Chi]S))/
       1848 - (Sqrt[Pi/1365]*\[Nu]^2*(\[Delta]*(481 + 456*\[Nu])*\[Chi]A + 
         (481 - 4824*\[Nu] + 12624*\[Nu]^2)*\[Chi]S))/198 - 
      (e^2*Sqrt[Pi/1365]*\[Nu]^2*(\[Delta]*(13115 + 798*\[Nu])*\[Chi]A + 
         (13115 - 125198*\[Nu] + 291816*\[Nu]^2)*\[Chi]S))/231) + 
    SO^2*x^8*\[Epsilon]^6*((Sqrt[Pi/1365]*\[Nu]^2*
        (\[Kappa]S*(835 - 5786*\[Nu] + 7224*\[Nu]^2) + 835*\[Chi]A^2 - 
         7456*\[Nu]*\[Chi]A^2 + 14448*\[Nu]^2*\[Chi]A^2 + 835*\[Chi]S^2 - 
         4116*\[Nu]*\[Chi]S^2 - 2016*\[Nu]^2*\[Chi]S^2 - 
         \[Delta]*(-835 + 4116*\[Nu])*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/462 + 
      (e^2*Sqrt[Pi/1365]*\[Nu]^2*(\[Kappa]S*(2063 - 13786*\[Nu] + 
           17208*\[Nu]^2) + 2063*\[Chi]A^2 - 17912*\[Nu]*\[Chi]A^2 + 
         34416*\[Nu]^2*\[Chi]A^2 + 2063*\[Chi]S^2 - 9660*\[Nu]*\[Chi]S^2 - 
         4224*\[Nu]^2*\[Chi]S^2 - \[Delta]*(-2063 + 9660*\[Nu])*
          (\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/44 + 
      (e^6*Sqrt[(5*Pi)/273]*\[Nu]^2*(\[Kappa]S*(2360501 - 15223858*\[Nu] + 
           18840528*\[Nu]^2) + 2360501*\[Chi]A^2 - 19944860*\[Nu]*\[Chi]A^2 + 
         37681056*\[Nu]^2*\[Chi]A^2 + 2360501*\[Chi]S^2 - 
         10502856*\[Nu]*\[Chi]S^2 - 4330368*\[Nu]^2*\[Chi]S^2 - 
         13*\[Delta]*(-181577 + 807912*\[Nu])*(\[Kappa]A + 
           2*\[Chi]A*\[Chi]S)))/7392 + (e^4*Sqrt[Pi/1365]*\[Nu]^2*
        (\[Kappa]S*(1336351 - 8728598*\[Nu] + 10839696*\[Nu]^2) + 
         1336351*\[Chi]A^2 - 11401300*\[Nu]*\[Chi]A^2 + 
         21679392*\[Nu]^2*\[Chi]A^2 + 1336351*\[Chi]S^2 - 
         6055896*\[Nu]*\[Chi]S^2 - 2544192*\[Nu]^2*\[Chi]S^2 - 
         \[Delta]*(-1336351 + 6055896*\[Nu])*(\[Kappa]A + 
           2*\[Chi]A*\[Chi]S)))/3696)
 
hdotDCMem[8, 0] = x^7*\[Epsilon]^4*
    (-(Sqrt[Pi/595]*\[Nu]^2*(75601 - 452070*\[Nu] + 733320*\[Nu]^2))/694980 - 
     (e^2*Sqrt[Pi/595]*\[Nu]^2*(241751 - 1461492*\[Nu] + 2269512*\[Nu]^2))/
      92664 - (e^4*Sqrt[Pi/595]*\[Nu]^2*(2266811 - 13678584*\[Nu] + 
        20742624*\[Nu]^2))/123552 - (e^6*Sqrt[(7*Pi)/85]*\[Nu]^2*
       (7716547 - 46470636*\[Nu] + 69582096*\[Nu]^2))/741312)


(* ::Subsection::Closed:: *)
(*DC memory for eccentric orbits*)


HlmDCMem[2, 0] = ((5*(e^(12/19) - ei^(12/19)))/(14*Sqrt[6]*ei^(12/19)) + 
      (5*(43199*e^2 - 361*e^(26/19)*ei^(12/19) - 42838*ei^2))/
       (525616*Sqrt[6]*(ei/e)^(12/19)) + 
      (5*(7855261648*e^4 - 68679167*e^(64/19)*ei^(12/19) - 
         9110443136*e^2*ei^2 + 1323860655*ei^4))/(97150656512*Sqrt[6]*
        (ei/e)^(12/19)))*x + x^2*\[Epsilon]^2*
     ((5*(e^(12/19) - ei^(12/19))*(19*ei^(12/19)*(-4075 + 5628*\[Nu]) + 
         e^(12/19)*(-145417 + 239316*\[Nu])))/(1072512*Sqrt[6]*ei^(24/19)) + 
      (5*e^(12/19)*(42*e^2*ei^(12/19)*(1082158253 - 1430751140*\[Nu]) + 
         3598392*ei^(50/19)*(-2833 + 5516*\[Nu]) + 302393*e^(50/19)*
          (-145417 + 239316*\[Nu]) + 109744*e^(26/19)*ei^(24/19)*
          (-22585 + 1216992*\[Nu]) - 9*e^(12/19)*ei^2*(-1243916559 + 
           18409137292*\[Nu])))/(140932366848*Sqrt[6]*ei^(24/19)) + 
      (5*e^(12/19)*(-120471319605*ei^(88/19)*(-2833 + 5516*\[Nu]) + 
         94403934716*e^(88/19)*(-145417 + 239316*\[Nu]) + 
         9595712*e^2*ei^(50/19)*(-1082158253 + 1430751140*\[Nu]) - 
         4147104*e^(50/19)*ei^2*(-1243916559 + 18409137292*\[Nu]) + 
         2476099*e^(64/19)*ei^(24/19)*(-2933732575 + 27297870964*\[Nu]) - 
         448*e^(12/19)*ei^4*(-23500911051568 + 15843764110629*\[Nu]) - 
         364*e^4*ei^(12/19)*(-42164487546687 + 54413587634588*\[Nu])))/
       (28219545498353664*Sqrt[6]*ei^(24/19))) + SO*x^(7/2)*\[Epsilon]^5*
     ((560*e^(24/19)*ei^(18/19)*(-145417 + 239316*\[Nu])*
         (-157*\[Delta]*\[Chi]A + (-157 + 110*\[Nu])*\[Chi]S) + 
        105*e^(30/19)*ei^(12/19)*(-2833 + 5516*\[Nu])*
         (-99391*\[Delta]*\[Chi]A + (-99391 + 22844*\[Nu])*\[Chi]S) + 
        2*e^(42/19)*(\[Delta]*(7875278603 + 32255910204*\[Nu])*\[Chi]A + 
          (7875278603 + 99877088*\[Nu] - 7901609520*\[Nu]^2)*\[Chi]S) + 
        1344*e^(12/19)*ei^(30/19)*(4*\[Delta]*(-5883581 + 3034822*\[Nu])*
           \[Chi]A + (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*
           \[Chi]S) - 5415*ei^(42/19)*(\[Delta]*(4888427 + 410172*\[Nu])*
           \[Chi]A + (4888427 - 7631920*\[Nu] + 74256*\[Nu]^2)*\[Chi]S))/
       (10270374912*Sqrt[6]*ei^(42/19)) + 
      (e^(12/19)*(5760*e^(12/19)*ei^(56/19)*(-1243916559 + 18409137292*\[Nu])*
          (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) + 
         420*e^(56/19)*ei^(12/19)*(-724653277 + 1072804096*\[Nu])*
          (-99391*\[Delta]*\[Chi]A + (-99391 + 22844*\[Nu])*\[Chi]S) + 
         728*e^(50/19)*ei^(18/19)*(-145417 + 239316*\[Nu])*
          (-113175949*\[Delta]*\[Chi]A + (-113175949 + 41301776*\[Nu])*
            \[Chi]S) + 1365*e^(18/19)*ei^(50/19)*(-2833 + 5516*\[Nu])*
          (1171697071*\[Delta]*\[Chi]A + (1171697071 + 6870205636*\[Nu])*
            \[Chi]S) + 168*e^2*ei^(30/19)*(\[Delta]*(-118659872446897 + 
             113510738272204*\[Nu])*\[Chi]A + (-118659872446897 + 
             286545005193108*\[Nu] - 61303880674400*\[Nu]^2)*\[Chi]S) + 
         1209572*e^(68/19)*(\[Delta]*(7875278603 + 32255910204*\[Nu])*
            \[Chi]A + (7875278603 + 99877088*\[Nu] - 7901609520*\[Nu]^2)*
            \[Chi]S) - 230297088*ei^(68/19)*
          (4*\[Delta]*(-5883581 + 3034822*\[Nu])*\[Chi]A + 
           (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*\[Chi]S) - 
         1954815*e^(26/19)*ei^(42/19)*(\[Delta]*(3593809211 + 
             51319330860*\[Nu])*\[Chi]A + (3593809211 - 4349305600*\[Nu] + 
             75051987024*\[Nu]^2)*\[Chi]S) + 504*e^(30/19)*ei^2*
          (7*\[Delta]*(-6961707660527 + 23826828158692*\[Nu])*\[Chi]A + 
           (-48731953623689 + 56400792295184*\[Nu] + 212302080351088*\[Nu]^2)*
            \[Chi]S)))/(1542363822784512*Sqrt[6]*ei^(42/19)) + 
      (e^(12/19)*(28302144*e^(50/19)*ei^(56/19)*(-1243916559 + 
           18409137292*\[Nu])*(113175949*\[Delta]*\[Chi]A + 
           (113175949 - 41301776*\[Nu])*\[Chi]S) + 812779520*e^(12/19)*
          ei^(94/19)*(-23500911051568 + 15843764110629*\[Nu])*
          (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) + 
         5159245*e^(94/19)*ei^(12/19)*(-96101931245055 + 132465768779708*
            \[Nu])*(-99391*\[Delta]*\[Chi]A + (-99391 + 22844*\[Nu])*
            \[Chi]S) + 20636980*e^(56/19)*ei^(50/19)*(-724653277 + 
           1072804096*\[Nu])*(1171697071*\[Delta]*\[Chi]A + 
           (1171697071 + 6870205636*\[Nu])*\[Chi]S) + 688278864*e^(88/19)*
          ei^(18/19)*(-145417 + 239316*\[Nu])*(-189611322769*\[Delta]*
            \[Chi]A + (-189611322769 + 49774532976*\[Nu])*\[Chi]S) + 
         1206660*e^(18/19)*ei^(88/19)*(-2833 + 5516*\[Nu])*
          (399181778038271*\[Delta]*\[Chi]A + (399181778038271 + 
             196461331296836*\[Nu])*\[Chi]S) + 91*e^(30/19)*ei^4*
          (\[Delta]*(-67851876903305430734671 + 8143200327587771321148*\[Nu])*
            \[Chi]A + (-67851876903305430734671 + 280784879015624996422712*
              \[Nu] - 169173767949134972320752*\[Nu]^2)*\[Chi]S) + 
         24752*e^4*ei^(30/19)*(\[Delta]*(-881417005050492936489 + 
             1269739664742855690748*\[Nu])*\[Chi]A + 
           (-881417005050492936489 + 2871898607753603335636*\[Nu] - 
             610400986757547833120*\[Nu]^2)*\[Chi]S) - 108805778368*e^2*
          ei^(68/19)*(\[Delta]*(-118659872446897 + 113510738272204*\[Nu])*
            \[Chi]A + (-118659872446897 + 286545005193108*\[Nu] - 
             61303880674400*\[Nu]^2)*\[Chi]S) + 1662935054833416*e^(106/19)*
          (\[Delta]*(7875278603 + 32255910204*\[Nu])*\[Chi]A + 
           (7875278603 + 99877088*\[Nu] - 7901609520*\[Nu]^2)*\[Chi]S) + 
         21856388688017520*ei^(106/19)*(4*\[Delta]*(-5883581 + 3034822*\[Nu])*
            \[Chi]A + (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*
            \[Chi]S) + 1152088300272*e^(68/19)*ei^2*
          (7*\[Delta]*(-6961707660527 + 23826828158692*\[Nu])*\[Chi]A + 
           (-48731953623689 + 56400792295184*\[Nu] + 212302080351088*\[Nu]^2)*
            \[Chi]S) - 163719665880*e^(64/19)*ei^(42/19)*
          (\[Delta]*(-99287898431985 + 1565517614766014*\[Nu])*\[Chi]A + 
           5*(-19857579686397 + 6993230143212*\[Nu] + 441867124045984*
              \[Nu]^2)*\[Chi]S)))/(875469182646356888518656*Sqrt[6]*
        ei^(42/19))) + \[Epsilon]^3*
     (((1885*(-e^(30/19) + e^(12/19)*ei^(18/19))*Pi)/(3192*Sqrt[6]*
          ei^(30/19)) - (5*e^(12/19)*(25055420*e^(56/19) - 
           27470894*e^2*ei^(18/19) - 925965*e^(26/19)*ei^(30/19) - 
           6596977*e^(18/19)*ei^2 + 9938416*ei^(56/19))*Pi)/
         (73747968*Sqrt[6]*ei^(30/19)) + 
        (5*e^(12/19)*(-9269831689358240*e^(94/19) + 9995224240679288*e^4*
            ei^(18/19) + 530147328125031*e^(64/19)*ei^(30/19) + 
           3801670677702820*e^(56/19)*ei^2 - 6279394966669792*e^2*
            ei^(56/19) + 889287733075248*e^(18/19)*ei^4 + 
           332896676445645*ei^(94/19))*Pi)/(14774283239718912*Sqrt[6]*
          ei^(30/19)))*x^(5/2) + SO*x^(5/2)*
       ((e^(30/19)*(99391*\[Delta]*\[Chi]A + (99391 - 22844*\[Nu])*\[Chi]S) - 
          160*e^(12/19)*ei^(18/19)*(157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*
             \[Chi]S) - 57*ei^(30/19)*(1303*\[Delta]*\[Chi]A + 
            (1303 - 92*\[Nu])*\[Chi]S))/(76608*Sqrt[6]*ei^(30/19)) + 
        (e^(12/19)*(-208*e^2*ei^(18/19)*(92307509*\[Delta]*\[Chi]A + 
             (92307509 - 26680576*\[Nu])*\[Chi]S) - 308655*e^(26/19)*
            ei^(30/19)*(41497*\[Delta]*\[Chi]A + (41497 - 313124*\[Nu])*
              \[Chi]S) + 431990*e^(56/19)*(99391*\[Delta]*\[Chi]A + 
             (99391 - 22844*\[Nu])*\[Chi]S) + 27416320*ei^(56/19)*
            (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) - 
           13*e^(18/19)*ei^2*(1171697071*\[Delta]*\[Chi]A + 
             (1171697071 + 6870205636*\[Nu])*\[Chi]S)))/(11504683008*Sqrt[6]*
          ei^(30/19)) + (e^(12/19)*(-23920*e^4*ei^(18/19)*
            (6334364140217*\[Delta]*\[Chi]A + (6334364140217 - 1267078023544*
                \[Nu])*\[Chi]S) - 29081782755*e^(64/19)*ei^(30/19)*
            (1225909*\[Delta]*\[Chi]A + (1225909 - 36705500*\[Nu])*\[Chi]S) + 
           914334272*e^2*ei^(56/19)*(92307509*\[Delta]*\[Chi]A + 
             (92307509 - 26680576*\[Nu])*\[Chi]S) + 3073551621140*e^(94/19)*
            (99391*\[Delta]*\[Chi]A + (99391 - 22844*\[Nu])*\[Chi]S) - 
           17660301137700*ei^(94/19)*(157*\[Delta]*\[Chi]A + 
             (157 - 110*\[Nu])*\[Chi]S) - 144068665*e^(56/19)*ei^2*
            (1171697071*\[Delta]*\[Chi]A + (1171697071 + 6870205636*\[Nu])*
              \[Chi]S) - 78*e^(18/19)*ei^4*(399181778038271*\[Delta]*
              \[Chi]A + (399181778038271 + 196461331296836*\[Nu])*\[Chi]S)))/
         (44322849719156736*Sqrt[6]*ei^(30/19)))) + 
    \[Epsilon]^4*(x^3*((5*(e^(12/19) - ei^(12/19))*
          (361*ei^(24/19)*(-151877213 - 187208280*\[Nu] + 39054960*\[Nu]^2) - 
           (e*ei)^(12/19)*(31176362099 + 92144559528*\[Nu] + 
             14676980976*\[Nu]^2) + e^(24/19)*(50392977379 - 
             385204834728*\[Nu] + 246696296112*\[Nu]^2)))/
         (338922372096*Sqrt[6]*ei^(36/19)) + 
        (5*e^(12/19)*(-6597052*ei^(62/19)*(-358353209 + 372157128*\[Nu] + 
             435997296*\[Nu]^2) + 302393*e^(62/19)*(50392977379 - 
             385204834728*\[Nu] + 246696296112*\[Nu]^2) + 260642*e^(26/19)*
            ei^(36/19)*(1023823302491 - 4774321617240*\[Nu] + 
             2523333592080*\[Nu]^2) - 1617*e^(50/19)*ei^(12/19)*
            (27565325190597 - 84987315736504*\[Nu] + 65207541903504*
              \[Nu]^2) + 1188*e^(12/19)*ei^(50/19)*(3524015611647 - 
             59014529687680*\[Nu] + 101544801302672*\[Nu]^2) + 
           132*e^2*ei^(24/19)*(-77776941858793 - 57290524699128*\[Nu] + 
             226944607584240*\[Nu]^2) - 32*e^(24/19)*ei^2*(7306273592457739 - 
             40736695045546269*\[Nu] + 24205030226070876*\[Nu]^2)))/
         (29690503588601856*Sqrt[6]*ei^(36/19)) + 
        (e^(12/19)*(12*e^(24/19)*ei^4*(15782029726452124264043 + 
             43183078773123411060120*\[Nu] - 21685168802452643835984*
              \[Nu]^2) + 6625922578275*ei^(100/19)*(-358353209 + 
             372157128*\[Nu] + 435997296*\[Nu]^2) + 3877164058040*e^(100/19)*
            (50392977379 - 385204834728*\[Nu] + 246696296112*\[Nu]^2) - 
           663536640*e^(62/19)*ei^2*(7306273592457739 - 40736695045546269*
              \[Nu] + 24205030226070876*\[Nu]^2) + 329321167*e^(64/19)*
            ei^(36/19)*(15334587311436637 - 76415710261605480*\[Nu] + 
             39628755350153424*\[Nu]^2) + 1774080*e^(12/19)*ei^(88/19)*
            (66578081009092144 - 174516409085861045*\[Nu] + 87394202834229564*
              \[Nu]^2) - 900900*e^(88/19)*ei^(12/19)*(741588616321469453 - 
             2221703929192049176*\[Nu] + 1647787690719105936*\[Nu]^2) + 
           31680*e^4*ei^(24/19)*(-8259865818436428113 - 2803675687047024972*
              \[Nu] + 16063790724563400000*\[Nu]^2) + 
           10560*e^2*(-85676*ei^(62/19)*(-77776941858793 - 57290524699128*
                \[Nu] + 226944607584240*\[Nu]^2) + 63*ei^2*(e*ei)^(12/19)*
              (235797495882898419 - 3828581745994024968*\[Nu] + 
               5016023129984824048*\[Nu]^2))))/(35670408535374978613248*
          Sqrt[6]*ei^(36/19))) + SO^2*x^3*
       ((5*(2*e^(12/19)*ei^(24/19)*(178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 
             356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 
             382*\[Delta]*\[Chi]A*\[Chi]S + 191*\[Chi]S^2 - 
             52*\[Nu]*\[Chi]S^2) + 19*ei^(36/19)*(32*\[Delta]*\[Kappa]A + 
             32*\[Kappa]S - 64*\[Kappa]S*\[Nu] + 33*\[Chi]A^2 - 
             128*\[Nu]*\[Chi]A^2 + 66*\[Delta]*\[Chi]A*\[Chi]S + 
             33*\[Chi]S^2 - 4*\[Nu]*\[Chi]S^2) + e^(36/19)*
            (964*\[Kappa]S*(-1 + 2*\[Nu]) - 1009*\[Chi]A^2 + 
             3856*\[Nu]*\[Chi]A^2 - 1009*\[Chi]S^2 + 180*\[Nu]*\[Chi]S^2 - 
             2*\[Delta]*(482*\[Kappa]A + 1009*\[Chi]A*\[Chi]S))))/
         (8512*Sqrt[6]*ei^(36/19)) - 
        (5*e^(12/19)*(599732*ei^(62/19)*(178*\[Delta]*\[Kappa]A + 
             178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 
             712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*\[Chi]S + 
             191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) - 13*e^(24/19)*ei^2*
            (10216736*\[Kappa]S*(-1 + 2*\[Nu]) + 50378283*\[Chi]A^2 + 
             40866944*\[Nu]*\[Chi]A^2 + 50378283*\[Chi]S^2 - 
             242380076*\[Nu]*\[Chi]S^2 - 2*\[Delta]*(5108368*\[Kappa]A - 
               50378283*\[Chi]A*\[Chi]S)) + 192052*e^(26/19)*ei^(36/19)*
            (2986*\[Kappa]S*(-1 + 2*\[Nu]) + 953*\[Chi]A^2 + 
             11944*\[Nu]*\[Chi]A^2 + 953*\[Chi]S^2 - 15756*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(-2986*\[Kappa]A + 1906*\[Chi]A*\[Chi]S)) + 
           907179*e^(62/19)*(-964*\[Kappa]S*(-1 + 2*\[Nu]) + 1009*\[Chi]A^2 - 
             3856*\[Nu]*\[Chi]A^2 + 1009*\[Chi]S^2 - 180*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(964*\[Kappa]A + 2018*\[Chi]A*\[Chi]S)) - 
           156*e^2*ei^(24/19)*(-3465533*\[Kappa]S*(-1 + 2*\[Nu]) + 
             3576925*\[Chi]A^2 - 13862132*\[Nu]*\[Chi]A^2 + 
             3576925*\[Chi]S^2 - 445568*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(3465533*\[Kappa]A + 7153850*\[Chi]A*\[Chi]S))))/
         (2237021696*Sqrt[6]*ei^(36/19)) + (Sqrt[3/2]*e^(12/19)*
          (15445040975*ei^(100/19)*(178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 
             356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 
             382*\[Delta]*\[Chi]A*\[Chi]S + 191*\[Chi]S^2 - 
             52*\[Nu]*\[Chi]S^2) + 6911840*e^(62/19)*ei^2*
            (10216736*\[Kappa]S*(-1 + 2*\[Nu]) + 50378283*\[Chi]A^2 + 
             40866944*\[Nu]*\[Chi]A^2 + 50378283*\[Chi]S^2 - 
             242380076*\[Nu]*\[Chi]S^2 - 2*\[Delta]*(5108368*\[Kappa]A - 
               50378283*\[Chi]A*\[Chi]S)) + 51998079*e^(64/19)*ei^(36/19)*
            (-4021854*\[Kappa]S*(-1 + 2*\[Nu]) - 3794671*\[Chi]A^2 - 
             16087416*\[Nu]*\[Chi]A^2 - 3794671*\[Chi]S^2 + 
             31266100*\[Nu]*\[Chi]S^2 + \[Delta]*(4021854*\[Kappa]A - 7589342*
                \[Chi]A*\[Chi]S)) - 298243389080*e^(100/19)*
            (-964*\[Kappa]S*(-1 + 2*\[Nu]) + 1009*\[Chi]A^2 - 
             3856*\[Nu]*\[Chi]A^2 + 1009*\[Chi]S^2 - 180*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(964*\[Kappa]A + 2018*\[Chi]A*\[Chi]S)) - 
           27416320*e^2*ei^(62/19)*(-3465533*\[Kappa]S*(-1 + 2*\[Nu]) + 
             3576925*\[Chi]A^2 - 13862132*\[Nu]*\[Chi]A^2 + 
             3576925*\[Chi]S^2 - 445568*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(3465533*\[Kappa]A + 7153850*\[Chi]A*\[Chi]S)) + 
           80*e^4*ei^(24/19)*(-2530412457698*\[Kappa]S*(-1 + 2*\[Nu]) + 
             2583790126097*\[Chi]A^2 - 10121649830792*\[Nu]*\[Chi]A^2 + 
             2583790126097*\[Chi]S^2 - 213510673596*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(2530412457698*\[Kappa]A + 5167580252194*\[Chi]A*
                \[Chi]S)) + 8*e^(24/19)*ei^4*(-4853042174758*\[Kappa]S*
              (-1 + 2*\[Nu]) + 4806182243753*\[Chi]A^2 - 19412168699032*\[Nu]*
              \[Chi]A^2 + 4806182243753*\[Chi]S^2 + 187439724020*\[Nu]*
              \[Chi]S^2 + \[Delta]*(4853042174758*\[Kappa]A + 9612364487506*
                \[Chi]A*\[Chi]S))))/(206736597057536*ei^(36/19)))) + 
    \[Epsilon]^6*
     (SO*x^4*((5*Pi*(-7*e^(48/19)*(6762971*\[Delta]*\[Chi]A + 
             (6762971 - 1994668*\[Nu])*\[Chi]S) + 377*e^(30/19)*ei^(18/19)*
            (124511*\[Delta]*\[Chi]A + (124511 - 40444*\[Nu])*\[Chi]S) + 
           149454*ei^(48/19)*(\[Delta]*\[Chi]A + \[Chi]S - 4*\[Nu]*\[Chi]S) + 
           8*e^(12/19)*ei^(36/19)*(31337*\[Delta]*\[Chi]A + 
             (31337 + 235316*\[Nu])*\[Chi]S)))/(34933248*Sqrt[6]*
          ei^(48/19)) - (e^(12/19)*Pi*(-13*e^(36/19)*ei^2*
            (5132086523719*\[Delta]*\[Chi]A + (5132086523719 - 7005190158236*
                \[Nu])*\[Chi]S) - 455*e^(56/19)*ei^(18/19)*
            (408268760993*\[Delta]*\[Chi]A + (408268760993 - 111244601092*
                \[Nu])*\[Chi]S) + 24191440*e^(74/19)*
            (6762971*\[Delta]*\[Chi]A + (6762971 - 1994668*\[Nu])*\[Chi]S) + 
           6854080*ei^(74/19)*(31337*\[Delta]*\[Chi]A + 
             (31337 + 235316*\[Nu])*\[Chi]S) - 1172889*e^(26/19)*ei^(48/19)*
            (11122477*\[Delta]*\[Chi]A + (11122477 + 246889772*\[Nu])*
              \[Chi]S) + 1040*e^2*ei^(36/19)*(65005615159*\[Delta]*\[Chi]A + 
             (65005615159 + 28798683196*\[Nu])*\[Chi]S) + 
           65*e^(18/19)*ei^(56/19)*(524587826887*\[Delta]*\[Chi]A + 
             (524587826887 + 2532014127172*\[Nu])*\[Chi]S)))/
         (5246135451648*Sqrt[6]*ei^(48/19)) - 
        (e^(12/19)*Pi*(-41860*e^(94/19)*ei^(18/19)*(680665042775389337*
              \[Delta]*\[Chi]A + (680665042775389337 - 167827572913243108*
                \[Nu])*\[Chi]S) - 3227138096*e^(74/19)*ei^2*
            (5132086523719*\[Delta]*\[Chi]A + (5132086523719 - 7005190158236*
                \[Nu])*\[Chi]S) + 3535661359739600*e^(112/19)*
            (6762971*\[Delta]*\[Chi]A + (6762971 - 1994668*\[Nu])*\[Chi]S) - 
           61811053981950*ei^(112/19)*(31337*\[Delta]*\[Chi]A + 
             (31337 + 235316*\[Nu])*\[Chi]S) - 12278974941*e^(64/19)*
            ei^(48/19)*(270023886961*\[Delta]*\[Chi]A + 
             (270023886961 + 6264760409036*\[Nu])*\[Chi]S) + 
           5460*e^(18/19)*ei^(94/19)*(186291062139098407*\[Delta]*\[Chi]A + 
             (186291062139098407 + 48983447376271972*\[Nu])*\[Chi]S) + 
           43056*e^(36/19)*ei^4*(56117733091629833*\[Delta]*\[Chi]A + 
             (56117733091629833 + 92077887474718076*\[Nu])*\[Chi]S) + 
           33488*e^4*ei^(36/19)*(467358287836566959*\[Delta]*\[Chi]A + 
             (467358287836566959 + 252002971654484168*\[Nu])*\[Chi]S) + 
           23345*e^2*(-2741632*ei^(74/19)*(65005615159*\[Delta]*\[Chi]A + 
               65005615159*\[Chi]S + 28798683196*\[Nu]*\[Chi]S) + 
             13*ei^2*(e*ei)^(18/19)*(31424619221651087*\[Delta]*\[Chi]A + 
               31424619221651087*\[Chi]S + 143440553890184612*\[Nu]*
                \[Chi]S))))/(282957072607096602624*Sqrt[6]*ei^(48/19))) + 
      SO^2*x^4*((-5*(-18*(e*ei)^(24/19)*(-145417 + 239316*\[Nu])*
            (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
             191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
              \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) - 
           e^(48/19)*(-252*\[Delta]*\[Kappa]A*(3243271 + 3227408*\[Nu]) + 
             252*\[Kappa]S*(-3243271 + 3259134*\[Nu] + 17101656*\[Nu]^2) - 
             4816146012*\[Chi]A^2 + 5914074257*\[Delta]^2*\[Chi]A^2 + 
             17498024628*\[Nu]*\[Chi]A^2 + 8619234624*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(1097928245 + 67532948*\[Nu])*\[Chi]A*\[Chi]S + 
             1097928245*\[Chi]S^2 + 1901625316*\[Nu]*\[Chi]S^2 - 
             958689760*\[Nu]^2*\[Chi]S^2) + 3249*ei^(48/19)*
            (8*\[Delta]*\[Kappa]A*(-44803 + 25242*\[Nu]) + 
             8*\[Kappa]S*(-44803 + 114848*\[Nu] + 10416*\[Nu]^2) - 
             1234422*\[Chi]A^2 + 1030673*\[Delta]^2*\[Chi]A^2 + 
             4987820*\[Nu]*\[Chi]A^2 + 166656*\[Nu]^2*\[Chi]A^2 + 
             14*\[Delta]*(-29107 + 79316*\[Nu])*\[Chi]A*\[Chi]S - 
             203749*\[Chi]S^2 + 1060292*\[Nu]*\[Chi]S^2 - 382144*\[Nu]^2*
              \[Chi]S^2) + 112*e^(30/19)*ei^(18/19)*
            (15604387*\[Delta]^2*\[Chi]A^2 + 2*\[Delta]*(15604387 - 7259759*
                \[Nu])*\[Chi]A*\[Chi]S + (15604387 - 14519518*\[Nu] + 2512840*
                \[Nu]^2)*\[Chi]S^2) - 324*e^(36/19)*ei^(12/19)*
            (-2833 + 5516*\[Nu])*(-964*\[Kappa]S*(-1 + 2*\[Nu]) + 
             1009*\[Chi]A^2 - 3856*\[Nu]*\[Chi]A^2 + 1009*\[Chi]S^2 - 
             180*\[Nu]*\[Chi]S^2 + \[Delta]*(964*\[Kappa]A + 2018*\[Chi]A*
                \[Chi]S)) + 8*e^(12/19)*ei^(36/19)*
            (-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 116610*\[Nu]^2) - 
             278950701*\[Chi]A^2 + 102215792*\[Delta]^2*\[Chi]A^2 + 
             1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*\[Chi]A^2 - 
             176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
             68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*
                (-124448 + 126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*
                \[Chi]A*\[Chi]S))))/(2934392832*Sqrt[6]*ei^(48/19)) + 
        (e^(12/19)*(-810*e^(12/19)*ei^(62/19)*(-1243916559 + 
             18409137292*\[Nu])*(178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 
             356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 
             382*\[Delta]*\[Chi]A*\[Chi]S + 191*\[Chi]S^2 - 
             52*\[Nu]*\[Chi]S^2) + 3023930*e^(74/19)*
            (-252*\[Delta]*\[Kappa]A*(3243271 + 3227408*\[Nu]) + 
             252*\[Kappa]S*(-3243271 + 3259134*\[Nu] + 17101656*\[Nu]^2) - 
             4816146012*\[Chi]A^2 + 5914074257*\[Delta]^2*\[Chi]A^2 + 
             17498024628*\[Nu]*\[Chi]A^2 + 8619234624*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(1097928245 + 67532948*\[Nu])*\[Chi]A*\[Chi]S + 
             1097928245*\[Chi]S^2 + 1901625316*\[Nu]*\[Chi]S^2 - 
             958689760*\[Nu]^2*\[Chi]S^2) + e^(36/19)*ei^2*
            (36*\[Delta]*\[Kappa]A*(2199940776583 + 45382625640716*\[Nu]) + 
             36*\[Kappa]S*(2199940776583 + 40982744087550*\[Nu] + 
               124668786014136*\[Nu]^2) - 40704717297445158*\[Chi]A^2 + 
             42217477632040165*\[Delta]^2*\[Chi]A^2 + 160552218323235072*
              \[Nu]*\[Chi]A^2 + 8976152593017792*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(1512760334595007 + 436752419096266*\[Nu])*\[Chi]A*
              \[Chi]S + 1512760334595007*\[Chi]S^2 + 3140155704738092*\[Nu]*
              \[Chi]S^2 + 2296792545445216*\[Nu]^2*\[Chi]S^2) + 
           6370*e^(18/19)*ei^(56/19)*(183956440147*\[Delta]^2*\[Chi]A^2 + 
             2*\[Delta]*(183956440147 + 474867803521*\[Nu])*\[Chi]A*\[Chi]S + 
             (183956440147 + 949735607042*\[Nu] - 755722619960*\[Nu]^2)*
              \[Chi]S^2) - 637*e^(56/19)*ei^(18/19)*
            (12285738307079*\[Delta]^2*\[Chi]A^2 + 2*\[Delta]*
              (12285738307079 - 3827691681826*\[Nu])*\[Chi]A*\[Chi]S + 
             (12285738307079 - 7655383363652*\[Nu] + 1110501117344*\[Nu]^2)*
              \[Chi]S^2) - 10530*e^(24/19)*ei^(50/19)*(-2833 + 5516*\[Nu])*
            (10216736*\[Kappa]S*(-1 + 2*\[Nu]) + 50378283*\[Chi]A^2 + 
             40866944*\[Nu]*\[Chi]A^2 + 50378283*\[Chi]S^2 - 
             242380076*\[Nu]*\[Chi]S^2 - 2*\[Delta]*(5108368*\[Kappa]A - 
               50378283*\[Chi]A*\[Chi]S)) + 2835*e^(62/19)*ei^(12/19)*
            (-1571689321 + 2383893876*\[Nu])*(-964*\[Kappa]S*(-1 + 2*\[Nu]) + 
             1009*\[Chi]A^2 - 3856*\[Nu]*\[Chi]A^2 + 1009*\[Chi]S^2 - 
             180*\[Nu]*\[Chi]S^2 + \[Delta]*(964*\[Kappa]A + 2018*\[Chi]A*
                \[Chi]S)) + 585*e^(50/19)*ei^(24/19)*(-145417 + 239316*\[Nu])*
            (-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 25904401*\[Chi]A^2 - 
             99734624*\[Nu]*\[Chi]A^2 + 25904401*\[Chi]S^2 - 
             3882980*\[Nu]*\[Chi]S^2 + \[Delta]*(24933656*\[Kappa]A + 
               51808802*\[Chi]A*\[Chi]S)) + 5997320*ei^(74/19)*
            (-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 116610*\[Nu]^2) - 
             278950701*\[Chi]A^2 + 102215792*\[Delta]^2*\[Chi]A^2 + 
             1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*\[Chi]A^2 - 
             176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
             68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*
                (-124448 + 126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*
                \[Chi]A*\[Chi]S)) - 2345778*e^(26/19)*ei^(48/19)*
            (4*\[Kappa]S*(-438867911 + 1043712010*\[Nu] + 802045608*
                \[Nu]^2) - 14355607086*\[Chi]A^2 + 13638069109*\[Delta]^2*
              \[Chi]A^2 + 57157773104*\[Nu]*\[Chi]A^2 + 6416364864*\[Nu]^2*
              \[Chi]A^2 - 717537977*\[Chi]S^2 + 2293324302*\[Nu]*\[Chi]S^2 + 
             3402549976*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(\[Kappa]A*
                (-877735822 + 331952376*\[Nu]) + (-717537977 + 1014334531*
                  \[Nu])*\[Chi]A*\[Chi]S)) - 10*e^2*ei^(36/19)*
            (-18*\[Kappa]S*(30936092472295 - 97020967578282*\[Nu] + 
               36669175275096*\[Nu]^2) - 2827049041226343*\[Chi]A^2 + 
             2206814975744288*\[Delta]^2*\[Chi]A^2 + 11721265129529220*\[Nu]*
              \[Chi]A^2 - 1320090309903456*\[Nu]^2*\[Chi]A^2 - 
             620234065482055*\[Chi]S^2 + 1521486691850320*\[Nu]*\[Chi]S^2 - 
             141962774214352*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(9*\[Kappa]A*
                (-30936092472295 + 35148782633692*\[Nu]) + 
               (-620234065482055 + 967277828237084*\[Nu])*\[Chi]A*\[Chi]S))))/
         (385590955696128*Sqrt[6]*ei^(48/19)) + 
        (e^(12/19)*(-564762240*e^(12/19)*ei^(100/19)*(-23500911051568 + 
             15843764110629*\[Nu])*(178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 
             356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 
             382*\[Delta]*\[Chi]A*\[Chi]S + 191*\[Chi]S^2 - 
             52*\[Nu]*\[Chi]S^2) - 16744*e^4*ei^(36/19)*
            (522*\[Delta]*\[Kappa]A*(-5919500807162614699 + 
               7069747731165174748*\[Nu]) - 522*\[Kappa]S*
              (5919500807162614699 - 18908749345490404146*\[Nu] + 
               7669345784813238456*\[Nu]^2) - 18928431938173283618721*
              \[Chi]A^2 + 16007908125685966389172*\[Delta]^2*\[Chi]A^2 + 
             78468871437192708323052*\[Nu]*\[Chi]A^2 - 8006796999345020948064*
              \[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(-2920523812487317229549 + 
               5349738482100710957140*\[Nu])*\[Chi]A*\[Chi]S - 
             2920523812487317229549*\[Chi]S^2 + 7944333279701848066112*\[Nu]*
              \[Chi]S^2 - 384534241038124715888*\[Nu]^2*\[Chi]S^2) + 
           22981798838307400*e^(112/19)*(-252*\[Delta]*\[Kappa]A*
              (3243271 + 3227408*\[Nu]) + 252*\[Kappa]S*(-3243271 + 3259134*
                \[Nu] + 17101656*\[Nu]^2) - 4816146012*\[Chi]A^2 + 
             5914074257*\[Delta]^2*\[Chi]A^2 + 17498024628*\[Nu]*\[Chi]A^2 + 
             8619234624*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(1097928245 + 67532948*
                \[Nu])*\[Chi]A*\[Chi]S + 1097928245*\[Chi]S^2 + 
             1901625316*\[Nu]*\[Chi]S^2 - 958689760*\[Nu]^2*\[Chi]S^2) + 
           12908552384*e^(74/19)*ei^2*(36*\[Delta]*\[Kappa]A*
              (2199940776583 + 45382625640716*\[Nu]) + 36*\[Kappa]S*
              (2199940776583 + 40982744087550*\[Nu] + 124668786014136*
                \[Nu]^2) - 40704717297445158*\[Chi]A^2 + 42217477632040165*
              \[Delta]^2*\[Chi]A^2 + 160552218323235072*\[Nu]*\[Chi]A^2 + 
             8976152593017792*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*
              (1512760334595007 + 436752419096266*\[Nu])*\[Chi]A*\[Chi]S + 
             1512760334595007*\[Chi]S^2 + 3140155704738092*\[Nu]*\[Chi]S^2 + 
             2296792545445216*\[Nu]^2*\[Chi]S^2) - 1794*e^(36/19)*ei^4*
            (1044*\[Delta]*\[Kappa]A*(-7365045929568369393 + 
               5742829720995328276*\[Nu]) - 1044*\[Kappa]S*
              (7365045929568369393 - 20472921580132067062*\[Nu] + 
               1517701378547771592*\[Nu]^2) - 36820157636844044300469*
              \[Chi]A^2 + 30945138675512339580484*\[Delta]^2*\[Chi]A^2 + 
             153005155284657572553108*\[Nu]*\[Chi]A^2 - 
             3168960478407747084096*\[Nu]^2*\[Chi]A^2 + 86*\[Delta]*
              (-136628347937946621395 + 186569357031292881772*\[Nu])*\[Chi]A*
              \[Chi]S - 5875018961331704719985*\[Chi]S^2 + 
             10320439967409792481160*\[Nu]*\[Chi]S^2 + 5836439695605526530544*
              \[Nu]^2*\[Chi]S^2) + 27824160*e^(18/19)*ei^(94/19)*
            (62671539152008547*\[Delta]^2*\[Chi]A^2 - 2*\[Delta]*
              (-62671539152008547 + 6532783285303279*\[Nu])*\[Chi]A*\[Chi]S + 
             (62671539152008547 - 13065566570606558*\[Nu] - 21610746442651960*
                \[Nu]^2)*\[Chi]S^2) - 21331856*e^(94/19)*ei^(18/19)*
            (3231434976186023027*\[Delta]^2*\[Chi]A^2 + 2*\[Delta]*
              (3231434976186023027 - 836498434314498718*\[Nu])*\[Chi]A*
              \[Chi]S + (3231434976186023027 - 1672996868628997436*\[Nu] + 
               213816389149906112*\[Nu]^2)*\[Chi]S^2) - 393316560*e^(62/19)*
            ei^(50/19)*(-1571689321 + 2383893876*\[Nu])*
            (10216736*\[Kappa]S*(-1 + 2*\[Nu]) + 50378283*\[Chi]A^2 + 
             40866944*\[Nu]*\[Chi]A^2 + 50378283*\[Chi]S^2 - 
             242380076*\[Nu]*\[Chi]S^2 - 2*\[Delta]*(5108368*\[Kappa]A - 
               50378283*\[Chi]A*\[Chi]S)) + 344151990*e^(100/19)*ei^(12/19)*
            (-117334502622439 + 164817749118204*\[Nu])*
            (-964*\[Kappa]S*(-1 + 2*\[Nu]) + 1009*\[Chi]A^2 - 
             3856*\[Nu]*\[Chi]A^2 + 1009*\[Chi]S^2 - 180*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(964*\[Kappa]A + 2018*\[Chi]A*\[Chi]S)) - 
           112376160*e^(50/19)*ei^(62/19)*(-1243916559 + 18409137292*\[Nu])*
            (-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 25904401*\[Chi]A^2 - 
             99734624*\[Nu]*\[Chi]A^2 + 25904401*\[Chi]S^2 - 
             3882980*\[Nu]*\[Chi]S^2 + \[Delta]*(24933656*\[Kappa]A + 
               51808802*\[Chi]A*\[Chi]S)) - 117994968*e^(24/19)*ei^(88/19)*
            (-2833 + 5516*\[Nu])*(-4853042174758*\[Kappa]S*(-1 + 2*\[Nu]) + 
             4806182243753*\[Chi]A^2 - 19412168699032*\[Nu]*\[Chi]A^2 + 
             4806182243753*\[Chi]S^2 + 187439724020*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(4853042174758*\[Kappa]A + 9612364487506*\[Chi]A*
                \[Chi]S)) - 2812402956178725*ei^(112/19)*
            (-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 116610*\[Nu]^2) - 
             278950701*\[Chi]A^2 + 102215792*\[Delta]^2*\[Chi]A^2 + 
             1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*\[Chi]A^2 - 
             176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
             68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*
                (-124448 + 126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*
                \[Chi]A*\[Chi]S)) - 12278974941*e^(64/19)*ei^(48/19)*
            (368*\[Kappa]S*(-5800115804177 + 9056885427979*\[Nu] + 
               28442618518314*\[Nu]^2) - 29636377823388555*\[Chi]A^2 + 
             31181399186580188*\[Delta]^2*\[Chi]A^2 + 119503932749870596*
              \[Nu]*\[Chi]A^2 + 20933767229479104*\[Nu]^2*\[Chi]A^2 + 
             1545021363191633*\[Chi]S^2 - 1650033761251200*\[Nu]*\[Chi]S^2 + 
             13396573285480112*\[Nu]^2*\[Chi]S^2 - 2*\[Delta]*
              (184*\[Kappa]A*(5800115804177 + 2543346180375*\[Nu]) + 
               (-1545021363191633 + 345806152467412*\[Nu])*\[Chi]A*
                \[Chi]S)) + 43701840*e^(88/19)*ei^(24/19)*
            (-145417 + 239316*\[Nu])*(-2948986828885*\[Kappa]S*
              (-1 + 2*\[Nu]) + \[Delta]*(2948986828885*\[Kappa]A + 
               6058122836776*\[Chi]A*\[Chi]S) - 4*((-757265354597 + 
                 2948986828885*\[Nu])*\[Chi]A^2 + (-757265354597 + 
                 80074589503*\[Nu])*\[Chi]S^2)) + 18676*e^2*
            (8281*ei^2*(e*ei)^(18/19)*(144833672963114999*\[Delta]^2*
                \[Chi]A^2 + 2*\[Delta]*(144833672963114999 + 
                 396134150578580894*\[Nu])*\[Chi]A*\[Chi]S + 
               (144833672963114999 + 792268301157161788*\[Nu] - 
                 333977019574551136*\[Nu]^2)*\[Chi]S^2) + 1713520*ei^(74/19)*
              (-18*\[Kappa]S*(30936092472295 - 97020967578282*\[Nu] + 
                 36669175275096*\[Nu]^2) - 2827049041226343*\[Chi]A^2 + 
               2206814975744288*\[Delta]^2*\[Chi]A^2 + 11721265129529220*
                \[Nu]*\[Chi]A^2 - 1320090309903456*\[Nu]^2*\[Chi]A^2 - 
               620234065482055*\[Chi]S^2 + 1521486691850320*\[Nu]*\[Chi]S^2 - 
               141962774214352*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*
                (9*\[Kappa]A*(-30936092472295 + 35148782633692*\[Nu]) + 
                 (-620234065482055 + 967277828237084*\[Nu])*\[Chi]A*
                  \[Chi]S)))))/(1081461931504323215228928*Sqrt[6]*
          ei^(48/19))))
 
HlmDCMem[4, 0] = ((e^(12/19) - ei^(12/19))/(504*Sqrt[2]*ei^(12/19)) + 
      (86398*e^2 + 1805*e^(26/19)*ei^(12/19) - 88203*ei^2)/
       (37844352*Sqrt[2]*(ei/e)^(12/19)) + 
      (15710523296*e^4 + 343395835*e^(64/19)*ei^(12/19) - 
        18758308416*e^2*ei^2 + 2704389285*ei^4)/(6994847268864*Sqrt[2]*
        (ei/e)^(12/19)))*x + x^2*\[Epsilon]^2*
     (((e^(12/19) - ei^(12/19))*(361*ei^(12/19)*(-9479 + 40124*\[Nu]) + 
         e^(12/19)*(-3920527 + 15455580*\[Nu])))/(283143168*Sqrt[2]*
        ei^(24/19)) + (e^(12/19)*(3*e^(12/19)*ei^2*(1859009125901 - 
           8124219720708*\[Nu]) + 27166524*ei^(50/19)*(-2833 + 5516*\[Nu]) + 
         302393*e^(50/19)*(-3920527 + 15455580*\[Nu]) + 
         54872*e^(26/19)*ei^(24/19)*(-84703117 + 364299110*\[Nu]) - 
         308*e^2*ei^(12/19)*(-1082158253 + 1430751140*\[Nu])))/
       (37206144847872*Sqrt[2]*ei^(24/19)) + 
      (e^(12/19)*(-2707093674285*ei^(88/19)*(-2833 + 5516*\[Nu]) + 
         283211804148*e^(88/19)*(-3920527 + 15455580*\[Nu]) + 
         217332192*e^2*ei^(50/19)*(-1082158253 + 1430751140*\[Nu]) - 
         4147104*e^(50/19)*ei^2*(-1859009125901 + 8124219720708*\[Nu]) + 
         2476099*e^(64/19)*ei^(24/19)*(-2906651517879 + 12249998813780*
            \[Nu]) - 8008*e^4*ei^(12/19)*(-42164487546687 + 
           54413587634588*\[Nu]) - 84*e^(12/19)*ei^4*(-5807829430373107 + 
           10448813703644316*\[Nu])))/(22349880034696101888*Sqrt[2]*
        ei^(24/19))) + SO*x^(7/2)*\[Epsilon]^5*
     ((70*e^(24/19)*ei^(18/19)*(-3920527 + 15455580*\[Nu])*
         (-157*\[Delta]*\[Chi]A + (-157 + 110*\[Nu])*\[Chi]S) - 
        1925*e^(30/19)*ei^(12/19)*(-2833 + 5516*\[Nu])*
         (5189*\[Delta]*\[Chi]A + (5189 + 10748*\[Nu])*\[Chi]S) + 
        1232*e^(12/19)*ei^(30/19)*(4*\[Delta]*(-5883581 + 3034822*\[Nu])*
           \[Chi]A + (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*
           \[Chi]S) - 21660*ei^(42/19)*(\[Delta]*(-1319709 + 12539422*\[Nu])*
           \[Chi]A + (-1319709 + 6548872*\[Nu] + 5249776*\[Nu]^2)*\[Chi]S) + 
        e^(42/19)*(\[Delta]*(-70975542727 + 481603460604*\[Nu])*\[Chi]A + 
          (-70975542727 + 289947696068*\[Nu] + 119614397280*\[Nu]^2)*
           \[Chi]S))/(1694611860480*Sqrt[2]*ei^(42/19)) + 
      (e^(12/19)*(240*e^(12/19)*ei^(56/19)*(-1859009125901 + 
           8124219720708*\[Nu])*(157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*
            \[Chi]S) - 7700*e^(56/19)*ei^(12/19)*(-724653277 + 
           1072804096*\[Nu])*(5189*\[Delta]*\[Chi]A + (5189 + 10748*\[Nu])*
            \[Chi]S) + 91*e^(50/19)*ei^(18/19)*(-3920527 + 15455580*\[Nu])*
          (-113175949*\[Delta]*\[Chi]A + (-113175949 + 41301776*\[Nu])*
            \[Chi]S) + 60060*e^(18/19)*ei^(50/19)*(-2833 + 5516*\[Nu])*
          (31277492*\[Delta]*\[Chi]A + (31277492 + 1656666347*\[Nu])*
            \[Chi]S) + 154*e^2*ei^(30/19)*(\[Delta]*(-118659872446897 + 
             113510738272204*\[Nu])*\[Chi]A + (-118659872446897 + 
             286545005193108*\[Nu] - 61303880674400*\[Nu]^2)*\[Chi]S) - 
         217332192*ei^(68/19)*(4*\[Delta]*(-5883581 + 3034822*\[Nu])*
            \[Chi]A + (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*
            \[Chi]S) + 604786*e^(68/19)*(\[Delta]*(-70975542727 + 
             481603460604*\[Nu])*\[Chi]A + (-70975542727 + 289947696068*
              \[Nu] + 119614397280*\[Nu]^2)*\[Chi]S) - 3909630*e^(26/19)*
          ei^(42/19)*(\[Delta]*(-193129614957 + 884416890364*\[Nu])*\[Chi]A + 
           (-193129614957 + 429487218226*\[Nu] + 833412262648*\[Nu]^2)*
            \[Chi]S) + 21*e^(30/19)*ei^2*(7*\[Delta]*(-4713900887677313 + 
             20661171666495828*\[Nu])*\[Chi]A + (-32997306213741191 + 
             73396240249492084*\[Nu] + 137612708364899040*\[Nu]^2)*\[Chi]S)))/
       (254490030759444480*Sqrt[2]*ei^(42/19)) + 
      (e^(12/19)*(3537768*e^(50/19)*ei^(56/19)*(-1859009125901 + 
           8124219720708*\[Nu])*(113175949*\[Delta]*\[Chi]A + 
           (113175949 - 41301776*\[Nu])*\[Chi]S) + 19049520*e^(12/19)*
          ei^(94/19)*(-5807829430373107 + 10448813703644316*\[Nu])*
          (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) - 
         283758475*e^(94/19)*ei^(12/19)*(-96101931245055 + 
           132465768779708*\[Nu])*(5189*\[Delta]*\[Chi]A + 
           (5189 + 10748*\[Nu])*\[Chi]S) + 2724081360*e^(56/19)*ei^(50/19)*
          (-724653277 + 1072804096*\[Nu])*(31277492*\[Delta]*\[Chi]A + 
           (31277492 + 1656666347*\[Nu])*\[Chi]S) + 258104574*e^(88/19)*
          ei^(18/19)*(-3920527 + 15455580*\[Nu])*(-189611322769*\[Delta]*
            \[Chi]A + (-189611322769 + 49774532976*\[Nu])*\[Chi]S) + 
         13273260*e^(18/19)*ei^(88/19)*(-2833 + 5516*\[Nu])*
          (70548776005597*\[Delta]*\[Chi]A + (70548776005597 + 
             205639257894652*\[Nu])*\[Chi]S) + 68068*e^4*ei^(30/19)*
          (\[Delta]*(-881417005050492936489 + 1269739664742855690748*\[Nu])*
            \[Chi]A + (-881417005050492936489 + 2871898607753603335636*
              \[Nu] - 610400986757547833120*\[Nu]^2)*\[Chi]S) - 
         308041215636*e^2*ei^(68/19)*(\[Delta]*(-118659872446897 + 
             113510738272204*\[Nu])*\[Chi]A + (-118659872446897 + 
             286545005193108*\[Nu] - 61303880674400*\[Nu]^2)*\[Chi]S) + 
         61391470345435230*ei^(106/19)*(4*\[Delta]*(-5883581 + 3034822*\[Nu])*
            \[Chi]A + (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*
            \[Chi]S) + 2494402582250124*e^(106/19)*
          (\[Delta]*(-70975542727 + 481603460604*\[Nu])*\[Chi]A + 
           (-70975542727 + 289947696068*\[Nu] + 119614397280*\[Nu]^2)*
            \[Chi]S) + 144011037534*e^(68/19)*ei^2*
          (7*\[Delta]*(-4713900887677313 + 20661171666495828*\[Nu])*\[Chi]A + 
           (-32997306213741191 + 73396240249492084*\[Nu] + 137612708364899040*
              \[Nu]^2)*\[Chi]S) - 40929916470*e^(64/19)*ei^(42/19)*
          (8*\[Delta]*(-16557657366135909 + 73963894327673623*\[Nu])*
            \[Chi]A + (-132461258929087272 + 283298459287770341*\[Nu] + 
             551036274021854348*\[Nu]^2)*\[Chi]S) - 273*e^(30/19)*ei^4*
          (\[Delta]*(-91579519882877550261559 + 1070089805632513857206652*
              \[Nu])*\[Chi]A + (-91579519882877550261559 - 
             1886006275288386965668792*\[Nu] + 3593907680873728598492112*
              \[Nu]^2)*\[Chi]S)))/(433357245409946659816734720*Sqrt[2]*
        ei^(42/19))) + \[Epsilon]^3*
     (((377*(-e^(30/19) + e^(12/19)*ei^(18/19))*Pi)/(114912*Sqrt[2]*
          ei^(30/19)) - (e^(12/19)*(50110840*e^(56/19) - 
           54941788*e^2*ei^(18/19) + 4629825*e^(26/19)*ei^(30/19) - 
           20261973*e^(18/19)*ei^2 + 20463096*ei^(56/19))*Pi)/
         (5309853696*Sqrt[2]*ei^(30/19)) + 
        (e^(12/19)*(-18539663378716480*e^(94/19) + 19990448481358576*e^4*
            ei^(18/19) - 2650736640625155*e^(64/19)*ei^(30/19) + 
           11676461601504180*e^(56/19)*ei^2 - 12929209445939952*e^2*
            ei^(56/19) + 1772656357202016*e^(18/19)*ei^4 + 
           680043025216815*ei^(94/19))*Pi)/(1063748393259761664*Sqrt[2]*
          ei^(30/19)))*x^(5/2) + SO*x^(5/2)*
       ((-8*e^(12/19)*ei^(18/19)*(157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*
             \[Chi]S) - 171*ei^(30/19)*(23*\[Delta]*\[Chi]A + 
            (23 + 68*\[Nu])*\[Chi]S) + e^(30/19)*(5189*\[Delta]*\[Chi]A + 
            (5189 + 10748*\[Nu])*\[Chi]S))/(689472*Sqrt[2]*ei^(30/19)) + 
        (e^(12/19)*(-26*e^2*ei^(18/19)*(92307509*\[Delta]*\[Chi]A + 
             (92307509 - 26680576*\[Nu])*\[Chi]S) - 308655*e^(26/19)*
            ei^(30/19)*(4271*\[Delta]*\[Chi]A + (4271 - 380058*\[Nu])*
              \[Chi]S) + 3528120*ei^(56/19)*(157*\[Delta]*\[Chi]A + 
             (157 - 110*\[Nu])*\[Chi]S) + 1079975*e^(56/19)*
            (5189*\[Delta]*\[Chi]A + (5189 + 10748*\[Nu])*\[Chi]S) - 
           78*e^(18/19)*ei^2*(31277492*\[Delta]*\[Chi]A + 
             (31277492 + 1656666347*\[Nu])*\[Chi]S)))/(258855367680*Sqrt[2]*
          ei^(30/19)) + (e^(12/19)*(-11960*e^4*ei^(18/19)*
            (6334364140217*\[Delta]*\[Chi]A + (6334364140217 - 1267078023544*
                \[Nu])*\[Chi]S) - 16156545975*e^(64/19)*ei^(30/19)*
            (408755*\[Delta]*\[Chi]A + (408755 - 335858562*\[Nu])*\[Chi]S) + 
           470651208*e^2*ei^(56/19)*(92307509*\[Delta]*\[Chi]A + 
             (92307509 - 26680576*\[Nu])*\[Chi]S) - 9019138265475*ei^(94/19)*
            (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) + 
           30735516211400*e^(94/19)*(5189*\[Delta]*\[Chi]A + 
             (5189 + 10748*\[Nu])*\[Chi]S) - 3457647960*e^(56/19)*ei^2*
            (31277492*\[Delta]*\[Chi]A + (31277492 + 1656666347*\[Nu])*
              \[Chi]S) - 156*e^(18/19)*ei^4*(70548776005597*\[Delta]*
              \[Chi]A + (70548776005597 + 205639257894652*\[Nu])*\[Chi]S)))/
         (3989056474724106240*Sqrt[2]*ei^(30/19)))) + 
    \[Epsilon]^4*(x^3*(((e^(12/19) - ei^(12/19))*
          (361*ei^(24/19)*(24215523937 - 140432459328*\[Nu] + 
             53657768304*\[Nu]^2) + (e*ei)^(12/19)*(9356738247901 - 
             51334739449056*\[Nu] + 18622282997808*\[Nu]^2) + 
           e^(24/19)*(17153749047583 - 97253461569600*\[Nu] + 
             78469874452368*\[Nu]^2)))/(317231340281856*Sqrt[2]*ei^(36/19)) + 
        (e^(12/19)*(-13583262*ei^(62/19)*(-358353209 + 372157128*\[Nu] + 
             435997296*\[Nu]^2) + 23261*e^(62/19)*(17153749047583 - 
             97253461569600*\[Nu] + 78469874452368*\[Nu]^2) + 
           260642*e^(26/19)*ei^(36/19)*(43165398517289 - 207702668397672*
              \[Nu] + 180915951308400*\[Nu]^2) + 264*e^2*ei^(24/19)*
            (-77776941858793 - 57290524699128*\[Nu] + 226944607584240*
              \[Nu]^2) - 441*e^(50/19)*ei^(12/19)*(743177219125107 - 
             3998012304886168*\[Nu] + 4211253658313520*\[Nu]^2) + 
           108*e^(12/19)*ei^(50/19)*(5266572853677533 - 33270208807235680*
              \[Nu] + 44813195979425328*\[Nu]^2) - 12*e^(24/19)*ei^2*
            (989592697031477957 - 4854040883091542142*\[Nu] + 
             4334686378660131480*\[Nu]^2)))/(2137716258379333632*Sqrt[2]*
          ei^(36/19)) + (e^(12/19)*(13535468371425*ei^(100/19)*
            (-358353209 + 372157128*\[Nu] + 435997296*\[Nu]^2) + 
           298243389080*e^(100/19)*(17153749047583 - 97253461569600*\[Nu] + 
             78469874452368*\[Nu]^2) + 329321167*e^(64/19)*ei^(36/19)*
            (700626752131304143 - 3371572907848553400*\[Nu] + 
             2927172440774116080*\[Nu]^2) - 248826240*e^(62/19)*ei^2*
            (989592697031477957 - 4854040883091542142*\[Nu] + 
             4334686378660131480*\[Nu]^2) + 63360*e^4*ei^(24/19)*
            (-8259865818436428113 - 2803675687047024972*\[Nu] + 
             16063790724563400000*\[Nu]^2) + 30240*e^(12/19)*ei^(88/19)*
            (16453580776247012131 - 61637476360362405440*\[Nu] + 
             57635656389302047056*\[Nu]^2) - 245700*e^(88/19)*ei^(12/19)*
            (19993660941849726443 - 105813828266912291512*\[Nu] + 
             106417934767940293680*\[Nu]^2) - 12*e^(24/19)*ei^4*
            (510435231187790895672403 - 2602527967770865254371160*\[Nu] + 
             1573115806240005997897200*\[Nu]^2) + 2880*e^2*
            (-646822*ei^(62/19)*(-77776941858793 - 57290524699128*\[Nu] + 
               226944607584240*\[Nu]^2) + 21*ei^2*(e*ei)^(12/19)*
              (352394775629730672441 - 2046564335232903537272*\[Nu] + 
               2213643876177691727952*\[Nu]^2))))/(12841347072734992300769280*
          Sqrt[2]*ei^(36/19))) + SO^2*x^3*
       ((-2*e^(36/19)*(723*\[Delta]*\[Kappa]A + 723*\[Kappa]S - 
            1446*\[Kappa]S*\[Nu] + 790*\[Chi]A^2 - 2892*\[Nu]*\[Chi]A^2 + 
            1580*\[Delta]*\[Chi]A*\[Chi]S + 790*\[Chi]S^2 - 
            268*\[Nu]*\[Chi]S^2) + 3*e^(12/19)*ei^(24/19)*
           (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
            191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
             \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) + 
          19*ei^(36/19)*(48*\[Delta]*\[Kappa]A + 48*\[Kappa]S - 
            96*\[Kappa]S*\[Nu] + 53*\[Chi]A^2 - 192*\[Nu]*\[Chi]A^2 + 
            106*\[Delta]*\[Chi]A*\[Chi]S + 53*\[Chi]S^2 - 
            20*\[Nu]*\[Chi]S^2))/(459648*Sqrt[2]*ei^(36/19)) + 
        (e^(12/19)*(-617421*ei^(62/19)*(178*\[Delta]*\[Kappa]A + 
             178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 
             712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*\[Chi]S + 
             191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) - 1209572*e^(62/19)*
            (-723*\[Kappa]S*(-1 + 2*\[Nu]) + 790*\[Chi]A^2 - 
             2892*\[Nu]*\[Chi]A^2 + 790*\[Chi]S^2 - 268*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(723*\[Kappa]A + 1580*\[Chi]A*\[Chi]S)) - 
           48013*e^(26/19)*ei^(36/19)*(-6564*\[Kappa]S*(-1 + 2*\[Nu]) + 
             59539*\[Chi]A^2 - 26256*\[Nu]*\[Chi]A^2 + 59539*\[Chi]S^2 - 
             211900*\[Nu]*\[Chi]S^2 + 2*\[Delta]*(3282*\[Kappa]A + 59539*
                \[Chi]A*\[Chi]S)) + 156*e^2*ei^(24/19)*
            (-3465533*\[Kappa]S*(-1 + 2*\[Nu]) + 3576925*\[Chi]A^2 - 
             13862132*\[Nu]*\[Chi]A^2 + 3576925*\[Chi]S^2 - 
             445568*\[Nu]*\[Chi]S^2 + \[Delta]*(3465533*\[Kappa]A + 7153850*
                \[Chi]A*\[Chi]S)) + 78*e^(24/19)*ei^2*
            (-9730201*\[Kappa]S*(-1 + 2*\[Nu]) + 43258141*\[Chi]A^2 - 
             38920804*\[Nu]*\[Chi]A^2 + 43258141*\[Chi]S^2 - 
             134111760*\[Nu]*\[Chi]S^2 + \[Delta]*(9730201*\[Kappa]A + 
               86516282*\[Chi]A*\[Chi]S))))/(80532781056*Sqrt[2]*
          ei^(36/19)) + (e^(12/19)*(18930724995*ei^(100/19)*
            (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
             191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
              \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) - 
           477189422528*e^(100/19)*(-723*\[Kappa]S*(-1 + 2*\[Nu]) + 
             790*\[Chi]A^2 - 2892*\[Nu]*\[Chi]A^2 + 790*\[Chi]S^2 - 
             268*\[Nu]*\[Chi]S^2 + \[Delta]*(723*\[Kappa]A + 1580*\[Chi]A*
                \[Chi]S)) - 33869952*e^2*ei^(62/19)*
            (-3465533*\[Kappa]S*(-1 + 2*\[Nu]) + 3576925*\[Chi]A^2 - 
             13862132*\[Nu]*\[Chi]A^2 + 3576925*\[Chi]S^2 - 
             445568*\[Nu]*\[Chi]S^2 + \[Delta]*(3465533*\[Kappa]A + 7153850*
                \[Chi]A*\[Chi]S)) + 49765248*e^(62/19)*ei^2*
            (-9730201*\[Kappa]S*(-1 + 2*\[Nu]) + 43258141*\[Chi]A^2 - 
             38920804*\[Nu]*\[Chi]A^2 + 43258141*\[Chi]S^2 - 
             134111760*\[Nu]*\[Chi]S^2 + \[Delta]*(9730201*\[Kappa]A + 
               86516282*\[Chi]A*\[Chi]S)) - 17332693*e^(64/19)*ei^(36/19)*
            (-18305358*\[Kappa]S*(-1 + 2*\[Nu]) + 112913425*\[Chi]A^2 - 
             73221432*\[Nu]*\[Chi]A^2 + 112913425*\[Chi]S^2 - 
             378432268*\[Nu]*\[Chi]S^2 + 2*\[Delta]*(9152679*\[Kappa]A + 
               112913425*\[Chi]A*\[Chi]S)) + 24*e^(24/19)*ei^4*
            (-2047975353262*\[Kappa]S*(-1 + 2*\[Nu]) + 2117163228955*
              \[Chi]A^2 - 8191901413048*\[Nu]*\[Chi]A^2 + 2117163228955*
              \[Chi]S^2 - 276751502772*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(2047975353262*\[Kappa]A + 4234326457910*\[Chi]A*
                \[Chi]S)) + 96*e^4*ei^(24/19)*(-2530412457698*\[Kappa]S*
              (-1 + 2*\[Nu]) + 2583790126097*\[Chi]A^2 - 10121649830792*\[Nu]*
              \[Chi]A^2 + 2583790126097*\[Chi]S^2 - 213510673596*\[Nu]*
              \[Chi]S^2 + \[Delta]*(2530412457698*\[Kappa]A + 5167580252194*
                \[Chi]A*\[Chi]S))))/(14885034988142592*Sqrt[2]*
          ei^(36/19)))) + \[Epsilon]^6*
     (SO*x^4*((Pi*(162450*ei^(48/19)*(\[Delta]*\[Chi]A + \[Chi]S - 
             4*\[Nu]*\[Chi]S) + 1885*e^(30/19)*ei^(18/19)*
            (6445*\[Delta]*\[Chi]A + (6445 + 9868*\[Nu])*\[Chi]S) + 
           2*e^(12/19)*ei^(36/19)*(31337*\[Delta]*\[Chi]A + 
             (31337 + 235316*\[Nu])*\[Chi]S) - 7*e^(48/19)*
            (1767707*\[Delta]*\[Chi]A + (1767707 + 2631716*\[Nu])*\[Chi]S)))/
         (314399232*Sqrt[2]*ei^(48/19)) - 
        (e^(12/19)*Pi*(-39*e^(36/19)*ei^2*(60799062235*\[Delta]*\[Chi]A + 
             (60799062235 - 1681766505644*\[Nu])*\[Chi]S) + 
           352812*ei^(74/19)*(31337*\[Delta]*\[Chi]A + (31337 + 235316*\[Nu])*
              \[Chi]S) + 4838288*e^(74/19)*(1767707*\[Delta]*\[Chi]A + 
             (1767707 + 2631716*\[Nu])*\[Chi]S) - 4691556*e^(26/19)*
            ei^(48/19)*(476311*\[Delta]*\[Chi]A + (476311 + 34775336*\[Nu])*
              \[Chi]S) + 52*e^2*ei^(36/19)*(65005615159*\[Delta]*\[Chi]A + 
             (65005615159 + 28798683196*\[Nu])*\[Chi]S) - 
           91*e^(56/19)*ei^(18/19)*(105398452691*\[Delta]*\[Chi]A + 
             (105398452691 + 152682705716*\[Nu])*\[Chi]S) + 
           78*e^(18/19)*ei^(56/19)*(28885111903*\[Delta]*\[Chi]A + 
             (28885111903 + 1245411730588*\[Nu])*\[Chi]S)))/
         (47215219064832*Sqrt[2]*ei^(48/19)) - 
        (e^(12/19)*Pi*(-193628285760*e^(74/19)*ei^2*
            (60799062235*\[Delta]*\[Chi]A + (60799062235 - 1681766505644*
                \[Nu])*\[Chi]S) - 63133967858325*ei^(112/19)*
            (31337*\[Delta]*\[Chi]A + (31337 + 235316*\[Nu])*\[Chi]S) + 
           14142645438958400*e^(112/19)*(1767707*\[Delta]*\[Chi]A + 
             (1767707 + 2631716*\[Nu])*\[Chi]S) - 61394874705*e^(64/19)*
            ei^(48/19)*(186875094583*\[Delta]*\[Chi]A + 
             (186875094583 + 13970833203308*\[Nu])*\[Chi]S) + 
           196560*e^(18/19)*ei^(94/19)*(3946331805444981*\[Delta]*\[Chi]A + 
             (3946331805444981 + 7919583360996556*\[Nu])*\[Chi]S) + 
           330096*e^(36/19)*ei^4*(8388764562288913*\[Delta]*\[Chi]A + 
             (8388764562288913 + 50345437127831536*\[Nu])*\[Chi]S) - 
           167440*e^(94/19)*ei^(18/19)*(175515568865658611*\[Delta]*\[Chi]A + 
             (175515568865658611 + 247878349535282036*\[Nu])*\[Chi]S) + 
           33488*e^4*ei^(36/19)*(467358287836566959*\[Delta]*\[Chi]A + 
             (467358287836566959 + 252002971654484168*\[Nu])*\[Chi]S) + 
           280140*e^2*(-235208*ei^(74/19)*(65005615159*\[Delta]*\[Chi]A + 
               65005615159*\[Chi]S + 28798683196*\[Nu]*\[Chi]S) + 
             13*ei^2*(e*ei)^(18/19)*(3493708141898143*\[Delta]*\[Chi]A + 
               3493708141898143*\[Chi]S + 140502049022356708*\[Nu]*
                \[Chi]S))))/(50932273069277388472320*Sqrt[2]*ei^(48/19))) + 
      SO^2*x^4*((108*(e*ei)^(24/19)*(-3920527 + 15455580*\[Nu])*
           (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
            191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
             \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) + 
          e^(48/19)*(-504*\[Delta]*\[Kappa]A*(-133053743 + 928541756*\[Nu]) + 
            504*\[Kappa]S*(133053743 - 1194649242*\[Nu] + 2313892152*\[Nu]^
                2) - 105124988001*\[Chi]A^2 + 272322040144*\[Delta]^2*
             \[Chi]A^2 - 139352970876*\[Nu]*\[Chi]A^2 + 2332403289216*\[Nu]^2*
             \[Chi]A^2 + 2*\[Delta]*(167197052143 - 251207045852*\[Nu])*
             \[Chi]A*\[Chi]S + 167197052143*\[Chi]S^2 + 57438831176*\[Nu]*
             \[Chi]S^2 - 75366630416*\[Nu]^2*\[Chi]S^2) - 
          16245*ei^(48/19)*(-96*\[Delta]*\[Kappa]A*(1979 + 89740*\[Nu]) + 
            96*\[Kappa]S*(-1979 - 85782*\[Nu] + 296072*\[Nu]^2) - 
            8023897*\[Chi]A^2 + 9605904*\[Delta]^2*\[Chi]A^2 + 
            19042564*\[Nu]*\[Chi]A^2 + 56845824*\[Nu]^2*\[Chi]A^2 + 
            14*\[Delta]*(226001 + 180188*\[Nu])*\[Chi]A*\[Chi]S + 
            1582007*\[Chi]S^2 + 15575656*\[Nu]*\[Chi]S^2 - 
            3065104*\[Nu]^2*\[Chi]S^2) - 98560*e^(30/19)*ei^(18/19)*
           (814673*\[Delta]^2*\[Chi]A^2 + 2*\[Delta]*(814673 + 558323*\[Nu])*
             \[Chi]A*\[Chi]S + (814673 + 1116646*\[Nu] - 1182280*\[Nu]^2)*
             \[Chi]S^2) + 19008*e^(36/19)*ei^(12/19)*(-2833 + 5516*\[Nu])*
           (-723*\[Kappa]S*(-1 + 2*\[Nu]) + 790*\[Chi]A^2 - 
            2892*\[Nu]*\[Chi]A^2 + 790*\[Chi]S^2 - 268*\[Nu]*\[Chi]S^2 + 
            \[Delta]*(723*\[Kappa]A + 1580*\[Chi]A*\[Chi]S)) - 
          352*e^(12/19)*ei^(36/19)*(-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 
              116610*\[Nu]^2) - 278950701*\[Chi]A^2 + 102215792*\[Delta]^2*
             \[Chi]A^2 + 1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*
             \[Chi]A^2 - 176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
            68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*(-124448 + 
                126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*\[Chi]A*
               \[Chi]S)))/(4648078245888*Sqrt[2]*ei^(48/19)) + 
        (e^(12/19)*(-810*e^(12/19)*ei^(62/19)*(-1859009125901 + 
             8124219720708*\[Nu])*(178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 
             356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 
             382*\[Delta]*\[Chi]A*\[Chi]S + 191*\[Chi]S^2 - 
             52*\[Nu]*\[Chi]S^2) + 1511965*e^(74/19)*
            (-504*\[Delta]*\[Kappa]A*(-133053743 + 928541756*\[Nu]) + 
             504*\[Kappa]S*(133053743 - 1194649242*\[Nu] + 2313892152*
                \[Nu]^2) - 105124988001*\[Chi]A^2 + 272322040144*\[Delta]^2*
              \[Chi]A^2 - 139352970876*\[Nu]*\[Chi]A^2 + 2332403289216*
              \[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(167197052143 - 251207045852*
                \[Nu])*\[Chi]A*\[Chi]S + 167197052143*\[Chi]S^2 + 
             57438831176*\[Nu]*\[Chi]S^2 - 75366630416*\[Nu]^2*\[Chi]S^2) - 
           280280*e^(56/19)*ei^(18/19)*(641413166941*\[Delta]^2*\[Chi]A^2 + 
             2*\[Delta]*(641413166941 + 538156238674*\[Nu])*\[Chi]A*\[Chi]S + 
             (641413166941 + 1076312477348*\[Nu] - 522485817248*\[Nu]^2)*
              \[Chi]S^2) + 6726720*e^(18/19)*ei^(56/19)*
            (4910566244*\[Delta]^2*\[Chi]A^2 + \[Delta]*(9821132488 + 
               256656092359*\[Nu])*\[Chi]A*\[Chi]S + 
             (4910566244 + 256656092359*\[Nu] - 182233298170*\[Nu]^2)*
              \[Chi]S^2) + 83160*e^(62/19)*ei^(12/19)*(-1571689321 + 
             2383893876*\[Nu])*(-723*\[Kappa]S*(-1 + 2*\[Nu]) + 
             790*\[Chi]A^2 - 2892*\[Nu]*\[Chi]A^2 + 790*\[Chi]S^2 - 
             268*\[Nu]*\[Chi]S^2 + \[Delta]*(723*\[Kappa]A + 1580*\[Chi]A*
                \[Chi]S)) + 1755*e^(50/19)*ei^(24/19)*(-3920527 + 
             15455580*\[Nu])*(-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 
             25904401*\[Chi]A^2 - 99734624*\[Nu]*\[Chi]A^2 + 
             25904401*\[Chi]S^2 - 3882980*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(24933656*\[Kappa]A + 51808802*\[Chi]A*\[Chi]S)) - 
           1389960*e^(24/19)*ei^(50/19)*(-2833 + 5516*\[Nu])*
            (-9730201*\[Kappa]S*(-1 + 2*\[Nu]) + 43258141*\[Chi]A^2 - 
             38920804*\[Nu]*\[Chi]A^2 + 43258141*\[Chi]S^2 - 
             134111760*\[Nu]*\[Chi]S^2 + \[Delta]*(9730201*\[Kappa]A + 
               86516282*\[Chi]A*\[Chi]S)) + 135832620*ei^(74/19)*
            (-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 116610*\[Nu]^2) - 
             278950701*\[Chi]A^2 + 102215792*\[Delta]^2*\[Chi]A^2 + 
             1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*\[Chi]A^2 - 
             176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
             68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*
                (-124448 + 126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*
                \[Chi]A*\[Chi]S)) - 4691556*e^(26/19)*ei^(48/19)*
            (30*\[Kappa]S*(6701235367 - 28461012890*\[Nu] + 49884658824*
                \[Nu]^2) + 75645263203*\[Chi]A^2 + 164941658112*\[Delta]^2*
              \[Chi]A^2 - 954564658552*\[Nu]*\[Chi]A^2 + 2993079529440*
              \[Nu]^2*\[Chi]A^2 + 240586921315*\[Chi]S^2 + 110028617900*\[Nu]*
              \[Chi]S^2 + 406681327360*\[Nu]^2*\[Chi]S^2 - 
             10*\[Delta]*(3*\[Kappa]A*(-6701235367 + 15058542156*\[Nu]) + 
               (-48117384263 + 54195498784*\[Nu])*\[Chi]A*\[Chi]S)) - 
           220*e^2*ei^(36/19)*(-18*\[Kappa]S*(30936092472295 - 97020967578282*
                \[Nu] + 36669175275096*\[Nu]^2) - 2827049041226343*
              \[Chi]A^2 + 2206814975744288*\[Delta]^2*\[Chi]A^2 + 
             11721265129529220*\[Nu]*\[Chi]A^2 - 1320090309903456*\[Nu]^2*
              \[Chi]A^2 - 620234065482055*\[Chi]S^2 + 1521486691850320*\[Nu]*
              \[Chi]S^2 - 141962774214352*\[Nu]^2*\[Chi]S^2 + 
             2*\[Delta]*(9*\[Kappa]A*(-30936092472295 + 35148782633692*
                  \[Nu]) + (-620234065482055 + 967277828237084*\[Nu])*\[Chi]A*
                \[Chi]S)) + 12*e^(36/19)*ei^2*(45*\[Kappa]S*
              (1288935001706461 - 4177294463778654*\[Nu] + 7918493972374440*
                \[Nu]^2) - 20556501925266276*\[Chi]A^2 + 81704163869614286*
              \[Delta]^2*\[Chi]A^2 - 79734899027872416*\[Nu]*\[Chi]A^2 + 
             712664457513699600*\[Nu]^2*\[Chi]A^2 + 61147661944348010*
              \[Chi]S^2 + 46845577848720400*\[Nu]*\[Chi]S^2 + 
             155611899223682720*\[Nu]^2*\[Chi]S^2 - 5*\[Delta]*
              (9*\[Kappa]A*(-1288935001706461 + 1599424460365732*\[Nu]) + 4*
                (-6114766194434801 + 5755766444010856*\[Nu])*\[Chi]A*
                \[Chi]S))))/(1526940184556666880*Sqrt[2]*ei^(48/19)) + 
        (e^(12/19)*(-105892920*e^(12/19)*ei^(100/19)*(-5807829430373107 + 
             10448813703644316*\[Nu])*(178*\[Delta]*\[Kappa]A + 
             178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 
             712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*\[Chi]S + 
             191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) - 368368*e^4*ei^(36/19)*
            (522*\[Delta]*\[Kappa]A*(-5919500807162614699 + 
               7069747731165174748*\[Nu]) - 522*\[Kappa]S*
              (5919500807162614699 - 18908749345490404146*\[Nu] + 
               7669345784813238456*\[Nu]^2) - 18928431938173283618721*
              \[Chi]A^2 + 16007908125685966389172*\[Delta]^2*\[Chi]A^2 + 
             78468871437192708323052*\[Nu]*\[Chi]A^2 - 8006796999345020948064*
              \[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(-2920523812487317229549 + 
               5349738482100710957140*\[Nu])*\[Chi]A*\[Chi]S - 
             2920523812487317229549*\[Chi]S^2 + 7944333279701848066112*\[Nu]*
              \[Chi]S^2 - 384534241038124715888*\[Nu]^2*\[Chi]S^2) + 
           11490899419153700*e^(112/19)*(-504*\[Delta]*\[Kappa]A*
              (-133053743 + 928541756*\[Nu]) + 504*\[Kappa]S*
              (133053743 - 1194649242*\[Nu] + 2313892152*\[Nu]^2) - 
             105124988001*\[Chi]A^2 + 272322040144*\[Delta]^2*\[Chi]A^2 - 
             139352970876*\[Nu]*\[Chi]A^2 + 2332403289216*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(167197052143 - 251207045852*\[Nu])*\[Chi]A*\[Chi]S + 
             167197052143*\[Chi]S^2 + 57438831176*\[Nu]*\[Chi]S^2 - 
             75366630416*\[Nu]^2*\[Chi]S^2) - 9386016640*e^(94/19)*ei^(18/19)*
            (168706584010919233*\[Delta]^2*\[Chi]A^2 + 2*\[Delta]*
              (168706584010919233 + 150437240685157342*\[Nu])*\[Chi]A*
              \[Chi]S + (168706584010919233 + 300874481370314684*\[Nu] - 
               100599656390439104*\[Nu]^2)*\[Chi]S^2) + 2448526080*e^(18/19)*
            ei^(94/19)*(11076157832878729*\[Delta]^2*\[Chi]A^2 + 
             2*\[Delta]*(11076157832878729 + 12262499064422347*\[Nu])*\[Chi]A*
              \[Chi]S + (11076157832878729 + 24524998128844694*\[Nu] - 
               22620318368411720*\[Nu]^2)*\[Chi]S^2) + 10095125040*e^(100/19)*
            ei^(12/19)*(-117334502622439 + 164817749118204*\[Nu])*
            (-723*\[Kappa]S*(-1 + 2*\[Nu]) + 790*\[Chi]A^2 - 
             2892*\[Nu]*\[Chi]A^2 + 790*\[Chi]S^2 - 268*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(723*\[Kappa]A + 1580*\[Chi]A*\[Chi]S)) - 
           112376160*e^(50/19)*ei^(62/19)*(-1859009125901 + 
             8124219720708*\[Nu])*(-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 
             25904401*\[Chi]A^2 - 99734624*\[Nu]*\[Chi]A^2 + 
             25904401*\[Chi]S^2 - 3882980*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(24933656*\[Kappa]A + 51808802*\[Chi]A*\[Chi]S)) - 
           51917785920*e^(62/19)*ei^(50/19)*(-1571689321 + 2383893876*\[Nu])*
            (-9730201*\[Kappa]S*(-1 + 2*\[Nu]) + 43258141*\[Chi]A^2 - 
             38920804*\[Nu]*\[Chi]A^2 + 43258141*\[Chi]S^2 - 
             134111760*\[Nu]*\[Chi]S^2 + \[Delta]*(9730201*\[Kappa]A + 
               86516282*\[Chi]A*\[Chi]S)) - 6489723240*e^(24/19)*ei^(88/19)*
            (-2833 + 5516*\[Nu])*(-2047975353262*\[Kappa]S*(-1 + 2*\[Nu]) + 
             2117163228955*\[Chi]A^2 - 8191901413048*\[Nu]*\[Chi]A^2 + 
             2117163228955*\[Chi]S^2 - 276751502772*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(2047975353262*\[Kappa]A + 4234326457910*\[Chi]A*
                \[Chi]S)) - 63197101826183325*ei^(112/19)*
            (-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 116610*\[Nu]^2) - 
             278950701*\[Chi]A^2 + 102215792*\[Delta]^2*\[Chi]A^2 + 
             1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*\[Chi]A^2 - 
             176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
             68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*
                (-124448 + 126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*
                \[Chi]A*\[Chi]S)) + 154902628608*e^(74/19)*ei^2*
            (45*\[Kappa]S*(1288935001706461 - 4177294463778654*\[Nu] + 
               7918493972374440*\[Nu]^2) - 20556501925266276*\[Chi]A^2 + 
             81704163869614286*\[Delta]^2*\[Chi]A^2 - 79734899027872416*\[Nu]*
              \[Chi]A^2 + 712664457513699600*\[Nu]^2*\[Chi]A^2 + 
             61147661944348010*\[Chi]S^2 + 46845577848720400*\[Nu]*
              \[Chi]S^2 + 155611899223682720*\[Nu]^2*\[Chi]S^2 - 
             5*\[Delta]*(9*\[Kappa]A*(-1288935001706461 + 1599424460365732*
                  \[Nu]) + 4*(-6114766194434801 + 5755766444010856*\[Nu])*
                \[Chi]A*\[Chi]S)) - 12278974941*e^(64/19)*ei^(48/19)*
            (5520*\[Kappa]S*(209905475372017 - 941420628564083*\[Nu] + 
               1512613179868518*\[Nu]^2) + 661395331763843437*\[Chi]A^2 + 
             764458173288617808*\[Delta]^2*\[Chi]A^2 - 6139672279226353948*
              \[Nu]*\[Chi]A^2 + 16699249505748438720*\[Nu]^2*\[Chi]A^2 + 
             1425853505052461245*\[Chi]S^2 + 260421436089885200*\[Nu]*
              \[Chi]S^2 + 2082778093318402480*\[Nu]^2*\[Chi]S^2 - 
             10*\[Delta]*(552*\[Kappa]A*(-209905475372017 + 521609677820049*
                  \[Nu]) + (-285170701010492249 + 323366951608109500*\[Nu])*
                \[Chi]A*\[Chi]S)) - 3588*e^(36/19)*ei^4*
            (1044*\[Kappa]S*(12865074844809000047 + 203738894213893801302*
                \[Nu] + 319207505993389781208*\[Nu]^2) - 
             269552471721024277111089*\[Chi]A^2 + 326073634913301114343504*
              \[Delta]^2*\[Chi]A^2 + 1115745822132711025804788*\[Nu]*
              \[Chi]A^2 + 666505272514197863162304*\[Nu]^2*\[Chi]A^2 + 
             56521163192276837232415*\[Chi]S^2 + 43212084697850551775000*
              \[Nu]*\[Chi]S^2 + 1270802617793555061375664*\[Nu]^2*\[Chi]S^2 + 
             2*\[Delta]*(522*\[Kappa]A*(12865074844809000047 + 
                 229469043903511801396*\[Nu]) + (56521163192276837232415 + 
                 40374009973232234567716*\[Nu])*\[Chi]A*\[Chi]S)) + 
           131105520*e^(88/19)*ei^(24/19)*(-3920527 + 15455580*\[Nu])*
            (-2948986828885*\[Kappa]S*(-1 + 2*\[Nu]) + \[Delta]*
              (2948986828885*\[Kappa]A + 6058122836776*\[Chi]A*\[Chi]S) - 
             4*((-757265354597 + 2948986828885*\[Nu])*\[Chi]A^2 + 
               (-757265354597 + 80074589503*\[Nu])*\[Chi]S^2)) + 
           4930464*e^2*(33124*ei^2*(e*ei)^(18/19)*(3866216072016148*
                \[Delta]^2*\[Chi]A^2 + \[Delta]*(7732432144032296 + 
                 203260333927841651*\[Nu])*\[Chi]A*\[Chi]S + 
               (3866216072016148 + 203260333927841651*\[Nu] - 
                 80534487366910472*\[Nu]^2)*\[Chi]S^2) + 147005*ei^(74/19)*
              (-18*\[Kappa]S*(30936092472295 - 97020967578282*\[Nu] + 
                 36669175275096*\[Nu]^2) - 2827049041226343*\[Chi]A^2 + 
               2206814975744288*\[Delta]^2*\[Chi]A^2 + 11721265129529220*
                \[Nu]*\[Chi]A^2 - 1320090309903456*\[Nu]^2*\[Chi]A^2 - 
               620234065482055*\[Chi]S^2 + 1521486691850320*\[Nu]*\[Chi]S^2 - 
               141962774214352*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*
                (9*\[Kappa]A*(-30936092472295 + 35148782633692*\[Nu]) + 
                 (-620234065482055 + 967277828237084*\[Nu])*\[Chi]A*
                  \[Chi]S)))))/(4282589248757119932306554880*Sqrt[2]*
          ei^(48/19))))
 
HlmDCMem[6, 0] = x^2*\[Epsilon]^2*
     ((5*(e^(24/19) - ei^(24/19))*(-839 + 3612*\[Nu]))/
       (1419264*Sqrt[273]*ei^(24/19)) + 
      (5*(3*ei^2*(48655927 - 184821756*\[Nu]) + 23261*e^2*
          (-839 + 3612*\[Nu]) + 722*e^(14/19)*ei^(24/19)*
          (-175141 + 651588*\[Nu])))/(14345920512*Sqrt[273]*(ei/e)^(24/19)) + 
      (5*(49*ei^4*(16680917087 - 53534704476*\[Nu]) + 1815460283*e^4*
          (-839 + 3612*\[Nu]) + 1303210*e^(52/19)*ei^(24/19)*
          (-12361241 + 45993108*\[Nu]) - 345592*e^2*ei^2*
          (-48655927 + 184821756*\[Nu])))/(718137652936704*Sqrt[273]*
        (ei/e)^(24/19))) + x^3*\[Epsilon]^4*
     ((-90*e^(24/19)*ei^(12/19)*(2376887 - 14860720*\[Nu] + 
          19923792*\[Nu]^2) - 19*ei^(36/19)*(45661561 - 255563280*\[Nu] + 
          363772080*\[Nu]^2) + e^(36/19)*(1081489489 - 6193167120*\[Nu] + 
          8704810800*\[Nu]^2))/(40772616192*Sqrt[273]*ei^(36/19)) + 
      (e^(24/19)*(302393*e^(50/19)*(1081489489 - 6193167120*\[Nu] + 
           8704810800*\[Nu]^2) + 1426672*e^(14/19)*ei^(36/19)*
          (13522442272 - 77062595925*\[Nu] + 115417999620*\[Nu]^2) - 
         735*e^2*ei^(12/19)*(159041293899 - 913298900408*\[Nu] + 
           984178414128*\[Nu]^2) + 2340*ei^(50/19)*(137842241191 - 
           791986128080*\[Nu] + 1019476806096*\[Nu]^2) - 
         168*e^(12/19)*ei^2*(118004643485102 - 672605899323435*\[Nu] + 
           1005703016654700*\[Nu]^2)))/(3571789905395712*Sqrt[273]*
        ei^(36/19)) + (e^(24/19)*(10651549610*e^(88/19)*
          (1081489489 - 6193167120*\[Nu] + 8704810800*\[Nu]^2) + 
         12600*ei^(88/19)*(47257038107471 - 243675756432400*\[Nu] + 
           295297429889616*\[Nu]^2) - 9570240*e^(50/19)*ei^2*
          (118004643485102 - 672605899323435*\[Nu] + 1005703016654700*
            \[Nu]^2) + 4952198*e^(52/19)*ei^(36/19)*(225971503921877 - 
           1288535252520960*\[Nu] + 1930944652565040*\[Nu]^2) - 
         1125*e^4*ei^(12/19)*(4278680271864451 - 24197109215592152*\[Nu] + 
           24870084486107952*\[Nu]^2) + 3600*e^2*ei^(50/19)*
          (9223243844976507 - 48292418443744984*\[Nu] + 50359242169554864*
            \[Nu]^2) - 3*e^(12/19)*ei^4*(10076003557682657627 - 
           53249972801790776760*\[Nu] + 62451968965269614640*\[Nu]^2)))/
       (58944738553044664320*Sqrt[273]*ei^(36/19))) + 
    SO*x^(7/2)*\[Epsilon]^5*
     ((5*(14*e^(24/19)*ei^(18/19)*(-839 + 3612*\[Nu])*
          (-157*\[Delta]*\[Chi]A + (-157 + 110*\[Nu])*\[Chi]S) - 
         57*ei^(42/19)*(\[Delta]*(-54403 + 446460*\[Nu])*\[Chi]A + 
           (-54403 + 66704*\[Nu] + 785904*\[Nu]^2)*\[Chi]S) + 
         e^(42/19)*(\[Delta]*(-4945093 + 33387396*\[Nu])*\[Chi]A + 
           (-4945093 + 13033364*\[Nu] + 39234048*\[Nu]^2)*\[Chi]S)))/
       (1698859008*Sqrt[273]*ei^(42/19)) + 
      (e^(24/19)*(240*ei^(56/19)*(-48655927 + 184821756*\[Nu])*
          (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) + 
         7*e^2*ei^(18/19)*(-839 + 3612*\[Nu])*(-113175949*\[Delta]*\[Chi]A + 
           (-113175949 + 41301776*\[Nu])*\[Chi]S) + 232610*e^(56/19)*
          (\[Delta]*(-4945093 + 33387396*\[Nu])*\[Chi]A + 
           (-4945093 + 13033364*\[Nu] + 39234048*\[Nu]^2)*\[Chi]S) - 
         823080*e^(14/19)*ei^(42/19)*(2*\[Delta]*(-16976455 + 
             116358648*\[Nu])*\[Chi]A + (-33952910 + 12168211*\[Nu] + 
             444154116*\[Nu]^2)*\[Chi]S) + 21*e^(18/19)*ei^2*
          (\[Delta]*(-1220333814247 + 8556007939836*\[Nu])*\[Chi]A + 
           (-1220333814247 + 87582195044*\[Nu] + 17156340864096*\[Nu]^2)*
            \[Chi]S)))/(19625219260416*Sqrt[273]*ei^(42/19)) + 
      (e^(24/19)*(25636*e^2*ei^(56/19)*(-48655927 + 184821756*\[Nu])*
          (113175949*\[Delta]*\[Chi]A + (113175949 - 41301776*\[Nu])*
            \[Chi]S) + 966280*ei^(94/19)*(-16680917087 + 53534704476*\[Nu])*
          (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) + 
         143871*e^4*ei^(18/19)*(-839 + 3612*\[Nu])*
          (-189611322769*\[Delta]*\[Chi]A + (-189611322769 + 
             49774532976*\[Nu])*\[Chi]S) + 6952069627230*e^(94/19)*
          (\[Delta]*(-4945093 + 33387396*\[Nu])*\[Chi]A + 
           (-4945093 + 13033364*\[Nu] + 39234048*\[Nu]^2)*\[Chi]S) - 
         4308412260*e^(52/19)*ei^(42/19)*(9*\[Delta]*(-36076614005 + 
             241522072148*\[Nu])*\[Chi]A + (-324689526045 + 
             177345933110*\[Nu] + 4018876732248*\[Nu]^2)*\[Chi]S) + 
         1043558243*e^(56/19)*ei^2*(\[Delta]*(-1220333814247 + 
             8556007939836*\[Nu])*\[Chi]A + (-1220333814247 + 
             87582195044*\[Nu] + 17156340864096*\[Nu]^2)*\[Chi]S) - 
         182*e^(18/19)*ei^4*(\[Delta]*(-163648399752181339 + 
             1326592427149021452*\[Nu])*\[Chi]A + (-163648399752181339 - 
             488821146525940972*\[Nu] + 3768709409950686432*\[Nu]^2)*
            \[Chi]S)))/(242164634222091829248*Sqrt[273]*ei^(42/19))) + 
    SO^2*x^4*\[Epsilon]^6*((5*(-e^(24/19) + ei^(24/19))*
        (e^(24/19)*(24*\[Delta]*\[Kappa]A*(-46297 + 179508*\[Nu]) - 
           24*\[Kappa]S*(46297 - 272102*\[Nu] + 397320*\[Nu]^2) - 
           1170697*\[Chi]A^2 + 9009156*\[Nu]*\[Chi]A^2 - 19071360*\[Nu]^2*
            \[Chi]A^2 + 2*\[Delta]*(-1170697 + 4564644*\[Nu])*\[Chi]A*
            \[Chi]S - 1170697*\[Chi]S^2 + 4802920*\[Nu]*\[Chi]S^2 - 
           2864400*\[Nu]^2*\[Chi]S^2) + 19*ei^(24/19)*
          (32*\[Delta]*\[Kappa]A*(-845 + 2856*\[Nu]) - 32*\[Kappa]S*
            (845 - 4546*\[Nu] + 7224*\[Nu]^2) - 27879*\[Chi]A^2 + 
           203164*\[Nu]*\[Chi]A^2 - 462336*\[Nu]^2*\[Chi]A^2 + 
           6*\[Delta]*(-9293 + 31668*\[Nu])*\[Chi]A*\[Chi]S - 
           27879*\[Chi]S^2 + 98360*\[Nu]*\[Chi]S^2 - 111216*\[Nu]^2*
            \[Chi]S^2)))/(862912512*Sqrt[273]*ei^(48/19)) + 
      (e^(24/19)*(-30*ei^(62/19)*(-48655927 + 184821756*\[Nu])*
          (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
           191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
            \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) + 
         54872*e^(14/19)*ei^(48/19)*(\[Delta]*\[Kappa]A*(-20323301 + 
             54435108*\[Nu]) + \[Kappa]S*(-20323301 + 95081710*\[Nu] - 
             140283528*\[Nu]^2) - 22066269*\[Chi]A^2 + 142387076*\[Nu]*
            \[Chi]A^2 - 280567056*\[Nu]^2*\[Chi]A^2 + 6*\[Delta]*
            (-7355423 + 20364624*\[Nu])*\[Chi]A*\[Chi]S - 
           22066269*\[Chi]S^2 + 68065744*\[Nu]*\[Chi]S^2 - 
           89461680*\[Nu]^2*\[Chi]S^2) + 116305*e^(62/19)*
          (-24*\[Delta]*\[Kappa]A*(-46297 + 179508*\[Nu]) + 
           24*\[Kappa]S*(46297 - 272102*\[Nu] + 397320*\[Nu]^2) + 
           1170697*\[Chi]A^2 - 9009156*\[Nu]*\[Chi]A^2 + 19071360*\[Nu]^2*
            \[Chi]A^2 + 2*\[Delta]*(1170697 - 4564644*\[Nu])*\[Chi]A*
            \[Chi]S + 1170697*\[Chi]S^2 - 4802920*\[Nu]*\[Chi]S^2 + 
           2864400*\[Nu]^2*\[Chi]S^2) + 12*e^(24/19)*ei^2*
          (\[Delta]*\[Kappa]A*(69227038931 - 162437718828*\[Nu]) + 
           \[Kappa]S*(69227038931 - 300891796690*\[Nu] + 459608148888*
              \[Nu]^2) + 75377741539*\[Chi]A^2 - 462762892736*\[Nu]*
            \[Chi]A^2 + 919216297776*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(75377741539 - 185854737012*\[Nu])*\[Chi]A*\[Chi]S + 
           75377741539*\[Chi]S^2 - 210457547444*\[Nu]*\[Chi]S^2 + 
           363133495200*\[Nu]^2*\[Chi]S^2) + 5*e^2*ei^(24/19)*
          (-839 + 3612*\[Nu])*(-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 
           25904401*\[Chi]A^2 - 99734624*\[Nu]*\[Chi]A^2 + 
           25904401*\[Chi]S^2 - 3882980*\[Nu]*\[Chi]S^2 + 
           \[Delta]*(24933656*\[Kappa]A + 51808802*\[Chi]A*\[Chi]S))))/
       (4361159835648*Sqrt[273]*ei^(48/19)) + 
      (e^(24/19)*(-10290*ei^(100/19)*(-16680917087 + 53534704476*\[Nu])*
          (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
           191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
            \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) + 
         15*e^(24/19)*ei^4*(32*\[Delta]*\[Kappa]A*(-400388252426881 + 
             802766187587868*\[Nu]) - 32*\[Kappa]S*(400388252426881 - 
             1603542692441630*\[Nu] + 3159441545844864*\[Nu]^2) - 
           13810688311230753*\[Chi]A^2 + 80248459845737732*\[Nu]*\[Chi]A^2 - 
           202204258934071296*\[Nu]^2*\[Chi]A^2 + 6*\[Delta]*
            (-4603562770410251 + 9666254511698988*\[Nu])*\[Chi]A*\[Chi]S - 
           13810688311230753*\[Chi]S^2 + 32991820469379208*\[Nu]*\[Chi]S^2 - 
           112691169051964944*\[Nu]^2*\[Chi]S^2) + 4952198*e^(52/19)*
          ei^(48/19)*(16*\[Delta]*\[Kappa]A*(-72845532337 + 
             205322523756*\[Nu]) - 16*\[Kappa]S*(72845532337 - 
             351013588430*\[Nu] + 503800610376*\[Nu]^2) - 
           1252030791783*\[Chi]A^2 + 8278732353532*\[Nu]*\[Chi]A^2 - 
           16121619532032*\[Nu]^2*\[Chi]A^2 + 6*\[Delta]*(-417343597261 + 
             1205539427988*\[Nu])*\[Chi]A*\[Chi]S - 1252030791783*\[Chi]S^2 + 
           3962627381528*\[Nu]*\[Chi]S^2 - 4306809627120*\[Nu]^2*\[Chi]S^2) + 
         331302601175*e^(100/19)*(-24*\[Delta]*\[Kappa]A*
            (-46297 + 179508*\[Nu]) + 24*\[Kappa]S*(46297 - 272102*\[Nu] + 
             397320*\[Nu]^2) + 1170697*\[Chi]A^2 - 9009156*\[Nu]*\[Chi]A^2 + 
           19071360*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(1170697 - 4564644*\[Nu])*
            \[Chi]A*\[Chi]S + 1170697*\[Chi]S^2 - 4802920*\[Nu]*\[Chi]S^2 + 
           2864400*\[Nu]^2*\[Chi]S^2) + 58059456*e^(62/19)*ei^2*
          (\[Delta]*\[Kappa]A*(69227038931 - 162437718828*\[Nu]) + 
           \[Kappa]S*(69227038931 - 300891796690*\[Nu] + 459608148888*
              \[Nu]^2) + 75377741539*\[Chi]A^2 - 462762892736*\[Nu]*
            \[Chi]A^2 + 919216297776*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(75377741539 - 185854737012*\[Nu])*\[Chi]A*\[Chi]S + 
           75377741539*\[Chi]S^2 - 210457547444*\[Nu]*\[Chi]S^2 + 
           363133495200*\[Nu]^2*\[Chi]S^2) - 1560*e^2*ei^(62/19)*
          (-48655927 + 184821756*\[Nu])*(-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 
           25904401*\[Chi]A^2 - 99734624*\[Nu]*\[Chi]A^2 + 
           25904401*\[Chi]S^2 - 3882980*\[Nu]*\[Chi]S^2 + 
           \[Delta]*(24933656*\[Kappa]A + 51808802*\[Chi]A*\[Chi]S)) + 
         140*e^4*ei^(24/19)*(-839 + 3612*\[Nu])*(-2948986828885*\[Kappa]S*
            (-1 + 2*\[Nu]) + \[Delta]*(2948986828885*\[Kappa]A + 
             6058122836776*\[Chi]A*\[Chi]S) - 
           4*((-757265354597 + 2948986828885*\[Nu])*\[Chi]A^2 + 
             (-757265354597 + 80074589503*\[Nu])*\[Chi]S^2))))/
       (4584590776347918336*Sqrt[273]*ei^(48/19)))
 
HlmDCMem[8, 0] = x^3*\[Epsilon]^4*
    (((e^(36/19) - ei^(36/19))*(75601 - 452070*\[Nu] + 733320*\[Nu]^2))/
      (213497856*Sqrt[119]*ei^(36/19)) + 
     (3323*e^2*(75601 - 452070*\[Nu] + 733320*\[Nu]^2) + 
       722*e^(2/19)*ei^(36/19)*(31645823 - 192093570*\[Nu] + 
         293380920*\[Nu]^2) - 3*ei^2*(7699835443 - 46731262050*\[Nu] + 
         71419282200*\[Nu]^2))/(205527269376*Sqrt[119]*(ei/e)^(36/19)) + 
     (409674985*e^4*(75601 - 452070*\[Nu] + 733320*\[Nu]^2) - 
       598140*e^2*ei^2*(7699835443 - 46731262050*\[Nu] + 
         71419282200*\[Nu]^2) + 130321*e^(40/19)*ei^(36/19)*
        (36088208969 - 218933563950*\[Nu] + 333741227400*\[Nu]^2) - 
       6*ei^4*(21407291285669 - 127500945901650*\[Nu] + 
         179197316814600*\[Nu]^2))/(11871255079157760*Sqrt[119]*
       (ei/e)^(36/19)))
 
HlmDCMem[l_, m_] := HlmDCMem[l, m] = 
    hideEO[(Collect[#1, {\[Epsilon], SO, x, EO}, Simplify] & )[
      showEO[cute[4][x*HDCMem[l, m] /. SO -> SO*\[Epsilon]]]]]
"DC memory expanded to e^6"
HlmDCMem[2, 0] = ((5*(e^(12/19) - ei^(12/19)))/(14*Sqrt[6]*ei^(12/19)) + 
      (5*(43199*e^2 - 361*e^(26/19)*ei^(12/19) - 42838*ei^2))/
       (525616*Sqrt[6]*(ei/e)^(12/19)) + 
      (5*(7855261648*e^4 - 68679167*e^(64/19)*ei^(12/19) - 
         9110443136*e^2*ei^2 + 1323860655*ei^4))/(97150656512*Sqrt[6]*
        (ei/e)^(12/19)) + (5*(387027635574528*e^6 - 3353289260037*e^(102/19)*
          ei^(12/19) - 440043298008416*e^4*ei^2 + 74786212261605*e^2*ei^4 - 
         18417260567680*ei^6))/(4769708632113152*Sqrt[6]*(ei/e)^(12/19)))*x + 
    x^2*\[Epsilon]^2*((5*(e^(12/19) - ei^(12/19))*
        (19*ei^(12/19)*(-4075 + 5628*\[Nu]) + e^(12/19)*
          (-145417 + 239316*\[Nu])))/(1072512*Sqrt[6]*ei^(24/19)) + 
      (5*e^(12/19)*(42*e^2*ei^(12/19)*(1082158253 - 1430751140*\[Nu]) + 
         3598392*ei^(50/19)*(-2833 + 5516*\[Nu]) + 302393*e^(50/19)*
          (-145417 + 239316*\[Nu]) + 109744*e^(26/19)*ei^(24/19)*
          (-22585 + 1216992*\[Nu]) - 9*e^(12/19)*ei^2*(-1243916559 + 
           18409137292*\[Nu])))/(140932366848*Sqrt[6]*ei^(24/19)) + 
      (5*e^(12/19)*(-120471319605*ei^(88/19)*(-2833 + 5516*\[Nu]) + 
         94403934716*e^(88/19)*(-145417 + 239316*\[Nu]) + 
         9595712*e^2*ei^(50/19)*(-1082158253 + 1430751140*\[Nu]) - 
         4147104*e^(50/19)*ei^2*(-1243916559 + 18409137292*\[Nu]) + 
         2476099*e^(64/19)*ei^(24/19)*(-2933732575 + 27297870964*\[Nu]) - 
         448*e^(12/19)*ei^4*(-23500911051568 + 15843764110629*\[Nu]) - 
         364*e^4*ei^(12/19)*(-42164487546687 + 54413587634588*\[Nu])))/
       (28219545498353664*Sqrt[6]*ei^(24/19)) + 
      (e^(12/19)*(150837364049299200*ei^(126/19)*(-2833 + 5516*\[Nu]) + 
         567006818558882520*e^(126/19)*(-145417 + 239316*\[Nu]) - 
         7089273807525*e^2*ei^(88/19)*(-1082158253 + 1430751140*\[Nu]) - 
         28570157825040*e^(88/19)*ei^2*(-1243916559 + 18409137292*\[Nu]) - 
         4555434240*e^(50/19)*ei^4*(-23500911051568 + 15843764110629*\[Nu]) + 
         1835179920*e^4*ei^(50/19)*(-42164487546687 + 54413587634588*\[Nu]) + 
         893871739*e^(102/19)*ei^(24/19)*(-68594457788131 + 
           531584534617500*\[Nu]) - 65520*e^6*ei^(12/19)*
          (-1481083213734247525 + 1878694147208499348*\[Nu]) + 
         9464*e^(12/19)*ei^6*(-2719168921326953119 + 2093118343804025580*
            \[Nu])))/(24938402504169086779392*Sqrt[6]*ei^(24/19))) + 
    SO*x^(7/2)*\[Epsilon]^5*
     ((560*e^(24/19)*ei^(18/19)*(-145417 + 239316*\[Nu])*
         (-157*\[Delta]*\[Chi]A + (-157 + 110*\[Nu])*\[Chi]S) + 
        105*e^(30/19)*ei^(12/19)*(-2833 + 5516*\[Nu])*
         (-99391*\[Delta]*\[Chi]A + (-99391 + 22844*\[Nu])*\[Chi]S) + 
        2*e^(42/19)*(\[Delta]*(7875278603 + 32255910204*\[Nu])*\[Chi]A + 
          (7875278603 + 99877088*\[Nu] - 7901609520*\[Nu]^2)*\[Chi]S) + 
        1344*e^(12/19)*ei^(30/19)*(4*\[Delta]*(-5883581 + 3034822*\[Nu])*
           \[Chi]A + (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*
           \[Chi]S) - 5415*ei^(42/19)*(\[Delta]*(4888427 + 410172*\[Nu])*
           \[Chi]A + (4888427 - 7631920*\[Nu] + 74256*\[Nu]^2)*\[Chi]S))/
       (10270374912*Sqrt[6]*ei^(42/19)) + 
      (e^(12/19)*(5760*e^(12/19)*ei^(56/19)*(-1243916559 + 18409137292*\[Nu])*
          (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) + 
         420*e^(56/19)*ei^(12/19)*(-724653277 + 1072804096*\[Nu])*
          (-99391*\[Delta]*\[Chi]A + (-99391 + 22844*\[Nu])*\[Chi]S) + 
         728*e^(50/19)*ei^(18/19)*(-145417 + 239316*\[Nu])*
          (-113175949*\[Delta]*\[Chi]A + (-113175949 + 41301776*\[Nu])*
            \[Chi]S) + 1365*e^(18/19)*ei^(50/19)*(-2833 + 5516*\[Nu])*
          (1171697071*\[Delta]*\[Chi]A + (1171697071 + 6870205636*\[Nu])*
            \[Chi]S) + 168*e^2*ei^(30/19)*(\[Delta]*(-118659872446897 + 
             113510738272204*\[Nu])*\[Chi]A + (-118659872446897 + 
             286545005193108*\[Nu] - 61303880674400*\[Nu]^2)*\[Chi]S) + 
         1209572*e^(68/19)*(\[Delta]*(7875278603 + 32255910204*\[Nu])*
            \[Chi]A + (7875278603 + 99877088*\[Nu] - 7901609520*\[Nu]^2)*
            \[Chi]S) - 230297088*ei^(68/19)*
          (4*\[Delta]*(-5883581 + 3034822*\[Nu])*\[Chi]A + 
           (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*\[Chi]S) - 
         1954815*e^(26/19)*ei^(42/19)*(\[Delta]*(3593809211 + 
             51319330860*\[Nu])*\[Chi]A + (3593809211 - 4349305600*\[Nu] + 
             75051987024*\[Nu]^2)*\[Chi]S) + 504*e^(30/19)*ei^2*
          (7*\[Delta]*(-6961707660527 + 23826828158692*\[Nu])*\[Chi]A + 
           (-48731953623689 + 56400792295184*\[Nu] + 212302080351088*\[Nu]^2)*
            \[Chi]S)))/(1542363822784512*Sqrt[6]*ei^(42/19)) + 
      (e^(12/19)*(28302144*e^(50/19)*ei^(56/19)*(-1243916559 + 
           18409137292*\[Nu])*(113175949*\[Delta]*\[Chi]A + 
           (113175949 - 41301776*\[Nu])*\[Chi]S) + 812779520*e^(12/19)*
          ei^(94/19)*(-23500911051568 + 15843764110629*\[Nu])*
          (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) + 
         5159245*e^(94/19)*ei^(12/19)*(-96101931245055 + 132465768779708*
            \[Nu])*(-99391*\[Delta]*\[Chi]A + (-99391 + 22844*\[Nu])*
            \[Chi]S) + 20636980*e^(56/19)*ei^(50/19)*(-724653277 + 
           1072804096*\[Nu])*(1171697071*\[Delta]*\[Chi]A + 
           (1171697071 + 6870205636*\[Nu])*\[Chi]S) + 688278864*e^(88/19)*
          ei^(18/19)*(-145417 + 239316*\[Nu])*(-189611322769*\[Delta]*
            \[Chi]A + (-189611322769 + 49774532976*\[Nu])*\[Chi]S) + 
         1206660*e^(18/19)*ei^(88/19)*(-2833 + 5516*\[Nu])*
          (399181778038271*\[Delta]*\[Chi]A + (399181778038271 + 
             196461331296836*\[Nu])*\[Chi]S) + 91*e^(30/19)*ei^4*
          (\[Delta]*(-67851876903305430734671 + 8143200327587771321148*\[Nu])*
            \[Chi]A + (-67851876903305430734671 + 280784879015624996422712*
              \[Nu] - 169173767949134972320752*\[Nu]^2)*\[Chi]S) + 
         24752*e^4*ei^(30/19)*(\[Delta]*(-881417005050492936489 + 
             1269739664742855690748*\[Nu])*\[Chi]A + 
           (-881417005050492936489 + 2871898607753603335636*\[Nu] - 
             610400986757547833120*\[Nu]^2)*\[Chi]S) - 108805778368*e^2*
          ei^(68/19)*(\[Delta]*(-118659872446897 + 113510738272204*\[Nu])*
            \[Chi]A + (-118659872446897 + 286545005193108*\[Nu] - 
             61303880674400*\[Nu]^2)*\[Chi]S) + 1662935054833416*e^(106/19)*
          (\[Delta]*(7875278603 + 32255910204*\[Nu])*\[Chi]A + 
           (7875278603 + 99877088*\[Nu] - 7901609520*\[Nu]^2)*\[Chi]S) + 
         21856388688017520*ei^(106/19)*(4*\[Delta]*(-5883581 + 3034822*\[Nu])*
            \[Chi]A + (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*
            \[Chi]S) + 1152088300272*e^(68/19)*ei^2*
          (7*\[Delta]*(-6961707660527 + 23826828158692*\[Nu])*\[Chi]A + 
           (-48731953623689 + 56400792295184*\[Nu] + 212302080351088*\[Nu]^2)*
            \[Chi]S) - 163719665880*e^(64/19)*ei^(42/19)*
          (\[Delta]*(-99287898431985 + 1565517614766014*\[Nu])*\[Chi]A + 
           5*(-19857579686397 + 6993230143212*\[Nu] + 441867124045984*
              \[Nu]^2)*\[Chi]S)))/(875469182646356888518656*Sqrt[6]*
        ei^(42/19)) + (e^(12/19)*(4901150592*e^(88/19)*ei^(56/19)*
          (-1243916559 + 18409137292*\[Nu])*(189611322769*\[Delta]*\[Chi]A + 
           (189611322769 - 49774532976*\[Nu])*\[Chi]S) + 
         731501568*e^(50/19)*ei^(94/19)*(-23500911051568 + 
           15843764110629*\[Nu])*(113175949*\[Delta]*\[Chi]A + 
           (113175949 - 41301776*\[Nu])*\[Chi]S) + 403999232*e^(12/19)*
          ei^(132/19)*(-2719168921326953119 + 2093118343804025580*\[Nu])*
          (-157*\[Delta]*\[Chi]A + (-157 + 110*\[Nu])*\[Chi]S) + 
         781170390*e^(132/19)*ei^(12/19)*(-137799513225092197 + 
           183838168203252052*\[Nu])*(-99391*\[Delta]*\[Chi]A + 
           (-99391 + 22844*\[Nu])*\[Chi]S) + 46433205*e^(94/19)*ei^(50/19)*
          (-96101931245055 + 132465768779708*\[Nu])*
          (1171697071*\[Delta]*\[Chi]A + (1171697071 + 6870205636*\[Nu])*
            \[Chi]S) + 3341520*e^(56/19)*ei^(88/19)*(-724653277 + 
           1072804096*\[Nu])*(399181778038271*\[Delta]*\[Chi]A + 
           (399181778038271 + 196461331296836*\[Nu])*\[Chi]S) - 
         118942200*e^(18/19)*ei^(126/19)*(-2833 + 5516*\[Nu])*
          (371951973273312503*\[Delta]*\[Chi]A + (371951973273312503 + 
             78429579789035624*\[Nu])*\[Chi]S) + 5550636*e^(126/19)*
          ei^(18/19)*(-145417 + 239316*\[Nu])*(-5038160345051358953*\[Delta]*
            \[Chi]A + (-5038160345051358953 + 1065846237840519832*\[Nu])*
            \[Chi]S) - 13*e^(30/19)*ei^6*
          (7*\[Delta]*(-25119013415045732523094873211 + 
             12028278989858245813112009388*\[Nu])*\[Chi]A + 
           (-175833093905320127661664112477 + 169281618499208901308220940984*
              \[Nu] - 17359652238135706399900241424*\[Nu]^2)*\[Chi]S) + 
         936*e^6*ei^(30/19)*(7*\[Delta]*(-532202728837787049264046867 + 
             1052272546464970879844248804*\[Nu])*\[Chi]A + 
           (-3725419101864509344848328069 + 15569484378055550886037435956*
              \[Nu] - 3304589732988741434471430560*\[Nu]^2)*\[Chi]S) + 
         38101518*e^(68/19)*ei^4*(\[Delta]*(-67851876903305430734671 + 
             8143200327587771321148*\[Nu])*\[Chi]A + 
           (-67851876903305430734671 + 280784879015624996422712*\[Nu] - 
             169173767949134972320752*\[Nu]^2)*\[Chi]S) - 
         2936287872*e^4*ei^(68/19)*(\[Delta]*(-881417005050492936489 + 
             1269739664742855690748*\[Nu])*\[Chi]A + 
           (-881417005050492936489 + 2871898607753603335636*\[Nu] - 
             610400986757547833120*\[Nu]^2)*\[Chi]S) + 1891418251847670*e^2*
          ei^(106/19)*(\[Delta]*(-118659872446897 + 113510738272204*\[Nu])*
            \[Chi]A + (-118659872446897 + 286545005193108*\[Nu] - 
             61303880674400*\[Nu]^2)*\[Chi]S) + 331547386563424279248*
          e^(144/19)*(\[Delta]*(7875278603 + 32255910204*\[Nu])*\[Chi]A + 
           (7875278603 + 99877088*\[Nu] - 7901609520*\[Nu]^2)*\[Chi]S) - 
         643894539653648424960*ei^(144/19)*
          (4*\[Delta]*(-5883581 + 3034822*\[Nu])*\[Chi]A + 
           (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*\[Chi]S) + 
         290118208027860576*e^(106/19)*ei^2*
          (7*\[Delta]*(-6961707660527 + 23826828158692*\[Nu])*\[Chi]A + 
           (-48731953623689 + 56400792295184*\[Nu] + 212302080351088*\[Nu]^2)*
            \[Chi]S) - 985046656378*e^(102/19)*ei^(42/19)*
          (\[Delta]*(-7276770611998778183 + 68408938744352000340*\[Nu])*
            \[Chi]A + (-7276770611998778183 + 3788968440211405420*\[Nu] + 
             95499743546993341920*\[Nu]^2)*\[Chi]S)))/
       (91020779981376432985507627008*Sqrt[6]*ei^(42/19))) + 
    \[Epsilon]^3*(((1885*(-e^(30/19) + e^(12/19)*ei^(18/19))*Pi)/
         (3192*Sqrt[6]*ei^(30/19)) - (5*e^(12/19)*(25055420*e^(56/19) - 
           27470894*e^2*ei^(18/19) - 925965*e^(26/19)*ei^(30/19) - 
           6596977*e^(18/19)*ei^2 + 9938416*ei^(56/19))*Pi)/
         (73747968*Sqrt[6]*ei^(30/19)) + 
        (5*e^(12/19)*(-9269831689358240*e^(94/19) + 9995224240679288*e^4*
            ei^(18/19) + 530147328125031*e^(64/19)*ei^(30/19) + 
           3801670677702820*e^(56/19)*ei^2 - 6279394966669792*e^2*
            ei^(56/19) + 889287733075248*e^(18/19)*ei^4 + 
           332896676445645*ei^(94/19))*Pi)/(14774283239718912*Sqrt[6]*
          ei^(30/19)) - (5*e^(12/19)*(59012693533231186607040*e^(132/19) - 
           63234058261818572623222*e^6*ei^(18/19) - 4028842426764406598817*
            e^(102/19)*ei^(30/19) - 28954341276045863517840*e^(94/19)*ei^2 + 
           47033422621517632044864*e^4*ei^(56/19) - 10549718199122305301280*
            e^(56/19)*ei^4 - 4329912525514387226415*e^2*ei^(94/19) + 
           4661736916809304025590*e^(18/19)*ei^6 + 389019617707412590080*
            ei^(132/19))*Pi)/(60930089634728135098368*Sqrt[6]*ei^(30/19)))*
       x^(5/2) + SO*x^(5/2)*
       ((e^(30/19)*(99391*\[Delta]*\[Chi]A + (99391 - 22844*\[Nu])*\[Chi]S) - 
          160*e^(12/19)*ei^(18/19)*(157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*
             \[Chi]S) - 57*ei^(30/19)*(1303*\[Delta]*\[Chi]A + 
            (1303 - 92*\[Nu])*\[Chi]S))/(76608*Sqrt[6]*ei^(30/19)) + 
        (e^(12/19)*(-208*e^2*ei^(18/19)*(92307509*\[Delta]*\[Chi]A + 
             (92307509 - 26680576*\[Nu])*\[Chi]S) - 308655*e^(26/19)*
            ei^(30/19)*(41497*\[Delta]*\[Chi]A + (41497 - 313124*\[Nu])*
              \[Chi]S) + 431990*e^(56/19)*(99391*\[Delta]*\[Chi]A + 
             (99391 - 22844*\[Nu])*\[Chi]S) + 27416320*ei^(56/19)*
            (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) - 
           13*e^(18/19)*ei^2*(1171697071*\[Delta]*\[Chi]A + 
             (1171697071 + 6870205636*\[Nu])*\[Chi]S)))/(11504683008*Sqrt[6]*
          ei^(30/19)) + (e^(12/19)*(-23920*e^4*ei^(18/19)*
            (6334364140217*\[Delta]*\[Chi]A + (6334364140217 - 1267078023544*
                \[Nu])*\[Chi]S) - 29081782755*e^(64/19)*ei^(30/19)*
            (1225909*\[Delta]*\[Chi]A + (1225909 - 36705500*\[Nu])*\[Chi]S) + 
           914334272*e^2*ei^(56/19)*(92307509*\[Delta]*\[Chi]A + 
             (92307509 - 26680576*\[Nu])*\[Chi]S) + 3073551621140*e^(94/19)*
            (99391*\[Delta]*\[Chi]A + (99391 - 22844*\[Nu])*\[Chi]S) - 
           17660301137700*ei^(94/19)*(157*\[Delta]*\[Chi]A + 
             (157 - 110*\[Nu])*\[Chi]S) - 144068665*e^(56/19)*ei^2*
            (1171697071*\[Delta]*\[Chi]A + (1171697071 + 6870205636*\[Nu])*
              \[Chi]S) - 78*e^(18/19)*ei^4*(399181778038271*\[Delta]*
              \[Chi]A + (399181778038271 + 196461331296836*\[Nu])*\[Chi]S)))/
         (44322849719156736*Sqrt[6]*ei^(30/19)) + 
        (5*e^(12/19)*(-71162*e^6*ei^(18/19)*(477199721661585533*\[Delta]*
              \[Chi]A + (477199721661585533 - 76702103546620216*\[Nu])*
              \[Chi]S) - 77766841293*e^(102/19)*ei^(30/19)*
            (18182549065*\[Delta]*\[Chi]A + (18182549065 - 3475155620576*
                \[Nu])*\[Chi]S) + 3751923392*e^4*ei^(56/19)*
            (6334364140217*\[Delta]*\[Chi]A + (6334364140217 - 1267078023544*
                \[Nu])*\[Chi]S) - 21015758353863*e^2*ei^(94/19)*
            (92307509*\[Delta]*\[Chi]A + (92307509 - 26680576*\[Nu])*
              \[Chi]S) + 652218098289469348*e^(132/19)*
            (99391*\[Delta]*\[Chi]A + (99391 - 22844*\[Nu])*\[Chi]S) + 
           687921516723983360*ei^(132/19)*(157*\[Delta]*\[Chi]A + 
             (157 - 110*\[Nu])*\[Chi]S) - 36575264291566*e^(94/19)*ei^2*
            (1171697071*\[Delta]*\[Chi]A + (1171697071 + 6870205636*\[Nu])*
              \[Chi]S) - 30844086*e^(56/19)*ei^4*(399181778038271*\[Delta]*
              \[Chi]A + (399181778038271 + 196461331296836*\[Nu])*\[Chi]S) + 
           10166*e^(18/19)*ei^6*(371951973273312503*\[Delta]*\[Chi]A + 
             (371951973273312503 + 78429579789035624*\[Nu])*\[Chi]S)))/
         (30465044817364067549184*Sqrt[6]*ei^(30/19)))) + 
    \[Epsilon]^4*(x^3*((5*(e^(12/19) - ei^(12/19))*
          (361*ei^(24/19)*(-151877213 - 187208280*\[Nu] + 39054960*\[Nu]^2) - 
           (e*ei)^(12/19)*(31176362099 + 92144559528*\[Nu] + 
             14676980976*\[Nu]^2) + e^(24/19)*(50392977379 - 
             385204834728*\[Nu] + 246696296112*\[Nu]^2)))/
         (338922372096*Sqrt[6]*ei^(36/19)) + 
        (5*e^(12/19)*(-6597052*ei^(62/19)*(-358353209 + 372157128*\[Nu] + 
             435997296*\[Nu]^2) + 302393*e^(62/19)*(50392977379 - 
             385204834728*\[Nu] + 246696296112*\[Nu]^2) + 260642*e^(26/19)*
            ei^(36/19)*(1023823302491 - 4774321617240*\[Nu] + 
             2523333592080*\[Nu]^2) - 1617*e^(50/19)*ei^(12/19)*
            (27565325190597 - 84987315736504*\[Nu] + 65207541903504*
              \[Nu]^2) + 1188*e^(12/19)*ei^(50/19)*(3524015611647 - 
             59014529687680*\[Nu] + 101544801302672*\[Nu]^2) + 
           132*e^2*ei^(24/19)*(-77776941858793 - 57290524699128*\[Nu] + 
             226944607584240*\[Nu]^2) - 32*e^(24/19)*ei^2*(7306273592457739 - 
             40736695045546269*\[Nu] + 24205030226070876*\[Nu]^2)))/
         (29690503588601856*Sqrt[6]*ei^(36/19)) + 
        (e^(12/19)*(-1797478588254148800*ei^(138/19)*(-358353209 + 
             372157128*\[Nu] + 435997296*\[Nu]^2) + 6421601321727355400*
            e^(138/19)*(50392977379 - 385204834728*\[Nu] + 246696296112*
              \[Nu]^2) + 144823736353725*e^2*ei^(100/19)*(-77776941858793 - 
             57290524699128*\[Nu] + 226944607584240*\[Nu]^2) - 
           1355899659154560*e^(100/19)*ei^2*(7306273592457739 - 
             40736695045546269*\[Nu] + 24205030226070876*\[Nu]^2) - 
           34606249920*e^4*ei^(62/19)*(-8259865818436428113 - 
             2803675687047024972*\[Nu] + 16063790724563400000*\[Nu]^2) + 
           16983563041*e^(102/19)*ei^(36/19)*(588408877431430301353 - 
             3010857804240096320328*\[Nu] + 1547516078643914175888*\[Nu]^2) + 
           158336640*e^(50/19)*ei^(88/19)*(4454845412928108978288 - 
             9406754868738873532081*\[Nu] + 4317023985663613391076*\[Nu]^2) - 
           8120112*e^(12/19)*ei^(126/19)*(7703405554119258186127 - 
             20928740038036277872544*\[Nu] + 11545640784423005099280*
              \[Nu]^2) - 39656682*e^(62/19)*ei^4*(-15782029726452124264043 - 
             43183078773123411060120*\[Nu] + 21685168802452643835984*
              \[Nu]^2) + 59073300*e^(88/19)*ei^(50/19)*
            (6343648677996200717931 - 102446648821941436404392*\[Nu] + 
             126754374243743232309232*\[Nu]^2) - 1171170*e^(126/19)*
            ei^(12/19)*(1019041603402746469825345 - 3015202618891145351861752*
              \[Nu] + 2202211984888771992152016*\[Nu]^2) + 
           1144*e^6*ei^(24/19)*(-497250243303835832331909187 - 
             133308866202460814417262984*\[Nu] + 849484252524317627474122320*
              \[Nu]^2) - 4*e^(24/19)*ei^6*(141655716761937010685634349403 - 
             141460070026150986951177051624*\[Nu] + 
             33344875954750833932647734000*\[Nu]^2)))/
         (34149850360329014024922464256*Sqrt[6]*ei^(36/19)) + 
        (e^(12/19)*(12*e^(24/19)*ei^4*(15782029726452124264043 + 
             43183078773123411060120*\[Nu] - 21685168802452643835984*
              \[Nu]^2) + 6625922578275*ei^(100/19)*(-358353209 + 
             372157128*\[Nu] + 435997296*\[Nu]^2) + 3877164058040*e^(100/19)*
            (50392977379 - 385204834728*\[Nu] + 246696296112*\[Nu]^2) - 
           663536640*e^(62/19)*ei^2*(7306273592457739 - 40736695045546269*
              \[Nu] + 24205030226070876*\[Nu]^2) + 329321167*e^(64/19)*
            ei^(36/19)*(15334587311436637 - 76415710261605480*\[Nu] + 
             39628755350153424*\[Nu]^2) + 1774080*e^(12/19)*ei^(88/19)*
            (66578081009092144 - 174516409085861045*\[Nu] + 87394202834229564*
              \[Nu]^2) - 900900*e^(88/19)*ei^(12/19)*(741588616321469453 - 
             2221703929192049176*\[Nu] + 1647787690719105936*\[Nu]^2) + 
           31680*e^4*ei^(24/19)*(-8259865818436428113 - 2803675687047024972*
              \[Nu] + 16063790724563400000*\[Nu]^2) + 
           10560*e^2*(-85676*ei^(62/19)*(-77776941858793 - 57290524699128*
                \[Nu] + 226944607584240*\[Nu]^2) + 63*ei^2*(e*ei)^(12/19)*
              (235797495882898419 - 3828581745994024968*\[Nu] + 
               5016023129984824048*\[Nu]^2))))/(35670408535374978613248*
          Sqrt[6]*ei^(36/19))) + SO^2*x^3*
       ((5*(2*e^(12/19)*ei^(24/19)*(178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 
             356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 
             382*\[Delta]*\[Chi]A*\[Chi]S + 191*\[Chi]S^2 - 
             52*\[Nu]*\[Chi]S^2) + 19*ei^(36/19)*(32*\[Delta]*\[Kappa]A + 
             32*\[Kappa]S - 64*\[Kappa]S*\[Nu] + 33*\[Chi]A^2 - 
             128*\[Nu]*\[Chi]A^2 + 66*\[Delta]*\[Chi]A*\[Chi]S + 
             33*\[Chi]S^2 - 4*\[Nu]*\[Chi]S^2) + e^(36/19)*
            (964*\[Kappa]S*(-1 + 2*\[Nu]) - 1009*\[Chi]A^2 + 
             3856*\[Nu]*\[Chi]A^2 - 1009*\[Chi]S^2 + 180*\[Nu]*\[Chi]S^2 - 
             2*\[Delta]*(482*\[Kappa]A + 1009*\[Chi]A*\[Chi]S))))/
         (8512*Sqrt[6]*ei^(36/19)) - 
        (5*e^(12/19)*(599732*ei^(62/19)*(178*\[Delta]*\[Kappa]A + 
             178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 
             712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*\[Chi]S + 
             191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) - 13*e^(24/19)*ei^2*
            (10216736*\[Kappa]S*(-1 + 2*\[Nu]) + 50378283*\[Chi]A^2 + 
             40866944*\[Nu]*\[Chi]A^2 + 50378283*\[Chi]S^2 - 
             242380076*\[Nu]*\[Chi]S^2 - 2*\[Delta]*(5108368*\[Kappa]A - 
               50378283*\[Chi]A*\[Chi]S)) + 192052*e^(26/19)*ei^(36/19)*
            (2986*\[Kappa]S*(-1 + 2*\[Nu]) + 953*\[Chi]A^2 + 
             11944*\[Nu]*\[Chi]A^2 + 953*\[Chi]S^2 - 15756*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(-2986*\[Kappa]A + 1906*\[Chi]A*\[Chi]S)) + 
           907179*e^(62/19)*(-964*\[Kappa]S*(-1 + 2*\[Nu]) + 1009*\[Chi]A^2 - 
             3856*\[Nu]*\[Chi]A^2 + 1009*\[Chi]S^2 - 180*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(964*\[Kappa]A + 2018*\[Chi]A*\[Chi]S)) - 
           156*e^2*ei^(24/19)*(-3465533*\[Kappa]S*(-1 + 2*\[Nu]) + 
             3576925*\[Chi]A^2 - 13862132*\[Nu]*\[Chi]A^2 + 
             3576925*\[Chi]S^2 - 445568*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(3465533*\[Kappa]A + 7153850*\[Chi]A*\[Chi]S))))/
         (2237021696*Sqrt[6]*ei^(36/19)) + (Sqrt[3/2]*e^(12/19)*
          (15445040975*ei^(100/19)*(178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 
             356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 
             382*\[Delta]*\[Chi]A*\[Chi]S + 191*\[Chi]S^2 - 
             52*\[Nu]*\[Chi]S^2) + 6911840*e^(62/19)*ei^2*
            (10216736*\[Kappa]S*(-1 + 2*\[Nu]) + 50378283*\[Chi]A^2 + 
             40866944*\[Nu]*\[Chi]A^2 + 50378283*\[Chi]S^2 - 
             242380076*\[Nu]*\[Chi]S^2 - 2*\[Delta]*(5108368*\[Kappa]A - 
               50378283*\[Chi]A*\[Chi]S)) + 51998079*e^(64/19)*ei^(36/19)*
            (-4021854*\[Kappa]S*(-1 + 2*\[Nu]) - 3794671*\[Chi]A^2 - 
             16087416*\[Nu]*\[Chi]A^2 - 3794671*\[Chi]S^2 + 
             31266100*\[Nu]*\[Chi]S^2 + \[Delta]*(4021854*\[Kappa]A - 7589342*
                \[Chi]A*\[Chi]S)) - 298243389080*e^(100/19)*
            (-964*\[Kappa]S*(-1 + 2*\[Nu]) + 1009*\[Chi]A^2 - 
             3856*\[Nu]*\[Chi]A^2 + 1009*\[Chi]S^2 - 180*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(964*\[Kappa]A + 2018*\[Chi]A*\[Chi]S)) - 
           27416320*e^2*ei^(62/19)*(-3465533*\[Kappa]S*(-1 + 2*\[Nu]) + 
             3576925*\[Chi]A^2 - 13862132*\[Nu]*\[Chi]A^2 + 
             3576925*\[Chi]S^2 - 445568*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(3465533*\[Kappa]A + 7153850*\[Chi]A*\[Chi]S)) + 
           80*e^4*ei^(24/19)*(-2530412457698*\[Kappa]S*(-1 + 2*\[Nu]) + 
             2583790126097*\[Chi]A^2 - 10121649830792*\[Nu]*\[Chi]A^2 + 
             2583790126097*\[Chi]S^2 - 213510673596*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(2530412457698*\[Kappa]A + 5167580252194*\[Chi]A*
                \[Chi]S)) + 8*e^(24/19)*ei^4*(-4853042174758*\[Kappa]S*
              (-1 + 2*\[Nu]) + 4806182243753*\[Chi]A^2 - 19412168699032*\[Nu]*
              \[Chi]A^2 + 4806182243753*\[Chi]S^2 + 187439724020*\[Nu]*
              \[Chi]S^2 + \[Delta]*(4853042174758*\[Kappa]A + 9612364487506*
                \[Chi]A*\[Chi]S))))/(206736597057536*ei^(36/19)) + 
        (e^(12/19)*(-12569780337441600*ei^(138/19)*(178*\[Delta]*\[Kappa]A + 
             178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 
             712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*\[Chi]S + 
             191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) + 18771306519*e^(102/19)*
            ei^(36/19)*(-59274639867*\[Kappa]S*(-1 + 2*\[Nu]) - 
             73851835093*\[Chi]A^2 - 237098559468*\[Nu]*\[Chi]A^2 - 
             73851835093*\[Chi]S^2 + 532505899840*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(59274639867*\[Kappa]A - 147703670186*\[Chi]A*
                \[Chi]S)) + 42371864348580*e^(100/19)*ei^2*
            (10216736*\[Kappa]S*(-1 + 2*\[Nu]) + 50378283*\[Chi]A^2 + 
             40866944*\[Nu]*\[Chi]A^2 + 50378283*\[Chi]S^2 - 
             242380076*\[Nu]*\[Chi]S^2 - 2*\[Delta]*(5108368*\[Kappa]A - 
               50378283*\[Chi]A*\[Chi]S)) - 1481907997321697400*e^(138/19)*
            (-964*\[Kappa]S*(-1 + 2*\[Nu]) + 1009*\[Chi]A^2 - 
             3856*\[Nu]*\[Chi]A^2 + 1009*\[Chi]S^2 - 180*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(964*\[Kappa]A + 2018*\[Chi]A*\[Chi]S)) + 
           13165794213975*e^2*ei^(100/19)*(-3465533*\[Kappa]S*
              (-1 + 2*\[Nu]) + 3576925*\[Chi]A^2 - 13862132*\[Nu]*\[Chi]A^2 + 
             3576925*\[Chi]S^2 - 445568*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(3465533*\[Kappa]A + 7153850*\[Chi]A*\[Chi]S)) - 
           262168560*e^4*ei^(62/19)*(-2530412457698*\[Kappa]S*
              (-1 + 2*\[Nu]) + 2583790126097*\[Chi]A^2 - 10121649830792*\[Nu]*
              \[Chi]A^2 + 2583790126097*\[Chi]S^2 - 213510673596*\[Nu]*
              \[Chi]S^2 + \[Delta]*(2530412457698*\[Kappa]A + 5167580252194*
                \[Chi]A*\[Chi]S)) + 79313364*e^(62/19)*ei^4*
            (-4853042174758*\[Kappa]S*(-1 + 2*\[Nu]) + 4806182243753*
              \[Chi]A^2 - 19412168699032*\[Nu]*\[Chi]A^2 + 4806182243753*
              \[Chi]S^2 + 187439724020*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(4853042174758*\[Kappa]A + 9612364487506*\[Chi]A*
                \[Chi]S)) + 7072*e^6*ei^(24/19)*(-151862337793381386*
              \[Kappa]S*(-1 + 2*\[Nu]) + 154258072076881051*\[Chi]A^2 - 
             607449351173525544*\[Nu]*\[Chi]A^2 + 154258072076881051*
              \[Chi]S^2 - 9582937133998660*\[Nu]*\[Chi]S^2 + 
             2*\[Delta]*(75931168896690693*\[Kappa]A + 154258072076881051*
                \[Chi]A*\[Chi]S)) - 272*e^(24/19)*ei^6*
            (-331165474334851236*\[Kappa]S*(-1 + 2*\[Nu]) + 
             340045705042503011*\[Chi]A^2 - 1324661897339404944*\[Nu]*
              \[Chi]A^2 + 340045705042503011*\[Chi]S^2 - 35520922830607100*
              \[Nu]*\[Chi]S^2 + \[Delta]*(331165474334851236*\[Kappa]A + 
               680091410085006022*\[Chi]A*\[Chi]S))))/(197923829398167355392*
          Sqrt[6]*ei^(36/19)))) + \[Epsilon]^6*
     (SO*x^4*((5*Pi*(-7*e^(48/19)*(6762971*\[Delta]*\[Chi]A + 
             (6762971 - 1994668*\[Nu])*\[Chi]S) + 377*e^(30/19)*ei^(18/19)*
            (124511*\[Delta]*\[Chi]A + (124511 - 40444*\[Nu])*\[Chi]S) + 
           149454*ei^(48/19)*(\[Delta]*\[Chi]A + \[Chi]S - 4*\[Nu]*\[Chi]S) + 
           8*e^(12/19)*ei^(36/19)*(31337*\[Delta]*\[Chi]A + 
             (31337 + 235316*\[Nu])*\[Chi]S)))/(34933248*Sqrt[6]*
          ei^(48/19)) - (e^(12/19)*Pi*(-13*e^(36/19)*ei^2*
            (5132086523719*\[Delta]*\[Chi]A + (5132086523719 - 7005190158236*
                \[Nu])*\[Chi]S) - 455*e^(56/19)*ei^(18/19)*
            (408268760993*\[Delta]*\[Chi]A + (408268760993 - 111244601092*
                \[Nu])*\[Chi]S) + 24191440*e^(74/19)*
            (6762971*\[Delta]*\[Chi]A + (6762971 - 1994668*\[Nu])*\[Chi]S) + 
           6854080*ei^(74/19)*(31337*\[Delta]*\[Chi]A + 
             (31337 + 235316*\[Nu])*\[Chi]S) - 1172889*e^(26/19)*ei^(48/19)*
            (11122477*\[Delta]*\[Chi]A + (11122477 + 246889772*\[Nu])*
              \[Chi]S) + 1040*e^2*ei^(36/19)*(65005615159*\[Delta]*\[Chi]A + 
             (65005615159 + 28798683196*\[Nu])*\[Chi]S) + 
           65*e^(18/19)*ei^(56/19)*(524587826887*\[Delta]*\[Chi]A + 
             (524587826887 + 2532014127172*\[Nu])*\[Chi]S)))/
         (5246135451648*Sqrt[6]*ei^(48/19)) - 
        (e^(12/19)*Pi*(-9784775*e^(132/19)*ei^(18/19)*
            (81412557579663315503749*\[Delta]*\[Chi]A + 
             (81412557579663315503749 - 18796986832204616606612*\[Nu])*
              \[Chi]S) - 22365200*e^(18/19)*ei^(132/19)*
            (564666014634414440503*\[Delta]*\[Chi]A + 
             (564666014634414440503 - 12179690576134339556*\[Nu])*\[Chi]S) - 
           113343201303652320*e^(112/19)*ei^2*(5132086523719*\[Delta]*
              \[Chi]A + (5132086523719 - 7005190158236*\[Nu])*\[Chi]S) + 
           96519641913542076028800*e^(150/19)*(6762971*\[Delta]*\[Chi]A + 
             (6762971 - 1994668*\[Nu])*\[Chi]S) + 227014100518914508800*
            ei^(150/19)*(31337*\[Delta]*\[Chi]A + (31337 + 235316*\[Nu])*
              \[Chi]S) + 138704005135495800*e^2*ei^(112/19)*
            (65005615159*\[Delta]*\[Chi]A + (65005615159 + 28798683196*\[Nu])*
              \[Chi]S) - 118205598765360*e^(102/19)*ei^(48/19)*
            (994492804490683*\[Delta]*\[Chi]A + (994492804490683 + 
               23425727829215456*\[Nu])*\[Chi]S) + 2568483979776*e^(74/19)*
            ei^4*(56117733091629833*\[Delta]*\[Chi]A + 
             (56117733091629833 + 92077887474718076*\[Nu])*\[Chi]S) - 
           495253887744*e^4*ei^(74/19)*(467358287836566959*\[Delta]*\[Chi]A + 
             (467358287836566959 + 252002971654484168*\[Nu])*\[Chi]S) + 
           1174173000*e^(94/19)*ei^(56/19)*(297124719644282766185*\[Delta]*
              \[Chi]A + (297124719644282766185 + 1319984854731693937916*
                \[Nu])*\[Chi]S) + 6126120*e^(56/19)*ei^(94/19)*
            (11302030209576524265191*\[Delta]*\[Chi]A + 
             (11302030209576524265191 + 3066739024847758928996*\[Nu])*
              \[Chi]S) + 481712*e^6*ei^(36/19)*(1103353378842183123126601*
              \[Delta]*\[Chi]A + (1103353378842183123126601 + 
               638529972502296236631652*\[Nu])*\[Chi]S) - 2737*e^(36/19)*ei^6*
            (5698758298498097478817957*\[Delta]*\[Chi]A + 
             (5698758298498097478817957 + 3149038754181755012665324*\[Nu])*
              \[Chi]S)))/(3667503955293555907840966656*Sqrt[6]*ei^(48/19)) - 
        (e^(12/19)*Pi*(-41860*e^(94/19)*ei^(18/19)*(680665042775389337*
              \[Delta]*\[Chi]A + (680665042775389337 - 167827572913243108*
                \[Nu])*\[Chi]S) - 3227138096*e^(74/19)*ei^2*
            (5132086523719*\[Delta]*\[Chi]A + (5132086523719 - 7005190158236*
                \[Nu])*\[Chi]S) + 3535661359739600*e^(112/19)*
            (6762971*\[Delta]*\[Chi]A + (6762971 - 1994668*\[Nu])*\[Chi]S) - 
           61811053981950*ei^(112/19)*(31337*\[Delta]*\[Chi]A + 
             (31337 + 235316*\[Nu])*\[Chi]S) - 12278974941*e^(64/19)*
            ei^(48/19)*(270023886961*\[Delta]*\[Chi]A + 
             (270023886961 + 6264760409036*\[Nu])*\[Chi]S) + 
           5460*e^(18/19)*ei^(94/19)*(186291062139098407*\[Delta]*\[Chi]A + 
             (186291062139098407 + 48983447376271972*\[Nu])*\[Chi]S) + 
           43056*e^(36/19)*ei^4*(56117733091629833*\[Delta]*\[Chi]A + 
             (56117733091629833 + 92077887474718076*\[Nu])*\[Chi]S) + 
           33488*e^4*ei^(36/19)*(467358287836566959*\[Delta]*\[Chi]A + 
             (467358287836566959 + 252002971654484168*\[Nu])*\[Chi]S) + 
           23345*e^2*(-2741632*ei^(74/19)*(65005615159*\[Delta]*\[Chi]A + 
               65005615159*\[Chi]S + 28798683196*\[Nu]*\[Chi]S) + 
             13*ei^2*(e*ei)^(18/19)*(31424619221651087*\[Delta]*\[Chi]A + 
               31424619221651087*\[Chi]S + 143440553890184612*\[Nu]*
                \[Chi]S))))/(282957072607096602624*Sqrt[6]*ei^(48/19))) + 
      SO^2*x^4*((-5*(-18*(e*ei)^(24/19)*(-145417 + 239316*\[Nu])*
            (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
             191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
              \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) - 
           e^(48/19)*(-252*\[Delta]*\[Kappa]A*(3243271 + 3227408*\[Nu]) + 
             252*\[Kappa]S*(-3243271 + 3259134*\[Nu] + 17101656*\[Nu]^2) - 
             4816146012*\[Chi]A^2 + 5914074257*\[Delta]^2*\[Chi]A^2 + 
             17498024628*\[Nu]*\[Chi]A^2 + 8619234624*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(1097928245 + 67532948*\[Nu])*\[Chi]A*\[Chi]S + 
             1097928245*\[Chi]S^2 + 1901625316*\[Nu]*\[Chi]S^2 - 
             958689760*\[Nu]^2*\[Chi]S^2) + 3249*ei^(48/19)*
            (8*\[Delta]*\[Kappa]A*(-44803 + 25242*\[Nu]) + 
             8*\[Kappa]S*(-44803 + 114848*\[Nu] + 10416*\[Nu]^2) - 
             1234422*\[Chi]A^2 + 1030673*\[Delta]^2*\[Chi]A^2 + 
             4987820*\[Nu]*\[Chi]A^2 + 166656*\[Nu]^2*\[Chi]A^2 + 
             14*\[Delta]*(-29107 + 79316*\[Nu])*\[Chi]A*\[Chi]S - 
             203749*\[Chi]S^2 + 1060292*\[Nu]*\[Chi]S^2 - 382144*\[Nu]^2*
              \[Chi]S^2) + 112*e^(30/19)*ei^(18/19)*
            (15604387*\[Delta]^2*\[Chi]A^2 + 2*\[Delta]*(15604387 - 7259759*
                \[Nu])*\[Chi]A*\[Chi]S + (15604387 - 14519518*\[Nu] + 2512840*
                \[Nu]^2)*\[Chi]S^2) - 324*e^(36/19)*ei^(12/19)*
            (-2833 + 5516*\[Nu])*(-964*\[Kappa]S*(-1 + 2*\[Nu]) + 
             1009*\[Chi]A^2 - 3856*\[Nu]*\[Chi]A^2 + 1009*\[Chi]S^2 - 
             180*\[Nu]*\[Chi]S^2 + \[Delta]*(964*\[Kappa]A + 2018*\[Chi]A*
                \[Chi]S)) + 8*e^(12/19)*ei^(36/19)*
            (-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 116610*\[Nu]^2) - 
             278950701*\[Chi]A^2 + 102215792*\[Delta]^2*\[Chi]A^2 + 
             1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*\[Chi]A^2 - 
             176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
             68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*
                (-124448 + 126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*
                \[Chi]A*\[Chi]S))))/(2934392832*Sqrt[6]*ei^(48/19)) + 
        (e^(12/19)*(-810*e^(12/19)*ei^(62/19)*(-1243916559 + 
             18409137292*\[Nu])*(178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 
             356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 
             382*\[Delta]*\[Chi]A*\[Chi]S + 191*\[Chi]S^2 - 
             52*\[Nu]*\[Chi]S^2) + 3023930*e^(74/19)*
            (-252*\[Delta]*\[Kappa]A*(3243271 + 3227408*\[Nu]) + 
             252*\[Kappa]S*(-3243271 + 3259134*\[Nu] + 17101656*\[Nu]^2) - 
             4816146012*\[Chi]A^2 + 5914074257*\[Delta]^2*\[Chi]A^2 + 
             17498024628*\[Nu]*\[Chi]A^2 + 8619234624*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(1097928245 + 67532948*\[Nu])*\[Chi]A*\[Chi]S + 
             1097928245*\[Chi]S^2 + 1901625316*\[Nu]*\[Chi]S^2 - 
             958689760*\[Nu]^2*\[Chi]S^2) + e^(36/19)*ei^2*
            (36*\[Delta]*\[Kappa]A*(2199940776583 + 45382625640716*\[Nu]) + 
             36*\[Kappa]S*(2199940776583 + 40982744087550*\[Nu] + 
               124668786014136*\[Nu]^2) - 40704717297445158*\[Chi]A^2 + 
             42217477632040165*\[Delta]^2*\[Chi]A^2 + 160552218323235072*
              \[Nu]*\[Chi]A^2 + 8976152593017792*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(1512760334595007 + 436752419096266*\[Nu])*\[Chi]A*
              \[Chi]S + 1512760334595007*\[Chi]S^2 + 3140155704738092*\[Nu]*
              \[Chi]S^2 + 2296792545445216*\[Nu]^2*\[Chi]S^2) + 
           6370*e^(18/19)*ei^(56/19)*(183956440147*\[Delta]^2*\[Chi]A^2 + 
             2*\[Delta]*(183956440147 + 474867803521*\[Nu])*\[Chi]A*\[Chi]S + 
             (183956440147 + 949735607042*\[Nu] - 755722619960*\[Nu]^2)*
              \[Chi]S^2) - 637*e^(56/19)*ei^(18/19)*
            (12285738307079*\[Delta]^2*\[Chi]A^2 + 2*\[Delta]*
              (12285738307079 - 3827691681826*\[Nu])*\[Chi]A*\[Chi]S + 
             (12285738307079 - 7655383363652*\[Nu] + 1110501117344*\[Nu]^2)*
              \[Chi]S^2) - 10530*e^(24/19)*ei^(50/19)*(-2833 + 5516*\[Nu])*
            (10216736*\[Kappa]S*(-1 + 2*\[Nu]) + 50378283*\[Chi]A^2 + 
             40866944*\[Nu]*\[Chi]A^2 + 50378283*\[Chi]S^2 - 
             242380076*\[Nu]*\[Chi]S^2 - 2*\[Delta]*(5108368*\[Kappa]A - 
               50378283*\[Chi]A*\[Chi]S)) + 2835*e^(62/19)*ei^(12/19)*
            (-1571689321 + 2383893876*\[Nu])*(-964*\[Kappa]S*(-1 + 2*\[Nu]) + 
             1009*\[Chi]A^2 - 3856*\[Nu]*\[Chi]A^2 + 1009*\[Chi]S^2 - 
             180*\[Nu]*\[Chi]S^2 + \[Delta]*(964*\[Kappa]A + 2018*\[Chi]A*
                \[Chi]S)) + 585*e^(50/19)*ei^(24/19)*(-145417 + 239316*\[Nu])*
            (-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 25904401*\[Chi]A^2 - 
             99734624*\[Nu]*\[Chi]A^2 + 25904401*\[Chi]S^2 - 
             3882980*\[Nu]*\[Chi]S^2 + \[Delta]*(24933656*\[Kappa]A + 
               51808802*\[Chi]A*\[Chi]S)) + 5997320*ei^(74/19)*
            (-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 116610*\[Nu]^2) - 
             278950701*\[Chi]A^2 + 102215792*\[Delta]^2*\[Chi]A^2 + 
             1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*\[Chi]A^2 - 
             176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
             68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*
                (-124448 + 126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*
                \[Chi]A*\[Chi]S)) - 2345778*e^(26/19)*ei^(48/19)*
            (4*\[Kappa]S*(-438867911 + 1043712010*\[Nu] + 802045608*
                \[Nu]^2) - 14355607086*\[Chi]A^2 + 13638069109*\[Delta]^2*
              \[Chi]A^2 + 57157773104*\[Nu]*\[Chi]A^2 + 6416364864*\[Nu]^2*
              \[Chi]A^2 - 717537977*\[Chi]S^2 + 2293324302*\[Nu]*\[Chi]S^2 + 
             3402549976*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(\[Kappa]A*
                (-877735822 + 331952376*\[Nu]) + (-717537977 + 1014334531*
                  \[Nu])*\[Chi]A*\[Chi]S)) - 10*e^2*ei^(36/19)*
            (-18*\[Kappa]S*(30936092472295 - 97020967578282*\[Nu] + 
               36669175275096*\[Nu]^2) - 2827049041226343*\[Chi]A^2 + 
             2206814975744288*\[Delta]^2*\[Chi]A^2 + 11721265129529220*\[Nu]*
              \[Chi]A^2 - 1320090309903456*\[Nu]^2*\[Chi]A^2 - 
             620234065482055*\[Chi]S^2 + 1521486691850320*\[Nu]*\[Chi]S^2 - 
             141962774214352*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(9*\[Kappa]A*
                (-30936092472295 + 35148782633692*\[Nu]) + 
               (-620234065482055 + 967277828237084*\[Nu])*\[Chi]A*\[Chi]S))))/
         (385590955696128*Sqrt[6]*ei^(48/19)) + 
        (e^(12/19)*(2916369456*e^(12/19)*ei^(138/19)*(-2719168921326953119 + 
             2093118343804025580*\[Nu])*(178*\[Delta]*\[Kappa]A + 
             178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 
             712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*\[Chi]S + 
             191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) + 161*e^(36/19)*ei^6*
            (8352*\[Delta]*\[Kappa]A*(-16199191005975146527063767 + 
               10672817768043867807307300*\[Nu]) - 8352*\[Kappa]S*
              (16199191005975146527063767 - 43071199779994160861434834*
                \[Nu] + 13127449498976183773791480*\[Nu]^2) - 
             271703657349649599005194133916*\[Chi]A^2 + 
             79027448227137743995475355269*\[Delta]^2*\[Chi]A^2 + 
             1149991954855717089462215205552*\[Nu]*\[Chi]A^2 - 
             219280916430898173757412881920*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-192676209122511855009718778647 + 
               122289055338280436394188580362*\[Nu])*\[Chi]A*\[Chi]S - 
             192676209122511855009718778647*\[Chi]S^2 + 
             181400785219442179346938490836*\[Nu]*\[Chi]S^2 - 
             14149960848178490011498153696*\[Nu]^2*\[Chi]S^2) - 
           28336*e^6*ei^(36/19)*(33408*\[Delta]*\[Kappa]A*
              (-138222143243089778236681 + 169677614533377548972647*\[Nu]) - 
             33408*\[Kappa]S*(138222143243089778236681 - 
               446121901019557105446009*\[Nu] + 187902755346189154082994*
                \[Nu]^2) - 30921703483728105055475662749*\[Chi]A^2 + 
             27011676587914586707090685734*\[Delta]^2*\[Chi]A^2 + 
             128203625765713966307445860916*\[Nu]*\[Chi]A^2 - 
             12554910501210974519209327104*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-3910026895813518348384977015 + 
               7969183702013203055213796514*\[Nu])*\[Chi]A*\[Chi]S - 
             3910026895813518348384977015*\[Chi]S^2 + 
             11421555573224860024884383108*\[Nu]*\[Chi]S^2 - 
             230627307042699097000726256*\[Nu]^2*\[Chi]S^2) + 
           20635578656*e^4*ei^(74/19)*(522*\[Delta]*\[Kappa]A*
              (-5919500807162614699 + 7069747731165174748*\[Nu]) - 
             522*\[Kappa]S*(5919500807162614699 - 18908749345490404146*
                \[Nu] + 7669345784813238456*\[Nu]^2) - 
             18928431938173283618721*\[Chi]A^2 + 16007908125685966389172*
              \[Delta]^2*\[Chi]A^2 + 78468871437192708323052*\[Nu]*
              \[Chi]A^2 - 8006796999345020948064*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-2920523812487317229549 + 5349738482100710957140*
                \[Nu])*\[Chi]A*\[Chi]S - 2920523812487317229549*\[Chi]S^2 + 
             7944333279701848066112*\[Nu]*\[Chi]S^2 - 384534241038124715888*
              \[Nu]^2*\[Chi]S^2) + 52281472703168624515600*e^(150/19)*
            (-252*\[Delta]*\[Kappa]A*(3243271 + 3227408*\[Nu]) + 
             252*\[Kappa]S*(-3243271 + 3259134*\[Nu] + 17101656*\[Nu]^2) - 
             4816146012*\[Chi]A^2 + 5914074257*\[Delta]^2*\[Chi]A^2 + 
             17498024628*\[Nu]*\[Chi]A^2 + 8619234624*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(1097928245 + 67532948*\[Nu])*\[Chi]A*\[Chi]S + 
             1097928245*\[Chi]S^2 + 1901625316*\[Nu]*\[Chi]S^2 - 
             958689760*\[Nu]^2*\[Chi]S^2) + 37781067101217440*e^(112/19)*ei^2*
            (36*\[Delta]*\[Kappa]A*(2199940776583 + 45382625640716*\[Nu]) + 
             36*\[Kappa]S*(2199940776583 + 40982744087550*\[Nu] + 
               124668786014136*\[Nu]^2) - 40704717297445158*\[Chi]A^2 + 
             42217477632040165*\[Delta]^2*\[Chi]A^2 + 160552218323235072*
              \[Nu]*\[Chi]A^2 + 8976152593017792*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(1512760334595007 + 436752419096266*\[Nu])*\[Chi]A*
              \[Chi]S + 1512760334595007*\[Chi]S^2 + 3140155704738092*\[Nu]*
              \[Chi]S^2 + 2296792545445216*\[Nu]^2*\[Chi]S^2) - 
           8918347152*e^(74/19)*ei^4*(1044*\[Delta]*\[Kappa]A*
              (-7365045929568369393 + 5742829720995328276*\[Nu]) - 
             1044*\[Kappa]S*(7365045929568369393 - 20472921580132067062*
                \[Nu] + 1517701378547771592*\[Nu]^2) - 
             36820157636844044300469*\[Chi]A^2 + 30945138675512339580484*
              \[Delta]^2*\[Chi]A^2 + 153005155284657572553108*\[Nu]*
              \[Chi]A^2 - 3168960478407747084096*\[Nu]^2*\[Chi]A^2 + 
             86*\[Delta]*(-136628347937946621395 + 186569357031292881772*
                \[Nu])*\[Chi]A*\[Chi]S - 5875018961331704719985*\[Chi]S^2 + 
             10320439967409792481160*\[Nu]*\[Chi]S^2 + 5836439695605526530544*
              \[Nu]^2*\[Chi]S^2) + 1994528536*e^(94/19)*ei^(56/19)*
            (38094625234921853405987*\[Delta]^2*\[Chi]A^2 + 
             2*\[Delta]*(38094625234921853405987 + 106199811219640248726242*
                \[Nu])*\[Chi]A*\[Chi]S + (38094625234921853405987 + 
               212399622439280497452484*\[Nu] - 64304086929034066692928*
                \[Nu]^2)*\[Chi]S^2) + 260155896*e^(56/19)*ei^(94/19)*
            (49342927045031166777799*\[Delta]^2*\[Chi]A^2 + 
             2*\[Delta]*(49342927045031166777799 + 2439721838610957441694*
                \[Nu])*\[Chi]A*\[Chi]S + (49342927045031166777799 + 
               4879443677221914883388*\[Nu] - 9550452106462359242336*\[Nu]^2)*
              \[Chi]S^2) - 28493264800*e^(18/19)*ei^(132/19)*
            (58396459803910062971*\[Delta]^2*\[Chi]A^2 + 2*\[Delta]*
              (58396459803910062971 - 14300636516592891181*\[Nu])*\[Chi]A*
              \[Chi]S + (58396459803910062971 - 28601273033185782362*\[Nu] - 
               8627253776793918640*\[Nu]^2)*\[Chi]S^2) - 
           1246580335*e^(132/19)*ei^(18/19)*(136020958334009543602051*
              \[Delta]^2*\[Chi]A^2 + 2*\[Delta]*(136020958334009543602051 - 
               31496851652822095944842*\[Nu])*\[Chi]A*\[Chi]S + 
             (136020958334009543602051 - 62993703305644191889684*\[Nu] + 
               7292971603193993485600*\[Nu]^2)*\[Chi]S^2) - 
           18387549180*e^(100/19)*ei^(50/19)*(-117334502622439 + 
             164817749118204*\[Nu])*(10216736*\[Kappa]S*(-1 + 2*\[Nu]) + 
             50378283*\[Chi]A^2 + 40866944*\[Nu]*\[Chi]A^2 + 
             50378283*\[Chi]S^2 - 242380076*\[Nu]*\[Chi]S^2 - 
             2*\[Delta]*(5108368*\[Kappa]A - 50378283*\[Chi]A*\[Chi]S)) + 
           60570750240*e^(138/19)*ei^(12/19)*(-1663373241307010153 + 
             2252417164528514620*\[Nu])*(-964*\[Kappa]S*(-1 + 2*\[Nu]) + 
             1009*\[Chi]A^2 - 3856*\[Nu]*\[Chi]A^2 + 1009*\[Chi]S^2 - 
             180*\[Nu]*\[Chi]S^2 + \[Delta]*(964*\[Kappa]A + 2018*\[Chi]A*
                \[Chi]S)) - 30174439680*e^(50/19)*ei^(100/19)*
            (-23500911051568 + 15843764110629*\[Nu])*
            (-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 25904401*\[Chi]A^2 - 
             99734624*\[Nu]*\[Chi]A^2 + 25904401*\[Chi]S^2 - 
             3882980*\[Nu]*\[Chi]S^2 + \[Delta]*(24933656*\[Kappa]A + 
               51808802*\[Chi]A*\[Chi]S)) - 1697312232*e^(62/19)*ei^(88/19)*
            (-1571689321 + 2383893876*\[Nu])*(-4853042174758*\[Kappa]S*
              (-1 + 2*\[Nu]) + 4806182243753*\[Chi]A^2 - 19412168699032*\[Nu]*
              \[Chi]A^2 + 4806182243753*\[Chi]S^2 + 187439724020*\[Nu]*
              \[Chi]S^2 + \[Delta]*(4853042174758*\[Kappa]A + 9612364487506*
                \[Chi]A*\[Chi]S)) + 1508721984*e^(24/19)*ei^(126/19)*
            (-2833 + 5516*\[Nu])*(-331165474334851236*\[Kappa]S*
              (-1 + 2*\[Nu]) + 340045705042503011*\[Chi]A^2 - 
             1324661897339404944*\[Nu]*\[Chi]A^2 + 340045705042503011*
              \[Chi]S^2 - 35520922830607100*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(331165474334851236*\[Kappa]A + 680091410085006022*
                \[Chi]A*\[Chi]S)) + 272408136*e^(126/19)*ei^(24/19)*
            (-145417 + 239316*\[Nu])*(-1165556009138332944*\[Kappa]S*
              (-1 + 2*\[Nu]) + 1189984978416247669*\[Chi]A^2 - 
             4662224036553331776*\[Nu]*\[Chi]A^2 + 1189984978416247669*
              \[Chi]S^2 - 97715877111658900*\[Nu]*\[Chi]S^2 + 
             2*\[Delta]*(582778004569166472*\[Kappa]A + 1189984978416247669*
                \[Chi]A*\[Chi]S)) + 860761797800884179200*ei^(150/19)*
            (-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 116610*\[Nu]^2) - 
             278950701*\[Chi]A^2 + 102215792*\[Delta]^2*\[Chi]A^2 + 
             1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*\[Chi]A^2 - 
             176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
             68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*
                (-124448 + 126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*
                \[Chi]A*\[Chi]S)) - 5779333547312325*e^2*ei^(112/19)*
            (-18*\[Kappa]S*(30936092472295 - 97020967578282*\[Nu] + 
               36669175275096*\[Nu]^2) - 2827049041226343*\[Chi]A^2 + 
             2206814975744288*\[Delta]^2*\[Chi]A^2 + 11721265129529220*\[Nu]*
              \[Chi]A^2 - 1320090309903456*\[Nu]^2*\[Chi]A^2 - 
             620234065482055*\[Chi]S^2 + 1521486691850320*\[Nu]*\[Chi]S^2 - 
             141962774214352*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(9*\[Kappa]A*
                (-30936092472295 + 35148782633692*\[Nu]) + 
               (-620234065482055 + 967277828237084*\[Nu])*\[Chi]A*\[Chi]S)) - 
           1477569984567*e^(102/19)*ei^(48/19)*
            (46*\[Kappa]S*(-712537575477233533 + 355569756646147118*\[Nu] + 
               6177242842172456376*\[Nu]^2) - 678842802091667076779*
              \[Chi]A^2 + 751787167076801473800*\[Delta]^2*\[Chi]A^2 + 
             2760035566056729353524*\[Nu]*\[Chi]A^2 + 568306341479865986592*
              \[Nu]^2*\[Chi]A^2 + 72944364985134397021*\[Chi]S^2 - 
             107042776711812379664*\[Nu]*\[Chi]S^2 + 382570925764065748080*
              \[Nu]^2*\[Chi]S^2 - 2*\[Delta]*(23*\[Kappa]A*
                (712537575477233533 + 1069505394308319948*\[Nu]) + 
               (-72944364985134397021 + 31189209510875666628*\[Nu])*\[Chi]A*
                \[Chi]S)) - 3232975680*e^(88/19)*ei^(62/19)*
            (-1243916559 + 18409137292*\[Nu])*(-2948986828885*\[Kappa]S*
              (-1 + 2*\[Nu]) + \[Delta]*(2948986828885*\[Kappa]A + 
               6058122836776*\[Chi]A*\[Chi]S) - 4*((-757265354597 + 
                 2948986828885*\[Nu])*\[Chi]A^2 + (-757265354597 + 
                 80074589503*\[Nu])*\[Chi]S^2))))/
         (1168100009760997556647347879936*Sqrt[6]*ei^(48/19)) + 
        (e^(12/19)*(-564762240*e^(12/19)*ei^(100/19)*(-23500911051568 + 
             15843764110629*\[Nu])*(178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 
             356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 
             382*\[Delta]*\[Chi]A*\[Chi]S + 191*\[Chi]S^2 - 
             52*\[Nu]*\[Chi]S^2) - 16744*e^4*ei^(36/19)*
            (522*\[Delta]*\[Kappa]A*(-5919500807162614699 + 
               7069747731165174748*\[Nu]) - 522*\[Kappa]S*
              (5919500807162614699 - 18908749345490404146*\[Nu] + 
               7669345784813238456*\[Nu]^2) - 18928431938173283618721*
              \[Chi]A^2 + 16007908125685966389172*\[Delta]^2*\[Chi]A^2 + 
             78468871437192708323052*\[Nu]*\[Chi]A^2 - 8006796999345020948064*
              \[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(-2920523812487317229549 + 
               5349738482100710957140*\[Nu])*\[Chi]A*\[Chi]S - 
             2920523812487317229549*\[Chi]S^2 + 7944333279701848066112*\[Nu]*
              \[Chi]S^2 - 384534241038124715888*\[Nu]^2*\[Chi]S^2) + 
           22981798838307400*e^(112/19)*(-252*\[Delta]*\[Kappa]A*
              (3243271 + 3227408*\[Nu]) + 252*\[Kappa]S*(-3243271 + 3259134*
                \[Nu] + 17101656*\[Nu]^2) - 4816146012*\[Chi]A^2 + 
             5914074257*\[Delta]^2*\[Chi]A^2 + 17498024628*\[Nu]*\[Chi]A^2 + 
             8619234624*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(1097928245 + 67532948*
                \[Nu])*\[Chi]A*\[Chi]S + 1097928245*\[Chi]S^2 + 
             1901625316*\[Nu]*\[Chi]S^2 - 958689760*\[Nu]^2*\[Chi]S^2) + 
           12908552384*e^(74/19)*ei^2*(36*\[Delta]*\[Kappa]A*
              (2199940776583 + 45382625640716*\[Nu]) + 36*\[Kappa]S*
              (2199940776583 + 40982744087550*\[Nu] + 124668786014136*
                \[Nu]^2) - 40704717297445158*\[Chi]A^2 + 42217477632040165*
              \[Delta]^2*\[Chi]A^2 + 160552218323235072*\[Nu]*\[Chi]A^2 + 
             8976152593017792*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*
              (1512760334595007 + 436752419096266*\[Nu])*\[Chi]A*\[Chi]S + 
             1512760334595007*\[Chi]S^2 + 3140155704738092*\[Nu]*\[Chi]S^2 + 
             2296792545445216*\[Nu]^2*\[Chi]S^2) - 1794*e^(36/19)*ei^4*
            (1044*\[Delta]*\[Kappa]A*(-7365045929568369393 + 
               5742829720995328276*\[Nu]) - 1044*\[Kappa]S*
              (7365045929568369393 - 20472921580132067062*\[Nu] + 
               1517701378547771592*\[Nu]^2) - 36820157636844044300469*
              \[Chi]A^2 + 30945138675512339580484*\[Delta]^2*\[Chi]A^2 + 
             153005155284657572553108*\[Nu]*\[Chi]A^2 - 
             3168960478407747084096*\[Nu]^2*\[Chi]A^2 + 86*\[Delta]*
              (-136628347937946621395 + 186569357031292881772*\[Nu])*\[Chi]A*
              \[Chi]S - 5875018961331704719985*\[Chi]S^2 + 
             10320439967409792481160*\[Nu]*\[Chi]S^2 + 5836439695605526530544*
              \[Nu]^2*\[Chi]S^2) + 27824160*e^(18/19)*ei^(94/19)*
            (62671539152008547*\[Delta]^2*\[Chi]A^2 - 2*\[Delta]*
              (-62671539152008547 + 6532783285303279*\[Nu])*\[Chi]A*\[Chi]S + 
             (62671539152008547 - 13065566570606558*\[Nu] - 21610746442651960*
                \[Nu]^2)*\[Chi]S^2) - 21331856*e^(94/19)*ei^(18/19)*
            (3231434976186023027*\[Delta]^2*\[Chi]A^2 + 2*\[Delta]*
              (3231434976186023027 - 836498434314498718*\[Nu])*\[Chi]A*
              \[Chi]S + (3231434976186023027 - 1672996868628997436*\[Nu] + 
               213816389149906112*\[Nu]^2)*\[Chi]S^2) - 393316560*e^(62/19)*
            ei^(50/19)*(-1571689321 + 2383893876*\[Nu])*
            (10216736*\[Kappa]S*(-1 + 2*\[Nu]) + 50378283*\[Chi]A^2 + 
             40866944*\[Nu]*\[Chi]A^2 + 50378283*\[Chi]S^2 - 
             242380076*\[Nu]*\[Chi]S^2 - 2*\[Delta]*(5108368*\[Kappa]A - 
               50378283*\[Chi]A*\[Chi]S)) + 344151990*e^(100/19)*ei^(12/19)*
            (-117334502622439 + 164817749118204*\[Nu])*
            (-964*\[Kappa]S*(-1 + 2*\[Nu]) + 1009*\[Chi]A^2 - 
             3856*\[Nu]*\[Chi]A^2 + 1009*\[Chi]S^2 - 180*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(964*\[Kappa]A + 2018*\[Chi]A*\[Chi]S)) - 
           112376160*e^(50/19)*ei^(62/19)*(-1243916559 + 18409137292*\[Nu])*
            (-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 25904401*\[Chi]A^2 - 
             99734624*\[Nu]*\[Chi]A^2 + 25904401*\[Chi]S^2 - 
             3882980*\[Nu]*\[Chi]S^2 + \[Delta]*(24933656*\[Kappa]A + 
               51808802*\[Chi]A*\[Chi]S)) - 117994968*e^(24/19)*ei^(88/19)*
            (-2833 + 5516*\[Nu])*(-4853042174758*\[Kappa]S*(-1 + 2*\[Nu]) + 
             4806182243753*\[Chi]A^2 - 19412168699032*\[Nu]*\[Chi]A^2 + 
             4806182243753*\[Chi]S^2 + 187439724020*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(4853042174758*\[Kappa]A + 9612364487506*\[Chi]A*
                \[Chi]S)) - 2812402956178725*ei^(112/19)*
            (-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 116610*\[Nu]^2) - 
             278950701*\[Chi]A^2 + 102215792*\[Delta]^2*\[Chi]A^2 + 
             1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*\[Chi]A^2 - 
             176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
             68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*
                (-124448 + 126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*
                \[Chi]A*\[Chi]S)) - 12278974941*e^(64/19)*ei^(48/19)*
            (368*\[Kappa]S*(-5800115804177 + 9056885427979*\[Nu] + 
               28442618518314*\[Nu]^2) - 29636377823388555*\[Chi]A^2 + 
             31181399186580188*\[Delta]^2*\[Chi]A^2 + 119503932749870596*
              \[Nu]*\[Chi]A^2 + 20933767229479104*\[Nu]^2*\[Chi]A^2 + 
             1545021363191633*\[Chi]S^2 - 1650033761251200*\[Nu]*\[Chi]S^2 + 
             13396573285480112*\[Nu]^2*\[Chi]S^2 - 2*\[Delta]*
              (184*\[Kappa]A*(5800115804177 + 2543346180375*\[Nu]) + 
               (-1545021363191633 + 345806152467412*\[Nu])*\[Chi]A*
                \[Chi]S)) + 43701840*e^(88/19)*ei^(24/19)*
            (-145417 + 239316*\[Nu])*(-2948986828885*\[Kappa]S*
              (-1 + 2*\[Nu]) + \[Delta]*(2948986828885*\[Kappa]A + 
               6058122836776*\[Chi]A*\[Chi]S) - 4*((-757265354597 + 
                 2948986828885*\[Nu])*\[Chi]A^2 + (-757265354597 + 
                 80074589503*\[Nu])*\[Chi]S^2)) + 18676*e^2*
            (8281*ei^2*(e*ei)^(18/19)*(144833672963114999*\[Delta]^2*
                \[Chi]A^2 + 2*\[Delta]*(144833672963114999 + 
                 396134150578580894*\[Nu])*\[Chi]A*\[Chi]S + 
               (144833672963114999 + 792268301157161788*\[Nu] - 
                 333977019574551136*\[Nu]^2)*\[Chi]S^2) + 1713520*ei^(74/19)*
              (-18*\[Kappa]S*(30936092472295 - 97020967578282*\[Nu] + 
                 36669175275096*\[Nu]^2) - 2827049041226343*\[Chi]A^2 + 
               2206814975744288*\[Delta]^2*\[Chi]A^2 + 11721265129529220*
                \[Nu]*\[Chi]A^2 - 1320090309903456*\[Nu]^2*\[Chi]A^2 - 
               620234065482055*\[Chi]S^2 + 1521486691850320*\[Nu]*\[Chi]S^2 - 
               141962774214352*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*
                (9*\[Kappa]A*(-30936092472295 + 35148782633692*\[Nu]) + 
                 (-620234065482055 + 967277828237084*\[Nu])*\[Chi]A*
                  \[Chi]S)))))/(1081461931504323215228928*Sqrt[6]*
          ei^(48/19))))
 
HlmDCMem[4, 0] = ((e^(12/19) - ei^(12/19))/(504*Sqrt[2]*ei^(12/19)) + 
      (86398*e^2 + 1805*e^(26/19)*ei^(12/19) - 88203*ei^2)/
       (37844352*Sqrt[2]*(ei/e)^(12/19)) + 
      (15710523296*e^4 + 343395835*e^(64/19)*ei^(12/19) - 
        18758308416*e^2*ei^2 + 2704389285*ei^4)/(6994847268864*Sqrt[2]*
        (ei/e)^(12/19)) + (258018423716352*e^6 + 5588815433395*e^(102/19)*
         ei^(12/19) - 302014870086032*e^4*ei^2 + 50924551699645*e^2*ei^4 - 
        12516920763360*ei^6)/(114473007170715648*Sqrt[2]*(ei/e)^(12/19)))*x + 
    x^2*\[Epsilon]^2*(((e^(12/19) - ei^(12/19))*
        (361*ei^(12/19)*(-9479 + 40124*\[Nu]) + e^(12/19)*
          (-3920527 + 15455580*\[Nu])))/(283143168*Sqrt[2]*ei^(24/19)) + 
      (e^(12/19)*(3*e^(12/19)*ei^2*(1859009125901 - 8124219720708*\[Nu]) + 
         27166524*ei^(50/19)*(-2833 + 5516*\[Nu]) + 302393*e^(50/19)*
          (-3920527 + 15455580*\[Nu]) + 54872*e^(26/19)*ei^(24/19)*
          (-84703117 + 364299110*\[Nu]) - 308*e^2*ei^(12/19)*
          (-1082158253 + 1430751140*\[Nu])))/(37206144847872*Sqrt[2]*
        ei^(24/19)) + (e^(12/19)*(-2707093674285*ei^(88/19)*
          (-2833 + 5516*\[Nu]) + 283211804148*e^(88/19)*
          (-3920527 + 15455580*\[Nu]) + 217332192*e^2*ei^(50/19)*
          (-1082158253 + 1430751140*\[Nu]) - 4147104*e^(50/19)*ei^2*
          (-1859009125901 + 8124219720708*\[Nu]) + 2476099*e^(64/19)*
          ei^(24/19)*(-2906651517879 + 12249998813780*\[Nu]) - 
         8008*e^4*ei^(12/19)*(-42164487546687 + 54413587634588*\[Nu]) - 
         84*e^(12/19)*ei^4*(-5807829430373107 + 10448813703644316*\[Nu])))/
       (22349880034696101888*Sqrt[2]*ei^(24/19)) + 
      (e^(12/19)*(75176626104740160*ei^(126/19)*(-2833 + 5516*\[Nu]) + 
         37800454570592168*e^(126/19)*(-3920527 + 15455580*\[Nu]) - 
         3540045574065*e^2*ei^(88/19)*(-1082158253 + 1430751140*\[Nu]) - 
         634892396112*e^(88/19)*ei^2*(-1859009125901 + 8124219720708*\[Nu]) + 
         923661816*e^4*ei^(50/19)*(-42164487546687 + 54413587634588*\[Nu]) + 
         893871739*e^(102/19)*ei^(24/19)*(-1263862590860999 + 
           5300436799281260*\[Nu]) - 18980976*e^(50/19)*ei^4*
          (-5807829430373107 + 10448813703644316*\[Nu]) - 
         32032*e^6*ei^(12/19)*(-1481083213734247525 + 1878694147208499348*
            \[Nu]) + 14560*e^(12/19)*ei^6*(-1696052141317888257 + 
           3349191345990271600*\[Nu])))/(2194579420366879636586496*Sqrt[2]*
        ei^(24/19))) + SO*x^(7/2)*\[Epsilon]^5*
     ((70*e^(24/19)*ei^(18/19)*(-3920527 + 15455580*\[Nu])*
         (-157*\[Delta]*\[Chi]A + (-157 + 110*\[Nu])*\[Chi]S) - 
        1925*e^(30/19)*ei^(12/19)*(-2833 + 5516*\[Nu])*
         (5189*\[Delta]*\[Chi]A + (5189 + 10748*\[Nu])*\[Chi]S) + 
        1232*e^(12/19)*ei^(30/19)*(4*\[Delta]*(-5883581 + 3034822*\[Nu])*
           \[Chi]A + (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*
           \[Chi]S) - 21660*ei^(42/19)*(\[Delta]*(-1319709 + 12539422*\[Nu])*
           \[Chi]A + (-1319709 + 6548872*\[Nu] + 5249776*\[Nu]^2)*\[Chi]S) + 
        e^(42/19)*(\[Delta]*(-70975542727 + 481603460604*\[Nu])*\[Chi]A + 
          (-70975542727 + 289947696068*\[Nu] + 119614397280*\[Nu]^2)*
           \[Chi]S))/(1694611860480*Sqrt[2]*ei^(42/19)) + 
      (e^(12/19)*(240*e^(12/19)*ei^(56/19)*(-1859009125901 + 
           8124219720708*\[Nu])*(157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*
            \[Chi]S) - 7700*e^(56/19)*ei^(12/19)*(-724653277 + 
           1072804096*\[Nu])*(5189*\[Delta]*\[Chi]A + (5189 + 10748*\[Nu])*
            \[Chi]S) + 91*e^(50/19)*ei^(18/19)*(-3920527 + 15455580*\[Nu])*
          (-113175949*\[Delta]*\[Chi]A + (-113175949 + 41301776*\[Nu])*
            \[Chi]S) + 60060*e^(18/19)*ei^(50/19)*(-2833 + 5516*\[Nu])*
          (31277492*\[Delta]*\[Chi]A + (31277492 + 1656666347*\[Nu])*
            \[Chi]S) + 154*e^2*ei^(30/19)*(\[Delta]*(-118659872446897 + 
             113510738272204*\[Nu])*\[Chi]A + (-118659872446897 + 
             286545005193108*\[Nu] - 61303880674400*\[Nu]^2)*\[Chi]S) - 
         217332192*ei^(68/19)*(4*\[Delta]*(-5883581 + 3034822*\[Nu])*
            \[Chi]A + (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*
            \[Chi]S) + 604786*e^(68/19)*(\[Delta]*(-70975542727 + 
             481603460604*\[Nu])*\[Chi]A + (-70975542727 + 289947696068*
              \[Nu] + 119614397280*\[Nu]^2)*\[Chi]S) - 3909630*e^(26/19)*
          ei^(42/19)*(\[Delta]*(-193129614957 + 884416890364*\[Nu])*\[Chi]A + 
           (-193129614957 + 429487218226*\[Nu] + 833412262648*\[Nu]^2)*
            \[Chi]S) + 21*e^(30/19)*ei^2*(7*\[Delta]*(-4713900887677313 + 
             20661171666495828*\[Nu])*\[Chi]A + (-32997306213741191 + 
             73396240249492084*\[Nu] + 137612708364899040*\[Nu]^2)*\[Chi]S)))/
       (254490030759444480*Sqrt[2]*ei^(42/19)) + 
      (e^(12/19)*(3537768*e^(50/19)*ei^(56/19)*(-1859009125901 + 
           8124219720708*\[Nu])*(113175949*\[Delta]*\[Chi]A + 
           (113175949 - 41301776*\[Nu])*\[Chi]S) + 19049520*e^(12/19)*
          ei^(94/19)*(-5807829430373107 + 10448813703644316*\[Nu])*
          (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) - 
         283758475*e^(94/19)*ei^(12/19)*(-96101931245055 + 
           132465768779708*\[Nu])*(5189*\[Delta]*\[Chi]A + 
           (5189 + 10748*\[Nu])*\[Chi]S) + 2724081360*e^(56/19)*ei^(50/19)*
          (-724653277 + 1072804096*\[Nu])*(31277492*\[Delta]*\[Chi]A + 
           (31277492 + 1656666347*\[Nu])*\[Chi]S) + 258104574*e^(88/19)*
          ei^(18/19)*(-3920527 + 15455580*\[Nu])*(-189611322769*\[Delta]*
            \[Chi]A + (-189611322769 + 49774532976*\[Nu])*\[Chi]S) + 
         13273260*e^(18/19)*ei^(88/19)*(-2833 + 5516*\[Nu])*
          (70548776005597*\[Delta]*\[Chi]A + (70548776005597 + 
             205639257894652*\[Nu])*\[Chi]S) + 68068*e^4*ei^(30/19)*
          (\[Delta]*(-881417005050492936489 + 1269739664742855690748*\[Nu])*
            \[Chi]A + (-881417005050492936489 + 2871898607753603335636*
              \[Nu] - 610400986757547833120*\[Nu]^2)*\[Chi]S) - 
         308041215636*e^2*ei^(68/19)*(\[Delta]*(-118659872446897 + 
             113510738272204*\[Nu])*\[Chi]A + (-118659872446897 + 
             286545005193108*\[Nu] - 61303880674400*\[Nu]^2)*\[Chi]S) + 
         61391470345435230*ei^(106/19)*(4*\[Delta]*(-5883581 + 3034822*\[Nu])*
            \[Chi]A + (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*
            \[Chi]S) + 2494402582250124*e^(106/19)*
          (\[Delta]*(-70975542727 + 481603460604*\[Nu])*\[Chi]A + 
           (-70975542727 + 289947696068*\[Nu] + 119614397280*\[Nu]^2)*
            \[Chi]S) + 144011037534*e^(68/19)*ei^2*
          (7*\[Delta]*(-4713900887677313 + 20661171666495828*\[Nu])*\[Chi]A + 
           (-32997306213741191 + 73396240249492084*\[Nu] + 137612708364899040*
              \[Nu]^2)*\[Chi]S) - 40929916470*e^(64/19)*ei^(42/19)*
          (8*\[Delta]*(-16557657366135909 + 73963894327673623*\[Nu])*
            \[Chi]A + (-132461258929087272 + 283298459287770341*\[Nu] + 
             551036274021854348*\[Nu]^2)*\[Chi]S) - 273*e^(30/19)*ei^4*
          (\[Delta]*(-91579519882877550261559 + 1070089805632513857206652*
              \[Nu])*\[Chi]A + (-91579519882877550261559 - 
             1886006275288386965668792*\[Nu] + 3593907680873728598492112*
              \[Nu]^2)*\[Chi]S)))/(433357245409946659816734720*Sqrt[2]*
        ei^(42/19)) + (e^(12/19)*(272286144*e^(88/19)*ei^(56/19)*
          (-1859009125901 + 8124219720708*\[Nu])*
          (189611322769*\[Delta]*\[Chi]A + (189611322769 - 49774532976*\[Nu])*
            \[Chi]S) + 7619808*e^(50/19)*ei^(94/19)*(-5807829430373107 + 
           10448813703644316*\[Nu])*(113175949*\[Delta]*\[Chi]A + 
           (113175949 - 41301776*\[Nu])*\[Chi]S) + 1553843200*e^(12/19)*
          ei^(132/19)*(-1696052141317888257 + 3349191345990271600*\[Nu])*
          (-157*\[Delta]*\[Chi]A + (-157 + 110*\[Nu])*\[Chi]S) - 
         19095276200*e^(132/19)*ei^(12/19)*(-137799513225092197 + 
           183838168203252052*\[Nu])*(5189*\[Delta]*\[Chi]A + 
           (5189 + 10748*\[Nu])*\[Chi]S) + 2724081360*e^(94/19)*ei^(50/19)*
          (-96101931245055 + 132465768779708*\[Nu])*
          (31277492*\[Delta]*\[Chi]A + (31277492 + 1656666347*\[Nu])*
            \[Chi]S) + 16336320*e^(56/19)*ei^(88/19)*(-724653277 + 
           1072804096*\[Nu])*(70548776005597*\[Delta]*\[Chi]A + 
           (70548776005597 + 205639257894652*\[Nu])*\[Chi]S) - 
         436121400*e^(18/19)*ei^(126/19)*(-2833 + 5516*\[Nu])*
          (110801570912666765*\[Delta]*\[Chi]A + (110801570912666765 + 
             267721377018031016*\[Nu])*\[Chi]S) + 925106*e^(126/19)*
          ei^(18/19)*(-3920527 + 15455580*\[Nu])*
          (-5038160345051358953*\[Delta]*\[Chi]A + 
           (-5038160345051358953 + 1065846237840519832*\[Nu])*\[Chi]S) - 
         26*e^(30/19)*ei^6*(7*\[Delta]*(-17430637311335814129825495207 + 
             8869318287655361435797614236*\[Nu])*\[Chi]A + 
           (-122014461179350698908778466449 + 267812540155742803963588981748*
              \[Nu] - 331688234075643406273569519808*\[Nu]^2)*\[Chi]S) + 
         1144*e^6*ei^(30/19)*(7*\[Delta]*(-532202728837787049264046867 + 
             1052272546464970879844248804*\[Nu])*\[Chi]A + 
           (-3725419101864509344848328069 + 15569484378055550886037435956*
              \[Nu] - 3304589732988741434471430560*\[Nu]^2)*\[Chi]S) - 
         3694647264*e^4*ei^(68/19)*(\[Delta]*(-881417005050492936489 + 
             1269739664742855690748*\[Nu])*\[Chi]A + 
           (-881417005050492936489 + 2871898607753603335636*\[Nu] - 
             610400986757547833120*\[Nu]^2)*\[Chi]S) + 2361210397901355*e^2*
          ei^(106/19)*(\[Delta]*(-118659872446897 + 113510738272204*\[Nu])*
            \[Chi]A + (-118659872446897 + 286545005193108*\[Nu] - 
             61303880674400*\[Nu]^2)*\[Chi]S) - 802284953789786987520*
          ei^(144/19)*(4*\[Delta]*(-5883581 + 3034822*\[Nu])*\[Chi]A + 
           (-23534324 + 39309761*\[Nu] - 8755460*\[Nu]^2)*\[Chi]S) + 
         221031591042282852832*e^(144/19)*(\[Delta]*(-70975542727 + 
             481603460604*\[Nu])*\[Chi]A + (-70975542727 + 289947696068*
              \[Nu] + 119614397280*\[Nu]^2)*\[Chi]S) + 16117678223770032*
          e^(106/19)*ei^2*(7*\[Delta]*(-4713900887677313 + 20661171666495828*
              \[Nu])*\[Chi]A + (-32997306213741191 + 73396240249492084*
              \[Nu] + 137612708364899040*\[Nu]^2)*\[Chi]S) - 
         2462616640945*e^(102/19)*ei^(42/19)*
          (\[Delta]*(-251156394462857060739 + 1110679241082372140452*\[Nu])*
            \[Chi]A + (-251156394462857060739 + 531445907574370510588*\[Nu] + 
             1029096221279830054624*\[Nu]^2)*\[Chi]S) - 50802024*e^(68/19)*
          ei^4*(\[Delta]*(-91579519882877550261559 + 
             1070089805632513857206652*\[Nu])*\[Chi]A + 
           (-91579519882877550261559 - 1886006275288386965668792*\[Nu] + 
             3593907680873728598492112*\[Nu]^2)*\[Chi]S)))/
       (20024571595902815256811677941760*Sqrt[2]*ei^(42/19))) + 
    \[Epsilon]^3*(((377*(-e^(30/19) + e^(12/19)*ei^(18/19))*Pi)/
         (114912*Sqrt[2]*ei^(30/19)) - (e^(12/19)*(50110840*e^(56/19) - 
           54941788*e^2*ei^(18/19) + 4629825*e^(26/19)*ei^(30/19) - 
           20261973*e^(18/19)*ei^2 + 20463096*ei^(56/19))*Pi)/
         (5309853696*Sqrt[2]*ei^(30/19)) + 
        (e^(12/19)*(-18539663378716480*e^(94/19) + 19990448481358576*e^4*
            ei^(18/19) - 2650736640625155*e^(64/19)*ei^(30/19) + 
           11676461601504180*e^(56/19)*ei^2 - 12929209445939952*e^2*
            ei^(56/19) + 1772656357202016*e^(18/19)*ei^4 + 
           680043025216815*ei^(94/19))*Pi)/(1063748393259761664*Sqrt[2]*
          ei^(30/19)) - (e^(12/19)*(118025387066462373214080*e^(132/19) - 
           126468116523637145246444*e^6*ei^(18/19) + 20144212133822032994085*
            e^(102/19)*ei^(30/19) - 88930442105228930366160*e^(94/19)*ei^2 + 
           96841331889577470919584*e^4*ei^(56/19) - 21029217357686808029760*
            e^(56/19)*ei^4 - 8845167348060810772005*e^2*ei^(94/19) + 
           9468844165436914242140*e^(18/19)*ei^6 + 793168079314903044480*
            ei^(132/19))*Pi)/(4386966453700425727082496*Sqrt[2]*ei^(30/19)))*
       x^(5/2) + SO*x^(5/2)*((-8*e^(12/19)*ei^(18/19)*(157*\[Delta]*\[Chi]A + 
            (157 - 110*\[Nu])*\[Chi]S) - 171*ei^(30/19)*
           (23*\[Delta]*\[Chi]A + (23 + 68*\[Nu])*\[Chi]S) + 
          e^(30/19)*(5189*\[Delta]*\[Chi]A + (5189 + 10748*\[Nu])*\[Chi]S))/
         (689472*Sqrt[2]*ei^(30/19)) + 
        (e^(12/19)*(-26*e^2*ei^(18/19)*(92307509*\[Delta]*\[Chi]A + 
             (92307509 - 26680576*\[Nu])*\[Chi]S) - 308655*e^(26/19)*
            ei^(30/19)*(4271*\[Delta]*\[Chi]A + (4271 - 380058*\[Nu])*
              \[Chi]S) + 3528120*ei^(56/19)*(157*\[Delta]*\[Chi]A + 
             (157 - 110*\[Nu])*\[Chi]S) + 1079975*e^(56/19)*
            (5189*\[Delta]*\[Chi]A + (5189 + 10748*\[Nu])*\[Chi]S) - 
           78*e^(18/19)*ei^2*(31277492*\[Delta]*\[Chi]A + 
             (31277492 + 1656666347*\[Nu])*\[Chi]S)))/(258855367680*Sqrt[2]*
          ei^(30/19)) + (e^(12/19)*(-11960*e^4*ei^(18/19)*
            (6334364140217*\[Delta]*\[Chi]A + (6334364140217 - 1267078023544*
                \[Nu])*\[Chi]S) - 16156545975*e^(64/19)*ei^(30/19)*
            (408755*\[Delta]*\[Chi]A + (408755 - 335858562*\[Nu])*\[Chi]S) + 
           470651208*e^2*ei^(56/19)*(92307509*\[Delta]*\[Chi]A + 
             (92307509 - 26680576*\[Nu])*\[Chi]S) - 9019138265475*ei^(94/19)*
            (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) + 
           30735516211400*e^(94/19)*(5189*\[Delta]*\[Chi]A + 
             (5189 + 10748*\[Nu])*\[Chi]S) - 3457647960*e^(56/19)*ei^2*
            (31277492*\[Delta]*\[Chi]A + (31277492 + 1656666347*\[Nu])*
              \[Chi]S) - 156*e^(18/19)*ei^4*(70548776005597*\[Delta]*
              \[Chi]A + (70548776005597 + 205639257894652*\[Nu])*\[Chi]S)))/
         (3989056474724106240*Sqrt[2]*ei^(30/19)) + 
        (e^(12/19)*(-142324*e^6*ei^(18/19)*(477199721661585533*\[Delta]*
              \[Chi]A + (477199721661585533 - 76702103546620216*\[Nu])*
              \[Chi]S) + 7725171552*e^4*ei^(56/19)*(6334364140217*\[Delta]*
              \[Chi]A + (6334364140217 - 1267078023544*\[Nu])*\[Chi]S) - 
           42931098143661*e^2*ei^(94/19)*(92307509*\[Delta]*\[Chi]A + 
             (92307509 - 26680576*\[Nu])*\[Chi]S) + 1402596073059068160*
            ei^(132/19)*(157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) + 
           26088723931578773920*e^(132/19)*(5189*\[Delta]*\[Chi]A + 
             (5189 + 10748*\[Nu])*\[Chi]S) - 3511225371990336*e^(94/19)*ei^2*
            (31277492*\[Delta]*\[Chi]A + (31277492 + 1656666347*\[Nu])*
              \[Chi]S) + 388834206465*e^(102/19)*ei^(30/19)*
            (20119529201*\[Delta]*\[Chi]A + (20119529201 + 14321836977536*
                \[Nu])*\[Chi]S) - 246752688*e^(56/19)*ei^4*
            (70548776005597*\[Delta]*\[Chi]A + (70548776005597 + 
               205639257894652*\[Nu])*\[Chi]S) + 60996*e^(18/19)*ei^6*
            (110801570912666765*\[Delta]*\[Chi]A + (110801570912666765 + 
               267721377018031016*\[Nu])*\[Chi]S)))/
         (2193483226850212863541248*Sqrt[2]*ei^(30/19)))) + 
    \[Epsilon]^4*(x^3*(((e^(12/19) - ei^(12/19))*
          (361*ei^(24/19)*(24215523937 - 140432459328*\[Nu] + 
             53657768304*\[Nu]^2) + (e*ei)^(12/19)*(9356738247901 - 
             51334739449056*\[Nu] + 18622282997808*\[Nu]^2) + 
           e^(24/19)*(17153749047583 - 97253461569600*\[Nu] + 
             78469874452368*\[Nu]^2)))/(317231340281856*Sqrt[2]*ei^(36/19)) + 
        (e^(12/19)*(-13583262*ei^(62/19)*(-358353209 + 372157128*\[Nu] + 
             435997296*\[Nu]^2) + 23261*e^(62/19)*(17153749047583 - 
             97253461569600*\[Nu] + 78469874452368*\[Nu]^2) + 
           260642*e^(26/19)*ei^(36/19)*(43165398517289 - 207702668397672*
              \[Nu] + 180915951308400*\[Nu]^2) + 264*e^2*ei^(24/19)*
            (-77776941858793 - 57290524699128*\[Nu] + 226944607584240*
              \[Nu]^2) - 441*e^(50/19)*ei^(12/19)*(743177219125107 - 
             3998012304886168*\[Nu] + 4211253658313520*\[Nu]^2) + 
           108*e^(12/19)*ei^(50/19)*(5266572853677533 - 33270208807235680*
              \[Nu] + 44813195979425328*\[Nu]^2) - 12*e^(24/19)*ei^2*
            (989592697031477957 - 4854040883091542142*\[Nu] + 
             4334686378660131480*\[Nu]^2)))/(2137716258379333632*Sqrt[2]*
          ei^(36/19)) + (e^(12/19)*(-3664860522606082800*ei^(138/19)*
            (-358353209 + 372157128*\[Nu] + 435997296*\[Nu]^2) + 
           493969332440565800*e^(138/19)*(17153749047583 - 97253461569600*
              \[Nu] + 78469874452368*\[Nu]^2) + 295846665832575*e^2*
            ei^(100/19)*(-77776941858793 - 57290524699128*\[Nu] + 
             226944607584240*\[Nu]^2) - 508462372182960*e^(100/19)*ei^2*
            (989592697031477957 - 4854040883091542142*\[Nu] + 
             4334686378660131480*\[Nu]^2) - 71253911520*e^4*ei^(62/19)*
            (-8259865818436428113 - 2803675687047024972*\[Nu] + 
             16063790724563400000*\[Nu]^2) - 51105600*e^(12/19)*ei^(126/19)*
            (4804915716353577432081 - 18843682694699911068412*\[Nu] + 
             18474139464482338145600*\[Nu]^2) + 16983563041*e^(102/19)*
            ei^(36/19)*(27796008920359825953091 - 133664618078343013891320*
              \[Nu] + 115642064393010874096560*\[Nu]^2) - 39656682*e^(62/19)*
            ei^4*(510435231187790895672403 - 2602527967770865254371160*
              \[Nu] + 1573115806240005997897200*\[Nu]^2) + 
           2698920*e^(50/19)*ei^(88/19)*(1100935288857247994770887 - 
             3563169034924665329533864*\[Nu] + 2847036794122806504723504*
              \[Nu]^2) + 5370300*e^(88/19)*ei^(50/19)*
            (9480459680820722404130209 - 54231424467663871918439768*\[Nu] + 
             55938546743552688123858768*\[Nu]^2) - 319410*e^(126/19)*
            ei^(12/19)*(27473955041458422392876695 - 
             144385493290682670414243352*\[Nu] + 142223936173959144505444080*
              \[Nu]^2) + 2288*e^6*ei^(24/19)*(-497250243303835832331909187 - 
             133308866202460814417262984*\[Nu] + 849484252524317627474122320*
              \[Nu]^2) + 4*e^(24/19)*ei^6*(-357063787763001955513370276971 + 
             183754669839655106920004215848*\[Nu] + 
             262523015186195756815360246800*\[Nu]^2)))/
         (12293946129718445048972087132160*Sqrt[2]*ei^(36/19)) + 
        (e^(12/19)*(13535468371425*ei^(100/19)*(-358353209 + 
             372157128*\[Nu] + 435997296*\[Nu]^2) + 298243389080*e^(100/19)*
            (17153749047583 - 97253461569600*\[Nu] + 78469874452368*
              \[Nu]^2) + 329321167*e^(64/19)*ei^(36/19)*(700626752131304143 - 
             3371572907848553400*\[Nu] + 2927172440774116080*\[Nu]^2) - 
           248826240*e^(62/19)*ei^2*(989592697031477957 - 4854040883091542142*
              \[Nu] + 4334686378660131480*\[Nu]^2) + 63360*e^4*ei^(24/19)*
            (-8259865818436428113 - 2803675687047024972*\[Nu] + 
             16063790724563400000*\[Nu]^2) + 30240*e^(12/19)*ei^(88/19)*
            (16453580776247012131 - 61637476360362405440*\[Nu] + 
             57635656389302047056*\[Nu]^2) - 245700*e^(88/19)*ei^(12/19)*
            (19993660941849726443 - 105813828266912291512*\[Nu] + 
             106417934767940293680*\[Nu]^2) - 12*e^(24/19)*ei^4*
            (510435231187790895672403 - 2602527967770865254371160*\[Nu] + 
             1573115806240005997897200*\[Nu]^2) + 2880*e^2*
            (-646822*ei^(62/19)*(-77776941858793 - 57290524699128*\[Nu] + 
               226944607584240*\[Nu]^2) + 21*ei^2*(e*ei)^(12/19)*
              (352394775629730672441 - 2046564335232903537272*\[Nu] + 
               2213643876177691727952*\[Nu]^2))))/(12841347072734992300769280*
          Sqrt[2]*ei^(36/19))) + SO^2*x^3*
       ((-2*e^(36/19)*(723*\[Delta]*\[Kappa]A + 723*\[Kappa]S - 
            1446*\[Kappa]S*\[Nu] + 790*\[Chi]A^2 - 2892*\[Nu]*\[Chi]A^2 + 
            1580*\[Delta]*\[Chi]A*\[Chi]S + 790*\[Chi]S^2 - 
            268*\[Nu]*\[Chi]S^2) + 3*e^(12/19)*ei^(24/19)*
           (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
            191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
             \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) + 
          19*ei^(36/19)*(48*\[Delta]*\[Kappa]A + 48*\[Kappa]S - 
            96*\[Kappa]S*\[Nu] + 53*\[Chi]A^2 - 192*\[Nu]*\[Chi]A^2 + 
            106*\[Delta]*\[Chi]A*\[Chi]S + 53*\[Chi]S^2 - 
            20*\[Nu]*\[Chi]S^2))/(459648*Sqrt[2]*ei^(36/19)) + 
        (e^(12/19)*(-617421*ei^(62/19)*(178*\[Delta]*\[Kappa]A + 
             178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 
             712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*\[Chi]S + 
             191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) - 1209572*e^(62/19)*
            (-723*\[Kappa]S*(-1 + 2*\[Nu]) + 790*\[Chi]A^2 - 
             2892*\[Nu]*\[Chi]A^2 + 790*\[Chi]S^2 - 268*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(723*\[Kappa]A + 1580*\[Chi]A*\[Chi]S)) - 
           48013*e^(26/19)*ei^(36/19)*(-6564*\[Kappa]S*(-1 + 2*\[Nu]) + 
             59539*\[Chi]A^2 - 26256*\[Nu]*\[Chi]A^2 + 59539*\[Chi]S^2 - 
             211900*\[Nu]*\[Chi]S^2 + 2*\[Delta]*(3282*\[Kappa]A + 59539*
                \[Chi]A*\[Chi]S)) + 156*e^2*ei^(24/19)*
            (-3465533*\[Kappa]S*(-1 + 2*\[Nu]) + 3576925*\[Chi]A^2 - 
             13862132*\[Nu]*\[Chi]A^2 + 3576925*\[Chi]S^2 - 
             445568*\[Nu]*\[Chi]S^2 + \[Delta]*(3465533*\[Kappa]A + 7153850*
                \[Chi]A*\[Chi]S)) + 78*e^(24/19)*ei^2*
            (-9730201*\[Kappa]S*(-1 + 2*\[Nu]) + 43258141*\[Chi]A^2 - 
             38920804*\[Nu]*\[Chi]A^2 + 43258141*\[Chi]S^2 - 
             134111760*\[Nu]*\[Chi]S^2 + \[Delta]*(9730201*\[Kappa]A + 
               86516282*\[Chi]A*\[Chi]S))))/(80532781056*Sqrt[2]*
          ei^(36/19)) + (e^(12/19)*(18930724995*ei^(100/19)*
            (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
             191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
              \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) - 
           477189422528*e^(100/19)*(-723*\[Kappa]S*(-1 + 2*\[Nu]) + 
             790*\[Chi]A^2 - 2892*\[Nu]*\[Chi]A^2 + 790*\[Chi]S^2 - 
             268*\[Nu]*\[Chi]S^2 + \[Delta]*(723*\[Kappa]A + 1580*\[Chi]A*
                \[Chi]S)) - 33869952*e^2*ei^(62/19)*
            (-3465533*\[Kappa]S*(-1 + 2*\[Nu]) + 3576925*\[Chi]A^2 - 
             13862132*\[Nu]*\[Chi]A^2 + 3576925*\[Chi]S^2 - 
             445568*\[Nu]*\[Chi]S^2 + \[Delta]*(3465533*\[Kappa]A + 7153850*
                \[Chi]A*\[Chi]S)) + 49765248*e^(62/19)*ei^2*
            (-9730201*\[Kappa]S*(-1 + 2*\[Nu]) + 43258141*\[Chi]A^2 - 
             38920804*\[Nu]*\[Chi]A^2 + 43258141*\[Chi]S^2 - 
             134111760*\[Nu]*\[Chi]S^2 + \[Delta]*(9730201*\[Kappa]A + 
               86516282*\[Chi]A*\[Chi]S)) - 17332693*e^(64/19)*ei^(36/19)*
            (-18305358*\[Kappa]S*(-1 + 2*\[Nu]) + 112913425*\[Chi]A^2 - 
             73221432*\[Nu]*\[Chi]A^2 + 112913425*\[Chi]S^2 - 
             378432268*\[Nu]*\[Chi]S^2 + 2*\[Delta]*(9152679*\[Kappa]A + 
               112913425*\[Chi]A*\[Chi]S)) + 24*e^(24/19)*ei^4*
            (-2047975353262*\[Kappa]S*(-1 + 2*\[Nu]) + 2117163228955*
              \[Chi]A^2 - 8191901413048*\[Nu]*\[Chi]A^2 + 2117163228955*
              \[Chi]S^2 - 276751502772*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(2047975353262*\[Kappa]A + 4234326457910*\[Chi]A*
                \[Chi]S)) + 96*e^4*ei^(24/19)*(-2530412457698*\[Kappa]S*
              (-1 + 2*\[Nu]) + 2583790126097*\[Chi]A^2 - 10121649830792*\[Nu]*
              \[Chi]A^2 + 2583790126097*\[Chi]S^2 - 213510673596*\[Nu]*
              \[Chi]S^2 + \[Delta]*(2530412457698*\[Kappa]A + 5167580252194*
                \[Chi]A*\[Chi]S))))/(14885034988142592*Sqrt[2]*ei^(36/19)) + 
        (e^(12/19)*(-25628395262979600*ei^(138/19)*(178*\[Delta]*\[Kappa]A + 
             178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 
             712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*\[Chi]S + 
             191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) - 3951754659524526400*
            e^(138/19)*(-723*\[Kappa]S*(-1 + 2*\[Nu]) + 790*\[Chi]A^2 - 
             2892*\[Nu]*\[Chi]A^2 + 790*\[Chi]S^2 - 268*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(723*\[Kappa]A + 1580*\[Chi]A*\[Chi]S)) + 
           26895151439325*e^2*ei^(100/19)*(-3465533*\[Kappa]S*
              (-1 + 2*\[Nu]) + 3576925*\[Chi]A^2 - 13862132*\[Nu]*\[Chi]A^2 + 
             3576925*\[Chi]S^2 - 445568*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(3465533*\[Kappa]A + 7153850*\[Chi]A*\[Chi]S)) + 
           508462372182960*e^(100/19)*ei^2*(-9730201*\[Kappa]S*
              (-1 + 2*\[Nu]) + 43258141*\[Chi]A^2 - 38920804*\[Nu]*
              \[Chi]A^2 + 43258141*\[Chi]S^2 - 134111760*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(9730201*\[Kappa]A + 86516282*\[Chi]A*\[Chi]S)) - 
           31285510865*e^(102/19)*ei^(36/19)*(-114560183793*\[Kappa]S*
              (-1 + 2*\[Nu]) + 651662338765*\[Chi]A^2 - 458240735172*\[Nu]*
              \[Chi]A^2 + 651662338765*\[Chi]S^2 - 2148408619888*\[Nu]*
              \[Chi]S^2 + \[Delta]*(114560183793*\[Kappa]A + 1303324677530*
                \[Chi]A*\[Chi]S)) + 396566820*e^(62/19)*ei^4*
            (-2047975353262*\[Kappa]S*(-1 + 2*\[Nu]) + 2117163228955*
              \[Chi]A^2 - 8191901413048*\[Nu]*\[Chi]A^2 + 2117163228955*
              \[Chi]S^2 - 276751502772*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(2047975353262*\[Kappa]A + 4234326457910*\[Chi]A*
                \[Chi]S)) - 539802360*e^4*ei^(62/19)*
            (-2530412457698*\[Kappa]S*(-1 + 2*\[Nu]) + 2583790126097*
              \[Chi]A^2 - 10121649830792*\[Nu]*\[Chi]A^2 + 2583790126097*
              \[Chi]S^2 - 213510673596*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(2530412457698*\[Kappa]A + 5167580252194*\[Chi]A*
                \[Chi]S)) + 14144*e^6*ei^(24/19)*(-151862337793381386*
              \[Kappa]S*(-1 + 2*\[Nu]) + 154258072076881051*\[Chi]A^2 - 
             607449351173525544*\[Nu]*\[Chi]A^2 + 154258072076881051*
              \[Chi]S^2 - 9582937133998660*\[Nu]*\[Chi]S^2 + 
             2*\[Delta]*(75931168896690693*\[Kappa]A + 154258072076881051*
                \[Chi]A*\[Chi]S)) - 1088*e^(24/19)*ei^6*
            (-173777694647922993*\[Kappa]S*(-1 + 2*\[Nu]) + 
             187183269443850718*\[Chi]A^2 - 695110778591691972*\[Nu]*
              \[Chi]A^2 + 187183269443850718*\[Chi]S^2 - 53622299183710900*
              \[Nu]*\[Chi]S^2 + \[Delta]*(173777694647922993*\[Kappa]A + 
               374366538887701436*\[Chi]A*\[Chi]S))))/
         (71252578583340247941120*Sqrt[2]*ei^(36/19)))) + 
    \[Epsilon]^6*
     (SO*x^4*((Pi*(162450*ei^(48/19)*(\[Delta]*\[Chi]A + \[Chi]S - 
             4*\[Nu]*\[Chi]S) + 1885*e^(30/19)*ei^(18/19)*
            (6445*\[Delta]*\[Chi]A + (6445 + 9868*\[Nu])*\[Chi]S) + 
           2*e^(12/19)*ei^(36/19)*(31337*\[Delta]*\[Chi]A + 
             (31337 + 235316*\[Nu])*\[Chi]S) - 7*e^(48/19)*
            (1767707*\[Delta]*\[Chi]A + (1767707 + 2631716*\[Nu])*\[Chi]S)))/
         (314399232*Sqrt[2]*ei^(48/19)) - 
        (e^(12/19)*Pi*(-39*e^(36/19)*ei^2*(60799062235*\[Delta]*\[Chi]A + 
             (60799062235 - 1681766505644*\[Nu])*\[Chi]S) + 
           352812*ei^(74/19)*(31337*\[Delta]*\[Chi]A + (31337 + 235316*\[Nu])*
              \[Chi]S) + 4838288*e^(74/19)*(1767707*\[Delta]*\[Chi]A + 
             (1767707 + 2631716*\[Nu])*\[Chi]S) - 4691556*e^(26/19)*
            ei^(48/19)*(476311*\[Delta]*\[Chi]A + (476311 + 34775336*\[Nu])*
              \[Chi]S) + 52*e^2*ei^(36/19)*(65005615159*\[Delta]*\[Chi]A + 
             (65005615159 + 28798683196*\[Nu])*\[Chi]S) - 
           91*e^(56/19)*ei^(18/19)*(105398452691*\[Delta]*\[Chi]A + 
             (105398452691 + 152682705716*\[Nu])*\[Chi]S) + 
           78*e^(18/19)*ei^(56/19)*(28885111903*\[Delta]*\[Chi]A + 
             (28885111903 + 1245411730588*\[Nu])*\[Chi]S)))/
         (47215219064832*Sqrt[2]*ei^(48/19)) - 
        (e^(12/19)*Pi*(-2737*e^(36/19)*ei^6*(1409952906044235180650719*
              \[Delta]*\[Chi]A + (1409952906044235180650719 - 
               266156392897263956221052*\[Nu])*\[Chi]S) - 1700148019554784800*
            e^(112/19)*ei^2*(60799062235*\[Delta]*\[Chi]A + 
             (60799062235 - 1681766505644*\[Nu])*\[Chi]S) + 
           57857088013686561600*ei^(150/19)*(31337*\[Delta]*\[Chi]A + 
             (31337 + 235316*\[Nu])*\[Chi]S) + 96519641913542076028800*
            e^(150/19)*(1767707*\[Delta]*\[Chi]A + (1767707 + 2631716*\[Nu])*
              \[Chi]S) + 35418155968520325*e^2*ei^(112/19)*
            (65005615159*\[Delta]*\[Chi]A + (65005615159 + 28798683196*\[Nu])*
              \[Chi]S) + 4922927627904*e^(74/19)*ei^4*
            (8388764562288913*\[Delta]*\[Chi]A + (8388764562288913 + 
               50345437127831536*\[Nu])*\[Chi]S) - 36939249614175*e^(102/19)*
            ei^(48/19)*(2774130035414197*\[Delta]*\[Chi]A + 
             (2774130035414197 + 208064117014088756*\[Nu])*\[Chi]S) - 
           127465330608*e^4*ei^(74/19)*(467358287836566959*\[Delta]*\[Chi]A + 
             (467358287836566959 + 252002971654484168*\[Nu])*\[Chi]S) + 
           55135080*e^(56/19)*ei^(94/19)*(244623349738220609297*\[Delta]*
              \[Chi]A + (244623349738220609297 + 454896098654954459292*\[Nu])*
              \[Chi]S) - 5591300*e^(18/19)*ei^(132/19)*
            (522183109407641892175*\[Delta]*\[Chi]A + 
             (522183109407641892175 + 805922124431695215388*\[Nu])*\[Chi]S) - 
           9784775*e^(132/19)*ei^(18/19)*(20978511342636109717999*\[Delta]*
              \[Chi]A + (20978511342636109717999 + 29184542459338607836420*
                \[Nu])*\[Chi]S) + 140900760*e^(94/19)*ei^(56/19)*
            (829523661780365350411*\[Delta]*\[Chi]A + 
             (829523661780365350411 + 32243072397146635524676*\[Nu])*
              \[Chi]S) + 120428*e^6*ei^(36/19)*(1103353378842183123126601*
              \[Delta]*\[Chi]A + (1103353378842183123126601 + 
               638529972502296236631652*\[Nu])*\[Chi]S)))/
         (165037677988210015852843499520*Sqrt[2]*ei^(48/19)) - 
        (e^(12/19)*Pi*(-193628285760*e^(74/19)*ei^2*
            (60799062235*\[Delta]*\[Chi]A + (60799062235 - 1681766505644*
                \[Nu])*\[Chi]S) - 63133967858325*ei^(112/19)*
            (31337*\[Delta]*\[Chi]A + (31337 + 235316*\[Nu])*\[Chi]S) + 
           14142645438958400*e^(112/19)*(1767707*\[Delta]*\[Chi]A + 
             (1767707 + 2631716*\[Nu])*\[Chi]S) - 61394874705*e^(64/19)*
            ei^(48/19)*(186875094583*\[Delta]*\[Chi]A + 
             (186875094583 + 13970833203308*\[Nu])*\[Chi]S) + 
           196560*e^(18/19)*ei^(94/19)*(3946331805444981*\[Delta]*\[Chi]A + 
             (3946331805444981 + 7919583360996556*\[Nu])*\[Chi]S) + 
           330096*e^(36/19)*ei^4*(8388764562288913*\[Delta]*\[Chi]A + 
             (8388764562288913 + 50345437127831536*\[Nu])*\[Chi]S) - 
           167440*e^(94/19)*ei^(18/19)*(175515568865658611*\[Delta]*\[Chi]A + 
             (175515568865658611 + 247878349535282036*\[Nu])*\[Chi]S) + 
           33488*e^4*ei^(36/19)*(467358287836566959*\[Delta]*\[Chi]A + 
             (467358287836566959 + 252002971654484168*\[Nu])*\[Chi]S) + 
           280140*e^2*(-235208*ei^(74/19)*(65005615159*\[Delta]*\[Chi]A + 
               65005615159*\[Chi]S + 28798683196*\[Nu]*\[Chi]S) + 
             13*ei^2*(e*ei)^(18/19)*(3493708141898143*\[Delta]*\[Chi]A + 
               3493708141898143*\[Chi]S + 140502049022356708*\[Nu]*
                \[Chi]S))))/(50932273069277388472320*Sqrt[2]*ei^(48/19))) + 
      SO^2*x^4*((108*(e*ei)^(24/19)*(-3920527 + 15455580*\[Nu])*
           (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
            191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
             \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) + 
          e^(48/19)*(-504*\[Delta]*\[Kappa]A*(-133053743 + 928541756*\[Nu]) + 
            504*\[Kappa]S*(133053743 - 1194649242*\[Nu] + 2313892152*\[Nu]^
                2) - 105124988001*\[Chi]A^2 + 272322040144*\[Delta]^2*
             \[Chi]A^2 - 139352970876*\[Nu]*\[Chi]A^2 + 2332403289216*\[Nu]^2*
             \[Chi]A^2 + 2*\[Delta]*(167197052143 - 251207045852*\[Nu])*
             \[Chi]A*\[Chi]S + 167197052143*\[Chi]S^2 + 57438831176*\[Nu]*
             \[Chi]S^2 - 75366630416*\[Nu]^2*\[Chi]S^2) - 
          16245*ei^(48/19)*(-96*\[Delta]*\[Kappa]A*(1979 + 89740*\[Nu]) + 
            96*\[Kappa]S*(-1979 - 85782*\[Nu] + 296072*\[Nu]^2) - 
            8023897*\[Chi]A^2 + 9605904*\[Delta]^2*\[Chi]A^2 + 
            19042564*\[Nu]*\[Chi]A^2 + 56845824*\[Nu]^2*\[Chi]A^2 + 
            14*\[Delta]*(226001 + 180188*\[Nu])*\[Chi]A*\[Chi]S + 
            1582007*\[Chi]S^2 + 15575656*\[Nu]*\[Chi]S^2 - 
            3065104*\[Nu]^2*\[Chi]S^2) - 98560*e^(30/19)*ei^(18/19)*
           (814673*\[Delta]^2*\[Chi]A^2 + 2*\[Delta]*(814673 + 558323*\[Nu])*
             \[Chi]A*\[Chi]S + (814673 + 1116646*\[Nu] - 1182280*\[Nu]^2)*
             \[Chi]S^2) + 19008*e^(36/19)*ei^(12/19)*(-2833 + 5516*\[Nu])*
           (-723*\[Kappa]S*(-1 + 2*\[Nu]) + 790*\[Chi]A^2 - 
            2892*\[Nu]*\[Chi]A^2 + 790*\[Chi]S^2 - 268*\[Nu]*\[Chi]S^2 + 
            \[Delta]*(723*\[Kappa]A + 1580*\[Chi]A*\[Chi]S)) - 
          352*e^(12/19)*ei^(36/19)*(-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 
              116610*\[Nu]^2) - 278950701*\[Chi]A^2 + 102215792*\[Delta]^2*
             \[Chi]A^2 + 1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*
             \[Chi]A^2 - 176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
            68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*(-124448 + 
                126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*\[Chi]A*
               \[Chi]S)))/(4648078245888*Sqrt[2]*ei^(48/19)) + 
        (e^(12/19)*(-810*e^(12/19)*ei^(62/19)*(-1859009125901 + 
             8124219720708*\[Nu])*(178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 
             356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 
             382*\[Delta]*\[Chi]A*\[Chi]S + 191*\[Chi]S^2 - 
             52*\[Nu]*\[Chi]S^2) + 1511965*e^(74/19)*
            (-504*\[Delta]*\[Kappa]A*(-133053743 + 928541756*\[Nu]) + 
             504*\[Kappa]S*(133053743 - 1194649242*\[Nu] + 2313892152*
                \[Nu]^2) - 105124988001*\[Chi]A^2 + 272322040144*\[Delta]^2*
              \[Chi]A^2 - 139352970876*\[Nu]*\[Chi]A^2 + 2332403289216*
              \[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(167197052143 - 251207045852*
                \[Nu])*\[Chi]A*\[Chi]S + 167197052143*\[Chi]S^2 + 
             57438831176*\[Nu]*\[Chi]S^2 - 75366630416*\[Nu]^2*\[Chi]S^2) - 
           280280*e^(56/19)*ei^(18/19)*(641413166941*\[Delta]^2*\[Chi]A^2 + 
             2*\[Delta]*(641413166941 + 538156238674*\[Nu])*\[Chi]A*\[Chi]S + 
             (641413166941 + 1076312477348*\[Nu] - 522485817248*\[Nu]^2)*
              \[Chi]S^2) + 6726720*e^(18/19)*ei^(56/19)*
            (4910566244*\[Delta]^2*\[Chi]A^2 + \[Delta]*(9821132488 + 
               256656092359*\[Nu])*\[Chi]A*\[Chi]S + 
             (4910566244 + 256656092359*\[Nu] - 182233298170*\[Nu]^2)*
              \[Chi]S^2) + 83160*e^(62/19)*ei^(12/19)*(-1571689321 + 
             2383893876*\[Nu])*(-723*\[Kappa]S*(-1 + 2*\[Nu]) + 
             790*\[Chi]A^2 - 2892*\[Nu]*\[Chi]A^2 + 790*\[Chi]S^2 - 
             268*\[Nu]*\[Chi]S^2 + \[Delta]*(723*\[Kappa]A + 1580*\[Chi]A*
                \[Chi]S)) + 1755*e^(50/19)*ei^(24/19)*(-3920527 + 
             15455580*\[Nu])*(-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 
             25904401*\[Chi]A^2 - 99734624*\[Nu]*\[Chi]A^2 + 
             25904401*\[Chi]S^2 - 3882980*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(24933656*\[Kappa]A + 51808802*\[Chi]A*\[Chi]S)) - 
           1389960*e^(24/19)*ei^(50/19)*(-2833 + 5516*\[Nu])*
            (-9730201*\[Kappa]S*(-1 + 2*\[Nu]) + 43258141*\[Chi]A^2 - 
             38920804*\[Nu]*\[Chi]A^2 + 43258141*\[Chi]S^2 - 
             134111760*\[Nu]*\[Chi]S^2 + \[Delta]*(9730201*\[Kappa]A + 
               86516282*\[Chi]A*\[Chi]S)) + 135832620*ei^(74/19)*
            (-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 116610*\[Nu]^2) - 
             278950701*\[Chi]A^2 + 102215792*\[Delta]^2*\[Chi]A^2 + 
             1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*\[Chi]A^2 - 
             176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
             68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*
                (-124448 + 126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*
                \[Chi]A*\[Chi]S)) - 4691556*e^(26/19)*ei^(48/19)*
            (30*\[Kappa]S*(6701235367 - 28461012890*\[Nu] + 49884658824*
                \[Nu]^2) + 75645263203*\[Chi]A^2 + 164941658112*\[Delta]^2*
              \[Chi]A^2 - 954564658552*\[Nu]*\[Chi]A^2 + 2993079529440*
              \[Nu]^2*\[Chi]A^2 + 240586921315*\[Chi]S^2 + 110028617900*\[Nu]*
              \[Chi]S^2 + 406681327360*\[Nu]^2*\[Chi]S^2 - 
             10*\[Delta]*(3*\[Kappa]A*(-6701235367 + 15058542156*\[Nu]) + 
               (-48117384263 + 54195498784*\[Nu])*\[Chi]A*\[Chi]S)) - 
           220*e^2*ei^(36/19)*(-18*\[Kappa]S*(30936092472295 - 97020967578282*
                \[Nu] + 36669175275096*\[Nu]^2) - 2827049041226343*
              \[Chi]A^2 + 2206814975744288*\[Delta]^2*\[Chi]A^2 + 
             11721265129529220*\[Nu]*\[Chi]A^2 - 1320090309903456*\[Nu]^2*
              \[Chi]A^2 - 620234065482055*\[Chi]S^2 + 1521486691850320*\[Nu]*
              \[Chi]S^2 - 141962774214352*\[Nu]^2*\[Chi]S^2 + 
             2*\[Delta]*(9*\[Kappa]A*(-30936092472295 + 35148782633692*
                  \[Nu]) + (-620234065482055 + 967277828237084*\[Nu])*\[Chi]A*
                \[Chi]S)) + 12*e^(36/19)*ei^2*(45*\[Kappa]S*
              (1288935001706461 - 4177294463778654*\[Nu] + 7918493972374440*
                \[Nu]^2) - 20556501925266276*\[Chi]A^2 + 81704163869614286*
              \[Delta]^2*\[Chi]A^2 - 79734899027872416*\[Nu]*\[Chi]A^2 + 
             712664457513699600*\[Nu]^2*\[Chi]A^2 + 61147661944348010*
              \[Chi]S^2 + 46845577848720400*\[Nu]*\[Chi]S^2 + 
             155611899223682720*\[Nu]^2*\[Chi]S^2 - 5*\[Delta]*
              (9*\[Kappa]A*(-1288935001706461 + 1599424460365732*\[Nu]) + 4*
                (-6114766194434801 + 5755766444010856*\[Nu])*\[Chi]A*
                \[Chi]S))))/(1526940184556666880*Sqrt[2]*ei^(48/19)) + 
        (e^(12/19)*(201902500800*e^(12/19)*ei^(138/19)*
            (-1696052141317888257 + 3349191345990271600*\[Nu])*
            (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
             191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
              \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) - 
           623392*e^6*ei^(36/19)*(33408*\[Delta]*\[Kappa]A*
              (-138222143243089778236681 + 169677614533377548972647*\[Nu]) - 
             33408*\[Kappa]S*(138222143243089778236681 - 
               446121901019557105446009*\[Nu] + 187902755346189154082994*
                \[Nu]^2) - 30921703483728105055475662749*\[Chi]A^2 + 
             27011676587914586707090685734*\[Delta]^2*\[Chi]A^2 + 
             128203625765713966307445860916*\[Nu]*\[Chi]A^2 - 
             12554910501210974519209327104*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-3910026895813518348384977015 + 
               7969183702013203055213796514*\[Nu])*\[Chi]A*\[Chi]S - 
             3910026895813518348384977015*\[Chi]S^2 + 
             11421555573224860024884383108*\[Nu]*\[Chi]S^2 - 
             230627307042699097000726256*\[Nu]^2*\[Chi]S^2) + 
           467372878896*e^4*ei^(74/19)*(522*\[Delta]*\[Kappa]A*
              (-5919500807162614699 + 7069747731165174748*\[Nu]) - 
             522*\[Kappa]S*(5919500807162614699 - 18908749345490404146*
                \[Nu] + 7669345784813238456*\[Nu]^2) - 
             18928431938173283618721*\[Chi]A^2 + 16007908125685966389172*
              \[Delta]^2*\[Chi]A^2 + 78468871437192708323052*\[Nu]*
              \[Chi]A^2 - 8006796999345020948064*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-2920523812487317229549 + 5349738482100710957140*
                \[Nu])*\[Chi]A*\[Chi]S - 2920523812487317229549*\[Chi]S^2 + 
             7944333279701848066112*\[Nu]*\[Chi]S^2 - 384534241038124715888*
              \[Nu]^2*\[Chi]S^2) + 26140736351584312257800*e^(150/19)*
            (-504*\[Delta]*\[Kappa]A*(-133053743 + 928541756*\[Nu]) + 
             504*\[Kappa]S*(133053743 - 1194649242*\[Nu] + 2313892152*
                \[Nu]^2) - 105124988001*\[Chi]A^2 + 272322040144*\[Delta]^2*
              \[Chi]A^2 - 139352970876*\[Nu]*\[Chi]A^2 + 2332403289216*
              \[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(167197052143 - 251207045852*
                \[Nu])*\[Chi]A*\[Chi]S + 167197052143*\[Chi]S^2 + 
             57438831176*\[Nu]*\[Chi]S^2 - 75366630416*\[Nu]^2*\[Chi]S^2) - 
           22163549768505*e^(102/19)*ei^(48/19)*(-138*\[Delta]*\[Kappa]A*
              (-14498752044329278759 + 37701813984195651164*\[Nu]) + 
             138*\[Kappa]S*(14498752044329278759 - 66699318072854208682*
                \[Nu] + 103884308553457687640*\[Nu]^2) + 
             1256983824005985698827*\[Chi]A^2 + 1237683276840027265248*
              \[Delta]^2*\[Chi]A^2 - 10930197089807692178932*\[Nu]*
              \[Chi]A^2 + 28672069160754321788640*\[Nu]^2*\[Chi]A^2 + 
             10*\[Delta]*(498933420169202592815 - 564150551657420161084*
                \[Nu])*\[Chi]A*\[Chi]S + 2494667100846012964075*\[Chi]S^2 + 
             260756277209547772784*\[Nu]*\[Chi]S^2 + 3430102131980339561872*
              \[Nu]^2*\[Chi]S^2) + 3864*e^(36/19)*ei^6*
            (522*\[Delta]*\[Kappa]A*(-306774235889116393472294781 + 
               590981029664619833905078868*\[Nu]) - 522*\[Kappa]S*
              (306774235889116393472294781 - 1204529501442852620849668430*
                \[Nu] + 713711612540245317793330536*\[Nu]^2) - 
             287442032170821380844717949560*\[Chi]A^2 + 
             81396770844384560434563464771*\[Delta]^2*\[Chi]A^2 + 
             1362616475368189466448541378968*\[Nu]*\[Chi]A^2 - 
             745114923492016111776237079584*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-206045261326436820410154484789 + 
               225689763317072717532279234734*\[Nu])*\[Chi]A*\[Chi]S - 
             206045261326436820410154484789*\[Chi]S^2 + 
             238531179949241491994888888740*\[Nu]*\[Chi]S^2 + 
             281445769123207109379502334336*\[Nu]^2*\[Chi]S^2) + 
           2106222134016*e^(94/19)*ei^(56/19)*(1016904766187869381924*
              \[Delta]^2*\[Chi]A^2 + \[Delta]*(2033809532375738763848 + 
               53569363462188306088943*\[Nu])*\[Chi]A*\[Chi]S + 
             (1016904766187869381924 + 53569363462188306088943*\[Nu] - 
               15506146749330476009456*\[Nu]^2)*\[Chi]S^2) + 
           22893718848*e^(56/19)*ei^(94/19)*(8720546124794990115893*
              \[Delta]^2*\[Chi]A^2 + 2*\[Delta]*(8720546124794990115893 + 
               10994779897934329223858*\[Nu])*\[Chi]A*\[Chi]S + 
             (8720546124794990115893 + 21989559795868658447716*\[Nu] - 
               9996612925135791413152*\[Nu]^2)*\[Chi]S^2) - 
           548495347400*e^(132/19)*ei^(18/19)*(7101374901099450873329*
              \[Delta]^2*\[Chi]A^2 + 2*\[Delta]*(7101374901099450873329 + 
               6526258624833549912314*\[Nu])*\[Chi]A*\[Chi]S + 
             (7101374901099450873329 + 13052517249667099824628*\[Nu] - 
               3431310575692919015200*\[Nu]^2)*\[Chi]S^2) - 
           1880555476800*e^(18/19)*ei^(132/19)*(17395846633288682105*
              \[Delta]^2*\[Chi]A^2 + 2*\[Delta]*(17395846633288682105 + 
               14922041695718762681*\[Nu])*\[Chi]A*\[Chi]S + 
             (17395846633288682105 + 29844083391437525362*\[Nu] - 
               29449351471983411760*\[Nu]^2)*\[Chi]S^2) + 
           1776742007040*e^(138/19)*ei^(12/19)*(-1663373241307010153 + 
             2252417164528514620*\[Nu])*(-723*\[Kappa]S*(-1 + 2*\[Nu]) + 
             790*\[Chi]A^2 - 2892*\[Nu]*\[Chi]A^2 + 790*\[Chi]S^2 - 
             268*\[Nu]*\[Chi]S^2 + \[Delta]*(723*\[Kappa]A + 1580*\[Chi]A*
                \[Chi]S)) - 5657707440*e^(50/19)*ei^(100/19)*
            (-5807829430373107 + 10448813703644316*\[Nu])*
            (-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 25904401*\[Chi]A^2 - 
             99734624*\[Nu]*\[Chi]A^2 + 25904401*\[Chi]S^2 - 
             3882980*\[Nu]*\[Chi]S^2 + \[Delta]*(24933656*\[Kappa]A + 
               51808802*\[Chi]A*\[Chi]S)) - 2427156491760*e^(100/19)*
            ei^(50/19)*(-117334502622439 + 164817749118204*\[Nu])*
            (-9730201*\[Kappa]S*(-1 + 2*\[Nu]) + 43258141*\[Chi]A^2 - 
             38920804*\[Nu]*\[Chi]A^2 + 43258141*\[Chi]S^2 - 
             134111760*\[Nu]*\[Chi]S^2 + \[Delta]*(9730201*\[Kappa]A + 
               86516282*\[Chi]A*\[Chi]S)) - 93352172760*e^(62/19)*ei^(88/19)*
            (-1571689321 + 2383893876*\[Nu])*(-2047975353262*\[Kappa]S*
              (-1 + 2*\[Nu]) + 2117163228955*\[Chi]A^2 - 8191901413048*\[Nu]*
              \[Chi]A^2 + 2117163228955*\[Chi]S^2 - 276751502772*\[Nu]*
              \[Chi]S^2 + \[Delta]*(2047975353262*\[Kappa]A + 4234326457910*
                \[Chi]A*\[Chi]S)) + 66383767296*e^(24/19)*ei^(126/19)*
            (-2833 + 5516*\[Nu])*(-173777694647922993*\[Kappa]S*
              (-1 + 2*\[Nu]) + 187183269443850718*\[Chi]A^2 - 
             695110778591691972*\[Nu]*\[Chi]A^2 + 187183269443850718*
              \[Chi]S^2 - 53622299183710900*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(173777694647922993*\[Kappa]A + 374366538887701436*
                \[Chi]A*\[Chi]S)) + 817224408*e^(126/19)*ei^(24/19)*
            (-3920527 + 15455580*\[Nu])*(-1165556009138332944*\[Kappa]S*
              (-1 + 2*\[Nu]) + 1189984978416247669*\[Chi]A^2 - 
             4662224036553331776*\[Nu]*\[Chi]A^2 + 1189984978416247669*
              \[Chi]S^2 - 97715877111658900*\[Nu]*\[Chi]S^2 + 
             2*\[Delta]*(582778004569166472*\[Kappa]A + 1189984978416247669*
                \[Chi]A*\[Chi]S)) + 19304981700566749387200*ei^(150/19)*
            (-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 116610*\[Nu]^2) - 
             278950701*\[Chi]A^2 + 102215792*\[Delta]^2*\[Chi]A^2 + 
             1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*\[Chi]A^2 - 
             176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
             68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*
                (-124448 + 126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*
                \[Chi]A*\[Chi]S)) - 129866571884574525*e^2*ei^(112/19)*
            (-18*\[Kappa]S*(30936092472295 - 97020967578282*\[Nu] + 
               36669175275096*\[Nu]^2) - 2827049041226343*\[Chi]A^2 + 
             2206814975744288*\[Delta]^2*\[Chi]A^2 + 11721265129529220*\[Nu]*
              \[Chi]A^2 - 1320090309903456*\[Nu]^2*\[Chi]A^2 - 
             620234065482055*\[Chi]S^2 + 1521486691850320*\[Nu]*\[Chi]S^2 - 
             141962774214352*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(9*\[Kappa]A*
                (-30936092472295 + 35148782633692*\[Nu]) + 
               (-620234065482055 + 967277828237084*\[Nu])*\[Chi]A*\[Chi]S)) + 
           453372805214609280*e^(112/19)*ei^2*(45*\[Kappa]S*
              (1288935001706461 - 4177294463778654*\[Nu] + 7918493972374440*
                \[Nu]^2) - 20556501925266276*\[Chi]A^2 + 81704163869614286*
              \[Delta]^2*\[Chi]A^2 - 79734899027872416*\[Nu]*\[Chi]A^2 + 
             712664457513699600*\[Nu]^2*\[Chi]A^2 + 61147661944348010*
              \[Chi]S^2 + 46845577848720400*\[Nu]*\[Chi]S^2 + 
             155611899223682720*\[Nu]^2*\[Chi]S^2 - 5*\[Delta]*
              (9*\[Kappa]A*(-1288935001706461 + 1599424460365732*\[Nu]) + 4*
                (-6114766194434801 + 5755766444010856*\[Nu])*\[Chi]A*
                \[Chi]S)) - 17836694304*e^(74/19)*ei^4*
            (1044*\[Kappa]S*(12865074844809000047 + 203738894213893801302*
                \[Nu] + 319207505993389781208*\[Nu]^2) - 
             269552471721024277111089*\[Chi]A^2 + 326073634913301114343504*
              \[Delta]^2*\[Chi]A^2 + 1115745822132711025804788*\[Nu]*
              \[Chi]A^2 + 666505272514197863162304*\[Nu]^2*\[Chi]A^2 + 
             56521163192276837232415*\[Chi]S^2 + 43212084697850551775000*
              \[Nu]*\[Chi]S^2 + 1270802617793555061375664*\[Nu]^2*\[Chi]S^2 + 
             2*\[Delta]*(522*\[Kappa]A*(12865074844809000047 + 
                 229469043903511801396*\[Nu]) + (56521163192276837232415 + 
                 40374009973232234567716*\[Nu])*\[Chi]A*\[Chi]S)) - 
           3232975680*e^(88/19)*ei^(62/19)*(-1859009125901 + 
             8124219720708*\[Nu])*(-2948986828885*\[Kappa]S*(-1 + 2*\[Nu]) + 
             \[Delta]*(2948986828885*\[Kappa]A + 6058122836776*\[Chi]A*
                \[Chi]S) - 4*((-757265354597 + 2948986828885*\[Nu])*
                \[Chi]A^2 + (-757265354597 + 80074589503*\[Nu])*\[Chi]S^2))))/
         (4625676038653550324323497604546560*Sqrt[2]*ei^(48/19)) + 
        (e^(12/19)*(-105892920*e^(12/19)*ei^(100/19)*(-5807829430373107 + 
             10448813703644316*\[Nu])*(178*\[Delta]*\[Kappa]A + 
             178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 
             712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*\[Chi]S + 
             191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) - 368368*e^4*ei^(36/19)*
            (522*\[Delta]*\[Kappa]A*(-5919500807162614699 + 
               7069747731165174748*\[Nu]) - 522*\[Kappa]S*
              (5919500807162614699 - 18908749345490404146*\[Nu] + 
               7669345784813238456*\[Nu]^2) - 18928431938173283618721*
              \[Chi]A^2 + 16007908125685966389172*\[Delta]^2*\[Chi]A^2 + 
             78468871437192708323052*\[Nu]*\[Chi]A^2 - 8006796999345020948064*
              \[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(-2920523812487317229549 + 
               5349738482100710957140*\[Nu])*\[Chi]A*\[Chi]S - 
             2920523812487317229549*\[Chi]S^2 + 7944333279701848066112*\[Nu]*
              \[Chi]S^2 - 384534241038124715888*\[Nu]^2*\[Chi]S^2) + 
           11490899419153700*e^(112/19)*(-504*\[Delta]*\[Kappa]A*
              (-133053743 + 928541756*\[Nu]) + 504*\[Kappa]S*
              (133053743 - 1194649242*\[Nu] + 2313892152*\[Nu]^2) - 
             105124988001*\[Chi]A^2 + 272322040144*\[Delta]^2*\[Chi]A^2 - 
             139352970876*\[Nu]*\[Chi]A^2 + 2332403289216*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(167197052143 - 251207045852*\[Nu])*\[Chi]A*\[Chi]S + 
             167197052143*\[Chi]S^2 + 57438831176*\[Nu]*\[Chi]S^2 - 
             75366630416*\[Nu]^2*\[Chi]S^2) - 9386016640*e^(94/19)*ei^(18/19)*
            (168706584010919233*\[Delta]^2*\[Chi]A^2 + 2*\[Delta]*
              (168706584010919233 + 150437240685157342*\[Nu])*\[Chi]A*
              \[Chi]S + (168706584010919233 + 300874481370314684*\[Nu] - 
               100599656390439104*\[Nu]^2)*\[Chi]S^2) + 2448526080*e^(18/19)*
            ei^(94/19)*(11076157832878729*\[Delta]^2*\[Chi]A^2 + 
             2*\[Delta]*(11076157832878729 + 12262499064422347*\[Nu])*\[Chi]A*
              \[Chi]S + (11076157832878729 + 24524998128844694*\[Nu] - 
               22620318368411720*\[Nu]^2)*\[Chi]S^2) + 10095125040*e^(100/19)*
            ei^(12/19)*(-117334502622439 + 164817749118204*\[Nu])*
            (-723*\[Kappa]S*(-1 + 2*\[Nu]) + 790*\[Chi]A^2 - 
             2892*\[Nu]*\[Chi]A^2 + 790*\[Chi]S^2 - 268*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(723*\[Kappa]A + 1580*\[Chi]A*\[Chi]S)) - 
           112376160*e^(50/19)*ei^(62/19)*(-1859009125901 + 
             8124219720708*\[Nu])*(-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 
             25904401*\[Chi]A^2 - 99734624*\[Nu]*\[Chi]A^2 + 
             25904401*\[Chi]S^2 - 3882980*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(24933656*\[Kappa]A + 51808802*\[Chi]A*\[Chi]S)) - 
           51917785920*e^(62/19)*ei^(50/19)*(-1571689321 + 2383893876*\[Nu])*
            (-9730201*\[Kappa]S*(-1 + 2*\[Nu]) + 43258141*\[Chi]A^2 - 
             38920804*\[Nu]*\[Chi]A^2 + 43258141*\[Chi]S^2 - 
             134111760*\[Nu]*\[Chi]S^2 + \[Delta]*(9730201*\[Kappa]A + 
               86516282*\[Chi]A*\[Chi]S)) - 6489723240*e^(24/19)*ei^(88/19)*
            (-2833 + 5516*\[Nu])*(-2047975353262*\[Kappa]S*(-1 + 2*\[Nu]) + 
             2117163228955*\[Chi]A^2 - 8191901413048*\[Nu]*\[Chi]A^2 + 
             2117163228955*\[Chi]S^2 - 276751502772*\[Nu]*\[Chi]S^2 + 
             \[Delta]*(2047975353262*\[Kappa]A + 4234326457910*\[Chi]A*
                \[Chi]S)) - 63197101826183325*ei^(112/19)*
            (-1008*\[Kappa]S*(124448 - 375411*\[Nu] + 116610*\[Nu]^2) - 
             278950701*\[Chi]A^2 + 102215792*\[Delta]^2*\[Chi]A^2 + 
             1165211292*\[Nu]*\[Chi]A^2 - 235085760*\[Nu]^2*\[Chi]A^2 - 
             176734909*\[Chi]S^2 + 376286320*\[Nu]*\[Chi]S^2 - 
             68029360*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(504*\[Kappa]A*
                (-124448 + 126515*\[Nu]) + (-176734909 + 212847404*\[Nu])*
                \[Chi]A*\[Chi]S)) + 154902628608*e^(74/19)*ei^2*
            (45*\[Kappa]S*(1288935001706461 - 4177294463778654*\[Nu] + 
               7918493972374440*\[Nu]^2) - 20556501925266276*\[Chi]A^2 + 
             81704163869614286*\[Delta]^2*\[Chi]A^2 - 79734899027872416*\[Nu]*
              \[Chi]A^2 + 712664457513699600*\[Nu]^2*\[Chi]A^2 + 
             61147661944348010*\[Chi]S^2 + 46845577848720400*\[Nu]*
              \[Chi]S^2 + 155611899223682720*\[Nu]^2*\[Chi]S^2 - 
             5*\[Delta]*(9*\[Kappa]A*(-1288935001706461 + 1599424460365732*
                  \[Nu]) + 4*(-6114766194434801 + 5755766444010856*\[Nu])*
                \[Chi]A*\[Chi]S)) - 12278974941*e^(64/19)*ei^(48/19)*
            (5520*\[Kappa]S*(209905475372017 - 941420628564083*\[Nu] + 
               1512613179868518*\[Nu]^2) + 661395331763843437*\[Chi]A^2 + 
             764458173288617808*\[Delta]^2*\[Chi]A^2 - 6139672279226353948*
              \[Nu]*\[Chi]A^2 + 16699249505748438720*\[Nu]^2*\[Chi]A^2 + 
             1425853505052461245*\[Chi]S^2 + 260421436089885200*\[Nu]*
              \[Chi]S^2 + 2082778093318402480*\[Nu]^2*\[Chi]S^2 - 
             10*\[Delta]*(552*\[Kappa]A*(-209905475372017 + 521609677820049*
                  \[Nu]) + (-285170701010492249 + 323366951608109500*\[Nu])*
                \[Chi]A*\[Chi]S)) - 3588*e^(36/19)*ei^4*
            (1044*\[Kappa]S*(12865074844809000047 + 203738894213893801302*
                \[Nu] + 319207505993389781208*\[Nu]^2) - 
             269552471721024277111089*\[Chi]A^2 + 326073634913301114343504*
              \[Delta]^2*\[Chi]A^2 + 1115745822132711025804788*\[Nu]*
              \[Chi]A^2 + 666505272514197863162304*\[Nu]^2*\[Chi]A^2 + 
             56521163192276837232415*\[Chi]S^2 + 43212084697850551775000*
              \[Nu]*\[Chi]S^2 + 1270802617793555061375664*\[Nu]^2*\[Chi]S^2 + 
             2*\[Delta]*(522*\[Kappa]A*(12865074844809000047 + 
                 229469043903511801396*\[Nu]) + (56521163192276837232415 + 
                 40374009973232234567716*\[Nu])*\[Chi]A*\[Chi]S)) + 
           131105520*e^(88/19)*ei^(24/19)*(-3920527 + 15455580*\[Nu])*
            (-2948986828885*\[Kappa]S*(-1 + 2*\[Nu]) + \[Delta]*
              (2948986828885*\[Kappa]A + 6058122836776*\[Chi]A*\[Chi]S) - 
             4*((-757265354597 + 2948986828885*\[Nu])*\[Chi]A^2 + 
               (-757265354597 + 80074589503*\[Nu])*\[Chi]S^2)) + 
           4930464*e^2*(33124*ei^2*(e*ei)^(18/19)*(3866216072016148*
                \[Delta]^2*\[Chi]A^2 + \[Delta]*(7732432144032296 + 
                 203260333927841651*\[Nu])*\[Chi]A*\[Chi]S + 
               (3866216072016148 + 203260333927841651*\[Nu] - 
                 80534487366910472*\[Nu]^2)*\[Chi]S^2) + 147005*ei^(74/19)*
              (-18*\[Kappa]S*(30936092472295 - 97020967578282*\[Nu] + 
                 36669175275096*\[Nu]^2) - 2827049041226343*\[Chi]A^2 + 
               2206814975744288*\[Delta]^2*\[Chi]A^2 + 11721265129529220*
                \[Nu]*\[Chi]A^2 - 1320090309903456*\[Nu]^2*\[Chi]A^2 - 
               620234065482055*\[Chi]S^2 + 1521486691850320*\[Nu]*\[Chi]S^2 - 
               141962774214352*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*
                (9*\[Kappa]A*(-30936092472295 + 35148782633692*\[Nu]) + 
                 (-620234065482055 + 967277828237084*\[Nu])*\[Chi]A*
                  \[Chi]S)))))/(4282589248757119932306554880*Sqrt[2]*
          ei^(48/19))))
 
HlmDCMem[6, 0] = x^2*\[Epsilon]^2*
     ((5*(e^(24/19) - ei^(24/19))*(-839 + 3612*\[Nu]))/
       (1419264*Sqrt[273]*ei^(24/19)) + 
      (5*(3*ei^2*(48655927 - 184821756*\[Nu]) + 23261*e^2*
          (-839 + 3612*\[Nu]) + 722*e^(14/19)*ei^(24/19)*
          (-175141 + 651588*\[Nu])))/(14345920512*Sqrt[273]*(ei/e)^(24/19)) + 
      (5*(49*ei^4*(16680917087 - 53534704476*\[Nu]) + 1815460283*e^4*
          (-839 + 3612*\[Nu]) + 1303210*e^(52/19)*ei^(24/19)*
          (-12361241 + 45993108*\[Nu]) - 345592*e^2*ei^2*
          (-48655927 + 184821756*\[Nu])))/(718137652936704*Sqrt[273]*
        (ei/e)^(24/19)) + (5*(64141042823403*e^6*(-839 + 3612*\[Nu]) - 
         14004979326*e^4*ei^2*(-48655927 + 184821756*\[Nu]) + 
         94091762*e^(90/19)*ei^(24/19)*(-7092285107 + 26389323996*\[Nu]) - 
         2930886*e^2*ei^4*(-16680917087 + 53534704476*\[Nu]) + 
         91*ei^6*(-100826162933963 + 336632397661644*\[Nu])))/
       (18665833875130810368*Sqrt[273]*(ei/e)^(24/19))) + 
    x^3*\[Epsilon]^4*((-90*e^(24/19)*ei^(12/19)*(2376887 - 14860720*\[Nu] + 
          19923792*\[Nu]^2) - 19*ei^(36/19)*(45661561 - 255563280*\[Nu] + 
          363772080*\[Nu]^2) + e^(36/19)*(1081489489 - 6193167120*\[Nu] + 
          8704810800*\[Nu]^2))/(40772616192*Sqrt[273]*ei^(36/19)) + 
      (e^(24/19)*(302393*e^(50/19)*(1081489489 - 6193167120*\[Nu] + 
           8704810800*\[Nu]^2) + 1426672*e^(14/19)*ei^(36/19)*
          (13522442272 - 77062595925*\[Nu] + 115417999620*\[Nu]^2) - 
         735*e^2*ei^(12/19)*(159041293899 - 913298900408*\[Nu] + 
           984178414128*\[Nu]^2) + 2340*ei^(50/19)*(137842241191 - 
           791986128080*\[Nu] + 1019476806096*\[Nu]^2) - 
         168*e^(12/19)*ei^2*(118004643485102 - 672605899323435*\[Nu] + 
           1005703016654700*\[Nu]^2)))/(3571789905395712*Sqrt[273]*
        ei^(36/19)) + (e^(24/19)*(10651549610*e^(88/19)*
          (1081489489 - 6193167120*\[Nu] + 8704810800*\[Nu]^2) + 
         12600*ei^(88/19)*(47257038107471 - 243675756432400*\[Nu] + 
           295297429889616*\[Nu]^2) - 9570240*e^(50/19)*ei^2*
          (118004643485102 - 672605899323435*\[Nu] + 1005703016654700*
            \[Nu]^2) + 4952198*e^(52/19)*ei^(36/19)*(225971503921877 - 
           1288535252520960*\[Nu] + 1930944652565040*\[Nu]^2) - 
         1125*e^4*ei^(12/19)*(4278680271864451 - 24197109215592152*\[Nu] + 
           24870084486107952*\[Nu]^2) + 3600*e^2*ei^(50/19)*
          (9223243844976507 - 48292418443744984*\[Nu] + 50359242169554864*
            \[Nu]^2) - 3*e^(12/19)*ei^4*(10076003557682657627 - 
           53249972801790776760*\[Nu] + 62451968965269614640*\[Nu]^2)))/
       (58944738553044664320*Sqrt[273]*ei^(36/19)) + 
      (e^(24/19)*(246984666220282900*e^(126/19)*(1081489489 - 
           6193167120*\[Nu] + 8704810800*\[Nu]^2) - 273787431175440*e^(88/19)*
          ei^2*(118004643485102 - 672605899323435*\[Nu] + 1005703016654700*
            \[Nu]^2) - 12066600*ei^(126/19)*(285640519591917179 - 
           1509836697319177360*\[Nu] + 1856864305501628304*\[Nu]^2) + 
         30391639126*e^(90/19)*ei^(36/19)*(1065585375488561129 - 
           6073335158886574320*\[Nu] + 9089722103689871280*\[Nu]^2) + 
         15743700*e^2*ei^(88/19)*(3162043667387864067 - 14693194487619523544*
            \[Nu] + 14586849543743306544*\[Nu]^2) - 138798387*e^(50/19)*ei^4*
          (10076003557682657627 - 53249972801790776760*\[Nu] + 
           62451968965269614640*\[Nu]^2) + 4475250*e^4*ei^(50/19)*
          (248132485058613685043 - 1277558337046140175096*\[Nu] + 
           1272572726077195264176*\[Nu]^2) - 20475*e^6*ei^(12/19)*
          (5879476988625155849615 - 33032453934069423808184*\[Nu] + 
           33238018725944961622512*\[Nu]^2) + 7350*e^(12/19)*ei^6*
          (2491914028617678683843 - 24444060361294631556888*\[Nu] + 
           53577100465496674243632*\[Nu]^2)))/(790048591332076669171138560*
        Sqrt[273]*ei^(36/19))) + SO*x^(7/2)*\[Epsilon]^5*
     ((5*(14*e^(24/19)*ei^(18/19)*(-839 + 3612*\[Nu])*
          (-157*\[Delta]*\[Chi]A + (-157 + 110*\[Nu])*\[Chi]S) - 
         57*ei^(42/19)*(\[Delta]*(-54403 + 446460*\[Nu])*\[Chi]A + 
           (-54403 + 66704*\[Nu] + 785904*\[Nu]^2)*\[Chi]S) + 
         e^(42/19)*(\[Delta]*(-4945093 + 33387396*\[Nu])*\[Chi]A + 
           (-4945093 + 13033364*\[Nu] + 39234048*\[Nu]^2)*\[Chi]S)))/
       (1698859008*Sqrt[273]*ei^(42/19)) + 
      (e^(24/19)*(240*ei^(56/19)*(-48655927 + 184821756*\[Nu])*
          (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) + 
         7*e^2*ei^(18/19)*(-839 + 3612*\[Nu])*(-113175949*\[Delta]*\[Chi]A + 
           (-113175949 + 41301776*\[Nu])*\[Chi]S) + 232610*e^(56/19)*
          (\[Delta]*(-4945093 + 33387396*\[Nu])*\[Chi]A + 
           (-4945093 + 13033364*\[Nu] + 39234048*\[Nu]^2)*\[Chi]S) - 
         823080*e^(14/19)*ei^(42/19)*(2*\[Delta]*(-16976455 + 
             116358648*\[Nu])*\[Chi]A + (-33952910 + 12168211*\[Nu] + 
             444154116*\[Nu]^2)*\[Chi]S) + 21*e^(18/19)*ei^2*
          (\[Delta]*(-1220333814247 + 8556007939836*\[Nu])*\[Chi]A + 
           (-1220333814247 + 87582195044*\[Nu] + 17156340864096*\[Nu]^2)*
            \[Chi]S)))/(19625219260416*Sqrt[273]*ei^(42/19)) + 
      (e^(24/19)*(25636*e^2*ei^(56/19)*(-48655927 + 184821756*\[Nu])*
          (113175949*\[Delta]*\[Chi]A + (113175949 - 41301776*\[Nu])*
            \[Chi]S) + 966280*ei^(94/19)*(-16680917087 + 53534704476*\[Nu])*
          (157*\[Delta]*\[Chi]A + (157 - 110*\[Nu])*\[Chi]S) + 
         143871*e^4*ei^(18/19)*(-839 + 3612*\[Nu])*
          (-189611322769*\[Delta]*\[Chi]A + (-189611322769 + 
             49774532976*\[Nu])*\[Chi]S) + 6952069627230*e^(94/19)*
          (\[Delta]*(-4945093 + 33387396*\[Nu])*\[Chi]A + 
           (-4945093 + 13033364*\[Nu] + 39234048*\[Nu]^2)*\[Chi]S) - 
         4308412260*e^(52/19)*ei^(42/19)*(9*\[Delta]*(-36076614005 + 
             241522072148*\[Nu])*\[Chi]A + (-324689526045 + 
             177345933110*\[Nu] + 4018876732248*\[Nu]^2)*\[Chi]S) + 
         1043558243*e^(56/19)*ei^2*(\[Delta]*(-1220333814247 + 
             8556007939836*\[Nu])*\[Chi]A + (-1220333814247 + 
             87582195044*\[Nu] + 17156340864096*\[Nu]^2)*\[Chi]S) - 
         182*e^(18/19)*ei^4*(\[Delta]*(-163648399752181339 + 
             1326592427149021452*\[Nu])*\[Chi]A + (-163648399752181339 - 
             488821146525940972*\[Nu] + 3768709409950686432*\[Nu]^2)*
            \[Chi]S)))/(242164634222091829248*Sqrt[273]*ei^(42/19)) + 
      (e^(24/19)*(17757792*e^4*ei^(56/19)*(-48655927 + 184821756*\[Nu])*
          (189611322769*\[Delta]*\[Chi]A + (189611322769 - 49774532976*\[Nu])*
            \[Chi]S) + 3478608*e^2*ei^(94/19)*(-16680917087 + 
           53534704476*\[Nu])*(113175949*\[Delta]*\[Chi]A + 
           (113175949 - 41301776*\[Nu])*\[Chi]S) + 28712320*ei^(132/19)*
          (-100826162933963 + 336632397661644*\[Nu])*(-157*\[Delta]*\[Chi]A + 
           (-157 + 110*\[Nu])*\[Chi]S) + 4641*e^6*ei^(18/19)*
          (-839 + 3612*\[Nu])*(-5038160345051358953*\[Delta]*\[Chi]A + 
           (-5038160345051358953 + 1065846237840519832*\[Nu])*\[Chi]S) + 
         5544270678318131760*e^(132/19)*(\[Delta]*(-4945093 + 33387396*\[Nu])*
            \[Chi]A + (-4945093 + 13033364*\[Nu] + 39234048*\[Nu]^2)*
            \[Chi]S) + 1051152927637176*e^(94/19)*ei^2*
          (\[Delta]*(-1220333814247 + 8556007939836*\[Nu])*\[Chi]A + 
           (-1220333814247 + 87582195044*\[Nu] + 17156340864096*\[Nu]^2)*
            \[Chi]S) - 1036891217240*e^(90/19)*ei^(42/19)*
          (5*\[Delta]*(-272121296956445 + 1802900172157476*\[Nu])*\[Chi]A + 
           (-1360606484782225 + 864467120228516*\[Nu] + 16407914694007872*
              \[Nu]^2)*\[Chi]S) - 304812144*e^(56/19)*ei^4*
          (\[Delta]*(-163648399752181339 + 1326592427149021452*\[Nu])*
            \[Chi]A + (-163648399752181339 - 488821146525940972*\[Nu] + 
             3768709409950686432*\[Nu]^2)*\[Chi]S) + 1547*e^(18/19)*ei^6*
          (\[Delta]*(-119662226893163888573 + 9841403996830518016884*\[Nu])*
            \[Chi]A + (-119662226893163888573 - 9896134462812045205844*
              \[Nu] + 40485734239865689499424*\[Nu]^2)*\[Chi]S)))/
       (100709490763209773213024256*Sqrt[273]*ei^(42/19))) + 
    SO^2*x^4*\[Epsilon]^6*((5*(-e^(24/19) + ei^(24/19))*
        (e^(24/19)*(24*\[Delta]*\[Kappa]A*(-46297 + 179508*\[Nu]) - 
           24*\[Kappa]S*(46297 - 272102*\[Nu] + 397320*\[Nu]^2) - 
           1170697*\[Chi]A^2 + 9009156*\[Nu]*\[Chi]A^2 - 19071360*\[Nu]^2*
            \[Chi]A^2 + 2*\[Delta]*(-1170697 + 4564644*\[Nu])*\[Chi]A*
            \[Chi]S - 1170697*\[Chi]S^2 + 4802920*\[Nu]*\[Chi]S^2 - 
           2864400*\[Nu]^2*\[Chi]S^2) + 19*ei^(24/19)*
          (32*\[Delta]*\[Kappa]A*(-845 + 2856*\[Nu]) - 32*\[Kappa]S*
            (845 - 4546*\[Nu] + 7224*\[Nu]^2) - 27879*\[Chi]A^2 + 
           203164*\[Nu]*\[Chi]A^2 - 462336*\[Nu]^2*\[Chi]A^2 + 
           6*\[Delta]*(-9293 + 31668*\[Nu])*\[Chi]A*\[Chi]S - 
           27879*\[Chi]S^2 + 98360*\[Nu]*\[Chi]S^2 - 111216*\[Nu]^2*
            \[Chi]S^2)))/(862912512*Sqrt[273]*ei^(48/19)) + 
      (e^(24/19)*(-30*ei^(62/19)*(-48655927 + 184821756*\[Nu])*
          (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
           191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
            \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) + 
         54872*e^(14/19)*ei^(48/19)*(\[Delta]*\[Kappa]A*(-20323301 + 
             54435108*\[Nu]) + \[Kappa]S*(-20323301 + 95081710*\[Nu] - 
             140283528*\[Nu]^2) - 22066269*\[Chi]A^2 + 142387076*\[Nu]*
            \[Chi]A^2 - 280567056*\[Nu]^2*\[Chi]A^2 + 6*\[Delta]*
            (-7355423 + 20364624*\[Nu])*\[Chi]A*\[Chi]S - 
           22066269*\[Chi]S^2 + 68065744*\[Nu]*\[Chi]S^2 - 
           89461680*\[Nu]^2*\[Chi]S^2) + 116305*e^(62/19)*
          (-24*\[Delta]*\[Kappa]A*(-46297 + 179508*\[Nu]) + 
           24*\[Kappa]S*(46297 - 272102*\[Nu] + 397320*\[Nu]^2) + 
           1170697*\[Chi]A^2 - 9009156*\[Nu]*\[Chi]A^2 + 19071360*\[Nu]^2*
            \[Chi]A^2 + 2*\[Delta]*(1170697 - 4564644*\[Nu])*\[Chi]A*
            \[Chi]S + 1170697*\[Chi]S^2 - 4802920*\[Nu]*\[Chi]S^2 + 
           2864400*\[Nu]^2*\[Chi]S^2) + 12*e^(24/19)*ei^2*
          (\[Delta]*\[Kappa]A*(69227038931 - 162437718828*\[Nu]) + 
           \[Kappa]S*(69227038931 - 300891796690*\[Nu] + 459608148888*
              \[Nu]^2) + 75377741539*\[Chi]A^2 - 462762892736*\[Nu]*
            \[Chi]A^2 + 919216297776*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(75377741539 - 185854737012*\[Nu])*\[Chi]A*\[Chi]S + 
           75377741539*\[Chi]S^2 - 210457547444*\[Nu]*\[Chi]S^2 + 
           363133495200*\[Nu]^2*\[Chi]S^2) + 5*e^2*ei^(24/19)*
          (-839 + 3612*\[Nu])*(-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 
           25904401*\[Chi]A^2 - 99734624*\[Nu]*\[Chi]A^2 + 
           25904401*\[Chi]S^2 - 3882980*\[Nu]*\[Chi]S^2 + 
           \[Delta]*(24933656*\[Kappa]A + 51808802*\[Chi]A*\[Chi]S))))/
       (4361159835648*Sqrt[273]*ei^(48/19)) + 
      (e^(24/19)*(-10290*ei^(100/19)*(-16680917087 + 53534704476*\[Nu])*
          (178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 356*\[Kappa]S*\[Nu] + 
           191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 382*\[Delta]*\[Chi]A*
            \[Chi]S + 191*\[Chi]S^2 - 52*\[Nu]*\[Chi]S^2) + 
         15*e^(24/19)*ei^4*(32*\[Delta]*\[Kappa]A*(-400388252426881 + 
             802766187587868*\[Nu]) - 32*\[Kappa]S*(400388252426881 - 
             1603542692441630*\[Nu] + 3159441545844864*\[Nu]^2) - 
           13810688311230753*\[Chi]A^2 + 80248459845737732*\[Nu]*\[Chi]A^2 - 
           202204258934071296*\[Nu]^2*\[Chi]A^2 + 6*\[Delta]*
            (-4603562770410251 + 9666254511698988*\[Nu])*\[Chi]A*\[Chi]S - 
           13810688311230753*\[Chi]S^2 + 32991820469379208*\[Nu]*\[Chi]S^2 - 
           112691169051964944*\[Nu]^2*\[Chi]S^2) + 4952198*e^(52/19)*
          ei^(48/19)*(16*\[Delta]*\[Kappa]A*(-72845532337 + 
             205322523756*\[Nu]) - 16*\[Kappa]S*(72845532337 - 
             351013588430*\[Nu] + 503800610376*\[Nu]^2) - 
           1252030791783*\[Chi]A^2 + 8278732353532*\[Nu]*\[Chi]A^2 - 
           16121619532032*\[Nu]^2*\[Chi]A^2 + 6*\[Delta]*(-417343597261 + 
             1205539427988*\[Nu])*\[Chi]A*\[Chi]S - 1252030791783*\[Chi]S^2 + 
           3962627381528*\[Nu]*\[Chi]S^2 - 4306809627120*\[Nu]^2*\[Chi]S^2) + 
         331302601175*e^(100/19)*(-24*\[Delta]*\[Kappa]A*
            (-46297 + 179508*\[Nu]) + 24*\[Kappa]S*(46297 - 272102*\[Nu] + 
             397320*\[Nu]^2) + 1170697*\[Chi]A^2 - 9009156*\[Nu]*\[Chi]A^2 + 
           19071360*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(1170697 - 4564644*\[Nu])*
            \[Chi]A*\[Chi]S + 1170697*\[Chi]S^2 - 4802920*\[Nu]*\[Chi]S^2 + 
           2864400*\[Nu]^2*\[Chi]S^2) + 58059456*e^(62/19)*ei^2*
          (\[Delta]*\[Kappa]A*(69227038931 - 162437718828*\[Nu]) + 
           \[Kappa]S*(69227038931 - 300891796690*\[Nu] + 459608148888*
              \[Nu]^2) + 75377741539*\[Chi]A^2 - 462762892736*\[Nu]*
            \[Chi]A^2 + 919216297776*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(75377741539 - 185854737012*\[Nu])*\[Chi]A*\[Chi]S + 
           75377741539*\[Chi]S^2 - 210457547444*\[Nu]*\[Chi]S^2 + 
           363133495200*\[Nu]^2*\[Chi]S^2) - 1560*e^2*ei^(62/19)*
          (-48655927 + 184821756*\[Nu])*(-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 
           25904401*\[Chi]A^2 - 99734624*\[Nu]*\[Chi]A^2 + 
           25904401*\[Chi]S^2 - 3882980*\[Nu]*\[Chi]S^2 + 
           \[Delta]*(24933656*\[Kappa]A + 51808802*\[Chi]A*\[Chi]S)) + 
         140*e^4*ei^(24/19)*(-839 + 3612*\[Nu])*(-2948986828885*\[Kappa]S*
            (-1 + 2*\[Nu]) + \[Delta]*(2948986828885*\[Kappa]A + 
             6058122836776*\[Chi]A*\[Chi]S) - 
           4*((-757265354597 + 2948986828885*\[Nu])*\[Chi]A^2 + 
             (-757265354597 + 80074589503*\[Nu])*\[Chi]S^2))))/
       (4584590776347918336*Sqrt[273]*ei^(48/19)) + 
      (e^(24/19)*(70070*ei^(138/19)*(-100826162933963 + 336632397661644*
            \[Nu])*(178*\[Delta]*\[Kappa]A + 178*\[Kappa]S - 
           356*\[Kappa]S*\[Nu] + 191*\[Chi]A^2 - 712*\[Nu]*\[Chi]A^2 + 
           382*\[Delta]*\[Chi]A*\[Chi]S + 191*\[Chi]S^2 - 
           52*\[Nu]*\[Chi]S^2) + 6579540*e^(62/19)*ei^4*
          (32*\[Delta]*\[Kappa]A*(-400388252426881 + 802766187587868*\[Nu]) - 
           32*\[Kappa]S*(400388252426881 - 1603542692441630*\[Nu] + 
             3159441545844864*\[Nu]^2) - 13810688311230753*\[Chi]A^2 + 
           80248459845737732*\[Nu]*\[Chi]A^2 - 202204258934071296*\[Nu]^2*
            \[Chi]A^2 + 6*\[Delta]*(-4603562770410251 + 9666254511698988*
              \[Nu])*\[Chi]A*\[Chi]S - 13810688311230753*\[Chi]S^2 + 
           32991820469379208*\[Nu]*\[Chi]S^2 - 112691169051964944*\[Nu]^2*
            \[Chi]S^2) + 66501428913444525*e^(138/19)*
          (-24*\[Delta]*\[Kappa]A*(-46297 + 179508*\[Nu]) + 
           24*\[Kappa]S*(46297 - 272102*\[Nu] + 397320*\[Nu]^2) + 
           1170697*\[Chi]A^2 - 9009156*\[Nu]*\[Chi]A^2 + 19071360*\[Nu]^2*
            \[Chi]A^2 + 2*\[Delta]*(1170697 - 4564644*\[Nu])*\[Chi]A*
            \[Chi]S + 1170697*\[Chi]S^2 - 4802920*\[Nu]*\[Chi]S^2 + 
           2864400*\[Nu]^2*\[Chi]S^2) + 14993809150320*e^(100/19)*ei^2*
          (\[Delta]*\[Kappa]A*(69227038931 - 162437718828*\[Nu]) + 
           \[Kappa]S*(69227038931 - 300891796690*\[Nu] + 459608148888*
              \[Nu]^2) + 75377741539*\[Chi]A^2 - 462762892736*\[Nu]*
            \[Chi]A^2 + 919216297776*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(75377741539 - 185854737012*\[Nu])*\[Chi]A*\[Chi]S + 
           75377741539*\[Chi]S^2 - 210457547444*\[Nu]*\[Chi]S^2 + 
           363133495200*\[Nu]^2*\[Chi]S^2) - 35754869560*e^(90/19)*ei^(48/19)*
          (\[Delta]*\[Kappa]A*(42913583736551 - 123729947141388*\[Nu]) + 
           \[Kappa]S*(42913583736551 - 209557114614490*\[Nu] + 
             296988447329112*\[Nu]^2) + 45845832221419*\[Chi]A^2 - 
           306652095373076*\[Nu]*\[Chi]A^2 + 593976894658224*\[Nu]^2*
            \[Chi]A^2 + 2*\[Delta]*(45845832221419 - 134997760426872*\[Nu])*
            \[Chi]A*\[Chi]S + 45845832221419*\[Chi]S^2 - 146726754366344*
            \[Nu]*\[Chi]S^2 + 144128359234608*\[Nu]^2*\[Chi]S^2) + 
         98*e^(24/19)*ei^6*(6*\[Delta]*\[Kappa]A*(-8541870653698479251 + 
             40518251799001665708*\[Nu]) - 6*\[Kappa]S*(8541870653698479251 - 
             57601993106398624210*\[Nu] + 52526063141904828696*\[Nu]^2) - 
           44288597791391875171*\[Chi]A^2 + 422929023658235362692*\[Nu]*
            \[Chi]A^2 - 630312757702857944352*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(-44288597791391875171 + 217924127969471860668*\[Nu])*
            \[Chi]A*\[Chi]S - 44288597791391875171*\[Chi]S^2 + 
           190073623446275859328*\[Nu]*\[Chi]S^2 + 442866816771334566960*
            \[Nu]^2*\[Chi]S^2) - 48510*e^2*ei^(100/19)*(-16680917087 + 
           53534704476*\[Nu])*(-24933656*\[Kappa]S*(-1 + 2*\[Nu]) + 
           25904401*\[Chi]A^2 - 99734624*\[Nu]*\[Chi]A^2 + 
           25904401*\[Chi]S^2 - 3882980*\[Nu]*\[Chi]S^2 + 
           \[Delta]*(24933656*\[Kappa]A + 51808802*\[Chi]A*\[Chi]S)) + 
         77*e^6*ei^(24/19)*(-839 + 3612*\[Nu])*(-1165556009138332944*
            \[Kappa]S*(-1 + 2*\[Nu]) + 1189984978416247669*\[Chi]A^2 - 
           4662224036553331776*\[Nu]*\[Chi]A^2 + 1189984978416247669*
            \[Chi]S^2 - 97715877111658900*\[Nu]*\[Chi]S^2 + 
           2*\[Delta]*(582778004569166472*\[Kappa]A + 1189984978416247669*
              \[Chi]A*\[Chi]S)) - 3960*e^4*ei^(62/19)*(-48655927 + 
           184821756*\[Nu])*(-2948986828885*\[Kappa]S*(-1 + 2*\[Nu]) + 
           \[Delta]*(2948986828885*\[Kappa]A + 6058122836776*\[Chi]A*
              \[Chi]S) - 4*((-757265354597 + 2948986828885*\[Nu])*\[Chi]A^2 + 
             (-757265354597 + 80074589503*\[Nu])*\[Chi]S^2))))/
       (436929839349062009094144*Sqrt[273]*ei^(48/19)))
 
HlmDCMem[8, 0] = x^3*\[Epsilon]^4*
    (((e^(36/19) - ei^(36/19))*(75601 - 452070*\[Nu] + 733320*\[Nu]^2))/
      (213497856*Sqrt[119]*ei^(36/19)) + 
     (3323*e^2*(75601 - 452070*\[Nu] + 733320*\[Nu]^2) + 
       722*e^(2/19)*ei^(36/19)*(31645823 - 192093570*\[Nu] + 
         293380920*\[Nu]^2) - 3*ei^2*(7699835443 - 46731262050*\[Nu] + 
         71419282200*\[Nu]^2))/(205527269376*Sqrt[119]*(ei/e)^(36/19)) + 
     (409674985*e^4*(75601 - 452070*\[Nu] + 733320*\[Nu]^2) - 
       598140*e^2*ei^2*(7699835443 - 46731262050*\[Nu] + 
         71419282200*\[Nu]^2) + 130321*e^(40/19)*ei^(36/19)*
        (36088208969 - 218933563950*\[Nu] + 333741227400*\[Nu]^2) - 
       6*ei^4*(21407291285669 - 127500945901650*\[Nu] + 
         179197316814600*\[Nu]^2))/(11871255079157760*Sqrt[119]*
       (ei/e)^(36/19)) + (79826976800350*e^6*(75601 - 452070*\[Nu] + 
         733320*\[Nu]^2) - 143795919735*e^4*ei^2*(7699835443 - 
         46731262050*\[Nu] + 71419282200*\[Nu]^2) - 2332746*e^2*ei^4*
        (21407291285669 - 127500945901650*\[Nu] + 179197316814600*\[Nu]^2) + 
       47045881*e^(78/19)*ei^(36/19)*(24432057125909 - 148165211424150*
          \[Nu] + 225585603066600*\[Nu]^2) + 
       1500*ei^6*(1120027365343529 - 7027894668537558*\[Nu] + 
         10940692640515848*\[Nu]^2))/(1337083202075696824320*Sqrt[119]*
       (ei/e)^(36/19)))


(* ::Subsection::Closed:: *)
(*DC memory for circular orbits*)


HlmDCMemCirc[2, 0] = (-5*x)/(14*Sqrt[6]) - 
    (5*x^2*\[Epsilon]^2*(-4075 + 5628*\[Nu]))/(56448*Sqrt[6]) + 
    (x^3*\[Epsilon]^4*(759386065 + 936041400*\[Nu] - 195274800*\[Nu]^2))/
     (938843136*Sqrt[6]) + 
    SO*(x^4*\[Epsilon]^6*((115*Pi*\[Delta]*\[Chi]A)/(5376*Sqrt[6]) - 
        (115*Pi*(-1 + 4*\[Nu])*\[Chi]S)/(5376*Sqrt[6])) + 
      x^(5/2)*\[Epsilon]^3*((-1303*\[Delta]*\[Chi]A)/(1344*Sqrt[6]) + 
        ((-1303 + 92*\[Nu])*\[Chi]S)/(1344*Sqrt[6])) + 
      x^(7/2)*\[Epsilon]^5*((-5*\[Delta]*(4888427 + 410172*\[Nu])*\[Chi]A)/
         (9483264*Sqrt[6]) - (5*(4888427 - 7631920*\[Nu] + 74256*\[Nu]^2)*
          \[Chi]S)/(9483264*Sqrt[6]))) + 
    SO^2*(x^3*\[Epsilon]^4*((5*\[Delta]*\[Kappa]A)/(14*Sqrt[6]) - 
        (5*\[Kappa]S*(-1 + 2*\[Nu]))/(14*Sqrt[6]) + 
        (5*(33 - 128*\[Nu])*\[Chi]A^2)/(448*Sqrt[6]) + 
        (55*Sqrt[3/2]*\[Delta]*\[Chi]A*\[Chi]S)/224 - 
        (5*(-33 + 4*\[Nu])*\[Chi]S^2)/(448*Sqrt[6])) + 
      x^4*\[Epsilon]^6*((-5*\[Delta]*\[Kappa]A*(-44803 + 25242*\[Nu]))/
         (112896*Sqrt[6]) - (5*\[Kappa]S*(-44803 + 114848*\[Nu] + 
           10416*\[Nu]^2))/(112896*Sqrt[6]) - 
        (5*(-203749 + 865128*\[Nu] + 166656*\[Nu]^2)*\[Chi]A^2)/
         (903168*Sqrt[6]) - (5*\[Delta]*(-29107 + 79316*\[Nu])*\[Chi]A*
          \[Chi]S)/(64512*Sqrt[6]) + (5*(203749 - 1060292*\[Nu] + 
           382144*\[Nu]^2)*\[Chi]S^2)/(903168*Sqrt[6])))
 
HlmDCMemCirc[4, 0] = -x/(504*Sqrt[2]) - (19*x^2*\[Epsilon]^2*(-9479 + 40124*\[Nu]))/
     (14902272*Sqrt[2]) - (x^3*\[Epsilon]^4*(24215523937 - 
       140432459328*\[Nu] + 53657768304*\[Nu]^2))/(878757175296*Sqrt[2]) + 
    SO*(x^4*\[Epsilon]^6*((25*Pi*\[Delta]*\[Chi]A)/(48384*Sqrt[2]) - 
        (25*Pi*(-1 + 4*\[Nu])*\[Chi]S)/(48384*Sqrt[2])) + 
      x^(5/2)*\[Epsilon]^3*((-23*\[Delta]*\[Chi]A)/(4032*Sqrt[2]) - 
        ((23 + 68*\[Nu])*\[Chi]S)/(4032*Sqrt[2])) + x^(7/2)*\[Epsilon]^5*
       ((\[Delta]*(1319709 - 12539422*\[Nu])*\[Chi]A)/(78236928*Sqrt[2]) - 
        ((-1319709 + 6548872*\[Nu] + 5249776*\[Nu]^2)*\[Chi]S)/
         (78236928*Sqrt[2]))) + 
    SO^2*(x^3*\[Epsilon]^4*((\[Delta]*\[Kappa]A)/(504*Sqrt[2]) + 
        (\[Kappa]S*(1 - 2*\[Nu]))/(504*Sqrt[2]) + 
        ((53 - 192*\[Nu])*\[Chi]A^2)/(24192*Sqrt[2]) + 
        (53*\[Delta]*\[Chi]A*\[Chi]S)/(12096*Sqrt[2]) + 
        ((53 - 20*\[Nu])*\[Chi]S^2)/(24192*Sqrt[2])) + 
      x^4*\[Epsilon]^6*((5*\[Delta]*\[Kappa]A*(1979 + 89740*\[Nu]))/
         (14902272*Sqrt[2]) - (5*\[Kappa]S*(-1979 - 85782*\[Nu] + 
           296072*\[Nu]^2))/(14902272*Sqrt[2]) - 
        (5*(1582007 - 19381052*\[Nu] + 56845824*\[Nu]^2)*\[Chi]A^2)/
         (1430618112*Sqrt[2]) - (5*\[Delta]*(226001 + 180188*\[Nu])*\[Chi]A*
          \[Chi]S)/(102187008*Sqrt[2]) + 
        (5*(-1582007 - 15575656*\[Nu] + 3065104*\[Nu]^2)*\[Chi]S^2)/
         (1430618112*Sqrt[2])))
 
HlmDCMemCirc[6, 0] = (-5*x^2*\[Epsilon]^2*(-839 + 3612*\[Nu]))/
     (1419264*Sqrt[273]) - (x^3*\[Epsilon]^4*(45661561 - 255563280*\[Nu] + 
       363772080*\[Nu]^2))/(2145927168*Sqrt[273]) + 
    SO*x^(7/2)*\[Epsilon]^5*((-5*\[Delta]*(-54403 + 446460*\[Nu])*\[Chi]A)/
       (29804544*Sqrt[273]) - (5*(-54403 + 66704*\[Nu] + 785904*\[Nu]^2)*
        \[Chi]S)/(29804544*Sqrt[273])) + SO^2*x^4*\[Epsilon]^6*
     ((5*\[Delta]*\[Kappa]A*(-845 + 2856*\[Nu]))/(1419264*Sqrt[273]) - 
      (5*\[Kappa]S*(845 - 4546*\[Nu] + 7224*\[Nu]^2))/(1419264*Sqrt[273]) - 
      (5*(27879 - 203164*\[Nu] + 462336*\[Nu]^2)*\[Chi]A^2)/
       (45416448*Sqrt[273]) + (5*\[Delta]*(-9293 + 31668*\[Nu])*\[Chi]A*
        \[Chi]S)/(7569408*Sqrt[273]) - 
      (5*(27879 - 98360*\[Nu] + 111216*\[Nu]^2)*\[Chi]S^2)/
       (45416448*Sqrt[273]))
 
HlmDCMemCirc[8, 0] = -(x^3*\[Epsilon]^4*(75601 - 452070*\[Nu] + 733320*\[Nu]^2))/
    (213497856*Sqrt[119])


(* ::Section::Closed:: *)
(*Phase redefinition for the full waveform*)


(* ::Text:: *)
(*This section includes the full waveform modes, expanded in eccentricity to O(e^6), and expressed in terms of (x,e,\[Xi],\[Psi]), where \[Psi] and \[Xi] are chosen to absorb the arbitrary gauge constant in the tail part, as explained in Sec. IV.D of the paper. Therefore, we only include the modes that have tail contributions in the spin part at 3PN and the nonspinning part at 2PN. We provide the oscillatory memory contribution separately for all modes, to facilitate comparison with the literature. The DC memory is the same as above, since it is not affected by the phase redefinition.*)
(**)
(*h[l, m] = 8 M \[Nu] / R x Sqrt[\[Pi]/5] E^(-i m \[Psi]) H\[Psi][l, m]*)
(*H\[Psi][l, m] = H\[Psi]InstTail[l, m] + H\[Psi]OscMem[l, m]*)


(* ::Subsection::Closed:: *)
(*relations between \[Psi], \[Xi], and \[ScriptL]*)


(* ::Subsubsection::Closed:: *)
(*\[Psi] in terms of \[Xi]*)


\[Psi]\[LetterSpace]\[Xi] = e*(I/E^(I*\[Xi]) - I*E^(I*\[Xi])) + 
    e^2*(((5*I)/8)/E^((2*I)*\[Xi]) - ((5*I)/8)*E^((2*I)*\[Xi])) + 
    e^3*((-I/8)/E^(I*\[Xi]) + (I/8)*E^(I*\[Xi]) + 
      ((13*I)/24)/E^((3*I)*\[Xi]) - ((13*I)/24)*E^((3*I)*\[Xi])) + 
    e^4*(((-11*I)/48)/E^((2*I)*\[Xi]) + ((11*I)/48)*E^((2*I)*\[Xi]) + 
      ((103*I)/192)/E^((4*I)*\[Xi]) - ((103*I)/192)*E^((4*I)*\[Xi])) + 
    e^5*(((5*I)/192)/E^(I*\[Xi]) - ((5*I)/192)*E^(I*\[Xi]) - 
      ((43*I)/128)/E^((3*I)*\[Xi]) + ((43*I)/128)*E^((3*I)*\[Xi]) + 
      ((1097*I)/1920)/E^((5*I)*\[Xi]) - ((1097*I)/1920)*E^((5*I)*\[Xi])) + 
    e^6*(((17*I)/384)/E^((2*I)*\[Xi]) - ((17*I)/384)*E^((2*I)*\[Xi]) - 
      ((451*I)/960)/E^((4*I)*\[Xi]) + ((451*I)/960)*E^((4*I)*\[Xi]) + 
      ((1223*I)/1920)/E^((6*I)*\[Xi]) - ((1223*I)/1920)*E^((6*I)*\[Xi])) + 
    \[Xi] + x*\[Epsilon]^2*(e*(((-I/2)*(-10 + \[Nu]))/E^(I*\[Xi]) + 
        (I/2)*E^(I*\[Xi])*(-10 + \[Nu])) + 
      e^3*(((-I/16)*(-46 + \[Nu]))/E^(I*\[Xi]) + (I/16)*E^(I*\[Xi])*
         (-46 + \[Nu]) - ((I/16)*(-62 + 9*\[Nu]))/E^((3*I)*\[Xi]) + 
        (I/16)*E^((3*I)*\[Xi])*(-62 + 9*\[Nu])) + 
      e^5*(((-I/384)*(-1202 + 41*\[Nu]))/E^(I*\[Xi]) + (I/384)*E^(I*\[Xi])*
         (-1202 + 41*\[Nu]) + ((I/256)*(-94 + 63*\[Nu]))/E^((3*I)*\[Xi]) - 
        (I/256)*E^((3*I)*\[Xi])*(-94 + 63*\[Nu]) - 
        ((I/3840)*(-19082 + 3125*\[Nu]))/E^((5*I)*\[Xi]) + 
        (I/3840)*E^((5*I)*\[Xi])*(-19082 + 3125*\[Nu])) + 3*\[Xi] + 
      e^2*(((-I/8)*(-31 + 4*\[Nu]))/E^((2*I)*\[Xi]) + (I/8)*E^((2*I)*\[Xi])*
         (-31 + 4*\[Nu]) + 3*\[Xi]) + 
      e^4*(((821*I)/192 - ((2*I)/3)*\[Nu])/E^((4*I)*\[Xi]) + 
        ((I/48)*(41 + 4*\[Nu]))/E^((2*I)*\[Xi]) - (I/48)*E^((2*I)*\[Xi])*
         (41 + 4*\[Nu]) + (I/192)*E^((4*I)*\[Xi])*(-821 + 128*\[Nu]) + 
        3*\[Xi]) + e^6*(((I/640)*(3815 - 648*\[Nu]))/E^((6*I)*\[Xi]) + 
        (I/15)*E^((4*I)*\[Xi])*(25 - 7*\[Nu]) + ((I/15)*(-25 + 7*\[Nu]))/
         E^((4*I)*\[Xi]) - ((I/384)*(-635 + 32*\[Nu]))/E^((2*I)*\[Xi]) + 
        (I/384)*E^((2*I)*\[Xi])*(-635 + 32*\[Nu]) + (I/640)*E^((6*I)*\[Xi])*
         (-3815 + 648*\[Nu]) + 3*\[Xi])) + x^2*\[Epsilon]^4*
     (e*(((I/24)*(624 - 235*\[Nu] + \[Nu]^2))/E^(I*\[Xi]) - 
        (I/24)*E^(I*\[Xi])*(624 - 235*\[Nu] + \[Nu]^2)) + 
      e^3*(((I/192)*(3756 - 1205*\[Nu] + 3*\[Nu]^2))/E^((3*I)*\[Xi]) - 
        (I/192)*E^((3*I)*\[Xi])*(3756 - 1205*\[Nu] + 3*\[Nu]^2) + 
        ((I/192)*(10140 - 4453*\[Nu] + 67*\[Nu]^2))/E^(I*\[Xi]) - 
        (I/192)*E^(I*\[Xi])*(10140 - 4453*\[Nu] + 67*\[Nu]^2)) + 
      e^5*(((I/4608)*(359304 - 166127*\[Nu] + 1445*\[Nu]^2))/E^(I*\[Xi]) - 
        (I/4608)*E^(I*\[Xi])*(359304 - 166127*\[Nu] + 1445*\[Nu]^2) + 
        ((I/3072)*(84240 - 31343*\[Nu] + 1797*\[Nu]^2))/E^((3*I)*\[Xi]) - 
        (I/3072)*E^((3*I)*\[Xi])*(84240 - 31343*\[Nu] + 1797*\[Nu]^2) - 
        ((I/46080)*(-1038048 + 338987*\[Nu] + 6055*\[Nu]^2))/
         E^((5*I)*\[Xi]) + (I/46080)*E^((5*I)*\[Xi])*
         (-1038048 + 338987*\[Nu] + 6055*\[Nu]^2)) + (27/2 - 7*\[Nu])*\[Xi] + 
      e^6*(((I/4608)*(207777 - 93122*\[Nu] + 1028*\[Nu]^2))/E^((2*I)*\[Xi]) - 
        (I/4608)*E^((2*I)*\[Xi])*(207777 - 93122*\[Nu] + 1028*\[Nu]^2) - 
        ((I/7680)*(-193707 + 66566*\[Nu] + 2172*\[Nu]^2))/E^((6*I)*\[Xi]) + 
        (I/7680)*E^((6*I)*\[Xi])*(-193707 + 66566*\[Nu] + 2172*\[Nu]^2) + 
        ((I/2880)*(71157 - 23291*\[Nu] + 2423*\[Nu]^2))/E^((4*I)*\[Xi]) - 
        (I/2880)*E^((4*I)*\[Xi])*(71157 - 23291*\[Nu] + 2423*\[Nu]^2) + 
        ((369 - 190*\[Nu])*\[Xi])/4) + 
      e^2*(((I/48)*(969 - 326*\[Nu] + 2*\[Nu]^2))/E^((2*I)*\[Xi]) - 
        (I/48)*E^((2*I)*\[Xi])*(969 - 326*\[Nu] + 2*\[Nu]^2) + 
        ((159 - 82*\[Nu])*\[Xi])/4) + 
      e^4*(((-I/288)*(-5925 + 1888*\[Nu] + 11*\[Nu]^2))/E^((4*I)*\[Xi]) + 
        (I/288)*E^((4*I)*\[Xi])*(-5925 + 1888*\[Nu] + 11*\[Nu]^2) + 
        ((I/144)*(4821 - 1952*\[Nu] + 62*\[Nu]^2))/E^((2*I)*\[Xi]) - 
        (I/144)*E^((2*I)*\[Xi])*(4821 - 1952*\[Nu] + 62*\[Nu]^2) + 
        (66 - 34*\[Nu])*\[Xi])) + 
    SO*(x^(3/2)*\[Epsilon]^3*(2*\[Xi]*(-2*\[Delta]*\[Chi]A + 
          (-2 + \[Nu])*\[Chi]S) + 
        e^5*(((-I/1920)*(7513*\[Delta]*\[Chi]A + (7513 - 2194*\[Nu])*
              \[Chi]S))/E^((5*I)*\[Xi]) + (I/1920)*E^((5*I)*\[Xi])*
           (7513*\[Delta]*\[Chi]A + (7513 - 2194*\[Nu])*\[Chi]S) - 
          ((I/192)*(1441*\[Delta]*\[Chi]A + (1441 - 658*\[Nu])*\[Chi]S))/
           E^(I*\[Xi]) + (I/192)*E^(I*\[Xi])*(1441*\[Delta]*\[Chi]A + 
            (1441 - 658*\[Nu])*\[Chi]S) - ((I/128)*(253*\[Delta]*\[Chi]A + 
             (253 - 122*\[Nu])*\[Chi]S))/E^((3*I)*\[Xi]) + 
          (I/128)*E^((3*I)*\[Xi])*(253*\[Delta]*\[Chi]A + (253 - 122*\[Nu])*
             \[Chi]S)) + e^3*(((-I/24)*(79*\[Delta]*\[Chi]A + 
             (79 - 26*\[Nu])*\[Chi]S))/E^((3*I)*\[Xi]) + 
          (I/24)*E^((3*I)*\[Xi])*(79*\[Delta]*\[Chi]A + (79 - 26*\[Nu])*
             \[Chi]S) - ((I/8)*(49*\[Delta]*\[Chi]A + (49 - 22*\[Nu])*
              \[Chi]S))/E^(I*\[Xi]) + (I/8)*E^(I*\[Xi])*
           (49*\[Delta]*\[Chi]A + (49 - 22*\[Nu])*\[Chi]S)) + 
        e*(((-I)*(5*\[Delta]*\[Chi]A + (5 - 2*\[Nu])*\[Chi]S))/E^(I*\[Xi]) + 
          I*E^(I*\[Xi])*(5*\[Delta]*\[Chi]A + (5 - 2*\[Nu])*\[Chi]S)) + 
        e^6*(((-I/960)*(4390*\[Delta]*\[Chi]A + (4390 - 1223*\[Nu])*\[Chi]S))/
           E^((6*I)*\[Xi]) + (I/960)*E^((6*I)*\[Xi])*(4390*\[Delta]*\[Chi]A + 
            (4390 - 1223*\[Nu])*\[Chi]S) - ((I/960)*(1030*\[Delta]*\[Chi]A + 
             (1030 - 643*\[Nu])*\[Chi]S))/E^((4*I)*\[Xi]) + 
          (I/960)*E^((4*I)*\[Xi])*(1030*\[Delta]*\[Chi]A + (1030 - 643*\[Nu])*
             \[Chi]S) - ((I/192)*(758*\[Delta]*\[Chi]A + (758 - 335*\[Nu])*
              \[Chi]S))/E^((2*I)*\[Xi]) + (I/192)*E^((2*I)*\[Xi])*
           (758*\[Delta]*\[Chi]A + (758 - 335*\[Nu])*\[Chi]S) - 
          (35*\[Xi]*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S))/8) + 
        e^4*(((-I/96)*(334*\[Delta]*\[Chi]A + (334 - 103*\[Nu])*\[Chi]S))/
           E^((4*I)*\[Xi]) + (I/96)*E^((4*I)*\[Xi])*(334*\[Delta]*\[Chi]A + 
            (334 - 103*\[Nu])*\[Chi]S) - ((I/12)*(38*\[Delta]*\[Chi]A + 
             (38 - 17*\[Nu])*\[Chi]S))/E^((2*I)*\[Xi]) + 
          (I/12)*E^((2*I)*\[Xi])*(38*\[Delta]*\[Chi]A + (38 - 17*\[Nu])*
             \[Chi]S) - (15*\[Xi]*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*
              \[Chi]S))/4) + e^2*(((-I/4)*(14*\[Delta]*\[Chi]A + 
             (14 - 5*\[Nu])*\[Chi]S))/E^((2*I)*\[Xi]) + (I/4)*E^((2*I)*\[Xi])*
           (14*\[Delta]*\[Chi]A + (14 - 5*\[Nu])*\[Chi]S) + 
          3*\[Xi]*(-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))) + 
      x^(5/2)*\[Epsilon]^5*((\[Xi]*(17*\[Delta]*(-4 + \[Nu])*\[Chi]A + 
           (-68 + 81*\[Nu] - 4*\[Nu]^2)*\[Chi]S))/2 + 
        e^2*((\[Xi]*(\[Delta]*(-460 + 141*\[Nu])*\[Chi]A + 
             (-460 + 521*\[Nu] - 48*\[Nu]^2)*\[Chi]S))/4 + 
          ((I/48)*(\[Delta]*(-1988 + 365*\[Nu])*\[Chi]A + 
             (-1988 + 1737*\[Nu] - 28*\[Nu]^2)*\[Chi]S))/E^((2*I)*\[Xi]) - 
          (I/48)*E^((2*I)*\[Xi])*(\[Delta]*(-1988 + 365*\[Nu])*\[Chi]A + 
            (-1988 + 1737*\[Nu] - 28*\[Nu]^2)*\[Chi]S)) + 
        e*(((I/6)*(\[Delta]*(-346 + 73*\[Nu])*\[Chi]A + 
             (-346 + 351*\[Nu] - 14*\[Nu]^2)*\[Chi]S))/E^(I*\[Xi]) - 
          (I/6)*E^(I*\[Xi])*(\[Delta]*(-346 + 73*\[Nu])*\[Chi]A + 
            (-346 + 351*\[Nu] - 14*\[Nu]^2)*\[Chi]S)) + 
        e^3*(((I/48)*(\[Delta]*(-7414 + 2083*\[Nu])*\[Chi]A + 
             (-7414 + 7893*\[Nu] - 650*\[Nu]^2)*\[Chi]S))/E^(I*\[Xi]) - 
          (I/48)*E^(I*\[Xi])*(\[Delta]*(-7414 + 2083*\[Nu])*\[Chi]A + 
            (-7414 + 7893*\[Nu] - 650*\[Nu]^2)*\[Chi]S) + 
          ((I/48)*(\[Delta]*(-1790 + 279*\[Nu])*\[Chi]A + 
             (-1790 + 1297*\[Nu] + 30*\[Nu]^2)*\[Chi]S))/E^((3*I)*\[Xi]) - 
          (I/48)*E^((3*I)*\[Xi])*(\[Delta]*(-1790 + 279*\[Nu])*\[Chi]A + 
            (-1790 + 1297*\[Nu] + 30*\[Nu]^2)*\[Chi]S)) + 
        e^4*(((I/144)*(\[Delta]*(-14528 + 4013*\[Nu])*\[Chi]A + 
             (-14528 + 15069*\[Nu] - 1210*\[Nu]^2)*\[Chi]S))/
           E^((2*I)*\[Xi]) - (I/144)*E^((2*I)*\[Xi])*
           (\[Delta]*(-14528 + 4013*\[Nu])*\[Chi]A + 
            (-14528 + 15069*\[Nu] - 1210*\[Nu]^2)*\[Chi]S) + 
          (5*\[Xi]*(\[Delta]*(-716 + 231*\[Nu])*\[Chi]A + 
             (-716 + 799*\[Nu] - 84*\[Nu]^2)*\[Chi]S))/16 + 
          ((I/1152)*(5*\[Delta]*(-8348 + 1037*\[Nu])*\[Chi]A + 
             (-41740 + 22929*\[Nu] + 2272*\[Nu]^2)*\[Chi]S))/
           E^((4*I)*\[Xi]) - (I/1152)*E^((4*I)*\[Xi])*
           (5*\[Delta]*(-8348 + 1037*\[Nu])*\[Chi]A + 
            (-41740 + 22929*\[Nu] + 2272*\[Nu]^2)*\[Chi]S)) + 
        e^6*(((I/11520)*(\[Delta]*(-972164 + 279971*\[Nu])*\[Chi]A + 
             (-972164 + 1047711*\[Nu] - 93208*\[Nu]^2)*\[Chi]S))/
           E^((4*I)*\[Xi]) - (I/11520)*E^((4*I)*\[Xi])*
           (\[Delta]*(-972164 + 279971*\[Nu])*\[Chi]A + 
            (-972164 + 1047711*\[Nu] - 93208*\[Nu]^2)*\[Chi]S) + 
          ((I/2304)*(5*\[Delta]*(-74696 + 21899*\[Nu])*\[Chi]A + 
             (-373480 + 390543*\[Nu] - 35654*\[Nu]^2)*\[Chi]S))/
           E^((2*I)*\[Xi]) - (I/2304)*E^((2*I)*\[Xi])*
           (5*\[Delta]*(-74696 + 21899*\[Nu])*\[Chi]A + 
            (-373480 + 390543*\[Nu] - 35654*\[Nu]^2)*\[Chi]S) + 
          (35*\[Xi]*(\[Delta]*(-324 + 107*\[Nu])*\[Chi]A + 
             (-324 + 359*\[Nu] - 40*\[Nu]^2)*\[Chi]S))/32 + 
          ((I/3840)*(\[Delta]*(-139048 + 4597*\[Nu])*\[Chi]A + 
             (-139048 + 4077*\[Nu] + 23494*\[Nu]^2)*\[Chi]S))/
           E^((6*I)*\[Xi]) - (I/3840)*E^((6*I)*\[Xi])*
           (\[Delta]*(-139048 + 4597*\[Nu])*\[Chi]A + 
            (-139048 + 4077*\[Nu] + 23494*\[Nu]^2)*\[Chi]S)) + 
        e^5*(((I/1152)*(\[Delta]*(-311882 + 93569*\[Nu])*\[Chi]A + 
             (-311882 + 333543*\[Nu] - 31582*\[Nu]^2)*\[Chi]S))/E^(I*\[Xi]) - 
          (I/1152)*E^(I*\[Xi])*(\[Delta]*(-311882 + 93569*\[Nu])*\[Chi]A + 
            (-311882 + 333543*\[Nu] - 31582*\[Nu]^2)*\[Chi]S) + 
          ((I/768)*(\[Delta]*(-66722 + 18613*\[Nu])*\[Chi]A + 
             (-66722 + 69395*\[Nu] - 5726*\[Nu]^2)*\[Chi]S))/
           E^((3*I)*\[Xi]) - (I/768)*E^((3*I)*\[Xi])*
           (\[Delta]*(-66722 + 18613*\[Nu])*\[Chi]A + 
            (-66722 + 69395*\[Nu] - 5726*\[Nu]^2)*\[Chi]S) + 
          ((I/11520)*(\[Delta]*(-416954 + 35489*\[Nu])*\[Chi]A + 
             (-416954 + 136839*\[Nu] + 42890*\[Nu]^2)*\[Chi]S))/
           E^((5*I)*\[Xi]) - (I/11520)*E^((5*I)*\[Xi])*
           (\[Delta]*(-416954 + 35489*\[Nu])*\[Chi]A + 
            (-416954 + 136839*\[Nu] + 42890*\[Nu]^2)*\[Chi]S)))) + 
    SO^2*(x^2*\[Epsilon]^4*((3*\[Xi]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
           \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
           \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/2 + 
        e*(((2*I)*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/E^(I*\[Xi]) - (2*I)*E^(I*\[Xi])*
           (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))) + 
        e^3*((((27*I)/8)*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/E^(I*\[Xi]) - ((27*I)/8)*E^(I*\[Xi])*
           (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)) + 
          (((11*I)/8)*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/E^((3*I)*\[Xi]) - ((11*I)/8)*E^((3*I)*\[Xi])*
           (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))) + 
        e^5*((((29*I)/6)*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/E^(I*\[Xi]) - ((29*I)/6)*E^(I*\[Xi])*
           (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)) + 
          (((23*I)/16)*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/E^((3*I)*\[Xi]) - ((23*I)/16)*E^((3*I)*\[Xi])*
           (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)) + 
          (((401*I)/240)*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/E^((5*I)*\[Xi]) - ((401*I)/240)*E^((5*I)*\[Xi])*
           (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))) + 
        e^2*((((23*I)/16)*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/E^((2*I)*\[Xi]) - ((23*I)/16)*E^((2*I)*\[Xi])*
           (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)) + 
          3*\[Xi]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
            4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
              2*\[Chi]A*\[Chi]S))) + 
        e^4*((((187*I)/96)*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/E^((2*I)*\[Xi]) - ((187*I)/96)*E^((2*I)*\[Xi])*
           (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)) + 
          (((565*I)/384)*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/E^((4*I)*\[Xi]) - ((565*I)/384)*E^((4*I)*\[Xi])*
           (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)) + 
          (9*\[Xi]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/2) + 
        e^6*((((689*I)/256)*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/E^((2*I)*\[Xi]) - ((689*I)/256)*E^((2*I)*\[Xi])*
           (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)) + 
          (((707*I)/640)*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/E^((4*I)*\[Xi]) - ((707*I)/640)*E^((4*I)*\[Xi])*
           (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)) + 
          (((2519*I)/1280)*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/E^((6*I)*\[Xi]) - ((2519*I)/1280)*E^((6*I)*\[Xi])*
           (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)) + 
          6*\[Xi]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
            4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
              2*\[Chi]A*\[Chi]S)))) + x^3*\[Epsilon]^6*
       ((\[Xi]*(\[Delta]*\[Kappa]A*(39 - 17*\[Nu]) + \[Kappa]S*
            (39 - 95*\[Nu] + 10*\[Nu]^2) + 67*\[Chi]A^2 - 
           279*\[Nu]*\[Chi]A^2 + 20*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(67 - 55*\[Nu])*\[Chi]A*\[Chi]S + 67*\[Chi]S^2 - 
           99*\[Nu]*\[Chi]S^2 + 28*\[Nu]^2*\[Chi]S^2))/2 + 
        e^3*(((-I/16)*(\[Delta]*\[Kappa]A*(-1742 + 743*\[Nu]) + 
             \[Kappa]S*(-1742 + 4227*\[Nu] - 626*\[Nu]^2) - 2898*\[Chi]A^2 + 
             12081*\[Nu]*\[Chi]A^2 - 1252*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-2898 + 2021*\[Nu])*\[Chi]A*\[Chi]S - 
             2898*\[Chi]S^2 + 3553*\[Nu]*\[Chi]S^2 - 836*\[Nu]^2*\[Chi]S^2))/
           E^(I*\[Xi]) + (I/16)*E^(I*\[Xi])*(\[Delta]*\[Kappa]A*
             (-1742 + 743*\[Nu]) + \[Kappa]S*(-1742 + 4227*\[Nu] - 
              626*\[Nu]^2) - 2898*\[Chi]A^2 + 12081*\[Nu]*\[Chi]A^2 - 
            1252*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(-2898 + 2021*\[Nu])*\[Chi]A*
             \[Chi]S - 2898*\[Chi]S^2 + 3553*\[Nu]*\[Chi]S^2 - 
            836*\[Nu]^2*\[Chi]S^2) - ((I/48)*(7*\[Delta]*\[Kappa]A*
              (-142 + 63*\[Nu]) + \[Kappa]S*(-994 + 2429*\[Nu] - 142*
                \[Nu]^2) - 1902*\[Chi]A^2 + 7747*\[Nu]*\[Chi]A^2 - 
             284*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(-1902 + 919*\[Nu])*\[Chi]A*
              \[Chi]S - 1902*\[Chi]S^2 + 1699*\[Nu]*\[Chi]S^2 - 
             332*\[Nu]^2*\[Chi]S^2))/E^((3*I)*\[Xi]) + (I/48)*E^((3*I)*\[Xi])*
           (7*\[Delta]*\[Kappa]A*(-142 + 63*\[Nu]) + \[Kappa]S*
             (-994 + 2429*\[Nu] - 142*\[Nu]^2) - 1902*\[Chi]A^2 + 
            7747*\[Nu]*\[Chi]A^2 - 284*\[Nu]^2*\[Chi]A^2 + 
            2*\[Delta]*(-1902 + 919*\[Nu])*\[Chi]A*\[Chi]S - 1902*\[Chi]S^2 + 
            1699*\[Nu]*\[Chi]S^2 - 332*\[Nu]^2*\[Chi]S^2)) + 
        e*(((-I/6)*(\[Delta]*\[Kappa]A*(-198 + 83*\[Nu]) + 
             \[Kappa]S*(-198 + 479*\[Nu] - 46*\[Nu]^2) - 
             2*(177 - 730*\[Nu] + 46*\[Nu]^2)*\[Chi]A^2 + 
             4*\[Delta]*(-177 + 121*\[Nu])*\[Chi]A*\[Chi]S - 
             2*(177 - 220*\[Nu] + 54*\[Nu]^2)*\[Chi]S^2))/E^(I*\[Xi]) + 
          (I/6)*E^(I*\[Xi])*(\[Delta]*\[Kappa]A*(-198 + 83*\[Nu]) + 
            \[Kappa]S*(-198 + 479*\[Nu] - 46*\[Nu]^2) - 
            2*(177 - 730*\[Nu] + 46*\[Nu]^2)*\[Chi]A^2 + 
            4*\[Delta]*(-177 + 121*\[Nu])*\[Chi]A*\[Chi]S - 
            2*(177 - 220*\[Nu] + 54*\[Nu]^2)*\[Chi]S^2)) + 
        e^4*(((-I/288)*(\[Delta]*\[Kappa]A*(-20949 + 8906*\[Nu]) + 
             \[Kappa]S*(-20949 + 50804*\[Nu] - 7528*\[Nu]^2) - 
             34953*\[Chi]A^2 + 145463*\[Nu]*\[Chi]A^2 - 15056*\[Nu]^2*
              \[Chi]A^2 + 2*\[Delta]*(-34953 + 23333*\[Nu])*\[Chi]A*\[Chi]S - 
             34953*\[Chi]S^2 + 41015*\[Nu]*\[Chi]S^2 - 9132*\[Nu]^2*
              \[Chi]S^2))/E^((2*I)*\[Xi]) + (I/288)*E^((2*I)*\[Xi])*
           (\[Delta]*\[Kappa]A*(-20949 + 8906*\[Nu]) + \[Kappa]S*
             (-20949 + 50804*\[Nu] - 7528*\[Nu]^2) - 34953*\[Chi]A^2 + 
            145463*\[Nu]*\[Chi]A^2 - 15056*\[Nu]^2*\[Chi]A^2 + 
            2*\[Delta]*(-34953 + 23333*\[Nu])*\[Chi]A*\[Chi]S - 
            34953*\[Chi]S^2 + 41015*\[Nu]*\[Chi]S^2 - 9132*\[Nu]^2*
             \[Chi]S^2) - ((I/1152)*(\[Delta]*\[Kappa]A*(-22533 + 10544*
                \[Nu]) + \[Kappa]S*(-22533 + 55610*\[Nu] - 1660*\[Nu]^2) - 
             45321*\[Chi]A^2 + 183263*\[Nu]*\[Chi]A^2 - 3320*\[Nu]^2*
              \[Chi]A^2 + 2*\[Delta]*(-45321 + 16601*\[Nu])*\[Chi]A*\[Chi]S - 
             45321*\[Chi]S^2 + 31223*\[Nu]*\[Chi]S^2 - 4884*\[Nu]^2*
              \[Chi]S^2))/E^((4*I)*\[Xi]) + (I/1152)*E^((4*I)*\[Xi])*
           (\[Delta]*\[Kappa]A*(-22533 + 10544*\[Nu]) + 
            \[Kappa]S*(-22533 + 55610*\[Nu] - 1660*\[Nu]^2) - 
            45321*\[Chi]A^2 + 183263*\[Nu]*\[Chi]A^2 - 3320*\[Nu]^2*
             \[Chi]A^2 + 2*\[Delta]*(-45321 + 16601*\[Nu])*\[Chi]A*\[Chi]S - 
            45321*\[Chi]S^2 + 31223*\[Nu]*\[Chi]S^2 - 4884*\[Nu]^2*
             \[Chi]S^2) - (3*\[Xi]*(\[Delta]*\[Kappa]A*(-234 + 103*\[Nu]) + 
             \[Kappa]S*(-234 + 571*\[Nu] - 98*\[Nu]^2) - 
             2*(185 - 778*\[Nu] + 98*\[Nu]^2)*\[Chi]A^2 + 20*\[Delta]*
              (-37 + 27*\[Nu])*\[Chi]A*\[Chi]S - 2*(185 - 232*\[Nu] + 56*
                \[Nu]^2)*\[Chi]S^2))/4) + 
        e^2*(((-I/48)*(\[Delta]*\[Kappa]A*(-1125 + 481*\[Nu]) + 
             \[Kappa]S*(-1125 + 2731*\[Nu] - 218*\[Nu]^2) - 2073*\[Chi]A^2 + 
             8497*\[Nu]*\[Chi]A^2 - 436*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-2073 + 1213*\[Nu])*\[Chi]A*\[Chi]S - 
             2073*\[Chi]S^2 + 2221*\[Nu]*\[Chi]S^2 - 492*\[Nu]^2*\[Chi]S^2))/
           E^((2*I)*\[Xi]) + (I/48)*E^((2*I)*\[Xi])*
           (\[Delta]*\[Kappa]A*(-1125 + 481*\[Nu]) + \[Kappa]S*
             (-1125 + 2731*\[Nu] - 218*\[Nu]^2) - 2073*\[Chi]A^2 + 
            8497*\[Nu]*\[Chi]A^2 - 436*\[Nu]^2*\[Chi]A^2 + 
            2*\[Delta]*(-2073 + 1213*\[Nu])*\[Chi]A*\[Chi]S - 
            2073*\[Chi]S^2 + 2221*\[Nu]*\[Chi]S^2 - 492*\[Nu]^2*\[Chi]S^2) - 
          (\[Xi]*(\[Delta]*\[Kappa]A*(-312 + 137*\[Nu]) + 
             \[Kappa]S*(-312 + 761*\[Nu] - 118*\[Nu]^2) + 
             8*\[Delta]*(-126 + 95*\[Nu])*\[Chi]A*\[Chi]S - 
             2*((252 - 1057*\[Nu] + 118*\[Nu]^2)*\[Chi]A^2 + (252 - 
                 331*\[Nu] + 84*\[Nu]^2)*\[Chi]S^2)))/4) + 
        e^6*(((-I/7680)*(\[Delta]*\[Kappa]A*(-134940 + 75941*\[Nu]) + 
             \[Kappa]S*(-134940 + 345821*\[Nu] + 24482*\[Nu]^2) - 
             321340*\[Chi]A^2 + 1273182*\[Nu]*\[Chi]A^2 + 48964*\[Nu]^2*
              \[Chi]A^2 + 8*\[Delta]*(-80335 + 3328*\[Nu])*\[Chi]A*\[Chi]S - 
             321340*\[Chi]S^2 + 38802*\[Nu]*\[Chi]S^2 + 20816*\[Nu]^2*
              \[Chi]S^2))/E^((6*I)*\[Xi]) + (I/7680)*E^((6*I)*\[Xi])*
           (\[Delta]*\[Kappa]A*(-134940 + 75941*\[Nu]) + 
            \[Kappa]S*(-134940 + 345821*\[Nu] + 24482*\[Nu]^2) - 
            321340*\[Chi]A^2 + 1273182*\[Nu]*\[Chi]A^2 + 48964*\[Nu]^2*
             \[Chi]A^2 + 8*\[Delta]*(-80335 + 3328*\[Nu])*\[Chi]A*\[Chi]S - 
            321340*\[Chi]S^2 + 38802*\[Nu]*\[Chi]S^2 + 20816*\[Nu]^2*
             \[Chi]S^2) - ((I/1920)*(\[Delta]*\[Kappa]A*(-120570 + 50659*
                \[Nu]) + \[Kappa]S*(-120570 + 291799*\[Nu] - 46402*\[Nu]^2) - 
             2*(98605 - 410764*\[Nu] + 46402*\[Nu]^2)*\[Chi]A^2 + 
             4*\[Delta]*(-98605 + 65139*\[Nu])*\[Chi]A*\[Chi]S - 
             2*(98605 - 113934*\[Nu] + 23418*\[Nu]^2)*\[Chi]S^2))/
           E^((4*I)*\[Xi]) + (I/1920)*E^((4*I)*\[Xi])*
           (\[Delta]*\[Kappa]A*(-120570 + 50659*\[Nu]) + 
            \[Kappa]S*(-120570 + 291799*\[Nu] - 46402*\[Nu]^2) - 
            2*(98605 - 410764*\[Nu] + 46402*\[Nu]^2)*\[Chi]A^2 + 
            4*\[Delta]*(-98605 + 65139*\[Nu])*\[Chi]A*\[Chi]S - 
            2*(98605 - 113934*\[Nu] + 23418*\[Nu]^2)*\[Chi]S^2) - 
          ((I/1536)*(\[Delta]*\[Kappa]A*(-210672 + 90317*\[Nu]) + 
             \[Kappa]S*(-210672 + 511661*\[Nu] - 83102*\[Nu]^2) - 
             2*(172304 - 719669*\[Nu] + 83102*\[Nu]^2)*\[Chi]A^2 + 
             8*\[Delta]*(-86152 + 58065*\[Nu])*\[Chi]A*\[Chi]S - 
             2*(172304 - 201807*\[Nu] + 44976*\[Nu]^2)*\[Chi]S^2))/
           E^((2*I)*\[Xi]) + (I/1536)*E^((2*I)*\[Xi])*
           (\[Delta]*\[Kappa]A*(-210672 + 90317*\[Nu]) + 
            \[Kappa]S*(-210672 + 511661*\[Nu] - 83102*\[Nu]^2) - 
            2*(172304 - 719669*\[Nu] + 83102*\[Nu]^2)*\[Chi]A^2 + 
            8*\[Delta]*(-86152 + 58065*\[Nu])*\[Chi]A*\[Chi]S - 
            2*(172304 - 201807*\[Nu] + 44976*\[Nu]^2)*\[Chi]S^2) - 
          (\[Xi]*(\[Delta]*\[Kappa]A*(-624 + 275*\[Nu]) + 
             \[Kappa]S*(-624 + 1523*\[Nu] - 274*\[Nu]^2) + 
             8*\[Delta]*(-244 + 175*\[Nu])*\[Chi]A*\[Chi]S - 
             2*((488 - 2055*\[Nu] + 274*\[Nu]^2)*\[Chi]A^2 + (488 - 
                 597*\[Nu] + 140*\[Nu]^2)*\[Chi]S^2)))/2) + 
        e^5*(((-I/11520)*(\[Delta]*\[Kappa]A*(-215850 + 108877*\[Nu]) + 
             \[Kappa]S*(-215850 + 540577*\[Nu] + 5566*\[Nu]^2) + 
             2*(-232515 + 931901*\[Nu] + 5566*\[Nu]^2)*\[Chi]A^2 + 
             4*\[Delta]*(-232515 + 52076*\[Nu])*\[Chi]A*\[Chi]S - 
             2*(232515 - 102311*\[Nu] + 7110*\[Nu]^2)*\[Chi]S^2))/
           E^((5*I)*\[Xi]) + (I/11520)*E^((5*I)*\[Xi])*
           (\[Delta]*\[Kappa]A*(-215850 + 108877*\[Nu]) + 
            \[Kappa]S*(-215850 + 540577*\[Nu] + 5566*\[Nu]^2) + 
            2*(-232515 + 931901*\[Nu] + 5566*\[Nu]^2)*\[Chi]A^2 + 
            4*\[Delta]*(-232515 + 52076*\[Nu])*\[Chi]A*\[Chi]S - 
            2*(232515 - 102311*\[Nu] + 7110*\[Nu]^2)*\[Chi]S^2) - 
          ((I/1152)*(\[Delta]*\[Kappa]A*(-255126 + 109783*\[Nu]) + 
             \[Kappa]S*(-255126 + 620035*\[Nu] - 101726*\[Nu]^2) - 
             2*(207369 - 867581*\[Nu] + 101726*\[Nu]^2)*\[Chi]A^2 + 
             4*\[Delta]*(-207369 + 144050*\[Nu])*\[Chi]A*\[Chi]S - 
             2*(207369 - 249995*\[Nu] + 57654*\[Nu]^2)*\[Chi]S^2))/
           E^(I*\[Xi]) + (I/1152)*E^(I*\[Xi])*(\[Delta]*\[Kappa]A*
             (-255126 + 109783*\[Nu]) + \[Kappa]S*(-255126 + 620035*\[Nu] - 
              101726*\[Nu]^2) - 2*(207369 - 867581*\[Nu] + 101726*\[Nu]^2)*
             \[Chi]A^2 + 4*\[Delta]*(-207369 + 144050*\[Nu])*\[Chi]A*
             \[Chi]S - 2*(207369 - 249995*\[Nu] + 57654*\[Nu]^2)*\[Chi]S^2) - 
          ((I/768)*(\[Delta]*\[Kappa]A*(-48954 + 20749*\[Nu]) + 
             \[Kappa]S*(-48954 + 118657*\[Nu] - 18002*\[Nu]^2) + 
             4*\[Delta]*(-40611 + 26576*\[Nu])*\[Chi]A*\[Chi]S - 
             2*((40611 - 168977*\[Nu] + 18002*\[Nu]^2)*\[Chi]A^2 + 
               (40611 - 46619*\[Nu] + 9918*\[Nu]^2)*\[Chi]S^2)))/
           E^((3*I)*\[Xi]) + (I/768)*E^((3*I)*\[Xi])*
           (\[Delta]*\[Kappa]A*(-48954 + 20749*\[Nu]) + 
            \[Kappa]S*(-48954 + 118657*\[Nu] - 18002*\[Nu]^2) + 
            4*\[Delta]*(-40611 + 26576*\[Nu])*\[Chi]A*\[Chi]S - 
            2*((40611 - 168977*\[Nu] + 18002*\[Nu]^2)*\[Chi]A^2 + 
              (40611 - 46619*\[Nu] + 9918*\[Nu]^2)*\[Chi]S^2)))))


(* ::Subsubsection::Closed:: *)
(*\[Xi] in terms of \[ScriptL]*)


\[Xi]\[LetterSpace]\[ScriptL] = \[ScriptL] + x^(3/2)*\[Epsilon]^3*
     (11/6 - 2*\[Gamma]E - 4*Log[2] - 2*Log[b] - 3*Log[x]) + 
    x^(5/2)*\[Epsilon]^5*(e^2*(-11/2 + 6*\[Gamma]E + Log[4096] + 6*Log[b] + 
        9*Log[x]) + e^4*(-11/2 + 6*\[Gamma]E + Log[4096] + 6*Log[b] + 
        9*Log[x]) + e^6*(-11/2 + 6*\[Gamma]E + Log[4096] + 6*Log[b] + 
        9*Log[x]) + ((6 + \[Nu])*(-11 + 12*\[Gamma]E + 24*Log[2] + 
         12*Log[b] + 18*Log[x]))/12) + SO*x^3*\[Epsilon]^6*
     (-((2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*(-11 + 12*\[Gamma]E + 
          24*Log[2] + 12*Log[b] + 18*Log[x]))/3 - 
      (e^2*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*(-11 + 12*\[Gamma]E + 
         24*Log[2] + 12*Log[b] + 18*Log[x]))/2 - 
      (5*e^4*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*
        (-11 + 12*\[Gamma]E + 24*Log[2] + 12*Log[b] + 18*Log[x]))/8 - 
      (35*e^6*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*
        (-11 + 12*\[Gamma]E + 24*Log[2] + 12*Log[b] + 18*Log[x]))/48) + 
    SO^2*x^(7/2)*\[Epsilon]^7*
     (((\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
         \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*
        (-11 + 12*\[Gamma]E + 24*Log[2] + 12*Log[b] + 18*Log[x]))/4 + 
      (e^2*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
         \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*
        (-11 + 12*\[Gamma]E + 24*Log[2] + 12*Log[b] + 18*Log[x]))/2 + 
      (3*e^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
         \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*
        (-11 + 12*\[Gamma]E + 24*Log[2] + 12*Log[b] + 18*Log[x]))/4 + 
      e^6*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
        \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*
       (-11 + 12*\[Gamma]E + 24*Log[2] + 12*Log[b] + 18*Log[x]))


(* ::Subsubsection::Closed:: *)
(*\[ScriptL] in terms of \[Xi]*)


\[ScriptL]\[LetterSpace]\[Xi] = \[Xi] + x^(3/2)*\[Epsilon]^3*
     (-11/6 + 2*\[Gamma]E + Log[16] + 2*Log[b] + 3*Log[x]) + 
    x^(5/2)*\[Epsilon]^5*(e^2*(11/2 - 6*\[Gamma]E - 12*Log[2] - 6*Log[b] - 
        9*Log[x]) + e^4*(11/2 - 6*\[Gamma]E - 12*Log[2] - 6*Log[b] - 
        9*Log[x]) + e^6*(11/2 - 6*\[Gamma]E - 12*Log[2] - 6*Log[b] - 
        9*Log[x]) - ((6 + \[Nu])*(-11 + 12*\[Gamma]E + 24*Log[2] + 
         12*Log[b] + 18*Log[x]))/12) + SO*x^3*\[Epsilon]^6*
     (((2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*(-11 + 12*\[Gamma]E + 
         24*Log[2] + 12*Log[b] + 18*Log[x]))/3 + 
      (e^2*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*(-11 + 12*\[Gamma]E + 
         24*Log[2] + 12*Log[b] + 18*Log[x]))/2 + 
      (5*e^4*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*
        (-11 + 12*\[Gamma]E + 24*Log[2] + 12*Log[b] + 18*Log[x]))/8 + 
      (35*e^6*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*
        (-11 + 12*\[Gamma]E + 24*Log[2] + 12*Log[b] + 18*Log[x]))/48) + 
    SO^2*x^(7/2)*\[Epsilon]^7*
     (-((\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
          \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*
         (-11 + 12*\[Gamma]E + 24*Log[2] + 12*Log[b] + 18*Log[x]))/4 - 
      (e^2*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
         \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*
        (-11 + 12*\[Gamma]E + 24*Log[2] + 12*Log[b] + 18*Log[x]))/2 - 
      (3*e^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
         \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*
        (-11 + 12*\[Gamma]E + 24*Log[2] + 12*Log[b] + 18*Log[x]))/4 - 
      e^6*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
        \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*
       (-11 + 12*\[Gamma]E + 24*Log[2] + 12*Log[b] + 18*Log[x]))


(* ::Subsection::Closed:: *)
(*instantaneous + tail*)


H\[Psi]InstTail[2, 0] = e*(1/(2*Sqrt[6]*E^(I*\[Xi])) + 
      E^(I*\[Xi])/(2*Sqrt[6])) + e^2*(1/(2*Sqrt[6]*E^((2*I)*\[Xi])) + 
      E^((2*I)*\[Xi])/(2*Sqrt[6])) + e^3*(-1/(16*Sqrt[6]*E^(I*\[Xi])) - 
      E^(I*\[Xi])/(16*Sqrt[6]) + (3*Sqrt[3/2])/(16*E^((3*I)*\[Xi])) + 
      (3*Sqrt[3/2]*E^((3*I)*\[Xi]))/16) + 
    e^4*(-1/(6*Sqrt[6]*E^((2*I)*\[Xi])) - E^((2*I)*\[Xi])/(6*Sqrt[6]) + 
      Sqrt[2/3]/(3*E^((4*I)*\[Xi])) + (Sqrt[2/3]*E^((4*I)*\[Xi]))/3) + 
    e^5*(1/(384*Sqrt[6]*E^(I*\[Xi])) + E^(I*\[Xi])/(384*Sqrt[6]) - 
      (27*Sqrt[3/2])/(256*E^((3*I)*\[Xi])) - (27*Sqrt[3/2]*E^((3*I)*\[Xi]))/
       256 + 625/(768*Sqrt[6]*E^((5*I)*\[Xi])) + (625*E^((5*I)*\[Xi]))/
       (768*Sqrt[6])) + e^6*(1/(48*Sqrt[6]*E^((2*I)*\[Xi])) + 
      E^((2*I)*\[Xi])/(48*Sqrt[6]) - (4*Sqrt[2/3])/(15*E^((4*I)*\[Xi])) - 
      (4*Sqrt[2/3]*E^((4*I)*\[Xi]))/15 + (27*Sqrt[3/2])/
       (80*E^((6*I)*\[Xi])) + (27*Sqrt[3/2]*E^((6*I)*\[Xi]))/80) + 
    x*\[Epsilon]^2*(e^4*((-93 + \[Nu])/(9*Sqrt[6]*E^((4*I)*\[Xi])) + 
        (E^((4*I)*\[Xi])*(-93 + \[Nu]))/(9*Sqrt[6]) - 
        (-57 + \[Nu])/(36*Sqrt[6]*E^((2*I)*\[Xi])) - 
        (E^((2*I)*\[Xi])*(-57 + \[Nu]))/(36*Sqrt[6])) + 
      e^3*((Sqrt[3/2]*(-495 + \[Nu]))/(224*E^((3*I)*\[Xi])) + 
        (Sqrt[3/2]*E^((3*I)*\[Xi])*(-495 + \[Nu]))/224 - 
        (177 + \[Nu])/(672*Sqrt[6]*E^(I*\[Xi])) - (E^(I*\[Xi])*(177 + \[Nu]))/
         (672*Sqrt[6])) + e^2*(-(339 + 5*\[Nu])/(84*Sqrt[6]*
          E^((2*I)*\[Xi])) - (E^((2*I)*\[Xi])*(339 + 5*\[Nu]))/
         (84*Sqrt[6])) + e*(-(183 + 11*\[Nu])/(84*Sqrt[6]*E^(I*\[Xi])) - 
        (E^(I*\[Xi])*(183 + 11*\[Nu]))/(84*Sqrt[6])) + 
      e^6*((-2307 + 19*\[Nu])/(2016*Sqrt[6]*E^((2*I)*\[Xi])) + 
        (E^((2*I)*\[Xi])*(-2307 + 19*\[Nu]))/(2016*Sqrt[6]) + 
        (9*Sqrt[3/2]*(-963 + 19*\[Nu]))/(1120*E^((6*I)*\[Xi])) + 
        (9*Sqrt[3/2]*E^((6*I)*\[Xi])*(-963 + 19*\[Nu]))/1120 - 
        (2*Sqrt[2/3]*(-858 + 19*\[Nu]))/(315*E^((4*I)*\[Xi])) - 
        (2*Sqrt[2/3]*E^((4*I)*\[Xi])*(-858 + 19*\[Nu]))/315) + 
      e^5*((-14919 + 13*\[Nu])/(16128*Sqrt[6]*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*(-14919 + 13*\[Nu]))/(16128*Sqrt[6]) + 
        (625*(-807 + 13*\[Nu]))/(32256*Sqrt[6]*E^((5*I)*\[Xi])) + 
        (625*E^((5*I)*\[Xi])*(-807 + 13*\[Nu]))/(32256*Sqrt[6]) - 
        (3*Sqrt[3/2]*(-1973 + 39*\[Nu]))/(3584*E^((3*I)*\[Xi])) - 
        (3*Sqrt[3/2]*E^((3*I)*\[Xi])*(-1973 + 39*\[Nu]))/3584)) + 
    x^2*\[Epsilon]^4*(e^2*((-2357 + (3389 - 265*\[Nu])*\[Nu])/
         (336*Sqrt[6]*E^((2*I)*\[Xi])) + 
        (E^((2*I)*\[Xi])*(-2357 + (3389 - 265*\[Nu])*\[Nu]))/(336*Sqrt[6])) + 
      e*((-856 + (4979 - 95*\[Nu])*\[Nu])/(1008*Sqrt[6]*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*(-856 + (4979 - 95*\[Nu])*\[Nu]))/(1008*Sqrt[6])) + 
      e^3*(-(Sqrt[3/2]*(5544 + \[Nu]*(-5145 + 709*\[Nu])))/
         (896*E^((3*I)*\[Xi])) - (Sqrt[3/2]*E^((3*I)*\[Xi])*
          (5544 + \[Nu]*(-5145 + 709*\[Nu])))/896 + 
        (-11016 + \[Nu]*(24711 + 709*\[Nu]))/(2688*Sqrt[6]*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*(-11016 + \[Nu]*(24711 + 709*\[Nu])))/(2688*Sqrt[6])) + 
      e^4*((-28594 + (20405 - 4091*\[Nu])*\[Nu])/(756*Sqrt[6]*
          E^((4*I)*\[Xi])) + (E^((4*I)*\[Xi])*(-28594 + (20405 - 4091*\[Nu])*
            \[Nu]))/(756*Sqrt[6]) + (6535 + \[Nu]*(30841 + 4091*\[Nu]))/
         (3024*Sqrt[6]*E^((2*I)*\[Xi])) + 
        (E^((2*I)*\[Xi])*(6535 + \[Nu]*(30841 + 4091*\[Nu])))/
         (3024*Sqrt[6])) + e^6*((-72872 + 5*(18719 - 661*\[Nu])*\[Nu])/
         (8064*Sqrt[6]*E^((2*I)*\[Xi])) + 
        (E^((2*I)*\[Xi])*(-72872 + 5*(18719 - 661*\[Nu])*\[Nu]))/
         (8064*Sqrt[6]) - (9*Sqrt[3/2]*(19312 + \[Nu]*(-9579 + 3305*\[Nu])))/
         (4480*E^((6*I)*\[Xi])) - (9*Sqrt[3/2]*E^((6*I)*\[Xi])*
          (19312 + \[Nu]*(-9579 + 3305*\[Nu])))/4480 + 
        (50333 + 2*\[Nu]*(7057 + 6610*\[Nu]))/(1260*Sqrt[6]*
          E^((4*I)*\[Xi])) + (E^((4*I)*\[Xi])*(50333 + 
           2*\[Nu]*(7057 + 6610*\[Nu])))/(1260*Sqrt[6])) + 
      e^5*((-717696 + (819745 - 2229*\[Nu])*\[Nu])/(64512*Sqrt[6]*
          E^(I*\[Xi])) + (E^(I*\[Xi])*(-717696 + (819745 - 2229*\[Nu])*
            \[Nu]))/(64512*Sqrt[6]) + 
        (3*Sqrt[3/2]*(23128 + \[Nu]*(17645 + 6687*\[Nu])))/
         (14336*E^((3*I)*\[Xi])) + (3*Sqrt[3/2]*E^((3*I)*\[Xi])*
          (23128 + \[Nu]*(17645 + 6687*\[Nu])))/14336 - 
        (25*(353880 + \[Nu]*(-207337 + 55725*\[Nu])))/(129024*Sqrt[6]*
          E^((5*I)*\[Xi])) - (25*E^((5*I)*\[Xi])*(353880 + 
           \[Nu]*(-207337 + 55725*\[Nu])))/(129024*Sqrt[6]))) + 
    SO^2*(x^2*\[Epsilon]^4*
       (e*(-(Sqrt[3/2]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
              4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
                2*\[Chi]A*\[Chi]S)))/(4*E^(I*\[Xi])) - 
          (Sqrt[3/2]*E^(I*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/4) + 
        e^2*(-((\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S))/(Sqrt[6]*E^((2*I)*\[Xi]))) - 
          (E^((2*I)*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/Sqrt[6]) + 
        e^3*((-35*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(32*Sqrt[6]*E^(I*\[Xi])) - 
          (35*E^(I*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(32*Sqrt[6]) - (15*Sqrt[3/2]*(\[Kappa]S - 
             2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/
           (32*E^((3*I)*\[Xi])) - (15*Sqrt[3/2]*E^((3*I)*\[Xi])*
            (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
             \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/32) + 
        e^4*(-((\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S))/(Sqrt[6]*E^((2*I)*\[Xi]))) - 
          (E^((2*I)*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/Sqrt[6] - (Sqrt[2/3]*(\[Kappa]S - 
             2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/E^((4*I)*\[Xi]) - 
          Sqrt[2/3]*E^((4*I)*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
            \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
            \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))) + 
        e^5*((-1183*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(768*Sqrt[6]*E^(I*\[Xi])) - 
          (1183*E^(I*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(768*Sqrt[6]) - (147*Sqrt[3/2]*
            (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
             \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/
           (512*E^((3*I)*\[Xi])) - (147*Sqrt[3/2]*E^((3*I)*\[Xi])*
            (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
             \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/512 - 
          (4375*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(1536*Sqrt[6]*E^((5*I)*\[Xi])) - 
          (4375*E^((5*I)*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(1536*Sqrt[6])) + 
        e^6*((-17*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(12*Sqrt[6]*E^((2*I)*\[Xi])) - 
          (17*E^((2*I)*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(12*Sqrt[6]) - (4*Sqrt[2/3]*(\[Kappa]S - 
             2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/
           (15*E^((4*I)*\[Xi])) - (4*Sqrt[2/3]*E^((4*I)*\[Xi])*
            (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
             \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/15 - 
          (27*Sqrt[3/2]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(20*E^((6*I)*\[Xi])) - 
          (27*Sqrt[3/2]*E^((6*I)*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/20)) + 
      x^3*\[Epsilon]^6*(e^2*((3*\[Delta]*\[Kappa]A*(-911 + 631*\[Nu]) - 
            3*\[Kappa]S*(911 + \[Nu]*(-2453 + 958*\[Nu])) + 
            (-3125 + (18971 - 5748*\[Nu])*\[Nu])*\[Chi]A^2 + 
            10*\[Delta]*(-625 + 1283*\[Nu])*\[Chi]A*\[Chi]S + 
            (-3125 + (6359 - 4928*\[Nu])*\[Nu])*\[Chi]S^2)/
           (252*Sqrt[6]*E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*
            (3*\[Delta]*\[Kappa]A*(-911 + 631*\[Nu]) - 3*\[Kappa]S*
              (911 + \[Nu]*(-2453 + 958*\[Nu])) + 
             (-3125 + (18971 - 5748*\[Nu])*\[Nu])*\[Chi]A^2 + 
             10*\[Delta]*(-625 + 1283*\[Nu])*\[Chi]A*\[Chi]S + 
             (-3125 + (6359 - 4928*\[Nu])*\[Nu])*\[Chi]S^2))/(252*Sqrt[6])) + 
        e*((3*\[Delta]*\[Kappa]A*(-637 + 671*\[Nu]) - 3*\[Kappa]S*
             (637 + \[Nu]*(-1945 + 686*\[Nu])) + 
            (-1015 + (8173 - 4116*\[Nu])*\[Nu])*\[Chi]A^2 + 
            2*\[Delta]*(-1015 + 4883*\[Nu])*\[Chi]A*\[Chi]S + 
            (-1015 + (5653 - 3472*\[Nu])*\[Nu])*\[Chi]S^2)/
           (504*Sqrt[6]*E^(I*\[Xi])) + (E^(I*\[Xi])*
            (3*\[Delta]*\[Kappa]A*(-637 + 671*\[Nu]) - 3*\[Kappa]S*
              (637 + \[Nu]*(-1945 + 686*\[Nu])) + 
             (-1015 + (8173 - 4116*\[Nu])*\[Nu])*\[Chi]A^2 + 
             2*\[Delta]*(-1015 + 4883*\[Nu])*\[Chi]A*\[Chi]S + 
             (-1015 + (5653 - 3472*\[Nu])*\[Nu])*\[Chi]S^2))/(504*Sqrt[6])) + 
        e^3*((-70287*\[Kappa]S + 9*\[Kappa]S*(21593 - 5726*\[Nu])*\[Nu] + 
            \[Delta]*\[Kappa]A*(-70287 + 53763*\[Nu]) + 
            (282845 - 103068*\[Nu])*\[Nu]*\[Chi]A^2 + 2*\[Delta]*
             (-60095 + 101251*\[Nu])*\[Chi]A*\[Chi]S + (160037 - 41216*\[Nu])*
             \[Nu]*\[Chi]S^2 - 60095*(\[Chi]A^2 + \[Chi]S^2))/
           (4032*Sqrt[6]*E^(I*\[Xi])) + (E^(I*\[Xi])*(-70287*\[Kappa]S + 
             9*\[Kappa]S*(21593 - 5726*\[Nu])*\[Nu] + \[Delta]*\[Kappa]A*
              (-70287 + 53763*\[Nu]) + (282845 - 103068*\[Nu])*\[Nu]*
              \[Chi]A^2 + 2*\[Delta]*(-60095 + 101251*\[Nu])*\[Chi]A*
              \[Chi]S + (160037 - 41216*\[Nu])*\[Nu]*\[Chi]S^2 - 
             60095*(\[Chi]A^2 + \[Chi]S^2)))/(4032*Sqrt[6]) + 
          (-3*\[Kappa]S*(3475 + \[Nu]*(-8917 + 3686*\[Nu])) + 
            (82363 - 22116*\[Nu])*\[Nu]*\[Chi]A^2 + 7*(2965 - 2752*\[Nu])*
             \[Nu]*\[Chi]S^2 + \[Delta]*(3*\[Kappa]A*(-3475 + 1967*\[Nu]) + 
              2*(-13897 + 23765*\[Nu])*\[Chi]A*\[Chi]S) - 
            13897*(\[Chi]A^2 + \[Chi]S^2))/(448*Sqrt[6]*E^((3*I)*\[Xi])) + 
          (E^((3*I)*\[Xi])*(-3*\[Kappa]S*(3475 + \[Nu]*(-8917 + 
                 3686*\[Nu])) + (82363 - 22116*\[Nu])*\[Nu]*\[Chi]A^2 + 
             7*(2965 - 2752*\[Nu])*\[Nu]*\[Chi]S^2 + \[Delta]*
              (3*\[Kappa]A*(-3475 + 1967*\[Nu]) + 2*(-13897 + 23765*\[Nu])*
                \[Chi]A*\[Chi]S) - 13897*(\[Chi]A^2 + \[Chi]S^2)))/
           (448*Sqrt[6])) + 
        e^5*(-(Sqrt[3/2]*(67911*\[Kappa]S + \[Delta]*\[Kappa]A*(67911 - 
                51221*\[Nu]) + \[Kappa]S*\[Nu]*(-187043 + 40666*\[Nu]) + 
              \[Nu]*(-228981 + 81332*\[Nu])*\[Chi]A^2 + 2*\[Delta]*(61639 - 
                87523*\[Nu])*\[Chi]A*\[Chi]S + \[Nu]*(-192621 + 11536*\[Nu])*
               \[Chi]S^2 + 61639*(\[Chi]A^2 + \[Chi]S^2)))/
           (7168*E^((3*I)*\[Xi])) - (Sqrt[3/2]*E^((3*I)*\[Xi])*
            (67911*\[Kappa]S + \[Delta]*\[Kappa]A*(67911 - 51221*\[Nu]) + 
             \[Kappa]S*\[Nu]*(-187043 + 40666*\[Nu]) + 
             \[Nu]*(-228981 + 81332*\[Nu])*\[Chi]A^2 + 2*\[Delta]*
              (61639 - 87523*\[Nu])*\[Chi]A*\[Chi]S + 
             \[Nu]*(-192621 + 11536*\[Nu])*\[Chi]S^2 + 
             61639*(\[Chi]A^2 + \[Chi]S^2)))/7168 + 
          (3*\[Delta]*\[Kappa]A*(-1305361 + 888575*\[Nu]) - 
            3*\[Kappa]S*(1305361 + \[Nu]*(-3499297 + 908606*\[Nu])) + 
            2*\[Delta]*(-3863443 + 5068895*\[Nu])*\[Chi]A*\[Chi]S - 
            3863443*(\[Chi]A^2 + \[Chi]S^2) + 
            \[Nu]*((17635489 - 5451636*\[Nu])*\[Chi]A^2 + (7956073 - 
                1942192*\[Nu])*\[Chi]S^2))/(96768*Sqrt[6]*E^(I*\[Xi])) + 
          (E^(I*\[Xi])*(3*\[Delta]*\[Kappa]A*(-1305361 + 888575*\[Nu]) - 
             3*\[Kappa]S*(1305361 + \[Nu]*(-3499297 + 908606*\[Nu])) + 
             2*\[Delta]*(-3863443 + 5068895*\[Nu])*\[Chi]A*\[Chi]S - 
             3863443*(\[Chi]A^2 + \[Chi]S^2) + \[Nu]*((17635489 - 
                 5451636*\[Nu])*\[Chi]A^2 + (7956073 - 1942192*\[Nu])*
                \[Chi]S^2)))/(96768*Sqrt[6]) + 
          (25*(3*\[Delta]*\[Kappa]A*(-191689 + 88355*\[Nu]) - 
             3*\[Kappa]S*(191689 + \[Nu]*(-471733 + 215606*\[Nu])) + 
             2*\[Delta]*(-849691 + 1315079*\[Nu])*\[Chi]A*\[Chi]S - 
             849691*(\[Chi]A^2 + \[Chi]S^2) + \[Nu]*((5124001 - 1293636*
                  \[Nu])*\[Chi]A^2 + (904921 - 1156624*\[Nu])*\[Chi]S^2)))/
           (193536*Sqrt[6]*E^((5*I)*\[Xi])) + (25*E^((5*I)*\[Xi])*
            (3*\[Delta]*\[Kappa]A*(-191689 + 88355*\[Nu]) - 
             3*\[Kappa]S*(191689 + \[Nu]*(-471733 + 215606*\[Nu])) + 
             2*\[Delta]*(-849691 + 1315079*\[Nu])*\[Chi]A*\[Chi]S - 
             849691*(\[Chi]A^2 + \[Chi]S^2) + \[Nu]*((5124001 - 1293636*
                  \[Nu])*\[Chi]A^2 + (904921 - 1156624*\[Nu])*\[Chi]S^2)))/
           (193536*Sqrt[6])) + 
        e^6*((15*\[Delta]*\[Kappa]A*(-68578 + 47545*\[Nu]) - 
            3*\[Kappa]S*(342890 + 3*\[Nu]*(-307835 + 84578*\[Nu])) + 
            52*(90362 - 29277*\[Nu])*\[Nu]*\[Chi]A^2 + 44*\[Delta]*
             (-46241 + 64723*\[Nu])*\[Chi]A*\[Chi]S + 
            4*(554549 - 140798*\[Nu])*\[Nu]*\[Chi]S^2 - 
            1017302*(\[Chi]A^2 + \[Chi]S^2))/(24192*Sqrt[6]*
            E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*(15*\[Delta]*\[Kappa]A*
              (-68578 + 47545*\[Nu]) - 3*\[Kappa]S*(342890 + 3*\[Nu]*
                (-307835 + 84578*\[Nu])) + 52*(90362 - 29277*\[Nu])*\[Nu]*
              \[Chi]A^2 + 44*\[Delta]*(-46241 + 64723*\[Nu])*\[Chi]A*
              \[Chi]S + 4*(554549 - 140798*\[Nu])*\[Nu]*\[Chi]S^2 - 
             1017302*(\[Chi]A^2 + \[Chi]S^2)))/(24192*Sqrt[6]) + 
          (-30756*\[Kappa]S + 9*\[Kappa]S*(9619 - 1114*\[Nu])*\[Nu] + 
            3*\[Delta]*\[Kappa]A*(-10252 + 8353*\[Nu]) + 
            4*(11287 - 5013*\[Nu])*\[Nu]*\[Chi]A^2 + 4*\[Delta]*
             (-12242 + 13555*\[Nu])*\[Chi]A*\[Chi]S + 
            64*\[Nu]*(1672 + 287*\[Nu])*\[Chi]S^2 - 
            24484*(\[Chi]A^2 + \[Chi]S^2))/(945*Sqrt[6]*E^((4*I)*\[Xi])) + 
          (E^((4*I)*\[Xi])*(-30756*\[Kappa]S + 9*\[Kappa]S*(9619 - 1114*
                \[Nu])*\[Nu] + 3*\[Delta]*\[Kappa]A*(-10252 + 8353*\[Nu]) + 
             4*(11287 - 5013*\[Nu])*\[Nu]*\[Chi]A^2 + 4*\[Delta]*
              (-12242 + 13555*\[Nu])*\[Chi]A*\[Chi]S + 64*\[Nu]*
              (1672 + 287*\[Nu])*\[Chi]S^2 - 24484*(\[Chi]A^2 + \[Chi]S^2)))/
           (945*Sqrt[6]) + (3*Sqrt[3/2]*(\[Delta]*\[Kappa]A*(-60494 + 26111*
                \[Nu]) + \[Kappa]S*(-60494 + 3*(49033 - 23478*\[Nu])*\[Nu]) + 
             4*\[Delta]*(-45843 + 70025*\[Nu])*\[Chi]A*\[Chi]S - 
             91686*(\[Chi]A^2 + \[Chi]S^2) + 4*\[Nu]*((140338 - 35217*\[Nu])*
                \[Chi]A^2 + (21373 - 31822*\[Nu])*\[Chi]S^2)))/
           (4480*E^((6*I)*\[Xi])) + (3*Sqrt[3/2]*E^((6*I)*\[Xi])*
            (\[Delta]*\[Kappa]A*(-60494 + 26111*\[Nu]) + 
             \[Kappa]S*(-60494 + 3*(49033 - 23478*\[Nu])*\[Nu]) + 
             4*\[Delta]*(-45843 + 70025*\[Nu])*\[Chi]A*\[Chi]S - 
             91686*(\[Chi]A^2 + \[Chi]S^2) + 4*\[Nu]*((140338 - 35217*\[Nu])*
                \[Chi]A^2 + (21373 - 31822*\[Nu])*\[Chi]S^2)))/4480) + 
        e^4*((3*\[Delta]*\[Kappa]A*(-10898 + 5467*\[Nu]) + 
            \[Kappa]S*(-32694 + 81789*\[Nu] - 35598*\[Nu]^2) + 
            4*(69269 - 17799*\[Nu])*\[Nu]*\[Chi]A^2 + 4*\[Delta]*
             (-23263 + 37121*\[Nu])*\[Chi]A*\[Chi]S + 56*(1027 - 1123*\[Nu])*
             \[Nu]*\[Chi]S^2 - 46526*(\[Chi]A^2 + \[Chi]S^2))/
           (756*Sqrt[6]*E^((4*I)*\[Xi])) + (E^((4*I)*\[Xi])*
            (3*\[Delta]*\[Kappa]A*(-10898 + 5467*\[Nu]) + 
             \[Kappa]S*(-32694 + 81789*\[Nu] - 35598*\[Nu]^2) + 
             4*(69269 - 17799*\[Nu])*\[Nu]*\[Chi]A^2 + 4*\[Delta]*
              (-23263 + 37121*\[Nu])*\[Chi]A*\[Chi]S + 56*(1027 - 1123*\[Nu])*
              \[Nu]*\[Chi]S^2 - 46526*(\[Chi]A^2 + \[Chi]S^2)))/
           (756*Sqrt[6]) + (3*\[Delta]*\[Kappa]A*(-5587 + 4205*\[Nu]) - 
            3*\[Kappa]S*(5587 + \[Nu]*(-15379 + 4010*\[Nu])) + 
            2*\[Delta]*(-15193 + 24347*\[Nu])*\[Chi]A*\[Chi]S - 
            15193*(\[Chi]A^2 + \[Chi]S^2) + \[Nu]*((67171 - 24060*\[Nu])*
               \[Chi]A^2 + (42295 - 8344*\[Nu])*\[Chi]S^2))/
           (756*Sqrt[6]*E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*
            (3*\[Delta]*\[Kappa]A*(-5587 + 4205*\[Nu]) - 3*\[Kappa]S*
              (5587 + \[Nu]*(-15379 + 4010*\[Nu])) + 2*\[Delta]*
              (-15193 + 24347*\[Nu])*\[Chi]A*\[Chi]S - 
             15193*(\[Chi]A^2 + \[Chi]S^2) + \[Nu]*((67171 - 24060*\[Nu])*
                \[Chi]A^2 + (42295 - 8344*\[Nu])*\[Chi]S^2)))/
           (756*Sqrt[6])))) + x^(3/2)*\[Epsilon]^3*
     (e^2*(Pi/(Sqrt[6]*E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*Pi)/Sqrt[6]) + 
      e^3*((9*Sqrt[3/2]*E^((3*I)*\[Xi])*(Pi - (2*I)*Log[3/2]))/16 + 
        (9*Sqrt[3/2]*(Pi + (2*I)*Log[3/2]))/(16*E^((3*I)*\[Xi])) - 
        (Pi - (2*I)*Log[2])/(16*Sqrt[6]*E^(I*\[Xi])) - 
        (E^(I*\[Xi])*(Pi + (2*I)*Log[2]))/(16*Sqrt[6])) + 
      e*((Pi - (2*I)*Log[2])/(2*Sqrt[6]*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*(Pi + (2*I)*Log[2]))/(2*Sqrt[6])) + 
      e^4*(-Pi/(3*Sqrt[6]*E^((2*I)*\[Xi])) - (E^((2*I)*\[Xi])*Pi)/
         (3*Sqrt[6]) + (4*Sqrt[2/3]*E^((4*I)*\[Xi])*(Pi - (2*I)*Log[2]))/3 + 
        (4*Sqrt[2/3]*(Pi + (2*I)*Log[2]))/(3*E^((4*I)*\[Xi]))) + 
      e^5*((-81*Sqrt[3/2]*E^((3*I)*\[Xi])*(Pi - (2*I)*Log[3/2]))/256 - 
        (81*Sqrt[3/2]*(Pi + (2*I)*Log[3/2]))/(256*E^((3*I)*\[Xi])) + 
        (Pi - (2*I)*Log[2])/(384*Sqrt[6]*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*(Pi + (2*I)*Log[2]))/(384*Sqrt[6]) + 
        (3125*E^((5*I)*\[Xi])*(Pi - (2*I)*Log[5/2]))/(768*Sqrt[6]) + 
        (3125*(Pi + (2*I)*Log[5/2]))/(768*Sqrt[6]*E^((5*I)*\[Xi]))) + 
      e^6*(Pi/(24*Sqrt[6]*E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*Pi)/
         (24*Sqrt[6]) - (16*Sqrt[2/3]*E^((4*I)*\[Xi])*(Pi - (2*I)*Log[2]))/
         15 - (16*Sqrt[2/3]*(Pi + (2*I)*Log[2]))/(15*E^((4*I)*\[Xi])) + 
        (81*Sqrt[3/2]*E^((6*I)*\[Xi])*(Pi - (2*I)*Log[3]))/40 + 
        (81*Sqrt[3/2]*(Pi + (2*I)*Log[3]))/(40*E^((6*I)*\[Xi])))) + 
    SO*(x^(3/2)*\[Epsilon]^3*
       (e^3*((47*\[Delta]*\[Chi]A + (47 - 25*\[Nu])*\[Chi]S)/
           (24*Sqrt[6]*E^(I*\[Xi])) + (E^(I*\[Xi])*(47*\[Delta]*\[Chi]A + 
             (47 - 25*\[Nu])*\[Chi]S))/(24*Sqrt[6]) + 
          (Sqrt[3/2]*(13*\[Delta]*\[Chi]A + (13 - 11*\[Nu])*\[Chi]S))/
           (8*E^((3*I)*\[Xi])) + (Sqrt[3/2]*E^((3*I)*\[Xi])*
            (13*\[Delta]*\[Chi]A + (13 - 11*\[Nu])*\[Chi]S))/8) + 
        e^6*((40*\[Delta]*\[Chi]A + (40 - 23*\[Nu])*\[Chi]S)/
           (18*Sqrt[6]*E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*
            (40*\[Delta]*\[Chi]A + (40 - 23*\[Nu])*\[Chi]S))/(18*Sqrt[6]) + 
          (9*Sqrt[3/2]*(11*\[Delta]*\[Chi]A + (11 - 10*\[Nu])*\[Chi]S))/
           (20*E^((6*I)*\[Xi])) + (9*Sqrt[3/2]*E^((6*I)*\[Xi])*
            (11*\[Delta]*\[Chi]A + (11 - 10*\[Nu])*\[Chi]S))/20 - 
          (14*Sqrt[2/3]*(4*\[Delta]*\[Chi]A + (4 - 5*\[Nu])*\[Chi]S))/
           (45*E^((4*I)*\[Xi])) - (14*Sqrt[2/3]*E^((4*I)*\[Xi])*
            (4*\[Delta]*\[Chi]A + (4 - 5*\[Nu])*\[Chi]S))/45) + 
        e*((7*\[Delta]*\[Chi]A + (7 - 5*\[Nu])*\[Chi]S)/
           (3*Sqrt[6]*E^(I*\[Xi])) + (E^(I*\[Xi])*(7*\[Delta]*\[Chi]A + 
             (7 - 5*\[Nu])*\[Chi]S))/(3*Sqrt[6])) + 
        e^2*((Sqrt[2/3]*(5*\[Delta]*\[Chi]A + (5 - 4*\[Nu])*\[Chi]S))/
           (3*E^((2*I)*\[Xi])) + (Sqrt[2/3]*E^((2*I)*\[Xi])*
            (5*\[Delta]*\[Chi]A + (5 - 4*\[Nu])*\[Chi]S))/3) + 
        e^4*((4*Sqrt[2/3]*(8*\[Delta]*\[Chi]A + (8 - 7*\[Nu])*\[Chi]S))/
           (9*E^((4*I)*\[Xi])) + (4*Sqrt[2/3]*E^((4*I)*\[Xi])*
            (8*\[Delta]*\[Chi]A + (8 - 7*\[Nu])*\[Chi]S))/9 + 
          (11*\[Delta]*\[Chi]A + (11 - 4*\[Nu])*\[Chi]S)/
           (9*Sqrt[6]*E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*
            (11*\[Delta]*\[Chi]A + (11 - 4*\[Nu])*\[Chi]S))/(9*Sqrt[6])) + 
        e^5*((1423*\[Delta]*\[Chi]A + (1423 - 773*\[Nu])*\[Chi]S)/
           (576*Sqrt[6]*E^(I*\[Xi])) + (E^(I*\[Xi])*(1423*\[Delta]*\[Chi]A + 
             (1423 - 773*\[Nu])*\[Chi]S))/(576*Sqrt[6]) + 
          (625*(19*\[Delta]*\[Chi]A + (19 - 17*\[Nu])*\[Chi]S))/
           (1152*Sqrt[6]*E^((5*I)*\[Xi])) + (625*E^((5*I)*\[Xi])*
            (19*\[Delta]*\[Chi]A + (19 - 17*\[Nu])*\[Chi]S))/(1152*Sqrt[6]) - 
          (3*Sqrt[3/2]*(\[Delta]*\[Chi]A + \[Chi]S - 11*\[Nu]*\[Chi]S))/
           (128*E^((3*I)*\[Xi])) - (3*Sqrt[3/2]*E^((3*I)*\[Xi])*
            (\[Delta]*\[Chi]A + \[Chi]S - 11*\[Nu]*\[Chi]S))/128)) + 
      x^(5/2)*\[Epsilon]^5*(e^2*((\[Delta]*(471 - 1043*\[Nu])*\[Chi]A + 
            (471 + \[Nu]*(-662 + 835*\[Nu]))*\[Chi]S)/(63*Sqrt[6]*
            E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*(\[Delta]*(471 - 1043*\[Nu])*
              \[Chi]A + (471 + \[Nu]*(-662 + 835*\[Nu]))*\[Chi]S))/
           (63*Sqrt[6])) + e*((\[Delta]*(204 - 1691*\[Nu])*\[Chi]A + 
            (204 + \[Nu]*(-1829 + 964*\[Nu]))*\[Chi]S)/(252*Sqrt[6]*
            E^(I*\[Xi])) + (E^(I*\[Xi])*(\[Delta]*(204 - 1691*\[Nu])*
              \[Chi]A + (204 + \[Nu]*(-1829 + 964*\[Nu]))*\[Chi]S))/
           (252*Sqrt[6])) + e^3*((\[Delta]*(4608 - 7589*\[Nu])*\[Chi]A + 
            (4608 + \[Nu]*(-2387 + 7012*\[Nu]))*\[Chi]S)/
           (224*Sqrt[6]*E^((3*I)*\[Xi])) + (E^((3*I)*\[Xi])*
            (\[Delta]*(4608 - 7589*\[Nu])*\[Chi]A + 
             (4608 + \[Nu]*(-2387 + 7012*\[Nu]))*\[Chi]S))/(224*Sqrt[6]) + 
          (\[Delta]*(32592 - 29143*\[Nu])*\[Chi]A + 
            (32592 + \[Nu]*(-78289 + 9308*\[Nu]))*\[Chi]S)/
           (2016*Sqrt[6]*E^(I*\[Xi])) + (E^(I*\[Xi])*
            (\[Delta]*(32592 - 29143*\[Nu])*\[Chi]A + 
             (32592 + \[Nu]*(-78289 + 9308*\[Nu]))*\[Chi]S))/
           (2016*Sqrt[6])) + e^4*((\[Delta]*(5973 - 5342*\[Nu])*\[Chi]A + 
            (5973 + \[Nu]*(-21581 + 229*\[Nu]))*\[Chi]S)/
           (378*Sqrt[6]*E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*
            (\[Delta]*(5973 - 5342*\[Nu])*\[Chi]A + 
             (5973 + \[Nu]*(-21581 + 229*\[Nu]))*\[Chi]S))/(378*Sqrt[6]) + 
          (2*\[Delta]*(3999 - 5908*\[Nu])*\[Chi]A + 
            (7998 + \[Nu]*(-542 + 11917*\[Nu]))*\[Chi]S)/
           (189*Sqrt[6]*E^((4*I)*\[Xi])) + (E^((4*I)*\[Xi])*
            (2*\[Delta]*(3999 - 5908*\[Nu])*\[Chi]A + 
             (7998 + \[Nu]*(-542 + 11917*\[Nu]))*\[Chi]S))/(189*Sqrt[6])) + 
        e^6*((3*Sqrt[3/2]*(2*\[Delta]*(8031 - 11231*\[Nu])*\[Chi]A + 
             (16062 + \[Nu]*(8852 + 25307*\[Nu]))*\[Chi]S))/
           (1120*E^((6*I)*\[Xi])) + (3*Sqrt[3/2]*E^((6*I)*\[Xi])*
            (2*\[Delta]*(8031 - 11231*\[Nu])*\[Chi]A + 
             (16062 + \[Nu]*(8852 + 25307*\[Nu]))*\[Chi]S))/1120 + 
          (2*\[Delta]*(97509 - 85513*\[Nu])*\[Chi]A + 
            (195018 + \[Nu]*(-390476 + 74239*\[Nu]))*\[Chi]S)/
           (6048*Sqrt[6]*E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*
            (2*\[Delta]*(97509 - 85513*\[Nu])*\[Chi]A + 
             (195018 + \[Nu]*(-390476 + 74239*\[Nu]))*\[Chi]S))/
           (6048*Sqrt[6]) + (2*\[Delta]*(8877 + 12218*\[Nu])*\[Chi]A + 
            (17754 - \[Nu]*(322166 + 96041*\[Nu]))*\[Chi]S)/
           (1890*Sqrt[6]*E^((4*I)*\[Xi])) + (E^((4*I)*\[Xi])*
            (2*\[Delta]*(8877 + 12218*\[Nu])*\[Chi]A + 
             (17754 - \[Nu]*(322166 + 96041*\[Nu]))*\[Chi]S))/
           (1890*Sqrt[6])) + 
        e^5*((Sqrt[3/2]*(\[Delta]*(18404 - 9339*\[Nu])*\[Chi]A + 
             (18404 - \[Nu]*(114949 + 16644*\[Nu]))*\[Chi]S))/
           (3584*E^((3*I)*\[Xi])) + (Sqrt[3/2]*E^((3*I)*\[Xi])*
            (\[Delta]*(18404 - 9339*\[Nu])*\[Chi]A + 
             (18404 - \[Nu]*(114949 + 16644*\[Nu]))*\[Chi]S))/3584 + 
          (25*(\[Delta]*(296436 - 419687*\[Nu])*\[Chi]A + 
             (296436 + \[Nu]*(78727 + 450556*\[Nu]))*\[Chi]S))/
           (96768*Sqrt[6]*E^((5*I)*\[Xi])) + (25*E^((5*I)*\[Xi])*
            (\[Delta]*(296436 - 419687*\[Nu])*\[Chi]A + 
             (296436 + \[Nu]*(78727 + 450556*\[Nu]))*\[Chi]S))/
           (96768*Sqrt[6]) + (\[Delta]*(1791012 - 1307579*\[Nu])*\[Chi]A + 
            (1791012 + \[Nu]*(-3238229 + 510868*\[Nu]))*\[Chi]S)/
           (48384*Sqrt[6]*E^(I*\[Xi])) + 
          (E^(I*\[Xi])*(\[Delta]*(1791012 - 1307579*\[Nu])*\[Chi]A + 
             (1791012 + \[Nu]*(-3238229 + 510868*\[Nu]))*\[Chi]S))/
           (48384*Sqrt[6]))) + x^3*\[Epsilon]^6*
       (e^2*((Sqrt[2/3]*(4*(3*I + 4*Pi)*\[Delta]*\[Chi]A + Pi*(16 - 11*\[Nu])*
              \[Chi]S - (6*I)*(-2 + \[Nu])*\[Chi]S))/(3*E^((2*I)*\[Xi])) + 
          (Sqrt[2/3]*E^((2*I)*\[Xi])*(4*(-3*I + 4*Pi)*\[Delta]*\[Chi]A + 
             Pi*(16 - 11*\[Nu])*\[Chi]S + (6*I)*(-2 + \[Nu])*\[Chi]S))/3) + 
        e^3*((3*Sqrt[3/2]*(19*Pi*\[Delta]*\[Chi]A + Pi*(19 - 14*\[Nu])*
              \[Chi]S + (4*I)*\[Delta]*\[Chi]A*(3 + 19*ArcCoth[5]) + 
             (2*I)*\[Chi]S*(6 + 38*ArcCoth[5] - \[Nu]*(3 + 28*ArcCoth[5]))))/
           (8*E^((3*I)*\[Xi])) + (3*Sqrt[3/2]*E^((3*I)*\[Xi])*
            (19*Pi*\[Delta]*\[Chi]A + Pi*(19 - 14*\[Nu])*\[Chi]S - 
             (4*I)*\[Delta]*\[Chi]A*(3 + 19*ArcCoth[5]) + (2*I)*\[Chi]S*
              (-6 - 38*ArcCoth[5] + \[Nu]*(3 + 28*ArcCoth[5]))))/8 + 
          (E^(I*\[Xi])*(113*Pi*\[Delta]*\[Chi]A + Pi*(113 - 58*\[Nu])*
              \[Chi]S + (2*I)*\[Delta]*\[Chi]A*(-66 + 113*Log[2]) - 
             (2*I)*\[Chi]S*(66 - 113*Log[2] + \[Nu]*(-33 + 58*Log[2]))))/
           (24*Sqrt[6]) + (113*Pi*\[Delta]*\[Chi]A + Pi*(113 - 58*\[Nu])*
             \[Chi]S + (2*I)*\[Delta]*\[Chi]A*(66 - 113*Log[2]) + 
            (2*I)*\[Chi]S*(66 - 113*Log[2] + \[Nu]*(-33 + 58*Log[2])))/
           (24*Sqrt[6]*E^(I*\[Xi]))) + 
        e^6*((2*(102*I + 91*Pi)*\[Delta]*\[Chi]A + Pi*(182 - 97*\[Nu])*
             \[Chi]S - (102*I)*(-2 + \[Nu])*\[Chi]S)/(18*Sqrt[6]*
            E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*(2*(-102*I + 91*Pi)*\[Delta]*
              \[Chi]A + Pi*(182 - 97*\[Nu])*\[Chi]S + (102*I)*(-2 + \[Nu])*
              \[Chi]S))/(18*Sqrt[6]) + (27*Sqrt[3/2]*
            (Pi*(28 - 23*\[Nu])*\[Chi]S + 4*\[Delta]*\[Chi]A*
              (3*I + 7*Pi + (14*I)*Log[3]) - (2*I)*\[Chi]S*(-6 - 28*Log[3] + 
               \[Nu]*(3 + 23*Log[3]))))/(20*E^((6*I)*\[Xi])) + 
          (27*Sqrt[3/2]*E^((6*I)*\[Xi])*(28*Pi*\[Delta]*\[Chi]A + 
             Pi*(28 - 23*\[Nu])*\[Chi]S - (4*I)*\[Delta]*\[Chi]A*
              (3 + 14*Log[3]) + (2*I)*\[Chi]S*(-6 - 28*Log[3] + \[Nu]*
                (3 + 23*Log[3]))))/20 - (28*Sqrt[2/3]*
            (2*\[Delta]*\[Chi]A*(Pi + (2*I)*(-3 + Log[2])) + 
             \[Chi]S*(Pi*(2 - 7*\[Nu]) + (2*I)*(-6 + \[Nu]*(3 - 7*Log[2]) + 
                 Log[4]))))/(45*E^((4*I)*\[Xi])) + 
          (28*Sqrt[2/3]*E^((4*I)*\[Xi])*(Pi*(-2*\[Delta]*\[Chi]A + 
               (-2 + 7*\[Nu])*\[Chi]S) + (2*I)*(\[Delta]*\[Chi]A*
                (-6 + Log[4]) + \[Chi]S*(-6 + \[Nu]*(3 - 7*Log[2]) + 
                 Log[4]))))/45) + 
        e^5*((E^(I*\[Xi])*(Pi*(3373 - 1748*\[Nu])*\[Chi]S + \[Delta]*\[Chi]A*
              (-3900*I + 3373*Pi + (6746*I)*Log[2]) + (2*I)*\[Chi]S*
              (975*(-2 + \[Nu]) + (3373 - 1748*\[Nu])*Log[2])))/
           (576*Sqrt[6]) + (Pi*(3373 - 1748*\[Nu])*\[Chi]S + 
            \[Delta]*\[Chi]A*(3900*I + 3373*Pi - (6746*I)*Log[2]) + 
            (2*I)*\[Chi]S*(-975*(-2 + \[Nu]) + (-3373 + 1748*\[Nu])*Log[2]))/
           (576*Sqrt[6]*E^(I*\[Xi])) + (3125*(25*Pi*\[Delta]*\[Chi]A + 
             5*Pi*(5 - 4*\[Nu])*\[Chi]S + (2*I)*\[Delta]*\[Chi]A*
              (6 + 25*Log[5/2]) - (2*I)*\[Chi]S*(3*(-2 + \[Nu]) + 5*
                (-5 + 4*\[Nu])*Log[5/2])))/(1152*Sqrt[6]*E^((5*I)*\[Xi])) + 
          (3125*E^((5*I)*\[Xi])*(25*Pi*\[Delta]*\[Chi]A + 5*Pi*(5 - 4*\[Nu])*
              \[Chi]S - (2*I)*\[Delta]*\[Chi]A*(6 + 25*Log[5/2]) + 
             (2*I)*\[Chi]S*(3*(-2 + \[Nu]) + 5*(-5 + 4*\[Nu])*Log[5/2])))/
           (1152*Sqrt[6]) + (9*Sqrt[3/2]*(29*Pi*\[Delta]*\[Chi]A + 
             Pi*(29 - 4*\[Nu])*\[Chi]S + (4*I)*\[Delta]*\[Chi]A*
              (15 + 29*ArcCoth[5]) + (2*I)*\[Chi]S*(30 + 58*ArcCoth[5] - 
               \[Nu]*(15 + Log[81/16]))))/(128*E^((3*I)*\[Xi])) + 
          (9*Sqrt[3/2]*E^((3*I)*\[Xi])*(29*Pi*\[Delta]*\[Chi]A + 
             Pi*(29 - 4*\[Nu])*\[Chi]S - (4*I)*\[Delta]*\[Chi]A*
              (15 + 29*ArcCoth[5]) + (2*I)*\[Chi]S*(-30 - 58*ArcCoth[5] + 
               \[Nu]*(15 + Log[81/16]))))/128) + 
        e^4*((4*(21*I + 16*Pi)*\[Delta]*\[Chi]A + Pi*(64 - 29*\[Nu])*
             \[Chi]S - (42*I)*(-2 + \[Nu])*\[Chi]S)/(9*Sqrt[6]*
            E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*(4*(-21*I + 16*Pi)*\[Delta]*
              \[Chi]A + Pi*(64 - 29*\[Nu])*\[Chi]S + (42*I)*(-2 + \[Nu])*
              \[Chi]S))/(9*Sqrt[6]) + (8*Sqrt[2/3]*E^((4*I)*\[Xi])*
            (22*Pi*\[Delta]*\[Chi]A + Pi*(22 - 17*\[Nu])*\[Chi]S + 
             (2*I)*\[Chi]S*(-6 - 22*Log[2] + \[Nu]*(3 + 17*Log[2])) - 
             (4*I)*\[Delta]*\[Chi]A*(3 + Log[2048])))/9 + 
          (8*Sqrt[2/3]*(22*Pi*\[Delta]*\[Chi]A + Pi*(22 - 17*\[Nu])*\[Chi]S + 
             (2*I)*\[Chi]S*(6 + 22*Log[2] - \[Nu]*(3 + 17*Log[2])) + 
             (4*I)*\[Delta]*\[Chi]A*(3 + Log[2048])))/(9*E^((4*I)*\[Xi]))) + 
        e*((13*Pi*\[Delta]*\[Chi]A + Pi*(13 - 8*\[Nu])*\[Chi]S + 
            (2*I)*\[Delta]*\[Chi]A*(6 - 13*Log[2]) + (2*I)*\[Chi]S*
             (6 - 13*Log[2] + \[Nu]*(-3 + Log[256])))/(3*Sqrt[6]*
            E^(I*\[Xi])) + (E^(I*\[Xi])*(13*Pi*\[Delta]*\[Chi]A + 
             Pi*(13 - 8*\[Nu])*\[Chi]S - (2*I)*\[Chi]S*(6 - 13*Log[2] + \[Nu]*
                (-3 + Log[256])) + (2*I)*\[Delta]*\[Chi]A*(-6 + Log[8192])))/
           (3*Sqrt[6]))))
 
H\[Psi]InstTail[2, 1] = 
   Sqrt[x]*((I/3)*\[Delta] + e*(((I/3)*\[Delta])/E^(I*\[Xi]) + 
        (I/3)*E^(I*\[Xi])*\[Delta]) + 
      e^2*((((5*I)/12)*\[Delta])/E^((2*I)*\[Xi]) + ((5*I)/12)*E^((2*I)*\[Xi])*
         \[Delta]) + e^3*(((-I/24)*\[Delta])/E^(I*\[Xi]) - 
        (I/24)*E^(I*\[Xi])*\[Delta] + (((13*I)/24)*\[Delta])/
         E^((3*I)*\[Xi]) + ((13*I)/24)*E^((3*I)*\[Xi])*\[Delta]) + 
      e^4*((((-11*I)/72)*\[Delta])/E^((2*I)*\[Xi]) - 
        ((11*I)/72)*E^((2*I)*\[Xi])*\[Delta] + (((103*I)/144)*\[Delta])/
         E^((4*I)*\[Xi]) + ((103*I)/144)*E^((4*I)*\[Xi])*\[Delta]) + 
      e^5*((((5*I)/576)*\[Delta])/E^(I*\[Xi]) + ((5*I)/576)*E^(I*\[Xi])*
         \[Delta] - (((43*I)/128)*\[Delta])/E^((3*I)*\[Xi]) - 
        ((43*I)/128)*E^((3*I)*\[Xi])*\[Delta] + (((1097*I)/1152)*\[Delta])/
         E^((5*I)*\[Xi]) + ((1097*I)/1152)*E^((5*I)*\[Xi])*\[Delta]) + 
      e^6*((((17*I)/576)*\[Delta])/E^((2*I)*\[Xi]) + 
        ((17*I)/576)*E^((2*I)*\[Xi])*\[Delta] - (((451*I)/720)*\[Delta])/
         E^((4*I)*\[Xi]) - ((451*I)/720)*E^((4*I)*\[Xi])*\[Delta] + 
        (((1223*I)/960)*\[Delta])/E^((6*I)*\[Xi]) + ((1223*I)/960)*
         E^((6*I)*\[Xi])*\[Delta]))*\[Epsilon] + x^(3/2)*\[Epsilon]^3*
     ((I/84)*\[Delta]*(-17 + 20*\[Nu]) + 
      e*(((I/84)*\[Delta]*(-43 + 2*\[Nu]))/E^(I*\[Xi]) + 
        (I/28)*E^(I*\[Xi])*\[Delta]*(41 + 2*\[Nu])) + 
      e^2*(((3*I)/56)*\[Delta]*(-9 + 4*\[Nu]) - (I/168)*E^((2*I)*\[Xi])*
         \[Delta]*(-643 + 26*\[Nu]) - ((I/168)*\[Delta]*(187 + 46*\[Nu]))/
         E^((2*I)*\[Xi])) + e^3*((I/672)*E^(I*\[Xi])*\[Delta]*
         (-835 + 186*\[Nu]) + ((I/672)*\[Delta]*(-5 + 206*\[Nu]))/
         E^(I*\[Xi]) - (I/672)*E^((3*I)*\[Xi])*\[Delta]*(-5065 + 358*\[Nu]) - 
        ((I/672)*\[Delta]*(1409 + 514*\[Nu]))/E^((3*I)*\[Xi])) + 
      e^4*(((3*I)/56)*\[Delta]*(-9 + 4*\[Nu]) + 
        ((I/48)*\[Delta]*(31 + 32*\[Nu]))/E^((2*I)*\[Xi]) - 
        (I/448)*E^((4*I)*\[Xi])*\[Delta]*(-5963 + 516*\[Nu]) + 
        (I/1008)*E^((2*I)*\[Xi])*\[Delta]*(-3665 + 568*\[Nu]) - 
        ((I/4032)*\[Delta]*(14725 + 6292*\[Nu]))/E^((4*I)*\[Xi])) + 
      e^5*((I/16128)*E^(I*\[Xi])*\[Delta]*(-5813 + 2806*\[Nu]) + 
        ((I/16128)*\[Delta]*(-4651 + 2834*\[Nu]))/E^(I*\[Xi]) + 
        (I/3584)*E^((3*I)*\[Xi])*\[Delta]*(-31467 + 4274*\[Nu]) + 
        ((I/3584)*\[Delta]*(7211 + 5206*\[Nu]))/E^((3*I)*\[Xi]) - 
        (I/32256)*E^((5*I)*\[Xi])*\[Delta]*(-715031 + 69010*\[Nu]) - 
        ((I/32256)*\[Delta]*(195479 + 90950*\[Nu]))/E^((5*I)*\[Xi])) + 
      e^6*(((3*I)/56)*\[Delta]*(-9 + 4*\[Nu]) + 
        ((I/576)*\[Delta]*(-226 + 39*\[Nu]))/E^((2*I)*\[Xi]) + 
        (I/4032)*E^((2*I)*\[Xi])*\[Delta]*(991 + 335*\[Nu]) - 
        ((I/960)*\[Delta]*(9343 + 4583*\[Nu]))/E^((6*I)*\[Xi]) + 
        ((I/5760)*\[Delta]*(26867 + 17164*\[Nu]))/E^((4*I)*\[Xi]) - 
        (I/6720)*E^((6*I)*\[Xi])*\[Delta]*(-239126 + 24743*\[Nu]) + 
        (I/40320)*E^((4*I)*\[Xi])*\[Delta]*(-752819 + 97476*\[Nu]))) + 
    SO^3*x^3*\[Epsilon]^6*
     ((-I/4)*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
          (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
          \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*
         \[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3) + 
      e*((((-5*I)/4)*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 2*\[Kappa]S*
                \[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
           (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3))/
         E^(I*\[Xi]) - ((5*I)/4)*E^(I*\[Xi])*
         (\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
            (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
            \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^
                2))*\[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + 
          \[Delta]*\[Chi]S^3)) + 
      e^3*((((-77*I)/32)*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 2*\[Kappa]S*
                \[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
           (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3))/
         E^(I*\[Xi]) - ((77*I)/32)*E^(I*\[Xi])*
         (\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
            (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
            \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^
                2))*\[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + 
          \[Delta]*\[Chi]S^3) - (((139*I)/32)*
          (\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
             \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*
                \[Chi]A^2))*\[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + 
           \[Delta]*\[Chi]S^3))/E^((3*I)*\[Xi]) - ((139*I)/32)*
         E^((3*I)*\[Xi])*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
            2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
          (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 
              2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
          (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3)) + 
      e^4*(((-65*I)/32)*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
            2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
          (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 
              2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
          (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3) - 
        (((169*I)/48)*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 2*\[Kappa]S*
                \[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
           (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3))/
         E^((2*I)*\[Xi]) - ((169*I)/48)*E^((2*I)*\[Xi])*
         (\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
            (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
            \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^
                2))*\[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + 
          \[Delta]*\[Chi]S^3) - (((1361*I)/192)*
          (\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
             \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*
                \[Chi]A^2))*\[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + 
           \[Delta]*\[Chi]S^3))/E^((4*I)*\[Xi]) - ((1361*I)/192)*
         E^((4*I)*\[Xi])*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
            2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
          (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 
              2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
          (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3)) + 
      e^5*((((-2977*I)/768)*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 2*\[Kappa]S*
                \[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
           (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3))/
         E^(I*\[Xi]) - ((2977*I)/768)*E^(I*\[Xi])*
         (\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
            (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
            \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^
                2))*\[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + 
          \[Delta]*\[Chi]S^3) - (((2287*I)/512)*
          (\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
             \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*
                \[Chi]A^2))*\[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + 
           \[Delta]*\[Chi]S^3))/E^((3*I)*\[Xi]) - ((2287*I)/512)*
         E^((3*I)*\[Xi])*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
            2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
          (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 
              2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
          (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3) - 
        (((17137*I)/1536)*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 2*\[Kappa]S*
                \[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
           (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3))/
         E^((5*I)*\[Xi]) - ((17137*I)/1536)*E^((5*I)*\[Xi])*
         (\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
            (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
            \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^
                2))*\[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + 
          \[Delta]*\[Chi]S^3)) + 
      e^6*(((-105*I)/32)*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
            2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
          (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 
              2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
          (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3) - 
        (((2023*I)/384)*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 2*\[Kappa]S*
                \[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
           (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3))/
         E^((2*I)*\[Xi]) - ((2023*I)/384)*E^((2*I)*\[Xi])*
         (\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
            (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
            \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^
                2))*\[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + 
          \[Delta]*\[Chi]S^3) - (((949*I)/192)*
          (\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
             \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*
                \[Chi]A^2))*\[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + 
           \[Delta]*\[Chi]S^3))/E^((4*I)*\[Xi]) - ((949*I)/192)*
         E^((4*I)*\[Xi])*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
            2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
          (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 
              2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
          (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3) - 
        (((2195*I)/128)*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 2*\[Kappa]S*
                \[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
           (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3))/
         E^((6*I)*\[Xi]) - ((2195*I)/128)*E^((6*I)*\[Xi])*
         (\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
            (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
            \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^
                2))*\[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + 
          \[Delta]*\[Chi]S^3)) + 
      e^2*((((-5*I)/2)*(\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           (\[Kappa]A - 4*\[Kappa]A*\[Nu] + \[Delta]*(\[Kappa]S - 2*\[Kappa]S*
                \[Nu] + (3 - 4*\[Nu])*\[Chi]A^2))*\[Chi]S + 
           (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + \[Delta]*\[Chi]S^3))/
         E^((2*I)*\[Xi]) - ((5*I)/2)*E^((2*I)*\[Xi])*
         (\[Chi]A*(\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
            (1 - 4*\[Nu])*\[Chi]A^2) + (\[Kappa]A - 4*\[Kappa]A*\[Nu] + 
            \[Delta]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + (3 - 4*\[Nu])*\[Chi]A^
                2))*\[Chi]S + (3 - 8*\[Nu])*\[Chi]A*\[Chi]S^2 + 
          \[Delta]*\[Chi]S^3) - I*(\[Delta]*\[Kappa]A*\[Chi]A + \[Chi]A^3 - 
          4*\[Nu]*\[Chi]A^3 + \[Kappa]S*(\[Chi]A - 2*\[Nu]*\[Chi]A) + 
          \[Kappa]A*\[Chi]S - 4*\[Kappa]A*\[Nu]*\[Chi]S + 
          3*\[Chi]A*\[Chi]S^2 - 8*\[Nu]*\[Chi]A*\[Chi]S^2 + 
          \[Delta]*\[Chi]S*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
            (3 - 4*\[Nu])*\[Chi]A^2 + \[Chi]S^2)))) + 
    SO^2*x^(5/2)*\[Epsilon]^5*((-I/6)*(\[Delta]*\[Kappa]S*(-1 + 6*\[Nu]) + 
        \[Kappa]A*(-1 + 8*\[Nu]) + 6*\[Delta]*(-1 + 2*\[Nu])*\[Chi]A^2 + 
        4*(-3 + 10*\[Nu])*\[Chi]A*\[Chi]S + 2*\[Delta]*(-3 + 2*\[Nu])*
         \[Chi]S^2) + e*((-I/6)*E^(I*\[Xi])*(\[Delta]*\[Kappa]S*
           (-5 + 14*\[Nu]) + \[Kappa]A*(-5 + 24*\[Nu]) + 
          \[Delta]*(-25 + 28*\[Nu])*\[Chi]A^2 + 2*(-25 + 69*\[Nu])*\[Chi]A*
           \[Chi]S + 5*\[Delta]*(-5 + 2*\[Nu])*\[Chi]S^2) - 
        ((I/6)*(\[Delta]*\[Kappa]S*(-1 + 14*\[Nu]) + \[Kappa]A*
            (-1 + 16*\[Nu]) + 7*\[Delta]*(-3 + 4*\[Nu])*\[Chi]A^2 + 
           2*(-21 + 61*\[Nu])*\[Chi]A*\[Chi]S + \[Delta]*(-21 + 10*\[Nu])*
            \[Chi]S^2))/E^(I*\[Xi])) + 
      e^5*((-I/2304)*E^((5*I)*\[Xi])*(-19497*\[Delta]*\[Kappa]S + 
          \[Kappa]A*(-19497 + 77656*\[Nu]) + 2*(-82477 + 215185*\[Nu])*
           \[Chi]A*\[Chi]S + \[Delta]*(38662*\[Nu]*(\[Kappa]S + 
              2*\[Chi]A^2) + 23138*\[Nu]*\[Chi]S^2 - 
            82477*(\[Chi]A^2 + \[Chi]S^2))) - 
        ((I/2304)*(2443*\[Delta]*\[Kappa]S + \[Kappa]A*(2443 + 33776*\[Nu]) + 
           2*(-60537 + 171305*\[Nu])*\[Chi]A*\[Chi]S + 
           \[Delta]*(38662*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
             23138*\[Nu]*\[Chi]S^2 - 60537*(\[Chi]A^2 + \[Chi]S^2))))/
         E^((5*I)*\[Xi]) - (I/1152)*E^(I*\[Xi])*(-1301*\[Delta]*\[Kappa]S + 
          \[Kappa]A*(-1301 + 8600*\[Nu]) + 2*(-8937 + 25709*\[Nu])*\[Chi]A*
           \[Chi]S + \[Delta]*(5998*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
            3674*\[Nu]*\[Chi]S^2 - 8937*(\[Chi]A^2 + \[Chi]S^2))) - 
        ((I/1152)*(-1041*\[Delta]*\[Kappa]S + \[Kappa]A*
            (-1041 + 8080*\[Nu]) + 2*(-8677 + 25189*\[Nu])*\[Chi]A*\[Chi]S + 
           \[Delta]*(5998*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
             3674*\[Nu]*\[Chi]S^2 - 8677*(\[Chi]A^2 + \[Chi]S^2))))/
         E^(I*\[Xi]) - ((I/768)*(-895*\[Delta]*\[Kappa]S + 
           \[Kappa]A*(-895 + 5288*\[Nu]) + 2*(-5723 + 16247*\[Nu])*\[Chi]A*
            \[Chi]S + \[Delta]*(3498*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
             2606*\[Nu]*\[Chi]S^2 - 5723*(\[Chi]A^2 + \[Chi]S^2))))/
         E^((3*I)*\[Xi]) - (I/768)*E^((3*I)*\[Xi])*
         (\[Kappa]A*(-595 + 4688*\[Nu]) + 2*(-5423 + 15647*\[Nu])*\[Chi]A*
           \[Chi]S + \[Delta]*(-595*\[Kappa]S + 3498*\[Kappa]S*\[Nu] + 
            6996*\[Nu]*\[Chi]A^2 + 2606*\[Nu]*\[Chi]S^2 - 
            5423*(\[Chi]A^2 + \[Chi]S^2)))) + 
      e^3*((-I/48)*E^(I*\[Xi])*(-41*\[Delta]*\[Kappa]S + 
          \[Kappa]A*(-41 + 256*\[Nu]) + 2*(-269 + 769*\[Nu])*\[Chi]A*
           \[Chi]S - 269*\[Delta]*(\[Chi]A^2 + \[Chi]S^2) + 
          6*\[Delta]*\[Nu]*(29*\[Kappa]S + 58*\[Chi]A^2 + 19*\[Chi]S^2)) - 
        ((I/48)*(-29*\[Delta]*\[Kappa]S + 29*\[Kappa]A*(-1 + 8*\[Nu]) + 
           2*(-257 + 745*\[Nu])*\[Chi]A*\[Chi]S - 257*\[Delta]*
            (\[Chi]A^2 + \[Chi]S^2) + 6*\[Delta]*\[Nu]*(29*\[Kappa]S + 
             58*\[Chi]A^2 + 19*\[Chi]S^2)))/E^(I*\[Xi]) - 
        (I/48)*E^((3*I)*\[Xi])*(\[Kappa]A*(-151 + 632*\[Nu]) + 
          2*(-675 + 1787*\[Nu])*\[Chi]A*\[Chi]S + 
          \[Delta]*(-151*\[Kappa]S + 330*\[Kappa]S*\[Nu] + 
            660*\[Nu]*\[Chi]A^2 + 214*\[Nu]*\[Chi]S^2 - 
            675*(\[Chi]A^2 + \[Chi]S^2))) - 
        ((I/48)*(5*\[Delta]*\[Kappa]S + 5*\[Kappa]A*(1 + 64*\[Nu]) + 
           2*(-519 + 1475*\[Nu])*\[Chi]A*\[Chi]S + \[Delta]*
            (330*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 214*\[Nu]*\[Chi]S^2 - 
             519*(\[Chi]A^2 + \[Chi]S^2))))/E^((3*I)*\[Xi])) + 
      e^4*(((-2*I)/3)*(\[Delta]*\[Kappa]S*(-1 + 6*\[Nu]) + 
          \[Kappa]A*(-1 + 8*\[Nu]) + \[Delta]*(-7 + 12*\[Nu])*\[Chi]A^2 + 
          (-14 + 43*\[Nu])*\[Chi]A*\[Chi]S + \[Delta]*(-7 + 3*\[Nu])*
           \[Chi]S^2) - (I/288)*E^((4*I)*\[Xi])*(-1517*\[Delta]*\[Kappa]S + 
          \[Kappa]A*(-1517 + 6168*\[Nu]) + 4*(-3287 + 8628*\[Nu])*\[Chi]A*
           \[Chi]S + 2*\[Delta]*(1567*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
            974*\[Nu]*\[Chi]S^2 - 3287*(\[Chi]A^2 + \[Chi]S^2))) - 
        ((I/288)*(131*\[Delta]*\[Kappa]S + \[Kappa]A*(131 + 2872*\[Nu]) + 
           4*(-2463 + 6980*\[Nu])*\[Chi]A*\[Chi]S + 2*\[Delta]*
            (1567*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 974*\[Nu]*\[Chi]S^2 - 
             2463*(\[Chi]A^2 + \[Chi]S^2))))/E^((4*I)*\[Xi]) - 
        (I/36)*E^((2*I)*\[Xi])*(-37*\[Delta]*\[Kappa]S + 
          \[Kappa]A*(-37 + 232*\[Nu]) + 85*(-6 + 17*\[Nu])*\[Chi]A*\[Chi]S + 
          \[Delta]*(158*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
            109*\[Nu]*\[Chi]S^2 - 255*(\[Chi]A^2 + \[Chi]S^2))) - 
        ((I/36)*(-29*\[Delta]*\[Kappa]S + \[Kappa]A*(-29 + 216*\[Nu]) + 
           (-494 + 1413*\[Nu])*\[Chi]A*\[Chi]S + \[Delta]*
            (158*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 109*\[Nu]*\[Chi]S^2 - 
             247*(\[Chi]A^2 + \[Chi]S^2))))/E^((2*I)*\[Xi])) + 
      e^6*((-I/288)*E^((2*I)*\[Xi])*(-425*\[Delta]*\[Kappa]S + 
          \[Kappa]A*(-425 + 2638*\[Nu]) + (-5668 + 16049*\[Nu])*\[Chi]A*
           \[Chi]S - 2834*\[Delta]*(\[Chi]A^2 + \[Chi]S^2) + 
          3*\[Delta]*\[Nu]*(596*(\[Kappa]S + 2*\[Chi]A^2) + 379*\[Chi]S^2)) - 
        (I/480)*E^((6*I)*\[Xi])*(\[Kappa]A*(-6350 + 24898*\[Nu]) + 
          (-52680 + 136811*\[Nu])*\[Chi]A*\[Chi]S + 
          \[Delta]*(-6350*\[Kappa]S + 12198*\[Kappa]S*\[Nu] + 
            24396*\[Nu]*\[Chi]A^2 + 7055*\[Nu]*\[Chi]S^2 - 
            26340*(\[Chi]A^2 + \[Chi]S^2))) - 
        ((I/480)*(38*\[Kappa]A*(26 + 269*\[Nu]) + (-38004 + 107459*\[Nu])*
            \[Chi]A*\[Chi]S + \[Delta]*(38*(26*\[Kappa]S + 321*\[Kappa]S*
                \[Nu] + 642*\[Nu]*\[Chi]A^2) + 7055*\[Nu]*\[Chi]S^2 - 
             19002*(\[Chi]A^2 + \[Chi]S^2))))/E^((6*I)*\[Xi]) - 
        ((I/2880)*(-5275*\[Delta]*\[Kappa]S + \[Kappa]A*
            (-5275 + 20552*\[Nu]) + 4*(-8975 + 25261*\[Nu])*\[Chi]A*\[Chi]S + 
           2*\[Delta]*(5001*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
             4620*\[Nu]*\[Chi]S^2 - 8975*(\[Chi]A^2 + \[Chi]S^2))))/
         E^((4*I)*\[Xi]) - (I/2880)*E^((4*I)*\[Xi])*(917*\[Delta]*\[Kappa]S + 
          \[Kappa]A*(917 + 8168*\[Nu]) + 4*(-5879 + 19069*\[Nu])*\[Chi]A*
           \[Chi]S + 2*\[Delta]*(5001*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
            4620*\[Nu]*\[Chi]S^2 - 5879*(\[Chi]A^2 + \[Chi]S^2))) - 
        ((I/288)*(\[Kappa]A*(-299 + 2386*\[Nu]) + (-5416 + 15545*\[Nu])*
            \[Chi]A*\[Chi]S + \[Delta]*(-299*\[Kappa]S + 1788*\[Kappa]S*
              \[Nu] + 3576*\[Nu]*\[Chi]A^2 + 1137*\[Nu]*\[Chi]S^2 - 
             2708*(\[Chi]A^2 + \[Chi]S^2))))/E^((2*I)*\[Xi]) - 
        (I/12)*(11*\[Kappa]A*(-1 + 8*\[Nu]) + 4*(-39 + 119*\[Nu])*\[Chi]A*
           \[Chi]S + \[Delta]*(-11*\[Kappa]S + 66*\[Kappa]S*\[Nu] + 
            132*\[Nu]*\[Chi]A^2 + 32*\[Nu]*\[Chi]S^2 - 
            78*(\[Chi]A^2 + \[Chi]S^2)))) + 
      e^2*((-I/12)*(5*\[Kappa]A*(-1 + 8*\[Nu]) + 4*(-17 + 53*\[Nu])*\[Chi]A*
           \[Chi]S + \[Delta]*(5*\[Kappa]S*(-1 + 6*\[Nu]) + 
            (-34 + 60*\[Nu])*\[Chi]A^2 + 2*(-17 + 8*\[Nu])*\[Chi]S^2)) - 
        (I/12)*E^((2*I)*\[Xi])*(\[Kappa]A*(-21 + 92*\[Nu]) + 
          2*(-98 + 263*\[Nu])*\[Chi]A*\[Chi]S + \[Delta]*(-21*\[Kappa]S + 
            50*\[Kappa]S*\[Nu] + 100*\[Nu]*\[Chi]A^2 + 34*\[Nu]*\[Chi]S^2 - 
            98*(\[Chi]A^2 + \[Chi]S^2))) - 
        ((I/12)*(\[Kappa]A*(-1 + 52*\[Nu]) + 2*(-78 + 223*\[Nu])*\[Chi]A*
            \[Chi]S + \[Delta]*(-\[Kappa]S + 50*\[Kappa]S*\[Nu] + 
             100*\[Nu]*\[Chi]A^2 + 34*\[Nu]*\[Chi]S^2 - 
             78*(\[Chi]A^2 + \[Chi]S^2))))/E^((2*I)*\[Xi]))) + 
    x^2*\[Epsilon]^4*
     (e^3*((E^(I*\[Xi])*\[Delta]*(19 + (34*I)*Pi - 68*Log[2]))/48 + 
        (E^((3*I)*\[Xi])*\[Delta]*(21 + (86*I)*Pi + 116*Log[2]))/144 + 
        (\[Delta]*(17 + (34*I)*Pi + 356*Log[2] - 324*Log[3]))/
         (48*E^(I*\[Xi])) + (\[Delta]*(255 + (510*I)*Pi - 3028*Log[2] + 
           972*Log[3]))/(144*E^((3*I)*\[Xi]))) + 
      e^4*((E^((2*I)*\[Xi])*\[Delta]*(36 + (49*I)*Pi - 128*Log[2]))/72 + 
        (\[Delta]*(42 + (83*I)*Pi + 808*Log[2] - 729*Log[3]))/96 + 
        (E^((4*I)*\[Xi])*\[Delta]*(50 + (525*I)*Pi + 328*Log[2] + 
           243*Log[3]))/576 + (\[Delta]*(14 + (28*I)*Pi - 2232*Log[2] + 
           1215*Log[3]))/(72*E^((2*I)*\[Xi])) + 
        (\[Delta]*(6*(283 + (566*I)*Pi + 3884*Log[2] + 81*Log[3]) - 
           15625*Log[5]))/(576*E^((4*I)*\[Xi]))) + 
      e^5*((E^((3*I)*\[Xi])*\[Delta]*(521 + (414*I)*Pi - (4388*Log[2])/3 - 
           324*Log[3]))/768 + (E^((5*I)*\[Xi])*\[Delta]*
          (-415 + (16382*I)*Pi + 17604*Log[2] + 4860*Log[3]))/11520 + 
        (\[Delta]*(559 + (1110*I)*Pi - 39892*Log[2] + 21384*Log[3]))/
         (1152*E^(I*\[Xi])) + E^(I*\[Xi])*
         ((\[Delta]*(549 + (1082*I)*Pi + 11668*Log[2]))/1152 - 
          9*\[Delta]*Log[3]) + (\[Delta]*(-759 - (1518*I)*Pi + 
           176100*Log[2] - 23328*Log[3] - 62500*Log[5]))/
         (2304*E^((3*I)*\[Xi])) + (\[Delta]*(54435 + (108870*I)*Pi - 
           299012*Log[2] - 553392*Log[3] + 312500*Log[5]))/
         (11520*E^((5*I)*\[Xi]))) + 
      e^6*((E^((2*I)*\[Xi])*\[Delta]*(1156 + (2353*I)*Pi + 29328*Log[2] - 
           25434*Log[3]))/2304 + (E^((4*I)*\[Xi])*\[Delta]*
          (2834 + (509*I)*Pi - 9272*Log[2] - 1701*Log[3]))/2880 + 
        (\[Delta]*((887*I)*Pi + 9*(50 - 4024*Log[2] + 2187*Log[3])))/864 + 
        (\[Delta]*(2584 + (5141*I)*Pi + 407904*Log[2] - 70713*Log[3] - 
           140625*Log[5]))/(4608*E^((2*I)*\[Xi])) + 
        (E^((6*I)*\[Xi])*\[Delta]*(-18072 + (153793*I)*Pi + 46496*Log[2] + 
           32805*Log[3] + 78125*Log[5]))/69120 + 
        (\[Delta]*(-4382 - (8764*I)*Pi - 301816*Log[2] - 144018*Log[3] + 
           234375*Log[5]))/(2880*E^((4*I)*\[Xi])) + 
        (\[Delta]*(510264 + (1020528*I)*Pi + 5576416*Log[2] + 
           3392037*Log[3] + 234375*Log[5] - 5764801*Log[7]))/
         (69120*E^((6*I)*\[Xi]))) + (\[Delta]*(1 + (2*I)*Pi + Log[16]))/6 + 
      e*((\[Delta]*(3 + (6*I)*Pi - 4*Log[2]))/(6*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*\[Delta]*(1 + (2*I)*Pi + Log[16]))/6) + 
      e^2*((\[Delta]*(1 + (2*I)*Pi - 4*Log[2]))/3 + 
        (\[Delta]*(1 + (2*I)*Pi + (20*Log[2])/3 - (27*Log[3])/4))/
         E^((2*I)*\[Xi]) + (E^((2*I)*\[Xi])*\[Delta]*(2 + (5*I)*Pi + 
           Log[256]))/12)) + 
    SO*(x*\[Epsilon]^2*((-I/2)*(\[Chi]A + \[Delta]*\[Chi]S) + 
        e*(((-I/2)*(\[Chi]A + \[Delta]*\[Chi]S))/E^(I*\[Xi]) - 
          (I/2)*E^(I*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)) + 
        e^2*((-I/4)*(\[Chi]A + \[Delta]*\[Chi]S) - 
          (((5*I)/8)*(\[Chi]A + \[Delta]*\[Chi]S))/E^((2*I)*\[Xi]) - 
          ((5*I)/8)*E^((2*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)) + 
        e^3*((((-3*I)/16)*(\[Chi]A + \[Delta]*\[Chi]S))/E^(I*\[Xi]) - 
          ((3*I)/16)*E^(I*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S) - 
          (((13*I)/16)*(\[Chi]A + \[Delta]*\[Chi]S))/E^((3*I)*\[Xi]) - 
          ((13*I)/16)*E^((3*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)) + 
        e^4*(((-3*I)/16)*(\[Chi]A + \[Delta]*\[Chi]S) - 
          ((I/12)*(\[Chi]A + \[Delta]*\[Chi]S))/E^((2*I)*\[Xi]) - 
          (I/12)*E^((2*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S) - 
          (((103*I)/96)*(\[Chi]A + \[Delta]*\[Chi]S))/E^((4*I)*\[Xi]) - 
          ((103*I)/96)*E^((4*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)) + 
        e^5*((((-65*I)/384)*(\[Chi]A + \[Delta]*\[Chi]S))/E^(I*\[Xi]) - 
          ((65*I)/384)*E^(I*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S) + 
          (((25*I)/256)*(\[Chi]A + \[Delta]*\[Chi]S))/E^((3*I)*\[Xi]) + 
          ((25*I)/256)*E^((3*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S) - 
          (((1097*I)/768)*(\[Chi]A + \[Delta]*\[Chi]S))/E^((5*I)*\[Xi]) - 
          ((1097*I)/768)*E^((5*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)) + 
        e^6*(((-5*I)/32)*(\[Chi]A + \[Delta]*\[Chi]S) - 
          (((21*I)/128)*(\[Chi]A + \[Delta]*\[Chi]S))/E^((2*I)*\[Xi]) - 
          ((21*I)/128)*E^((2*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S) + 
          (((129*I)/320)*(\[Chi]A + \[Delta]*\[Chi]S))/E^((4*I)*\[Xi]) + 
          ((129*I)/320)*E^((4*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S) - 
          (((1223*I)/640)*(\[Chi]A + \[Delta]*\[Chi]S))/E^((6*I)*\[Xi]) - 
          ((1223*I)/640)*E^((6*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S))) + 
      x^2*\[Epsilon]^4*((I/42)*((-7 + 205*\[Nu])*\[Chi]A + 
          \[Delta]*(-7 + 33*\[Nu])*\[Chi]S) + 
        e*((I/84)*E^(I*\[Xi])*((-364 + 991*\[Nu])*\[Chi]A + 
            13*\[Delta]*(-28 + 11*\[Nu])*\[Chi]S) + 
          ((I/84)*((-70 + 825*\[Nu])*\[Chi]A + \[Delta]*(-70 + 169*\[Nu])*
              \[Chi]S))/E^(I*\[Xi])) + 
        e^2*((I/24)*(9*(-3 + 26*\[Nu])*\[Chi]A + \[Delta]*(-27 + 50*\[Nu])*
             \[Chi]S) + (I/336)*E^((2*I)*\[Xi])*((-3353 + 7180*\[Nu])*
             \[Chi]A + \[Delta]*(-3353 + 996*\[Nu])*\[Chi]S) + 
          ((I/336)*((-413 + 5520*\[Nu])*\[Chi]A + \[Delta]*(-413 + 1256*
                \[Nu])*\[Chi]S))/E^((2*I)*\[Xi])) + 
        e^3*((I/96)*E^(I*\[Xi])*((-176 + 1131*\[Nu])*\[Chi]A + 
            \[Delta]*(-176 + 235*\[Nu])*\[Chi]S) + 
          ((I/672)*((-1526 + 8083*\[Nu])*\[Chi]A + \[Delta]*(-1526 + 1619*
                \[Nu])*\[Chi]S))/E^(I*\[Xi]) + (I/672)*E^((3*I)*\[Xi])*
           ((-12530 + 23909*\[Nu])*\[Chi]A + \[Delta]*(-12530 + 3237*\[Nu])*
             \[Chi]S) + ((I/672)*((-1064 + 17435*\[Nu])*\[Chi]A + 
             \[Delta]*(-1064 + 4251*\[Nu])*\[Chi]S))/E^((3*I)*\[Xi])) + 
        e^4*((I/144)*E^((2*I)*\[Xi])*((10 + 1533*\[Nu])*\[Chi]A + 
            5*\[Delta]*(2 + 73*\[Nu])*\[Chi]S) + 
          (I/112)*((-182 + 1433*\[Nu])*\[Chi]A + \[Delta]*(-182 + 317*\[Nu])*
             \[Chi]S) + ((I/1008)*((-3164 + 12557*\[Nu])*\[Chi]A + 
             \[Delta]*(-3164 + 2269*\[Nu])*\[Chi]S))/E^((2*I)*\[Xi]) + 
          (I/2016)*E^((4*I)*\[Xi])*((-64358 + 114281*\[Nu])*\[Chi]A + 
            \[Delta]*(-64358 + 15205*\[Nu])*\[Chi]S) + 
          ((I/2016)*((-3794 + 80085*\[Nu])*\[Chi]A + \[Delta]*
              (-3794 + 20561*\[Nu])*\[Chi]S))/E^((4*I)*\[Xi])) + 
        e^6*(((5*I)/1344)*((-539 + 4094*\[Nu])*\[Chi]A + 
            \[Delta]*(-539 + 918*\[Nu])*\[Chi]S) + (I/5760)*E^((4*I)*\[Xi])*
           ((116855 - 75630*\[Nu])*\[Chi]A + \[Delta]*(116855 + 1946*\[Nu])*
             \[Chi]S) + (I/16128)*E^((2*I)*\[Xi])*((-58793 + 272248*\[Nu])*
             \[Chi]A + \[Delta]*(-58793 + 55520*\[Nu])*\[Chi]S) + 
          ((I/16128)*((-48797 + 266604*\[Nu])*\[Chi]A + 
             \[Delta]*(-48797 + 56404*\[Nu])*\[Chi]S))/E^((2*I)*\[Xi]) + 
          ((I/40320)*((-242767 + 69518*\[Nu])*\[Chi]A - 
             \[Delta]*(242767 + 80186*\[Nu])*\[Chi]S))/E^((4*I)*\[Xi]) + 
          ((I/3840)*(25*(-307 + 13564*\[Nu])*\[Chi]A + \[Delta]*
              (-7675 + 93876*\[Nu])*\[Chi]S))/E^((6*I)*\[Xi]) + 
          (I/26880)*E^((6*I)*\[Xi])*((-2211097 + 3591808*\[Nu])*\[Chi]A + 
            \[Delta]*(-2211097 + 466344*\[Nu])*\[Chi]S)) + 
        e^5*(((I/10752)*(3*(-15400 + 35383*\[Nu])*\[Chi]A + 
             55*\[Delta]*(-840 + 211*\[Nu])*\[Chi]S))/E^((3*I)*\[Xi]) + 
          (I/10752)*E^((3*I)*\[Xi])*((67578 + 41907*\[Nu])*\[Chi]A + 
            \[Delta]*(67578 + 21667*\[Nu])*\[Chi]S) + (I/16128)*E^(I*\[Xi])*
           ((-44072 + 241299*\[Nu])*\[Chi]A + \[Delta]*(-44072 + 50531*\[Nu])*
             \[Chi]S) + ((I/16128)*((-42602 + 240469*\[Nu])*\[Chi]A + 
             \[Delta]*(-42602 + 50661*\[Nu])*\[Chi]S))/E^(I*\[Xi]) + 
          (I/32256)*E^((5*I)*\[Xi])*((-1679048 + 2834155*\[Nu])*\[Chi]A + 
            \[Delta]*(-1679048 + 372027*\[Nu])*\[Chi]S) + 
          ((I/32256)*((-66458 + 1923645*\[Nu])*\[Chi]A + 
             \[Delta]*(-66458 + 514637*\[Nu])*\[Chi]S))/E^((5*I)*\[Xi]))) + 
      x^3*\[Epsilon]^6*((-I/1512)*((1248 + 7*\[Nu]*(-1193 + 780*\[Nu]))*
           \[Chi]A + \[Delta]*(1248 - \[Nu]*(7027 + 716*\[Nu]))*\[Chi]S) + 
        e*((-I/3024)*E^(I*\[Xi])*((42696 + \[Nu]*(-189343 + 53019*\[Nu]))*
             \[Chi]A + \[Delta]*(42696 + \[Nu]*(-76991 + 731*\[Nu]))*
             \[Chi]S) - ((I/3024)*((28950 + \[Nu]*(-124517 + 32469*\[Nu]))*
              \[Chi]A + \[Delta]*(28950 - \[Nu]*(51829 + 4811*\[Nu]))*
              \[Chi]S))/E^(I*\[Xi])) + 
        e^2*((-I/3024)*(5*(4860 + \[Nu]*(-38519 + 18210*\[Nu]))*\[Chi]A + 
            \[Delta]*(24300 + \[Nu]*(-85007 + 3722*\[Nu]))*\[Chi]S) - 
          (I/6048)*E^((2*I)*\[Xi])*((249798 + \[Nu]*(-974975 + 264384*\[Nu]))*
             \[Chi]A + \[Delta]*(249798 + \[Nu]*(-346987 + 11680*\[Nu]))*
             \[Chi]S) - ((I/6048)*((66066 + \[Nu]*(-411503 + 104172*\[Nu]))*
              \[Chi]A + \[Delta]*(66066 - \[Nu]*(183451 + 27764*\[Nu]))*
              \[Chi]S))/E^((2*I)*\[Xi])) + 
        e^3*(((-I/8064)*((228112 + \[Nu]*(-1119767 + 420643*\[Nu]))*\[Chi]A + 
             \[Delta]*(228112 + \[Nu]*(-432999 + 16643*\[Nu]))*\[Chi]S))/
           E^(I*\[Xi]) - (I/24192)*E^(I*\[Xi])*
           ((595698 + \[Nu]*(-3488023 + 1351479*\[Nu]))*\[Chi]A + 
            \[Delta]*(595698 + \[Nu]*(-1492151 + 79031*\[Nu]))*\[Chi]S) - 
          ((I/8064)*((43162 + \[Nu]*(-760019 + 201939*\[Nu]))*\[Chi]A + 
             \[Delta]*(43162 - \[Nu]*(403443 + 85421*\[Nu]))*\[Chi]S))/
           E^((3*I)*\[Xi]) - (I/24192)*E^((3*I)*\[Xi])*
           ((2290512 + \[Nu]*(-8264491 + 2278119*\[Nu]))*\[Chi]A + 
            \[Delta]*(2290512 + \[Nu]*(-2723387 + 138311*\[Nu]))*\[Chi]S)) + 
        e^4*((-I/24192)*((423801 + 2*\[Nu]*(-1796779 + 854634*\[Nu]))*
             \[Chi]A + \[Delta]*(423801 + 2*\[Nu]*(-741659 + 54554*\[Nu]))*
             \[Chi]S) - ((I/36288)*((1655415 - 7532846*\[Nu] + 2747934*
                \[Nu]^2)*\[Chi]A + \[Delta]*(1655415 + 2*\[Nu]*(-1329455 + 
                 77471*\[Nu]))*\[Chi]S))/E^((2*I)*\[Xi]) - 
          (I/36288)*E^((2*I)*\[Xi])*((1017327 - 6949330*\[Nu] + 
              2699202*\[Nu]^2)*\[Chi]A + \[Delta]*(1017327 + 
              2*\[Nu]*(-1568089 + 87217*\[Nu]))*\[Chi]S) - 
          (I/145152)*E^((4*I)*\[Xi])*((27838341 + 46*\[Nu]*(-2076415 + 
                587142*\[Nu]))*\[Chi]A + \[Delta]*(27838341 + 
              2*\[Nu]*(-14964961 + 960298*\[Nu]))*\[Chi]S) + 
          ((I/145152)*((1687827 + 2*(8567173 - 2513658*\[Nu])*\[Nu])*
              \[Chi]A + \[Delta]*(1687827 + 2*\[Nu]*(5885261 + 1569382*
                  \[Nu]))*\[Chi]S))/E^((4*I)*\[Xi])) + 
        e^6*(((-5*I)/16128)*((92239 - 817726*\[Nu] + 391296*\[Nu]^2)*
             \[Chi]A + \[Delta]*(92239 - 328606*\[Nu] + 28448*\[Nu]^2)*
             \[Chi]S) - (I/1451520)*E^((4*I)*\[Xi])*
           ((-109081707 + 20064658*\[Nu] + 30469416*\[Nu]^2)*\[Chi]A - 
            \[Delta]*(109081707 + 94903774*\[Nu] + 711800*\[Nu]^2)*\[Chi]S) - 
          ((I/580608)*((36374415 + 2*\[Nu]*(-95728547 + 37351842*\[Nu]))*
              \[Chi]A + \[Delta]*(36374415 - 73971398*\[Nu] + 3901940*
                \[Nu]^2)*\[Chi]S))/E^((2*I)*\[Xi]) - (I/967680)*
           E^((6*I)*\[Xi])*((620804751 - 2012219990*\[Nu] + 598058556*\[Nu]^
                2)*\[Chi]A + \[Delta]*(620804751 - 593629270*\[Nu] + 
              49822636*\[Nu]^2)*\[Chi]S) + 
          ((I/967680)*((111728577 + 2*(64258391 - 28839138*\[Nu])*\[Nu])*
              \[Chi]A + \[Delta]*(111728577 + 578*\[Nu]*(351943 + 
                 123110*\[Nu]))*\[Chi]S))/E^((6*I)*\[Xi]) - 
          (I/193536)*E^((2*I)*\[Xi])*((11696101 - 68432906*\[Nu] + 
              27241108*\[Nu]^2)*\[Chi]A + \[Delta]*(11696101 + 
              2*\[Nu]*(-13907637 + 976594*\[Nu]))*\[Chi]S) - 
          ((I/483840)*((57491327 - 201233010*\[Nu] + 71166152*\[Nu]^2)*
              \[Chi]A + \[Delta]*(57491327 + 2*\[Nu]*(-22526705 + 
                 5993716*\[Nu]))*\[Chi]S))/E^((4*I)*\[Xi])) + 
        e^5*(((-I/129024)*((9281544 + \[Nu]*(-37994347 + 13693495*\[Nu]))*
              \[Chi]A + \[Delta]*(9281544 + \[Nu]*(-11469227 + 1334391*
                  \[Nu]))*\[Chi]S))/E^((3*I)*\[Xi]) - 
          (I/387072)*E^((3*I)*\[Xi])*((2253114 + \[Nu]*(-66220891 + 
                28375755*\[Nu]))*\[Chi]A + \[Delta]*(2253114 + 
              \[Nu]*(-37101419 + 1773803*\[Nu]))*\[Chi]S) - 
          ((I/580608)*((26349270 + \[Nu]*(-145845541 + 58927509*\[Nu]))*
              \[Chi]A + \[Delta]*(26349270 + \[Nu]*(-57079253 + 3341237*
                  \[Nu]))*\[Chi]S))/E^(I*\[Xi]) - (I/580608)*E^(I*\[Xi])*
           ((24749400 + \[Nu]*(-150704383 + 61824891*\[Nu]))*\[Chi]A + 
            \[Delta]*(24749400 + \[Nu]*(-61746431 + 4200443*\[Nu]))*
             \[Chi]S) - (I/1161216)*E^((5*I)*\[Xi])*
           ((418011720 + \[Nu]*(-1386769855 + 402400539*\[Nu]))*\[Chi]A + 
            \[Delta]*(418011720 + \[Nu]*(-419607935 + 31466459*\[Nu]))*
             \[Chi]S) + ((I/1161216)*((55481370 + (155801125 - 53640309*
                  \[Nu])*\[Nu])*\[Chi]A + \[Delta]*(55481370 + \[Nu]*
                (151827605 + 47513131*\[Nu]))*\[Chi]S))/E^((5*I)*\[Xi]))) + 
      x^(5/2)*\[Epsilon]^5*(e^3*((E^(I*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)*
            (-23 - (42*I)*Pi + 52*Log[2]))/32 - (I/96)*E^((3*I)*\[Xi])*
           (\[Chi]A + \[Delta]*\[Chi]S)*(86*Pi - I*(21 + 116*Log[2])) + 
          ((\[Chi]A + \[Delta]*\[Chi]S)*(-255 - (510*I)*Pi + 3028*Log[2] - 
             972*Log[3]))/(96*E^((3*I)*\[Xi])) + 
          ((\[Chi]A + \[Delta]*\[Chi]S)*(-29 - (58*I)*Pi - 340*Log[2] + 
             324*Log[3]))/(32*E^(I*\[Xi]))) + 
        e^4*((E^((2*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)*(-21 - (32*I)*Pi + 
             52*Log[2]))/24 + ((\[Chi]A + \[Delta]*\[Chi]S)*
            (-25 - (50*I)*Pi + 996*Log[2] - 486*Log[3]))/
           (24*E^((2*I)*\[Xi])) + (E^((4*I)*\[Xi])*(\[Chi]A + 
             \[Delta]*\[Chi]S)*(-50 - (525*I)*Pi - 328*Log[2] - 243*Log[3]))/
           384 + ((\[Chi]A + \[Delta]*\[Chi]S)*(-64 - (127*I)*Pi - 
             768*Log[2] + 729*Log[3]))/64 + ((\[Chi]A + \[Delta]*\[Chi]S)*
            ((-3396*I)*Pi - 6*(283 + 3884*Log[2] + 81*Log[3]) + 
             15625*Log[5]))/(384*E^((4*I)*\[Xi]))) + 
        e^5*(((\[Chi]A + \[Delta]*\[Chi]S)*(-979 - (1950*I)*Pi + 
             35908*Log[2] - 17496*Log[3]))/(768*E^(I*\[Xi])) + 
          (E^((5*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)*(415 - (16382*I)*Pi - 
             17604*Log[2] - 4860*Log[3]))/7680 + 
          (E^((3*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)*
            (-1731 - (1930*I)*Pi + 3460*Log[2] + 972*Log[3]))/1536 + 
          (E^(I*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)*(-849 - (1634*I)*Pi - 
             11140*Log[2] + 10368*Log[3]))/768 + 
          ((\[Chi]A + \[Delta]*\[Chi]S)*(-54435 - (108870*I)*Pi + 
             299012*Log[2] + 553392*Log[3] - 312500*Log[5]))/
           (7680*E^((5*I)*\[Xi])) + ((\[Chi]A + \[Delta]*\[Chi]S)*
            (-1281 - (2562*I)*Pi - 151876*Log[2] + 15552*Log[3] + 
             62500*Log[5]))/(1536*E^((3*I)*\[Xi]))) + 
        e^6*(((\[Chi]A + \[Delta]*\[Chi]S)*((-3133*I)*Pi + 
             9*(-176 + 7296*Log[2] - 3645*Log[3])))/1152 + 
          (E^((4*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)*
            (-5918 - (3643*I)*Pi + 16904*Log[2] + 2187*Log[3]))/3840 + 
          (E^((2*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)*
            (-1876 - (3497*I)*Pi - 27856*Log[2] + 25434*Log[3]))/1536 + 
          ((\[Chi]A + \[Delta]*\[Chi]S)*(274 + (548*I)*Pi + 487112*Log[2] + 
             285606*Log[3] - 390625*Log[5]))/(3840*E^((4*I)*\[Xi])) + 
          (E^((6*I)*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)*
            (18072 - (153793*I)*Pi - 46496*Log[2] - 32805*Log[3] - 
             78125*Log[5]))/46080 + ((\[Chi]A + \[Delta]*\[Chi]S)*
            (-4760 - (9493*I)*Pi - 348000*Log[2] + 43497*Log[3] + 
             140625*Log[5]))/(3072*E^((2*I)*\[Xi])) + 
          ((\[Chi]A + \[Delta]*\[Chi]S)*(-510264 - (1020528*I)*Pi - 
             5576416*Log[2] - 3392037*Log[3] - 234375*Log[5] + 
             5764801*Log[7]))/(46080*E^((6*I)*\[Xi]))) - 
        (I/4)*(\[Chi]A + \[Delta]*\[Chi]S)*(2*Pi - I*(1 + Log[16])) + 
        e*(((\[Chi]A + \[Delta]*\[Chi]S)*(-3 - (6*I)*Pi + Log[16]))/
           (4*E^(I*\[Xi])) - (I/4)*E^(I*\[Xi])*(\[Chi]A + \[Delta]*\[Chi]S)*
           (2*Pi - I*(1 + Log[16]))) + 
        e^2*(((\[Chi]A + \[Delta]*\[Chi]S)*(-12 - (24*I)*Pi - 80*Log[2] + 
             81*Log[3]))/(8*E^((2*I)*\[Xi])) - (I/8)*E^((2*I)*\[Xi])*
           (\[Chi]A + \[Delta]*\[Chi]S)*(5*Pi - (2*I)*(1 + Log[16])) + 
          ((\[Chi]A + \[Delta]*\[Chi]S)*(-5 - (10*I)*Pi + Log[4096]))/8)))
 
H\[Psi]InstTail[2, 2] = 1 + e*(1/(4*E^(I*\[Xi])) + (5*E^(I*\[Xi]))/4) + 
    e^2*(-1/2 + 1/(4*E^((2*I)*\[Xi])) + (7*E^((2*I)*\[Xi]))/4) + 
    e^3*(-5/(32*E^(I*\[Xi])) - (33*E^(I*\[Xi]))/32 + 9/(32*E^((3*I)*\[Xi])) + 
      (77*E^((3*I)*\[Xi]))/32) + e^4*(-1/8 - 1/(6*E^((2*I)*\[Xi])) - 
      (11*E^((2*I)*\[Xi]))/6 + 1/(3*E^((4*I)*\[Xi])) + 
      (79*E^((4*I)*\[Xi]))/24) + e^5*(-47/(768*E^(I*\[Xi])) - 
      (11*E^(I*\[Xi]))/768 - 117/(512*E^((3*I)*\[Xi])) - 
      (1585*E^((3*I)*\[Xi]))/512 + 625/(1536*E^((5*I)*\[Xi])) + 
      (6901*E^((5*I)*\[Xi]))/1536) + e^6*(-1/16 - 1/(32*E^((2*I)*\[Xi])) + 
      E^((2*I)*\[Xi])/3 - 1/(3*E^((4*I)*\[Xi])) - (403*E^((4*I)*\[Xi]))/80 + 
      81/(160*E^((6*I)*\[Xi])) + (49*E^((6*I)*\[Xi]))/8) + 
    x*\[Epsilon]^2*((-107 + 55*\[Nu])/42 + 
      e*((E^(I*\[Xi])*(-31 + 35*\[Nu]))/24 + (-257 + 169*\[Nu])/
         (168*E^(I*\[Xi]))) + e^2*((-221 + 89*\[Nu])/84 + 
        (-289 + 181*\[Nu])/(168*E^((2*I)*\[Xi])) + 
        (E^((2*I)*\[Xi])*(-71 + 283*\[Nu]))/168) + 
      e^4*((-19/36 + (5*\[Nu])/9)/E^((2*I)*\[Xi]) + 
        (E^((4*I)*\[Xi])*(160 + 157*\[Nu]))/63 + (-379 + 215*\[Nu])/
         (144*E^((4*I)*\[Xi])) + (-549 + 233*\[Nu])/336 + 
        (E^((2*I)*\[Xi])*(-1367 + 328*\[Nu]))/252) + 
      e^3*((5*(-223 + 179*\[Nu]))/(1344*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*(-685 + 237*\[Nu]))/192 + (-2813 + 1677*\[Nu])/
         (1344*E^((3*I)*\[Xi])) + (E^((3*I)*\[Xi])*(1027 + 2729*\[Nu]))/
         1344) + e^6*((-877 + 377*\[Nu])/672 + 
        (5*E^((2*I)*\[Xi])*(-299 + 379*\[Nu]))/4032 + (-590 + 447*\[Nu])/
         (1008*E^((2*I)*\[Xi])) + (-1832 + 939*\[Nu])/(420*E^((6*I)*\[Xi])) + 
        (4873 + 2755*\[Nu])/(10080*E^((4*I)*\[Xi])) + 
        (E^((4*I)*\[Xi])*(-69071 + 7167*\[Nu]))/5040 + 
        (E^((6*I)*\[Xi])*(62233 + 25835*\[Nu]))/6720) + 
      e^5*((-2815 + 9411*\[Nu])/(21504*E^((3*I)*\[Xi])) + 
        (-20249 + 16465*\[Nu])/(32256*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*(-38305 + 21661*\[Nu]))/32256 + 
        (E^((3*I)*\[Xi])*(-183655 + 29423*\[Nu]))/21504 + 
        (-217369 + 117265*\[Nu])/(64512*E^((5*I)*\[Xi])) + 
        (E^((5*I)*\[Xi])*(337087 + 199165*\[Nu]))/64512)) + 
    x^2*\[Epsilon]^4*((-2173 + \[Nu]*(-7483 + 2047*\[Nu]))/1512 + 
      e*((E^(I*\[Xi])*(-17240 + \[Nu]*(-4965 + 2597*\[Nu])))/2016 + 
        (-34168 + \[Nu]*(-35131 + 2947*\[Nu]))/(6048*E^(I*\[Xi]))) + 
      e^2*((-34399 + \[Nu]*(-22897 + 829*\[Nu]))/3024 + 
        (E^((2*I)*\[Xi])*(-48099 + \[Nu]*(14467 + 1881*\[Nu])))/2016 + 
        (-77575 - \[Nu]*(42697 + 4523*\[Nu]))/(6048*E^((2*I)*\[Xi]))) + 
      e^4*((E^((4*I)*\[Xi])*(-1258321 + (775017 - 35033*\[Nu])*\[Nu]))/
         12096 + (-51565 - \[Nu]*(49063 + 605*\[Nu]))/4032 + 
        (-21593 + \[Nu]*(-59167 + 6779*\[Nu]))/(6048*E^((2*I)*\[Xi])) + 
        (E^((2*I)*\[Xi])*(-41117 + \[Nu]*(-462659 + 24575*\[Nu])))/18144 + 
        (-348401 - \[Nu]*(123302 + 62557*\[Nu]))/(9072*E^((4*I)*\[Xi]))) + 
      e^3*(E^((3*I)*\[Xi])*(-1432/27 + ((1320281 - 11753*\[Nu])*\[Nu])/
           48384) + (E^(I*\[Xi])*(-546032 + \[Nu]*(-632117 + 14285*\[Nu])))/
         48384 + (-189328 + \[Nu]*(-519889 + 15457*\[Nu]))/
         (48384*E^(I*\[Xi])) + (-1122208 - \[Nu]*(459619 + 144821*\[Nu]))/
         (48384*E^((3*I)*\[Xi]))) + 
      e^6*((E^((6*I)*\[Xi])*(-78904240 + (55886825 - 4196279*\[Nu])*\[Nu]))/
         241920 + (-339815 - \[Nu]*(372953 + 9979*\[Nu]))/24192 + 
        (-641191 - 2*\[Nu]*(1166348 + 30157*\[Nu]))/
         (145152*E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*(-2184614 - 
           \[Nu]*(2387339 + 197119*\[Nu])))/145152 + 
        (1054745 + \[Nu]*(-226585 + 866989*\[Nu]))/(90720*E^((4*I)*\[Xi])) + 
        (-22510693 - 2*\[Nu]*(3707906 + 2900389*\[Nu]))/
         (241920*E^((6*I)*\[Xi])) + (E^((4*I)*\[Xi])*(36131317 + 
           \[Nu]*(-44839757 + 4173977*\[Nu])))/362880) + 
      e^5*((E^((5*I)*\[Xi])*(-146291832 + (98415547 - 6238443*\[Nu])*\[Nu]))/
         774144 + (-491/189 - (\[Nu]*(6330297 + 31631*\[Nu]))/387072)/
         E^(I*\[Xi]) + (E^(I*\[Xi])*(-13378576 - 
           \[Nu]*(20348767 + 694625*\[Nu])))/1161216 + 
        (492184 + \[Nu]*(-5939153 + 2830337*\[Nu]))/
         (774144*E^((3*I)*\[Xi])) + (E^((3*I)*\[Xi])*(21738664 + 
           \[Nu]*(-43930853 + 3420269*\[Nu])))/774144 + 
        (-140886712 - \[Nu]*(46898683 + 31148909*\[Nu]))/
         (2322432*E^((5*I)*\[Xi])))) + 
    SO^2*(x^2*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
        4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
          2*\[Chi]A*\[Chi]S) + e*((15*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/(8*E^(I*\[Xi])) + 
          (15*E^(I*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/8) + 
        e^2*((3*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/2 + (21*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/(8*E^((2*I)*\[Xi])) + 
          (23*E^((2*I)*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/8) + 
        e^3*((155*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(64*E^(I*\[Xi])) + 
          (131*E^(I*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/64 + (237*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/
           (64*E^((3*I)*\[Xi])) + (277*E^((3*I)*\[Xi])*(\[Kappa]S - 
             2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/64) + 
        e^4*((15*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/8 + (5*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/(2*E^((2*I)*\[Xi])) + 
          2*E^((2*I)*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
            4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
              2*\[Chi]A*\[Chi]S)) + (167*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/
           (32*E^((4*I)*\[Xi])) + (205*E^((4*I)*\[Xi])*(\[Kappa]S - 
             2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/32) + 
        e^5*((4691*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(1536*E^(I*\[Xi])) + 
          (1313*E^(I*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/512 + (2421*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/
           (1024*E^((3*I)*\[Xi])) + (1605*E^((3*I)*\[Xi])*
            (\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
             \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/1024 + 
          (7513*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(1024*E^((5*I)*\[Xi])) + 
          (28763*E^((5*I)*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/3072) + 
        e^6*((35*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/16 + (1217*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/
           (384*E^((2*I)*\[Xi])) + (1019*E^((2*I)*\[Xi])*(\[Kappa]S - 
             2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/384 + 
          (341*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/(192*E^((4*I)*\[Xi])) + 
          (71*E^((4*I)*\[Xi])*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
             4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*
                \[Chi]S)))/192 + (1317*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/
           (128*E^((6*I)*\[Xi])) + (1735*E^((6*I)*\[Xi])*(\[Kappa]S - 
             2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
             \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/128)) + 
      x^3*\[Epsilon]^6*((9*\[Delta]*\[Kappa]A*(4 - 7*\[Nu]) - 
          3*\[Kappa]S*(-12 + \[Nu]*(45 + 68*\[Nu])) + 
          (8 + (73 - 408*\[Nu])*\[Nu])*\[Chi]A^2 + 4*\[Delta]*(4 - 91*\[Nu])*
           \[Chi]A*\[Chi]S + (8 + 7*\[Nu]*(-67 + 24*\[Nu]))*\[Chi]S^2)/63 + 
        e*((E^(I*\[Xi])*(7539*\[Kappa]S + \[Delta]*\[Kappa]A*(7539 - 4365*
                \[Nu]) - 3*\[Kappa]S*\[Nu]*(6481 + 2506*\[Nu]) - 
             \[Nu]*(45637 + 15036*\[Nu])*\[Chi]A^2 + 2*\[Delta]*
              (12355 - 14179*\[Nu])*\[Chi]A*\[Chi]S + 
             \[Nu]*(-32141 + 8400*\[Nu])*\[Chi]S^2 + 
             12355*(\[Chi]A^2 + \[Chi]S^2)))/1008 + 
          (8139*\[Kappa]S + \[Delta]*\[Kappa]A*(8139 - 5013*\[Nu]) - 
            3*\[Kappa]S*\[Nu]*(7097 + 2090*\[Nu]) - 
            \[Nu]*(62573 + 12540*\[Nu])*\[Chi]A^2 + 2*\[Delta]*
             (16763 - 15611*\[Nu])*\[Chi]A*\[Chi]S + 
            \[Nu]*(-35701 + 8400*\[Nu])*\[Chi]S^2 + 
            16763*(\[Chi]A^2 + \[Chi]S^2))/(1008*E^(I*\[Xi]))) + 
        e^2*((792*\[Kappa]S + \[Delta]*\[Kappa]A*(792 - 483*\[Nu]) - 
            3*\[Kappa]S*\[Nu]*(689 + 222*\[Nu]) - 4*\[Nu]*(547 + 333*\[Nu])*
             \[Chi]A^2 + 4*\[Delta]*(326 - 665*\[Nu])*\[Chi]A*\[Chi]S + 
            280*\[Nu]*(-11 + 3*\[Nu])*\[Chi]S^2 + 
            652*(\[Chi]A^2 + \[Chi]S^2))/126 + 
          (E^((2*I)*\[Xi])*(8865*\[Kappa]S - 3*\[Kappa]S*\[Nu]*
              (7575 + 6518*\[Nu]) - \[Nu]*(63611 + 39108*\[Nu])*\[Chi]A^2 + 
             \[Nu]*(-46519 + 11592*\[Nu])*\[Chi]S^2 + 5*\[Delta]*
              (\[Kappa]A*(1773 - 999*\[Nu]) + 2*(3733 - 3547*\[Nu])*\[Chi]A*
                \[Chi]S) + 18665*(\[Chi]A^2 + \[Chi]S^2)))/1008 + 
          (12675*\[Kappa]S + \[Delta]*\[Kappa]A*(12675 - 7929*\[Nu]) - 
            3*\[Kappa]S*\[Nu]*(11093 + 3938*\[Nu]) - 11*\[Nu]*
             (8947 + 2148*\[Nu])*\[Chi]A^2 + 2*\[Delta]*(27515 - 19941*\[Nu])*
             \[Chi]A*\[Chi]S + \[Nu]*(-51525 + 9016*\[Nu])*\[Chi]S^2 + 
            27515*(\[Chi]A^2 + \[Chi]S^2))/(1008*E^((2*I)*\[Xi]))) + 
        e^3*((E^((3*I)*\[Xi])*(21237*\[Kappa]S + 3*\[Delta]*\[Kappa]A*
              (7079 - 7703*\[Nu]) - 3*\[Kappa]S*\[Nu]*(21861 + 113458*
                \[Nu]) - \[Nu]*(294607 + 680748*\[Nu])*\[Chi]A^2 + 
             2*\[Delta]*(123829 - 90729*\[Nu])*\[Chi]A*\[Chi]S + 
             \[Nu]*(-382167 + 94304*\[Nu])*\[Chi]S^2 + 
             123829*(\[Chi]A^2 + \[Chi]S^2)))/8064 + 
          (E^(I*\[Xi])*(200451*\[Kappa]S + 3*\[Delta]*\[Kappa]A*
              (66817 - 35401*\[Nu]) - 3*\[Kappa]S*\[Nu]*(169035 + 2494*
                \[Nu]) - 3*\[Nu]*(331131 + 4988*\[Nu])*\[Chi]A^2 + 
             2*\[Delta]*(248163 - 300775*\[Nu])*\[Chi]A*\[Chi]S + 
             \[Nu]*(-600809 + 142240*\[Nu])*\[Chi]S^2 + 
             248163*(\[Chi]A^2 + \[Chi]S^2)))/8064 + 
          (155181*\[Kappa]S + 3*\[Delta]*\[Kappa]A*(51727 - 33935*\[Nu]) - 
            3*\[Kappa]S*\[Nu]*(137389 + 53106*\[Nu]) - 
            3*\[Nu]*(386485 + 106212*\[Nu])*\[Chi]A^2 + 2*\[Delta]*
             (337293 - 194905*\[Nu])*\[Chi]A*\[Chi]S + 
            \[Nu]*(-579527 + 69664*\[Nu])*\[Chi]S^2 + 
            337293*(\[Chi]A^2 + \[Chi]S^2))/(8064*E^((3*I)*\[Xi])) + 
          (258555*\[Kappa]S + 3*\[Delta]*\[Kappa]A*(86185 - 44321*\[Nu]) - 
            3*\[Kappa]S*\[Nu]*(216691 + 5566*\[Nu]) - 
            \[Nu]*(1683713 + 33396*\[Nu])*\[Chi]A^2 + 2*\[Delta]*
             (423419 - 391767*\[Nu])*\[Chi]A*\[Chi]S + 
            \[Nu]*(-793497 + 183904*\[Nu])*\[Chi]S^2 + 
            423419*(\[Chi]A^2 + \[Chi]S^2))/(8064*E^(I*\[Xi]))) + 
        e^6*((5*(\[Kappa]S*(155 - 2*\[Nu]*(191 + 14*\[Nu])) + 
             (161 - \[Nu]*(619 + 56*\[Nu]))*\[Chi]A^2 + 
             (161 + \[Nu]*(-437 + 104*\[Nu]))*\[Chi]S^2 + 
             \[Delta]*(\[Kappa]A*(155 - 72*\[Nu]) + 2*(161 - 206*\[Nu])*
                \[Chi]A*\[Chi]S)))/24 + (E^((6*I)*\[Xi])*
            (3*\[Delta]*\[Kappa]A*(-3535889 + 1248900*\[Nu]) - 
             3*\[Kappa]S*(3535889 - 8320678*\[Nu] + 6928292*\[Nu]^2) + 
             (57851191 - 41569752*\[Nu])*\[Nu]*\[Chi]A^2 + 
             2*\[Delta]*(-11295235 + 12778737*\[Nu])*\[Chi]A*\[Chi]S + 
             (12887223 - 2656640*\[Nu])*\[Nu]*\[Chi]S^2 - 
             11295235*(\[Chi]A^2 + \[Chi]S^2)))/80640 + 
          (E^((2*I)*\[Xi])*(3*\[Kappa]S*(984163 + 64*\[Nu]*(-38218 + 
                 1073*\[Nu])) + \[Nu]*(-15547577 + 412032*\[Nu])*\[Chi]A^2 + 
             5*\[Nu]*(-1671597 + 386960*\[Nu])*\[Chi]S^2 + 
             \[Delta]*(3*\[Kappa]A*(984163 - 477626*\[Nu]) + 2*(3848489 - 
                 4255803*\[Nu])*\[Chi]A*\[Chi]S) + 3848489*(\[Chi]A^2 + 
               \[Chi]S^2)))/48384 + (3681327*\[Kappa]S + \[Delta]*\[Kappa]A*
             (3681327 - 1801938*\[Nu]) + 48*\[Kappa]S*\[Nu]*
             (-190929 + 2246*\[Nu]) + 3*\[Nu]*(-7944837 + 71872*\[Nu])*
             \[Chi]A^2 + 2*\[Delta]*(5966799 - 5166061*\[Nu])*\[Chi]A*
             \[Chi]S + \[Nu]*(-10364807 + 2235184*\[Nu])*\[Chi]S^2 + 
            5966799*(\[Chi]A^2 + \[Chi]S^2))/(48384*E^((2*I)*\[Xi])) + 
          (\[Delta]*\[Kappa]A*(5344971 - 4057356*\[Nu]) - 
            3*\[Kappa]S*(-1781657 + 4915766*\[Nu] + 2000316*\[Nu]^2) - 
            3*\[Nu]*(11166485 + 4000632*\[Nu])*\[Chi]A^2 + 
            2*\[Delta]*(10875195 - 2834393*\[Nu])*\[Chi]A*\[Chi]S - 
            \[Nu]*(15670111 + 456064*\[Nu])*\[Chi]S^2 + 
            10875195*(\[Chi]A^2 + \[Chi]S^2))/(80640*E^((6*I)*\[Xi])) + 
          (7121193*\[Kappa]S + \[Delta]*\[Kappa]A*(7121193 - 3039645*\[Nu]) - 
            3*\[Kappa]S*\[Nu]*(5760677 + 45482*\[Nu]) - 
            \[Nu]*(51701959 + 272892*\[Nu])*\[Chi]A^2 + 2*\[Delta]*
             (12674545 - 11183823*\[Nu])*\[Chi]A*\[Chi]S + 
            7*\[Nu]*(-3051981 + 665480*\[Nu])*\[Chi]S^2 + 
            12674545*(\[Chi]A^2 + \[Chi]S^2))/(120960*E^((4*I)*\[Xi])) + 
          (E^((4*I)*\[Xi])*(3*\[Kappa]S*(5106557 + \[Nu]*(-12370131 + 
                 2139514*\[Nu])) + 57*\[Nu]*(-1551505 + 225212*\[Nu])*
              \[Chi]A^2 + \[Nu]*(-38532661 + 7884296*\[Nu])*\[Chi]S^2 + 
             \[Delta]*(\[Kappa]A*(15319671 - 6471051*\[Nu]) + 2*(21176655 - 
                 21130913*\[Nu])*\[Chi]A*\[Chi]S) + 21176655*
              (\[Chi]A^2 + \[Chi]S^2)))/120960) + 
        e^5*((E^(I*\[Xi])*(8743791*\[Kappa]S + \[Delta]*\[Kappa]A*
              (8743791 - 4567869*\[Nu]) + 3*\[Kappa]S*\[Nu]*
              (-7351817 + 267526*\[Nu]) + \[Nu]*(-41041985 + 1605156*\[Nu])*
              \[Chi]A^2 + 2*\[Delta]*(10103135 - 12151431*\[Nu])*\[Chi]A*
              \[Chi]S + \[Nu]*(-23673417 + 5501104*\[Nu])*\[Chi]S^2 + 
             10103135*(\[Chi]A^2 + \[Chi]S^2)))/193536 + 
          (6379185*\[Kappa]S + \[Delta]*\[Kappa]A*(6379185 - 2994351*\[Nu]) - 
            3*\[Kappa]S*\[Nu]*(5250907 + 192750*\[Nu]) - 
            15*\[Nu]*(2986361 + 77100*\[Nu])*\[Chi]A^2 + 
            2*\[Delta]*(11249505 - 9532729*\[Nu])*\[Chi]A*\[Chi]S + 
            \[Nu]*(-19268063 + 4015984*\[Nu])*\[Chi]S^2 + 
            11249505*(\[Chi]A^2 + \[Chi]S^2))/(129024*E^((3*I)*\[Xi])) + 
          (E^((3*I)*\[Xi])*(9421017*\[Kappa]S + \[Delta]*\[Kappa]A*
              (9421017 - 4170999*\[Nu]) + 3*\[Kappa]S*\[Nu]*
              (-7671011 + 736338*\[Nu]) + \[Nu]*(-53971423 + 4418028*\[Nu])*
              \[Chi]A^2 + 2*\[Delta]*(13167529 - 13432209*\[Nu])*\[Chi]A*
              \[Chi]S + 7*\[Nu]*(-3651873 + 780560*\[Nu])*\[Chi]S^2 + 
             13167529*(\[Chi]A^2 + \[Chi]S^2)))/129024 + 
          (3*\[Kappa]S*(4275213 + \[Nu]*(-10683537 + 366886*\[Nu])) + 
            \[Nu]*(-79554041 + 2201316*\[Nu])*\[Chi]A^2 + 
            \[Nu]*(-35435137 + 8082480*\[Nu])*\[Chi]S^2 + 
            \[Delta]*(9*\[Kappa]A*(1425071 - 711037*\[Nu]) + 2*(19721591 - 
                18051407*\[Nu])*\[Chi]A*\[Chi]S) + 19721591*
             (\[Chi]A^2 + \[Chi]S^2))/(193536*E^(I*\[Xi])) + 
          (17080167*\[Kappa]S + 3*\[Delta]*\[Kappa]A*(5693389 - 
              4133499*\[Nu]) - 3*\[Kappa]S*\[Nu]*(15520277 + 6322178*\[Nu]) - 
            \[Nu]*(113866865 + 37933068*\[Nu])*\[Chi]A^2 + 
            2*\[Delta]*(35685719 - 12586095*\[Nu])*\[Chi]A*\[Chi]S + 
            \[Nu]*(-54048201 + 1002064*\[Nu])*\[Chi]S^2 + 
            35685719*(\[Chi]A^2 + \[Chi]S^2))/(387072*E^((5*I)*\[Xi])) + 
          (E^((5*I)*\[Xi])*(9*\[Delta]*\[Kappa]A*(-2421353 + 801135*\[Nu]) - 
             3*\[Kappa]S*(7264059 + \[Nu]*(-16931523 + 19240706*\[Nu])) + 
             2*\[Delta]*(-19117729 + 24061945*\[Nu])*\[Chi]A*\[Chi]S - 
             19117729*(\[Chi]A^2 + \[Chi]S^2) + \[Nu]*((111490759 - 
                 115444236*\[Nu])*\[Chi]A^2 + (13104047 - 2343600*\[Nu])*
                \[Chi]S^2)))/387072) + 
        e^4*((-3*\[Kappa]S*(-5802 + \[Nu]*(14509 + 2022*\[Nu])) - 
            2*(-8479 + \[Nu]*(31669 + 6066*\[Nu]))*\[Chi]A^2 + 
            2*(8479 + 7*\[Nu]*(-3827 + 948*\[Nu]))*\[Chi]S^2 + 
            \[Delta]*(\[Kappa]A*(17406 - 8715*\[Nu]) + 4*(8479 - 12271*\[Nu])*
               \[Chi]A*\[Chi]S))/1008 + (E^((4*I)*\[Xi])*
            (3*\[Delta]*\[Kappa]A*(-65203 + 15759*\[Nu]) - 
             3*\[Kappa]S*(65203 + \[Nu]*(-146165 + 331438*\[Nu])) + 
             (770087 - 1988628*\[Nu])*\[Nu]*\[Chi]A^2 + 62*\[Delta]*
              (-1391 + 4517*\[Nu])*\[Chi]A*\[Chi]S + \[Nu]*(-317549 + 84504*
                \[Nu])*\[Chi]S^2 - 43121*(\[Chi]A^2 + \[Chi]S^2)))/12096 + 
          (3*\[Delta]*\[Kappa]A*(40510 - 20013*\[Nu]) - 3*\[Kappa]S*
             (-40510 + \[Nu]*(101033 + 4334*\[Nu])) - 
            4*\[Nu]*(207623 + 6501*\[Nu])*\[Chi]A^2 + 44*\[Delta]*
             (9551 - 8225*\[Nu])*\[Chi]A*\[Chi]S + 56*\[Nu]*
             (-6641 + 1422*\[Nu])*\[Chi]S^2 + 210122*(\[Chi]A^2 + \[Chi]S^2))/
           (3024*E^((2*I)*\[Xi])) + (353505*\[Kappa]S + \[Delta]*\[Kappa]A*
             (353505 - 244305*\[Nu]) - 3*\[Kappa]S*\[Nu]*
             (317105 + 127418*\[Nu]) - \[Nu]*(2502287 + 764508*\[Nu])*
             \[Chi]A^2 + 2*\[Delta]*(756089 - 345539*\[Nu])*\[Chi]A*\[Chi]S + 
            \[Nu]*(-1213147 + 82824*\[Nu])*\[Chi]S^2 + 
            756089*(\[Chi]A^2 + \[Chi]S^2))/(12096*E^((4*I)*\[Xi])) + 
          (E^((2*I)*\[Xi])*(3*\[Kappa]S*(43526 + \[Nu]*(-107623 + 
                 3038*\[Nu])) + \[Delta]*(\[Kappa]A*(130578 - 61713*\[Nu]) + 
               4*(89873 - 96275*\[Nu])*\[Chi]A*\[Chi]S) + 
             179746*(\[Chi]A^2 + \[Chi]S^2) + 4*\[Nu]*((-181372 + 4557*\[Nu])*
                \[Chi]A^2 + (-94649 + 21252*\[Nu])*\[Chi]S^2)))/3024))) + 
    SO*(x^(3/2)*\[Epsilon]^3*((-4*(\[Delta]*\[Chi]A + \[Chi]S - 
           \[Nu]*\[Chi]S))/3 + e*((-14*\[Delta]*\[Chi]A + (-14 + 5*\[Nu])*
             \[Chi]S)/(6*E^(I*\[Xi])) + (E^(I*\[Xi])*(-16*\[Delta]*\[Chi]A + 
             (-16 + 13*\[Nu])*\[Chi]S))/6) + 
        e^2*((-4*(\[Delta]*\[Chi]A + \[Chi]S - \[Nu]*\[Chi]S))/3 + 
          (-17*\[Delta]*\[Chi]A + (-17 + 3*\[Nu])*\[Chi]S)/
           (6*E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*(-23*\[Delta]*\[Chi]A + 
             (-23 + 17*\[Nu])*\[Chi]S))/6) + 
        e^3*((-5*E^(I*\[Xi])*(18*\[Delta]*\[Chi]A + (18 - 17*\[Nu])*\[Chi]S))/
           48 + (-174*\[Delta]*\[Chi]A + (-174 + 7*\[Nu])*\[Chi]S)/
           (48*E^((3*I)*\[Xi])) + (-116*\[Delta]*\[Chi]A + (-116 + 69*\[Nu])*
             \[Chi]S)/(48*E^(I*\[Xi])) + (E^((3*I)*\[Xi])*
            (-268*\[Delta]*\[Chi]A + (-268 + 183*\[Nu])*\[Chi]S))/48) + 
        e^4*((-4*(\[Delta]*\[Chi]A + \[Chi]S - \[Nu]*\[Chi]S))/3 + 
          (-341*\[Delta]*\[Chi]A - (341 + 25*\[Nu])*\[Chi]S)/
           (72*E^((4*I)*\[Xi])) + (E^((2*I)*\[Xi])*(-41*\[Delta]*\[Chi]A + 
             (-41 + 53*\[Nu])*\[Chi]S))/36 + (-79*\[Delta]*\[Chi]A + 
            (-79 + 55*\[Nu])*\[Chi]S)/(36*E^((2*I)*\[Xi])) + 
          (E^((4*I)*\[Xi])*(-583*\[Delta]*\[Chi]A + (-583 + 373*\[Nu])*
              \[Chi]S))/72) + 
        e^6*((-2*(26*\[Delta]*\[Chi]A + (26 - 63*\[Nu])*\[Chi]S))/
           (45*E^((4*I)*\[Xi])) + (2*E^((4*I)*\[Xi])*(66*\[Delta]*\[Chi]A + 
             (66 - 5*\[Nu])*\[Chi]S))/45 - (4*(\[Delta]*\[Chi]A + \[Chi]S - 
             \[Nu]*\[Chi]S))/3 - (11*(363*\[Delta]*\[Chi]A + 
             (363 + 95*\[Nu])*\[Chi]S))/(480*E^((6*I)*\[Xi])) + 
          (-711*\[Delta]*\[Chi]A + (-711 + 385*\[Nu])*\[Chi]S)/
           (288*E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*(-605*\[Delta]*\[Chi]A + 
             (-605 + 531*\[Nu])*\[Chi]S))/288 + 
          (E^((6*I)*\[Xi])*(-8027*\[Delta]*\[Chi]A + (-8027 + 4641*\[Nu])*
              \[Chi]S))/480) + 
        e^5*((-7*(434*\[Delta]*\[Chi]A + (434 - 251*\[Nu])*\[Chi]S))/
           (1152*E^(I*\[Xi])) + (E^((3*I)*\[Xi])*(206*\[Delta]*\[Chi]A + 
             (206 + 687*\[Nu])*\[Chi]S))/768 + (-1416*\[Delta]*\[Chi]A + 
            (-1416 + 1495*\[Nu])*\[Chi]S)/(768*E^((3*I)*\[Xi])) + 
          (E^(I*\[Xi])*(-2576*\[Delta]*\[Chi]A + (-2576 + 2229*\[Nu])*
              \[Chi]S))/1152 + (-14414*\[Delta]*\[Chi]A - 
            (14414 + 2487*\[Nu])*\[Chi]S)/(2304*E^((5*I)*\[Xi])) + 
          (E^((5*I)*\[Xi])*(-26888*\[Delta]*\[Chi]A + (-26888 + 16289*\[Nu])*
              \[Chi]S))/2304)) + x^(5/2)*\[Epsilon]^5*
       ((-16*\[Delta]*(10 + 7*\[Nu])*\[Chi]A + 2*(5 + 4*\[Nu])*
           (-16 + 33*\[Nu])*\[Chi]S)/63 + 
        e*((-(\[Delta]*(6506 + 1925*\[Nu])*\[Chi]A) + 
            (-6506 + \[Nu]*(9077 + 4368*\[Nu]))*\[Chi]S)/(504*E^(I*\[Xi])) + 
          (E^(I*\[Xi])*(-(\[Delta]*(2950 + 3013*\[Nu])*\[Chi]A) + 
             (-2950 + \[Nu]*(6949 + 4392*\[Nu]))*\[Chi]S))/504) + 
        e^2*((-(\[Delta]*(1048 + 281*\[Nu])*\[Chi]A) + 
            (-1048 + 3*\[Nu]*(643 + 336*\[Nu]))*\[Chi]S)/126 + 
          (-(\[Delta]*(9094 + 4835*\[Nu])*\[Chi]A) + 
            (-9094 + \[Nu]*(12459 + 8966*\[Nu]))*\[Chi]S)/
           (504*E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*
            (\[Delta]*(1678 - 10025*\[Nu])*\[Chi]A + 
             (1678 + \[Nu]*(3905 + 9282*\[Nu]))*\[Chi]S))/504) + 
        e^3*((E^(I*\[Xi])*(\[Delta]*(-101566 + 18183*\[Nu])*\[Chi]A + 
             (-101566 + \[Nu]*(185153 + 30368*\[Nu]))*\[Chi]S))/4032 + 
          (11*\[Delta]*(-13534 + 773*\[Nu])*\[Chi]A + 
            (-148874 + \[Nu]*(208193 + 41960*\[Nu]))*\[Chi]S)/
           (4032*E^(I*\[Xi])) + (-(\[Delta]*(101594 + 74931*\[Nu])*\[Chi]A) + 
            (-101594 + \[Nu]*(139819 + 133232*\[Nu]))*\[Chi]S)/
           (4032*E^((3*I)*\[Xi])) + (E^((3*I)*\[Xi])*
            (\[Delta]*(126562 - 188227*\[Nu])*\[Chi]A + 
             (126562 + \[Nu]*(-76149 + 141256*\[Nu]))*\[Chi]S))/4032) + 
        e^4*((-2*\[Delta]*(2343 + 149*\[Nu])*\[Chi]A + 
            (-4686 + \[Nu]*(8420 + 2787*\[Nu]))*\[Chi]S)/252 + 
          (-(\[Delta]*(213794 + 197291*\[Nu])*\[Chi]A) + 
            (-213794 + 47*\[Nu]*(6373 + 7434*\[Nu]))*\[Chi]S)/
           (6048*E^((4*I)*\[Xi])) + (E^((2*I)*\[Xi])*
            (\[Delta]*(-150074 + 40069*\[Nu])*\[Chi]A + 
             (-150074 + \[Nu]*(238271 + 20646*\[Nu]))*\[Chi]S))/3024 + 
          (\[Delta]*(-129454 + 4079*\[Nu])*\[Chi]A + 
            (-129454 + \[Nu]*(173197 + 35778*\[Nu]))*\[Chi]S)/
           (3024*E^((2*I)*\[Xi])) + (E^((4*I)*\[Xi])*
            (\[Delta]*(545042 - 566845*\[Nu])*\[Chi]A + 
             (545042 + \[Nu]*(-476027 + 375282*\[Nu]))*\[Chi]S))/6048) + 
        e^6*((E^((4*I)*\[Xi])*(\[Delta]*(-6484648 + 2921469*\[Nu])*\[Chi]A + 
             (-6484648 + (8214683 - 616264*\[Nu])*\[Nu])*\[Chi]S))/30240 + 
          (\[Delta]*(-1638560 + 445297*\[Nu])*\[Chi]A + 
            (-1638560 + (1911479 - 318508*\[Nu])*\[Nu])*\[Chi]S)/
           (30240*E^((4*I)*\[Xi])) + (2*\[Delta]*(-7591 + 71*\[Nu])*\[Chi]A + 
            (-15182 + \[Nu]*(26804 + 7011*\[Nu]))*\[Chi]S)/504 + 
          (E^((2*I)*\[Xi])*(\[Delta]*(-1121786 + 212269*\[Nu])*\[Chi]A + 
             (-1121786 + \[Nu]*(2122667 + 241760*\[Nu]))*\[Chi]S))/24192 + 
          (7*\[Delta]*(-240202 + 18993*\[Nu])*\[Chi]A + 
            (-1681414 + \[Nu]*(2343689 + 425312*\[Nu]))*\[Chi]S)/
           (24192*E^((2*I)*\[Xi])) + (-(\[Delta]*(2843146 + 3523827*\[Nu])*
              \[Chi]A) + (-2843146 + \[Nu]*(4219211 + 6378352*\[Nu]))*
             \[Chi]S)/(40320*E^((6*I)*\[Xi])) + 
          (E^((6*I)*\[Xi])*(\[Delta]*(15560570 - 12084401*\[Nu])*\[Chi]A + 
             (15560570 + \[Nu]*(-15789807 + 6886784*\[Nu]))*\[Chi]S))/
           40320) + e^5*((E^((3*I)*\[Xi])*(\[Delta]*(-6641638 + 2484809*
                \[Nu])*\[Chi]A + (-6641638 + (9249991 - 15704*\[Nu])*\[Nu])*
              \[Chi]S))/64512 + (\[Delta]*(-6241442 + 844507*\[Nu])*\[Chi]A + 
            (-6241442 + 7*\[Nu]*(1264307 + 176976*\[Nu]))*\[Chi]S)/
           (96768*E^(I*\[Xi])) + (\[Delta]*(-3162010 + 293769*\[Nu])*
             \[Chi]A + (-3162010 + \[Nu]*(4015895 + 452912*\[Nu]))*\[Chi]S)/
           (64512*E^((3*I)*\[Xi])) + (E^(I*\[Xi])*
            (\[Delta]*(-3824158 + 907771*\[Nu])*\[Chi]A + 
             (-3824158 + \[Nu]*(7078005 + 789608*\[Nu]))*\[Chi]S))/96768 + 
          (-(\[Delta]*(9643306 + 10508705*\[Nu])*\[Chi]A) + 
            (-9643306 + \[Nu]*(13869681 + 18771944*\[Nu]))*\[Chi]S)/
           (193536*E^((5*I)*\[Xi])) + (E^((5*I)*\[Xi])*
            (\[Delta]*(38415658 - 33357473*\[Nu])*\[Chi]A + 
             (38415658 + \[Nu]*(-37214623 + 20267616*\[Nu]))*\[Chi]S))/
           193536)) + x^3*\[Epsilon]^6*
       ((-4*((I + 2*Pi)*\[Delta]*\[Chi]A + (I + 2*Pi - 2*(I + Pi)*\[Nu])*
            \[Chi]S))/3 + e*((\[Delta]*\[Chi]A*(167*I - 128*Pi - 
              (1728*I)*ArcCoth[5]) + \[Chi]S*(167*I - (10*I)*\[Nu] + 
              8*Pi*(-16 + 7*\[Nu]) + (432*I)*(-2 + \[Nu])*Log[3/2]))/
           (12*E^(I*\[Xi])) + (E^(I*\[Xi])*(4*Pi*(-37 + 22*\[Nu])*\[Chi]S + 
             \[Delta]*\[Chi]A*(43*I - 148*Pi - (312*I)*Log[2]) + 
             I*\[Chi]S*(43 + 22*\[Nu] + 24*(-13 + 6*\[Nu])*Log[2])))/12) + 
        e^2*((E^((2*I)*\[Xi])*(Pi*(-175 + 88*\[Nu])*\[Chi]S + 
             \[Delta]*\[Chi]A*(49*I - 175*Pi - (402*I)*Log[2]) + 
             I*\[Chi]S*(49 + 10*\[Nu] + 6*(-67 + 30*\[Nu])*Log[2])))/6 + 
          (-47*Pi*\[Delta]*\[Chi]A + Pi*(-47 + 26*\[Nu])*\[Chi]S + 
            I*\[Delta]*\[Chi]A*(56 + 1038*Log[2] - 837*Log[3]) - 
            I*\[Chi]S*(-56 - 1038*Log[2] + \[Nu]*(4 + 468*Log[2] - 
                378*Log[3]) + 837*Log[3]))/3 + (-131*Pi*\[Delta]*\[Chi]A + 
            Pi*(-131 + 36*\[Nu])*\[Chi]S + I*\[Delta]*\[Chi]A*
             (115 - 3274*Log[2] + 1674*Log[3]) + I*\[Chi]S*(115 + 22*\[Nu] - 
              3274*Log[2] + 1428*\[Nu]*Log[2] + 54*(31 - 14*\[Nu])*Log[3]))/
           (6*E^((2*I)*\[Xi]))) + 
        e^3*((E^((3*I)*\[Xi])*(\[Delta]*\[Chi]A*(4473*I - 16228*Pi - 
               (35752*I)*Log[2]) + \[Chi]S*((9*I)*(497 + 26*\[Nu]) + 4*Pi*
                (-4057 + 1879*\[Nu]) + (8*I)*(-4469 + 1973*\[Nu])*Log[2])))/
           288 + (E^(I*\[Xi])*(-2480*Pi*\[Delta]*\[Chi]A + 
             20*Pi*(-124 + 65*\[Nu])*\[Chi]S + I*\[Delta]*\[Chi]A*
              (3153 + 72128*Log[2] - 57456*Log[3]) - I*\[Chi]S*
              (-3153 - 72128*Log[2] + 30*\[Nu]*(17 + 1036*Log[2] - 
                 828*Log[3]) + 57456*Log[3])))/96 + 
          (-2572*Pi*\[Delta]*\[Chi]A + 4*Pi*(-643 + 303*\[Nu])*\[Chi]S + 
            I*\[Delta]*\[Chi]A*(3485 - 202952*Log[2] + 117288*Log[3]) + 
            I*\[Chi]S*(3485 - 526*\[Nu] - 202952*Log[2] + 83400*\[Nu]*Log[
                2] - 648*(-181 + 74*\[Nu])*Log[3]))/(96*E^(I*\[Xi])) + 
          (-11184*Pi*\[Delta]*\[Chi]A + 12*Pi*(-932 + 155*\[Nu])*\[Chi]S + 
            I*\[Delta]*\[Chi]A*(7197 + 513856*Log[2] - 66096*Log[3] - 
              202000*Log[5]) - I*\[Chi]S*(-7197 - 3642*\[Nu] - 
              513856*Log[2] + 206536*\[Nu]*Log[2] + 66096*Log[3] - 
              29160*\[Nu]*Log[3] + 2000*(101 - 38*\[Nu])*Log[5]))/
           (288*E^((3*I)*\[Xi]))) + 
        e^4*((E^((4*I)*\[Xi])*(-6949*Pi*\[Delta]*\[Chi]A + 
             Pi*(-6949 + 3076*\[Nu])*\[Chi]S - (4*I)*\[Delta]*\[Chi]A*
              (-482 + 3677*Log[2]) + (4*I)*\[Chi]S*(482 - 19*\[Nu] + 
               (-3677 + 1598*\[Nu])*Log[2])))/72 + 
          (E^((2*I)*\[Xi])*(-2793*Pi*\[Delta]*\[Chi]A + 
             3*Pi*(-931 + 458*\[Nu])*\[Chi]S + I*\[Delta]*\[Chi]A*
              (3314 + 96488*Log[2] - 77085*Log[3]) - I*\[Chi]S*
              (-3314 - 96488*Log[2] + \[Nu]*(580 + 40448*Log[2] - 
                 32346*Log[3]) + 77085*Log[3])))/72 + 
          (-205*Pi*\[Delta]*\[Chi]A + 5*Pi*(-41 + 20*\[Nu])*\[Chi]S + 
            I*\[Delta]*\[Chi]A*(259 - 26342*Log[2] + 15579*Log[3]) + 
            I*\[Chi]S*(259 - 32*\[Nu] - 26342*Log[2] + 10548*\[Nu]*Log[2] + 
              27*(577 - 230*\[Nu])*Log[3]))/6 + (-2716*Pi*\[Delta]*\[Chi]A + 
            4*Pi*(-679 + 313*\[Nu])*\[Chi]S + (2*I)*\[Delta]*\[Chi]A*
             (1607 + 248620*Log[2] - 54756*Log[3]) - (2*I)*\[Chi]S*
             (-1607 + 190*\[Nu] - 248620*Log[2] + 95200*\[Nu]*Log[2] - 
              4212*(-13 + 5*\[Nu])*Log[3]) - (125*I)*(1183*\[Delta]*\[Chi]A + 
              (1183 - 454*\[Nu])*\[Chi]S)*Log[5])/(72*E^((2*I)*\[Xi])) + 
          (-4654*Pi*\[Delta]*\[Chi]A + Pi*(-4654 + 376*\[Nu])*\[Chi]S + 
            I*\[Chi]S*(2260 - 191876*Log[2] - 112779*Log[3] + 
              \[Nu]*(2068 + 75416*Log[2] + 38772*Log[3] - 56750*Log[5]) + 
              147875*Log[5]) - I*\[Delta]*\[Chi]A*(191876*Log[2] + 
              112779*Log[3] - 5*(452 + 29575*Log[5])))/
           (72*E^((4*I)*\[Xi]))) + 
        e^5*((E^(I*\[Xi])*(-100420*Pi*\[Delta]*\[Chi]A + 
             4*Pi*(-25105 + 12306*\[Nu])*\[Chi]S + I*\[Delta]*\[Chi]A*
              (140555 - 17877912*Log[2] + 10683360*Log[3]) + 
             I*\[Chi]S*(140555 - 23674*\[Nu] - 17877912*Log[2] + 7030272*
                \[Nu]*Log[2] + 216*(49460 - 19363*\[Nu])*Log[3])))/2304 + 
          (E^((5*I)*\[Xi])*(-3538812*Pi*\[Delta]*\[Chi]A + 
             12*Pi*(-294901 + 127124*\[Nu])*\[Chi]S - I*\[Delta]*\[Chi]A*
              (-1006495 + 7201112*Log[2] + 130248*Log[3]) + 
             I*\[Chi]S*(1006495 - 105170*\[Nu] - 7201112*Log[2] + 3124928*
                \[Nu]*Log[2] + 1944*(-67 + 8*\[Nu])*Log[3])))/23040 + 
          (E^((3*I)*\[Xi])*(-263968*Pi*\[Delta]*\[Chi]A + 
             16*Pi*(-16498 + 7609*\[Nu])*\[Chi]S + I*\[Delta]*\[Chi]A*
              (288099 + 10173920*Log[2] - 8110368*Log[3]) - 
             I*\[Chi]S*(-288099 + 53154*\[Nu] - 10173920*Log[2] + 4170080*
                \[Nu]*Log[2] - 3888*(-2086 + 855*\[Nu])*Log[3])))/4608 + 
          (-110480*Pi*\[Delta]*\[Chi]A + 16*Pi*(-6905 + 3107*\[Nu])*\[Chi]S + 
            I*\[Delta]*\[Chi]A*(144839 + 33053696*Log[2] - 7799328*Log[3] - 
              9302000*Log[5]) - I*\[Chi]S*(-144839 + 21946*\[Nu] - 
              33053696*Log[2] + 12396704*\[Nu]*Log[2] + 7799328*Log[3] - 
              2857680*\[Nu]*Log[3] + 1000*(9302 - 3551*\[Nu])*Log[5]))/
           (2304*E^(I*\[Xi])) + (-220164*Pi*\[Delta]*\[Chi]A + 
            12*Pi*(-18347 + 9440*\[Nu])*\[Chi]S + I*\[Chi]S*
             (255255 - 39522*\[Nu] - 55608152*Log[2] + 19956800*\[Nu]*Log[
                2] - 17765568*Log[3] + 6367248*\[Nu]*Log[3] - 
              5000*(-7063 + 2524*\[Nu])*Log[5]) - I*\[Delta]*\[Chi]A*
             (55608152*Log[2] + 17765568*Log[3] - 35*(7293 + 1009000*
                 Log[5])))/(4608*E^((3*I)*\[Xi])) + 
          (\[Delta]*\[Chi]A*(-2383520*Pi + I*(876235 + 103795872*Log[2] + 
                99839520*Log[3] - 31620000*Log[5] - 72518432*Log[7])) + 
            \[Chi]S*(80*Pi*(-29794 + 249*\[Nu]) + I*(876235 + 1289710*\[Nu] + 
                103795872*Log[2] + 99839520*Log[3] - 31620000*Log[5] - 
                72518432*Log[7]) - (48*I)*\[Nu]*(746426*Log[2] + 
                748935*Log[3] - 266875*Log[5] - 480886*Log[7])))/
           (23040*E^((5*I)*\[Xi]))) + 
        e^6*((E^((2*I)*\[Xi])*(-58835*Pi*\[Delta]*\[Chi]A + 
             Pi*(-58835 + 28902*\[Nu])*\[Chi]S + (7*I)*\[Delta]*\[Chi]A*
              (12292 - 2103704*Log[2] + 1265949*Log[3]) + 
             I*\[Chi]S*(86044 - 14725928*Log[2] + 2*\[Nu]*(-7348 + 
                 2851080*Log[2] - 1708695*Log[3]) + 8861643*Log[3])))/1152 + 
          (E^((6*I)*\[Xi])*(-4024913*Pi*\[Delta]*\[Chi]A + 
             Pi*(-4024913 + 1708922*\[Nu])*\[Chi]S - (9*I)*\[Delta]*\[Chi]A*
              (-132172 + 911784*Log[2] + 28593*Log[3]) + (3*I)*\[Chi]S*
              (396516 - 2735352*Log[2] + 8*\[Nu]*(-7707 + 145046*Log[2]) + 
               243*(-353 + 66*\[Nu])*Log[3])))/17280 + 
          (E^((4*I)*\[Xi])*(\[Delta]*\[Chi]A*(-231080*Pi + (3*I)*
                (78576 + 3331000*Log[2] - 2634615*Log[3])) + 
             \[Chi]S*(4*Pi*(-57770 + 24989*\[Nu]) - (8*I)*(-29466 + 
                 5490*\[Nu] - 1249125*Log[2] + 503854*\[Nu]*Log[2]) + (27*I)*
                (-292735 + 117986*\[Nu])*Log[3])))/2880 + 
          (20*Pi*(-3539 + 1510*\[Nu])*\[Chi]S + (4*I)*\[Chi]S*
             (20697 - 7531874*Log[2] + \[Nu]*(-2262 + 2623268*Log[2])) + 
            (54*I)*(-149683 + 55874*\[Nu])*\[Chi]S*Log[3] - 
            2*\[Delta]*\[Chi]A*(35390*Pi + I*(-41394 + 15063748*Log[2] + 
                4041441*Log[3])) + (625*I)*(29171*\[Delta]*\[Chi]A + 
              (29171 - 10358*\[Nu])*\[Chi]S)*Log[5])/(1152*E^((2*I)*\[Xi])) + 
          (-98177*Pi*\[Delta]*\[Chi]A + Pi*(-98177 + 45410*\[Nu])*\[Chi]S + 
            I*\[Delta]*\[Chi]A*(123552 + 44079056*Log[2] - 10906569*Log[3] - 
              11930875*Log[5]) - I*\[Chi]S*(-123552 - 44079056*Log[2] + 
              10906569*Log[3] + 2*\[Nu]*(8208 + 8164336*Log[2] - 
                1966113*Log[3] - 2255375*Log[5]) + 11930875*Log[5]))/1728 + 
          (-310184*Pi*\[Delta]*\[Chi]A + 1528*Pi*(-203 + 138*\[Nu])*\[Chi]S + 
            I*\[Delta]*\[Chi]A*(401824 + 106489936*Log[2] + 91977201*Log[3] - 
              51084375*Log[5] - 48612361*Log[7]) - I*\[Chi]S*
             (-401824 - 106489936*Log[2] - 91977201*Log[3] + 51084375*Log[
                5] + 2*\[Nu]*(57232 + 18072432*Log[2] + 15415677*Log[3] - 
                8634375*Log[5] - 8241261*Log[7]) + 48612361*Log[7]))/
           (5760*E^((4*I)*\[Xi])) + (-2787228*Pi*\[Delta]*\[Chi]A - 
            36*Pi*(77423 + 4294*\[Nu])*\[Chi]S - I*\[Delta]*\[Chi]A*
             (-763092 + 350349512*Log[2] + 47601513*Log[3] + 916250*Log[5] - 
              145837083*Log[7]) + I*\[Chi]S*(763092 - 350349512*Log[2] - 
              47601513*Log[3] - 916250*Log[5] + 2*\[Nu]*(869220 + 
                54990664*Log[2] + 9267777*Log[3] + 96250*Log[5] - 
                24723783*Log[7]) + 145837083*Log[7]))/
           (17280*E^((6*I)*\[Xi]))))) + x^(3/2)*\[Epsilon]^3*
     (2*Pi + e*((11*Pi + (54*I)*Log[3/2])/(4*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*(13*Pi + (6*I)*Log[2]))/4) + 
      e^3*((E^((3*I)*\[Xi])*(703*Pi + (454*I)*Log[2]))/96 + 
        (E^(I*\[Xi])*(59*Pi - (2*I)*(811*Log[2] - 702*Log[3])))/32 + 
        (61*Pi + (4838*I)*Log[2] - (2754*I)*Log[3])/(32*E^(I*\[Xi])) + 
        (555*Pi - (2*I)*(6683*Log[2] - 486*Log[3] - 3125*Log[5]))/
         (96*E^((3*I)*\[Xi]))) + 
      e^4*(2*Pi + (5*E^((4*I)*\[Xi])*(25*Pi + (17*I)*Log[2]))/12 + 
        (I/4)*(1054*Log[2] - 621*Log[3]) + E^((2*I)*\[Xi])*
         ((11*Pi)/8 - ((467*I)/6)*Log[2] + ((531*I)/8)*Log[3]) + 
        (199*Pi + I*(3674*Log[2] + 2943*Log[3] - 3125*Log[5]))/
         (24*E^((4*I)*\[Xi])) + (34*Pi - (10220*I)*Log[2] + (2106*I)*Log[3] + 
          (3125*I)*Log[5])/(24*E^((2*I)*\[Xi]))) + 
      e^6*((E^((4*I)*\[Xi])*(-1708*Pi - (3*I)*(53344*Log[2] - 44289*Log[3])))/
         960 + (E^((2*I)*\[Xi])*(689*Pi + (241592*I)*Log[2] - 
           (147501*I)*Log[3]))/384 + (E^((6*I)*\[Xi])*(115751*Pi + 
           (243*I)*(344*Log[2] + 17*Log[3])))/5760 + 
        (796*Pi + I*(533240*Log[2] + 112482*Log[3] - 303125*Log[5]))/
         (384*E^((2*I)*\[Xi])) + (1145*Pi - (734864*I)*Log[2] + 
          (200961*I)*Log[3] + (184375*I)*Log[5])/576 + 
        (96516*Pi + I*(6639176*Log[2] + 527067*Log[3] + 31250*Log[5] - 
            2470629*Log[7]))/(5760*E^((6*I)*\[Xi])) + 
        ((-32*Pi)/15 - (I/1920)*(1817152*Log[2] + 1614087*Log[3] - 
            890625*Log[5] - 823543*Log[7]))/E^((4*I)*\[Xi])) + 
      e^5*((E^((3*I)*\[Xi])*(469*Pi - (2*I)*(88579*Log[2] - 74520*Log[3])))/
         1536 + (E^(I*\[Xi])*(1465*Pi + (6*I)*(53185*Log[2] - 31986*Log[3])))/
         768 + (E^((5*I)*\[Xi])*(111771*Pi + (76406*I)*Log[2] + 
           (2754*I)*Log[3]))/7680 + 
        (357*Pi + (2*I)*(517483*Log[2] + 169452*Log[3] - 334375*Log[5]))/
         (1536*E^((3*I)*\[Xi])) + (1547*Pi - (2*I)*(301531*Log[2] - 
            75735*Log[3] - 81250*Log[5]))/(768*E^(I*\[Xi])) + 
        (90755*Pi - (2*I)*(1002963*Log[2] + 930960*Log[3] - 187500*Log[5] - 
            823543*Log[7]))/(7680*E^((5*I)*\[Xi]))) + 
      e^2*(2*Pi + E^((2*I)*\[Xi])*(5*Pi + (3*I)*Log[2]) + 
        (4*Pi + (59*I)*Log[2] - (27*I)*Log[3])/E^((2*I)*\[Xi]) + 
        (3*I)*Log[19683/1024]))
 
H\[Psi]InstTail[3, 0] = 
   x*\[Epsilon]^2*(e*((1 - 3*\[Nu])/(4*Sqrt[42]*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*(-1 + 3*\[Nu]))/(4*Sqrt[42])) + 
      e^2*((1 - 3*\[Nu])/(2*Sqrt[42]*E^((2*I)*\[Xi])) + 
        (E^((2*I)*\[Xi])*(-1 + 3*\[Nu]))/(2*Sqrt[42])) + 
      e^3*((5*(-1 + 3*\[Nu]))/(32*Sqrt[42]*E^(I*\[Xi])) - 
        (5*E^(I*\[Xi])*(-1 + 3*\[Nu]))/(32*Sqrt[42]) - 
        (9*Sqrt[3/14]*(-1 + 3*\[Nu]))/(32*E^((3*I)*\[Xi])) + 
        (9*Sqrt[3/14]*E^((3*I)*\[Xi])*(-1 + 3*\[Nu]))/32) + 
      e^4*((5*E^((2*I)*\[Xi])*(1 - 3*\[Nu]))/(12*Sqrt[42]) + 
        (2*Sqrt[2/21]*(1 - 3*\[Nu]))/(3*E^((4*I)*\[Xi])) + 
        (5*(-1 + 3*\[Nu]))/(12*Sqrt[42]*E^((2*I)*\[Xi])) + 
        (2*Sqrt[2/21]*E^((4*I)*\[Xi])*(-1 + 3*\[Nu]))/3) + 
      e^5*((3125*(1 - 3*\[Nu]))/(1536*Sqrt[42]*E^((5*I)*\[Xi])) + 
        (11*(-1 + 3*\[Nu]))/(768*Sqrt[42]*E^(I*\[Xi])) - 
        (11*E^(I*\[Xi])*(-1 + 3*\[Nu]))/(768*Sqrt[42]) + 
        (153*Sqrt[3/14]*(-1 + 3*\[Nu]))/(512*E^((3*I)*\[Xi])) - 
        (153*Sqrt[3/14]*E^((3*I)*\[Xi])*(-1 + 3*\[Nu]))/512 + 
        (3125*E^((5*I)*\[Xi])*(-1 + 3*\[Nu]))/(1536*Sqrt[42])) + 
      e^6*((1 - 3*\[Nu])/(24*Sqrt[42]*E^((2*I)*\[Xi])) + 
        (E^((2*I)*\[Xi])*(-1 + 3*\[Nu]))/(24*Sqrt[42]) + 
        (13*Sqrt[2/21]*(-1 + 3*\[Nu]))/(15*E^((4*I)*\[Xi])) - 
        (13*Sqrt[2/21]*E^((4*I)*\[Xi])*(-1 + 3*\[Nu]))/15 - 
        (81*Sqrt[3/14]*(-1 + 3*\[Nu]))/(80*E^((6*I)*\[Xi])) + 
        (81*Sqrt[3/14]*E^((6*I)*\[Xi])*(-1 + 3*\[Nu]))/80)) + 
    x^2*\[Epsilon]^4*(e*((-112 + (341 - 23*\[Nu])*\[Nu])/
         (72*Sqrt[42]*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*(112 + \[Nu]*(-341 + 23*\[Nu])))/(72*Sqrt[42])) + 
      e^2*((-191 + (583 - 25*\[Nu])*\[Nu])/(36*Sqrt[42]*E^((2*I)*\[Xi])) + 
        (E^((2*I)*\[Xi])*(191 + \[Nu]*(-583 + 25*\[Nu])))/(36*Sqrt[42])) + 
      e^6*((-845 + (2611 - 145*\[Nu])*\[Nu])/(432*Sqrt[42]*E^((2*I)*\[Xi])) - 
        (27*Sqrt[3/14]*(169 + 11*(-47 + \[Nu])*\[Nu]))/
         (160*E^((6*I)*\[Xi])) + (27*Sqrt[3/14]*E^((6*I)*\[Xi])*
          (169 + 11*(-47 + \[Nu])*\[Nu]))/160 + 
        (4901 - \[Nu]*(14863 + 131*\[Nu]))/(135*Sqrt[42]*E^((4*I)*\[Xi])) + 
        (E^((4*I)*\[Xi])*(-4901 + \[Nu]*(14863 + 131*\[Nu])))/
         (135*Sqrt[42]) + (E^((2*I)*\[Xi])*(845 + \[Nu]*(-2611 + 145*\[Nu])))/
         (432*Sqrt[42])) + 
      e^4*(-(Sqrt[2/21]*(349 + \[Nu]*(-1067 + 29*\[Nu])))/
         (27*E^((4*I)*\[Xi])) + (Sqrt[2/21]*E^((4*I)*\[Xi])*
          (349 + \[Nu]*(-1067 + 29*\[Nu])))/27 + 
        (731 - \[Nu]*(2155 + 191*\[Nu]))/(216*Sqrt[42]*E^((2*I)*\[Xi])) + 
        (E^((2*I)*\[Xi])*(-731 + \[Nu]*(2155 + 191*\[Nu])))/(216*Sqrt[42])) + 
      e^3*((-2 + (115 - 313*\[Nu])*\[Nu])/(576*Sqrt[42]*E^(I*\[Xi])) - 
        (3*Sqrt[3/14]*(90 + \[Nu]*(-275 + 9*\[Nu])))/(64*E^((3*I)*\[Xi])) + 
        (3*Sqrt[3/14]*E^((3*I)*\[Xi])*(90 + \[Nu]*(-275 + 9*\[Nu])))/64 + 
        (E^(I*\[Xi])*(2 + \[Nu]*(-115 + 313*\[Nu])))/(576*Sqrt[42])) + 
      e^5*((-7460 + (23761 - 3691*\[Nu])*\[Nu])/(13824*Sqrt[42]*
          E^(I*\[Xi])) - (3125*(428 + \[Nu]*(-1309 + 31*\[Nu])))/
         (27648*Sqrt[42]*E^((5*I)*\[Xi])) + (3125*E^((5*I)*\[Xi])*
          (428 + \[Nu]*(-1309 + 31*\[Nu])))/(27648*Sqrt[42]) - 
        (3*Sqrt[3/14]*(-1524 + \[Nu]*(4591 + 123*\[Nu])))/
         (1024*E^((3*I)*\[Xi])) + (3*Sqrt[3/14]*E^((3*I)*\[Xi])*
          (-1524 + \[Nu]*(4591 + 123*\[Nu])))/1024 + 
        (E^(I*\[Xi])*(7460 + \[Nu]*(-23761 + 3691*\[Nu])))/
         (13824*Sqrt[42]))) + SO^2*x^3*\[Epsilon]^6*
     (e*((3*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
            2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
          2*\[Delta]*(-3 + 35*\[Nu])*\[Chi]A*\[Chi]S + 
          (-3 + (61 - 32*\[Nu])*\[Nu])*\[Chi]S^2)/(6*Sqrt[42]*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*(-3*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(3 - 35*\[Nu])*\[Chi]A*\[Chi]S + 
           (3 + \[Nu]*(-61 + 32*\[Nu]))*\[Chi]S^2))/(6*Sqrt[42])) + 
      e^5*((381*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
            2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
          2*\[Delta]*(-381 + 4516*\[Nu])*\[Chi]A*\[Chi]S + 
          (-381 + 23*(343 - 152*\[Nu])*\[Nu])*\[Chi]S^2)/
         (576*Sqrt[42]*E^(I*\[Xi])) + 
        (3125*(3*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(-3 + 34*\[Nu])*\[Chi]A*\[Chi]S + 
           (-3 + (59 - 40*\[Nu])*\[Nu])*\[Chi]S^2))/(1152*Sqrt[42]*
          E^((5*I)*\[Xi])) + (9*Sqrt[3/14]*(3*(-1 + 3*\[Nu])*
            (\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             (1 - 4*\[Nu])*\[Chi]A^2) + 2*\[Delta]*(-3 + 38*\[Nu])*\[Chi]A*
            \[Chi]S + (-3 + (67 - 8*\[Nu])*\[Nu])*\[Chi]S^2))/
         (128*E^((3*I)*\[Xi])) + (9*Sqrt[3/14]*E^((3*I)*\[Xi])*
          (-3*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(3 - 38*\[Nu])*\[Chi]A*\[Chi]S + 
           (3 + \[Nu]*(-67 + 8*\[Nu]))*\[Chi]S^2))/128 + 
        (3125*E^((5*I)*\[Xi])*(-3*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + 
             \[Kappa]S - 2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(3 - 34*\[Nu])*\[Chi]A*\[Chi]S + 
           (3 + \[Nu]*(-59 + 40*\[Nu]))*\[Chi]S^2))/(1152*Sqrt[42]) + 
        (E^(I*\[Xi])*(-381*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(381 - 4516*\[Nu])*\[Chi]A*\[Chi]S + 
           (381 + 23*\[Nu]*(-343 + 152*\[Nu]))*\[Chi]S^2))/(576*Sqrt[42])) + 
      e^2*((15*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
            2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
          2*\[Delta]*(-15 + 173*\[Nu])*\[Chi]A*\[Chi]S + 
          (-15 + (301 - 176*\[Nu])*\[Nu])*\[Chi]S^2)/(12*Sqrt[42]*
          E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*(-15*(-1 + 3*\[Nu])*
            (\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             (1 - 4*\[Nu])*\[Chi]A^2) + 2*\[Delta]*(15 - 173*\[Nu])*\[Chi]A*
            \[Chi]S + (15 + \[Nu]*(-301 + 176*\[Nu]))*\[Chi]S^2))/
         (12*Sqrt[42])) + 
      e^4*((57*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
            2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
          2*\[Delta]*(-57 + 683*\[Nu])*\[Chi]A*\[Chi]S + 
          (-57 + (1195 - 464*\[Nu])*\[Nu])*\[Chi]S^2)/
         (72*Sqrt[42]*E^((2*I)*\[Xi])) + 
        (Sqrt[2/21]*(21*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(-21 + 239*\[Nu])*\[Chi]A*\[Chi]S + 
           (-21 + (415 - 272*\[Nu])*\[Nu])*\[Chi]S^2))/(9*E^((4*I)*\[Xi])) + 
        (Sqrt[2/21]*E^((4*I)*\[Xi])*(-21*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + 
             \[Kappa]S - 2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(21 - 239*\[Nu])*\[Chi]A*\[Chi]S + 
           (21 + \[Nu]*(-415 + 272*\[Nu]))*\[Chi]S^2))/9 + 
        (E^((2*I)*\[Xi])*(-57*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + 
             \[Kappa]S - 2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(57 - 683*\[Nu])*\[Chi]A*\[Chi]S + 
           (57 + \[Nu]*(-1195 + 464*\[Nu]))*\[Chi]S^2))/(72*Sqrt[42])) + 
      e^3*((51*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
            2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
          2*\[Delta]*(-51 + 605*\[Nu])*\[Chi]A*\[Chi]S + 
          (-51 + (1057 - 464*\[Nu])*\[Nu])*\[Chi]S^2)/
         (96*Sqrt[42]*E^(I*\[Xi])) + 
        (3*Sqrt[3/14]*(9*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(-9 + 103*\[Nu])*\[Chi]A*\[Chi]S + 
           (-9 + (179 - 112*\[Nu])*\[Nu])*\[Chi]S^2))/(32*E^((3*I)*\[Xi])) + 
        (3*Sqrt[3/14]*E^((3*I)*\[Xi])*(-9*(-1 + 3*\[Nu])*
            (\[Delta]*\[Kappa]A + \[Kappa]S - 2*\[Kappa]S*\[Nu] + 
             (1 - 4*\[Nu])*\[Chi]A^2) + 2*\[Delta]*(9 - 103*\[Nu])*\[Chi]A*
            \[Chi]S + (9 + \[Nu]*(-179 + 112*\[Nu]))*\[Chi]S^2))/32 + 
        (E^(I*\[Xi])*(-51*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(51 - 605*\[Nu])*\[Chi]A*\[Chi]S + 
           (51 + \[Nu]*(-1057 + 464*\[Nu]))*\[Chi]S^2))/(96*Sqrt[42])) + 
      e^6*((165*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
            2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
          2*\[Delta]*(-165 + 1951*\[Nu])*\[Chi]A*\[Chi]S + 
          (-165 + (3407 - 1552*\[Nu])*\[Nu])*\[Chi]S^2)/
         (144*Sqrt[42]*E^((2*I)*\[Xi])) + (Sqrt[7/6]*E^((4*I)*\[Xi])*
          (3*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(-3 + 25*\[Nu])*\[Chi]A*\[Chi]S + 
           (-3 + (41 - 112*\[Nu])*\[Nu])*\[Chi]S^2))/45 + 
        (27*Sqrt[3/14]*(27*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(-27 + 305*\[Nu])*\[Chi]A*\[Chi]S + 
           (-27 + 23*(23 - 16*\[Nu])*\[Nu])*\[Chi]S^2))/
         (160*E^((6*I)*\[Xi])) + (27*Sqrt[3/14]*E^((6*I)*\[Xi])*
          (-27*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(27 - 305*\[Nu])*\[Chi]A*\[Chi]S + 
           (27 + 23*\[Nu]*(-23 + 16*\[Nu]))*\[Chi]S^2))/160 + 
        (Sqrt[7/6]*(-3*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + \[Kappa]S - 
             2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(3 - 25*\[Nu])*\[Chi]A*\[Chi]S + 
           (3 + \[Nu]*(-41 + 112*\[Nu]))*\[Chi]S^2))/(45*E^((4*I)*\[Xi])) + 
        (E^((2*I)*\[Xi])*(-165*(-1 + 3*\[Nu])*(\[Delta]*\[Kappa]A + 
             \[Kappa]S - 2*\[Kappa]S*\[Nu] + (1 - 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(165 - 1951*\[Nu])*\[Chi]A*\[Chi]S + 
           (165 + \[Nu]*(-3407 + 1552*\[Nu]))*\[Chi]S^2))/(144*Sqrt[42]))) + 
    SO*(x^(3/2)*\[Epsilon]^3*(e*((\[Nu]*\[Chi]S)/(Sqrt[42]*E^(I*\[Xi])) - 
          (E^(I*\[Xi])*\[Nu]*\[Chi]S)/Sqrt[42]) + 
        e^2*((Sqrt[2/21]*\[Nu]*\[Chi]S)/E^((2*I)*\[Xi]) - 
          Sqrt[2/21]*E^((2*I)*\[Xi])*\[Nu]*\[Chi]S) + 
        e^3*(-(\[Nu]*\[Chi]S)/(8*Sqrt[42]*E^(I*\[Xi])) + 
          (E^(I*\[Xi])*\[Nu]*\[Chi]S)/(8*Sqrt[42]) + 
          (9*Sqrt[3/14]*\[Nu]*\[Chi]S)/(8*E^((3*I)*\[Xi])) - 
          (9*Sqrt[3/14]*E^((3*I)*\[Xi])*\[Nu]*\[Chi]S)/8) + 
        e^4*(-(Sqrt[2/21]*\[Nu]*\[Chi]S)/(3*E^((2*I)*\[Xi])) + 
          (Sqrt[2/21]*E^((2*I)*\[Xi])*\[Nu]*\[Chi]S)/3 + 
          (8*Sqrt[2/21]*\[Nu]*\[Chi]S)/(3*E^((4*I)*\[Xi])) - 
          (8*Sqrt[2/21]*E^((4*I)*\[Xi])*\[Nu]*\[Chi]S)/3) + 
        e^5*((\[Nu]*\[Chi]S)/(192*Sqrt[42]*E^(I*\[Xi])) - 
          (E^(I*\[Xi])*\[Nu]*\[Chi]S)/(192*Sqrt[42]) - 
          (81*Sqrt[3/14]*\[Nu]*\[Chi]S)/(128*E^((3*I)*\[Xi])) + 
          (81*Sqrt[3/14]*E^((3*I)*\[Xi])*\[Nu]*\[Chi]S)/128 + 
          (3125*\[Nu]*\[Chi]S)/(384*Sqrt[42]*E^((5*I)*\[Xi])) - 
          (3125*E^((5*I)*\[Xi])*\[Nu]*\[Chi]S)/(384*Sqrt[42])) + 
        e^6*((\[Nu]*\[Chi]S)/(12*Sqrt[42]*E^((2*I)*\[Xi])) - 
          (E^((2*I)*\[Xi])*\[Nu]*\[Chi]S)/(12*Sqrt[42]) - 
          (32*Sqrt[2/21]*\[Nu]*\[Chi]S)/(15*E^((4*I)*\[Xi])) + 
          (32*Sqrt[2/21]*E^((4*I)*\[Xi])*\[Nu]*\[Chi]S)/15 + 
          (81*Sqrt[3/14]*\[Nu]*\[Chi]S)/(20*E^((6*I)*\[Xi])) - 
          (81*Sqrt[3/14]*E^((6*I)*\[Xi])*\[Nu]*\[Chi]S)/20)) + 
      x^(5/2)*\[Epsilon]^5*(e^4*((8*\[Delta]*(-1 + \[Nu])*\[Chi]A + 
            (-8 + (65 - 33*\[Nu])*\[Nu])*\[Chi]S)/(9*Sqrt[42]*
            E^((2*I)*\[Xi])) - (4*Sqrt[2/21]*(\[Delta]*(-23 + 53*\[Nu])*
              \[Chi]A + (-23 + 7*(29 - 6*\[Nu])*\[Nu])*\[Chi]S))/
           (9*E^((4*I)*\[Xi])) + (4*Sqrt[2/21]*E^((4*I)*\[Xi])*
            (\[Delta]*(-23 + 53*\[Nu])*\[Chi]A + (-23 + 7*(29 - 6*\[Nu])*
                \[Nu])*\[Chi]S))/9 + (E^((2*I)*\[Xi])*
            (-8*\[Delta]*(-1 + \[Nu])*\[Chi]A + (8 + \[Nu]*(-65 + 33*\[Nu]))*
              \[Chi]S))/(9*Sqrt[42])) + 
        e^3*((-3*Sqrt[3/14]*(\[Delta]*(-74 + 177*\[Nu])*\[Chi]A + 
             (-74 + (649 - 136*\[Nu])*\[Nu])*\[Chi]S))/(64*E^((3*I)*\[Xi])) + 
          (3*Sqrt[3/14]*E^((3*I)*\[Xi])*(\[Delta]*(-74 + 177*\[Nu])*\[Chi]A + 
             (-74 + (649 - 136*\[Nu])*\[Nu])*\[Chi]S))/64 + 
          (\[Delta]*(86 - 303*\[Nu])*\[Chi]A + (86 - \[Nu]*(823 + 40*\[Nu]))*
             \[Chi]S)/(192*Sqrt[42]*E^(I*\[Xi])) + 
          (E^(I*\[Xi])*(\[Delta]*(-86 + 303*\[Nu])*\[Chi]A + 
             (-86 + \[Nu]*(823 + 40*\[Nu]))*\[Chi]S))/(192*Sqrt[42])) + 
        e^2*((E^((2*I)*\[Xi])*(\[Delta]*(-28 + 71*\[Nu])*\[Chi]A + 
             (-28 + (243 - 52*\[Nu])*\[Nu])*\[Chi]S))/(6*Sqrt[42]) + 
          (\[Delta]*(28 - 71*\[Nu])*\[Chi]A + (28 + \[Nu]*(-243 + 52*\[Nu]))*
             \[Chi]S)/(6*Sqrt[42]*E^((2*I)*\[Xi]))) + 
        e*((E^(I*\[Xi])*(\[Delta]*(-38 + 107*\[Nu])*\[Chi]A + 
             (-38 + (323 - 72*\[Nu])*\[Nu])*\[Chi]S))/(24*Sqrt[42]) + 
          (\[Delta]*(38 - 107*\[Nu])*\[Chi]A + (38 + \[Nu]*(-323 + 72*\[Nu]))*
             \[Chi]S)/(24*Sqrt[42]*E^(I*\[Xi]))) + 
        e^6*((4*Sqrt[2/21]*(\[Delta]*(-103 + 207*\[Nu])*\[Chi]A + 
             (-103 + (908 - 217*\[Nu])*\[Nu])*\[Chi]S))/
           (45*E^((4*I)*\[Xi])) - (4*Sqrt[2/21]*E^((4*I)*\[Xi])*
            (\[Delta]*(-103 + 207*\[Nu])*\[Chi]A + (-103 + (908 - 217*\[Nu])*
                \[Nu])*\[Chi]S))/45 + (E^((2*I)*\[Xi])*
            (\[Delta]*(-224 + 621*\[Nu])*\[Chi]A + (-224 + (2041 - 212*\[Nu])*
                \[Nu])*\[Chi]S))/(144*Sqrt[42]) - 
          (27*Sqrt[3/14]*(\[Delta]*(-64 + 141*\[Nu])*\[Chi]A + 
             (-64 + (569 - 116*\[Nu])*\[Nu])*\[Chi]S))/(80*E^((6*I)*\[Xi])) + 
          (27*Sqrt[3/14]*E^((6*I)*\[Xi])*(\[Delta]*(-64 + 141*\[Nu])*
              \[Chi]A + (-64 + (569 - 116*\[Nu])*\[Nu])*\[Chi]S))/80 + 
          (\[Delta]*(224 - 621*\[Nu])*\[Chi]A + 
            (224 + \[Nu]*(-2041 + 212*\[Nu]))*\[Chi]S)/(144*Sqrt[42]*
            E^((2*I)*\[Xi]))) + 
        e^5*((E^(I*\[Xi])*(\[Delta]*(-3470 + 10327*\[Nu])*\[Chi]A + 
             (-3470 + (31887 - 2216*\[Nu])*\[Nu])*\[Chi]S))/(4608*Sqrt[42]) + 
          (3*Sqrt[3/14]*(\[Delta]*(-670 + 1263*\[Nu])*\[Chi]A + 
             (-670 + (5831 - 1608*\[Nu])*\[Nu])*\[Chi]S))/
           (1024*E^((3*I)*\[Xi])) - (3*Sqrt[3/14]*E^((3*I)*\[Xi])*
            (\[Delta]*(-670 + 1263*\[Nu])*\[Chi]A + 
             (-670 + (5831 - 1608*\[Nu])*\[Nu])*\[Chi]S))/1024 + 
          (3125*E^((5*I)*\[Xi])*(\[Delta]*(-110 + 247*\[Nu])*\[Chi]A - 
             5*(22 + 5*\[Nu]*(-39 + 8*\[Nu]))*\[Chi]S))/(9216*Sqrt[42]) + 
          (3125*(\[Delta]*(110 - 247*\[Nu])*\[Chi]A + 
             5*(22 + 5*\[Nu]*(-39 + 8*\[Nu]))*\[Chi]S))/(9216*Sqrt[42]*
            E^((5*I)*\[Xi])) + (\[Delta]*(3470 - 10327*\[Nu])*\[Chi]A + 
            (3470 + \[Nu]*(-31887 + 2216*\[Nu]))*\[Chi]S)/
           (4608*Sqrt[42]*E^(I*\[Xi])))) + x^3*\[Epsilon]^6*
       (e^2*(I*Sqrt[2/21]*E^((2*I)*\[Xi])*(\[Delta]*\[Chi]A + \[Chi]S + 
            (-5 + (2*I)*Pi)*\[Nu]*\[Chi]S) + 
          (Sqrt[2/21]*(I*\[Delta]*\[Chi]A + (I + (-5*I + 2*Pi)*\[Nu])*
              \[Chi]S))/E^((2*I)*\[Xi])) + 
        e^4*(((-I/3)*Sqrt[2/21]*(\[Delta]*\[Chi]A + \[Chi]S + 
             (-5 - (2*I)*Pi)*\[Nu]*\[Chi]S))/E^((2*I)*\[Xi]) + 
          (Sqrt[2/21]*E^((2*I)*\[Xi])*((-I)*\[Delta]*\[Chi]A + 
             (-I + (5*I + 2*Pi)*\[Nu])*\[Chi]S))/3 + 
          (16*Sqrt[2/21]*(I*\[Delta]*\[Chi]A + \[Chi]S*(I + 2*Pi*\[Nu] + I*
                \[Nu]*(-5 + Log[16]))))/(3*E^((4*I)*\[Xi])) + 
          ((16*I)/3)*Sqrt[2/21]*E^((4*I)*\[Xi])*(\[Delta]*\[Chi]A + \[Chi]S + 
            \[Nu]*\[Chi]S*(-5 + (2*I)*Pi + Log[16]))) + 
        e*(((I/2)*E^(I*\[Xi])*(\[Delta]*\[Chi]A + \[Chi]S + 
             \[Nu]*\[Chi]S*(-5 + (2*I)*Pi - 4*Log[2])))/Sqrt[42] + 
          (I*\[Delta]*\[Chi]A + \[Chi]S*(I + 2*Pi*\[Nu] - I*\[Nu]*(5 + 
                Log[16])))/(2*Sqrt[42]*E^(I*\[Xi]))) + 
        e^3*(((-I/16)*(\[Delta]*\[Chi]A + \[Chi]S + \[Nu]*\[Chi]S*
              (-5 - (2*I)*Pi - 4*Log[2])))/(Sqrt[42]*E^(I*\[Xi])) + 
          (27*Sqrt[3/14]*(2*Pi*\[Nu]*\[Chi]S + I*(\[Delta]*\[Chi]A + 
               \[Chi]S + \[Nu]*\[Chi]S*(-5 + Log[81/16]))))/
           (16*E^((3*I)*\[Xi])) + ((27*I)/16)*Sqrt[3/14]*E^((3*I)*\[Xi])*
           (\[Delta]*\[Chi]A + \[Chi]S + \[Nu]*\[Chi]S*(-5 + (2*I)*Pi + 
              Log[81/16])) + (E^(I*\[Xi])*((-I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(-I + 2*Pi*\[Nu] + I*\[Nu]*(5 + Log[16]))))/
           (16*Sqrt[42])) + e^5*(((I/384)*E^(I*\[Xi])*(\[Delta]*\[Chi]A + 
             \[Chi]S + \[Nu]*\[Chi]S*(-5 + (2*I)*Pi - 4*Log[2])))/Sqrt[42] + 
          (243*Sqrt[3/14]*E^((3*I)*\[Xi])*(2*Pi*\[Nu]*\[Chi]S - 
             I*(\[Delta]*\[Chi]A + \[Chi]S + \[Nu]*\[Chi]S*(-5 + 
                 Log[81/16]))))/256 - (((243*I)/256)*Sqrt[3/14]*
            (\[Delta]*\[Chi]A + \[Chi]S + \[Nu]*\[Chi]S*(-5 - (2*I)*Pi + Log[
                81/16])))/E^((3*I)*\[Xi]) + (I*\[Delta]*\[Chi]A + 
            \[Chi]S*(I + 2*Pi*\[Nu] - I*\[Nu]*(5 + Log[16])))/
           (384*Sqrt[42]*E^(I*\[Xi])) + (15625*(I*\[Delta]*\[Chi]A + 
             \[Chi]S*(I + 2*Pi*\[Nu] + I*\[Nu]*(-5 + Log[625/16]))))/
           (768*Sqrt[42]*E^((5*I)*\[Xi])) + (((15625*I)/768)*E^((5*I)*\[Xi])*
            (\[Delta]*\[Chi]A + \[Chi]S + \[Nu]*\[Chi]S*(-5 + (2*I)*Pi + Log[
                625/16])))/Sqrt[42]) + 
        e^6*(((I/12)*E^((2*I)*\[Xi])*(\[Delta]*\[Chi]A + \[Chi]S + 
             (-5 + (2*I)*Pi)*\[Nu]*\[Chi]S))/Sqrt[42] + 
          (I*\[Delta]*\[Chi]A + (I + (-5*I + 2*Pi)*\[Nu])*\[Chi]S)/
           (12*Sqrt[42]*E^((2*I)*\[Xi])) + (64*Sqrt[2/21]*E^((4*I)*\[Xi])*
            ((-I)*\[Delta]*\[Chi]A + \[Chi]S*(-I + \[Nu]*(5*I + 2*Pi - 
                 (4*I)*Log[2]))))/15 - (((64*I)/15)*Sqrt[2/21]*
            (\[Delta]*\[Chi]A + \[Chi]S + \[Nu]*\[Chi]S*(-5 - (2*I)*Pi + Log[
                16])))/E^((4*I)*\[Xi]) + (243*Sqrt[3/14]*
            (I*\[Delta]*\[Chi]A + \[Chi]S*(I + 2*Pi*\[Nu] + I*\[Nu]*
                (-5 + Log[81]))))/(20*E^((6*I)*\[Xi])) + 
          ((243*I)/20)*Sqrt[3/14]*E^((6*I)*\[Xi])*(\[Delta]*\[Chi]A + 
            \[Chi]S + \[Nu]*\[Chi]S*(-5 + (2*I)*Pi + Log[81])))))
 
H\[Psi]InstTail[3, 1] = Sqrt[x]*(((I/12)*\[Delta])/Sqrt[14] + 
      e*((((-5*I)/12)*\[Delta])/(Sqrt[14]*E^(I*\[Xi])) + 
        ((I/12)*E^(I*\[Xi])*\[Delta])/Sqrt[14]) + 
      e^2*(((-I/4)*\[Delta])/Sqrt[14] - (((25*I)/48)*\[Delta])/
         (Sqrt[14]*E^((2*I)*\[Xi])) - ((I/48)*E^((2*I)*\[Xi])*\[Delta])/
         Sqrt[14]) + e^3*(((-I/96)*\[Delta])/(Sqrt[14]*E^(I*\[Xi])) - 
        (((19*I)/96)*E^(I*\[Xi])*\[Delta])/Sqrt[14] - (((65*I)/96)*\[Delta])/
         (Sqrt[14]*E^((3*I)*\[Xi])) - (((11*I)/96)*E^((3*I)*\[Xi])*\[Delta])/
         Sqrt[14]) + e^4*(((-I/16)*\[Delta])/Sqrt[14] + 
        (((43*I)/288)*\[Delta])/(Sqrt[14]*E^((2*I)*\[Xi])) - 
        (((53*I)/288)*E^((2*I)*\[Xi])*\[Delta])/Sqrt[14] - 
        (((515*I)/576)*\[Delta])/(Sqrt[14]*E^((4*I)*\[Xi])) - 
        (((131*I)/576)*E^((4*I)*\[Xi])*\[Delta])/Sqrt[14]) + 
      e^5*((((-73*I)/2304)*\[Delta])/(Sqrt[14]*E^(I*\[Xi])) - 
        (((43*I)/2304)*E^(I*\[Xi])*\[Delta])/Sqrt[14] + 
        (((197*I)/512)*\[Delta])/(Sqrt[14]*E^((3*I)*\[Xi])) - 
        (((73*I)/512)*E^((3*I)*\[Xi])*\[Delta])/Sqrt[14] - 
        (((5485*I)/4608)*\[Delta])/(Sqrt[14]*E^((5*I)*\[Xi])) - 
        (((1735*I)/4608)*E^((5*I)*\[Xi])*\[Delta])/Sqrt[14]) + 
      e^6*(((-I/32)*\[Delta])/Sqrt[14] - (((109*I)/2304)*\[Delta])/
         (Sqrt[14]*E^((2*I)*\[Xi])) + ((5*I)/2304)*Sqrt[7/2]*E^((2*I)*\[Xi])*
         \[Delta] + (((2159*I)/2880)*\[Delta])/(Sqrt[14]*E^((4*I)*\[Xi])) - 
        (((29*I)/576)*E^((4*I)*\[Xi])*\[Delta])/Sqrt[14] - 
        (((1223*I)/768)*\[Delta])/(Sqrt[14]*E^((6*I)*\[Xi])) - 
        (((2227*I)/3840)*E^((6*I)*\[Xi])*\[Delta])/Sqrt[14]))*\[Epsilon] + 
    x^(3/2)*\[Epsilon]^3*(e^2*(((I/288)*\[Delta]*(1301 - 340*\[Nu]))/
         (Sqrt[14]*E^((2*I)*\[Xi])) + ((I/144)*\[Delta]*(79 - 26*\[Nu]))/
         Sqrt[14] + ((I/288)*E^((2*I)*\[Xi])*\[Delta]*(421 - 8*\[Nu]))/
         Sqrt[14]) + e*(((I/72)*\[Delta]*(154 - 41*\[Nu]))/
         (Sqrt[14]*E^(I*\[Xi])) - (I/72)*Sqrt[7/2]*E^(I*\[Xi])*\[Delta]*
         (-2 + \[Nu])) - ((I/18)*\[Delta]*(4 + \[Nu]))/Sqrt[14] + 
      e^3*(((I/576)*\[Delta]*(4664 - 1183*\[Nu]))/
         (Sqrt[14]*E^((3*I)*\[Xi])) + ((I/576)*\[Delta]*(848 - 139*\[Nu]))/
         (Sqrt[14]*E^(I*\[Xi])) + ((I/576)*E^((3*I)*\[Xi])*\[Delta]*
          (2144 + 23*\[Nu]))/Sqrt[14] - ((I/576)*E^(I*\[Xi])*\[Delta]*
          (56 + 173*\[Nu]))/Sqrt[14]) + 
      e^4*(((I/864)*\[Delta]*(11674 - 2885*\[Nu]))/
         (Sqrt[14]*E^((4*I)*\[Xi])) + ((I/288)*\[Delta]*(347 - 88*\[Nu]))/
         Sqrt[14] + ((I/864)*E^((4*I)*\[Xi])*\[Delta]*(6494 + 101*\[Nu]))/
         Sqrt[14] + ((I/1728)*\[Delta]*(-109 + 158*\[Nu]))/
         (Sqrt[14]*E^((2*I)*\[Xi])) - ((I/1728)*E^((2*I)*\[Xi])*\[Delta]*
          (2621 + 602*\[Nu]))/Sqrt[14]) + 
      e^5*(((I/27648)*\[Delta]*(598310 - 144661*\[Nu]))/
         (Sqrt[14]*E^((5*I)*\[Xi])) + ((I/13824)*E^(I*\[Xi])*\[Delta]*
          (14954 - 5407*\[Nu]))/Sqrt[14] + ((I/13824)*\[Delta]*
          (28654 - 5345*\[Nu]))/(Sqrt[14]*E^(I*\[Xi])) + 
        ((I/9216)*Sqrt[7/2]*\[Delta]*(-4894 + 1205*\[Nu]))/E^((3*I)*\[Xi]) - 
        ((I/9216)*E^((3*I)*\[Xi])*\[Delta]*(46390 + 3859*\[Nu]))/Sqrt[14] + 
        ((I/27648)*E^((5*I)*\[Xi])*\[Delta]*(377410 + 5749*\[Nu]))/
         Sqrt[14]) + e^6*(((I/23040)*\[Delta]*(777397 - 184508*\[Nu]))/
         (Sqrt[14]*E^((6*I)*\[Xi])) + ((I/13824)*\[Delta]*
          (32431 - 6548*\[Nu]))/(Sqrt[14]*E^((2*I)*\[Xi])) + 
        ((I/13824)*E^((2*I)*\[Xi])*\[Delta]*(25799 - 5824*\[Nu]))/Sqrt[14] - 
        (((5*I)/576)*\[Delta]*(-161 + 40*\[Nu]))/Sqrt[14] + 
        ((I/23040)*E^((6*I)*\[Xi])*\[Delta]*(537277 + 7264*\[Nu]))/Sqrt[14] - 
        ((I/34560)*E^((4*I)*\[Xi])*\[Delta]*(422443 + 17758*\[Nu]))/
         Sqrt[14] + ((I/34560)*\[Delta]*(-384443 + 88954*\[Nu]))/
         (Sqrt[14]*E^((4*I)*\[Xi])))) + 
    SO*(x^2*\[Epsilon]^4*(((I/24)*((-4 + 11*\[Nu])*\[Chi]A + 
           \[Delta]*(-4 + 13*\[Nu])*\[Chi]S))/Sqrt[14] + 
        e*(((I/24)*E^(I*\[Xi])*((-14 + 51*\[Nu])*\[Chi]A + 
             \[Delta]*(-14 + 17*\[Nu])*\[Chi]S))/Sqrt[14] + 
          ((I/24)*((-50 + 181*\[Nu])*\[Chi]A + \[Delta]*(-50 + 63*\[Nu])*
              \[Chi]S))/(Sqrt[14]*E^(I*\[Xi]))) + 
        e^2*(((I/48)*E^((2*I)*\[Xi])*((-110 + 431*\[Nu])*\[Chi]A + 
             5*\[Delta]*(-22 + 17*\[Nu])*\[Chi]S))/Sqrt[14] + 
          ((I/24)*((-24 + 83*\[Nu])*\[Chi]A + \[Delta]*(-24 + 37*\[Nu])*
              \[Chi]S))/Sqrt[14] + ((I/48)*((-146 + 505*\[Nu])*\[Chi]A + 
             \[Delta]*(-146 + 243*\[Nu])*\[Chi]S))/(Sqrt[14]*
            E^((2*I)*\[Xi]))) + 
        e^3*(((I/192)*E^(I*\[Xi])*((-58 + 105*\[Nu])*\[Chi]A + 
             \[Delta]*(-58 + 227*\[Nu])*\[Chi]S))/Sqrt[14] + 
          ((I/192)*((-502 + 1895*\[Nu])*\[Chi]A + \[Delta]*(-502 + 421*\[Nu])*
              \[Chi]S))/(Sqrt[14]*E^(I*\[Xi])) + ((I/192)*E^((3*I)*\[Xi])*
            ((-922 + 3665*\[Nu])*\[Chi]A + \[Delta]*(-922 + 643*\[Nu])*
              \[Chi]S))/Sqrt[14] + ((I/192)*((-886 + 2975*\[Nu])*\[Chi]A + 
             \[Delta]*(-886 + 1717*\[Nu])*\[Chi]S))/(Sqrt[14]*
            E^((3*I)*\[Xi]))) + e^4*(((-I/144)*E^((2*I)*\[Xi])*
            ((11 + 80*\[Nu])*\[Chi]A + \[Delta]*(11 - 174*\[Nu])*\[Chi]S))/
           Sqrt[14] + ((I/64)*((-100 + 353*\[Nu])*\[Chi]A + 
             5*\[Delta]*(-20 + 27*\[Nu])*\[Chi]S))/Sqrt[14] + 
          ((I/144)*((-317 + 1221*\[Nu])*\[Chi]A + \[Delta]*(-317 + 173*\[Nu])*
              \[Chi]S))/(Sqrt[14]*E^((2*I)*\[Xi])) + 
          ((I/1152)*E^((4*I)*\[Xi])*((-9956 + 39855*\[Nu])*\[Chi]A + 
             \[Delta]*(-9956 + 6617*\[Nu])*\[Chi]S))/Sqrt[14] + 
          ((I/1152)*((-8084 + 26599*\[Nu])*\[Chi]A + \[Delta]*
              (-8084 + 17217*\[Nu])*\[Chi]S))/(Sqrt[14]*E^((4*I)*\[Xi]))) + 
        e^6*(((I/96)*((-194 + 691*\[Nu])*\[Chi]A + \[Delta]*(-194 + 245*
                \[Nu])*\[Chi]S))/Sqrt[14] - ((I/2880)*E^((4*I)*\[Xi])*
            ((-15730 + 68073*\[Nu])*\[Chi]A + 13*\[Delta]*(-1210 + 323*\[Nu])*
              \[Chi]S))/Sqrt[14] + ((I/1152)*((-3835 + 14241*\[Nu])*\[Chi]A + 
             5*\[Delta]*(-767 + 745*\[Nu])*\[Chi]S))/(Sqrt[14]*
            E^((2*I)*\[Xi])) + ((I/1152)*E^((2*I)*\[Xi])*
            ((-1633 + 5552*\[Nu])*\[Chi]A + \[Delta]*(-1633 + 2386*\[Nu])*
              \[Chi]S))/Sqrt[14] - ((I/2880)*((-2602 + 2933*\[Nu])*\[Chi]A + 
             \[Delta]*(-2602 + 22891*\[Nu])*\[Chi]S))/
           (Sqrt[14]*E^((4*I)*\[Xi])) + ((I/1920)*E^((6*I)*\[Xi])*
            ((-44813 + 180742*\[Nu])*\[Chi]A + \[Delta]*(-44813 + 28364*
                \[Nu])*\[Chi]S))/Sqrt[14] + 
          ((I/1920)*((-30575 + 98107*\[Nu])*\[Chi]A + \[Delta]*
              (-30575 + 72611*\[Nu])*\[Chi]S))/(Sqrt[14]*E^((6*I)*\[Xi]))) + 
        e^5*(((I/3072)*E^((3*I)*\[Xi])*((4698 - 22489*\[Nu])*\[Chi]A + 
             \[Delta]*(4698 + 1573*\[Nu])*\[Chi]S))/Sqrt[14] + 
          ((I/3072)*((-3834 + 17057*\[Nu])*\[Chi]A - \[Delta]*
              (3834 + 4997*\[Nu])*\[Chi]S))/(Sqrt[14]*E^((3*I)*\[Xi])) + 
          ((I/4608)*E^(I*\[Xi])*(47*(-74 + 217*\[Nu])*\[Chi]A + 
             \[Delta]*(-3478 + 7629*\[Nu])*\[Chi]S))/Sqrt[14] + 
          ((I/4608)*((-15754 + 59233*\[Nu])*\[Chi]A + \[Delta]*
              (-15754 + 13907*\[Nu])*\[Chi]S))/(Sqrt[14]*E^(I*\[Xi])) + 
          ((I/9216)*E^((5*I)*\[Xi])*((-133606 + 537199*\[Nu])*\[Chi]A + 
             \[Delta]*(-133606 + 86261*\[Nu])*\[Chi]S))/Sqrt[14] + 
          ((I/9216)*((-97786 + 317129*\[Nu])*\[Chi]A + \[Delta]*
              (-97786 + 221931*\[Nu])*\[Chi]S))/(Sqrt[14]*
            E^((5*I)*\[Xi])))) + x^3*\[Epsilon]^6*
       (((-I/216)*((-70 + \[Nu]*(-59 + 931*\[Nu]))*\[Chi]A + 
           \[Delta]*(-70 + 9*\[Nu]*(11 + 5*\[Nu]))*\[Chi]S))/Sqrt[14] + 
        e*(((-I/432)*E^(I*\[Xi])*((1816 + \[Nu]*(-8086 + 4193*\[Nu]))*
              \[Chi]A + \[Delta]*(1816 - 3042*\[Nu] + 339*\[Nu]^2)*\[Chi]S))/
           Sqrt[14] - ((I/432)*((-1336 + \[Nu]*(-62 + 13327*\[Nu]))*\[Chi]A + 
             \[Delta]*(-1336 + 3*\[Nu]*(-30 + 607*\[Nu]))*\[Chi]S))/
           (Sqrt[14]*E^(I*\[Xi]))) + 
        e^2*(((-I/864)*((-496 + \[Nu]*(-6125 + 27472*\[Nu]))*\[Chi]A + 
             \[Delta]*(-496 - 3843*\[Nu] + 4416*\[Nu]^2)*\[Chi]S))/Sqrt[14] - 
          ((I/1728)*E^((2*I)*\[Xi])*((24488 + \[Nu]*(-109813 + 53794*\[Nu]))*
              \[Chi]A + \[Delta]*(24488 + \[Nu]*(-34507 + 4918*\[Nu]))*
              \[Chi]S))/Sqrt[14] - 
          ((I/1728)*((-10520 + \[Nu]*(-20713 + 145814*\[Nu]))*\[Chi]A + 
             \[Delta]*(-10520 + 3*\[Nu]*(4131 + 9670*\[Nu]))*\[Chi]S))/
           (Sqrt[14]*E^((2*I)*\[Xi]))) + 
        e^3*(((-I/3456)*E^((3*I)*\[Xi])*((127084 + \[Nu]*(-573184 + 
                 257669*\[Nu]))*\[Chi]A + \[Delta]*(127084 + 9*\[Nu]*
                (-15536 + 2823*\[Nu]))*\[Chi]S))/Sqrt[14] - 
          ((I/3456)*E^(I*\[Xi])*((21404 + \[Nu]*(-117940 + 142069*\[Nu]))*
              \[Chi]A + \[Delta]*(21404 + \[Nu]*(-60748 + 23119*\[Nu]))*
              \[Chi]S))/Sqrt[14] - 
          ((I/3456)*((48100 + \[Nu]*(-288024 + 322763*\[Nu]))*\[Chi]A + 
             \[Delta]*(48100 + \[Nu]*(-119992 + 43921*\[Nu]))*\[Chi]S))/
           (Sqrt[14]*E^(I*\[Xi])) - 
          ((I/3456)*((-23852 + \[Nu]*(-195828 + 655931*\[Nu]))*\[Chi]A + 
             \[Delta]*(-23852 + \[Nu]*(78068 + 152545*\[Nu]))*\[Chi]S))/
           (Sqrt[14]*E^((3*I)*\[Xi]))) + 
        e^4*(((-I/3456)*((26876 + \[Nu]*(-192959 + 300106*\[Nu]))*\[Chi]A + 
             \[Delta]*(26876 - 79497*\[Nu] + 51270*\[Nu]^2)*\[Chi]S))/
           Sqrt[14] - ((I/5184)*E^((2*I)*\[Xi])*
            ((76522 + \[Nu]*(-376273 + 389342*\[Nu]))*\[Chi]A + 
             \[Delta]*(76522 + 3*\[Nu]*(-82797 + 19322*\[Nu]))*\[Chi]S))/
           Sqrt[14] - ((I/5184)*((58790 + \[Nu]*(-395483 + 583318*\[Nu]))*
              \[Chi]A + \[Delta]*(58790 + 3*\[Nu]*(-80223 + 24490*\[Nu]))*
              \[Chi]S))/(Sqrt[14]*E^((2*I)*\[Xi])) - 
          ((I/20736)*E^((4*I)*\[Xi])*((1688392 + \[Nu]*(-7648975 + 
                 3229844*\[Nu]))*\[Chi]A + \[Delta]*(1688392 + 3*\[Nu]*
                (-494051 + 110092*\[Nu]))*\[Chi]S))/Sqrt[14] - 
          ((I/20736)*((-103120 + \[Nu]*(-3160295 + 7853824*\[Nu]))*\[Chi]A + 
             \[Delta]*(-103120 + 3*\[Nu]*(382309 + 668904*\[Nu]))*\[Chi]S))/
           (Sqrt[14]*E^((4*I)*\[Xi]))) + 
        e^5*(((-I/55296)*E^((3*I)*\[Xi])*((702160 + \[Nu]*(-3386206 + 
                 5445809*\[Nu]))*\[Chi]A + \[Delta]*(702160 - 4843386*\[Nu] + 
               844851*\[Nu]^2)*\[Chi]S))/Sqrt[14] - 
          ((I/55296)*((641872 + 5*\[Nu]*(-666294 + 1116787*\[Nu]))*\[Chi]A + 
             7*\[Delta]*(91696 + \[Nu]*(-658534 + 16411*\[Nu]))*\[Chi]S))/
           (Sqrt[14]*E^((3*I)*\[Xi])) - ((I/82944)*E^(I*\[Xi])*
            ((999424 + \[Nu]*(-5885098 + 7109825*\[Nu]))*\[Chi]A + 
             \[Delta]*(999424 + 3*\[Nu]*(-781706 + 435537*\[Nu]))*\[Chi]S))/
           Sqrt[14] - ((I/82944)*((3241376 + 331*\[Nu]*(-52702 + 
                 47245*\[Nu]))*\[Chi]A + \[Delta]*(3241376 + \[Nu]*
                (-6393214 + 2235613*\[Nu]))*\[Chi]S))/
           (Sqrt[14]*E^(I*\[Xi])) - ((I/165888)*E^((5*I)*\[Xi])*
            ((26996240 + \[Nu]*(-122763778 + 49638337*\[Nu]))*\[Chi]A + 
             \[Delta]*(26996240 + \[Nu]*(-19262566 + 5164771*\[Nu]))*
              \[Chi]S))/Sqrt[14] - ((I/165888)*
            ((189520 + \[Nu]*(-54823378 + 116363087*\[Nu]))*\[Chi]A + 
             \[Delta]*(189520 + 3*\[Nu]*(6647758 + 10595727*\[Nu]))*\[Chi]S))/
           (Sqrt[14]*E^((5*I)*\[Xi]))) + 
        e^6*(((-I/1152)*((24572 + \[Nu]*(-150145 + 184368*\[Nu]))*\[Chi]A + 
             \[Delta]*(24572 + \[Nu]*(-57539 + 31696*\[Nu]))*\[Chi]S))/
           Sqrt[14] - ((I/20736)*E^((2*I)*\[Xi])*
            ((621907 + \[Nu]*(-3206184 + 2854955*\[Nu]))*\[Chi]A + 
             \[Delta]*(621907 + 4*\[Nu]*(-299821 + 116938*\[Nu]))*\[Chi]S))/
           Sqrt[14] - ((I/34560)*E^((6*I)*\[Xi])*
            ((10488521 + 14*\[Nu]*(-3418369 + 1340612*\[Nu]))*\[Chi]A + 
             \[Delta]*(10488521 + 3*\[Nu]*(-2045282 + 655949*\[Nu]))*
              \[Chi]S))/Sqrt[14] - ((I/20736)*
            ((688049 + \[Nu]*(-4138573 + 4585507*\[Nu]))*\[Chi]A + 
             \[Delta]*(688049 + \[Nu]*(-1531903 + 725764*\[Nu]))*\[Chi]S))/
           (Sqrt[14]*E^((2*I)*\[Xi])) - ((I/103680)*E^((4*I)*\[Xi])*
            ((-1469794 + \[Nu]*(6949097 + 9369433*\[Nu]))*\[Chi]A + 
             \[Delta]*(-1469794 + \[Nu]*(-15025045 + 1763959*\[Nu]))*
              \[Chi]S))/Sqrt[14] - ((I/103680)*
            ((1343242 + \[Nu]*(341793 + 79991*\[Nu]))*\[Chi]A + 
             \[Delta]*(1343242 - \[Nu]*(17521837 + 4631351*\[Nu]))*\[Chi]S))/
           (Sqrt[14]*E^((4*I)*\[Xi])) - 
          ((I/34560)*((475403 + \[Nu]*(-22116559 + 42693334*\[Nu]))*\[Chi]A + 
             \[Delta]*(475403 + \[Nu]*(8354795 + 12282157*\[Nu]))*\[Chi]S))/
           (Sqrt[14]*E^((6*I)*\[Xi]))))) + SO^2*x^(5/2)*\[Epsilon]^5*
     (((I/24)*(\[Kappa]A*(-5 + 4*\[Nu]) + 2*(-5 + 4*\[Nu])*\[Chi]A*\[Chi]S - 
         \[Delta]*(\[Kappa]S*(5 + 6*\[Nu]) + (5 + 12*\[Nu])*\[Chi]A^2 + 
           5*\[Chi]S^2)))/Sqrt[14] + 
      e^2*(((-I/48)*E^((2*I)*\[Xi])*(-27*\[Delta]*\[Kappa]S + 
           \[Kappa]A*(-27 + 140*\[Nu]) + 86*\[Delta]*\[Nu]*
            (\[Kappa]S + 2*\[Chi]A^2) + 2*(-27 + 140*\[Nu])*\[Chi]A*\[Chi]S - 
           27*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/Sqrt[14] - 
        ((I/48)*(19*\[Delta]*\[Kappa]S + \[Kappa]A*(19 + 4*\[Nu]) + 
           42*\[Delta]*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 2*(19 + 4*\[Nu])*
            \[Chi]A*\[Chi]S + 19*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/
         Sqrt[14] - ((I/48)*(47*\[Delta]*\[Kappa]S + \[Kappa]A*
            (47 + 4*\[Nu]) + 98*\[Delta]*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
           2*(47 + 4*\[Nu])*\[Chi]A*\[Chi]S + 47*\[Delta]*
            (\[Chi]A^2 + \[Chi]S^2)))/(Sqrt[14]*E^((2*I)*\[Xi]))) + 
      e^3*(((-I/192)*E^((3*I)*\[Xi])*(-283*\[Delta]*\[Kappa]S + 
           \[Kappa]A*(-283 + 1244*\[Nu]) + 678*\[Delta]*\[Nu]*
            (\[Kappa]S + 2*\[Chi]A^2) + 2*(-283 + 1244*\[Nu])*\[Chi]A*
            \[Chi]S - 283*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/Sqrt[14] - 
        ((I/192)*(-41*\[Delta]*\[Kappa]S + \[Kappa]A*(-41 + 628*\[Nu]) + 
           546*\[Delta]*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 2*(-41 + 628*\[Nu])*
            \[Chi]A*\[Chi]S - 41*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/
         (Sqrt[14]*E^(I*\[Xi])) - ((I/192)*E^(I*\[Xi])*
          (127*\[Delta]*\[Kappa]S + \[Kappa]A*(127 - 140*\[Nu]) + 
           114*\[Delta]*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 2*(127 - 140*\[Nu])*
            \[Chi]A*\[Chi]S + 127*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/
         Sqrt[14] - ((I/192)*(413*\[Delta]*\[Kappa]S + 
           \[Kappa]A*(413 - 292*\[Nu]) + 534*\[Delta]*\[Nu]*
            (\[Kappa]S + 2*\[Chi]A^2) + 2*(413 - 292*\[Nu])*\[Chi]A*\[Chi]S + 
           413*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/
         (Sqrt[14]*E^((3*I)*\[Xi]))) + 
      e^4*(((-I/1152)*E^((4*I)*\[Xi])*(-3347*\[Delta]*\[Kappa]S + 
           \[Kappa]A*(-3347 + 13788*\[Nu]) + 7094*\[Delta]*\[Nu]*
            (\[Kappa]S + 2*\[Chi]A^2) + 2*(-3347 + 13788*\[Nu])*\[Chi]A*
            \[Chi]S - 3347*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/Sqrt[14] - 
        ((I/288)*(-19*\[Delta]*\[Kappa]S + \[Kappa]A*(-19 + 876*\[Nu]) + 
           838*\[Delta]*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 2*(-19 + 876*\[Nu])*
            \[Chi]A*\[Chi]S - 19*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/
         (Sqrt[14]*E^((2*I)*\[Xi])) - ((I/48)*(25*\[Delta]*\[Kappa]S + 
           \[Kappa]A*(25 + 28*\[Nu]) + 78*\[Delta]*\[Nu]*(\[Kappa]S + 
             2*\[Chi]A^2) + 2*(25 + 28*\[Nu])*\[Chi]A*\[Chi]S + 
           25*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/Sqrt[14] - 
        ((I/288)*E^((2*I)*\[Xi])*(181*\[Delta]*\[Kappa]S + 
           \[Kappa]A*(181 - 52*\[Nu]) + 310*\[Delta]*\[Nu]*
            (\[Kappa]S + 2*\[Chi]A^2) + 2*(181 - 52*\[Nu])*\[Chi]A*\[Chi]S + 
           181*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/Sqrt[14] - 
        ((I/1152)*(4565*\[Delta]*\[Kappa]S + \[Kappa]A*(4565 - 4676*\[Nu]) + 
           4454*\[Delta]*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
           2*(4565 - 4676*\[Nu])*\[Chi]A*\[Chi]S + 4565*\[Delta]*
            (\[Chi]A^2 + \[Chi]S^2)))/(Sqrt[14]*E^((4*I)*\[Xi]))) + 
      e^6*(((I/11520)*E^((4*I)*\[Xi])*(-25985*\[Delta]*\[Kappa]S + 
           \[Kappa]A*(-25985 + 56596*\[Nu]) + 4626*\[Delta]*\[Nu]*
            (\[Kappa]S + 2*\[Chi]A^2) + 2*(-25985 + 56596*\[Nu])*\[Chi]A*
            \[Chi]S - 25985*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/Sqrt[14] - 
        ((I/11520)*(-17287*\[Delta]*\[Kappa]S + \[Kappa]A*
            (-17287 + 66956*\[Nu]) + 32382*\[Delta]*\[Nu]*
            (\[Kappa]S + 2*\[Chi]A^2) + 2*(-17287 + 66956*\[Nu])*\[Chi]A*
            \[Chi]S - 17287*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/
         (Sqrt[14]*E^((4*I)*\[Xi])) - ((I/1920)*E^((6*I)*\[Xi])*
          (-16433*\[Delta]*\[Kappa]S + \[Kappa]A*(-16433 + 63724*\[Nu]) + 
           30858*\[Delta]*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
           2*(-16433 + 63724*\[Nu])*\[Chi]A*\[Chi]S - 16433*\[Delta]*
            (\[Chi]A^2 + \[Chi]S^2)))/Sqrt[14] - 
        ((I/1152)*(-59*\[Delta]*\[Kappa]S + \[Kappa]A*(-59 + 4996*\[Nu]) + 
           4878*\[Delta]*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
           2*(-59 + 4996*\[Nu])*\[Chi]A*\[Chi]S - 59*\[Delta]*
            (\[Chi]A^2 + \[Chi]S^2)))/(Sqrt[14]*E^((2*I)*\[Xi])) - 
        ((I/96)*(59*\[Delta]*\[Kappa]S + \[Kappa]A*(59 + 116*\[Nu]) + 
           234*\[Delta]*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 2*(59 + 116*\[Nu])*
            \[Chi]A*\[Chi]S + 59*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/
         Sqrt[14] - ((I/288)*E^((2*I)*\[Xi])*(199*\[Delta]*\[Kappa]S + 
           \[Kappa]A*(199 + 142*\[Nu]) + 540*\[Delta]*\[Nu]*
            (\[Kappa]S + 2*\[Chi]A^2) + 2*(199 + 142*\[Nu])*\[Chi]A*\[Chi]S + 
           199*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/Sqrt[14] - 
        ((I/960)*(10505*\[Delta]*\[Kappa]S + \[Kappa]A*
            (10505 - 13672*\[Nu]) + 7338*\[Delta]*\[Nu]*(\[Kappa]S + 
             2*\[Chi]A^2) + 2*(10505 - 13672*\[Nu])*\[Chi]A*\[Chi]S + 
           10505*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/
         (Sqrt[14]*E^((6*I)*\[Xi]))) + 
      e^5*(((-I/9216)*E^((5*I)*\[Xi])*(-47319*\[Delta]*\[Kappa]S + 
           \[Kappa]A*(-47319 + 187948*\[Nu]) + 93310*\[Delta]*\[Nu]*
            (\[Kappa]S + 2*\[Chi]A^2) + 2*(-47319 + 187948*\[Nu])*\[Chi]A*
            \[Chi]S - 47319*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/Sqrt[14] - 
        ((I/4608)*(-1935*\[Delta]*\[Kappa]S + \[Kappa]A*
            (-1935 + 23404*\[Nu]) + 19534*\[Delta]*\[Nu]*(\[Kappa]S + 
             2*\[Chi]A^2) + 2*(-1935 + 23404*\[Nu])*\[Chi]A*\[Chi]S - 
           1935*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/(Sqrt[14]*E^(I*\[Xi])) - 
        ((I/3072)*(-1153*\[Delta]*\[Kappa]S + \[Kappa]A*
            (-1153 + 11444*\[Nu]) + 9138*\[Delta]*\[Nu]*(\[Kappa]S + 
             2*\[Chi]A^2) + 2*(-1153 + 11444*\[Nu])*\[Chi]A*\[Chi]S - 
           1153*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/
         (Sqrt[14]*E^((3*I)*\[Xi])) - ((I/3072)*E^((3*I)*\[Xi])*
          (3251*\[Delta]*\[Kappa]S + \[Kappa]A*(3251 - 3772*\[Nu]) + 
           2730*\[Delta]*\[Nu]*(\[Kappa]S + 2*\[Chi]A^2) + 
           2*(3251 - 3772*\[Nu])*\[Chi]A*\[Chi]S + 3251*\[Delta]*
            (\[Chi]A^2 + \[Chi]S^2)))/Sqrt[14] - 
        ((I/4608)*E^(I*\[Xi])*(4501*\[Delta]*\[Kappa]S + 
           \[Kappa]A*(4501 - 4420*\[Nu]) + 4582*\[Delta]*\[Nu]*
            (\[Kappa]S + 2*\[Chi]A^2) + 2*(4501 - 4420*\[Nu])*\[Chi]A*
            \[Chi]S + 4501*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/Sqrt[14] - 
        ((I/9216)*(62101*\[Delta]*\[Kappa]S + \[Kappa]A*
            (62101 - 74212*\[Nu]) + 49990*\[Delta]*\[Nu]*(\[Kappa]S + 
             2*\[Chi]A^2) + 2*(62101 - 74212*\[Nu])*\[Chi]A*\[Chi]S + 
           62101*\[Delta]*(\[Chi]A^2 + \[Chi]S^2)))/
         (Sqrt[14]*E^((5*I)*\[Xi]))) + 
      e*(((-I/24)*(\[Kappa]A*(5 + 28*\[Nu]) + \[Delta]*\[Kappa]S*
            (5 + 38*\[Nu]) + \[Delta]*(5 + 76*\[Nu])*\[Chi]A^2 + 
           2*(5 + 28*\[Nu])*\[Chi]A*\[Chi]S + 5*\[Delta]*\[Chi]S^2))/
         (Sqrt[14]*E^(I*\[Xi])) - ((I/24)*E^(I*\[Xi])*(\[Kappa]A + 
           12*\[Kappa]A*\[Nu] + 2*(1 + 12*\[Nu])*\[Chi]A*\[Chi]S + 
           \[Delta]*(\[Kappa]S + 14*\[Kappa]S*\[Nu] + \[Chi]A^2 + 
             28*\[Nu]*\[Chi]A^2 + \[Chi]S^2)))/Sqrt[14])) + 
    x^2*\[Epsilon]^4*
     (e^3*((E^(I*\[Xi])*\[Delta]*(-539 - (275*I)*Pi - 50*Log[2]))/
         (480*Sqrt[14]) + (E^((3*I)*\[Xi])*\[Delta]*(903 - (325*I)*Pi + 
           650*Log[2]))/(1440*Sqrt[14]) + 
        (\[Delta]*(-7287 - (5205*I)*Pi + 32480*Log[2] - 7290*Log[6]))/
         (1440*Sqrt[14]*E^((3*I)*\[Xi])) + 
        (\[Delta]*(-469 - (335*I)*Pi - 5020*Log[2] + 2430*Log[6]))/
         (480*Sqrt[14]*E^(I*\[Xi]))) + 
      e^4*((E^((2*I)*\[Xi])*\[Delta]*(-2352 - (745*I)*Pi - 160*Log[2]))/
         (1440*Sqrt[14]) + (\[Delta]*((-72540*I)*Pi - 
           18*(5642 + 24060*Log[2] + 405*Log[3]) + 296875*Log[5]))/
         (11520*Sqrt[14]*E^((4*I)*\[Xi])) + 
        (\[Delta]*(-812 - (580*I)*Pi + 59175*Log[2] - 20655*Log[6]))/
         (1440*Sqrt[14]*E^((2*I)*\[Xi])) + (E^((4*I)*\[Xi])*\[Delta]*
          (13804 - (6735*I)*Pi + 23105*Log[2] - 8505*Log[6]))/
         (11520*Sqrt[14]) + (\[Delta]*(-2436 - (1685*I)*Pi - 22735*Log[2] + 
           10935*Log[6]))/(1920*Sqrt[14])) + 
      e^5*((E^((3*I)*\[Xi])*\[Delta]*((-1479*I)*Pi - 
           7*(1653 + 622*Log[2] - 486*Log[3])))/(4608*Sqrt[14]) + 
        (E^((5*I)*\[Xi])*\[Delta]*(48979 - (27209*I)*Pi + 5202*Log[2] - 
           17010*Log[3]))/(23040*Sqrt[14]) + 
        (\[Delta]*(-334971 - (239265*I)*Pi + 567374*Log[2] + 1110024*Log[3] - 
           593750*Log[5]))/(23040*Sqrt[14]*E^((5*I)*\[Xi])) + 
        (\[Delta]*(4683 + (3345*I)*Pi - 331614*Log[2] + 42768*Log[3] + 
           118750*Log[5]))/(4608*Sqrt[14]*E^((3*I)*\[Xi])) + 
        (E^(I*\[Xi])*\[Delta]*(-3297 - (2171*I)*Pi - 32230*Log[2] + 
           15552*Log[6]))/(2304*Sqrt[14]) + 
        (\[Delta]*((-2361*I)*Pi + 104386*Log[2] - 37*(91 + 972*Log[6])))/
         (2304*Sqrt[14]*E^(I*\[Xi]))) + 
      e^6*((\[Delta]*((-1919*I)*Pi + 9*(-308 + 6856*Log[2] - 3645*Log[3])))/
         (1728*Sqrt[14]) + (E^((4*I)*\[Xi])*\[Delta]*(-57722 + (2350*I)*Pi - 
           1900*Log[2] + 13365*Log[3]))/(14400*Sqrt[14]) + 
        (E^((2*I)*\[Xi])*\[Delta]*(-70336 - (47965*I)*Pi - 408960*Log[2] + 
           383940*Log[3]))/(46080*Sqrt[14]) + (E^((6*I)*\[Xi])*\[Delta]*
          (4963392 - (2957035*I)*Pi + 3076480*Log[2] - 1148175*Log[3] - 
           2421875*Log[5]))/(1382400*Sqrt[14]) + 
        (\[Delta]*(68558 + (48970*I)*Pi + 1505380*Log[2] + 718065*Log[3] - 
           1171875*Log[5]))/(14400*Sqrt[14]*E^((4*I)*\[Xi])) + 
        (\[Delta]*(-159712 - (112595*I)*Pi - 7690560*Log[2] + 
           1332855*Log[3] + 2671875*Log[5]))/(92160*Sqrt[14]*
          E^((2*I)*\[Xi])) + (\[Delta]*(-32020128 - (22871520*I)*Pi - 
           115917760*Log[2] - 67676715*Log[3] - 4453125*Log[5] + 
           119413735*Log[7]))/(1382400*Sqrt[14]*E^((6*I)*\[Xi]))) + 
      (\[Delta]*(7 + (5*I)*Pi + Log[1024]))/(60*Sqrt[14]) + 
      e*((\[Delta]*(-63 - (45*I)*Pi - 10*Log[2]))/(60*Sqrt[14]*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*\[Delta]*(7 + (5*I)*Pi + Log[1024]))/(60*Sqrt[14])) + 
      e^2*((E^((2*I)*\[Xi])*\[Delta]*(14 - I*Pi + 20*Log[2]))/(48*Sqrt[14]) + 
        (\[Delta]*(-126 - (90*I)*Pi - 487*Log[2] + 243*Log[6]))/
         (48*Sqrt[14]*E^((2*I)*\[Xi])) + 
        (\[Delta]*(-49 - (35*I)*Pi + Log[1024]))/(60*Sqrt[14])))
 
H\[Psi]InstTail[3, 2] = 
   x*\[Epsilon]^2*(e^2*((Sqrt[5/7]*(1 - 3*\[Nu]))/6 + 
        (2*Sqrt[5/7]*(1 - 3*\[Nu]))/(3*E^((2*I)*\[Xi])) + 
        (5*Sqrt[5/7]*E^((2*I)*\[Xi])*(1 - 3*\[Nu]))/6) + 
      (Sqrt[5/7]*(1 - 3*\[Nu]))/3 + 
      e*((-11*Sqrt[5/7]*(-1 + 3*\[Nu]))/(24*E^(I*\[Xi])) - 
        (13*Sqrt[5/7]*E^(I*\[Xi])*(-1 + 3*\[Nu]))/24) + 
      e^3*((Sqrt[35]*E^(I*\[Xi])*(1 - 3*\[Nu]))/192 - 
        (17*Sqrt[5/7]*(-1 + 3*\[Nu]))/(192*E^(I*\[Xi])) - 
        (185*Sqrt[5/7]*(-1 + 3*\[Nu]))/(192*E^((3*I)*\[Xi])) - 
        (239*Sqrt[5/7]*E^((3*I)*\[Xi])*(-1 + 3*\[Nu]))/192) + 
      e^4*((Sqrt[5/7]*(1 - 3*\[Nu]))/8 + (Sqrt[35]*(-1 + 3*\[Nu]))/
         (72*E^((2*I)*\[Xi])) + (17*Sqrt[5/7]*E^((2*I)*\[Xi])*(-1 + 3*\[Nu]))/
         72 - (199*Sqrt[5/7]*(-1 + 3*\[Nu]))/(144*E^((4*I)*\[Xi])) - 
        (263*Sqrt[5/7]*E^((4*I)*\[Xi])*(-1 + 3*\[Nu]))/144) + 
      e^5*((-551*Sqrt[5/7]*(-1 + 3*\[Nu]))/(4608*E^(I*\[Xi])) - 
        (529*Sqrt[5/7]*E^(I*\[Xi])*(-1 + 3*\[Nu]))/4608 + 
        (1361*Sqrt[5/7]*(-1 + 3*\[Nu]))/(3072*E^((3*I)*\[Xi])) + 
        (2279*Sqrt[5/7]*E^((3*I)*\[Xi])*(-1 + 3*\[Nu]))/3072 - 
        (2593*Sqrt[35]*(-1 + 3*\[Nu]))/(9216*E^((5*I)*\[Xi])) - 
        (24401*Sqrt[5/7]*E^((5*I)*\[Xi])*(-1 + 3*\[Nu]))/9216) + 
      e^6*((-5*Sqrt[5/7]*(-1 + 3*\[Nu]))/48 - (83*Sqrt[5/7]*(-1 + 3*\[Nu]))/
         (576*E^((2*I)*\[Xi])) - (13*Sqrt[35]*E^((2*I)*\[Xi])*(-1 + 3*\[Nu]))/
         576 + (1507*(-1 + 3*\[Nu]))/(288*Sqrt[35]*E^((4*I)*\[Xi])) + 
        (2339*E^((4*I)*\[Xi])*(-1 + 3*\[Nu]))/(288*Sqrt[35]) - 
        (383*Sqrt[7/5]*(-1 + 3*\[Nu]))/(192*E^((6*I)*\[Xi])) - 
        (3653*E^((6*I)*\[Xi])*(-1 + 3*\[Nu]))/(192*Sqrt[35]))) + 
    x^2*\[Epsilon]^4*((-193 + 5*(145 - 73*\[Nu])*\[Nu])/(54*Sqrt[35]) + 
      e^2*((-374 + 5*(278 - 71*\[Nu])*\[Nu])/(54*Sqrt[35]) + 
        (-3637 + 5*(2461 - 61*\[Nu])*\[Nu])/(216*Sqrt[35]*E^((2*I)*\[Xi])) + 
        (E^((2*I)*\[Xi])*(971 - 5*\[Nu]*(531 + 97*\[Nu])))/(72*Sqrt[35])) + 
      e*((-3626 + 5*(2549 - 479*\[Nu])*\[Nu])/(432*Sqrt[35]*E^(I*\[Xi])) + 
        (E^(I*\[Xi])*(298 - 5*\[Nu]*(83 + 215*\[Nu])))/(144*Sqrt[35])) + 
      e^3*((E^(I*\[Xi])*(-41776 + 5*(29641 - 5971*\[Nu])*\[Nu]))/
         (3456*Sqrt[35]) + (-24836 + 5*(18971 - 5705*\[Nu])*\[Nu])/
         (3456*Sqrt[35]*E^(I*\[Xi])) + (E^((3*I)*\[Xi])*
          (119644 - 5*\[Nu]*(69691 + 2183*\[Nu])))/(3456*Sqrt[35]) + 
        (-106376 + 5*\[Nu]*(70199 + 5635*\[Nu]))/(3456*Sqrt[35]*
          E^((3*I)*\[Xi]))) + 
      e^4*((E^((2*I)*\[Xi])*(-34669 + 5*(23317 - 3607*\[Nu])*\[Nu]))/
         (1296*Sqrt[35]) + (-1303 + 5*(967 - 211*\[Nu])*\[Nu])/
         (144*Sqrt[35]) - (1773 - 8345*\[Nu] + 6675*\[Nu]^2)/
         (432*Sqrt[35]*E^((2*I)*\[Xi])) + 
        (E^((4*I)*\[Xi])*(61589 + 5*\[Nu]*(-36781 + 929*\[Nu])))/
         (864*Sqrt[35]) + (-137833 + 5*\[Nu]*(89329 + 13883*\[Nu]))/
         (2592*Sqrt[35]*E^((4*I)*\[Xi]))) + 
      e^6*((-59263 + 5*(43069 - 5638*\[Nu])*\[Nu])/(5184*Sqrt[35]*
          E^((2*I)*\[Xi])) + (E^((2*I)*\[Xi])*(-38093 + 
           35*(4289 - 1004*\[Nu])*\[Nu]))/(5184*Sqrt[35]) - 
        (E^((4*I)*\[Xi])*(1680173 - 5385385*\[Nu] + 668260*\[Nu]^2))/
         (12960*Sqrt[35]) - (Sqrt[5/7]*(929 + 5*\[Nu]*(-689 + 140*\[Nu])))/
         432 + (E^((6*I)*\[Xi])*(1996652 + 5*\[Nu]*(-1219622 + 90815*\[Nu])))/
         (8640*Sqrt[35]) + (428987 - 5*\[Nu]*(232547 + 191020*\[Nu]))/
         (12960*Sqrt[35]*E^((4*I)*\[Xi])) + 
        (-1236808 + 5*\[Nu]*(781882 + 201863*\[Nu]))/(8640*Sqrt[35]*
          E^((6*I)*\[Xi]))) + 
      e^5*(-(Sqrt[5/7]*(53642 + \[Nu]*(-199331 + 40489*\[Nu])))/
         (27648*E^(I*\[Xi])) + (Sqrt[5/7]*(76198 - 
           5*\[Nu]*(27559 + 74699*\[Nu])))/(55296*E^((3*I)*\[Xi])) - 
        (Sqrt[5/7]*E^(I*\[Xi])*(165694 + \[Nu]*(-614195 + 129929*\[Nu])))/
         82944 + (E^((5*I)*\[Xi])*(7317094 + 5*\[Nu]*(-4430175 + 
             245629*\[Nu])))/(55296*Sqrt[35]) - 
        (Sqrt[5/7]*E^((3*I)*\[Xi])*(668426 + \[Nu]*(-2177383 + 
             287941*\[Nu])))/55296 + 
        (-14672318 + 5*\[Nu]*(9378425 + 1994437*\[Nu]))/
         (165888*Sqrt[35]*E^((5*I)*\[Xi])))) + SO^2*x^3*\[Epsilon]^6*
     ((-2*Sqrt[5/7]*(3*\[Delta]*\[Kappa]A*(-1 + 2*\[Nu]) - 
         3*(\[Kappa]S + \[Chi]A^2) - 6*\[Nu]*(-2 + 3*\[Nu])*
          (\[Kappa]S + 2*\[Chi]A^2) + 2*\[Delta]*(-3 + 16*\[Nu])*\[Chi]A*
          \[Chi]S + (-3 + 4*(5 - 4*\[Nu])*\[Nu])*\[Chi]S^2))/9 + 
      e*((Sqrt[5/7]*E^(I*\[Xi])*(3*\[Delta]*\[Kappa]A*(40 - 91*\[Nu]) + 
           120*(\[Kappa]S + \[Chi]A^2) + 9*\[Nu]*
            (\[Kappa]S*(-57 + 80*\[Nu]) + (-103 + 160*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(120 - 743*\[Nu])*\[Chi]A*\[Chi]S + 
           (120 + \[Nu]*(-1039 + 584*\[Nu]))*\[Chi]S^2))/72 + 
        (Sqrt[5/7]*(3*\[Delta]*\[Kappa]A*(44 - 83*\[Nu]) + 
           132*(\[Kappa]S + \[Chi]A^2) + 9*\[Nu]*
            (\[Kappa]S*(-57 + 88*\[Nu]) + (-119 + 176*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(132 - 799*\[Nu])*\[Chi]A*\[Chi]S + 
           (132 + \[Nu]*(-1055 + 616*\[Nu]))*\[Chi]S^2))/(72*E^(I*\[Xi]))) + 
      e^2*((-2*Sqrt[5/7]*(12*\[Delta]*\[Kappa]A*(-1 + 2*\[Nu]) - 
           12*(\[Kappa]S + \[Chi]A^2) - 24*\[Nu]*(-2 + 3*\[Nu])*
            (\[Kappa]S + 2*\[Chi]A^2) + \[Delta]*(-24 + 143*\[Nu])*\[Chi]A*
            \[Chi]S - (-3 + 2*\[Nu])*(-4 + 29*\[Nu])*\[Chi]S^2))/9 + 
        (Sqrt[5/7]*E^((2*I)*\[Xi])*(3*(\[Delta]*\[Kappa]A*(79 - 191*\[Nu]) + 
             \[Kappa]S*(79 + \[Nu]*(-349 + 474*\[Nu])) + 
             (79 + \[Nu]*(-599 + 948*\[Nu]))*\[Chi]A^2) + 
           2*\[Delta]*(237 - 1549*\[Nu])*\[Chi]A*\[Chi]S + 
           (237 + \[Nu]*(-2249 + 1072*\[Nu]))*\[Chi]S^2))/72 + 
        (Sqrt[5/7]*(3*\[Delta]*\[Kappa]A*(89 - 161*\[Nu]) + 
           267*(\[Kappa]S + \[Chi]A^2) + 9*\[Nu]*
            (\[Kappa]S*(-113 + 178*\[Nu]) + (-243 + 356*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(267 - 1643*\[Nu])*\[Chi]A*\[Chi]S + 
           (267 + \[Nu]*(-2167 + 1136*\[Nu]))*\[Chi]S^2))/
         (72*E^((2*I)*\[Xi]))) + 
      e^4*((5*Sqrt[5/7]*(39*(\[Kappa]S + \[Delta]*(\[Kappa]A - 2*\[Kappa]A*
                \[Nu]) + \[Chi]A^2 + 2*\[Nu]*(-2 + 3*\[Nu])*(\[Kappa]S + 2*
                \[Chi]A^2)) + 2*\[Delta]*(39 - 238*\[Nu])*\[Chi]A*\[Chi]S + 
           (39 + 8*\[Nu]*(-40 + 23*\[Nu]))*\[Chi]S^2))/36 + 
        (Sqrt[5/7]*E^((4*I)*\[Xi])*(3*\[Delta]*\[Kappa]A*(715 - 1849*\[Nu]) + 
           2145*(\[Kappa]S + \[Chi]A^2) + 9*\[Nu]*
            (\[Kappa]S*(-1093 + 1430*\[Nu]) + (-1767 + 2860*\[Nu])*
              \[Chi]A^2) + 10*\[Delta]*(429 - 2941*\[Nu])*\[Chi]A*\[Chi]S + 
           (2145 + \[Nu]*(-22087 + 8624*\[Nu]))*\[Chi]S^2))/216 + 
        (Sqrt[5/7]*(3*\[Delta]*\[Kappa]A*(827 - 1415*\[Nu]) + 
           2481*(\[Kappa]S + \[Chi]A^2) + 9*\[Nu]*
            (\[Kappa]S*(-1023 + 1654*\[Nu]) + (-2285 + 3308*\[Nu])*
              \[Chi]A^2) + 2*\[Delta]*(2481 - 15295*\[Nu])*\[Chi]A*\[Chi]S + 
           (2481 + \[Nu]*(-19949 + 9280*\[Nu]))*\[Chi]S^2))/
         (216*E^((4*I)*\[Xi])) + (Sqrt[5/7]*E^((2*I)*\[Xi])*
          (3*\[Delta]*\[Kappa]A*(761 - 1565*\[Nu]) + 
           2283*(\[Kappa]S + \[Chi]A^2) + 9*\[Nu]*
            (\[Kappa]S*(-1029 + 1522*\[Nu]) + (-2015 + 3044*\[Nu])*
              \[Chi]A^2) + 2*\[Delta]*(2283 - 14131*\[Nu])*\[Chi]A*\[Chi]S + 
           (2283 + \[Nu]*(-19259 + 11104*\[Nu]))*\[Chi]S^2))/432 + 
        (Sqrt[5/7]*(3*\[Delta]*\[Kappa]A*(799 - 1579*\[Nu]) + 
           2397*(\[Kappa]S + \[Chi]A^2) + 9*\[Nu]*
            (\[Kappa]S*(-1059 + 1598*\[Nu]) + (-2137 + 3196*\[Nu])*
              \[Chi]A^2) + 2*\[Delta]*(2397 - 15077*\[Nu])*\[Chi]A*\[Chi]S + 
           (2397 + \[Nu]*(-20509 + 11552*\[Nu]))*\[Chi]S^2))/
         (432*E^((2*I)*\[Xi]))) + 
      e^6*((-5*Sqrt[35]*(3*\[Delta]*\[Kappa]A*(-1 + 2*\[Nu]) - 
           3*(\[Kappa]S + \[Chi]A^2) - 6*\[Nu]*(-2 + 3*\[Nu])*
            (\[Kappa]S + 2*\[Chi]A^2) + \[Delta]*(-6 + 37*\[Nu])*\[Chi]A*
            \[Chi]S + (-3 + (25 - 14*\[Nu])*\[Nu])*\[Chi]S^2))/12 + 
        (E^((6*I)*\[Xi])*(3*\[Delta]*\[Kappa]A*(24655 - 66376*\[Nu]) + 
           73965*(\[Kappa]S + \[Chi]A^2) + 18*\[Nu]*
            (\[Kappa]S*(-19281 + 24655*\[Nu]) + (-30029 + 49310*\[Nu])*
              \[Chi]A^2) + 2*\[Delta]*(73965 - 518668*\[Nu])*\[Chi]A*
            \[Chi]S + (73965 - 792674*\[Nu] + 270064*\[Nu]^2)*\[Chi]S^2))/
         (576*Sqrt[35]) + (33*(\[Delta]*\[Kappa]A*(2639 - 4348*\[Nu]) + 
            2639*(\[Kappa]S + \[Chi]A^2) + 2*\[Nu]*(\[Kappa]S*(-4813 + 
                7917*\[Nu]) + (-11021 + 15834*\[Nu])*\[Chi]A^2)) + 
          26*\[Delta]*(6699 - 41068*\[Nu])*\[Chi]A*\[Chi]S + 
          (87087 - 688730*\[Nu] + 296896*\[Nu]^2)*\[Chi]S^2)/
         (576*Sqrt[35]*E^((6*I)*\[Xi])) + 
        (3*(\[Delta]*\[Kappa]A*(2285 - 5054*\[Nu]) + 
            2285*(\[Kappa]S + \[Chi]A^2) + 6*\[Nu]*(\[Kappa]S*(-1604 + 
                2285*\[Nu]) + 2*(-1483 + 2285*\[Nu])*\[Chi]A^2)) + 
          \[Delta]*(13710 - 90709*\[Nu])*\[Chi]A*\[Chi]S + 
          (6855 + \[Nu]*(-64741 + 40766*\[Nu]))*\[Chi]S^2)/
         (216*Sqrt[35]*E^((4*I)*\[Xi])) + 
        (E^((4*I)*\[Xi])*(3*(\[Delta]*\[Kappa]A*(2341 - 3942*\[Nu]) + 
             2*\[Kappa]S*\[Nu]*(-4312 + 7023*\[Nu]) + 
             4*\[Nu]*(-4867 + 7023*\[Nu])*\[Chi]A^2 + 
             2341*(\[Kappa]S + \[Chi]A^2)) + \[Delta]*(14046 - 82757*\[Nu])*
            \[Chi]A*\[Chi]S + (7023 + \[Nu]*(-52445 + 40894*\[Nu]))*
            \[Chi]S^2))/(216*Sqrt[35]) + (Sqrt[5/7]*E^((2*I)*\[Xi])*
          (3*\[Delta]*\[Kappa]A*(5020 - 10339*\[Nu]) + 
           15060*(\[Kappa]S + \[Chi]A^2) + 9*\[Nu]*
            (\[Kappa]S*(-6793 + 10040*\[Nu]) + (-13287 + 20080*\[Nu])*
              \[Chi]A^2) + 2*\[Delta]*(15060 - 93733*\[Nu])*\[Chi]A*\[Chi]S + 
           (15060 + \[Nu]*(-128123 + 70336*\[Nu]))*\[Chi]S^2))/1728 + 
        (Sqrt[5/7]*(3*(\[Delta]*\[Kappa]A*(5240 - 10329*\[Nu]) + 
             \[Kappa]S*\[Nu]*(-20809 + 31440*\[Nu]) + 
             \[Nu]*(-42071 + 62880*\[Nu])*\[Chi]A^2 + 
             5240*(\[Kappa]S + \[Chi]A^2)) + 2*\[Delta]*(15720 - 98723*\[Nu])*
            \[Chi]A*\[Chi]S + (15720 + \[Nu]*(-134113 + 73328*\[Nu]))*
            \[Chi]S^2))/(1728*E^((2*I)*\[Xi]))) + 
      e^5*((Sqrt[5/7]*E^((3*I)*\[Xi])*
          (3*(\[Delta]*\[Kappa]A*(19522 - 38367*\[Nu]) + \[Kappa]S*\[Nu]*
              (-77411 + 117132*\[Nu]) + \[Nu]*(-156853 + 234264*\[Nu])*
              \[Chi]A^2 + 19522*(\[Kappa]S + \[Chi]A^2)) + 
           2*\[Delta]*(58566 - 359815*\[Nu])*\[Chi]A*\[Chi]S + 
           (58566 + \[Nu]*(-483335 + 299512*\[Nu]))*\[Chi]S^2))/9216 + 
        (Sqrt[5/7]*(3*(\[Delta]*\[Kappa]A*(20170 - 40931*\[Nu]) + 
             \[Kappa]S*\[Nu]*(-81271 + 121020*\[Nu]) + 
             \[Nu]*(-160769 + 242040*\[Nu])*\[Chi]A^2 + 
             20170*(\[Kappa]S + \[Chi]A^2)) + 2*\[Delta]*
            (60510 - 387043*\[Nu])*\[Chi]A*\[Chi]S + 
           (60510 + \[Nu]*(-533819 + 307672*\[Nu]))*\[Chi]S^2))/
         (9216*E^((3*I)*\[Xi])) + (Sqrt[5/7]*E^(I*\[Xi])*
          (3*\[Delta]*\[Kappa]A*(32174 - 65453*\[Nu]) + 
           96522*(\[Kappa]S + \[Chi]A^2) + 9*\[Nu]*
            (\[Kappa]S*(-43267 + 64348*\[Nu]) + (-85429 + 128696*\[Nu])*
              \[Chi]A^2) + 2*\[Delta]*(96522 - 595325*\[Nu])*\[Chi]A*
            \[Chi]S + (96522 + \[Nu]*(-807877 + 452456*\[Nu]))*\[Chi]S^2))/
         13824 + (Sqrt[5/7]*(3*(\[Delta]*\[Kappa]A*(33190 - 65921*\[Nu]) + 
             \[Kappa]S*\[Nu]*(-132301 + 199140*\[Nu]) + 
             \[Nu]*(-265979 + 398280*\[Nu])*\[Chi]A^2 + 
             33190*(\[Kappa]S + \[Chi]A^2)) + 2*\[Delta]*
            (99570 - 620617*\[Nu])*\[Chi]A*\[Chi]S + 
           (99570 + \[Nu]*(-841577 + 468040*\[Nu]))*\[Chi]S^2))/
         (13824*E^(I*\[Xi])) + (Sqrt[5/7]*E^((5*I)*\[Xi])*
          (3*(\[Delta]*\[Kappa]A*(149062 - 394109*\[Nu]) + \[Kappa]S*\[Nu]*
              (-692233 + 894372*\[Nu]) + \[Nu]*(-1096511 + 1788744*\[Nu])*
              \[Chi]A^2 + 149062*(\[Kappa]S + \[Chi]A^2)) + 
           2*\[Delta]*(447186 - 3105757*\[Nu])*\[Chi]A*\[Chi]S + 
           (447186 + \[Nu]*(-4710725 + 1709992*\[Nu]))*\[Chi]S^2))/27648 + 
        (Sqrt[5/7]*(21*(\[Delta]*\[Kappa]A*(24866 - 41687*\[Nu]) + 
             24866*(\[Kappa]S + \[Chi]A^2) + 3*\[Nu]*(\[Kappa]S*(-30473 + 
                 49732*\[Nu]) + (-68991 + 99464*\[Nu])*\[Chi]A^2)) + 
           2*\[Delta]*(522186 - 3211097*\[Nu])*\[Chi]A*\[Chi]S + 
           (522186 + \[Nu]*(-4164505 + 1858952*\[Nu]))*\[Chi]S^2))/
         (27648*E^((5*I)*\[Xi]))) + 
      e^3*((Sqrt[5/7]*(3*(\[Delta]*\[Kappa]A*(791 - 1562*\[Nu]) + 
             791*(\[Kappa]S + \[Chi]A^2) + 6*\[Nu]*(\[Kappa]S*(-524 + 
                 791*\[Nu]) + 2*(-529 + 791*\[Nu])*\[Chi]A^2)) + 
           2*\[Delta]*(2373 - 14696*\[Nu])*\[Chi]A*\[Chi]S + 
           (2373 + 8*\[Nu]*(-2480 + 1417*\[Nu]))*\[Chi]S^2))/
         (576*E^(I*\[Xi])) + (Sqrt[5/7]*E^(I*\[Xi])*
          (3*\[Delta]*\[Kappa]A*(757 - 1560*\[Nu]) + 
           2271*(\[Kappa]S + \[Chi]A^2) + 6*\[Nu]*
            (\[Kappa]S*(-1537 + 2271*\[Nu]) + (-3005 + 4542*\[Nu])*
              \[Chi]A^2) + 2*\[Delta]*(2271 - 13906*\[Nu])*\[Chi]A*\[Chi]S + 
           (2271 + 2*\[Nu]*(-9433 + 5444*\[Nu]))*\[Chi]S^2))/576 + 
        (Sqrt[5/7]*(3*\[Delta]*\[Kappa]A*(1291 - 2264*\[Nu]) + 
           3873*(\[Kappa]S + \[Chi]A^2) + 6*\[Nu]*
            (\[Kappa]S*(-2423 + 3873*\[Nu]) + (-5323 + 7746*\[Nu])*
              \[Chi]A^2) + 2*\[Delta]*(3873 - 23902*\[Nu])*\[Chi]A*\[Chi]S + 
           (3873 + 2*\[Nu]*(-15679 + 7676*\[Nu]))*\[Chi]S^2))/
         (576*E^((3*I)*\[Xi])) + (Sqrt[5/7]*E^((3*I)*\[Xi])*
          (6*\[Kappa]S*\[Nu]*(-2548 + 3387*\[Nu]) + 
           12*\[Nu]*(-2113 + 3387*\[Nu])*\[Chi]A^2 + 
           3387*(\[Kappa]S + \[Chi]A^2) + 
           (3387 + 8*\[Nu]*(-4222 + 1799*\[Nu]))*\[Chi]S^2 + 
           \[Delta]*(\[Kappa]A*(3387 - 8514*\[Nu]) + 2*(3387 - 22792*\[Nu])*
              \[Chi]A*\[Chi]S)))/576)) + 
    SO*(x^(3/2)*\[Epsilon]^3*((4*Sqrt[5/7]*\[Nu]*\[Chi]S)/3 + 
        e*((11*Sqrt[5/7]*\[Nu]*\[Chi]S)/(6*E^(I*\[Xi])) + 
          (13*Sqrt[5/7]*E^(I*\[Xi])*\[Nu]*\[Chi]S)/6) + 
        e^2*((4*Sqrt[5/7]*\[Nu]*\[Chi]S)/3 + (8*Sqrt[5/7]*\[Nu]*\[Chi]S)/
           (3*E^((2*I)*\[Xi])) + (10*Sqrt[5/7]*E^((2*I)*\[Xi])*\[Nu]*\[Chi]S)/
           3) + e^3*((61*Sqrt[5/7]*\[Nu]*\[Chi]S)/(48*E^(I*\[Xi])) + 
          (59*Sqrt[5/7]*E^(I*\[Xi])*\[Nu]*\[Chi]S)/48 + 
          (185*Sqrt[5/7]*\[Nu]*\[Chi]S)/(48*E^((3*I)*\[Xi])) + 
          (239*Sqrt[5/7]*E^((3*I)*\[Xi])*\[Nu]*\[Chi]S)/48) + 
        e^4*((4*Sqrt[5/7]*\[Nu]*\[Chi]S)/3 + (17*Sqrt[5/7]*\[Nu]*\[Chi]S)/
           (18*E^((2*I)*\[Xi])) + (13*Sqrt[5/7]*E^((2*I)*\[Xi])*\[Nu]*
            \[Chi]S)/18 + (199*Sqrt[5/7]*\[Nu]*\[Chi]S)/
           (36*E^((4*I)*\[Xi])) + (263*Sqrt[5/7]*E^((4*I)*\[Xi])*\[Nu]*
            \[Chi]S)/36) + e^5*((221*Sqrt[35]*\[Nu]*\[Chi]S)/
           (1152*E^(I*\[Xi])) + (1549*Sqrt[5/7]*E^(I*\[Xi])*\[Nu]*\[Chi]S)/
           1152 + (17*Sqrt[35]*\[Nu]*\[Chi]S)/(768*E^((3*I)*\[Xi])) - 
          (367*Sqrt[5/7]*E^((3*I)*\[Xi])*\[Nu]*\[Chi]S)/768 + 
          (2593*Sqrt[35]*\[Nu]*\[Chi]S)/(2304*E^((5*I)*\[Xi])) + 
          (24401*Sqrt[5/7]*E^((5*I)*\[Xi])*\[Nu]*\[Chi]S)/2304) + 
        e^6*((4*Sqrt[5/7]*\[Nu]*\[Chi]S)/3 + (199*Sqrt[5/7]*\[Nu]*\[Chi]S)/
           (144*E^((2*I)*\[Xi])) + (29*Sqrt[35]*E^((2*I)*\[Xi])*\[Nu]*
            \[Chi]S)/144 - (64*\[Nu]*\[Chi]S)/(9*Sqrt[35]*E^((4*I)*\[Xi])) - 
          (128*E^((4*I)*\[Xi])*\[Nu]*\[Chi]S)/(9*Sqrt[35]) + 
          (383*Sqrt[7/5]*\[Nu]*\[Chi]S)/(48*E^((6*I)*\[Xi])) + 
          (3653*E^((6*I)*\[Xi])*\[Nu]*\[Chi]S)/(48*Sqrt[35]))) + 
      x^(5/2)*\[Epsilon]^5*((Sqrt[5/7]*(\[Delta]*(-8 + 43*\[Nu])*\[Chi]A - 
           (8 + \[Nu]*(3 + 26*\[Nu]))*\[Chi]S))/9 + 
        e*((Sqrt[5/7]*E^(I*\[Xi])*(\[Delta]*(-394 + 1713*\[Nu])*\[Chi]A + 
             (-394 + (1481 - 1064*\[Nu])*\[Nu])*\[Chi]S))/144 + 
          (Sqrt[5/7]*(\[Delta]*(-182 + 1639*\[Nu])*\[Chi]A - 
             (182 + \[Nu]*(81 + 1208*\[Nu]))*\[Chi]S))/(144*E^(I*\[Xi]))) + 
        e^2*((Sqrt[5/7]*E^((2*I)*\[Xi])*(\[Delta]*(-424 + 1689*\[Nu])*
              \[Chi]A + (-424 + (2141 - 1040*\[Nu])*\[Nu])*\[Chi]S))/72 + 
          (Sqrt[5/7]*(\[Delta]*(-88 + 571*\[Nu])*\[Chi]A + 
             (-88 + (187 - 452*\[Nu])*\[Nu])*\[Chi]S))/36 + 
          (Sqrt[5/7]*(\[Delta]*(-104 + 1541*\[Nu])*\[Chi]A - 
             (104 + \[Nu]*(199 + 1248*\[Nu]))*\[Chi]S))/
           (72*E^((2*I)*\[Xi]))) + 
        e^4*((Sqrt[5/7]*E^((2*I)*\[Xi])*(\[Delta]*(-1432 + 10763*\[Nu])*
              \[Chi]A + (-1432 + (2103 - 9496*\[Nu])*\[Nu])*\[Chi]S))/432 + 
          (Sqrt[5/7]*E^((4*I)*\[Xi])*(4*\[Delta]*(-1049 + 3813*\[Nu])*
              \[Chi]A + (-4196 + (25879 - 9169*\[Nu])*\[Nu])*\[Chi]S))/216 + 
          (Sqrt[5/7]*(\[Delta]*(-1880 + 11211*\[Nu])*\[Chi]A + 
             (-1880 + (5479 - 8968*\[Nu])*\[Nu])*\[Chi]S))/
           (432*E^((2*I)*\[Xi])) + (Sqrt[5/7]*(\[Delta]*(-72 + 485*\[Nu])*
              \[Chi]A + (-72 + (193 - 400*\[Nu])*\[Nu])*\[Chi]S))/18 + 
          (Sqrt[5/7]*(4*\[Delta]*(-31 + 3275*\[Nu])*\[Chi]A - 
             (124 + \[Nu]*(3633 + 11857*\[Nu]))*\[Chi]S))/
           (216*E^((4*I)*\[Xi]))) + 
        e^3*((Sqrt[5/7]*E^((3*I)*\[Xi])*(\[Delta]*(-12770 + 48217*\[Nu])*
              \[Chi]A + (-12770 + (73185 - 29336*\[Nu])*\[Nu])*\[Chi]S))/
           1152 + (Sqrt[5/7]*(\[Delta]*(-3686 + 24827*\[Nu])*\[Chi]A + 
             (-3686 + (9203 - 20008*\[Nu])*\[Nu])*\[Chi]S))/
           (1152*E^(I*\[Xi])) + (Sqrt[5/7]*E^(I*\[Xi])*
            (\[Delta]*(-3802 + 24509*\[Nu])*\[Chi]A + 
             (-3802 + 279*(35 - 72*\[Nu])*\[Nu])*\[Chi]S))/1152 + 
          (Sqrt[5/7]*(\[Delta]*(-1502 + 42511*\[Nu])*\[Chi]A - 
             (1502 + \[Nu]*(8793 + 36680*\[Nu]))*\[Chi]S))/
           (1152*E^((3*I)*\[Xi]))) + 
        e^5*((Sqrt[5/7]*E^((5*I)*\[Xi])*(\[Delta]*(-1799786 + 6353481*\[Nu])*
              \[Chi]A + (-1799786 + (11658049 - 3776392*\[Nu])*\[Nu])*
              \[Chi]S))/55296 + (Sqrt[5/7]*(\[Delta]*(-130142 + 905427*\[Nu])*
              \[Chi]A + (-130142 + (370171 - 750808*\[Nu])*\[Nu])*\[Chi]S))/
           (27648*E^(I*\[Xi])) + (Sqrt[5/7]*E^(I*\[Xi])*
            (\[Delta]*(-136162 + 897733*\[Nu])*\[Chi]A + 
             (-136162 + (405757 - 746376*\[Nu])*\[Nu])*\[Chi]S))/27648 + 
          (Sqrt[5/7]*(\[Delta]*(-117670 + 490111*\[Nu])*\[Chi]A + 
             (-117670 + (408967 - 364280*\[Nu])*\[Nu])*\[Chi]S))/
           (18432*E^((3*I)*\[Xi])) + (Sqrt[5/7]*E^((3*I)*\[Xi])*
            (\[Delta]*(-22810 + 427273*\[Nu])*\[Chi]A - 
             (22810 + \[Nu]*(280799 + 451112*\[Nu]))*\[Chi]S))/18432 + 
          (Sqrt[5/7]*(\[Delta]*(65834 + 5344111*\[Nu])*\[Chi]A + 
             (65834 - 7*\[Nu]*(260303 + 718056*\[Nu]))*\[Chi]S))/
           (55296*E^((5*I)*\[Xi]))) + 
        e^6*((E^((6*I)*\[Xi])*(\[Delta]*(-304544 + 1049983*\[Nu])*\[Chi]A + 
             (-304544 + (2044131 - 617504*\[Nu])*\[Nu])*\[Chi]S))/
           (1152*Sqrt[35]) + (Sqrt[5/7]*(5*\[Delta]*(-3664 + 26681*\[Nu])*
              \[Chi]A + (-18320 + 3*(17523 - 37264*\[Nu])*\[Nu])*\[Chi]S))/
           (3456*E^((2*I)*\[Xi])) + (Sqrt[5/7]*E^((2*I)*\[Xi])*
            (\[Delta]*(-20704 + 132613*\[Nu])*\[Chi]A + 
             (-20704 + 7*(9751 - 15728*\[Nu])*\[Nu])*\[Chi]S))/3456 + 
          (Sqrt[5/7]*(\[Delta]*(-200 + 1369*\[Nu])*\[Chi]A + 
             (-200 + (585 - 1148*\[Nu])*\[Nu])*\[Chi]S))/36 + 
          (\[Delta]*(-87464 + 162623*\[Nu])*\[Chi]A + 
            (-87464 + 523*(677 - 148*\[Nu])*\[Nu])*\[Chi]S)/
           (1728*Sqrt[35]*E^((4*I)*\[Xi])) + (E^((4*I)*\[Xi])*
            (\[Delta]*(44408 + 83999*\[Nu])*\[Chi]A + 
             (44408 - 3*\[Nu]*(198787 + 62836*\[Nu]))*\[Chi]S))/
           (1728*Sqrt[35]) + (\[Delta]*(27088 + 868039*\[Nu])*\[Chi]A + 
            (27088 - \[Nu]*(346101 + 843008*\[Nu]))*\[Chi]S)/
           (1152*Sqrt[35]*E^((6*I)*\[Xi])))) + x^3*\[Epsilon]^6*
       ((4*Sqrt[5/7]*(I*\[Delta]*\[Chi]A + (I + (-5*I + 2*Pi)*\[Nu])*
            \[Chi]S))/3 + e*((Sqrt[5/7]*((49*I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(49*I + \[Nu]*(-245*I + 98*Pi + (648*I)*ArcCoth[5]))))/
           (12*E^(I*\[Xi])) + (Sqrt[5/7]*E^(I*\[Xi])*
            ((29*I)*\[Delta]*\[Chi]A + \[Chi]S*(29*I + \[Nu]*(-145*I + 
                 58*Pi + (12*I)*Log[2]))))/12) + 
        e^2*((Sqrt[5/7]*E^((2*I)*\[Xi])*((23*I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(23*I + \[Nu]*(-115*I + 46*Pi + (12*I)*Log[2]))))/6 + 
          (Sqrt[5/7]*((53*I)*\[Delta]*\[Chi]A + \[Chi]S*(53*I + 106*Pi*
                \[Nu] + I*\[Nu]*(-265 + 836*Log[2] - 324*Log[3]))))/
           (6*E^((2*I)*\[Xi])) + (2*Sqrt[5/7]*((8*I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(8*I + 16*Pi*\[Nu] - I*\[Nu]*(40 + 84*Log[2] - 
                 81*Log[3]))))/3) + 
        e^4*((Sqrt[5/7]*E^((4*I)*\[Xi])*(611*Pi*\[Nu]*\[Chi]S + 
             (8*I)*(37*(\[Delta]*\[Chi]A + \[Chi]S - 5*\[Nu]*\[Chi]S) + 23*
                \[Nu]*\[Chi]S*Log[2])))/36 + (Sqrt[5/7]*E^((2*I)*\[Xi])*
            ((359*I)*\[Delta]*\[Chi]A + \[Chi]S*(359*I + 711*Pi*\[Nu] - I*
                \[Nu]*(1795 + 5068*Log[2] - 4779*Log[3]))))/36 + 
          (Sqrt[5/7]*((65*I)*\[Delta]*\[Chi]A + \[Chi]S*(65*I + 130*Pi*
                \[Nu] + I*\[Nu]*(-325 + 3544*Log[2] - 1863*Log[3]))))/6 + 
          (Sqrt[5/7]*(2132*Pi*\[Nu]*\[Chi]S + I*(1066*\[Delta]*\[Chi]A + 
               \[Chi]S*(1066 + \[Nu]*(-5330 + 17848*Log[2] + 17577*Log[3] - 
                   15625*Log[5])))))/(36*E^((4*I)*\[Xi])) + 
          (Sqrt[5/7]*(818*Pi*\[Nu]*\[Chi]S + I*(409*\[Delta]*\[Chi]A + 
               \[Chi]S*(409 + \[Nu]*(-2045 - 41908*Log[2] + 6318*Log[3] + 
                   15625*Log[5])))))/(36*E^((2*I)*\[Xi]))) + 
        e^3*((Sqrt[5/7]*E^((3*I)*\[Xi])*((1647*I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(1647*I + \[Nu]*(-8235*I + 3322*Pi + (964*I)*Log[2]))))/
           288 + (Sqrt[5/7]*((811*I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(811*I + 1622*Pi*\[Nu] + I*\[Nu]*(-4055 + 32980*Log[2] - 
                 16524*Log[3]))))/(96*E^(I*\[Xi])) + 
          (Sqrt[5/7]*E^(I*\[Xi])*((711*I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(711*I + 1422*Pi*\[Nu] - I*\[Nu]*(3555 + 8860*Log[2] - 
                 8424*Log[3]))))/96 + (Sqrt[5/7]*(9654*Pi*\[Nu]*\[Chi]S + 
             I*(4827*\[Delta]*\[Chi]A + \[Chi]S*(4827 + \[Nu]*(-24135 - 
                   117484*Log[2] + 5832*Log[3] + 62500*Log[5])))))/
           (288*E^((3*I)*\[Xi]))) + 
        e^5*((Sqrt[5/7]*E^(I*\[Xi])*((31069*I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(31069*I + 62306*Pi*\[Nu] + I*\[Nu]*(-155345 + 
                 2131212*Log[2] - 1151496*Log[3]))))/2304 + 
          (Sqrt[5/7]*E^((3*I)*\[Xi])*((60597*I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(60597*I + 116518*Pi*\[Nu] - I*\[Nu]*(302985 + 
                 953684*Log[2] - 894240*Log[3]))))/4608 + 
          (E^((5*I)*\[Xi])*((265385*I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(265385*I + 569934*Pi*\[Nu] + I*\[Nu]*(-1326925 + 
                 189404*Log[2] - 16524*Log[3]))))/(4608*Sqrt[35]) + 
          (Sqrt[5/7]*(117474*Pi*\[Nu]*\[Chi]S + I*(58737*\[Delta]*\[Chi]A + 
               \[Chi]S*(58737 + \[Nu]*(-293685 + 9547964*Log[2] + 4272912*
                    Log[3] - 6687500*Log[5])))))/(4608*E^((3*I)*\[Xi])) + 
          (Sqrt[5/7]*(67298*Pi*\[Nu]*\[Chi]S + I*(33649*\[Delta]*\[Chi]A + 
               \[Chi]S*(33649 + \[Nu]*(-168245 - 4844612*Log[2] + 
                   908820*Log[3] + 1625000*Log[5])))))/(2304*E^(I*\[Xi])) + 
          ((1151885*I)*\[Delta]*\[Chi]A + \[Chi]S*(1151885*I + 
              2303770*Pi*\[Nu] - I*\[Nu]*(5759425 + 26671284*Log[2] + 
                22368960*Log[3] - 3750000*Log[5] - 23059204*Log[7])))/
           (4608*Sqrt[35]*E^((5*I)*\[Xi]))) + 
        e^6*((Sqrt[5/7]*E^((2*I)*\[Xi])*((9442*I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(9442*I + 19079*Pi*\[Nu] + I*\[Nu]*(-47210 + 
                 801656*Log[2] - 442503*Log[3]))))/576 + 
          (E^((4*I)*\[Xi])*((24648*I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(24648*I + 44659*Pi*\[Nu] - (4*I)*\[Nu]*(30810 + 
                 107888*Log[2] - 100683*Log[3]))))/(288*Sqrt[35]) + 
          (E^((6*I)*\[Xi])*((136602*I)*\[Delta]*\[Chi]A + 
             \[Chi]S*(136602*I + 310351*Pi*\[Nu] + (3*I)*\[Nu]*(-227670 + 
                 31432*Log[2] - 4131*Log[3]))))/(1728*Sqrt[35]) + 
          (Sqrt[5/7]*(21564*Pi*\[Nu]*\[Chi]S + I*(10782*\[Delta]*\[Chi]A + 
               \[Chi]S*(10782 + \[Nu]*(-53910 + 2396968*Log[2] + 792342*
                    Log[3] - 1515625*Log[5])))))/(576*E^((2*I)*\[Xi])) + 
          (Sqrt[5/7]*(30247*Pi*\[Nu]*\[Chi]S + I*(15120*\[Delta]*\[Chi]A + 
               \[Chi]S*(15120 + \[Nu]*(-75600 - 2920288*Log[2] + 602883*
                    Log[3] + 921875*Log[5])))))/864 + 
          ((706662*I)*\[Delta]*\[Chi]A + \[Chi]S*(706662*I + 
              1413324*Pi*\[Nu] + I*\[Nu]*(-3533310 + 50716744*Log[2] + 
                3155841*Log[3] + 156250*Log[5] - 17294403*Log[7])))/
           (1728*Sqrt[35]*E^((6*I)*\[Xi])) + ((28816*I)*\[Delta]*\[Chi]A + 
            \[Chi]S*(28816*I + 57632*Pi*\[Nu] - I*\[Nu]*(144080 + 
                10634560*Log[2] + 9671157*Log[3] - 4453125*Log[5] - 
                5764801*Log[7])))/(576*Sqrt[35]*E^((4*I)*\[Xi])))))
 
H\[Psi]InstTail[3, 3] = Sqrt[x]*(((-3*I)/4)*Sqrt[15/14]*\[Delta] + 
      e*((((-5*I)/4)*Sqrt[5/42]*\[Delta])/E^(I*\[Xi]) - 
        ((5*I)/4)*Sqrt[15/14]*E^(I*\[Xi])*\[Delta]) + 
      e^2*((I/4)*Sqrt[15/14]*\[Delta] - (((23*I)/16)*Sqrt[5/42]*\[Delta])/
         E^((2*I)*\[Xi]) - ((95*I)/16)*Sqrt[5/42]*E^((2*I)*\[Xi])*\[Delta]) + 
      e^3*(((I/32)*Sqrt[35/6]*\[Delta])/E^(I*\[Xi]) + ((85*I)/32)*Sqrt[5/42]*
         E^(I*\[Xi])*\[Delta] - (((19*I)/32)*Sqrt[15/14]*\[Delta])/
         E^((3*I)*\[Xi]) - ((97*I)/32)*Sqrt[15/14]*E^((3*I)*\[Xi])*
         \[Delta]) + e^4*((I/16)*Sqrt[15/14]*\[Delta] + 
        (((5*I)/32)*Sqrt[15/14]*\[Delta])/E^((2*I)*\[Xi]) + 
        ((589*I)/96)*Sqrt[5/42]*E^((2*I)*\[Xi])*\[Delta] - 
        (((437*I)/192)*Sqrt[5/42]*\[Delta])/E^((4*I)*\[Xi]) - 
        ((871*I)/64)*Sqrt[5/42]*E^((4*I)*\[Xi])*\[Delta]) + 
      e^5*((((55*I)/768)*Sqrt[5/42]*\[Delta])/E^(I*\[Xi]) - 
        ((155*I)/768)*Sqrt[5/42]*E^(I*\[Xi])*\[Delta] + 
        (((165*I)/512)*Sqrt[15/14]*\[Delta])/E^((3*I)*\[Xi]) + 
        ((2055*I)/512)*Sqrt[15/14]*E^((3*I)*\[Xi])*\[Delta] - 
        (((4541*I)/1536)*Sqrt[5/42]*\[Delta])/E^((5*I)*\[Xi]) - 
        ((30791*I)/1536)*Sqrt[5/42]*E^((5*I)*\[Xi])*\[Delta]) + 
      e^6*((I/32)*Sqrt[15/14]*\[Delta] - (((11*I)/768)*Sqrt[5/42]*\[Delta])/
         E^((2*I)*\[Xi]) - ((173*I)/768)*Sqrt[35/6]*E^((2*I)*\[Xi])*
         \[Delta] + (((341*I)/192)*Sqrt[5/42]*\[Delta])/E^((4*I)*\[Xi]) + 
        ((4181*I)/192)*Sqrt[5/42]*E^((4*I)*\[Xi])*\[Delta] - 
        (((331*I)/256)*Sqrt[15/14]*\[Delta])/E^((6*I)*\[Xi]) - 
        ((2491*I)/256)*Sqrt[15/14]*E^((6*I)*\[Xi])*\[Delta]))*\[Epsilon] + 
    x^(3/2)*\[Epsilon]^3*(((-3*I)/2)*Sqrt[15/14]*\[Delta]*(-2 + \[Nu]) + 
      e*(((-I/72)*Sqrt[5/42]*\[Delta]*(-662 + 295*\[Nu]))/E^(I*\[Xi]) - 
        (I/72)*Sqrt[5/42]*E^(I*\[Xi])*\[Delta]*(-682 + 521*\[Nu])) + 
      e^3*((I/64)*Sqrt[5/42]*E^((3*I)*\[Xi])*\[Delta]*(344 - 1015*\[Nu]) - 
        ((I/64)*Sqrt[15/14]*\[Delta]*(-376 + 123*\[Nu]))/E^((3*I)*\[Xi]) - 
        ((I/576)*Sqrt[5/42]*\[Delta]*(-4360 + 2333*\[Nu]))/E^(I*\[Xi]) - 
        (I/576)*Sqrt[5/42]*E^(I*\[Xi])*\[Delta]*(-9992 + 2683*\[Nu])) + 
      e^2*(((I/288)*Sqrt[5/42]*\[Delta]*(3625 - 1376*\[Nu]))/
         E^((2*I)*\[Xi]) - (I/48)*Sqrt[5/42]*\[Delta]*(-545 + 214*\[Nu]) - 
        (I/288)*Sqrt[5/42]*E^((2*I)*\[Xi])*\[Delta]*(-2489 + 3148*\[Nu])) + 
      e^4*((I/96)*Sqrt[5/42]*\[Delta]*(901 - 392*\[Nu]) - 
        ((I/864)*Sqrt[5/42]*\[Delta]*(-21416 + 6055*\[Nu]))/E^((4*I)*\[Xi]) - 
        ((I/1728)*Sqrt[5/42]*\[Delta]*(-10399 + 7202*\[Nu]))/
         E^((2*I)*\[Xi]) - (I/1728)*Sqrt[5/42]*E^((2*I)*\[Xi])*\[Delta]*
         (-48815 + 7402*\[Nu]) - (I/864)*Sqrt[5/42]*E^((4*I)*\[Xi])*\[Delta]*
         (1676 + 19457*\[Nu])) + 
      e^6*(((I/512)*\[Delta]*(125369 - 26200*\[Nu]))/
         (Sqrt[210]*E^((6*I)*\[Xi])) - ((19*I)/192)*Sqrt[5/42]*\[Delta]*
         (-89 + 40*\[Nu]) + ((I/6912)*E^((4*I)*\[Xi])*\[Delta]*
          (2862001 + 36850*\[Nu]))/Sqrt[210] - (I/13824)*Sqrt[5/42]*
         E^((2*I)*\[Xi])*\[Delta]*(-47035 + 46364*\[Nu]) - 
        ((I/13824)*Sqrt[5/42]*\[Delta]*(-106547 + 50920*\[Nu]))/
         E^((2*I)*\[Xi]) - ((I/512)*E^((6*I)*\[Xi])*\[Delta]*
          (102191 + 112052*\[Nu]))/Sqrt[210] - 
        ((I/6912)*\[Delta]*(131839 + 166438*\[Nu]))/
         (Sqrt[210]*E^((4*I)*\[Xi]))) + 
      e^5*(((-29*I)/13824)*Sqrt[5/42]*E^(I*\[Xi])*\[Delta]*
         (-4022 + 1813*\[Nu]) - (I/3072)*Sqrt[5/42]*E^((3*I)*\[Xi])*\[Delta]*
         (-147298 + 8111*\[Nu]) - ((I/3072)*Sqrt[5/42]*\[Delta]*
          (-8350 + 13553*\[Nu]))/E^((3*I)*\[Xi]) - 
        ((I/13824)*Sqrt[5/42]*\[Delta]*(-106322 + 53215*\[Nu]))/E^(I*\[Xi]) - 
        ((I/27648)*Sqrt[5/42]*\[Delta]*(-963850 + 235019*\[Nu]))/
         E^((5*I)*\[Xi]) - (I/27648)*Sqrt[5/42]*E^((5*I)*\[Xi])*\[Delta]*
         (436330 + 872149*\[Nu]))) + SO^2*x^(5/2)*\[Epsilon]^5*
     (((9*I)/8)*Sqrt[15/14]*(\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + 
        \[Kappa]A*(-1 + 4*\[Nu]) + \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 
        2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S - \[Delta]*\[Chi]S^2) + 
      e*((((59*I)/8)*Sqrt[5/42]*(\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + 
           \[Kappa]A*(-1 + 4*\[Nu]) + \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 
           2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S - \[Delta]*\[Chi]S^2))/
         E^(I*\[Xi]) + ((3*I)/8)*Sqrt[105/2]*E^(I*\[Xi])*
         (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
          \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
           \[Chi]S - \[Delta]*\[Chi]S^2)) + 
      e^2*(((47*I)/16)*Sqrt[15/14]*(\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + 
          \[Kappa]A*(-1 + 4*\[Nu]) + \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 
          2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S - \[Delta]*\[Chi]S^2) + 
        (((191*I)/16)*Sqrt[5/42]*(\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + 
           \[Kappa]A*(-1 + 4*\[Nu]) + \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 
           2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S - \[Delta]*\[Chi]S^2))/
         E^((2*I)*\[Xi]) + ((79*I)/16)*Sqrt[15/14]*E^((2*I)*\[Xi])*
         (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
          \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
           \[Chi]S - \[Delta]*\[Chi]S^2)) + 
      e^3*((((829*I)/64)*Sqrt[5/42]*(\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + 
           \[Kappa]A*(-1 + 4*\[Nu]) + \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 
           2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S - \[Delta]*\[Chi]S^2))/
         E^(I*\[Xi]) + ((107*I)/64)*Sqrt[35/6]*E^(I*\[Xi])*
         (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
          \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
           \[Chi]S - \[Delta]*\[Chi]S^2) + (((397*I)/64)*Sqrt[15/14]*
          (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
           \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
            \[Chi]S - \[Delta]*\[Chi]S^2))/E^((3*I)*\[Xi]) + 
        ((541*I)/64)*Sqrt[15/14]*E^((3*I)*\[Xi])*
         (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
          \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
           \[Chi]S - \[Delta]*\[Chi]S^2)) + 
      e^4*(((11*I)/16)*Sqrt[105/2]*(\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + 
          \[Kappa]A*(-1 + 4*\[Nu]) + \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 
          2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S - \[Delta]*\[Chi]S^2) + 
        (((163*I)/32)*Sqrt[15/14]*(\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + 
           \[Kappa]A*(-1 + 4*\[Nu]) + \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 
           2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S - \[Delta]*\[Chi]S^2))/
         E^((2*I)*\[Xi]) + ((1315*I)/96)*Sqrt[5/42]*E^((2*I)*\[Xi])*
         (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
          \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
           \[Chi]S - \[Delta]*\[Chi]S^2) + (((10883*I)/384)*Sqrt[5/42]*
          (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
           \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
            \[Chi]S - \[Delta]*\[Chi]S^2))/E^((4*I)*\[Xi]) + 
        ((1763*I)/128)*Sqrt[15/14]*E^((4*I)*\[Xi])*
         (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
          \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
           \[Chi]S - \[Delta]*\[Chi]S^2)) + 
      e^5*((((3263*I)/512)*Sqrt[15/14]*(\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + 
           \[Kappa]A*(-1 + 4*\[Nu]) + \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 
           2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S - \[Delta]*\[Chi]S^2))/
         E^(I*\[Xi]) + ((26563*I)/1536)*Sqrt[5/42]*E^(I*\[Xi])*
         (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
          \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
           \[Chi]S - \[Delta]*\[Chi]S^2) + (((805*I)/1024)*Sqrt[105/2]*
          (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
           \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
            \[Chi]S - \[Delta]*\[Chi]S^2))/E^((3*I)*\[Xi]) + 
        ((625*I)/1024)*Sqrt[105/2]*E^((3*I)*\[Xi])*
         (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
          \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
           \[Chi]S - \[Delta]*\[Chi]S^2) + (((130483*I)/3072)*Sqrt[5/42]*
          (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
           \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
            \[Chi]S - \[Delta]*\[Chi]S^2))/E^((5*I)*\[Xi]) + 
        ((66821*I)/1024)*Sqrt[5/42]*E^((5*I)*\[Xi])*
         (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
          \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
           \[Chi]S - \[Delta]*\[Chi]S^2)) + 
      e^6*(((215*I)/32)*Sqrt[15/14]*(\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + 
          \[Kappa]A*(-1 + 4*\[Nu]) + \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 
          2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S - \[Delta]*\[Chi]S^2) + 
        (((4195*I)/192)*Sqrt[5/42]*(\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + 
           \[Kappa]A*(-1 + 4*\[Nu]) + \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 
           2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S - \[Delta]*\[Chi]S^2))/
         E^((2*I)*\[Xi]) + ((7897*I)/384)*Sqrt[5/42]*E^((2*I)*\[Xi])*
         (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
          \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
           \[Chi]S - \[Delta]*\[Chi]S^2) + 
        (((58319*I)/768)*(\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + 
           \[Kappa]A*(-1 + 4*\[Nu]) + \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 
           2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S - \[Delta]*\[Chi]S^2))/
         (Sqrt[210]*E^((4*I)*\[Xi])) + (((23287*I)/768)*E^((4*I)*\[Xi])*
          (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
           \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
            \[Chi]S - \[Delta]*\[Chi]S^2))/Sqrt[210] + 
        (((13421*I)/128)*Sqrt[3/70]*(\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + 
           \[Kappa]A*(-1 + 4*\[Nu]) + \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 
           2*(-1 + 4*\[Nu])*\[Chi]A*\[Chi]S - \[Delta]*\[Chi]S^2))/
         E^((6*I)*\[Xi]) + ((2689*I)/16)*Sqrt[3/70]*E^((6*I)*\[Xi])*
         (\[Delta]*\[Kappa]S*(-1 + 2*\[Nu]) + \[Kappa]A*(-1 + 4*\[Nu]) + 
          \[Delta]*(-1 + 4*\[Nu])*\[Chi]A^2 + 2*(-1 + 4*\[Nu])*\[Chi]A*
           \[Chi]S - \[Delta]*\[Chi]S^2))) + 
    SO*(x^2*\[Epsilon]^4*(((-3*I)/8)*Sqrt[15/14]*((-4 + 19*\[Nu])*\[Chi]A + 
          \[Delta]*(-4 + 5*\[Nu])*\[Chi]S) + 
        e*(((-I/8)*Sqrt[5/42]*((-70 + 327*\[Nu])*\[Chi]A + 
             5*\[Delta]*(-14 + 9*\[Nu])*\[Chi]S))/E^(I*\[Xi]) - 
          (I/8)*Sqrt[15/14]*E^(I*\[Xi])*((-30 + 139*\[Nu])*\[Chi]A + 
            3*\[Delta]*(-10 + 11*\[Nu])*\[Chi]S)) + 
        e^2*((-I/8)*Sqrt[15/14]*((-24 + 113*\[Nu])*\[Chi]A + 
            \[Delta]*(-24 + 23*\[Nu])*\[Chi]S) - 
          ((I/16)*Sqrt[5/42]*((-198 + 947*\[Nu])*\[Chi]A + 
             \[Delta]*(-198 + 97*\[Nu])*\[Chi]S))/E^((2*I)*\[Xi]) - 
          (I/16)*Sqrt[5/42]*E^((2*I)*\[Xi])*((-346 + 1589*\[Nu])*\[Chi]A + 
            \[Delta]*(-346 + 359*\[Nu])*\[Chi]S)) + 
        e^3*(((-I/64)*Sqrt[15/14]*((-374 + 1823*\[Nu])*\[Chi]A + 
             \[Delta]*(-374 + 141*\[Nu])*\[Chi]S))/E^((3*I)*\[Xi]) - 
          (I/64)*Sqrt[5/42]*E^(I*\[Xi])*((-622 + 2939*\[Nu])*\[Chi]A + 
            \[Delta]*(-622 + 609*\[Nu])*\[Chi]S) - 
          ((I/64)*Sqrt[5/42]*((-802 + 3669*\[Nu])*\[Chi]A + 
             \[Delta]*(-802 + 615*\[Nu])*\[Chi]S))/E^(I*\[Xi]) - 
          (I/64)*Sqrt[15/14]*E^((3*I)*\[Xi])*((-794 + 3633*\[Nu])*\[Chi]A + 
            \[Delta]*(-794 + 795*\[Nu])*\[Chi]S)) + 
        e^4*(((-I/16)*Sqrt[5/42]*(3*(-71 + 322*\[Nu])*\[Chi]A + 
             \[Delta]*(-213 + 172*\[Nu])*\[Chi]S))/E^((2*I)*\[Xi]) - 
          (I/64)*Sqrt[15/14]*((-260 + 1217*\[Nu])*\[Chi]A + 
            \[Delta]*(-260 + 231*\[Nu])*\[Chi]S) - (I/48)*Sqrt[5/42]*
           E^((2*I)*\[Xi])*((-361 + 1731*\[Nu])*\[Chi]A + 
            \[Delta]*(-361 + 375*\[Nu])*\[Chi]S) - 
          ((I/384)*Sqrt[5/42]*((-9508 + 47103*\[Nu])*\[Chi]A + 
             \[Delta]*(-9508 + 2697*\[Nu])*\[Chi]S))/E^((4*I)*\[Xi]) - 
          (I/128)*Sqrt[5/42]*E^((4*I)*\[Xi])*((-7772 + 35485*\[Nu])*\[Chi]A + 
            \[Delta]*(-7772 + 7579*\[Nu])*\[Chi]S)) + 
        e^6*((-I/32)*Sqrt[15/14]*((-158 + 737*\[Nu])*\[Chi]A + 
            \[Delta]*(-158 + 135*\[Nu])*\[Chi]S) - 
          ((I/128)*Sqrt[3/70]*((-10435 + 53128*\[Nu])*\[Chi]A + 
             5*\[Delta]*(-2087 + 270*\[Nu])*\[Chi]S))/E^((6*I)*\[Xi]) + 
          ((I/192)*E^((4*I)*\[Xi])*((-22690 + 100301*\[Nu])*\[Chi]A + 
             5*\[Delta]*(-4538 + 3591*\[Nu])*\[Chi]S))/Sqrt[210] - 
          ((I/384)*Sqrt[5/42]*((-6349 + 29210*\[Nu])*\[Chi]A + 
             \[Delta]*(-6349 + 4608*\[Nu])*\[Chi]S))/E^((2*I)*\[Xi]) - 
          (I/384)*Sqrt[5/42]*E^((2*I)*\[Xi])*(7*(-841 + 3921*\[Nu])*\[Chi]A + 
            \[Delta]*(-5887 + 5247*\[Nu])*\[Chi]S) - 
          ((I/192)*((-10550 + 41679*\[Nu])*\[Chi]A + \[Delta]*
              (-10550 + 15633*\[Nu])*\[Chi]S))/(Sqrt[210]*E^((4*I)*\[Xi])) - 
          (I/128)*Sqrt[3/70]*E^((6*I)*\[Xi])*((-31585 + 143843*\[Nu])*
             \[Chi]A + \[Delta]*(-31585 + 29631*\[Nu])*\[Chi]S)) + 
        e^5*((I/1024)*Sqrt[15/14]*E^((3*I)*\[Xi])*
           ((-586 + 1889*\[Nu])*\[Chi]A - \[Delta]*(586 + 5*\[Nu])*\[Chi]S) - 
          ((I/3072)*Sqrt[5/42]*((-107070 + 538099*\[Nu])*\[Chi]A + 
             5*\[Delta]*(-21414 + 4333*\[Nu])*\[Chi]S))/E^((5*I)*\[Xi]) - 
          ((I/1024)*Sqrt[15/14]*(5*(-898 + 3941*\[Nu])*\[Chi]A + 
             \[Delta]*(-4490 + 4379*\[Nu])*\[Chi]S))/E^((3*I)*\[Xi]) - 
          ((I/1536)*Sqrt[5/42]*(5*(-4822 + 22079*\[Nu])*\[Chi]A + 
             \[Delta]*(-24110 + 18073*\[Nu])*\[Chi]S))/E^(I*\[Xi]) - 
          (I/1536)*Sqrt[5/42]*E^(I*\[Xi])*((-20178 + 94717*\[Nu])*\[Chi]A + 
            \[Delta]*(-20178 + 18119*\[Nu])*\[Chi]S) - (I/3072)*Sqrt[5/42]*
           E^((5*I)*\[Xi])*((-294530 + 1342789*\[Nu])*\[Chi]A + 
            5*\[Delta]*(-58906 + 56243*\[Nu])*\[Chi]S))) + 
      x^3*\[Epsilon]^6*((-I/8)*Sqrt[3/70]*((10 + \[Nu]*(-279 + 407*\[Nu]))*
           \[Chi]A + \[Delta]*(10 + \[Nu]*(-1 + 241*\[Nu]))*\[Chi]S) + 
        e*((-I/432)*Sqrt[5/42]*E^(I*\[Xi])*
           ((-4648 + \[Nu]*(6234 + 42751*\[Nu]))*\[Chi]A + 
            \[Delta]*(-4648 + \[Nu]*(17806 + 18749*\[Nu]))*\[Chi]S) - 
          ((I/432)*((-55480 + \[Nu]*(135498 + 172981*\[Nu]))*\[Chi]A + 
             \[Delta]*(-55480 + \[Nu]*(101662 + 82943*\[Nu]))*\[Chi]S))/
           (Sqrt[210]*E^(I*\[Xi]))) + 
        e^3*(((-I/3456)*((-1654340 + \[Nu]*(5757252 + 1624409*\[Nu]))*
              \[Chi]A + \[Delta]*(-1654340 + 13*\[Nu]*(231596 + 
                 107479*\[Nu]))*\[Chi]S))/(Sqrt[210]*E^(I*\[Xi])) - 
          (I/3456)*Sqrt[5/42]*E^(I*\[Xi])*
           ((-206060 + \[Nu]*(685200 + 275579*\[Nu]))*\[Chi]A + 
            \[Delta]*(-206060 + \[Nu]*(507104 + 260209*\[Nu]))*\[Chi]S) - 
          ((I/384)*((-76340 + \[Nu]*(4472 + 703649*\[Nu]))*\[Chi]A + 
             \[Delta]*(-76340 + \[Nu]*(236248 + 291907*\[Nu]))*\[Chi]S))/
           (Sqrt[210]*E^((3*I)*\[Xi])) - ((I/384)*E^((3*I)*\[Xi])*
            ((145620 + \[Nu]*(-759596 + 1244703*\[Nu]))*\[Chi]A + 
             \[Delta]*(145620 + \[Nu]*(143916 + 320669*\[Nu]))*\[Chi]S))/
           Sqrt[210]) + 
        e^2*(((-I/288)*((-28720 + \[Nu]*(55881 + 129112*\[Nu]))*\[Chi]A + 
             \[Delta]*(-28720 + \[Nu]*(88039 + 78536*\[Nu]))*\[Chi]S))/
           Sqrt[210] - ((I/1728)*((-291320 + \[Nu]*(441723 + 1617266*\[Nu]))*
              \[Chi]A + \[Delta]*(-291320 + \[Nu]*(670277 + 686518*\[Nu]))*
              \[Chi]S))/(Sqrt[210]*E^((2*I)*\[Xi])) - 
          ((I/1728)*E^((2*I)*\[Xi])*((46280 + \[Nu]*(-619017 + 2385526*
                  \[Nu]))*\[Chi]A + \[Delta]*(46280 + \[Nu]*(625417 + 
                 766658*\[Nu]))*\[Chi]S))/Sqrt[210]) + 
        e^6*(((-I/20736)*Sqrt[5/42]*((-4682981 + \[Nu]*(16869918 + 
                 3699239*\[Nu]))*\[Chi]A + \[Delta]*(-4682981 + 8643974*
                \[Nu] + 3658576*\[Nu]^2)*\[Chi]S))/E^((2*I)*\[Xi]) - 
          (I/1152)*Sqrt[5/42]*((-178916 + \[Nu]*(624555 + 193256*\[Nu]))*
             \[Chi]A + \[Delta]*(-178916 + \[Nu]*(397697 + 186472*\[Nu]))*
             \[Chi]S) - ((I/768)*E^((6*I)*\[Xi])*
            ((4981051 + \[Nu]*(-21969889 + 18256946*\[Nu]))*\[Chi]A + 
             \[Delta]*(4981051 + \[Nu]*(-2351491 + 3250301*\[Nu]))*\[Chi]S))/
           Sqrt[210] + ((I/103680)*E^((4*I)*\[Xi])*
            ((368971870 + \[Nu]*(-1534510977 + 654848981*\[Nu]))*\[Chi]A + 
             \[Delta]*(368971870 - \[Nu]*(488509363 + 7403957*\[Nu]))*
              \[Chi]S))/Sqrt[210] - 
          ((I/20736)*((-21223762 + 13*\[Nu]*(6388761 + 316981*\[Nu]))*
              \[Chi]A + \[Delta]*(-21223762 + \[Nu]*(30065407 + 7597343*
                  \[Nu]))*\[Chi]S))/(Sqrt[210]*E^((4*I)*\[Xi])) - 
          ((I/3840)*((-724515 + 4*\[Nu]*(-2944324 + 8935197*\[Nu]))*\[Chi]A + 
             \[Delta]*(-724515 + \[Nu]*(9354476 + 15531939*\[Nu]))*\[Chi]S))/
           (Sqrt[210]*E^((6*I)*\[Xi])) - ((I/20736)*E^((2*I)*\[Xi])*
            ((-17556275 + 229*\[Nu]*(281553 + 52991*\[Nu]))*\[Chi]A + 
             \[Delta]*(-17556275 + \[Nu]*(41461103 + 15853372*\[Nu]))*
              \[Chi]S))/Sqrt[210]) + 
        e^4*(((-I/1152)*((-438220 + \[Nu]*(1404459 + 760318*\[Nu]))*\[Chi]A + 
             \[Delta]*(-438220 + \[Nu]*(1053181 + 595154*\[Nu]))*\[Chi]S))/
           Sqrt[210] + ((I/5184)*E^((2*I)*\[Xi])*
            ((3546770 + \[Nu]*(-13502049 + 786322*\[Nu]))*\[Chi]A + 
             \[Delta]*(3546770 - \[Nu]*(7213511 + 2288134*\[Nu]))*\[Chi]S))/
           Sqrt[210] - ((I/5184)*((-3239710 + \[Nu]*(11208291 + 3686422*
                  \[Nu]))*\[Chi]A + \[Delta]*(-3239710 + \[Nu]*(5767909 + 
                 2788406*\[Nu]))*\[Chi]S))/(Sqrt[210]*E^((2*I)*\[Xi])) - 
          ((I/20736)*((-4513480 + \[Nu]*(-10518063 + 68334124*\[Nu]))*
              \[Chi]A + \[Delta]*(-4513480 + \[Nu]*(20158903 + 28544972*
                  \[Nu]))*\[Chi]S))/(Sqrt[210]*E^((4*I)*\[Xi])) - 
          ((I/20736)*E^((4*I)*\[Xi])*((25978240 + \[Nu]*(-121447239 + 
                 140225672*\[Nu]))*\[Chi]A + \[Delta]*(25978240 + \[Nu]*
                (1236239 + 30844336*\[Nu]))*\[Chi]S))/Sqrt[210]) + 
        e^5*(((-I/18432)*Sqrt[5/42]*((-2946496 + \[Nu]*(10572078 + 
                 2676331*\[Nu]))*\[Chi]A + \[Delta]*(-2946496 + \[Nu]*
                (4906762 + 2132657*\[Nu]))*\[Chi]S))/E^((3*I)*\[Xi]) + 
          ((I/18432)*E^((3*I)*\[Xi])*((28904000 + \[Nu]*(-116901462 + 
                 34996231*\[Nu]))*\[Chi]A + \[Delta]*(28904000 - \[Nu]*
                (46793858 + 7244827*\[Nu]))*\[Chi]S))/Sqrt[210] - 
          ((I/82944)*E^(I*\[Xi])*((-50927120 + \[Nu]*(178530438 + 
                 54569531*\[Nu]))*\[Chi]A + \[Delta]*(-50927120 + \[Nu]*
                (118132562 + 53494513*\[Nu]))*\[Chi]S))/Sqrt[210] - 
          ((I/82944)*((-78875440 + \[Nu]*(289901478 + 43096981*\[Nu]))*
              \[Chi]A + \[Delta]*(-78875440 + \[Nu]*(145338802 + 
                 56136863*\[Nu]))*\[Chi]S))/(Sqrt[210]*E^(I*\[Xi])) - 
          ((I/165888)*((-36195680 + \[Nu]*(-241332306 + 934499093*\[Nu]))*
              \[Chi]A + \[Delta]*(-36195680 + \[Nu]*(255202826 + 397369039*
                  \[Nu]))*\[Chi]S))/(Sqrt[210]*E^((5*I)*\[Xi])) - 
          ((I/165888)*E^((5*I)*\[Xi])*((508125920 + \[Nu]*(-2284554786 + 
                 2162436283*\[Nu]))*\[Chi]A + \[Delta]*(508125920 + \[Nu]*
                (-150286694 + 422468609*\[Nu]))*\[Chi]S))/Sqrt[210]))) + 
    x^2*\[Epsilon]^4*
     (e*((\[Delta]*(-329 - (235*I)*Pi + 2900*Log[2] - 810*Log[6]))/
         (4*Sqrt[210]*E^(I*\[Xi])) + (Sqrt[3/70]*E^(I*\[Xi])*\[Delta]*
          (-133 - (95*I)*Pi - 540*Log[2] + 270*Log[6]))/4) + 
      e^2*((Sqrt[3/70]*\[Delta]*(-119 - (85*I)*Pi + 3980*Log[2] - 
           1350*Log[6]))/4 + (Sqrt[5/42]*\[Delta]*(-434 - (310*I)*Pi - 
           7331*Log[2] + 3125*Log[5] + 567*Log[6]))/(16*E^((2*I)*\[Xi])) + 
        (Sqrt[5/42]*E^((2*I)*\[Xi])*\[Delta]*(-574 - (410*I)*Pi - 
           2773*Log[2] + 1377*Log[6]))/16) + 
      e^3*((Sqrt[3/70]*E^((3*I)*\[Xi])*\[Delta]*(-3199 - (2285*I)*Pi - 
           8650*Log[2] + 8460*Log[3]))/32 + (Sqrt[3/70]*\[Delta]*
          (-2289 - (1635*I)*Pi + 40570*Log[2] + 25560*Log[3] - 31250*Log[5]))/
         (32*E^((3*I)*\[Xi])) + (E^(I*\[Xi])*\[Delta]*(-3157 - (2255*I)*Pi + 
           227850*Log[2] - 81000*Log[6]))/(32*Sqrt[210]) + 
        (\[Delta]*(-3227 - (2305*I)*Pi - 347950*Log[2] + 93750*Log[5] + 
           56700*Log[6]))/(32*Sqrt[210]*E^(I*\[Xi]))) + 
      e^5*((E^((5*I)*\[Xi])*\[Delta]*((-819431*I)*Pi - 
           7*(164669 + 499038*Log[2] - 483408*Log[3])))/(1536*Sqrt[210]) + 
        (Sqrt[5/42]*E^((3*I)*\[Xi])*\[Delta]*(-1911 - (1515*I)*Pi + 
           1611094*Log[2] - 940896*Log[3]))/512 + 
        (Sqrt[5/42]*\[Delta]*(-19537 - (13955*I)*Pi + 9203386*Log[2] + 
           1521342*Log[3] - 4950000*Log[5]))/(768*E^(I*\[Xi])) + 
        (Sqrt[5/42]*E^(I*\[Xi])*\[Delta]*(-19607 - (14005*I)*Pi - 
           6983706*Log[2] + 1902042*Log[3] + 1762500*Log[5]))/768 + 
        (\[Delta]*(-768733 - (549095*I)*Pi + 58992818*Log[2] + 
           9796464*Log[3] - 750000*Log[5] - 24706290*Log[7]))/
         (1536*Sqrt[210]*E^((5*I)*\[Xi])) + (Sqrt[5/42]*\[Delta]*
          (-7329 - (5235*I)*Pi - 4313590*Log[2] - 3339792*Log[3] + 
           2175000*Log[5] + 1647086*Log[7]))/(512*E^((3*I)*\[Xi]))) + 
      e^4*((E^((2*I)*\[Xi])*\[Delta]*(164*(-49 - (35*I)*Pi + 5210*Log[2]) - 
           488025*Log[3]))/(96*Sqrt[210]) + (E^((4*I)*\[Xi])*\[Delta]*
          (-123116 - (87815*I)*Pi - 356360*Log[2] + 346545*Log[3]))/
         (256*Sqrt[210]) + (Sqrt[3/70]*\[Delta]*(-1064 - (760*I)*Pi + 
           230160*Log[2] + 65025*Log[3] - 140625*Log[5]))/
         (32*E^((2*I)*\[Xi])) + (Sqrt[3/70]*\[Delta]*(-4956 - (3540*I)*Pi - 
           963880*Log[2] + 236790*Log[3] + 265625*Log[5]))/128 + 
        (\[Delta]*(-253988 - (181420*I)*Pi - 6207320*Log[2] - 
           5600745*Log[3] + 1968750*Log[5] + 4117715*Log[7]))/
         (768*Sqrt[210]*E^((4*I)*\[Xi]))) + 
      e^6*((E^((6*I)*\[Xi])*\[Delta]*((-4158895*I)*Pi - 
           368*(15981 + 50030*Log[2] - 48330*Log[3])))/(5120*Sqrt[210]) + 
        (E^((4*I)*\[Xi])*\[Delta]*(133574 + (88295*I)*Pi + 25260980*Log[2] - 
           14998770*Log[3]))/(960*Sqrt[210]) + 
        (Sqrt[5/42]*\[Delta]*((-1260*I)*Pi + 1637416*Log[2] + 
           252*(-7 + 708*Log[3]) - 821875*Log[5]))/64 + 
        (E^((2*I)*\[Xi])*\[Delta]*(-873376 - (615965*I)*Pi - 
           507244480*Log[2] + 147272985*Log[3] + 120328125*Log[5]))/
         (6144*Sqrt[210]) + (\[Delta]*(7294 + (5210*I)*Pi + 
           146465380*Log[2] + 50751360*Log[3] - 12890625*Log[5] - 
           70001155*Log[7]))/(960*Sqrt[210]*E^((4*I)*\[Xi])) + 
        (\[Delta]*(-7652736 - (5466240*I)*Pi - 1030088960*Log[2] + 
           418021965*Log[3] - 390625*Log[5] + 144120025*Log[7]))/
         (10240*Sqrt[210]*E^((6*I)*\[Xi])) + 
        (\[Delta]*(-854336 - (610240*I)*Pi - 756648320*Log[2] - 
           477855045*Log[3] + 400578125*Log[5] + 210003465*Log[7]))/
         (6144*Sqrt[210]*E^((2*I)*\[Xi]))) + 
      (9*Sqrt[3/70]*\[Delta]*(-7 - (5*I)*Pi + Log[59049/1024]))/4)


(* ::Subsection::Closed:: *)
(*oscillatory memory*)


H\[Psi]OscMem[2, 0] = x^(5/2)*\[Epsilon]^5*
    ((((16*I)/7)*Sqrt[6]*e*(-1 + E^((2*I)*\[Xi]))*\[Nu])/E^(I*\[Xi]) + 
     (((647*I)/42)*e^2*(-1 + E^((4*I)*\[Xi]))*\[Nu])/
      (Sqrt[6]*E^((2*I)*\[Xi])) + (((80*I)/63)*Sqrt[2/3]*e^3*
       (-8 - 21*E^((2*I)*\[Xi]) + 21*E^((4*I)*\[Xi]) + 8*E^((6*I)*\[Xi]))*
       \[Nu])/E^((3*I)*\[Xi]) + ((I/336)*e^4*(-9413 - 14264*E^((2*I)*\[Xi]) + 
        14264*E^((6*I)*\[Xi]) + 9413*E^((8*I)*\[Xi]))*\[Nu])/
      (Sqrt[6]*E^((4*I)*\[Xi])) + 
     ((I/1260)*e^5*(-49471 - 51925*E^((2*I)*\[Xi]) - 161350*E^((4*I)*\[Xi]) + 
        161350*E^((6*I)*\[Xi]) + 51925*E^((8*I)*\[Xi]) + 
        49471*E^((10*I)*\[Xi]))*\[Nu])/(Sqrt[6]*E^((5*I)*\[Xi])) + 
     ((I/5040)*e^6*(-279128 - 207983*E^((2*I)*\[Xi]) - 
        454160*E^((4*I)*\[Xi]) + 454160*E^((8*I)*\[Xi]) + 
        207983*E^((10*I)*\[Xi]) + 279128*E^((12*I)*\[Xi]))*\[Nu])/
      (Sqrt[6]*E^((6*I)*\[Xi])))
 
H\[Psi]OscMem[2, 2] = x^(3/2)*\[Epsilon]^3*
     (((-13*I)/252)*e^2*E^((2*I)*\[Xi])*\[Nu] - ((13*I)/126)*e^3*E^(I*\[Xi])*
       (-1 + E^((2*I)*\[Xi]))*\[Nu] - (I/1008)*e^4*(39 - 70*E^((2*I)*\[Xi]) + 
        169*E^((4*I)*\[Xi]))*\[Nu] - 
      ((I/3024)*e^5*(13 - 555*E^((2*I)*\[Xi]) - 225*E^((4*I)*\[Xi]) + 
         767*E^((6*I)*\[Xi]))*\[Nu])/E^(I*\[Xi]) - 
      ((I/12096)*e^6*(26 + 1320*E^((2*I)*\[Xi]) - 1509*E^((4*I)*\[Xi]) - 
         1352*E^((6*I)*\[Xi]) + 4485*E^((8*I)*\[Xi]))*\[Nu])/
       E^((2*I)*\[Xi])) + SO*x^2*\[Epsilon]^4*
     (((13*I)/378)*e^2*E^((2*I)*\[Xi])*\[Nu]*(-2*\[Delta]*\[Chi]A + 
        (-2 + \[Nu])*\[Chi]S) + ((13*I)/189)*e^3*E^(I*\[Xi])*
       (-1 + E^((2*I)*\[Xi]))*\[Nu]*(-2*\[Delta]*\[Chi]A + 
        (-2 + \[Nu])*\[Chi]S) + (I/1512)*e^4*(39 - 44*E^((2*I)*\[Xi]) + 
        169*E^((4*I)*\[Xi]))*\[Nu]*(-2*\[Delta]*\[Chi]A + 
        (-2 + \[Nu])*\[Chi]S) + ((I/4536)*e^5*(13 - 711*E^((2*I)*\[Xi]) - 
         69*E^((4*I)*\[Xi]) + 767*E^((6*I)*\[Xi]))*\[Nu]*
        (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/E^(I*\[Xi]) + 
      ((I/18144)*e^6*(26 + 1554*E^((2*I)*\[Xi]) - 1695*E^((4*I)*\[Xi]) - 
         338*E^((6*I)*\[Xi]) + 4485*E^((8*I)*\[Xi]))*\[Nu]*
        (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/E^((2*I)*\[Xi])) + 
    \[Epsilon]^5*(x^(5/2)*((((-4*I)/63)*e*(-1 + 3*E^((2*I)*\[Xi]))*\[Nu])/
         E^(I*\[Xi]) - ((I/1512)*e^3*\[Nu]*(-342 - 162*E^((2*I)*\[Xi]) + 
           E^((4*I)*\[Xi])*(7449 - 12032*\[Nu]) + E^((6*I)*\[Xi])*
            (-6465 + 12032*\[Nu])))/E^((3*I)*\[Xi]) - 
        ((I/3024)*e^2*\[Nu]*(-390 - 456*E^((2*I)*\[Xi]) + 
           E^((4*I)*\[Xi])*(-6273 + 12110*\[Nu])))/E^((2*I)*\[Xi]) - 
        ((I/12096)*e^4*\[Nu]*(-4498 - 1476*E^((2*I)*\[Xi]) + 
           3*E^((4*I)*\[Xi])*(-8623 + 11902*\[Nu]) + E^((6*I)*\[Xi])*
            (-9670 + 21724*\[Nu]) + E^((8*I)*\[Xi])*(-84705 + 155558*\[Nu])))/
         E^((4*I)*\[Xi]) - ((I/36288)*e^5*\[Nu]*(-21372 - 
           4212*E^((2*I)*\[Xi]) + E^((6*I)*\[Xi])*(614925 - 1023564*\[Nu]) + 
           93*E^((8*I)*\[Xi])*(-1825 + 3324*\[Nu]) + E^((4*I)*\[Xi])*
            (-13635 + 11876*\[Nu]) + E^((10*I)*\[Xi])*(-385821 + 
             702556*\[Nu])))/E^((5*I)*\[Xi]) - 
        ((I/725760)*e^6*\[Nu]*(-660582 - 43776*E^((2*I)*\[Xi]) + 
           20*E^((4*I)*\[Xi])*(-10497 + 5899*\[Nu]) - 60*E^((8*I)*\[Xi])*
            (-87781 + 146170*\[Nu]) + 60*E^((6*I)*\[Xi])*
            (-105495 + 163828*\[Nu]) + 4*E^((10*I)*\[Xi])*
            (-1466547 + 2634340*\[Nu]) + 3*E^((12*I)*\[Xi])*
            (-3769863 + 6816610*\[Nu])))/E^((6*I)*\[Xi])) + 
      SO^2*x^(5/2)*(((13*I)/4536)*e^2*E^((2*I)*\[Xi])*\[Nu]*
         (9*\[Delta]*\[Kappa]A + \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 
          92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 
          23*\[Chi]S^2 + 32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2) + 
        (I/18144)*e^4*(39 - 18*E^((2*I)*\[Xi]) + 169*E^((4*I)*\[Xi]))*\[Nu]*
         (9*\[Delta]*\[Kappa]A + \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 
          92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 
          23*\[Chi]S^2 + 32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2) + 
        ((I/54432)*e^5*(13 - 867*E^((2*I)*\[Xi]) + 87*E^((4*I)*\[Xi]) + 
           767*E^((6*I)*\[Xi]))*\[Nu]*(9*\[Delta]*\[Kappa]A + 
           \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 92*\[Nu]*\[Chi]A^2 + 
           2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 23*\[Chi]S^2 + 
           32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2))/E^(I*\[Xi]) + 
        ((I/217728)*e^6*(26 + 1788*E^((2*I)*\[Xi]) - 1725*E^((4*I)*\[Xi]) + 
           676*E^((6*I)*\[Xi]) + 4485*E^((8*I)*\[Xi]))*\[Nu]*
          (9*\[Delta]*\[Kappa]A + \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 
           92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 
           23*\[Chi]S^2 + 32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2))/
         E^((2*I)*\[Xi]) - ((13*I)/2268)*e^3*E^(I*\[Xi])*
         (-1 + E^((2*I)*\[Xi]))*\[Nu]*(-9*\[Delta]*\[Kappa]A + 
          9*\[Kappa]S*(-1 + 2*\[Nu]) + 23*\[Chi]A^2 - 92*\[Nu]*\[Chi]A^2 + 
          2*\[Delta]*(23 - 16*\[Nu])*\[Chi]A*\[Chi]S + 23*\[Chi]S^2 - 
          32*\[Nu]*\[Chi]S^2 + 8*\[Nu]^2*\[Chi]S^2))) + 
    \[Epsilon]^6*(x^3*(((-29*I)/126)*e^2*E^((2*I)*\[Xi])*Pi*\[Nu] - 
        ((29*I)/63)*e^3*E^(I*\[Xi])*(-1 + E^((2*I)*\[Xi]))*Pi*\[Nu] - 
        (I/1512)*e^4*(261 + 473*E^((2*I)*\[Xi]) + 1131*E^((4*I)*\[Xi]))*Pi*
         \[Nu] - ((I/1512)*e^5*(29 - 3121*E^((2*I)*\[Xi]) + 
           1381*E^((4*I)*\[Xi]) + 1711*E^((6*I)*\[Xi]))*Pi*\[Nu])/
         E^(I*\[Xi]) - ((I/24192)*e^6*(232 + 23076*E^((2*I)*\[Xi]) - 
           13897*E^((4*I)*\[Xi]) + 36892*E^((6*I)*\[Xi]) + 
           40020*E^((8*I)*\[Xi]))*Pi*\[Nu])/E^((2*I)*\[Xi])) + 
      SO*x^3*((I/3024)*e^3*E^(I*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
         (\[Delta]*(20017 - 32172*\[Nu])*\[Chi]A + 
          (20017 - 58604*\[Nu] + 16320*\[Nu]^2)*\[Chi]S) + 
        (I/6048)*e^2*E^((2*I)*\[Xi])*\[Nu]*(\[Delta]*(20537 - 32380*\[Nu])*
           \[Chi]A + (20537 - 59228*\[Nu] + 16424*\[Nu]^2)*\[Chi]S) + 
        (I/72576)*e^4*\[Nu]*(-(\[Delta]*(39*E^((4*I)*\[Xi])*(-19577 + 31996*
                \[Nu]) + 4*E^((2*I)*\[Xi])*(-71809 + 90264*\[Nu]) + 
             3*(-57451 + 95476*\[Nu]))*\[Chi]A) + 
          (39*E^((4*I)*\[Xi])*(19577 - 58076*\[Nu] + 16232*\[Nu]^2) + 
            4*E^((2*I)*\[Xi])*(71809 - 173864*\[Nu] + 45240*\[Nu]^2) + 
            3*(57451 - 172692*\[Nu] + 48440*\[Nu]^2))*\[Chi]S) + 
        ((I/72576)*e^5*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
          (-(\[Delta]*(18977 - 31756*\[Nu] + 10*E^((2*I)*\[Xi])*(-198933 + 
                307628*\[Nu]) + E^((4*I)*\[Xi])*(-1132123 + 1878596*\[Nu]))*
             \[Chi]A) + (-18977 + 57356*\[Nu] - 16112*\[Nu]^2 + 
             2*E^((2*I)*\[Xi])*(994665 - 2823964*\[Nu] + 778160*\[Nu]^2) + 
             E^((4*I)*\[Xi])*(1132123 - 3398980*\[Nu] + 953104*\[Nu]^2))*
            \[Chi]S))/E^(I*\[Xi]) - ((I/290304)*e^6*
          (\[Delta]*\[Nu]*(-36914 + E^((4*I)*\[Xi])*(2109756 - 3763332*
                \[Nu]) + 63096*\[Nu] + 18*E^((2*I)*\[Xi])*(-200561 + 322628*
                \[Nu]) + 3*E^((8*I)*\[Xi])*(-2166235 + 3645492*\[Nu]) + 
             E^((6*I)*\[Xi])*(-5679950 + 8033592*\[Nu]))*\[Chi]A - 
           \[Nu]*(36914 - 113464*\[Nu] + 32016*\[Nu]^2 - 12*E^((4*I)*\[Xi])*
              (175813 - 559523*\[Nu] + 159204*\[Nu]^2) + 18*E^((2*I)*\[Xi])*
              (200561 - 587332*\[Nu] + 163336*\[Nu]^2) + 3*E^((8*I)*\[Xi])*
              (2166235 - 6576596*\[Nu] + 1849656*\[Nu]^2) + 2*E^((6*I)*\[Xi])*
              (2839975 - 7520924*\[Nu] + 2023608*\[Nu]^2))*\[Chi]S))/
         E^((2*I)*\[Xi])) + SO^3*x^3*(((-13*I)/3402)*e^2*E^((2*I)*\[Xi])*
         \[Nu]*(9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
             \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
            (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
             \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
            3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
            (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)) - 
        ((13*I)/1701)*e^3*E^(I*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
         (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
             \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
            (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
             \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
            3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
            (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)) - 
        (I/13608)*e^4*(39 + 8*E^((2*I)*\[Xi]) + 169*E^((4*I)*\[Xi]))*\[Nu]*
         (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
             \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
            (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
             \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
            3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
            (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)) - 
        ((I/40824)*e^5*(13 - 1023*E^((2*I)*\[Xi]) + 243*E^((4*I)*\[Xi]) + 
           767*E^((6*I)*\[Xi]))*\[Nu]*(9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + 
             \[Delta]*(-2 + \[Nu])*\[Chi]S) + 2*\[Delta]*\[Chi]A*
            (9*\[Kappa]S*(-1 + 2*\[Nu]) + (7 - 28*\[Nu])*\[Chi]A^2 + 
             3*(7 - 13*\[Nu] + 4*\[Nu]^2)*\[Chi]S^2) + 
           \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/
         E^(I*\[Xi]) - ((I/163296)*e^6*(26 + 2022*E^((2*I)*\[Xi]) - 
           1599*E^((4*I)*\[Xi]) + 1690*E^((6*I)*\[Xi]) + 
           4485*E^((8*I)*\[Xi]))*\[Nu]*(9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + 
             \[Delta]*(-2 + \[Nu])*\[Chi]S) + 2*\[Delta]*\[Chi]A*
            (9*\[Kappa]S*(-1 + 2*\[Nu]) + (7 - 28*\[Nu])*\[Chi]A^2 + 
             3*(7 - 13*\[Nu] + 4*\[Nu]^2)*\[Chi]S^2) + 
           \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/
         E^((2*I)*\[Xi])))
 
H\[Psi]OscMem[3, 1] = x^2*\[Epsilon]^4*
     ((-22*Sqrt[14]*e*E^(I*\[Xi])*\[Delta]*\[Nu])/135 - 
      (22*Sqrt[14]*e^2*(-1 + E^((2*I)*\[Xi]))*\[Delta]*\[Nu])/135 - 
      (e^3*(-308 + 8735*E^((2*I)*\[Xi]) + 2772*E^((4*I)*\[Xi]))*\[Delta]*
        \[Nu])/(1080*Sqrt[14]*E^(I*\[Xi])) - 
      (e^4*(-616 - 33597*E^((2*I)*\[Xi]) + 24357*E^((4*I)*\[Xi]) + 
         9856*E^((6*I)*\[Xi]))*\[Delta]*\[Nu])/(3240*Sqrt[14]*
        E^((2*I)*\[Xi])) - (e^5*(-4158 - 32365*E^((2*I)*\[Xi]) + 
         374448*E^((4*I)*\[Xi]) + 202581*E^((6*I)*\[Xi]) + 
         96250*E^((8*I)*\[Xi]))*\[Delta]*\[Nu])/(25920*Sqrt[14]*
        E^((3*I)*\[Xi])) - (e^6*(-9856 - 52915*E^((2*I)*\[Xi]) - 
         1591890*E^((4*I)*\[Xi]) + 804325*E^((6*I)*\[Xi]) + 
         550960*E^((8*I)*\[Xi]) + 299376*E^((10*I)*\[Xi]))*\[Delta]*\[Nu])/
       (64800*Sqrt[14]*E^((4*I)*\[Xi]))) + SO*x^(5/2)*\[Epsilon]^5*
     ((e*E^(I*\[Xi])*\[Nu]*((-5729 + 19712*\[Nu])*\[Chi]A + 
         \[Delta]*(-5729 + 2464*\[Nu])*\[Chi]S))/(1620*Sqrt[14]) + 
      (e^2*(-1 + E^((2*I)*\[Xi]))*\[Nu]*((-5729 + 19712*\[Nu])*\[Chi]A + 
         \[Delta]*(-5729 + 2464*\[Nu])*\[Chi]S))/(1620*Sqrt[14]) + 
      (e^3*\[Nu]*((5729 - 19712*\[Nu] + 9*E^((4*I)*\[Xi])*
            (-5729 + 19712*\[Nu]) + 4*E^((2*I)*\[Xi])*(-46249 + 
             159472*\[Nu]))*\[Chi]A + \[Delta]*(5729 - 2464*\[Nu] + 
           9*E^((4*I)*\[Xi])*(-5729 + 2464*\[Nu]) + 4*E^((2*I)*\[Xi])*
            (-46249 + 19934*\[Nu]))*\[Chi]S))/(12960*Sqrt[14]*E^(I*\[Xi])) + 
      (e^5*\[Nu]*((154683 - 532224*\[Nu] + 8640*E^((6*I)*\[Xi])*
            (-1013 + 3494*\[Nu]) + 625*E^((8*I)*\[Xi])*
            (-5729 + 19712*\[Nu]) - 128*E^((2*I)*\[Xi])*
            (-10462 + 36061*\[Nu]) + 54*E^((4*I)*\[Xi])*(-344511 + 
             1188608*\[Nu]))*\[Chi]A + \[Delta]*(154683 - 66528*\[Nu] + 
           2160*E^((6*I)*\[Xi])*(-4052 + 1747*\[Nu]) + 625*E^((8*I)*\[Xi])*
            (-5729 + 2464*\[Nu]) - 16*E^((2*I)*\[Xi])*(-83696 + 
             36061*\[Nu]) + 54*E^((4*I)*\[Xi])*(-344511 + 148576*\[Nu]))*
          \[Chi]S))/(622080*Sqrt[14]*E^((3*I)*\[Xi])) + 
      (e^4*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
        ((-5729 + 19712*\[Nu] + 16*E^((4*I)*\[Xi])*(-5729 + 19712*\[Nu]) + 
           E^((2*I)*\[Xi])*(-351971 + 1213088*\[Nu]))*\[Chi]A + 
         \[Delta]*(-5729 + 2464*\[Nu] + 16*E^((4*I)*\[Xi])*
            (-5729 + 2464*\[Nu]) + E^((2*I)*\[Xi])*(-351971 + 151636*\[Nu]))*
          \[Chi]S))/(19440*Sqrt[14]*E^((2*I)*\[Xi])) + 
      (e^6*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
        ((-91664 + 315392*\[Nu] + 486*E^((8*I)*\[Xi])*(-5729 + 19712*\[Nu]) + 
           27*E^((2*I)*\[Xi])*(-23707 + 81696*\[Nu]) + 3*E^((4*I)*\[Xi])*
            (-6347263 + 21890464*\[Nu]) + E^((6*I)*\[Xi])*
            (-8809174 + 30364672*\[Nu]))*\[Chi]A + 
         \[Delta]*(-91664 + 39424*\[Nu] + 486*E^((8*I)*\[Xi])*
            (-5729 + 2464*\[Nu]) + 27*E^((2*I)*\[Xi])*(-23707 + 
             10212*\[Nu]) + 3*E^((4*I)*\[Xi])*(-6347263 + 2736308*\[Nu]) + 
           E^((6*I)*\[Xi])*(-8809174 + 3795584*\[Nu]))*\[Chi]S))/
       (388800*Sqrt[14]*E^((4*I)*\[Xi]))) + 
    \[Epsilon]^6*(x^3*((-121*\[Delta]*\[Nu])/(45*Sqrt[14]) - 
        (e^2*\[Delta]*\[Nu]*(65714 + E^((2*I)*\[Xi])*(335463 - 6616*\[Nu]) + 
           E^((4*I)*\[Xi])*(-167097 + 6616*\[Nu])))/(11880*Sqrt[14]*
          E^((2*I)*\[Xi])) - (e*\[Delta]*\[Nu]*(39732 + E^((2*I)*\[Xi])*
            (-138607 + 20168*\[Nu])))/(11880*Sqrt[14]*E^(I*\[Xi])) + 
        (e^3*\[Delta]*\[Nu]*(-848870 + E^((2*I)*\[Xi])*(-1158395 + 
             20168*\[Nu]) + E^((6*I)*\[Xi])*(1870811 + 35320*\[Nu]) + 
           E^((4*I)*\[Xi])*(2343930 + 563612*\[Nu])))/(95040*Sqrt[14]*
          E^((3*I)*\[Xi])) + (e^4*\[Delta]*\[Nu]*(-2010855 + 
           4*E^((2*I)*\[Xi])*(-497155 + 5042*\[Nu]) + 4*E^((8*I)*\[Xi])*
            (1031179 + 51460*\[Nu]) + 6*E^((6*I)*\[Xi])*
            (420121 + 247072*\[Nu]) - 3*E^((4*I)*\[Xi])*(3983461 + 
             569480*\[Nu])))/(142560*Sqrt[14]*E^((4*I)*\[Xi])) + 
        (e^6*\[Delta]*\[Nu]*(-95514177 - 10*E^((4*I)*\[Xi])*
            (7780026 + 151307*\[Nu]) + E^((2*I)*\[Xi])*(-52452347 + 
             322688*\[Nu]) + 3*E^((12*I)*\[Xi])*(61570609 + 4525184*\[Nu]) - 
           30*E^((6*I)*\[Xi])*(16558956 + 5511281*\[Nu]) + 
           10*E^((8*I)*\[Xi])*(8495720 + 9338059*\[Nu]) + E^((10*I)*\[Xi])*
            (14922117 + 59572670*\[Nu])))/(2851200*Sqrt[14]*
          E^((6*I)*\[Xi])) + (e^5*\[Delta]*\[Nu]*(-99854216 - 
           8*E^((4*I)*\[Xi])*(15366264 + 430285*\[Nu]) + E^((2*I)*\[Xi])*
            (-74658597 + 544536*\[Nu]) + 24*E^((8*I)*\[Xi])*
            (2416139 + 2878583*\[Nu]) + E^((10*I)*\[Xi])*(197166351 + 
             12981176*\[Nu]) + 2*E^((6*I)*\[Xi])*(106584163 + 
             61386984*\[Nu])))/(4561920*Sqrt[14]*E^((5*I)*\[Xi]))) + 
      SO^2*x^3*(-(e*E^(I*\[Xi])*\[Nu]*(2772*\[Kappa]A*(-1 + 4*\[Nu]) + 
            (17372 - 73737*\[Nu] + 39424*\[Nu]^2)*\[Chi]A*\[Chi]S + 
            \[Delta]*(2772*\[Kappa]S*(-1 + 2*\[Nu]) + (8686 - 28336*\[Nu])*
               \[Chi]A^2 + (8686 - 10657*\[Nu] + 2464*\[Nu]^2)*\[Chi]S^2)))/
         (2430*Sqrt[14]) - (e^2*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
          (2772*\[Kappa]A*(-1 + 4*\[Nu]) + (17372 - 73737*\[Nu] + 
             39424*\[Nu]^2)*\[Chi]A*\[Chi]S + \[Delta]*
            (2772*\[Kappa]S*(-1 + 2*\[Nu]) + (8686 - 28336*\[Nu])*\[Chi]A^2 + 
             (8686 - 10657*\[Nu] + 2464*\[Nu]^2)*\[Chi]S^2)))/
         (2430*Sqrt[14]) - (e^4*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
          (9*(616 + 41605*E^((2*I)*\[Xi]) + 9856*E^((4*I)*\[Xi]))*\[Kappa]A*
            (-1 + 4*\[Nu]) + 2*(17372 - 73737*\[Nu] + 39424*\[Nu]^2 + 
             16*E^((4*I)*\[Xi])*(17372 - 73737*\[Nu] + 39424*\[Nu]^2) + 
             5*E^((2*I)*\[Xi])*(234187 - 994977*\[Nu] + 532544*\[Nu]^2))*
            \[Chi]A*\[Chi]S + \[Delta]*(9*(616 + 41605*E^((2*I)*\[Xi]) + 9856*
                E^((4*I)*\[Xi]))*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (17372 - 56672*\[Nu] - 64*E^((4*I)*\[Xi])*(-4343 + 
                 14168*\[Nu]) - 5*E^((2*I)*\[Xi])*(-234187 + 765532*\[Nu]))*
              \[Chi]A^2 + (17372 - 21314*\[Nu] + 4928*\[Nu]^2 + 32*
                E^((4*I)*\[Xi])*(8686 - 10657*\[Nu] + 2464*\[Nu]^2) + 5*
                E^((2*I)*\[Xi])*(234187 - 287674*\[Nu] + 66568*\[Nu]^2))*
              \[Chi]S^2)))/(58320*Sqrt[14]*E^((2*I)*\[Xi])) - 
        (e^3*\[Nu]*(9*(-308 + 11199*E^((2*I)*\[Xi]) + 2772*E^((4*I)*\[Xi]))*
            \[Kappa]A*(-1 + 4*\[Nu]) + (-17372 + 73737*\[Nu] - 
             39424*\[Nu]^2 + 9*E^((4*I)*\[Xi])*(17372 - 73737*\[Nu] + 39424*
                \[Nu]^2) + 6*E^((2*I)*\[Xi])*(105011 - 446256*\[Nu] + 238912*
                \[Nu]^2))*\[Chi]A*\[Chi]S + \[Delta]*
            (9*(-308 + 11199*E^((2*I)*\[Xi]) + 2772*E^((4*I)*\[Xi]))*
              \[Kappa]S*(-1 + 2*\[Nu]) + (-8686 + E^((2*I)*\[Xi])*
                (315033 - 1030308*\[Nu]) + E^((4*I)*\[Xi])*(78174 - 
                 255024*\[Nu]) + 28336*\[Nu])*\[Chi]A^2 + 
             (-8686 + 10657*\[Nu] - 2464*\[Nu]^2 + 9*E^((4*I)*\[Xi])*
                (8686 - 10657*\[Nu] + 2464*\[Nu]^2) + E^((2*I)*\[Xi])*
                (315033 - 387096*\[Nu] + 89592*\[Nu]^2))*\[Chi]S^2)))/
         (19440*Sqrt[14]*E^(I*\[Xi])) - 
        (e^5*\[Nu]*(18*(-4158 - 39757*E^((2*I)*\[Xi]) + 
             643224*E^((4*I)*\[Xi]) + 269109*E^((6*I)*\[Xi]) + 
             96250*E^((8*I)*\[Xi]))*\[Kappa]A*(-1 + 4*\[Nu]) + 
           (-27*(17372 - 73737*\[Nu] + 39424*\[Nu]^2) + 625*E^((8*I)*\[Xi])*
              (17372 - 73737*\[Nu] + 39424*\[Nu]^2) + 108*E^((6*I)*\[Xi])*
              (280289 - 1191294*\[Nu] + 637888*\[Nu]^2) - 4*E^((2*I)*\[Xi])*
              (1118819 - 4753674*\[Nu] + 2544448*\[Nu]^2) + 6*E^((4*I)*\[Xi])*
              (12052972 - 51240087*\[Nu] + 27444224*\[Nu]^2))*\[Chi]A*
            \[Chi]S + \[Delta]*(18*(-4158 - 39757*E^((2*I)*\[Xi]) + 643224*
                E^((4*I)*\[Xi]) + 269109*E^((6*I)*\[Xi]) + 96250*
                E^((8*I)*\[Xi]))*\[Kappa]S*(-1 + 2*\[Nu]) - 
             2*(117261 + E^((2*I)*\[Xi])*(1118819 - 3657644*\[Nu]) - 382536*
                \[Nu] + 625*E^((8*I)*\[Xi])*(-4343 + 14168*\[Nu]) + 27*
                E^((6*I)*\[Xi])*(-280289 + 916964*\[Nu]) + 6*E^((4*I)*\[Xi])*
                (-3013243 + 9862768*\[Nu]))*\[Chi]A^2 + 
             (E^((2*I)*\[Xi])*(-2237638 + 2748856*\[Nu] - 636112*\[Nu]^2) - 
               27*(8686 - 10657*\[Nu] + 2464*\[Nu]^2) + 625*E^((8*I)*\[Xi])*
                (8686 - 10657*\[Nu] + 2464*\[Nu]^2) + 54*E^((6*I)*\[Xi])*
                (280289 - 344468*\[Nu] + 79736*\[Nu]^2) + 6*E^((4*I)*\[Xi])*
                (6026486 - 7408607*\[Nu] + 1715264*\[Nu]^2))*\[Chi]S^2)))/
         (933120*Sqrt[14]*E^((3*I)*\[Xi])) - (e^6*(-1 + E^((2*I)*\[Xi]))*
          \[Nu]*(9*(9856 + 75091*E^((2*I)*\[Xi]) + 2486761*E^((4*I)*\[Xi]) + 
             1047456*E^((6*I)*\[Xi]) + 299376*E^((8*I)*\[Xi]))*\[Kappa]A*
            (-1 + 4*\[Nu]) + 2*(16*(17372 - 73737*\[Nu] + 39424*\[Nu]^2) + 
             486*E^((8*I)*\[Xi])*(17372 - 73737*\[Nu] + 39424*\[Nu]^2) + 
             18*E^((6*I)*\[Xi])*(1637564 - 6957819*\[Nu] + 3724288*\[Nu]^2) + 
             E^((2*I)*\[Xi])*(2113697 - 8979687*\[Nu] + 4805824*\[Nu]^2) + 
             E^((4*I)*\[Xi])*(69927587 - 297217377*\[Nu] + 159152704*
                \[Nu]^2))*\[Chi]A*\[Chi]S + \[Delta]*
            (9*(9856 + 75091*E^((2*I)*\[Xi]) + 2486761*E^((4*I)*\[Xi]) + 
               1047456*E^((6*I)*\[Xi]) + 299376*E^((8*I)*\[Xi]))*\[Kappa]S*
              (-1 + 2*\[Nu]) + (277952 + E^((4*I)*\[Xi])*(69927587 - 
                 228782012*\[Nu]) + E^((2*I)*\[Xi])*(2113697 - 6908372*
                  \[Nu]) - 906752*\[Nu] - 1944*E^((8*I)*\[Xi])*(-4343 + 
                 14168*\[Nu]) - 72*E^((6*I)*\[Xi])*(-409391 + 1338416*\[Nu]))*
              \[Chi]A^2 + (32*(8686 - 10657*\[Nu] + 2464*\[Nu]^2) + 972*
                E^((8*I)*\[Xi])*(8686 - 10657*\[Nu] + 2464*\[Nu]^2) + 36*
                E^((6*I)*\[Xi])*(818782 - 1005859*\[Nu] + 232768*\[Nu]^2) + 
               E^((2*I)*\[Xi])*(2113697 - 2596214*\[Nu] + 600728*\[Nu]^2) + 
               E^((4*I)*\[Xi])*(69927587 - 85942394*\[Nu] + 19894088*
                  \[Nu]^2))*\[Chi]S^2)))/(1166400*Sqrt[14]*E^((4*I)*\[Xi]))))
 
H\[Psi]OscMem[3, 3] = x^2*\[Epsilon]^4*
     ((-269*e^3*E^((3*I)*\[Xi])*\[Delta]*\[Nu])/(648*Sqrt[210]) - 
      (269*e^4*E^((2*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Delta]*\[Nu])/
       (216*Sqrt[210]) - (e^5*E^(I*\[Xi])*(5649 - 12754*E^((2*I)*\[Xi]) + 
         13719*E^((4*I)*\[Xi]))*\[Delta]*\[Nu])/(5184*Sqrt[210]) - 
      (e^6*(-538 - 1044*E^((2*I)*\[Xi]) - 11061*E^((4*I)*\[Xi]) + 
         12643*E^((6*I)*\[Xi]))*\[Delta]*\[Nu])/(2592*Sqrt[210])) + 
    SO*x^(5/2)*\[Epsilon]^5*((e^3*E^((3*I)*\[Xi])*\[Nu]*
        (4*(-121 + 538*\[Nu])*\[Chi]A + \[Delta]*(-484 + 269*\[Nu])*\[Chi]S))/
       (972*Sqrt[210]) + (e^4*E^((2*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
        (4*(-121 + 538*\[Nu])*\[Chi]A + \[Delta]*(-484 + 269*\[Nu])*\[Chi]S))/
       (324*Sqrt[210]) + (e^5*E^(I*\[Xi])*\[Nu]*
        ((E^((2*I)*\[Xi])*(42041 - 186848*\[Nu]) + 168*(-121 + 538*\[Nu]) + 
           408*E^((4*I)*\[Xi])*(-121 + 538*\[Nu]))*\[Chi]A + 
         \[Delta]*(E^((2*I)*\[Xi])*(42041 - 23356*\[Nu]) + 
           42*(-484 + 269*\[Nu]) + 102*E^((4*I)*\[Xi])*(-484 + 269*\[Nu]))*
          \[Chi]S))/(15552*Sqrt[210]) + (e^6*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
        ((32*(-121 + 538*\[Nu]) + 752*E^((4*I)*\[Xi])*(-121 + 538*\[Nu]) + 
           E^((2*I)*\[Xi])*(-22949 + 102272*\[Nu]))*\[Chi]A + 
         \[Delta]*(8*(-484 + 269*\[Nu]) + 188*E^((4*I)*\[Xi])*
            (-484 + 269*\[Nu]) + E^((2*I)*\[Xi])*(-22949 + 12784*\[Nu]))*
          \[Chi]S))/(15552*Sqrt[210])) + \[Epsilon]^6*
     (x^3*((11*\[Delta]*\[Nu])/(27*Sqrt[210]) + 
        (e*(33 + E^((2*I)*\[Xi]))*\[Delta]*\[Nu])/(18*Sqrt[210]*
          E^(I*\[Xi])) - (e^2*(-713 - 570*E^((2*I)*\[Xi]) + 
           595*E^((4*I)*\[Xi]))*\[Delta]*\[Nu])/(180*Sqrt[210]*
          E^((2*I)*\[Xi])) - (e^4*\[Delta]*\[Nu]*(-612744 - 
           303083*E^((2*I)*\[Xi]) - 275649*E^((4*I)*\[Xi]) - 
           5*E^((6*I)*\[Xi])*(-419447 + 349080*\[Nu]) + 5*E^((8*I)*\[Xi])*
            (-346099 + 349080*\[Nu])))/(47520*Sqrt[210]*E^((4*I)*\[Xi])) - 
        (e^3*\[Delta]*\[Nu]*(-1056132 - 662013*E^((2*I)*\[Xi]) - 
           1219680*E^((4*I)*\[Xi]) + 5*E^((6*I)*\[Xi])*(-268351 + 
             354998*\[Nu])))/(142560*Sqrt[210]*E^((3*I)*\[Xi])) - 
        (e^5*\[Delta]*\[Nu]*(-24493326 - 8776944*E^((2*I)*\[Xi]) - 
           9541686*E^((4*I)*\[Xi]) + E^((8*I)*\[Xi])*(71473275 - 
             58127180*\[Nu]) + 30*E^((10*I)*\[Xi])*(-3130156 + 
             2922795*\[Nu]) + 6*E^((6*I)*\[Xi])*(-12171484 + 5975745*\[Nu])))/
         (1140480*Sqrt[210]*E^((5*I)*\[Xi])) - 
        (e^6*\[Delta]*\[Nu]*(-39587427 - 8702001*E^((2*I)*\[Xi]) - 
           12173502*E^((4*I)*\[Xi]) + E^((8*I)*\[Xi])*(98457378 - 
             82937820*\[Nu]) + E^((10*I)*\[Xi])*(83203062 - 69887220*\[Nu]) + 
           E^((6*I)*\[Xi])*(1148696 - 6567340*\[Nu]) + 2*E^((12*I)*\[Xi])*
            (-88612207 + 79696190*\[Nu])))/(1140480*Sqrt[210]*
          E^((6*I)*\[Xi]))) + SO^2*x^3*
       (-(e^3*E^((3*I)*\[Xi])*\[Nu]*(2421*\[Kappa]A*(-1 + 4*\[Nu]) + 
            2*(5323 - 27108*\[Nu] + 17216*\[Nu]^2)*\[Chi]A*\[Chi]S + 
            \[Delta]*(2421*\[Kappa]S*(-1 + 2*\[Nu]) + (5323 - 24748*\[Nu])*
               \[Chi]A^2 + (5323 - 8176*\[Nu] + 2152*\[Nu]^2)*\[Chi]S^2)))/
         (11664*Sqrt[210]) - (e^4*E^((2*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*
          \[Nu]*(2421*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(5323 - 27108*\[Nu] + 
             17216*\[Nu]^2)*\[Chi]A*\[Chi]S + \[Delta]*
            (2421*\[Kappa]S*(-1 + 2*\[Nu]) + (5323 - 24748*\[Nu])*\[Chi]A^2 + 
             (5323 - 8176*\[Nu] + 2152*\[Nu]^2)*\[Chi]S^2)))/
         (3888*Sqrt[210]) - (e^6*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
          (9*(538 + 4810*E^((2*I)*\[Xi]) + 12643*E^((4*I)*\[Xi]))*\[Kappa]A*
            (-1 + 4*\[Nu]) + 2*(10646 - 54216*\[Nu] + 34432*\[Nu]^2 + 
             47*E^((4*I)*\[Xi])*(5323 - 27108*\[Nu] + 17216*\[Nu]^2) + 
             5*E^((2*I)*\[Xi])*(18994 - 96849*\[Nu] + 61568*\[Nu]^2))*\[Chi]A*
            \[Chi]S + \[Delta]*(9*(538 + 4810*E^((2*I)*\[Xi]) + 12643*
                E^((4*I)*\[Xi]))*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (10646 + E^((2*I)*\[Xi])*(94970 - 442520*\[Nu]) - 49496*\[Nu] - 
               47*E^((4*I)*\[Xi])*(-5323 + 24748*\[Nu]))*\[Chi]A^2 + 
             (2*(5323 - 8176*\[Nu] + 2152*\[Nu]^2) + 47*E^((4*I)*\[Xi])*
                (5323 - 8176*\[Nu] + 2152*\[Nu]^2) + 10*E^((2*I)*\[Xi])*
                (9497 - 14609*\[Nu] + 3848*\[Nu]^2))*\[Chi]S^2)))/
         (46656*Sqrt[210]) - (e^5*E^(I*\[Xi])*\[Nu]*
          (9*(1883 - 3534*E^((2*I)*\[Xi]) + 4573*E^((4*I)*\[Xi]))*\[Kappa]A*
            (-1 + 4*\[Nu]) + 2*(7*(5323 - 27108*\[Nu] + 17216*\[Nu]^2) + 
             17*E^((4*I)*\[Xi])*(5323 - 27108*\[Nu] + 17216*\[Nu]^2) - 
             6*E^((2*I)*\[Xi])*(11663 - 59373*\[Nu] + 37696*\[Nu]^2))*\[Chi]A*
            \[Chi]S + \[Delta]*(9*(1883 - 3534*E^((2*I)*\[Xi]) + 4573*
                E^((4*I)*\[Xi]))*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (37261 + E^((4*I)*\[Xi])*(90491 - 420716*\[Nu]) - 173236*\[Nu] + 
               6*E^((2*I)*\[Xi])*(-11663 + 54188*\[Nu]))*\[Chi]A^2 + 
             (7*(5323 - 8176*\[Nu] + 2152*\[Nu]^2) + 17*E^((4*I)*\[Xi])*
                (5323 - 8176*\[Nu] + 2152*\[Nu]^2) - 6*E^((2*I)*\[Xi])*
                (11663 - 17906*\[Nu] + 4712*\[Nu]^2))*\[Chi]S^2)))/
         (31104*Sqrt[210])))
 
H\[Psi]OscMem[4, 0] = x^(5/2)*\[Epsilon]^5*
    ((((4*I)/105)*Sqrt[2]*e*(-1 + E^((2*I)*\[Xi]))*\[Nu])/E^(I*\[Xi]) + 
     (((143*I)/1680)*e^2*(-1 + E^((4*I)*\[Xi]))*\[Nu])/
      (Sqrt[2]*E^((2*I)*\[Xi])) + ((I/1890)*e^3*(-211 - 567*E^((2*I)*\[Xi]) + 
        567*E^((4*I)*\[Xi]) + 211*E^((6*I)*\[Xi]))*\[Nu])/
      (Sqrt[2]*E^((3*I)*\[Xi])) + 
     ((I/8064)*e^4*(-1235 - 1928*E^((2*I)*\[Xi]) + 1928*E^((6*I)*\[Xi]) + 
        1235*E^((8*I)*\[Xi]))*\[Nu])/(Sqrt[2]*E^((4*I)*\[Xi])) + 
     ((I/9450)*e^5*(-2019 - 2200*E^((2*I)*\[Xi]) - 6825*E^((4*I)*\[Xi]) + 
        6825*E^((6*I)*\[Xi]) + 2200*E^((8*I)*\[Xi]) + 2019*E^((10*I)*\[Xi]))*
       \[Nu])/(Sqrt[2]*E^((5*I)*\[Xi])) + 
     ((I/67200)*e^6*(-20168 - 15773*E^((2*I)*\[Xi]) - 34160*E^((4*I)*\[Xi]) + 
        34160*E^((8*I)*\[Xi]) + 15773*E^((10*I)*\[Xi]) + 
        20168*E^((12*I)*\[Xi]))*\[Nu])/(Sqrt[2]*E^((6*I)*\[Xi])))
 
H\[Psi]OscMem[4, 2] = x^(3/2)*\[Epsilon]^3*
     ((((-13*I)/3024)*e^2*E^((2*I)*\[Xi])*\[Nu])/Sqrt[5] - 
      (((13*I)/1512)*e^3*E^(I*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu])/Sqrt[5] - 
      ((I/12096)*e^4*(39 - 70*E^((2*I)*\[Xi]) + 169*E^((4*I)*\[Xi]))*\[Nu])/
       Sqrt[5] - ((I/36288)*e^5*(13 - 555*E^((2*I)*\[Xi]) - 
         225*E^((4*I)*\[Xi]) + 767*E^((6*I)*\[Xi]))*\[Nu])/
       (Sqrt[5]*E^(I*\[Xi])) - ((I/145152)*e^6*(26 + 1320*E^((2*I)*\[Xi]) - 
         1509*E^((4*I)*\[Xi]) - 1352*E^((6*I)*\[Xi]) + 4485*E^((8*I)*\[Xi]))*
        \[Nu])/(Sqrt[5]*E^((2*I)*\[Xi]))) + SO*x^2*\[Epsilon]^4*
     ((((13*I)/4536)*e^2*E^((2*I)*\[Xi])*\[Nu]*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu])*\[Chi]S))/Sqrt[5] + (((13*I)/2268)*e^3*E^(I*\[Xi])*
        (-1 + E^((2*I)*\[Xi]))*\[Nu]*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu])*\[Chi]S))/Sqrt[5] + 
      ((I/18144)*e^4*(39 - 44*E^((2*I)*\[Xi]) + 169*E^((4*I)*\[Xi]))*\[Nu]*
        (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/Sqrt[5] + 
      ((I/54432)*e^5*(13 - 711*E^((2*I)*\[Xi]) - 69*E^((4*I)*\[Xi]) + 
         767*E^((6*I)*\[Xi]))*\[Nu]*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu])*\[Chi]S))/(Sqrt[5]*E^(I*\[Xi])) + 
      ((I/217728)*e^6*(26 + 1554*E^((2*I)*\[Xi]) - 1695*E^((4*I)*\[Xi]) - 
         338*E^((6*I)*\[Xi]) + 4485*E^((8*I)*\[Xi]))*\[Nu]*
        (-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/
       (Sqrt[5]*E^((2*I)*\[Xi]))) + \[Epsilon]^5*
     (x^(5/2)*(((-I/189)*e*(-1 + 3*E^((2*I)*\[Xi]))*\[Nu])/
         (Sqrt[5]*E^(I*\[Xi])) - ((I/199584)*e^3*\[Nu]*
          (-3762 - 1782*E^((2*I)*\[Xi]) + E^((4*I)*\[Xi])*
            (275055 - 710146*\[Nu]) + E^((6*I)*\[Xi])*(-264231 + 
             710146*\[Nu])))/(Sqrt[5]*E^((3*I)*\[Xi])) - 
        ((I/399168)*e^2*\[Nu]*(-4290 - 5016*E^((2*I)*\[Xi]) + 
           E^((4*I)*\[Xi])*(-262119 + 711004*\[Nu])))/
         (Sqrt[5]*E^((2*I)*\[Xi])) - ((I/1596672)*e^4*\[Nu]*
          (-49478 - 16236*E^((2*I)*\[Xi]) + 3*E^((4*I)*\[Xi])*
            (-287969 + 708716*\[Nu]) + 2*E^((6*I)*\[Xi])*
            (-243221 + 726760*\[Nu]) + E^((8*I)*\[Xi])*(-3442263 + 
             9222460*\[Nu])))/(Sqrt[5]*E^((4*I)*\[Xi])) - 
        ((I/4790016)*e^5*\[Nu]*(-235092 - 46332*E^((2*I)*\[Xi]) + 
           E^((6*I)*\[Xi])*(23528307 - 61881090*\[Nu]) + 99*E^((8*I)*\[Xi])*
            (-71153 + 195502*\[Nu]) + E^((4*I)*\[Xi])*(-343101 + 
             708430*\[Nu]) + E^((10*I)*\[Xi])*(-15637875 + 41817962*\[Nu])))/
         (Sqrt[5]*E^((5*I)*\[Xi])) - ((I/95800320)*e^6*\[Nu]*
          (-7266402 - 481536*E^((2*I)*\[Xi]) + 20*E^((4*I)*\[Xi])*
            (-212025 + 353786*\[Nu]) - 240*E^((8*I)*\[Xi])*
            (-844467 + 2213690*\[Nu]) + 60*E^((6*I)*\[Xi])*
            (-3859449 + 9935450*\[Nu]) + 4*E^((10*I)*\[Xi])*
            (-60248397 + 164745050*\[Nu]) + 3*E^((12*I)*\[Xi])*
            (-152510193 + 407214260*\[Nu])))/(Sqrt[5]*E^((6*I)*\[Xi]))) + 
      SO^2*x^(5/2)*((((13*I)/54432)*e^2*E^((2*I)*\[Xi])*\[Nu]*
          (9*\[Delta]*\[Kappa]A + \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 
           92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 
           23*\[Chi]S^2 + 32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2))/
         Sqrt[5] + ((I/217728)*e^4*(39 - 18*E^((2*I)*\[Xi]) + 
           169*E^((4*I)*\[Xi]))*\[Nu]*(9*\[Delta]*\[Kappa]A + 
           \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 92*\[Nu]*\[Chi]A^2 + 
           2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 23*\[Chi]S^2 + 
           32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2))/Sqrt[5] + 
        ((I/653184)*e^5*(13 - 867*E^((2*I)*\[Xi]) + 87*E^((4*I)*\[Xi]) + 
           767*E^((6*I)*\[Xi]))*\[Nu]*(9*\[Delta]*\[Kappa]A + 
           \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 92*\[Nu]*\[Chi]A^2 + 
           2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 23*\[Chi]S^2 + 
           32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2))/(Sqrt[5]*E^(I*\[Xi])) + 
        ((I/2612736)*e^6*(26 + 1788*E^((2*I)*\[Xi]) - 1725*E^((4*I)*\[Xi]) + 
           676*E^((6*I)*\[Xi]) + 4485*E^((8*I)*\[Xi]))*\[Nu]*
          (9*\[Delta]*\[Kappa]A + \[Kappa]S*(9 - 18*\[Nu]) - 23*\[Chi]A^2 + 
           92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(-23 + 16*\[Nu])*\[Chi]A*\[Chi]S - 
           23*\[Chi]S^2 + 32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2))/
         (Sqrt[5]*E^((2*I)*\[Xi])) - (((13*I)/27216)*e^3*E^(I*\[Xi])*
          (-1 + E^((2*I)*\[Xi]))*\[Nu]*(-9*\[Delta]*\[Kappa]A + 
           9*\[Kappa]S*(-1 + 2*\[Nu]) + 23*\[Chi]A^2 - 92*\[Nu]*\[Chi]A^2 + 
           2*\[Delta]*(23 - 16*\[Nu])*\[Chi]A*\[Chi]S + 23*\[Chi]S^2 - 
           32*\[Nu]*\[Chi]S^2 + 8*\[Nu]^2*\[Chi]S^2))/Sqrt[5])) + 
    \[Epsilon]^6*(x^3*((((-29*I)/1512)*e^2*E^((2*I)*\[Xi])*Pi*\[Nu])/
         Sqrt[5] - (((29*I)/756)*e^3*E^(I*\[Xi])*(-1 + E^((2*I)*\[Xi]))*Pi*
          \[Nu])/Sqrt[5] - ((I/18144)*e^4*(261 + 473*E^((2*I)*\[Xi]) + 
           1131*E^((4*I)*\[Xi]))*Pi*\[Nu])/Sqrt[5] - 
        ((I/18144)*e^5*(29 - 3121*E^((2*I)*\[Xi]) + 1381*E^((4*I)*\[Xi]) + 
           1711*E^((6*I)*\[Xi]))*Pi*\[Nu])/(Sqrt[5]*E^(I*\[Xi])) - 
        ((I/290304)*e^6*(232 + 23076*E^((2*I)*\[Xi]) - 
           13897*E^((4*I)*\[Xi]) + 36892*E^((6*I)*\[Xi]) + 
           40020*E^((8*I)*\[Xi]))*Pi*\[Nu])/(Sqrt[5]*E^((2*I)*\[Xi]))) + 
      SO*x^3*(((I/199584)*e^2*E^((2*I)*\[Xi])*\[Nu]*
          (\[Delta]*(206184 - 474241*\[Nu])*\[Chi]A + 
           4*(51546 - 177424*\[Nu] + 59441*\[Nu]^2)*\[Chi]S))/Sqrt[5] + 
        ((I/99792)*e^3*E^(I*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
          (\[Delta]*(204754 - 473669*\[Nu])*\[Chi]A + 
           2*(102377 - 353990*\[Nu] + 118739*\[Nu]^2)*\[Chi]S))/Sqrt[5] + 
        ((I/798336)*e^4*\[Nu]*(-(\[Delta]*(-607112 + 1418147*\[Nu] + 
              8*E^((2*I)*\[Xi])*(-106733 + 238883*\[Nu]) + 13*E^((4*I)*\[Xi])*
               (-203544 + 473185*\[Nu]))*\[Chi]A) + 
           4*(151778 - 528840*\[Nu] + 177751*\[Nu]^2 + 13*E^((4*I)*\[Xi])*
              (50886 - 176632*\[Nu] + 59309*\[Nu]^2) + E^((2*I)*\[Xi])*
              (213466 - 716433*\[Nu] + 238982*\[Nu]^2))*\[Chi]S))/Sqrt[5] + 
        ((I/2395008)*e^5*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
          (-(\[Delta]*(201894 - 472525*\[Nu] + 10*E^((2*I)*\[Xi])*(-2009238 + 
                4644485*\[Nu]) + E^((4*I)*\[Xi])*(-11946066 + 27892703*
                 \[Nu]))*\[Chi]A) + 2*(-100947 + 352274*\[Nu] - 
             118453*\[Nu]^2 + 2*E^((2*I)*\[Xi])*(5023095 - 17341448*\[Nu] + 
               5818105*\[Nu]^2) + E^((4*I)*\[Xi])*(5973033 - 20804758*\[Nu] + 
               6992159*\[Nu]^2))*\[Chi]S))/(Sqrt[5]*E^(I*\[Xi])) + 
        ((I/9580032)*e^6*\[Nu]*(-(\[Delta]*(-400928 + E^((4*I)*\[Xi])*(
                23784192 - 56532531*\[Nu]) + 943906*\[Nu] + 6*E^((2*I)*\[Xi])*
               (-6283156 + 14661725*\[Nu]) + 3*E^((8*I)*\[Xi])*(-23173480 + 
                54322643*\[Nu]) + 2*E^((6*I)*\[Xi])*(-27260420 + 
                62273449*\[Nu]))*\[Chi]A) + 
           (8*(50116 - 175708*\[Nu] + 59155*\[Nu]^2) + 12*E^((2*I)*\[Xi])*
              (3141578 - 10929155*\[Nu] + 3673772*\[Nu]^2) + 
             12*E^((8*I)*\[Xi])*(5793370 - 20242456*\[Nu] + 6808831*
                \[Nu]^2) - 3*E^((4*I)*\[Xi])*(7928064 - 28014157*\[Nu] + 
               9448472*\[Nu]^2) + 4*E^((6*I)*\[Xi])*(13630210 - 46557223*
                \[Nu] + 15589276*\[Nu]^2))*\[Chi]S))/
         (Sqrt[5]*E^((2*I)*\[Xi]))) + SO^3*x^3*
       ((((-13*I)/40824)*e^2*E^((2*I)*\[Xi])*\[Nu]*
          (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
              \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
              \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/Sqrt[5] - 
        (((13*I)/20412)*e^3*E^(I*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
          (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
              \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
              \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/Sqrt[5] - 
        ((I/163296)*e^4*(39 + 8*E^((2*I)*\[Xi]) + 169*E^((4*I)*\[Xi]))*\[Nu]*
          (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
              \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
              \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/Sqrt[5] - 
        ((I/489888)*e^5*(13 - 1023*E^((2*I)*\[Xi]) + 243*E^((4*I)*\[Xi]) + 
           767*E^((6*I)*\[Xi]))*\[Nu]*(9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + 
             \[Delta]*(-2 + \[Nu])*\[Chi]S) + 2*\[Delta]*\[Chi]A*
            (9*\[Kappa]S*(-1 + 2*\[Nu]) + (7 - 28*\[Nu])*\[Chi]A^2 + 
             3*(7 - 13*\[Nu] + 4*\[Nu]^2)*\[Chi]S^2) + 
           \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/
         (Sqrt[5]*E^(I*\[Xi])) - ((I/1959552)*e^6*
          (26 + 2022*E^((2*I)*\[Xi]) - 1599*E^((4*I)*\[Xi]) + 
           1690*E^((6*I)*\[Xi]) + 4485*E^((8*I)*\[Xi]))*\[Nu]*
          (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
              \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
              \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
             3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
             (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)))/
         (Sqrt[5]*E^((2*I)*\[Xi]))))
 
H\[Psi]OscMem[4, 4] = x^(3/2)*\[Epsilon]^3*
     (((5*I)/6912)*Sqrt[5/7]*e^4*E^((4*I)*\[Xi])*\[Nu] + 
      ((5*I)/1728)*Sqrt[5/7]*e^5*E^((3*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*
       \[Nu] + ((5*I)/13824)*Sqrt[5/7]*e^6*E^((2*I)*\[Xi])*
       (11 - 27*E^((2*I)*\[Xi]) + 21*E^((4*I)*\[Xi]))*\[Nu]) + 
    SO*x^2*\[Epsilon]^4*(((5*I)/10368)*Sqrt[5/7]*e^4*E^((4*I)*\[Xi])*\[Nu]*
       (2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S) + 
      ((5*I)/2592)*Sqrt[5/7]*e^5*E^((3*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
       (2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S) - 
      ((5*I)/20736)*Sqrt[5/7]*e^6*E^((2*I)*\[Xi])*(11 - 26*E^((2*I)*\[Xi]) + 
        21*E^((4*I)*\[Xi]))*\[Nu]*(-2*\[Delta]*\[Chi]A + 
        (-2 + \[Nu])*\[Chi]S)) + \[Epsilon]^5*
     (x^(5/2)*(((I/9)*\[Nu])/Sqrt[35] + ((I/45)*e*(7 + 15*E^((2*I)*\[Xi]))*
          \[Nu])/(Sqrt[35]*E^(I*\[Xi])) + 
        ((I/4320)*e^2*(1037 + 993*E^((2*I)*\[Xi]) + 3255*E^((4*I)*\[Xi]))*
          \[Nu])/(Sqrt[35]*E^((2*I)*\[Xi])) + 
        ((I/1080)*e^3*(394 + 317*E^((2*I)*\[Xi]) + 239*E^((4*I)*\[Xi]) + 
           1650*E^((6*I)*\[Xi]))*\[Nu])/(Sqrt[35]*E^((3*I)*\[Xi])) - 
        ((I/4561920)*e^4*\[Nu]*(-2492886 - 1581360*E^((2*I)*\[Xi]) - 
           1709400*E^((4*I)*\[Xi]) + 426096*E^((6*I)*\[Xi]) + 
           5*E^((8*I)*\[Xi])*(-2754285 + 288332*\[Nu])))/
         (Sqrt[35]*E^((4*I)*\[Xi])) - ((I/1140480)*e^5*\[Nu]*
          (-922900 - 436524*E^((2*I)*\[Xi]) - 495264*E^((4*I)*\[Xi]) - 
           590216*E^((6*I)*\[Xi]) + E^((8*I)*\[Xi])*(1792359 - 
             1449910*\[Nu]) + 5*E^((10*I)*\[Xi])*(-1303539 + 289982*\[Nu])))/
         (Sqrt[35]*E^((5*I)*\[Xi])) - ((I/45619200)*e^6*\[Nu]*
          (-54176298 - 17196410*E^((2*I)*\[Xi]) - 22597630*E^((4*I)*\[Xi]) - 
           23166000*E^((6*I)*\[Xi]) + 55*E^((8*I)*\[Xi])*(-1249333 + 
             1459660*\[Nu]) + 15*E^((12*I)*\[Xi])*(-31510139 + 
             10201620*\[Nu]) - 4*E^((10*I)*\[Xi])*(-55303447 + 
             45587775*\[Nu])))/(Sqrt[35]*E^((6*I)*\[Xi]))) + 
      SO^2*x^(5/2)*(((5*I)/124416)*Sqrt[5/7]*e^4*E^((4*I)*\[Xi])*\[Nu]*
         (-9*\[Delta]*\[Kappa]A + 9*\[Kappa]S*(-1 + 2*\[Nu]) + 23*\[Chi]A^2 - 
          92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(23 - 16*\[Nu])*\[Chi]A*\[Chi]S + 
          23*\[Chi]S^2 - 32*\[Nu]*\[Chi]S^2 + 8*\[Nu]^2*\[Chi]S^2) + 
        ((5*I)/31104)*Sqrt[5/7]*e^5*E^((3*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*
         \[Nu]*(-9*\[Delta]*\[Kappa]A + 9*\[Kappa]S*(-1 + 2*\[Nu]) + 
          23*\[Chi]A^2 - 92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(23 - 16*\[Nu])*
           \[Chi]A*\[Chi]S + 23*\[Chi]S^2 - 32*\[Nu]*\[Chi]S^2 + 
          8*\[Nu]^2*\[Chi]S^2) + ((5*I)/248832)*Sqrt[5/7]*e^6*E^((2*I)*\[Xi])*
         (11 - 25*E^((2*I)*\[Xi]) + 21*E^((4*I)*\[Xi]))*\[Nu]*
         (-9*\[Delta]*\[Kappa]A + 9*\[Kappa]S*(-1 + 2*\[Nu]) + 23*\[Chi]A^2 - 
          92*\[Nu]*\[Chi]A^2 + 2*\[Delta]*(23 - 16*\[Nu])*\[Chi]A*\[Chi]S + 
          23*\[Chi]S^2 - 32*\[Nu]*\[Chi]S^2 + 8*\[Nu]^2*\[Chi]S^2))) + 
    \[Epsilon]^6*(x^3*((((19*I)/1152)*e^4*E^((4*I)*\[Xi])*Pi*\[Nu])/
         Sqrt[35] + (((19*I)/288)*e^5*E^((3*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*
          Pi*\[Nu])/Sqrt[35] + ((I/414720)*e^6*E^((2*I)*\[Xi])*
          (37620 - 76793*E^((2*I)*\[Xi]) + 71820*E^((4*I)*\[Xi]))*Pi*\[Nu])/
         Sqrt[35]) + SO*x^3*(((I/152064)*e^4*E^((4*I)*\[Xi])*\[Nu]*
          (\[Delta]*(20960 - 63921*\[Nu])*\[Chi]A + 
           4*(5240 - 17008*\[Nu] + 7887*\[Nu]^2)*\[Chi]S))/Sqrt[35] + 
        ((I/114048)*e^5*E^((3*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
          (\[Delta]*(65630 - 192863*\[Nu])*\[Chi]A + 
           2*(32815 - 103698*\[Nu] + 47597*\[Nu]^2)*\[Chi]S))/Sqrt[35] + 
        ((I/912384)*e^6*E^((2*I)*\[Xi])*\[Nu]*
          (-(\[Delta]*(33*(-22960 + 64721*\[Nu]) - 8*E^((2*I)*\[Xi])*(
                -208460 + 581701*\[Nu]) + E^((4*I)*\[Xi])*(-1430480 + 
                4071023*\[Nu]))*\[Chi]A) + 
           (E^((2*I)*\[Xi])*(-1667680 + 5076631*\[Nu] - 2297104*\[Nu]^2) + 
             132*(5740 - 17608*\[Nu] + 7987*\[Nu]^2) + 4*E^((4*I)*\[Xi])*
              (357620 - 1104504*\[Nu] + 502381*\[Nu]^2))*\[Chi]S))/
         Sqrt[35]) + SO^3*x^3*(((5*I)/93312)*Sqrt[5/7]*e^4*E^((4*I)*\[Xi])*
         \[Nu]*(9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
             \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
            (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
             \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
            3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
            (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)) + 
        ((5*I)/23328)*Sqrt[5/7]*e^5*E^((3*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*
         \[Nu]*(9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
             \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
            (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
             \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
            3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
            (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2)) + 
        ((5*I)/186624)*Sqrt[5/7]*e^6*E^((2*I)*\[Xi])*
         (11 - 24*E^((2*I)*\[Xi]) + 21*E^((4*I)*\[Xi]))*\[Nu]*
         (9*\[Kappa]A*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
             \[Chi]S) + 2*\[Delta]*\[Chi]A*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
            (7 - 28*\[Nu])*\[Chi]A^2 + 3*(7 - 13*\[Nu] + 4*\[Nu]^2)*
             \[Chi]S^2) + \[Chi]S*(-9*\[Kappa]S*(2 - 5*\[Nu] + 2*\[Nu]^2) + 
            3*(14 - 69*\[Nu] + 52*\[Nu]^2)*\[Chi]A^2 + 
            (14 - 39*\[Nu] + 24*\[Nu]^2 - 4*\[Nu]^3)*\[Chi]S^2))))
 
H\[Psi]OscMem[5, 1] = x^2*\[Epsilon]^4*((-43*e*E^(I*\[Xi])*\[Delta]*\[Nu])/
       (216*Sqrt[385]) - (43*e^2*(-1 + E^((2*I)*\[Xi]))*\[Delta]*\[Nu])/
       (216*Sqrt[385]) - (e^3*(-43 + 1300*E^((2*I)*\[Xi]) + 
         387*E^((4*I)*\[Xi]))*\[Delta]*\[Nu])/(1728*Sqrt[385]*E^(I*\[Xi])) - 
      (e^4*(-43 - 2466*E^((2*I)*\[Xi]) + 1821*E^((4*I)*\[Xi]) + 
         688*E^((6*I)*\[Xi]))*\[Delta]*\[Nu])/(2592*Sqrt[385]*
        E^((2*I)*\[Xi])) - (e^5*(-1161 - 9520*E^((2*I)*\[Xi]) + 
         112926*E^((4*I)*\[Xi]) + 60912*E^((6*I)*\[Xi]) + 
         26875*E^((8*I)*\[Xi]))*\[Delta]*\[Nu])/(82944*Sqrt[385]*
        E^((3*I)*\[Xi])) - (e^6*(-688 - 3895*E^((2*I)*\[Xi]) - 
         118770*E^((4*I)*\[Xi]) + 60775*E^((6*I)*\[Xi]) + 
         41680*E^((8*I)*\[Xi]) + 20898*E^((10*I)*\[Xi]))*\[Delta]*\[Nu])/
       (51840*Sqrt[385]*E^((4*I)*\[Xi]))) + SO*x^(5/2)*\[Epsilon]^5*
     ((43*e*E^(I*\[Xi])*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*
          \[Chi]S))/(324*Sqrt[385]) + (43*e^2*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
        ((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*\[Chi]S))/
       (324*Sqrt[385]) + (e^3*(-43 + 1472*E^((2*I)*\[Xi]) + 
         387*E^((4*I)*\[Xi]))*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(2592*Sqrt[385]*E^(I*\[Xi])) + 
      (e^4*(-43 - 2724*E^((2*I)*\[Xi]) + 2079*E^((4*I)*\[Xi]) + 
         688*E^((6*I)*\[Xi]))*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(3888*Sqrt[385]*E^((2*I)*\[Xi])) + 
      (e^5*(-1161 - 10552*E^((2*I)*\[Xi]) + 150318*E^((4*I)*\[Xi]) + 
         70200*E^((6*I)*\[Xi]) + 26875*E^((8*I)*\[Xi]))*\[Nu]*
        ((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*\[Chi]S))/
       (124416*Sqrt[385]*E^((3*I)*\[Xi])) + 
      (e^6*(-688 - 4325*E^((2*I)*\[Xi]) - 147300*E^((4*I)*\[Xi]) + 
         82855*E^((6*I)*\[Xi]) + 48560*E^((8*I)*\[Xi]) + 
         20898*E^((10*I)*\[Xi]))*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(77760*Sqrt[385]*E^((4*I)*\[Xi]))) + 
    \[Epsilon]^6*(x^3*((-13*\[Delta]*\[Nu])/(63*Sqrt[385]) - 
        (e^2*\[Delta]*\[Nu]*(15587 + 56*E^((4*I)*\[Xi])*
            (-2355 + 5156*\[Nu]) - 7*E^((2*I)*\[Xi])*(-24963 + 41248*\[Nu])))/
         (39312*Sqrt[385]*E^((2*I)*\[Xi])) - 
        (e*\[Delta]*\[Nu]*(9789 + 7*E^((2*I)*\[Xi])*(-17254 + 41807*\[Nu])))/
         (39312*Sqrt[385]*E^(I*\[Xi])) - (e^4*\[Delta]*\[Nu]*
          (449319 + E^((2*I)*\[Xi])*(611242 - 292649*\[Nu]) + 
           E^((8*I)*\[Xi])*(-2648878 + 4529777*\[Nu]) + 3*E^((6*I)*\[Xi])*
            (-2701778 + 7717255*\[Nu]) - 3*E^((4*I)*\[Xi])*
            (-3994993 + 9129631*\[Nu])))/(471744*Sqrt[385]*E^((4*I)*\[Xi])) - 
        (e^3*\[Delta]*\[Nu]*(194987 + E^((2*I)*\[Xi])*(392582 - 
             292649*\[Nu]) + E^((6*I)*\[Xi])*(-1324070 + 2571233*\[Nu]) + 
           E^((4*I)*\[Xi])*(-5742651 + 16138108*\[Nu])))/
         (314496*Sqrt[385]*E^((3*I)*\[Xi])) - 
        (e^6*\[Delta]*\[Nu]*(40778088 + E^((4*I)*\[Xi])*(70260975 - 
             89235230*\[Nu]) + E^((2*I)*\[Xi])*(31691461 - 9364768*\[Nu]) - 
           840*E^((6*I)*\[Xi])*(-1839552 + 4585619*\[Nu]) + 
           3*E^((12*I)*\[Xi])*(-68188307 + 90318326*\[Nu]) + 
           10*E^((8*I)*\[Xi])*(-89386741 + 254900422*\[Nu]) + 
           E^((10*I)*\[Xi])*(-375664593 + 1130560760*\[Nu])))/
         (18869760*Sqrt[385]*E^((6*I)*\[Xi])) - 
        (e^5*\[Delta]*\[Nu]*(21778913 + E^((4*I)*\[Xi])*(67666026 - 
             108284176*\[Nu]) + E^((2*I)*\[Xi])*(22257906 - 7901523*\[Nu]) + 
           E^((10*I)*\[Xi])*(-116720742 + 175517881*\[Nu]) + 
           3*E^((8*I)*\[Xi])*(-91270681 + 265880048*\[Nu]) + 
           E^((6*I)*\[Xi])*(-798034940 + 2243839626*\[Nu])))/
         (15095808*Sqrt[385]*E^((5*I)*\[Xi]))) + 
      SO^2*x^3*((-43*e*E^(I*\[Xi])*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(3888*Sqrt[385]) - 
        (43*e^2*(-1 + E^((2*I)*\[Xi]))*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(3888*Sqrt[385]) - 
        (e^3*(-43 + 1644*E^((2*I)*\[Xi]) + 387*E^((4*I)*\[Xi]))*\[Nu]*
          (9*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(23 - 108*\[Nu] + 64*\[Nu]^2)*
            \[Chi]A*\[Chi]S + \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (23 - 92*\[Nu])*\[Chi]A^2 + (23 - 32*\[Nu] + 8*\[Nu]^2)*
              \[Chi]S^2)))/(31104*Sqrt[385]*E^(I*\[Xi])) - 
        (e^4*(-43 - 2982*E^((2*I)*\[Xi]) + 2337*E^((4*I)*\[Xi]) + 
           688*E^((6*I)*\[Xi]))*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(46656*Sqrt[385]*
          E^((2*I)*\[Xi])) - (e^5*(-1161 - 11584*E^((2*I)*\[Xi]) + 
           191838*E^((4*I)*\[Xi]) + 79488*E^((6*I)*\[Xi]) + 
           26875*E^((8*I)*\[Xi]))*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(1492992*Sqrt[385]*
          E^((3*I)*\[Xi])) - (e^6*(-688 - 4755*E^((2*I)*\[Xi]) - 
           178410*E^((4*I)*\[Xi]) + 107515*E^((6*I)*\[Xi]) + 
           55440*E^((8*I)*\[Xi]) + 20898*E^((10*I)*\[Xi]))*\[Nu]*
          (9*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(23 - 108*\[Nu] + 64*\[Nu]^2)*
            \[Chi]A*\[Chi]S + \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (23 - 92*\[Nu])*\[Chi]A^2 + (23 - 32*\[Nu] + 8*\[Nu]^2)*
              \[Chi]S^2)))/(933120*Sqrt[385]*E^((4*I)*\[Xi]))))
 
H\[Psi]OscMem[5, 3] = x^2*\[Epsilon]^4*
     ((-31*e^3*E^((3*I)*\[Xi])*\[Delta]*\[Nu])/(1296*Sqrt[330]) - 
      (31*e^4*E^((2*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Delta]*\[Nu])/
       (432*Sqrt[330]) - (e^5*E^(I*\[Xi])*(651 - 1481*E^((2*I)*\[Xi]) + 
         1581*E^((4*I)*\[Xi]))*\[Delta]*\[Nu])/(10368*Sqrt[330]) - 
      (e^6*(-124 - 207*E^((2*I)*\[Xi]) - 2583*E^((4*I)*\[Xi]) + 
         2914*E^((6*I)*\[Xi]))*\[Delta]*\[Nu])/(10368*Sqrt[330])) + 
    SO*x^(5/2)*\[Epsilon]^5*((31*e^3*E^((3*I)*\[Xi])*\[Nu]*
        ((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*\[Chi]S))/
       (1944*Sqrt[330]) + (31*e^4*E^((2*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*
        \[Nu]*((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*\[Chi]S))/
       (648*Sqrt[330]) + (e^5*E^(I*\[Xi])*(651 - 1357*E^((2*I)*\[Xi]) + 
         1581*E^((4*I)*\[Xi]))*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(15552*Sqrt[330]) + 
      (e^6*(-124 - 579*E^((2*I)*\[Xi]) - 2211*E^((4*I)*\[Xi]) + 
         2914*E^((6*I)*\[Xi]))*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(15552*Sqrt[330])) + 
    \[Epsilon]^6*(x^3*(-(\[Delta]*\[Nu])/(189*Sqrt[330]) - 
        (e*(-123 + 134*E^((2*I)*\[Xi]))*\[Delta]*\[Nu])/
         (2016*Sqrt[330]*E^(I*\[Xi])) - (e^2*(-1531 - 1215*E^((2*I)*\[Xi]) + 
           3290*E^((4*I)*\[Xi]))*\[Delta]*\[Nu])/(10080*Sqrt[330]*
          E^((2*I)*\[Xi])) - (e^4*\[Delta]*\[Nu]*(-848172 - 
           406939*E^((2*I)*\[Xi]) - 348582*E^((4*I)*\[Xi]) + 
           50*E^((8*I)*\[Xi])*(-276347 + 485436*\[Nu]) - 5*E^((6*I)*\[Xi])*
            (-3050263 + 4854360*\[Nu])))/(1572480*Sqrt[330]*
          E^((4*I)*\[Xi])) - (e^3*\[Delta]*\[Nu]*(-2835027 - 
           1712178*E^((2*I)*\[Xi]) - 3948165*E^((4*I)*\[Xi]) + 
           10*E^((6*I)*\[Xi])*(-2515391 + 4865644*\[Nu])))/
         (9434880*Sqrt[330]*E^((3*I)*\[Xi])) - 
        (e^6*\[Delta]*\[Nu]*(-56661579 - 12008997*E^((2*I)*\[Xi]) - 
           16806426*E^((4*I)*\[Xi]) + E^((8*I)*\[Xi])*(728924394 - 
             1160358360*\[Nu]) + E^((6*I)*\[Xi])*(53182588 - 
             96297320*\[Nu]) - 12*E^((10*I)*\[Xi])*(-51887827 + 
             84664930*\[Nu]) + 8*E^((12*I)*\[Xi])*(-168374936 + 
             284079355*\[Nu])))/(37739520*Sqrt[330]*E^((6*I)*\[Xi])) - 
        (e^5*\[Delta]*\[Nu]*(-138314397 - 48121398*E^((2*I)*\[Xi]) - 
           51649962*E^((4*I)*\[Xi]) + 30*E^((10*I)*\[Xi])*
            (-96411197 + 164709720*\[Nu]) + 12*E^((6*I)*\[Xi])*
            (-122612419 + 169394820*\[Nu]) - 5*E^((8*I)*\[Xi])*
            (-412828047 + 659176112*\[Nu])))/(150958080*Sqrt[330]*
          E^((5*I)*\[Xi]))) + SO^2*x^3*
       ((-31*e^3*E^((3*I)*\[Xi])*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(23328*Sqrt[330]) - 
        (31*e^4*E^((2*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
          (9*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(23 - 108*\[Nu] + 64*\[Nu]^2)*
            \[Chi]A*\[Chi]S + \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (23 - 92*\[Nu])*\[Chi]A^2 + (23 - 32*\[Nu] + 8*\[Nu]^2)*
              \[Chi]S^2)))/(7776*Sqrt[330]) - 
        (e^5*E^(I*\[Xi])*(217 - 411*E^((2*I)*\[Xi]) + 527*E^((4*I)*\[Xi]))*
          \[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(23 - 108*\[Nu] + 64*\[Nu]^2)*
            \[Chi]A*\[Chi]S + \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (23 - 92*\[Nu])*\[Chi]A^2 + (23 - 32*\[Nu] + 8*\[Nu]^2)*
              \[Chi]S^2)))/(62208*Sqrt[330]) - 
        (e^6*(-124 - 951*E^((2*I)*\[Xi]) - 1839*E^((4*I)*\[Xi]) + 
           2914*E^((6*I)*\[Xi]))*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(186624*Sqrt[330])))
 
H\[Psi]OscMem[5, 5] = x^2*\[Epsilon]^4*
     ((5*e^5*E^((5*I)*\[Xi])*\[Delta]*\[Nu])/(1152*Sqrt[66]) + 
      (25*e^6*E^((4*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Delta]*\[Nu])/
       (1152*Sqrt[66])) + SO*x^(5/2)*\[Epsilon]^5*
     ((-5*e^5*E^((5*I)*\[Xi])*\[Nu]*((-2 + 8*\[Nu])*\[Chi]A + 
         \[Delta]*(-2 + \[Nu])*\[Chi]S))/(1728*Sqrt[66]) - 
      (25*e^6*E^((4*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
        ((-2 + 8*\[Nu])*\[Chi]A + \[Delta]*(-2 + \[Nu])*\[Chi]S))/
       (1728*Sqrt[66])) + \[Epsilon]^6*
     (x^3*((3*Sqrt[3/22]*\[Delta]*\[Nu])/35 + (e*(914 + 1773*E^((2*I)*\[Xi]))*
          \[Delta]*\[Nu])/(2016*Sqrt[66]*E^(I*\[Xi])) + 
        (e^2*(7870 + 8909*E^((2*I)*\[Xi]) + 21845*E^((4*I)*\[Xi]))*\[Delta]*
          \[Nu])/(10080*Sqrt[66]*E^((2*I)*\[Xi])) + 
        (e^3*(20724 + 19659*E^((2*I)*\[Xi]) + 19794*E^((4*I)*\[Xi]) + 
           74653*E^((6*I)*\[Xi]))*\[Delta]*\[Nu])/(16128*Sqrt[66]*
          E^((3*I)*\[Xi])) + (e^4*(248125 + 191760*E^((2*I)*\[Xi]) + 
           213882*E^((4*I)*\[Xi]) + 121805*E^((6*I)*\[Xi]) + 
           1105260*E^((8*I)*\[Xi]))*\[Delta]*\[Nu])/(120960*Sqrt[66]*
          E^((4*I)*\[Xi])) - (e^6*\[Delta]*\[Nu]*(-20655258 - 
           9236019*E^((2*I)*\[Xi]) - 11150165*E^((4*I)*\[Xi]) - 
           11854128*E^((6*I)*\[Xi]) - 14895270*E^((8*I)*\[Xi]) - 
           68*E^((10*I)*\[Xi])*(-451411 + 58240*\[Nu]) + 4*E^((12*I)*\[Xi])*
            (-32905349 + 990080*\[Nu])))/(4193280*Sqrt[66]*E^((6*I)*\[Xi])) - 
        (e^5*\[Delta]*\[Nu]*(-161250310 - 97871085*E^((2*I)*\[Xi]) - 
           110328400*E^((4*I)*\[Xi]) - 121370210*E^((6*I)*\[Xi]) + 
           46703670*E^((8*I)*\[Xi]) + E^((10*I)*\[Xi])*(-866689133 + 
             9395568*\[Nu])))/(50319360*Sqrt[66]*E^((5*I)*\[Xi]))) + 
      SO^2*x^3*((5*e^5*E^((5*I)*\[Xi])*\[Nu]*(9*\[Kappa]A*(-1 + 4*\[Nu]) + 
           2*(23 - 108*\[Nu] + 64*\[Nu]^2)*\[Chi]A*\[Chi]S + 
           \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + (23 - 92*\[Nu])*\[Chi]A^2 + 
             (23 - 32*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)))/(20736*Sqrt[66]) + 
        (25*e^6*E^((4*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
          (9*\[Kappa]A*(-1 + 4*\[Nu]) + 2*(23 - 108*\[Nu] + 64*\[Nu]^2)*
            \[Chi]A*\[Chi]S + \[Delta]*(9*\[Kappa]S*(-1 + 2*\[Nu]) + 
             (23 - 92*\[Nu])*\[Chi]A^2 + (23 - 32*\[Nu] + 8*\[Nu]^2)*
              \[Chi]S^2)))/(20736*Sqrt[66])))
 
H\[Psi]OscMem[6, 0] = 0
 
H\[Psi]OscMem[6, 2] = x^(5/2)*\[Epsilon]^5*
     (((-I/59136)*e^2*E^((2*I)*\[Xi])*\[Nu]*(-2783 + 8904*\[Nu]))/Sqrt[65] - 
      ((I/29568)*e^3*E^(I*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
        (-2783 + 8904*\[Nu]))/Sqrt[65] - 
      ((I/2128896)*e^5*(-1 + E^((2*I)*\[Xi]))*\[Nu]*(8349 - 26712*\[Nu] + 
         177*E^((4*I)*\[Xi])*(-2783 + 8904*\[Nu]) + 98*E^((2*I)*\[Xi])*
          (-7711 + 24888*\[Nu])))/(Sqrt[65]*E^(I*\[Xi])) - 
      ((I/2128896)*e^4*\[Nu]*(27*(-2783 + 8904*\[Nu]) + 117*E^((4*I)*\[Xi])*
          (-2783 + 8904*\[Nu]) + E^((2*I)*\[Xi])*(-68926 + 231168*\[Nu])))/
       Sqrt[65] - ((I/8515584)*e^6*\[Nu]*
        (E^((4*I)*\[Xi])*(1282418 - 4119864*\[Nu]) + 6*(-2783 + 8904*\[Nu]) + 
         1035*E^((8*I)*\[Xi])*(-2783 + 8904*\[Nu]) + 24*E^((2*I)*\[Xi])*
          (-60797 + 195846*\[Nu]) + 8*E^((6*I)*\[Xi])*
          (-222629 + 729582*\[Nu])))/(Sqrt[65]*E^((2*I)*\[Xi]))) + 
    SO*x^3*\[Epsilon]^6*(((I/88704)*e^2*E^((2*I)*\[Xi])*\[Nu]*
        (-2783 + 8904*\[Nu])*(-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/
       Sqrt[65] + ((I/44352)*e^3*E^(I*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
        (-2783 + 8904*\[Nu])*(-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/
       Sqrt[65] + ((I/3193344)*e^4*\[Nu]*(27*(-2783 + 8904*\[Nu]) + 
         117*E^((4*I)*\[Xi])*(-2783 + 8904*\[Nu]) + 20*E^((2*I)*\[Xi])*
          (-5951 + 19572*\[Nu]))*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu])*\[Chi]S))/Sqrt[65] + 
      ((I/3193344)*e^5*(-1 + E^((2*I)*\[Xi]))*\[Nu]*(8349 - 26712*\[Nu] + 
         177*E^((4*I)*\[Xi])*(-2783 + 8904*\[Nu]) + E^((2*I)*\[Xi])*
          (-855866 + 2759568*\[Nu]))*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu])*\[Chi]S))/(Sqrt[65]*E^(I*\[Xi])) + 
      ((I/12773376)*e^6*\[Nu]*(6*(-2783 + 8904*\[Nu]) + 1035*E^((8*I)*\[Xi])*
          (-2783 + 8904*\[Nu]) - 36*E^((4*I)*\[Xi])*(-27619 + 88242*\[Nu]) + 
         30*E^((2*I)*\[Xi])*(-53647 + 172704*\[Nu]) + E^((6*I)*\[Xi])*
          (-2432254 + 7920192*\[Nu]))*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu])*\[Chi]S))/(Sqrt[65]*E^((2*I)*\[Xi])))
 
H\[Psi]OscMem[6, 4] = x^(5/2)*\[Epsilon]^5*
     (((-I/101376)*e^4*E^((4*I)*\[Xi])*\[Nu]*(-703 + 2244*\[Nu]))/Sqrt[78] - 
      ((I/25344)*e^5*E^((3*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
        (-703 + 2244*\[Nu]))/Sqrt[78] - ((I/202752)*e^6*E^((2*I)*\[Xi])*\[Nu]*
        (-7733 + E^((2*I)*\[Xi])*(17656 - 56388*\[Nu]) + 24684*\[Nu] + 
         21*E^((4*I)*\[Xi])*(-703 + 2244*\[Nu])))/Sqrt[78]) + 
    SO*x^3*\[Epsilon]^6*(((I/152064)*e^4*E^((4*I)*\[Xi])*\[Nu]*
        (-703 + 2244*\[Nu])*(-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/
       Sqrt[78] + ((I/38016)*e^5*E^((3*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*\[Nu]*
        (-703 + 2244*\[Nu])*(-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/
       Sqrt[78] + ((I/304128)*e^6*E^((2*I)*\[Xi])*\[Nu]*
        (-7733 + E^((2*I)*\[Xi])*(16953 - 54144*\[Nu]) + 24684*\[Nu] + 
         21*E^((4*I)*\[Xi])*(-703 + 2244*\[Nu]))*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu])*\[Chi]S))/Sqrt[78])
 
H\[Psi]OscMem[6, 6] = (((25*I)/110592)*e^6*E^((6*I)*\[Xi])*x^(5/2)*
      \[Epsilon]^5*\[Nu]*(-19 + 64*\[Nu]))/Sqrt[143] - 
    (((25*I)/165888)*e^6*E^((6*I)*\[Xi])*SO*x^3*\[Epsilon]^6*\[Nu]*
      (-19 + 64*\[Nu])*(-2*\[Delta]*\[Chi]A + (-2 + \[Nu])*\[Chi]S))/Sqrt[143]
 
H\[Psi]OscMem[7, 1] = x^3*\[Epsilon]^6*
    ((-5*e*E^(I*\[Xi])*\[Delta]*\[Nu]*(-5023 + 16296*\[Nu]))/
      (5189184*Sqrt[2]) - (5*e^2*(-1 + E^((2*I)*\[Xi]))*\[Delta]*\[Nu]*
       (-5023 + 16296*\[Nu]))/(5189184*Sqrt[2]) - 
     (5*e^4*(-1 + E^((2*I)*\[Xi]))*\[Delta]*\[Nu]*(-5023 + 16296*\[Nu] + 
        16*E^((4*I)*\[Xi])*(-5023 + 16296*\[Nu]) + 70*E^((2*I)*\[Xi])*
         (-8143 + 24972*\[Nu])))/(62270208*Sqrt[2]*E^((2*I)*\[Xi])) - 
     (5*e^3*\[Delta]*\[Nu]*(5023 - 16296*\[Nu] + 9*E^((4*I)*\[Xi])*
         (-5023 + 16296*\[Nu]) + 18*E^((2*I)*\[Xi])*(-18693 + 56896*\[Nu])))/
      (41513472*Sqrt[2]*E^(I*\[Xi])) - (e^6*(-1 + E^((2*I)*\[Xi]))*\[Delta]*
       \[Nu]*(8*(-5023 + 16296*\[Nu]) + 243*E^((8*I)*\[Xi])*
         (-5023 + 16296*\[Nu]) + 7*E^((2*I)*\[Xi])*(-71207 + 218964*\[Nu]) + 
        28*E^((4*I)*\[Xi])*(-779033 + 2338491*\[Nu]) + 
        3*E^((6*I)*\[Xi])*(-2449103 + 7495656*\[Nu])))/
      (124540416*Sqrt[2]*E^((4*I)*\[Xi])) - 
     (5*e^5*\[Delta]*\[Nu]*(135621 + E^((2*I)*\[Xi])*
         (2219764 - 6796608*\[Nu]) - 439992*\[Nu] + 625*E^((8*I)*\[Xi])*
         (-5023 + 16296*\[Nu]) + 108*E^((6*I)*\[Xi])*
         (-158191 + 479472*\[Nu]) + 18*E^((4*I)*\[Xi])*
         (-2796539 + 8327368*\[Nu])))/(1992646656*Sqrt[2]*E^((3*I)*\[Xi])))
 
H\[Psi]OscMem[7, 3] = x^3*\[Epsilon]^6*
    ((-5*e^3*E^((3*I)*\[Xi])*\[Delta]*\[Nu]*(-12539 + 28896*\[Nu]))/
      (6918912*Sqrt[6]) - (5*e^4*E^((2*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*
       \[Delta]*\[Nu]*(-12539 + 28896*\[Nu]))/(2306304*Sqrt[6]) - 
     (5*e^6*(-1 + E^((2*I)*\[Xi]))*\[Delta]*\[Nu]*(-25078 + 57792*\[Nu] + 
        70*E^((2*I)*\[Xi])*(-4915 + 11532*\[Nu]) + 47*E^((4*I)*\[Xi])*
         (-12539 + 28896*\[Nu])))/(27675648*Sqrt[6]) - 
     (5*e^5*E^(I*\[Xi])*\[Delta]*\[Nu]*(-263319 + E^((2*I)*\[Xi])*
         (414302 - 945168*\[Nu]) + 606816*\[Nu] + 51*E^((4*I)*\[Xi])*
         (-12539 + 28896*\[Nu])))/(55351296*Sqrt[6]))
 
H\[Psi]OscMem[7, 5] = x^3*\[Epsilon]^6*
    (-(e^5*E^((5*I)*\[Xi])*\[Delta]*\[Nu]*(-6547 + 14664*\[Nu]))/
      (1797120*Sqrt[66]) - (e^6*E^((4*I)*\[Xi])*(-1 + E^((2*I)*\[Xi]))*
       \[Delta]*\[Nu]*(-6547 + 14664*\[Nu]))/(359424*Sqrt[66]))
 
H\[Psi]OscMem[7, 7] = 0
