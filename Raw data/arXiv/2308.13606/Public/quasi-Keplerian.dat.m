(* ::Package:: *)

(* ::Chapter:: *)
(*Supplementary file to "Spin effects in gravitational waveforms and fluxes for binaries on eccentric orbits to the third post-Newtonian order"*)


(* ::Text:: *)
(*by Quentin Henry and Mohammed Khalil*)


(* ::Text:: *)
(*This file contains relations that enter the quasi-Keplerian parametrization (QKP), in harmonic coordinates using the covariant spin-supplementary condition. Aligned-spin contributions are included to 3PN, and nonspinning contributions to 2PN.*)
(**)
(*The 3PN nonspinning QKP was derived in [Memmesheimer, Gopakumar, and Sch\[ADoubleDot]fer arXiv: gr-qc/0407049] in harmonic and ADM coordinates. The QKP for aligned spins was derived in [Tessmer, Hartung, and Sch\[ADoubleDot]fer arXiv: 1003.2735  and 1207.6961] to 3.5PN in ADM coordinates using the Newton-Wigner SSC.*)


(*
Notation:
We use geometric units with c = G = 1
\[Nu] = m1 m2 / M^2 is the symmetric mass ratio 
\[Delta] = (m1 - m2) / M is the antisymmetric mass ratio 
M = m1 + m2 is the total mass

\[Epsilon] = 1 is a PN counting parameter
SO = 1 is a spin order counting parameter

et, er, and e\[Phi] are the time, radial, and azimuthal eccentricities, respectively
\[ScriptL] is the mean anomaly
u is the eccentric anomaly
n is the mean motion
K is the periastron advance
rd and \[Phi]d denote dr/dt and d\[Phi]/dt
x = (M \[CapitalOmega])^(2/3) = (n \[CapitalKappa])^(2/3), where \[CapitalOmega] is the azimuthal orbital frequency

\[Chi]1 = S1 / m1^2, and \[Chi]2 = S2 / m2^2 are the dimensionless spins
\[Chi]S = 1/2 (\[Chi]1 + \[Chi]2)
\[Chi]A = 1/2 (\[Chi]1 - \[Chi]2)
\[Kappa]1 and \[Kappa]2 are the spin quadrupole constants, which equal 1 for black holes
\[Kappa]S = 1/2 ((\[Kappa]1 - 1) \[Chi]1^2 + (\[Kappa]2 - 1) \[Chi]2^2)
\[Kappa]A = 1/2 ((\[Kappa]1 - 1) \[Chi]1^2 - (\[Kappa]2 - 1) \[Chi]2^2)

We use dimensionless variables, given by
mE = -E / \[Mu], L = L / (M \[Mu]), r = r / M, \[Phi]d = \[Phi]d M
Note that mE and L are denoted \tilde{E} and h, respectively, in the paper.
*)


(* ::Section::Closed:: *)
(*energy and orbital angular momentum*)


Energy = (-2 + r*rd^2 + r^3*\[Phi]d^2)/(2*r) + 
    (\[Epsilon]^2*(4/r^2 + (4*rd^2*(3 + 2*\[Nu]))/r + 
       4*r*(3 + \[Nu])*\[Phi]d^2 - 3*(-1 + 3*\[Nu])*(rd^2 + r^2*\[Phi]d^2)^
         2))/8 + (\[Epsilon]^4*(-8 - 60*\[Nu] + 
       4*r*rd^2*(9 + 7*\[Nu] + 8*\[Nu]^2) + 2*r^3*(14 - 55*\[Nu] + 4*\[Nu]^2)*
        \[Phi]d^2 + 5*r^3*(1 - 7*\[Nu] + 13*\[Nu]^2)*(rd^2 + r^2*\[Phi]d^2)^
         3 - 2*r^2*(3*rd^4*(-7 + 8*\[Nu] + 16*\[Nu]^2) + 
         2*r^2*rd^2*(-21 + 22*\[Nu] + 42*\[Nu]^2)*\[Phi]d^2 + 
         r^4*(-21 + 23*\[Nu] + 27*\[Nu]^2)*\[Phi]d^4)))/(16*r^3) + 
    SO*((2*\[Epsilon]^3*\[Nu]*\[Phi]d*\[Chi]S)/r + 
      (\[Epsilon]^5*\[Phi]d*(\[Nu]*(\[Delta]*\[Chi]A + \[Chi]S + 
           4*\[Nu]*\[Chi]S) + 2*r^3*\[Phi]d^2*(\[Delta]*(-2 + \[Nu])*
            \[Chi]A + (-2 + 6*\[Nu] - 7*\[Nu]^2)*\[Chi]S) - 
         r*rd^2*(\[Delta]*(4 + \[Nu])*\[Chi]A + (4 - 9*\[Nu] + 8*\[Nu]^2)*
            \[Chi]S)))/(2*r^2)) + 
    SO^2*(-(\[Epsilon]^4*(\[Delta]*\[Kappa]A + \[Kappa]S - 
          2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
          2*\[Delta]*\[Chi]A*\[Chi]S + \[Chi]S^2))/(2*r^3) + 
      (\[Epsilon]^6*(r*rd^2*(\[Kappa]S*(3 - 6*\[Nu] - 8*\[Nu]^2) + 
           3*\[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 - 16*\[Nu]^2*\[Chi]A^2 + 
           3*\[Chi]S^2 - 8*\[Nu]*\[Chi]S^2 + 3*\[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S)) + r^3*\[Phi]d^2*
          (\[Kappa]S*(3 - 3*\[Nu] - 2*\[Nu]^2) + 3*\[Chi]A^2 - 
           13*\[Nu]*\[Chi]A^2 - 4*\[Nu]^2*\[Chi]A^2 + 3*\[Chi]S^2 - 
           17*\[Nu]*\[Chi]S^2 + 3*\[Delta]*(\[Kappa]A + \[Kappa]A*\[Nu] + 
             2*\[Chi]A*\[Chi]S - 6*\[Nu]*\[Chi]A*\[Chi]S)) + 
         2*(\[Kappa]S + \[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
           \[Chi]S^2 + 4*\[Nu]*\[Chi]S^2 + 4*\[Nu]^2*\[Chi]S^2 + 
           \[Delta]*(\[Kappa]A + 3*\[Kappa]A*\[Nu] + 2*(1 + 2*\[Nu])*\[Chi]A*
              \[Chi]S))))/(4*r^4))


AngMtm = r^2*\[Phi]d - (r*\[Epsilon]^2*\[Phi]d*(-2*(3 + \[Nu]) + 
       r*rd^2*(-1 + 3*\[Nu]) + r^3*(-1 + 3*\[Nu])*\[Phi]d^2))/2 + 
    (\[Epsilon]^4*\[Phi]d*(28 - 82*\[Nu] + 8*\[Nu]^2 - 
       4*r*rd^2*(-7 + 12*\[Nu] + 14*\[Nu]^2) - 
       4*r^3*(-7 + 10*\[Nu] + 9*\[Nu]^2)*\[Phi]d^2 + 
       3*r^2*(1 - 7*\[Nu] + 13*\[Nu]^2)*(rd^2 + r^2*\[Phi]d^2)^2))/8 + 
    (SO^2*\[Epsilon]^6*\[Phi]d*(3*\[Delta]*\[Kappa]A*(1 + \[Nu]) + 
       \[Kappa]S*(3 - 3*\[Nu] - 2*\[Nu]^2) + 3*\[Chi]A^2 - 
       13*\[Nu]*\[Chi]A^2 - 4*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(3 - 5*\[Nu])*
        \[Chi]A*\[Chi]S + 3*\[Chi]S^2 - 9*\[Nu]*\[Chi]S^2 + 
       8*\[Nu]^2*\[Chi]S^2))/(2*r) + 
    SO*((\[Epsilon]^3*(-2*\[Delta]*\[Chi]A + 
         (-2 + \[Nu]*(2 + r*rd^2 + r^3*\[Phi]d^2))*\[Chi]S))/r + 
      (\[Epsilon]^5*(\[Nu]*(rd^2 + r^2*\[Phi]d^2)^2*(\[Delta]*\[Chi]A + 
           (7 - 22*\[Nu])*\[Chi]S) + (4*\[Nu]*(\[Delta]*\[Chi]A + 
            (3 + 4*\[Nu])*\[Chi]S))/r^2 + 4*r*\[Phi]d^2*
          (3*\[Delta]*(-2 + \[Nu])*\[Chi]A + (-6 + 21*\[Nu] - 10*\[Nu]^2)*
            \[Chi]S) + (4*rd^2*(\[Delta]*(-2 + \[Nu])*\[Chi]A + 
            (-2 + 11*\[Nu] + 8*\[Nu]^2)*\[Chi]S))/r))/8)


(* ::Section::Closed:: *)
(*ar, er, et, and e\[Phi] as functions of E and L*)


ar\[LetterSpace]EL = 1/(2*mE) + (\[Epsilon]^2*(-7 + \[Nu]))/4 + 
    \[Epsilon]^4*((-4 + 7*\[Nu])/L^2 + (mE*(1 + \[Nu]^2))/8) + 
    SO*((2*\[Epsilon]^3*(\[Delta]*\[Chi]A + \[Chi]S - \[Nu]*\[Chi]S))/L + 
      (\[Epsilon]^5*(\[Delta]*(16 + 2*L^2*mE*(-1 + \[Nu]) - 5*\[Nu])*
          \[Chi]A + (16 - 2*L^2*mE*(-1 + \[Nu])^2 - 21*\[Nu] + 2*\[Nu]^2)*
          \[Chi]S))/L^3) + 
    SO^2*(-(\[Epsilon]^4*(\[Delta]*\[Kappa]A + \[Kappa]S - 
          2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
          2*\[Delta]*\[Chi]A*\[Chi]S + \[Chi]S^2))/(2*L^2) + 
      (\[Epsilon]^6*(\[Delta]*\[Kappa]A*(-18 + 4*L^2*mE + 5*\[Nu]) + 
         \[Kappa]S*(-18 + 4*L^2*mE*(-1 + \[Nu])^2 + 41*\[Nu] - 6*\[Nu]^2) - 
         42*\[Chi]A^2 + 4*L^2*mE*\[Chi]A^2 + 175*\[Nu]*\[Chi]A^2 - 
         20*L^2*mE*\[Nu]*\[Chi]A^2 - 12*\[Nu]^2*\[Chi]A^2 + 
         8*L^2*mE*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(L^2*mE*(4 - 8*\[Nu]) + 
           7*(-6 + 5*\[Nu]))*\[Chi]A*\[Chi]S - 42*\[Chi]S^2 + 
         4*L^2*mE*\[Chi]S^2 + 63*\[Nu]*\[Chi]S^2 - 12*L^2*mE*\[Nu]*
          \[Chi]S^2 - 24*\[Nu]^2*\[Chi]S^2 + 8*L^2*mE*\[Nu]^2*\[Chi]S^2))/
       (2*L^4))


er\[LetterSpace]EL = Sqrt[1 - 2*L^2*mE] + 
    (mE*\[Epsilon]^2*(-2*(-6 + \[Nu]) + 5*L^2*mE*(-3 + \[Nu])))/
     (2*Sqrt[1 - 2*L^2*mE]) + (mE*\[Epsilon]^4*(32*(4 - 7*\[Nu]) + 
       8*L^2*mE*(-35 + 99*\[Nu]) - 4*L^4*mE^2*(50 + 148*\[Nu] + \[Nu]^2) + 
       L^6*mE^3*(415 - 210*\[Nu] + 7*\[Nu]^2)))/
     (8*L^2*(1 - 2*L^2*mE)^(3/2)) + 
    SO*((4*mE*\[Epsilon]^3*(2*(-1 + L^2*mE)*\[Delta]*\[Chi]A + 
         L^2*mE*(2 - 3*\[Nu])*\[Chi]S + 2*(-1 + \[Nu])*\[Chi]S))/
       (L*Sqrt[1 - 2*L^2*mE]) - (mE*\[Epsilon]^5*
        (4*\[Delta]*(16 - 5*\[Nu])*\[Chi]A + 4*(16 - 21*\[Nu] + 2*\[Nu]^2)*
          \[Chi]S + 2*L^2*mE*(5*\[Delta]*(-18 + 7*\[Nu])*\[Chi]A + 
           (-90 + 125*\[Nu] - 14*\[Nu]^2)*\[Chi]S) - 
         2*L^6*mE^3*(5*\[Delta]*(-10 + \[Nu])*\[Chi]A + 
           (-50 + 70*\[Nu] - 7*\[Nu]^2)*\[Chi]S) + 3*L^4*mE^2*
          (\[Delta]*(12 - 19*\[Nu])*\[Chi]A + (12 - 29*\[Nu] + 6*\[Nu]^2)*
            \[Chi]S)))/(L^3*(1 - 2*L^2*mE)^(3/2))) + 
    SO^2*((-2*mE*(-1 + L^2*mE)*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
         \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
         \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/
       (L^2*Sqrt[1 - 2*L^2*mE]) - (mE*\[Epsilon]^6*
        (\[Kappa]S*(-36 + 82*\[Nu] - 12*\[Nu]^2 + L^4*mE^2*(-87 + 197*\[Nu] - 
             46*\[Nu]^2) + L^6*mE^3*(-13 + 33*\[Nu] + 2*\[Nu]^2) + 
           L^2*mE*(121 - 275*\[Nu] + 46*\[Nu]^2)) - 84*\[Chi]A^2 + 
         257*L^2*mE*\[Chi]A^2 - 119*L^4*mE^2*\[Chi]A^2 - 
         77*L^6*mE^3*\[Chi]A^2 + 350*\[Nu]*\[Chi]A^2 - 1083*L^2*mE*\[Nu]*
          \[Chi]A^2 + 535*L^4*mE^2*\[Nu]*\[Chi]A^2 + 299*L^6*mE^3*\[Nu]*
          \[Chi]A^2 - 24*\[Nu]^2*\[Chi]A^2 + 92*L^2*mE*\[Nu]^2*\[Chi]A^2 - 
         92*L^4*mE^2*\[Nu]^2*\[Chi]A^2 + 4*L^6*mE^3*\[Nu]^2*\[Chi]A^2 - 
         84*\[Chi]S^2 + 257*L^2*mE*\[Chi]S^2 - 119*L^4*mE^2*\[Chi]S^2 - 
         77*L^6*mE^3*\[Chi]S^2 + 126*\[Nu]*\[Chi]S^2 - 399*L^2*mE*\[Nu]*
          \[Chi]S^2 + 219*L^4*mE^2*\[Nu]*\[Chi]S^2 + 119*L^6*mE^3*\[Nu]*
          \[Chi]S^2 - 48*\[Nu]^2*\[Chi]S^2 + 152*L^2*mE*\[Nu]^2*\[Chi]S^2 - 
         76*L^4*mE^2*\[Nu]^2*\[Chi]S^2 - 64*L^6*mE^3*\[Nu]^2*\[Chi]S^2 + 
         \[Delta]*(\[Kappa]A*(-36 + 11*L^2*mE*(11 - 3*\[Nu]) + 10*\[Nu] + 
             L^6*mE^3*(-13 + 7*\[Nu]) + L^4*mE^2*(-87 + 23*\[Nu])) + 
           2*(-84 + L^2*mE*(257 - 227*\[Nu]) + 70*\[Nu] + 11*L^6*mE^3*
              (-7 + 5*\[Nu]) + L^4*mE^2*(-119 + 139*\[Nu]))*\[Chi]A*
            \[Chi]S)))/(L^4*(1 - 2*L^2*mE)^(3/2)))


et\[LetterSpace]EL = Sqrt[1 - 2*L^2*mE] - 
    (mE*\[Epsilon]^2*(4 - 4*\[Nu] + L^2*mE*(-17 + 7*\[Nu])))/
     (2*Sqrt[1 - 2*L^2*mE]) + (mE*\[Epsilon]^4*(64 - 112*\[Nu] + 
       24*Sqrt[2]*L*Sqrt[mE]*(-5 + 2*\[Nu]) - 96*Sqrt[2]*L^3*mE^(3/2)*
        (-5 + 2*\[Nu]) + 96*Sqrt[2]*L^5*mE^(5/2)*(-5 + 2*\[Nu]) + 
       8*L^2*mE*(-15 + 50*\[Nu] + 3*\[Nu]^2) - 
       4*L^4*mE^2*(90 + 73*\[Nu] + 22*\[Nu]^2) + 
       L^6*mE^3*(607 - 138*\[Nu] + 79*\[Nu]^2)))/
     (8*L^2*(1 - 2*L^2*mE)^(3/2)) + 
    SO*((-4*mE*\[Epsilon]^3*(\[Delta]*\[Chi]A + \[Chi]S + 
         (-1 + L^2*mE)*\[Nu]*\[Chi]S))/(L*Sqrt[1 - 2*L^2*mE]) - 
      (mE*\[Epsilon]^5*(\[Delta]*(32 + 3*L^4*mE^2*(30 - 17*\[Nu]) + 
           8*Sqrt[2]*L*Sqrt[mE]*(-3 + \[Nu]) - 32*Sqrt[2]*L^3*mE^(3/2)*
            (-3 + \[Nu]) + 32*Sqrt[2]*L^5*mE^(5/2)*(-3 + \[Nu]) - 10*\[Nu] + 
           2*L^6*mE^3*\[Nu] + 2*L^2*mE*(-59 + 22*\[Nu]))*\[Chi]A + 
         (32 - 42*\[Nu] + 2*L^6*mE^3*(18 - 11*\[Nu])*\[Nu] + 4*\[Nu]^2 - 
           4*Sqrt[2]*L*Sqrt[mE]*(6 - 8*\[Nu] + \[Nu]^2) + 
           16*Sqrt[2]*L^3*mE^(3/2)*(6 - 8*\[Nu] + \[Nu]^2) - 
           16*Sqrt[2]*L^5*mE^(5/2)*(6 - 8*\[Nu] + \[Nu]^2) - 
           2*L^2*mE*(59 - 81*\[Nu] + 13*\[Nu]^2) + L^4*mE^2*
            (90 - 167*\[Nu] + 48*\[Nu]^2))*\[Chi]S))/
       (L^3*(1 - 2*L^2*mE)^(3/2))) + 
    SO^2*((mE*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
         4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
           2*\[Chi]A*\[Chi]S)))/(L^2*Sqrt[1 - 2*L^2*mE]) - 
      (mE*\[Epsilon]^6*(\[Kappa]S*(-36 + 82*\[Nu] - 12*\[Nu]^2 + 
           L^4*mE^2*(-97 + 225*\[Nu] - 62*\[Nu]^2) + 2*Sqrt[2]*L*Sqrt[mE]*
            (14 - 33*\[Nu] + 6*\[Nu]^2) - 8*Sqrt[2]*L^3*mE^(3/2)*
            (14 - 33*\[Nu] + 6*\[Nu]^2) + 8*Sqrt[2]*L^5*mE^(5/2)*
            (14 - 33*\[Nu] + 6*\[Nu]^2) + L^2*mE*(125 - 285*\[Nu] + 
             54*\[Nu]^2)) - 84*\[Chi]A^2 + 44*Sqrt[2]*L*Sqrt[mE]*\[Chi]A^2 + 
         269*L^2*mE*\[Chi]A^2 - 176*Sqrt[2]*L^3*mE^(3/2)*\[Chi]A^2 - 
         161*L^4*mE^2*\[Chi]A^2 + 176*Sqrt[2]*L^5*mE^(5/2)*\[Chi]A^2 + 
         350*\[Nu]*\[Chi]A^2 - 184*Sqrt[2]*L*Sqrt[mE]*\[Nu]*\[Chi]A^2 - 
         1131*L^2*mE*\[Nu]*\[Chi]A^2 + 736*Sqrt[2]*L^3*mE^(3/2)*\[Nu]*
          \[Chi]A^2 + 699*L^4*mE^2*\[Nu]*\[Chi]A^2 - 736*Sqrt[2]*L^5*mE^(5/2)*
          \[Nu]*\[Chi]A^2 - 24*\[Nu]^2*\[Chi]A^2 + 24*Sqrt[2]*L*Sqrt[mE]*
          \[Nu]^2*\[Chi]A^2 + 108*L^2*mE*\[Nu]^2*\[Chi]A^2 - 
         96*Sqrt[2]*L^3*mE^(3/2)*\[Nu]^2*\[Chi]A^2 - 124*L^4*mE^2*\[Nu]^2*
          \[Chi]A^2 + 96*Sqrt[2]*L^5*mE^(5/2)*\[Nu]^2*\[Chi]A^2 - 
         84*\[Chi]S^2 + 44*Sqrt[2]*L*Sqrt[mE]*\[Chi]S^2 + 
         269*L^2*mE*\[Chi]S^2 - 176*Sqrt[2]*L^3*mE^(3/2)*\[Chi]S^2 - 
         161*L^4*mE^2*\[Chi]S^2 + 176*Sqrt[2]*L^5*mE^(5/2)*\[Chi]S^2 + 
         126*\[Nu]*\[Chi]S^2 - 64*Sqrt[2]*L*Sqrt[mE]*\[Nu]*\[Chi]S^2 - 
         407*L^2*mE*\[Nu]*\[Chi]S^2 + 256*Sqrt[2]*L^3*mE^(3/2)*\[Nu]*
          \[Chi]S^2 + 279*L^4*mE^2*\[Nu]*\[Chi]S^2 - 256*Sqrt[2]*L^5*mE^(5/2)*
          \[Nu]*\[Chi]S^2 - 48*\[Nu]^2*\[Chi]S^2 + 16*Sqrt[2]*L*Sqrt[mE]*
          \[Nu]^2*\[Chi]S^2 + 160*L^2*mE*\[Nu]^2*\[Chi]S^2 - 
         64*Sqrt[2]*L^3*mE^(3/2)*\[Nu]^2*\[Chi]S^2 - 120*L^4*mE^2*\[Nu]^2*
          \[Chi]S^2 + 64*Sqrt[2]*L^5*mE^(5/2)*\[Nu]^2*\[Chi]S^2 + 
         \[Delta]*(\[Kappa]A*(-36 + 2*Sqrt[2]*L*Sqrt[mE]*(14 - 5*\[Nu]) + 
             10*\[Nu] + 8*Sqrt[2]*L^3*mE^(3/2)*(-14 + 5*\[Nu]) - 
             8*Sqrt[2]*L^5*mE^(5/2)*(-14 + 5*\[Nu]) - 5*L^2*mE*
              (-25 + 7*\[Nu]) + L^4*mE^2*(-97 + 31*\[Nu])) + 
           2*(-84 + L^2*mE*(269 - 231*\[Nu]) + 4*Sqrt[2]*L*Sqrt[mE]*
              (11 - 9*\[Nu]) + 70*\[Nu] + 16*Sqrt[2]*L^3*mE^(3/2)*
              (-11 + 9*\[Nu]) - 16*Sqrt[2]*L^5*mE^(5/2)*(-11 + 9*\[Nu]) + 
             L^4*mE^2*(-161 + 167*\[Nu]))*\[Chi]A*\[Chi]S)))/
       (2*L^4*(1 - 2*L^2*mE)^(3/2)))


e\[Phi]\[LetterSpace]EL = Sqrt[1 - 2*L^2*mE] + 
    (mE*\[Epsilon]^2*(12 + L^2*mE*(-15 + \[Nu])))/(2*Sqrt[1 - 2*L^2*mE]) + 
    (mE*\[Epsilon]^4*(416 - 91*\[Nu] - 15*\[Nu]^2 - 
       12*L^4*mE^2*(-20 + 5*\[Nu] + 7*\[Nu]^2) + 
       2*L^6*mE^3*(415 - 94*\[Nu] + 11*\[Nu]^2) + 
       2*L^2*mE*(-600 + 125*\[Nu] + 33*\[Nu]^2)))/
     (16*L^2*(1 - 2*L^2*mE)^(3/2)) + 
    SO*((4*mE*(-1 + L^2*mE)*\[Epsilon]^3*(2*\[Delta]*\[Chi]A - 
         (-2 + \[Nu])*\[Chi]S))/(L*Sqrt[1 - 2*L^2*mE]) - 
      (mE*\[Epsilon]^5*(\[Delta]*(384 - 67*\[Nu] + 16*L^6*mE^3*(25 + \[Nu]) - 
           8*L^4*mE^2*(-82 + 33*\[Nu]) + 2*L^2*mE*(-616 + 129*\[Nu]))*
          \[Chi]A + (384 - 419*\[Nu] + 6*\[Nu]^2 - 8*L^6*mE^3*
            (-50 + 13*\[Nu] + \[Nu]^2) - 14*L^2*mE*(88 - 107*\[Nu] + 
             2*\[Nu]^2) + 8*L^4*mE^2*(82 - 155*\[Nu] + 5*\[Nu]^2))*\[Chi]S))/
       (4*L^3*(1 - 2*L^2*mE)^(3/2))) + 
    SO^2*(-((mE*(-3 + 4*L^2*mE)*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
          \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
          \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/
        (L^2*Sqrt[1 - 2*L^2*mE])) + 
      (mE*\[Epsilon]^6*(\[Kappa]S*(528 - 1261*\[Nu] + 102*\[Nu]^2 - 
           8*L^6*mE^3*(-26 + 41*\[Nu] + 18*\[Nu]^2) - 
           6*L^2*mE*(314 - 763*\[Nu] + 74*\[Nu]^2) + 4*L^4*mE^2*
            (379 - 970*\[Nu] + 140*\[Nu]^2)) + 1104*\[Chi]A^2 - 
         3740*L^2*mE*\[Chi]A^2 + 2540*L^4*mE^2*\[Chi]A^2 + 
         720*L^6*mE^3*\[Chi]A^2 - 4505*\[Nu]*\[Chi]A^2 + 
         15314*L^2*mE*\[Nu]*\[Chi]A^2 - 10544*L^4*mE^2*\[Nu]*\[Chi]A^2 - 
         2824*L^6*mE^3*\[Nu]*\[Chi]A^2 + 204*\[Nu]^2*\[Chi]A^2 - 
         888*L^2*mE*\[Nu]^2*\[Chi]A^2 + 1120*L^4*mE^2*\[Nu]^2*\[Chi]A^2 - 
         288*L^6*mE^3*\[Nu]^2*\[Chi]A^2 + 1104*\[Chi]S^2 - 
         3740*L^2*mE*\[Chi]S^2 + 2540*L^4*mE^2*\[Chi]S^2 + 
         720*L^6*mE^3*\[Chi]S^2 - 1369*\[Nu]*\[Chi]S^2 + 
         5026*L^2*mE*\[Nu]*\[Chi]S^2 - 4480*L^4*mE^2*\[Nu]*\[Chi]S^2 + 
         56*L^6*mE^3*\[Nu]*\[Chi]S^2 + 320*\[Nu]^2*\[Chi]S^2 - 
         1152*L^2*mE*\[Nu]^2*\[Chi]S^2 + 992*L^4*mE^2*\[Nu]^2*\[Chi]S^2 + 
         \[Delta]*(\[Kappa]A*(528 + 4*L^4*mE^2*(379 - 212*\[Nu]) - 
             205*\[Nu] + 8*L^6*mE^3*(26 + 11*\[Nu]) + 6*L^2*mE*
              (-314 + 135*\[Nu])) + 2*(1104 - 729*\[Nu] + 
             8*L^6*mE^3*(90 + 7*\[Nu]) + 10*L^2*mE*(-374 + 269*\[Nu]) - 
             4*L^4*mE^2*(-635 + 608*\[Nu]))*\[Chi]A*\[Chi]S)))/
       (8*L^4*(1 - 2*L^2*mE)^(3/2)))


(* ::Section::Closed:: *)
(*relations between er, et, and e\[Phi]*)


(* ::Text:: *)
(*et(er, ar)*)


et\[LetterSpace]er = er + (er*\[Epsilon]^2*(-8 + 3*\[Nu]))/(2*ar) + 
    (er*\[Epsilon]^4*(-96 + 60*Sqrt[1 - er^2] + (11 - 24*Sqrt[1 - er^2])*
        \[Nu] - 15*\[Nu]^2 + er^2*(128 - 67*\[Nu] + 15*\[Nu]^2)))/
     (8*ar^2*(-1 + er^2)) + 
    SO*((2*er*\[Epsilon]^3*(\[Delta]*\[Chi]A + \[Chi]S - \[Nu]*\[Chi]S))/
       Sqrt[-(ar^3*(-1 + er^2))] - ((-1 + er)*er*(1 + er)*\[Epsilon]^5*
        (-2*Sqrt[ar]*(-1 + er^2)*(2*\[Delta]*(-3 + \[Nu])*\[Chi]A - 
           (6 - 8*\[Nu] + \[Nu]^2)*\[Chi]S) + Sqrt[ar - ar*er^2]*
          ((-1 + 3*er^2)*\[Delta]*(-5 + 2*\[Nu])*\[Chi]A + 
           (5 - 7*\[Nu] + 5*\[Nu]^2 - 3*er^2*(5 - 7*\[Nu] + 2*\[Nu]^2))*
            \[Chi]S)))/(ar^3*(-1 + er^2)^3)) + 
    SO^2*((er*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
         4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
           2*\[Chi]A*\[Chi]S)))/(2*ar^2*(-1 + er^2)) - 
      (er*\[Epsilon]^6*(\[Kappa]S*(4 + 14*Sqrt[1 - er^2] - 
           (10 + 33*Sqrt[1 - er^2])*\[Nu] + (-2 + 6*Sqrt[1 - er^2])*\[Nu]^2 + 
           2*er^2*(4 - 9*\[Nu] + 3*\[Nu]^2)) + \[Delta]*
          (\[Kappa]A*(4 + 14*Sqrt[1 - er^2] - 2*er^2*(-4 + \[Nu]) - 2*\[Nu] - 
             5*Sqrt[1 - er^2]*\[Nu]) - 4*(-4 - 11*Sqrt[1 - er^2] + 5*\[Nu] + 
             9*Sqrt[1 - er^2]*\[Nu] + er^2*(-6 + 4*\[Nu]))*\[Chi]A*\[Chi]S) + 
         2*((4 + 11*Sqrt[1 - er^2] - (17 + 46*Sqrt[1 - er^2])*\[Nu] + 
             (-2 + 6*Sqrt[1 - er^2])*\[Nu]^2 + er^2*(6 - 26*\[Nu] + 6*
                \[Nu]^2))*\[Chi]A^2 + (4 + 11*Sqrt[1 - er^2] - 
             (9 + 16*Sqrt[1 - er^2])*\[Nu] + (2 + 4*Sqrt[1 - er^2])*\[Nu]^2 + 
             2*er^2*(3 - 3*\[Nu] + \[Nu]^2))*\[Chi]S^2)))/
       (2*ar^3*(-1 + er^2)^2))


(* ::Text:: *)
(*e\[Phi](er, ar)*)


e\[Phi]\[LetterSpace]er = er + (er*\[Epsilon]^2*\[Nu])/(2*ar) - 
    (er*\[Epsilon]^4*(160 + (328 + 29*er^2)*\[Nu] - 15*er^2*\[Nu]^2))/
     (32*ar^2*(-1 + er^2)) + 
    SO*((2*er*Sqrt[ar - ar*er^2]*\[Epsilon]^3*\[Nu]*\[Chi]S)/
       (ar^2*(-1 + er^2)) - (er*Sqrt[ar - ar*er^2]*\[Epsilon]^5*
        (\[Delta]*(128 + (4 + 9*er^2)*\[Nu])*\[Chi]A + 
         (128 + (-132 + er^2)*\[Nu] - 18*er^2*\[Nu]^2)*\[Chi]S))/
       (8*ar^3*(-1 + er^2)^2)) + 
    SO^2*(-(er*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
          4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
            2*\[Chi]A*\[Chi]S)))/(2*ar^2*(-1 + er^2)) - 
      (er*\[Epsilon]^6*(\[Kappa]S*(-192 + (492 + 9*er^2)*\[Nu] + 
           2*(4 + er^2)*\[Nu]^2) - 384*\[Chi]A^2 + 1516*\[Nu]*\[Chi]A^2 - 
         11*er^2*\[Nu]*\[Chi]A^2 + 16*\[Nu]^2*\[Chi]A^2 + 
         4*er^2*\[Nu]^2*\[Chi]A^2 - 384*\[Chi]S^2 + 460*\[Nu]*\[Chi]S^2 + 
         21*er^2*\[Nu]*\[Chi]S^2 - 96*\[Nu]^2*\[Chi]S^2 + 
         32*er^2*\[Nu]^2*\[Chi]S^2 + \[Delta]*
          (3*\[Kappa]A*(-64 + 3*(12 + er^2)*\[Nu]) + 
           2*(-384 + 5*(44 + er^2)*\[Nu])*\[Chi]A*\[Chi]S)))/
       (16*ar^3*(-1 + er^2)^2))


(* ::Text:: *)
(*et / er as a function of (E,L)*)


et\[LetterSpace]erEL = er (et\[LetterSpace]erEL = 1 + mE*\[Epsilon]^2*(-8 + 3*\[Nu]) + 
    (mE*\[Epsilon]^4*(-8 + 14*\[Nu] + 3*Sqrt[2]*L*Sqrt[mE]*(-5 + 2*\[Nu]) + 
       L^2*mE*(36 - 19*\[Nu] + 6*\[Nu]^2)))/L^2 + 
    SO*((4*mE*\[Epsilon]^3*(\[Delta]*\[Chi]A + \[Chi]S - \[Nu]*\[Chi]S))/L + 
      (2*mE*\[Epsilon]^5*(\[Delta]*(16 - 4*Sqrt[2]*L*Sqrt[mE]*(-3 + \[Nu]) - 
           5*\[Nu] + L^2*mE*(-11 + 7*\[Nu]))*\[Chi]A + 
         (16 - 21*\[Nu] + 2*\[Nu]^2 + L^2*mE*(-11 + 18*\[Nu] - 7*\[Nu]^2) + 
           2*Sqrt[2]*L*Sqrt[mE]*(6 - 8*\[Nu] + \[Nu]^2))*\[Chi]S))/L^3) + 
    SO^2*(-((mE*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
          4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
            2*\[Chi]A*\[Chi]S)))/L^2) - 
      (mE*\[Epsilon]^6*(\[Kappa]S*(36 - 82*\[Nu] + 12*\[Nu]^2 + 
           L^2*mE*(-17 + 39*\[Nu] - 18*\[Nu]^2) + 2*Sqrt[2]*L*Sqrt[mE]*
            (14 - 33*\[Nu] + 6*\[Nu]^2)) + 84*\[Chi]A^2 + 
         44*Sqrt[2]*L*Sqrt[mE]*\[Chi]A^2 - 17*L^2*mE*\[Chi]A^2 - 
         350*\[Nu]*\[Chi]A^2 - 184*Sqrt[2]*L*Sqrt[mE]*\[Nu]*\[Chi]A^2 + 
         81*L^2*mE*\[Nu]*\[Chi]A^2 + 24*\[Nu]^2*\[Chi]A^2 + 
         24*Sqrt[2]*L*Sqrt[mE]*\[Nu]^2*\[Chi]A^2 - 36*L^2*mE*\[Nu]^2*
          \[Chi]A^2 + 84*\[Chi]S^2 + 44*Sqrt[2]*L*Sqrt[mE]*\[Chi]S^2 - 
         17*L^2*mE*\[Chi]S^2 - 126*\[Nu]*\[Chi]S^2 - 64*Sqrt[2]*L*Sqrt[mE]*
          \[Nu]*\[Chi]S^2 + 29*L^2*mE*\[Nu]*\[Chi]S^2 + 
         48*\[Nu]^2*\[Chi]S^2 + 16*Sqrt[2]*L*Sqrt[mE]*\[Nu]^2*\[Chi]S^2 - 
         16*L^2*mE*\[Nu]^2*\[Chi]S^2 + \[Delta]*
          (\[Kappa]A*(36 + 2*Sqrt[2]*L*Sqrt[mE]*(14 - 5*\[Nu]) - 10*\[Nu] + 
             L^2*mE*(-17 + 5*\[Nu])) + 2*(84 + 4*Sqrt[2]*L*Sqrt[mE]*
              (11 - 9*\[Nu]) - 70*\[Nu] + L^2*mE*(-17 + 21*\[Nu]))*\[Chi]A*
            \[Chi]S)))/(2*L^4)))


(* ::Text:: *)
(*e\[Phi] / er as a function of (E,L)*)


e\[Phi]\[LetterSpace]erEL = er (1 + mE*\[Epsilon]^2*\[Nu] + 
    (mE*\[Epsilon]^4*(160 + (357 - 2*L^2*mE)*\[Nu] + (-15 + 22*L^2*mE)*
        \[Nu]^2))/(16*L^2) + SO*((-4*mE*\[Epsilon]^3*\[Nu]*\[Chi]S)/L + 
      (mE*\[Epsilon]^5*(\[Delta]*(-128 + (-13 + 2*L^2*mE)*\[Nu])*\[Chi]A + 
         (-128 + (83 - 38*L^2*mE)*\[Nu] + (26 - 28*L^2*mE)*\[Nu]^2)*\[Chi]S))/
       (4*L^3)) + SO^2*((mE*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
         \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
         \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/L^2 + 
      (mE*\[Epsilon]^6*(\[Kappa]S*(240 - 605*\[Nu] + 6*\[Nu]^2 - 
           2*L^2*mE*(2 - 23*\[Nu] + 18*\[Nu]^2)) + 432*\[Chi]A^2 - 
         4*L^2*mE*\[Chi]A^2 - 1705*\[Nu]*\[Chi]A^2 + 14*L^2*mE*\[Nu]*
          \[Chi]A^2 + 12*\[Nu]^2*\[Chi]A^2 - 72*L^2*mE*\[Nu]^2*\[Chi]A^2 + 
         432*\[Chi]S^2 - 4*L^2*mE*\[Chi]S^2 - 361*\[Nu]*\[Chi]S^2 + 
         126*L^2*mE*\[Nu]*\[Chi]S^2 - 64*\[Nu]^2*\[Chi]S^2 + 
         64*L^2*mE*\[Nu]^2*\[Chi]S^2 + \[Delta]*
          (\[Kappa]A*(240 - 125*\[Nu] + 2*L^2*mE*(-2 + 19*\[Nu])) + 
           2*(432 - 169*\[Nu] + 2*L^2*mE*(-2 + 31*\[Nu]))*\[Chi]A*\[Chi]S)))/
       (8*L^4)))


(* ::Section::Closed:: *)
(*mean motion and periastron advance*)


(* ::Text:: *)
(*mean motion n, and periastron advance K are provided as functions of E and L, and as functions of x and et*)


n\[LetterSpace]EL = 2*Sqrt[2]*mE^(3/2) + 
    (mE^(5/2)*\[Epsilon]^2*(-15 + \[Nu]))/Sqrt[2] + 
    \[Epsilon]^4*((12*mE^3*(-5 + 2*\[Nu]))/L + 
      (mE^(7/2)*(555 + 30*\[Nu] + 11*\[Nu]^2))/(8*Sqrt[2])) + 
    (32*mE^3*Sqrt[1 - 2*L^2*mE]*(-1 - Sqrt[1 - 2*L^2*mE] + 
       L^2*mE*(2 + Sqrt[1 - 2*L^2*mE]))*SO*\[Epsilon]^5*
      (2*\[Delta]*(-3 + \[Nu])*\[Chi]A - (6 - 8*\[Nu] + \[Nu]^2)*\[Chi]S))/
     (L^2*(1 - 2*L^2*mE + Sqrt[1 - 2*L^2*mE])^2) - 
    (8*mE^3*Sqrt[1 - 2*L^2*mE]*(-1 - Sqrt[1 - 2*L^2*mE] + 
       L^2*mE*(2 + Sqrt[1 - 2*L^2*mE]))*SO^2*\[Epsilon]^6*
      (\[Delta]*\[Kappa]A*(-14 + 5*\[Nu]) + \[Kappa]S*(-14 + 33*\[Nu] - 
         6*\[Nu]^2) + 4*\[Delta]*(-11 + 9*\[Nu])*\[Chi]A*\[Chi]S - 
       2*((11 - 46*\[Nu] + 6*\[Nu]^2)*\[Chi]A^2 + (11 - 16*\[Nu] + 4*\[Nu]^2)*
          \[Chi]S^2)))/(L^3*(1 - 2*L^2*mE + Sqrt[1 - 2*L^2*mE])^2)


n\[LetterSpace]xe = x^(3/2) + (3*x^(5/2)*\[Epsilon]^2)/(-1 + et^2) + 
    (x^(7/2)*\[Epsilon]^4*(-18 + 28*\[Nu] + et^2*(-51 + 26*\[Nu])))/
     (4*(-1 + et^2)^2) + 
    SO*((2*x^3*\[Epsilon]^3*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S))/
       (1 - et^2)^(3/2) - (x^4*\[Epsilon]^5*
        (\[Delta]*(-20 + 17*\[Nu] + 4*et^2*(-15 + 7*\[Nu]))*\[Chi]A - 
         (20 - 57*\[Nu] + 4*\[Nu]^2 + 2*et^2*(30 - 29*\[Nu] + 7*\[Nu]^2))*
          \[Chi]S))/(2*(1 - et^2)^(5/2))) + 
    SO^2*((-3*x^(7/2)*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
         \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
         \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/(2*(-1 + et^2)^2) + 
      (x^(9/2)*\[Epsilon]^6*(\[Kappa]S*(42 - 118*\[Nu] + 20*\[Nu]^2 + 
           et^2*(78 - 191*\[Nu] + 58*\[Nu]^2)) + 
         \[Delta]*(\[Kappa]A*(42 + et^2*(78 - 35*\[Nu]) - 34*\[Nu]) + 
           4*(17 + et^2*(51 - 25*\[Nu]) - 39*\[Nu])*\[Chi]A*\[Chi]S) + 
         2*((17 - 79*\[Nu] + 20*\[Nu]^2 + et^2*(51 - 220*\[Nu] + 58*\[Nu]^2))*
            \[Chi]A^2 + (17 + et^2*(51 - 34*\[Nu]) - 67*\[Nu] + 20*\[Nu]^2)*
            \[Chi]S^2)))/(4*(-1 + et^2)^3))


K\[LetterSpace]EL = 1 + (3*\[Epsilon]^2)/L^2 + 
    (\[Epsilon]^4*(105 - 30*\[Nu] + 6*L^2*mE*(-5 + 2*\[Nu])))/(4*L^4) + 
    SO*((\[Epsilon]^3*(-4*\[Delta]*\[Chi]A + 2*(-2 + \[Nu])*\[Chi]S))/L^3 + 
      (\[Epsilon]^5*(21*\[Delta]*(-8 + \[Nu])*\[Chi]A - 16*L^2*mE*\[Delta]*
          (-3 + \[Nu])*\[Chi]A + 8*L^2*mE*(6 - 8*\[Nu] + \[Nu]^2)*\[Chi]S - 
         3*(56 - 49*\[Nu] + 2*\[Nu]^2)*\[Chi]S))/(2*L^5)) + 
    SO^2*((3*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
         4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
           2*\[Chi]A*\[Chi]S)))/(2*L^4) + 
      (3*\[Epsilon]^6*(5*\[Kappa]S*(12 - 27*\[Nu] + 2*\[Nu]^2) - 
         2*L^2*mE*\[Kappa]S*(14 - 33*\[Nu] + 6*\[Nu]^2) - 
         2*(-5*(14 - 57*\[Nu] + 2*\[Nu]^2) + 2*L^2*mE*(11 - 46*\[Nu] + 
             6*\[Nu]^2))*\[Chi]A^2 - 2*(-70 + 75*\[Nu] - 20*\[Nu]^2 + 
           2*L^2*mE*(11 - 16*\[Nu] + 4*\[Nu]^2))*\[Chi]S^2 + 
         \[Delta]*(\[Kappa]A*(-15*(-4 + \[Nu]) + 2*L^2*mE*(-14 + 5*\[Nu])) + 
           8*(35 - 20*\[Nu] + L^2*mE*(-11 + 9*\[Nu]))*\[Chi]A*\[Chi]S)))/
       (4*L^6))


K\[LetterSpace]xe = 1 - (3*x*\[Epsilon]^2)/(-1 + et^2) - 
    (x^2*\[Epsilon]^4*(-54 + 28*\[Nu] + et^2*(-51 + 26*\[Nu])))/
     (4*(-1 + et^2)^2) + SO*((-(x/(-1 + et^2)))^(3/2)*\[Epsilon]^3*
       (-4*\[Delta]*\[Chi]A + 2*(-2 + \[Nu])*\[Chi]S) + 
      ((-(x/(-1 + et^2)))^(5/2)*\[Epsilon]^5*
        (\[Delta]*(17*(-4 + \[Nu]) + 4*et^2*(-15 + 7*\[Nu]))*\[Chi]A - 
         (68 - 81*\[Nu] + 4*\[Nu]^2 + 2*et^2*(30 - 29*\[Nu] + 7*\[Nu]^2))*
          \[Chi]S))/2) + 
    SO^2*((3*x^2*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
         4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
           2*\[Chi]A*\[Chi]S)))/(2*(-1 + et^2)^2) - 
      (x^3*\[Epsilon]^6*(\[Kappa]S*(78 - 190*\[Nu] + 20*\[Nu]^2 + 
           et^2*(78 - 191*\[Nu] + 58*\[Nu]^2)) + 
         \[Delta]*(\[Kappa]A*(78 + et^2*(78 - 35*\[Nu]) - 34*\[Nu]) - 
           4*(-67 + 55*\[Nu] + et^2*(-51 + 25*\[Nu]))*\[Chi]A*\[Chi]S) + 
         2*((67 - 279*\[Nu] + 20*\[Nu]^2 + et^2*(51 - 220*\[Nu] + 58*
                \[Nu]^2))*\[Chi]A^2 + (67 + et^2*(51 - 34*\[Nu]) - 99*\[Nu] + 
             28*\[Nu]^2)*\[Chi]S^2)))/(4*(-1 + et^2)^3))


(* ::Section::Closed:: *)
(*relations for \[CapitalDelta]\[Phi]/K*)


(* ::Text:: *)
(*2 \[Pi] (\[Phi] - \[Phi]0) / \[CapitalPhi] = \[CapitalDelta]\[Phi] / K = v + \[ScriptCapitalG]2 Sin[2 v] + \[ScriptCapitalG]3 Sin[3 v]*)


\[ScriptCapitalG]2 = (er^2*\[Epsilon]^4*(2 + 38*\[Nu] - 6*\[Nu]^2))/
     (16*ar^2*(-1 + er)^2*(1 + er)^2) + 
    (er^2*SO*\[Epsilon]^5*\[Nu]*(-5*\[Delta]*\[Chi]A + 
       (-7 + 6*\[Nu])*\[Chi]S))/(4*ar^2*(-1 + er)^2*(1 + er)^2*
      Sqrt[ar - ar*er^2]) + (er^2*SO^2*\[Epsilon]^6*
      (\[Delta]*\[Kappa]A*(-6 + 7*\[Nu]) + \[Kappa]S*(-6 + 19*\[Nu] + 
         6*\[Nu]^2) + 2*(-3 + 7*\[Nu] + 6*\[Nu]^2)*\[Chi]A^2 - 
       12*\[Delta]*(1 + 2*\[Nu])*\[Chi]A*\[Chi]S + 
       2*(-3 - 7*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2))/(8*ar^3*(-1 + er)^3*
      (1 + er)^3)
\[ScriptCapitalG]3 = (er^3*\[Epsilon]^4*(1 - 3*\[Nu])*\[Nu])/
     (32*ar^2*(-1 + er)^2*(1 + er)^2) + 
    (er^3*SO*\[Epsilon]^5*\[Nu]*(-(\[Delta]*\[Chi]A) + 
       (-1 + 2*\[Nu])*\[Chi]S))/(8*ar^2*(-1 + er)^2*(1 + er)^2*
      Sqrt[ar - ar*er^2]) + (er^3*SO^2*\[Epsilon]^6*\[Nu]*
      (\[Kappa]S + 2*\[Kappa]S*\[Nu] - 3*\[Chi]A^2 + 4*\[Nu]*\[Chi]A^2 - 
       3*\[Chi]S^2 + \[Delta]*(\[Kappa]A - 6*\[Chi]A*\[Chi]S)))/
     (16*ar^3*(-1 + er)^3*(1 + er)^3)


(* ::Text:: *)
(*\[ScriptCapitalG]2 and \[ScriptCapitalG]3 in terms of energy and angular momentum*)


\[ScriptCapitalG]2\[LetterSpace]EL = 
   ((-1 + 2*L^2*mE)*\[Epsilon]^4*(-1 - 19*\[Nu] + 3*\[Nu]^2))/(8*L^4) - 
    ((-1 + 2*L^2*mE)*SO*\[Epsilon]^5*\[Nu]*(-5*\[Delta]*\[Chi]A + 
       (-7 + 6*\[Nu])*\[Chi]S))/(4*L^5) + 
    ((-1 + 2*L^2*mE)*SO^2*\[Epsilon]^6*(\[Delta]*\[Kappa]A*(-6 + 7*\[Nu]) + 
       \[Kappa]S*(-6 + 19*\[Nu] + 6*\[Nu]^2) + 2*(-3 + 7*\[Nu] + 6*\[Nu]^2)*
        \[Chi]A^2 - 12*\[Delta]*(1 + 2*\[Nu])*\[Chi]A*\[Chi]S + 
       2*(-3 - 7*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2))/(8*L^6)
\[ScriptCapitalG]3\[LetterSpace]EL = 
   -((1 - 2*L^2*mE)^(3/2)*\[Epsilon]^4*\[Nu]*(-1 + 3*\[Nu]))/(32*L^4) - 
    ((1 - 2*L^2*mE)^(3/2)*SO*\[Epsilon]^5*\[Nu]*(\[Delta]*\[Chi]A + \[Chi]S - 
       2*\[Nu]*\[Chi]S))/(8*L^5) - ((1 - 2*L^2*mE)^(3/2)*SO^2*\[Epsilon]^6*
      \[Nu]*(\[Kappa]S + 2*\[Kappa]S*\[Nu] - 3*\[Chi]A^2 + 
       4*\[Nu]*\[Chi]A^2 - 3*\[Chi]S^2 + \[Delta]*(\[Kappa]A - 
         6*\[Chi]A*\[Chi]S)))/(16*L^6)


(* ::Section::Closed:: *)
(*relations for the mean anomaly \[ScriptL]*)


(* ::Text:: *)
(*\[ScriptL] = n (t - t0) = u - et Sin[u] + \[ScriptCapitalF]vu (v - u) + \[ScriptCapitalF]v Sin[v]*)


\[ScriptCapitalF]vu = (-3*\[Epsilon]^4*(-5 + 2*\[Nu]))/
     (2*ar^2*Sqrt[1 - er^2]) + 
    (SO*\[Epsilon]^5*(-4*\[Delta]*(-3 + \[Nu])*\[Chi]A + 
       2*(6 + (-8 + \[Nu])*\[Nu])*\[Chi]S))/(ar^(5/2)*(-1 + er^2)) - 
    (SO^2*\[Epsilon]^6*(\[Delta]*\[Kappa]A*(-14 + 5*\[Nu]) + 
       \[Kappa]S*(-14 + 33*\[Nu] - 6*\[Nu]^2) - 2*(11 - 46*\[Nu] + 6*\[Nu]^2)*
        \[Chi]A^2 + 4*\[Delta]*(-11 + 9*\[Nu])*\[Chi]A*\[Chi]S - 
       2*(11 + 4*(-4 + \[Nu])*\[Nu])*\[Chi]S^2))/(2*ar^3*(1 - er^2)^(3/2))
\[ScriptCapitalF]v = -(er*\[Epsilon]^4*(-15 + \[Nu])*\[Nu])/
     (8*ar^2*Sqrt[1 - er^2]) + (er*SO*\[Epsilon]^5*
      (\[Delta]*(4 + \[Nu])*\[Chi]A + (4 - 3*\[Nu] - 2*\[Nu]^2)*\[Chi]S))/
     (2*ar^(5/2)*(-1 + er^2)) - (er*SO^2*\[Epsilon]^6*
      (\[Delta]*\[Kappa]A*(-8 + 3*\[Nu]) + \[Kappa]S*(-8 + 19*\[Nu] - 
         2*\[Nu]^2) - 8*\[Chi]A^2 + 31*\[Nu]*\[Chi]A^2 - 
       4*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(-8 + 7*\[Nu])*\[Chi]A*\[Chi]S - 
       8*\[Chi]S^2 + 15*\[Nu]*\[Chi]S^2))/(4*ar^3*(1 - er^2)^(3/2))


(* ::Text:: *)
(*\[ScriptCapitalF]vu and \[ScriptCapitalF]v in terms of energy and angular momentum*)


\[ScriptCapitalF]vu\[LetterSpace]EL = 
   (-3*Sqrt[2]*mE^(3/2)*\[Epsilon]^4*(-5 + 2*\[Nu]))/L - 
    (2*Sqrt[2]*mE^(3/2)*SO*\[Epsilon]^5*(-4*\[Delta]*(-3 + \[Nu])*\[Chi]A + 
       2*(6 + (-8 + \[Nu])*\[Nu])*\[Chi]S))/L^2 - 
    (Sqrt[2]*mE^(3/2)*SO^2*\[Epsilon]^6*(\[Delta]*\[Kappa]A*(-14 + 5*\[Nu]) + 
       \[Kappa]S*(-14 + 33*\[Nu] - 6*\[Nu]^2) - 2*(11 - 46*\[Nu] + 6*\[Nu]^2)*
        \[Chi]A^2 + 4*\[Delta]*(-11 + 9*\[Nu])*\[Chi]A*\[Chi]S - 
       2*(11 + 4*(-4 + \[Nu])*\[Nu])*\[Chi]S^2))/L^3
\[ScriptCapitalF]v\[LetterSpace]EL = 
   -(mE^(3/2)*Sqrt[1/2 - L^2*mE]*\[Epsilon]^4*(-15 + \[Nu])*\[Nu])/(2*L) - 
    (mE^(3/2)*Sqrt[2 - 4*L^2*mE]*SO*\[Epsilon]^5*
      (\[Delta]*(4 + \[Nu])*\[Chi]A + (4 - 3*\[Nu] - 2*\[Nu]^2)*\[Chi]S))/
     L^2 - (mE^(3/2)*Sqrt[1/2 - L^2*mE]*SO^2*\[Epsilon]^6*
      (\[Delta]*\[Kappa]A*(-8 + 3*\[Nu]) + \[Kappa]S*(-8 + 19*\[Nu] - 
         2*\[Nu]^2) - 8*\[Chi]A^2 + 31*\[Nu]*\[Chi]A^2 - 
       4*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(-8 + 7*\[Nu])*\[Chi]A*\[Chi]S - 
       8*\[Chi]S^2 + 15*\[Nu]*\[Chi]S^2))/L^3


(* ::Text:: *)
(*\[ScriptL](x,et,u)*)


\[ScriptL]\[LetterSpace]xeu = u - et*Sin[u] + 
    (x^2*\[Epsilon]^4*(12*u*(-5 + 2*\[Nu]) - 24*(-5 + 2*\[Nu])*
        ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]] - 
       et*(-15 + \[Nu])*\[Nu]*Sin[2*ArcTan[Sqrt[-((1 + et)/(-1 + et))]*
            Tan[u/2]]]))/(8*Sqrt[1 - et^2]) + 
    (SO*x^(5/2)*\[Epsilon]^5*(8*u*\[Delta]*(-3 + \[Nu])*\[Chi]A - 
       4*u*(6 - 8*\[Nu] + \[Nu]^2)*\[Chi]S - 
       8*(2*\[Delta]*(-3 + \[Nu])*\[Chi]A - (6 - 8*\[Nu] + \[Nu]^2)*\[Chi]S)*
        ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]] + 
       et*(\[Delta]*(4 + \[Nu])*\[Chi]A + (4 - 3*\[Nu] - 2*\[Nu]^2)*\[Chi]S)*
        Sin[2*ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]]))/
     (2*(-1 + et^2)) + (SO^2*x^3*\[Epsilon]^6*
      (2*u*(-14*\[Kappa]S + 33*\[Kappa]S*\[Nu] - 6*\[Kappa]S*\[Nu]^2 + 
         \[Delta]*\[Kappa]A*(-14 + 5*\[Nu]) - 22*\[Chi]A^2 + 
         92*\[Nu]*\[Chi]A^2 - 12*\[Nu]^2*\[Chi]A^2 + 
         4*\[Delta]*(-11 + 9*\[Nu])*\[Chi]A*\[Chi]S - 22*\[Chi]S^2 + 
         32*\[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2) + 
       4*(\[Delta]*\[Kappa]A*(14 - 5*\[Nu]) + \[Kappa]S*(14 - 33*\[Nu] + 
           6*\[Nu]^2) + 22*\[Chi]A^2 - 92*\[Nu]*\[Chi]A^2 + 
         12*\[Nu]^2*\[Chi]A^2 + 4*\[Delta]*(11 - 9*\[Nu])*\[Chi]A*\[Chi]S + 
         22*\[Chi]S^2 - 32*\[Nu]*\[Chi]S^2 + 8*\[Nu]^2*\[Chi]S^2)*
        ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]] + 
       et*(\[Delta]*\[Kappa]A*(8 - 3*\[Nu]) + \[Kappa]S*(8 - 19*\[Nu] + 
           2*\[Nu]^2) + 8*\[Chi]A^2 - 31*\[Nu]*\[Chi]A^2 + 
         4*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(8 - 7*\[Nu])*\[Chi]A*\[Chi]S + 
         8*\[Chi]S^2 - 15*\[Nu]*\[Chi]S^2)*
        Sin[2*ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]]))/
     (4*(1 - et^2)^(3/2))


(* ::Section::Closed:: *)
(*u(x, et, \[ScriptL]) expanded to O(e^8)*)


u\[LetterSpace]xe\[ScriptL] = \[ScriptL] + et*Sin[\[ScriptL]] + 
    et^2*Cos[\[ScriptL]]*Sin[\[ScriptL]] + x^2*\[Epsilon]^4*
     ((et*(-60 + 9*\[Nu] + \[Nu]^2)*Sin[\[ScriptL]])/8 + 
      (3*et^2*(-50 + 5*\[Nu] + \[Nu]^2)*Cos[\[ScriptL]]*Sin[\[ScriptL]])/8 + 
      (et^3*(-440 + 71*\[Nu] + 7*\[Nu]^2 + (-760 + 49*\[Nu] + 17*\[Nu]^2)*
          Cos[2*\[ScriptL]])*Sin[\[ScriptL]])/32 + 
      (et^4*Cos[\[ScriptL]]*(-11*(-15 + \[Nu])*\[Nu] + 
         (-2925 + 105*\[Nu] + 71*\[Nu]^2)*Cos[2*\[ScriptL]])*Sin[\[ScriptL]])/
       48 + (et^5*(29*(-1660 + 289*\[Nu] + 25*\[Nu]^2) + 
         4*(-21520 + 3133*\[Nu] + 365*\[Nu]^2)*Cos[2*\[ScriptL]] + 
         (-100980 + 1167*\[Nu] + 2615*\[Nu]^2)*Cos[4*\[ScriptL]])*
        Sin[\[ScriptL]])/2560 + (et^7*(-29855040 + 5456586*\[Nu] + 
         432362*\[Nu]^2 + (-55153080 + 9023697*\[Nu] + 869169*\[Nu]^2)*
          Cos[2*\[ScriptL]] + 6*(-8245920 + 1476513*\[Nu] + 121457*\[Nu]^2)*
          Cos[4*\[ScriptL]] - 86873160*Cos[6*\[ScriptL]] - 
         2534241*\[Nu]*Cos[6*\[ScriptL]] + 2485567*\[Nu]^2*Cos[6*\[ScriptL]])*
        Sin[\[ScriptL]])/1290240 + (et^6*(-17600 + 185*\[Nu] + 457*\[Nu]^2 + 
         (9550 + 3920*\[Nu] - 516*\[Nu]^2)*Cos[2*\[ScriptL]] + 
         (-32900 - 325*\[Nu] + 899*\[Nu]^2)*Cos[4*\[ScriptL]])*
        Sin[2*\[ScriptL]])/640 + (et^8*(2551640 + 771604*\[Nu] - 
         119484*\[Nu]^2 + 81*(-127295 - 7612*\[Nu] + 3902*\[Nu]^2)*
          Cos[2*\[ScriptL]] + (5418840 + 1405644*\[Nu] - 238212*\[Nu]^2)*
          Cos[4*\[ScriptL]] - 9519335*Cos[6*\[ScriptL]] - 
         445576*\[Nu]*Cos[6*\[ScriptL]] + 283554*\[Nu]^2*Cos[6*\[ScriptL]])*
        Sin[2*\[ScriptL]])/107520) + 
    (et^3*(-Sin[\[ScriptL]] + 3*Sin[3*\[ScriptL]]))/8 + 
    (et^4*(-Sin[2*\[ScriptL]] + 2*Sin[4*\[ScriptL]]))/6 + 
    (et^5*(2*Sin[\[ScriptL]] - 81*Sin[3*\[ScriptL]] + 125*Sin[5*\[ScriptL]]))/
     384 + (et^6*(5*Sin[2*\[ScriptL]] - 64*Sin[4*\[ScriptL]] + 
       81*Sin[6*\[ScriptL]]))/240 + 
    (et^7*(-5*Sin[\[ScriptL]] + 2187*Sin[3*\[ScriptL]] - 
       15625*Sin[5*\[ScriptL]] + 16807*Sin[7*\[ScriptL]]))/46080 + 
    SO*x^(5/2)*\[Epsilon]^5*((et*(-7*\[Delta]*(-4 + \[Nu])*\[Chi]A + 
         (28 - 35*\[Nu] + 2*\[Nu]^2)*\[Chi]S)*Sin[\[ScriptL]])/2 - 
      (et^2*(\[Delta]*(-72 + 17*\[Nu])*\[Chi]A + (-72 + 89*\[Nu] - 4*\[Nu]^2)*
          \[Chi]S)*Cos[\[ScriptL]]*Sin[\[ScriptL]])/2 - 
      (et^4*Cos[\[ScriptL]]*(2*(\[Delta]*(-86 + 31*\[Nu])*\[Chi]A + 
           (-86 + 117*\[Nu] - 17*\[Nu]^2)*\[Chi]S) + 
         (\[Delta]*(-1454 + 319*\[Nu])*\[Chi]A + (-1454 + 1773*\[Nu] - 
             53*\[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]])*Sin[\[ScriptL]])/12 - 
      (et^3*(\[Delta]*(-780 + 197*\[Nu])*\[Chi]A + 
         (-780 + 977*\[Nu] - 58*\[Nu]^2)*\[Chi]S + 
         (\[Delta]*(-1116 + 253*\[Nu])*\[Chi]A + (-1116 + 1369*\[Nu] - 
             50*\[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]])*Sin[\[ScriptL]])/24 - 
      (et^6*Cos[\[ScriptL]]*(-29404*\[Delta]*\[Chi]A + 6544*\[Delta]*\[Nu]*
          \[Chi]A - 29404*\[Chi]S + 35948*\[Nu]*\[Chi]S - 
         1178*\[Nu]^2*\[Chi]S + 2*(\[Delta]*(1556 + 459*\[Nu])*\[Chi]A + 
           (1556 - 1097*\[Nu] - 858*\[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]] + 
         (\[Delta]*(-50268 + 10463*\[Nu])*\[Chi]A + 
           (-50268 + 60731*\[Nu] - 1186*\[Nu]^2)*\[Chi]S)*Cos[4*\[ScriptL]])*
        Sin[\[ScriptL]])/240 - (et^8*(1862880*\[Delta]*\[Chi]A - 
         65048*\[Delta]*\[Nu]*\[Chi]A + 1862880*\[Chi]S - 
         1927928*\[Nu]*\[Chi]S - 324848*\[Nu]^2*\[Chi]S + 
         3*(2*\[Delta]*(-2752515 + 564544*\[Nu])*\[Chi]A + 
           (-5505030 + 6634118*\[Nu] - 110737*\[Nu]^2)*\[Chi]S)*
          Cos[2*\[ScriptL]] - 24*(\[Delta]*(-214110 + 23917*\[Nu])*\[Chi]A + 
           (-214110 + 238027*\[Nu] + 18547*\[Nu]^2)*\[Chi]S)*
          Cos[4*\[ScriptL]] - 14825850*\[Delta]*\[Chi]A*Cos[6*\[ScriptL]] + 
         2957072*\[Delta]*\[Nu]*\[Chi]A*Cos[6*\[ScriptL]] - 
         14825850*\[Chi]S*Cos[6*\[ScriptL]] + 17782922*\[Nu]*\[Chi]S*
          Cos[6*\[ScriptL]] - 202543*\[Nu]^2*\[Chi]S*Cos[6*\[ScriptL]])*
        Sin[2*\[ScriptL]])/80640 + 
      (et^5*(-30*(\[Delta]*(-1220 + 381*\[Nu])*\[Chi]A + 
           (-1220 + 1601*\[Nu] - 174*\[Nu]^2)*\[Chi]S)*Sin[\[ScriptL]] - 
         25*(\[Delta]*(-516 + 305*\[Nu])*\[Chi]A + 
           (-516 + 821*\[Nu] - 238*\[Nu]^2)*\[Chi]S)*Sin[3*\[ScriptL]] + 
         3*(\[Delta]*(50852 - 10849*\[Nu])*\[Chi]A + 
           (50852 - 61701*\[Nu] + 1502*\[Nu]^2)*\[Chi]S)*Sin[5*\[ScriptL]]))/
       3840 + (et^7*(-35*(\[Delta]*(-176852 + 55719*\[Nu])*\[Chi]A + 
           (-176852 + 232571*\[Nu] - 25782*\[Nu]^2)*\[Chi]S)*
          Sin[\[ScriptL]] - 21*(7*\[Delta]*(-39396 + 11443*\[Nu])*\[Chi]A + 
           (-275772 + 355873*\[Nu] - 32450*\[Nu]^2)*\[Chi]S)*
          Sin[3*\[ScriptL]] - 9171652*\[Delta]*\[Chi]A*Sin[5*\[ScriptL]] + 
         495579*\[Delta]*\[Nu]*\[Chi]A*Sin[5*\[ScriptL]] - 
         9171652*\[Chi]S*Sin[5*\[ScriptL]] + 9667231*\[Nu]*\[Chi]S*
          Sin[5*\[ScriptL]] + 1398978*\[Nu]^2*\[Chi]S*Sin[5*\[ScriptL]] + 
         44691532*\[Delta]*\[Chi]A*Sin[7*\[ScriptL]] - 9097521*\[Delta]*\[Nu]*
          \[Chi]A*Sin[7*\[ScriptL]] + 44691532*\[Chi]S*Sin[7*\[ScriptL]] - 
         53789053*\[Nu]*\[Chi]S*Sin[7*\[ScriptL]] + 820410*\[Nu]^2*\[Chi]S*
          Sin[7*\[ScriptL]]))/645120) + SO^2*x^3*\[Epsilon]^6*
     ((et*(\[Delta]*\[Kappa]A*(-36 + 13*\[Nu]) + \[Kappa]S*
          (-36 + 85*\[Nu] - 14*\[Nu]^2) - 52*\[Chi]A^2 + 
         215*\[Nu]*\[Chi]A^2 - 28*\[Nu]^2*\[Chi]A^2 + 
         2*\[Delta]*(-52 + 43*\[Nu])*\[Chi]A*\[Chi]S - 52*\[Chi]S^2 + 
         79*\[Nu]*\[Chi]S^2 - 16*\[Nu]^2*\[Chi]S^2)*Sin[\[ScriptL]])/4 + 
      (et^2*(\[Delta]*\[Kappa]A*(-94 + 34*\[Nu]) + 
         \[Kappa]S*(-94 + 222*\[Nu] - 36*\[Nu]^2) - 134*\[Chi]A^2 + 
         553*\[Nu]*\[Chi]A^2 - 72*\[Nu]^2*\[Chi]A^2 + 
         2*\[Delta]*(-134 + 111*\[Nu])*\[Chi]A*\[Chi]S - 134*\[Chi]S^2 + 
         205*\[Nu]*\[Chi]S^2 - 40*\[Nu]^2*\[Chi]S^2)*Cos[\[ScriptL]]*
        Sin[\[ScriptL]])/4 + (et^4*Cos[\[ScriptL]]*
        (2*(\[Delta]*\[Kappa]A*(-476 + 171*\[Nu]) + \[Kappa]S*
            (-476 + 1123*\[Nu] - 194*\[Nu]^2) - 716*\[Chi]A^2 + 
           2977*\[Nu]*\[Chi]A^2 - 388*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(-716 + 589*\[Nu])*\[Chi]A*\[Chi]S - 716*\[Chi]S^2 + 
           1065*\[Nu]*\[Chi]S^2 - 240*\[Nu]^2*\[Chi]S^2) + 
         (-3866*\[Kappa]S + 9133*\[Kappa]S*\[Nu] - 1454*\[Kappa]S*\[Nu]^2 + 
           \[Delta]*\[Kappa]A*(-3866 + 1401*\[Nu]) - 5426*\[Chi]A^2 + 
           22342*\[Nu]*\[Chi]A^2 - 2908*\[Nu]^2*\[Chi]A^2 + 
           4*\[Delta]*(-2713 + 2252*\[Nu])*\[Chi]A*\[Chi]S - 5426*\[Chi]S^2 + 
           8370*\[Nu]*\[Chi]S^2 - 1560*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]])*
        Sin[\[ScriptL]])/48 - (et^3*(1216*\[Kappa]S + \[Delta]*\[Kappa]A*
          (1216 - 439*\[Nu]) - 2871*\[Kappa]S*\[Nu] + 474*\[Kappa]S*\[Nu]^2 + 
         1760*\[Chi]A^2 - 7279*\[Nu]*\[Chi]A^2 + 948*\[Nu]^2*\[Chi]A^2 + 
         10*\[Delta]*(352 - 291*\[Nu])*\[Chi]A*\[Chi]S + 1760*\[Chi]S^2 - 
         2671*\[Nu]*\[Chi]S^2 + 544*\[Nu]^2*\[Chi]S^2 + 
         (\[Delta]*\[Kappa]A*(1472 - 533*\[Nu]) + \[Kappa]S*
            (1472 - 3477*\[Nu] + 558*\[Nu]^2) + 2080*\[Chi]A^2 - 
           8573*\[Nu]*\[Chi]A^2 + 1116*\[Nu]^2*\[Chi]A^2 + 
           10*\[Delta]*(416 - 345*\[Nu])*\[Chi]A*\[Chi]S + 2080*\[Chi]S^2 - 
           3197*\[Nu]*\[Chi]S^2 + 608*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]])*
        Sin[\[ScriptL]])/48 + (et^5*(-60692*\[Delta]*\[Kappa]A - 
         60692*\[Kappa]S + 21905*\[Delta]*\[Kappa]A*\[Nu] + 
         143289*\[Kappa]S*\[Nu] - 23718*\[Kappa]S*\[Nu]^2 - 88036*\[Chi]A^2 + 
         364211*\[Nu]*\[Chi]A^2 - 47436*\[Nu]^2*\[Chi]A^2 - 
         176072*\[Delta]*\[Chi]A*\[Chi]S + 145518*\[Delta]*\[Nu]*\[Chi]A*
          \[Chi]S - 88036*\[Chi]S^2 + 133451*\[Nu]*\[Chi]S^2 - 
         27344*\[Nu]^2*\[Chi]S^2 + 12*(\[Delta]*\[Kappa]A*
            (-7592 + 2745*\[Nu]) + \[Kappa]S*(-7592 + 17929*\[Nu] - 
             2918*\[Nu]^2) - 10856*\[Chi]A^2 + 44821*\[Nu]*\[Chi]A^2 - 
           5836*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(-10856 + 8989*\[Nu])*\[Chi]A*
            \[Chi]S - 10856*\[Chi]S^2 + 16581*\[Nu]*\[Chi]S^2 - 
           3264*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]] + 
         (\[Delta]*\[Kappa]A*(-68044 + 24675*\[Nu]) + 
           \[Kappa]S*(-68044 + 160763*\[Nu] - 25426*\[Nu]^2) - 
           94972*\[Chi]A^2 + 390737*\[Nu]*\[Chi]A^2 - 50852*\[Nu]^2*
            \[Chi]A^2 + 2*\[Delta]*(-94972 + 78893*\[Nu])*\[Chi]A*\[Chi]S - 
           94972*\[Chi]S^2 + 146937*\[Nu]*\[Chi]S^2 - 26928*\[Nu]^2*
            \[Chi]S^2)*Cos[4*\[ScriptL]])*Sin[\[ScriptL]])/1280 + 
      (et^6*Cos[\[ScriptL]]*(-45128*\[Delta]*\[Kappa]A - 45128*\[Kappa]S + 
         16333*\[Delta]*\[Kappa]A*\[Nu] + 106589*\[Kappa]S*\[Nu] - 
         17182*\[Kappa]S*\[Nu]^2 - 64008*\[Chi]A^2 + 263961*\[Nu]*\[Chi]A^2 - 
         34364*\[Nu]^2*\[Chi]A^2 - 128016*\[Delta]*\[Chi]A*\[Chi]S + 
         106114*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S - 64008*\[Chi]S^2 + 
         98185*\[Nu]*\[Chi]S^2 - 18880*\[Nu]^2*\[Chi]S^2 + 
         (-12906*\[Kappa]S + 30403*\[Kappa]S*\[Nu] - 5714*\[Kappa]S*\[Nu]^2 + 
           \[Delta]*\[Kappa]A*(-12906 + 4591*\[Nu]) - 20866*\[Chi]A^2 + 
           87572*\[Nu]*\[Chi]A^2 - 11428*\[Nu]^2*\[Chi]A^2 + 
           4*\[Delta]*(-10433 + 8507*\[Nu])*\[Chi]A*\[Chi]S - 
           20866*\[Chi]S^2 + 29920*\[Nu]*\[Chi]S^2 - 7960*\[Nu]^2*\[Chi]S^2)*
          Cos[2*\[ScriptL]] + (\[Delta]*\[Kappa]A*(-67636 + 24541*\[Nu]) + 
           \[Kappa]S*(-67636 + 159813*\[Nu] - 25134*\[Nu]^2) - 
           93956*\[Chi]A^2 + 386287*\[Nu]*\[Chi]A^2 - 50268*\[Nu]^2*
            \[Chi]A^2 + 2*\[Delta]*(-93956 + 78099*\[Nu])*\[Chi]A*\[Chi]S - 
           93956*\[Chi]S^2 + 145735*\[Nu]*\[Chi]S^2 - 26320*\[Nu]^2*
            \[Chi]S^2)*Cos[4*\[ScriptL]])*Sin[\[ScriptL]])/480 + 
      (et^8*(-2283184*\[Delta]*\[Kappa]A - 2283184*\[Kappa]S + 
         785536*\[Delta]*\[Kappa]A*\[Nu] + 5351904*\[Kappa]S*\[Nu] - 
         1277376*\[Kappa]S*\[Nu]^2 - 4544240*\[Chi]A^2 + 
         19516696*\[Nu]*\[Chi]A^2 - 2554752*\[Nu]^2*\[Chi]A^2 - 
         9088480*\[Delta]*\[Chi]A*\[Chi]S + 7245840*\[Delta]*\[Nu]*\[Chi]A*
          \[Chi]S - 4544240*\[Chi]S^2 + 5906104*\[Nu]*\[Chi]S^2 - 
         2261056*\[Nu]^2*\[Chi]S^2 + 3*(-16125966*\[Kappa]S + 
           38099461*\[Kappa]S*\[Nu] - 6028574*\[Kappa]S*\[Nu]^2 + 
           \[Delta]*\[Kappa]A*(-16125966 + 5847529*\[Nu]) - 
           22516630*\[Chi]A^2 + 92644064*\[Nu]*\[Chi]A^2 - 
           12057148*\[Nu]^2*\[Chi]A^2 + 20*\[Delta]*(-2251663 + 
             1870351*\[Nu])*\[Chi]A*\[Chi]S - 22516630*\[Chi]S^2 + 
           34829476*\[Nu]*\[Chi]S^2 - 6390664*\[Nu]^2*\[Chi]S^2)*
          Cos[2*\[ScriptL]] + 72*(\[Delta]*\[Kappa]A*(53906 - 20290*\[Nu]) + 
           2*\[Kappa]S*(26953 - 64051*\[Nu] + 6362*\[Nu]^2) + 
           51498*\[Chi]A^2 - 197523*\[Nu]*\[Chi]A^2 + 25448*\[Nu]^2*
            \[Chi]A^2 + 2*\[Delta]*(51498 - 45437*\[Nu])*\[Chi]A*\[Chi]S + 
           51498*\[Chi]S^2 - 99343*\[Nu]*\[Chi]S^2 - 2408*\[Nu]^2*\[Chi]S^2)*
          Cos[4*\[ScriptL]] - 40264730*\[Delta]*\[Kappa]A*Cos[6*\[ScriptL]] - 
         40264730*\[Kappa]S*Cos[6*\[ScriptL]] + 14623307*\[Delta]*\[Kappa]A*
          \[Nu]*Cos[6*\[ScriptL]] + 95152767*\[Kappa]S*\[Nu]*
          Cos[6*\[ScriptL]] - 14825850*\[Kappa]S*\[Nu]^2*Cos[6*\[ScriptL]] - 
         55495666*\[Chi]A^2*Cos[6*\[ScriptL]] + 227896808*\[Nu]*\[Chi]A^2*
          Cos[6*\[ScriptL]] - 29651700*\[Nu]^2*\[Chi]A^2*Cos[6*\[ScriptL]] - 
         110991332*\[Delta]*\[Chi]A*\[Chi]S*Cos[6*\[ScriptL]] + 
         92357748*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S*Cos[6*\[ScriptL]] - 
         55495666*\[Chi]S^2*Cos[6*\[ScriptL]] + 86443604*\[Nu]*\[Chi]S^2*
          Cos[6*\[ScriptL]] - 15230936*\[Nu]^2*\[Chi]S^2*Cos[6*\[ScriptL]])*
        Sin[2*\[ScriptL]])/322560 + 
      (et^7*(35*(\[Delta]*\[Kappa]A*(-495776 + 177847*\[Nu]) + 
           \[Kappa]S*(-495776 + 1169399*\[Nu] - 204634*\[Nu]^2) - 
           753984*\[Chi]A^2 + 3139551*\[Nu]*\[Chi]A^2 - 409268*\[Nu]^2*
            \[Chi]A^2 + 2*\[Delta]*(-753984 + 619391*\[Nu])*\[Chi]A*\[Chi]S - 
           753984*\[Chi]S^2 + 1115167*\[Nu]*\[Chi]S^2 - 258208*\[Nu]^2*
            \[Chi]S^2)*Sin[\[ScriptL]] + 
         63*(\[Delta]*\[Kappa]A*(-244368 + 87827*\[Nu]) + 
           \[Kappa]S*(-244368 + 576563*\[Nu] - 99202*\[Nu]^2) - 
           366320*\[Chi]A^2 + 1522387*\[Nu]*\[Chi]A^2 - 198404*\[Nu]^2*
            \[Chi]A^2 + 10*\[Delta]*(-73264 + 60295*\[Nu])*\[Chi]A*\[Chi]S - 
           366320*\[Chi]S^2 + 545843*\[Nu]*\[Chi]S^2 - 121952*\[Nu]^2*
            \[Chi]S^2)*Sin[3*\[ScriptL]] - 2787344*\[Delta]*\[Kappa]A*
          Sin[5*\[ScriptL]] - 2787344*\[Kappa]S*Sin[5*\[ScriptL]] + 
         932785*\[Delta]*\[Kappa]A*\[Nu]*Sin[5*\[ScriptL]] + 
         6507473*\[Kappa]S*\[Nu]*Sin[5*\[ScriptL]] - 1821526*\[Kappa]S*
          \[Nu]^2*Sin[5*\[ScriptL]] - 6386352*\[Chi]A^2*Sin[5*\[ScriptL]] + 
         27783777*\[Nu]*\[Chi]A^2*Sin[5*\[ScriptL]] - 3643052*\[Nu]^2*
          \[Chi]A^2*Sin[5*\[ScriptL]] - 12772704*\[Delta]*\[Chi]A*\[Chi]S*
          Sin[5*\[ScriptL]] + 10051426*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S*
          Sin[5*\[ScriptL]] - 6386352*\[Chi]S^2*Sin[5*\[ScriptL]] + 
         7813057*\[Nu]*\[Chi]S^2*Sin[5*\[ScriptL]] - 3599008*\[Nu]^2*
          \[Chi]S^2*Sin[5*\[ScriptL]] - 60425344*\[Delta]*\[Kappa]A*
          Sin[7*\[ScriptL]] - 60425344*\[Kappa]S*Sin[7*\[ScriptL]] + 
         21935561*\[Delta]*\[Kappa]A*\[Nu]*Sin[7*\[ScriptL]] + 
         142786249*\[Kappa]S*\[Nu]*Sin[7*\[ScriptL]] - 22345766*\[Kappa]S*
          \[Nu]^2*Sin[7*\[ScriptL]] - 83591520*\[Chi]A^2*Sin[7*\[ScriptL]] + 
         343463601*\[Nu]*\[Chi]A^2*Sin[7*\[ScriptL]] - 44691532*\[Nu]^2*
          \[Chi]A^2*Sin[7*\[ScriptL]] - 167183040*\[Delta]*\[Chi]A*\[Chi]S*
          Sin[7*\[ScriptL]] + 139045730*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S*
          Sin[7*\[ScriptL]] - 83591520*\[Chi]S^2*Sin[7*\[ScriptL]] + 
         129948209*\[Nu]*\[Chi]S^2*Sin[7*\[ScriptL]] - 23166176*\[Nu]^2*
          \[Chi]S^2*Sin[7*\[ScriptL]]))/1290240) + 
    (et^8*(-7*Sin[2*\[ScriptL]] + 448*Sin[4*\[ScriptL]] - 
       2187*Sin[6*\[ScriptL]] + 2048*Sin[8*\[ScriptL]]))/5040


(* ::Section::Closed:: *)
(*r, dr/dt, d\[Phi]/dt, and \[Phi] as functions of x, et and u*)


(* ::Subsection::Closed:: *)
(*r*)


r\[LetterSpace]xeu = (1 - et*Cos[u])/x + 
    (\[Epsilon]^2*(6 + 2*et^2*(-9 + \[Nu]) - 2*\[Nu] + 
       et*(18 - 7*\[Nu] + et^2*(-6 + 7*\[Nu]))*Cos[u]))/(6*(-1 + et^2)) + 
    (x*\[Epsilon]^4*(2*(180 + 99*\[Nu] + 4*\[Nu]^2 + 
         et^2*(378 - 438*\[Nu] - 8*\[Nu]^2) + et^4*(36 + 15*\[Nu] + 
           4*\[Nu]^2)) - et*(648 - 567*\[Nu] + 35*\[Nu]^2 + 
         et^2*(468 + 150*\[Nu] - 70*\[Nu]^2) + et^4*(72 - 231*\[Nu] + 
           35*\[Nu]^2))*Cos[u] + 36*(1 - et^2)^(3/2)*(-5 + 2*\[Nu])*
        (2 + et*Cos[u])))/(72*(-1 + et^2)^2) + 
    SO^2*((x*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
         4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
           2*\[Chi]A*\[Chi]S))*(1 + et^2 - 2*et*Cos[u]))/(2*(-1 + et^2)^2) + 
      (x^2*\[Epsilon]^6*(99*\[Kappa]S + 306*et^2*\[Kappa]S + 
         45*et^4*\[Kappa]S - 252*\[Kappa]S*\[Nu] - 750*et^2*\[Kappa]S*\[Nu] - 
         105*et^4*\[Kappa]S*\[Nu] + 18*\[Kappa]S*\[Nu]^2 + 
         168*et^2*\[Kappa]S*\[Nu]^2 + 48*et^4*\[Kappa]S*\[Nu]^2 - 
         3*\[Delta]*\[Kappa]A*(-33 + 5*et^4*(-3 + \[Nu]) + 18*\[Nu] + 
           2*et^2*(-51 + 23*\[Nu])) + 139*\[Chi]A^2 + 474*et^2*\[Chi]A^2 + 
         45*et^4*\[Chi]A^2 - 574*\[Nu]*\[Chi]A^2 - 2007*et^2*\[Nu]*
          \[Chi]A^2 - 213*et^4*\[Nu]*\[Chi]A^2 + 36*\[Nu]^2*\[Chi]A^2 + 
         336*et^2*\[Nu]^2*\[Chi]A^2 + 96*et^4*\[Nu]^2*\[Chi]A^2 + 
         2*\[Delta]*(139 + et^2*(474 - 309*\[Nu]) + et^4*(45 - 33*\[Nu]) - 
           122*\[Nu])*\[Chi]A*\[Chi]S + 139*\[Chi]S^2 + 474*et^2*\[Chi]S^2 + 
         45*et^4*\[Chi]S^2 - 226*\[Nu]*\[Chi]S^2 - 507*et^2*\[Nu]*\[Chi]S^2 - 
         33*et^4*\[Nu]*\[Chi]S^2 + 40*\[Nu]^2*\[Chi]S^2 + 
         120*et^2*\[Nu]^2*\[Chi]S^2 + 
         et*(-3*\[Kappa]S*(72 - 178*\[Nu] + 20*\[Nu]^2 + 
             et^2*(78 - 191*\[Nu] + 58*\[Nu]^2)) + \[Delta]*
            (3*\[Kappa]A*(-72 + 34*\[Nu] + et^2*(-78 + 35*\[Nu])) + 
             4*(-176 + 157*\[Nu] + 3*et^2*(-51 + 25*\[Nu]))*\[Chi]A*
              \[Chi]S) - 2*((176 - 737*\[Nu] + 60*\[Nu]^2 + 3*et^2*
                (51 - 220*\[Nu] + 58*\[Nu]^2))*\[Chi]A^2 + 
             (176 - 281*\[Nu] + 80*\[Nu]^2 - 51*et^2*(-3 + 2*\[Nu]))*
              \[Chi]S^2))*Cos[u] + 3*(1 - et^2)^(3/2)*
          (\[Delta]*\[Kappa]A*(-14 + 5*\[Nu]) + \[Kappa]S*(-14 + 33*\[Nu] - 
             6*\[Nu]^2) + 4*\[Delta]*(-11 + 9*\[Nu])*\[Chi]A*\[Chi]S - 
           2*((11 - 46*\[Nu] + 6*\[Nu]^2)*\[Chi]A^2 + 
             (11 - 16*\[Nu] + 4*\[Nu]^2)*\[Chi]S^2))*(2 + et*Cos[u])))/
       (18*(1 - et^2)^3)) + 
    SO*((-2*(-(x/(-1 + et^2)))^(3/2)*\[Epsilon]^3*
        (\[Delta]*(\[Chi]A + 3*et^2*\[Chi]A) + \[Chi]S - 
         3*et^2*(-1 + \[Nu])*\[Chi]S + \[Nu]*\[Chi]S - 4*et*\[Delta]*\[Chi]A*
          Cos[u] + 2*et*(-2 + \[Nu])*\[Chi]S*Cos[u]))/(3*x) + 
      (x^(3/2)*\[Epsilon]^5*(2*(-1 + et^2)^2*(2*\[Delta]*(-3 + \[Nu])*
            \[Chi]A - (6 - 8*\[Nu] + \[Nu]^2)*\[Chi]S)*(2 + et*Cos[u]) + 
         Sqrt[1 - et^2]*(\[Delta]*(et^2*(84 - 33*\[Nu]) + 
             et^4*(12 - 8*\[Nu]) - 4*(-6 + \[Nu]))*\[Chi]A + 
           (24 - 28*\[Nu] + 4*et^4*(3 - 5*\[Nu] + 2*\[Nu]^2) + 
             et^2*(84 - 87*\[Nu] + 10*\[Nu]^2))*\[Chi]S + 
           et*(\[Delta]*(-60 + 17*\[Nu] + 4*et^2*(-15 + 7*\[Nu]))*\[Chi]A - 
             (60 - 77*\[Nu] + 4*\[Nu]^2 + 2*et^2*(30 - 29*\[Nu] + 7*\[Nu]^2))*
              \[Chi]S)*Cos[u])))/(3*(-1 + et^2)^3))


(* ::Subsection::Closed:: *)
(*dr/dt*)


rd\[LetterSpace]xeu = (et*Sqrt[x]*Sin[u])/(1 - et*Cos[u]) + 
    (et*x^(3/2)*\[Epsilon]^2*(-7*\[Nu] + et^2*(-6 + 7*\[Nu]))*Sin[u])/
     (6*(-1 + et^2)*(-1 + et*Cos[u])) + 
    (et*x^(5/2)*\[Epsilon]^4*(-864 + 405*\[Nu] - 9*et^6*(-15 + \[Nu])*\[Nu] + 
       35*\[Nu]^2 + et^4*(-468 - 285*\[Nu] + 53*\[Nu]^2) - 
       et^2*(-846 - 69*\[Nu] + 79*\[Nu]^2) - 
       6*et*(-252 + 153*\[Nu] + 16*\[Nu]^2 + et^2*(63 + 66*\[Nu] - 
           32*\[Nu]^2) + et^4*(-54 - 57*\[Nu] + 16*\[Nu]^2))*Cos[u] + 
       3*et^2*(-324 + 189*\[Nu] + 35*\[Nu]^2 + et^2*(-234 + 366*\[Nu] - 
           70*\[Nu]^2) + et^4*(72 - 231*\[Nu] + 35*\[Nu]^2))*Cos[u]^2 - 
       et^3*(-324 + 189*\[Nu] + 35*\[Nu]^2 + et^2*(-234 + 366*\[Nu] - 
           70*\[Nu]^2) + et^4*(72 - 231*\[Nu] + 35*\[Nu]^2))*Cos[u]^3 + 
       36*(1 - et^2)^(3/2)*(-5 + 2*\[Nu])*(-4 + et*Cos[u])*
        (-1 + et*Cos[u])^2)*Sin[u])/(72*(-1 + et^2)^2*(-1 + et*Cos[u])^4) + 
    SO^2*((et*x^(5/2)*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
         \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
         \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*Sin[u])/
       (2*(-1 + et^2)^2*(-1 + et*Cos[u])) + 
      (et*x^(7/2)*\[Epsilon]^6*(468*\[Kappa]S - 396*et^2*\[Kappa]S + 
         396*et^4*\[Kappa]S - 72*et^6*\[Kappa]S - 1191*\[Kappa]S*\[Nu] + 
         957*et^2*\[Kappa]S*\[Nu] - 936*et^4*\[Kappa]S*\[Nu] + 
         171*et^6*\[Kappa]S*\[Nu] + 294*\[Kappa]S*\[Nu]^2 - 
         186*et^2*\[Kappa]S*\[Nu]^2 + 144*et^4*\[Kappa]S*\[Nu]^2 - 
         18*et^6*\[Kappa]S*\[Nu]^2 + 3*\[Delta]*\[Kappa]A*
          (156 + et^4*(132 - 48*\[Nu]) - 85*\[Nu] + 3*et^6*(-8 + 3*\[Nu]) + 
           11*et^2*(-12 + 5*\[Nu])) + 652*\[Chi]A^2 - 612*et^2*\[Chi]A^2 + 
         540*et^4*\[Chi]A^2 - 72*et^6*\[Chi]A^2 - 2809*\[Nu]*\[Chi]A^2 + 
         2550*et^2*\[Nu]*\[Chi]A^2 - 2214*et^4*\[Nu]*\[Chi]A^2 + 
         279*et^6*\[Nu]*\[Chi]A^2 + 588*\[Nu]^2*\[Chi]A^2 - 
         372*et^2*\[Nu]^2*\[Chi]A^2 + 288*et^4*\[Nu]^2*\[Chi]A^2 - 
         36*et^6*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(652 - 653*\[Nu] - 
           90*et^4*(-6 + 5*\[Nu]) + 9*et^6*(-8 + 7*\[Nu]) + 
           12*et^2*(-51 + 52*\[Nu]))*\[Chi]A*\[Chi]S + 652*\[Chi]S^2 - 
         612*et^2*\[Chi]S^2 + 540*et^4*\[Chi]S^2 - 72*et^6*\[Chi]S^2 - 
         1105*\[Nu]*\[Chi]S^2 + 1146*et^2*\[Nu]*\[Chi]S^2 - 
         846*et^4*\[Nu]*\[Chi]S^2 + 135*et^6*\[Nu]*\[Chi]S^2 + 
         280*\[Nu]^2*\[Chi]S^2 - 288*et^2*\[Nu]^2*\[Chi]S^2 + 
         144*et^4*\[Nu]^2*\[Chi]S^2 + 
         3*et*(-3*\[Kappa]S*(92 - 246*\[Nu] + 72*\[Nu]^2 - 
             4*et^2*(-5 + 10*\[Nu] + \[Nu]^2) + et^4*(20 - 47*\[Nu] + 10*
                \[Nu]^2)) - 364*\[Chi]A^2 - 36*et^2*\[Chi]A^2 - 
           108*et^4*\[Chi]A^2 + 1612*\[Nu]*\[Chi]A^2 + 123*et^2*\[Nu]*
            \[Chi]A^2 + 459*et^4*\[Nu]*\[Chi]A^2 - 432*\[Nu]^2*\[Chi]A^2 + 
           24*et^2*\[Nu]^2*\[Chi]A^2 - 60*et^4*\[Nu]^2*\[Chi]A^2 - 
           364*\[Chi]S^2 - 36*et^2*\[Chi]S^2 - 108*et^4*\[Chi]S^2 + 
           676*\[Nu]*\[Chi]S^2 - 153*et^2*\[Nu]*\[Chi]S^2 + 
           147*et^4*\[Nu]*\[Chi]S^2 - 184*\[Nu]^2*\[Chi]S^2 + 
           96*et^2*\[Nu]^2*\[Chi]S^2 - 48*et^4*\[Nu]^2*\[Chi]S^2 + 
           \[Delta]*(3*\[Kappa]A*(-92 - 20*et^2 + 62*\[Nu] + et^4*
                (-20 + 7*\[Nu])) + 2*(52*(-7 + 8*\[Nu]) + 3*et^4*
                (-36 + 29*\[Nu]) - 3*et^2*(12 + 29*\[Nu]))*\[Chi]A*\[Chi]S))*
          Cos[u] - 3*et^2*(-3*\[Kappa]S*(72 - 199*\[Nu] + 62*\[Nu]^2 + 
             2*et^2*(30 - 67*\[Nu] + 8*\[Nu]^2)) - 256*\[Chi]A^2 - 
           252*et^2*\[Chi]A^2 + 1153*\[Nu]*\[Chi]A^2 + 1041*et^2*\[Nu]*
            \[Chi]A^2 - 372*\[Nu]^2*\[Chi]A^2 - 96*et^2*\[Nu]^2*\[Chi]A^2 - 
           256*\[Chi]S^2 - 252*et^2*\[Chi]S^2 + 529*\[Nu]*\[Chi]S^2 + 
           141*et^2*\[Nu]*\[Chi]S^2 - 136*\[Nu]^2*\[Chi]S^2 + 
           \[Delta]*(3*\[Kappa]A*(-72 + 55*\[Nu] + 2*et^2*(-30 + 7*\[Nu])) + 
             2*(-256 + 329*\[Nu] + 3*et^2*(-84 + 29*\[Nu]))*\[Chi]A*\[Chi]S))*
          Cos[u]^2 + et^3*(-3*\[Kappa]S*(72 - 199*\[Nu] + 62*\[Nu]^2 + 
             2*et^2*(30 - 67*\[Nu] + 8*\[Nu]^2)) - 256*\[Chi]A^2 - 
           252*et^2*\[Chi]A^2 + 1153*\[Nu]*\[Chi]A^2 + 1041*et^2*\[Nu]*
            \[Chi]A^2 - 372*\[Nu]^2*\[Chi]A^2 - 96*et^2*\[Nu]^2*\[Chi]A^2 - 
           256*\[Chi]S^2 - 252*et^2*\[Chi]S^2 + 529*\[Nu]*\[Chi]S^2 + 
           141*et^2*\[Nu]*\[Chi]S^2 - 136*\[Nu]^2*\[Chi]S^2 + 
           \[Delta]*(3*\[Kappa]A*(-72 + 55*\[Nu] + 2*et^2*(-30 + 7*\[Nu])) + 
             2*(-256 + 329*\[Nu] + 3*et^2*(-84 + 29*\[Nu]))*\[Chi]A*\[Chi]S))*
          Cos[u]^3 - 6*(1 - et^2)^(3/2)*(\[Delta]*\[Kappa]A*(-14 + 5*\[Nu]) + 
           \[Kappa]S*(-14 + 33*\[Nu] - 6*\[Nu]^2) + 4*\[Delta]*
            (-11 + 9*\[Nu])*\[Chi]A*\[Chi]S - 
           2*((11 - 46*\[Nu] + 6*\[Nu]^2)*\[Chi]A^2 + 
             (11 - 16*\[Nu] + 4*\[Nu]^2)*\[Chi]S^2))*(-4 + et*Cos[u])*
          (-1 + et*Cos[u])^2)*Sin[u])/(36*(-1 + et^2)^3*
        (-1 + et*Cos[u])^4)) + 
    SO*((et*x^(3/2)*\[Epsilon]^3*(3*Sqrt[x - et^2*x]*(\[Delta]*\[Chi]A + 
           \[Chi]S - \[Nu]*\[Chi]S) + Sqrt[-(x/(-1 + et^2))]*
          (\[Delta]*(\[Chi]A + 3*et^2*\[Chi]A) + (1 - 3*et^2*(-1 + \[Nu]) + 
             \[Nu])*\[Chi]S))*Sin[u])/(3*(-1 + et^2)*(-1 + et*Cos[u])) + 
      (et*x^3*\[Epsilon]^5*(16*(-1 + et^2)^2*(2*\[Delta]*(-3 + \[Nu])*
            \[Chi]A - (6 - 8*\[Nu] + \[Nu]^2)*\[Chi]S)*(-4 + et*Cos[u])*
          (-1 + et*Cos[u])^2 + 3*Sqrt[1 - et^2]*(-176*\[Delta]*\[Chi]A + 
           40*et^2*\[Delta]*\[Chi]A - 200*et^4*\[Delta]*\[Chi]A + 
           16*et^6*\[Delta]*\[Chi]A + 92*\[Delta]*\[Nu]*\[Chi]A + 
           30*et^2*\[Delta]*\[Nu]*\[Chi]A + 24*et^4*\[Delta]*\[Nu]*\[Chi]A + 
           4*et^6*\[Delta]*\[Nu]*\[Chi]A - 176*\[Chi]S + 40*et^2*\[Chi]S - 
           200*et^4*\[Chi]S + 16*et^6*\[Chi]S + 268*\[Nu]*\[Chi]S - 
           34*et^2*\[Nu]*\[Chi]S + 188*et^4*\[Nu]*\[Chi]S - 
           12*et^6*\[Nu]*\[Chi]S - 40*\[Nu]^2*\[Chi]S - 12*et^2*\[Nu]^2*
            \[Chi]S - 8*et^6*\[Nu]^2*\[Chi]S + 
           et*(\[Delta]*(320 - 216*\[Nu] - 4*et^4*(-29 + 9*\[Nu]) + et^2*
                (44 + 27*\[Nu]))*\[Chi]A + (et^2*(44 + 55*\[Nu] - 
                 30*\[Nu]^2) + 8*(40 - 67*\[Nu] + 12*\[Nu]^2) + 2*et^4*
                (58 - 67*\[Nu] + 12*\[Nu]^2))*\[Chi]S)*Cos[u] - 
           6*et^2*(\[Delta]*(20 + 12*et^2 - 15*\[Nu])*\[Chi]A + 
             (20 - 6*et^2*(-2 + \[Nu]) - 35*\[Nu] + 6*\[Nu]^2)*\[Chi]S)*
            Cos[2*u] + 20*et^3*\[Delta]*\[Chi]A*Cos[3*u] + 
           12*et^5*\[Delta]*\[Chi]A*Cos[3*u] - 15*et^3*\[Delta]*\[Nu]*\[Chi]A*
            Cos[3*u] + 20*et^3*\[Chi]S*Cos[3*u] + 12*et^5*\[Chi]S*Cos[3*u] - 
           35*et^3*\[Nu]*\[Chi]S*Cos[3*u] - 6*et^5*\[Nu]*\[Chi]S*Cos[3*u] + 
           6*et^3*\[Nu]^2*\[Chi]S*Cos[3*u]))*Sin[u])/(24*(-1 + et^2)^3*
        (-1 + et*Cos[u])^4))


(* ::Subsection::Closed:: *)
(*d\[Phi]/dt*)


\[Phi]d\[LetterSpace]xeu = 
   -((et*x^(5/2)*\[Epsilon]^2*(-4 + \[Nu])*(et - Cos[u]))/
      (Sqrt[1 - et^2]*(-1 + et*Cos[u])^3)) + Sqrt[-((-1 + et^2)*x^3)]/
     (-1 + et*Cos[u])^2 + (x^(7/2)*\[Epsilon]^4*
      (-18*(-1 + et^2)*(-5 + 2*\[Nu])*(-1 + et*Cos[u])^2*
        (1 - 2*et^2 + et*Cos[u]) + Sqrt[1 - et^2]*(90 - 36*\[Nu] - 
         6*et^6*\[Nu]*(3 + 2*\[Nu]) + et^2*(75 + 50*\[Nu] - 2*\[Nu]^2) + 
         et^4*(-60 - 26*\[Nu] + 20*\[Nu]^2) - et*(246 - 67*\[Nu] + \[Nu]^2 + 
           et^4*(-12 - 97*\[Nu] + \[Nu]^2) + et^2*(81 + 74*\[Nu] + 
             16*\[Nu]^2))*Cos[u] + et^2*(114 - 35*\[Nu] + 5*\[Nu]^2 + 
           et^2*(153 - 38*\[Nu] - 4*\[Nu]^2) + et^4*(48 - 17*\[Nu] + 
             17*\[Nu]^2))*Cos[u]^2 - et^3*(-2*(21 + 11*\[Nu] + 4*\[Nu]^2) + 
           et^2*(147 - 8*\[Nu] + 14*\[Nu]^2))*Cos[u]^3)))/
     (12*(-1 + et^2)^2*(-1 + et*Cos[u])^5) + 
    SO*((2*et*x^3*\[Epsilon]^3*(\[Delta]*\[Chi]A + \[Chi]S)*(et - Cos[u]))/
       ((-1 + et^2)*(-1 + et*Cos[u])^3) + 
      (x^4*\[Epsilon]^5*(-(\[Delta]*(et^4*(36 - 54*\[Nu]) - 24*(-3 + \[Nu]) + 
            6*et^6*(-2 + \[Nu]) + et^2*(56 + 31*\[Nu]))*\[Chi]A) + 
         (-12*(6 - 8*\[Nu] + \[Nu]^2) + 6*et^6*(2 - 3*\[Nu] + 2*\[Nu]^2) - 
           6*et^4*(6 - 27*\[Nu] + 8*\[Nu]^2) + et^2*(-56 - 57*\[Nu] + 
             50*\[Nu]^2))*\[Chi]S - et*(\[Delta]*(-224 + 44*\[Nu] + 
             6*et^4*(-6 + 11*\[Nu]) + et^2*(-196 + 13*\[Nu]))*\[Chi]A + 
           (-4*(56 - 60*\[Nu] + \[Nu]^2) - 6*et^4*(6 - 23*\[Nu] + 4*
                \[Nu]^2) + et^2*(-196 + 171*\[Nu] + 34*\[Nu]^2))*\[Chi]S)*
          Cos[u] + et^2*(\[Delta]*(8*(-29 + 5*\[Nu]) + 
             6*et^4*(-6 + 5*\[Nu]) + et^2*(-188 + 53*\[Nu]))*\[Chi]A + 
           (-4*(58 - 81*\[Nu] + 2*\[Nu]^2) - 6*et^4*(6 - 5*\[Nu] + 2*
                \[Nu]^2) + et^2*(-188 + 195*\[Nu] + 26*\[Nu]^2))*\[Chi]S)*
          Cos[u]^2 - et^3*(\[Delta]*(-68 + 14*\[Nu] + 3*et^2*(-28 + 9*\[Nu]))*
            \[Chi]A + (-68 + 162*\[Nu] - 4*\[Nu]^2 + 3*et^2*(-28 + 7*\[Nu] + 
               2*\[Nu]^2))*\[Chi]S)*Cos[u]^3 - 12*Sqrt[1 - et^2]*
          (2*\[Delta]*(-3 + \[Nu])*\[Chi]A - (6 - 8*\[Nu] + \[Nu]^2)*\[Chi]S)*
          (-1 + et*Cos[u])^2*(1 - 2*et^2 + et*Cos[u])))/
       (6*(-1 + et^2)^2*(-1 + et*Cos[u])^5)) + 
    SO^2*((et*x^(7/2)*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
         \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
         \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*(et - Cos[u]))/
       ((1 - et^2)^(3/2)*(-1 + et*Cos[u])^3) + 
      (x^4*\[Epsilon]^6*(6*(-1 + et^2)*Sqrt[x]*(\[Delta]*\[Kappa]A*
            (-14 + 5*\[Nu]) + \[Kappa]S*(-14 + 33*\[Nu] - 6*\[Nu]^2) + 
           4*\[Delta]*(-11 + 9*\[Nu])*\[Chi]A*\[Chi]S - 
           2*((11 - 46*\[Nu] + 6*\[Nu]^2)*\[Chi]A^2 + 
             (11 - 16*\[Nu] + 4*\[Nu]^2)*\[Chi]S^2))*(-1 + et*Cos[u])^2*
          (1 - 2*et^2 + et*Cos[u]) - 24*Sqrt[x - et^2*x]*(-1 + et*Cos[u])^2*
          (-((1 + 4*et^2 + 3*et^4)*(-1 + 4*\[Nu])*\[Chi]A^2) + 
           2*\[Delta]*(1 + et^4*(3 - 6*\[Nu]) - et^2*(-4 + \[Nu]) + \[Nu])*
            \[Chi]A*\[Chi]S + ((1 + \[Nu])^2 + 3*et^4*(1 - 4*\[Nu] + 3*
                \[Nu]^2) - 2*et^2*(-2 + \[Nu] + 3*\[Nu]^2))*\[Chi]S^2 + 
           et*(-1 + 3*et^2)*((1 + 3*et^2)*(-1 + 4*\[Nu])*\[Chi]A^2 + 
             2*\[Delta]*(-1 + 3*et^2*(-1 + \[Nu]))*\[Chi]A*\[Chi]S - 
             (-1 + 3*et^2*(-1 + \[Nu]) - \[Nu])*(-1 + \[Nu])*\[Chi]S^2)*
            Cos[u]) + Sqrt[x - et^2*x]*(-84*\[Kappa]S - 114*et^2*\[Kappa]S - 
           72*et^4*\[Kappa]S + 24*et^6*\[Kappa]S + 198*\[Kappa]S*\[Nu] + 
           311*et^2*\[Kappa]S*\[Nu] + 160*et^4*\[Kappa]S*\[Nu] - 
           48*et^6*\[Kappa]S*\[Nu] - 36*\[Kappa]S*\[Nu]^2 + 
           38*et^2*\[Kappa]S*\[Nu]^2 - 128*et^4*\[Kappa]S*\[Nu]^2 + 
           24*et^6*\[Kappa]S*\[Nu]^2 + \[Delta]*\[Kappa]A*(-84 + 24*et^6 + 
             30*\[Nu] + 8*et^4*(-9 + 2*\[Nu]) + et^2*(-114 + 83*\[Nu])) - 
           108*\[Chi]A^2 - 138*et^2*\[Chi]A^2 + 24*et^6*\[Chi]A^2 + 
           456*\[Nu]*\[Chi]A^2 + 512*et^2*\[Nu]*\[Chi]A^2 + 
           112*et^4*\[Nu]*\[Chi]A^2 - 120*et^6*\[Nu]*\[Chi]A^2 - 
           72*\[Nu]^2*\[Chi]A^2 + 76*et^2*\[Nu]^2*\[Chi]A^2 - 
           256*et^4*\[Nu]^2*\[Chi]A^2 + 48*et^6*\[Nu]^2*\[Chi]A^2 + 
           4*\[Delta]*(-54 + et^6*(12 - 24*\[Nu]) + et^2*(-69 + \[Nu]) + 
             66*\[Nu] + 56*et^4*\[Nu])*\[Chi]A*\[Chi]S - 108*\[Chi]S^2 - 
           138*et^2*\[Chi]S^2 + 24*et^6*\[Chi]S^2 + 240*\[Nu]*\[Chi]S^2 + 
           44*et^2*\[Nu]*\[Chi]S^2 + 112*et^4*\[Nu]*\[Chi]S^2 - 
           72*et^6*\[Nu]*\[Chi]S^2 - 24*\[Nu]^2*\[Chi]S^2 - 
           96*et^2*\[Nu]^2*\[Chi]S^2 + 72*et^4*\[Nu]^2*\[Chi]S^2 - 
           et*(-(\[Kappa]S*(300 - 740*\[Nu] + 64*\[Nu]^2 + 2*et^4*
                 (36 - 97*\[Nu] + 74*\[Nu]^2) + et^2*(366 - 929*\[Nu] + 
                  94*\[Nu]^2))) + \[Delta]*(\[Kappa]A*(20*(-15 + 7*\[Nu]) + 
                 et^4*(-72 + 50*\[Nu]) + et^2*(-366 + 197*\[Nu])) + 4*
                (-258 + et^4*(144 - 167*\[Nu]) + 211*\[Nu] + 
                 et^2*(-219 + 253*\[Nu]))*\[Chi]A*\[Chi]S) - 
             2*((258 - 1057*\[Nu] + 64*\[Nu]^2 + et^2*(219 - 910*\[Nu] + 
                   94*\[Nu]^2) + et^4*(-144 + 527*\[Nu] + 148*\[Nu]^2))*
                \[Chi]A^2 + (258 - 397*\[Nu] + 36*\[Nu]^2 + et^4*(-144 + 
                   383*\[Nu] - 252*\[Nu]^2) + et^2*(219 - 472*\[Nu] + 
                   288*\[Nu]^2))*\[Chi]S^2))*Cos[u] + 
           et^2*(-(\[Kappa]S*(348 - 868*\[Nu] + 56*\[Nu]^2 + 2*et^4*
                 (36 - 89*\[Nu] + 46*\[Nu]^2) + et^2*(318 - 817*\[Nu] + 
                  158*\[Nu]^2))) + \[Delta]*(\[Kappa]A*
                (et^4*(-72 + 34*\[Nu]) + 4*(-87 + 43*\[Nu]) + 
                 et^2*(-318 + 181*\[Nu])) + 4*(-354 + et^4*(216 - 
                   259*\[Nu]) + 323*\[Nu] + et^2*(-195 + 233*\[Nu]))*\[Chi]A*
                \[Chi]S) - 2*((354 - 1457*\[Nu] + 56*\[Nu]^2 + 
                 et^4*(-216 + 835*\[Nu] + 92*\[Nu]^2) + et^2*(195 - 
                   818*\[Nu] + 158*\[Nu]^2))*\[Chi]A^2 + (354 - 605*\[Nu] + 
                 108*\[Nu]^2 + et^4*(-216 + 547*\[Nu] - 324*\[Nu]^2) + 
                 et^2*(195 - 428*\[Nu] + 288*\[Nu]^2))*\[Chi]S^2))*Cos[u]^2 + 
           et^3*(\[Kappa]S*(108 - 278*\[Nu] + 4*\[Nu]^2 + et^2*
                (138 - 343*\[Nu] + 98*\[Nu]^2)) + \[Delta]*(\[Kappa]A*
                (108 + et^2*(138 - 67*\[Nu]) - 62*\[Nu]) + 4*(138 + 
                 et^2*(81 - 53*\[Nu]) + 108*et^4*(-1 + \[Nu]) - 154*\[Nu])*
                \[Chi]A*\[Chi]S) + 2*((138 - 568*\[Nu] + 4*\[Nu]^2 + 
                 108*et^4*(-1 + 4*\[Nu]) + et^2*(81 - 344*\[Nu] + 
                   98*\[Nu]^2))*\[Chi]A^2 + (138 - 108*et^4*(-1 + \[Nu])^2 - 
                 292*\[Nu] + 84*\[Nu]^2 + et^2*(81 - 86*\[Nu] + 48*\[Nu]^2))*
                \[Chi]S^2))*Cos[u]^3)))/(12*(-1 + et^2)^3*(-1 + et*Cos[u])^5))


(* ::Subsection::Closed:: *)
(*\[Phi]*)


\[Phi]\[LetterSpace]xeu = 2*ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]] - 
    (x*\[Epsilon]^2*(6*ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]*
        (-1 + et*Cos[u]) + et*Sqrt[1 - et^2]*(-4 + \[Nu])*Sin[u]))/
     ((-1 + et^2)*(-1 + et*Cos[u])) + 
    (x^2*\[Epsilon]^4*(-2*et*Sqrt[1 - et^2]*
        (et^2*(384 - 275*\[Nu] - 7*\[Nu]^2) + 4*(-600 + 89*\[Nu] + \[Nu]^2) + 
         et*(1632 + 28*\[Nu] - 52*\[Nu]^2 + et^2*(384 - 109*\[Nu] + 
             55*\[Nu]^2))*Cos[u])*Sin[u] - 3*(-1 + et*Cos[u])*
        (32*(-54 + 28*\[Nu] + et^2*(-51 + 26*\[Nu]))*
          ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]*(-1 + et*Cos[u]) + 
         et*(96*(-1 + et^2)*(-5 + 2*\[Nu])*Sin[u] + 2*et*(-1 + et*Cos[u])*
            (8*(-1 - 19*\[Nu] + 3*\[Nu]^2)*Cos[2*ArcTan[
                 Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]] + 
             et*\[Nu]*(-1 + 3*\[Nu])*(1 + 2*Cos[4*ArcTan[
                   Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]]))*
            Sin[2*ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]]))))/
     (192*(-1 + et^2)^2*(-1 + et*Cos[u])^2) + 
    SO^2*((x^2*\[Epsilon]^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
         4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
           2*\[Chi]A*\[Chi]S))*(3*ArcTan[Sqrt[-((1 + et)/(-1 + et))]*
            Tan[u/2]]*(-1 + et*Cos[u]) + (-1 + et)*et*Sqrt[(1 + et)/(1 - et)]*
          Sin[u]))/((-1 + et^2)^2*(-1 + et*Cos[u])) + 
      (x^3*\[Epsilon]^6*(-2*et*Sqrt[1 - et^2]*(1392*\[Kappa]S + 
           96*et^2*\[Kappa]S + \[Delta]*\[Kappa]A*(1392 + 
             et^2*(96 - 115*\[Nu]) - 548*\[Nu]) - 3332*\[Kappa]S*\[Nu] - 
           307*et^2*\[Kappa]S*\[Nu] + 280*\[Kappa]S*\[Nu]^2 + 
           218*et^2*\[Kappa]S*\[Nu]^2 + 2544*\[Chi]A^2 - 
           10340*\[Nu]*\[Chi]A^2 - 103*et^2*\[Nu]*\[Chi]A^2 + 
           560*\[Nu]^2*\[Chi]A^2 + 436*et^2*\[Nu]^2*\[Chi]A^2 - 
           2*\[Delta]*(-2544 + (1316 + 103*et^2)*\[Nu])*\[Chi]A*\[Chi]S + 
           2544*\[Chi]S^2 - 2468*\[Nu]*\[Chi]S^2 - 103*et^2*\[Nu]*\[Chi]S^2 + 
           384*\[Nu]^2*\[Chi]S^2 + et*(-(\[Kappa]S*(4*(300 - 725*\[Nu] + 
                  46*\[Nu]^2) + et^2*(288 - 739*\[Nu] + 314*\[Nu]^2))) - 
             2256*\[Chi]A^2 - 288*et^2*\[Chi]A^2 + 9140*\[Nu]*\[Chi]A^2 + 
             1303*et^2*\[Nu]*\[Chi]A^2 - 368*\[Nu]^2*\[Chi]A^2 - 
             628*et^2*\[Nu]^2*\[Chi]A^2 - 2256*\[Chi]S^2 - 
             288*et^2*\[Chi]S^2 + 2420*\[Nu]*\[Chi]S^2 + 151*et^2*\[Nu]*
              \[Chi]S^2 - 384*\[Nu]^2*\[Chi]S^2 + \[Delta]*(100*\[Kappa]A*
                (-12 + 5*\[Nu]) + et^2*\[Kappa]A*(-288 + 163*\[Nu]) + 2*et^2*
                (-288 + 151*\[Nu])*\[Chi]A*\[Chi]S + 8*(-564 + 317*\[Nu])*
                \[Chi]A*\[Chi]S))*Cos[u])*Sin[u] + 3*(-1 + et*Cos[u])*
          (16*(-(\[Kappa]S*(78 - 190*\[Nu] + 20*\[Nu]^2 + et^2*
                 (78 - 191*\[Nu] + 58*\[Nu]^2))) + \[Delta]*(\[Kappa]A*
                (-78 + 34*\[Nu] + et^2*(-78 + 35*\[Nu])) + 4*(-67 + 
                 55*\[Nu] + et^2*(-51 + 25*\[Nu]))*\[Chi]A*\[Chi]S) - 
             2*((67 - 279*\[Nu] + 20*\[Nu]^2 + et^2*(51 - 220*\[Nu] + 
                   58*\[Nu]^2))*\[Chi]A^2 + (67 + et^2*(51 - 34*\[Nu]) - 
                 99*\[Nu] + 28*\[Nu]^2)*\[Chi]S^2))*
            ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]*(-1 + et*Cos[u]) + 
           et*(16*(-1 + et^2)*(-14*\[Kappa]S + 33*\[Kappa]S*\[Nu] - 6*
                \[Kappa]S*\[Nu]^2 + \[Delta]*\[Kappa]A*(-14 + 5*\[Nu]) - 22*
                \[Chi]A^2 + 92*\[Nu]*\[Chi]A^2 - 12*\[Nu]^2*\[Chi]A^2 + 4*
                \[Delta]*(-11 + 9*\[Nu])*\[Chi]A*\[Chi]S - 22*\[Chi]S^2 + 32*
                \[Nu]*\[Chi]S^2 - 8*\[Nu]^2*\[Chi]S^2)*Sin[u] + 
             2*et*(-1 + et*Cos[u])*(4*(\[Delta]*\[Kappa]A*(-6 + 7*\[Nu]) + 
                 \[Kappa]S*(-6 + 19*\[Nu] + 6*\[Nu]^2) + 2*(-3 + 7*\[Nu] + 
                   6*\[Nu]^2)*\[Chi]A^2 - 12*\[Delta]*(1 + 2*\[Nu])*\[Chi]A*
                  \[Chi]S + 2*(-3 - 7*\[Nu] + 8*\[Nu]^2)*\[Chi]S^2)*
                Cos[2*ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]] + et*
                \[Nu]*(\[Kappa]S + 2*\[Kappa]S*\[Nu] - 3*\[Chi]A^2 + 
                 4*\[Nu]*\[Chi]A^2 - 3*\[Chi]S^2 + \[Delta]*(\[Kappa]A - 
                   6*\[Chi]A*\[Chi]S))*(1 + 2*Cos[4*ArcTan[
                     Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]]))*
              Sin[2*ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]]))))/
       (96*(-1 + et^2)^3*(-1 + et*Cos[u])^2)) + 
    SO*(-2*(-(x/(-1 + et^2)))^(3/2)*\[Epsilon]^3*
       ((4*\[Delta]*\[Chi]A - 2*(-2 + \[Nu])*\[Chi]S)*
         ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]] - 
        (et*Sqrt[1 - et^2]*(\[Delta]*\[Chi]A + \[Chi]S)*Sin[u])/
         (-1 + et*Cos[u])) + \[Epsilon]^5*
       ((et*x^(5/2)*(3 + et^2*(-9 + \[Nu]) - \[Nu])*(\[Delta]*\[Chi]A + 
           \[Chi]S)*Sin[u])/((-1 + et^2)^2*(-1 + et*Cos[u])) - 
        (et*x^(5/2)*(\[Delta]*(98 - 13*\[Nu])*\[Chi]A + (98 - 69*\[Nu])*
            \[Chi]S + et*(\[Delta]*(-82 + 9*\[Nu])*\[Chi]A + 
             (-82 + 65*\[Nu])*\[Chi]S)*Cos[u])*Sin[u])/
         (2*(-1 + et)^2*(1 + et)^2*(-1 + et*Cos[u])^2) + 
        (et*Sqrt[1 - et^2]*x^(5/2)*((16*(-2*\[Delta]*(-3 + \[Nu])*\[Chi]A + 
              (6 - 8*\[Nu] + \[Nu]^2)*\[Chi]S)*(-1 + et*Cos[u]))/
            (-1 + et^2)^2 + (8*\[Nu]^2*\[Chi]S + et^2*(3*\[Delta]*
                (-40 + 13*\[Nu])*\[Chi]A + (-120 + 103*\[Nu] - 14*\[Nu]^2)*
                \[Chi]S) + (-8*et*\[Nu]^2*\[Chi]S + et^3*
                (\[Delta]*(184 - 55*\[Nu])*\[Chi]A + (184 - 119*\[Nu] + 
                   14*\[Nu]^2)*\[Chi]S))*Cos[u])/(1 - et^2)^(5/2))*Sin[u])/
         (8*(-1 + et*Cos[u])^2) + ((-(x/(-1 + et^2)))^(5/2)*
          ((16*et*Sqrt[1 - et^2]*(-4 + \[Nu])*(\[Delta]*(\[Chi]A + 
                3*et^2*\[Chi]A) + (1 - 3*et^2*(-1 + \[Nu]) + \[Nu])*\[Chi]S)*
             Sin[u])/(-1 + et*Cos[u]) + 
           3*(8*(\[Delta]*(17*(-4 + \[Nu]) + 4*et^2*(-15 + 7*\[Nu]))*
                \[Chi]A - (68 - 81*\[Nu] + 4*\[Nu]^2 + 2*et^2*(30 - 
                   29*\[Nu] + 7*\[Nu]^2))*\[Chi]S)*ArcTan[Sqrt[
                 -((1 + et)/(-1 + et))]*Tan[u/2]] + et^2*\[Nu]*
              (-2*(5*\[Delta]*\[Chi]A + (7 - 6*\[Nu])*\[Chi]S)*
                Sin[4*ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]] - et*
                (\[Delta]*\[Chi]A + \[Chi]S - 2*\[Nu]*\[Chi]S)*
                Sin[6*ArcTan[Sqrt[-((1 + et)/(-1 + et))]*Tan[u/2]]]))))/24))


(* ::Section::Closed:: *)
(*r, dr/dt, d\[Phi]/dt, and \[Phi] as functions of x, et and \[ScriptL],  expanded to O(e^8)*)


(* ::Subsection::Closed:: *)
(*r*)


r\[LetterSpace]xel = (1 - et*Cos[\[ScriptL]] + et^2*Sin[\[ScriptL]]^2 + 
      (3*et^3*Cos[\[ScriptL]]*Sin[\[ScriptL]]^2)/2 + 
      (2*et^4*(1 + 2*Cos[2*\[ScriptL]])*Sin[\[ScriptL]]^2)/3 + 
      (5*et^5*(23*Cos[\[ScriptL]] + 25*Cos[3*\[ScriptL]])*Sin[\[ScriptL]]^2)/
       96 + (et^6*(11 + 22*Cos[2*\[ScriptL]] + 27*Cos[4*\[ScriptL]])*
        Sin[\[ScriptL]]^2)/20 + (7*et^7*(1682*Cos[\[ScriptL]] + 
         1677*Cos[3*\[ScriptL]] + 2401*Cos[5*\[ScriptL]])*Sin[\[ScriptL]]^2)/
       11520 + (et^8*(151 + 302*Cos[2*\[ScriptL]] + 295*Cos[4*\[ScriptL]] + 
         512*Cos[6*\[ScriptL]])*Sin[\[ScriptL]]^2)/315)/x + 
    \[Epsilon]^2*((-3 + \[Nu])/3 + (et*(-18 + 7*\[Nu])*Cos[\[ScriptL]])/6 + 
      (et^4*(54 - 7*\[Nu]*Cos[2*\[ScriptL]] + (-18 + 7*\[Nu])*
          Cos[4*\[ScriptL]]))/18 + (et^5*Cos[\[ScriptL]]*
        (-3006 + 945*\[Nu] - 4*(-738 + 455*\[Nu])*Cos[2*\[ScriptL]] + 
         125*(-18 + 7*\[Nu])*Cos[4*\[ScriptL]]))/1152 + 
      (et^7*Cos[\[ScriptL]]*(130860 - 153370*\[Nu] + 45*(-13802 + 6811*\[Nu])*
          Cos[2*\[ScriptL]] + (516276 - 270774*\[Nu])*Cos[4*\[ScriptL]] - 
         302526*Cos[6*\[ScriptL]] + 117649*\[Nu]*Cos[6*\[ScriptL]]))/138240 + 
      et^6*(3 + ((-50 + 7*\[Nu])*Cos[2*\[ScriptL]])/96 - 
        ((-8 + 7*\[Nu])*Cos[4*\[ScriptL]])/15 + 
        (9*(-18 + 7*\[Nu])*Cos[6*\[ScriptL]])/160) + 
      (et^8*(22680 - 7*(477 + 7*\[Nu])*Cos[2*\[ScriptL]] + 
         112*(-27 + 14*\[Nu])*Cos[4*\[ScriptL]] + 8019*Cos[6*\[ScriptL]] - 
         5103*\[Nu]*Cos[6*\[ScriptL]] - 9216*Cos[8*\[ScriptL]] + 
         3584*\[Nu]*Cos[8*\[ScriptL]]))/7560 + 
      et^2*(2 + (3 - (7*\[Nu])/6)*Sin[\[ScriptL]]^2) - 
      (et^3*Cos[\[ScriptL]]*(8 + (-18 + 7*\[Nu])*Sin[\[ScriptL]]^2))/4) + 
    x*\[Epsilon]^4*((\[Nu]*(171 + 4*\[Nu]))/36 - 
      (et*(828 - 639*\[Nu] + 35*\[Nu]^2)*Cos[\[ScriptL]])/72 - 
      (et^3*Cos[\[ScriptL]]*(2904 - 865*\[Nu] - 59*\[Nu]^2 + 
         (-432 - 495*\[Nu] + 59*\[Nu]^2)*Cos[2*\[ScriptL]]))/96 + 
      (et^5*Cos[\[ScriptL]]*(-9*(17652 - 30529*\[Nu] + 1437*\[Nu]^2) + 
         4*(-173304 - 3303*\[Nu] + 6451*\[Nu]^2)*Cos[2*\[ScriptL]] + 
         (272340 + 56979*\[Nu] - 12871*\[Nu]^2)*Cos[4*\[ScriptL]]))/13824 + 
      (et^4*(27*(2892 - 1351*\[Nu] + \[Nu]^2) + 
         4*(-6876 + 1755*\[Nu] + 151*\[Nu]^2)*Cos[2*\[ScriptL]] + 
         (10116 + 3681*\[Nu] - 631*\[Nu]^2)*Cos[4*\[ScriptL]]))/1728 - 
      (et^7*Cos[\[ScriptL]]*(141550560 - 23479830*\[Nu] - 3265970*\[Nu]^2 + 
         45*(-2848848 - 792947*\[Nu] + 144239*\[Nu]^2)*Cos[2*\[ScriptL]] - 
         18*(-8539728 - 843611*\[Nu] + 330551*\[Nu]^2)*Cos[4*\[ScriptL]] - 
         70385904*Cos[6*\[ScriptL]] - 9072153*\[Nu]*Cos[6*\[ScriptL]] + 
         2725133*\[Nu]^2*Cos[6*\[ScriptL]]))/1658880 + 
      (et^6*(90*(9032 - 4447*\[Nu] + \[Nu]^2) - 
         5*(34422 - 19669*\[Nu] + 139*\[Nu]^2)*Cos[2*\[ScriptL]] + 
         2*(-106464 - 565*\[Nu] + 3835*\[Nu]^2)*Cos[4*\[ScriptL]] + 
         9*(18942 + 2895*\[Nu] - 785*\[Nu]^2)*Cos[6*\[ScriptL]]))/11520 + 
      (et^8*(277886700 - 140397705*\[Nu] + 14175*\[Nu]^2 + 
         28*(-1955061 + 1015992*\[Nu] + 1702*\[Nu]^2)*Cos[2*\[ScriptL]] - 
         28*(214884 - 622053*\[Nu] + 35699*\[Nu]^2)*Cos[4*\[ScriptL]] - 
         104119668*Cos[6*\[ScriptL]] - 11889504*\[Nu]*Cos[6*\[ScriptL]] + 
         4159512*\[Nu]^2*Cos[6*\[ScriptL]] + 86235228*Cos[8*\[ScriptL]] + 
         10273869*\[Nu]*Cos[8*\[ScriptL]] - 3221771*\[Nu]^2*
          Cos[8*\[ScriptL]]))/2903040 + et^2*(18 - (17*\[Nu])/3 + 
        (4 - (31*\[Nu])/4 + (11*\[Nu]^2)/18)*Sin[\[ScriptL]]^2)) + 
    SO^2*(x^2*\[Epsilon]^6*((-3*\[Delta]*\[Kappa]A*(-5 + 8*\[Nu]) - 
          3*\[Kappa]S*(-5 + 18*\[Nu] + 6*\[Nu]^2) + 7*\[Chi]A^2 - 
          22*\[Nu]*\[Chi]A^2 - 36*\[Nu]^2*\[Chi]A^2 + 14*\[Delta]*
           (1 - 2*\[Nu])*\[Chi]A*\[Chi]S + 7*\[Chi]S^2 - 34*\[Nu]*\[Chi]S^2 - 
          8*\[Nu]^2*\[Chi]S^2)/18 + 
        (et*(3*\[Delta]*\[Kappa]A*(-86 + 39*\[Nu]) - 3*\[Kappa]S*
            (86 - 211*\[Nu] + 26*\[Nu]^2) - 2*(209 - 875*\[Nu] + 78*\[Nu]^2)*
            \[Chi]A^2 + 4*\[Delta]*(-209 + 184*\[Nu])*\[Chi]A*\[Chi]S - 
           2*(209 - 329*\[Nu] + 92*\[Nu]^2)*\[Chi]S^2)*Cos[\[ScriptL]])/18 + 
        (et^2*(2100*\[Kappa]S - 5337*\[Kappa]S*\[Nu] + 702*\[Kappa]S*
            \[Nu]^2 + 3140*\[Chi]A^2 - 13169*\[Nu]*\[Chi]A^2 + 
           1404*\[Nu]^2*\[Chi]A^2 + 3140*\[Chi]S^2 - 4193*\[Nu]*\[Chi]S^2 + 
           896*\[Nu]^2*\[Chi]S^2 + \[Delta]*(-3*\[Kappa]A*(-700 + 379*
                \[Nu]) + 2*(3140 - 2401*\[Nu])*\[Chi]A*\[Chi]S) + 
           (3*\[Delta]*\[Kappa]A*(-64 + 39*\[Nu]) - 3*\[Kappa]S*
              (64 - 167*\[Nu] + 10*\[Nu]^2) - 368*\[Chi]A^2 + 
             1565*\[Nu]*\[Chi]A^2 - 60*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-368 + 349*\[Nu])*\[Chi]A*\[Chi]S - 368*\[Chi]S^2 + 
             605*\[Nu]*\[Chi]S^2 - 224*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]]))/
         72 - (et^3*Cos[\[ScriptL]]*(696*\[Kappa]S + \[Delta]*\[Kappa]A*
            (696 - 301*\[Nu]) - 1693*\[Kappa]S*\[Nu] + 290*\[Kappa]S*
            \[Nu]^2 + 1044*\[Chi]A^2 - 4381*\[Nu]*\[Chi]A^2 + 
           580*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(1044 - 829*\[Nu])*\[Chi]A*
            \[Chi]S + 1044*\[Chi]S^2 - 1453*\[Nu]*\[Chi]S^2 + 
           336*\[Nu]^2*\[Chi]S^2 + (-6*\[Kappa]S*(11 - 24*\[Nu] + 6*
                \[Nu]^2) - 70*\[Chi]A^2 + 277*\[Nu]*\[Chi]A^2 - 
             72*\[Nu]^2*\[Chi]A^2 - 70*\[Chi]S^2 + 97*\[Nu]*\[Chi]S^2 + 
             8*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(\[Kappa]A*(-33 + 6*\[Nu]) + 
               (-70 + 47*\[Nu])*\[Chi]A*\[Chi]S))*Cos[2*\[ScriptL]]))/12 + 
        (et^4*(9*(\[Delta]*\[Kappa]A*(9276 - 4617*\[Nu]) + 
             3*\[Kappa]S*(3092 - 7723*\[Nu] + 1298*\[Nu]^2) + 
             13804*\[Chi]A^2 - 58183*\[Nu]*\[Chi]A^2 + 7788*\[Nu]^2*
              \[Chi]A^2 + 2*\[Delta]*(13804 - 10211*\[Nu])*\[Chi]A*\[Chi]S + 
             13804*\[Chi]S^2 - 17455*\[Nu]*\[Chi]S^2 + 3760*\[Nu]^2*
              \[Chi]S^2) + 8*(30*\[Delta]*\[Kappa]A*(-95 + 42*\[Nu]) - 
             6*\[Kappa]S*(475 - 1160*\[Nu] + 196*\[Nu]^2) - 4258*\[Chi]A^2 + 
             17875*\[Nu]*\[Chi]A^2 - 2352*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-4258 + 3377*\[Nu])*\[Chi]A*\[Chi]S - 
             4258*\[Chi]S^2 + 5911*\[Nu]*\[Chi]S^2 - 1360*\[Nu]^2*\[Chi]S^2)*
            Cos[2*\[ScriptL]] + (\[Delta]*\[Kappa]A*(6492 - 1971*\[Nu]) + 
             3*\[Kappa]S*(2164 - 4985*\[Nu] + 934*\[Nu]^2) + 8396*\[Chi]A^2 - 
             34217*\[Nu]*\[Chi]A^2 + 5604*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(8396 - 6613*\[Nu])*\[Chi]A*\[Chi]S + 
             8396*\[Chi]S^2 - 12593*\[Nu]*\[Chi]S^2 + 1520*\[Nu]^2*\[Chi]S^2)*
            Cos[4*\[ScriptL]]))/864 - (et^5*Cos[\[ScriptL]]*
          (9*(-3*\[Delta]*\[Kappa]A*(-9586 + 4541*\[Nu]) + 
             3*\[Kappa]S*(9586 - 23713*\[Nu] + 4382*\[Nu]^2) + 
             2*(22315 - 94705*\[Nu] + 13146*\[Nu]^2)*\[Chi]A^2 + 
             4*\[Delta]*(22315 - 17024*\[Nu])*\[Chi]A*\[Chi]S + 
             2*(22315 - 28603*\[Nu] + 7348*\[Nu]^2)*\[Chi]S^2) + 
           8*(-3*\[Delta]*\[Kappa]A*(-9556 + 3759*\[Nu]) + 
             3*\[Kappa]S*(9556 - 22871*\[Nu] + 3910*\[Nu]^2) + 
             40424*\[Chi]A^2 - 167783*\[Nu]*\[Chi]A^2 + 23460*\[Nu]^2*
              \[Chi]A^2 + 2*\[Delta]*(40424 - 32203*\[Nu])*\[Chi]A*\[Chi]S + 
             40424*\[Chi]S^2 - 58319*\[Nu]*\[Chi]S^2 + 10976*\[Nu]^2*
              \[Chi]S^2)*Cos[2*\[ScriptL]] + 
           (3*\[Delta]*\[Kappa]A*(-29810 + 9813*\[Nu]) - 3*\[Kappa]S*
              (29810 - 69433*\[Nu] + 12110*\[Nu]^2) - 
             2*(59771 - 244613*\[Nu] + 36330*\[Nu]^2)*\[Chi]A^2 + 
             4*\[Delta]*(-59771 + 48244*\[Nu])*\[Chi]A*\[Chi]S - 
             2*(59771 - 90959*\[Nu] + 13556*\[Nu]^2)*\[Chi]S^2)*
            Cos[4*\[ScriptL]]))/3456 - (et^7*Cos[\[ScriptL]]*
          (41613540*\[Delta]*\[Kappa]A + 41613540*\[Kappa]S - 
           17936820*\[Delta]*\[Kappa]A*\[Nu] - 101163900*\[Kappa]S*\[Nu] + 
           18428640*\[Kappa]S*\[Nu]^2 + 60986620*\[Chi]A^2 - 
           256075390*\[Nu]*\[Chi]A^2 + 36857280*\[Nu]^2*\[Chi]A^2 + 
           121973240*\[Delta]*\[Chi]A*\[Chi]S - 93926020*\[Delta]*\[Nu]*
            \[Chi]A*\[Chi]S + 60986620*\[Chi]S^2 - 81797110*\[Nu]*\[Chi]S^2 + 
           17520400*\[Nu]^2*\[Chi]S^2 - 45*(\[Delta]*\[Kappa]A*
              (221586 - 62020*\[Nu]) + \[Kappa]S*(221586 - 505192*\[Nu] + 
               70820*\[Nu]^2) + 269286*\[Chi]A^2 - 1069717*\[Nu]*\[Chi]A^2 + 
             141640*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(269286 - 236287*\[Nu])*
              \[Chi]A*\[Chi]S + 269286*\[Chi]S^2 - 480001*\[Nu]*\[Chi]S^2 + 
             49848*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]] + 
           18*(\[Delta]*\[Kappa]A*(1245246 - 465118*\[Nu]) + 
             2*\[Kappa]S*(622623 - 1477805*\[Nu] + 243856*\[Nu]^2) + 
             1706226*\[Chi]A^2 - 7030369*\[Nu]*\[Chi]A^2 + 975424*\[Nu]^2*
              \[Chi]A^2 + 2*\[Delta]*(1706226 - 1381327*\[Nu])*\[Chi]A*
              \[Chi]S + 1706226*\[Chi]S^2 - 2557189*\[Nu]*\[Chi]S^2 + 
             435192*\[Nu]^2*\[Chi]S^2)*Cos[4*\[ScriptL]] - 
           11940918*\[Delta]*\[Kappa]A*Cos[6*\[ScriptL]] - 
           11940918*\[Kappa]S*Cos[6*\[ScriptL]] + 4132044*\[Delta]*\[Kappa]A*
            \[Nu]*Cos[6*\[ScriptL]] + 28013880*\[Kappa]S*\[Nu]*
            Cos[6*\[ScriptL]] - 4628076*\[Kappa]S*\[Nu]^2*Cos[6*\[ScriptL]] - 
           16216498*\[Chi]A^2*Cos[6*\[ScriptL]] + 66486367*\[Nu]*\[Chi]A^2*
            Cos[6*\[ScriptL]] - 9256152*\[Nu]^2*\[Chi]A^2*Cos[6*\[ScriptL]] - 
           32432996*\[Delta]*\[Chi]A*\[Chi]S*Cos[6*\[ScriptL]] + 
           26584762*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S*Cos[6*\[ScriptL]] - 
           16216498*\[Chi]S^2*Cos[6*\[ScriptL]] + 24964387*\[Nu]*\[Chi]S^2*
            Cos[6*\[ScriptL]] - 4073896*\[Nu]^2*\[Chi]S^2*Cos[6*\[ScriptL]]))/
         207360 + (et^6*(1165200*\[Delta]*\[Kappa]A + 1165200*\[Kappa]S - 
           563910*\[Delta]*\[Kappa]A*\[Nu] - 2894310*\[Kappa]S*\[Nu] + 
           525060*\[Kappa]S*\[Nu]^2 + 1727120*\[Chi]A^2 - 
           7294430*\[Nu]*\[Chi]A^2 + 1050120*\[Nu]^2*\[Chi]A^2 + 
           3454240*\[Delta]*\[Chi]A*\[Chi]S - 2521820*\[Delta]*\[Nu]*\[Chi]A*
            \[Chi]S + 1727120*\[Chi]S^2 - 2135870*\[Nu]*\[Chi]S^2 + 
           458240*\[Nu]^2*\[Chi]S^2 + 5*(2*\[Delta]*\[Kappa]A*
              (-27207 + 12173*\[Nu]) - 2*\[Kappa]S*(27207 - 66587*\[Nu] + 
               12106*\[Nu]^2) - 81558*\[Chi]A^2 + 343801*\[Nu]*\[Chi]A^2 - 
             48424*\[Nu]^2*\[Chi]A^2 + 46*\[Delta]*(-3546 + 2729*\[Nu])*
              \[Chi]A*\[Chi]S - 81558*\[Chi]S^2 + 107965*\[Nu]*\[Chi]S^2 - 
             24936*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]] + 
           2*(\[Delta]*\[Kappa]A*(-59184 + 24887*\[Nu]) + 
             \[Kappa]S*(-59184 + 143255*\[Nu] - 23962*\[Nu]^2) - 
             83856*\[Chi]A^2 + 349223*\[Nu]*\[Chi]A^2 - 47924*\[Nu]^2*
              \[Chi]A^2 + 2*\[Delta]*(-83856 + 66455*\[Nu])*\[Chi]A*\[Chi]S - 
             83856*\[Chi]S^2 + 119111*\[Nu]*\[Chi]S^2 - 23136*\[Nu]^2*
              \[Chi]S^2)*Cos[4*\[ScriptL]] + 113958*\[Delta]*\[Kappa]A*
            Cos[6*\[ScriptL]] + 113958*\[Kappa]S*Cos[6*\[ScriptL]] - 
           38754*\[Delta]*\[Kappa]A*\[Nu]*Cos[6*\[ScriptL]] - 
           266670*\[Kappa]S*\[Nu]*Cos[6*\[ScriptL]] + 44964*\[Kappa]S*\[Nu]^2*
            Cos[6*\[ScriptL]] + 154062*\[Chi]A^2*Cos[6*\[ScriptL]] - 
           631341*\[Nu]*\[Chi]A^2*Cos[6*\[ScriptL]] + 89928*\[Nu]^2*\[Chi]A^2*
            Cos[6*\[ScriptL]] + 308124*\[Delta]*\[Chi]A*\[Chi]S*
            Cos[6*\[ScriptL]] - 251190*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S*
            Cos[6*\[ScriptL]] + 154062*\[Chi]S^2*Cos[6*\[ScriptL]] - 
           236097*\[Nu]*\[Chi]S^2*Cos[6*\[ScriptL]] + 37512*\[Nu]^2*\[Chi]S^2*
            Cos[6*\[ScriptL]]))/5760 + (et^8*(502063380*\[Delta]*\[Kappa]A + 
           502063380*\[Kappa]S - 239664285*\[Delta]*\[Kappa]A*\[Nu] - 
           1243791045*\[Kappa]S*\[Nu] + 234312750*\[Kappa]S*\[Nu]^2 + 
           742284900*\[Chi]A^2 - 3138543135*\[Nu]*\[Chi]A^2 + 
           468625500*\[Nu]^2*\[Chi]A^2 + 1484569800*\[Delta]*\[Chi]A*
            \[Chi]S - 1076015430*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S + 
           742284900*\[Chi]S^2 - 906611895*\[Nu]*\[Chi]S^2 + 
           193772880*\[Nu]^2*\[Chi]S^2 - 14*(\[Delta]*\[Kappa]A*
              (7742022 - 3478383*\[Nu]) + 3*\[Kappa]S*(2580674 - 6320809*
                \[Nu] + 1188182*\[Nu]^2) + 2*(5780087 - 24401681*\[Nu] + 
               3564546*\[Nu]^2)*\[Chi]A^2 + 4*\[Delta]*(5780087 - 4376392*
                \[Nu])*\[Chi]A*\[Chi]S + 2*(5780087 - 7471451*\[Nu] + 1701332*
                \[Nu]^2)*\[Chi]S^2)*Cos[2*\[ScriptL]] - 
           28*(-15*\[Delta]*\[Kappa]A*(-73460 + 33381*\[Nu]) + 
             15*\[Kappa]S*(73460 - 180301*\[Nu] + 33758*\[Nu]^2) + 
             1655420*\[Chi]A^2 - 6995417*\[Nu]*\[Chi]A^2 + 1012740*\[Nu]^2*
              \[Chi]A^2 + 2*\[Delta]*(1655420 - 1254709*\[Nu])*\[Chi]A*
              \[Chi]S + 1655420*\[Chi]S^2 - 2135681*\[Nu]*\[Chi]S^2 + 
             496688*\[Nu]^2*\[Chi]S^2)*Cos[4*\[ScriptL]] - 
           50516460*\[Delta]*\[Kappa]A*Cos[6*\[ScriptL]] - 
           50516460*\[Kappa]S*Cos[6*\[ScriptL]] + 19785870*\[Delta]*\[Kappa]A*
            \[Nu]*Cos[6*\[ScriptL]] + 120818790*\[Kappa]S*\[Nu]*
            Cos[6*\[ScriptL]] - 19603620*\[Kappa]S*\[Nu]^2*
            Cos[6*\[ScriptL]] - 68822460*\[Chi]A^2*Cos[6*\[ScriptL]] + 
           283966884*\[Nu]*\[Chi]A^2*Cos[6*\[ScriptL]] - 39207240*\[Nu]^2*
            \[Chi]A^2*Cos[6*\[ScriptL]] - 137644920*\[Delta]*\[Chi]A*\[Chi]S*
            Cos[6*\[ScriptL]] + 110828736*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S*
            Cos[6*\[ScriptL]] - 68822460*\[Chi]S^2*Cos[6*\[ScriptL]] + 
           102151692*\[Nu]*\[Chi]S^2*Cos[6*\[ScriptL]] - 17186256*\[Nu]^2*
            \[Chi]S^2*Cos[6*\[ScriptL]] + 59306388*\[Delta]*\[Kappa]A*
            Cos[8*\[ScriptL]] + 59306388*\[Kappa]S*Cos[8*\[ScriptL]] - 
           20744667*\[Delta]*\[Kappa]A*\[Nu]*Cos[8*\[ScriptL]] - 
           139357443*\[Kappa]S*\[Nu]*Cos[8*\[ScriptL]] + 22700994*\[Kappa]S*
            \[Nu]^2*Cos[8*\[ScriptL]] + 80677156*\[Chi]A^2*
            Cos[8*\[ScriptL]] - 330811213*\[Nu]*\[Chi]A^2*Cos[8*\[ScriptL]] + 
           45401988*\[Nu]^2*\[Chi]A^2*Cos[8*\[ScriptL]] + 161354312*\[Delta]*
            \[Chi]A*\[Chi]S*Cos[8*\[ScriptL]] - 132720802*\[Delta]*\[Nu]*
            \[Chi]A*\[Chi]S*Cos[8*\[ScriptL]] + 80677156*\[Chi]S^2*
            Cos[8*\[ScriptL]] - 124618213*\[Nu]*\[Chi]S^2*Cos[8*\[ScriptL]] + 
           20584336*\[Nu]^2*\[Chi]S^2*Cos[8*\[ScriptL]]))/1451520) + 
      x*\[Epsilon]^4*((\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
          4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
            2*\[Chi]A*\[Chi]S))/2 - et*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
          \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
          \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*Cos[\[ScriptL]] + 
        (et^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(21 - 4*Cos[2*\[ScriptL]] - 
           2*Cos[4*\[ScriptL]]))/6 - (et^5*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
           \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
           \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*(874*Cos[\[ScriptL]] + 
           153*Cos[3*\[ScriptL]] + 125*Cos[5*\[ScriptL]]))/384 - 
        (et^6*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(-1200 + 215*Cos[2*\[ScriptL]] + 
           64*Cos[4*\[ScriptL]] + 81*Cos[6*\[ScriptL]]))/240 - 
        (et^7*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(134845*Cos[\[ScriptL]] + 
           24543*Cos[3*\[ScriptL]] + 8125*Cos[5*\[ScriptL]] + 
           16807*Cos[7*\[ScriptL]]))/46080 - 
        (et^8*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(-16380 + 2821*Cos[2*\[ScriptL]] + 
           952*Cos[4*\[ScriptL]] + 243*Cos[6*\[ScriptL]] + 
           1024*Cos[8*\[ScriptL]]))/2520 + 
        (et^2*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(3 + 2*Sin[\[ScriptL]]^2))/2 + 
        (et^3*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*Cos[\[ScriptL]]*(-4 + 3*Sin[\[ScriptL]]^2))/
         2)) + SO*(Sqrt[x]*\[Epsilon]^3*
       ((-2*(\[Delta]*\[Chi]A + \[Chi]S + \[Nu]*\[Chi]S))/3 + 
        (4*et*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*Cos[\[ScriptL]])/
         3 + 2*et^3*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*
         Cos[\[ScriptL]]^3 + et^4*((-25*\[Delta]*\[Chi]A - 25*\[Chi]S + 
            11*\[Nu]*\[Chi]S)/4 + (5*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*
              \[Chi]S)*Cos[2*\[ScriptL]])/9 + 
          (4*(2*\[Delta]*\[Chi]A + 2*\[Chi]S - \[Nu]*\[Chi]S)*
            Cos[4*\[ScriptL]])/9) + (et^5*(2*\[Delta]*\[Chi]A - 
           (-2 + \[Nu])*\[Chi]S)*(514*Cos[\[ScriptL]] + 
           81*Cos[3*\[ScriptL]] + 125*Cos[5*\[ScriptL]]))/288 + 
        (et^6*(-925*\[Delta]*\[Chi]A - 925*\[Chi]S + 425*\[Nu]*\[Chi]S + 
           80*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*Cos[2*\[ScriptL]] + 
           16*(2*\[Delta]*\[Chi]A + 2*\[Chi]S - \[Nu]*\[Chi]S)*
            Cos[4*\[ScriptL]] + 108*\[Delta]*\[Chi]A*Cos[6*\[ScriptL]] + 
           108*\[Chi]S*Cos[6*\[ScriptL]] - 54*\[Nu]*\[Chi]S*
            Cos[6*\[ScriptL]]))/120 + (et^7*(2*\[Delta]*\[Chi]A - 
           (-2 + \[Nu])*\[Chi]S)*(70165*Cos[\[ScriptL]] + 
           13203*Cos[3*\[ScriptL]] + 625*Cos[5*\[ScriptL]] + 
           16807*Cos[7*\[ScriptL]]))/34560 + 
        (et^8*(-540225*\[Delta]*\[Chi]A - 540225*\[Chi]S + 
           253575*\[Nu]*\[Chi]S + 44912*(2*\[Delta]*\[Chi]A - 
             (-2 + \[Nu])*\[Chi]S)*Cos[2*\[ScriptL]] + 
           16352*(2*\[Delta]*\[Chi]A + 2*\[Chi]S - \[Nu]*\[Chi]S)*
            Cos[4*\[ScriptL]] - 11664*\[Delta]*\[Chi]A*Cos[6*\[ScriptL]] - 
           11664*\[Chi]S*Cos[6*\[ScriptL]] + 5832*\[Nu]*\[Chi]S*
            Cos[6*\[ScriptL]] + 65536*\[Delta]*\[Chi]A*Cos[8*\[ScriptL]] + 
           65536*\[Chi]S*Cos[8*\[ScriptL]] - 32768*\[Nu]*\[Chi]S*
            Cos[8*\[ScriptL]]))/60480 + et^2*(-3*\[Delta]*\[Chi]A + 
          (-3 + \[Nu])*\[Chi]S + ((-8*\[Delta]*\[Chi]A + 4*(-2 + \[Nu])*
              \[Chi]S)*Sin[\[ScriptL]]^2)/3)) + x^(3/2)*\[Epsilon]^5*
       ((-4*\[Nu]*(\[Delta]*\[Chi]A + \[Chi]S - \[Nu]*\[Chi]S))/3 + 
        et*(\[Delta]*(24 - 7*\[Nu])*\[Chi]A + (24 - 31*\[Nu] + 2*\[Nu]^2)*
           \[Chi]S)*Cos[\[ScriptL]] + (et^3*Cos[\[ScriptL]]*
          (\[Delta]*(972 - 307*\[Nu])*\[Chi]A + (972 - 1159*\[Nu] + 
             104*\[Nu]^2)*\[Chi]S + (-84*\[Delta]*\[Chi]A + 9*\[Delta]*\[Nu]*
              \[Chi]A - 84*\[Chi]S + 93*\[Nu]*\[Chi]S)*Cos[2*\[ScriptL]]))/
         12 + (et^5*Cos[\[ScriptL]]*(14424*\[Delta]*\[Chi]A - 
           5929*\[Delta]*\[Nu]*\[Chi]A + 14424*\[Chi]S - 
           16993*\[Nu]*\[Chi]S + 2582*\[Nu]^2*\[Chi]S - 
           4*(\[Delta]*(-4996 + 1251*\[Nu])*\[Chi]A + 
             (-4996 + 5887*\[Nu] - 276*\[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]] + 
           (\[Delta]*(-7240 + 1437*\[Nu])*\[Chi]A + (-7240 + 8677*\[Nu] - 198*
                \[Nu]^2)*\[Chi]S)*Cos[4*\[ScriptL]]))/192 + 
        (et^4*(-6228*\[Delta]*\[Chi]A + 2085*\[Delta]*\[Nu]*\[Chi]A - 
           6228*\[Chi]S + 6873*\[Nu]*\[Chi]S - 630*\[Nu]^2*\[Chi]S - 
           6*(\[Delta]*(-316 + 97*\[Nu])*\[Chi]A + (-316 + 373*\[Nu] - 30*
                \[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]] + 
           (\[Delta]*(-516 + 97*\[Nu])*\[Chi]A + (-516 + 613*\[Nu] - 14*
                \[Nu]^2)*\[Chi]S)*Cos[4*\[ScriptL]]))/48 + 
        (et^7*Cos[\[ScriptL]]*(30*(\[Delta]*(185332 - 55855*\[Nu])*\[Chi]A + 
             (185332 - 216587*\[Nu] + 17212*\[Nu]^2)*\[Chi]S) + 
           5*(\[Delta]*(-516228 + 67357*\[Nu])*\[Chi]A + 
             (-516228 + 607585*\[Nu] + 22432*\[Nu]^2)*\[Chi]S)*
            Cos[2*\[ScriptL]] + (\[Delta]*(4098216 - 877550*\[Nu])*\[Chi]A + 
             2*(2049108 - 2412883*\[Nu] + 50012*\[Nu]^2)*\[Chi]S)*
            Cos[4*\[ScriptL]] + 3*(\[Delta]*(-648292 + 129365*\[Nu])*
              \[Chi]A + (-648292 + 777657*\[Nu] - 13728*\[Nu]^2)*\[Chi]S)*
            Cos[6*\[ScriptL]]))/23040 + 
        (et^6*(30*(13*\[Delta]*(-1764 + 613*\[Nu])*\[Chi]A + 
             (-22932 + 25501*\[Nu] - 2682*\[Nu]^2)*\[Chi]S) - 
           5*(\[Delta]*(-31344 + 10823*\[Nu])*\[Chi]A + 
             (-31344 + 36887*\[Nu] - 4048*\[Nu]^2)*\[Chi]S)*
            Cos[2*\[ScriptL]] - 4*(\[Delta]*(-26232 + 6587*\[Nu])*\[Chi]A + 
             (-26232 + 30419*\[Nu] - 1228*\[Nu]^2)*\[Chi]S)*
            Cos[4*\[ScriptL]] + 9*(\[Delta]*(-9272 + 1857*\[Nu])*\[Chi]A + 
             (-9272 + 11129*\[Nu] - 228*\[Nu]^2)*\[Chi]S)*Cos[6*\[ScriptL]]))/
         2880 + (et^8*(-29638980*\[Delta]*\[Chi]A + 10493385*\[Delta]*\[Nu]*
            \[Chi]A - 29638980*\[Chi]S + 33076365*\[Nu]*\[Chi]S - 
           3705030*\[Nu]^2*\[Chi]S - 14*(\[Delta]*(-441508 + 154899*\[Nu])*
              \[Chi]A + (-441508 + 514807*\[Nu] - 58686*\[Nu]^2)*\[Chi]S)*
            Cos[2*\[ScriptL]] - 14*(5*\[Delta]*(-23276 + 9415*\[Nu])*
              \[Chi]A + (-116380 + 138495*\[Nu] - 21442*\[Nu]^2)*\[Chi]S)*
            Cos[4*\[ScriptL]] + 5249880*\[Delta]*\[Chi]A*Cos[6*\[ScriptL]] - 
           1103490*\[Delta]*\[Nu]*\[Chi]A*Cos[6*\[ScriptL]] + 
           5249880*\[Chi]S*Cos[6*\[ScriptL]] - 6081210*\[Nu]*\[Chi]S*
            Cos[6*\[ScriptL]] + 71172*\[Nu]^2*\[Chi]S*Cos[6*\[ScriptL]] - 
           4826212*\[Delta]*\[Chi]A*Cos[8*\[ScriptL]] + 954901*\[Delta]*\[Nu]*
            \[Chi]A*Cos[8*\[ScriptL]] - 4826212*\[Chi]S*Cos[8*\[ScriptL]] + 
           5781113*\[Nu]*\[Chi]S*Cos[8*\[ScriptL]] - 85214*\[Nu]^2*\[Chi]S*
            Cos[8*\[ScriptL]]))/80640 + 
        (et^2*(10*\[Delta]*(-24 + 7*\[Nu])*\[Chi]A - 
           2*(120 - 125*\[Nu] + 6*\[Nu]^2)*\[Chi]S + 
           3*(\[Delta]*(-20 + 7*\[Nu])*\[Chi]A + (-20 + 27*\[Nu] - 2*\[Nu]^2)*
              \[Chi]S)*Sin[\[ScriptL]]^2))/6))


(* ::Subsection::Closed:: *)
(*dr/dt*)


rd\[LetterSpace]xel = x^(3/2)*\[Epsilon]^2*((-7*et*\[Nu]*Sin[\[ScriptL]])/6 - 
      (7*et^2*\[Nu]*Cos[\[ScriptL]]*Sin[\[ScriptL]])/3 - 
      (et^3*(8 + 7*\[Nu] + 21*\[Nu]*Cos[2*\[ScriptL]])*Sin[\[ScriptL]])/8 - 
      (et^5*(4*(648 + 385*\[Nu])*Cos[2*\[ScriptL]] + 
         7*(288 + 115*\[Nu] + 625*\[Nu]*Cos[4*\[ScriptL]]))*Sin[\[ScriptL]])/
       1152 - (et^7*(324720 + 82418*\[Nu] + 3*(156480 + 55027*\[Nu])*
          Cos[2*\[ScriptL]] + 6*(75000 + 9653*\[Nu])*Cos[4*\[ScriptL]] + 
         823543*\[Nu]*Cos[6*\[ScriptL]])*Sin[\[ScriptL]])/138240 - 
      (et^4*(9 - 7*\[Nu] + 28*\[Nu]*Cos[2*\[ScriptL]])*Sin[2*\[ScriptL]])/9 - 
      (et^6*(40 + 301*\[Nu] - 64*(-5 + 7*\[Nu])*Cos[2*\[ScriptL]] + 
         567*\[Nu]*Cos[4*\[ScriptL]])*Sin[2*\[ScriptL]])/120 - 
      (et^8*(1341 - 2194*\[Nu] + 96*(-3 + 52*\[Nu])*Cos[2*\[ScriptL]] - 
         2187*(-1 + 2*\[Nu])*Cos[4*\[ScriptL]] + 4096*\[Nu]*
          Cos[6*\[ScriptL]])*Sin[2*\[ScriptL]])/540) + 
    x^(5/2)*\[Epsilon]^4*((et*(-144 + 117*\[Nu] + 35*\[Nu]^2)*
        Sin[\[ScriptL]])/72 + (et^3*(-2460 + 1205*\[Nu] + 59*\[Nu]^2 + 
         3*(-1404 + 261*\[Nu] + 59*\[Nu]^2)*Cos[2*\[ScriptL]])*
        Sin[\[ScriptL]])/96 + (et^5*(-866304 + 377919*\[Nu] + 12809*\[Nu]^2 + 
         4*(-363996 + 124155*\[Nu] + 6389*\[Nu]^2)*Cos[2*\[ScriptL]] + 
         5*(-393840 + 37521*\[Nu] + 12871*\[Nu]^2)*Cos[4*\[ScriptL]])*
        Sin[\[ScriptL]])/13824 + (et^7*(-183011256 + 75468438*\[Nu] + 
         2266666*\[Nu]^2 + 3*(-105410124 + 36327417*\[Nu] + 1524839*\[Nu]^2)*
          Cos[2*\[ScriptL]] + 6*(-47091276 + 17152383*\[Nu] + 492001*\[Nu]^2)*
          Cos[4*\[ScriptL]] - 607056156*Cos[6*\[ScriptL]] + 
         25437573*\[Nu]*Cos[6*\[ScriptL]] + 19075931*\[Nu]^2*
          Cos[6*\[ScriptL]])*Sin[\[ScriptL]])/1658880 + 
      (et^2*(-342 + 99*\[Nu] + 22*\[Nu]^2)*Sin[2*\[ScriptL]])/36 + 
      (et^4*(882 + 1917*\[Nu] - 151*\[Nu]^2 + (-17892 + 2367*\[Nu] + 
           631*\[Nu]^2)*Cos[2*\[ScriptL]])*Sin[2*\[ScriptL]])/216 + 
      (et^6*(-355494 + 38681*\[Nu] + 10945*\[Nu]^2 + 
         (268032 + 89812*\[Nu] - 15340*\[Nu]^2)*Cos[2*\[ScriptL]] + 
         27*(-24774 + 1641*\[Nu] + 785*\[Nu]^2)*Cos[4*\[ScriptL]])*
        Sin[2*\[ScriptL]])/2880 + (et^8*(32234553 + 8347680*\[Nu] - 
         1565774*\[Nu]^2 + 3*(-36776268 + 107535*\[Nu] + 1240519*\[Nu]^2)*
          Cos[2*\[ScriptL]] - 243*(-280209 - 54840*\[Nu] + 12838*\[Nu]^2)*
          Cos[4*\[ScriptL]] - 102160476*Cos[6*\[ScriptL]] + 
         2112435*\[Nu]*Cos[6*\[ScriptL]] + 3221771*\[Nu]^2*Cos[6*\[ScriptL]])*
        Sin[2*\[ScriptL]])/181440) + Sqrt[x]*(et*Sin[\[ScriptL]] + 
      (et^8*(302*Cos[\[ScriptL]] + 309*Cos[3*\[ScriptL]] - 
         139*Cos[5*\[ScriptL]] + 2048*Cos[7*\[ScriptL]])*Sin[\[ScriptL]])/
       315 + et^2*Sin[2*\[ScriptL]] - 
      (3*et^3*(Sin[\[ScriptL]] - 3*Sin[3*\[ScriptL]]))/8 - 
      (2*et^4*(Sin[2*\[ScriptL]] - 2*Sin[4*\[ScriptL]]))/3 + 
      (5*et^5*(2*Sin[\[ScriptL]] - 81*Sin[3*\[ScriptL]] + 
         125*Sin[5*\[ScriptL]]))/384 + 
      (et^6*(5*Sin[2*\[ScriptL]] - 64*Sin[4*\[ScriptL]] + 
         81*Sin[6*\[ScriptL]]))/40 + 
      (7*et^7*(-5*Sin[\[ScriptL]] + 2187*Sin[3*\[ScriptL]] - 
         15625*Sin[5*\[ScriptL]] + 16807*Sin[7*\[ScriptL]]))/46080) + 
    SO^2*(x^(7/2)*\[Epsilon]^6*((et*(3*\[Delta]*\[Kappa]A*(-44 + 45*\[Nu]) - 
           3*\[Kappa]S*(44 - 133*\[Nu] + 50*\[Nu]^2) - 124*\[Chi]A^2 + 
           601*\[Nu]*\[Chi]A^2 - 300*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(-124 + 221*\[Nu])*\[Chi]A*\[Chi]S - 124*\[Chi]S^2 + 
           337*\[Nu]*\[Chi]S^2 - 88*\[Nu]^2*\[Chi]S^2)*Sin[\[ScriptL]])/36 - 
        (et^3*(1848*\[Kappa]S - 4769*\[Kappa]S*\[Nu] + 1186*\[Kappa]S*
            \[Nu]^2 + 2336*\[Chi]A^2 - 10105*\[Nu]*\[Chi]A^2 + 
           2372*\[Nu]^2*\[Chi]A^2 + 2336*\[Chi]S^2 - 3961*\[Nu]*\[Chi]S^2 + 
           872*\[Nu]^2*\[Chi]S^2 + \[Delta]*(\[Kappa]A*(1848 - 1073*\[Nu]) + 
             2*(2336 - 2361*\[Nu])*\[Chi]A*\[Chi]S) + 
           3*(\[Delta]*\[Kappa]A*(912 - 417*\[Nu]) + 3*\[Kappa]S*
              (304 - 747*\[Nu] + 150*\[Nu]^2) + 1240*\[Chi]A^2 - 
             5209*\[Nu]*\[Chi]A^2 + 900*\[Nu]^2*\[Chi]A^2 + 
             10*\[Delta]*(248 - 229*\[Nu])*\[Chi]A*\[Chi]S + 1240*\[Chi]S^2 - 
             2041*\[Nu]*\[Chi]S^2 + 424*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]])*
          Sin[\[ScriptL]])/48 + (et^5*(-847932*\[Delta]*\[Kappa]A - 
           847932*\[Kappa]S + 437805*\[Delta]*\[Kappa]A*\[Nu] + 
           2133669*\[Kappa]S*\[Nu] - 483666*\[Kappa]S*\[Nu]^2 - 
           1112548*\[Chi]A^2 + 4753291*\[Nu]*\[Chi]A^2 - 967332*\[Nu]^2*
            \[Chi]A^2 - 2225096*\[Delta]*\[Chi]A*\[Chi]S + 
           2066254*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S - 1112548*\[Chi]S^2 + 
           1763155*\[Nu]*\[Chi]S^2 - 375880*\[Nu]^2*\[Chi]S^2 + 
           4*(147*\[Delta]*\[Kappa]A*(-2240 + 1017*\[Nu]) - 
             3*\[Kappa]S*(109760 - 269353*\[Nu] + 54434*\[Nu]^2) - 
             450232*\[Chi]A^2 + 1895335*\[Nu]*\[Chi]A^2 - 326604*\[Nu]^2*
              \[Chi]A^2 + 2*\[Delta]*(-450232 + 404015*\[Nu])*\[Chi]A*
              \[Chi]S - 450232*\[Chi]S^2 + 713623*\[Nu]*\[Chi]S^2 - 
             149560*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]] + 
           5*(3*\[Delta]*\[Kappa]A*(-86620 + 35001*\[Nu]) - 
             3*\[Kappa]S*(86620 - 208241*\[Nu] + 36970*\[Nu]^2) - 
             359084*\[Chi]A^2 + 1491077*\[Nu]*\[Chi]A^2 - 221820*\[Nu]^2*
              \[Chi]A^2 + 2*\[Delta]*(-359084 + 312601*\[Nu])*\[Chi]A*
              \[Chi]S - 359084*\[Chi]S^2 + 570461*\[Nu]*\[Chi]S^2 - 
             111224*\[Nu]^2*\[Chi]S^2)*Cos[4*\[ScriptL]])*Sin[\[ScriptL]])/
         6912 - (et^6*Cos[\[ScriptL]]*(254790*\[Delta]*\[Kappa]A + 
           254790*\[Kappa]S - 111560*\[Delta]*\[Kappa]A*\[Nu] - 
           621140*\[Kappa]S*\[Nu] + 119152*\[Kappa]S*\[Nu]^2 + 
           343998*\[Chi]A^2 - 1439915*\[Nu]*\[Chi]A^2 + 238304*\[Nu]^2*
            \[Chi]A^2 + 687996*\[Delta]*\[Chi]A*\[Chi]S - 
           606034*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S + 343998*\[Chi]S^2 - 
           542111*\[Nu]*\[Chi]S^2 + 106104*\[Nu]^2*\[Chi]S^2 + 
           4*(5*\[Delta]*\[Kappa]A*(-2160 + 91*\[Nu]) + \[Kappa]S*
              (-10800 + 22055*\[Nu] + 1286*\[Nu]^2) - 14736*\[Chi]A^2 + 
             54215*\[Nu]*\[Chi]A^2 + 2572*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-14736 + 10199*\[Nu])*\[Chi]A*\[Chi]S - 
             14736*\[Chi]S^2 + 25127*\[Nu]*\[Chi]S^2 - 1248*\[Nu]^2*
              \[Chi]S^2)*Cos[2*\[ScriptL]] + 
           27*(\[Delta]*\[Kappa]A*(16550 - 6520*\[Nu]) + 2*\[Kappa]S*
              (8275 - 19810*\[Nu] + 3416*\[Nu]^2) + 22878*\[Chi]A^2 - 
             94755*\[Nu]*\[Chi]A^2 + 13664*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(22878 - 19697*\[Nu])*\[Chi]A*\[Chi]S + 
             22878*\[Chi]S^2 - 36151*\[Nu]*\[Chi]S^2 + 6904*\[Nu]^2*
              \[Chi]S^2)*Cos[4*\[ScriptL]])*Sin[\[ScriptL]])/720 + 
        (et^7*(-224322816*\[Delta]*\[Kappa]A - 224322816*\[Kappa]S + 
           110069478*\[Delta]*\[Kappa]A*\[Nu] + 558715110*\[Kappa]S*\[Nu] - 
           121629612*\[Kappa]S*\[Nu]^2 - 298715504*\[Chi]A^2 + 
           1270334030*\[Nu]*\[Chi]A^2 - 243259224*\[Nu]^2*\[Chi]A^2 - 
           597431008*\[Delta]*\[Chi]A*\[Chi]S + 535819772*\[Delta]*\[Nu]*
            \[Chi]A*\[Chi]S - 298715504*\[Chi]S^2 + 460347758*\[Nu]*
            \[Chi]S^2 - 96889712*\[Nu]^2*\[Chi]S^2 - 
           3*(\[Delta]*\[Kappa]A*(120533424 - 54155157*\[Nu]) + 
             3*\[Kappa]S*(40177808 - 98407335*\[Nu] + 19793166*\[Nu]^2) + 
             165156376*\[Chi]A^2 - 695411965*\[Nu]*\[Chi]A^2 + 
             118758996*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(165156376 - 145472189*
                \[Nu])*\[Chi]A*\[Chi]S + 165156376*\[Chi]S^2 - 
             256157917*\[Nu]*\[Chi]S^2 + 53499208*\[Nu]^2*\[Chi]S^2)*
            Cos[2*\[ScriptL]] + 6*(3*\[Delta]*\[Kappa]A*(-18035392 + 7548541*
                \[Nu]) - 3*\[Kappa]S*(18035392 - 43619325*\[Nu] + 8202954*
                \[Nu]^2) - 75482264*\[Chi]A^2 + 315524075*\[Nu]*\[Chi]A^2 - 
             49217724*\[Nu]^2*\[Chi]A^2 + 22*\[Delta]*(-6862024 + 5954921*
                \[Nu])*\[Chi]A*\[Chi]S - 75482264*\[Chi]S^2 + 
             117413243*\[Nu]*\[Chi]S^2 - 24208472*\[Nu]^2*\[Chi]S^2)*
            Cos[4*\[ScriptL]] - 410582256*\[Delta]*\[Kappa]A*
            Cos[6*\[ScriptL]] - 410582256*\[Kappa]S*Cos[6*\[ScriptL]] + 
           159109713*\[Delta]*\[Kappa]A*\[Nu]*Cos[6*\[ScriptL]] + 
           980274225*\[Kappa]S*\[Nu]*Cos[6*\[ScriptL]] - 165586722*\[Kappa]S*
            \[Nu]^2*Cos[6*\[ScriptL]] - 567004984*\[Chi]A^2*
            Cos[6*\[ScriptL]] + 2344096825*\[Nu]*\[Chi]A^2*
            Cos[6*\[ScriptL]] - 331173444*\[Nu]^2*\[Chi]A^2*
            Cos[6*\[ScriptL]] - 1134009968*\[Delta]*\[Chi]A*\[Chi]S*
            Cos[6*\[ScriptL]] + 969553522*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S*
            Cos[6*\[ScriptL]] - 567004984*\[Chi]S^2*Cos[6*\[ScriptL]] + 
           893476633*\[Nu]*\[Chi]S^2*Cos[6*\[ScriptL]] - 167717032*\[Nu]^2*
            \[Chi]S^2*Cos[6*\[ScriptL]])*Sin[\[ScriptL]])/829440 + 
        (et^2*(3*\[Delta]*\[Kappa]A*(-38 + 21*\[Nu]) - 3*\[Kappa]S*
            (38 - 97*\[Nu] + 23*\[Nu]^2) - 2*(74 - 317*\[Nu] + 69*\[Nu]^2)*
            \[Chi]A^2 + 8*\[Delta]*(-37 + 38*\[Nu])*\[Chi]A*\[Chi]S - 
           2*(74 - 131*\[Nu] + 29*\[Nu]^2)*\[Chi]S^2)*Sin[2*\[ScriptL]])/9 - 
        (et^4*(780*\[Kappa]S + \[Delta]*\[Kappa]A*(780 - 819*\[Nu]) - 
           2379*\[Kappa]S*\[Nu] + 978*\[Kappa]S*\[Nu]^2 + 796*\[Chi]A^2 - 
           4000*\[Nu]*\[Chi]A^2 + 1956*\[Nu]^2*\[Chi]A^2 + 
           4*\[Delta]*(398 - 553*\[Nu])*\[Chi]A*\[Chi]S + 796*\[Chi]S^2 - 
           1396*\[Nu]*\[Chi]S^2 + 472*\[Nu]^2*\[Chi]S^2 + 
           (\[Delta]*\[Kappa]A*(11676 - 4923*\[Nu]) + 3*\[Kappa]S*
              (3892 - 9425*\[Nu] + 1750*\[Nu]^2) + 16076*\[Chi]A^2 - 
             67025*\[Nu]*\[Chi]A^2 + 10500*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(16076 - 14269*\[Nu])*\[Chi]A*\[Chi]S + 
             16076*\[Chi]S^2 - 25817*\[Nu]*\[Chi]S^2 + 5168*\[Nu]^2*
              \[Chi]S^2)*Cos[2*\[ScriptL]])*Sin[2*\[ScriptL]])/108 + 
        (et^8*(15612006*\[Delta]*\[Kappa]A + 15612006*\[Kappa]S - 
           2525427*\[Delta]*\[Kappa]A*\[Nu] - 33749439*\[Kappa]S*\[Nu] + 
           1074234*\[Kappa]S*\[Nu]^2 + 21908750*\[Chi]A^2 - 
           84341858*\[Nu]*\[Chi]A^2 + 2148468*\[Nu]^2*\[Chi]A^2 + 
           43817500*\[Delta]*\[Chi]A*\[Chi]S - 33894992*\[Delta]*\[Nu]*
            \[Chi]A*\[Chi]S + 21908750*\[Chi]S^2 - 37188134*\[Nu]*\[Chi]S^2 + 
           4189640*\[Nu]^2*\[Chi]S^2 + 6*(3*\[Delta]*\[Kappa]A*
              (-8345356 + 3344697*\[Nu]) - 3*\[Kappa]S*(8345356 - 20035409*
                \[Nu] + 3494134*\[Nu]^2) - 34201460*\[Chi]A^2 + 
             141810149*\[Nu]*\[Chi]A^2 - 20964804*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-34201460 + 29339833*\[Nu])*\[Chi]A*\[Chi]S - 
             34201460*\[Chi]S^2 + 53675357*\[Nu]*\[Chi]S^2 - 
             10030160*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]] + 
           243*(\[Delta]*\[Kappa]A*(187414 - 61115*\[Nu]) + 
             \[Kappa]S*(187414 - 435943*\[Nu] + 54698*\[Nu]^2) + 
             2*(124287 - 500385*\[Nu] + 54698*\[Nu]^2)*\[Chi]A^2 + 
             4*\[Delta]*(124287 - 103100*\[Nu])*\[Chi]A*\[Chi]S + 
             2*(124287 - 202963*\[Nu] + 28420*\[Nu]^2)*\[Chi]S^2)*
            Cos[4*\[ScriptL]] - 139846440*\[Delta]*\[Kappa]A*
            Cos[6*\[ScriptL]] - 139846440*\[Kappa]S*Cos[6*\[ScriptL]] + 
           53580726*\[Delta]*\[Kappa]A*\[Nu]*Cos[6*\[ScriptL]] + 
           333273606*\[Kappa]S*\[Nu]*Cos[6*\[ScriptL]] - 55428996*\[Kappa]S*
            \[Nu]^2*Cos[6*\[ScriptL]] - 192811592*\[Chi]A^2*
            Cos[6*\[ScriptL]] + 796003994*\[Nu]*\[Chi]A^2*Cos[6*\[ScriptL]] - 
           110857992*\[Nu]^2*\[Chi]A^2*Cos[6*\[ScriptL]] - 
           385623184*\[Delta]*\[Chi]A*\[Chi]S*Cos[6*\[ScriptL]] + 
           328159556*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S*Cos[6*\[ScriptL]] - 
           192811592*\[Chi]S^2*Cos[6*\[ScriptL]] + 303401930*\[Nu]*\[Chi]S^2*
            Cos[6*\[ScriptL]] - 56110880*\[Nu]^2*\[Chi]S^2*Cos[6*\[ScriptL]])*
          Sin[2*\[ScriptL]])/181440) + x^(5/2)*\[Epsilon]^4*
       (-(et*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
            \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*
           Sin[\[ScriptL]])/2 - et^2*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
          \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
          \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*Cos[\[ScriptL]]*
         Sin[\[ScriptL]] - (et^8*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(9550*Cos[\[ScriptL]] + 
           6729*Cos[3*\[ScriptL]] + 4825*Cos[5*\[ScriptL]] + 
           4096*Cos[7*\[ScriptL]])*Sin[\[ScriptL]])/1260 - 
        (et^3*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(13*Sin[\[ScriptL]] + 9*Sin[3*\[ScriptL]]))/
         16 - (2*et^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(Sin[2*\[ScriptL]] + Sin[4*\[ScriptL]]))/3 - 
        (et^5*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(874*Sin[\[ScriptL]] + 
           459*Sin[3*\[ScriptL]] + 625*Sin[5*\[ScriptL]]))/768 - 
        (et^6*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(215*Sin[2*\[ScriptL]] + 
           128*Sin[4*\[ScriptL]] + 243*Sin[6*\[ScriptL]]))/240 - 
        (et^7*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(134845*Sin[\[ScriptL]] + 
           73629*Sin[3*\[ScriptL]] + 40625*Sin[5*\[ScriptL]] + 
           117649*Sin[7*\[ScriptL]]))/92160)) + 
    SO*(x^3*\[Epsilon]^5*(-(et*(\[Delta]*(-36 + 37*\[Nu])*\[Chi]A + 
            (-36 + 73*\[Nu] - 14*\[Nu]^2)*\[Chi]S)*Sin[\[ScriptL]])/6 - 
        (et^3*(\[Delta]*(-1356 + 673*\[Nu])*\[Chi]A + 
           (-1356 + 1957*\[Nu] - 242*\[Nu]^2)*\[Chi]S + 
           9*(\[Delta]*(-236 + 85*\[Nu])*\[Chi]A + (-236 + 321*\[Nu] - 26*
                \[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]])*Sin[\[ScriptL]])/24 - 
        (et^5*(-181452*\[Delta]*\[Chi]A + 77695*\[Delta]*\[Nu]*\[Chi]A - 
           181452*\[Chi]S + 247915*\[Nu]*\[Chi]S - 27674*\[Nu]^2*\[Chi]S + 
           4*(\[Delta]*(-72300 + 26449*\[Nu])*\[Chi]A + 
             (-72300 + 96805*\[Nu] - 8882*\[Nu]^2)*\[Chi]S)*
            Cos[2*\[ScriptL]] + 5*(\[Delta]*(-65940 + 18497*\[Nu])*\[Chi]A + 
             (-65940 + 84437*\[Nu] - 4438*\[Nu]^2)*\[Chi]S)*
            Cos[4*\[ScriptL]])*Sin[\[ScriptL]])/1152 + 
        (et^6*Cos[\[ScriptL]]*(186072*\[Delta]*\[Chi]A - 56389*\[Delta]*\[Nu]*
            \[Chi]A + 186072*\[Chi]S - 238501*\[Nu]*\[Chi]S + 
           14276*\[Nu]^2*\[Chi]S + 8*(\[Delta]*(-10632 + 239*\[Nu])*\[Chi]A + 
             (-10632 + 11591*\[Nu] + 1244*\[Nu]^2)*\[Chi]S)*
            Cos[2*\[ScriptL]] - 27*(\[Delta]*(-12512 + 3279*\[Nu])*\[Chi]A + 
             (-12512 + 15791*\[Nu] - 696*\[Nu]^2)*\[Chi]S)*Cos[4*\[ScriptL]])*
          Sin[\[ScriptL]])/360 - (et^7*(-43042392*\[Delta]*\[Chi]A + 
           17307386*\[Delta]*\[Nu]*\[Chi]A - 43042392*\[Chi]S + 
           57509378*\[Nu]*\[Chi]S - 6156292*\[Nu]^2*\[Chi]S + 
           3*(\[Delta]*(-23687388 + 8537609*\[Nu])*\[Chi]A + 
             (-23687388 + 31288997*\[Nu] - 2913058*\[Nu]^2)*\[Chi]S)*
            Cos[2*\[ScriptL]] + 6*(\[Delta]*(-10670652 + 3547081*\[Nu])*
              \[Chi]A + (-10670652 + 13992733*\[Nu] - 1192442*\[Nu]^2)*
              \[Chi]S)*Cos[4*\[ScriptL]] - 102861612*\[Delta]*\[Chi]A*
            Cos[6*\[ScriptL]] + 25594261*\[Delta]*\[Nu]*\[Chi]A*
            Cos[6*\[ScriptL]] - 102861612*\[Chi]S*Cos[6*\[ScriptL]] + 
           128455873*\[Nu]*\[Chi]S*Cos[6*\[ScriptL]] - 4788602*\[Nu]^2*
            \[Chi]S*Cos[6*\[ScriptL]])*Sin[\[ScriptL]])/138240 - 
        (et^2*(\[Delta]*(-60 + 29*\[Nu])*\[Chi]A + 
           (-60 + 89*\[Nu] - 10*\[Nu]^2)*\[Chi]S)*Sin[2*\[ScriptL]])/3 - 
        (et^4*(2*\[Delta]*(-18 + 83*\[Nu])*\[Chi]A - 
           4*(9 - 37*\[Nu] + 23*\[Nu]^2)*\[Chi]S + 
           (\[Delta]*(-2988 + 923*\[Nu])*\[Chi]A + (-2988 + 3911*\[Nu] - 250*
                \[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]])*Sin[2*\[ScriptL]])/18 - 
        (et^8*(7603908*\[Delta]*\[Chi]A - 532949*\[Delta]*\[Nu]*\[Chi]A + 
           7603908*\[Chi]S - 8577605*\[Nu]*\[Chi]S - 635006*\[Nu]^2*\[Chi]S + 
           48*(\[Delta]*(-777288 + 194699*\[Nu])*\[Chi]A + 
             (-777288 + 965435*\[Nu] - 33724*\[Nu]^2)*\[Chi]S)*
            Cos[2*\[ScriptL]] - 243*(\[Delta]*(-72236 + 10517*\[Nu])*
              \[Chi]A + (-72236 + 84265*\[Nu] + 2750*\[Nu]^2)*\[Chi]S)*
            Cos[4*\[ScriptL]] - 34855512*\[Delta]*\[Chi]A*Cos[6*\[ScriptL]] + 
           8318078*\[Delta]*\[Nu]*\[Chi]A*Cos[6*\[ScriptL]] - 
           34855512*\[Chi]S*Cos[6*\[ScriptL]] + 43173590*\[Nu]*\[Chi]S*
            Cos[6*\[ScriptL]] - 1363252*\[Nu]^2*\[Chi]S*Cos[6*\[ScriptL]])*
          Sin[2*\[ScriptL]])/30240) + x^2*\[Epsilon]^3*
       ((2*et*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*Sin[\[ScriptL]])/
         3 + 3*et^3*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*
         Cos[\[ScriptL]]^2*Sin[\[ScriptL]] + 
        (et^6*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*
          (137*Cos[\[ScriptL]] + 97*Cos[3*\[ScriptL]] + 81*Cos[5*\[ScriptL]])*
          Sin[\[ScriptL]])/30 + (2*et^2*(2*\[Delta]*\[Chi]A - 
           (-2 + \[Nu])*\[Chi]S)*Sin[2*\[ScriptL]])/3 + 
        (et^4*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*
          (5*Sin[2*\[ScriptL]] + 8*Sin[4*\[ScriptL]]))/9 + 
        (et^5*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*
          (514*Sin[\[ScriptL]] + 243*Sin[3*\[ScriptL]] + 
           625*Sin[5*\[ScriptL]]))/576 + 
        (et^7*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*
          (70165*Sin[\[ScriptL]] + 39609*Sin[3*\[ScriptL]] + 
           3125*Sin[5*\[ScriptL]] + 117649*Sin[7*\[ScriptL]]))/69120 + 
        (et^8*(2*\[Delta]*\[Chi]A - (-2 + \[Nu])*\[Chi]S)*
          (5614*Sin[2*\[ScriptL]] + 4088*Sin[4*\[ScriptL]] - 
           2187*Sin[6*\[ScriptL]] + 16384*Sin[8*\[ScriptL]]))/7560))


(* ::Subsection::Closed:: *)
(*d\[Phi]/dt*)


\[Phi]d\[LetterSpace]xel = x^(5/2)*\[Epsilon]^2*
     (-(et*(-4 + \[Nu])*Cos[\[ScriptL]]) - 2*et^2*(-4 + \[Nu])*
       Cos[2*\[ScriptL]] - (et^3*(-4 + \[Nu])*(Cos[\[ScriptL]] + 
         27*Cos[3*\[ScriptL]]))/8 + (et^4*(-4 + \[Nu])*(Cos[2*\[ScriptL]] - 
         16*Cos[4*\[ScriptL]]))/3 - (et^5*(-4 + \[Nu])*(82*Cos[\[ScriptL]] - 
         567*Cos[3*\[ScriptL]] + 3125*Cos[5*\[ScriptL]]))/384 - 
      (et^6*(-4 + \[Nu])*(20*Cos[2*\[ScriptL]] - 224*Cos[4*\[ScriptL]] + 
         729*Cos[6*\[ScriptL]]))/60 - 
      (et^7*(-4 + \[Nu])*(8485*Cos[\[ScriptL]] + 31347*Cos[3*\[ScriptL]] - 
         359375*Cos[5*\[ScriptL]] + 823543*Cos[7*\[ScriptL]]))/46080 - 
      (et^8*(-4 + \[Nu])*(574*Cos[2*\[ScriptL]] + 4144*Cos[4*\[ScriptL]] - 
         37179*Cos[6*\[ScriptL]] + 65536*Cos[8*\[ScriptL]]))/2520) + 
    x^(3/2)*(1 + 2*et*Cos[\[ScriptL]] + (5*et^2*Cos[2*\[ScriptL]])/2 + 
      (et^3*(-Cos[\[ScriptL]] + 13*Cos[3*\[ScriptL]]))/4 + 
      (et^4*(-22*Cos[2*\[ScriptL]] + 103*Cos[4*\[ScriptL]]))/24 + 
      (et^5*(10*Cos[\[ScriptL]] - 387*Cos[3*\[ScriptL]] + 
         1097*Cos[5*\[ScriptL]]))/192 + 
      (et^6*(85*Cos[2*\[ScriptL]] - 1804*Cos[4*\[ScriptL]] + 
         3669*Cos[6*\[ScriptL]]))/480 + 
      (et^7*(107*Cos[\[ScriptL]] + 2565*Cos[3*\[ScriptL]] - 
         29785*Cos[5*\[ScriptL]] + 47273*Cos[7*\[ScriptL]]))/4608 + 
      (et^8*(602*Cos[2*\[ScriptL]] + 57722*Cos[4*\[ScriptL]] - 
         427302*Cos[6*\[ScriptL]] + 556403*Cos[8*\[ScriptL]]))/40320) + 
    x^(7/2)*\[Epsilon]^4*((et*(156 - 31*\[Nu] + \[Nu]^2)*Cos[\[ScriptL]])/
       12 + (et^2*(138 - 22*\[Nu] + \[Nu]^2)*Cos[2*\[ScriptL]])/6 + 
      (et^3*Cos[\[ScriptL]]*(-816 + 49*\[Nu] + 29*\[Nu]^2 + 
         9*(352 - 51*\[Nu] + \[Nu]^2)*Cos[2*\[ScriptL]]))/48 + 
      (et^4*(4*(969 - 158*\[Nu] + 62*\[Nu]^2)*Cos[2*\[ScriptL]] + 
         (6141 - 922*\[Nu] - 44*\[Nu]^2)*Cos[4*\[ScriptL]]))/144 - 
      (et^5*Cos[\[ScriptL]]*(-64116 + 25341*\[Nu] + 9669*\[Nu]^2 - 
         4*(-12336 + 7085*\[Nu] + 5557*\[Nu]^2)*Cos[2*\[ScriptL]] + 
         (-232620 + 42191*\[Nu] + 6055*\[Nu]^2)*Cos[4*\[ScriptL]]))/2304 - 
      (et^7*Cos[\[ScriptL]]*(-8641680 - 11398330*\[Nu] - 5512610*\[Nu]^2 + 
         15*(463536 + 1717811*\[Nu] + 708247*\[Nu]^2)*Cos[2*\[ScriptL]] - 
         6*(3933576 + 3279881*\[Nu] + 1827277*\[Nu]^2)*Cos[4*\[ScriptL]] - 
         26169744*Cos[6*\[ScriptL]] + 12586771*\[Nu]*Cos[6*\[ScriptL]] + 
         4028087*\[Nu]^2*Cos[6*\[ScriptL]]))/276480 + 
      (et^6*(5*(29271 - 8510*\[Nu] + 1028*\[Nu]^2)*Cos[2*\[ScriptL]] + 
         4*(84387 + 2554*\[Nu] + 9692*\[Nu]^2)*Cos[4*\[ScriptL]] - 
         27*(-11451 + 2998*\[Nu] + 724*\[Nu]^2)*Cos[6*\[ScriptL]]))/5760 + 
      (et^8*(14*(479487 - 130750*\[Nu] + 19876*\[Nu]^2)*Cos[2*\[ScriptL]] - 
         112*(-61521 + 39632*\[Nu] + 2053*\[Nu]^2)*Cos[4*\[ScriptL]] + 
         35241318*Cos[6*\[ScriptL]] + 8547444*\[Nu]*Cos[6*\[ScriptL]] + 
         5544936*\[Nu]^2*Cos[6*\[ScriptL]] + 5647047*Cos[8*\[ScriptL]] - 
         9193670*\[Nu]*Cos[8*\[ScriptL]] - 3423544*\[Nu]^2*
          Cos[8*\[ScriptL]]))/241920) + 
    SO*(x^3*\[Epsilon]^3*(-2*et*(\[Delta]*\[Chi]A + \[Chi]S)*
         Cos[\[ScriptL]] - 4*et^2*(\[Delta]*\[Chi]A + \[Chi]S)*
         Cos[2*\[ScriptL]] - (et^3*(\[Delta]*\[Chi]A + \[Chi]S)*
          (5*Cos[\[ScriptL]] + 27*Cos[3*\[ScriptL]]))/4 - 
        (4*et^4*(\[Delta]*\[Chi]A + \[Chi]S)*(Cos[2*\[ScriptL]] + 
           8*Cos[4*\[ScriptL]]))/3 - (et^5*(\[Delta]*\[Chi]A + \[Chi]S)*
          (250*Cos[\[ScriptL]] + 81*Cos[3*\[ScriptL]] + 
           3125*Cos[5*\[ScriptL]]))/192 - (et^6*(\[Delta]*\[Chi]A + \[Chi]S)*
          (55*Cos[2*\[ScriptL]] - 64*Cos[4*\[ScriptL]] + 
           729*Cos[6*\[ScriptL]]))/30 - (et^7*(\[Delta]*\[Chi]A + \[Chi]S)*
          (29965*Cos[\[ScriptL]] + 55647*Cos[3*\[ScriptL]] - 
           171875*Cos[5*\[ScriptL]] + 823543*Cos[7*\[ScriptL]]))/23040 - 
        (et^8*(\[Delta]*\[Chi]A + \[Chi]S)*(1127*Cos[2*\[ScriptL]] + 
           2240*Cos[4*\[ScriptL]] - 10935*Cos[6*\[ScriptL]] + 
           32768*Cos[8*\[ScriptL]]))/630) + x^4*\[Epsilon]^5*
       ((2*et*(\[Delta]*(-38 + 5*\[Nu])*\[Chi]A + 
           2*(-19 + 18*\[Nu] + \[Nu]^2)*\[Chi]S)*Cos[\[ScriptL]])/3 + 
        (et^2*(\[Delta]*(-220 + 7*\[Nu])*\[Chi]A + 
           (-220 + 117*\[Nu] + 40*\[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]])/6 - 
        (et^3*Cos[\[ScriptL]]*(\[Delta]*(76 - 109*\[Nu])*\[Chi]A + 
           2*(38 - 207*\[Nu] + 52*\[Nu]^2)*\[Chi]S + 
           3*(\[Delta]*(156 + 25*\[Nu])*\[Chi]A + 4*(39 + 10*\[Nu] - 17*
                \[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]]))/6 + 
        (et^4*((\[Delta]*(-3188 + 737*\[Nu])*\[Chi]A + 
             (-3188 + 3531*\[Nu] - 4*\[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]] - 
           (\[Delta]*(958 + 785*\[Nu])*\[Chi]A + (958 + 2829*\[Nu] - 1261*
                \[Nu]^2)*\[Chi]S)*Cos[4*\[ScriptL]]))/36 + 
        (et^5*Cos[\[ScriptL]]*(3*\[Delta]*(7454 - 6395*\[Nu])*\[Chi]A + 
           6*(3727 - 13674*\[Nu] + 3209*\[Nu]^2)*\[Chi]S + 
           4*(\[Delta]*(-21386 + 12089*\[Nu])*\[Chi]A - 
             2*(10693 - 26475*\[Nu] + 4862*\[Nu]^2)*\[Chi]S)*
            Cos[2*\[ScriptL]] + (\[Delta]*(6542 - 28979*\[Nu])*\[Chi]A + 
             2*(3271 - 60558*\[Nu] + 18701*\[Nu]^2)*\[Chi]S)*
            Cos[4*\[ScriptL]]))/288 - (et^7*Cos[\[ScriptL]]*
          (18940720*\[Delta]*\[Chi]A - 9483970*\[Delta]*\[Nu]*\[Chi]A + 
           18940720*\[Chi]S - 44126220*\[Nu]*\[Chi]S + 7893560*\[Nu]^2*
            \[Chi]S + 15*(\[Delta]*(-2110640 + 1159931*\[Nu])*\[Chi]A - 
             28*(75380 - 192618*\[Nu] + 37435*\[Nu]^2)*\[Chi]S)*
            Cos[2*\[ScriptL]] - 6*(\[Delta]*(-6892280 + 3274853*\[Nu])*
              \[Chi]A - 2*(3446140 - 7615179*\[Nu] + 1298030*\[Nu]^2)*
              \[Chi]S)*Cos[4*\[ScriptL]] - 16953520*\[Delta]*\[Chi]A*
            Cos[6*\[ScriptL]] + 12425323*\[Delta]*\[Nu]*\[Chi]A*
            Cos[6*\[ScriptL]] - 16953520*\[Chi]S*Cos[6*\[ScriptL]] + 
           57580968*\[Nu]*\[Chi]S*Cos[6*\[ScriptL]] - 13034900*\[Nu]^2*
            \[Chi]S*Cos[6*\[ScriptL]]))/34560 + 
        (et^6*(5*(\[Delta]*(-31562 + 7517*\[Nu])*\[Chi]A + 
             (-31562 + 35433*\[Nu] + 71*\[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]] + 
           16*(2*\[Delta]*(-9697 + 2611*\[Nu])*\[Chi]A + 
             (-19394 + 24264*\[Nu] - 1099*\[Nu]^2)*\[Chi]S)*
            Cos[4*\[ScriptL]] - 27*(\[Delta]*(-4942 + 5291*\[Nu])*\[Chi]A + 
             (-4942 + 23607*\[Nu] - 6027*\[Nu]^2)*\[Chi]S)*
            Cos[6*\[ScriptL]]))/1440 + 
        (et^8*(7*(\[Delta]*(-1181546 + 291701*\[Nu])*\[Chi]A + 
             (-1181546 + 1354209*\[Nu] - 2209*\[Nu]^2)*\[Chi]S)*
            Cos[2*\[ScriptL]] + 14*(\[Delta]*(-675082 + 148687*\[Nu])*
              \[Chi]A + (-675082 + 679719*\[Nu] + 28579*\[Nu]^2)*\[Chi]S)*
            Cos[4*\[ScriptL]] - 36320562*\[Delta]*\[Chi]A*Cos[6*\[ScriptL]] + 
           11589237*\[Delta]*\[Nu]*\[Chi]A*Cos[6*\[ScriptL]] - 
           36320562*\[Chi]S*Cos[6*\[ScriptL]] + 55002969*\[Nu]*\[Chi]S*
            Cos[6*\[ScriptL]] - 5428701*\[Nu]^2*\[Chi]S*Cos[6*\[ScriptL]] + 
           30992722*\[Delta]*\[Chi]A*Cos[8*\[ScriptL]] - 18680692*\[Delta]*
            \[Nu]*\[Chi]A*Cos[8*\[ScriptL]] + 30992722*\[Chi]S*
            Cos[8*\[ScriptL]] - 88793418*\[Nu]*\[Chi]S*Cos[8*\[ScriptL]] + 
           18472823*\[Nu]^2*\[Chi]S*Cos[8*\[ScriptL]]))/60480)) + 
    SO^2*(x^(7/2)*\[Epsilon]^4*(et*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
          \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
          \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*Cos[\[ScriptL]] + 
        2*et^2*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
          4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
            2*\[Chi]A*\[Chi]S))*Cos[2*\[ScriptL]] + 
        (9*et^3*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(Cos[\[ScriptL]] + 3*Cos[3*\[ScriptL]]))/8 + 
        (et^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(5*Cos[2*\[ScriptL]] + 
           16*Cos[4*\[ScriptL]]))/3 + (et^5*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
           \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
           \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*(514*Cos[\[ScriptL]] + 
           729*Cos[3*\[ScriptL]] + 3125*Cos[5*\[ScriptL]]))/384 + 
        (et^6*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(40*Cos[2*\[ScriptL]] + 
           32*Cos[4*\[ScriptL]] + 243*Cos[6*\[ScriptL]]))/20 + 
        (et^7*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(70165*Cos[\[ScriptL]] + 
           118827*Cos[3*\[ScriptL]] + 15625*Cos[5*\[ScriptL]] + 
           823543*Cos[7*\[ScriptL]]))/46080 + 
        (et^8*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(5614*Cos[2*\[ScriptL]] + 
           8176*Cos[4*\[ScriptL]] - 6561*Cos[6*\[ScriptL]] + 
           65536*Cos[8*\[ScriptL]]))/2520) + x^(9/2)*\[Epsilon]^6*
       (-(et*(\[Delta]*\[Kappa]A*(-108 + 55*\[Nu]) + \[Kappa]S*
             (-108 + 271*\[Nu] - 14*\[Nu]^2) - 204*\[Chi]A^2 + 
            829*\[Nu]*\[Chi]A^2 - 28*\[Nu]^2*\[Chi]A^2 + 
            2*\[Delta]*(-204 + 133*\[Nu])*\[Chi]A*\[Chi]S - 204*\[Chi]S^2 + 
            253*\[Nu]*\[Chi]S^2 - 48*\[Nu]^2*\[Chi]S^2)*Cos[\[ScriptL]])/6 - 
        (et^2*(\[Delta]*\[Kappa]A*(-162 + 95*\[Nu]) + 
           \[Kappa]S*(-162 + 419*\[Nu] + 2*\[Nu]^2) - 330*\[Chi]A^2 + 
           1322*\[Nu]*\[Chi]A^2 + 4*\[Nu]^2*\[Chi]A^2 + 
           4*\[Delta]*(-165 + 76*\[Nu])*\[Chi]A*\[Chi]S - 330*\[Chi]S^2 + 
           302*\[Nu]*\[Chi]S^2 - 36*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]])/6 - 
        (et^3*Cos[\[ScriptL]]*(-136*\[Kappa]S + 279*\[Kappa]S*\[Nu] - 
           182*\[Kappa]S*\[Nu]^2 + \[Delta]*\[Kappa]A*(-136 + 7*\[Nu]) - 
           88*\[Chi]A^2 + 479*\[Nu]*\[Chi]A^2 - 364*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(-88 + 387*\[Nu])*\[Chi]A*\[Chi]S - 88*\[Chi]S^2 + 
           647*\[Nu]*\[Chi]S^2 - 224*\[Nu]^2*\[Chi]S^2 + 
           (\[Delta]*\[Kappa]A*(-488 + 359*\[Nu]) + \[Kappa]S*
              (-488 + 1335*\[Nu] + 138*\[Nu]^2) - 1144*\[Chi]A^2 + 
             4487*\[Nu]*\[Chi]A^2 + 276*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(-1144 + 219*\[Nu])*\[Chi]A*\[Chi]S - 
             1144*\[Chi]S^2 + 527*\[Nu]*\[Chi]S^2 + 64*\[Nu]^2*\[Chi]S^2)*
            Cos[2*\[ScriptL]]))/8 + 
        (et^4*(4*(\[Delta]*\[Kappa]A*(2910 - 1453*\[Nu]) + 
             \[Kappa]S*(2910 - 7273*\[Nu] + 878*\[Nu]^2) + 
             2*(2499 - 10298*\[Nu] + 878*\[Nu]^2)*\[Chi]A^2 + 
             4*\[Delta]*(2499 - 1781*\[Nu])*\[Chi]A*\[Chi]S + 
             2*(2499 - 3260*\[Nu] + 624*\[Nu]^2)*\[Chi]S^2)*
            Cos[2*\[ScriptL]] + (\[Delta]*\[Kappa]A*(3570 - 4139*\[Nu]) + 
             \[Kappa]S*(3570 - 11279*\[Nu] - 3734*\[Nu]^2) + 
             11562*\[Chi]A^2 - 43676*\[Nu]*\[Chi]A^2 - 7468*\[Nu]^2*
              \[Chi]A^2 + 4*\[Delta]*(5781 + 1541*\[Nu])*\[Chi]A*\[Chi]S + 
             11562*\[Chi]S^2 + 3592*\[Nu]*\[Chi]S^2 - 3768*\[Nu]^2*\[Chi]S^2)*
            Cos[4*\[ScriptL]]))/144 - (et^5*Cos[\[ScriptL]]*
          (3*(11*\[Delta]*\[Kappa]A*(1356 + 445*\[Nu]) + 
             \[Kappa]S*(14916 - 24937*\[Nu] + 27074*\[Nu]^2) - 
             2652*\[Chi]A^2 - 7483*\[Nu]*\[Chi]A^2 + 54148*\[Nu]^2*
              \[Chi]A^2 - 34*\[Delta]*(156 + 3011*\[Nu])*\[Chi]A*\[Chi]S - 
             2652*\[Chi]S^2 - 84283*\[Nu]*\[Chi]S^2 + 27792*\[Nu]^2*
              \[Chi]S^2) + 4*(\[Delta]*\[Kappa]A*(-68016 + 15133*\[Nu]) + 
             \[Kappa]S*(-68016 + 151165*\[Nu] - 55778*\[Nu]^2) - 
             74304*\[Chi]A^2 + 335269*\[Nu]*\[Chi]A^2 - 111556*\[Nu]^2*
              \[Chi]A^2 + 2*\[Delta]*(-74304 + 137473*\[Nu])*\[Chi]A*
              \[Chi]S - 74304*\[Chi]S^2 + 236893*\[Nu]*\[Chi]S^2 - 
             65856*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]] + 
           (\[Delta]*\[Kappa]A*(-6252 + 77735*\[Nu]) + \[Kappa]S*
              (-6252 + 90239*\[Nu] + 133202*\[Nu]^2) - 170316*\[Chi]A^2 + 
             588101*\[Nu]*\[Chi]A^2 + 266404*\[Nu]^2*\[Chi]A^2 - 
             2*\[Delta]*(170316 + 219331*\[Nu])*\[Chi]A*\[Chi]S - 
             170316*\[Chi]S^2 - 345499*\[Nu]*\[Chi]S^2 + 155856*\[Nu]^2*
              \[Chi]S^2)*Cos[4*\[ScriptL]]))/1152 - 
        (et^7*Cos[\[ScriptL]]*(-49742400*\[Delta]*\[Kappa]A - 
           49742400*\[Kappa]S + 8454650*\[Delta]*\[Kappa]A*\[Nu] + 
           107939450*\[Kappa]S*\[Nu] - 41695300*\[Kappa]S*\[Nu]^2 - 
           50490720*\[Chi]A^2 + 230399330*\[Nu]*\[Chi]A^2 - 
           83390600*\[Nu]^2*\[Chi]A^2 - 100981440*\[Delta]*\[Chi]A*\[Chi]S + 
           208684180*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S - 50490720*\[Chi]S^2 + 
           180247730*\[Nu]*\[Chi]S^2 - 49109760*\[Nu]^2*\[Chi]S^2 + 
           15*(\[Delta]*\[Kappa]A*(4516200 - 81155*\[Nu]) + 
             5*\[Kappa]S*(903240 - 1822711*\[Nu] + 963446*\[Nu]^2) + 
             3125304*\[Chi]A^2 - 15762011*\[Nu]*\[Chi]A^2 + 
             9634460*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(3125304 - 11073263*
                \[Nu])*\[Chi]A*\[Chi]S + 3125304*\[Chi]S^2 - 
             18885731*\[Nu]*\[Chi]S^2 + 5401152*\[Nu]^2*\[Chi]S^2)*
            Cos[2*\[ScriptL]] - 6*(\[Delta]*\[Kappa]A*(19181952 - 4210025*
                \[Nu]) + \[Kappa]S*(19181952 - 42573929*\[Nu] + 14617498*
                \[Nu]^2) + 21295440*\[Chi]A^2 - 95083037*\[Nu]*\[Chi]A^2 + 
             29234996*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(21295440 - 37339937*
                \[Nu])*\[Chi]A*\[Chi]S + 21295440*\[Chi]S^2 - 
             64778597*\[Nu]*\[Chi]S^2 + 16949376*\[Nu]^2*\[Chi]S^2)*
            Cos[4*\[ScriptL]] + 40032552*\[Delta]*\[Kappa]A*
            Cos[6*\[ScriptL]] + 40032552*\[Kappa]S*Cos[6*\[ScriptL]] + 
           8434525*\[Delta]*\[Kappa]A*\[Nu]*Cos[6*\[ScriptL]] - 
           71630579*\[Kappa]S*\[Nu]*Cos[6*\[ScriptL]] + 56528398*\[Kappa]S*
            \[Nu]^2*Cos[6*\[ScriptL]] + 12186360*\[Chi]A^2*
            Cos[6*\[ScriptL]] - 88442267*\[Nu]*\[Chi]A^2*Cos[6*\[ScriptL]] + 
           113056796*\[Nu]^2*\[Chi]A^2*Cos[6*\[ScriptL]] + 
           24372720*\[Delta]*\[Chi]A*\[Chi]S*Cos[6*\[ScriptL]] - 
           259359454*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S*Cos[6*\[ScriptL]] + 
           12186360*\[Chi]S^2*Cos[6*\[ScriptL]] - 219662627*\[Nu]*\[Chi]S^2*
            Cos[6*\[ScriptL]] + 71516736*\[Nu]^2*\[Chi]S^2*
            Cos[6*\[ScriptL]]))/138240 + 
        (et^6*(5*(\[Delta]*\[Kappa]A*(48738 - 24371*\[Nu]) + 
             \[Kappa]S*(48738 - 121847*\[Nu] + 16058*\[Nu]^2) + 
             83290*\[Chi]A^2 - 344332*\[Nu]*\[Chi]A^2 + 32116*\[Nu]^2*
              \[Chi]A^2 + 4*\[Delta]*(41645 - 30907*\[Nu])*\[Chi]A*\[Chi]S + 
             83290*\[Chi]S^2 - 112456*\[Nu]*\[Chi]S^2 + 23144*\[Nu]^2*
              \[Chi]S^2)*Cos[2*\[ScriptL]] + 
           4*(\[Delta]*\[Kappa]A*(85686 - 41105*\[Nu]) + 
             \[Kappa]S*(85686 - 212477*\[Nu] + 31070*\[Nu]^2) + 
             2*(69359 - 287354*\[Nu] + 31070*\[Nu]^2)*\[Chi]A^2 + 
             4*\[Delta]*(69359 - 51409*\[Nu])*\[Chi]A*\[Chi]S + 
             2*(69359 - 92900*\[Nu] + 16204*\[Nu]^2)*\[Chi]S^2)*
            Cos[4*\[ScriptL]] - 9*(89*\[Delta]*\[Kappa]A*(114 + 85*\[Nu]) + 
             \[Kappa]S*(10146 - 12727*\[Nu] + 24090*\[Nu]^2) - 
             8422*\[Chi]A^2 + 16772*\[Nu]*\[Chi]A^2 + 48180*\[Nu]^2*
              \[Chi]A^2 - 4*\[Delta]*(4211 + 24739*\[Nu])*\[Chi]A*\[Chi]S - 
             8422*\[Chi]S^2 - 82040*\[Nu]*\[Chi]S^2 + 29768*\[Nu]^2*
              \[Chi]S^2)*Cos[6*\[ScriptL]]))/1920 + 
        (et^8*(-112*(-387726*\[Kappa]S + 969109*\[Kappa]S*\[Nu] - 
             135134*\[Kappa]S*\[Nu]^2 + \[Delta]*\[Kappa]A*(-387726 + 193657*
                \[Nu]) - 658518*\[Chi]A^2 + 2728252*\[Nu]*\[Chi]A^2 - 
             270268*\[Nu]^2*\[Chi]A^2 + 28*\[Delta]*(-47037 + 35819*\[Nu])*
              \[Chi]A*\[Chi]S - 658518*\[Chi]S^2 + 908752*\[Nu]*\[Chi]S^2 - 
             193272*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]] + 
           56*(\[Delta]*\[Kappa]A*(916086 - 475481*\[Nu]) + 
             \[Kappa]S*(916086 - 2307653*\[Nu] + 291046*\[Nu]^2) + 
             2*(787839 - 3248764*\[Nu] + 291046*\[Nu]^2)*\[Chi]A^2 + 
             4*\[Delta]*(787839 - 541999*\[Nu])*\[Chi]A*\[Chi]S + 
             2*(787839 - 986590*\[Nu] + 186360*\[Nu]^2)*\[Chi]S^2)*
            Cos[4*\[ScriptL]] + 104022144*\[Delta]*\[Kappa]A*
            Cos[6*\[ScriptL]] + 104022144*\[Kappa]S*Cos[6*\[ScriptL]] - 
           42669504*\[Delta]*\[Kappa]A*\[Nu]*Cos[6*\[ScriptL]] - 
           250713792*\[Kappa]S*\[Nu]*Cos[6*\[ScriptL]] + 49019904*\[Kappa]S*
            \[Nu]^2*Cos[6*\[ScriptL]] + 151466112*\[Chi]A^2*
            Cos[6*\[ScriptL]] - 636934752*\[Nu]*\[Chi]A^2*Cos[6*\[ScriptL]] + 
           98039808*\[Nu]^2*\[Chi]A^2*Cos[6*\[ScriptL]] + 302932224*\[Delta]*
            \[Chi]A*\[Chi]S*Cos[6*\[ScriptL]] - 279422784*\[Delta]*\[Nu]*
            \[Chi]A*\[Chi]S*Cos[6*\[ScriptL]] + 151466112*\[Chi]S^2*
            Cos[6*\[ScriptL]] - 248352480*\[Nu]*\[Chi]S^2*Cos[6*\[ScriptL]] + 
           46604160*\[Nu]^2*\[Chi]S^2*Cos[6*\[ScriptL]] - 76973682*\[Delta]*
            \[Kappa]A*Cos[8*\[ScriptL]] - 76973682*\[Kappa]S*
            Cos[8*\[ScriptL]] - 3208741*\[Delta]*\[Kappa]A*\[Nu]*
            Cos[8*\[ScriptL]] + 150738623*\[Kappa]S*\[Nu]*Cos[8*\[ScriptL]] - 
           85240858*\[Kappa]S*\[Nu]^2*Cos[8*\[ScriptL]] - 50818026*\[Chi]A^2*
            Cos[8*\[ScriptL]] + 263038964*\[Nu]*\[Chi]A^2*Cos[8*\[ScriptL]] - 
           170481716*\[Nu]^2*\[Chi]A^2*Cos[8*\[ScriptL]] - 
           101636052*\[Delta]*\[Chi]A*\[Chi]S*Cos[8*\[ScriptL]] + 
           417808444*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S*Cos[8*\[ScriptL]] - 
           50818026*\[Chi]S^2*Cos[8*\[ScriptL]] + 358041584*\[Nu]*\[Chi]S^2*
            Cos[8*\[ScriptL]] - 108922584*\[Nu]^2*\[Chi]S^2*
            Cos[8*\[ScriptL]]))/241920))


(* ::Subsection::Closed:: *)
(*\[Phi]*)


\[Phi]\[LetterSpace]xel = \[ScriptL] + 2*et*Sin[\[ScriptL]] + 
    (5*et^2*Cos[\[ScriptL]/2]*(-Sin[\[ScriptL]/2] + Sin[(3*\[ScriptL])/2]))/
     2 + (et^3*Cos[\[ScriptL]/2]*(10*Sin[\[ScriptL]/2] + 
       13*(-Sin[(3*\[ScriptL])/2] + Sin[(5*\[ScriptL])/2])))/6 + 
    (et^4*Cos[\[ScriptL]/2]*(-59*Sin[\[ScriptL]/2] + 
       59*Sin[(3*\[ScriptL])/2] + 103*(-Sin[(5*\[ScriptL])/2] + 
         Sin[(7*\[ScriptL])/2])))/48 + 
    (et^5*Cos[\[ScriptL]/2]*(502*Sin[\[ScriptL]/2] - 
       452*Sin[(3*\[ScriptL])/2] + 452*Sin[(5*\[ScriptL])/2] - 
       1097*Sin[(7*\[ScriptL])/2] + 1097*Sin[(9*\[ScriptL])/2]))/480 + 
    (et^6*Cos[\[ScriptL]/2]*(-406*Sin[\[ScriptL]/2] + 
       406*Sin[(3*\[ScriptL])/2] - 321*Sin[(5*\[ScriptL])/2] + 
       321*Sin[(7*\[ScriptL])/2] - 1223*Sin[(9*\[ScriptL])/2] + 
       1223*Sin[(11*\[ScriptL])/2]))/480 + 
    (et^7*Cos[\[ScriptL]/2]*(12308*Sin[\[ScriptL]/2] - 
       11559*Sin[(3*\[ScriptL])/2] + 11559*Sin[(5*\[ScriptL])/2] - 
       5574*Sin[(7*\[ScriptL])/2] + 5574*Sin[(9*\[ScriptL])/2] - 
       47273*Sin[(11*\[ScriptL])/2] + 47273*Sin[(13*\[ScriptL])/2]))/16128 + 
    (et^8*Cos[\[ScriptL]/2]*(-104519*Sin[\[ScriptL]/2] + 
       104519*Sin[(3*\[ScriptL])/2] - 102111*Sin[(5*\[ScriptL])/2] + 
       102111*Sin[(7*\[ScriptL])/2] + 13333*Sin[(9*\[ScriptL])/2] - 
       13333*Sin[(11*\[ScriptL])/2] - 556403*Sin[(13*\[ScriptL])/2] + 
       556403*Sin[(15*\[ScriptL])/2]))/161280 + 
    x*\[Epsilon]^2*(3*\[ScriptL] - et*(-10 + \[Nu])*Sin[\[ScriptL]] - 
      (et^3*(-54 + 5*\[Nu] + (-62 + 9*\[Nu])*Cos[2*\[ScriptL]])*
        Sin[\[ScriptL]])/4 - (et^5*(-14846 + 1295*\[Nu] + 
         4*(-4418 + 545*\[Nu])*Cos[2*\[ScriptL]] + (-19082 + 3125*\[Nu])*
          Cos[4*\[ScriptL]])*Sin[\[ScriptL]])/960 - 
      (et^7*(-2714420 + 226478*\[Nu] + 9*(-380550 + 43729*\[Nu])*
          Cos[2*\[ScriptL]] + 18*(-141926 + 17801*\[Nu])*Cos[4*\[ScriptL]] - 
         4712362*Cos[6*\[ScriptL]] + 823543*\[Nu]*Cos[6*\[ScriptL]])*
        Sin[\[ScriptL]])/161280 + et^2*(3*\[ScriptL] + 
        (31/4 - \[Nu])*Sin[2*\[ScriptL]]) + 
      (et^4*(288*\[ScriptL] + 4*(41 + 4*\[Nu])*Sin[2*\[ScriptL]] + 
         (821 - 128*\[Nu])*Sin[4*\[ScriptL]]))/96 + 
      (et^6*(2880*\[ScriptL] - 5*(-635 + 32*\[Nu])*Sin[2*\[ScriptL]] + 
         128*(-25 + 7*\[Nu])*Sin[4*\[ScriptL]] + 3*(3815 - 648*\[Nu])*
          Sin[6*\[ScriptL]]))/960 + 
      (et^8*(967680*\[ScriptL] - 56*(-17963 + 656*\[Nu])*Sin[2*\[ScriptL]] - 
         28*(-35921 + 4736*\[Nu])*Sin[4*\[ScriptL]] - 
         3649032*Sin[6*\[ScriptL]] + 793152*\[Nu]*Sin[6*\[ScriptL]] + 
         5863513*Sin[8*\[ScriptL]] - 1048576*\[Nu]*Sin[8*\[ScriptL]]))/
       322560) + x^2*\[Epsilon]^4*(\[ScriptL]*(27/2 - 7*\[Nu]) + 
      (et*(624 - 235*\[Nu] + \[Nu]^2)*Sin[\[ScriptL]])/12 + 
      (et^3*(6948 - 2829*\[Nu] + 35*\[Nu]^2 + (3756 - 1205*\[Nu] + 3*\[Nu]^2)*
          Cos[2*\[ScriptL]])*Sin[\[ScriptL]])/48 + 
      (et^5*(2947344 - 1235201*\[Nu] + 17675*\[Nu]^2 + 
         4*(575412 - 202283*\[Nu] + 5225*\[Nu]^2)*Cos[2*\[ScriptL]] + 
         (1038048 - 338987*\[Nu] - 6055*\[Nu]^2)*Cos[4*\[ScriptL]])*
        Sin[\[ScriptL]])/11520 + (et^2*(6*\[ScriptL]*(159 - 82*\[Nu]) + 
         (969 - 326*\[Nu] + 2*\[Nu]^2)*Sin[2*\[ScriptL]]))/24 + 
      (et^4*(288*\[ScriptL]*(33 - 17*\[Nu]) + 2*(4821 - 1952*\[Nu] + 
           62*\[Nu]^2)*Sin[2*\[ScriptL]] + (5925 - 1888*\[Nu] - 11*\[Nu]^2)*
          Sin[4*\[ScriptL]]))/144 + 
      (et^6*(5*(207777 - 93122*\[Nu] + 1028*\[Nu]^2)*Sin[2*\[ScriptL]] + 
         8*(71157 - 23291*\[Nu] + 2423*\[Nu]^2)*Sin[4*\[ScriptL]] - 
         3*(960*\[ScriptL]*(-369 + 190*\[Nu]) + (-193707 + 66566*\[Nu] + 
             2172*\[Nu]^2)*Sin[6*\[ScriptL]])))/11520 + 
      (et^7*(35*(22929900 - 10801553*\[Nu] + 80303*\[Nu]^2)*Sin[\[ScriptL]] + 
         63*(4156332 - 1880513*\[Nu] + 12591*\[Nu]^2)*Sin[3*\[ScriptL]] + 
         181095684*Sin[5*\[ScriptL]] - 46504451*\[Nu]*Sin[5*\[ScriptL]] + 
         9709805*\[Nu]^2*Sin[5*\[ScriptL]] + 221342196*Sin[7*\[ScriptL]] - 
         81943639*\[Nu]*Sin[7*\[ScriptL]] - 4028087*\[Nu]^2*
          Sin[7*\[ScriptL]]))/3870720 + 
      (et^8*(114670080*\[ScriptL] - 59028480*\[ScriptL]*\[Nu] + 
         28*(4006317 - 1824214*\[Nu] + 19876*\[Nu]^2)*Sin[2*\[ScriptL]] - 
         14*(-3811935 + 1848502*\[Nu] + 16424*\[Nu]^2)*Sin[4*\[ScriptL]] + 
         44384220*Sin[6*\[ScriptL]] - 6349608*\[Nu]*Sin[6*\[ScriptL]] + 
         3696624*\[Nu]^2*Sin[6*\[ScriptL]] + 63106581*Sin[8*\[ScriptL]] - 
         25718482*\[Nu]*Sin[8*\[ScriptL]] - 1711772*\[Nu]^2*
          Sin[8*\[ScriptL]]))/967680) + 
    SO*(x^(3/2)*\[Epsilon]^3*(-4*\[ScriptL]*\[Delta]*\[Chi]A + 
        2*\[ScriptL]*(-2 + \[Nu])*\[Chi]S - 2*et*(5*\[Delta]*\[Chi]A + 
          (5 - 2*\[Nu])*\[Chi]S)*Sin[\[ScriptL]] - 
        (et^3*(113*\[Delta]*\[Chi]A + (113 - 46*\[Nu])*\[Chi]S + 
           (79*\[Delta]*\[Chi]A + (79 - 26*\[Nu])*\[Chi]S)*Cos[2*\[ScriptL]])*
          Sin[\[ScriptL]])/6 - 
        (et^5*(11*(1169*\[Delta]*\[Chi]A + 1169*\[Chi]S - 
             482*\[Nu]*\[Chi]S) + 4*(2827*\[Delta]*\[Chi]A + 
             (2827 - 1006*\[Nu])*\[Chi]S)*Cos[2*\[ScriptL]] + 
           (7513*\[Delta]*\[Chi]A + 7513*\[Chi]S - 2194*\[Nu]*\[Chi]S)*
            Cos[4*\[ScriptL]])*Sin[\[ScriptL]])/480 - 
        (et^7*(2752546*\[Delta]*\[Chi]A + 2752546*\[Chi]S - 
           1145644*\[Nu]*\[Chi]S + 3*(903319*\[Delta]*\[Chi]A + 
             (903319 - 332866*\[Nu])*\[Chi]S)*Cos[2*\[ScriptL]] + 
           6*(300029*\[Delta]*\[Chi]A + 300029*\[Chi]S - 101438*\[Nu]*
              \[Chi]S)*Cos[4*\[ScriptL]] + 1769003*\[Delta]*\[Chi]A*
            Cos[6*\[ScriptL]] + 1769003*\[Chi]S*Cos[6*\[ScriptL]] - 
           472730*\[Nu]*\[Chi]S*Cos[6*\[ScriptL]])*Sin[\[ScriptL]])/80640 + 
        et^2*(-6*\[ScriptL]*\[Delta]*\[Chi]A + 3*\[ScriptL]*(-2 + \[Nu])*
           \[Chi]S + 5*\[Nu]*\[Chi]S*Cos[\[ScriptL]]*Sin[\[ScriptL]] - 
          7*(\[Delta]*\[Chi]A + \[Chi]S)*Sin[2*\[ScriptL]]) + 
        (et^4*(-180*\[ScriptL]*(2*\[Delta]*\[Chi]A + 2*\[Chi]S - 
             \[Nu]*\[Chi]S) - 8*(38*\[Delta]*\[Chi]A + (38 - 17*\[Nu])*
              \[Chi]S)*Sin[2*\[ScriptL]] + (-334*\[Delta]*\[Chi]A - 
             334*\[Chi]S + 103*\[Nu]*\[Chi]S)*Sin[4*\[ScriptL]]))/48 + 
        (et^6*(-4200*\[ScriptL]*\[Delta]*\[Chi]A - 4200*\[ScriptL]*\[Chi]S + 
           2100*\[ScriptL]*\[Nu]*\[Chi]S - 5*(758*\[Delta]*\[Chi]A + 
             (758 - 335*\[Nu])*\[Chi]S)*Sin[2*\[ScriptL]] + 
           (-1030*\[Delta]*\[Chi]A - 1030*\[Chi]S + 643*\[Nu]*\[Chi]S)*
            Sin[4*\[ScriptL]] - 4390*\[Delta]*\[Chi]A*Sin[6*\[ScriptL]] - 
           4390*\[Chi]S*Sin[6*\[ScriptL]] + 1223*\[Nu]*\[Chi]S*
            Sin[6*\[ScriptL]]))/480 + 
        (et^8*(-1587600*\[ScriptL]*\[Delta]*\[Chi]A - 1587600*\[ScriptL]*
            \[Chi]S + 793800*\[ScriptL]*\[Nu]*\[Chi]S - 
           448*(3224*\[Delta]*\[Chi]A + (3224 - 1451*\[Nu])*\[Chi]S)*
            Sin[2*\[ScriptL]] - 56*(13622*\[Delta]*\[Chi]A + 13622*\[Chi]S - 
             5531*\[Nu]*\[Chi]S)*Sin[4*\[ScriptL]] + 373248*\[Delta]*\[Chi]A*
            Sin[6*\[ScriptL]] + 373248*\[Chi]S*Sin[6*\[ScriptL]] + 
           46656*\[Nu]*\[Chi]S*Sin[6*\[ScriptL]] - 2161382*\[Delta]*\[Chi]A*
            Sin[8*\[ScriptL]] - 2161382*\[Chi]S*Sin[8*\[ScriptL]] + 
           556403*\[Nu]*\[Chi]S*Sin[8*\[ScriptL]]))/161280) + 
      x^(5/2)*\[Epsilon]^5*((\[ScriptL]*(17*\[Delta]*(-4 + \[Nu])*\[Chi]A + 
           (-68 + 81*\[Nu] - 4*\[Nu]^2)*\[Chi]S))/2 + 
        (et*(\[Delta]*(-346 + 73*\[Nu])*\[Chi]A + 
           (-346 + 351*\[Nu] - 14*\[Nu]^2)*\[Chi]S)*Sin[\[ScriptL]])/3 + 
        (et^3*(\[Delta]*(-4602 + 1181*\[Nu])*\[Chi]A + 
           (-4602 + 4595*\[Nu] - 310*\[Nu]^2)*\[Chi]S + 
           (\[Delta]*(-1790 + 279*\[Nu])*\[Chi]A + (-1790 + 1297*\[Nu] + 30*
                \[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]])*Sin[\[ScriptL]])/12 + 
        (et^5*(-2268302*\[Delta]*\[Chi]A + 625187*\[Delta]*\[Nu]*\[Chi]A - 
           2268302*\[Chi]S + 2256597*\[Nu]*\[Chi]S - 179410*\[Nu]^2*\[Chi]S + 
           4*(\[Delta]*(-354446 + 78671*\[Nu])*\[Chi]A + 
             (-354446 + 294441*\[Nu] - 10750*\[Nu]^2)*\[Chi]S)*
            Cos[2*\[ScriptL]] + (\[Delta]*(-416954 + 35489*\[Nu])*\[Chi]A + 
             (-416954 + 136839*\[Nu] + 42890*\[Nu]^2)*\[Chi]S)*
            Cos[4*\[ScriptL]])*Sin[\[ScriptL]])/2880 - 
        (et^7*(636608492*\[Delta]*\[Chi]A - 182167190*\[Delta]*\[Nu]*
            \[Chi]A + 636608492*\[Chi]S - 632951466*\[Nu]*\[Chi]S + 
           54659668*\[Nu]^2*\[Chi]S - 9*(\[Delta]*(-53948806 + 13337795*
                \[Nu])*\[Chi]A + (-53948806 + 46886693*\[Nu] - 2642514*
                \[Nu]^2)*\[Chi]S)*Cos[2*\[ScriptL]] - 
           18*(\[Delta]*(-13354762 + 2726165*\[Nu])*\[Chi]A + 
             (-13354762 + 9453675*\[Nu] - 72086*\[Nu]^2)*\[Chi]S)*
            Cos[4*\[ScriptL]] + 68665258*\[Delta]*\[Chi]A*Cos[6*\[ScriptL]] + 
           2913515*\[Delta]*\[Nu]*\[Chi]A*Cos[6*\[ScriptL]] + 
           68665258*\[Chi]S*Cos[6*\[ScriptL]] + 28077693*\[Nu]*\[Chi]S*
            Cos[6*\[ScriptL]] - 18292162*\[Nu]^2*\[Chi]S*Cos[6*\[ScriptL]])*
          Sin[\[ScriptL]])/483840 + 
        (et^2*(6*\[ScriptL]*(\[Delta]*(-460 + 141*\[Nu])*\[Chi]A + 
             (-460 + 521*\[Nu] - 48*\[Nu]^2)*\[Chi]S) + 
           (\[Delta]*(-1988 + 365*\[Nu])*\[Chi]A + (-1988 + 1737*\[Nu] - 28*
                \[Nu]^2)*\[Chi]S)*Sin[2*\[ScriptL]]))/24 + 
        (et^4*(180*\[ScriptL]*(\[Delta]*(-716 + 231*\[Nu])*\[Chi]A + 
             (-716 + 799*\[Nu] - 84*\[Nu]^2)*\[Chi]S) + 
           8*(\[Delta]*(-14528 + 4013*\[Nu])*\[Chi]A + 
             (-14528 + 15069*\[Nu] - 1210*\[Nu]^2)*\[Chi]S)*
            Sin[2*\[ScriptL]] + (5*\[Delta]*(-8348 + 1037*\[Nu])*\[Chi]A + 
             (-41740 + 22929*\[Nu] + 2272*\[Nu]^2)*\[Chi]S)*
            Sin[4*\[ScriptL]]))/576 + 
        (et^6*(6300*\[ScriptL]*(\[Delta]*(-324 + 107*\[Nu])*\[Chi]A + 
             (-324 + 359*\[Nu] - 40*\[Nu]^2)*\[Chi]S) + 
           5*(5*\[Delta]*(-74696 + 21899*\[Nu])*\[Chi]A + 
             (-373480 + 390543*\[Nu] - 35654*\[Nu]^2)*\[Chi]S)*
            Sin[2*\[ScriptL]] + (\[Delta]*(-972164 + 279971*\[Nu])*\[Chi]A + 
             (-972164 + 1047711*\[Nu] - 93208*\[Nu]^2)*\[Chi]S)*
            Sin[4*\[ScriptL]] + 3*(\[Delta]*(-139048 + 4597*\[Nu])*\[Chi]A + 
             (-139048 + 4077*\[Nu] + 23494*\[Nu]^2)*\[Chi]S)*
            Sin[6*\[ScriptL]]))/5760 + 
        (et^8*(-974786400*\[ScriptL]*\[Delta]*\[Chi]A + 326251800*\[ScriptL]*
            \[Delta]*\[Nu]*\[Chi]A - 974786400*\[ScriptL]*\[Chi]S + 
           1075599000*\[ScriptL]*\[Nu]*\[Chi]S - 123832800*\[ScriptL]*\[Nu]^2*
            \[Chi]S + 112*(2*\[Delta]*(-4047709 + 1228816*\[Nu])*\[Chi]A + 
             (-8095418 + 8526822*\[Nu] - 835303*\[Nu]^2)*\[Chi]S)*
            Sin[2*\[ScriptL]] + 56*(\[Delta]*(-7482008 + 2143481*\[Nu])*
              \[Chi]A + (-7482008 + 7516005*\[Nu] - 655714*\[Nu]^2)*\[Chi]S)*
            Sin[4*\[ScriptL]] - 385658208*\[Delta]*\[Chi]A*
            Sin[6*\[ScriptL]] + 124142976*\[Delta]*\[Nu]*\[Chi]A*
            Sin[6*\[ScriptL]] - 385658208*\[Chi]S*Sin[6*\[ScriptL]] + 
           489480480*\[Nu]*\[Chi]S*Sin[6*\[ScriptL]] - 53942544*\[Nu]^2*
            \[Chi]S*Sin[6*\[ScriptL]] - 127947356*\[Delta]*\[Chi]A*
            Sin[8*\[ScriptL]] - 21180391*\[Delta]*\[Nu]*\[Chi]A*
            Sin[8*\[ScriptL]] - 127947356*\[Chi]S*Sin[8*\[ScriptL]] - 
           144470271*\[Nu]*\[Chi]S*Sin[8*\[ScriptL]] + 54631544*\[Nu]^2*
            \[Chi]S*Sin[8*\[ScriptL]]))/1935360)) + 
    SO^2*(x^2*\[Epsilon]^4*((3*\[ScriptL]*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + 
           \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + 
           \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S)))/2 + 
        4*et*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 4*\[Nu]*\[Chi]A^2 + 
          \[Chi]S^2 + \[Delta]*(\[Kappa]A + 2*\[Chi]A*\[Chi]S))*
         Sin[\[ScriptL]] + (et^2*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(24*\[ScriptL] + 23*Sin[2*\[ScriptL]]))/8 + 
        (et^3*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*Cos[\[ScriptL]/2]*(38*Sin[\[ScriptL]/2] + 
           11*(-Sin[(3*\[ScriptL])/2] + Sin[(5*\[ScriptL])/2])))/2 + 
        (et^4*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(864*\[ScriptL] + 748*Sin[2*\[ScriptL]] + 
           565*Sin[4*\[ScriptL]]))/192 + 
        (et^5*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*Cos[\[ScriptL]/2]*(1906*Sin[\[ScriptL]/2] - 
           746*Sin[(3*\[ScriptL])/2] + 746*Sin[(5*\[ScriptL])/2] - 
           401*Sin[(7*\[ScriptL])/2] + 401*Sin[(9*\[ScriptL])/2]))/60 + 
        (et^6*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(3840*\[ScriptL] + 3445*Sin[2*\[ScriptL]] + 
           1414*Sin[4*\[ScriptL]] + 2519*Sin[6*\[ScriptL]]))/640 + 
        (et^7*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*Cos[\[ScriptL]/2]*
          (3692396*Sin[\[ScriptL]/2] - 1662081*Sin[(3*\[ScriptL])/2] + 
           1662081*Sin[(5*\[ScriptL])/2] - 1017402*Sin[(7*\[ScriptL])/2] + 
           1017402*Sin[(9*\[ScriptL])/2] - 766319*Sin[(11*\[ScriptL])/2] + 
           766319*Sin[(13*\[ScriptL])/2]))/80640 + 
        (et^8*(\[Kappa]S - 2*\[Kappa]S*\[Nu] + \[Chi]A^2 - 
           4*\[Nu]*\[Chi]A^2 + \[Chi]S^2 + \[Delta]*(\[Kappa]A + 
             2*\[Chi]A*\[Chi]S))*(4838400*\[ScriptL] + 
           4405016*Sin[2*\[ScriptL]] + 2165884*Sin[4*\[ScriptL]] + 
           476424*Sin[6*\[ScriptL]] + 3766361*Sin[8*\[ScriptL]]))/645120) + 
      x^3*\[Epsilon]^6*(-(\[ScriptL]*(\[Delta]*\[Kappa]A*(-39 + 17*\[Nu]) + 
            \[Kappa]S*(-39 + 95*\[Nu] - 10*\[Nu]^2) - 67*\[Chi]A^2 + 
            279*\[Nu]*\[Chi]A^2 - 20*\[Nu]^2*\[Chi]A^2 + 
            2*\[Delta]*(-67 + 55*\[Nu])*\[Chi]A*\[Chi]S - 67*\[Chi]S^2 + 
            99*\[Nu]*\[Chi]S^2 - 28*\[Nu]^2*\[Chi]S^2))/2 - 
        (et*(\[Delta]*\[Kappa]A*(-198 + 83*\[Nu]) + \[Kappa]S*
            (-198 + 479*\[Nu] - 46*\[Nu]^2) - 2*(177 - 730*\[Nu] + 
             46*\[Nu]^2)*\[Chi]A^2 + 4*\[Delta]*(-177 + 121*\[Nu])*\[Chi]A*
            \[Chi]S - 2*(177 - 220*\[Nu] + 54*\[Nu]^2)*\[Chi]S^2)*
          Sin[\[ScriptL]])/3 + (et^3*(3110*\[Kappa]S - 7555*\[Kappa]S*\[Nu] + 
           1010*\[Kappa]S*\[Nu]^2 - 5*\[Delta]*\[Kappa]A*(-622 + 267*\[Nu]) + 
           5298*\[Chi]A^2 - 21995*\[Nu]*\[Chi]A^2 + 2020*\[Nu]^2*\[Chi]A^2 + 
           2*\[Delta]*(5298 - 3491*\[Nu])*\[Chi]A*\[Chi]S + 5298*\[Chi]S^2 - 
           6179*\[Nu]*\[Chi]S^2 + 1420*\[Nu]^2*\[Chi]S^2 + 
           (\[Delta]*\[Kappa]A*(994 - 441*\[Nu]) + \[Kappa]S*
              (994 - 2429*\[Nu] + 142*\[Nu]^2) + 1902*\[Chi]A^2 - 
             7747*\[Nu]*\[Chi]A^2 + 284*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(1902 - 919*\[Nu])*\[Chi]A*\[Chi]S + 1902*\[Chi]S^2 - 
             1699*\[Nu]*\[Chi]S^2 + 332*\[Nu]^2*\[Chi]S^2)*Cos[2*\[ScriptL]])*
          Sin[\[ScriptL]])/12 - (et^5*(-1750710*\[Delta]*\[Kappa]A - 
           1750710*\[Kappa]S + 758971*\[Delta]*\[Kappa]A*\[Nu] + 
           4260391*\[Kappa]S*\[Nu] - 640862*\[Kappa]S*\[Nu]^2 - 
           2915370*\[Chi]A^2 + 12142366*\[Nu]*\[Chi]A^2 - 
           1281724*\[Nu]^2*\[Chi]A^2 - 5830740*\[Delta]*\[Chi]A*\[Chi]S + 
           3782432*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S - 2915370*\[Chi]S^2 + 
           3301546*\[Nu]*\[Chi]S^2 - 732420*\[Nu]^2*\[Chi]S^2 + 
           8*(-2*\[Kappa]S*(59385 - 145027*\[Nu] + 16529*\[Nu]^2) - 
             210420*\[Chi]A^2 + 866639*\[Nu]*\[Chi]A^2 - 66116*\[Nu]^2*
              \[Chi]A^2 - 210420*\[Chi]S^2 + 200399*\[Nu]*\[Chi]S^2 - 
             38970*\[Nu]^2*\[Chi]S^2 + 2*\[Delta]*(\[Kappa]A*(-59385 + 
                 26257*\[Nu]) + 7*(-30060 + 16097*\[Nu])*\[Chi]A*\[Chi]S))*
            Cos[2*\[ScriptL]] + (\[Delta]*\[Kappa]A*(-215850 + 108877*
                \[Nu]) + \[Kappa]S*(-215850 + 540577*\[Nu] + 5566*\[Nu]^2) + 
             2*(-232515 + 931901*\[Nu] + 5566*\[Nu]^2)*\[Chi]A^2 + 
             4*\[Delta]*(-232515 + 52076*\[Nu])*\[Chi]A*\[Chi]S - 
             2*(232515 - 102311*\[Nu] + 7110*\[Nu]^2)*\[Chi]S^2)*
            Cos[4*\[ScriptL]])*Sin[\[ScriptL]])/2880 - 
        (et^7*(-547054260*\[Delta]*\[Kappa]A - 547054260*\[Kappa]S + 
           238427530*\[Delta]*\[Kappa]A*\[Nu] + 1332536050*\[Kappa]S*\[Nu] - 
           213069308*\[Kappa]S*\[Nu]^2 - 898941564*\[Chi]A^2 + 
           3751264438*\[Nu]*\[Chi]A^2 - 426138616*\[Nu]^2*\[Chi]A^2 - 
           1797883128*\[Delta]*\[Chi]A*\[Chi]S + 1156897340*\[Delta]*\[Nu]*
            \[Chi]A*\[Chi]S - 898941564*\[Chi]S^2 + 1001399158*\[Nu]*
            \[Chi]S^2 - 217405320*\[Nu]^2*\[Chi]S^2 - 
           3*(\[Delta]*\[Kappa]A*(124613970 - 55118225*\[Nu]) + 
             \[Kappa]S*(124613970 - 304346165*\[Nu] + 41246182*\[Nu]^2) + 
             214024806*\[Chi]A^2 - 885760667*\[Nu]*\[Chi]A^2 + 
             82492364*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(214024806 - 119056055*
                \[Nu])*\[Chi]A*\[Chi]S + 214024806*\[Chi]S^2 - 
             208450667*\[Nu]*\[Chi]S^2 + 40514100*\[Nu]^2*\[Chi]S^2)*
            Cos[2*\[ScriptL]] + 6*(29*\[Delta]*\[Kappa]A*(-914010 + 421061*
                \[Nu]) + \[Kappa]S*(-26506290 + 65223349*\[Nu] - 6558278*
                \[Nu]^2) - 48273030*\[Chi]A^2 + 197779999*\[Nu]*\[Chi]A^2 - 
             13116556*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(-48273030 + 20874787*
                \[Nu])*\[Chi]A*\[Chi]S - 48273030*\[Chi]S^2 + 
             37061695*\[Nu]*\[Chi]S^2 - 5691828*\[Nu]^2*\[Chi]S^2)*
            Cos[4*\[ScriptL]] - 29874090*\[Delta]*\[Kappa]A*
            Cos[6*\[ScriptL]] - 29874090*\[Kappa]S*Cos[6*\[ScriptL]] + 
           19977821*\[Delta]*\[Kappa]A*\[Nu]*Cos[6*\[ScriptL]] + 
           79726001*\[Kappa]S*\[Nu]*Cos[6*\[ScriptL]] + 13761362*\[Kappa]S*
            \[Nu]^2*Cos[6*\[ScriptL]] - 83416878*\[Chi]A^2*
            Cos[6*\[ScriptL]] + 325325087*\[Nu]*\[Chi]A^2*Cos[6*\[ScriptL]] + 
           27522724*\[Nu]^2*\[Chi]A^2*Cos[6*\[ScriptL]] - 166833756*\[Delta]*
            \[Chi]A*\[Chi]S*Cos[6*\[ScriptL]] - 34384874*\[Delta]*\[Nu]*
            \[Chi]A*\[Chi]S*Cos[6*\[ScriptL]] - 83416878*\[Chi]S^2*
            Cos[6*\[ScriptL]] - 26042449*\[Nu]*\[Chi]S^2*Cos[6*\[ScriptL]] + 
           15903708*\[Nu]^2*\[Chi]S^2*Cos[6*\[ScriptL]])*Sin[\[ScriptL]])/
         483840 + (et^2*(6*\[ScriptL]*(\[Delta]*\[Kappa]A*(312 - 137*\[Nu]) + 
             \[Kappa]S*(312 - 761*\[Nu] + 118*\[Nu]^2) + 504*\[Chi]A^2 - 
             2114*\[Nu]*\[Chi]A^2 + 236*\[Nu]^2*\[Chi]A^2 + 
             8*\[Delta]*(126 - 95*\[Nu])*\[Chi]A*\[Chi]S + 504*\[Chi]S^2 - 
             662*\[Nu]*\[Chi]S^2 + 168*\[Nu]^2*\[Chi]S^2) + 
           (\[Delta]*\[Kappa]A*(1125 - 481*\[Nu]) + \[Kappa]S*
              (1125 - 2731*\[Nu] + 218*\[Nu]^2) + 2073*\[Chi]A^2 - 
             8497*\[Nu]*\[Chi]A^2 + 436*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(2073 - 1213*\[Nu])*\[Chi]A*\[Chi]S + 
             2073*\[Chi]S^2 - 2221*\[Nu]*\[Chi]S^2 + 492*\[Nu]^2*\[Chi]S^2)*
            Sin[2*\[ScriptL]]))/24 + 
        (et^4*(432*\[ScriptL]*(\[Delta]*\[Kappa]A*(234 - 103*\[Nu]) + 
             \[Kappa]S*(234 - 571*\[Nu] + 98*\[Nu]^2) + 
             2*(185 - 778*\[Nu] + 98*\[Nu]^2)*\[Chi]A^2 + 20*\[Delta]*
              (37 - 27*\[Nu])*\[Chi]A*\[Chi]S + 2*(185 - 232*\[Nu] + 56*
                \[Nu]^2)*\[Chi]S^2) + 4*(\[Delta]*\[Kappa]A*(20949 - 8906*
                \[Nu]) + \[Kappa]S*(20949 - 50804*\[Nu] + 7528*\[Nu]^2) + 
             34953*\[Chi]A^2 - 145463*\[Nu]*\[Chi]A^2 + 15056*\[Nu]^2*
              \[Chi]A^2 + 2*\[Delta]*(34953 - 23333*\[Nu])*\[Chi]A*\[Chi]S + 
             34953*\[Chi]S^2 - 41015*\[Nu]*\[Chi]S^2 + 9132*\[Nu]^2*
              \[Chi]S^2)*Sin[2*\[ScriptL]] + 
           (\[Delta]*\[Kappa]A*(22533 - 10544*\[Nu]) + \[Kappa]S*
              (22533 - 55610*\[Nu] + 1660*\[Nu]^2) + 45321*\[Chi]A^2 - 
             183263*\[Nu]*\[Chi]A^2 + 3320*\[Nu]^2*\[Chi]A^2 + 
             2*\[Delta]*(45321 - 16601*\[Nu])*\[Chi]A*\[Chi]S + 
             45321*\[Chi]S^2 - 31223*\[Nu]*\[Chi]S^2 + 4884*\[Nu]^2*
              \[Chi]S^2)*Sin[4*\[ScriptL]]))/576 + 
        et^6*(312*\[ScriptL]*\[Delta]*\[Kappa]A + 312*\[ScriptL]*\[Kappa]S - 
          (275*\[ScriptL]*\[Delta]*\[Kappa]A*\[Nu])/2 - 
          (1523*\[ScriptL]*\[Kappa]S*\[Nu])/2 + 137*\[ScriptL]*\[Kappa]S*
           \[Nu]^2 + 488*\[ScriptL]*\[Chi]A^2 - 2055*\[ScriptL]*\[Nu]*
           \[Chi]A^2 + 274*\[ScriptL]*\[Nu]^2*\[Chi]A^2 + 
          976*\[ScriptL]*\[Delta]*\[Chi]A*\[Chi]S - 700*\[ScriptL]*\[Delta]*
           \[Nu]*\[Chi]A*\[Chi]S + 488*\[ScriptL]*\[Chi]S^2 - 
          597*\[ScriptL]*\[Nu]*\[Chi]S^2 + 140*\[ScriptL]*\[Nu]^2*\[Chi]S^2 - 
          ((\[Delta]*\[Kappa]A*(-210672 + 90317*\[Nu]) + 
             \[Kappa]S*(-210672 + 511661*\[Nu] - 83102*\[Nu]^2) - 
             2*(172304 - 719669*\[Nu] + 83102*\[Nu]^2)*\[Chi]A^2 + 
             8*\[Delta]*(-86152 + 58065*\[Nu])*\[Chi]A*\[Chi]S - 
             2*(172304 - 201807*\[Nu] + 44976*\[Nu]^2)*\[Chi]S^2)*
            Sin[2*\[ScriptL]])/768 - 
          ((\[Delta]*\[Kappa]A*(-120570 + 50659*\[Nu]) + 
             \[Kappa]S*(-120570 + 291799*\[Nu] - 46402*\[Nu]^2) - 
             2*(98605 - 410764*\[Nu] + 46402*\[Nu]^2)*\[Chi]A^2 + 
             4*\[Delta]*(-98605 + 65139*\[Nu])*\[Chi]A*\[Chi]S - 
             2*(98605 - 113934*\[Nu] + 23418*\[Nu]^2)*\[Chi]S^2)*
            Sin[4*\[ScriptL]])/960 + (2249*\[Delta]*\[Kappa]A*
            Sin[6*\[ScriptL]])/64 + (2249*\[Kappa]S*Sin[6*\[ScriptL]])/64 - 
          (75941*\[Delta]*\[Kappa]A*\[Nu]*Sin[6*\[ScriptL]])/3840 - 
          (345821*\[Kappa]S*\[Nu]*Sin[6*\[ScriptL]])/3840 - 
          (12241*\[Kappa]S*\[Nu]^2*Sin[6*\[ScriptL]])/1920 + 
          (16067*\[Chi]A^2*Sin[6*\[ScriptL]])/192 - 
          (212197*\[Nu]*\[Chi]A^2*Sin[6*\[ScriptL]])/640 - 
          (12241*\[Nu]^2*\[Chi]A^2*Sin[6*\[ScriptL]])/960 + 
          (16067*\[Delta]*\[Chi]A*\[Chi]S*Sin[6*\[ScriptL]])/96 - 
          (104*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S*Sin[6*\[ScriptL]])/15 + 
          (16067*\[Chi]S^2*Sin[6*\[ScriptL]])/192 - 
          (6467*\[Nu]*\[Chi]S^2*Sin[6*\[ScriptL]])/640 - 
          (1301*\[Nu]^2*\[Chi]S^2*Sin[6*\[ScriptL]])/240) + 
        (et^8*(943488000*\[ScriptL]*\[Delta]*\[Kappa]A + 943488000*\[ScriptL]*
            \[Kappa]S - 416102400*\[ScriptL]*\[Delta]*\[Kappa]A*\[Nu] - 
           2303078400*\[ScriptL]*\[Kappa]S*\[Nu] + 425779200*\[ScriptL]*
            \[Kappa]S*\[Nu]^2 + 1466035200*\[ScriptL]*\[Chi]A^2 - 
           6178636800*\[ScriptL]*\[Nu]*\[Chi]A^2 + 851558400*\[ScriptL]*
            \[Nu]^2*\[Chi]A^2 + 2932070400*\[ScriptL]*\[Delta]*\[Chi]A*
            \[Chi]S - 2080512000*\[ScriptL]*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S + 
           1466035200*\[ScriptL]*\[Chi]S^2 - 1766016000*\[ScriptL]*\[Nu]*
            \[Chi]S^2 + 406425600*\[ScriptL]*\[Nu]^2*\[Chi]S^2 + 
           56*(\[Delta]*\[Kappa]A*(15271743 - 6577238*\[Nu]) + 
             \[Kappa]S*(15271743 - 37120724*\[Nu] + 6336040*\[Nu]^2) + 
             24667419*\[Chi]A^2 - 103250849*\[Nu]*\[Chi]A^2 + 
             12672080*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(24667419 - 16699379*
                \[Nu])*\[Chi]A*\[Chi]S + 24667419*\[Chi]S^2 - 
             28817585*\[Nu]*\[Chi]S^2 + 6395748*\[Nu]^2*\[Chi]S^2)*
            Sin[2*\[ScriptL]] + 28*(\[Delta]*\[Kappa]A*(13679031 - 5880533*
                \[Nu]) + \[Kappa]S*(13679031 - 33238595*\[Nu] + 5350402*
                \[Nu]^2) + 22517715*\[Chi]A^2 - 93885983*\[Nu]*\[Chi]A^2 + 
             10700804*\[Nu]^2*\[Chi]A^2 + 2*\[Delta]*(22517715 - 14479559*
                \[Nu])*\[Chi]A*\[Chi]S + 22517715*\[Chi]S^2 - 
             25143995*\[Nu]*\[Chi]S^2 + 5339004*\[Nu]^2*\[Chi]S^2)*
            Sin[4*\[ScriptL]] + 292065048*\[Delta]*\[Kappa]A*
            Sin[6*\[ScriptL]] + 292065048*\[Kappa]S*Sin[6*\[ScriptL]] - 
           116900784*\[Delta]*\[Kappa]A*\[Nu]*Sin[6*\[ScriptL]] - 
           701030880*\[Kappa]S*\[Nu]*Sin[6*\[ScriptL]] + 130239936*\[Kappa]S*
            \[Nu]^2*Sin[6*\[ScriptL]] + 450447480*\[Chi]A^2*
            Sin[6*\[ScriptL]] - 1889441064*\[Nu]*\[Chi]A^2*
            Sin[6*\[ScriptL]] + 260479872*\[Nu]^2*\[Chi]A^2*
            Sin[6*\[ScriptL]] + 900894960*\[Delta]*\[Chi]A*\[Chi]S*
            Sin[6*\[ScriptL]] - 674563824*\[Delta]*\[Nu]*\[Chi]A*\[Chi]S*
            Sin[6*\[ScriptL]] + 450447480*\[Chi]S^2*Sin[6*\[ScriptL]] - 
           586912680*\[Nu]*\[Chi]S^2*Sin[6*\[ScriptL]] + 117834912*\[Nu]^2*
            \[Chi]S^2*Sin[6*\[ScriptL]] + 44748573*\[Delta]*\[Kappa]A*
            Sin[8*\[ScriptL]] + 44748573*\[Kappa]S*Sin[8*\[ScriptL]] - 
           41022478*\[Delta]*\[Kappa]A*\[Nu]*Sin[8*\[ScriptL]] - 
           130519624*\[Kappa]S*\[Nu]*Sin[8*\[ScriptL]] - 49674400*\[Kappa]S*
            \[Nu]^2*Sin[8*\[ScriptL]] + 167973729*\[Chi]A^2*
            Sin[8*\[ScriptL]] - 639926539*\[Nu]*\[Chi]A^2*Sin[8*\[ScriptL]] - 
           99348800*\[Nu]^2*\[Chi]A^2*Sin[8*\[ScriptL]] + 335947458*\[Delta]*
            \[Chi]A*\[Chi]S*Sin[8*\[ScriptL]] + 190155262*\[Delta]*\[Nu]*
            \[Chi]A*\[Chi]S*Sin[8*\[ScriptL]] + 167973729*\[Chi]S^2*
            Sin[8*\[ScriptL]] + 158186885*\[Nu]*\[Chi]S^2*Sin[8*\[ScriptL]] - 
           62184732*\[Nu]^2*\[Chi]S^2*Sin[8*\[ScriptL]]))/1935360))
