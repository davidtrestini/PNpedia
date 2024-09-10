(* ::Package:: *)

(* ::Chapter:: *)
(*Supplementary file to "Complete gravitational-waveform amplitude modes for quasi-circular compact binaries to the 3.5PN order"*)


(* ::Text:: *)
(*Author: Quentin Henry*)
(*This file contains the spin-weighted spherical harmonics modes including non-spinning and spinning contributions to the 3.5PN order, for quasi-circular orbits and aligned spins (See Sec. III. A. in the paper)*)


(*
Notation:
M = m1 + m2 is the total mass
\[Nu] = m1 m2 / M is the symmetric mass ratio 
\[Delta] = (m1 - m2) / M is the antisymmetric mass ratio 
x = (G M \[Omega] / c^3)^(2/3), where \[Omega] is the orbital frequency
\[Psi] is a phase defined in Eq. (3.1)
R is the distance from the observer to the source
S = S1 + S2
\[CapitalSigma] = M/m2 S2 - M/m1 S1
(\[Kappa]1, \[Kappa]2, \[Lambda]1, \[Lambda]2) are the spin multipole constants, which equal 1 for black holes
\[Kappa]p = \[Kappa]1 + \[Kappa]2, \[Kappa]m = \[Kappa]1 - \[Kappa]2
\[Lambda]p = \[Lambda]1 + \[Lambda]2, \[Lambda]m = \[Lambda]1 - \[Lambda]2  
*)


h[2, 0] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(-5/(14*Sqrt[6]) + 
      (20375*x)/(56448*Sqrt[6]) - (335*x*\[Nu])/(672*Sqrt[6])))/(c^2*R)
 
h[2, 1] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((I/3)*Sqrt[x]*\[Delta] - 
      ((17*I)/84)*x^(3/2)*\[Delta] + (x^2*\[Delta])/6 + 
      (I/3)*Pi*x^2*\[Delta] - (((43*I)/21)*S*x^2*\[Delta])/(G*M^2) - 
      ((43*I)/378)*x^(5/2)*\[Delta] + (I*S^2*x^(5/2)*\[Delta])/(G^2*M^4) - 
      (17*x^3*\[Delta])/168 - ((17*I)/84)*Pi*x^3*\[Delta] - 
      (((331*I)/756)*S*x^3*\[Delta])/(G*M^2) + ((15223771*I)/4365900)*x^(7/2)*
       \[Delta] - ((214*I)/315)*EulerGamma*x^(7/2)*\[Delta] - 
      (109*Pi*x^(7/2)*\[Delta])/630 + (I/18)*Pi^2*x^(7/2)*\[Delta] - 
      (181*S*x^(7/2)*\[Delta])/(210*G*M^2) - 
      (((43*I)/21)*Pi*S*x^(7/2)*\[Delta])/(G*M^2) + 
      (((41*I)/42)*S^2*x^(7/2)*\[Delta])/(G^2*M^4) - 
      ((I/3)*S^2*x^(5/2)*\[Kappa]m)/(G^2*M^4) + 
      (((23*I)/48)*S^2*x^(7/2)*\[Kappa]m)/(G^2*M^4) + 
      ((I/2)*S^2*x^(5/2)*\[Delta]*\[Kappa]p)/(G^2*M^4) + 
      (((47*I)/336)*S^2*x^(7/2)*\[Delta]*\[Kappa]p)/(G^2*M^4) + 
      ((5*I)/21)*x^(3/2)*\[Delta]*\[Nu] - ((509*I)/378)*x^(5/2)*\[Delta]*
       \[Nu] + (353*x^3*\[Delta]*\[Nu])/84 + (I/14)*Pi*x^3*\[Delta]*\[Nu] + 
      (((386*I)/189)*S*x^3*\[Delta]*\[Nu])/(G*M^2) - 
      ((102119*I)/7128)*x^(7/2)*\[Delta]*\[Nu] + ((205*I)/384)*Pi^2*x^(7/2)*
       \[Delta]*\[Nu] - ((I/7)*S^2*x^(7/2)*\[Delta]*\[Nu])/(G^2*M^4) - 
      (((191*I)/72)*S^2*x^(7/2)*\[Kappa]m*\[Nu])/(G^2*M^4) - 
      ((I/14)*S^2*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu])/(G^2*M^4) + 
      ((79*I)/504)*x^(5/2)*\[Delta]*\[Nu]^2 - ((4211*I)/24948)*x^(7/2)*
       \[Delta]*\[Nu]^2 + ((2263*I)/24948)*x^(7/2)*\[Delta]*\[Nu]^3 + 
      ((I/2)*x*\[CapitalSigma])/(G*M^2) - (((79*I)/42)*x^2*\[CapitalSigma])/
       (G*M^2) + (x^(5/2)*\[CapitalSigma])/(4*G*M^2) + 
      ((I/2)*Pi*x^(5/2)*\[CapitalSigma])/(G*M^2) - 
      ((I/3)*S*x^(5/2)*\[CapitalSigma])/(G^2*M^4) + 
      (((293*I)/756)*x^3*\[CapitalSigma])/(G*M^2) + 
      ((I/2)*S^2*x^3*\[CapitalSigma])/(G^3*M^6) - 
      (79*x^(7/2)*\[CapitalSigma])/(84*G*M^2) - 
      (((79*I)/42)*Pi*x^(7/2)*\[CapitalSigma])/(G*M^2) - 
      (((29*I)/21)*S*x^(7/2)*\[CapitalSigma])/(G^2*M^4) - 
      (((5*I)/6)*S*x^(5/2)*\[Delta]*\[Kappa]m*\[CapitalSigma])/(G^2*M^4) + 
      (((19*I)/56)*S*x^(7/2)*\[Delta]*\[Kappa]m*\[CapitalSigma])/(G^2*M^4) + 
      (((5*I)/6)*S*x^(5/2)*\[Kappa]p*\[CapitalSigma])/(G^2*M^4) + 
      ((I/4)*S^2*x^3*\[Kappa]p*\[CapitalSigma])/(G^3*M^6) - 
      (((19*I)/56)*S*x^(7/2)*\[Kappa]p*\[CapitalSigma])/(G^2*M^4) + 
      (((139*I)/42)*x^2*\[Nu]*\[CapitalSigma])/(G*M^2) - 
      ((4*I)*S*x^(5/2)*\[Nu]*\[CapitalSigma])/(G^2*M^4) - 
      (((2615*I)/1512)*x^3*\[Nu]*\[CapitalSigma])/(G*M^2) + 
      (1951*x^(7/2)*\[Nu]*\[CapitalSigma])/(280*G*M^2) + 
      (((257*I)/84)*Pi*x^(7/2)*\[Nu]*\[CapitalSigma])/(G*M^2) + 
      (((100*I)/21)*S*x^(7/2)*\[Nu]*\[CapitalSigma])/(G^2*M^4) - 
      (((1301*I)/504)*S*x^(7/2)*\[Delta]*\[Kappa]m*\[Nu]*\[CapitalSigma])/
       (G^2*M^4) - ((2*I)*S*x^(5/2)*\[Kappa]p*\[Nu]*\[CapitalSigma])/
       (G^2*M^4) + (((1019*I)/504)*S*x^(7/2)*\[Kappa]p*\[Nu]*\[CapitalSigma])/
       (G^2*M^4) - (((1723*I)/378)*x^3*\[Nu]^2*\[CapitalSigma])/(G*M^2) + 
      (((4*I)/7)*S*x^(7/2)*\[Nu]^2*\[CapitalSigma])/(G^2*M^4) + 
      (((2*I)/7)*S*x^(7/2)*\[Kappa]p*\[Nu]^2*\[CapitalSigma])/(G^2*M^4) - 
      ((I/2)*x^(5/2)*\[Delta]*\[CapitalSigma]^2)/(G^2*M^4) + 
      ((I/2)*S*x^3*\[Delta]*\[CapitalSigma]^2)/(G^3*M^6) - 
      (((6*I)/7)*x^(7/2)*\[Delta]*\[CapitalSigma]^2)/(G^2*M^4) - 
      (((5*I)/12)*x^(5/2)*\[Kappa]m*\[CapitalSigma]^2)/(G^2*M^4) - 
      ((I/4)*S*x^3*\[Kappa]m*\[CapitalSigma]^2)/(G^3*M^6) + 
      (((19*I)/112)*x^(7/2)*\[Kappa]m*\[CapitalSigma]^2)/(G^2*M^4) + 
      (((5*I)/12)*x^(5/2)*\[Delta]*\[Kappa]p*\[CapitalSigma]^2)/(G^2*M^4) + 
      ((I/4)*S*x^3*\[Delta]*\[Kappa]p*\[CapitalSigma]^2)/(G^3*M^6) - 
      (((19*I)/112)*x^(7/2)*\[Delta]*\[Kappa]p*\[CapitalSigma]^2)/(G^2*M^4) - 
      (I*x^(5/2)*\[Delta]*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) + 
      (((59*I)/21)*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) + 
      (((4*I)/3)*x^(5/2)*\[Kappa]m*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) - 
      (((751*I)/504)*x^(7/2)*\[Kappa]m*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) - 
      ((I/2)*x^(5/2)*\[Delta]*\[Kappa]p*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) + 
      (((145*I)/126)*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu]*\[CapitalSigma]^2)/
       (G^2*M^4) + ((I/7)*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma]^2)/
       (G^2*M^4) + (((1265*I)/504)*x^(7/2)*\[Kappa]m*\[Nu]^2*
        \[CapitalSigma]^2)/(G^2*M^4) + ((I/14)*x^(7/2)*\[Delta]*\[Kappa]p*
        \[Nu]^2*\[CapitalSigma]^2)/(G^2*M^4) - 
      ((I/8)*x^3*\[Delta]*\[Kappa]m*\[CapitalSigma]^3)/(G^3*M^6) + 
      ((I/8)*x^3*\[Kappa]p*\[CapitalSigma]^3)/(G^3*M^6) - 
      ((I/2)*x^3*\[Nu]*\[CapitalSigma]^3)/(G^3*M^6) - 
      ((I/4)*x^3*\[Kappa]p*\[Nu]*\[CapitalSigma]^3)/(G^3*M^6) + 
      (2*x^2*\[Delta]*Log[2])/3 - (17*x^3*\[Delta]*Log[2])/42 - 
      ((319*I)/315)*x^(7/2)*\[Delta]*Log[2] + (2*Pi*x^(7/2)*\[Delta]*Log[2])/
       3 - (86*S*x^(7/2)*\[Delta]*Log[2])/(21*G*M^2) + 
      (x^3*\[Delta]*\[Nu]*Log[2])/7 + (x^(5/2)*\[CapitalSigma]*Log[2])/
       (G*M^2) - (79*x^(7/2)*\[CapitalSigma]*Log[2])/(21*G*M^2) + 
      (257*x^(7/2)*\[Nu]*\[CapitalSigma]*Log[2])/(42*G*M^2) - 
      ((2*I)/3)*x^(7/2)*\[Delta]*Log[2]^2 - ((107*I)/315)*x^(7/2)*\[Delta]*
       Log[x]))/(c^2*E^(I*\[Psi])*R)
 
h[2, 2] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(1 - (107*x)/42 + 2*Pi*x^(3/2) - 
      (2*S*x^(3/2))/(G*M^2) - (2173*x^2)/1512 + (2*S^2*x^2)/(G^2*M^4) - 
      (107*Pi*x^(5/2))/21 - (163*S*x^(5/2))/(63*G*M^2) + 
      (27027409*x^3)/646800 - (856*EulerGamma*x^3)/105 + 
      ((428*I)/105)*Pi*x^3 + (2*Pi^2*x^3)/3 - (((4*I)/3)*S*x^3)/(G*M^2) - 
      (4*Pi*S*x^3)/(G*M^2) - (404*S^2*x^3)/(63*G^2*M^4) - 
      (2173*Pi*x^(7/2))/756 + (1061*S*x^(7/2))/(84*G*M^2) + 
      (4*Pi*S^2*x^(7/2))/(G^2*M^4) + (32*S^3*x^(7/2))/(3*G^3*M^6) + 
      (55*S^2*x^3*\[Delta]*\[Kappa]m)/(42*G^2*M^4) + 
      (S^2*x^2*\[Kappa]p)/(G^2*M^4) - (31*S^2*x^3*\[Kappa]p)/(42*G^2*M^4) + 
      (2*Pi*S^2*x^(7/2)*\[Kappa]p)/(G^2*M^4) - (2*S^3*x^(7/2)*\[Kappa]p)/
       (3*G^3*M^6) - (2*S^3*x^(7/2)*\[Lambda]p)/(G^3*M^6) + (55*x*\[Nu])/42 - 
      (1069*x^2*\[Nu])/216 - (24*I)*x^(5/2)*\[Nu] + (34*Pi*x^(5/2)*\[Nu])/
       21 - (92*S*x^(5/2)*\[Nu])/(63*G*M^2) - (278185*x^3*\[Nu])/33264 + 
      (41*Pi^2*x^3*\[Nu])/96 + (68*S^2*x^3*\[Nu])/(21*G^2*M^4) + 
      ((14333*I)/162)*x^(7/2)*\[Nu] - (2495*Pi*x^(7/2)*\[Nu])/378 + 
      (4043*S*x^(7/2)*\[Nu])/(84*G*M^2) + (34*S^2*x^3*\[Kappa]p*\[Nu])/
       (21*G^2*M^4) + (2047*x^2*\[Nu]^2)/1512 - (20261*x^3*\[Nu]^2)/2772 - 
      ((4066*I)/945)*x^(7/2)*\[Nu]^2 + (40*Pi*x^(7/2)*\[Nu]^2)/27 + 
      (499*S*x^(7/2)*\[Nu]^2)/(84*G*M^2) + (114635*x^3*\[Nu]^3)/99792 - 
      (2*x^(3/2)*\[Delta]*\[CapitalSigma])/(3*G*M^2) + 
      (2*S*x^2*\[Delta]*\[CapitalSigma])/(G^2*M^4) - 
      (x^(5/2)*\[Delta]*\[CapitalSigma])/(21*G*M^2) - 
      (4*Pi*x^3*\[Delta]*\[CapitalSigma])/(3*G*M^2) - 
      (481*S*x^3*\[Delta]*\[CapitalSigma])/(63*G^2*M^4) + 
      (3931*x^(7/2)*\[Delta]*\[CapitalSigma])/(756*G*M^2) + 
      (4*Pi*S*x^(7/2)*\[Delta]*\[CapitalSigma])/(G^2*M^4) + 
      (52*S^2*x^(7/2)*\[Delta]*\[CapitalSigma])/(3*G^3*M^6) - 
      (S*x^2*\[Kappa]m*\[CapitalSigma])/(G^2*M^4) + 
      (43*S*x^3*\[Kappa]m*\[CapitalSigma])/(21*G^2*M^4) - 
      (2*Pi*S*x^(7/2)*\[Kappa]m*\[CapitalSigma])/(G^2*M^4) - 
      (7*S^2*x^(7/2)*\[Kappa]m*\[CapitalSigma])/(3*G^3*M^6) + 
      (S*x^2*\[Delta]*\[Kappa]p*\[CapitalSigma])/(G^2*M^4) - 
      (43*S*x^3*\[Delta]*\[Kappa]p*\[CapitalSigma])/(21*G^2*M^4) + 
      (2*Pi*S*x^(7/2)*\[Delta]*\[Kappa]p*\[CapitalSigma])/(G^2*M^4) - 
      (S^2*x^(7/2)*\[Delta]*\[Kappa]p*\[CapitalSigma])/(3*G^3*M^6) + 
      (3*S^2*x^(7/2)*\[Lambda]m*\[CapitalSigma])/(G^3*M^6) - 
      (3*S^2*x^(7/2)*\[Delta]*\[Lambda]p*\[CapitalSigma])/(G^3*M^6) + 
      (20*x^(5/2)*\[Delta]*\[Nu]*\[CapitalSigma])/(63*G*M^2) + 
      (68*S*x^3*\[Delta]*\[Nu]*\[CapitalSigma])/(21*G^2*M^4) + 
      (7813*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma])/(378*G*M^2) - 
      (48*S*x^3*\[Kappa]m*\[Nu]*\[CapitalSigma])/(7*G^2*M^4) + 
      (34*S*x^3*\[Delta]*\[Kappa]p*\[Nu]*\[CapitalSigma])/(21*G^2*M^4) + 
      (1025*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma])/(252*G*M^2) - 
      (5*x^3*\[CapitalSigma]^2)/(3*G^2*M^4) + 
      (20*S*x^(7/2)*\[CapitalSigma]^2)/(3*G^3*M^6) - 
      (x^2*\[Delta]*\[Kappa]m*\[CapitalSigma]^2)/(2*G^2*M^4) + 
      (43*x^3*\[Delta]*\[Kappa]m*\[CapitalSigma]^2)/(42*G^2*M^4) - 
      (Pi*x^(7/2)*\[Delta]*\[Kappa]m*\[CapitalSigma]^2)/(G^2*M^4) - 
      (3*S*x^(7/2)*\[Delta]*\[Kappa]m*\[CapitalSigma]^2)/(G^3*M^6) + 
      (x^2*\[Kappa]p*\[CapitalSigma]^2)/(2*G^2*M^4) - 
      (43*x^3*\[Kappa]p*\[CapitalSigma]^2)/(42*G^2*M^4) + 
      (Pi*x^(7/2)*\[Kappa]p*\[CapitalSigma]^2)/(G^2*M^4) + 
      (3*S*x^(7/2)*\[Kappa]p*\[CapitalSigma]^2)/(G^3*M^6) + 
      (3*S*x^(7/2)*\[Delta]*\[Lambda]m*\[CapitalSigma]^2)/(G^3*M^6) - 
      (3*S*x^(7/2)*\[Lambda]p*\[CapitalSigma]^2)/(G^3*M^6) - 
      (2*x^2*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) + 
      (172*x^3*\[Nu]*\[CapitalSigma]^2)/(21*G^2*M^4) - 
      (4*Pi*x^(7/2)*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) - 
      (112*S*x^(7/2)*\[Nu]*\[CapitalSigma]^2)/(3*G^3*M^6) - 
      (89*x^3*\[Delta]*\[Kappa]m*\[Nu]*\[CapitalSigma]^2)/(42*G^2*M^4) - 
      (x^2*\[Kappa]p*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) + 
      (25*x^3*\[Kappa]p*\[Nu]*\[CapitalSigma]^2)/(6*G^2*M^4) - 
      (2*Pi*x^(7/2)*\[Kappa]p*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) - 
      (2*S*x^(7/2)*\[Kappa]p*\[Nu]*\[CapitalSigma]^2)/(3*G^3*M^6) + 
      (6*S*x^(7/2)*\[Lambda]p*\[Nu]*\[CapitalSigma]^2)/(G^3*M^6) - 
      (68*x^3*\[Nu]^2*\[CapitalSigma]^2)/(21*G^2*M^4) - 
      (34*x^3*\[Kappa]p*\[Nu]^2*\[CapitalSigma]^2)/(21*G^2*M^4) - 
      (5*x^(7/2)*\[Kappa]m*\[CapitalSigma]^3)/(3*G^3*M^6) + 
      (5*x^(7/2)*\[Delta]*\[Kappa]p*\[CapitalSigma]^3)/(3*G^3*M^6) + 
      (x^(7/2)*\[Lambda]m*\[CapitalSigma]^3)/(G^3*M^6) - 
      (x^(7/2)*\[Delta]*\[Lambda]p*\[CapitalSigma]^3)/(G^3*M^6) - 
      (20*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma]^3)/(3*G^3*M^6) + 
      (11*x^(7/2)*\[Kappa]m*\[Nu]*\[CapitalSigma]^3)/(3*G^3*M^6) - 
      (x^(7/2)*\[Delta]*\[Kappa]p*\[Nu]*\[CapitalSigma]^3)/(3*G^3*M^6) - 
      (3*x^(7/2)*\[Lambda]m*\[Nu]*\[CapitalSigma]^3)/(G^3*M^6) + 
      (x^(7/2)*\[Delta]*\[Lambda]p*\[Nu]*\[CapitalSigma]^3)/(G^3*M^6) - 
      (1712*x^3*Log[2])/105 - (428*x^3*Log[x])/105))/(c^2*E^((2*I)*\[Psi])*R)
 
h[3, 0] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((-2*I)/5)*Sqrt[6/7]*x^(5/2)*\[Nu] + 
      (((5017*I)/540)*x^(7/2)*\[Nu])/Sqrt[42] + (((5*I)/9)*x^(7/2)*\[Nu]^2)/
       Sqrt[42]))/(c^2*R)
 
h[3, 1] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((I/12)*Sqrt[x]*\[Delta])/Sqrt[14] - 
      (I/9)*Sqrt[2/7]*x^(3/2)*\[Delta] + (Sqrt[7/2]*x^2*\[Delta])/60 + 
      ((I/12)*Pi*x^2*\[Delta])/Sqrt[14] + ((I/24)*S*x^2*\[Delta])/
       (Sqrt[14]*G*M^2) + (((607*I)/2376)*x^(5/2)*\[Delta])/Sqrt[14] + 
      ((I/4)*S^2*x^(5/2)*\[Delta])/(Sqrt[14]*G^2*M^4) - 
      (Sqrt[14]*x^3*\[Delta])/45 - (I/9)*Sqrt[2/7]*Pi*x^3*\[Delta] - 
      (((79*I)/216)*S*x^3*\[Delta])/(Sqrt[14]*G*M^2) + 
      (((10753397*I)/18162144)*x^(7/2)*\[Delta])/Sqrt[14] - 
      (((13*I)/126)*EulerGamma*x^(7/2)*\[Delta])/Sqrt[14] + 
      (41*Pi*x^(7/2)*\[Delta])/(630*Sqrt[14]) + 
      ((I/72)*Pi^2*x^(7/2)*\[Delta])/Sqrt[14] - (47*S*x^(7/2)*\[Delta])/
       (240*Sqrt[14]*G*M^2) + ((I/24)*Pi*S*x^(7/2)*\[Delta])/
       (Sqrt[14]*G*M^2) - (((149*I)/72)*S^2*x^(7/2)*\[Delta])/
       (Sqrt[14]*G^2*M^4) - ((I/3)*S^2*x^(5/2)*\[Kappa]m)/
       (Sqrt[14]*G^2*M^4) + (((13*I)/16)*S^2*x^(7/2)*\[Kappa]m)/
       (Sqrt[14]*G^2*M^4) + ((I/8)*S^2*x^(5/2)*\[Delta]*\[Kappa]p)/
       (Sqrt[14]*G^2*M^4) - (((53*I)/144)*S^2*x^(7/2)*\[Delta]*\[Kappa]p)/
       (Sqrt[14]*G^2*M^4) - ((I/18)*x^(3/2)*\[Delta]*\[Nu])/Sqrt[14] - 
      ((17*I)/297)*Sqrt[2/7]*x^(5/2)*\[Delta]*\[Nu] + 
      (x^3*\[Delta]*\[Nu])/(180*Sqrt[14]) - (I/72)*Sqrt[7/2]*Pi*x^3*\[Delta]*
       \[Nu] + (((443*I)/216)*S*x^3*\[Delta]*\[Nu])/(Sqrt[14]*G*M^2) - 
      (((1738843*I)/1853280)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[14] + 
      (((41*I)/768)*Pi^2*x^(7/2)*\[Delta]*\[Nu])/Sqrt[14] - 
      (((11*I)/18)*S^2*x^(7/2)*\[Delta]*\[Nu])/(Sqrt[14]*G^2*M^4) - 
      (((2*I)/9)*Sqrt[2/7]*S^2*x^(7/2)*\[Kappa]m*\[Nu])/(G^2*M^4) - 
      (((11*I)/36)*S^2*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu])/(Sqrt[14]*G^2*M^4) - 
      (((247*I)/2376)*x^(5/2)*\[Delta]*\[Nu]^2)/Sqrt[14] + 
      (((327059*I)/370656)*x^(7/2)*\[Delta]*\[Nu]^2)/Sqrt[14] - 
      (((17525*I)/185328)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[14] + 
      (((5*I)/24)*x^2*\[CapitalSigma])/(Sqrt[14]*G*M^2) + 
      ((I/4)*S*x^(5/2)*\[CapitalSigma])/(Sqrt[14]*G^2*M^4) - 
      (((149*I)/216)*x^3*\[CapitalSigma])/(Sqrt[14]*G*M^2) + 
      (Sqrt[7/2]*x^(7/2)*\[CapitalSigma])/(24*G*M^2) + 
      (((5*I)/24)*Pi*x^(7/2)*\[CapitalSigma])/(Sqrt[14]*G*M^2) - 
      (((115*I)/36)*S*x^(7/2)*\[CapitalSigma])/(Sqrt[14]*G^2*M^4) - 
      (((11*I)/24)*S*x^(5/2)*\[Delta]*\[Kappa]m*\[CapitalSigma])/
       (Sqrt[14]*G^2*M^4) + (((85*I)/72)*S*x^(7/2)*\[Delta]*\[Kappa]m*
        \[CapitalSigma])/(Sqrt[14]*G^2*M^4) + 
      (((11*I)/24)*S*x^(5/2)*\[Kappa]p*\[CapitalSigma])/(Sqrt[14]*G^2*M^4) - 
      (((85*I)/72)*S*x^(7/2)*\[Kappa]p*\[CapitalSigma])/(Sqrt[14]*G^2*M^4) - 
      (((5*I)/8)*x^2*\[Nu]*\[CapitalSigma])/(Sqrt[14]*G*M^2) - 
      (I*S*x^(5/2)*\[Nu]*\[CapitalSigma])/(Sqrt[14]*G^2*M^4) + 
      (((25*I)/54)*Sqrt[7/2]*x^3*\[Nu]*\[CapitalSigma])/(G*M^2) - 
      (11*x^(7/2)*\[Nu]*\[CapitalSigma])/(240*Sqrt[14]*G*M^2) - 
      (((5*I)/8)*Pi*x^(7/2)*\[Nu]*\[CapitalSigma])/(Sqrt[14]*G*M^2) + 
      ((5*I)*Sqrt[2/7]*S*x^(7/2)*\[Nu]*\[CapitalSigma])/(G^2*M^4) - 
      (((5*I)/36)*S*x^(7/2)*\[Delta]*\[Kappa]m*\[Nu]*\[CapitalSigma])/
       (Sqrt[14]*G^2*M^4) - ((I/2)*S*x^(5/2)*\[Kappa]p*\[Nu]*\[CapitalSigma])/
       (Sqrt[14]*G^2*M^4) + (((29*I)/18)*S*x^(7/2)*\[Kappa]p*\[Nu]*
        \[CapitalSigma])/(Sqrt[14]*G^2*M^4) - 
      (((841*I)/216)*x^3*\[Nu]^2*\[CapitalSigma])/(Sqrt[14]*G*M^2) + 
      (((11*I)/9)*Sqrt[2/7]*S*x^(7/2)*\[Nu]^2*\[CapitalSigma])/(G^2*M^4) + 
      (((11*I)/9)*S*x^(7/2)*\[Kappa]p*\[Nu]^2*\[CapitalSigma])/
       (Sqrt[14]*G^2*M^4) - (((9*I)/8)*x^(7/2)*\[Delta]*\[CapitalSigma]^2)/
       (Sqrt[14]*G^2*M^4) - (((11*I)/48)*x^(5/2)*\[Kappa]m*\[CapitalSigma]^2)/
       (Sqrt[14]*G^2*M^4) + (((85*I)/144)*x^(7/2)*\[Kappa]m*
        \[CapitalSigma]^2)/(Sqrt[14]*G^2*M^4) + 
      (((11*I)/48)*x^(5/2)*\[Delta]*\[Kappa]p*\[CapitalSigma]^2)/
       (Sqrt[14]*G^2*M^4) - (((85*I)/144)*x^(7/2)*\[Delta]*\[Kappa]p*
        \[CapitalSigma]^2)/(Sqrt[14]*G^2*M^4) - 
      ((I/4)*x^(5/2)*\[Delta]*\[Nu]*\[CapitalSigma]^2)/(Sqrt[14]*G^2*M^4) + 
      (((25*I)/9)*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma]^2)/
       (Sqrt[14]*G^2*M^4) + ((I/12)*Sqrt[7/2]*x^(5/2)*\[Kappa]m*\[Nu]*
        \[CapitalSigma]^2)/(G^2*M^4) - (((233*I)/144)*x^(7/2)*\[Kappa]m*\[Nu]*
        \[CapitalSigma]^2)/(Sqrt[14]*G^2*M^4) - 
      ((I/8)*x^(5/2)*\[Delta]*\[Kappa]p*\[Nu]*\[CapitalSigma]^2)/
       (Sqrt[14]*G^2*M^4) + ((I/16)*Sqrt[7/2]*x^(7/2)*\[Delta]*\[Kappa]p*
        \[Nu]*\[CapitalSigma]^2)/(G^2*M^4) + 
      (((11*I)/18)*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma]^2)/
       (Sqrt[14]*G^2*M^4) - ((I/6)*x^(7/2)*\[Kappa]m*\[Nu]^2*
        \[CapitalSigma]^2)/(Sqrt[14]*G^2*M^4) + 
      (((11*I)/36)*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu]^2*\[CapitalSigma]^2)/
       (Sqrt[14]*G^2*M^4) + (x^2*\[Delta]*Log[2])/(6*Sqrt[14]) - 
      (2*Sqrt[2/7]*x^3*\[Delta]*Log[2])/9 - ((53*I)/315)*Sqrt[2/7]*x^(7/2)*
       \[Delta]*Log[2] + (Pi*x^(7/2)*\[Delta]*Log[2])/(6*Sqrt[14]) + 
      (S*x^(7/2)*\[Delta]*Log[2])/(12*Sqrt[14]*G*M^2) - 
      (Sqrt[7/2]*x^3*\[Delta]*\[Nu]*Log[2])/36 + 
      (5*x^(7/2)*\[CapitalSigma]*Log[2])/(12*Sqrt[14]*G*M^2) - 
      (5*x^(7/2)*\[Nu]*\[CapitalSigma]*Log[2])/(4*Sqrt[14]*G*M^2) - 
      ((I/6)*x^(7/2)*\[Delta]*Log[2]^2)/Sqrt[14] - 
      (((13*I)/252)*x^(7/2)*\[Delta]*Log[x])/Sqrt[14]))/(c^2*E^(I*\[Psi])*R)
 
h[3, 2] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((Sqrt[5/7]*x)/3 + 
      (2*Sqrt[5/7]*S*x^(3/2))/(3*G*M^2) - (193*x^2)/(54*Sqrt[35]) - 
      I*Sqrt[5/7]*x^(5/2) + (2*Sqrt[5/7]*Pi*x^(5/2))/3 - 
      (13*Sqrt[5/7]*S*x^(5/2))/(3*G*M^2) - (1451*x^3)/(2376*Sqrt[35]) - 
      (((2*I)/3)*Sqrt[5/7]*S*x^3)/(G*M^2) + (4*Sqrt[5/7]*Pi*S*x^3)/
       (3*G*M^2) - (8*Sqrt[5/7]*S^2*x^3)/(9*G^2*M^4) + 
      (((193*I)/18)*x^(7/2))/Sqrt[35] - (193*Pi*x^(7/2))/(27*Sqrt[35]) + 
      (4859*S*x^(7/2))/(396*Sqrt[35]*G*M^2) + (4*Sqrt[5/7]*S^3*x^(7/2))/
       (3*G^3*M^6) - (Sqrt[5/7]*S^2*x^3*\[Delta]*\[Kappa]m)/(3*G^2*M^4) + 
      (Sqrt[5/7]*S^2*x^3*\[Kappa]p)/(G^2*M^4) + 
      (2*Sqrt[5/7]*S^3*x^(7/2)*\[Kappa]p)/(3*G^3*M^6) - Sqrt[5/7]*x*\[Nu] + 
      (145*Sqrt[5/7]*x^2*\[Nu])/54 + ((22*I)*x^(5/2)*\[Nu])/Sqrt[35] - 
      2*Sqrt[5/7]*Pi*x^(5/2)*\[Nu] + (73*Sqrt[5/7]*S*x^(5/2)*\[Nu])/
       (9*G*M^2) - (17387*x^3*\[Nu])/(2376*Sqrt[35]) - 
      (4*Sqrt[5/7]*S^2*x^3*\[Nu])/(G^2*M^4) - 
      (((258929*I)/3240)*x^(7/2)*\[Nu])/Sqrt[35] + 
      (136*Sqrt[5/7]*Pi*x^(7/2)*\[Nu])/27 - (15413*Sqrt[5/7]*S*x^(7/2)*\[Nu])/
       (1188*G*M^2) - (2*Sqrt[5/7]*S^2*x^3*\[Kappa]p*\[Nu])/(G^2*M^4) - 
      (73*Sqrt[5/7]*x^2*\[Nu]^2)/54 + (5557*x^3*\[Nu]^2)/(132*Sqrt[35]) + 
      (((33751*I)/270)*x^(7/2)*\[Nu]^2)/Sqrt[35] - 
      (46*Sqrt[5/7]*Pi*x^(7/2)*\[Nu]^2)/27 - 
      (419*Sqrt[5/7]*S*x^(7/2)*\[Nu]^2)/(132*G*M^2) - 
      (763*Sqrt[7/5]*x^3*\[Nu]^3)/792 + (2*Sqrt[5/7]*x^(3/2)*\[Delta]*
        \[CapitalSigma])/(3*G*M^2) - (31*Sqrt[5/7]*x^(5/2)*\[Delta]*
        \[CapitalSigma])/(9*G*M^2) - ((2*I)*Sqrt[5/7]*x^3*\[Delta]*
        \[CapitalSigma])/(G*M^2) + (4*Sqrt[5/7]*Pi*x^3*\[Delta]*
        \[CapitalSigma])/(3*G*M^2) - (20*Sqrt[5/7]*S*x^3*\[Delta]*
        \[CapitalSigma])/(9*G^2*M^4) + (19241*x^(7/2)*\[Delta]*
        \[CapitalSigma])/(1188*Sqrt[35]*G*M^2) + 
      (8*Sqrt[5/7]*S^2*x^(7/2)*\[Delta]*\[CapitalSigma])/(3*G^3*M^6) - 
      (4*Sqrt[5/7]*S*x^3*\[Kappa]m*\[CapitalSigma])/(3*G^2*M^4) - 
      (2*Sqrt[5/7]*S^2*x^(7/2)*\[Kappa]m*\[CapitalSigma])/(3*G^3*M^6) + 
      (4*Sqrt[5/7]*S*x^3*\[Delta]*\[Kappa]p*\[CapitalSigma])/(3*G^2*M^4) + 
      (4*Sqrt[5/7]*S^2*x^(7/2)*\[Delta]*\[Kappa]p*\[CapitalSigma])/
       (3*G^3*M^6) + (10*Sqrt[5/7]*x^(5/2)*\[Delta]*\[Nu]*\[CapitalSigma])/
       (3*G*M^2) - (4*Sqrt[5/7]*S*x^3*\[Delta]*\[Nu]*\[CapitalSigma])/
       (G^2*M^4) - (1616*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma])/
       (33*Sqrt[35]*G*M^2) + (10*Sqrt[5/7]*S*x^3*\[Kappa]m*\[Nu]*
        \[CapitalSigma])/(3*G^2*M^4) - (2*Sqrt[5/7]*S*x^3*\[Delta]*\[Kappa]p*
        \[Nu]*\[CapitalSigma])/(G^2*M^4) - (16153*x^(7/2)*\[Delta]*\[Nu]^2*
        \[CapitalSigma])/(1188*Sqrt[35]*G*M^2) - 
      (4*Sqrt[5/7]*x^3*\[CapitalSigma]^2)/(3*G^2*M^4) + 
      (4*Sqrt[5/7]*S*x^(7/2)*\[CapitalSigma]^2)/(3*G^3*M^6) - 
      (2*Sqrt[5/7]*x^3*\[Delta]*\[Kappa]m*\[CapitalSigma]^2)/(3*G^2*M^4) - 
      (Sqrt[5/7]*S*x^(7/2)*\[Delta]*\[Kappa]m*\[CapitalSigma]^2)/(G^3*M^6) + 
      (2*Sqrt[5/7]*x^3*\[Kappa]p*\[CapitalSigma]^2)/(3*G^2*M^4) + 
      (Sqrt[5/7]*S*x^(7/2)*\[Kappa]p*\[CapitalSigma]^2)/(G^3*M^6) + 
      (8*Sqrt[5/7]*x^3*\[Nu]*\[CapitalSigma]^2)/(3*G^2*M^4) - 
      (20*Sqrt[5/7]*S*x^(7/2)*\[Nu]*\[CapitalSigma]^2)/(3*G^3*M^6) + 
      (4*Sqrt[5/7]*x^3*\[Delta]*\[Kappa]m*\[Nu]*\[CapitalSigma]^2)/
       (3*G^2*M^4) - (8*Sqrt[5/7]*x^3*\[Kappa]p*\[Nu]*\[CapitalSigma]^2)/
       (3*G^2*M^4) - (10*Sqrt[5/7]*S*x^(7/2)*\[Kappa]p*\[Nu]*
        \[CapitalSigma]^2)/(3*G^3*M^6) + (4*Sqrt[5/7]*x^3*\[Nu]^2*
        \[CapitalSigma]^2)/(G^2*M^4) + (2*Sqrt[5/7]*x^3*\[Kappa]p*\[Nu]^2*
        \[CapitalSigma]^2)/(G^2*M^4) - (Sqrt[5/7]*x^(7/2)*\[Kappa]m*
        \[CapitalSigma]^3)/(3*G^3*M^6) + (Sqrt[5/7]*x^(7/2)*\[Delta]*
        \[Kappa]p*\[CapitalSigma]^3)/(3*G^3*M^6) - 
      (4*Sqrt[5/7]*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma]^3)/(3*G^3*M^6) + 
      (4*Sqrt[5/7]*x^(7/2)*\[Kappa]m*\[Nu]*\[CapitalSigma]^3)/(3*G^3*M^6) - 
      (2*Sqrt[5/7]*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu]*\[CapitalSigma]^3)/
       (3*G^3*M^6)))/(c^2*E^((2*I)*\[Psi])*R)
 
h[3, 3] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((-3*I)/4)*Sqrt[15/14]*Sqrt[x]*
       \[Delta] + (3*I)*Sqrt[15/14]*x^(3/2)*\[Delta] - 
      (9*Sqrt[21/10]*x^2*\[Delta])/4 - ((9*I)/4)*Sqrt[15/14]*Pi*x^2*
       \[Delta] + (((3*I)/8)*Sqrt[105/2]*S*x^2*\[Delta])/(G*M^2) - 
      ((369*I)/88)*Sqrt[3/70]*x^(5/2)*\[Delta] - 
      (((9*I)/4)*Sqrt[15/14]*S^2*x^(5/2)*\[Delta])/(G^2*M^4) + 
      9*Sqrt[21/10]*x^3*\[Delta] + (9*I)*Sqrt[15/14]*Pi*x^3*\[Delta] - 
      (((139*I)/8)*Sqrt[3/70]*S*x^3*\[Delta])/(G*M^2) - 
      ((58164441*I)/224224)*Sqrt[3/70]*x^(7/2)*\[Delta] + 
      ((117*I)/14)*Sqrt[15/14]*EulerGamma*x^(7/2)*\[Delta] - 
      (369*Sqrt[3/70]*Pi*x^(7/2)*\[Delta])/14 - ((9*I)/8)*Sqrt[15/14]*Pi^2*
       x^(7/2)*\[Delta] + (639*Sqrt[3/70]*S*x^(7/2)*\[Delta])/(16*G*M^2) + 
      (((9*I)/8)*Sqrt[105/2]*Pi*S*x^(7/2)*\[Delta])/(G*M^2) + 
      (((69*I)/8)*Sqrt[15/14]*S^2*x^(7/2)*\[Delta])/(G^2*M^4) - 
      (((3*I)/16)*Sqrt[105/2]*S^2*x^(7/2)*\[Kappa]m)/(G^2*M^4) - 
      (((9*I)/8)*Sqrt[15/14]*S^2*x^(5/2)*\[Delta]*\[Kappa]p)/(G^2*M^4) + 
      (((45*I)/16)*Sqrt[15/14]*S^2*x^(7/2)*\[Delta]*\[Kappa]p)/(G^2*M^4) - 
      ((3*I)/2)*Sqrt[15/14]*x^(3/2)*\[Delta]*\[Nu] + 
      ((919*I)/22)*Sqrt[3/70]*x^(5/2)*\[Delta]*\[Nu] - 
      (48103*x^3*\[Delta]*\[Nu])/(108*Sqrt[210]) - ((27*I)/8)*Sqrt[15/14]*Pi*
       x^3*\[Delta]*\[Nu] + (((83*I)/8)*Sqrt[3/70]*S*x^3*\[Delta]*\[Nu])/
       (G*M^2) + ((7055*I)/4576)*Sqrt[15/14]*x^(7/2)*\[Delta]*\[Nu] - 
      ((123*I)/256)*Sqrt[15/14]*Pi^2*x^(7/2)*\[Delta]*\[Nu] - 
      (((9*I)/2)*Sqrt[15/14]*S^2*x^(7/2)*\[Delta]*\[Nu])/(G^2*M^4) + 
      ((3*I)*Sqrt[30/7]*S^2*x^(7/2)*\[Kappa]m*\[Nu])/(G^2*M^4) - 
      (((9*I)/4)*Sqrt[15/14]*S^2*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu])/
       (G^2*M^4) - ((887*I)/88)*Sqrt[3/70]*x^(5/2)*\[Delta]*\[Nu]^2 + 
      ((318841*I)/4576)*Sqrt[3/70]*x^(7/2)*\[Delta]*\[Nu]^2 - 
      ((24711*I)/2288)*Sqrt[3/70]*x^(7/2)*\[Delta]*\[Nu]^3 + 
      (((9*I)/8)*Sqrt[15/14]*x^2*\[CapitalSigma])/(G*M^2) - 
      (((9*I)/4)*Sqrt[15/14]*S*x^(5/2)*\[CapitalSigma])/(G^2*M^4) - 
      (((129*I)/8)*Sqrt[3/70]*x^3*\[CapitalSigma])/(G*M^2) + 
      (27*Sqrt[21/10]*x^(7/2)*\[CapitalSigma])/(8*G*M^2) + 
      (((27*I)/8)*Sqrt[15/14]*Pi*x^(7/2)*\[CapitalSigma])/(G*M^2) + 
      (((39*I)/4)*Sqrt[15/14]*S*x^(7/2)*\[CapitalSigma])/(G^2*M^4) + 
      (((9*I)/8)*Sqrt[15/14]*S*x^(5/2)*\[Delta]*\[Kappa]m*\[CapitalSigma])/
       (G^2*M^4) - (((33*I)/8)*Sqrt[15/14]*S*x^(7/2)*\[Delta]*\[Kappa]m*
        \[CapitalSigma])/(G^2*M^4) - (((9*I)/8)*Sqrt[15/14]*S*x^(5/2)*
        \[Kappa]p*\[CapitalSigma])/(G^2*M^4) + 
      (((33*I)/8)*Sqrt[15/14]*S*x^(7/2)*\[Kappa]p*\[CapitalSigma])/
       (G^2*M^4) - (((27*I)/8)*Sqrt[15/14]*x^2*\[Nu]*\[CapitalSigma])/
       (G*M^2) + ((9*I)*Sqrt[15/14]*S*x^(5/2)*\[Nu]*\[CapitalSigma])/
       (G^2*M^4) + ((9*I)*Sqrt[15/14]*x^3*\[Nu]*\[CapitalSigma])/(G*M^2) - 
      (8797*x^(7/2)*\[Nu]*\[CapitalSigma])/(48*Sqrt[210]*G*M^2) - 
      (((81*I)/8)*Sqrt[15/14]*Pi*x^(7/2)*\[Nu]*\[CapitalSigma])/(G*M^2) - 
      ((24*I)*Sqrt[30/7]*S*x^(7/2)*\[Nu]*\[CapitalSigma])/(G^2*M^4) + 
      (((33*I)/4)*Sqrt[15/14]*S*x^(7/2)*\[Delta]*\[Kappa]m*\[Nu]*
        \[CapitalSigma])/(G^2*M^4) + (((9*I)/2)*Sqrt[15/14]*S*x^(5/2)*
        \[Kappa]p*\[Nu]*\[CapitalSigma])/(G^2*M^4) - 
      (((39*I)/2)*Sqrt[15/14]*S*x^(7/2)*\[Kappa]p*\[Nu]*\[CapitalSigma])/
       (G^2*M^4) + (((15*I)/8)*Sqrt[15/14]*x^3*\[Nu]^2*\[CapitalSigma])/
       (G*M^2) + ((9*I)*Sqrt[30/7]*S*x^(7/2)*\[Nu]^2*\[CapitalSigma])/
       (G^2*M^4) + ((9*I)*Sqrt[15/14]*S*x^(7/2)*\[Kappa]p*\[Nu]^2*
        \[CapitalSigma])/(G^2*M^4) + (((9*I)/8)*Sqrt[15/14]*x^(7/2)*\[Delta]*
        \[CapitalSigma]^2)/(G^2*M^4) + (((9*I)/16)*Sqrt[15/14]*x^(5/2)*
        \[Kappa]m*\[CapitalSigma]^2)/(G^2*M^4) - 
      (((33*I)/16)*Sqrt[15/14]*x^(7/2)*\[Kappa]m*\[CapitalSigma]^2)/
       (G^2*M^4) - (((9*I)/16)*Sqrt[15/14]*x^(5/2)*\[Delta]*\[Kappa]p*
        \[CapitalSigma]^2)/(G^2*M^4) + (((33*I)/16)*Sqrt[15/14]*x^(7/2)*
        \[Delta]*\[Kappa]p*\[CapitalSigma]^2)/(G^2*M^4) + 
      (((9*I)/4)*Sqrt[15/14]*x^(5/2)*\[Delta]*\[Nu]*\[CapitalSigma]^2)/
       (G^2*M^4) - ((6*I)*Sqrt[30/7]*x^(7/2)*\[Delta]*\[Nu]*
        \[CapitalSigma]^2)/(G^2*M^4) - (((9*I)/4)*Sqrt[15/14]*x^(5/2)*
        \[Kappa]m*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) + 
      (((177*I)/16)*Sqrt[15/14]*x^(7/2)*\[Kappa]m*\[Nu]*\[CapitalSigma]^2)/
       (G^2*M^4) + (((9*I)/8)*Sqrt[15/14]*x^(5/2)*\[Delta]*\[Kappa]p*\[Nu]*
        \[CapitalSigma]^2)/(G^2*M^4) - (((111*I)/16)*Sqrt[15/14]*x^(7/2)*
        \[Delta]*\[Kappa]p*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) + 
      (((9*I)/2)*Sqrt[15/14]*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma]^2)/
       (G^2*M^4) - (((3*I)/2)*Sqrt[105/2]*x^(7/2)*\[Kappa]m*\[Nu]^2*
        \[CapitalSigma]^2)/(G^2*M^4) + (((9*I)/4)*Sqrt[15/14]*x^(7/2)*
        \[Delta]*\[Kappa]p*\[Nu]^2*\[CapitalSigma]^2)/(G^2*M^4) - 
      (9*Sqrt[15/14]*x^2*\[Delta]*Log[2])/2 + 9*Sqrt[30/7]*x^3*\[Delta]*
       Log[2] + ((477*I)/7)*Sqrt[6/35]*x^(7/2)*\[Delta]*Log[2] - 
      (27*Sqrt[15/14]*Pi*x^(7/2)*\[Delta]*Log[2])/2 + 
      (9*Sqrt[105/2]*S*x^(7/2)*\[Delta]*Log[2])/(4*G*M^2) - 
      (27*Sqrt[15/14]*x^3*\[Delta]*\[Nu]*Log[2])/4 + 
      (27*Sqrt[15/14]*x^(7/2)*\[CapitalSigma]*Log[2])/(4*G*M^2) - 
      (81*Sqrt[15/14]*x^(7/2)*\[Nu]*\[CapitalSigma]*Log[2])/(4*G*M^2) + 
      ((27*I)/2)*Sqrt[15/14]*x^(7/2)*\[Delta]*Log[2]^2 + 
      (9*Sqrt[15/14]*x^2*\[Delta]*Log[3])/2 - 9*Sqrt[30/7]*x^3*\[Delta]*
       Log[3] - ((369*I)/7)*Sqrt[3/70]*x^(7/2)*\[Delta]*Log[3] + 
      (27*Sqrt[15/14]*Pi*x^(7/2)*\[Delta]*Log[3])/2 - 
      (9*Sqrt[105/2]*S*x^(7/2)*\[Delta]*Log[3])/(4*G*M^2) + 
      (27*Sqrt[15/14]*x^3*\[Delta]*\[Nu]*Log[3])/4 - 
      (27*Sqrt[15/14]*x^(7/2)*\[CapitalSigma]*Log[3])/(4*G*M^2) + 
      (81*Sqrt[15/14]*x^(7/2)*\[Nu]*\[CapitalSigma]*Log[3])/(4*G*M^2) - 
      (27*I)*Sqrt[15/14]*x^(7/2)*\[Delta]*Log[2]*Log[3] + 
      ((27*I)/2)*Sqrt[15/14]*x^(7/2)*\[Delta]*Log[3]^2 + 
      ((117*I)/28)*Sqrt[15/14]*x^(7/2)*\[Delta]*Log[x]))/
    (c^2*E^((3*I)*\[Psi])*R)
 
h[4, 0] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(-1/504*1/Sqrt[2] + 
      (180101*x)/(14902272*Sqrt[2]) - (27227*x*\[Nu])/(532224*Sqrt[2])))/
    (c^2*R)
 
h[4, 1] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((I/84)*x^(3/2)*\[Delta])/Sqrt[10] + 
      ((I/168)*Sqrt[5/2]*S*x^2*\[Delta])/(G*M^2) - 
      (((101*I)/2772)*x^(5/2)*\[Delta])/Sqrt[10] + (4*Sqrt[2/5]*x^3*\[Delta])/
       315 + ((I/84)*Pi*x^3*\[Delta])/Sqrt[10] - 
      (((1147*I)/5544)*S*x^3*\[Delta])/(Sqrt[10]*G*M^2) + 
      (((21491*I)/630630)*x^(7/2)*\[Delta])/Sqrt[10] + 
      (53*S*x^(7/2)*\[Delta])/(1008*Sqrt[10]*G*M^2) + 
      ((I/168)*Sqrt[5/2]*Pi*S*x^(7/2)*\[Delta])/(G*M^2) - 
      ((I/56)*Sqrt[5/2]*S^2*x^(7/2)*\[Delta])/(G^2*M^4) - 
      ((I/24)*S^2*x^(7/2)*\[Kappa]m)/(Sqrt[10]*G^2*M^4) + 
      ((I/14)*S^2*x^(7/2)*\[Delta]*\[Kappa]p)/(Sqrt[10]*G^2*M^4) - 
      ((I/42)*x^(3/2)*\[Delta]*\[Nu])/Sqrt[10] + 
      (((337*I)/3696)*x^(5/2)*\[Delta]*\[Nu])/Sqrt[10] - 
      (1661*x^3*\[Delta]*\[Nu])/(2520*Sqrt[10]) - 
      ((I/42)*Pi*x^3*\[Delta]*\[Nu])/Sqrt[10] + 
      (((1139*I)/5544)*S*x^3*\[Delta]*\[Nu])/(Sqrt[10]*G*M^2) - 
      (((73427*I)/617760)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[10] - 
      ((I/42)*Sqrt[5/2]*S^2*x^(7/2)*\[Delta]*\[Nu])/(G^2*M^4) + 
      ((I/12)*S^2*x^(7/2)*\[Kappa]m*\[Nu])/(Sqrt[10]*G^2*M^4) - 
      ((I/84)*Sqrt[5/2]*S^2*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu])/(G^2*M^4) - 
      (((83*I)/2772)*x^(5/2)*\[Delta]*\[Nu]^2)/Sqrt[10] + 
      (((196957*I)/864864)*x^(7/2)*\[Delta]*\[Nu]^2)/Sqrt[10] - 
      ((239*I)/48048)*Sqrt[5/2]*x^(7/2)*\[Delta]*\[Nu]^3 + 
      ((I/168)*Sqrt[5/2]*x^2*\[CapitalSigma])/(G*M^2) - 
      (((103*I)/616)*x^3*\[CapitalSigma])/(Sqrt[10]*G*M^2) + 
      (2*Sqrt[2/5]*x^(7/2)*\[CapitalSigma])/(63*G*M^2) + 
      ((I/168)*Sqrt[5/2]*Pi*x^(7/2)*\[CapitalSigma])/(G*M^2) - 
      ((I/28)*Sqrt[5/2]*S*x^(7/2)*\[CapitalSigma])/(G^2*M^4) - 
      (((19*I)/168)*S*x^(7/2)*\[Delta]*\[Kappa]m*\[CapitalSigma])/
       (Sqrt[10]*G^2*M^4) + (((19*I)/168)*S*x^(7/2)*\[Kappa]p*
        \[CapitalSigma])/(Sqrt[10]*G^2*M^4) - 
      ((I/56)*Sqrt[5/2]*x^2*\[Nu]*\[CapitalSigma])/(G*M^2) + 
      (((29*I)/231)*Sqrt[5/2]*x^3*\[Nu]*\[CapitalSigma])/(G*M^2) - 
      (181*x^(7/2)*\[Nu]*\[CapitalSigma])/(1008*Sqrt[10]*G*M^2) - 
      ((I/56)*Sqrt[5/2]*Pi*x^(7/2)*\[Nu]*\[CapitalSigma])/(G*M^2) + 
      ((I/14)*Sqrt[5/2]*S*x^(7/2)*\[Nu]*\[CapitalSigma])/(G^2*M^4) + 
      ((I/7)*S*x^(7/2)*\[Delta]*\[Kappa]m*\[Nu]*\[CapitalSigma])/
       (Sqrt[10]*G^2*M^4) - (((3*I)/7)*S*x^(7/2)*\[Kappa]p*\[Nu]*
        \[CapitalSigma])/(Sqrt[10]*G^2*M^4) - 
      (((37*I)/616)*Sqrt[5/2]*x^3*\[Nu]^2*\[CapitalSigma])/(G*M^2) + 
      ((I/21)*Sqrt[10]*S*x^(7/2)*\[Nu]^2*\[CapitalSigma])/(G^2*M^4) + 
      ((I/21)*Sqrt[5/2]*S*x^(7/2)*\[Kappa]p*\[Nu]^2*\[CapitalSigma])/
       (G^2*M^4) - ((I/56)*Sqrt[5/2]*x^(7/2)*\[Delta]*\[CapitalSigma]^2)/
       (G^2*M^4) - (((19*I)/336)*x^(7/2)*\[Kappa]m*\[CapitalSigma]^2)/
       (Sqrt[10]*G^2*M^4) + (((19*I)/336)*x^(7/2)*\[Delta]*\[Kappa]p*
        \[CapitalSigma]^2)/(Sqrt[10]*G^2*M^4) + 
      ((I/42)*Sqrt[5/2]*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) + 
      (((43*I)/168)*x^(7/2)*\[Kappa]m*\[Nu]*\[CapitalSigma]^2)/
       (Sqrt[10]*G^2*M^4) - ((I/7)*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu]*
        \[CapitalSigma]^2)/(Sqrt[10]*G^2*M^4) + 
      ((I/42)*Sqrt[5/2]*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma]^2)/
       (G^2*M^4) - (((17*I)/84)*x^(7/2)*\[Kappa]m*\[Nu]^2*\[CapitalSigma]^2)/
       (Sqrt[10]*G^2*M^4) + ((I/84)*Sqrt[5/2]*x^(7/2)*\[Delta]*\[Kappa]p*
        \[Nu]^2*\[CapitalSigma]^2)/(G^2*M^4) + (x^3*\[Delta]*Log[2])/
       (42*Sqrt[10]) + (Sqrt[5/2]*S*x^(7/2)*\[Delta]*Log[2])/(84*G*M^2) - 
      (x^3*\[Delta]*\[Nu]*Log[2])/(21*Sqrt[10]) + 
      (Sqrt[5/2]*x^(7/2)*\[CapitalSigma]*Log[2])/(84*G*M^2) - 
      (Sqrt[5/2]*x^(7/2)*\[Nu]*\[CapitalSigma]*Log[2])/(28*G*M^2)))/
    (c^2*E^(I*\[Psi])*R)
 
h[4, 2] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((Sqrt[5]*x)/63 - 
      (437*x^2)/(1386*Sqrt[5]) - ((I/3)*x^(5/2))/Sqrt[5] + 
      (2*Sqrt[5]*Pi*x^(5/2))/63 - (4*S*x^(5/2))/(189*Sqrt[5]*G*M^2) + 
      (346013*x^3)/(840840*Sqrt[5]) + (4*Sqrt[5]*S^2*x^3)/(63*G^2*M^4) + 
      (((437*I)/330)*x^(7/2))/Sqrt[5] - (437*Pi*x^(7/2))/(693*Sqrt[5]) - 
      (86*S*x^(7/2))/(231*Sqrt[5]*G*M^2) - 
      (Sqrt[5]*S^2*x^3*\[Delta]*\[Kappa]m)/(21*G^2*M^4) + 
      (5*Sqrt[5]*S^2*x^3*\[Kappa]p)/(63*G^2*M^4) - (Sqrt[5]*x*\[Nu])/21 + 
      (115*Sqrt[5]*x^2*\[Nu])/594 + (((4*I)/3)*x^(5/2)*\[Nu])/Sqrt[5] - 
      (2*Sqrt[5]*Pi*x^(5/2)*\[Nu])/21 + (4*S*x^(5/2)*\[Nu])/
       (63*Sqrt[5]*G*M^2) - (606751*x^3*\[Nu])/(360360*Sqrt[5]) - 
      (4*Sqrt[5]*S^2*x^3*\[Nu])/(21*G^2*M^4) - 
      (((83029*I)/11088)*x^(7/2)*\[Nu])/Sqrt[5] + 
      (772*Sqrt[5]*Pi*x^(7/2)*\[Nu])/2079 + (6653*S*x^(7/2)*\[Nu])/
       (2079*Sqrt[5]*G*M^2) - (2*Sqrt[5]*S^2*x^3*\[Kappa]p*\[Nu])/
       (21*G^2*M^4) - (19*Sqrt[5]*x^2*\[Nu]^2)/1386 + 
      (400453*x^3*\[Nu]^2)/(324324*Sqrt[5]) + 
      (((31027*I)/4620)*x^(7/2)*\[Nu]^2)/Sqrt[5] + 
      (2*Sqrt[5]*Pi*x^(7/2)*\[Nu]^2)/99 - (1387*S*x^(7/2)*\[Nu]^2)/
       (231*Sqrt[5]*G*M^2) + (25783*x^3*\[Nu]^3)/(216216*Sqrt[5]) + 
      (4*x^(5/2)*\[Delta]*\[CapitalSigma])/(21*Sqrt[5]*G*M^2) + 
      (4*Sqrt[5]*S*x^3*\[Delta]*\[CapitalSigma])/(63*G^2*M^4) - 
      (626*x^(7/2)*\[Delta]*\[CapitalSigma])/(693*Sqrt[5]*G*M^2) - 
      (8*Sqrt[5]*S*x^3*\[Kappa]m*\[CapitalSigma])/(63*G^2*M^4) + 
      (8*Sqrt[5]*S*x^3*\[Delta]*\[Kappa]p*\[CapitalSigma])/(63*G^2*M^4) - 
      (8*x^(5/2)*\[Delta]*\[Nu]*\[CapitalSigma])/(21*Sqrt[5]*G*M^2) - 
      (4*Sqrt[5]*S*x^3*\[Delta]*\[Nu]*\[CapitalSigma])/(21*G^2*M^4) + 
      (6698*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma])/(2079*Sqrt[5]*G*M^2) + 
      (2*Sqrt[5]*S*x^3*\[Kappa]m*\[Nu]*\[CapitalSigma])/(7*G^2*M^4) - 
      (2*Sqrt[5]*S*x^3*\[Delta]*\[Kappa]p*\[Nu]*\[CapitalSigma])/
       (21*G^2*M^4) - (145*Sqrt[5]*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma])/
       (231*G*M^2) - (4*Sqrt[5]*x^3*\[Delta]*\[Kappa]m*\[CapitalSigma]^2)/
       (63*G^2*M^4) + (4*Sqrt[5]*x^3*\[Kappa]p*\[CapitalSigma]^2)/
       (63*G^2*M^4) - (4*Sqrt[5]*x^3*\[Nu]*\[CapitalSigma]^2)/(63*G^2*M^4) + 
      (2*Sqrt[5]*x^3*\[Delta]*\[Kappa]m*\[Nu]*\[CapitalSigma]^2)/
       (21*G^2*M^4) - (2*Sqrt[5]*x^3*\[Kappa]p*\[Nu]*\[CapitalSigma]^2)/
       (9*G^2*M^4) + (4*Sqrt[5]*x^3*\[Nu]^2*\[CapitalSigma]^2)/(21*G^2*M^4) + 
      (2*Sqrt[5]*x^3*\[Kappa]p*\[Nu]^2*\[CapitalSigma]^2)/(21*G^2*M^4)))/
    (c^2*E^((2*I)*\[Psi])*R)
 
h[4, 3] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((((-9*I)/4)*x^(3/2)*\[Delta])/Sqrt[70] - 
      (((9*I)/8)*Sqrt[5/14]*S*x^2*\[Delta])/(G*M^2) + 
      (((351*I)/44)*x^(5/2)*\[Delta])/Sqrt[70] - (36*Sqrt[2/35]*x^3*\[Delta])/
       5 - (((27*I)/4)*Pi*x^3*\[Delta])/Sqrt[70] + 
      (((3909*I)/88)*S*x^3*\[Delta])/(Sqrt[70]*G*M^2) - 
      (((32427*I)/10010)*x^(7/2)*\[Delta])/Sqrt[70] - 
      (477*S*x^(7/2)*\[Delta])/(16*Sqrt[70]*G*M^2) - 
      (((27*I)/8)*Sqrt[5/14]*Pi*S*x^(7/2)*\[Delta])/(G*M^2) + 
      (((27*I)/8)*Sqrt[5/14]*S^2*x^(7/2)*\[Delta])/(G^2*M^4) + 
      (((27*I)/8)*S^2*x^(7/2)*\[Kappa]m)/(Sqrt[70]*G^2*M^4) - 
      ((9*I)*S^2*x^(7/2)*\[Delta]*\[Kappa]p)/(Sqrt[70]*G^2*M^4) + 
      (((9*I)/2)*x^(3/2)*\[Delta]*\[Nu])/Sqrt[70] - ((543*I)/176)*Sqrt[7/10]*
       x^(5/2)*\[Delta]*\[Nu] + (16301*x^3*\[Delta]*\[Nu])/(360*Sqrt[70]) + 
      (((27*I)/2)*Pi*x^3*\[Delta]*\[Nu])/Sqrt[70] - 
      (((4353*I)/88)*S*x^3*\[Delta]*\[Nu])/(Sqrt[70]*G*M^2) + 
      (((745821*I)/22880)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[70] + 
      (((9*I)/2)*Sqrt[5/14]*S^2*x^(7/2)*\[Delta]*\[Nu])/(G^2*M^4) - 
      (((27*I)/4)*S^2*x^(7/2)*\[Kappa]m*\[Nu])/(Sqrt[70]*G^2*M^4) + 
      (((9*I)/4)*Sqrt[5/14]*S^2*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu])/(G^2*M^4) + 
      (((393*I)/44)*x^(5/2)*\[Delta]*\[Nu]^2)/Sqrt[70] - 
      ((44931*I)/4576)*Sqrt[7/10]*x^(7/2)*\[Delta]*\[Nu]^2 + 
      (((26883*I)/2288)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[70] - 
      (((9*I)/8)*Sqrt[5/14]*x^2*\[CapitalSigma])/(G*M^2) + 
      (((3249*I)/88)*x^3*\[CapitalSigma])/(Sqrt[70]*G*M^2) - 
      (18*Sqrt[2/35]*x^(7/2)*\[CapitalSigma])/(G*M^2) - 
      (((27*I)/8)*Sqrt[5/14]*Pi*x^(7/2)*\[CapitalSigma])/(G*M^2) + 
      (((27*I)/4)*Sqrt[5/14]*S*x^(7/2)*\[CapitalSigma])/(G^2*M^4) + 
      (((99*I)/8)*S*x^(7/2)*\[Delta]*\[Kappa]m*\[CapitalSigma])/
       (Sqrt[70]*G^2*M^4) - (((99*I)/8)*S*x^(7/2)*\[Kappa]p*\[CapitalSigma])/
       (Sqrt[70]*G^2*M^4) + (((27*I)/8)*Sqrt[5/14]*x^2*\[Nu]*\[CapitalSigma])/
       (G*M^2) - (((639*I)/22)*Sqrt[5/14]*x^3*\[Nu]*\[CapitalSigma])/
       (G*M^2) + (6007*x^(7/2)*\[Nu]*\[CapitalSigma])/(48*Sqrt[70]*G*M^2) + 
      (((81*I)/8)*Sqrt[5/14]*Pi*x^(7/2)*\[Nu]*\[CapitalSigma])/(G*M^2) - 
      (((27*I)/2)*Sqrt[5/14]*S*x^(7/2)*\[Nu]*\[CapitalSigma])/(G^2*M^4) - 
      ((9*I)*Sqrt[2/35]*S*x^(7/2)*\[Delta]*\[Kappa]m*\[Nu]*\[CapitalSigma])/
       (G^2*M^4) + ((27*I)*Sqrt[2/35]*S*x^(7/2)*\[Kappa]p*\[Nu]*
        \[CapitalSigma])/(G^2*M^4) + (((1467*I)/88)*Sqrt[5/14]*x^3*\[Nu]^2*
        \[CapitalSigma])/(G*M^2) - ((9*I)*Sqrt[10/7]*S*x^(7/2)*\[Nu]^2*
        \[CapitalSigma])/(G^2*M^4) - ((9*I)*Sqrt[5/14]*S*x^(7/2)*\[Kappa]p*
        \[Nu]^2*\[CapitalSigma])/(G^2*M^4) + 
      (((27*I)/8)*Sqrt[5/14]*x^(7/2)*\[Delta]*\[CapitalSigma]^2)/(G^2*M^4) + 
      (((99*I)/16)*x^(7/2)*\[Kappa]m*\[CapitalSigma]^2)/(Sqrt[70]*G^2*M^4) - 
      (((99*I)/16)*x^(7/2)*\[Delta]*\[Kappa]p*\[CapitalSigma]^2)/
       (Sqrt[70]*G^2*M^4) - (((9*I)/2)*Sqrt[5/14]*x^(7/2)*\[Delta]*\[Nu]*
        \[CapitalSigma]^2)/(G^2*M^4) - (((243*I)/8)*x^(7/2)*\[Kappa]m*\[Nu]*
        \[CapitalSigma]^2)/(Sqrt[70]*G^2*M^4) + 
      ((9*I)*Sqrt[2/35]*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu]*\[CapitalSigma]^2)/
       (G^2*M^4) - (((9*I)/2)*Sqrt[5/14]*x^(7/2)*\[Delta]*\[Nu]^2*
        \[CapitalSigma]^2)/(G^2*M^4) + (((117*I)/4)*x^(7/2)*\[Kappa]m*\[Nu]^2*
        \[CapitalSigma]^2)/(Sqrt[70]*G^2*M^4) - 
      (((9*I)/4)*Sqrt[5/14]*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu]^2*
        \[CapitalSigma]^2)/(G^2*M^4) - (27*x^3*\[Delta]*Log[2])/
       (2*Sqrt[70]) - (27*Sqrt[5/14]*S*x^(7/2)*\[Delta]*Log[2])/(4*G*M^2) + 
      (27*x^3*\[Delta]*\[Nu]*Log[2])/Sqrt[70] - 
      (27*Sqrt[5/14]*x^(7/2)*\[CapitalSigma]*Log[2])/(4*G*M^2) + 
      (81*Sqrt[5/14]*x^(7/2)*\[Nu]*\[CapitalSigma]*Log[2])/(4*G*M^2) + 
      (27*x^3*\[Delta]*Log[3])/(2*Sqrt[70]) + 
      (27*Sqrt[5/14]*S*x^(7/2)*\[Delta]*Log[3])/(4*G*M^2) - 
      (27*x^3*\[Delta]*\[Nu]*Log[3])/Sqrt[70] + 
      (27*Sqrt[5/14]*x^(7/2)*\[CapitalSigma]*Log[3])/(4*G*M^2) - 
      (81*Sqrt[5/14]*x^(7/2)*\[Nu]*\[CapitalSigma]*Log[3])/(4*G*M^2)))/
    (c^2*E^((3*I)*\[Psi])*R)
 
h[4, 4] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((-8*Sqrt[5/7]*x)/9 + 
      (2372*x^2)/(99*Sqrt[35]) + ((16*I)/3)*Sqrt[7/5]*x^(5/2) - 
      (32*Sqrt[5/7]*Pi*x^(5/2))/9 + (608*S*x^(5/2))/(27*Sqrt[35]*G*M^2) - 
      (1068671*x^3)/(45045*Sqrt[35]) - (32*Sqrt[5/7]*S^2*x^3)/(9*G^2*M^4) - 
      ((4744*I)/165)*Sqrt[7/5]*x^(7/2) + (9488*Pi*x^(7/2))/(99*Sqrt[35]) - 
      (6992*S*x^(7/2))/(99*Sqrt[35]*G*M^2) - (16*Sqrt[5/7]*S^2*x^3*\[Kappa]p)/
       (9*G^2*M^4) + (8*Sqrt[5/7]*x*\[Nu])/3 - (5092*Sqrt[5/7]*x^2*\[Nu])/
       297 - (((1193*I)/9)*x^(5/2)*\[Nu])/Sqrt[35] + 
      (32*Sqrt[5/7]*Pi*x^(5/2)*\[Nu])/3 - (608*S*x^(5/2)*\[Nu])/
       (9*Sqrt[35]*G*M^2) + (1088119*x^3*\[Nu])/(6435*Sqrt[35]) + 
      (32*Sqrt[5/7]*S^2*x^3*\[Nu])/(3*G^2*M^4) + 
      (((31525499*I)/31680)*x^(7/2)*\[Nu])/Sqrt[35] - 
      (19840*Sqrt[5/7]*Pi*x^(7/2)*\[Nu])/297 + (80504*S*x^(7/2)*\[Nu])/
       (297*Sqrt[35]*G*M^2) + (16*Sqrt[5/7]*S^2*x^3*\[Kappa]p*\[Nu])/
       (3*G^2*M^4) + (100*Sqrt[35]*x^2*\[Nu]^2)/99 - 
      (293758*x^3*\[Nu]^2)/(1053*Sqrt[35]) - 
      (((4096237*I)/4752)*x^(7/2)*\[Nu]^2)/Sqrt[35] + 
      (2272*Sqrt[5/7]*Pi*x^(7/2)*\[Nu]^2)/99 - (7768*S*x^(7/2)*\[Nu]^2)/
       (99*Sqrt[35]*G*M^2) + (226097*x^3*\[Nu]^3)/(3861*Sqrt[35]) + 
      (32*x^(5/2)*\[Delta]*\[CapitalSigma])/(3*Sqrt[35]*G*M^2) - 
      (32*Sqrt[5/7]*S*x^3*\[Delta]*\[CapitalSigma])/(9*G^2*M^4) - 
      (544*x^(7/2)*\[Delta]*\[CapitalSigma])/(11*Sqrt[35]*G*M^2) + 
      (16*Sqrt[5/7]*S*x^3*\[Kappa]m*\[CapitalSigma])/(9*G^2*M^4) - 
      (16*Sqrt[5/7]*S*x^3*\[Delta]*\[Kappa]p*\[CapitalSigma])/(9*G^2*M^4) - 
      (64*x^(5/2)*\[Delta]*\[Nu]*\[CapitalSigma])/(3*Sqrt[35]*G*M^2) + 
      (32*Sqrt[5/7]*S*x^3*\[Delta]*\[Nu]*\[CapitalSigma])/(3*G^2*M^4) + 
      (6928*Sqrt[5/7]*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma])/(297*G*M^2) - 
      (16*Sqrt[5/7]*S*x^3*\[Kappa]m*\[Nu]*\[CapitalSigma])/(3*G^2*M^4) + 
      (16*Sqrt[5/7]*S*x^3*\[Delta]*\[Kappa]p*\[Nu]*\[CapitalSigma])/
       (3*G^2*M^4) + (536*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma])/
       (99*Sqrt[35]*G*M^2) + (8*Sqrt[5/7]*x^3*\[Delta]*\[Kappa]m*
        \[CapitalSigma]^2)/(9*G^2*M^4) - (8*Sqrt[5/7]*x^3*\[Kappa]p*
        \[CapitalSigma]^2)/(9*G^2*M^4) + (32*Sqrt[5/7]*x^3*\[Nu]*
        \[CapitalSigma]^2)/(9*G^2*M^4) - (8*Sqrt[5/7]*x^3*\[Delta]*\[Kappa]m*
        \[Nu]*\[CapitalSigma]^2)/(3*G^2*M^4) + 
      (40*Sqrt[5/7]*x^3*\[Kappa]p*\[Nu]*\[CapitalSigma]^2)/(9*G^2*M^4) - 
      (32*Sqrt[5/7]*x^3*\[Nu]^2*\[CapitalSigma]^2)/(3*G^2*M^4) - 
      (16*Sqrt[5/7]*x^3*\[Kappa]p*\[Nu]^2*\[CapitalSigma]^2)/(3*G^2*M^4) - 
      ((64*I)/9)*Sqrt[5/7]*x^(5/2)*Log[2] + (((18976*I)/99)*x^(7/2)*Log[2])/
       Sqrt[35] + ((64*I)/3)*Sqrt[5/7]*x^(5/2)*\[Nu]*Log[2] - 
      ((39680*I)/297)*Sqrt[5/7]*x^(7/2)*\[Nu]*Log[2] + 
      ((4544*I)/99)*Sqrt[5/7]*x^(7/2)*\[Nu]^2*Log[2]))/
    (c^2*E^((4*I)*\[Psi])*R)
 
h[5, 0] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((((4117*I)/7560)*x^(7/2)*\[Nu])/
       Sqrt[462] - (((257*I)/90)*x^(7/2)*\[Nu]^2)/Sqrt[462]))/(c^2*R)
 
h[5, 1] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((I/288)*x^(3/2)*\[Delta])/Sqrt[385] - 
      (((179*I)/11232)*x^(5/2)*\[Delta])/Sqrt[385] + 
      (181*x^3*\[Delta])/(20160*Sqrt[385]) + ((I/288)*Pi*x^3*\[Delta])/
       Sqrt[385] + ((I/216)*S*x^3*\[Delta])/(Sqrt[385]*G*M^2) + 
      (((5023*I)/168480)*x^(7/2)*\[Delta])/Sqrt[385] + 
      ((I/288)*Sqrt[5/77]*S^2*x^(7/2)*\[Delta])/(G^2*M^4) - 
      ((I/48)*S^2*x^(7/2)*\[Kappa]m)/(Sqrt[385]*G^2*M^4) + 
      (((17*I)/576)*S^2*x^(7/2)*\[Delta]*\[Kappa]p)/(Sqrt[385]*G^2*M^4) - 
      ((I/144)*x^(3/2)*\[Delta]*\[Nu])/Sqrt[385] + (I/351)*Sqrt[11/35]*
       x^(5/2)*\[Delta]*\[Nu] - (313*x^3*\[Delta]*\[Nu])/(720*Sqrt[385]) - 
      ((I/144)*Pi*x^3*\[Delta]*\[Nu])/Sqrt[385] - 
      ((I/108)*S*x^3*\[Delta]*\[Nu])/(Sqrt[385]*G*M^2) - 
      (((49447*I)/673920)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[385] - 
      ((I/144)*Sqrt[5/77]*S^2*x^(7/2)*\[Delta]*\[Nu])/(G^2*M^4) + 
      ((I/24)*S^2*x^(7/2)*\[Kappa]m*\[Nu])/(Sqrt[385]*G^2*M^4) - 
      ((I/288)*Sqrt[5/77]*S^2*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu])/(G^2*M^4) - 
      ((I/2808)*x^(5/2)*\[Delta]*\[Nu]^2)/Sqrt[385] + 
      (((17*I)/648)*x^(7/2)*\[Delta]*\[Nu]^2)/Sqrt[385] + 
      ((41*I)/33696)*Sqrt[7/55]*x^(7/2)*\[Delta]*\[Nu]^3 + 
      ((I/432)*Sqrt[7/55]*x^3*\[CapitalSigma])/(G*M^2) + 
      ((I/288)*Sqrt[5/77]*S*x^(7/2)*\[CapitalSigma])/(G^2*M^4) - 
      (((29*I)/576)*S*x^(7/2)*\[Delta]*\[Kappa]m*\[CapitalSigma])/
       (Sqrt[385]*G^2*M^4) + (((29*I)/576)*S*x^(7/2)*\[Kappa]p*
        \[CapitalSigma])/(Sqrt[385]*G^2*M^4) - 
      ((I/432)*Sqrt[35/11]*x^3*\[Nu]*\[CapitalSigma])/(G*M^2) - 
      ((I/48)*Sqrt[5/77]*S*x^(7/2)*\[Nu]*\[CapitalSigma])/(G^2*M^4) + 
      (((17*I)/288)*S*x^(7/2)*\[Delta]*\[Kappa]m*\[Nu]*\[CapitalSigma])/
       (Sqrt[385]*G^2*M^4) - (((17*I)/96)*S*x^(7/2)*\[Kappa]p*\[Nu]*
        \[CapitalSigma])/(Sqrt[385]*G^2*M^4) + 
      ((I/432)*Sqrt[35/11]*x^3*\[Nu]^2*\[CapitalSigma])/(G*M^2) + 
      ((I/36)*Sqrt[5/77]*S*x^(7/2)*\[Nu]^2*\[CapitalSigma])/(G^2*M^4) + 
      ((I/72)*Sqrt[5/77]*S*x^(7/2)*\[Kappa]p*\[Nu]^2*\[CapitalSigma])/
       (G^2*M^4) - (((29*I)/1152)*x^(7/2)*\[Kappa]m*\[CapitalSigma]^2)/
       (Sqrt[385]*G^2*M^4) + (((29*I)/1152)*x^(7/2)*\[Delta]*\[Kappa]p*
        \[CapitalSigma]^2)/(Sqrt[385]*G^2*M^4) - 
      ((I/288)*Sqrt[5/77]*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma]^2)/
       (G^2*M^4) + ((I/64)*Sqrt[7/55]*x^(7/2)*\[Kappa]m*\[Nu]*
        \[CapitalSigma]^2)/(G^2*M^4) - (((17*I)/288)*x^(7/2)*\[Delta]*
        \[Kappa]p*\[Nu]*\[CapitalSigma]^2)/(Sqrt[385]*G^2*M^4) + 
      ((I/144)*Sqrt[5/77]*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma]^2)/
       (G^2*M^4) - ((I/144)*Sqrt[11/35]*x^(7/2)*\[Kappa]m*\[Nu]^2*
        \[CapitalSigma]^2)/(G^2*M^4) + ((I/288)*Sqrt[5/77]*x^(7/2)*\[Delta]*
        \[Kappa]p*\[Nu]^2*\[CapitalSigma]^2)/(G^2*M^4) + 
      (x^3*\[Delta]*Log[2])/(144*Sqrt[385]) - (x^3*\[Delta]*\[Nu]*Log[2])/
       (72*Sqrt[385])))/(c^2*E^(I*\[Psi])*R)
 
h[5, 2] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((2*x^2)/(27*Sqrt[55]) + 
      (2*S*x^(5/2))/(9*Sqrt[55]*G*M^2) - (3911*x^3)/(12285*Sqrt[55]) - 
      (((52*I)/135)*x^(7/2))/Sqrt[55] + (4*Pi*x^(7/2))/(27*Sqrt[55]) - 
      (71*S*x^(7/2))/(39*Sqrt[55]*G*M^2) - (2*Sqrt[5/11]*x^2*\[Nu])/27 - 
      (2*S*x^(5/2)*\[Nu])/(3*Sqrt[55]*G*M^2) + (3079*x^3*\[Nu])/
       (1755*Sqrt[55]) + (((16237*I)/4536)*x^(7/2)*\[Nu])/Sqrt[55] - 
      (4*Sqrt[5/11]*Pi*x^(7/2)*\[Nu])/27 + (229*Sqrt[11/5]*S*x^(7/2)*\[Nu])/
       (351*G*M^2) + (2*Sqrt[5/11]*x^2*\[Nu]^2)/27 - 
      (826*x^3*\[Nu]^2)/(351*Sqrt[55]) - (((1861*I)/270)*x^(7/2)*\[Nu]^2)/
       Sqrt[55] + (4*Sqrt[5/11]*Pi*x^(7/2)*\[Nu]^2)/27 - 
      (493*S*x^(7/2)*\[Nu]^2)/(117*Sqrt[55]*G*M^2) + 
      (7*Sqrt[11/5]*x^3*\[Nu]^3)/117 + (2*x^(5/2)*\[Delta]*\[CapitalSigma])/
       (9*Sqrt[55]*G*M^2) - (107*Sqrt[5/11]*x^(7/2)*\[Delta]*\[CapitalSigma])/
       (351*G*M^2) - (4*x^(5/2)*\[Delta]*\[Nu]*\[CapitalSigma])/
       (9*Sqrt[55]*G*M^2) + (488*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma])/
       (117*Sqrt[55]*G*M^2) - (113*Sqrt[5/11]*x^(7/2)*\[Delta]*\[Nu]^2*
        \[CapitalSigma])/(351*G*M^2)))/(c^2*E^((2*I)*\[Psi])*R)
 
h[5, 3] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((-9*I)/32)*Sqrt[3/110]*x^(3/2)*
       \[Delta] + ((621*I)/416)*Sqrt[3/110]*x^(5/2)*\[Delta] - 
      (4887*Sqrt[3/110]*x^3*\[Delta])/2240 - ((27*I)/32)*Sqrt[3/110]*Pi*x^3*
       \[Delta] + (((3*I)/8)*Sqrt[3/110]*S*x^3*\[Delta])/(G*M^2) - 
      ((3399*I)/14560)*Sqrt[33/10]*x^(7/2)*\[Delta] - 
      (((9*I)/32)*Sqrt[15/22]*S^2*x^(7/2)*\[Delta])/(G^2*M^4) + 
      (((9*I)/8)*Sqrt[3/110]*S^2*x^(7/2)*\[Kappa]m)/(G^2*M^4) - 
      (((117*I)/64)*Sqrt[3/110]*S^2*x^(7/2)*\[Delta]*\[Kappa]p)/(G^2*M^4) + 
      ((9*I)/16)*Sqrt[3/110]*x^(3/2)*\[Delta]*\[Nu] - 
      ((87*I)/26)*Sqrt[3/110]*x^(5/2)*\[Delta]*\[Nu] + 
      (41851*x^3*\[Delta]*\[Nu])/(2160*Sqrt[330]) + ((27*I)/16)*Sqrt[3/110]*
       Pi*x^3*\[Delta]*\[Nu] - (((3*I)/4)*Sqrt[3/110]*S*x^3*\[Delta]*\[Nu])/
       (G*M^2) + ((46611*I)/58240)*Sqrt[33/10]*x^(7/2)*\[Delta]*\[Nu] + 
      (((9*I)/16)*Sqrt[15/22]*S^2*x^(7/2)*\[Delta]*\[Nu])/(G^2*M^4) - 
      (((9*I)/4)*Sqrt[3/110]*S^2*x^(7/2)*\[Kappa]m*\[Nu])/(G^2*M^4) + 
      (((9*I)/32)*Sqrt[15/22]*S^2*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu])/
       (G^2*M^4) + ((3*I)/52)*Sqrt[33/10]*x^(5/2)*\[Delta]*\[Nu]^2 - 
      ((4887*I)/728)*Sqrt[3/110]*x^(7/2)*\[Delta]*\[Nu]^2 + 
      ((219*I)/2912)*Sqrt[15/22]*x^(7/2)*\[Delta]*\[Nu]^3 - 
      (((9*I)/16)*Sqrt[3/110]*x^3*\[CapitalSigma])/(G*M^2) - 
      (((9*I)/32)*Sqrt[15/22]*S*x^(7/2)*\[CapitalSigma])/(G^2*M^4) + 
      (((189*I)/64)*Sqrt[3/110]*S*x^(7/2)*\[Delta]*\[Kappa]m*\[CapitalSigma])/
       (G^2*M^4) - (((189*I)/64)*Sqrt[3/110]*S*x^(7/2)*\[Kappa]p*
        \[CapitalSigma])/(G^2*M^4) + (((9*I)/16)*Sqrt[15/22]*x^3*\[Nu]*
        \[CapitalSigma])/(G*M^2) + (((27*I)/16)*Sqrt[15/22]*S*x^(7/2)*\[Nu]*
        \[CapitalSigma])/(G^2*M^4) - (((117*I)/32)*Sqrt[3/110]*S*x^(7/2)*
        \[Delta]*\[Kappa]m*\[Nu]*\[CapitalSigma])/(G^2*M^4) + 
      (((351*I)/32)*Sqrt[3/110]*S*x^(7/2)*\[Kappa]p*\[Nu]*\[CapitalSigma])/
       (G^2*M^4) - (((9*I)/16)*Sqrt[15/22]*x^3*\[Nu]^2*\[CapitalSigma])/
       (G*M^2) - (((9*I)/4)*Sqrt[15/22]*S*x^(7/2)*\[Nu]^2*\[CapitalSigma])/
       (G^2*M^4) - (((9*I)/8)*Sqrt[15/22]*S*x^(7/2)*\[Kappa]p*\[Nu]^2*
        \[CapitalSigma])/(G^2*M^4) + (((189*I)/128)*Sqrt[3/110]*x^(7/2)*
        \[Kappa]m*\[CapitalSigma]^2)/(G^2*M^4) - 
      (((189*I)/128)*Sqrt[3/110]*x^(7/2)*\[Delta]*\[Kappa]p*
        \[CapitalSigma]^2)/(G^2*M^4) + (((9*I)/32)*Sqrt[15/22]*x^(7/2)*
        \[Delta]*\[Nu]*\[CapitalSigma]^2)/(G^2*M^4) - 
      (((423*I)/64)*Sqrt[3/110]*x^(7/2)*\[Kappa]m*\[Nu]*\[CapitalSigma]^2)/
       (G^2*M^4) + (((117*I)/32)*Sqrt[3/110]*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu]*
        \[CapitalSigma]^2)/(G^2*M^4) - (((9*I)/16)*Sqrt[15/22]*x^(7/2)*
        \[Delta]*\[Nu]^2*\[CapitalSigma]^2)/(G^2*M^4) + 
      (((81*I)/16)*Sqrt[3/110]*x^(7/2)*\[Kappa]m*\[Nu]^2*\[CapitalSigma]^2)/
       (G^2*M^4) - (((9*I)/32)*Sqrt[15/22]*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu]^2*
        \[CapitalSigma]^2)/(G^2*M^4) - (27*Sqrt[3/110]*x^3*\[Delta]*Log[2])/
       16 + (27*Sqrt[3/110]*x^3*\[Delta]*\[Nu]*Log[2])/8 + 
      (27*Sqrt[3/110]*x^3*\[Delta]*Log[3])/16 - 
      (27*Sqrt[3/110]*x^3*\[Delta]*\[Nu]*Log[3])/8))/(c^2*E^((3*I)*\[Psi])*R)
 
h[5, 4] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((-32*x^2)/(9*Sqrt[165]) - 
      (32*S*x^(5/2))/(3*Sqrt[165]*G*M^2) + (71216*x^3)/(4095*Sqrt[165]) + 
      (((1664*I)/45)*x^(7/2))/Sqrt[165] - (128*Pi*x^(7/2))/(9*Sqrt[165]) + 
      (3856*S*x^(7/2))/(39*Sqrt[165]*G*M^2) + (32*Sqrt[5/33]*x^2*\[Nu])/9 + 
      (32*S*x^(5/2)*\[Nu])/(Sqrt[165]*G*M^2) - (5264*Sqrt[11/15]*x^3*\[Nu])/
       585 - (((3351011*I)/15120)*x^(7/2)*\[Nu])/Sqrt[165] + 
      (128*Sqrt[5/33]*Pi*x^(7/2)*\[Nu])/9 - (47024*S*x^(7/2)*\[Nu])/
       (117*Sqrt[165]*G*M^2) - (32*Sqrt[5/33]*x^2*\[Nu]^2)/9 + 
      (16672*x^3*\[Nu]^2)/(117*Sqrt[165]) + ((331*I)/12)*Sqrt[11/15]*x^(7/2)*
       \[Nu]^2 - (128*Sqrt[5/33]*Pi*x^(7/2)*\[Nu]^2)/9 + 
      (3376*S*x^(7/2)*\[Nu]^2)/(13*Sqrt[165]*G*M^2) - 
      (1808*x^3*\[Nu]^3)/(39*Sqrt[165]) - 
      (32*x^(5/2)*\[Delta]*\[CapitalSigma])/(3*Sqrt[165]*G*M^2) + 
      (9904*x^(7/2)*\[Delta]*\[CapitalSigma])/(117*Sqrt[165]*G*M^2) + 
      (64*x^(5/2)*\[Delta]*\[Nu]*\[CapitalSigma])/(3*Sqrt[165]*G*M^2) - 
      (640*Sqrt[5/33]*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma])/(13*G*M^2) + 
      (13072*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma])/
       (117*Sqrt[165]*G*M^2) - (((256*I)/9)*x^(7/2)*Log[2])/Sqrt[165] + 
      ((256*I)/9)*Sqrt[5/33]*x^(7/2)*\[Nu]*Log[2] - ((256*I)/9)*Sqrt[5/33]*
       x^(7/2)*\[Nu]^2*Log[2]))/(c^2*E^((4*I)*\[Psi])*R)
 
h[5, 5] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((((625*I)/96)*x^(3/2)*\[Delta])/
       Sqrt[66] - (((164375*I)/3744)*x^(5/2)*\[Delta])/Sqrt[66] + 
      (113125*x^3*\[Delta])/(1344*Sqrt[66]) + (((3125*I)/96)*Pi*x^3*\[Delta])/
       Sqrt[66] - (((3125*I)/72)*S*x^3*\[Delta])/(Sqrt[66]*G*M^2) + 
      ((521875*I)/78624)*Sqrt[11/6]*x^(7/2)*\[Delta] + 
      (((3125*I)/96)*S^2*x^(7/2)*\[Delta])/(Sqrt[66]*G^2*M^4) + 
      (((3125*I)/192)*S^2*x^(7/2)*\[Delta]*\[Kappa]p)/(Sqrt[66]*G^2*M^4) - 
      (((625*I)/48)*x^(3/2)*\[Delta]*\[Nu])/Sqrt[66] + 
      (((26875*I)/234)*x^(5/2)*\[Delta]*\[Nu])/Sqrt[66] - 
      (17639*x^3*\[Delta]*\[Nu])/(80*Sqrt[66]) - 
      (((3125*I)/48)*Pi*x^3*\[Delta]*\[Nu])/Sqrt[66] + 
      (((3125*I)/36)*S*x^3*\[Delta]*\[Nu])/(Sqrt[66]*G*M^2) - 
      (((117978125*I)/314496)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[66] - 
      (((3125*I)/48)*S^2*x^(7/2)*\[Delta]*\[Nu])/(Sqrt[66]*G^2*M^4) - 
      (((3125*I)/96)*S^2*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu])/
       (Sqrt[66]*G^2*M^4) - ((2500*I)/117)*Sqrt[2/33]*x^(5/2)*\[Delta]*
       \[Nu]^2 + ((773125*I)/19656)*Sqrt[11/6]*x^(7/2)*\[Delta]*\[Nu]^2 - 
      (((6604375*I)/78624)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[66] - 
      (((3125*I)/144)*x^3*\[CapitalSigma])/(Sqrt[66]*G*M^2) + 
      (((3125*I)/96)*S*x^(7/2)*\[CapitalSigma])/(Sqrt[66]*G^2*M^4) - 
      (((3125*I)/192)*S*x^(7/2)*\[Delta]*\[Kappa]m*\[CapitalSigma])/
       (Sqrt[66]*G^2*M^4) + (((3125*I)/192)*S*x^(7/2)*\[Kappa]p*
        \[CapitalSigma])/(Sqrt[66]*G^2*M^4) + 
      (((15625*I)/144)*x^3*\[Nu]*\[CapitalSigma])/(Sqrt[66]*G*M^2) - 
      (((3125*I)/16)*S*x^(7/2)*\[Nu]*\[CapitalSigma])/(Sqrt[66]*G^2*M^4) + 
      (((3125*I)/96)*S*x^(7/2)*\[Delta]*\[Kappa]m*\[Nu]*\[CapitalSigma])/
       (Sqrt[66]*G^2*M^4) - (((3125*I)/32)*S*x^(7/2)*\[Kappa]p*\[Nu]*
        \[CapitalSigma])/(Sqrt[66]*G^2*M^4) - 
      (((15625*I)/144)*x^3*\[Nu]^2*\[CapitalSigma])/(Sqrt[66]*G*M^2) + 
      (((3125*I)/12)*S*x^(7/2)*\[Nu]^2*\[CapitalSigma])/(Sqrt[66]*G^2*M^4) + 
      (((3125*I)/24)*S*x^(7/2)*\[Kappa]p*\[Nu]^2*\[CapitalSigma])/
       (Sqrt[66]*G^2*M^4) - (((3125*I)/384)*x^(7/2)*\[Kappa]m*
        \[CapitalSigma]^2)/(Sqrt[66]*G^2*M^4) + 
      (((3125*I)/384)*x^(7/2)*\[Delta]*\[Kappa]p*\[CapitalSigma]^2)/
       (Sqrt[66]*G^2*M^4) - (((3125*I)/96)*x^(7/2)*\[Delta]*\[Nu]*
        \[CapitalSigma]^2)/(Sqrt[66]*G^2*M^4) + 
      (((3125*I)/64)*x^(7/2)*\[Kappa]m*\[Nu]*\[CapitalSigma]^2)/
       (Sqrt[66]*G^2*M^4) - (((3125*I)/96)*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu]*
        \[CapitalSigma]^2)/(Sqrt[66]*G^2*M^4) + 
      (((3125*I)/48)*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma]^2)/
       (Sqrt[66]*G^2*M^4) - (((3125*I)/48)*x^(7/2)*\[Kappa]m*\[Nu]^2*
        \[CapitalSigma]^2)/(Sqrt[66]*G^2*M^4) + 
      (((3125*I)/96)*x^(7/2)*\[Delta]*\[Kappa]p*\[Nu]^2*\[CapitalSigma]^2)/
       (Sqrt[66]*G^2*M^4) + (3125*x^3*\[Delta]*Log[2])/(48*Sqrt[66]) - 
      (3125*x^3*\[Delta]*\[Nu]*Log[2])/(24*Sqrt[66]) - 
      (3125*x^3*\[Delta]*Log[5])/(48*Sqrt[66]) + 
      (3125*x^3*\[Delta]*\[Nu]*Log[5])/(24*Sqrt[66])))/
    (c^2*E^((5*I)*\[Psi])*R)
 
h[6, 0] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((4195*x)/(1419264*Sqrt[273]) - 
      (215*x*\[Nu])/(16896*Sqrt[273])))/(c^2*R)
 
h[6, 1] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((I/8316)*x^(5/2)*\[Delta])/Sqrt[26] + 
      ((I/2376)*S*x^3*\[Delta])/(Sqrt[26]*G*M^2) - 
      (((125*I)/199584)*x^(7/2)*\[Delta])/Sqrt[26] - 
      ((I/2079)*x^(5/2)*\[Delta]*\[Nu])/Sqrt[26] - 
      ((I/1188)*S*x^3*\[Delta]*\[Nu])/(Sqrt[26]*G*M^2) + 
      (((277*I)/99792)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[26] + 
      ((I/2772)*x^(5/2)*\[Delta]*\[Nu]^2)/Sqrt[26] - 
      (((289*I)/99792)*x^(7/2)*\[Delta]*\[Nu]^2)/Sqrt[26] + 
      (((17*I)/24948)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[26] + 
      ((I/2376)*x^3*\[CapitalSigma])/(Sqrt[26]*G*M^2) - 
      (((5*I)/2376)*x^3*\[Nu]*\[CapitalSigma])/(Sqrt[26]*G*M^2) + 
      (((5*I)/2376)*x^3*\[Nu]^2*\[CapitalSigma])/(Sqrt[26]*G*M^2)))/
    (c^2*E^(I*\[Psi])*R)
 
h[6, 2] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((2*x^2)/(297*Sqrt[65]) - 
      (3*x^3)/(77*Sqrt[65]) - (((83*I)/2079)*x^(7/2))/Sqrt[65] + 
      (4*Pi*x^(7/2))/(297*Sqrt[65]) + (4*S*x^(7/2))/(693*Sqrt[65]*G*M^2) - 
      (2*Sqrt[5/13]*x^2*\[Nu])/297 + (59*x^3*\[Nu])/(297*Sqrt[65]) + 
      (((799789*I)/1995840)*x^(7/2)*\[Nu])/Sqrt[65] - 
      (4*Sqrt[5/13]*Pi*x^(7/2)*\[Nu])/297 - (4*Sqrt[5/13]*S*x^(7/2)*\[Nu])/
       (693*G*M^2) + (2*Sqrt[5/13]*x^2*\[Nu]^2)/297 - 
      (64*x^3*\[Nu]^2)/(297*Sqrt[65]) - (((19193*I)/23760)*x^(7/2)*\[Nu]^2)/
       Sqrt[65] + (4*Sqrt[5/13]*Pi*x^(7/2)*\[Nu]^2)/297 + 
      (4*Sqrt[5/13]*S*x^(7/2)*\[Nu]^2)/(693*G*M^2) + 
      (7*x^3*\[Nu]^3)/(297*Sqrt[65]) + (68*x^(7/2)*\[Delta]*\[CapitalSigma])/
       (2079*Sqrt[65]*G*M^2) - (272*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma])/
       (2079*Sqrt[65]*G*M^2) + (68*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma])/
       (693*Sqrt[65]*G*M^2)))/(c^2*E^((2*I)*\[Psi])*R)
 
h[6, 3] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((((-81*I)/616)*x^(5/2)*\[Delta])/
       Sqrt[65] - (((81*I)/176)*S*x^3*\[Delta])/(Sqrt[65]*G*M^2) + 
      (((513*I)/704)*x^(7/2)*\[Delta])/Sqrt[65] + 
      (((81*I)/154)*x^(5/2)*\[Delta]*\[Nu])/Sqrt[65] + 
      (((81*I)/88)*S*x^3*\[Delta]*\[Nu])/(Sqrt[65]*G*M^2) - 
      (((1161*I)/352)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[65] - 
      (((243*I)/616)*x^(5/2)*\[Delta]*\[Nu]^2)/Sqrt[65] + 
      (((1269*I)/352)*x^(7/2)*\[Delta]*\[Nu]^2)/Sqrt[65] - 
      (((81*I)/88)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[65] - 
      (((81*I)/176)*x^3*\[CapitalSigma])/(Sqrt[65]*G*M^2) + 
      (((81*I)/176)*Sqrt[5/13]*x^3*\[Nu]*\[CapitalSigma])/(G*M^2) - 
      (((81*I)/176)*Sqrt[5/13]*x^3*\[Nu]^2*\[CapitalSigma])/(G*M^2)))/
    (c^2*E^((3*I)*\[Psi])*R)
 
h[6, 4] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((-128*Sqrt[2/39]*x^2)/495 + 
      (1984*Sqrt[2/39]*x^3)/1155 + ((10624*I)/3465)*Sqrt[2/39]*x^(7/2) - 
      (512*Sqrt[2/39]*Pi*x^(7/2))/495 + (256*Sqrt[2/39]*S*x^(7/2))/
       (385*G*M^2) + (128*Sqrt[2/39]*x^2*\[Nu])/99 - 
      (4544*Sqrt[2/39]*x^3*\[Nu])/495 - (((686443*I)/19008)*x^(7/2)*\[Nu])/
       Sqrt[78] + (512*Sqrt[2/39]*Pi*x^(7/2)*\[Nu])/99 - 
      (256*Sqrt[2/39]*S*x^(7/2)*\[Nu])/(77*G*M^2) - 
      (128*Sqrt[2/39]*x^2*\[Nu]^2)/99 + (512*Sqrt[2/39]*x^3*\[Nu]^2)/45 + 
      (((8497*I)/176)*x^(7/2)*\[Nu]^2)/Sqrt[78] - 
      (512*Sqrt[2/39]*Pi*x^(7/2)*\[Nu]^2)/99 + 
      (256*Sqrt[2/39]*S*x^(7/2)*\[Nu]^2)/(77*G*M^2) - 
      (1216*Sqrt[2/39]*x^3*\[Nu]^3)/495 - (256*Sqrt[2/39]*x^(7/2)*\[Delta]*
        \[CapitalSigma])/(693*G*M^2) + (1024*Sqrt[2/39]*x^(7/2)*\[Delta]*
        \[Nu]*\[CapitalSigma])/(693*G*M^2) - 
      (256*Sqrt[2/39]*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma])/(231*G*M^2) - 
      ((1024*I)/495)*Sqrt[2/39]*x^(7/2)*Log[2] + ((1024*I)/99)*Sqrt[2/39]*
       x^(7/2)*\[Nu]*Log[2] - ((1024*I)/99)*Sqrt[2/39]*x^(7/2)*\[Nu]^2*
       Log[2]))/(c^2*E^((4*I)*\[Psi])*R)
 
h[6, 5] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((((3125*I)/504)*x^(5/2)*\[Delta])/
       Sqrt[429] + (((3125*I)/144)*S*x^3*\[Delta])/(Sqrt[429]*G*M^2) - 
      (((465625*I)/12096)*x^(7/2)*\[Delta])/Sqrt[429] - 
      (((3125*I)/126)*x^(5/2)*\[Delta]*\[Nu])/Sqrt[429] - 
      (((3125*I)/72)*S*x^3*\[Delta]*\[Nu])/(Sqrt[429]*G*M^2) + 
      (((1090625*I)/6048)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[429] + 
      (((3125*I)/168)*x^(5/2)*\[Delta]*\[Nu]^2)/Sqrt[429] - 
      (((1278125*I)/6048)*x^(7/2)*\[Delta]*\[Nu]^2)/Sqrt[429] + 
      (((90625*I)/1512)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[429] + 
      (((3125*I)/144)*x^3*\[CapitalSigma])/(Sqrt[429]*G*M^2) - 
      (((15625*I)/144)*x^3*\[Nu]*\[CapitalSigma])/(Sqrt[429]*G*M^2) + 
      (((15625*I)/144)*x^3*\[Nu]^2*\[CapitalSigma])/(Sqrt[429]*G*M^2)))/
    (c^2*E^((5*I)*\[Psi])*R)
 
h[6, 6] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((54*x^2)/(5*Sqrt[143]) - 
      (3051*x^3)/(35*Sqrt[143]) - (((6723*I)/35)*x^(7/2))/Sqrt[143] + 
      (324*Pi*x^(7/2))/(5*Sqrt[143]) - (3132*S*x^(7/2))/
       (35*Sqrt[143]*G*M^2) - (54*x^2*\[Nu])/Sqrt[143] + 
      (189*Sqrt[13/11]*x^3*\[Nu])/5 + (((21787499*I)/20160)*x^(7/2)*\[Nu])/
       Sqrt[143] - (324*Pi*x^(7/2)*\[Nu])/Sqrt[143] + 
      (3132*S*x^(7/2)*\[Nu])/(7*Sqrt[143]*G*M^2) + (54*x^2*\[Nu]^2)/
       Sqrt[143] - (3456*x^3*\[Nu]^2)/(5*Sqrt[143]) - 
      (((323903*I)/240)*x^(7/2)*\[Nu]^2)/Sqrt[143] + (324*Pi*x^(7/2)*\[Nu]^2)/
       Sqrt[143] - (3132*S*x^(7/2)*\[Nu]^2)/(7*Sqrt[143]*G*M^2) + 
      (81*Sqrt[13/11]*x^3*\[Nu]^3)/5 - (324*x^(7/2)*\[Delta]*\[CapitalSigma])/
       (7*Sqrt[143]*G*M^2) + (1296*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma])/
       (7*Sqrt[143]*G*M^2) - (972*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma])/
       (7*Sqrt[143]*G*M^2) + (((648*I)/5)*x^(7/2)*Log[3])/Sqrt[143] - 
      ((648*I)*x^(7/2)*\[Nu]*Log[3])/Sqrt[143] + 
      ((648*I)*x^(7/2)*\[Nu]^2*Log[3])/Sqrt[143]))/(c^2*E^((6*I)*\[Psi])*R)
 
h[7, 1] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((I/864864)*x^(5/2)*\[Delta])/Sqrt[2] - 
      (((223*I)/29405376)*x^(7/2)*\[Delta])/Sqrt[2] - 
      ((I/216216)*x^(5/2)*\[Delta]*\[Nu])/Sqrt[2] + 
      (((1361*I)/44108064)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[2] + 
      ((I/288288)*x^(5/2)*\[Delta]*\[Nu]^2)/Sqrt[2] - 
      (((43*I)/1696464)*x^(7/2)*\[Delta]*\[Nu]^2)/Sqrt[2] + 
      (((19*I)/7351344)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[2]))/
    (c^2*E^(I*\[Psi])*R)
 
h[7, 2] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(x^3/(3003*Sqrt[3]) + 
      (4*S*x^(7/2))/(3003*Sqrt[3]*G*M^2) - (x^3*\[Nu])/(429*Sqrt[3]) - 
      (20*S*x^(7/2)*\[Nu])/(3003*Sqrt[3]*G*M^2) + (2*x^3*\[Nu]^2)/
       (429*Sqrt[3]) + (20*S*x^(7/2)*\[Nu]^2)/(3003*Sqrt[3]*G*M^2) - 
      (x^3*\[Nu]^3)/(429*Sqrt[3]) + (4*x^(7/2)*\[Delta]*\[CapitalSigma])/
       (3003*Sqrt[3]*G*M^2) - (16*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma])/
       (3003*Sqrt[3]*G*M^2) + (4*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma])/
       (1001*Sqrt[3]*G*M^2)))/(c^2*E^((2*I)*\[Psi])*R)
 
h[7, 3] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((-243*I)/160160)*Sqrt[3/2]*x^(5/2)*
       \[Delta] + ((58077*I)/5445440)*Sqrt[3/2]*x^(7/2)*\[Delta] + 
      ((243*I)/40040)*Sqrt[3/2]*x^(5/2)*\[Delta]*\[Nu] - 
      ((3483*I)/77792)*Sqrt[3/2]*x^(7/2)*\[Delta]*\[Nu] - 
      ((729*I)/160160)*Sqrt[3/2]*x^(5/2)*\[Delta]*\[Nu]^2 + 
      ((7857*I)/194480)*Sqrt[3/2]*x^(7/2)*\[Delta]*\[Nu]^2 - 
      ((243*I)/38896)*Sqrt[3/2]*x^(7/2)*\[Delta]*\[Nu]^3))/
    (c^2*E^((3*I)*\[Psi])*R)
 
h[7, 4] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((-128*Sqrt[2/33]*x^3)/1365 - 
      (512*Sqrt[2/33]*S*x^(7/2))/(1365*G*M^2) + (128*Sqrt[2/33]*x^3*\[Nu])/
       195 + (512*Sqrt[2/33]*S*x^(7/2)*\[Nu])/(273*G*M^2) - 
      (256*Sqrt[2/33]*x^3*\[Nu]^2)/195 - (512*Sqrt[2/33]*S*x^(7/2)*\[Nu]^2)/
       (273*G*M^2) + (128*Sqrt[2/33]*x^3*\[Nu]^3)/195 - 
      (512*Sqrt[2/33]*x^(7/2)*\[Delta]*\[CapitalSigma])/(1365*G*M^2) + 
      (2048*Sqrt[2/33]*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma])/(1365*G*M^2) - 
      (512*Sqrt[2/33]*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma])/(455*G*M^2)))/
    (c^2*E^((4*I)*\[Psi])*R)
 
h[7, 5] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((((15625*I)/26208)*x^(5/2)*\[Delta])/
       Sqrt[66] - (((4234375*I)/891072)*x^(7/2)*\[Delta])/Sqrt[66] - 
      (((15625*I)/6552)*x^(5/2)*\[Delta]*\[Nu])/Sqrt[66] + 
      ((2546875*I)/1336608)*Sqrt[11/6]*x^(7/2)*\[Delta]*\[Nu] + 
      (((15625*I)/8736)*x^(5/2)*\[Delta]*\[Nu]^2)/Sqrt[66] - 
      (((14359375*I)/668304)*x^(7/2)*\[Delta]*\[Nu]^2)/Sqrt[66] + 
      (((1046875*I)/222768)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[66]))/
    (c^2*E^((5*I)*\[Psi])*R)
 
h[7, 6] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((81*Sqrt[3/143]*x^3)/35 + 
      (324*Sqrt[3/143]*S*x^(7/2))/(35*G*M^2) - (81*Sqrt[3/143]*x^3*\[Nu])/5 - 
      (324*Sqrt[3/143]*S*x^(7/2)*\[Nu])/(7*G*M^2) + 
      (162*Sqrt[3/143]*x^3*\[Nu]^2)/5 + (324*Sqrt[3/143]*S*x^(7/2)*\[Nu]^2)/
       (7*G*M^2) - (81*Sqrt[3/143]*x^3*\[Nu]^3)/5 + 
      (324*Sqrt[3/143]*x^(7/2)*\[Delta]*\[CapitalSigma])/(35*G*M^2) - 
      (1296*Sqrt[3/143]*x^(7/2)*\[Delta]*\[Nu]*\[CapitalSigma])/(35*G*M^2) + 
      (972*Sqrt[3/143]*x^(7/2)*\[Delta]*\[Nu]^2*\[CapitalSigma])/(35*G*M^2)))/
    (c^2*E^((6*I)*\[Psi])*R)
 
h[7, 7] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((-16807*I)/1440)*Sqrt[7/858]*x^(5/2)*
       \[Delta] + ((487403*I)/48960)*Sqrt[77/78]*x^(7/2)*\[Delta] + 
      ((16807*I)/360)*Sqrt[7/858]*x^(5/2)*\[Delta]*\[Nu] - 
      ((7479115*I)/14688)*Sqrt[7/858]*x^(7/2)*\[Delta]*\[Nu] - 
      ((16807*I)/480)*Sqrt[7/858]*x^(5/2)*\[Delta]*\[Nu]^2 + 
      ((21496153*I)/36720)*Sqrt[7/858]*x^(7/2)*\[Delta]*\[Nu]^2 - 
      ((386561*I)/2448)*Sqrt[7/858]*x^(7/2)*\[Delta]*\[Nu]^3))/
    (c^2*E^((7*I)*\[Psi])*R)
 
h[8, 1] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((I/741312)*x^(7/2)*\[Delta])/
       Sqrt[238] - ((I/123552)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[238] + 
      (((5*I)/370656)*x^(7/2)*\[Delta]*\[Nu]^2)/Sqrt[238] - 
      ((I/185328)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[238]))/(c^2*E^(I*\[Psi])*R)
 
h[8, 2] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(x^3/(9009*Sqrt[85]) - 
      (x^3*\[Nu])/(1287*Sqrt[85]) + (2*x^3*\[Nu]^2)/(1287*Sqrt[85]) - 
      (x^3*\[Nu]^3)/(1287*Sqrt[85])))/(c^2*E^((2*I)*\[Psi])*R)
 
h[8, 3] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((-81*I)/5824)*Sqrt[3/1870]*x^(7/2)*
       \[Delta] + ((243*I)/2912)*Sqrt[3/1870]*x^(7/2)*\[Delta]*\[Nu] - 
      ((81*I)/2912)*Sqrt[15/374]*x^(7/2)*\[Delta]*\[Nu]^2 + 
      ((81*I)/1456)*Sqrt[3/1870]*x^(7/2)*\[Delta]*\[Nu]^3))/
    (c^2*E^((3*I)*\[Psi])*R)
 
h[8, 4] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((-128*Sqrt[2/187]*x^3)/4095 + 
      (128*Sqrt[2/187]*x^3*\[Nu])/585 - (256*Sqrt[2/187]*x^3*\[Nu]^2)/585 + 
      (128*Sqrt[2/187]*x^3*\[Nu]^3)/585))/(c^2*E^((4*I)*\[Psi])*R)
 
h[8, 5] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((((78125*I)/36288)*x^(7/2)*\[Delta])/
       Sqrt[4862] - (((78125*I)/6048)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[4862] + 
      (((390625*I)/18144)*x^(7/2)*\[Delta]*\[Nu]^2)/Sqrt[4862] - 
      (((78125*I)/9072)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[4862]))/
    (c^2*E^((5*I)*\[Psi])*R)
 
h[8, 6] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((243*Sqrt[3/17017]*x^3)/35 - 
      (243*Sqrt[3/17017]*x^3*\[Nu])/5 + (486*Sqrt[3/17017]*x^3*\[Nu]^2)/5 - 
      (243*Sqrt[3/17017]*x^3*\[Nu]^3)/5))/(c^2*E^((6*I)*\[Psi])*R)
 
h[8, 7] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((-117649*I)/5184)*Sqrt[7/24310]*x^(7/2)*
       \[Delta] + ((117649*I)/864)*Sqrt[7/24310]*x^(7/2)*\[Delta]*\[Nu] - 
      ((117649*I)/2592)*Sqrt[35/4862]*x^(7/2)*\[Delta]*\[Nu]^2 + 
      ((117649*I)/1296)*Sqrt[7/24310]*x^(7/2)*\[Delta]*\[Nu]^3))/
    (c^2*E^((7*I)*\[Psi])*R)
 
h[8, 8] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((-16384*Sqrt[2/85085]*x^3)/63 + 
      (16384*Sqrt[2/85085]*x^3*\[Nu])/9 - (32768*Sqrt[2/85085]*x^3*\[Nu]^2)/
       9 + (16384*Sqrt[2/85085]*x^3*\[Nu]^3)/9))/(c^2*E^((8*I)*\[Psi])*R)
 
h[9, 1] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((I/9165312)*x^(7/2)*\[Delta])/
       Sqrt[2090] - ((I/1527552)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[2090] + 
      (I/4582656)*Sqrt[5/418]*x^(7/2)*\[Delta]*\[Nu]^2 - 
      ((I/2291328)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[2090]))/(c^2*E^(I*\[Psi])*R)
 
h[9, 3] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*(((-81*I)/113152)*Sqrt[3/665]*x^(7/2)*
       \[Delta] + ((243*I)/56576)*Sqrt[3/665]*x^(7/2)*\[Delta]*\[Nu] - 
      ((81*I)/56576)*Sqrt[15/133]*x^(7/2)*\[Delta]*\[Nu]^2 + 
      ((81*I)/28288)*Sqrt[3/665]*x^(7/2)*\[Delta]*\[Nu]^3))/
    (c^2*E^((3*I)*\[Psi])*R)
 
h[9, 5] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((((390625*I)/4935168)*x^(7/2)*\[Delta])/
       Sqrt[247] - (((390625*I)/822528)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[247] + 
      (((1953125*I)/2467584)*x^(7/2)*\[Delta]*\[Nu]^2)/Sqrt[247] - 
      (((390625*I)/1233792)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[247]))/
    (c^2*E^((5*I)*\[Psi])*R)
 
h[9, 7] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*
     ((((-5764801*I)/1410048)*x^(7/2)*\[Delta])/Sqrt[1235] + 
      (((5764801*I)/235008)*x^(7/2)*\[Delta]*\[Nu])/Sqrt[1235] - 
      ((5764801*I)/705024)*Sqrt[5/247]*x^(7/2)*\[Delta]*\[Nu]^2 + 
      (((5764801*I)/352512)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[1235]))/
    (c^2*E^((7*I)*\[Psi])*R)
 
h[9, 9] = (8*G*M*Sqrt[Pi/5]*x*\[Nu]*((((1594323*I)/7168)*x^(7/2)*\[Delta])/
       Sqrt[20995] - (((4782969*I)/3584)*x^(7/2)*\[Delta]*\[Nu])/
       Sqrt[20995] + ((1594323*I)/3584)*Sqrt[5/4199]*x^(7/2)*\[Delta]*
       \[Nu]^2 - (((1594323*I)/1792)*x^(7/2)*\[Delta]*\[Nu]^3)/Sqrt[20995]))/
    (c^2*E^((9*I)*\[Psi])*R)
