(* ::Package:: *)

(* ::Chapter:: *)
(*Supplementary file to  "Gravitational-wave tails of memory " (David Trestini, Luc Blanchet)*)


(* ::Subtitle:: *)
(*Expression of the kernels needed for the expression the the crude radiative quadrupole moment in terms of the radiative canonical moments*)


(* ::Text:: *)
(*This file contains the explicit expressions for the kernels K and L that enter in the definitions of the functionals \[Psi] and \[Chi] as given by Eq. 5.9. These kernels were defined in an integral form in  Eqs. 4.26, 4.30b and 4.31b and their general structure is given in Eq. 5.5. Two examples of the explicit expression of the kernels were given in Eq. 5.6. Only the kernels needed for the radiative quadrupole were computed: these are determined by Eq. 5.11c and Table I. *)


(* ::Section:: *)
(*Notations*)


(* ::Text:: *)
(*We set c = 1.*)
(*Log[x] is the natural logarithm of x*)
(*r0 is the regularization constant associated to the finite part operator.*)
(*PolyLog[n, x] = Subscript[Li, n](x) for n >= 2  (note that for n = 1, Subscript[Li, 1](x) = - Log[1-x])*)
(*\[Rho] and \[Tau] are the two time-variables associated to the kernels*)
(*The kernels are denoted*)
(*	Kbar[k,m,l][\[Rho],\[Tau]] , see eq. 4.26a*)
(*	Lbar[k,m,l][\[Rho],\[Tau]] , see eq. 4.26b*)


(* ::Section::Closed:: *)
(*Kernels associated to the coefficients A of Eq. 5.11c*)


Kbar[1, 0, 2] = (3*\[Rho])/(2*\[Tau]) - Log[\[Rho]/(2*r0)]^2/4 + Log[\[Rho]/(2*r0)]*(-3/2 + Log[\[Tau]/(2*r0)]/2) + (3*Log[\[Tau]/(2*r0)])/2 + (-3/2 - (3*\[Rho]^2)/(2*\[Tau]^2) - (3*\[Rho])/\[Tau])*Log[1 + \[Tau]/\[Rho]] - 
 PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 1, 2] = 1 - \[Rho]/(2*\[Tau]) - Log[\[Rho]/(2*r0)]^2/4 + Log[\[Rho]/(2*r0)]*(-1/2 + Log[\[Tau]/(2*r0)]/2) + (3*Log[\[Tau]/(2*r0)])/2 + (-1/2 + \[Rho]^2/(2*\[Tau]^2))*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 2, 2] = 9/4 - Log[\[Rho]/(2*r0)]^2/4 + (3*Log[\[Tau]/(2*r0)])/2 + (Log[\[Rho]/(2*r0)]*Log[\[Tau]/(2*r0)])/2 - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 3, 2] = 13/4 + \[Tau]/(6*\[Rho]) - \[Tau]^2/(3*\[Rho]^2) - Log[\[Rho]/(2*r0)]^2/4 + Log[\[Rho]/(2*r0)]*(1/3 + \[Tau]^3/(3*\[Rho]^3) + Log[\[Tau]/(2*r0)]/2) + (3/2 - \[Tau]^3/(3*\[Rho]^3))*Log[\[Tau]/(2*r0)] + 
 (1/3 + \[Tau]^3/(3*\[Rho]^3))*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 4, 2] = 4 + (7*\[Tau])/(12*\[Rho]) - (35*\[Tau]^2)/(24*\[Rho]^2) - (7*\[Tau]^3)/(4*\[Rho]^3) - Log[\[Rho]/(2*r0)]^2/4 + 
 Log[\[Rho]/(2*r0)]*(7/12 + (7*\[Tau]^3)/(3*\[Rho]^3) + (7*\[Tau]^4)/(4*\[Rho]^4) + Log[\[Tau]/(2*r0)]/2) + (3/2 - (7*\[Tau]^3)/(3*\[Rho]^3) - (7*\[Tau]^4)/(4*\[Rho]^4))*Log[\[Tau]/(2*r0)] + 
 (7/12 + (7*\[Tau]^3)/(3*\[Rho]^3) + (7*\[Tau]^4)/(4*\[Rho]^4))*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 5, 2]  = 23/5 + (73*\[Tau])/(60*\[Rho]) - (463*\[Tau]^2)/(120*\[Rho]^2) - (243*\[Tau]^3)/(20*\[Rho]^3) - (36*\[Tau]^4)/(5*\[Rho]^4) - Log[\[Rho]/(2*r0)]^2/4 + 
 Log[\[Rho]/(2*r0)]*(47/60 + (28*\[Tau]^3)/(3*\[Rho]^3) + (63*\[Tau]^4)/(4*\[Rho]^4) + (36*\[Tau]^5)/(5*\[Rho]^5) + Log[\[Tau]/(2*r0)]/2) + 
 (3/2 - (28*\[Tau]^3)/(3*\[Rho]^3) - (63*\[Tau]^4)/(4*\[Rho]^4) - (36*\[Tau]^5)/(5*\[Rho]^5))*Log[\[Tau]/(2*r0)] + (47/60 + (28*\[Tau]^3)/(3*\[Rho]^3) + (63*\[Tau]^4)/(4*\[Rho]^4) + (36*\[Tau]^5)/(5*\[Rho]^5))*
  Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 6, 2]  = 51/10 + (41*\[Tau])/(20*\[Rho]) - (163*\[Tau]^2)/(20*\[Rho]^2) - (2899*\[Tau]^3)/(60*\[Rho]^3) - (1309*\[Tau]^4)/(20*\[Rho]^4) - (55*\[Tau]^5)/(2*\[Rho]^5) - Log[\[Rho]/(2*r0)]^2/4 + 
 Log[\[Rho]/(2*r0)]*(19/20 + (28*\[Tau]^3)/\[Rho]^3 + (315*\[Tau]^4)/(4*\[Rho]^4) + (396*\[Tau]^5)/(5*\[Rho]^5) + (55*\[Tau]^6)/(2*\[Rho]^6) + Log[\[Tau]/(2*r0)]/2) + 
 (3/2 - (28*\[Tau]^3)/\[Rho]^3 - (315*\[Tau]^4)/(4*\[Rho]^4) - (396*\[Tau]^5)/(5*\[Rho]^5) - (55*\[Tau]^6)/(2*\[Rho]^6))*Log[\[Tau]/(2*r0)] + 
 (19/20 + (28*\[Tau]^3)/\[Rho]^3 + (315*\[Tau]^4)/(4*\[Rho]^4) + (396*\[Tau]^5)/(5*\[Rho]^5) + (55*\[Tau]^6)/(2*\[Rho]^6))*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 0, 3]  = (10*\[Rho]^2)/(3*\[Tau]^2) + (35*\[Rho])/(6*\[Tau]) - Log[\[Rho]/(2*r0)]^2/4 + Log[\[Rho]/(2*r0)]*(-11/6 + Log[\[Tau]/(2*r0)]/2) + (11*Log[\[Tau]/(2*r0)])/6 + 
 (-11/6 - (10*\[Rho]^3)/(3*\[Tau]^3) - (15*\[Rho]^2)/(2*\[Tau]^2) - (6*\[Rho])/\[Tau])*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 1, 3] = 1 - (5*\[Rho]^2)/(3*\[Tau]^2) - (5*\[Rho])/(3*\[Tau]) - Log[\[Rho]/(2*r0)]^2/4 + Log[\[Rho]/(2*r0)]*(-5/6 + Log[\[Tau]/(2*r0)]/2) + (11*Log[\[Tau]/(2*r0)])/6 + 
 (-5/6 + (5*\[Rho]^3)/(3*\[Tau]^3) + (5*\[Rho]^2)/(2*\[Tau]^2))*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 2, 3] = 9/4 + \[Rho]^2/(3*\[Tau]^2) - \[Rho]/(6*\[Tau]) - Log[\[Rho]/(2*r0)]^2/4 + Log[\[Rho]/(2*r0)]*(-1/3 + Log[\[Tau]/(2*r0)]/2) + (11*Log[\[Tau]/(2*r0)])/6 + (-1/3 - \[Rho]^3/(3*\[Tau]^3))*Log[1 + \[Tau]/\[Rho]] - 
 PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 3, 3] = 121/36 - Log[\[Rho]/(2*r0)]^2/4 + (11*Log[\[Tau]/(2*r0)])/6 + (Log[\[Rho]/(2*r0)]*Log[\[Tau]/(2*r0)])/2 - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 4, 3] = 77/18 + \[Tau]/(12*\[Rho]) - \[Tau]^2/(8*\[Rho]^2) + \[Tau]^3/(4*\[Rho]^3) - Log[\[Rho]/(2*r0)]^2/4 + Log[\[Rho]/(2*r0)]*(1/4 - \[Tau]^4/(4*\[Rho]^4) + Log[\[Tau]/(2*r0)]/2) + 
 (11/6 + \[Tau]^4/(4*\[Rho]^4))*Log[\[Tau]/(2*r0)] + (1/4 - \[Tau]^4/(4*\[Rho]^4))*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 5, 3] = 451/90 + (3*\[Tau])/(10*\[Rho]) - (21*\[Tau]^2)/(40*\[Rho]^2) + (27*\[Tau]^3)/(20*\[Rho]^3) + (9*\[Tau]^4)/(5*\[Rho]^4) - Log[\[Rho]/(2*r0)]^2/4 + 
 Log[\[Rho]/(2*r0)]*(9/20 - (9*\[Tau]^4)/(4*\[Rho]^4) - (9*\[Tau]^5)/(5*\[Rho]^5) + Log[\[Tau]/(2*r0)]/2) + (11/6 + (9*\[Tau]^4)/(4*\[Rho]^4) + (9*\[Tau]^5)/(5*\[Rho]^5))*Log[\[Tau]/(2*r0)] + 
 (9/20 - (9*\[Tau]^4)/(4*\[Rho]^4) - (9*\[Tau]^5)/(5*\[Rho]^5))*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 6, 3] = 253/45 + (19*\[Tau])/(30*\[Rho]) - (79*\[Tau]^2)/(60*\[Rho]^2) + (793*\[Tau]^3)/(180*\[Rho]^3) + (913*\[Tau]^4)/(60*\[Rho]^4) + (55*\[Tau]^5)/(6*\[Rho]^5) - Log[\[Rho]/(2*r0)]^2/4 + 
 Log[\[Rho]/(2*r0)]*(37/60 - (45*\[Tau]^4)/(4*\[Rho]^4) - (99*\[Tau]^5)/(5*\[Rho]^5) - (55*\[Tau]^6)/(6*\[Rho]^6) + Log[\[Tau]/(2*r0)]/2) + 
 (11/6 + (45*\[Tau]^4)/(4*\[Rho]^4) + (99*\[Tau]^5)/(5*\[Rho]^5) + (55*\[Tau]^6)/(6*\[Rho]^6))*Log[\[Tau]/(2*r0)] + (37/60 - (45*\[Tau]^4)/(4*\[Rho]^4) - (99*\[Tau]^5)/(5*\[Rho]^5) - (55*\[Tau]^6)/(6*\[Rho]^6))*
  Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 0, 4] = (35*\[Rho]^3)/(4*\[Tau]^3) + (455*\[Rho]^2)/(24*\[Tau]^2) + (55*\[Rho])/(4*\[Tau]) - Log[\[Rho]/(2*r0)]^2/4 + Log[\[Rho]/(2*r0)]*(-25/12 + Log[\[Tau]/(2*r0)]/2) + (25*Log[\[Tau]/(2*r0)])/12 + 
 (-25/12 - (35*\[Rho]^4)/(4*\[Tau]^4) - (70*\[Rho]^3)/(3*\[Tau]^3) - (45*\[Rho]^2)/(2*\[Tau]^2) - (10*\[Rho])/\[Tau])*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 1, 4]= 1 - (21*\[Rho]^3)/(4*\[Tau]^3) - (217*\[Rho]^2)/(24*\[Tau]^2) - (41*\[Rho])/(12*\[Tau]) - Log[\[Rho]/(2*r0)]^2/4 + Log[\[Rho]/(2*r0)]*(-13/12 + Log[\[Tau]/(2*r0)]/2) + (25*Log[\[Tau]/(2*r0)])/12 + 
 (-13/12 + (21*\[Rho]^4)/(4*\[Tau]^4) + (35*\[Rho]^3)/(3*\[Tau]^3) + (15*\[Rho]^2)/(2*\[Tau]^2))*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 2, 4] = 9/4 + (7*\[Rho]^3)/(4*\[Tau]^3) + (35*\[Rho]^2)/(24*\[Tau]^2) - (7*\[Rho])/(12*\[Tau]) - Log[\[Rho]/(2*r0)]^2/4 + Log[\[Rho]/(2*r0)]*(-7/12 + Log[\[Tau]/(2*r0)]/2) + (25*Log[\[Tau]/(2*r0)])/12 + 
 (-7/12 - (7*\[Rho]^4)/(4*\[Tau]^4) - (7*\[Rho]^3)/(3*\[Tau]^3))*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 3, 4] = 121/36 - \[Rho]^3/(4*\[Tau]^3) + \[Rho]^2/(8*\[Tau]^2) - \[Rho]/(12*\[Tau]) - Log[\[Rho]/(2*r0)]^2/4 + Log[\[Rho]/(2*r0)]*(-1/4 + Log[\[Tau]/(2*r0)]/2) + (25*Log[\[Tau]/(2*r0)])/12 + 
 (-1/4 + \[Rho]^4/(4*\[Tau]^4))*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 4, 4] = 625/144 - Log[\[Rho]/(2*r0)]^2/4 + (25*Log[\[Tau]/(2*r0)])/12 + (Log[\[Rho]/(2*r0)]*Log[\[Tau]/(2*r0)])/2 - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 5, 4] = 745/144 + \[Tau]/(20*\[Rho]) - \[Tau]^2/(15*\[Rho]^2) + \[Tau]^3/(10*\[Rho]^3) - \[Tau]^4/(5*\[Rho]^4) - Log[\[Rho]/(2*r0)]^2/4 + Log[\[Rho]/(2*r0)]*(1/5 + \[Tau]^5/(5*\[Rho]^5) + Log[\[Tau]/(2*r0)]/2) + 
 (25/12 - \[Tau]^5/(5*\[Rho]^5))*Log[\[Tau]/(2*r0)] + (1/5 + \[Tau]^5/(5*\[Rho]^5))*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


Kbar[1, 6, 4] = 845/144 + (11*\[Tau])/(60*\[Rho]) - (11*\[Tau]^2)/(40*\[Rho]^2) + (22*\[Tau]^3)/(45*\[Rho]^3) - (77*\[Tau]^4)/(60*\[Rho]^4) - (11*\[Tau]^5)/(6*\[Rho]^5) - Log[\[Rho]/(2*r0)]^2/4 + 
 Log[\[Rho]/(2*r0)]*(11/30 + (11*\[Tau]^5)/(5*\[Rho]^5) + (11*\[Tau]^6)/(6*\[Rho]^6) + Log[\[Tau]/(2*r0)]/2) + (25/12 - (11*\[Tau]^5)/(5*\[Rho]^5) - (11*\[Tau]^6)/(6*\[Rho]^6))*Log[\[Tau]/(2*r0)] + 
 (11/30 + (11*\[Tau]^5)/(5*\[Rho]^5) + (11*\[Tau]^6)/(6*\[Rho]^6))*Log[1 + \[Tau]/\[Rho]] - PolyLog[2, -(\[Tau]/\[Rho])]/2


(* ::Section::Closed:: *)
(*Kernels associated to the coefficients B of Eq. 5.11c*)


Kbar[2, 0, 2] = (-2 - (2*\[Rho])/\[Tau] + (1 + (2*\[Rho]^2)/\[Tau]^2 + (3*\[Rho])/\[Tau])*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 1, 2] = (1/2 + \[Rho]/\[Tau] + (-(\[Rho]^2/\[Tau]^2) - \[Rho]/\[Tau])*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 2, 2] = (1/10 - \[Rho]/(5*\[Tau]) + \[Tau]/(10*\[Rho]) - \[Tau]^2/(5*\[Rho]^2) + (\[Tau]^3*Log[\[Rho]/(2*r0)])/(5*\[Rho]^3) - (\[Tau]^3*Log[\[Tau]/(2*r0)])/(5*\[Rho]^3) + (\[Rho]^2/(5*\[Tau]^2) + \[Tau]^3/(5*\[Rho]^3))*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 3, 2] = (\[Tau]/(6*\[Rho]) - \[Tau]^2/(2*\[Rho]^2) - \[Tau]^3/\[Rho]^3 + (\[Tau]^3/\[Rho]^3 + \[Tau]^4/\[Rho]^4)*Log[\[Rho]/(2*r0)] + (-(\[Tau]^3/\[Rho]^3) - \[Tau]^4/\[Rho]^4)*Log[\[Tau]/(2*r0)] + 
  (\[Tau]^3/\[Rho]^3 + \[Tau]^4/\[Rho]^4)*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 4, 2] = (\[Tau]/(6*\[Rho]) - (5*\[Tau]^2)/(6*\[Rho]^2) - (5*\[Tau]^3)/\[Rho]^3 - (4*\[Tau]^4)/\[Rho]^4 + ((3*\[Tau]^3)/\[Rho]^3 + (7*\[Tau]^4)/\[Rho]^4 + (4*\[Tau]^5)/\[Rho]^5)*Log[\[Rho]/(2*r0)] + 
  ((-3*\[Tau]^3)/\[Rho]^3 - (7*\[Tau]^4)/\[Rho]^4 - (4*\[Tau]^5)/\[Rho]^5)*Log[\[Tau]/(2*r0)] + ((3*\[Tau]^3)/\[Rho]^3 + (7*\[Tau]^4)/\[Rho]^4 + (4*\[Tau]^5)/\[Rho]^5)*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 5, 2]  = (\[Tau]/(6*\[Rho]) - (5*\[Tau]^2)/(4*\[Rho]^2) - (15*\[Tau]^3)/\[Rho]^3 - (57*\[Tau]^4)/(2*\[Rho]^4) - (15*\[Tau]^5)/\[Rho]^5 + ((7*\[Tau]^3)/\[Rho]^3 + (28*\[Tau]^4)/\[Rho]^4 + (36*\[Tau]^5)/\[Rho]^5 + (15*\[Tau]^6)/\[Rho]^6)*
   Log[\[Rho]/(2*r0)] + ((-7*\[Tau]^3)/\[Rho]^3 - (28*\[Tau]^4)/\[Rho]^4 - (36*\[Tau]^5)/\[Rho]^5 - (15*\[Tau]^6)/\[Rho]^6)*Log[\[Tau]/(2*r0)] + 
  ((7*\[Tau]^3)/\[Rho]^3 + (28*\[Tau]^4)/\[Rho]^4 + (36*\[Tau]^5)/\[Rho]^5 + (15*\[Tau]^6)/\[Rho]^6)*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 0, 3]  = (-8/3 - (5*\[Rho]^2)/\[Tau]^2 - (15*\[Rho])/(2*\[Tau]) + (1 + (5*\[Rho]^3)/\[Tau]^3 + (10*\[Rho]^2)/\[Tau]^2 + (6*\[Rho])/\[Tau])*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 1, 3] = (1/2 + (3*\[Rho]^2)/\[Tau]^2 + (7*\[Rho])/(2*\[Tau]) + ((-3*\[Rho]^3)/\[Tau]^3 - (5*\[Rho]^2)/\[Tau]^2 - (2*\[Rho])/\[Tau])*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 2, 3] = (1/6 - \[Rho]^2/\[Tau]^2 - \[Rho]/(2*\[Tau]) + (\[Rho]^3/\[Tau]^3 + \[Rho]^2/\[Tau]^2)*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 3, 3] = (1/21 + \[Rho]^2/(7*\[Tau]^2) - \[Rho]/(14*\[Tau]) + \[Tau]/(21*\[Rho]) - \[Tau]^2/(14*\[Rho]^2) + \[Tau]^3/(7*\[Rho]^3) - (\[Tau]^4*Log[\[Rho]/(2*r0)])/(7*\[Rho]^4) + (\[Tau]^4*Log[\[Tau]/(2*r0)])/(7*\[Rho]^4) + 
  (-1/7*\[Rho]^3/\[Tau]^3 - \[Tau]^4/(7*\[Rho]^4))*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 4, 3] = (\[Tau]/(12*\[Rho]) - \[Tau]^2/(6*\[Rho]^2) + \[Tau]^3/(2*\[Rho]^3) + \[Tau]^4/\[Rho]^4 + (-(\[Tau]^4/\[Rho]^4) - \[Tau]^5/\[Rho]^5)*Log[\[Rho]/(2*r0)] + (\[Tau]^4/\[Rho]^4 + \[Tau]^5/\[Rho]^5)*Log[\[Tau]/(2*r0)] + 
  (-(\[Tau]^4/\[Rho]^4) - \[Tau]^5/\[Rho]^5)*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 5, 3] = (\[Tau]/(12*\[Rho]) - \[Tau]^2/(4*\[Rho]^2) + (7*\[Tau]^3)/(6*\[Rho]^3) + (13*\[Tau]^4)/(2*\[Rho]^4) + (5*\[Tau]^5)/\[Rho]^5 + ((-4*\[Tau]^4)/\[Rho]^4 - (9*\[Tau]^5)/\[Rho]^5 - (5*\[Tau]^6)/\[Rho]^6)*Log[\[Rho]/(2*r0)] + 
  ((4*\[Tau]^4)/\[Rho]^4 + (9*\[Tau]^5)/\[Rho]^5 + (5*\[Tau]^6)/\[Rho]^6)*Log[\[Tau]/(2*r0)] + ((-4*\[Tau]^4)/\[Rho]^4 - (9*\[Tau]^5)/\[Rho]^5 - (5*\[Tau]^6)/\[Rho]^6)*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 0, 4]= (-19/6 - (14*\[Rho]^3)/\[Tau]^3 - (28*\[Rho]^2)/\[Tau]^2 - (103*\[Rho])/(6*\[Tau]) + (1 + (14*\[Rho]^4)/\[Tau]^4 + (35*\[Rho]^3)/\[Tau]^3 + (30*\[Rho]^2)/\[Tau]^2 + (10*\[Rho])/\[Tau])*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 1, 4] = (1/2 + (28*\[Rho]^3)/(3*\[Tau]^3) + (49*\[Rho]^2)/(3*\[Tau]^2) + (137*\[Rho])/(18*\[Tau]) + ((-28*\[Rho]^4)/(3*\[Tau]^4) - (21*\[Rho]^3)/\[Tau]^3 - (15*\[Rho]^2)/\[Tau]^2 - (10*\[Rho])/(3*\[Tau]))*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 2, 4] = (1/6 - (4*\[Rho]^3)/\[Tau]^3 - (5*\[Rho]^2)/\[Tau]^2 - (5*\[Rho])/(6*\[Tau]) + ((4*\[Rho]^4)/\[Tau]^4 + (7*\[Rho]^3)/\[Tau]^3 + (3*\[Rho]^2)/\[Tau]^2)*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 3, 4] = (1/12 + \[Rho]^3/\[Tau]^3 + \[Rho]^2/(2*\[Tau]^2) - \[Rho]/(6*\[Tau]) + (-(\[Rho]^4/\[Tau]^4) - \[Rho]^3/\[Tau]^3)*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 4, 4] = (1/36 - \[Rho]^3/(9*\[Tau]^3) + \[Rho]^2/(18*\[Tau]^2) - \[Rho]/(27*\[Tau]) + \[Tau]/(36*\[Rho]) - \[Tau]^2/(27*\[Rho]^2) + \[Tau]^3/(18*\[Rho]^3) - \[Tau]^4/(9*\[Rho]^4) + (\[Tau]^5*Log[\[Rho]/(2*r0)])/(9*\[Rho]^5) - 
  (\[Tau]^5*Log[\[Tau]/(2*r0)])/(9*\[Rho]^5) + (\[Rho]^4/(9*\[Tau]^4) + \[Tau]^5/(9*\[Rho]^5))*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Kbar[2, 5, 4] = (\[Tau]/(20*\[Rho]) - \[Tau]^2/(12*\[Rho]^2) + \[Tau]^3/(6*\[Rho]^3) - \[Tau]^4/(2*\[Rho]^4) - \[Tau]^5/\[Rho]^5 + (\[Tau]^5/\[Rho]^5 + \[Tau]^6/\[Rho]^6)*Log[\[Rho]/(2*r0)] + (-(\[Tau]^5/\[Rho]^5) - \[Tau]^6/\[Rho]^6)*Log[\[Tau]/(2*r0)] + 
  (\[Tau]^5/\[Rho]^5 + \[Tau]^6/\[Rho]^6)*Log[1 + \[Tau]/\[Rho]])/\[Tau]


(* ::Section::Closed:: *)
(*Kernels associated to the coefficients C of Eq. 5.11c*)


Lbar[1, 0, 2] = (9*\[Rho])/(4*\[Tau]) - Log[\[Rho]/(2*r0)]^3/6 + Log[\[Rho]/(2*r0)]^2*(-3/2 + Log[\[Tau]/(2*r0)]/4) + (9*Log[\[Tau]/(2*r0)])/4 + 
 Log[\[Rho]/(2*r0)]*(-9/4 + (3*\[Rho])/(2*\[Tau]) + (3*Log[\[Tau]/(2*r0)])/2) + (-9/4 - (3*\[Rho]^2)/(4*\[Tau]^2) - (3*\[Rho])/\[Tau] + (-3/2 - (3*\[Rho]^2)/(2*\[Tau]^2) - (3*\[Rho])/\[Tau])*Log[\[Rho]/(2*r0)])*
  Log[1 + \[Tau]/\[Rho]] + ((3*\[Rho]^2)/(2*\[Tau]^2) + (3*\[Rho])/\[Tau] - Log[\[Rho]/(2*r0)]/2)*PolyLog[2, -(\[Tau]/\[Rho])] - PolyLog[3, -(\[Tau]/\[Rho])]/2


Lbar[1, 1, 2] = 5/6 - (5*\[Rho])/(12*\[Tau]) - Log[\[Rho]/(2*r0)]^3/6 + Log[\[Rho]/(2*r0)]^2*(-1 + Log[\[Tau]/(2*r0)]/4) + (9/4 - \[Tau]/(6*\[Rho]))*Log[\[Tau]/(2*r0)] + 
 Log[\[Rho]/(2*r0)]*(-5/4 - \[Rho]/(2*\[Tau]) + \[Tau]/(6*\[Rho]) + (3*Log[\[Tau]/(2*r0)])/2) + (-5/4 - \[Rho]^2/(12*\[Tau]^2) - (3*\[Rho])/(2*\[Tau]) + \[Tau]/(6*\[Rho]) + (-1/2 + \[Rho]^2/(2*\[Tau]^2))*Log[\[Rho]/(2*r0)])*
  Log[1 + \[Tau]/\[Rho]] + (-1 - \[Rho]^2/(2*\[Tau]^2) - Log[\[Rho]/(2*r0)]/2)*PolyLog[2, -(\[Tau]/\[Rho])] - PolyLog[3, -(\[Tau]/\[Rho])]/2


Lbar[1, 2, 2] = 33/16 - \[Rho]/(8*\[Tau]) + \[Tau]/(8*\[Rho]) - Log[\[Rho]/(2*r0)]^3/6 + Log[\[Rho]/(2*r0)]^2*(-3/4 + Log[\[Tau]/(2*r0)]/4) + (9/4 - \[Tau]/(2*\[Rho]) + \[Tau]^2/(8*\[Rho]^2))*Log[\[Tau]/(2*r0)] + 
 Log[\[Rho]/(2*r0)]*(\[Tau]/(2*\[Rho]) - \[Tau]^2/(8*\[Rho]^2) + (3*Log[\[Tau]/(2*r0)])/2) + (\[Rho]^2/(8*\[Tau]^2) - \[Rho]/(2*\[Tau]) + \[Tau]/(2*\[Rho]) - \[Tau]^2/(8*\[Rho]^2))*Log[1 + \[Tau]/\[Rho]] + 
 (-3/2 - Log[\[Rho]/(2*r0)]/2)*PolyLog[2, -(\[Tau]/\[Rho])] - PolyLog[3, -(\[Tau]/\[Rho])]/2


Lbar[1, 3, 2]= 261/80 - \[Rho]/(40*\[Tau]) + (91*\[Tau])/(180*\[Rho]) + (43*\[Tau]^2)/(180*\[Rho]^2) - Log[\[Rho]/(2*r0)]^3/6 + Log[\[Rho]/(2*r0)]^2*(-7/12 + \[Tau]^3/(3*\[Rho]^3) + Log[\[Tau]/(2*r0)]/4) + 
 (9/4 - \[Tau]/\[Rho] + (5*\[Tau]^2)/(8*\[Rho]^2) + (43*\[Tau]^3)/(180*\[Rho]^3))*Log[\[Tau]/(2*r0)] + 
 Log[\[Rho]/(2*r0)]*(1 + (7*\[Tau])/(6*\[Rho]) - (23*\[Tau]^2)/(24*\[Rho]^2) - (43*\[Tau]^3)/(180*\[Rho]^3) + (3/2 - \[Tau]^3/(3*\[Rho]^3))*Log[\[Tau]/(2*r0)]) + 
 (10/9 + \[Rho]^2/(40*\[Tau]^2) - \[Rho]/(4*\[Tau]) + \[Tau]/\[Rho] - (5*\[Tau]^2)/(8*\[Rho]^2) - (43*\[Tau]^3)/(180*\[Rho]^3) + (1/3 + \[Tau]^3/(3*\[Rho]^3))*Log[\[Rho]/(2*r0)])*Log[1 + \[Tau]/\[Rho]] + 
 (-11/6 - Log[\[Rho]/(2*r0)]/2)*PolyLog[2, -(\[Tau]/\[Rho])] - PolyLog[3, -(\[Tau]/\[Rho])]/2


Lbar[1, 4, 2] = 1027/240 - \[Rho]/(120*\[Tau]) + (839*\[Tau])/(720*\[Rho]) + (427*\[Tau]^2)/(288*\[Rho]^2) + (91*\[Tau]^3)/(240*\[Rho]^3) - Log[\[Rho]/(2*r0)]^3/6 + 
 Log[\[Rho]/(2*r0)]^2*(-11/24 + (7*\[Tau]^3)/(3*\[Rho]^3) + (7*\[Tau]^4)/(4*\[Rho]^4) + Log[\[Tau]/(2*r0)]/4) + 
 (9/4 - (5*\[Tau])/(3*\[Rho]) + (15*\[Tau]^2)/(8*\[Rho]^2) + (301*\[Tau]^3)/(180*\[Rho]^3) + (91*\[Tau]^4)/(240*\[Rho]^4))*Log[\[Tau]/(2*r0)] + 
 Log[\[Rho]/(2*r0)]*(7/4 + (9*\[Tau])/(4*\[Rho]) - (10*\[Tau]^2)/(3*\[Rho]^2) - (154*\[Tau]^3)/(45*\[Rho]^3) - (91*\[Tau]^4)/(240*\[Rho]^4) + 
   (3/2 - (7*\[Tau]^3)/(3*\[Rho]^3) - (7*\[Tau]^4)/(4*\[Rho]^4))*Log[\[Tau]/(2*r0)]) + (301/144 + \[Rho]^2/(120*\[Tau]^2) - (3*\[Rho])/(20*\[Tau]) + (5*\[Tau])/(3*\[Rho]) - (15*\[Tau]^2)/(8*\[Rho]^2) - 
   (301*\[Tau]^3)/(180*\[Rho]^3) - (91*\[Tau]^4)/(240*\[Rho]^4) + (7/12 + (7*\[Tau]^3)/(3*\[Rho]^3) + (7*\[Tau]^4)/(4*\[Rho]^4))*Log[\[Rho]/(2*r0)])*Log[1 + \[Tau]/\[Rho]] + 
 (-25/12 - Log[\[Rho]/(2*r0)]/2)*PolyLog[2, -(\[Tau]/\[Rho])] - PolyLog[3, -(\[Tau]/\[Rho])]/2


Lbar[1, 5, 2]  = 2871/560 - \[Rho]/(280*\[Tau]) + (51287*\[Tau])/(25200*\[Rho]) + (260053*\[Tau]^2)/(50400*\[Rho]^2) + (8811*\[Tau]^3)/(2800*\[Rho]^3) + (93*\[Tau]^4)/(175*\[Rho]^4) - Log[\[Rho]/(2*r0)]^3/6 + 
 Log[\[Rho]/(2*r0)]^2*(-43/120 + (28*\[Tau]^3)/(3*\[Rho]^3) + (63*\[Tau]^4)/(4*\[Rho]^4) + (36*\[Tau]^5)/(5*\[Rho]^5) + Log[\[Tau]/(2*r0)]/4) + 
 (9/4 - (5*\[Tau])/(2*\[Rho]) + (35*\[Tau]^2)/(8*\[Rho]^2) + (301*\[Tau]^3)/(45*\[Rho]^3) + (273*\[Tau]^4)/(80*\[Rho]^4) + (93*\[Tau]^5)/(175*\[Rho]^5))*Log[\[Tau]/(2*r0)] + 
 Log[\[Rho]/(2*r0)]*(47/20 + (223*\[Tau])/(60*\[Rho]) - (247*\[Tau]^2)/(30*\[Rho]^2) - (3391*\[Tau]^3)/(180*\[Rho]^3) - (849*\[Tau]^4)/(80*\[Rho]^4) - (93*\[Tau]^5)/(175*\[Rho]^5) + 
   (3/2 - (28*\[Tau]^3)/(3*\[Rho]^3) - (63*\[Tau]^4)/(4*\[Rho]^4) - (36*\[Tau]^5)/(5*\[Rho]^5))*Log[\[Tau]/(2*r0)]) + 
 (10669/3600 + \[Rho]^2/(280*\[Tau]^2) - \[Rho]/(10*\[Tau]) + (5*\[Tau])/(2*\[Rho]) - (35*\[Tau]^2)/(8*\[Rho]^2) - (301*\[Tau]^3)/(45*\[Rho]^3) - (273*\[Tau]^4)/(80*\[Rho]^4) - (93*\[Tau]^5)/(175*\[Rho]^5) + 
   (47/60 + (28*\[Tau]^3)/(3*\[Rho]^3) + (63*\[Tau]^4)/(4*\[Rho]^4) + (36*\[Tau]^5)/(5*\[Rho]^5))*Log[\[Rho]/(2*r0)])*Log[1 + \[Tau]/\[Rho]] + (-137/60 - Log[\[Rho]/(2*r0)]/2)*PolyLog[2, -(\[Tau]/\[Rho])] - 
 PolyLog[3, -(\[Tau]/\[Rho])]/2


Lbar[1, 6, 2]  = 6549/1120 - \[Rho]/(560*\[Tau]) + (4239*\[Tau])/(1400*\[Rho]) + (149829*\[Tau]^2)/(11200*\[Rho]^2) + (179947*\[Tau]^3)/(12600*\[Rho]^3) + (13519*\[Tau]^4)/(2400*\[Rho]^4) + (143*\[Tau]^5)/(336*\[Rho]^5) - 
 Log[\[Rho]/(2*r0)]^3/6 + Log[\[Rho]/(2*r0)]^2*(-11/40 + (28*\[Tau]^3)/\[Rho]^3 + (315*\[Tau]^4)/(4*\[Rho]^4) + (396*\[Tau]^5)/(5*\[Rho]^5) + (55*\[Tau]^6)/(2*\[Rho]^6) + Log[\[Tau]/(2*r0)]/4) + 
 (9/4 - (7*\[Tau])/(2*\[Rho]) + (35*\[Tau]^2)/(4*\[Rho]^2) + (301*\[Tau]^3)/(15*\[Rho]^3) + (273*\[Tau]^4)/(16*\[Rho]^4) + (1023*\[Tau]^5)/(175*\[Rho]^5) + (143*\[Tau]^6)/(336*\[Rho]^6))*Log[\[Tau]/(2*r0)] + 
 Log[\[Rho]/(2*r0)]*(57/20 + (111*\[Tau])/(20*\[Rho]) - (169*\[Tau]^2)/(10*\[Rho]^2) - (4103*\[Tau]^3)/(60*\[Rho]^3) - (6601*\[Tau]^4)/(80*\[Rho]^4) - (11671*\[Tau]^5)/(350*\[Rho]^5) - (143*\[Tau]^6)/(336*\[Rho]^6) + 
   (3/2 - (28*\[Tau]^3)/\[Rho]^3 - (315*\[Tau]^4)/(4*\[Rho]^4) - (396*\[Tau]^5)/(5*\[Rho]^5) - (55*\[Tau]^6)/(2*\[Rho]^6))*Log[\[Tau]/(2*r0)]) + 
 (1501/400 + \[Rho]^2/(560*\[Tau]^2) - \[Rho]/(14*\[Tau]) + (7*\[Tau])/(2*\[Rho]) - (35*\[Tau]^2)/(4*\[Rho]^2) - (301*\[Tau]^3)/(15*\[Rho]^3) - (273*\[Tau]^4)/(16*\[Rho]^4) - (1023*\[Tau]^5)/(175*\[Rho]^5) - 
   (143*\[Tau]^6)/(336*\[Rho]^6) + (19/20 + (28*\[Tau]^3)/\[Rho]^3 + (315*\[Tau]^4)/(4*\[Rho]^4) + (396*\[Tau]^5)/(5*\[Rho]^5) + (55*\[Tau]^6)/(2*\[Rho]^6))*Log[\[Rho]/(2*r0)])*Log[1 + \[Tau]/\[Rho]] + 
 (-49/20 - Log[\[Rho]/(2*r0)]/2)*PolyLog[2, -(\[Tau]/\[Rho])] - PolyLog[3, -(\[Tau]/\[Rho])]/2


(* ::Section:: *)
(*Kernels associated to the coefficients D of Eq. 5.11c*)


Lbar[2, 0, 2] = (-11/3 - (8*\[Rho])/(3*\[Tau]) + (-2 - (2*\[Rho])/\[Tau] + \[Tau]/(6*\[Rho]))*Log[\[Rho]/(2*r0)] - (\[Tau]*Log[\[Tau]/(2*r0)])/(6*\[Rho]) + 
  (1 + (2*\[Rho]^2)/(3*\[Tau]^2) + (3*\[Rho])/(2*\[Tau]) + \[Tau]/(6*\[Rho]) + (1 + (2*\[Rho]^2)/\[Tau]^2 + (3*\[Rho])/\[Tau])*Log[\[Rho]/(2*r0)])*Log[1 + \[Tau]/\[Rho]] + 
  (-1 - (2*\[Rho]^2)/\[Tau]^2 - (3*\[Rho])/\[Tau])*PolyLog[2, -(\[Tau]/\[Rho])])/\[Tau]


Lbar[2, 1, 2] = (13/24 + (13*\[Rho])/(12*\[Tau]) + \[Tau]/(12*\[Rho]) + (1/2 + \[Rho]/\[Tau] + \[Tau]/(6*\[Rho]) - \[Tau]^2/(12*\[Rho]^2))*Log[\[Rho]/(2*r0)] + (-1/6*\[Tau]/\[Rho] + \[Tau]^2/(12*\[Rho]^2))*Log[\[Tau]/(2*r0)] + 
  (1/2 - \[Rho]^2/(12*\[Tau]^2) + \[Rho]/(6*\[Tau]) + \[Tau]/(6*\[Rho]) - \[Tau]^2/(12*\[Rho]^2) + (-(\[Rho]^2/\[Tau]^2) - \[Rho]/\[Tau])*Log[\[Rho]/(2*r0)])*Log[1 + \[Tau]/\[Rho]] + (\[Rho]^2/\[Tau]^2 + \[Rho]/\[Tau])*PolyLog[2, -(\[Tau]/\[Rho])])/
 \[Tau]


Lbar[2, 2, 2] = (137/600 - (17*\[Rho])/(300*\[Tau]) + (107*\[Tau])/(600*\[Rho]) + (43*\[Tau]^2)/(300*\[Rho]^2) + (\[Tau]^3*Log[\[Rho]/(2*r0)]^2)/(5*\[Rho]^3) + 
  (-1/6*\[Tau]/\[Rho] + \[Tau]^2/(4*\[Rho]^2) + (43*\[Tau]^3)/(300*\[Rho]^3))*Log[\[Tau]/(2*r0)] + Log[\[Rho]/(2*r0)]*(1/10 - \[Rho]/(5*\[Tau]) + (4*\[Tau])/(15*\[Rho]) - (9*\[Tau]^2)/(20*\[Rho]^2) - 
    (43*\[Tau]^3)/(300*\[Rho]^3) - (\[Tau]^3*Log[\[Tau]/(2*r0)])/(5*\[Rho]^3)) + (1/6 - (43*\[Rho]^2)/(300*\[Tau]^2) - \[Rho]/(4*\[Tau]) + \[Tau]/(6*\[Rho]) - \[Tau]^2/(4*\[Rho]^2) - (43*\[Tau]^3)/(300*\[Rho]^3) + 
    (\[Rho]^2/(5*\[Tau]^2) + \[Tau]^3/(5*\[Rho]^3))*Log[\[Rho]/(2*r0)])*Log[1 + \[Tau]/\[Rho]] - (\[Rho]^2*PolyLog[2, -(\[Tau]/\[Rho])])/(5*\[Tau]^2))/\[Tau]


Lbar[2, 3, 2] = (1/15 - \[Rho]/(30*\[Tau]) + (77*\[Tau])/(360*\[Rho]) + (73*\[Tau]^2)/(120*\[Rho]^2) + (13*\[Tau]^3)/(60*\[Rho]^3) + (\[Tau]^3/\[Rho]^3 + \[Tau]^4/\[Rho]^4)*Log[\[Rho]/(2*r0)]^2 + 
  (-1/6*\[Tau]/\[Rho] + \[Tau]^2/(2*\[Rho]^2) + (43*\[Tau]^3)/(60*\[Rho]^3) + (13*\[Tau]^4)/(60*\[Rho]^4))*Log[\[Tau]/(2*r0)] + 
  Log[\[Rho]/(2*r0)]*(\[Tau]/(3*\[Rho]) - \[Tau]^2/\[Rho]^2 - (103*\[Tau]^3)/(60*\[Rho]^3) - (13*\[Tau]^4)/(60*\[Rho]^4) + (-(\[Tau]^3/\[Rho]^3) - \[Tau]^4/\[Rho]^4)*Log[\[Tau]/(2*r0)]) + 
  (1/12 + \[Rho]^2/(30*\[Tau]^2) - \[Rho]/(20*\[Tau]) + \[Tau]/(6*\[Rho]) - \[Tau]^2/(2*\[Rho]^2) - (43*\[Tau]^3)/(60*\[Rho]^3) - (13*\[Tau]^4)/(60*\[Rho]^4) + (\[Tau]^3/\[Rho]^3 + \[Tau]^4/\[Rho]^4)*Log[\[Rho]/(2*r0)])*
   Log[1 + \[Tau]/\[Rho]])/\[Tau]


Lbar[2, 4, 2] = (2/105 - \[Rho]/(210*\[Tau]) + (479*\[Tau])/(2520*\[Rho]) + (751*\[Tau]^2)/(504*\[Rho]^2) + (115*\[Tau]^3)/(84*\[Rho]^3) + (31*\[Tau]^4)/(105*\[Rho]^4) + 
  ((3*\[Tau]^3)/\[Rho]^3 + (7*\[Tau]^4)/\[Rho]^4 + (4*\[Tau]^5)/\[Rho]^5)*Log[\[Rho]/(2*r0)]^2 + (-1/6*\[Tau]/\[Rho] + (5*\[Tau]^2)/(6*\[Rho]^2) + (43*\[Tau]^3)/(20*\[Rho]^3) + (91*\[Tau]^4)/(60*\[Rho]^4) + 
    (31*\[Tau]^5)/(105*\[Rho]^5))*Log[\[Tau]/(2*r0)] + Log[\[Rho]/(2*r0)]*(\[Tau]/(3*\[Rho]) - (5*\[Tau]^2)/(3*\[Rho]^2) - (143*\[Tau]^3)/(20*\[Rho]^3) - (331*\[Tau]^4)/(60*\[Rho]^4) - (31*\[Tau]^5)/(105*\[Rho]^5) + 
    ((-3*\[Tau]^3)/\[Rho]^3 - (7*\[Tau]^4)/\[Rho]^4 - (4*\[Tau]^5)/\[Rho]^5)*Log[\[Tau]/(2*r0)]) + (1/20 + \[Rho]^2/(210*\[Tau]^2) - \[Rho]/(60*\[Tau]) + \[Tau]/(6*\[Rho]) - (5*\[Tau]^2)/(6*\[Rho]^2) - (43*\[Tau]^3)/(20*\[Rho]^3) - 
    (91*\[Tau]^4)/(60*\[Rho]^4) - (31*\[Tau]^5)/(105*\[Rho]^5) + ((3*\[Tau]^3)/\[Rho]^3 + (7*\[Tau]^4)/\[Rho]^4 + (4*\[Tau]^5)/\[Rho]^5)*Log[\[Rho]/(2*r0)])*Log[1 + \[Tau]/\[Rho]])/\[Tau]


Lbar[2, 5, 2]  = (13/1680 - \[Rho]/(840*\[Tau]) + (46*\[Tau])/(315*\[Rho]) + (1889*\[Tau]^2)/(672*\[Rho]^2) + (809*\[Tau]^3)/(168*\[Rho]^3) + (1423*\[Tau]^4)/(560*\[Rho]^4) + (13*\[Tau]^5)/(56*\[Rho]^5) + 
  ((7*\[Tau]^3)/\[Rho]^3 + (28*\[Tau]^4)/\[Rho]^4 + (36*\[Tau]^5)/\[Rho]^5 + (15*\[Tau]^6)/\[Rho]^6)*Log[\[Rho]/(2*r0)]^2 + 
  (-1/6*\[Tau]/\[Rho] + (5*\[Tau]^2)/(4*\[Rho]^2) + (301*\[Tau]^3)/(60*\[Rho]^3) + (91*\[Tau]^4)/(15*\[Rho]^4) + (93*\[Tau]^5)/(35*\[Rho]^5) + (13*\[Tau]^6)/(56*\[Rho]^6))*Log[\[Tau]/(2*r0)] + 
  Log[\[Rho]/(2*r0)]*(\[Tau]/(3*\[Rho]) - (5*\[Tau]^2)/(2*\[Rho]^2) - (1201*\[Tau]^3)/(60*\[Rho]^3) - (1037*\[Tau]^4)/(30*\[Rho]^4) - (618*\[Tau]^5)/(35*\[Rho]^5) - (13*\[Tau]^6)/(56*\[Rho]^6) + 
    ((-7*\[Tau]^3)/\[Rho]^3 - (28*\[Tau]^4)/\[Rho]^4 - (36*\[Tau]^5)/\[Rho]^5 - (15*\[Tau]^6)/\[Rho]^6)*Log[\[Tau]/(2*r0)]) + 
  (1/30 + \[Rho]^2/(840*\[Tau]^2) - \[Rho]/(140*\[Tau]) + \[Tau]/(6*\[Rho]) - (5*\[Tau]^2)/(4*\[Rho]^2) - (301*\[Tau]^3)/(60*\[Rho]^3) - (91*\[Tau]^4)/(15*\[Rho]^4) - (93*\[Tau]^5)/(35*\[Rho]^5) - (13*\[Tau]^6)/(56*\[Rho]^6) + 
    ((7*\[Tau]^3)/\[Rho]^3 + (28*\[Tau]^4)/\[Rho]^4 + (36*\[Tau]^5)/\[Rho]^5 + (15*\[Tau]^6)/\[Rho]^6)*Log[\[Rho]/(2*r0)])*Log[1 + \[Tau]/\[Rho]])/\[Tau]




