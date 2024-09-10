(* ::Package:: *)

(* ::Chapter:: *)
(*Supplementary file to  "Gravitational-wave tails of memory " (David Trestini, Luc Blanchet)*)


(* ::Subtitle:: *)
(*Expression of the quadratic metric in the radiative construction*)


(* ::Text:: *)
(*Here, we give the relevant pieces of the quadratic metric in the radiative multipolar post-Minkowskian construction. We give it as a 3+1 projection for each relevant interaction, where h00 is the time-time component,  h0i is the time-space part and hij is the space-space part of the metric. The moments we consider are the ADM mass M, the ADM angular moment (or spin) Si and the canonical quadrupole moment of the radiative construction Mij (denoted with an overbar in the paper).*)
(**)
(*Note that these sources are associated the radiative construction and are given in terms of the radiative coordinates (unlike other computations, where the harmonic construction is used).*)


(* ::Section:: *)
(*Notations*)


(* ::Text:: *)
(*Time indices are designated as 0, and the only two possible free spatial indices are i and j, and are always upstairs. All other indices are dummy indices, and are downstairs when there is a minus sign in from, and upstairs otherwise (convention of the xTensor library).*)
(**)
(*Log[x] is the natural logarithm of x*)
(*r0 is the regularization constant associated to the finite part operator in the radiative construction.*)
(*G is the Newton's constant of gravity, and c is the speed of light*)
(*t is the time coordinate  and r is the radial coordinate (in radiative coordinates)*)
(**)
(*ni[i] is the unit 3-vectors pointing from the origin towards the field point*)
(*STFnL[a1, a2, ..., aN] is the symmetric trace-free part of the product of N unit 3-vectors (with the indices "a1", "a2", ..., "aN") pointing from the origin towards the field point, and defined in Footnote 1.*)
(*Metricdelta[a,b] is the two-index three-dimensional Kronecker symbol with indices "a" and "b"*)
(**)
(*M is the constant ADM mass, and SiVertBarj[a, b] is related to the constant ADM angular moment, and is defined in the Section III.A , right after Eq. 3.1, with index "a" before the vertical bar, and index "b" after. *)
(**)
(*barMij[a,b] is the canonical mass quadrupole moment associated to the radiative construction, defined by Eq. 2.4, with the two indices "a" and "b". Its N-th time derivative is denoted dt[N][barMij][a,b]. It is assumed to be evaluated at retarded time u = t - r, except when it is inside a head of IntegralOverQm or Antiderivative (introduced here-after).*)
(**)
(*Antiderivative[F], where F is some expression depending on barMij, is the antiderivative of F that cancels at negative infinity. More precisely, it is the integral over the variable s, spanning from 0 to positive infinity, of F(t-r-s). This type of term is illustrated by Eq. 3.13b.*)
(**)
(*IntegralOverBarQm[m][F], where F is some expression depending on barMij, is the integral over the variable  x, spanning from 1 to positive infinity, of Qbar[m](x,r) F(t-rx), where Qbar[m](x,r) is defined in Eq. 3.4. This type of term is illustrated by Eq. 3.12 and Eq. 3.13a.*)


(* ::Section:: *)
(*Formatting (for elegant display when using Mathematica)*)


(* ::Text:: *)
(*This section should be run so as to have an elegant and clear display when using Mathematica . It is however not necessary, as all notations are explained in the previous section . In our formatting, we display all indices downstairs, since they are space indices and their covariance is unimportant .*)


Format[dt[N_][barMij][L__]]:= StringReplace["\!\(\*SubsuperscriptBox[OverscriptBox[\(M\), \(_\)], \(Indices\), \((Number)\)]\)", {"Number"-> ToString[N, StandardForm],"Indices"-> ToString[Times@@(Flatten[{L}]/. Times[any_,-1]:> any), StandardForm]}]


Format[barMij[L__]]:= StringReplace["\!\(\*SubscriptBox[OverscriptBox[\(M\), \(_\)], \(Indices\)]\)", "Indices"-> ToString[Times@@(Flatten[{L}]/. Times[any_,-1]:> any), StandardForm]]


Format[SiVerticalBarj[i_,j_]]:= StringReplace["\!\(\*SubscriptBox[\(S\), \(IndexI | IndexJ\)]\)", {"IndexI"-> ToString[i/. Times[any_,-1]:> any, StandardForm],"IndexJ"-> ToString[j/. Times[any_,-1]:> any, StandardForm]}]


Format[Metricdelta[i_,j_]]:= StringReplace["\!\(\*SubscriptBox[\(\[Delta]\), \(IndexI\\\ IndexJ\)]\)", {"IndexI"-> ToString[i/. Times[any_,-1]:> any, StandardForm],"IndexJ"-> ToString[j/. Times[any_,-1]:> any, StandardForm]}]


Format[ni[i_]]:= StringReplace["\!\(\*SubscriptBox[n, \(index\)]\)", "index"-> ToString[i /. Times[any_,-1]:> any, StandardForm]] 


Format[STFnL[L__]]:=StringReplace["\!\(\*SubscriptBox[OverscriptBox[\(n\), \(^\)], \(Indices\)]\)", "Indices"-> ToString[Times@@(Flatten[{L}]/. Times[any_,-1]:> any), StandardForm]]


Format[IntegralOverBarQm[m_][expr_]]:= StringReplace["\!\(\*SuperscriptBox[SubscriptBox[\(\[Integral]\), \(1\)], \(\[Infinity]\)]\)dx \!\(\*SubscriptBox[OverscriptBox[\(Q\), \(_\)], \(IndexForQ\)]\)(x,r){Expression}(t-rx)", {"IndexForQ"-> ToString[m, StandardForm],"Expression"-> ToString[expr, StandardForm]}]


Format[Antiderivative[expr_]]:= StringReplace["\!\(\*SuperscriptBox[SubscriptBox[\(\[Integral]\), \(0\)], \(\[Infinity]\)]\)d\[Tau]{Expression}(t-r-\[Tau])", {"Expression"-> ToString[expr, StandardForm]}]


(* ::Section:: *)
(*Quadratic metric associated to the mass - mass interaction (eq. 3.2a)*)


h00MxM = (-3*M^2)/r^2


h0iMxM = 0


hijMxM = -((M^2*ni[i]*ni[j])/r^2)


(* ::Section:: *)
(*Quadratic metric associated to the mass - spin interaction (eq. 3.2b)*)


h00MxSi = 0


h0iMxSi = (2*M*ni[a]*SiVerticalBarj[i, -a])/r^3


hijMxSi = 0


(* ::Section:: *)
(*Quadratic metric associated to the mass - quadrupole interaction (eq. 3.3)*)


h00MxMij = (-21*M*barMij[-a, -b]*ni[a]*ni[b])/r^4 + 8*M*ni[a]*ni[b]*IntegralOverBarQm[2][dt[4][barMij][-a, -b]] - 
 (117*M*ni[a]*ni[b]*dt[1][barMij][-a, -b])/(5*r^3) + (23*M*ni[a]*ni[b]*dt[2][barMij][-a, -b])/(5*r^2) + 
 (118*M*ni[a]*ni[b]*dt[3][barMij][-a, -b])/(15*r)


h0iMxMij = 8*M*ni[a]*IntegralOverBarQm[1][dt[4][barMij][i, -a]] - (M*ni[a]*ni[b]*ni[i]*dt[1][barMij][-a, -b])/r^3 - 
 (5*M*ni[a]*dt[1][barMij][i, -a])/r^3 + (3*M*ni[a]*ni[b]*ni[i]*dt[2][barMij][-a, -b])/r^2 - 
 (107*M*ni[a]*dt[2][barMij][i, -a])/(15*r^2) + (M*ni[a]*ni[b]*ni[i]*dt[3][barMij][-a, -b])/r + 
 (43*M*ni[a]*dt[3][barMij][i, -a])/(15*r)


hijMxMij = -((M*barMij[i, j])/r^4) - (M*barMij[-a, -b]*Metricdelta[i, j]*ni[a]*ni[b])/(2*r^4) + (3*M*barMij[j, -a]*ni[a]*ni[i])/r^4 + 
 (3*M*barMij[i, -a]*ni[a]*ni[j])/r^4 - (15*M*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[j])/(2*r^4) + 
 8*M*IntegralOverBarQm[0][dt[4][barMij][i, j]] - (M*Metricdelta[i, j]*ni[a]*ni[b]*dt[1][barMij][-a, -b])/(2*r^3) - 
 (15*M*ni[a]*ni[b]*ni[i]*ni[j]*dt[1][barMij][-a, -b])/(2*r^3) + (3*M*ni[a]*ni[j]*dt[1][barMij][i, -a])/r^3 - 
 (M*dt[1][barMij][i, j])/r^3 + (3*M*ni[a]*ni[i]*dt[1][barMij][j, -a])/r^3 - 
 (2*M*Metricdelta[i, j]*ni[a]*ni[b]*dt[2][barMij][-a, -b])/r^2 - (3*M*ni[a]*ni[b]*ni[i]*ni[j]*dt[2][barMij][-a, -b])/r^2 + 
 (3*M*ni[a]*ni[j]*dt[2][barMij][i, -a])/r^2 - (4*M*dt[2][barMij][i, j])/r^2 + (3*M*ni[a]*ni[i]*dt[2][barMij][j, -a])/r^2 - 
 (M*Metricdelta[i, j]*ni[a]*ni[b]*dt[3][barMij][-a, -b])/(2*r) - (M*ni[a]*ni[b]*ni[i]*ni[j]*dt[3][barMij][-a, -b])/(2*r) + 
 (2*M*ni[a]*ni[j]*dt[3][barMij][i, -a])/r - (107*M*dt[3][barMij][i, j])/(15*r) + (2*M*ni[a]*ni[i]*dt[3][barMij][j, -a])/r


(* ::Section:: *)
(*Quadratic metric associated to the spin - quadrupole interaction (eq. 3.6)*)


h00SixMij = (-6*ni[a]*ni[b]*SiVerticalBarj[-a, k]*dt[1][barMij][-b, -k])/r^4 - (6*ni[a]*ni[b]*SiVerticalBarj[-a, k]*dt[2][barMij][-b, -k])/
  r^3 - (4*ni[a]*ni[b]*SiVerticalBarj[-a, k]*dt[3][barMij][-b, -k])/r^2 - 
 (4*ni[a]*ni[b]*SiVerticalBarj[-a, k]*dt[4][barMij][-b, -k])/(3*r)


h0iSixMij = (barMij[i, b]*ni[a]*SiVerticalBarj[-a, -b])/r^5 - (9*barMij[-a, k]*ni[a]*ni[b]*ni[i]*SiVerticalBarj[-b, -k])/r^5 + 
 (4*barMij[-a, b]*ni[a]*SiVerticalBarj[i, -b])/r^5 - (3*barMij[-a, -b]*ni[a]*ni[b]*ni[k]*SiVerticalBarj[i, -k])/(2*r^5) + 
 (4*ni[a]*SiVerticalBarj[i, b]*dt[1][barMij][-a, -b])/r^4 - (9*ni[a]*ni[b]*ni[i]*SiVerticalBarj[-a, k]*dt[1][barMij][-b, -k])/
  r^4 - (3*ni[a]*ni[b]*ni[k]*SiVerticalBarj[i, -a]*dt[1][barMij][-b, -k])/(2*r^4) + 
 (ni[a]*SiVerticalBarj[-a, b]*dt[1][barMij][i, -b])/r^4 + (5*ni[a]*SiVerticalBarj[i, b]*dt[2][barMij][-a, -b])/r^3 - 
 (3*ni[a]*ni[b]*ni[i]*SiVerticalBarj[-a, k]*dt[2][barMij][-b, -k])/r^3 - 
 (23*ni[a]*ni[b]*ni[k]*SiVerticalBarj[i, -a]*dt[2][barMij][-b, -k])/(2*r^3) + 
 (ni[a]*SiVerticalBarj[-a, b]*dt[2][barMij][i, -b])/r^3 + (5*ni[a]*SiVerticalBarj[i, b]*dt[3][barMij][-a, -b])/(3*r^2) - 
 (5*ni[a]*ni[b]*ni[k]*SiVerticalBarj[i, -a]*dt[3][barMij][-b, -k])/r^2 - (4*ni[a]*SiVerticalBarj[-a, b]*dt[3][barMij][i, -b])/
  (3*r^2) - (5*ni[a]*ni[b]*ni[k]*SiVerticalBarj[i, -a]*dt[4][barMij][-b, -k])/(6*r) - 
 (4*ni[a]*SiVerticalBarj[-a, b]*dt[4][barMij][i, -b])/(3*r)


hijSixMij = (-2*ni[a]*ni[j]*SiVerticalBarj[i, b]*dt[1][barMij][-a, -b])/r^4 - (2*ni[a]*ni[i]*SiVerticalBarj[j, b]*dt[1][barMij][-a, -b])/
  r^4 - (Metricdelta[i, j]*ni[a]*ni[b]*SiVerticalBarj[-a, k]*dt[1][barMij][-b, -k])/r^4 + 
 (9*ni[a]*ni[b]*ni[i]*ni[j]*SiVerticalBarj[-a, k]*dt[1][barMij][-b, -k])/r^4 + 
 (9*ni[a]*ni[b]*ni[j]*ni[k]*SiVerticalBarj[i, -a]*dt[1][barMij][-b, -k])/(2*r^4) + 
 (9*ni[a]*ni[b]*ni[i]*ni[k]*SiVerticalBarj[j, -a]*dt[1][barMij][-b, -k])/(2*r^4) + 
 (SiVerticalBarj[j, a]*dt[1][barMij][i, -a])/r^4 - (3*ni[a]*ni[j]*SiVerticalBarj[-a, b]*dt[1][barMij][i, -b])/r^4 - 
 (ni[a]*ni[b]*SiVerticalBarj[j, -a]*dt[1][barMij][i, -b])/r^4 + (SiVerticalBarj[i, a]*dt[1][barMij][j, -a])/r^4 - 
 (3*ni[a]*ni[i]*SiVerticalBarj[-a, b]*dt[1][barMij][j, -b])/r^4 - (ni[a]*ni[b]*SiVerticalBarj[i, -a]*dt[1][barMij][j, -b])/
  r^4 - (2*ni[a]*ni[j]*SiVerticalBarj[i, b]*dt[2][barMij][-a, -b])/r^3 - 
 (2*ni[a]*ni[i]*SiVerticalBarj[j, b]*dt[2][barMij][-a, -b])/r^3 - 
 (Metricdelta[i, j]*ni[a]*ni[b]*SiVerticalBarj[-a, k]*dt[2][barMij][-b, -k])/r^3 + 
 (9*ni[a]*ni[b]*ni[i]*ni[j]*SiVerticalBarj[-a, k]*dt[2][barMij][-b, -k])/r^3 + 
 (9*ni[a]*ni[b]*ni[j]*ni[k]*SiVerticalBarj[i, -a]*dt[2][barMij][-b, -k])/(2*r^3) + 
 (9*ni[a]*ni[b]*ni[i]*ni[k]*SiVerticalBarj[j, -a]*dt[2][barMij][-b, -k])/(2*r^3) + 
 (SiVerticalBarj[j, a]*dt[2][barMij][i, -a])/r^3 - (3*ni[a]*ni[j]*SiVerticalBarj[-a, b]*dt[2][barMij][i, -b])/r^3 - 
 (ni[a]*ni[b]*SiVerticalBarj[j, -a]*dt[2][barMij][i, -b])/r^3 + (SiVerticalBarj[i, a]*dt[2][barMij][j, -a])/r^3 - 
 (3*ni[a]*ni[i]*SiVerticalBarj[-a, b]*dt[2][barMij][j, -b])/r^3 - (ni[a]*ni[b]*SiVerticalBarj[i, -a]*dt[2][barMij][j, -b])/
  r^3 - (2*ni[a]*ni[j]*SiVerticalBarj[i, b]*dt[3][barMij][-a, -b])/(3*r^2) - 
 (2*ni[a]*ni[i]*SiVerticalBarj[j, b]*dt[3][barMij][-a, -b])/(3*r^2) + 
 (8*Metricdelta[i, j]*ni[a]*ni[b]*SiVerticalBarj[-a, k]*dt[3][barMij][-b, -k])/(3*r^2) + 
 (10*ni[a]*ni[b]*ni[i]*ni[j]*SiVerticalBarj[-a, k]*dt[3][barMij][-b, -k])/(3*r^2) + 
 (5*ni[a]*ni[b]*ni[j]*ni[k]*SiVerticalBarj[i, -a]*dt[3][barMij][-b, -k])/(3*r^2) + 
 (5*ni[a]*ni[b]*ni[i]*ni[k]*SiVerticalBarj[j, -a]*dt[3][barMij][-b, -k])/(3*r^2) + 
 (7*SiVerticalBarj[j, a]*dt[3][barMij][i, -a])/(3*r^2) - (13*ni[a]*ni[j]*SiVerticalBarj[-a, b]*dt[3][barMij][i, -b])/(3*r^2) - 
 (10*ni[a]*ni[b]*SiVerticalBarj[j, -a]*dt[3][barMij][i, -b])/(3*r^2) + (7*SiVerticalBarj[i, a]*dt[3][barMij][j, -a])/(3*r^2) - 
 (13*ni[a]*ni[i]*SiVerticalBarj[-a, b]*dt[3][barMij][j, -b])/(3*r^2) - 
 (10*ni[a]*ni[b]*SiVerticalBarj[i, -a]*dt[3][barMij][j, -b])/(3*r^2) + 
 (Metricdelta[i, j]*ni[a]*ni[b]*SiVerticalBarj[-a, k]*dt[4][barMij][-b, -k])/r + 
 (ni[a]*ni[b]*ni[i]*ni[j]*SiVerticalBarj[-a, k]*dt[4][barMij][-b, -k])/(3*r) + 
 (ni[a]*ni[b]*ni[j]*ni[k]*SiVerticalBarj[i, -a]*dt[4][barMij][-b, -k])/(6*r) + 
 (ni[a]*ni[b]*ni[i]*ni[k]*SiVerticalBarj[j, -a]*dt[4][barMij][-b, -k])/(6*r) - 
 (4*ni[a]*ni[j]*SiVerticalBarj[-a, b]*dt[4][barMij][i, -b])/(3*r) - (ni[a]*ni[b]*SiVerticalBarj[j, -a]*dt[4][barMij][i, -b])/
  r - (4*ni[a]*ni[i]*SiVerticalBarj[-a, b]*dt[4][barMij][j, -b])/(3*r) - 
 (ni[a]*ni[b]*SiVerticalBarj[i, -a]*dt[4][barMij][j, -b])/r


(* ::Section:: *)
(*Quadratic metric associated to the quadrupole - quadrupole interaction (not published in paper)*)


h00MijxMij = (4*Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][a, b]])/(5*r) - 
 (63*barMij[-a, -b]*barMij[-k, -l]*ni[a]*ni[b]*ni[k]*ni[l])/(4*r^6) - 
 (4*IntegralOverBarQm[0][dt[3][barMij][-a, -b]*dt[3][barMij][a, b]])/21 - 
 (8*ni[a]*ni[b]*IntegralOverBarQm[0][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/35 - 
 (4*IntegralOverBarQm[1][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][a, b]]])/(7*r) - 
 (24*ni[a]*ni[b]*IntegralOverBarQm[1][Antiderivative[dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]]])/(35*r) + 
 (16*ni[a]*ni[b]*IntegralOverBarQm[2][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/21 + 
 (2*ni[a]*ni[b]*ni[k]*ni[l]*IntegralOverBarQm[2][dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]])/21 - 
 (4*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][a, b]]])/(3*r) + 
 (56*ni[a]*ni[b]*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]]])/(15*r) + 
 (2*ni[a]*ni[b]*ni[k]*ni[l]*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]]])/(3*r) + 
 (2*IntegralOverBarQm[4][dt[3][barMij][-a, -b]*dt[3][barMij][a, b]])/11 - 
 (172*ni[a]*ni[b]*IntegralOverBarQm[4][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/385 - 
 (19*ni[a]*ni[b]*ni[k]*ni[l]*IntegralOverBarQm[4][dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]])/77 - 
 (2*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][a, b]]])/(21*r) + 
 (20*ni[a]*ni[b]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]]])/(21*r) - 
 (5*ni[a]*ni[b]*ni[k]*ni[l]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]]])/(3*r) + 
 (2*IntegralOverBarQm[6][dt[3][barMij][-a, -b]*dt[3][barMij][a, b]])/231 - 
 (20*ni[a]*ni[b]*IntegralOverBarQm[6][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/231 + 
 (5*ni[a]*ni[b]*ni[k]*ni[l]*IntegralOverBarQm[6][dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]])/33 - 
 (7*dt[1][barMij][-a, -b]*dt[1][barMij][a, b])/(4*r^4) + (29*ni[a]*ni[b]*dt[1][barMij][-a, k]*dt[1][barMij][-b, -k])/(2*r^4) - 
 (63*barMij[-a, -b]*ni[a]*ni[b]*ni[k]*ni[l]*dt[1][barMij][-k, -l])/(2*r^5) - 
 (291*ni[a]*ni[b]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[1][barMij][-k, -l])/(8*r^4) + 
 (5*barMij[a, b]*dt[2][barMij][-a, -b])/r^4 + (3*dt[1][barMij][a, b]*dt[2][barMij][-a, -b])/(2*r^3) + 
 (dt[2][barMij][-a, -b]*dt[2][barMij][a, b])/(2*r^2) - (42*barMij[-a, k]*ni[a]*ni[b]*dt[2][barMij][-b, -k])/r^4 - 
 (13*ni[a]*ni[b]*dt[1][barMij][-a, k]*dt[2][barMij][-b, -k])/r^3 - (3*ni[a]*ni[b]*dt[2][barMij][-a, k]*dt[2][barMij][-b, -k])/
  r^2 + (54*barMij[-a, -b]*ni[a]*ni[b]*ni[k]*ni[l]*dt[2][barMij][-k, -l])/r^4 + 
 (51*ni[a]*ni[b]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[2][barMij][-k, -l])/(4*r^3) + 
 (25*ni[a]*ni[b]*ni[k]*ni[l]*dt[2][barMij][-a, -b]*dt[2][barMij][-k, -l])/(4*r^2) + 
 (23*barMij[a, b]*dt[3][barMij][-a, -b])/(7*r^3) + (13*dt[1][barMij][a, b]*dt[3][barMij][-a, -b])/(21*r^2) - 
 (17*dt[2][barMij][a, b]*dt[3][barMij][-a, -b])/(42*r) - (174*barMij[-a, k]*ni[a]*ni[b]*dt[3][barMij][-b, -k])/(7*r^3) - 
 (172*ni[a]*ni[b]*dt[1][barMij][-a, k]*dt[3][barMij][-b, -k])/(21*r^2) + 
 (29*ni[a]*ni[b]*dt[2][barMij][-a, k]*dt[3][barMij][-b, -k])/(21*r) + 
 (69*barMij[-a, -b]*ni[a]*ni[b]*ni[k]*ni[l]*dt[3][barMij][-k, -l])/(2*r^3) + 
 (83*ni[a]*ni[b]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[3][barMij][-k, -l])/(6*r^2) + 
 (37*ni[a]*ni[b]*ni[k]*ni[l]*dt[2][barMij][-a, -b]*dt[3][barMij][-k, -l])/(12*r) + 
 (34*barMij[a, b]*dt[4][barMij][-a, -b])/(21*r^2) + (229*dt[1][barMij][a, b]*dt[4][barMij][-a, -b])/(210*r) - 
 (172*barMij[-a, k]*ni[a]*ni[b]*dt[4][barMij][-b, -k])/(21*r^2) - (47*ni[a]*ni[b]*dt[1][barMij][-a, k]*dt[4][barMij][-b, -k])/
  (21*r) + (25*barMij[-a, -b]*ni[a]*ni[b]*ni[k]*ni[l]*dt[4][barMij][-k, -l])/(3*r^2) + 
 (29*ni[a]*ni[b]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[4][barMij][-k, -l])/(12*r) + 
 (52*barMij[a, b]*dt[5][barMij][-a, -b])/(105*r) - (34*barMij[-a, k]*ni[a]*ni[b]*dt[5][barMij][-b, -k])/(21*r) + 
 (5*barMij[-a, -b]*ni[a]*ni[b]*ni[k]*ni[l]*dt[5][barMij][-k, -l])/(6*r)


h0iMijxMij = (-4*Antiderivative[dt[2][barMij][i, b]*dt[3][barMij][-a, -b]]*ni[a])/(5*r^2) + 
 (4*Antiderivative[dt[2][barMij][-a, b]*dt[3][barMij][i, -b]]*ni[a])/(5*r^2) + 
 (2*Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][a, b]]*ni[i])/(5*r) - 
 (24*ni[a]*ni[b]*ni[i]*IntegralOverBarQm[1][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/35 + 
 (24*ni[a]*IntegralOverBarQm[1][dt[3][barMij][-a, -b]*dt[3][barMij][i, b]])/35 - 
 (24*ni[a]*ni[b]*ni[i]*IntegralOverBarQm[2][Antiderivative[dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]]])/(7*r) + 
 (24*ni[a]*IntegralOverBarQm[2][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][i, b]]])/(7*r) + 
 (28*ni[a]*ni[b]*ni[i]*IntegralOverBarQm[3][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/45 - 
 (2*ni[a]*ni[b]*ni[k]*IntegralOverBarQm[3][dt[3][barMij][-b, -k]*dt[3][barMij][i, -a]])/9 - 
 (28*ni[a]*IntegralOverBarQm[3][dt[3][barMij][-a, -b]*dt[3][barMij][i, b]])/45 + 
 (2*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*IntegralOverBarQm[3][dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]])/9 - 
 (4*ni[a]*ni[b]*ni[i]*IntegralOverBarQm[4][Antiderivative[dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]]])/(7*r) - 
 (2*ni[a]*ni[b]*ni[k]*IntegralOverBarQm[4][Antiderivative[dt[3][barMij][-b, -k]*dt[3][barMij][i, -a]]])/r + 
 (4*ni[a]*IntegralOverBarQm[4][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][i, b]]])/(7*r) + 
 (2*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*IntegralOverBarQm[4][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]]])/r + 
 (4*ni[a]*ni[b]*ni[i]*IntegralOverBarQm[5][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/63 + 
 (2*ni[a]*ni[b]*ni[k]*IntegralOverBarQm[5][dt[3][barMij][-b, -k]*dt[3][barMij][i, -a]])/9 - 
 (4*ni[a]*IntegralOverBarQm[5][dt[3][barMij][-a, -b]*dt[3][barMij][i, b]])/63 - 
 (2*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*IntegralOverBarQm[5][dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]])/9 - 
 (3*barMij[i, b]*ni[a]*dt[1][barMij][-a, -b])/(2*r^5) + (barMij[a, b]*ni[i]*dt[1][barMij][-a, -b])/(2*r^5) + 
 (ni[i]*dt[1][barMij][-a, -b]*dt[1][barMij][a, b])/(2*r^4) - (9*barMij[-a, k]*ni[a]*ni[b]*ni[i]*dt[1][barMij][-b, -k])/
  (2*r^5) + (9*barMij[i, -a]*ni[a]*ni[b]*ni[k]*dt[1][barMij][-b, -k])/r^5 - 
 (9*ni[a]*ni[b]*ni[i]*dt[1][barMij][-a, k]*dt[1][barMij][-b, -k])/(2*r^4) - 
 (9*ni[a]*ni[b]*ni[k]*dt[1][barMij][-b, -k]*dt[1][barMij][i, -a])/(4*r^4) + (3*barMij[-a, b]*ni[a]*dt[1][barMij][i, -b])/r^5 + 
 (3*ni[a]*dt[1][barMij][-a, -b]*dt[1][barMij][i, b])/(2*r^4) - (45*barMij[-a, -b]*ni[a]*ni[b]*ni[k]*dt[1][barMij][i, -k])/
  (4*r^5) - (9*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*dt[1][barMij][-k, -l])/(2*r^5) - 
 (9*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[1][barMij][-k, -l])/(2*r^4) - 
 (3*barMij[i, b]*ni[a]*dt[2][barMij][-a, -b])/(2*r^4) + (barMij[a, b]*ni[i]*dt[2][barMij][-a, -b])/(2*r^4) + 
 (ni[i]*dt[1][barMij][a, b]*dt[2][barMij][-a, -b])/(2*r^3) - (4*ni[a]*dt[1][barMij][i, b]*dt[2][barMij][-a, -b])/r^3 - 
 (9*barMij[-a, k]*ni[a]*ni[b]*ni[i]*dt[2][barMij][-b, -k])/(2*r^4) + (9*barMij[i, -a]*ni[a]*ni[b]*ni[k]*dt[2][barMij][-b, -k])/
  r^4 - (7*ni[a]*ni[b]*ni[i]*dt[1][barMij][-a, k]*dt[2][barMij][-b, -k])/(2*r^3) + 
 (47*ni[a]*ni[b]*ni[k]*dt[1][barMij][i, -a]*dt[2][barMij][-b, -k])/(4*r^3) + 
 (ni[a]*ni[b]*ni[i]*dt[2][barMij][-a, k]*dt[2][barMij][-b, -k])/r^2 + 
 (7*ni[a]*ni[b]*ni[k]*dt[2][barMij][-b, -k]*dt[2][barMij][i, -a])/(2*r^2) + (3*barMij[-a, b]*ni[a]*dt[2][barMij][i, -b])/r^4 + 
 (9*ni[a]*dt[1][barMij][-a, b]*dt[2][barMij][i, -b])/(2*r^3) - (ni[a]*dt[2][barMij][-a, -b]*dt[2][barMij][i, b])/r^2 - 
 (45*barMij[-a, -b]*ni[a]*ni[b]*ni[k]*dt[2][barMij][i, -k])/(4*r^4) - 
 (21*ni[a]*ni[b]*ni[k]*dt[1][barMij][-a, -b]*dt[2][barMij][i, -k])/(2*r^3) - 
 (9*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*dt[2][barMij][-k, -l])/(2*r^4) - 
 (6*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[2][barMij][-k, -l])/r^3 - 
 (3*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*dt[2][barMij][-a, -b]*dt[2][barMij][-k, -l])/(2*r^2) - 
 (7*barMij[i, b]*ni[a]*dt[3][barMij][-a, -b])/(6*r^3) + (2*barMij[a, b]*ni[i]*dt[3][barMij][-a, -b])/(3*r^3) + 
 (51*ni[i]*dt[1][barMij][a, b]*dt[3][barMij][-a, -b])/(70*r^2) - (106*ni[a]*dt[1][barMij][i, b]*dt[3][barMij][-a, -b])/
  (105*r^2) - (17*ni[i]*dt[2][barMij][a, b]*dt[3][barMij][-a, -b])/(42*r) + 
 (44*ni[a]*dt[2][barMij][i, b]*dt[3][barMij][-a, -b])/(35*r) - (17*barMij[-a, k]*ni[a]*ni[b]*ni[i]*dt[3][barMij][-b, -k])/
  (6*r^3) + (29*barMij[i, -a]*ni[a]*ni[b]*ni[k]*dt[3][barMij][-b, -k])/(6*r^3) - 
 (ni[a]*ni[b]*ni[i]*dt[1][barMij][-a, k]*dt[3][barMij][-b, -k])/(3*r^2) + 
 (16*ni[a]*ni[b]*ni[k]*dt[1][barMij][i, -a]*dt[3][barMij][-b, -k])/(3*r^2) + 
 (28*ni[a]*ni[b]*ni[i]*dt[2][barMij][-a, k]*dt[3][barMij][-b, -k])/(15*r) + 
 (11*ni[a]*ni[b]*ni[k]*dt[2][barMij][i, -a]*dt[3][barMij][-b, -k])/(10*r) - 
 (8*barMij[-a, b]*ni[a]*dt[3][barMij][i, -b])/(3*r^3) - (274*ni[a]*dt[1][barMij][-a, b]*dt[3][barMij][i, -b])/(105*r^2) - 
 (61*ni[a]*dt[2][barMij][-a, b]*dt[3][barMij][i, -b])/(35*r) + (67*barMij[-a, -b]*ni[a]*ni[b]*ni[k]*dt[3][barMij][i, -k])/
  (12*r^3) + (10*ni[a]*ni[b]*ni[k]*dt[1][barMij][-a, -b]*dt[3][barMij][i, -k])/(3*r^2) + 
 (13*ni[a]*ni[b]*ni[k]*dt[2][barMij][-a, -b]*dt[3][barMij][i, -k])/(5*r) - 
 (7*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*dt[3][barMij][-k, -l])/(4*r^3) - 
 (7*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[3][barMij][-k, -l])/(4*r^2) - 
 (37*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*dt[2][barMij][-a, -b]*dt[3][barMij][-k, -l])/(60*r) - 
 (64*barMij[i, b]*ni[a]*dt[4][barMij][-a, -b])/(105*r^2) + (51*barMij[a, b]*ni[i]*dt[4][barMij][-a, -b])/(70*r^2) + 
 (229*ni[i]*dt[1][barMij][a, b]*dt[4][barMij][-a, -b])/(210*r) + (11*ni[a]*dt[1][barMij][i, b]*dt[4][barMij][-a, -b])/(35*r) - 
 (4*barMij[-a, k]*ni[a]*ni[b]*ni[i]*dt[4][barMij][-b, -k])/(3*r^2) + 
 (11*barMij[i, -a]*ni[a]*ni[b]*ni[k]*dt[4][barMij][-b, -k])/(6*r^2) - 
 (11*ni[a]*ni[b]*ni[i]*dt[1][barMij][-a, k]*dt[4][barMij][-b, -k])/(30*r) + 
 (23*ni[a]*ni[b]*ni[k]*dt[1][barMij][i, -a]*dt[4][barMij][-b, -k])/(20*r) - 
 (211*barMij[-a, b]*ni[a]*dt[4][barMij][i, -b])/(105*r^2) - (153*ni[a]*dt[1][barMij][-a, b]*dt[4][barMij][i, -b])/(70*r) + 
 (10*barMij[-a, -b]*ni[a]*ni[b]*ni[k]*dt[4][barMij][i, -k])/(3*r^2) + 
 (7*ni[a]*ni[b]*ni[k]*dt[1][barMij][-a, -b]*dt[4][barMij][i, -k])/(5*r) - 
 (barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*dt[4][barMij][-k, -l])/(4*r^2) - 
 (2*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[4][barMij][-k, -l])/(15*r) - 
 (17*barMij[i, b]*ni[a]*dt[5][barMij][-a, -b])/(70*r) + (52*barMij[a, b]*ni[i]*dt[5][barMij][-a, -b])/(105*r) - 
 (7*barMij[-a, k]*ni[a]*ni[b]*ni[i]*dt[5][barMij][-b, -k])/(30*r) + (3*barMij[i, -a]*ni[a]*ni[b]*ni[k]*dt[5][barMij][-b, -k])/
  (10*r) - (8*barMij[-a, b]*ni[a]*dt[5][barMij][i, -b])/(7*r) + (11*barMij[-a, -b]*ni[a]*ni[b]*ni[k]*dt[5][barMij][i, -k])/
  (20*r) - (barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[k]*ni[l]*dt[5][barMij][-k, -l])/(60*r)


hijMijxMij = (barMij[i, a]*barMij[j, -a])/(2*r^6) + (barMij[-a, -b]*barMij[a, b]*Metricdelta[i, j])/(4*r^6) - 
 (barMij[-a, -b]*barMij[i, j]*ni[a]*ni[b])/r^6 - (barMij[i, -a]*barMij[j, -b]*ni[a]*ni[b])/(2*r^6) - 
 (5*barMij[-a, k]*barMij[-b, -k]*Metricdelta[i, j]*ni[a]*ni[b])/(2*r^6) - (7*barMij[-a, -b]*barMij[j, b]*ni[a]*ni[i])/
  (2*r^6) - (7*barMij[-a, -b]*barMij[i, b]*ni[a]*ni[j])/(2*r^6) - (5*barMij[-a, -b]*barMij[a, b]*ni[i]*ni[j])/(4*r^6) + 
 (15*barMij[-a, k]*barMij[-b, -k]*ni[a]*ni[b]*ni[i]*ni[j])/r^6 + (6*barMij[-b, -k]*barMij[j, -a]*ni[a]*ni[b]*ni[i]*ni[k])/
  r^6 + (6*barMij[-b, -k]*barMij[i, -a]*ni[a]*ni[b]*ni[j]*ni[k])/r^6 + 
 (3*barMij[-a, -b]*barMij[-k, -l]*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l])/(2*r^6) - 
 (75*barMij[-a, -b]*barMij[-k, -l]*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l])/(4*r^6) - 
 (8*Metricdelta[i, j]*IntegralOverBarQm[0][dt[3][barMij][-a, -b]*dt[3][barMij][a, b]])/35 + 
 (44*ni[i]*ni[j]*IntegralOverBarQm[0][dt[3][barMij][-a, -b]*dt[3][barMij][a, b]])/105 + 
 (8*Metricdelta[i, j]*ni[a]*ni[b]*IntegralOverBarQm[0][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/35 - 
 (8*ni[a]*ni[j]*IntegralOverBarQm[0][dt[3][barMij][-a, -b]*dt[3][barMij][i, b]])/35 + 
 (16*IntegralOverBarQm[0][dt[3][barMij][i, a]*dt[3][barMij][j, -a]])/35 - 
 (8*ni[a]*ni[i]*IntegralOverBarQm[0][dt[3][barMij][-a, -b]*dt[3][barMij][j, b]])/35 - 
 (24*Metricdelta[i, j]*IntegralOverBarQm[1][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][a, b]]])/(35*r) + 
 (44*ni[i]*ni[j]*IntegralOverBarQm[1][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][a, b]]])/(35*r) + 
 (24*Metricdelta[i, j]*ni[a]*ni[b]*IntegralOverBarQm[1][Antiderivative[dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]]])/(35*r) - 
 (24*ni[a]*ni[j]*IntegralOverBarQm[1][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][i, b]]])/(35*r) + 
 (48*IntegralOverBarQm[1][Antiderivative[dt[3][barMij][i, a]*dt[3][barMij][j, -a]]])/(35*r) - 
 (24*ni[a]*ni[i]*IntegralOverBarQm[1][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][j, b]]])/(35*r) + 
 (20*Metricdelta[i, j]*IntegralOverBarQm[2][dt[3][barMij][-a, -b]*dt[3][barMij][a, b]])/63 - 
 (20*ni[i]*ni[j]*IntegralOverBarQm[2][dt[3][barMij][-a, -b]*dt[3][barMij][a, b]])/63 - 
 (16*Metricdelta[i, j]*ni[a]*ni[b]*IntegralOverBarQm[2][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/63 - 
 (32*ni[a]*ni[b]*ni[i]*ni[j]*IntegralOverBarQm[2][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/21 + 
 (2*ni[a]*ni[b]*ni[j]*ni[k]*IntegralOverBarQm[2][dt[3][barMij][-b, -k]*dt[3][barMij][i, -a]])/21 + 
 (52*ni[a]*ni[j]*IntegralOverBarQm[2][dt[3][barMij][-a, -b]*dt[3][barMij][i, b]])/63 - 
 (4*ni[a]*ni[b]*IntegralOverBarQm[2][dt[3][barMij][-a, -b]*dt[3][barMij][i, j]])/63 + 
 (2*ni[a]*ni[b]*ni[i]*ni[k]*IntegralOverBarQm[2][dt[3][barMij][-b, -k]*dt[3][barMij][j, -a]])/21 - 
 (40*IntegralOverBarQm[2][dt[3][barMij][i, a]*dt[3][barMij][j, -a]])/63 - 
 (8*ni[a]*ni[b]*IntegralOverBarQm[2][dt[3][barMij][i, -a]*dt[3][barMij][j, -b]])/63 + 
 (52*ni[a]*ni[i]*IntegralOverBarQm[2][dt[3][barMij][-a, -b]*dt[3][barMij][j, b]])/63 - 
 (2*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*IntegralOverBarQm[2][dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]])/21 + 
 (28*Metricdelta[i, j]*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][a, b]]])/(45*r) + 
 (32*ni[i]*ni[j]*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][a, b]]])/(45*r) - 
 (8*Metricdelta[i, j]*ni[a]*ni[b]*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]]])/(45*r) - 
 (32*ni[a]*ni[b]*ni[i]*ni[j]*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]]])/(3*r) + 
 (2*ni[a]*ni[b]*ni[j]*ni[k]*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][-b, -k]*dt[3][barMij][i, -a]]])/(3*r) + 
 (188*ni[a]*ni[j]*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][i, b]]])/(45*r) - 
 (4*ni[a]*ni[b]*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][i, j]]])/(9*r) + 
 (2*ni[a]*ni[b]*ni[i]*ni[k]*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][-b, -k]*dt[3][barMij][j, -a]]])/(3*r) - 
 (56*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][i, a]*dt[3][barMij][j, -a]]])/(45*r) - 
 (8*ni[a]*ni[b]*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][i, -a]*dt[3][barMij][j, -b]]])/(9*r) + 
 (188*ni[a]*ni[i]*IntegralOverBarQm[3][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][j, b]]])/(45*r) - 
 (2*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*IntegralOverBarQm[3][
    Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]]])/(3*r) - 
 (32*Metricdelta[i, j]*IntegralOverBarQm[4][dt[3][barMij][-a, -b]*dt[3][barMij][a, b]])/385 - 
 (38*ni[i]*ni[j]*IntegralOverBarQm[4][dt[3][barMij][-a, -b]*dt[3][barMij][a, b]])/385 - 
 (8*Metricdelta[i, j]*ni[a]*ni[b]*IntegralOverBarQm[4][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/385 + 
 (108*ni[a]*ni[b]*ni[i]*ni[j]*IntegralOverBarQm[4][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/77 - 
 (26*ni[a]*ni[b]*ni[j]*ni[k]*IntegralOverBarQm[4][dt[3][barMij][-b, -k]*dt[3][barMij][i, -a]])/77 - 
 (212*ni[a]*ni[j]*IntegralOverBarQm[4][dt[3][barMij][-a, -b]*dt[3][barMij][i, b]])/385 + 
 (8*ni[a]*ni[b]*IntegralOverBarQm[4][dt[3][barMij][-a, -b]*dt[3][barMij][i, j]])/77 - 
 (26*ni[a]*ni[b]*ni[i]*ni[k]*IntegralOverBarQm[4][dt[3][barMij][-b, -k]*dt[3][barMij][j, -a]])/77 + 
 (64*IntegralOverBarQm[4][dt[3][barMij][i, a]*dt[3][barMij][j, -a]])/385 + 
 (16*ni[a]*ni[b]*IntegralOverBarQm[4][dt[3][barMij][i, -a]*dt[3][barMij][j, -b]])/77 - 
 (212*ni[a]*ni[i]*IntegralOverBarQm[4][dt[3][barMij][-a, -b]*dt[3][barMij][j, b]])/385 + 
 (12*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*IntegralOverBarQm[4][dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]])/77 + 
 (5*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*IntegralOverBarQm[4][dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]])/11 + 
 (4*Metricdelta[i, j]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][a, b]]])/(63*r) + 
 (2*ni[i]*ni[j]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][a, b]]])/(63*r) - 
 (32*Metricdelta[i, j]*ni[a]*ni[b]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]]])/(63*r) - 
 (4*ni[a]*ni[b]*ni[i]*ni[j]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]]])/(3*r) - 
 (8*ni[a]*ni[b]*ni[j]*ni[k]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-b, -k]*dt[3][barMij][i, -a]]])/(3*r) + 
 (32*ni[a]*ni[j]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][i, b]]])/(63*r) + 
 (4*ni[a]*ni[b]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][i, j]]])/(9*r) - 
 (8*ni[a]*ni[b]*ni[i]*ni[k]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-b, -k]*dt[3][barMij][j, -a]]])/(3*r) - 
 (8*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][i, a]*dt[3][barMij][j, -a]]])/(63*r) + 
 (8*ni[a]*ni[b]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][i, -a]*dt[3][barMij][j, -b]]])/(9*r) + 
 (32*ni[a]*ni[i]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][j, b]]])/(63*r) + 
 (2*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*IntegralOverBarQm[5][
    Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]]])/(3*r) + 
 (5*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*IntegralOverBarQm[5][Antiderivative[dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]]])/r - 
 (4*Metricdelta[i, j]*IntegralOverBarQm[6][dt[3][barMij][-a, -b]*dt[3][barMij][a, b]])/693 - 
 (2*ni[i]*ni[j]*IntegralOverBarQm[6][dt[3][barMij][-a, -b]*dt[3][barMij][a, b]])/693 + 
 (32*Metricdelta[i, j]*ni[a]*ni[b]*IntegralOverBarQm[6][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/693 + 
 (4*ni[a]*ni[b]*ni[i]*ni[j]*IntegralOverBarQm[6][dt[3][barMij][-a, k]*dt[3][barMij][-b, -k]])/33 + 
 (8*ni[a]*ni[b]*ni[j]*ni[k]*IntegralOverBarQm[6][dt[3][barMij][-b, -k]*dt[3][barMij][i, -a]])/33 - 
 (32*ni[a]*ni[j]*IntegralOverBarQm[6][dt[3][barMij][-a, -b]*dt[3][barMij][i, b]])/693 - 
 (4*ni[a]*ni[b]*IntegralOverBarQm[6][dt[3][barMij][-a, -b]*dt[3][barMij][i, j]])/99 + 
 (8*ni[a]*ni[b]*ni[i]*ni[k]*IntegralOverBarQm[6][dt[3][barMij][-b, -k]*dt[3][barMij][j, -a]])/33 + 
 (8*IntegralOverBarQm[6][dt[3][barMij][i, a]*dt[3][barMij][j, -a]])/693 - 
 (8*ni[a]*ni[b]*IntegralOverBarQm[6][dt[3][barMij][i, -a]*dt[3][barMij][j, -b]])/99 - 
 (32*ni[a]*ni[i]*IntegralOverBarQm[6][dt[3][barMij][-a, -b]*dt[3][barMij][j, b]])/693 - 
 (2*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*IntegralOverBarQm[6][dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]])/33 - 
 (5*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*IntegralOverBarQm[6][dt[3][barMij][-a, -b]*dt[3][barMij][-k, -l]])/11 + 
 (barMij[a, b]*Metricdelta[i, j]*dt[1][barMij][-a, -b])/(2*r^5) - (barMij[i, j]*ni[a]*ni[b]*dt[1][barMij][-a, -b])/r^5 - 
 (7*barMij[j, b]*ni[a]*ni[i]*dt[1][barMij][-a, -b])/(2*r^5) - (7*barMij[i, b]*ni[a]*ni[j]*dt[1][barMij][-a, -b])/(2*r^5) - 
 (5*barMij[a, b]*ni[i]*ni[j]*dt[1][barMij][-a, -b])/(2*r^5) + (13*Metricdelta[i, j]*dt[1][barMij][-a, -b]*dt[1][barMij][a, b])/
  (12*r^4) - (7*ni[i]*ni[j]*dt[1][barMij][-a, -b]*dt[1][barMij][a, b])/(3*r^4) - 
 (5*barMij[-a, k]*Metricdelta[i, j]*ni[a]*ni[b]*dt[1][barMij][-b, -k])/r^5 + 
 (30*barMij[-a, k]*ni[a]*ni[b]*ni[i]*ni[j]*dt[1][barMij][-b, -k])/r^5 + 
 (6*barMij[j, -a]*ni[a]*ni[b]*ni[i]*ni[k]*dt[1][barMij][-b, -k])/r^5 + 
 (6*barMij[i, -a]*ni[a]*ni[b]*ni[j]*ni[k]*dt[1][barMij][-b, -k])/r^5 - 
 (37*Metricdelta[i, j]*ni[a]*ni[b]*dt[1][barMij][-a, k]*dt[1][barMij][-b, -k])/(6*r^4) + 
 (17*ni[a]*ni[b]*ni[i]*ni[j]*dt[1][barMij][-a, k]*dt[1][barMij][-b, -k])/r^4 + (barMij[j, a]*dt[1][barMij][i, -a])/(2*r^5) + 
 (5*ni[a]*ni[b]*ni[j]*ni[k]*dt[1][barMij][-b, -k]*dt[1][barMij][i, -a])/(4*r^4) - 
 (barMij[j, -a]*ni[a]*ni[b]*dt[1][barMij][i, -b])/(2*r^5) - (7*barMij[-a, b]*ni[a]*ni[j]*dt[1][barMij][i, -b])/(2*r^5) - 
 (17*ni[a]*ni[j]*dt[1][barMij][-a, -b]*dt[1][barMij][i, b])/(6*r^4) - (barMij[-a, -b]*ni[a]*ni[b]*dt[1][barMij][i, j])/r^5 + 
 (11*ni[a]*ni[b]*dt[1][barMij][-a, -b]*dt[1][barMij][i, j])/(6*r^4) + 
 (6*barMij[-a, -b]*ni[a]*ni[b]*ni[j]*ni[k]*dt[1][barMij][i, -k])/r^5 + (barMij[i, a]*dt[1][barMij][j, -a])/(2*r^5) + 
 (5*ni[a]*ni[b]*ni[i]*ni[k]*dt[1][barMij][-b, -k]*dt[1][barMij][j, -a])/(4*r^4) + 
 (dt[1][barMij][i, a]*dt[1][barMij][j, -a])/(3*r^4) - (barMij[i, -a]*ni[a]*ni[b]*dt[1][barMij][j, -b])/(2*r^5) - 
 (7*barMij[-a, b]*ni[a]*ni[i]*dt[1][barMij][j, -b])/(2*r^5) - (ni[a]*ni[b]*dt[1][barMij][i, -a]*dt[1][barMij][j, -b])/
  (3*r^4) - (17*ni[a]*ni[i]*dt[1][barMij][-a, -b]*dt[1][barMij][j, b])/(6*r^4) + 
 (6*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[k]*dt[1][barMij][j, -k])/r^5 + 
 (3*barMij[-a, -b]*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*dt[1][barMij][-k, -l])/r^5 - 
 (75*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*dt[1][barMij][-k, -l])/(2*r^5) + 
 (43*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[1][barMij][-k, -l])/(8*r^4) - 
 (35*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[1][barMij][-k, -l])/(2*r^4) - 
 (2*barMij[a, b]*Metricdelta[i, j]*dt[2][barMij][-a, -b])/(3*r^4) - (barMij[i, j]*ni[a]*ni[b]*dt[2][barMij][-a, -b])/(3*r^4) - 
 (5*barMij[j, b]*ni[a]*ni[i]*dt[2][barMij][-a, -b])/(3*r^4) - (5*barMij[i, b]*ni[a]*ni[j]*dt[2][barMij][-a, -b])/(3*r^4) - 
 (4*barMij[a, b]*ni[i]*ni[j]*dt[2][barMij][-a, -b])/(3*r^4) + (Metricdelta[i, j]*dt[1][barMij][a, b]*dt[2][barMij][-a, -b])/
  r^3 - (7*ni[i]*ni[j]*dt[1][barMij][a, b]*dt[2][barMij][-a, -b])/(2*r^3) - 
 (ni[a]*ni[j]*dt[1][barMij][i, b]*dt[2][barMij][-a, -b])/r^3 + (5*ni[a]*ni[b]*dt[1][barMij][i, j]*dt[2][barMij][-a, -b])/
  (2*r^3) - (ni[a]*ni[i]*dt[1][barMij][j, b]*dt[2][barMij][-a, -b])/r^3 + 
 (Metricdelta[i, j]*dt[2][barMij][-a, -b]*dt[2][barMij][a, b])/r^2 - (3*ni[i]*ni[j]*dt[2][barMij][-a, -b]*dt[2][barMij][a, b])/
  (2*r^2) + (10*barMij[-a, k]*Metricdelta[i, j]*ni[a]*ni[b]*dt[2][barMij][-b, -k])/(3*r^4) + 
 (14*barMij[-a, k]*ni[a]*ni[b]*ni[i]*ni[j]*dt[2][barMij][-b, -k])/r^4 + 
 (5*barMij[j, -a]*ni[a]*ni[b]*ni[i]*ni[k]*dt[2][barMij][-b, -k])/(2*r^4) + 
 (5*barMij[i, -a]*ni[a]*ni[b]*ni[j]*ni[k]*dt[2][barMij][-b, -k])/(2*r^4) - 
 (4*Metricdelta[i, j]*ni[a]*ni[b]*dt[1][barMij][-a, k]*dt[2][barMij][-b, -k])/r^3 + 
 (18*ni[a]*ni[b]*ni[i]*ni[j]*dt[1][barMij][-a, k]*dt[2][barMij][-b, -k])/r^3 - 
 (9*ni[a]*ni[b]*ni[j]*ni[k]*dt[1][barMij][i, -a]*dt[2][barMij][-b, -k])/(4*r^3) - 
 (9*ni[a]*ni[b]*ni[i]*ni[k]*dt[1][barMij][j, -a]*dt[2][barMij][-b, -k])/(4*r^3) - 
 (2*Metricdelta[i, j]*ni[a]*ni[b]*dt[2][barMij][-a, k]*dt[2][barMij][-b, -k])/r^2 + 
 (5*ni[a]*ni[b]*ni[i]*ni[j]*dt[2][barMij][-a, k]*dt[2][barMij][-b, -k])/r^2 + (7*barMij[j, a]*dt[2][barMij][i, -a])/(6*r^4) + 
 (dt[1][barMij][j, a]*dt[2][barMij][i, -a])/r^3 - (ni[a]*ni[b]*ni[j]*ni[k]*dt[2][barMij][-b, -k]*dt[2][barMij][i, -a])/
  (2*r^2) - (7*barMij[j, -a]*ni[a]*ni[b]*dt[2][barMij][i, -b])/(6*r^4) - (37*barMij[-a, b]*ni[a]*ni[j]*dt[2][barMij][i, -b])/
  (6*r^4) - (11*ni[a]*ni[j]*dt[1][barMij][-a, b]*dt[2][barMij][i, -b])/(2*r^3) - 
 (ni[a]*ni[b]*dt[1][barMij][j, -a]*dt[2][barMij][i, -b])/r^3 - (ni[a]*ni[j]*dt[2][barMij][-a, -b]*dt[2][barMij][i, b])/r^2 - 
 (29*barMij[-a, -b]*ni[a]*ni[b]*dt[2][barMij][i, j])/(6*r^4) - (2*ni[a]*ni[b]*dt[1][barMij][-a, -b]*dt[2][barMij][i, j])/r^3 + 
 (3*ni[a]*ni[b]*dt[2][barMij][-a, -b]*dt[2][barMij][i, j])/r^2 + 
 (37*barMij[-a, -b]*ni[a]*ni[b]*ni[j]*ni[k]*dt[2][barMij][i, -k])/(4*r^4) + 
 (9*ni[a]*ni[b]*ni[j]*ni[k]*dt[1][barMij][-a, -b]*dt[2][barMij][i, -k])/(2*r^3) + 
 (7*barMij[i, a]*dt[2][barMij][j, -a])/(6*r^4) + (dt[1][barMij][i, a]*dt[2][barMij][j, -a])/r^3 - 
 (ni[a]*ni[b]*ni[i]*ni[k]*dt[2][barMij][-b, -k]*dt[2][barMij][j, -a])/(2*r^2) + 
 (2*dt[2][barMij][i, a]*dt[2][barMij][j, -a])/r^2 - (7*barMij[i, -a]*ni[a]*ni[b]*dt[2][barMij][j, -b])/(6*r^4) - 
 (37*barMij[-a, b]*ni[a]*ni[i]*dt[2][barMij][j, -b])/(6*r^4) - (11*ni[a]*ni[i]*dt[1][barMij][-a, b]*dt[2][barMij][j, -b])/
  (2*r^3) - (ni[a]*ni[b]*dt[1][barMij][i, -a]*dt[2][barMij][j, -b])/r^3 - 
 (2*ni[a]*ni[b]*dt[2][barMij][i, -a]*dt[2][barMij][j, -b])/r^2 - (ni[a]*ni[i]*dt[2][barMij][-a, -b]*dt[2][barMij][j, b])/r^2 + 
 (37*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[k]*dt[2][barMij][j, -k])/(4*r^4) + 
 (9*ni[a]*ni[b]*ni[i]*ni[k]*dt[1][barMij][-a, -b]*dt[2][barMij][j, -k])/(2*r^3) - 
 (11*barMij[-a, -b]*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*dt[2][barMij][-k, -l])/(2*r^4) - 
 (16*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*dt[2][barMij][-k, -l])/r^4 + 
 (9*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[2][barMij][-k, -l])/(4*r^3) - 
 (27*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[2][barMij][-k, -l])/(2*r^3) + 
 (Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*dt[2][barMij][-a, -b]*dt[2][barMij][-k, -l])/r^2 - 
 (9*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*dt[2][barMij][-a, -b]*dt[2][barMij][-k, -l])/(4*r^2) - 
 (5*barMij[a, b]*Metricdelta[i, j]*dt[3][barMij][-a, -b])/(6*r^3) - (barMij[j, b]*ni[a]*ni[i]*dt[3][barMij][-a, -b])/(2*r^3) - 
 (barMij[i, b]*ni[a]*ni[j]*dt[3][barMij][-a, -b])/(2*r^3) - (barMij[a, b]*ni[i]*ni[j]*dt[3][barMij][-a, -b])/(2*r^3) - 
 (7*Metricdelta[i, j]*dt[1][barMij][a, b]*dt[3][barMij][-a, -b])/(10*r^2) - 
 (43*ni[i]*ni[j]*dt[1][barMij][a, b]*dt[3][barMij][-a, -b])/(30*r^2) - (ni[a]*ni[j]*dt[1][barMij][i, b]*dt[3][barMij][-a, -b])/
  (5*r^2) + (37*ni[a]*ni[b]*dt[1][barMij][i, j]*dt[3][barMij][-a, -b])/(30*r^2) - 
 (ni[a]*ni[i]*dt[1][barMij][j, b]*dt[3][barMij][-a, -b])/(5*r^2) + 
 (289*Metricdelta[i, j]*dt[2][barMij][a, b]*dt[3][barMij][-a, -b])/(315*r) - 
 (119*ni[i]*ni[j]*dt[2][barMij][a, b]*dt[3][barMij][-a, -b])/(90*r) + 
 (31*ni[a]*ni[j]*dt[2][barMij][i, b]*dt[3][barMij][-a, -b])/(45*r) + 
 (91*ni[a]*ni[b]*dt[2][barMij][i, j]*dt[3][barMij][-a, -b])/(90*r) + 
 (31*ni[a]*ni[i]*dt[2][barMij][j, b]*dt[3][barMij][-a, -b])/(45*r) + 
 (5*barMij[-a, k]*Metricdelta[i, j]*ni[a]*ni[b]*dt[3][barMij][-b, -k])/r^3 + 
 (4*barMij[-a, k]*ni[a]*ni[b]*ni[i]*ni[j]*dt[3][barMij][-b, -k])/r^3 + 
 (barMij[j, -a]*ni[a]*ni[b]*ni[i]*ni[k]*dt[3][barMij][-b, -k])/(2*r^3) + 
 (barMij[i, -a]*ni[a]*ni[b]*ni[j]*ni[k]*dt[3][barMij][-b, -k])/(2*r^3) + 
 (19*Metricdelta[i, j]*ni[a]*ni[b]*dt[1][barMij][-a, k]*dt[3][barMij][-b, -k])/(5*r^2) + 
 (77*ni[a]*ni[b]*ni[i]*ni[j]*dt[1][barMij][-a, k]*dt[3][barMij][-b, -k])/(15*r^2) - 
 (19*ni[a]*ni[b]*ni[j]*ni[k]*dt[1][barMij][i, -a]*dt[3][barMij][-b, -k])/(15*r^2) - 
 (19*ni[a]*ni[b]*ni[i]*ni[k]*dt[1][barMij][j, -a]*dt[3][barMij][-b, -k])/(15*r^2) + 
 (46*Metricdelta[i, j]*ni[a]*ni[b]*dt[2][barMij][-a, k]*dt[3][barMij][-b, -k])/(45*r) + 
 (37*ni[a]*ni[b]*ni[i]*ni[j]*dt[2][barMij][-a, k]*dt[3][barMij][-b, -k])/(15*r) - 
 (3*ni[a]*ni[b]*ni[j]*ni[k]*dt[2][barMij][i, -a]*dt[3][barMij][-b, -k])/(5*r) - 
 (3*ni[a]*ni[b]*ni[i]*ni[k]*dt[2][barMij][j, -a]*dt[3][barMij][-b, -k])/(5*r) + (barMij[j, a]*dt[3][barMij][i, -a])/r^3 + 
 (14*dt[1][barMij][j, a]*dt[3][barMij][i, -a])/(15*r^2) + (179*dt[2][barMij][j, a]*dt[3][barMij][i, -a])/(315*r) - 
 (barMij[j, -a]*ni[a]*ni[b]*dt[3][barMij][i, -b])/r^3 - (5*barMij[-a, b]*ni[a]*ni[j]*dt[3][barMij][i, -b])/r^3 - 
 (31*ni[a]*ni[j]*dt[1][barMij][-a, b]*dt[3][barMij][i, -b])/(5*r^2) + 
 (22*ni[a]*ni[b]*dt[1][barMij][j, -a]*dt[3][barMij][i, -b])/(15*r^2) - 
 (104*ni[a]*ni[j]*dt[2][barMij][-a, b]*dt[3][barMij][i, -b])/(45*r) + 
 (31*ni[a]*ni[b]*dt[2][barMij][j, -a]*dt[3][barMij][i, -b])/(45*r) - (9*barMij[-a, -b]*ni[a]*ni[b]*dt[3][barMij][i, j])/
  (2*r^3) - (113*ni[a]*ni[b]*dt[1][barMij][-a, -b]*dt[3][barMij][i, j])/(30*r^2) + 
 (91*ni[a]*ni[b]*dt[2][barMij][-a, -b]*dt[3][barMij][i, j])/(90*r) + 
 (29*barMij[-a, -b]*ni[a]*ni[b]*ni[j]*ni[k]*dt[3][barMij][i, -k])/(4*r^3) + 
 (71*ni[a]*ni[b]*ni[j]*ni[k]*dt[1][barMij][-a, -b]*dt[3][barMij][i, -k])/(15*r^2) + 
 (9*ni[a]*ni[b]*ni[j]*ni[k]*dt[2][barMij][-a, -b]*dt[3][barMij][i, -k])/(10*r) + (barMij[i, a]*dt[3][barMij][j, -a])/r^3 + 
 (14*dt[1][barMij][i, a]*dt[3][barMij][j, -a])/(15*r^2) + (179*dt[2][barMij][i, a]*dt[3][barMij][j, -a])/(315*r) - 
 (barMij[i, -a]*ni[a]*ni[b]*dt[3][barMij][j, -b])/r^3 - (5*barMij[-a, b]*ni[a]*ni[i]*dt[3][barMij][j, -b])/r^3 - 
 (31*ni[a]*ni[i]*dt[1][barMij][-a, b]*dt[3][barMij][j, -b])/(5*r^2) + 
 (22*ni[a]*ni[b]*dt[1][barMij][i, -a]*dt[3][barMij][j, -b])/(15*r^2) - 
 (104*ni[a]*ni[i]*dt[2][barMij][-a, b]*dt[3][barMij][j, -b])/(45*r) + 
 (31*ni[a]*ni[b]*dt[2][barMij][i, -a]*dt[3][barMij][j, -b])/(45*r) + 
 (29*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[k]*dt[3][barMij][j, -k])/(4*r^3) + 
 (71*ni[a]*ni[b]*ni[i]*ni[k]*dt[1][barMij][-a, -b]*dt[3][barMij][j, -k])/(15*r^2) + 
 (9*ni[a]*ni[b]*ni[i]*ni[k]*dt[2][barMij][-a, -b]*dt[3][barMij][j, -k])/(10*r) - 
 (13*barMij[-a, -b]*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*dt[3][barMij][-k, -l])/(2*r^3) - 
 (7*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*dt[3][barMij][-k, -l])/(2*r^3) - 
 (173*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[3][barMij][-k, -l])/(60*r^2) - 
 (49*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[3][barMij][-k, -l])/(20*r^2) - 
 (3*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*dt[2][barMij][-a, -b]*dt[3][barMij][-k, -l])/(10*r) - 
 (37*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*dt[2][barMij][-a, -b]*dt[3][barMij][-k, -l])/(60*r) - 
 (3*barMij[a, b]*Metricdelta[i, j]*dt[4][barMij][-a, -b])/(5*r^2) + (barMij[i, j]*ni[a]*ni[b]*dt[4][barMij][-a, -b])/
  (30*r^2) - (barMij[j, b]*ni[a]*ni[i]*dt[4][barMij][-a, -b])/(10*r^2) - (barMij[i, b]*ni[a]*ni[j]*dt[4][barMij][-a, -b])/
  (10*r^2) - (2*barMij[a, b]*ni[i]*ni[j]*dt[4][barMij][-a, -b])/(15*r^2) + 
 (452*Metricdelta[i, j]*dt[1][barMij][a, b]*dt[4][barMij][-a, -b])/(315*r) - 
 (31*ni[i]*ni[j]*dt[1][barMij][a, b]*dt[4][barMij][-a, -b])/(90*r) - (ni[a]*ni[j]*dt[1][barMij][i, b]*dt[4][barMij][-a, -b])/
  (45*r) + (29*ni[a]*ni[b]*dt[1][barMij][i, j]*dt[4][barMij][-a, -b])/(90*r) - 
 (ni[a]*ni[i]*dt[1][barMij][j, b]*dt[4][barMij][-a, -b])/(45*r) + 
 (12*barMij[-a, k]*Metricdelta[i, j]*ni[a]*ni[b]*dt[4][barMij][-b, -k])/(5*r^2) + 
 (11*barMij[-a, k]*ni[a]*ni[b]*ni[i]*ni[j]*dt[4][barMij][-b, -k])/(15*r^2) + 
 (barMij[j, -a]*ni[a]*ni[b]*ni[i]*ni[k]*dt[4][barMij][-b, -k])/(30*r^2) + 
 (barMij[i, -a]*ni[a]*ni[b]*ni[j]*ni[k]*dt[4][barMij][-b, -k])/(30*r^2) + 
 (74*Metricdelta[i, j]*ni[a]*ni[b]*dt[1][barMij][-a, k]*dt[4][barMij][-b, -k])/(45*r) + 
 (8*ni[a]*ni[b]*ni[i]*ni[j]*dt[1][barMij][-a, k]*dt[4][barMij][-b, -k])/(15*r) - 
 (3*ni[a]*ni[b]*ni[j]*ni[k]*dt[1][barMij][i, -a]*dt[4][barMij][-b, -k])/(20*r) - 
 (3*ni[a]*ni[b]*ni[i]*ni[k]*dt[1][barMij][j, -a]*dt[4][barMij][-b, -k])/(20*r) + 
 (19*barMij[j, a]*dt[4][barMij][i, -a])/(30*r^2) + (106*dt[1][barMij][j, a]*dt[4][barMij][i, -a])/(315*r) + 
 (17*barMij[j, -a]*ni[a]*ni[b]*dt[4][barMij][i, -b])/(30*r^2) - (31*barMij[-a, b]*ni[a]*ni[j]*dt[4][barMij][i, -b])/(10*r^2) - 
 (227*ni[a]*ni[j]*dt[1][barMij][-a, b]*dt[4][barMij][i, -b])/(90*r) + 
 (44*ni[a]*ni[b]*dt[1][barMij][j, -a]*dt[4][barMij][i, -b])/(45*r) + (barMij[-a, -b]*ni[a]*ni[b]*dt[4][barMij][i, j])/
  (30*r^2) - (8*ni[a]*ni[b]*dt[1][barMij][-a, -b]*dt[4][barMij][i, j])/(45*r) + 
 (23*barMij[-a, -b]*ni[a]*ni[b]*ni[j]*ni[k]*dt[4][barMij][i, -k])/(15*r^2) + 
 (3*ni[a]*ni[b]*ni[j]*ni[k]*dt[1][barMij][-a, -b]*dt[4][barMij][i, -k])/(5*r) + 
 (19*barMij[i, a]*dt[4][barMij][j, -a])/(30*r^2) + (106*dt[1][barMij][i, a]*dt[4][barMij][j, -a])/(315*r) + 
 (17*barMij[i, -a]*ni[a]*ni[b]*dt[4][barMij][j, -b])/(30*r^2) - (31*barMij[-a, b]*ni[a]*ni[i]*dt[4][barMij][j, -b])/(10*r^2) - 
 (227*ni[a]*ni[i]*dt[1][barMij][-a, b]*dt[4][barMij][j, -b])/(90*r) + 
 (44*ni[a]*ni[b]*dt[1][barMij][i, -a]*dt[4][barMij][j, -b])/(45*r) + 
 (23*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[k]*dt[4][barMij][j, -k])/(15*r^2) + 
 (3*ni[a]*ni[b]*ni[i]*ni[k]*dt[1][barMij][-a, -b]*dt[4][barMij][j, -k])/(5*r) - 
 (89*barMij[-a, -b]*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*dt[4][barMij][-k, -l])/(60*r^2) - 
 (7*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*dt[4][barMij][-k, -l])/(20*r^2) - 
 (9*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[4][barMij][-k, -l])/(20*r) - 
 (2*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*dt[1][barMij][-a, -b]*dt[4][barMij][-k, -l])/(15*r) + 
 (163*barMij[a, b]*Metricdelta[i, j]*dt[5][barMij][-a, -b])/(315*r) + (barMij[i, j]*ni[a]*ni[b]*dt[5][barMij][-a, -b])/
  (90*r) - (barMij[j, b]*ni[a]*ni[i]*dt[5][barMij][-a, -b])/(90*r) - (barMij[i, b]*ni[a]*ni[j]*dt[5][barMij][-a, -b])/(90*r) - 
 (barMij[a, b]*ni[i]*ni[j]*dt[5][barMij][-a, -b])/(45*r) + 
 (28*barMij[-a, k]*Metricdelta[i, j]*ni[a]*ni[b]*dt[5][barMij][-b, -k])/(45*r) + 
 (barMij[-a, k]*ni[a]*ni[b]*ni[i]*ni[j]*dt[5][barMij][-b, -k])/(15*r) - (73*barMij[j, a]*dt[5][barMij][i, -a])/(315*r) + 
 (13*barMij[j, -a]*ni[a]*ni[b]*dt[5][barMij][i, -b])/(45*r) - (41*barMij[-a, b]*ni[a]*ni[j]*dt[5][barMij][i, -b])/(45*r) + 
 (barMij[-a, -b]*ni[a]*ni[b]*dt[5][barMij][i, j])/(9*r) + (3*barMij[-a, -b]*ni[a]*ni[b]*ni[j]*ni[k]*dt[5][barMij][i, -k])/
  (20*r) - (73*barMij[i, a]*dt[5][barMij][j, -a])/(315*r) + (13*barMij[i, -a]*ni[a]*ni[b]*dt[5][barMij][j, -b])/(45*r) - 
 (41*barMij[-a, b]*ni[a]*ni[i]*dt[5][barMij][j, -b])/(45*r) + (3*barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[k]*dt[5][barMij][j, -k])/
  (20*r) - (3*barMij[-a, -b]*Metricdelta[i, j]*ni[a]*ni[b]*ni[k]*ni[l]*dt[5][barMij][-k, -l])/(20*r) - 
 (barMij[-a, -b]*ni[a]*ni[b]*ni[i]*ni[j]*ni[k]*ni[l]*dt[5][barMij][-k, -l])/(60*r)
