(* ::Package:: *)

(* ::Text:: *)
(*Ancillary file to the article "Gravitational Wave Flux and Quadrupole Modes from Quasi-Circular Non-Spinning Compact Binaries to the Fourth Post-Newtonian Order" and the letter "Gravitational-Wave Phasing of Compact Binary Systems to the Fourth-and-a-Half post-Newtonian Order"*)
(**)
(*Authors: L. Blanchet, G. Faye, Q. Henry, F. Larrouturou and D. Trestini*)
(*Type: Mathematica/xTensor program*)
(**)
(*This file contains the results associated with the gravitational-wave phase at 4.5PN and the (2,2) mode at 4PN, presented in the  paper arXiv:2304.11186 and the letter  arXiv:2304.11185*)
(**)


(*
Content :
FluxCanLoc4PN : local sector of the canonical flux on generic orbits, at 4PN, see eq. (6.3) of the paper.
FluxCirc45PN : gravitational flux on quasi-circular orbits, at 4.5PN, see eq. (6.11) of the paper or eq. (4) of the letter.
Freq4N : chirp frequency x in terms of the rescaled time, at 4.5PN, see eq. (6) of the letter.
Phase45PN : observable gravitational phase, at 4.5PN, see eq. (8) of the letter.
PhaseSPA45PN : gravitational phase in the stationary phase approximation, at 4.5PN, see eq. (9) of the letter.
H224PN : dominant quadrupole mode, at 4PN, see eq. (6.16) of the paper or eq. (11) of the letter.
*)


(*
Notations:
m = m1 + m2 is the total mass.
\[Nu] = m1 m2 / m^2 is the symmetric mass ratio.
xi[i] = r ni[i] is the position in the center-of-mass frame.
vi[i] is the velocity in the center-of-mass frame.
r0 is the scale associated with Hadamard regularization.
r0p is the scale associated with UV dimensional regularization, see footnote 10 of arXiv:2003.13672.

x = (\[Pi]Gmf/c^3)^(3/2) is defined in eq. (6.9) of the paper or eq. (1) of the letter.
\[Tau] is the rescaled time variable, see eq. (5) of the letter.
v = x^(1/2) is the SPA variable.
\[Psi]0 is the arbitrary constant associated to Phase45PN, see eq. (8) of the letter.
T0 and \[CapitalPsi]0 are the two arbitrary constants associated to PhaseSPA45PN, see eq. (9) of the letter.
*)


FluxCanLoc4PN = (-2683093592*G^8*m^9*\[Nu]^2)/(19864845*c^13*r^9) + 
     (300568*G^7*m^8*\[Nu]^2)/(10395*c^11*r^8) - (4048*G^6*m^7*\[Nu]^2)/
      (945*c^9*r^7) + (32*G^5*m^6*\[Nu]^2)/(105*c^7*r^6) + 
     (122461316*G^8*m^9*\[Nu]^3)/(361179*c^13*r^9) + 
     (82*G^8*m^9*Pi^2*\[Nu]^3)/(105*c^13*r^9) - (478784*G^7*m^8*\[Nu]^3)/
      (4455*c^11*r^8) + (608*G^6*m^7*\[Nu]^3)/(35*c^9*r^7) - 
     (128*G^5*m^6*\[Nu]^3)/(105*c^7*r^6) + (461586340*G^8*m^9*\[Nu]^4)/
      (567567*c^13*r^9) - (328*G^8*m^9*Pi^2*\[Nu]^4)/(105*c^13*r^9) - 
     (48304*G^7*m^8*\[Nu]^4)/(1485*c^11*r^8) - (128*G^6*m^7*\[Nu]^4)/
      (135*c^9*r^7) - (13688176*G^8*m^9*\[Nu]^5)/(405405*c^13*r^9) - 
     (27712*G^7*m^8*\[Nu]^5)/(31185*c^11*r^8) + (54272*G^8*m^9*\[Nu]^6)/
      (405405*c^13*r^9) + (3328*G^8*m^9*\[Nu]^2*Log[r/r0])/(2205*c^13*r^9) - 
     (13312*G^8*m^9*\[Nu]^3*Log[r/r0])/(2205*c^13*r^9) + 
     (1408*G^8*m^9*\[Nu]^3*Log[r/r0p])/(105*c^13*r^9) - 
     (5632*G^8*m^9*\[Nu]^4*Log[r/r0p])/(105*c^13*r^9) - 
     (4619968*G^7*m^8*\[Nu]^3*Scalar[ni[a]*vi[-a]])/(99225*c^12*r^8) + 
     (192*G^6*m^7*\[Nu]^3*Scalar[ni[a]*vi[-a]])/(175*c^10*r^7) + 
     (2074144*G^7*m^8*\[Nu]^4*Scalar[ni[a]*vi[-a]])/(99225*c^12*r^8) - 
     (350787478967*G^7*m^8*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2)/
      (33108075*c^13*r^8) + (25620808*G^6*m^7*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2)/
      (28875*c^11*r^7) - (212638*G^5*m^6*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2)/
      (945*c^9*r^6) + (5872*G^4*m^5*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2)/
      (105*c^7*r^5) - (88*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2)/
      (15*c^5*r^4) - (8223972079399*G^7*m^8*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2)/
      (297972675*c^13*r^8) + (7113791*G^7*m^8*Pi^2*\[Nu]^3*
       Scalar[ni[a]*vi[-a]]^2)/(5040*c^13*r^8) + 
     (24880796*G^6*m^7*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2)/(14175*c^11*r^7) - 
     (1763*G^6*m^7*Pi^2*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2)/(30*c^11*r^7) - 
     (6532*G^5*m^6*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2)/(315*c^9*r^6) - 
     (16*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2)/(7*c^7*r^5) - 
     (62489606891*G^7*m^8*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2)/
      (28378350*c^13*r^8) + (14323*G^7*m^8*Pi^2*\[Nu]^4*
       Scalar[ni[a]*vi[-a]]^2)/(105*c^13*r^8) + 
     (230308*G^6*m^7*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2)/(1155*c^11*r^7) - 
     (512*G^5*m^6*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2)/(45*c^9*r^6) + 
     (161364376*G^7*m^8*\[Nu]^5*Scalar[ni[a]*vi[-a]]^2)/(405405*c^13*r^8) - 
     (1376576*G^6*m^7*\[Nu]^5*Scalar[ni[a]*vi[-a]]^2)/(31185*c^11*r^7) - 
     (55496729*G^7*m^8*\[Nu]^6*Scalar[ni[a]*vi[-a]]^2)/(405405*c^13*r^8) + 
     (13872832*G^7*m^8*\[Nu]^2*Log[r/r0]*Scalar[ni[a]*vi[-a]]^2)/
      (11025*c^13*r^8) - (54784*G^6*m^7*\[Nu]^2*Log[r/r0]*
       Scalar[ni[a]*vi[-a]]^2)/(1575*c^11*r^7) - 
     (653824*G^7*m^8*\[Nu]^3*Log[r/r0]*Scalar[ni[a]*vi[-a]]^2)/
      (3675*c^13*r^8) + (52352*G^7*m^8*\[Nu]^3*Log[r/r0p]*
       Scalar[ni[a]*vi[-a]]^2)/(21*c^13*r^8) - 
     (2816*G^6*m^7*\[Nu]^3*Log[r/r0p]*Scalar[ni[a]*vi[-a]]^2)/(45*c^11*r^7) - 
     (216064*G^7*m^8*\[Nu]^4*Log[r/r0p]*Scalar[ni[a]*vi[-a]]^2)/
      (315*c^13*r^8) + (2582336*G^6*m^7*\[Nu]^3*Scalar[ni[a]*vi[-a]]^3)/
      (99225*c^12*r^7) + (182816*G^5*m^6*\[Nu]^3*Scalar[ni[a]*vi[-a]]^3)/
      (1575*c^10*r^6) + (14233216*G^6*m^7*\[Nu]^4*Scalar[ni[a]*vi[-a]]^3)/
      (99225*c^12*r^7) - (278448915058*G^6*m^7*\[Nu]^2*
       Scalar[ni[a]*vi[-a]]^4)/(4729725*c^13*r^7) + 
     (231255634*G^5*m^6*\[Nu]^2*Scalar[ni[a]*vi[-a]]^4)/(51975*c^11*r^6) - 
     (17872*G^4*m^5*\[Nu]^2*Scalar[ni[a]*vi[-a]]^4)/(63*c^9*r^5) + 
     (1374*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^4)/(35*c^7*r^4) - 
     (472247121109*G^6*m^7*\[Nu]^3*Scalar[ni[a]*vi[-a]]^4)/
      (11036025*c^13*r^7) + (940127*G^6*m^7*Pi^2*\[Nu]^3*
       Scalar[ni[a]*vi[-a]]^4)/(240*c^13*r^7) + 
     (170684*G^5*m^6*\[Nu]^3*Scalar[ni[a]*vi[-a]]^4)/(495*c^11*r^6) - 
     (697*G^5*m^6*Pi^2*\[Nu]^3*Scalar[ni[a]*vi[-a]]^4)/(5*c^11*r^6) + 
     (487768*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^4)/(945*c^9*r^5) - 
     (248*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^4)/(7*c^7*r^4) + 
     (5344151986*G^6*m^7*\[Nu]^4*Scalar[ni[a]*vi[-a]]^4)/(525525*c^13*r^7) - 
     (12013*G^6*m^7*Pi^2*\[Nu]^4*Scalar[ni[a]*vi[-a]]^4)/(28*c^13*r^7) + 
     (439742*G^5*m^6*\[Nu]^4*Scalar[ni[a]*vi[-a]]^4)/(297*c^11*r^6) - 
     (22864*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]^4)/(189*c^9*r^5) + 
     (438245876*G^6*m^7*\[Nu]^5*Scalar[ni[a]*vi[-a]]^4)/(135135*c^13*r^7) - 
     (3301324*G^5*m^6*\[Nu]^5*Scalar[ni[a]*vi[-a]]^4)/(10395*c^11*r^6) - 
     (33411722*G^6*m^7*\[Nu]^6*Scalar[ni[a]*vi[-a]]^4)/(45045*c^13*r^7) + 
     (39072*G^6*m^7*\[Nu]^2*Log[r/r0]*Scalar[ni[a]*vi[-a]]^4)/(5*c^13*r^7) - 
     (6848*G^5*m^6*\[Nu]^2*Log[r/r0]*Scalar[ni[a]*vi[-a]]^4)/(15*c^11*r^6) + 
     (779456*G^6*m^7*\[Nu]^3*Log[r/r0]*Scalar[ni[a]*vi[-a]]^4)/
      (1225*c^13*r^7) + (54208*G^6*m^7*\[Nu]^3*Log[r/r0p]*
       Scalar[ni[a]*vi[-a]]^4)/(15*c^13*r^7) + 
     (3872*G^6*m^7*\[Nu]^4*Log[r/r0p]*Scalar[ni[a]*vi[-a]]^4)/
      (105*c^13*r^7) + (8021704*G^5*m^6*\[Nu]^3*Scalar[ni[a]*vi[-a]]^5)/
      (735*c^12*r^6) - (17312*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^5)/
      (35*c^10*r^5) - (6133488*G^5*m^6*\[Nu]^4*Scalar[ni[a]*vi[-a]]^5)/
      (1225*c^12*r^6) - (702363704057*G^5*m^6*\[Nu]^2*Scalar[ni[a]*vi[-a]]^6)/
      (22072050*c^13*r^6) + (10953902*G^4*m^5*\[Nu]^2*Scalar[ni[a]*vi[-a]]^6)/
      (17325*c^11*r^5) - (10004*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^6)/
      (315*c^9*r^4) + (571186117207*G^5*m^6*\[Nu]^3*Scalar[ni[a]*vi[-a]]^6)/
      (8108100*c^13*r^6) + (238723*G^5*m^6*Pi^2*\[Nu]^3*
       Scalar[ni[a]*vi[-a]]^6)/(448*c^13*r^6) - 
     (268706*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^6)/(77*c^11*r^5) + 
     (80936*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^6)/(315*c^9*r^4) - 
     (511326653*G^5*m^6*\[Nu]^4*Scalar[ni[a]*vi[-a]]^6)/(90090*c^13*r^6) - 
     (6109*G^5*m^6*Pi^2*\[Nu]^4*Scalar[ni[a]*vi[-a]]^6)/(7*c^13*r^6) + 
     (10410038*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]^6)/(3465*c^11*r^5) - 
     (33616*G^3*m^4*\[Nu]^4*Scalar[ni[a]*vi[-a]]^6)/(315*c^9*r^4) + 
     (370344998*G^5*m^6*\[Nu]^5*Scalar[ni[a]*vi[-a]]^6)/(27027*c^13*r^6) - 
     (1720954*G^4*m^5*\[Nu]^5*Scalar[ni[a]*vi[-a]]^6)/(3465*c^11*r^5) - 
     (8735219*G^5*m^6*\[Nu]^6*Scalar[ni[a]*vi[-a]]^6)/(8190*c^13*r^6) + 
     (192240*G^5*m^6*\[Nu]^2*Log[r/r0]*Scalar[ni[a]*vi[-a]]^6)/
      (49*c^13*r^6) + (1152*G^5*m^6*\[Nu]^3*Log[r/r0]*Scalar[ni[a]*vi[-a]]^6)/
      (c^13*r^6) + (37622824*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^7)/
      (6615*c^12*r^5) - (93205696*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]^7)/
      (6615*c^12*r^5) - (19607293*G^4*m^5*\[Nu]^2*Scalar[ni[a]*vi[-a]]^8)/
      (42042*c^13*r^5) + (301585*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^8)/
      (1386*c^11*r^4) + (10047409627*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^8)/
      (1051050*c^13*r^5) - (65168*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^8)/
      (63*c^11*r^4) - (192617639*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]^8)/
      (6930*c^13*r^5) + (549844*G^3*m^4*\[Nu]^4*Scalar[ni[a]*vi[-a]]^8)/
      (693*c^11*r^4) + (425266136*G^4*m^5*\[Nu]^5*Scalar[ni[a]*vi[-a]]^8)/
      (45045*c^13*r^5) - (104816*G^3*m^4*\[Nu]^5*Scalar[ni[a]*vi[-a]]^8)/
      (693*c^11*r^4) + (4982053*G^4*m^5*\[Nu]^6*Scalar[ni[a]*vi[-a]]^8)/
      (30030*c^13*r^5) - (103711*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^10)/
      (385*c^13*r^4) + (8680813*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^10)/
      (3003*c^13*r^4) - (23488204*G^3*m^4*\[Nu]^4*Scalar[ni[a]*vi[-a]]^10)/
      (3003*c^13*r^4) + (5127872*G^3*m^4*\[Nu]^5*Scalar[ni[a]*vi[-a]]^10)/
      (3003*c^13*r^4) + (398872*G^3*m^4*\[Nu]^6*Scalar[ni[a]*vi[-a]]^10)/
      (3003*c^13*r^4) + (5232680272544*G^7*m^8*\[Nu]^2*Scalar[vi[-a]*vi[a]])/
      (496621125*c^13*r^8) - (787471744*G^6*m^7*\[Nu]^2*Scalar[vi[-a]*vi[a]])/
      (779625*c^11*r^7) + (562946*G^5*m^6*\[Nu]^2*Scalar[vi[-a]*vi[a]])/
      (2835*c^9*r^6) - (1088*G^4*m^5*\[Nu]^2*Scalar[vi[-a]*vi[a]])/
      (21*c^7*r^5) + (32*G^3*m^4*\[Nu]^2*Scalar[vi[-a]*vi[a]])/(5*c^5*r^4) + 
     (11310753492221*G^7*m^8*\[Nu]^3*Scalar[vi[-a]*vi[a]])/
      (496621125*c^13*r^8) - (834503*G^7*m^8*Pi^2*\[Nu]^3*
       Scalar[vi[-a]*vi[a]])/(840*c^13*r^8) - 
     (25425164*G^6*m^7*\[Nu]^3*Scalar[vi[-a]*vi[a]])/(14175*c^11*r^7) + 
     (451*G^6*m^7*Pi^2*\[Nu]^3*Scalar[vi[-a]*vi[a]])/(10*c^11*r^7) + 
     (18184*G^5*m^6*\[Nu]^3*Scalar[vi[-a]*vi[a]])/(315*c^9*r^6) + 
     (64*G^4*m^5*\[Nu]^3*Scalar[vi[-a]*vi[a]])/(21*c^7*r^5) + 
     (26171918846*G^7*m^8*\[Nu]^4*Scalar[vi[-a]*vi[a]])/(14189175*c^13*r^8) - 
     (21836*G^7*m^8*Pi^2*\[Nu]^4*Scalar[vi[-a]*vi[a]])/(105*c^13*r^8) + 
     (11920*G^6*m^7*\[Nu]^4*Scalar[vi[-a]*vi[a]])/(231*c^11*r^7) + 
     (416*G^5*m^6*\[Nu]^4*Scalar[vi[-a]*vi[a]])/(135*c^9*r^6) + 
     (18124*G^7*m^8*\[Nu]^5*Scalar[vi[-a]*vi[a]])/(99*c^13*r^8) - 
     (6728*G^6*m^7*\[Nu]^5*Scalar[vi[-a]*vi[a]])/(2835*c^11*r^7) - 
     (228626*G^7*m^8*\[Nu]^6*Scalar[vi[-a]*vi[a]])/(10395*c^13*r^8) - 
     (2659712*G^7*m^8*\[Nu]^2*Log[r/r0]*Scalar[vi[-a]*vi[a]])/
      (2205*c^13*r^8) + (27392*G^6*m^7*\[Nu]^2*Log[r/r0]*
       Scalar[vi[-a]*vi[a]])/(525*c^11*r^7) - 
     (279808*G^7*m^8*\[Nu]^3*Log[r/r0]*Scalar[vi[-a]*vi[a]])/
      (11025*c^13*r^8) - (289024*G^7*m^8*\[Nu]^3*Log[r/r0p]*
       Scalar[vi[-a]*vi[a]])/(105*c^13*r^8) + 
     (1408*G^6*m^7*\[Nu]^3*Log[r/r0p]*Scalar[vi[-a]*vi[a]])/(15*c^11*r^7) + 
     (49536*G^7*m^8*\[Nu]^4*Log[r/r0p]*Scalar[vi[-a]*vi[a]])/(35*c^13*r^8) - 
     (7338848*G^6*m^7*\[Nu]^3*Scalar[ni[a]*vi[-a]]*Scalar[vi[-a]*vi[a]])/
      (33075*c^12*r^7) - (4832*G^5*m^6*\[Nu]^3*Scalar[ni[a]*vi[-a]]*
       Scalar[vi[-a]*vi[a]])/(45*c^10*r^6) - 
     (721696*G^6*m^7*\[Nu]^4*Scalar[ni[a]*vi[-a]]*Scalar[vi[-a]*vi[a]])/
      (33075*c^12*r^7) + (4038478621136*G^6*m^7*\[Nu]^2*
       Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]])/(55180125*c^13*r^7) - 
     (142516036*G^5*m^6*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]])/
      (23625*c^11*r^6) + (39896*G^4*m^5*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(105*c^9*r^5) - 
     (5948*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]])/
      (105*c^7*r^4) + (7745487110642*G^6*m^7*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(165540375*c^13*r^7) - 
     (2800589*G^6*m^7*Pi^2*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(630*c^13*r^7) - 
     (301478*G^5*m^6*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]])/
      (385*c^11*r^6) + (861*G^5*m^6*Pi^2*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(5*c^11*r^6) - 
     (68104*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]])/
      (105*c^9*r^5) + (1856*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(35*c^7*r^4) - 
     (191341538624*G^6*m^7*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(14189175*c^13*r^7) + 
     (26937*G^6*m^7*Pi^2*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]])/
      (70*c^13*r^7) - (1812988*G^5*m^6*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(1155*c^11*r^6) + 
     (3464*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]])/
      (21*c^9*r^5) - (21106136*G^6*m^7*\[Nu]^5*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(12285*c^13*r^7) + 
     (3968648*G^5*m^6*\[Nu]^5*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]])/
      (10395*c^11*r^6) + (94375252*G^6*m^7*\[Nu]^6*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(135135*c^13*r^7) - 
     (107581888*G^6*m^7*\[Nu]^2*Log[r/r0]*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(11025*c^13*r^7) + 
     (109568*G^5*m^6*\[Nu]^2*Log[r/r0]*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(175*c^11*r^6) - 
     (613952*G^6*m^7*\[Nu]^3*Log[r/r0]*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(441*c^13*r^7) - 
     (1476992*G^6*m^7*\[Nu]^3*Log[r/r0p]*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(315*c^13*r^7) - 
     (8096*G^6*m^7*\[Nu]^4*Log[r/r0p]*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]])/(315*c^13*r^7) - 
     (203720096*G^5*m^6*\[Nu]^3*Scalar[ni[a]*vi[-a]]^3*Scalar[vi[-a]*vi[a]])/
      (11025*c^12*r^6) + (459584*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^3*
       Scalar[vi[-a]*vi[a]])/(525*c^10*r^5) + 
     (20621536*G^5*m^6*\[Nu]^4*Scalar[ni[a]*vi[-a]]^3*Scalar[vi[-a]*vi[a]])/
      (2205*c^12*r^6) + (841346328572*G^5*m^6*\[Nu]^2*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]])/(15049125*c^13*r^6) - 
     (1430314*G^4*m^5*\[Nu]^2*Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]])/
      (1155*c^11*r^5) + (8072*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]])/(105*c^9*r^4) - 
     (133918092919027*G^5*m^6*\[Nu]^3*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]])/(993242250*c^13*r^6) - 
     (918441*G^5*m^6*Pi^2*\[Nu]^3*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]])/(1120*c^13*r^6) + 
     (70316074*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]])/
      (10395*c^11*r^5) - (20276*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]])/(35*c^9*r^4) + 
     (36645503*G^5*m^6*\[Nu]^4*Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]])/
      (3861*c^13*r^6) + (11029*G^5*m^6*Pi^2*\[Nu]^4*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]])/(7*c^13*r^6) - 
     (7344286*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]])/
      (1155*c^11*r^5) + (10096*G^3*m^4*\[Nu]^4*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]])/(35*c^9*r^4) - 
     (346417333*G^5*m^6*\[Nu]^5*Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]])/
      (12285*c^13*r^6) + (2787280*G^4*m^5*\[Nu]^5*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]])/(2079*c^11*r^5) + 
     (12343106*G^5*m^6*\[Nu]^6*Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]])/
      (3861*c^13*r^6) - (8632816*G^5*m^6*\[Nu]^2*Log[r/r0]*
       Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]])/(1225*c^13*r^6) - 
     (341632*G^5*m^6*\[Nu]^3*Log[r/r0]*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]])/(147*c^13*r^6) - 
     (144416728*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^5*Scalar[vi[-a]*vi[a]])/
      (11025*c^12*r^5) + (384856048*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]^5*
       Scalar[vi[-a]*vi[a]])/(11025*c^12*r^5) + 
     (208616293*G^4*m^5*\[Nu]^2*Scalar[ni[a]*vi[-a]]^6*Scalar[vi[-a]*vi[a]])/
      (150150*c^13*r^5) - (2011958*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^6*
       Scalar[vi[-a]*vi[a]])/(3465*c^11*r^4) - 
     (1085899901*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^6*Scalar[vi[-a]*vi[a]])/
      (40950*c^13*r^5) + (10358396*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^6*
       Scalar[vi[-a]*vi[a]])/(3465*c^11*r^4) + 
     (1434092515*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]^6*Scalar[vi[-a]*vi[a]])/
      (18018*c^13*r^5) - (10577128*G^3*m^4*\[Nu]^4*Scalar[ni[a]*vi[-a]]^6*
       Scalar[vi[-a]*vi[a]])/(3465*c^11*r^4) - 
     (1632275384*G^4*m^5*\[Nu]^5*Scalar[ni[a]*vi[-a]]^6*Scalar[vi[-a]*vi[a]])/
      (45045*c^13*r^5) + (578912*G^3*m^4*\[Nu]^5*Scalar[ni[a]*vi[-a]]^6*
       Scalar[vi[-a]*vi[a]])/(693*c^11*r^4) + 
     (16632257*G^4*m^5*\[Nu]^6*Scalar[ni[a]*vi[-a]]^6*Scalar[vi[-a]*vi[a]])/
      (6930*c^13*r^5) + (10059068*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^8*
       Scalar[vi[-a]*vi[a]])/(9009*c^13*r^4) - 
     (103349531*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^8*Scalar[vi[-a]*vi[a]])/
      (9009*c^13*r^4) + (280637620*G^3*m^4*\[Nu]^4*Scalar[ni[a]*vi[-a]]^8*
       Scalar[vi[-a]*vi[a]])/(9009*c^13*r^4) - 
     (95166598*G^3*m^4*\[Nu]^5*Scalar[ni[a]*vi[-a]]^8*Scalar[vi[-a]*vi[a]])/
      (9009*c^13*r^4) + (5148040*G^3*m^4*\[Nu]^6*Scalar[ni[a]*vi[-a]]^8*
       Scalar[vi[-a]*vi[a]])/(9009*c^13*r^4) - 
     (264916580144*G^6*m^7*\[Nu]^2*Scalar[vi[-a]*vi[a]]^2)/
      (18393375*c^13*r^7) + (1197229882*G^5*m^6*\[Nu]^2*
       Scalar[vi[-a]*vi[a]]^2)/(779625*c^11*r^6) - 
     (3952*G^4*m^5*\[Nu]^2*Scalar[vi[-a]*vi[a]]^2)/(35*c^9*r^5) + 
     (314*G^3*m^4*\[Nu]^2*Scalar[vi[-a]*vi[a]]^2)/(21*c^7*r^4) - 
     (4013494967431*G^6*m^7*\[Nu]^3*Scalar[vi[-a]*vi[a]]^2)/
      (496621125*c^13*r^7) + (205997*G^6*m^7*Pi^2*\[Nu]^3*
       Scalar[vi[-a]*vi[a]]^2)/(240*c^13*r^7) + 
     (1276672*G^5*m^6*\[Nu]^3*Scalar[vi[-a]*vi[a]]^2)/(10395*c^11*r^6) - 
     (369*G^5*m^6*Pi^2*\[Nu]^3*Scalar[vi[-a]*vi[a]]^2)/(10*c^11*r^6) + 
     (41896*G^4*m^5*\[Nu]^3*Scalar[vi[-a]*vi[a]]^2)/(315*c^9*r^5) - 
     (568*G^3*m^4*\[Nu]^3*Scalar[vi[-a]*vi[a]]^2)/(35*c^7*r^4) + 
     (223158212*G^6*m^7*\[Nu]^4*Scalar[vi[-a]*vi[a]]^2)/(155925*c^13*r^7) - 
     (1517*G^6*m^7*Pi^2*\[Nu]^4*Scalar[vi[-a]*vi[a]]^2)/(28*c^13*r^7) + 
     (2601814*G^5*m^6*\[Nu]^4*Scalar[vi[-a]*vi[a]]^2)/(10395*c^11*r^6) - 
     (1592*G^4*m^5*\[Nu]^4*Scalar[vi[-a]*vi[a]]^2)/(45*c^9*r^5) + 
     (284332*G^6*m^7*\[Nu]^5*Scalar[vi[-a]*vi[a]]^2)/(36855*c^13*r^7) - 
     (1294264*G^5*m^6*\[Nu]^5*Scalar[vi[-a]*vi[a]]^2)/(31185*c^11*r^6) + 
     (1133260*G^6*m^7*\[Nu]^6*Scalar[vi[-a]*vi[a]]^2)/(81081*c^13*r^7) + 
     (7004192*G^6*m^7*\[Nu]^2*Log[r/r0]*Scalar[vi[-a]*vi[a]]^2)/
      (3675*c^13*r^7) - (27392*G^5*m^6*\[Nu]^2*Log[r/r0]*
       Scalar[vi[-a]*vi[a]]^2)/(175*c^11*r^6) + 
     (7309504*G^6*m^7*\[Nu]^3*Log[r/r0]*Scalar[vi[-a]*vi[a]]^2)/
      (11025*c^13*r^7) + (11264*G^6*m^7*\[Nu]^3*Log[r/r0p]*
       Scalar[vi[-a]*vi[a]]^2)/(15*c^13*r^7) + 
     (2816*G^6*m^7*\[Nu]^4*Log[r/r0p]*Scalar[vi[-a]*vi[a]]^2)/(15*c^13*r^7) + 
     (247407416*G^5*m^6*\[Nu]^3*Scalar[ni[a]*vi[-a]]*Scalar[vi[-a]*vi[a]]^2)/
      (33075*c^12*r^6) - (595232*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]*
       Scalar[vi[-a]*vi[a]]^2)/(1575*c^10*r^5) - 
     (140406992*G^5*m^6*\[Nu]^4*Scalar[ni[a]*vi[-a]]*Scalar[vi[-a]*vi[a]]^2)/
      (33075*c^12*r^6) - (330860601266*G^5*m^6*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^2)/(12733875*c^13*r^6) + 
     (54678*G^4*m^5*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]]^2)/
      (77*c^11*r^5) - (2292*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^2)/(35*c^9*r^4) + 
     (4410089814139*G^5*m^6*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^2)/(60196500*c^13*r^6) + 
     (108201*G^5*m^6*Pi^2*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^2)/(320*c^13*r^6) - 
     (405386*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]]^2)/
      (105*c^11*r^5) + (13704*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^2)/(35*c^9*r^4) - 
     (165373057*G^5*m^6*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^2)/(32175*c^13*r^6) - 
     (58917*G^5*m^6*Pi^2*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^2)/(70*c^13*r^6) + 
     (1488754*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]]^2)/
      (385*c^11*r^5) - (25168*G^3*m^4*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^2)/(105*c^9*r^4) + 
     (2249568533*G^5*m^6*\[Nu]^5*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^2)/(135135*c^13*r^6) - 
     (3724396*G^4*m^5*\[Nu]^5*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]]^2)/
      (3465*c^11*r^5) - (50572768*G^5*m^6*\[Nu]^6*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^2)/(19305*c^13*r^6) + 
     (4134416*G^5*m^6*\[Nu]^2*Log[r/r0]*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^2)/(1225*c^13*r^6) + 
     (1448512*G^5*m^6*\[Nu]^3*Log[r/r0]*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^2)/(1225*c^13*r^6) + 
     (105119848*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^3*
       Scalar[vi[-a]*vi[a]]^2)/(11025*c^12*r^5) - 
     (307952032*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]^3*
       Scalar[vi[-a]*vi[a]]^2)/(11025*c^12*r^5) - 
     (68948249*G^4*m^5*\[Nu]^2*Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]]^2)/
      (45045*c^13*r^5) + (204349*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]]^2)/(385*c^11*r^4) + 
     (3401965282*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]]^2)/(135135*c^13*r^5) - 
     (3522149*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]]^2)/
      (1155*c^11*r^4) - (10543843493*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]]^2)/(135135*c^13*r^5) + 
     (4709506*G^3*m^4*\[Nu]^4*Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]]^2)/
      (1155*c^11*r^4) + (6319709261*G^4*m^5*\[Nu]^5*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]]^2)/(135135*c^13*r^5) - 
     (1751152*G^3*m^4*\[Nu]^5*Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]]^2)/
      (1155*c^11*r^4) - (983746018*G^4*m^5*\[Nu]^6*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]]^2)/(135135*c^13*r^5) - 
     (75939146*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^6*Scalar[vi[-a]*vi[a]]^2)/
      (45045*c^13*r^4) + (763721359*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^6*
       Scalar[vi[-a]*vi[a]]^2)/(45045*c^13*r^4) - 
     (60988178*G^3*m^4*\[Nu]^4*Scalar[ni[a]*vi[-a]]^6*Scalar[vi[-a]*vi[a]]^2)/
      (1287*c^13*r^4) + (220021336*G^3*m^4*\[Nu]^5*Scalar[ni[a]*vi[-a]]^6*
       Scalar[vi[-a]*vi[a]]^2)/(9009*c^13*r^4) - 
     (187807952*G^3*m^4*\[Nu]^6*Scalar[ni[a]*vi[-a]]^6*
       Scalar[vi[-a]*vi[a]]^2)/(45045*c^13*r^4) + 
     (2338013724323*G^5*m^6*\[Nu]^2*Scalar[vi[-a]*vi[a]]^3)/
      (993242250*c^13*r^6) - (123338*G^4*m^5*\[Nu]^2*Scalar[vi[-a]*vi[a]]^3)/
      (1155*c^11*r^5) + (752*G^3*m^4*\[Nu]^2*Scalar[vi[-a]*vi[a]]^3)/
      (35*c^9*r^4) - (329208584321*G^5*m^6*\[Nu]^3*Scalar[vi[-a]*vi[a]]^3)/
      (36786750*c^13*r^6) - (3765*G^5*m^6*Pi^2*\[Nu]^3*
       Scalar[vi[-a]*vi[a]]^3)/(112*c^13*r^6) + 
     (190642*G^4*m^5*\[Nu]^3*Scalar[vi[-a]*vi[a]]^3)/(315*c^11*r^5) - 
     (21988*G^3*m^4*\[Nu]^3*Scalar[vi[-a]*vi[a]]^3)/(315*c^9*r^4) + 
     (1522068269*G^5*m^6*\[Nu]^4*Scalar[vi[-a]*vi[a]]^3)/(1351350*c^13*r^6) + 
     (8979*G^5*m^6*Pi^2*\[Nu]^4*Scalar[vi[-a]*vi[a]]^3)/(70*c^13*r^6) - 
     (1910026*G^4*m^5*\[Nu]^4*Scalar[vi[-a]*vi[a]]^3)/(3465*c^11*r^5) + 
     (3544*G^3*m^4*\[Nu]^4*Scalar[vi[-a]*vi[a]]^3)/(63*c^9*r^4) - 
     (144390232*G^5*m^6*\[Nu]^5*Scalar[vi[-a]*vi[a]]^3)/(57915*c^13*r^6) + 
     (151216*G^4*m^5*\[Nu]^5*Scalar[vi[-a]*vi[a]]^3)/(693*c^11*r^5) + 
     (51294577*G^5*m^6*\[Nu]^6*Scalar[vi[-a]*vi[a]]^3)/(115830*c^13*r^6) - 
     (1144624*G^5*m^6*\[Nu]^2*Log[r/r0]*Scalar[vi[-a]*vi[a]]^3)/
      (3675*c^13*r^6) - (8704*G^5*m^6*\[Nu]^3*Log[r/r0]*
       Scalar[vi[-a]*vi[a]]^3)/(3675*c^13*r^6) - 
     (70568888*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]*Scalar[vi[-a]*vi[a]]^3)/
      (33075*c^12*r^5) + (235304336*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]*
       Scalar[vi[-a]*vi[a]]^3)/(33075*c^12*r^5) + 
     (9683045*G^4*m^5*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]]^3)/
      (14014*c^13*r^5) - (62998*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^3)/(315*c^11*r^4) - 
     (1932606463*G^4*m^5*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^3)/(210210*c^13*r^5) + 
     (4479652*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]]^3)/
      (3465*c^11*r^4) + (34574507*G^4*m^5*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^3)/(1170*c^13*r^5) - 
     (357608*G^3*m^4*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]]^3)/
      (165*c^11*r^4) - (93552772*G^4*m^5*\[Nu]^5*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^3)/(4095*c^13*r^5) + 
     (247840*G^3*m^4*\[Nu]^5*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]]^3)/
      (231*c^11*r^4) + (26250403*G^4*m^5*\[Nu]^6*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^3)/(4290*c^13*r^5) + 
     (3389068*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]]^3)/
      (3003*c^13*r^4) - (31163159*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]]^3)/(2730*c^13*r^4) + 
     (101981699*G^3*m^4*\[Nu]^4*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]]^3)/(3003*c^13*r^4) - 
     (129936544*G^3*m^4*\[Nu]^5*Scalar[ni[a]*vi[-a]]^4*
       Scalar[vi[-a]*vi[a]]^3)/(5005*c^13*r^4) + 
     (37360016*G^3*m^4*\[Nu]^6*Scalar[ni[a]*vi[-a]]^4*Scalar[vi[-a]*vi[a]]^3)/
      (5005*c^13*r^4) - (28306423*G^4*m^5*\[Nu]^2*Scalar[vi[-a]*vi[a]]^4)/
      (450450*c^13*r^5) + (16063*G^3*m^4*\[Nu]^2*Scalar[vi[-a]*vi[a]]^4)/
      (462*c^11*r^4) + (284386423*G^4*m^5*\[Nu]^3*Scalar[vi[-a]*vi[a]]^4)/
      (286650*c^13*r^5) - (694427*G^3*m^4*\[Nu]^3*Scalar[vi[-a]*vi[a]]^4)/
      (3465*c^11*r^4) - (302988173*G^4*m^5*\[Nu]^4*Scalar[vi[-a]*vi[a]]^4)/
      (90090*c^13*r^5) + (241634*G^3*m^4*\[Nu]^4*Scalar[vi[-a]*vi[a]]^4)/
      (693*c^11*r^4) + (133726466*G^4*m^5*\[Nu]^5*Scalar[vi[-a]*vi[a]]^4)/
      (45045*c^13*r^5) - (54352*G^3*m^4*\[Nu]^5*Scalar[vi[-a]*vi[a]]^4)/
      (231*c^11*r^4) - (121150651*G^4*m^5*\[Nu]^6*Scalar[vi[-a]*vi[a]]^4)/
      (90090*c^13*r^5) - (2997671*G^3*m^4*\[Nu]^2*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^4)/(9009*c^13*r^4) + 
     (31225996*G^3*m^4*\[Nu]^3*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]]^4)/
      (9009*c^13*r^4) - (169896742*G^3*m^4*\[Nu]^4*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^4)/(15015*c^13*r^4) + 
     (79102168*G^3*m^4*\[Nu]^5*Scalar[ni[a]*vi[-a]]^2*Scalar[vi[-a]*vi[a]]^4)/
      (6435*c^13*r^4) - (25653168*G^3*m^4*\[Nu]^6*Scalar[ni[a]*vi[-a]]^2*
       Scalar[vi[-a]*vi[a]]^4)/(5005*c^13*r^4) + 
     (2074808*G^3*m^4*\[Nu]^2*Scalar[vi[-a]*vi[a]]^5)/(45045*c^13*r^4) - 
     (1002763*G^3*m^4*\[Nu]^3*Scalar[vi[-a]*vi[a]]^5)/(2310*c^13*r^4) + 
     (63955069*G^3*m^4*\[Nu]^4*Scalar[vi[-a]*vi[a]]^5)/(45045*c^13*r^4) - 
     (85146046*G^3*m^4*\[Nu]^5*Scalar[vi[-a]*vi[a]]^5)/(45045*c^13*r^4) + 
     (473072*G^3*m^4*\[Nu]^6*Scalar[vi[-a]*vi[a]]^5)/(429*c^13*r^4)
 
FluxCirc45PN = (32*c^5*x^5*\[Nu]^2)/(5*G) - (2494*c^5*x^6*\[Nu]^2)/(105*G) + 
     (128*c^5*Pi*x^(13/2)*\[Nu]^2)/(5*G) - (89422*c^5*x^7*\[Nu]^2)/(2835*G) - 
     (8191*c^5*Pi*x^(15/2)*\[Nu]^2)/(105*G) + (6643739519*c^5*x^8*\[Nu]^2)/
      (10914750*G) - (54784*c^5*EulerGamma*x^8*\[Nu]^2)/(525*G) + 
     (512*c^5*Pi^2*x^8*\[Nu]^2)/(15*G) - (13028*c^5*Pi*x^(17/2)*\[Nu]^2)/
      (63*G) - (323105549467*c^5*x^9*\[Nu]^2)/(496621125*G) + 
     (3721552*c^5*EulerGamma*x^9*\[Nu]^2)/(11025*G) - 
     (21904*c^5*Pi^2*x^9*\[Nu]^2)/(315*G) + 
     (265978667519*c^5*Pi*x^(19/2)*\[Nu]^2)/(116424000*G) - 
     (219136*c^5*EulerGamma*Pi*x^(19/2)*\[Nu]^2)/(525*G) - 
     (56*c^5*x^6*\[Nu]^3)/(3*G) + (37084*c^5*x^7*\[Nu]^3)/(315*G) - 
     (2332*c^5*Pi*x^(15/2)*\[Nu]^3)/(15*G) - (134543*c^5*x^8*\[Nu]^3)/
      (1215*G) + (82*c^5*Pi^2*x^8*\[Nu]^3)/(15*G) + 
     (42949*c^5*Pi*x^(17/2)*\[Nu]^3)/(54*G) - (1452202403629*c^5*x^9*\[Nu]^3)/
      (229209750*G) + (1327296*c^5*EulerGamma*x^9*\[Nu]^3)/(1225*G) - 
     (267127*c^5*Pi^2*x^9*\[Nu]^3)/(720*G) + 
     (2062241*c^5*Pi*x^(19/2)*\[Nu]^3)/(3465*G) + 
     (328*c^5*Pi^3*x^(19/2)*\[Nu]^3)/(15*G) + (208*c^5*x^7*\[Nu]^4)/(9*G) - 
     (188806*c^5*x^8*\[Nu]^4)/(945*G) + (77354*c^5*Pi*x^(17/2)*\[Nu]^4)/
      (189*G) + (2571400*c^5*x^9*\[Nu]^4)/(1701*G) - 
     (3157*c^5*Pi^2*x^9*\[Nu]^4)/(60*G) - (26622581*c^5*Pi*x^(19/2)*\[Nu]^4)/
      (9072*G) - (1240*c^5*x^8*\[Nu]^5)/(81*G) + (5500*c^5*x^9*\[Nu]^5)/
      (63*G) - (3719141*c^5*Pi*x^(19/2)*\[Nu]^5)/(5940*G) + 
     (16*c^5*x^9*\[Nu]^6)/(3*G) - (109568*c^5*x^8*\[Nu]^2*Log[2])/(525*G) + 
     (638896*c^5*x^9*\[Nu]^2*Log[2])/(735*G) + 
     (15329984*c^5*x^9*\[Nu]^3*Log[2])/(11025*G) - 
     (9477*c^5*x^9*\[Nu]^2*Log[3])/(49*G) + (37908*c^5*x^9*\[Nu]^3*Log[3])/
      (49*G) - (27392*c^5*x^8*\[Nu]^2*Log[x])/(525*G) + 
     (1860776*c^5*x^9*\[Nu]^2*Log[x])/(11025*G) + 
     (663648*c^5*x^9*\[Nu]^3*Log[x])/(1225*G) - 
     (109568*c^5*Pi*x^(19/2)*\[Nu]^2*Log[16*x])/(525*G)
 
Freq4N = (-9965202491753717*Pi)/(23073008910336000*\[Tau]^(11/8)) + 
     (107*EulerGamma*Pi)/(2400*\[Tau]^(11/8)) + 
     (23*Pi^3)/(2400*\[Tau]^(11/8)) + (8248609881163*Pi*\[Nu])/
      (10987147100160*\[Tau]^(11/8)) - (3157*Pi^3*\[Nu])/
      (122880*\[Tau]^(11/8)) - (3590973803*Pi*\[Nu]^2)/
      (83235962880*\[Tau]^(11/8)) - (520159*Pi*\[Nu]^3)/
      (6539968512*\[Tau]^(11/8)) - (113868647*Pi)/(1734082560*\[Tau]^(9/8)) - 
     (31821*Pi*\[Nu])/(573440*\[Tau]^(9/8)) + (294941*Pi*\[Nu]^2)/
      (15482880*\[Tau]^(9/8)) - 10052469856691/(24034384281600*\[Tau]) + 
     (107*EulerGamma)/(1680*\[Tau]) + Pi^2/(24*\[Tau]) + 
     (3147553127*\[Nu])/(3121348608*\[Tau]) - (451*Pi^2*\[Nu])/
      (12288*\[Tau]) - (15211*\[Nu]^2)/(1769472*\[Tau]) + 
     (25565*\[Nu]^3)/(1327104*\[Tau]) - (11891*Pi)/(215040*\[Tau]^(7/8)) + 
     (109*Pi*\[Nu])/(7680*\[Tau]^(7/8)) + 19583/(1016064*\[Tau]^(3/4)) + 
     (24401*\[Nu])/(774144*\[Tau]^(3/4)) + (31*\[Nu]^2)/(1152*\[Tau]^(3/4)) - 
     Pi/(20*\[Tau]^(5/8)) + 743/(16128*Sqrt[\[Tau]]) + 
     (11*\[Nu])/(192*Sqrt[\[Tau]]) + 1/(4*\[Tau]^(1/4)) - 
     (107*Pi*Log[\[Tau]/256])/(19200*\[Tau]^(11/8)) - 
     (107*Log[\[Tau]/256])/(13440*\[Tau]) - (2518977598355703073*Log[\[Tau]])/
      (15117435438052147200*\[Tau]^(5/4)) + (9203*EulerGamma*Log[\[Tau]])/
      (860160*\[Tau]^(5/4)) + (9049*Pi^2*Log[\[Tau]])/
      (1032192*\[Tau]^(5/4)) + (718143266031997*\[Nu]*Log[\[Tau]])/
      (2307300891033600*\[Tau]^(5/4)) + (244493*EulerGamma*\[Nu]*Log[\[Tau]])/
      (4515840*\[Tau]^(5/4)) - (65577*Pi^2*\[Nu]*Log[\[Tau]])/
      (7340032*\[Tau]^(5/4)) - (1502014727*\[Nu]^2*Log[\[Tau]])/
      (33294385152*\[Tau]^(5/4)) + (2255*Pi^2*\[Nu]^2*Log[\[Tau]])/
      (1572864*\[Tau]^(5/4)) - (258479*\[Nu]^3*Log[\[Tau]])/
      (132120576*\[Tau]^(5/4)) + (1195*\[Nu]^4*Log[\[Tau]])/
      (1048576*\[Tau]^(5/4)) + (14873*Log[2]*Log[\[Tau]])/
      (4515840*\[Tau]^(5/4)) + (15761*\[Nu]*Log[2]*Log[\[Tau]])/
      (188160*\[Tau]^(5/4)) + (47385*Log[3]*Log[\[Tau]])/
      (6422528*\[Tau]^(5/4)) - (47385*\[Nu]*Log[3]*Log[\[Tau]])/
      (1605632*\[Tau]^(5/4)) - (9203*Log[\[Tau]]^2)/(13762560*\[Tau]^(5/4)) - 
     (244493*\[Nu]*Log[\[Tau]]^2)/(72253440*\[Tau]^(5/4))
 
Phase45PN = -55/(384*x^(3/2)) - 27145/(32256*Sqrt[x]) + 
     (15737765635*Sqrt[x])/390168576 - (2255*Pi^2*Sqrt[x])/1536 - 
     (378515*Pi*x)/387072 + (680712846248317*x^(3/2))/10815472926720 + 
     (244493*EulerGamma*x^(3/2))/21168 - (109295*Pi^2*x^(3/2))/57344 - 
     (1492917260735*Pi*x^2)/34334834688 + (2255*Pi^3*x^2)/1536 - 
     1/(32*x^(5/2)*\[Nu]) - 3715/(32256*x^(3/2)*\[Nu]) + 
     (5*Pi)/(16*x*\[Nu]) - 15293365/(32514048*Sqrt[x]*\[Nu]) - 
     (12348611926451*Sqrt[x])/(600859607040*\[Nu]) + 
     (107*EulerGamma*Sqrt[x])/(42*\[Nu]) + (5*Pi^2*Sqrt[x])/(3*\[Nu]) - 
     (77096675*Pi*x)/(65028096*\[Nu]) - (2550713843998885153*x^(3/2))/
      (70862978615869440*\[Nu]) + (9203*EulerGamma*x^(3/2))/(4032*\[Nu]) + 
     (45245*Pi^2*x^(3/2))/(24192*\[Nu]) + (93098188434443*Pi*x^2)/
      (4806876856320*\[Nu]) - (107*EulerGamma*Pi*x^2)/(42*\[Nu]) - 
     (5*Pi^3*x^2)/(6*\[Nu]) - (3085*\[Nu])/(4608*Sqrt[x]) - 
     (76055*Sqrt[x]*\[Nu])/221184 + (74045*Pi*x*\[Nu])/193536 - 
     (7510073635*x^(3/2)*\[Nu])/780337152 + (11275*Pi^2*x^(3/2)*\[Nu])/
      36864 + (45293335*Pi*x^2*\[Nu])/32514048 + (127825*Sqrt[x]*\[Nu]^2)/
      165888 - (1292395*x^(3/2)*\[Nu]^2)/3096576 + (10323755*Pi*x^2*\[Nu]^2)/
      51093504 + (5975*x^(3/2)*\[Nu]^3)/24576 + \[Psi]0 - 
     (78975*x^(3/2)*Log[3])/12544 + (78975*x^(3/2)*Log[3])/(50176*\[Nu]) + 
     (622757*x^(3/2)*Log[16])/84672 + (107*Sqrt[x]*Log[16])/(84*\[Nu]) + 
     (252755*x^(3/2)*Log[16])/(338688*\[Nu]) - (107*Pi*x^2*Log[16])/
      (84*\[Nu]) + (65*Pi*Log[x])/512 + (244493*x^(3/2)*Log[x])/42336 - 
     (38645*Pi*Log[x])/(43008*\[Nu]) + (107*Sqrt[x]*Log[x])/(84*\[Nu]) + 
     (9203*x^(3/2)*Log[x])/(8064*\[Nu]) - (107*Pi*x^2*Log[x])/(84*\[Nu])
 
PhaseSPA45PN = 2*f*Pi*T0 + 55/(384*v^3) + 27145/(21504*v) - 
     (15737765635*v)/130056192 + (2255*Pi^2*v)/512 + (378515*Pi*v^2)/64512 - 
     (1492917260735*Pi*v^4)/5722472448 + (2255*Pi^3*v^4)/256 + 
     3/(128*v^5*\[Nu]) + 3715/(32256*v^3*\[Nu]) - (3*Pi)/(8*v^2*\[Nu]) + 
     15293365/(21676032*v*\[Nu]) + (11583231236531*v)/(200286535680*\[Nu]) - 
     (107*EulerGamma*v)/(14*\[Nu]) - (5*Pi^2*v)/\[Nu] + 
     (77096675*Pi*v^2)/(10838016*\[Nu]) + (105344279473163*Pi*v^4)/
      (801146142720*\[Nu]) - (107*EulerGamma*Pi*v^4)/(7*\[Nu]) - 
     (5*Pi^3*v^4)/\[Nu] + (3085*\[Nu])/(3072*v) + (76055*v*\[Nu])/73728 - 
     (74045*Pi*v^2*\[Nu])/32256 + (45293335*Pi*v^4*\[Nu])/5419008 - 
     (127825*v*\[Nu]^2)/55296 + (10323755*Pi*v^4*\[Nu]^2)/8515584 - 
     \[CapitalPsi]0 - (107*v*Log[2])/(7*\[Nu]) - (214*Pi*v^4*Log[2])/
      (7*\[Nu]) - (65*Pi*Log[v])/128 + (680712846248317*v^3*Log[v])/
      1802578821120 + (244493*EulerGamma*v^3*Log[v])/3528 - 
     (327885*Pi^2*v^3*Log[v])/28672 + (38645*Pi*Log[v])/(10752*\[Nu]) - 
     (107*v*Log[v])/(14*\[Nu]) - (2550713843998885153*v^3*Log[v])/
      (11810496435978240*\[Nu]) + (9203*EulerGamma*v^3*Log[v])/(672*\[Nu]) + 
     (45245*Pi^2*v^3*Log[v])/(4032*\[Nu]) - (107*Pi*v^4*Log[v])/(7*\[Nu]) - 
     (7510073635*v^3*\[Nu]*Log[v])/130056192 + (11275*Pi^2*v^3*\[Nu]*Log[v])/
      6144 - (1292395*v^3*\[Nu]^2*Log[v])/516096 + (5975*v^3*\[Nu]^3*Log[v])/
      4096 + (622757*v^3*Log[2]*Log[v])/3528 + (252755*v^3*Log[2]*Log[v])/
      (14112*\[Nu]) - (236925*v^3*Log[3]*Log[v])/6272 + 
     (236925*v^3*Log[3]*Log[v])/(25088*\[Nu]) + (244493*v^3*Log[v]^2)/7056 + 
     (9203*v^3*Log[v]^2)/(1344*\[Nu])
 
H224PN = 1 - (107*x)/42 + 2*Pi*x^(3/2) - (2173*x^2)/1512 - 
     (107*Pi*x^(5/2))/21 + (27027409*x^3)/646800 - (856*EulerGamma*x^3)/105 + 
     ((428*I)/105)*Pi*x^3 + (2*Pi^2*x^3)/3 - (2173*Pi*x^(7/2))/756 - 
     (846557506853*x^4)/12713500800 + (45796*EulerGamma*x^4)/2205 - 
     ((22898*I)/2205)*Pi*x^4 - (107*Pi^2*x^4)/63 + (55*x*\[Nu])/42 - 
     (1069*x^2*\[Nu])/216 - (24*I)*x^(5/2)*\[Nu] + (34*Pi*x^(5/2)*\[Nu])/21 - 
     (278185*x^3*\[Nu])/33264 + (41*Pi^2*x^3*\[Nu])/96 + 
     ((14333*I)/162)*x^(7/2)*\[Nu] - (2495*Pi*x^(7/2)*\[Nu])/378 - 
     (336005827477*x^4*\[Nu])/4237833600 + (15284*EulerGamma*x^4*\[Nu])/441 - 
     ((219314*I)/2205)*Pi*x^4*\[Nu] - (9755*Pi^2*x^4*\[Nu])/32256 + 
     (2047*x^2*\[Nu]^2)/1512 - (20261*x^3*\[Nu]^2)/2772 - 
     ((4066*I)/945)*x^(7/2)*\[Nu]^2 + (40*Pi*x^(7/2)*\[Nu]^2)/27 + 
     (256450291*x^4*\[Nu]^2)/7413120 - (1025*Pi^2*x^4*\[Nu]^2)/1008 + 
     (114635*x^3*\[Nu]^3)/99792 - (81579187*x^4*\[Nu]^3)/15567552 + 
     (26251249*x^4*\[Nu]^4)/31135104 - (428*x^3*Log[16*x])/105 + 
     (22898*x^4*Log[16*x])/2205 + (7642*x^4*\[Nu]*Log[16*x])/441
