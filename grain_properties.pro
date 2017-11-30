;==================================== SPDUST.2  ====================================;
;                                                                                   ;
;                               GRAIN_PROPERTIES.PRO                                ;
;                                                                                   ;
; Grain geometrical properties, dipole moments and size distribution function.      ;
; Refers to Section 3 of Ali-Haimoud, Hirata & Dickinson, 2009, MNRAS, 395, 1055.   ;
;                                                                                   ;      
; Author : Yacine Ali-Haimoud        yacine@tapir.caltech.edu                       ; 
;                                                                                   ;
; Revision history : Written December 2008                                          ;
;                    Revised March 2010 (Minor changes                              ;
;                    commented. Removed the linear grain option, not used).         ;
;                                                                                   ;
;===================================================================================;



pro physconst
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Stores the needed constants in cgs units

    pi    = 3.1415926535897932385d
    c     = 2.99792458d10           ; speed of light
    q     = 4.8032068d-10           ; elementary charge
    k     = 1.380650d-16            ; Boltzmann constant
    mp    = 1.6726231d-24           ; proton mass
    me    = 9.1093898d-28           ; electron mass
    h     = 6.6260688d-27           ; Planck constant
    debye = 1d-18                   ; 1 Debye in cgs units
    eV    = 1.60217653d-12          ; 1eV in ergs

end

;----------------------------------------------------------------------------

pro parameters
common grainparams, a2, d, rho, epsilon

; ---> Stores some grain geometrical properties

    a2      = 6d-8           ; radius under which grains are disklike
    d       = 3.35d-8        ; disk thickness (graphite interlayer separation)
    rho     = 2.24d          ; carbon density in g/cm^3
    epsilon = 0.01           ; charge centroid displacement

end

;----------------------------------------------------------------------------

function N_C, a
common cgsconst   , pi, c, q, k, mp, me,  h, debye, eV
common grainparams, a2, d, rho, epsilon

; ---> Number of carbon atoms in a grain

    return, floor(4d*pi/3d * a^3 *rho/(12d *mp)) + 1

end

;----------------------------------------------------------------------------

function N_H, a

; ---> Number of hydrogen atoms in a grain
; Eq (8) of Draine & Li, 2001, ApJ 551, 807

    Nc = N_C(a)

    if Nc LT 25 then begin
       return, floor(0.5 * Nc + 0.5)
    endif 

    if Nc LT 100 then begin
       return, floor(2.5 * sqrt(Nc) + 0.5)
    endif 
   
    return, floor(0.25 * Nc + 0.5)

end

;----------------------------------------------------------------------------

function Inertia, a
common cgsconst   , pi, c, q, k, mp, me,  h, debye, eV
common grainparams, a2, d, rho, epsilon

; ---> Largest moment of inertia

    Mass    = (12d * N_C(a) + N_H(a))*mp
    Isphere = 0.4d * Mass *a^2

    if (a LE a2) then begin
       return, 5d/3d * a/d * Isphere         ;DL98b (A8)
    endif 
   
    return, Isphere


end

;----------------------------------------------------------------------------

function asurf, a     
common grainparams, a2, d, rho, epsilon

; ---> Surface-equivalent radius, needed for grain charge calculation

    if (a LE a2) then begin
      b = sqrt(4d/3d * a^3/d)               ;DL98(A6)
      return, sqrt(b^2 /2d + b*d /2d)
    endif

    return, a

end

;----------------------------------------------------------------------------

function acx, a      
common grainparams, a2, d, rho, epsilon

; ---> Cylindrical excitation equivalent radius
; (AHD09) Eq. (4) and (65)
; MODIFIED FEB 10 : take the limit d-> 0 (thin disk) 

    if (a LE a2) then begin
          R = sqrt(4d/3d * a^3/d)  ; radius of the thin disk of equivalent "volume" (really, mass)
      return, (3d/8d)^0.25d * R ; (3d/8d *R^3 *(2d * d + R))^0.25d
    endif

    return, a

end

;----------------------------------------------------------------------------

function rms_dipole, a, Z2, beta                      
common cgsconst   , pi, c, q, k, mp, me,  h, debye, eV
common grainparams, a2, d, rho, epsilon

; Modified February 10. 
; Because we now need to account for both parallel and
; perpendicular components of the dipole moment separately, I will
; approximate the Z-averages by <mu(Z) * function(Z)> \approx
; mu(sqrt(<Z^2>)) * <function(Z)>.
; This is justified because the Z-dependent part is small compared to
; the intrinsic dipole part anyway.

   muZ     = epsilon *sqrt(Z2) *q *acx(a) ; charge induced part
   N_at    = N_C(a) + N_H(a) 
   return, sqrt(N_at *beta^2 + muZ^2)   
 
end

;----------------------------------------------------------------------------

function size_dist, a
common cgsconst,         pi, c, q, k, mp, me, h, debye, eV
common grainparams,      a2, d, rho, epsilon
common size_dist_params, bc, alpha_g, beta_g, at_g, ac_g, C_g

;  ---> Grain size distribution, using Weingartner & Draine, 2001a prescription. The line of their
; table 1 is given by the user in param_file.


; --- Parameters of the lognormal populations, see WD01a
  mc    = 12d*mp
  bci   = [0.75, 0.25]*bc
  a0i   = [3.5, 30.]*1d-8
  sigma = 0.4d
  amin  = 3.5d-8
;--- --- --- --- --- --- --- --- --- --- --- --- --- ---

  Bi   = 3d/(2d*pi)^1.5 * exp(-4.5d*sigma^2)/(rho*a0i^3*sigma) * bci *mc/(1d + erf(3d*sigma/sqrt(2d) + alog(a0i/amin)/(sigma * sqrt(2d))))
  D_a  = total(Bi/a *exp(-0.5d*(alog(a/a0i)/sigma)^2), /double) 

  if(beta_g GE 0) then begin
     F_a = 1. + beta_g *a/at_g
  endif else begin
     F_a = 1./(1. - beta_g *a/at_g)
  endelse

  cutoff = 1d 
  if (a GT at_g) then begin
     cutoff = exp(-((a - at_g)/ac_g)^3)
  endif

  return, D_a  + C_g/a * (a/at_g)^alpha_g * F_a * cutoff

end

;----------------------------------------------------------------------------



