;==================================== SPDUST.2 ==========================================;
;                                                                                        ;
;                               H2_PHOTOEMISSION.PRO                                     ;                      
;                                                                                        ;
; Rotational excitation rates through H2 formation and photoemission                     ;
; of electrons.                                                                          ;
; Refers to Sections 8-9 of Ali-Haimoud, Hirata & Dickinson, 2009, MNRAS, 395, 1055.     ; 
;                                                                                        ;         
; Author : Yacine Ali-Haimoud        yacine@tapir.caltech.edu                            ; 
;                                                                                        ;
; Revision history : Written December 2008                                               ;
;                    Modified November 2009 (changed aexc -> acx in function GH2)        ;
;                      the rest is left unchanged from SPDUST.1                          ;
;========================================================================================;



function GH2, env, a
common cgsconst

; ---> Excitation rate due to random H2 formation

  T     = env.T
  nh    = env.nh
  ax    = acx(a)
  hbar  = h/(2d*pi)
  y     = env.y
  gamma = env.gamma

  Ef = 0.2d * eV
  JJplus1 = 1d2

  return, gamma/4d * (1d - y)*Ef/(k*T) * (1d + JJplus1*hbar^2/(2d * mp * Ef * ax^2) )

end

;-----------------------------------------------------------------------------------------

function FGpeZ, env, a, Z
common cgsconst
common refr_indices, hnu_tab, la_tab

;--- parameter ---
  Nnu = 500
;-----------------

  T   = env.T
  nh = env.nh
  Chi = env.Chi
  as = asurf(a)

  Jpeisrf = Jpeisrf(a)

  Jpepos = Chi *Jpeisrf.Jpepos
  Jpeneg = Chi *Jpeisrf.Jpeneg

; --- Gpe --- 

  hnu_pet = max([1d-5, hnu_pet(Z,a)])
  hnu_pdt = max([1d-5, hnu_pdt(Z,a)])
  hnu_pet_tab = makelogtab(hnu_pet, 13.6d, Nnu)
  Dnu_over_nu_pet = DX_over_X(hnu_pet, 13.6d, Nnu)
  hnu_pdt_tab = makelogtab(hnu_pdt, 13.6d, Nnu) 
  Dnu_over_nu_pdt = DX_over_X(hnu_pdt, 13.6d, Nnu)

  Ytab = interpol(Y(Z,a), alog(hnu_tab), alog(hnu_pet_tab))
  Qtab = Qabs(a, Z, hnu_pet_tab)

  if( Z GE 0 ) then begin
     E_low  = -(Z + 1d) *q^2/a 
     E_high = (hnu_pet_tab - hnu_pet) *eV  
     Epe    = (0.5d* E_high * (E_high - 2d*E_low)/(E_high - 3d*E_low)) ; in erg
     second_term = 0d 
     third_term = Jpepos[Z]* (Z + 1d)*q^2/as
  endif else begin
     E_min   = E_min(Z,a)
     E_low  = E_min* eV
     E_high = (E_min + hnu_pet_tab - hnu_pet) * eV
     Epe    = 0.5d* (E_high + E_low) ; in erg
     second_term = Dnu_over_nu_pdt * total(sigma_pdt(hnu_pdt_tab, Z, a) * nu_uisrf(hnu_pdt_tab)/eV /hnu_pdt_tab $
                                           * (hnu_pdt_tab - hnu_pdt + E_min), /double )
     third_term = Jpeneg[-Z]* (Z + 1d)*q^2/as
  endelse

  first_term = c * pi* a^2 * Dnu_over_nu_pet * total(Ytab * Qtab * nu_uisrf(hnu_pet_tab)/eV /hnu_pet_tab * Epe , /double)

  Gpe =  me/(4d * nh * sqrt(8d *pi * mp *k*T) *as^2 *k*T ) * (first_term + second_term + third_term)

; --- Fpe ---

  if Z GE 0 then begin
     Jpe = Jpepos[Z]
  endif else begin
     Jpe = Jpeneg[-Z]
  endelse

  Fpe = me/mp * Jpe/(2d *pi *as^2 *nh *sqrt(2d *k *T/(pi *mp)))

  return, {Fpe: Fpe, Gpe: Gpe}

end

;-----------------------------------------------------------------------------------------

function FGpe_averaged, env, a, fZ

  NZ = N_elements(fZ[0,*])
  Fpe = 0d
  Gpe = 0d

  for i = 0, NZ - 1 do begin
     FGpe = FGpeZ(env, a, fZ[0,i]) 
     Fpe = Fpe + fZ[1,i] * FGpe.Fpe
     Gpe = Gpe + fZ[1,i] * FGpe.Gpe
  endfor

  return, {Fpe: Fpe, Gpe: Gpe}

end

;-----------------------------------------------------------------------------------------



