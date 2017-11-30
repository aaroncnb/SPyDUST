;=============================== SPDUST.2 =====================================;
;                                                                              ;
;                              PLASMADRAG.PRO                                  ;
;                                                                              ;
; Normalized rotational excitation rates for plasma drag.                      ;
; See section 6 of Ali-Haimoud, Hirata & Dickinson, 2009, MNRAS, 395, 1055     ;
; and section 5 of Silsbee, Ali-Haimoud & Hirata, 2010.                        ;             
;                                                                              ;
; Author : Yacine Ali-Haimoud        yacine@tapir.caltech.edu                  ;
;                                                                              ;
; Revision history : Written December 2008                                     ;
;                    Revised March 2010  (tumbling disklike grains)            ;
;==============================================================================;


pro compute_small_e
common cgsconst,   pi, c, q, k, mp, me, h, debye, eV 
common smalletabs, smalle_tab, Gamma_tab, Gamma_max

; ---> Computes and stores I(Zg < 0, e_1 <<1), as a function of 
; Gamma = omega*b/v *(e-1)

;-- creating arrays----
  gamma_min = 1d-20
  gamma_max = 10d^1.5
  Ngamma    = 300
  Gamma_tab = makelogtab(gamma_min, gamma_max, Ngamma)
  smalle_tab= dblarr(Ngamma)
;----------------------

  I = complex(0d, 1d, /double)

;----parameters---
  Ny = 20000L
  ymed = 1d
  ymax = 1d20                   ; sqrt(2d/e_1)
;-----------------

  u    = exp(I*pi/6d) * [maketab(0d, ymed, Ny), makelogtab(ymed, ymax, Ny)]
  Du   = exp(I*pi/6d) * [ymed/Ny + dblarr(Ny) , DX_over_X(ymed, ymax, Ny)*makelogtab(ymed, ymax, Ny)]
  fcos = (u^2 - 1d)/(u^2 + 1d)^2
  fsin = u/(u^2 + 1d)^2

  for ig = 0, Ngamma - 1 do begin
     Gamma = Gamma_tab[ig]
     time = Gamma * (u + u^3/3d)
     intcos = 4d *total(Real_part(exp(I*time) *fcos *Du), /double)^2
     intsin = 16d *total(Imaginary(exp(I*time) *fsin *Du), /double)^2
     smalle_tab[ig] = intcos + intsin
  endfor

  gamma_max = max(Gamma_tab)

  print, 'I(Zg<0, parabolic) stored'

end

;-----------------------------------------------------------------------------------------

pro  compute_int_plasma
common cgsconst,    pi, c, q, k, mp, me, h, debye, eV
common smalletabs,  smalle_tab, Gamma_tab, Gamma_max
common plasma_tabs, Ipos_tab, Ineg_tab, rot_tab, e_1tab, rot_min, rot_max, e_1min, e_1max

; ---> Computes and stores I(Zg, e-1), the reduced excitation
; coefficient for plasma drag

  compute_small_e               ; getting the values for small e_1, Zg < 0

;-- creating arrays--- --- --- --- --- --
  e_1min  = 1d-15
  e_1max  = 1d4
  Ne_1    = 100
  rot_min = 1d-7
  rot_max = 1d/e_1min * Gamma_max
  Nrot    = 100

  e_1small  = 1d-2              ; use normal expression for rot < rot_small
  rot_small = 1d3               ; use smalle expression for e_1 < e_1small and rot > rot_small

  rot_tab  = makelogtab(rot_min, rot_max, Nrot)
  e_1tab   = makelogtab(e_1min, e_1max, Ne_1)
  Ipos_tab = dblarr(Nrot, Ne_1)
  Ineg_tab = dblarr(Nrot, Ne_1)
;--- --- --- --- --- --- --- --- --- --- ---

  I = complex(0d, 1d, /double)

;----parameters---
  Ny = 10000L
  ymax = 1d50
  ymed = 800d
;-----------------

  y = [maketab(0d, ymed, Ny), makelogtab(ymed, ymax, Ny)]

;---- Zg > 0 -----------

  Dz = -I* [ymed/Ny + dblarr(Ny) , DX_over_X(ymed, ymax, Ny)*makelogtab(ymed, ymax, Ny)]
  log = 0.5d *alog(1d + 4d/y^2) + I *atan(2d/y)
  z = 1d - I*y

  for ie = 0, Ne_1 -1 do begin
     e_1 = e_1tab[ie]
     A = e_1/(e_1 + 2d)
     time = 1d/sqrt(e_1*(e_1 + 2d)) * (log + 2d*(e_1 + 1d) *z/(z^2 -1d))
     fcos = (z^2 - A)/(z^2 + A)^2
     fsin = z/(z^2 + A)^2
     for ir = 0, Nrot -1 do begin
        rot = rot_tab[ir]
        intcos = 4d *A *total(Real_part(exp(I*rot*time) *fcos *Dz), /double)^2
        intsin = 16d *A^2 *total(Imaginary(exp(I*rot*time) *fsin *Dz), /double)^2
        Ipos_tab[ir, ie] = intcos + intsin
     endfor
  endfor


  inde   = max(where(e_1tab LT e_1small))
  indrot = max(where(rot_tab LT rot_small))

;---- Zg < 0 non small e_1 ------------

  Dz = exp(I*pi/4d)* [ymed/Ny + dblarr(Ny) , DX_over_X(ymed, ymax, Ny)*makelogtab(ymed, ymax, Ny)]
  log = 0.5d *alog(1d + 4d/y^2 + 2d*sqrt(2d)/y) - I *atan(sqrt(2d)/(sqrt(2d)+y))
  z = 1d + exp(I*pi/4d)*y

  for ie = 0, Ne_1 - 1 do begin
     e_1 = e_1tab[ie]
     A = e_1/(e_1 + 2d)
     time = 1d/sqrt(e_1*(e_1 + 2d)) * (log - 2d*(e_1 + 1d) *z/(z^2 -1d))
     fcos = (1d - A*z^2)/(1d + A*z^2)^2
     fsin = z/(1 + A*z^2)^2
     for ir = 0, indrot - 1 do begin
        rot = rot_tab[ir]
        intcos = 4d * A *total(Real_part(exp(I*rot*time) *fcos *Dz), /double)^2
        intsin = 16d * A^2  *total(Imaginary(exp(I*rot*time) *fsin *Dz), /double)^2
        Ineg_tab[ir, ie] = intcos + intsin
     endfor
  endfor

;---- extending Zg < 0 for nearly parabolic case

  for ie = 0, inde - 1 do begin
     e_1 = e_1tab[ie]
     for ir = indrot, Nrot-1 do begin
        rot = rot_tab[ir]
        Gamma = rot *e_1
        if(Gamma LT Gamma_max) then begin
           Ineg_tab[ir, ie] = interpol(smalle_tab, Gamma_tab, Gamma)
        endif
     endfor
  endfor

;--- Some useful things for further interpolation

  rot_min  = min(rot_tab)   
  rot_max  = max(rot_tab)
  e_1min   = min(e_1tab)
  e_1max   = max(e_1tab)

  ind = where(Ipos_tab EQ 0, count)
  if(count NE 0) then begin
     Ipos_tab[ind] = 1d-30      ; the precision of those integrals is of order 10^-15 to 10^-20 anyway
  endif
  ind = where(Ineg_tab EQ 0, count)
  if(count NE 0) then begin
     Ineg_tab[ind] = 1d-30    
  endif

  print, 'I(rot, e, Zg <> 0) stored'
end

;-----------------------------------------------------------------------------------------

function int_plasma, rot_new, e_1_new, Zg
common plasma_tabs, Ipos_tab, Ineg_tab, rot_tab, e_1tab, rot_min, rot_max, e_1min, e_1max

; ---> Returns I(omega * b/v, e-1, Zg non zero) (interpolates the arrays
; calculated by compute_int_plasma)
; !!! rot_new and e_1_new must be arrays of the same dimension !!!
; returns an array of that dimension containing the interpolated integrals

  Nelem  = N_elements(rot_new)  ; = N_elements(e_1_new) as well !
  if(Nelem NE N_elements(e_1_new)) then begin
     print, 'error in int_plasma : e_1_new and rot_new should have the same dimensions'
     print, Nelem, '(rot_new)'
     print, N_elements(e_1_new), '(e_1_new)'
  endif
  result = dblarr(Nelem)

  ind = where(e_1_new GE e_1max and rot_new LT 100d, count)
  if (count NE 0) then begin
     rot = rot_new[ind]
     result[ind] = rot^2* (beselk(rot, 0)^2 + beselk(rot, 1)^2)
  endif

  ind = where(rot_new LE rot_min, count)
  if (count NE 0) then begin    ; note that there is some overlap in the region rot -> 0, e -> infinity, where I = 1
     e_1 = e_1_new[ind]
     result[ind] =  e_1*(e_1 + 2d)/(e_1 + 1d)^2
  endif

  ind = where((rot_new LT rot_max) and (rot_new GT rot_min) and (e_1_new LT e_1max) and (e_1_new GT e_1min), count) 
  if (count NE 0) then begin    ; interpolating the computed arrays
     rot = rot_new[ind]
     e_1 = e_1_new[ind]
     Drot_over_rot = alog(rot_tab[1]/rot_tab[0])
     irot = floor(alog(rot/min(rot_tab))/Drot_over_rot) ; rot_tab[irot[*]] < rot[*] < rot_tab[irot[*]+1]  
     De_over_e = alog(e_1tab[1]/e_1tab[0])
     ie = floor(alog(e_1/min(e_1tab))/De_over_e) ; e_1tab[ie[*]] < e_1[*] < e_1tab[ie[*]+1] 
     
     alpha = alog(rot/rot_tab[irot])/Drot_over_rot
     beta  = alog(e_1/e_1tab[ie])/De_over_e
     if (Zg GT 0) then begin
        result[ind] = exp(alpha*(beta* alog(Ipos_tab[irot+1, ie+1])+ (1d - beta)* alog(Ipos_tab[irot + 1, ie])) $
                          + (1d - alpha)*(beta* alog(Ipos_tab[irot, ie + 1])+(1d - beta) * alog(Ipos_tab[irot, ie])))
     endif else begin
        result[ind] = exp(alpha*(beta* alog(Ineg_tab[irot+1, ie+1])+ (1d - beta)* alog(Ineg_tab[irot + 1, ie])) $
                          + (1d - alpha)*(beta* alog(Ineg_tab[irot, ie + 1])+(1d - beta) * alog(Ineg_tab[irot, ie])))
     endelse
  endif

  return, result

end

;-----------------------------------------------------------------------------------------

function little_gp_charged, psi, Omega

  u_min = 1d-10
  u_max = 5d
  Nu    = 250L
  c_min = 1d-10
  c_max = 5d10
  Nc    = 250L

  u_1d      =  makelogtab(u_min, u_max, Nu)
  Du_over_u = DX_over_X(u_min, u_max, Nu)

  Ntot  = Nc*Nu                 ; work with arrays of Nc*Nu indices. index c*Nu+l corresponds to column c, line l
  u_arr = dblarr(Ntot)
  for col = 0, Nc-1 do begin
     u_arr[col*Nu + indgen(Nu)] = u_1d
  endfor
  c_arr     = dblarr(Ntot)
  Dc_over_c = dblarr(Ntot)      ; will contain Dc/c if no collision, 0 if collision 

; --- Positively charged grains ---

  if psi GT 0 then begin

     for iu = 0, Nu-1 do begin
        u = u_1d[iu]
        if u^2 LE psi then begin
           c1 = c_min
        endif else begin
           c1 = sqrt(1d - psi/u^2)
        endelse
        if (c1 LT c_max) then begin
           c_arr[iu + Nu * indgen(Nc)] = makelogtab(c1, c_max, Nc)
           Dc_over_c[iu + Nu * indgen(Nc)] = DX_over_X(c1, c_max, Nc)
        endif
     endfor

     stuff  = (2d/psi *c_arr *u_arr^2)^2
     e_1arr = sqrt(1d + stuff) - 1d
     ind    = where(stuff LT 1d-10, count)
     if(count NE 0) then begin
        e_1arr[ind] = 0.5d *stuff[ind]
     endif

     rot_arr = Omega * c_arr /u_arr
     Int_arr = int_plasma(rot_arr, e_1arr, 1d)    
     
     return,  2d *Du_over_u *total(u_arr^2*exp(-u_arr^2) * Int_arr * Dc_over_c , /double)
  endif

; --- Negatively charged grains ---

  if psi LT 0 then begin

     for iu = 0, Nu-1 do begin
        u = u_1d[iu]
        c1 = sqrt(1d - psi/u^2)
        if (c1 LT c_max) then begin
           c_arr[iu + Nu * indgen(Nc)] = makelogtab(c1, c_max, Nc)
           Dc_over_c[iu + Nu * indgen(Nc)] = DX_over_X(c1, c_max, Nc)
        endif
     endfor

     stuff  = (2d/psi *c_arr *u_arr^2)^2
     e_1arr = sqrt(1d + stuff) - 1d
     ind    = where(stuff LT 1d-10, count)
     if(count NE 0) then begin
        e_1arr[ind] = 0.5d *stuff[ind]
     endif

     rot_arr = Omega * c_arr /u_arr
     Int_arr = int_plasma(rot_arr, e_1arr, -1d)    
     
     return,  2d *Du_over_u *total(u_arr^2*exp(-u_arr^2) * Int_arr * Dc_over_c , /double)
  endif

end

;-----------------------------------------------------------------------------------------

function little_gp_neutral, phi, Omega

  u_min = 1d-10
  u_max = 5d
  Nu    = 250L
  X_max = 20d
  NX    = 250L

  u_tab     =  makelogtab(u_min, u_max, Nu)
  Du_over_u = DX_over_X(u_min, u_max, Nu)
  Int_c     = dblarr(Nu) 

  for iu = 0, Nu - 1 do begin
     u = u_tab[iu]
     X_min = Omega/u *sqrt(1d + phi/u)
     if X_min LT X_max then begin
        X_tab = makelogtab(X_min, X_max, NX)
        DX_over_X = DX_over_X(X_min, X_max, NX)
        Int_c[iu] = total(X_tab^2*(Beselk(X_tab, 0)^2 + Beselk(X_tab, 1)^2), /double) * DX_over_X
     endif
  endfor

  return,  2d *Du_over_u *total(u_tab^2*exp(-u_tab^2) * Int_c, /double)

end

;-----------------------------------------------------------------------------------------

pro set_up_gp_arrays
common gp_arrays, psi_min, psi_max, Npsi, phi_min, phi_max, Nphi, Omega_min, Omega_max, NOmega, gp_pos, gp_neg, gp_neutral

; ---> Sets up the values of the arrays over which little_gp is interpolated

  psi_min = 1d-5
  psi_max = 1d6
  Npsi    = 110

  phi_min = 1d-2
  phi_max = 5d2
  Nphi    = 300

  Omega_min = 1d-10
  Omega_max = 1d5
  NOmega    = 150


end

;-----------------------------------------------------------------------------------------

pro compute_little_gp
common gp_arrays, psi_min, psi_max, Npsi, phi_min, phi_max, Nphi, Omega_min, Omega_max, NOmega, gp_pos, gp_neg, gp_neutral
common path, SpDust_data_dir

; ---> Computes arrays for gp(psi, Omega,  +or-) and stores
; them in a file in the folder Data files

  compute_int_plasma

  set_up_gp_arrays


  psi_tab   = makelogtab(psi_min, psi_max, Npsi)
  phi_tab   = makelogtab(phi_min, phi_max, Nphi)
  Omega_tab = makelogtab(Omega_min, Omega_max, NOmega)

  gp_pos = dblarr(Npsi, NOmega)
  gp_neg = dblarr(Npsi, NOmega)
  gp_neutral = dblarr(Nphi, NOmega)

; --- charged grains ---

  for ipsi = 0, Npsi - 1 do begin
     psi = psi_tab[ipsi]
     for iOmega = 0, NOmega - 1 do begin
        Omega = Omega_tab[iOmega]
        gp_pos[ipsi, iOmega] = little_gp_charged(psi, Omega)
        gp_neg[ipsi, iOmega] = little_gp_charged(-psi, Omega)

     endfor
  endfor

  openw, 1, SpDust_data_dir + 'gp_pos_'+ strtrim(Npsi, 2)+ 'psi_' + strtrim(NOmega, 2) + 'Omega' 
  printf, 1, gp_pos
  close, 1

  openw, 1, SpDust_data_dir + 'gp_neg_' + strtrim(Npsi, 2)+ 'psi_' + strtrim(NOmega, 2) + 'Omega' 
  printf, 1, gp_neg
  close, 1

  for iphi = 0, Nphi - 1 do begin
     phi = phi_tab[iphi]
     for iOmega = 0, Nomega -1 do begin
        Omega = Omega_tab[iOmega]
        gp_neutral[iphi, iOmega] =  little_gp_neutral(phi, Omega)
     endfor
  endfor

  openw, 1, SpDust_data_dir + 'gp_neutral_' + strtrim(Nphi, 2)+ 'phi_' + strtrim(NOmega, 2) + 'Omega'
  printf, 1, gp_neutral
  close, 1


end

;-----------------------------------------------------------------------------------------

function little_gp_charged_interpol, psi_arr, Omega_arr
common gp_arrays, psi_min, psi_max, Npsi, phi_min, phi_max, Nphi, Omega_min, Omega_max, NOmega, gp_pos, gp_neg, gp_neutral
common warnings, warning_phi_min, warning_phi_max, warning_psi_min, warning_psi_max, warning_Omega_min, warning_Omega_max

; ---> Returns gp(psi_arr, Omega_arr) for given arrays of psi 
; values and Omega_values

; --- finding the index and weight for |psi| ---
; --- psi should always be within the array but just in case...

; Modified February 2010 : replaced the woarnings by common
; variables. Will be stated only once instead of many many times in
; case they are needed.


  psi_tab = makelogtab(psi_min, psi_max, Npsi)
  psi_indices = floor(Npsi *alog(abs(psi_arr)/psi_min)/alog(psi_max/psi_min) - 0.5d)

  ind = where(psi_indices LT 0, count)
  if count NE 0 then begin
;  print, '|psi| is less than psi_min in "little_gp_charged_interpol". You should set psi_min to a smaller value in "set_up_gp_arrays" and run "compute_little_gp".'
     warning_psi_min = 1
     psi_indices[ind] = 0       ; using minimal value
  endif

  ind = where(psi_indices GT Npsi - 2, count)
  if count NE 0 then begin
;  print, '|psi| is greater than psi_max in
;  "little_gp_charged_interpol". You should set psi_max to a greater
;  value in "set_up_gp_arrays" and run "compute_little_gp".'
     warning_psi_max = 1
     psi_indices[ind] = Npsi - 2 ; using maximal value
  endif

  psi_coeff = 1d - Npsi * alog(abs(psi_arr)/psi_tab[psi_indices])/alog(psi_max/psi_min)
  ind = where(psi_coeff GT 1d, count) ; for the eventual out-of-the-array values
  if count NE 0 then begin
     psi_coeff[ind] = 1d
  endif
  ind = where(psi_coeff LT 0d, count) ; for the eventual out-of-the-array values
  if count NE 0 then begin
     psi_coeff[ind] = 0d
  endif

; --- finding the indices and weights for Omega ---

  Omega_tab = makelogtab(Omega_min, Omega_max, NOmega)
  Omega_indices = floor(NOmega *alog(Omega_arr/Omega_min)/alog(Omega_max/Omega_min) - 0.5d)

  ind = where(Omega_indices LT 0, count)
  if count NE 0 then begin
;  print, 'Omega is less than Omega_min in "little_gp_charged_interpol" (plasmadrag.pro). You should set Omega_min to a smaller value in "set_up_gp_arrays" and run "compute_little_gp".'
     warning_Omega_min = 1
     Omega_indices[ind] = 0
  endif

  ind = where(Omega_indices GT NOmega - 2, count)
  if count NE 0 then begin
; print, 'Omega is greater than Omega_max in "little_gp_charged_interpol"(plasmadrag.pro). You should set Omega_max to a greater value in "set_up_gp_arrays" and run "compute_little_gp".'
     warning_Omega_max = 1
     Omega_indices[ind] = NOmega - 2
  endif

  Omega_coeff = 1d - NOmega * alog(Omega_arr/Omega_tab[Omega_indices])/alog(Omega_max/Omega_min)
  ind = where(Omega_coeff GT 1d, count)
  if count NE 0 then begin
     Omega_coeff[ind] = 1d
  endif
  ind = where(Omega_coeff LT 0d, count)
  if count NE 0 then begin
     Omega_coeff[ind] = 0d
  endif

; --- now using the calculated gp_pos, gp_neg ---

;   Omega_j+1 +       +(3)  +(4)
;             |
;             |           
;   Omega_j   +       +(1)  +(2)         
;             |     
;             |- - - -+ - - + - - -
;                   psi_i   psi_i+1
; i = psi_indices
; j = Omega_indices
; |psi| = psi_coeff * psi_i + (1 - psi_coeff) * psi_i+1
; Omega =  Omega_coeff * Omega_j + (1 - Omega_coeff) * Omega_j+1

; --- --- --- --- --- --- --- ---

  Npsi_arr   = N_elements(psi_arr)
  NOmega_arr = N_elements(Omega_arr)
  gp_charged = dblarr(Npsi_arr, NOmega_arr) 

  psi_indices   = matrix_multiply(psi_indices, 1 + intarr(NOmega_arr), /btranspose)
  Omega_indices = matrix_multiply(1 + intarr(Npsi_arr), Omega_indices, /btranspose)

  psi_coeff   = matrix_multiply(psi_coeff, 1d + dblarr(NOmega_arr), /btranspose)
  Omega_coeff = matrix_multiply(1d + dblarr(Npsi_arr), Omega_coeff, /btranspose)


; --- positively charged grains ---

  ind_pos = where(psi_arr GT 0 , Npos)
  if Npos NE 0 then begin
     psi_pos_indices     = psi_indices[ind_pos, *]
     Omega_pos_indices   = Omega_indices[ind_pos, *]
     psi_pos_coeff       = psi_coeff[ind_pos, *]
     Omega_pos_coeff     = Omega_coeff[ind_pos, *] 
     
     gp_charged[ind_pos, *] = exp( psi_pos_coeff * (Omega_pos_coeff *alog(gp_pos[psi_pos_indices, Omega_pos_indices]) $
                                                    +(1d - Omega_pos_coeff) * alog(gp_pos[psi_pos_indices, Omega_pos_indices + 1])) $
                                   + (1d - psi_pos_coeff) * (Omega_pos_coeff *alog(gp_pos[psi_pos_indices +1, Omega_pos_indices]) $
                                                             +(1d - Omega_pos_coeff) *alog(gp_pos[psi_pos_indices + 1, Omega_pos_indices + 1])) )        
  endif

; --- negatively charged grains ---

  ind_neg = where(psi_arr LT 0 , Nneg)
  if Nneg NE 0 then begin
     psi_neg_indices     = psi_indices[ind_neg, *]
     Omega_neg_indices   = Omega_indices[ind_neg, *]
     psi_neg_coeff       = psi_coeff[ind_neg, *]
     Omega_neg_coeff     = Omega_coeff[ind_neg,*] 
     
     gp_charged[ind_neg, *] = exp( psi_neg_coeff * (Omega_neg_coeff *alog(gp_neg[psi_neg_indices, Omega_neg_indices]) $
                                                    +(1d - Omega_neg_coeff) *alog(gp_neg[psi_neg_indices, Omega_neg_indices + 1])) $
                                   + (1d - psi_neg_coeff) * (Omega_neg_coeff *alog(gp_neg[psi_neg_indices +1, Omega_neg_indices]) $
                                                             +(1d - Omega_neg_coeff) *alog(gp_neg[psi_neg_indices + 1, Omega_neg_indices + 1])))        
  endif

  return, gp_charged

end 

;-----------------------------------------------------------------------------------------

function little_gp_neutral_interpol, phi, Omega_arr
common gp_arrays, psi_min, psi_max, Npsi, phi_min, phi_max, Nphi, Omega_min, Omega_max, NOmega, gp_pos, gp_neg, gp_neutral
common warnings, warning_phi_min, warning_phi_max, warning_psi_min, warning_psi_max, warning_Omega_min, warning_Omega_max

; ---> Returns gp(phi, Omega_arr) for a neutral grain, for
; one value of phi, and an array of Omega

; --- finding the index and weight for phi ---
; --- phi should always be within the array but just in case...

; Modified February 2010 : replaced the woarnings by common
; variables. Will be stated only once instead of many many times in
; case they are needed.

  phi_tab     = makelogtab(phi_min, phi_max, Nphi)
  phi_index   = floor(Nphi *alog(phi/phi_min)/alog(phi_max/phi_min) - 0.5d)

  if phi_index LT 0 then begin
;  print, 'phi is less than phi_min in "little_gp_neutral_interpol". You should set phi_min to a smaller value in "set_up_gp_arrays" and run "compute_little_gp".'
     warning_phi_min = 1
     phi_index = 0              ; using minimal value
  endif

  if phi_index GT Nphi - 2 then begin
;  print, 'phi is greater than phi_max in "little_gp_neutral_interpol". You should set phi_max to a greater value in "set_up_gp_arrays" and run "compute_little_gp".'
     warning_phi_max = 1
     phi_index = Nphi - 2       ; using minimal value
  endif

  phi_coeff = 1d - Nphi * alog(phi/phi_tab[phi_index])/alog(phi_max/phi_min)
  phi_coeff = max(0d, phi_coeff)
  phi_coeff = min(1d, phi_coeff)

  gp_neu_phi = reform(phi_coeff * gp_neutral[phi_index, * ] + (1d - phi_coeff) * gp_neutral[phi_index + 1, * ])

; --- finding the indices and weights for Omega ---

  Omega_tab = makelogtab(Omega_min, Omega_max, NOmega)
  Omega_indices = floor(NOmega *alog(Omega_arr/Omega_min)/alog(Omega_max/Omega_min) - 0.5d)

  ind = where(Omega_indices LT 0, count)
  if count NE 0 then begin
;  print, 'Omega is less than Omega_min in "little_gp_neutral_interpol". You should set Omega_min to a smaller value in "set_up_gp_arrays" and run "compute_little_gp".'
     warning_Omega_min = 1
     Omega_indices[ind] = 0
  endif

  ind = where(Omega_indices GT NOmega - 2, count)
  if count NE 0 then begin
     warning_Omega_max = 1
     Omega_indices[ind] = NOmega - 2
  endif

  Omega_coeff = 1d - NOmega * alog(Omega_arr/Omega_tab[Omega_indices])/alog(Omega_max/Omega_min)
  ind = where(Omega_coeff GT 1d, count)
  if count NE 0 then begin
     Omega_coeff[ind] = 1d
  endif
  ind = where(Omega_coeff LT 0d, count)
  if count NE 0 then begin
     Omega_coeff[ind] = 0d
  endif


  return, Omega_coeff * gp_neu_phi[Omega_indices] + (1d - Omega_coeff) * gp_neu_phi[Omega_indices + 1]

end

;----------------------------------------------------------------------------------------

function Gp_sphere_per_mu2_averaged, env, a, fZ, omega
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Returns G_p^(AHD09)(a, omega)/mu_ip^2 averaged over grain charges.
; Note that in principle < Gp > = < mu_ip(Zg)^2 * G_per_mu2>
; Because of the weak dependence on charge in mu_perp, we will take
; <Gp> \approx <mu_ip^2> * <G_per_mu2>  (averages over charge)

; Returned value : array of the same dimensions as omega.

; --- Added November 2009 ---

    T   = env.T
    acx = acx(a)
    xh  = env.xh
    xC  = env.xC

    Nomega = N_elements(omega)

    Zg_arr = reform(fZ[0,*])
    fZ_arr = reform(fZ[1,*])

    little_gp_H   = dblarr(Nomega)
    little_gp_C   = dblarr(Nomega)

    Omega_H = sqrt(mp/(2d *k *T)) * acx * omega
    Omega_C = sqrt(12d) * Omega_H

; --- Neutral grain contribution ---

      phi = sqrt(2d/(acx *k *T)) *q 
      little_gp_H = fZ[1,0] * little_gp_neutral_interpol(phi, Omega_H) 
      little_gp_C = fZ[1,0] * little_gp_neutral_interpol(phi, Omega_C)

; --- Charged grain contribution ---

    ind_charged = where(Zg_arr NE 0, N_charged)
    if N_charged NE 0 then begin
      psi_arr = Zg_arr[ind_charged] * q^2/(acx *k *T)
      little_gp_H = little_gp_H + matrix_multiply(little_gp_charged_interpol(psi_arr, Omega_H), fZ_arr[ind_charged], /atranspose)
      little_gp_C = little_gp_C + matrix_multiply(little_gp_charged_interpol(psi_arr, Omega_C), fZ_arr[ind_charged], /atranspose)
    endif

    Gp_over_mu2 = (q/(acx^2 *k *T))^2 *(xh *little_gp_H + xC *sqrt(12d) *little_gp_C)

    return, Gp_over_mu2

end

;---------------------------------------------------------------------------------------

function FGp_averaged, env, a, fZ, omega, mu_ip, mu_op, tumbling = tumbling

; Returns a structure {Fp, Gp}. 
; Each element is an array with dimensions [Nomega, Nmu].
; (mu_para and mu_perp must be of the same dimension Nmu).


; --- Added November 2009 ---


    if keyword_set(tumbling) then begin
; --- Disklike tumbling grain. Using the 2 point Gaussian integration method
       Gp_op = 2d/3d *Gp_sphere_per_mu2_averaged(env, a, fZ, 2d *omega)   
       Fp_op = 2d *Gp_op
       omegaG_plus  = (3d + sqrt(3d/5d))/2d * omega
       omegaG_minus = (3d - sqrt(3d/5d))/2d * omega
       Gp_ip = (Gp_sphere_per_mu2_averaged(env, a, fZ, omegaG_plus) $
                + Gp_sphere_per_mu2_averaged(env, a, fZ, omegaG_minus))/3d
  
       omegaF_plus  = (8d + sqrt(13d/3d))/5d * omega
       omegaF_minus = (8d - sqrt(13d/3d))/5d * omega 
       Fp_ip = 0.5d *(Gp_sphere_per_mu2_averaged(env, a, fZ, omegaF_plus) $
                      + Gp_sphere_per_mu2_averaged(env, a, fZ, omegaF_minus))
       
       Fp = matrix_multiply(Fp_ip, mu_ip^2, /btranspose) $
          + matrix_multiply(Fp_op, mu_op^2, /btranspose)
       Gp = matrix_multiply(Gp_ip, mu_ip^2, /btranspose) $
          + matrix_multiply(Gp_op, mu_op^2, /btranspose) 

    endif else begin
; --- Standard spherical grain with K = J
       Gp = Gp_sphere_per_mu2_averaged(env, a, fZ, omega)
       Gp = matrix_multiply(Gp, mu_ip^2, /btranspose)
       Fp = Gp
    endelse
    
     
    return, {Fp: Fp, Gp: Gp}

end

;-----------------------------------------------------------------------------------------

