; ================================ SPDUST.2 =====================================;
;                                                                                ;
;                               CHARGE_DIST.PRO                                  ;
;                                                                                ; 
; First computes the photoemission rate according to Weingartner & Draine,       ;   
; 2001b, ApJ 134, 263. All equation numbers refer to that paper,                 ;
; unless explicitly mentioned. Then compute collisional charging rates           ;
; according to Draine & Sutin, 1987, ApJ, 320, 803.                              ;
; Finally, gets the grain charge distribution function for given grain size      ;
; and environmental conditions.                                                  ;
;                                                                                ;      
; Author : Yacine Ali-Haimoud        yacine@tapir.caltech.edu                    ; 
;                                                                                ;
; Revision history: Written December 2008                                        ;
;                   Revised March 2010 (Minor change in function charge_dist)    ;
;================================================================================;



function IPv, Z, a              ; ---> Equation (2), in eV
common cgsconst, pi, c, q, k, mp, me, h, debye, eV
common paramjpe, W

    return, W + ((Z + 0.5d)*q^2/a + (Z + 2d)*q^2/a * (0.3d-8)/a)/eV

end

;---------------------------------------------------------------------------------------

function E_min, Z, a            ; ---> Equation (7), in eV
common cgsconst, pi, c, q, k, mp, me, h, debye, eV
common paramjpe, W

  if(Z LT -1) then begin
     return, -(Z + 1d)*q^2/a /(1d + (27d-8/a)^0.75d )/eV
  endif else begin
     return, 0d
  endelse

end

;---------------------------------------------------------------------------------------

function hnu_pet, Z, a          ; ---> Equation (6), in eV

  if(Z GE -1) then begin 
     return, max([0d, IPv(Z,a)])
  endif else begin
     return, max([0d, IPv(Z,a) + E_min(Z,a)])
  endelse

end 

;---------------------------------------------------------------------------------------

function theta, hnu_tab, Z, a   ; ---> Equation (9), in eV
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

  hnu_pet = hnu_pet(Z,a)

  if(Z GE 0) then begin
     return, hnu_tab - hnu_pet + (Z + 1d)*q^2/a /eV
  endif else begin
     return, hnu_tab - hnu_pet
  endelse

end

;---------------------------------------------------------------------------------------

function y2, hnu_tab, Z, a      ; ---> Equation (11) and preceding paragraph
  common cgsconst, pi, c, q, k, mp, me, h, debye, eV

  Nnu = N_elements(hnu_tab)
  y2_tab = dblarr(Nnu)

  if (Z GE 0) then begin
     E_low  = -(Z + 1d)*q^2/a /eV
     E_high = hnu_tab - hnu_pet(Z,a)
     ind = where(E_high GT max([0d, E_low]), count)
     if count NE 0  then begin
        E_high = E_high[ind] 
        y2_tab[ind] = E_high^2*(E_high - 3d*E_low)/(E_high - E_low)^3
     endif 
     
  endif else begin
     ind = where(hnu_tab GT hnu_pet(Z,a), count)
     if count NE 0 then begin
        y2_tab[ind] =  1d
     endif 
  endelse

  return, y2_tab

end

;--------------------------------------------------------------------------------------- 

pro l_a
common cgsconst,     pi, c, q, k, mp, me, h, debye, eV
common path,         SpDust_data_dir
common refr_indices, hnu_tab, la_tab

; ---> Stores la(hnu) needed for equation (15) of WD01b.    
; Uses tabulated values for Im(m_para) and
; Im(m_perp), taken from Draine website and given for a = 100A, T = 20K.
; We assume the changes due to grain size and temperature are negligible.

  fileperp   = SpDust_data_dir + 'perpendicular.out'
  filepara   = SpDust_data_dir + 'parallel.out'

  readcol, fileperp, lamb_perp, Imm_perp, comment = ';', format = 'D,X,X,X,D', /silent 
  readcol, filepara, lamb_para, Imm_para, comment = ';', format = 'D,X,X,X,D', /silent

; There are more wavelength values given for m_para to account for line effects,
; so we interpolate Im(m_perp) at those values.

  Imm_perp = exp( interpol(alog(Imm_perp), alog(lamb_perp), alog(lamb_para)) )

  lamb_tab = lamb_para*1d-4     ; in cm
  hnu_tab  = h*c/(lamb_tab)/eV

  la_tab = 3d*lamb_tab/(4d*pi)/(2d*Imm_perp + Imm_para) ; WD01 (15) 

  print, 'l_a computed'

end

;--------------------------------------------------------------------------------------- 

function fy1, x                 ; ---> Takes care of small beta case in eq (13)

  xmin = 1d-5
  ind = where(x LT xmin, count, complement = rest)
  result = dblarr(N_elements(x))

  if(count NE 0) then begin
     xsmall = x[ind]
     result[ind] = xsmall/3d 
  endif

  if(count NE N_elements(x)) then begin
     x = x[rest]
     result[rest] = (x^2- 2d*x + 2d - 2d*exp(-x))/x^2
  endif

  return, result

end

;--------------------------------------------------------------------------------------- 

function y1, a
common cgsconst,     pi, c, q, k, mp, me, h, debye, eV
common refr_indices, hnu_tab, la_tab

; ---> Returns y1(a, hnu) for the values stored in hnutab. Equation (13) in WD01b

  l_e = 10d-8                   ; electron escape length in cm

  beta  = a/la_tab
  alpha = a/la_tab + a/l_e

  return, fy1(alpha)/fy1(beta)

end

;--------------------------------------------------------------------------------------- 

function y0, theta              ; ---> Equation (16) in WD01b
common paramjpe, W

  return, 9d-3*(theta/W)^5/(1d + 3.7d-2*(theta/W)^5)

end

;--------------------------------------------------------------------------------------- 

function Y, Z, a                ; ---> Equation (12) for hnu values in hnu_tab
common refr_indices, hnu_tab, la_tab

  N_hnu = N_elements(hnu_tab)
  Y = dblarr(N_hnu)

  y0 = y0(theta(hnu_tab,Z,a))
  y1 = y1(a)

  y0y1 = y0*y1
  ind = where(y0y1 GT 1d, count)
  if count NE 0 then begin
     y0y1[ind] = 1d
  endif

  Y = y2(hnu_tab, Z, a) * y0y1

  return, Y

end

;--------------------------------------------------------------------------------------- 

function Planck_B, nu, T        ; ---> Planck function for nu in Hz
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

  return, 2d*h*nu^3/c^2 /(exp(h*nu/(k*T))- 1d)

end

;---------------------------------------------------------------------------------------

function nu_uisrf, hnu_tab
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Average interstellar radiation field spectrum nu*u_nu
; (Mathis, Mezger & Panagia 1983) see eq (31) of WD01b
; IMPORTANT NOTE : the ambiant radiation field is always supposed to
; be a multiple of u_ISRF. Of course you may change that but you will
; need to re-run some calculations for the Infrared damping and
; excitation coefficients (see readme file). 

  uisrf_tab = dblarr(N_elements(hnu_tab))

  ind1 = where((hnu_tab LT 13.6) and (hnu_tab GT 11.2), count1)
  if (count1 NE 0) then begin
     uisrf_tab[ind1] = 3.328d-9 * hnu_tab[ind1]^(-4.4172d)
  endif

  ind2 = where((hnu_tab LT 11.2) and (hnu_tab GT 9.26), count2)
  if(count2 NE 0) then begin
     uisrf_tab[ind2] = 8.463d-13/hnu_tab[ind2]
  endif

  ind3 = where((hnu_tab LT 9.26) and (hnu_tab GT 5.04), count3)
  if(count3 NE 0) then begin
     uisrf_tab[ind3] =  2.055d-14*hnu_tab[ind3]^0.6678d
  endif

  ind4 = where(hnu_tab LT 5.04, count4)
  if(count4 NE 0) then begin
     nu = hnu_tab[ind4] * eV /h
     uisrf_tab[ind4] =  4*pi*nu/c*(1.d-14*Planck_B(nu, 7500.d) + 1.65d-13*Planck_B(nu, 4000.d)+ 4d-13*Planck_B(nu, 3000.d))
  endif

  return, uisrf_tab

end

;---------------------------------------------------------------------------------------

pro readPAH
common cgsconst, pi, c, q, k, mp, me, h, debye, eV
common path,     SpDust_data_dir
common Qabstabs, Nrad, Nwav, a_tab, Qabs_hnu_tab, Q_abs_ion, Q_abs_neu

; ---> Reads absorption efficiencies in files given in Bruce T. Draine
; website. Thanks to Nate for his help.

  PAHion = SpDust_data_dir + 'PAHion_read.out'
  PAHneu = SpDust_data_dir + 'PAHneu_read.out' 

; --- Number of grain radii and wavelengths
  Nrad = 30
  Nwav = 1201L
; --- --- --- --- --- --- --- --- --- --- -

  readcol, PAHion, C1, C2_ion, format = 'D, X, D, X, X', comment = ';', /silent
  readcol, PAHneu,  C2_neu, format = 'X, X, D, X, X', comment = ';', /silent

  a_tab      = dblarr(Nrad)
  lambda_tab = dblarr(Nwav)
  Q_abs_neu  = dblarr(Nrad, Nwav)
  Q_abs_ion  = dblarr(Nrad, Nwav)

  for i = 0, Nrad - 1 do begin
     a_tab[i] = C1[(Nwav+1)*i]*1d-4
     for j = 0, Nwav - 1 do begin
        lambda_tab[j]  = C1[(Nwav+1)*i + j+1]
        Q_abs_neu[i,j] = C2_neu[(Nwav+1)*i + j+1]
        Q_abs_ion[i,j] = C2_ion[(Nwav+1)*i + j+1]
     endfor
  endfor

  Qabs_hnu_tab = h * c/(lambda_tab*1d-4) /eV

  print, 'readPAH computed'

end

;---------------------------------------------------------------------------------------

function Qabs, a, Z, hnu_tab
common cgsconst, pi, c, q, k, mp, me, h, debye, eV
common Qabstabs, Nrad, Nwav, a_tab, Qabs_hnu_tab, Q_abs_ion, Q_abs_neu

; ---> Interpolates the tables read in readPAH, returns the absorption efficiencies

; --- Max wavelength given by DL98 --- 
  hnu_min = min(Qabs_hnu_tab)
; --- --- --- --- --- --- --- --- ---

  N_nu = N_elements(hnu_tab)
  Qtab = dblarr(N_nu)

  ind  = where(hnu_tab GT hnu_min, count)
  if (count NE 0) then begin
     if (Z EQ 0) then begin    
        Qtab[ind] = reform(log_biinterplog(Q_abs_neu, a_tab, Qabs_hnu_tab, a, hnu_tab[ind]))   
     endif else begin
        Qtab[ind] = reform(log_biinterplog(Q_abs_ion, a_tab, Qabs_hnu_tab, a, hnu_tab[ind])) 
     endelse
  endif

; At long wavelengths, Qabs \propto nu^2

  ind = where(hnu_tab LE hnu_min, count)      

  if (count NE 0) then begin
     if (Z EQ 0) then begin    
        Qtab[ind] = (hnu_tab[ind]/hnu_min)^2 * (log_biinterplog(Q_abs_neu, a_tab, Qabs_hnu_tab, a, hnu_min))[0,0]
     endif else begin
        Qtab[ind] = (hnu_tab[ind]/hnu_min)^2 * (log_biinterplog(Q_abs_ion, a_tab, Qabs_hnu_tab, a, hnu_min))[0,0]
     endelse
  endif

  return, Qtab

end

;---------------------------------------------------------------------------------------

function EA, Z, a               ; ---> Equation (4), in eV  for Z < 0
common cgsconst, pi, c, q, k, mp, me, h, debye, eV
common paramjpe, W

  return, W + ( (Z - 0.5d)*q^2/a - q^2/a * 4d-8/(a + 7d-8) )/eV

end

;---------------------------------------------------------------------------------------

function hnu_pdt, Z, a          ; ---> Equation (18), in eV  for Z < 0

  return, max([0d, EA(Z+1,a) + E_min(Z,a)])

end

;---------------------------------------------------------------------------------------

function sigma_pdt, hnu_tab, Z, a ; ---> Equation (20), in cm^2, for Z < 0

  DeltaE = 3d                   ; in eV
  x      = (hnu_tab - hnu_pdt(Z,a))/DeltaE
  sigma  = dblarr(N_elements(hnu_tab))

  if( (Z GE 0) ) then begin
     return, sigma
  endif 

  ind = where(x GT 0, count)
  if (count NE 0) then begin 
     x = x[ind]
     sigma[ind] =  1.2d-17 * abs(Z) * x/(1.d + x^2/3.d)^2
  endif

  return, sigma

end

;---------------------------------------------------------------------------------------

function first_term, Z, a  
common cgsconst, pi, c, q, k, mp, me, h, debye, eV
common Qabstabs, Nrad, Nwav, a_tab, Qabs_hnu_tab, Q_abs_ion, Q_abs_neu
common refr_indices, hnu_tab, la_tab

; ---> First term in eq (25) for the standard interstellar radiation field (31)

  hnu_min = max([1d-5, hnu_pet(Z,a)])
  hnu_max = 13.6

;--Parameter------;
  Nnu     = 500                 
;-----------------;

  if ( hnu_min GT hnu_max) then begin
     return, 0d
  endif

  hnu    = makelogtab(hnu_min, hnu_max, Nnu)
  Dnu_over_nu   = Dx_over_x(hnu_min, hnu_max, Nnu)

  Ytab = interpol(Y(Z,a), alog(hnu_tab), alog(hnu))

  Qtab = Qabs(a, Z, hnu)
  utab = nu_uisrf(hnu)

  return, c*pi*a^2 * Dnu_over_nu * total( Ytab*Qtab* utab/(hnu*eV), /double)

end

;---------------------------------------------------------------------------------------

function second_term, Z,a   
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Second term in eq (25) for the standard interstellar radiation field (31), for Z < 0

  hnu_min = max([1d-5, hnu_pdt(Z,a)]) 
  hnu_max = 13.6d

;--Parameter------;
  Nnu     = 500                 ;
;-----------------;

  if((Z GE 0) or (hnu_pdt(Z,a) GT hnu_max)) then begin 
     return, 0d 
  endif 

  hnu         = makelogtab(hnu_min, hnu_max, Nnu)
  Dnu_over_nu = DX_over_X(hnu_min, hnu_max, Nnu)

  utab = nu_uisrf(hnu)
  sigmatab = sigma_pdt(hnu, Z,a)

  return, c * total( sigmatab* utab/(hnu*eV), /double) * Dnu_over_nu

end

;---------------------------------------------------------------------------------------

function Jpe, Z, a

  return, first_term(Z,a) + second_term(Z,a)

end

;---------------------------------------------------------------------------------------

function Zmax, a                ; ---> Equation (22)
common paramjpe, W

  ;---
  hnu_max = 13.6d               ; !!! change if needed !!!
  ;----
  aA = a/1d-8                   ; a in A   

  return, floor( ((hnu_max-W)/14.4d *aA + 0.5d - 0.3d/aA) / (1d + 0.3d/aA) )

end

;---------------------------------------------------------------------------------------

function Zmin, a                ; ---> Equations (23), (24)

  aA = a/1d-8                   ; a in A   
  U_ait = -(3.9d + 0.12d *aA + 2d/aA)

  return, floor(U_ait/14.4d*aA) + 1

end

;---------------------------------------------------------------------------------------

pro JPEisrf_calc
common paramjpe, W
common jpe_arrays, a_values, Jpe_pos_isrf, Jpe_neg_isrf

; ---> Computes and stores the photoemission rate Jpe(a, Z), for an
; array of grain radii and charges, for tha averge interstellar
; radiation field.

;---- PARAMETERS -----
  W = 4.4                       
  a_min = 3.5d-8                
  a_max = 1d-6                  
  Na = 30                       
;--------------------- 

  physconst
  l_a
  readpah

  a_values = makelogtab(a_min, a_max, Na)

  Z_min = Zmin(a_max)
  Z_max = Zmax(a_max)

  Jpe_neg_isrf = dblarr(Na, abs(Z_min) + 1)
  Jpe_pos_isrf = dblarr(Na, Z_max + 1)

  for i = 0, Na - 1 do begin
     a = a_values[i]
     for j = 0, abs(Zmin(a)) do begin 
        Z = -j
        Jpe_neg_isrf[i,j] = Jpe(Z,a) 
     endfor
     for j = 0, Zmax(a) do begin
        Z = j
        Jpe_pos_isrf[i,j] = Jpe(Z,a)
     endfor
  endfor

  print, 'Jpeisrf computed'

end

;---------------------------------------------------------------------------------------

function JPEisrf, a
common jpe_arrays, a_values, Jpe_pos_isrf, Jpe_neg_isrf

; ---> Interpolates the arrays calculated in JPEisrf_calc
; Returns Jpe(a, Z) where Z runs over the widest range of possible charges.

  Jpepos = dblarr(N_elements(Jpe_pos_isrf[0,*]))
  Jpeneg = dblarr(N_elements(Jpe_neg_isrf[0,*]))

  if (a LE min(a_values)) then begin
     Jpepos = Jpe_pos_isrf[0, *]
     Jpeneg = Jpe_neg_isrf[0, *] 
  endif else begin
     if (a GE max(a_values)) then begin
        N_a = N_elements(a_values)
        Jpepos = Jpe_pos_isrf[N_a -1, *]
        Jpeneg = Jpe_neg_isrf[N_a -1, *] 
     endif else begin
        ia = max(where(a_values LE a))
        alpha = alog(a/a_values[ia])/alog(a_values[ia+1]/a_values[ia])
        Jpepos = (1d - alpha) *Jpe_pos_isrf[ia, *] + alpha * Jpe_pos_isrf[ia+1, *]
        Jpeneg = (1d - alpha) *Jpe_neg_isrf[ia, *] + alpha * Jpe_neg_isrf[ia+1, *]
     endelse
  endelse

  Jpepos = reform(Jpepos)
  Jpeneg = reform(Jpeneg)

  Zmin = Zmin(a)
  Zmax = Zmax(a)

  Jpepos = Jpepos[0:Zmax]
  Jpeneg = Jpeneg[0:-Zmin]

  return, {Jpepos: Jpepos, Jpeneg: Jpeneg}

end

;---------------------------------------------------------------------------------------

function Jtilde, tau, nu       
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Returns the fuction Jtilde defined in Draine & Sutin, 1987, ApJ
; 320, 803. Equation numbers refer to this paper.

  if(nu EQ 0) then begin
     return, 1d + sqrt(pi/(2d*tau)) ; (3.3)
  endif

  if (nu LT 0) then begin
     return, (1d - nu/tau)*(1d + sqrt(2d/(tau - 2d*nu))) ; (3.4)
  endif

  if (nu GT 0) then begin
     ksi = 1d + 1d/sqrt(3d*nu)                                  
     theta = nu/ksi - 0.5d/(ksi^2*(ksi^2 - 1d))                 ; (2.4a)
     if (theta/tau LT 700d) then begin                          ; avoid floating underflow
        return, (1d + (4d*tau + 3d*nu)^(-0.5d))^2 * exp(-theta/tau) ; (3.5)
     endif else begin
        return, 0d
     endelse
  endif

end

;---------------------------------------------------------------------------------------

function se, Z, a              

; ---> Sticking coefficient for electrons. Note that s_i = 1 for ions
; Equations refer to WD01b

  l_e = 1d-7                    ; "electron escape length" (10 A) see below (28)
  Nc  = N_C(a)
  Zmin = Zmin(a)

  if ((Z EQ 0)or((Z LT 0)and(Z GT Zmin))) then begin
     return, 0.5d *(1d - exp(-a/l_e)) /(1d + exp(20d - Nc)) ; (28),(29)
  endif 

  if((Z LT 0) and (Z LE Zmin)) then begin 
     return, 0d                 ; (29)
  endif 

  if(Z GT 0) then begin
     return, 0.5d *(1d - exp(-a/l_e)) ; (30)
  endif

end

;----------------------------------------------------------------------------------------

function J_ion, env, Z, a
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

  nh = env.nh
  T  = env.T
  xh = env.xh
  xC = env.xC

  as  = asurf(a)
  tau = as *k *T/q^2
  nu  = Z                       ; assume singly charged ions

  return, nh *sqrt(8d *k *T/(pi *mp)) * pi *as^2 * (xh + xc/sqrt(12d)) *Jtilde(tau, nu)

end

;-----------------------------------------------------------------------------------------

function J_electron, env, Z, a
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

  nh = env.nh
  T  = env.T
  xh = env.xh
  xC = env.xC

  as  = asurf(a)
  tau = as *k *T/q^2
  nu  = -Z                      ; assume singly charged ions

  return, nh *(xh + xC) *se(Z, a) *sqrt(8d *k *T/(pi *me)) * pi *as^2 *Jtilde(tau, nu)

end

;-----------------------------------------------------------------------------------------

function charge_dist, env, a 
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Computes the charge distribution of grain of radius a in a
; given environment. Returns fZ s.t. fZ[0,*] = charges , fZ[1,*] =
; value of the distribution function.

;--- minimum value of distribution function ---;
  fmin = 1d-10                     
;--- --- --- --- --- --- --- --- --- --- --- --;

  Z_min = Zmin(a)
  Z_max = Zmax(a)

  fneg = dblarr(abs(Z_min) + 1)
  fpos = dblarr(Z_max + 1) 

  Chi = env.Chi        

  Ji_neg  = dblarr(abs(Z_min) + 1)
  Je_neg  = dblarr(abs(Z_min) + 1)
  Ji_pos  = dblarr(Z_max + 1)
  Je_pos  = dblarr(Z_max + 1)
  Jpe_ISRF = JPEisrf(a)


  for j = 0, abs(Z_min) do begin
     Z = -j 
     Ji_neg[j] = J_ion(env, Z, a)
     Je_neg[j] = J_electron(env, Z, a)
  endfor
  Jpe_neg = Chi *Jpe_ISRF.Jpeneg

  for Z = 0, Z_max do begin
     Ji_pos[Z] = J_ion(env, Z, a)
     Je_pos[Z] = J_electron(env, Z, a)
  endfor

  Jpe_pos = Chi *Jpe_ISRF.Jpepos 


; Now we compute the charge distribution function using DL98b (4)

  fpos[0] = 1d
  fneg[0] = 1d

  Z = -1
  while (Z GE Z_min)  do begin
     fneg[-Z] = fneg[-(Z+1)] *Je_neg[-(Z+1)] /(Ji_neg[-Z] + Jpe_neg[-Z])
     Z = Z - 1
     fneg = fneg/total(fneg)
  endwhile

  Z = 1
  while ( Z LE Z_max)  do begin
     fpos[Z] = fpos[Z-1] *(Ji_pos[Z-1] + Jpe_pos[Z-1]) /Je_pos[Z]
     Z = Z + 1
     fpos = fpos/total(fpos)
  endwhile

  if(fneg[0] GT 0d) then begin
     fneg = fpos[0]/fneg[0] * fneg
  endif else begin
     fpos = 0
  endelse

  norm = total(fneg, /double) + total(fpos, /double) - fpos[0]

  fneg = fneg /norm
  fpos = fpos /norm

  Zneg = - indgen(abs(Z_min) + 1)
  Zpos = indgen(Z_max + 1)

; Keep only the values for which f_Z > fmin 
; and make sure that fZ[0,0] = [0, f(0)] (changed Feb 10)

  ind = where(((fneg GT fmin) AND (Zneg NE 0)), NZneg)
  if (NZneg NE 0) then begin
     fneg = fneg[ind]
     Zneg = Zneg[ind]
  endif 

  ind = where(((fpos GT fmin) OR (Zpos EQ 0)), NZpos)
  fpos = fpos[ind]
  Zpos = Zpos[ind]

; Now putting everything in a single array

  fZ    = dblarr(2, NZpos + NZneg)

  fZ_array = fpos
  Z_array  = Zpos 

  if (NZneg NE 0) then begin
     Z_array  = [Z_array, Zneg]
     fZ_array = [fZ_array, fneg]
  endif

; renormalize
  fZ_array = fZ_array/total(fZ_array)

  fZ[0,*] = Z_array
  fZ[1,*] = fZ_array


  return, fZ

end

;---------------------------------------------------------------------------------------
