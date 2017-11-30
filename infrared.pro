;=================================== SPDUST.2 =========================================;
;                                                                                      ;
;                                   INFRARED.PRO                                       ; 
;                                                                                      ;
; Normalized rotational damping and excitation rates through infrared emission.        ;                                        
; Refers to Section 7 of Ali-Haimoud, Hirata & Dickinson, 2009, MNRAS, 395, 1055       ;           
; and Section 6 of Silsbee, Ali-Haimoud & Hirata, 2010.                                ;
; We use Draine & Li, 2001, ApJ 551, 807 (DL01) "thermal continuous"                   ;
; approximation to compute the infrared spectrum.                                      ;                 
;                                                                                      ;          
; Author : Yacine Ali-Haimoud        yacine@tapir.caltech.edu                          ; 
;                                                                                      ;
; Revision history : Written December 2008                                             ;
;                    Revised March 2010 (tumbling grains + correction to G_IR)         ;
;======================================================================================;



function f2, x

; ---> Returns f2(x) given by equation (10) in DL01 (note the typo : should
; be n instead of 1/n). x can be an array. 

  Nx    = N_elements(x)
  Ny    = 500
  y     = maketab(0d, 1d, Ny)
  Dy    = 1d/Ny

  integrand = dblarr(Ny, Nx)
  y_over_x = matrix_multiply(y , 1d/x, /btranspose)
  ind = where(y_over_x LT 500d, count)
  if count NE 0 then begin
     integrand[ind] = 1d/(exp(y_over_x[ind]) - 1d)
  endif
  integrand = matrix_multiply(y^2, 1 + dblarr(Nx), /btranspose)*integrand

  return, 2d * Dy * total(integrand, 1, /double)

end

;------------------------------------------------------------------------------------------------------------

function Energy_modes, a  
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Returns E_j for C-C out-of-plane and in-plane modes, and C-H
; bending modes, for a grain of radius a.

; Debye temperatures for out-of-plane and in-plane C-C modes
  Theta_op = 863d   
  Theta_ip = 2504d 

  Nc = N_C(a)
  Nm = (Nc - 2)*[1, 2]          ; number of o-p and i-p modes. needed for beta2

  if Nc LE 54 then begin
     beta2_tab = [0, 0]
  endif else begin 
     if Nc LE 102 then begin
        beta2_tab = (Nc - 54d)/52d/(2d*Nm - 1d)
     endif else begin
        beta2_tab = ((Nc - 2d)/52d * (102d/Nc)^(2d/3d) - 1d)/(2d*Nm - 1d)
     endelse
  endelse


; out-of-plane modes
  beta2 = beta2_tab[0]
  Nm    = Nc - 2
  NewNm = Nm
  Eop_tab = dblarr(NewNm + 1)    
  delta_tab = 0.5 + dblarr(NewNm + 1) ; eq (5)
  delta_tab[2] = 1                    ; eq (6)
  delta_tab[3] = 1                    ; eq (6)  
  Eop_tab[1:NewNm] = k *Theta_op *sqrt( (1d - beta2)/Nm *(1 + dindgen(NewNm) - delta_tab[1:NewNm]) + beta2 )

; in-plane modes
  beta2 = beta2_tab[1]
  Nm    = 2*(Nc - 2)
  NewNm = Nm
  Eip_tab = dblarr(NewNm + 1)      
  delta_tab = 0.5 + dblarr(NewNm + 1)
  delta_tab[2] = 1
  delta_tab[3] = 1
  Eip_tab[1:NewNm] = k *Theta_ip *sqrt( (1d - beta2)/Nm *(1 + dindgen(NewNm) - delta_tab[1:NewNm]) + beta2 )

; C-H stretching modes
  lambda_inverse = [886d, 1161d, 3030d]

  return, {op: Eop_tab[1:*], ip: Eip_tab[1:*], CH: h*c*lambda_inverse}

end

;------------------------------------------------------------------------------------------------------------

function EPAH, a, T, Energy_modes
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Returns the average energy E_PAH (a, T) in the thermal
; approximation, using the exact spectrum for NC < cutoff (equation (2))
; and the large mode limit for NC > cutoff (equation (33))
; Energy_modes = Energy_modes(a) are given to gain speed

  cutoff = 6d4
  NT     = N_elements(T)
  result = dblarr(NT)

; Debye temperatures for out-of-plane and in-plane C-C modes
  Theta_op = 863d   
  Theta_ip = 2504d 

; Number of carbon atoms
  Nc = N_C(a)

; Number of Hydrogen atoms (equation (8))
  NH = N_H(a)

; Energy of C-H bending and stretching modes

;Energy_modes = Energy_modes(a)
  E_Hmodes = dblarr(N_elements( Energy_modes.CH), NT)
  E_over_kT =  matrix_multiply( Energy_modes.CH, 1d/(k*T), /btranspose)

  ind = where(E_over_kT LT 600d, count) ; avoid floating overflow
  if (count NE 0) then begin
     E_Hmodes[ind] = 1d/(exp(E_over_kT[ind]) - 1d)
  endif
  E_Hmodes = matrix_multiply( Energy_modes.CH, 1 + dblarr(NT), /btranspose) * E_Hmodes

  E_Hmodes = NH * total(E_Hmodes, 1, /double)

  if(Nc GT cutoff) then begin   ; continuous limit for large NC, equation (33)
     return,(Nc - 2d)*k *(Theta_op *f2(T/Theta_op) + 2d * Theta_ip*f2(T/Theta_ip)) $
            +  E_Hmodes
  endif else begin              ; exact mode spectrum for small NC
     
     
                                ; op modes
     Eop_tab = Energy_modes.op
     Eop_temp = dblarr(N_elements(Eop_tab), NT) 

     E_over_kT = matrix_multiply(Eop_tab, 1d/(k*T), /btranspose)
     ind = where(E_over_kT LT 600d, count) ; avoid floating overflow
     if (count NE 0) then begin
        Eop_temp[ind] = 1d/(exp(E_over_kT[ind]) - 1d)
     endif

     Eop_temp = matrix_multiply(Eop_tab, 1 + dblarr(NT), /btranspose) * Eop_temp
     Eop_bar = total(Eop_temp, 1, /double)

                                ; ip modes
     Eip_tab = Energy_modes.ip
     Eip_temp = dblarr(N_elements(Eip_tab), NT)

     E_over_kT = matrix_multiply(Eip_tab, 1d/(k*T), /btranspose)
     ind = where(E_over_kT LT 600d, count) ; avoid floating overflow
     if (count NE 0) then begin
        Eip_temp[ind] = 1d/(exp(E_over_kT[ind]) - 1d)
     endif
     Eip_temp = matrix_multiply(Eip_tab, 1 + dblarr(NT), /btranspose) * Eip_temp
     Eip_bar = total(Eip_temp, 1, /double)
     
     return,   Eop_bar + Eip_bar + E_Hmodes

  endelse

end

;---------------------------------------------------------------------------------------------------------------

function Temp, a, Energy, Energy_modes
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Solve for T s.t. EPAH(a, T) = Energy. 
; Energy_modes = Energy_modes(a) are given to gain speed

  Tmin  = 1d
  Tmax  = 1d4
  NT    = 100
  T_tab = makelogtab(Tmin, Tmax, NT)
  E_tab = EPAH(a, T_tab, Energy_modes)

  Temperature = exp(interpol(alog(T_tab), alog(E_tab), alog(Energy)))
;energy_modes = Energy_modes(a)
  modes = [energy_modes.op, energy_modes.ip, energy_modes.CH]

  hbar_omega1 = min(modes, index)
  modes[index] = 2*max(modes)
  for j = 2, 20 do begin        ; find hbar_omega_20
     hbar_omega20 = min(modes, index)
     modes[index] = 2*max(modes)
  endfor

  ind = where(Energy LE hbar_omega20, count)
  if count NE 0 then begin
     Temperature[ind] = hbar_omega1/(k*alog(2d))
  endif

  return, Temperature

end

;---------------------------------------------------------------------------------------------------------------

function Energy_bins, a, Energy_modes, M, Energy_max     
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Returns an array E_tab containing the central, minimum and maximum
; values of M energy bins
; Energy_modes = Energy_modes(a) are given to gain speed.
; Retuns an array of 3* (M+1) elements : Ebin, Emin, Emax 
; Energy_max is the highest energy bin considered

  jmax = 10                     ; the ten first energy bins are calculated as in DL01

  modes = [ (Energy_modes.op)[0:jmax], (Energy_modes.ip)[0:jmax], Energy_modes.CH]

  hbar_omega = dblarr(2*jmax)

  for j = 1, 2*jmax-1 do begin
     hbar_omega[j] = min(modes, index)
     modes[index] = 2*max(modes)
  endfor 

; now the ten first energy bins according to appendix B of DL01

  E_tab       = dblarr(3, M + 1) ; center, min, max, of the energy bins

  E_tab[1, 1] = 1.5d * hbar_omega[1] - 0.5d * hbar_omega[2]
  for j = 1, 2 do begin
     E_tab[0,j]   = hbar_omega[j] 
     E_tab[1,j+1] = 0.5d *(hbar_omega[j] + hbar_omega[j+1])
     E_tab[2,j]   = E_tab[1,j+1]
  endfor
  for j = 3, jmax do begin
     E_tab[2,j] = 0.5d *(hbar_omega[2*j - 2] + hbar_omega[2*j - 1])
     E_tab[0,j] = 0.5d *(hbar_omega[2*j - 3] + hbar_omega[2*j - 2])  
     E_tab[1,j+1] = E_tab[2,j] 
  endfor

; next bins ( up to K ) are linearly spaced

  DeltaE   = E_tab[2, jmax] - E_tab[1, jmax] 

  K_tab = 12 + indgen(M - 11) 

  E_M_max = exp((M - K_tab + 1) *alog(E_tab[2, jmax] + (K_tab - 10) *DeltaE) - (M- K_tab) * alog(E_tab[2, jmax] + (K_tab - 11) *DeltaE))
; Highest energy bin. This is a decreasing function of K 

  if E_M_max[0] LT Energy_max then begin
     print, strtrim(M, 2)+ ' energy bins are not enough, using '+ strtrim(M+50, 2)+ ' bins instead (function "Energy_bins" in infrared.pro)'
     return, Energy_bins(a, Energy_modes, M + 50, Energy_max) 
  endif else begin
     index = max(where(E_M_max GE Energy_max))
  endelse

  BigK = K_tab[index] 

  E_tab[2,jmax+1:BigK] = E_tab[2, jmax] + (dindgen(BigK + 1)- jmax)[jmax+1:BigK] *DeltaE
  E_tab[1,jmax+2:BigK] = E_tab[2, jmax+1:BigK-1]
  E_tab[0,jmax+1:BigK] = 0.5d *(E_tab[1, jmax+1:BigK] + E_tab[2, jmax+1:BigK]) 
  if (BigK LT M) then begin     ; the last bins are logarithmically spaced
     E_tab[2,BigK+1:M] = E_tab[2, BigK] * exp((dindgen(M+1)- BigK)[BigK+1:M] * alog(E_tab[2,BigK]/E_tab[2,BigK-1]) )
     E_tab[1,BigK+1:M] = E_tab[2, BigK :M-1]
     E_tab[0,BigK+1:M] = sqrt(E_tab[2,BigK+1:M]*E_tab[1,BigK+1:M])
  endif

  return, E_tab

end

;---------------------------------------------------------------------------------------------------------------

function Btilde_ISRF_func, a, Z, Energy_modes, M, Energy_max   
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Computes the matrix Bji given in equation (53) for the average interstellar
; radiation field ISRF, for a grain of radius a and charge Z, using
; the thermal continuous approximation.

  Energies       = Energy_bins(a, Energy_modes, M, Energy_max)

  E_tab          = Energies[0,*]       ; energy bins
  M              = N_elements(E_tab) - 1 ; M initial may have been modified in energy_bins

  E_max          = Energies[2,*]
  E_min          = Energies[1,*]
  T_tab          = [0, Temp(a, E_tab[1:*], Energy_modes)] ; temperatures for the thermal approximation

  NEarr          = 500          ; number of elements in the integral over energies
  Energy_min     = h*c/(1d2)    ; corresponds to a wavelength of 100 cm. Basically zero energy.

  T_upward       = dblarr(M + 1, M + 1) ; upward transition rates
  T_downward     = dblarr(M + 1)        ; downward transition rates. T_downward[u] = T_{u-1, u}

  T_downward_td  = dblarr(M + 1, M + 1) ; thermal discrete downward rates

; --- upward transitions ---


; -- case u < M
  for u = 1, M do begin

     Eu      = E_tab[u]
     Emax_u  = E_max[u]
     Emin_u  = E_min[u]
     DeltaEu = Emax_u - Emin_u
     Tu      = T_tab[u]
                                ; case L = 0 : excitations from the ground state 
     W1 = Emin_u
     W4 = Emax_u
     E_array   = makelogtab(W1, W4, NEarr)
     DE_over_E = DX_over_X(W1, W4, NEarr) 
     Cabs      = pi*a^2*Qabs(a, Z, E_array/eV)
     T_upward[u , 0] = c /Eu * total(nu_uisrf(E_array/eV) * Cabs, /double) * DE_over_E
     
     T_downward_td[0 , u] = 8d *pi/(h^3*c^2)/Eu * (E_max[1] - E_min[1])/DeltaEu *total( E_array^4/(exp(E_array/(k*Tu)) - 1d) * Cabs, /double) * DE_over_E
     
                                ; case L > 0 
     for L = 1, u-1 do begin
        EL      = E_tab[L]
        DeltaEL = E_max[L] - E_min[L]
        W1 = Emin_u - E_max[L]
        W2 = min([Emin_u - E_min[L], Emax_u - E_max[L]])
        W3 = max([Emin_u - E_min[L], Emax_u - E_max[L]])
        W4 = Emax_u - E_min[L]

        if(max([W1, Energy_min]) LT W4) then begin
           E_array  = makelogtab(max([W1, Energy_min]), W4, NEarr)
           DE_over_E = DX_over_X(max([W1, Energy_min]), W4, NEarr)
           GuL      = dblarr(NEarr)
           ind1     = where((E_array GT W1) and (E_array LT W2), count1)
           if count1 NE 0 then begin
              GuL[ind1]= (E_array[ind1] - W1)/(DeltaEu * DeltaEL)
           endif
           ind2     = where((E_array GT W2) and (E_array LT W3), count2)
           if count2 NE 0 then begin 
              GuL[ind2]= min([DeltaEu, DeltaEL])/(DeltaEu * DeltaEL)
           endif
           ind3     = where((E_array GT W3) and (E_array LT W4), count3)
           if count3 NE 0 then begin
              GuL[ind3]= (W4 - E_array[ind3])/(DeltaEu * DeltaEL)
           endif
           u_arr    = nu_uisrf(E_array /eV) ; !!! carefull nu_uisrf = nu * u_nu !!! so E uE = nu u_nu
           Cabs     = pi *a^2 *Qabs(a, Z, E_array/eV)
           T_upward[u,L] = c * DeltaEu/(Eu - EL) *  total(GuL * Cabs * u_arr, /double) * DE_over_E
           
           T_downward_td[L,u] = 8d *pi/(h^3*c^2) * DeltaEL/(Eu - EL)  * total( GuL * Cabs * E_array^4/(exp(E_array/(k*Tu)) - 1d), /double) * DE_over_E 
           
        endif
     endfor
  endfor

; -- case u = M, L > 0 : include the upward transitions above the
;    highest level


  EM = E_tab[M]
  for L = 1, M-1 do begin
     EL = E_tab[L] 
     W1 = E_min[M] - E_max[L]
     Wc = E_min[M] - E_min[L]
     integral_1 = 0
     if(max([W1, Energy_min]) LT Wc) then begin
        E_array   = makelogtab(max([W1, Energy_min]), Wc, NEarr)
        DE_over_E = DX_over_X(max([W1, Energy_min]), Wc, NEarr)
        u_arr      = nu_uisrf(E_array /eV)
        Cabs       = pi *a^2 *Qabs(a, Z, E_array/eV)  
        integral_1 = total((E_array - W1)/(Wc - W1) * Cabs * u_arr, /double) * DE_over_E
     endif
     integral_2 = 0
     if (Wc LT 13.6 * eV) then begin ; assuming a radiation field that stops at 13.6 eV 
        E_array   = makelogtab(Wc, 13.6*eV, NEarr)
        DE_over_E = DX_over_X(Wc, 13.6*eV, NEarr)
        u_arr      = nu_uisrf(E_array /eV)
        Cabs       = pi *a^2 *Qabs(a, Z, E_array/eV) 
        integral_2 = total(Cabs * u_arr, /double) * DE_over_E
     endif
     T_upward[M,L] = c/(EM - EL) * (integral_1 + integral_2)
  endfor

; -- case u = M, L = 0
  EM = E_tab[M]
  Wc = E_min[M]
  integral_2 = 0
  if (Wc LT 13.6 * eV) then begin ; assuming a radiation field that stops at 13.6 eV 
     E_array  = makelogtab(Wc, 13.6*eV, NEarr)
     DE_over_E = DX_over_X(Wc, 13.6*eV, NEarr)
     u_arr      = nu_uisrf(E_array /eV)
     Cabs       = pi *a^2 *Qabs(a, Z, E_array/eV) 
     integral_2 = total(Cabs * u_arr, /double) * DE_over_E
  endif
  T_upward[M,0] = c/EM * integral_2

  for u = 2, M do begin         ; "intrabin" contribution to upward transitions
     DeltaEu_1 = E_max[u-1] - E_min[u-1]
     if( DeltaEu_1 GT Energy_min) then begin
        E_array = makelogtab(Energy_min, DeltaEu_1, NEarr)
        DE_overE = DX_over_X(Energy_min, DeltaEu_1, NEarr)
        Cabs = pi * a^2 * Qabs(a, Z, E_array/eV)
        T_upward[u, u-1] = T_upward[u, u-1] + c/(E_tab[u] - E_tab[u-1]) * total( (1d - E_array/DeltaEu_1) * Cabs *nu_uisrf(E_array/eV), /double) * DE_over_E
     endif 
  endfor


; --- downward transitions, "thermal continuous" cooling approximation


  for u = 1, M do begin
     T_downward[u] = 1d/(E_tab[u] - E_tab[u-1]) * total ( (E_tab[u] - E_tab[0:u-1]) * reform(T_downward_td[0:u-1, u]), /double)
     DeltaEu = E_max[u] - E_min[u]
     if(Energy_min LT DeltaEu) then begin ; "intrabin" contribution
        E_array = makelogtab(Energy_min, DeltaEu, NEarr)
        DE_over_E = DX_over_X(Energy_min, DeltaEu, NEarr)
        Cabs = pi * a^2 * Qabs(a, Z, E_array/eV) 
        T_downward[u] = T_downward[u] + 8d * pi /(h^3 * c^2)/(E_tab[u] - E_tab[u-1])  $
                        * total((1d - E_array/DeltaEu) * E_array^4 * Cabs/(exp(E_array/(k*T_tab[u])) - 1d), /double) * DE_over_E
     endif 
  endfor


; we now define Btilde_{j,i} = 1/T_downward[j] * SUM_{u=j..M}(T_upward[u,i]) for i < j
  Tinv   = reverse(T_upward, 1)
  Binv   = total(Tinv, 1, /cumulative)
  Btab   = reverse(Binv, 1) 
  T_downward      = matrix_multiply(T_downward, 1d + dblarr(M + 1), /btranspose)
  T_downward[0,*] = 1d + dblarr(M+1)
  Btilde  = Btab/T_downward

  return, Btilde

end

;---------------------------------------------------------------------------------------------------------------

function distribution, a, Z, Chi, Energy_modes, M, Energy_max 
common path, pathname
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Returns the distribution function for the energy levels.
; needs a previous run of Btilde_ISRF_calc to provide the files
; containing Btilde for different radii. 
; NOTE : This actually returns the distribution for the closest a - value
; of logtab(3.5E-8, 1E-6, 30)

  X_tab = 1d + dblarr(M + 1)
  Energy_max = Energy_max/ 3d

  while X_tab[M] GT 1d-14 do begin

     print, 'E_M = ' +strtrim(Energy_max/eV, 2)+' eV is not high enough in "distribution" (infrared.pro). Using '$
            +strtrim(3d*Energy_max/eV, 2)+' eV instead'
     
     Energy_max = Energy_max * 3d
     Btilde = Chi *Btilde_ISRF_func(a, Z, Energy_modes, M, Energy_max) 
     X_tab = 1d + dblarr(M + 1) ; M may have changed after using Btilde
     for j = 1, M do begin
        X_tab[j] =  total(reform(Btilde[j, 0 : j-1]) * X_tab[0 : j-1], /double)
        X_tab[0:j] = X_tab[0:j]/total(X_tab[0:j], /double)
     endfor
  endwhile

  Energy_bins = Energy_bins(a, Energy_modes, M, Energy_max)  

  return, {Energy_bins: Energy_bins, P_tab: X_tab}

end

;---------------------------------------------------------------------------------------------------------------

function IRemission, a, Z, Chi, nu_tab, Energy_modes, distribution
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

;  distribution = distribution(a, Z, Chi, Energy_modes, M, Energy_max)

  P_tab        = distribution.P_tab
  M            = N_elements(P_tab) - 1
  Energy_bins  = distribution.Energy_bins
  E_tab        = reform(Energy_bins[0,1:M])
  T_tab        = Temp(a, E_tab, Energy_modes)

  N_nu   = N_elements(nu_tab)
  F_nu   = dblarr(N_nu) 

  new_nu_tab = matrix_multiply(nu_tab, 1+dblarr(M), /btranspose)
  new_E_tab  = matrix_multiply(1+dblarr(N_nu), E_tab, /btranspose)
  new_P_tab  = matrix_multiply(1+dblarr(N_nu), P_tab[1:M], /btranspose)
  new_T_tab  = matrix_multiply(1+dblarr(N_nu), T_tab, /btranspose)

  ind = where((new_E_tab LT h*new_nu_tab) or (h*new_nu_tab/(k*new_T_tab) GT 700d), count) ; avoid floating overflow

  if (count NE 0) then begin
     new_P_tab[ind] = 0d
     new_nu_tab[ind] = 700d *k *new_T_tab[ind]/h ; just to avoid floating overflow     
     F_nu = total(new_P_tab/(exp(h*new_nu_tab/(k*new_T_tab)) - 1d), 2)
  endif


  F_nu = 2d*h*nu_tab^3/c^2* pi* a^2 * Qabs(a, Z, h*nu_tab/eV) *  F_nu 
  
  return, F_nu

end

;---------------------------------------------------------------------------------------------------------------

function FGIR_integrals, a, Z, Chi, Energy_modes, M, Energy_max 
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Returns \int F_nu/nu^2 dnu and \int F_nu/nu d nu

  distribution = distribution(a, Z, Chi, Energy_modes, M, Energy_max)
; now M and Energy_max may have a different (higher) value, changed in
; Energy_bins and distribution

  lambda_max = 1d-1      
  nu_min = c/lambda_max
  nu_max = Energy_max/h
  N_lambda = 1000L
  nu_tab      = makelogtab(nu_min, nu_max, N_lambda)
  Dnu_over_nu = DX_over_X(nu_min, nu_max, N_lambda)

  F_nu = IRemission(a, Z, Chi, nu_tab, Energy_modes, distribution)

  FIR_integral = total(F_nu/nu_tab, /double) * Dnu_over_nu
  GIR_integral = total(F_nu, /double) * Dnu_over_nu

  return, {FIR_integral: FIR_integral, GIR_integral: GIR_integral}

end

;---------------------------------------------------------------------------------------------------------------

pro set_up_IR_arrays
  common IR_arrays, a_min, a_max, Na, chi_min, chi_max, Nchi, FIR_integral_charged, $
     GIR_integral_charged, FIR_integral_neutral,GIR_integral_neutral, Tev_tab 

  a_min = 3.5d-8
  a_max = 1d-6
  Na    = 30
  a_tab = makelogtab(a_min, a_max, Na)

  chi_0 = 1d-5
  chi_N = 1d10
  Nchi  = 30
  chi_min = exp(alog(chi_0) - 1d/(2d*Nchi) *alog(chi_N/chi_0)) 
  chi_max = exp(alog(chi_N) - 1d/(2d*Nchi) *alog(chi_N/chi_0)) 
  chi_tab = makelogtab(chi_min, chi_max, Nchi) ; all this business to have an array [1E-5, 3.2E-5, 1E-4, 3.2E-4, etc...]
  
end

;---------------------------------------------------------------------------------------------------------------

pro compute_FGIR_integrals
common cgsconst
common path, SpDust_data_dir
common IR_arrays, a_min, a_max, Na, chi_min, chi_max, Nchi, FIR_integral_charged, $
                  GIR_integral_charged, FIR_integral_neutral,GIR_integral_neutral, Tev_tab 

  set_up_IR_arrays

  FIR_integral_charged = dblarr(Na, Nchi)
  GIR_integral_charged = dblarr(Na, Nchi)
  FIR_integral_neutral = dblarr(Na, Nchi)
  GIR_integral_neutral = dblarr(Na, Nchi)

  for ia = 0, Na - 1 do begin
     a = a_tab[ia]
     Energy_modes = Energy_modes(a)
     M = 100
     Energy_max = 13.6d *eV

     for ichi = 0, Nchi - 1 do begin
        Chi = chi_tab[ichi]
        Chi = chi_tab[ichi]
        FGIR_integral_charged = FGIR_integrals(a, 1, Chi, Energy_modes, M, Energy_max)
        FIR_integral_charged[ia, ichi] = FGIR_integral_charged.FIR_integral
        GIR_integral_charged[ia, ichi] = FGIR_integral_charged.GIR_integral
        FGIR_integral_neutral = FGIR_integrals(a, 0, Chi, Energy_modes, M, Energy_max)
        FIR_integral_neutral[ia, ichi] = FGIR_integral_neutral.FIR_integral
        GIR_integral_neutral[ia, ichi] = FGIR_integral_neutral.GIR_integral
     endfor

  endfor

  openw, 1, SpDust_data_dir+'FIR_integral_charged_'+strtrim(Na, 2)+'a_'+strtrim(Nchi, 2)+'chi'
  printf, 1, FIR_integral_charged
  close, 1
  openw, 1, SpDust_data_dir+'GIR_integral_charged_'+strtrim(Na, 2)+'a_'+strtrim(Nchi, 2)+'chi'
  printf, 1, GIR_integral_charged
  close, 1
  openw, 1, SpDust_data_dir+'FIR_integral_neutral_'+strtrim(Na, 2)+'a_'+strtrim(Nchi, 2)+'chi'
  printf, 1, FIR_integral_neutral
  close, 1
  openw, 1, SpDust_data_dir+'GIR_integral_neutral_'+strtrim(Na, 2)+'a_'+strtrim(Nchi, 2)+'chi'
  printf, 1, GIR_integral_neutral
  close, 1

end

;---------------------------------------------------------------------------------------------------------------

function FGIR_integrals_interpol, a, Z, Chi
common IR_arrays, a_min, a_max, Na, chi_min, chi_max, Nchi, FIR_integral_charged, $
                  GIR_integral_charged, FIR_integral_neutral,GIR_integral_neutral, Tev_tab  

  if Z NE 0 then begin
     FIR_integral_tab = FIR_integral_charged
     GIR_integral_tab = GIR_integral_charged
  endif else begin
     FIR_integral_tab = FIR_integral_neutral
     GIR_integral_tab = GIR_integral_neutral
  endelse

; --- finding the indices and coefficients alpha and beta for a and chi s.t. a = alpha *a_i + (1-alpha)*a_i+1 
  a_tab = makelogtab(a_min, a_max, Na)
  chi_tab = makelogtab(chi_min, chi_max, Nchi)

  if a LE min(a_tab) then begin ; approximate values less then atab[0] by atab[0] 
     ia = 0             
     alpha = 1d
  endif else begin
     if a GE max(a_tab) then begin ; approximate values greater then atab[Na-1] by atab[Na-1] 
        ia = Na - 2
        alpha = 0d
     endif else begin
        ia = max(where(a_tab LE a))
        alpha = 1d - Na * alog(a/a_tab[ia])/alog(a_max/a_min)
     endelse
  endelse

; --- case of low radiation field : FIR, GIR are linear in Chi --- 

  if chi LE min(chi_tab) then begin
     FIR_integral_0 = exp( alpha *alog(FIR_integral_tab[ia, 0]) + (1d - alpha) *alog(FIR_integral_tab[ia+1, 0]) )   
     FIR_integral = chi/chi_tab[0] *FIR_integral_0
     
     GIR_integral_0 = exp( alpha *alog(GIR_integral_tab[ia, 0]) + (1d - alpha) *alog(GIR_integral_tab[ia+1, 0]) )   
     GIR_integral = chi/chi_tab[0] *GIR_integral_0

     return, {FIR_integral: FIR_integral, GIR_integral: GIR_integral}

  endif 

; --- case of high radiation field : approximate FIR, GIR by a power-law ---

  if chi GE max(chi_tab) then begin
     print, 'you are using chi= '+strtrim(Chi, 2)+'! The code is written for chi < 1E10. It assumes a power-law in chi for higher values. If you wish to extend its validity to higher values, set up chi_max to a higher value in "set_up_FGIR_integrals" (infrared.pro), and run "compute_FGIR_integrals" and "@compileSpDust".'
     FIR_integral_last        = exp( alpha *alog(FIR_integral_tab[ia, Nchi-1]) + (1d - alpha) *alog(FIR_integral_tab[ia+1, Nchi -1]) )
     FIR_integral_before_last = exp( alpha *alog(FIR_integral_tab[ia, Nchi-2]) + (1d - alpha) *alog(FIR_integral_tab[ia+1, Nchi -2]) )  
     index = alog(FIR_integral_last/FIR_integral_before_last) / alog(chi_tab[Nchi-1]/chi_tab[Nchi -2])
     FIR_integral = (chi/chi_tab[Nchi-1])^index *FIR_integral_last

     GIR_integral_last        = exp( alpha *alog(GIR_integral_tab[ia, Nchi-1]) + (1d - alpha) *alog(GIR_integral_tab[ia+1, Nchi -1]) )
     GIR_integral_before_last = exp( alpha *alog(GIR_integral_tab[ia, Nchi-2]) + (1d - alpha) *alog(GIR_integral_tab[ia+1, Nchi -2]) )  
     index = alog(GIR_integral_last/GIR_integral_before_last) / alog(chi_tab[Nchi-1]/chi_tab[Nchi -2])
     GIR_integral = (chi/chi_tab[Nchi-1])^index *GIR_integral_last

     return, {FIR_integral: FIR_integral, GIR_integral: GIR_integral}

  endif 

; --- for the other cases, just interpolate ---

  ichi =  max(where(chi_tab LE chi))
  beta  = 1d - Nchi * alog(chi/chi_tab[ichi])/alog(chi_max/chi_min)

  FIR_integral = exp( alpha* (beta* alog(FIR_integral_tab[ia, ichi]) + (1d - beta) *alog(FIR_integral_tab[ia, ichi+1])) $
                      + (1d - alpha) *(beta* alog(FIR_integral_tab[ia+1, ichi]) + (1d - beta) *alog(FIR_integral_tab[ia+1, ichi+1]))) 

  GIR_integral = exp( alpha* (beta* alog(GIR_integral_tab[ia, ichi]) + (1d - beta) *alog(GIR_integral_tab[ia, ichi+1])) $
                      + (1d - alpha) *(beta* alog(GIR_integral_tab[ia+1, ichi]) + (1d - beta) *alog(GIR_integral_tab[ia+1, ichi+1])))


  return, {FIR_integral: FIR_integral, GIR_integral: GIR_integral}

end

;---------------------------------------------------------------------------------------------------------------

function FGIR, env, a, Zg
common cgsconst, pi, c, q, k, mp, me, h, debye, eV
common grainparams, a2, d, rho, epsilon

; --- Revised February 2010 to account for tumbling disklike grains --- 
; Also correcting a mistake in G_IR (each photon carries 2 hbar^2 !)

    Chi = env.Chi
    T   = env.T
    nh  = env.nh
    Inertia = Inertia(a)
    acx     = acx(a)
    tau_H   = 1d/(nh *mp *sqrt(2d* k*T/(pi*mp)) * 4d*pi*acx^4/(3d*Inertia))

    FGIR_integrals = FGIR_integrals_interpol(a, Zg, Chi)

    IntF = FGIR_integrals.FIR_integral
    IntG = FGIR_integrals.GIR_integral

    FIR = 2d * tau_H/(pi * Inertia) * IntF
    GIR = h * tau_H/(3d *pi * Inertia * k * T) * IntG    ; corrected (= 2 * old version)

; --- Addition : disklike grains, K randomized ---
; --- Note that this is always used, as during thermal spikes the
;     orientation is randomized 

    if a LT a2 then begin
       FIR = 5d/3d * FIR
    endif

return, {FIR: FIR, GIR: GIR} 

end

;---------------------------------------------------------------------------------------------------------------

function FGIR_averaged, env, a, fZ

    f0 = fZ[1,0]

    FGIR_neutral = FGIR(env, a, 0)
    FGIR_charged = FGIR(env, a, 1)

    FIR = f0 * FGIR_neutral.FIR + (1d - f0) * FGIR_charged.FIR
    GIR = f0 * FGIR_neutral.GIR + (1d - f0) * FGIR_charged.GIR

return, {FIR: FIR, GIR: GIR}

end

;---------------------------------------------------------------------------------------------------------------


