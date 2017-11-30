;======================================= SPDUST.2 =====================================;
;                                                                                      ;
;                                     COLLISIONS.PRO                                   ;
;                                                                                      ;
; Normalized damping and excitation rates for collisions with neutral                  ;
; species and ions.                                                                    ;
; Refers to Section 5 of Ali-Haimoud, Hirata & Dickinson, 2009, MNRAS, 395, 1055       ;
; and Section 7 of Silsbee, Ali-Haimoud & Hirata, 2010                                 ;           
;                                                                                      ;
; Author : Yacine Ali-Haimoud        yacine@tapir.caltech.edu                          ; 
;                                                                                      ;
; Revision history : Written December 2008                                             ;
;                    Revised March 2010 : different damping rates for disklike         ; 
;                    grains spinning around a non principal axis                       ;
;======================================================================================;



function Tev, a, Chi
common cgsconst

; Thermal spikes limit

    hnu_tab  = makelogtab(1d-2, 13.6d, 500)
    Qabs     = Qabs(a, 1, hnu_tab)
    nu_uisrf =  nu_uisrf(hnu_tab)

    E_photon = total(Qabs * nu_uisrf, /double)/ $
               total(Qabs * nu_uisrf/ hnu_tab, /double) *eV
    Energy_modes = Energy_modes(a)
    T_q = Temp(a, E_photon, Energy_modes)

; Steady emission limit. Take cross section of ionized PAH

    Dnu_over_nu = DX_over_X(1d-2, 13.6d, 500d)
    Qu_star = Chi *total(Qabs * nu_uisrf, /double) * Dnu_over_nu    ;Q_star * u_star in DL98b 

    hnu_0 = 1d-4 
    Qlambda_0 = Qabs(a, 1, hnu_0) * (c/(hnu_0 *eV/h))^2 

    zeta6  = pi^6/945d        ; zeta(6)
    gamma6 = 120d

    T_c = h*c/k*(Qu_star/(8d *pi *h *c *Qlambda_0 *gamma6 *zeta6))^(1d/6d)

    return, max([T_q, T_c])

end

;----------------------------------------------------------------------------------

pro compute_Tev
common path, SpDust_data_dir
common IR_arrays, a_min, a_max, Na, chi_min, chi_max, Nchi, FIR_integral_charged, GIR_integral_charged, FIR_integral_neutral,GIR_integral_neutral, Tev_tab 

    set_up_IR_arrays

    a_tab   = makelogtab(a_min, a_max, Na)
    chi_tab = makelogtab(chi_min, chi_max, Nchi)

    Tev_tab = dblarr(Na, Nchi)

    for ia = 0, Na -1 do begin
      a = a_tab[ia]
      for ichi = 0, Nchi-1 do begin
         Chi = chi_tab[ichi]
         Tev_tab[ia, ichi] = Tev(a, Chi)
     endfor
    endfor

    openw, 1, SpDust_data_dir+'Tev_'+strtrim(Na,2)+'a_' + strtrim(Nchi, 2)+'chi'
    printf, 1, Tev_tab
    close, 1

end

;------------------------------------------------------------------------------------

function Tev_interpol, a, Chi
common IR_arrays, a_min, a_max, Na, chi_min, chi_max, Nchi, FIR_integral_charged, GIR_integral_charged, FIR_integral_neutral,GIR_integral_neutral, Tev_tab 


; --- finding the indices and coefficients alpha and beta for a and chi s.t. a = alpha *a_i + (1-alpha)*a_i+1 
    a_tab = makelogtab(a_min, a_max, Na)
    chi_tab = makelogtab(chi_min, chi_max, Nchi)


   if a LE min(a_tab) then begin       ; approximate values less then atab[0] by atab[0] 
      ia = 0             
      alpha = 1d
    endif else begin
      if a GE max(a_tab) then begin  ; approximate values greater then atab[Na-1] by atab[Na-1] 
        ia = Na - 2
        alpha = 0d
      endif else begin
        ia = max(where(a_tab LE a))
        alpha = 1d - Na * alog(a/a_tab[ia])/alog(a_max/a_min)
      endelse
    endelse


    if chi LE min(chi_tab) then begin
      print, 'you are using chi= '+strtrim(Chi, 2)+'! The code is written for chi > 1E-6. If you wish to extend its validity to lower values, set up chi_min to a lower value in "set_up_Tev" (collisions.pro), and run "compute_Tev" and "@compileSpDust".'
      ichi = 0
      beta = 1d
    endif else begin

      if chi GE max(chi_tab) then begin
        print, 'you are using chi= '+strtrim(Chi, 2)+'! The code is written for chi < 1E8. If you wish to extend its validity to higher values, set up chi_max to a higher value in "set_up_Tev" (collisions.pro), and run "compute_Tev" and "@compileSpDust".'
        ichi = Nchi - 2
        beta = 0d
       endif else begin
          ichi = max(where(chi_tab LE chi))
          beta  = 1d - Nchi * alog(chi/chi_tab[ichi])/alog(chi_max/chi_min)
       endelse
    endelse

; --- now interpolating ---

    return, exp( alpha* (beta* alog(Tev_tab[ia, ichi]) + (1d - beta) *alog(Tev_tab[ia, ichi+1])) $
                 + (1d - alpha) *(beta *alog(Tev_tab[ia+1, ichi]) + (1d - beta) *alog(Tev_tab[ia+1, ichi+1]))) 

end 

;------------------------------------------------------------------------------------

function Tev_effective, env, a, get_info = get_info
common cgsconst,    pi, c, q, k, mp, me, h, debye, eV
common grainparams, a2, d, rho, epsilon

; Effective Evaporation temperature as in AHD09, Section 5.1.4.
; Modified February 2010 to include the get_info tag, which only
; returns the ratio R_{coll/abs}/Nsites 

; --- First check if the user wants a constant evaporation temperature
;     instead ---

    tag_names  = tag_names(env)
    Tev_present = STRCMP( tag_names, 'Tev', /FOLD_CASE )
    if total(Tev_present) NE 0 then begin
       return, env.Tev
    endif 

; --- Now the more physical model for Tev (section 5.1.4 of the paper) ---

    nh = env.nh
    T  = env.T
    Chi = env.Chi

    hnu_tab     = makelogtab(1d-2, 13.6d, 500)
    Dnu_over_nu = DX_over_X(1d-2, 13.6d, 500)
    Qabs     = Qabs(a, 1, hnu_tab)
    nu_uisrf =  nu_uisrf(hnu_tab)

; ratio of collision rate to absorption rate

    ratio  = nh *sqrt(8d *k *T/(pi *mp))/(Chi  * total( Qabs *nu_uisrf/(hnu_tab*eV), /double) *Dnu_over_nu *c)

; number of sites 

    if a LT a2 then begin
      Nsites = N_C(a)
    endif else begin
      Nsites = N_C(a) *3d *d/a
    endelse 
   
    if keyword_set(get_info) then begin
       return, ratio/Nsites
    endif

    if Nsites GT ratio then begin
      Tev = Tev_interpol(a, Chi) 
    endif else begin              ; in that case no more sticking collisions : the atoms bounce back
      Tev = T
    endelse

    return, Tev

end

;------------------------------------------------------------------------------------

function FGn , env, a, Tev, Zg_tab, tumbling = tumbling     
common cgsconst

; ---> Collisions with NEUTRAL impactors.
; Returns F_n and G_n for an array of grain charges (AHD09 Eq.(74)).
; Tev need to be calculated before hand for speed reasons.

; Modified February 2010 : In the case of disklike grains, Fn_in is non-zero.

;--- polarizabilities for neutral species : Hydrogen, Helium and H2 ---

    a0       = 0.52918d-8      ; Bohr radius
    alpha_H  = 4.5d *a0^3      ; J. Chem. Phys. 49, 4845 (1968); DOI:10.1063/1.1669968  
    alpha_He = 1.38d *a0^3     ; M A Thomas et al 1972 J. Phys. B: At. Mol. Phys. 5 L229-L232 
    alpha_H2 = 5.315d *a0^3    ; W.C. Marlow, Proc. Phys. Soc 1965, vol. 86 (ground state)
;--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

    T     = env.T   
    xh    = env.xh
    y     = env.y
    acx   = acx(a)

    e_n = sqrt(q^2/(2d*acx^4*k*T) * Zg_tab^2)
    e_e = sqrt(q^2/(2d*acx^4*k*Tev) * Zg_tab^2)

; Contribution of atomic Hydrogen 

    eps_n = e_n * sqrt(alpha_H)
    eps_e = e_e * sqrt(alpha_H)
    FnH_ev = (1d - xh - y)*(exp(-eps_n^2)+ sqrt(pi)* eps_n * erf(eps_n)) $ 
        /(exp(-eps_e^2)+ sqrt(pi)* eps_e * erf(eps_e)) *(exp(-eps_e^2) + 2d *eps_e^2)  
    GnH_in = 0.5d * (1d - xh - y)*(exp(-eps_n^2) + 2d *eps_n^2)

   if (keyword_set(tumbling)) then begin
      FnH_in = FnH_ev *2d/3d/(1d + sqrt(2d/3d)*eps_n)  ; (r_c/acx)^2 = sqrt(2/3) epsilon_n 
   endif else begin
      FnH_in = 0d
   endelse

; Contribution of neutral Helium. Assuming all Helium is neutral and
; XHe = 0.25 (mass fraction)
; nHe/nH * (mHe/mH)^(1/2) = 1/12 * 2 = 1/6

    eps_n = e_n * sqrt(alpha_He)
    eps_e = e_e * sqrt(alpha_He)
    FnHe_ev   = 1d/6d * (exp(-eps_n^2)+ sqrt(pi)* eps_n * erf(eps_n)) $ 
        /(exp(-eps_e^2)+ sqrt(pi)* eps_e * erf(eps_e)) *(exp(-eps_e^2) + 2d *eps_e^2)  
    GnHe_in = 1d/12d* (exp(-eps_n^2) + 2d *eps_n^2)

   if (keyword_set(tumbling)) then begin
      FnHe_in = FnHe_ev *2d/3d/(1d + sqrt(2d/3d)*eps_n)  ; (r_c/acx)^2 = sqrt(2/3) epsilon_n 
   endif else begin
      FnHe_in = 0d
   endelse

; Contribution of molecular Hydrogen 
; nH2/nH *sqrt(mH2/mH) = y/2 *sqrt(2) = y/sqrt(2) 

    eps_n = e_n * sqrt(alpha_H2)
    eps_e = e_e * sqrt(alpha_H2)
    FnH2_ev   = y/sqrt(2d) *(exp(-eps_n^2)+ sqrt(pi)* eps_n * erf(eps_n)) $ 
        /(exp(-eps_e^2)+ sqrt(pi)* eps_e * erf(eps_e)) *(exp(-eps_e^2) + 2d *eps_e^2)  
    GnH2_in = 0.5d *y/sqrt(2d) *(exp(-eps_n^2) + 2d *eps_n^2)
   
   if (keyword_set(tumbling)) then begin
      FnH2_in = FnH2_ev *2d/3d/(1d + sqrt(2d/3d)*eps_n)  ; (r_c/acx)^2 = sqrt(2/3) epsilon_n 
   endif else begin
      FnH2_in = 0d 
   endelse

; Sum of the three

    Fn_ev = FnH_ev + FnHe_ev + FnH2_ev
    Fn_in = FnH_in + FnHe_in + FnH2_in
    Fn = Fn_ev + Fn_in
    Gn = GnH_in + GnHe_in + GnH2_in + 0.5d * Tev/T * Fn_ev

    return, {Fn: Fn, Gn: Gn}

end

;----------------------------------------------------------------------------------

function FGn_averaged, env, a, Tev, fZ, tumbling = tumbling

; ---> Returns the actual Fn and Gn, averaged over grain charges,
; given the grain charge distribution function

; Modified February 2010 : includes the /tumbling option

    FGn = FGn(env, a, Tev, fZ[0,*], tumbling = tumbling)
    Fn = total(FGn.Fn * fZ[1,*], /double)
    Gn = total(FGn.Gn * fZ[1,*], /double)

    return, {Fn: Fn, Gn: Gn}

end

;----------------------------------------------------------------------------------

function g1g2, psi, mu_tilde

; CHANGED november 09. Put g1 and g2 together (and minor changes in array dimensions)
; ---> Returns a structure {g1(psi, mu_tilde), g2(psi, mutilde)} 
; where psi = double and mu_tilde = array(Ndipole)
; g1 is defined in AHD09 Eq.(90)
; g2 is defined in AHD09 Eq.(91)

    Ndipole = N_elements(mu_tilde)
    g1 = dblarr(Ndipole)
    g2 = dblarr(Ndipole)
 
    index = where(mu_tilde LE abs(psi) , count)
    if (count NE 0) then begin
        mu_i  = mu_tilde[index]
        if (psi LT 0d) then begin
             g1[index] = 1d - psi
             g2[index] = 1d - psi + 0.5d *psi^2 + 1d/6d *mu_i^2
        endif else begin
             if (psi GT 600d) then begin  ; to avoid floating underflow
                g1[index] = 0d
                g2[index] = 0d
             endif else begin 
                g1[index] = exp(-psi) *sinh(mu_i)/mu_i
                g2[index] = g1[index]
             endelse
        endelse
     endif

     index = where(mu_tilde GT abs(psi) , count)
     if count NE 0 then begin
          mu_i = mu_tilde[index]
          g1[index] = (1d - exp(-(psi + mu_i)) + mu_i - psi + 0.5d *(mu_i - psi)^2)/(2d *mu_i)
          g2[index] = g1[index] + (mu_i - psi)^3/(12d * mu_i)
     endif

     return, {g1: g1, g2: g2}

end

;----------------------------------------------------------------------------------

function h1, phi, mu_tilde
common cgsconst

; ---> Returns h1(phi, mu_tilde) as in AHD09 Eq.(103)

    mu = mu_tilde

    u_0       = 0.5d * (- phi + sqrt(phi^2 + 4d*mu)) 
    Dl98      = 1d + sqrt(pi)/2d * phi
    coeff_1   = mu/4 + phi^2/(4d*mu) + (1d - mu)/(2d*mu)
    coeff_erf = sqrt(pi) *phi *(3d - 2d*mu)/(8d *mu)
    coeff_exp = -(4d + phi^2 + phi *sqrt(phi^2 + 4d*mu))/(8d*mu)

    return, DL98 + coeff_1 + coeff_erf *erf(u_0) + coeff_exp *exp(-u_0^2) 

end

;----------------------------------------------------------------------------------

function h2, phi, mu_tilde
common cgsconst

; ---> Returns h2(phi, mu_tilde) as in AHD09 Eq.(104)

    mu = mu_tilde

    u_0       = 0.5d * (- phi + sqrt(phi^2 + 4d*mu)) 
    Dl98      = 1d + 3d*sqrt(pi)/4d * phi + 0.5d *phi^2
    coeff_1   = mu^2/12d + mu/4d + phi^2/(2d*mu) + (1d - mu)/(2d*mu) - phi^2/4d
    coeff_erf = sqrt(pi) *phi/(32d*mu) *(4d*mu^2 - 12d*mu + 15d + 2d*phi^2)
    coeff_exp = (phi^2*(2d*mu -9d) - 16d + (2d*mu -7d)*phi*sqrt(phi^2 + 4d*mu))/(32d*mu)

    return, DL98 + coeff_1 + coeff_erf *erf(u_0) + coeff_exp *exp(-u_0^2) 

end

;----------------------------------------------------------------------------------

function FGi, env, a, Tev, Zg, mu_tab
common cgsconst

; ---> Returns Fi and Gi of the same dimensions as mu_tab, for one
; value of the grain charge Zg. AHD09 Eq.(89), (100), (101), (102).

; Modified February 2010 : minor changes due to the minor changes in
; the format of g1, g2.

    T        = env.T   
    xh       = env.xh
    xM       = env.xC
    acx      = acx(a)
    mu_tilde = q * mu_tab/(acx^2 *k *T)

 
;--- polarizabilities for neutral species corresponding to the
;    incoming ions : Hydrogen and Carbon ---

    a0      = 0.53d-8       ; Bohr radius
    alpha_H = 4.5d *a0^3     ; J. Chem. Phys. 49, 4845 (1968); DOI:10.1063/1.1669968  
    alpha_C = 1.54d-24      ; Physical Review A, vol. 5, Issue 2, pp. 516-520
;-----------------------------------------

; ------ Neutral grains ------

    if Zg EQ 0 then begin
        phi   = sqrt(2d)*q/sqrt(acx*k*T)
        Fi    = (xh + xM *sqrt(12d)) *h1(phi, mu_tilde)
        Gi_in = 0.5d *(xh + xM *sqrt(12d))*h2(phi, mu_tilde)  
    endif else begin

; ------ Charged grains ------
        psi = Zg *q^2/(acx*k*T)
        e_i = sqrt(q^2/(2d * acx^4*k*Tev) * Zg^2)

; Contribution of H+, C+
        eps_i = e_i * sqrt([alpha_H, alpha_C])
        Fi    = [xh, xM] * (exp(-eps_i^2) + 2d*eps_i^2)/(exp(-eps_i^2) + sqrt(pi) *eps_i *erf(eps_i))
        g1g2  = g1g2(psi, mu_tilde)
        g1    = g1g2.g1
        g2    = g1g2.g2
        Fi    = total(Fi, /double) *g1
        Gi_in = 0.5d *(xh + xM *sqrt(12d)) *g2
    endelse
   
    Gi_ev = 0.5d *Tev/T *Fi

    return, {Fi: Fi, Gi: Gi_in + Gi_ev}

end

;----------------------------------------------------------------------------------

function FGi_averaged, env, a, Tev, mu_tab, fZ
; ---> Returns the actual Fi and Gi, averaged over grain charges,
; given the grain charge distribution function. This is an array of
; the same dimensions as mu_tab.

    Fi = 0d * mu_tab     
    Gi = Fi   
 
    NZg = N_elements(fZ[0,*])    
    for i = 0, NZg - 1 do begin
        FGi = FGi(env, a, Tev, fZ[0,i], mu_tab)
        Fi  = Fi + fZ[1,i] *FGi.Fi
        Gi  = Gi + fZ[1,i] *FGi.Gi   
    endfor   

return, {Fi: Fi, Gi: Gi}

end

;----------------------------------------------------------------------------------
