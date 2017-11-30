;================================= SPDUST.2 =======================================;
;                                                                                  ;
;                                EMISSIVITY.PRO                                    ;  
;                                                                                  ;
; Computation of the rotational distribution function and the emissivity.          ;
; Refers to Sections 2 and 4 of Ali-Haimoud, Hirata & Dickinson, MNRAS, 395. 1055, ;
; and Sections 3 and 4 of Silsbee, Ali-Haimoud & Hirata, 2010                      ;               
;                                                                                  ;
; Author : Yacine Ali-Haimoud        yacine@tapir.caltech.edu                      ;
;                                                                                  ;
; Revision history : Written December 2008                                         ;
;                    Revised March 2010 : modifications for grains                 ;
;                    spinning around a non-principal axis of inertia               ;
;                    May 2012: fixed a bug line 418                                ;
;==================================================================================;



function tau_H, env, a
common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Returns the characteristic damping time through collisions with
; neutral H atoms (AHD09 Eq. (24)).

    T       = env.T
    nh      = env.nh
    acx     = acx(a)
    Inertia = Inertia(a)

    return, 1d/(nh *mp *sqrt(2d* k*T/(pi*mp)) * 4d *pi *acx^4/(3d *Inertia))

end

;------------------------------------------------------------------------------------------------

function tau_ed_inv, env, a, mu_ip, mu_op, tumbling = tumbling 
  common cgsconst, pi, c, q, k, mp, me, h, debye, eV

; ---> Returns the INVERSE OF the characteristic damping time through electric dipole
; radiation, given mu_ip, mu_op.
; These two should be of the same dimensions, and so will be the output.
; ---- Modified February 2010 : tumbling disklike grains ----

   T       = env.T
   Inertia = Inertia(a)
   
   if keyword_set(tumbling) then begin  ; Tumbling disklike grains, SAH10 Eq. (47) 
        return, 3d *k *T *(82d/45d *mu_ip^2 + 32d/9d * mu_op^2)/(Inertia^2 *c^3)
   endif else begin  ; Spherical grains or diklike grains with K = J, AHD09 Eq. (29)
        return, 2d *k *T *mu_ip^2/(Inertia^2 *c^3)                             
   endelse

end 

;-------------------------------------------------------------------------------------------------

function f_rot, env, a, fZ, mu_ip, mu_op, tumbling = tumbling
common cgsconst,    pi, c, q, k, mp, me, h, debye, eV
common grainparams, a2, d, rho, epsilon

; ---> Returns the rotational distribution function f_a normalized s.t. 
; /int f_a(/omega)  4 pi omega^2 d/omega = 1. 
; The distribution function is calculated at frequencies centered on the approximate 
; peak frequency of the emitted power.
; The output is a structure {omega, f_a} containing the values of
; the frequencies at which the function is calculated, and the
; rotational distribution function f_a, array[Nmu, Nomega], where
; Nmu = N_elements(mu_ip) = N_elements(mu_op)

; Modified March 2010 : tumbling disklike grains

 ; --- Number of omega values for which the function is calculated ---
    Nomega = 1000L  

    Nmu = N_elements(mu_ip)

    T       = env.T 
    Inertia = Inertia(a)
   
; --- Characteristic timescales ---
    tau_H      = tau_H(env, a) 
    tau_ed_inv = tau_ed_inv(env, a, mu_ip, mu_op, tumbling = tumbling)      

; --- Evaporation temperature ---
    Tev = Tev_effective(env, a)

; --- F's and G's (except for plasma drag) ---
    FGn  = FGn_averaged(env, a, Tev, fZ, tumbling = tumbling)
    Fn   = FGn.Fn
    Gn   = FGn.Gn
 
    mu_tot = sqrt(mu_ip^2 + mu_op^2)
    FGi  = FGi_averaged(env, a, Tev, mu_tot, fZ)  
    Fi   = FGi.Fi
    Gi   = FGi.Gi

    FGIR = FGIR_averaged(env, a, fZ)   
    FIR  = FGIR.FIR
    GIR  = FGIR.GIR
 
    FGpe = FGpe_averaged(env, a, fZ) 
    Fpe  = FGpe.Fpe
    Gpe  = FGpe.Gpe  

    GH2  = GH2(env, a) 

; --- Array of omegas around the approximate peak ---

    omega_peak_th = sqrt(6d *k *T/Inertia)    ; peak of the spectrum if thermal rotation

; -- Peak frequency for the lowest and highest values of mu_ip, mu_op

    FGp = FGp_averaged(env, a, fZ, omega_peak_th, [min(mu_ip), max(mu_ip)], $ 
                             [min(mu_op), max(mu_op)], tumbling = tumbling)
    Fp_th  = reform(FGp.Fp)
    Gp_th  = reform(FGp.Gp)

    F_low    = Fn + min(Fi) + FIR + Fpe +       min(Fp_th)
    G_low    = Gn + min(Gi) + GIR + Gpe + GH2 + min(Gp_th)
    xi_low   = 8d * G_low/F_low^2 * tau_H *min(tau_ed_inv)      ; non-Gaussianity parameter AHD09 Eq.(167) 
   
    F_high  = Fn + max(Fi) + FIR + Fpe +       max(Fp_th)
    G_high  = Gn + max(Gi) + GIR + Gpe + GH2 + max(Gp_th)
    xi_high = 8d * G_high/F_high^2 * tau_H *max(tau_ed_inv)  

    omega_peak_low   = omega_peak_th *sqrt(2d * G_low/F_low/(1d + sqrt(1d + xi_low))) ; AHD09 Eq. (166)
    omega_peak_high  = omega_peak_th *sqrt(2d * G_high/F_high/(1d + sqrt(1d + xi_high)))

; -- Array omega --

    omega_min = 5d-3 * min([omega_peak_low, omega_peak_high])
    omega_max = 6d * max([omega_peak_low, omega_peak_high])
    omega     = makelogtab(omega_min, omega_max, Nomega) 
    Dln_omega = DX_over_X(omega_min, omega_max, Nomega) 

; --- Fp(omega), Gp(omega) ---

    FGp = FGp_averaged(env, a, fZ, omega, mu_ip, mu_op, tumbling = tumbling)
    Fp  = FGp.Fp
    Gp  = FGp.Gp

; --- Rotational distribution function, AHD09 Eq.(33) ---

    f_a       = dblarr(Nomega, Nmu)
   
    F             = Fn + FIR + Fpe +       matrix_multiply(1d + dblarr(Nomega), Fi, /btranspose) + Fp    
    G             = Gn + GIR + Gpe + GH2 + matrix_multiply(1d + dblarr(Nomega), Gi, /btranspose) + Gp 
    tau_ed_inv    = matrix_multiply(1d + dblarr(Nomega), tau_ed_inv, /btranspose) 
    omega_tab     = matrix_multiply(omega, 1d + dblarr(Nmu), /btranspose)

    X         = Inertia * omega_tab^2/(k *T)
    integrand = F/G *X + tau_H/(3d *G) *tau_ed_inv *X^2
    exponent  = total(integrand, 1, /cumulative, /double) *Dln_omega  
    norm      = 4d *pi *total(omega_tab^3 * exp(- exponent), 1, /double) *Dln_omega     ; = array[Nmu]
    norm      = matrix_multiply(1d + dblarr(Nomega), norm, /btranspose)

    f_a     = 1d/norm *exp(-exponent)    ; rotational distribution function s.t. /int f(omega)* 4pi omega^2 domega =1
    f_a     = transpose(f_a)             ; array[Nmu, Nomega] for faster integration over dipole values

 
    return, {omega: omega, f_a: f_a}

end

;-------------------------------------------------------------------------------------------------

function mu2_fa, env, a, fZ, mu_rms, ip, Ndipole, tumbling = tumbling
  common grainparams, a2, d, rho, epsilon
; Returns a [3, Nomega] array :
;      [omega, <mu_ip^2 fa(omega)>/<mu_ip^2>, 
;              <mu_op^2 fa(omega)>/<mu_op^2>]
;(if <mu_ip^2> or <mu_op^2> = 0, the corresponding component will be
;zero) 
; <mu_ip^2> = ip * mu_rms^2 
; <mu_op^2> = (1 - ip) * mu_rms^2 
; The average is made over a multi-gaussian distribution of dipole
; moments, with rms sqrt(ip)*mu_rms for the (2D) inplane component, and
; sqrt(1-ip)*mu_rms for the (1D) out-of plane component.
; Ndipole is the number of values used in the averaging, for each 
; component ip, op (so Ndipole^2 values used total).

; Note : I assume the rms of the dipole moment is given by AHD09
; Eq.(11) (instead of assuming mu = mu(intrinsic) + mu(charge), with mu(charge)
; single valued and mu(intrinsic) randomly distributed). That should make only
; a minor difference because mu(charge) << mu(intrinsic).

; BUG corrected from SpDust.1 : Dx_tab was wrong

  op = 1d - ip

  if (Ndipole EQ 1) then begin  ; in case the user does not want the averaging
     f_rot     = f_rot(env, a, fZ, mu_rms *sqrt(ip), mu_rms *sqrt(op), tumbling = tumbling)
     omega     = f_rot.omega
     f_a       = f_rot.f_a
     Nomega    = N_elements(omega)
     mu_ip2_fa = f_a
     mu_op2_fa = f_a
  endif else begin

     xmin  = 5d-3
     xmed  = 0.5d
     xmax  = 5d 
     x_tab  = [makelogtab(xmin, xmed, Ndipole/2), maketab(xmed, xmax, Ndipole/2)]
     Dx_tab = [DX_over_X(xmin, xmed, Ndipole/2) * makelogtab(xmin, xmed, Ndipole/2), $
               (xmax - xmed)/(Ndipole/2) + dblarr(Ndipole/2)]
     
     if (a LT a2) then begin  ; need a 2D gaussian for disc-like grains

        mu_ip    = sqrt(ip) *mu_rms *matrix_multiply(x_tab, 1d + dblarr(Ndipole), /btranspose) 
        mu_op    = sqrt(op) *mu_rms *matrix_multiply(1d + dblarr(Ndipole), x_tab ,/btranspose) 
        Dmu_ip   = matrix_multiply(Dx_tab, 1d + dblarr(Ndipole), /btranspose) 
        Dmu_op   = transpose(Dmu_ip)
        
        if (ip eq 0d) then begin
           Proba =  exp(- 0.5d *mu_op^2/mu_rms^2) *Dmu_op
        endif 
        if (op eq 0d) then begin
           Proba =  mu_ip/mu_rms * exp(- mu_ip^2/mu_rms^2) * Dmu_ip  
        endif
        if ((ip Ne 0d) and (op Ne 0d)) then begin   
           Proba = mu_ip/mu_rms * exp(- mu_ip^2/(ip * mu_rms^2)) *exp(- 0.5d *mu_op^2/(op *mu_rms^2)) *Dmu_ip *Dmu_op 
        endif
        
        Proba    = Proba/total(Proba, /double)
        ; Probability of cell [i,j] = Proba[i,j]

        ; Now make 1D arrays
        indices = indgen(Ndipole * Ndipole)
        mu_ip   = mu_ip[indices]    
        mu_op   = mu_op[indices]
        Proba   = Proba[indices]

        f_rot    = f_rot(env, a, fZ, mu_ip, mu_op, tumbling = tumbling)
        omega    = f_rot.omega
        f_a      = f_rot.f_a
        Nomega   = N_elements(omega)

        Proba    = matrix_multiply(Proba, 1d + dblarr(Nomega), /btranspose)
        mu_ip    = matrix_multiply(mu_ip, 1d + dblarr(Nomega), /btranspose)    
        mu_op    = matrix_multiply(mu_op, 1d + dblarr(Nomega), /btranspose) 
        
        mu_ip2_fa = total(mu_ip^2/mu_rms^2 *Proba * f_a, 1, /double)    
        mu_op2_fa = total(mu_op^2/mu_rms^2 *Proba * f_a, 1, /double)   
        if ((ip Ne 0d) and (op Ne 0d)) then begin    
           mu_ip2_fa = mu_ip2_fa /ip
           mu_op2_fa = mu_op2_fa /op
        endif
  
     endif else begin   ; for spherical grains, averaging on grain orientation first
                        ; -> mu_ip^2 should be replaced by 2/3 mu_tot^2 
                        ; BEFORE averaging over dipoles (so the above calculation does not apply) 
     
        ; mu_tot = mu_rms *x_tab
        ; Dmu    = Dx_tab
        Proba  = x_tab^2 * exp(-1.5d *x_tab^2) *Dx_tab
        Proba = Proba/total(Proba,/double)

        f_rot    = f_rot(env, a, fZ, sqrt(2d/3d) *mu_rms *x_tab, mu_rms/sqrt(3d) *x_tab, tumbling = tumbling)
        omega    = f_rot.omega
        f_a      = f_rot.f_a
        Nomega   = N_elements(omega)
        Proba    = matrix_multiply(Proba, 1d + dblarr(Nomega), /btranspose)
        x_tab    = matrix_multiply(x_tab, 1d + dblarr(Nomega), /btranspose)
        mu2_fa   = total(x_tab^2 *Proba *f_a, 1, /double)
        mu_ip2_fa = mu2_fa      
        mu_op2_fa = mu2_fa  

     endelse

  endelse    
  
  Result = dblarr(3, Nomega)
  Result[0,*] = omega
  Result[1,*] = mu_ip2_fa
  Result[2,*] = mu_op2_fa

  return, Result

end

;-------------------------------------------------------------------------------------------------

function fa_ip_fa_op, omega, mu_ip2_fa, mu_op2_fa 

; Annex functions fa_ip, fa_op defined as follows :
; omega^2 *fa_ip(omega) = 1/4 \int_{omega/3}^omega (3 - omega/x)^2 x^2 <mu_ip^2 fa(x)>/<mu_ip^2> dlogx
;                       + 1/2 \int_{omega}^{infty} (1 -(omega/x)^2)x^2 <mu_ip^2 fa(x)>/<mu_ip^2> dlogx
; fa_op(omega) = 1/8 <mu_op^2 fa(omega/2)>/<mu_op^2>
; mu_ip2_fa = <mu_ip^2 fa>/<mu_ip^2> and idem for
; mu_op2_fa. They are given as arrays of the same dimensions of omega
; (logarithmically spaced).
; Only used for tumbling grains    
; Written March 2010

    Nomega     = N_elements(omega) 
    Dln_omega  = sqrt(omega[1]/omega[0]) - sqrt(omega[0]/omega[1])

; -- values where interpolation at omega/3 is possible
    ind = where(omega/3d GT min(omega), count)

; -- interpolating fa_op at omega/2
    fa_op = interpol(mu_op2_fa, alog(omega), alog(omega[ind]/2d), /spline)/8d

; -- some useful integrals --
    int0 = reverse(total(reverse(mu_ip2_fa), /cumulative, /double)) * Dln_omega  
       ; = \int_{omega}^{infty} <mu_ip^2 fa(x)>/<mu_ip^2> dlogx
    int1 = reverse(total(reverse(omega *mu_ip2_fa), /cumulative, /double)) /omega * Dln_omega
       ; = 1/omega \int_{omega}^{infty} x <mu_ip^2 fa(x)>/<mu_ip^2> dlogx
    int2 = reverse(total(reverse(omega^2 *mu_ip2_fa), /cumulative, /double)) /omega^2 * Dln_omega
       ; = 1/omega^2 \int_{omega}^{infty} x^2 <mu_ip^2 fa(x)>/<mu_ip^2> dlogx

; -- interpolating them at omega/3
    int0_low = interpol(int0, alog(omega), alog(omega[ind]/3d), /spline)
           ; = \int_{omega/3}^{infty} <mu_ip^2 fa(x)>/<mu_ip^2> dlogx
    int1_low = interpol(int1, alog(omega), alog(omega[ind]/3d), /spline)
           ; = 3/omega \int_{omega/3}^{infty} x <mu_ip^2 fa(x)>/<mu_ip^2> dlogx
    int2_low = interpol(int2, alog(omega), alog(omega[ind]/3d), /spline)
           ; = 9/omega^2 \int_{omega/3}^{infty} x^2 <mu_ip^2 fa(x)>/<mu_ip^2> dlogx
   
    fa_ip = (int0_low - 2d *int1_low + int2_low - 3d * int0[ind] + 6d *int1[ind] - 7d *int2[ind])/4d
 
    result = dblarr(3, count)
    result[0,*] = omega[ind]
    result[1,*] = fa_ip
    result[2,*] = fa_op

    return, result

end

;-------------------------------------------------------------------------------------------------

function dP_dnu_dOmega, env, a, beta, ip, Ndipole, tumbling = tumbling  
common cgsconst,    pi, c, q, k, mp, me, h, debye, eV
common grainparams, a2, d, rho, epsilon

; ---> Returns the power radiated by a grain of radius a, per
; frequency interval, per steradian. This function returns an array
; Power_per_grain[2, Nnu] such that
; Power_per_grain[0,*] = nu
; Power_per_grain[1,*] = dP/dnu/dOmega


; --- Modified February 2010 to include the tumbling disklike grain option ---
 
    op = 1d - ip

    fZ         = charge_dist(env, a) 
    Z2         = total(fZ[0,*]^2 *fZ[1,*], /double)
    mu_rms     = rms_dipole(a, Z2, beta)
   
    mu2_fa     = mu2_fa(env, a, fZ, mu_rms, ip, Ndipole, tumbling = tumbling) 
    omega      = reform(mu2_fa[0,*])
    mu_ip2_fa  = reform(mu2_fa[1,*])
    mu_op2_fa  = reform(mu2_fa[2,*])   
    
    if (keyword_set(tumbling)) then begin
         fa_ip_fa_op = fa_ip_fa_op(omega, mu_ip2_fa, mu_op2_fa)
         nu          = reform(fa_ip_fa_op[0,*])/(2d *pi)
         fa_ip       = reform(fa_ip_fa_op[1,*])
         fa_op       = reform(fa_ip_fa_op[2,*])
         dPdnudOmega = 2d/(3d *c^3) *(2d *pi *nu)^6 * 2d *pi * mu_rms^2 *(ip *fa_ip + 2d/3d *op *fa_op)
    endif else begin   ; spherical / disklike with K = J case
         nu          = omega/(2d * pi)
         dPdnudOmega = 2d/(3d *c^3) *omega^6 *2d *pi* mu_rms^2 *ip *mu_ip2_fa
    endelse


; --- keep only the non zero values (for future log-interpolation) ---

     ind = where(dPdnudOmega GT 0d, count)
     if (count NE 0) then begin
        Non_zero_power =  dblarr(2, count) 
        Non_zero_power[0,*] = nu[ind]
        Non_zero_power[1,*] = dPdnudOmega[ind]
     endif else begin
        print, 'Power vanishes for a ='+strtrim(a, 2)
        return, [0d,0d]
     endelse   

     return, Non_zero_power

end

;-------------------------------------------------------------------------------------------------

function emissivity, env, beta, ip, Ndipole, nu_tab, tumbling = tumbling
common grainparams, a2, d, rho, epsilon

; ---> Returns j_nu/nH in cgs units (ergs/s/sr/H atom) for the given
; environment.

; --- Modified February 2010 to account for disklike tumbling dust ---


; --- Grain parameters ---
   a_min     = 3.5d-8
   a_max     = 3.5d-7 ;  used to be 1d-6 but large grains do not contribute near the peak
   Na        = 30
   a_tab     = makelogtab(a_min, a_max, Na)
   Da_over_a = DX_over_X(a_min, a_max, Na)

; --- get the emissivity --- 
   emiss = dblarr(N_elements(nu_tab))

   for ia = 0, Na-1 do begin
    a = a_tab[ia]  

  ; --- Use the tumbling case only for disklike grains ---  
     if (a LT a2) then begin
        power_per_grain = dP_dnu_dOmega(env, a, beta, ip, Ndipole, tumbling = tumbling)
     endif else begin
        power_per_grain = dP_dnu_dOmega(env, a, beta, 2d/3d, Ndipole)  
     endelse

     nu_tab_a        = power_per_grain[0,*]
     emiss_a         = power_per_grain[1,*]
     ind = where((nu_tab GT min(nu_tab_a)) AND (nu_tab LT max(nu_tab_a)), count) 
     if count NE 0 then begin 
       emiss[ind] = emiss[ind] + exp(interpol(alog(emiss_a), alog(nu_tab_a), alog(nu_tab[ind]), /spline)) * size_dist(a) * a * Da_over_a
     endif
   endfor

  return, emiss

end

;------------------------------------------------------------------------------------------------

pro read_gaunt_factor
common path, SpDust_data_dir
common gff_data, gamma2_tab, u_tab, gff_tab

; Using tabulated gaunt factors form Sutherland, 1998, MNRAS, 300, 321

  readcol, SpDust_data_dir+'gff.dat', gamma2, u,  gff, comment = ';', /silent

  Ngamma2 = 41
  Nu      = 81

  gamma2_tab = dblarr(Ngamma2)
  u_tab      = dblarr(Nu)
  gff_tab    = dblarr(Ngamma2, Nu)
  
  gamma2_tab = gamma2[Nu *indgen(Ngamma2)]
  u_tab      = u[0 : Nu-1]

  for i = 0, Ngamma2 -1 do begin
     gff_tab[i,*] = gff[i *Nu + indgen(Nu)]
  endfor

  print, 'gaunt factor stored'

end

;------------------------------------------------------------------------------------------------

function gaunt_factor, gamma2, u
common gff_data, gamma2_tab, u_tab, gff_tab

  Ngamma2 = N_elements(gamma2_tab)

  if(gamma2 GE max(gamma2_tab)) then begin
     index = Ngamma2-1
  endif else begin
     if (gamma2 LE min(gamma2_tab)) then begin
        index = 0
     endif else begin
        index = max(where(gamma2_tab LT gamma2))
        if alog(gamma2_tab[index + 1]/gamma2) LT alog(gamma2/gamma2_tab[index]) then begin
           index = index + 1
        endif
     endelse
  endelse

  gff_new  = reform(gff_tab[index,*])  

  return, interpol(gff_new, u_tab, u)

end

;------------------------------------------------------------------------------------------------

function free_free, env, nu_tab
common cgsconst

; ---> Returns j_nu/nH in cgs units (ergs/s/sr/Hatom) for free - free emission, for the
;                                              given environment.
; References : Radiative Processes in Astrophysics, Rybicki & Lightman. 

  nh = env.nh
  T  = env.T
  xh = env.xh
  xC = env.xC

  factor = 2d^5 *pi *q^6/(3d *me *c^3) *sqrt(2d *pi/(3d *k *me)) /(4d*pi)

  Ry = 13.6 *eV
  gamma2 = Ry/(k *T)            ; assuming Zion = 1
  u      = h*nu_tab/(k*T) 

  return, factor * (xh + xC)^2 *nh /sqrt(T) *exp(-h *nu_tab/(k*T)) *gaunt_factor(gamma2, u)

end


;------------------------------------------------------------------------------------------------

