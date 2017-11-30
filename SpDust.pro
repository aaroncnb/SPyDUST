;=========================== SPDUST.2 ===================================;
;                                                                        ;
;                            SPDUST.PRO                                  ;  
;                                                                        ;
; Main program. Computes the emissivity for environmental parameters     ;
; given as an input.                                                     ;
;                                                                        ;        
; Author : Yacine Ali-Haimoud        yacine@tapir.caltech.edu            ;
;                                                                        ;
; Revision history : Written December 2008                               ;
;                    Revised March 2010 : - disklike grains spinning     ;
;                    around a non-principal axis of inertia              ;
;                     - made somewhat more modular                       ;
;========================================================================;




function invalid_tags, environment
common cgsconst

; Checks that all tags are present in the provided environment. 
; Returns 0 if everything looks ok, 1 otherwise.
; Added February 2010

  tag_names  = tag_names(environment)

; --- checking that all the parameters are present and reading them ---

  invalid = 0

  nh_present = STRCMP( tag_names, 'nh', /FOLD_CASE )
  if total(nh_present) EQ 0 then begin
     print, 'invalid environment : please provide the number density nh'
     invalid = invalid + 1
  endif else begin
     nh = environment.nh
  endelse

  T_present = STRCMP( tag_names, 'T', /FOLD_CASE )
  if total(T_present) EQ 0 then begin
     print, 'invalid environment : please provide the temperature T'
     invalid = invalid + 1
  endif else begin
     T = environment.T
  endelse

  Chi_present = STRCMP( tag_names, 'Chi', /FOLD_CASE )
  if total(Chi_present) EQ 0 then begin
     print, 'invalid environment : please provide the radiation field intensity Chi'
     invalid = invalid + 1
  endif else begin
     Chi = environment.Chi
  endelse

  xh_present = STRCMP( tag_names, 'xh', /FOLD_CASE )
  if total(xh_present) EQ 0 then begin
     print, 'invalid environment : please provide the hydrogen ionisation fraction xh'
     invalid = invalid + 1
  endif else begin
     xh = environment.xh
  endelse

  xC_present = STRCMP( tag_names, 'xC', /FOLD_CASE )
  if total(xC_present) EQ 0 then begin
     print, 'invalid environment : please provide the carbon ionisation fraction xC = n(C+)/nH'
     invalid = invalid + 1
  endif else begin
     xC = environment.xC
  endelse

  y_present = STRCMP( tag_names, 'y', /FOLD_CASE )
  if total(y_present) EQ 0 then begin
     print, 'invalid environment : please provide the molecular hydrogen abundance y = 2*n(H2)/nH'
     invalid = invalid + 1
  endif else begin
     y = environment.y
  endelse

  gamma_present = STRCMP( tag_names, 'gamma', /FOLD_CASE )
  if total(gamma_present) EQ 0 then begin
    print, 'H2 formation efficiency gamma was not provided. gamma = 0 is assumed'
    gamma = 0d
  endif else begin
    gamma = environment.gamma
  endelse

   env = {nh: nh, T: T, Chi: Chi, xh: xh, xC: xC, y: y, gamma: gamma}

   Tev_present = STRCMP( tag_names, 'Tev', /FOLD_CASE )
   if total(Tev_present) NE 0 then begin
     print, 'Assuming a constant evaporation temperature Tev = '+strtrim(environment.Tev,2)
     env = {nh: nh, T: T, Chi: Chi, xh: xh, xC: xC, y: y, gamma: gamma, Tev: environment.Tev}
   endif

  return, {invalid:invalid, env:env}

end

;------------------------------------------------------------------------------------------------

pro SPDUST, environment, output_file, case1 = case1, $
            min_freq = min_freq , max_freq = max_freq, N_freq = n_freq, Ndipole = Ndipole, freefree = freefree, verbose = verbose
  common cgsconst,         pi, c, q, k, mp, me, h, debye, eV
  common size_dist_arrays, bc1e5_tab, alpha_tab, beta_tab, at_tab, ac_tab, C_tab
  common size_dist_params, bc, alpha_g, beta_g, at_g, ac_g, C_g
  common warnings, warning_phi_min, warning_phi_max, warning_psi_min, warning_psi_max, warning_Omega_min, warning_Omega_max

; ---> Computes the emissivity for environment parameters given in a
; structure of the form (example) :
; environment = {nh : 30d, T: 100d, Chi: 1d, xh: 1d-3, xC: 3d-4, y :
; 1d-1, gamma: 0d, dipole: 9.3d, line:7} (see comments below).
; Returns an array of
; frequencies vs emissivities (in Jy cm^2 sr^-1 per H atom) in
; output_file. 
; The user may change the range and number of frequencies
; used (min_freq and max_freq are given in GHz). The default is 200
; frequencies from 0.001 to 500 GHz.
; Ndipole is the number of dipole moments (default = 20) used in the
; numerical averaging over a gaussian distribution.
  
; If the keyword freefree is set, the free-free emissivity is also
; returned. (gaunt factor from Sutherland, 1998, MNRAS, 300, 321.

; If the keyword verbose is set, prints out some info.

; Addition (February 2010) : by default, the dislike grains are supposed to
; have a random orientation with respect to the angular momentum axis
; (as described by Silsbee, Ali-Haimoud & Hirata 2010).
; If the keyword /case1 is set, assumes that the disklike grains are
; rotating aaround their axis of main inertia. 
; the inplane keyword can force the <mu_op^2>/<mu^2> ratio (by default 2/3)
; to be = ip.
; Minor addition : warnings for plasma drag integrals (if the default
; interpolation bounds are not sufficient)
  
  tag_names  = tag_names(environment)
  invalid_tags = invalid_tags(environment)
  invalid = invalid_tags.invalid
  env = invalid_tags.env
 
  dipole_present = STRCMP( tag_names, 'dipole', /FOLD_CASE )
  beta_present   = STRCMP( tag_names, 'beta', /FOLD_CASE )
  if (total(dipole_present) + total(beta_present)  EQ 0) then begin
     print, 'invalid environment : please specify the dipole moment of the grains (either dipole(a = 1E-7 cm) or beta, in debye)'
     invalid = invalid + 1
  endif else begin
     if total(beta_present)  NE 0 then begin
        beta0 = environment.beta *debye
        mu_1d_7 = sqrt(N_C(1d-7) + N_H(1d-7))* environment.beta
     endif else begin
        mu_1d_7 = environment.dipole 
        beta0   = mu_1d_7/sqrt(N_C(1d-7) + N_H(1d-7)) *debye
     endelse
  endelse

  line_present = STRCMP( tag_names, 'line', /FOLD_CASE )
  if total(line_present) EQ 0 then begin
     print, 'invalid environment : please specify the grain size distribution parameters by providing the line number (starting at 1) of Weingartner & Draine, 2001a'
     invalid = invalid + 1
  endif else begin
     line = environment.line - 1
  endelse

  if invalid NE 0 then begin
     return
  endif

  bc      = 1d-5*bc1e5_tab[line]
  alpha_g = alpha_tab[line]
  beta_g  = beta_tab[line]
  at_g    = 1d-4*at_tab[line]
  ac_g    = 1d-4*ac_tab[line]
  C_g     = C_tab[line]


; --- calculating the emissivity ---
  
                                ; Number of dipole moments used in the averaging
                                ; (I suppressed the keyword accuracy-Boost)
  Ndip = 20
  if (keyword_set(Ndipole)) then begin
     Ndip = Ndipole 
  endif

  ip = 2d/3d
  ip_present = STRCMP( tag_names, 'inplane', /FOLD_CASE )
  if total(ip_present) NE 0 then begin
     print, 'assuming that <mu_ip^2>/<mu^2> = '+strtrim(environment.inplane,2) + ' for disklike grains'
     ip = environment.inplane
  endif
  

  GHz   = 1d9
  numin = 0.5 *GHz   ; used to be 0.01 but that was not very useful
  numax = 500 *GHz
  Nnu   = 200
  if(keyword_set(min_freq)) then begin
     numin = min_freq *GHz
  endif
  if(keyword_set(max_freq)) then begin
     numax = max_freq *GHz
  endif
  if(keyword_set(n_freq)) then begin
     Nnu = n_freq
  endif

  nu_tab = makelogtab(numin, numax, Nnu)

  warning_phi_min = 0
  warning_phi_max = 0
  warning_psi_min = 0
  warning_psi_max = 0
  warning_Omega_min = 0
  warning_Omega_max = 0

  if (keyword_set(case1)) then begin
     print, 'assuming that disklike grains spin around their axis of greatest inertia'
     jnu_per_H = emissivity(env, beta0, ip, Ndip, nu_tab)
  endif else begin
     jnu_per_H = emissivity(env, beta0, ip, Ndip, nu_tab, /tumbling)
  endelse

; --- warnings in case interpolation in plasma drag went wrong.--------
  if (warning_phi_min EQ 1 ) then begin
     print, 'phi is less than phi_min in "little_gp_neutral_interpol". You should set phi_min to a smaller value in "set_up_gp_arrays" and run "compute_little_gp".'
  endif
  if (warning_phi_max EQ 1 ) then begin
     print, 'phi is greater than phi_max in "little_gp_neutral_interpol". You should set phi_max to a greater value in "set_up_gp_arrays" and run "compute_little_gp".'
  endif
  if (warning_psi_min EQ 1 ) then begin
     print, '|psi| is less than psi_min in "little_gp_charged_interpol". You should set psi_min to a smaller value in "set_up_gp_arrays" and run "compute_little_gp".'
  endif
  if (warning_psi_max EQ 1 ) then begin
     print, '|psi| is greater than psi_max in "little_gp_charged_interpol". You should set psi_max to a greater value in "set_up_gp_arrays" and run "compute_little_gp".'
  endif
  if (warning_Omega_min EQ 1 ) then begin
     print, 'Omega is less than Omega_min in "little_gp_charged_interpol" (plasmadrag.pro). You should set Omega_min to a smaller value in "set_up_gp_arrays" and run "compute_little_gp".'
  endif
  if (warning_Omega_max EQ 1 ) then begin
     print, 'Omega is greater than Omega_max in "little_gp_charged_interpol"(plasmadrag.pro). You should set Omega_max to a greater value in "set_up_gp_arrays" and run "compute_little_gp".'
  endif
; --------------------------------------------------------------------

  Jy = 1d-23
  result = dblarr(2, Nnu)
  if keyword_set(freefree) then begin
     result = dblarr(3, Nnu)
     result[2,*] = free_free(env, nu_tab)/Jy
  endif
  result[0,*] = nu_tab/GHz
  result[1,*] = jnu_per_H/Jy


  openw, 1, output_file

  printf, 1, ';', '============================ SPDUST.2 ==============================='
  printf, 1, ';'
  printf, 1, ';', '    Rotational emission from a population of spinning dust grains,'
  printf, 1, ';', '    as described by Ali-Haimoud, Hirata & Dickinson, 2009, MNRAS, 395, 1055'
  printf, 1, ';', '    and in Silsbee, Ali-Haimoud & Hirata, 2010'
  printf, 1, ';', '    for an environment with parameters given below. '
  printf, 1, ';'
  printf, 1, ';', '    nH = ', strtrim(env.nh,2), ' cm^-3'
  printf, 1, ';', '    T = ', strtrim(env.T,2), ' K' 
  printf, 1, ';', '    Chi = ', strtrim(env.Chi,2)
  printf, 1, ';', '    xh = ', strtrim(env.xh,2)
  printf, 1, ';', '    xC = ', strtrim(env.xC,2)
  printf, 1, ';', '    mu(1E-7 cm) = ',  strtrim(mu_1d_7,2), ' debye (this corresponds to beta = ', strtrim(beta0/debye,2), ' debye)'
  printf, 1, ';'
  if (keyword_set (case1)) then begin
     printf, 1, ';', '    Diklike grains are asssumed to spin around their axis of greatest inertia'
  endif else begin
     printf, 1, ';', '    Diklike grains are asssumed to be randomly oriented with respect to their angular momentum vector.'
  endelse
  printf, 1, ';'
  printf, 1, ';', '====================================================================='
  printf, 1, ';'
  if  keyword_set(freefree) then begin
     printf, 1, ';', 'nu(GHz)       ', 'j_nu/nH(Jy sr-1 cm2/H)     ', 'j_nu/nH (free-free) (Jy sr-1 cm2/H)' 
  endif else begin
     printf, 1, ';', 'nu(GHz)       ', 'j_nu/nH(Jy sr-1 cm2/H)     '
  endelse
  printf, 1, ';'
  printf, 1, result
  close, 1

  if (keyword_set(verbose)) then begin
      get_info, environment, case1 = case1
  endif


end

;----------------------------------------------------------------------------------------------------

pro get_info, environment, case1 = case1

; This routine computes only the rotational distribution function for
; the smallest grains, with a = 3.5E-8 cm, and returns some
; information about the charge, dominant excitation and damping
; mechanisms, and known uncertainties.

; Written February 2010

  common cgsconst  

  tag_names  = tag_names(environment)
  invalid_tags = invalid_tags(environment)
  invalid = invalid_tags.invalid
  env = invalid_tags.env
 
  dipole_present = STRCMP( tag_names, 'dipole', /FOLD_CASE )
  beta_present   = STRCMP( tag_names, 'beta', /FOLD_CASE )
  if (total(dipole_present) + total(beta_present)  EQ 0) then begin
     print, 'invalid environment : please specify the dipole moment of the grains (either dipole(a = 1E-7 cm) or beta, in debye)'
     invalid = invalid + 1
  endif else begin
     if total(beta_present)  NE 0 then begin
        beta0 = environment.beta *debye
        mu_1d_7 = sqrt(N_C(1d-7) + N_H(1d-7))* environment.beta
     endif else begin
        mu_1d_7 = environment.dipole 
        beta0   = mu_1d_7/sqrt(N_C(1d-7) + N_H(1d-7)) *debye
     endelse
  endelse

  if (invalid NE 0) then begin
     return
  endif

  ip = 2d/3d
  ip_present = STRCMP( tag_names, 'inplane', /FOLD_CASE )
  if total(ip_present) NE 0 then begin
     print, 'assuming that <mu_ip^2>/<mu^2> = '+strtrim(environment.inplane,2) + ' for disklike grains'
     ip = environment.inplane
  endif

    a = 3.5d-8
    fZ = charge_dist(env, a)
    Zaverage = total(fZ[0,*] *fZ[1,*], /double)
    Z2 = total(fZ[0,*]^2 * fZ[1,*], /double)
    mu_rms = rms_dipole(a, Z2, beta0)

    T       = env.T 
    Inertia = Inertia(a)
    tau_H      = tau_H(env, a) 
    Tev = Tev_effective(env, a)

    FGi  = FGi_averaged(env, a, Tev, mu_rms, fZ)  
    Fi   = FGi.Fi
    Gi   = FGi.Gi

    FGIR = FGIR_averaged(env, a, fZ)   
    FIR  = FGIR.FIR
    GIR  = FGIR.GIR
 
    FGpe = FGpe_averaged(env, a, fZ) 
    Fpe  = FGpe.Fpe
    Gpe  = FGpe.Gpe  

    GH2  = GH2(env, a) 

    omega_peak_th = sqrt(6d *k *T/Inertia)    ; peak of the spectrum if thermal rotation

    if (keyword_set(case1)) then begin
       tau_ed_inv = tau_ed_inv(env, a, mu_rms *sqrt(ip), mu_rms * sqrt(1d - ip))   
       FGn  = FGn_averaged(env, a, Tev, fZ)
       FGp = FGp_averaged(env, a, fZ, omega_peak_th, mu_rms *sqrt(ip), mu_rms * sqrt(1d - ip))
    endif else begin
       tau_ed_inv = tau_ed_inv(env, a, mu_rms *sqrt(ip), mu_rms * sqrt(1d - ip), /tumbling)   
       FGn  = FGn_averaged(env, a, Tev, fZ, /tumbling)
       FGp = FGp_averaged(env, a, fZ, omega_peak_th, mu_rms *sqrt(ip), mu_rms * sqrt(1d - ip), /tumbling)
    endelse

    Fn   = FGn.Fn
    Gn   = FGn.Gn
    Fp  = FGp.Fp
    Gp  = FGp.Gp
 
    F = Fn + Fi + FIR + Fp + Fpe
    G = Gn + Gi + GIR + Gp + Gpe + GH2
    xi = 8d * G * tau_H/F^2 * tau_ed_inv

   ;--- Checking whether case 1 may actually be a better description
   
    hnu_tab     = makelogtab(1d-2, 13.6d, 500)
    Dnu_over_nu = DX_over_X(1d-2, 13.6d, 500)
    Qabs     = Qabs(a, 1, hnu_tab)
    nu_uisrf =  nu_uisrf(hnu_tab)
    abs_time = 1d/(env.Chi *pi *a^2 *total( Qabs *nu_uisrf/(hnu_tab*eV), /double) *Dnu_over_nu *c)
    Lchange_time = min([tau_H/F, sqrt(tau_H/(G *tau_ed_inv))])
   
    print, 'All values are provided for a grain volume-equivalent radius of 3.5E-8 cm'
    print, 'Average grain charge <Z> = ' + strtrim(Zaverage, 2)
    print, '---- Damping rates : ---- '
    print, 'Fn  = '+ strtrim(Fn, 2)
    print, 'Fi  = '+ strtrim(Fi, 2)
    print, 'FIR = '+ strtrim(FIR,2)
    print, 'Fp  = '+ strtrim(Fp, 2)
    print, 'Fpe = '+ strtrim(Fpe, 2)
    print, '---- Excitation rates : ---- '
    print, 'Gn  = '+ strtrim(Gn, 2)
    print, 'Gi  = '+ strtrim(Gi, 2)
    print, 'GIR = '+ strtrim(GIR,2)
    print, 'Gp  = '+ strtrim(Gp, 2)
    print, 'Gpe = '+ strtrim(Gpe, 2)
    print, 'GH2 = '+ strtrim(GH2, 2)
    print, '---- Non gaussianity parameter : ---- '
    print, 'xi = 8 G tau_H/(F^2 tau_ed) = ' + strtrim(xi, 2)
    print, '---- Timescales ----'
    print, 'Time between UV photon absorptions: tau_UV = '+ strtrim(abs_time, 2) + ' sec'
    print, 'Time to change angular momentum: tau_L = '+ strtrim(Lchange_time, 2) + 'sec'
    print, 'If tau_UV < tau_L, case 2 is appropriate. Otherwise, case 1 may be more appropriate'
    print, '---- Checking known uncertainties ----'
    print, ' - Evaporation temperature: R_{col/abs}/Nsites = ' + strtrim(Tev_effective(env, a, /get_info), 2)
    print, 'If this ratio is close to 1 then the evaporation temperature is uncertain'
    Tev_present = STRCMP( tag_names, 'Tev', /FOLD_CASE )
    if total(Tev_present) NE 0 then begin
       print, 'But you chose to use a constant evaporation temperature anyway'
    endif
    if NOT(keyword_set(case1)) then begin
       print, ' - Interpolation function for Fn^(in)(Zg <> 0) :' 
       print, 'Fn/F = ' + strtrim(Fn/F, 2)
       print, 'xi   = ' + strtrim(xi, 2)
       charged_fraction = 1d - fZ[1,0]
       print, 'fraction of charged grains: '+ strtrim(charged_fraction, 2) 
       ratio = 1.5d-8 * (env.T/8d3)^(-0.25d)/acx(a)
       print, 'rc/ acx = ' + strtrim(ratio, 2) 
       print, 'IF Fn is the dominant damping'
       print, '   and xi << 1 '
       print, '   and grains are mostly charged'
       print, '   and rc/acx is close to 1'
       print, 'THEN there could be a significant uncertainty in the spectrum'
    endif

end
