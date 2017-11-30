;========================== SPDUST.2 =========================================;
;                                                                             ;
;                     COMPILE_SPDUST.pro                                      ;
;                                                                             ;
; IDL batch script to compile all subroutines for SPDUST                      ;
; and set the Data Files directory.                                           ;
; Simply type "@compile_spdust" at the IDL command line.                      ;
; Sets up the path of the data files, then reads out table 1 from             ;
; WD01a to obtain the parameter of the size distribution function,            ;
; then compiles all the IDL routines, and                                     ;
; finally calculates arrays needed for computing spinning dust                ;
; emissivities.                                                               ;
;                                                                             ;             
; Author : Yacine Ali-Haimoud        yacine@tapir.caltech.edu                 ; 
;                                                                             ;
; Revision history : Written December 2008                                    ;
;                    March 2010 (new release): no changes                     ;
;==============================================================================

common path, SpDust_data_dir
common size_dist_arrays, bc1e5_tab, alpha_tab, beta_tab, at_tab, ac_tab, C_tab
common gp_arrays, psi_min, psi_max, Npsi, phi_min, phi_max, Nphi, Omega_min, Omega_max, NOmega, gp_pos, gp_neg, gp_neutral
common IR_arrays, a_min, a_max, Na, chi_min, chi_max, Nchi, FIR_integral_charged, GIR_integral_charged, FIR_integral_neutral,GIR_integral_neutral, Tev_tab  
common free_free_tab, b_tab

;--------------- MODIFY THE PATH OF THE DATA FILES HERE ---------------
   SpDust_data_dir = '/Applications/itt/idl/lib/SPDUST.2/Data_Files/'
;----------------------------------------------------------------------

size_dist_file  = SpDust_data_dir + 'sizedists_table1.out'
readcol, size_dist_file, bc1e5_tab, alpha_tab, beta_tab, at_tab, ac_tab, C_tab, format = 'X, D, X, D, D, D, D, D', /silent

.run subroutines
.run grain_properties
.run charge_dist
.run infrared
.run collisions
.run plasmadrag
.run H2_photoemission
.run emissivity
.run SpDust

physconst             ; storing the usual physical constants
parameters            ; storing dust grain parameters
Jpeisrf_calc          ; computing the photoelectric emission rate 

; --- settting up plasma drag arrays ---

set_up_gp_arrays     

gp_pos     = dblarr(Npsi, NOmega)
gp_neg     = dblarr(Npsi, NOmega) 
gp_neutral = dblarr(Nphi, NOmega)

openr, 1, SpDust_data_dir + 'gp_pos_'+ strtrim(Npsi, 2)+ 'psi_' + strtrim(NOmega, 2) + 'Omega' 
readf, 1, gp_pos
close, 1

openr, 1, SpDust_data_dir + 'gp_neg_'+ strtrim(Npsi, 2)+ 'psi_' + strtrim(NOmega, 2) + 'Omega' 
readf, 1, gp_neg
close, 1

openr, 1, SpDust_data_dir + 'gp_neutral_'+ strtrim(Nphi, 2)+ 'phi_' + strtrim(NOmega, 2) + 'Omega'
readf, 1, gp_neutral
close, 1

print, 'gp_arrays stored'


; --- setting up infrared emission arrays ---

set_up_IR_arrays
      
FIR_integral_charged = dblarr(Na, Nchi)
GIR_integral_charged = dblarr(Na, Nchi)
FIR_integral_neutral = dblarr(Na, Nchi)
GIR_integral_neutral = dblarr(Na, Nchi)

openr, 1, SpDust_data_dir+'FIR_integral_charged_'+strtrim(Na, 2)+'a_'+strtrim(Nchi, 2)+'chi'
readf, 1, FIR_integral_charged
close, 1

openr, 1, SpDust_data_dir+'GIR_integral_charged_'+strtrim(Na, 2)+'a_'+strtrim(Nchi, 2)+'chi'
readf, 1, GIR_integral_charged
close, 1

openr, 1, SpDust_data_dir+'FIR_integral_neutral_'+strtrim(Na, 2)+'a_'+strtrim(Nchi, 2)+'chi'
readf, 1, FIR_integral_neutral
close, 1

openr, 1, SpDust_data_dir+'GIR_integral_neutral_'+strtrim(Na, 2)+'a_'+strtrim(Nchi, 2)+'chi'
readf, 1, GIR_integral_neutral
close, 1

print, 'IR integrals stored'

; --- setting up T_ev array ---

Tev_tab = dblarr(Na, Nchi)

openr, 1, SpDust_data_dir+'Tev_'+strtrim(Na,2)+'a_' + strtrim(Nchi, 2)+'chi'
readf, 1, Tev_tab
close, 1

print, 'Tev array strored'

; --- Reading gaunt factor for free-free emission (Sutherland, 1998) ---

read_gaunt_factor
