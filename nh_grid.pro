;========================================================================;
;                          NH_GRID.PRO                                   ;  
; Example IDL routine to run SPDUST on a grid of densities               ;     
;                                                                        ;
; Last modified February 2010                                            ;             
; Author : Yacine Ali-Haimoud        yacine@tapir.caltech.edu            ;
;========================================================================;

pro nh_grid

; --- input parameters (the nh tag has to be defined but can be
;     assigned any value for now) ---

  input_params = {nh: 0d, T: 1d2, chi: 1d, xh: 1d-3, xc: 3d-4, y: 0d, gamma: 0d, dipole: 9.3d, line: 7}

; --- array of number densities ---

  nh_array = [1d, 3d, 10d, 30d, 100d]

; --- Path of the output files (we will add a suffix to distinguish
;     between them) ---

  output_file = '~/Desktop/cnm_nh_'

  for i = 0, 4 do begin
     nh = nh_array[i]
     input_params.nh = nh
     SPDUST, input_params, output_file + strtrim(floor(nh),2)
     print, nh
  endfor

; this will produce 5 files '~/Desktop/cnm_nh_1',
; '~/Desktop/cnm_nh_3', etc... containing the spinning dust
; emissivities for each density.

;--- now read the output files and plot them ---

  readcol, output_file + '1', nu, jnu
  plot, nu, jnu, /xlo, /ylog, xrange = [2d, 1d3], xstyle = 1, yrange = [5d-21, 5d-17], ystyle = 1, linestyle = 1

  for i = 1, 4 do begin
     nh = nh_array[i]
     readcol, output_file + strtrim(floor(nh),2), nu, jnu
     oplot, nu, jnu, linestyle = i+1
  endfor


end
