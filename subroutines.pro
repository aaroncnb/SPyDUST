;========================== SPDUST.2 ====================================;
;                                                                        ;
;                       SUBROUTINES.PRO                                  ;  
;                                                                        ;
; Miscelanious functions for array creation and interpolation            ;
;                                                                        ;             
; Author : Yacine Ali-Haimoud        yacine@tapir.caltech.edu            ;
;                                                                        ;
; Revision history : Written December 2008                               ;
;========================================================================;



function maketab, xmin, xmax, Npt

; ---> Creates a linearly spaced table

  return, xmin + (xmax - xmin)/Npt *(0.5d + dindgen(Npt))
end

;----------------------------------------------------------------------------

function makelogtab, xmin, xmax, Npt

; ---> Creates a logarithmically spaced table

  return, xmin *exp( 1d/Npt *alog(xmax/xmin) *(0.5d + dindgen(Npt)) )
end

;----------------------------------------------------------------------------

function DX_over_X, xmin, xmax, Npt

; ---> Returns Dx/x for a logarithmically spaced table created by makelogtab

  return, exp(0.5d /Npt *alog(xmax/xmin))- exp(-0.5d /Npt *alog(xmax/xmin))
end

;----------------------------------------------------------------------------

function biinterplog, f, X, Y, Xnew, Ynew

; ---> f is assumed to be an array of Nx*Ny elements
; interpolates f given at X, Y at Xnew, Ynew.
; X and Y are logarithmically spaced, f is linearly spaced

Nx = N_elements(X)
Ny = N_elements(Y)

Xind = dindgen(Nx)
Yind = dindgen(Ny)

Xindnew = interpol(Xind, alog(X), alog(Xnew))
Yindnew = interpol(Yind, alog(Y), alog(Ynew))

return, interpolate(f, Xindnew, Yindnew, /grid)

end

;----------------------------------------------------------------------------

function log_biinterplog, f, X, Y, Xnew, Ynew

; ---> Same as above, with f logarithmically sampled

return, exp(biinterplog( alog(f), X, Y, Xnew, Ynew) )

end

;----------------------------------------------------------------------------

