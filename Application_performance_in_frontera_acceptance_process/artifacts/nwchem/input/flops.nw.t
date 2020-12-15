restart gg-stack
memory stack 2000 mb heap 200 mb global 1500 mb noverify
 
title "GG stack"

scf
  noprint "final vectors analysis"
  maxiter 100
end

basis spherical noprint
  *  library aug-cc-pvtz
end

ccsd
  freeze core atomic
  thresh 1000.0
end

set ccsd:converged T

task ccsd(t) energy
