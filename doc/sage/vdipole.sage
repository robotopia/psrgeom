# Get dipole magnetic field
load('bdipole.sage')

det = var('det')
rho = var('rho')
c = var('c')

assume(det,'real')
assume(rho,'real')
assume(c,'real')
assume(c > 0)

Blen = B.norm()
Bnorm = B / Blen
Ph = vector([-y,x,0])
rho = sqrt(x^2+y^2)
Ph_dot_Bnorm = Ph.dot_product(Bnorm)
det = sqrt(Om^2*Ph_dot_Bnorm^2 - 4*(rho^2*Om^2-c^2))
VBpos = (-Om*Ph_dot_Bnorm + det)/2
VBneg = (-Om*Ph_dot_Bnorm - det)/2
Vpos = VBpos*Bnorm + Om*Ph
Vneg = VBneg*Bnorm + Om*Ph


