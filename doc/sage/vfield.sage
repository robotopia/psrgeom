load("bfield.sage")

# Variables
det = var('det')
Bnormx = var('Bnormx',latex_name='\\hat{B}_x')
Bnormy = var('Bnormy',latex_name='\\hat{B}_y')
Bnormz = var('Bnormz',latex_name='\\hat{B}_z')
assume(det,'real')

Ph = vector([-y,x,0])
rho = sqrt(x^2+y^2)

Bnorm = vector([Bnormx, Bnormy, Bnormy])

Ph_dot_Bnorm = Ph.dot_product(Bnorm)


det = sqrt(Om^2*Ph_dot_Bnorm^2 - 4*(rho^2*Om^2-c^2))

VBpos = (-Om*Ph_dot_Bnorm + det)/2
VBneg = (-Om*Ph_dot_Bnorm - det)/2

Vpos = VBpos*Bnorm + Om*Ph
Vneg = VBneg*Bnorm + Om*Ph


