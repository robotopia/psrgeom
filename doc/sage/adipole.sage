# Calculate acceleration field

load('vdipole.sage')

Apos = diff(Vpos,t) + Vpos[0]*diff(Vpos,x) + Vpos[1]*diff(Vpos,y) + Vpos[2]*diff(Vpos,z)
Aneg = diff(Vneg,t) + Vneg[0]*diff(Vneg,x) + Vneg[1]*diff(Vneg,y) + Vneg[2]*diff(Vneg,z)

Vposlen2 = var('Vposlen2')
assume(Vposlen2, 'real')
assume(Vposlen2 > 0)
Vposlen2 = Vpos.dot_product(Vpos)
Vposlen = sqrt(Vposlen2)

numvec = Vpos.cross_product(Apos)
numlen2 = var('numlen2')
assume(numlen2, 'real')
assume(numlen2 > 0)
numlen2 = numvec.dot_product(numvec)
num = sqrt(numlen2)

den = Vposlen^3

kappa = num / den
