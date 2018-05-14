# Variables
x = var('x')
y = var('y')
z = var('z')
r = var('r')
al = var('al',latex_name='\\alpha')
m0 = var('m0',latex_name='m_0')
Om = var('Om',latex_name='\\Omega')
t = var('t')

assume(x,'real')
assume(y,'real')
assume(z,'real')
assume(r,'real')
assume(al,'real')
assume(m0,'real')
assume(Om,'real')
assume(t,'real')

xhat = vector([1,0,0])
yhat = vector([0,1,0])
zhat = vector([0,0,1])

# Construct the magnetic axis, inclined by angle alpha and rotated by angle Om*t
m = m0 * vector([sin(al)*cos(Om*t), sin(al)*sin(Om*t), cos(al)])

rhat = vector([x/r,y/r,z/r])

# Now put in the standard dipole field aligned along the z-axis
B = 1/r^3 * (3*m.dot_product(rhat)*rhat - m)
