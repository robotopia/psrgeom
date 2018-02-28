# Variables
x = var('x')
y = var('y')
z = var('z')
r = var('r')
rL = var('rL',latex_name='r_L')
al = var('al',latex_name='\\alpha')
ph = var('ph',latex_name='\\varphi')
m0 = var('m0',latex_name='m_0')
Om = var('Om',latex_name='\\Omega')
t = var('t')
th = var('th',latex_name='\\theta')
c = var('c')
rmin = var('rmin',latex_name='r_\\ast')

assume(x,'real')
assume(y,'real')
assume(z,'real')
assume(r,'real')
assume(rL,'real')
assume(al,'real')
assume(ph,'real')
assume(m0,'real')
assume(Om,'real')
assume(t,'real')
assume(th,'real')
assume(c,'real')
assume(rmin,'real')

k = 1/rL
rhat = vector([x/r,y/r,z/r])
#ph = (r-rmin)/rL # leave commented to keep the "ph" as an abbreviation for "(r-rmin)/rL"
m0 = 1 # uncomment to make the resulting expressions slightly easier to read

xhat = vector([1,0,0])
yhat = vector([0,1,0])
zhat = vector([0,0,1])

mperp = m0*sin(al)
mz    = m0*cos(al)

Bpart1 = mperp*((xhat - rhat*(rhat.dot_product(xhat)))*(k^2*cos(ph)/r) + (3*rhat*(xhat.dot_product(rhat)) - xhat)*(k*r*sin(ph) + cos(ph))/r^3)
Bpart2 = mperp*((rhat*(rhat.dot_product(yhat)) - yhat)*(k^2*sin(ph)/r) + (3*rhat*(yhat.dot_product(rhat)) - yhat)*(k*r*cos(ph) - sin(ph))/r^3)
Bpart3 = 1/r^3 * (3*rhat*(rhat.dot_product(zhat)) - zhat)

B1 = Bpart1+Bpart2+Bpart3
B1x = B1[0]
B1y = B1[1]
B1z = B1[2]

A1 =  1/r^5
A2 =  sin(al)*cos(ph)/r^3/rL^2
A3 =  sin(al)*cos(ph)/r^5    + sin(al)*sin(ph)/r^4/rL
A4 = -sin(al)*sin(ph)/r^3/rL^2
A5 =  sin(al)*cos(ph)/r^4/rL - sin(al)*sin(ph)/r^5

B2x = A1*(3*x*z)     + A2*(r^2-x^2) + A3*(3*x^2-r^2) - A4*(x*y)     + A5*(3*x*y)
B2y = A1*(3*y*z)     - A2*(x*y)     + A3*(3*x*y)     + A4*(r^2-y^2) + A5*(3*y^2-r^2)
B2z = A1*(3*z^2-r^2) - A2*(x*z)     + A3*(3*x*z)     - A4*(y*z)     + A5*(3*y*z)

# Now, we have
# B1x == B2x
# B1y == B2y
# B1z == B2z
