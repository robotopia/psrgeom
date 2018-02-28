# Variables
x = var('x')
y = var('y')
z = var('z')
t = var('t')

c = var('c')
Om  = var('Om',latex_name='\\Omega')
det = var('det')
assume(det, 'real')

function('Bnormx')
function('Bnormy')
function('Bnormz')

Bnorm = vector([Bnormx(x,y,z,t),Bnormy(x,y,z,t),Bnormz(x,y,z,t)])

Ph = vector([-y,x,0])

PdotB = Ph.dot_product(Bnorm)

# Check:
#diff(PdotB,t) == Ph.dot_product(diff(Bnorm,t))
#diff(PdotB,x) == Ph.dot_product(diff(Bnorm,x)) + Bnormy(x,y,z,t)
#diff(PdotB,y) == Ph.dot_product(diff(Bnorm,y)) - Bnormx(x,y,z,t)
#diff(PdotB,z) == Ph.dot_product(diff(Bnorm,z))

rho   = sqrt(x^2+y^2)
det   = sqrt(Om^2*PdotB^2 - 4*(rho^2*Om^2-c^2))
VB    = (-Om*PdotB + det)/2

dVB_dt = diff(VB,t)

chi = Om^2 / (2*sqrt(Om^2*PdotB^2 - 4*(rho^2*Om^2 - c^2)))

# Check:
#diff(VB,t) == -Om*diff(PdotB,t)/2 + chi*(PdotB*diff(PdotB,t))
#diff(VB,x) == -Om*diff(PdotB,x)/2 + chi*(PdotB*diff(PdotB,x) - 4*x)
#diff(VB,y) == -Om*diff(PdotB,y)/2 + chi*(PdotB*diff(PdotB,y) - 4*y)
#diff(VB,z) == -Om*diff(PdotB,z)/2 + chi*(PdotB*diff(PdotB,z))

function('Bx')
function('By')
function('Bz')

B = vector([Bx(x,y,z,t),By(x,y,z,t),Bz(x,y,z,t)])
Bl = sqrt(B[0]^2 + B[1]^2 + B[2]^2)
Bn = B/Bl

# Check:
#(diff(Bn,t)[0]).factor() == ((1/Bl * (diff(B,t) - (Bn.dot_product(diff(B,t))*Bn)))[0]).factor()
#(diff(Bn,x)[0]).factor() == ((1/Bl * (diff(B,x) - (Bn.dot_product(diff(B,x))*Bn)))[0]).factor()
#(diff(Bn,y)[0]).factor() == ((1/Bl * (diff(B,y) - (Bn.dot_product(diff(B,y))*Bn)))[0]).factor()
#(diff(Bn,z)[0]).factor() == ((1/Bl * (diff(B,z) - (Bn.dot_product(diff(B,z))*Bn)))[0]).factor()

xp(x,y,z,t) = x*cos(Om*t) - y*sin(Om*t)
yp(x,y,z,t) = x*sin(Om*t) + y*cos(Om*t)
zp(x,y,z,t) = z

r = sqrt(x^2 + y^2 + z^2)
rL = var('rL',latex_name='r_L')
al = var('al',latex_name='\\alpha')
ph = var('ph',latex_name='\\varphi')

A1 =  1/r^5
A2 =  sin(al)*cos(ph)/r^3/rL^2
A3 =  sin(al)*cos(ph)/r^5    + sin(al)*sin(ph)/r^4/rL
A4 = -sin(al)*sin(ph)/r^3/rL^2
A5 =  sin(al)*cos(ph)/r^4/rL - sin(al)*sin(ph)/r^5

# Check:
#diff(A1,x) == -5*x/r^7
#diff(A2,x) == -3*x*sin(al)*cos(ph)/r^5/rL^2
#diff(A3,x) == -x*sin(al)/r^6 * (5*cos(ph)/r + 4*sin(ph)/rL)
#diff(A4,x) ==  3*x*sin(al)*sin(ph)/r^5/rL^2
#diff(A5,x) == -x*sin(al)/r^6 * (4*cos(ph)/rL - 5*sin(ph)/r)

function("A1 A2 A3 A4 A5")

B2x = A1(x,y,z)*(3*x*z)     + A2(x,y,z)*(r^2-x^2) + A3(x,y,z)*(3*x^2-r^2) - A4(x,y,z)*(x*y)     + A5(x,y,z)*(3*x*y)
B2y = A1(x,y,z)*(3*y*z)     - A2(x,y,z)*(x*y)     + A3(x,y,z)*(3*x*y)     + A4(x,y,z)*(r^2-y^2) + A5(x,y,z)*(3*y^2-r^2)
B2z = A1(x,y,z)*(3*z^2-r^2) - A2(x,y,z)*(x*z)     + A3(x,y,z)*(3*x*z)     - A4(x,y,z)*(y*z)     + A5(x,y,z)*(3*y*z)

# Check:
#diff(Bx2,x) == diff(A1(x,y,z),x)*3*x*z + A1(x,y,z)*3*z + diff(A2(x,y,z),x)*(r^2-x^2) ...

# Check:
#diff(B2x,x) == diff(A

