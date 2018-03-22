# Sage program for verifying the calculation of the acceleration field

###################################################################

x = var('x')
y = var('y')
z = var('z')
t = var('t')

c = 299792458
P = 0.5
al = n(30.0*pi/180.0)
rp = 1e4
rL = n(c*P/(2*pi))
Om = n(2*pi/P)

r = sqrt(x^2 + y^2 + z^2)
De = (r-rp)/rL


A1 =  1/r^5
A2 =  sin(al)*cos(De)/r^3/rL^2
A3 =  sin(al)*cos(De)/r^5    + sin(al)*sin(De)/r^4/rL
A4 = -sin(al)*sin(De)/r^3/rL^2
A5 =  sin(al)*cos(De)/r^4/rL - sin(al)*sin(De)/r^5

xp = x*cos(Om*t) - y*sin(Om*t)
yp = x*sin(Om*t) + y*cos(Om*t)
zp = z

Bx = A1*(3*xp*zp)    + A2*(r^2-xp^2) + A3*(3*xp^2-r^2) - A4*(xp*yp)    + A5*(3*xp*yp)
By = A1*(3*yp*zp)    - A2*(xp*yp)    + A3*(3*xp*yp)    + A4*(r^2-yp^2) + A5*(3*yp^2-r^2)
Bz = A1*(3*zp^2-r^2) - A2*(xp*zp)    + A3*(3*xp*zp)    - A4*(yp*zp)    + A5*(3*yp*zp)

Blen = sqrt(Bx^2 + By^2 + Bz^2)

Bnormx = Bx / Blen
Bnormy = By / Blen
Bnormz = Bz / Blen

phx = -y
phy =  x
phz =  0

rho = sqrt(x^2 + y^2)
ph_dot_Bnorm = phx*Bnormx + phy*Bnormy + phz*Bnormz

VBpos = -Om*ph_dot_Bnorm + sqrt(Om^2*ph_dot_Bnorm^2 - (rho^2*Om^2 - c^2))
VBneg = -Om*ph_dot_Bnorm - sqrt(Om^2*ph_dot_Bnorm^2 - (rho^2*Om^2 - c^2))

VxP = Om*phx + VBpos*Bnormx
VyP = Om*phy + VBpos*Bnormy
VzP = Om*phz + VBpos*Bnormz

VxN = Om*phx + VBneg*Bnormx
VyN = Om*phy + VBneg*Bnormy
VzN = Om*phz + VBneg*Bnormz

VlenP = sqrt(VxP^2 + VyP^2 + VzP^2)
VlenN = sqrt(VxN^2 + VyN^2 + VzN^2)

VPnormx = VxP / VlenP
VPnormy = VyP / VlenP
VPnormz = VzP / VlenP

VNnormx = VxN / VlenN
VNnormy = VyN / VlenN
VNnormz = VzN / VlenN

AxP = diff(VxP,t) + VxP*diff(VxP,x) + VyP*diff(VxP,y) + VzP*diff(VxP,z)
AyP = diff(VyP,t) + VxP*diff(VyP,x) + VyP*diff(VyP,y) + VzP*diff(VyP,z)
AzP = diff(VzP,t) + VxP*diff(VzP,x) + VyP*diff(VzP,y) + VzP*diff(VzP,z)

AxN = diff(VxN,t) + VxN*diff(VxN,x) + VyN*diff(VxN,y) + VzN*diff(VxN,z)
AyN = diff(VyN,t) + VxN*diff(VyN,x) + VyN*diff(VyN,y) + VzN*diff(VyN,z)
AzN = diff(VzN,t) + VxN*diff(VzN,x) + VyN*diff(VzN,y) + VzN*diff(VzN,z)

AlenP = sqrt(AxP^2 + AyP^2 + AzP^2)
AlenN = sqrt(AxN^2 + AyN^2 + AzN^2)

APnormx = AxP / AlenP
APnormy = AyP / AlenP
APnormz = AzP / AlenP

ANnormx = AxN / AlenN
ANnormy = AyN / AlenN
ANnormz = AzN / AlenN

VP_dot_AP = VPnormx*APnormx + VPnormy*APnormy + VPnormz*APnormz
VN_dot_AN = VNnormx*ANnormx + VNnormy*ANnormy + VNnormz*ANnormz
