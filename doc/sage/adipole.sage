load('bdipole.sage')

dr_dx = x/r
dr_dy = y/r
dr_dz = z/r

dBdx = diff(B,x) + diff(B,r)*dr_dx
dBdy = diff(B,y) + diff(B,r)*dr_dy
dBdz = diff(B,z) + diff(B,r)*dr_dz


# The following is what's written in my C code:
r7i = 1/r^7
rr = r^2

xx = x*x
yy = y*y
zz = z*z

xy = x*y
xz = x*z
yz = y*z

dBdx_x = 3*r7i*(
                 z*cos(al)*(  rr - 5*xx) +
                 x*sin(al)*(3*rr - 5*xx)
               )
dBdy_x = 3*r7i*y*(
                 cos(al)*(   - 5*xz) +
                 sin(al)*(rr - 5*xx)
               )
dBdz_x = 3*r7i*(
                 x*cos(al)*(rr - 5*zz) +
                 z*sin(al)*(rr - 5*xx)
               )
dBdx_y = dBdy_x
dBdy_y = 3*r7i*(
                 z*cos(al)*(rr - 5*yy) +
                 x*sin(al)*(rr - 5*yy)
               )
dBdz_y = 3*r7i*y*(
                 cos(al)*(rr - 5*zz) +
                 sin(al)*(   - 5*xz)
               )
dBdx_z = dBdz_x
dBdy_z = dBdz_y
dBdz_z = 3*r7i*(
                 z*cos(al)*(3*rr - 5*zz) +
                 x*sin(al)*(  rr - 5*zz)
               )

# Here is the solution for A for an aligned rotator (considering only the xz-plane)
load('vdipole.sage')

Vx = 3*x*z*sqrt(c^2-Om^2*x^2)/sqrt((x^2+4*z^2)*(x^2+z^2))
Vy = Om*x
Vz = (2*z^2 - x^2)*sqrt(c^2-Om^2*x^2)/sqrt((x^2+4*z^2)*(x^2+z^2))

V = vector([Vx, Vy, Vz])
dV_dx = diff(V,x)
dV_dy = diff(V,y)
dV_dz = diff(V,z)

A = Vx*dV_dx + Vy*dV_dy + Vz*dV_dz

Ax = (A[0].simplify_trig().factor())(x=R*sin(th)^3,z=R*sin(th)^2*cos(th)).simplify_trig().factor()
Ay = 3*Om*cos(th)*sin(th)*sqrt(c^2-Om^2*R^2*sin(th)^6)/sqrt(3*cos(th)^2+1)
Az = (A[2].simplify_trig().factor())(x=R*sin(th)^3,z=R*sin(th)^2*cos(th)).simplify_trig().factor()
Asimp = vector([Ax,Ay,Az])

# Re-acquire the velocity components from vdipole (where they've been simplified)
load('vdipole.sage')

# Calculate the curvature and find its minima
VxA = Vsimp.cross_product(Asimp)
VxAnorm = sqrt(VxA[0]^2 + VxA[1]^2 + VxA[2]^2)

# |V| = c, so but this is just a constant, irrelevant to finding the minimum, so
# we'll leave it out in the calculation of kappa
kappa = VxAnorm.simplify_trig()
dk_dth = (diff(kappa,th)).simplify_trig().factor()

num = dk_dth.numerator()
num2 = (num^2).factor().numerator()/cos(th)^2/9/R^2

# Isolate the only factor with interesting roots

num3 = 12*Om^4*R^4*sin(th)^16 - 40*Om^4*R^4*sin(th)^14 + 40*Om^4*R^4*sin(th)^12 - 33*Om^2*R^2*c^2*sin(th)^10 + 110*Om^2*R^2*c^2*sin(th)^8 - 128*Om^2*R^2*c^2*sin(th)^6 + 32*Om^2*R^2*c^2*sin(th)^4 - 6*c^4*sin(th)^4 + 20*c^4*sin(th)^2 - 8*c^4

# Express the field line constant as a fraction of the light cylinder radius
num4 = num3(R=r*c/Om).factor()/c^4

# Get the terms in r^2
Qa = num4.collect(r).operands()[0]/r^4
Qb = num4.collect(r).operands()[2]/r^2
Qc = num4.collect(r).operands()[1] + num4.collect(r).operands()[3] + num4.collect(r).operands()[4]
