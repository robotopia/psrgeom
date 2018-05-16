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
