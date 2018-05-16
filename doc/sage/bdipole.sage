# Create variables
x = var('x')    # Spatial coordinate 1
y = var('y')    # Spatial coordinate 2
z = var('z')    # Spatial coordinate 3
r = var('r')    # sqrt(x^2 + y^2 + z^2)
al = var('al')  # angle between rotation and magnetic axes
m0 = var('m0')  # "mu nought", the magnetic field strength

m0 = 1  # comment this if you want expressions to contain m0 explicitly

rhat = vector([x/r, y/r, z/r])

mu = m0*vector([sin(al), 0, cos(al)])

B = 1/r^3 * (3*mu.dot_product(rhat)*rhat - mu)
