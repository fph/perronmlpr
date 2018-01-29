# https://ask.sagemath.org/question/11070/find-algebraic-solutions-to-system-of-polynomial-equations/

import numpy as np
import time

starttime = time.time()

R = matrix(np.loadtxt('R32.txt'))
R = R.change_ring(QQ)
n = R.dimensions()[0]
v = vector([1/n]*n)
#x = vector(var(', '.join('x' + str(i) for i in range(1,n+1))))
Ring = PolynomialRing(QQ, 'x', n, order='lex'); x = vector(Ring.gens())

kronx = []
for xi in x:
	kronx = kronx + list(xi*x)

kronx = vector(kronx)


#def find_zeros(alpha):
"""
Finds all values of x[-1] for which the MLPR equation is satisfied.
TODO: Does not attempt to verify that these have the prescribed signs (for now).
"""

alpha = 0.7

eqns = 	alpha * R * kronx + (1-alpha)*v - x

#numsol = solve(list(eqns), x)

J = Ring.ideal(list(eqns))
gbasis = J.groebner_basis()
eq = gbasis[-1]
y = var('y')
poly = sum(eq.monomial_coefficient(x[-1]^i)* y^i for i in range(0, 2^n+1))
numsol = poly.roots(ring=RR)

for sol, mult in numsol:
	for J in gbasis[:-1]:
		solve()

#	return tuple(a[0] for a in numsol)


# zz = find_zeros(alpha = 99/100)

endtime = time.time()
print "Time elapsed: %s" % (endtime - starttime)
