"""Exact checks for the degree-one quaternion-order example.

All polynomial computations are over QQ. The script checks the three affine
smoothness ideals, the contact-conic identity, the full Hilbert polynomials,
the degree-two multiplication matrix, and the numerical Euler calculations.
The categorical interpretation and the nontriviality of the residue are proved
in the accompanying mathematical argument, not by a finite experiment here.
"""

import sympy as sp

u, v, w = sp.symbols('u v w')
q = sp.Matrix([[u, v, w**2], [v, w, u**2], [w**2, u**2, 2*v**3]])
delta = sp.expand(q.det())
print('Coefficient field: QQ')
print('Quadratic-form matrix:', q)
print('Quintic discriminant:', delta)
assert delta == -u**5 + 2*u**2*v*w**2 + 2*u*v**3*w - 2*v**5 - w**5

for chart_variable, affine_variables in [(w, (u,v)), (v, (u,w)), (u, (v,w))]:
    affine_equation = sp.expand(delta.subs(chart_variable, 1))
    singularity_ideal = [affine_equation] + [sp.diff(affine_equation, z) for z in affine_variables]
    basis = sp.groebner(singularity_ideal, *affine_variables, domain=sp.QQ)
    assert list(basis) == [1]
    print(f'Chart {chart_variable}=1 singularity ideal Groebner basis:', list(basis))

s, t = sp.symbols('s t')
restriction = sp.factor(delta.subs({u:s**2, v:s*t, w:t**2}))
assert sp.expand(restriction + (s**5-t**5)**2) == 0
assert sp.gcd(s**5-1, 5*s**4) == 1
print('On contact conic (u,v,w)=(s^2,st,t^2):', restriction)
print('The five contact roots are distinct; no root lies at t=0.')

m, n = sp.symbols('m n', integer=True)
chi_line = lambda degree: (degree+2)*(degree+1)/2
even = sp.expand(chi_line(m)+chi_line(m-1)+2*chi_line(m-2))
odd = sp.expand(2*chi_line(m)+chi_line(m-1)+chi_line(m-2))
assert even == 2*m**2+m+1
assert odd == 2*m**2+3*m+2
print('Even Hilbert polynomial:', even)
print('Odd Hilbert polynomial:', odd)
z = sp.symbols('z')
hilbert = 1/(1-z) + z/(1-z)**3
weighted = (1-z**6)/((1-z)**2*(1-z**2)*(1-z**3))
assert sp.cancel(hilbert-weighted) == 0
print('Full Hilbert series identity:', sp.factor(weighted))

# Ordered input: X^2, X*Y, Y*X, Y^2.
# Ordered output: u, v, w, e_X wedge e_Y.
multiplication = sp.Matrix([[1,0,0,0], [0,1,1,0], [0,0,0,1], [0,1,-1,0]])
assert multiplication.det() != 0
print('Degree-two multiplication matrix:', multiplication)
print('Degree-two multiplication determinant:', multiplication.det())

p = lambda k: 1+sp.Rational(1,2)*k*(k+1)
gram = sp.Matrix([[p(j-i) for j in range(3)] for i in range(3)])
assert gram.det() == 1
coordinates = sp.Matrix([(n-1)*(n-2)/2, n*(2-n), n*(n-1)/2])
assert all(sp.expand(x) == 0 for x in gram*coordinates-sp.Matrix([p(n-i) for i in range(3)]))
assert sp.expand(sum(coordinates)) == 1
assert sp.expand((coordinates.T*gram*coordinates)[0]) == 1
assert sp.expand((sp.Matrix([[-5,-3,-1]])*coordinates)[0]) == 2*n-5
print('Euler Gram matrix:', gram)
print('Gram determinant:', gram.det())
print('Coordinates of [T^n A] in [A],[TA],[T^2 A]:', coordinates.T)

x, y, z0 = sp.symbols('x y z0')
vector = sp.Matrix([x,y,z0])
self_euler = sp.expand((vector.T*gram*vector)[0].subs(x,1-y-z0))
assert sp.expand(self_euler-(1+y+4*z0-(y+2*z0)**2)) == 0
print('Self Euler for x+y+z=1:', self_euler)
print('All exact checks passed.')
