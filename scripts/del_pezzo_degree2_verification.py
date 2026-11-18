"""Exact auxiliary calculations for a degree-2 Clifford-order example.

Coefficient field QQ.  The mathematical argument interpreting the discriminant
cover, duality, AZ ampleness, and the numerical classification is separate.
"""
import sympy as s

a,b,c,d=s.symbols('a b c d')
q=s.Matrix([[b,a+c,a+d],[a+c,b,b],[a+d,b,c]])
delta=s.expand(q.det())
print('Symmetric matrix q:')
print(q)
print('Determinant:',delta)
coeff=s.Matrix([[s.expand(q[i,j]).coeff(z) for z in (a,b,c,d)]
                for i,j in ((0,0),(0,1),(0,2),(1,1),(1,2),(2,2))])
assert coeff.rank()==4
print('Span of the six entries in H0(O(1,1)): exact QQ-rank',coeff.rank())

x,y=s.symbols('x y')
for first in (0,1):
    for second in (0,1):
        sx=(1,x) if first==0 else (x,1)
        ty=(1,y) if second==0 else (y,1)
        substitutions={a:sx[0]*ty[0],b:sx[0]*ty[1],
                       c:sx[1]*ty[0],d:sx[1]*ty[1]}
        f=s.expand(delta.subs(substitutions, simultaneous=True))
        groebner=s.groebner([f,s.diff(f,x),s.diff(f,y)],x,y,domain=s.QQ)
        assert list(groebner)==[1]
        print('Chart s%d=t%d=1:'%(first,second),f)
        print('  singularity-ideal Groebner basis over QQ:',list(groebner))

m=s.symbols('m',integer=True,nonnegative=True)
even=s.expand((m+1)**2+3*m**2)
odd=s.expand(3*(m+1)**2+m**2)
assert s.expand(even-1-(2*m)*(2*m+1))==0
assert s.expand(odd-1-(2*m+1)*(2*m+2))==0
print('Even Hilbert dimension:',even)
print('Odd Hilbert dimension:',odd)
t=s.symbols('t')
standard=(1-t**4)/((1-t)**3*(1-t**2))
claimed=1/(1-t)+2*t/(1-t)**3
assert s.cancel(standard-claimed)==0
print('Full degree-2 Hilbert series identity:',s.factor(standard))

G=s.Matrix([[1,3,2,2],[1,1,0,0],[0,0,1,-3],[0,0,-3,1]])
assert G.det()==16
print('Euler Gram matrix for (A,P,A(1,0),A(0,1)):')
print(G)
print('Gram determinant:',G.det())
alpha,beta,z=s.symbols('alpha beta z')
vec=s.Matrix([1-(alpha+beta)/2,z,(alpha-z)/2,(beta-z)/2])
assert s.simplify(sum(vec)-1)==0
selfeuler=s.factor((vec.T*G*vec)[0])
assert selfeuler==1-2*alpha*beta+2*z
chi=s.simplify((G*vec)[0].subs(z,alpha*beta))
assert s.simplify(chi-(1+alpha*beta+(alpha+beta)/2))==0
print('Rank-one self Euler with normalized determinant (alpha,beta):',selfeuler)
print('After exceptionality z=alpha*beta, chi(A,F):',chi)
print('All assertions passed.')
