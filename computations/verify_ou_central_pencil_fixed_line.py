"""The fixed line and the moving conic in a central pencil on the Fermat cubic."""

from sympy import Matrix, Rational, binomial, cancel, expand, factor, hessian, symbols


def verify():
    u, w, x, lam = symbols('u w x lam')
    plane_section = expand(u**3+(x-u)**3+w**3+(lam*x-w)**3)
    conic = cancel(plane_section/x)
    matrix = hessian(conic, (u, w, x))/2
    coords = Matrix([u, w, x])
    assert expand((coords.T*matrix*coords)[0, 0]-conic) == 0
    assert factor(matrix.det()-Rational(9, 4)*lam*(1+lam**3)) == 0
    assert conic.subs(x, 0) == 3*u**2+3*lam*w**2
    # In k[u,w,x]/(x*q), the x-power torsion is (q), isomorphic to
    # k[u,w](-2). The saturation is the coordinate ring of the conic q.
    n = symbols('n', integer=True, positive=True)
    def hilbert(degree):
        return expand((binomial(n+2, 2)-binomial(n+2-degree, 2)).expand(func=True))
    assert hilbert(3) == 3*n and hilbert(2) == 2*n+1
    assert expand(hilbert(3)-hilbert(2)) == n-1
    print('PASS fixed-line pencil on the Fermat cubic, determinant', factor(matrix.det()))
    print('PASS full, moving and fixed Hilbert polynomials', hilbert(3), hilbert(2), n-1)


if __name__ == '__main__':
    verify()
