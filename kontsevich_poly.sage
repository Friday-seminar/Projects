import itertools
prime = 7
a, b, c = 6, 3, -1 # Kontsevich parameters


def check_affine_solutions(P, field):
    variabls = P.variables()
    all_affine_points = itertools.product(field, repeat=len(variabls))
    all_affine_values = map(lambda point: {k: v for k, v in zip(variabls, point)}, all_affine_points)
    solutions = filter(lambda x: P.subs(x) == 0 , all_affine_values)
    return len(list(solutions))

def Hecke_matrix(Px, field):
    q = field.cardinality() # field = F_q
    R.<x,y,z,t,w> = PolynomialRing(field, 5) # polynomial ring over F_q with the five variables x, y, z, t, w
    res = Matrix(ZZ, q, q)
    for yi, zi in itertools.product(range(q), repeat=2):
        res[yi, zi] = check_affine_solutions(Px.subs({y: yi, z: zi}), field)
    return res

def check_commutators(matrices, q):
    for j in range(q):
        for i in range(j):
            print('---'*q)
            print('[H_{}, H_{}] = '.format(i, j))
            print(matrices[i] * matrices[j] - matrices[j] * matrices[i], '\n')


def main(prime, a, b, c):
    field, q = GF(prime), prime 
    R.<x,y,z,t,w> = PolynomialRing(field, 5) # polynomial ring over F_q with the five variables x, y, z, t, w
    sigma1, sigma2, sigma3 = x + y + z, y*z + z*x + x*y, x*y*z # symmetric polys
    f = (a + sigma1) * t*t - (sigma2 - b) * t + (sigma3 + c) 
    P = f.discriminant(t) # Kontsevich poly

    my_equation = w^2 - P

    Hecke_operators = [ - Hecke_matrix( my_equation.subs({x:i}), field ) + matrix(ZZ, q, q, lambda i1, i2: 2) for i in range(q)]

    for i, op in enumerate(Hecke_operators):
        print('--'*q)
        print('H_{}= '.format(i))
        print(op)
    check_commutators(Hecke_operators, q)


print('Prime =', prime)
print('Kontsevich parameters: a = {}, b = {}, c = {}'.format(a,b,c))
print('Work in progress...')
main(prime, a, b, c)
print('Done!')
