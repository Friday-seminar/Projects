''' Niels thesis '''
''' Playing with Hecke algebras '''

# The goal is to compute alpha^x_{y, z}


import itertools
import numpy as np
prime = 5
t = 2 # Define the divisor D = {0, 1, t, inf}


def commutator(A, B):
    return A * B - B * A

def M0(z, t): # z != 0
    return t/z

def M1(z, t): # z != 1
    return (z - t) / (z - 1)

def Mt(z, t): # z != t
    return t * (z - 1) / (z - t)

def check_commutators(matrices, q):
    ''' Check that the Hecke matrices commute '''
    for j in range(q):
        for i in range(j):
            print('---'*q)
            print('[H_{}, H_{}] = '.format(i, j))
            print(commutator(matrices[i], matrices[j]))


def main(prime, t):
    field, field_elements, q = GF(prime), [i for i in GF(prime)], prime 
    R.<x,y,z,r> = PolynomialRing(field, 4) # F_q[x, y, z, r]
    t = field(t) # now t belongs to F_q
    D = {field(0) , field(1), field(t)} # the divisor w/o infty :sad

    equation_x_is_not_y = (y * r - x) * ( (y - 1) * (y - t) * r - (x - 1) * (x - t)) + z * (x - y) * (x - y) * r
    equation_x_is_y = z + (y * r - 1) * ( (y - 1) * (y - t) * r - (2 * y -(1 + t))) 

    alpha = np.zeros((q, q, q)) # consists of alpha^x_{z,y}
    for ix, iy, iz in itertools.product(field_elements,repeat=3):
        xx = field_elements[ix]
        yy = field_elements[iy]
        zz = field_elements[iz]
        if ix != iy: # x != y
            specialized_equation = equation_x_is_not_y.subs({x:xx, y:yy, z:zz})
            if specialized_equation in field:
                if specialized_equation == 0:
                    num_nonzero_roots = q - 1
                else:
                    num_nonzero_roots = 0
            else:
                specialized_equation = specialized_equation.polynomial(r)
                num_nonzero_roots = len(specialized_equation.roots())
                if specialized_equation.subs({r:0}) == 0:
                    num_nonzero_roots -= 1

            if xx in D and yy in D:
                if xx == field(0) and z == M0(y, t):
                    corrections = q
                elif xx == field(1) and z == M1(y, t):
                    corrections = q
                elif xx == field(t) and z == Mt(y, t):
                    corrections = q
                else:
                    corrections = 1 # non clear from the paper

            elif xx in D or yy in D:
                corrections = 1
            else:
                corrections = 2

            if zz in D:
                if zz == field(0) and yy == M0(x, t):
                    corrections += q
                elif zz == field(1) and yy == M1(x, t):
                    corrections += q
                elif zz == field(t) and yy == Mt(x, t):
                    corrections += q

            alpha[ix, iy, iz] = num_nonzero_roots - corrections
        else: # x = y
            specialized_equation = equation_x_is_y.subs({x:xx, y:yy, z:zz})
            if specialized_equation in field:
                if specialized_equation == 0:
                    num_roots = q
                else:
                    num_roots = 0
            else:
                specialized_equation = specialized_equation.polynomial(r)
                num_roots = len(specialized_equation.roots())

            alpha[ix, iy, iz] = num_roots - q + 1

    matrices = [Matrix(ZZ, ai) for ai in alpha]
    check_commutators(matrices, q)

print('Prime = {}, t = {}'.format(prime, t))
print('Work in progress...')
main(prime, t)
print('Done!')
