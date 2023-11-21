C(U,V) = U*V-V*U;	/* commutator of 2 matrices */
nr(k) = 2*issquare(k) - (k==0); /* number of solutions of w^2 = k */
X;Y;Z; u=varhigher("u"); /* set priority of variables, just for user */
P = (u^3+a*u^2+b*u+c) - (u-X)*(u-Y)*(u-Z);  /* Roubtsov's quadratic polynomial */
F = poldisc(P,u);			   /* Roubtsov discriminant = gen Kon */
FT = substvec(F,[a,b,c],[-1-T,T,0]); /* MK's poly in VR's parameters */
FT == (X*Y+Y*Z+X*Z-T)^2 + 4 *X*Y*Z * (1+T-(X+Y+Z)) /* check it is indeed MK poly */

/* Next function is Hecke operators from MK's paper 2007
   p is a prime,
   x,y,z,t - points on affine line (shall be projective?)
*/
H(x,y,z,t,p)  =  
{
/* next line reduces x,y,z,t modulo p */
x *= Mod(1,p); y *= Mod(1,p); z *= Mod(1,p); t*=Mod(1,p);
2 - nr(substvec(FT,[X,Y,Z,T],[x,y,z,t]))  \
		-  (x==y) * (1 + p*((x==0)+(x==1)+(x==t))) \
		+  p*(x!=0)*(x!=1)*(x!=t) * ( (z==0)*(x*y==t) + (z==1)*((1-x)*y==t-x) + (z==t)*(y*(t-x)==t*(1-x)) ) 
}
/* for x=1,...,p matrix M(t,p)[x] is Hecke operator from MK */
M(t=-1,p=3) = vector(p,x,matrix(p,p,y,z,H(x,y,z,t,p)));

W = M(2,11);
/* check that all matrices  W[n] and W[m] commpute modulo p */
matrix(11,11,n,m,C(W[n],W[m]) %11 != 0) == 0

