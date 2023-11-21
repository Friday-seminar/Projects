C(U,V) = U*V-V*U;
nr(k) = 2*issquare(k) - (k==0);
X;Y;Z;
u=varhigher("u"); 
P = (u^3+a*u^2+b*u+c) - (u-X)*(u-Y)*(u-Z);
F = poldisc(P,u);
FT = substvec(F,[a,b,c],[-1-T,T,0]);
FT == (X*Y+Y*Z+X*Z-T)^2 + 4 *X*Y*Z * (1+T-(X+Y+Z))


H(x,y,z,t,p)  =  
{
x *= Mod(1,p); y *= Mod(1,p); z *= Mod(1,p); t*=Mod(1,p);
2 - nr(substvec(FT,[X,Y,Z,T],[x,y,z,t]))  \
		-  (x==y) * (1 + p*((x==0)+(x==1)+(x==t))) \
		+  p*(x!=0)*(x!=1)*(x!=t) * ( (z==0)*(x*y==t) + (z==1)*((1-x)*y==t-x) + (z==t)*(y*(t-x)==t*(1-x)) ) 
}

M(t=-1,p=3) = vector(p,x,matrix(p,p,y,z,H(x,y,z,t,p)));

W = M(2,11);

matrix(11,11,n,m,C(W[n],W[m]) %11 != 0) == 0

