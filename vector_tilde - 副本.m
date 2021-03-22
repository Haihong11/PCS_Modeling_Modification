function skew=vector_tilde(vec)

skew        =zeros(3,3);
skew(1,2)   =-vec(3);
skew(1,3)   =vec(2);
skew(2,1)   =vec(3);
skew(2,3)   =-vec(1);
skew(3,1)   =-vec(2);
skew(3,2)   =vec(1);     %%skew-symmetric matrix   Lie algebra so(3)of SO(3)