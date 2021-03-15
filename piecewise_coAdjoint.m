function coAdjg=piecewise_coAdjoint(x,theta,xi)

coadjxi       =matrix_coadj(xi);

if theta==0
    coAdjg        =diag([1 1 1 1 1 1])+x*coadjxi;
else
    coAdjg        =diag([1 1 1 1 1 1])+((3*sin(x*theta)-x*theta*cos(x*theta))/(2*theta))*coadjxi+...
                 ((4-4*cos(x*theta)-x*theta*sin(x*theta))/(2*theta^2))*coadjxi^2+...
                 ((sin(x*theta)-x*theta*cos(x*theta))/(2*theta^3))*coadjxi^3+...
                 ((2-2*cos(x*theta)-x*theta*sin(x*theta))/(2*theta^4))*coadjxi^4;
end