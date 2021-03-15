function Adjg=piecewise_Adjoint(x,theta,xi)

adjxi       =matrix_adj(xi);

if theta==0
    Adjg        =diag([1 1 1 1 1 1])+x*adjxi;
else
    Adjg        =diag([1 1 1 1 1 1])+((3*sin(x*theta)-x*theta*cos(x*theta))/(2*theta))*adjxi+...
                 ((4-4*cos(x*theta)-x*theta*sin(x*theta))/(2*theta^2))*adjxi^2+...
                 ((sin(x*theta)-x*theta*cos(x*theta))/(2*theta^3))*adjxi^3+...
                 ((2-2*cos(x*theta)-x*theta*sin(x*theta))/(2*theta^4))*adjxi^4;
end

% eof