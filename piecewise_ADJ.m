function ADJg=piecewise_ADJ(x,theta,xi)  %%指数映射的切算子

adjxi       =matrix_adj(xi);

if theta==0
    ADJg        =x*diag([1 1 1 1 1 1])+((x^2)/2)*adjxi;
else
    ADJg        =x*diag([1 1 1 1 1 1])+((4-4*cos(x*theta)-x*theta*sin(x*theta))/(2*theta^2))*adjxi+...
                 ((4*x*theta-5*sin(x*theta)+x*theta*cos(x*theta))/(2*theta^3))*adjxi^2+...
                 ((2-2*cos(x*theta)-x*theta*sin(x*theta))/(2*theta^4))*adjxi^3+...
                 ((2*x*theta-3*sin(x*theta)+x*theta*cos(x*theta))/(2*theta^5))*adjxi^4;
end