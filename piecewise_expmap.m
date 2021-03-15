function g=piecewise_expmap(x,theta,xi) %指数映射

xihat          =vector_hat(xi);
if theta==0
    g           =diag([1 1 1 1])+x*xihat;
else
    g           =diag([1 1 1 1])+x*xihat+...
                 ((1-cos(x*theta))/(theta^2))*xihat^2+...
                 ((x*theta-sin(x*theta))/(theta^3))*xihat^3;
end