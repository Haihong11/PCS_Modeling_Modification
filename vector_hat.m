function se3=vector_hat(screw)

se3         =zeros(4,4);
se3(1:3,1:3)=vector_tilde(screw(1:3)); 
se3(1:3,4)  =screw(4:6); 