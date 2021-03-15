function coadj=matrix_coadj(screw)

coadj         =zeros(6,6);
coadj(1:3,1:3)=vector_tilde(screw(1:3));
coadj(1:3,4:6)=vector_tilde(screw(4:6));
coadj(4:6,4:6)=vector_tilde(screw(1:3));