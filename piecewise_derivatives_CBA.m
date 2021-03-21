function dz = piecewise_derivatives_CBA(t,z)
global gv
% Code for getting parameters.



num_piece = 2;
num_disc = 20; % This should be big enough.

M = zeros(6*num_piece);

for m = 1:num_piece
    for n = 1:num_piece
        
        M_mn = zeros(6);
        
        for i = max(m,n):N
            
            X = Li;
            
            for disc_index = 1:num_disc
                Sm = get_S(m); % Compute Sm.
                Sn = get_S(n);
                
                X = X + dx; % The actual position of disc.

                M_mn = M_mn + Sm' * Ma * Sn * dX;
            end
        end
        M(6*m-5:6*m, 6*n-5:6*n) = M_mn; % Fill in the integrated block.
    end
end

z1       =Xidot;
z2       =GIM^-1*(Tau+GM*matrix_Adjoint(g_r^-1)*Gra+ECL-(GCM1+GCM2)*Xidot);
dz       =[z1;z2];

end
