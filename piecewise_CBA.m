function dz = piecewise_CBA(t,z)
global gv


% global variable
L           =gv.L;
Eps         =gv.Eps;
Ipsi        =gv.Ipsi;
M           =gv.M;       %  M    =   ro_arm*diag([I J J A A A]); % mass matrix of the disc
xi_star     =gv.xi_star;
Gra         =gv.Gra;
dX          =gv.dX;      %  dX   =    L/(num_disc-1);
X           =gv.X;       %   X   =    linspace(0,L,num_disc);        % [m]
tact        =gv.tact;    %  加载时间
trel        =gv.trel;    %  松弛时间
Fax         =gv.Fax;
Fay         =gv.Fay;
Faz         =gv.Faz;
Famx        =gv.Famx;
Famy        =gv.Famy;
Famz        =gv.Famz;
Fpx         =gv.Fpx;
Fpy         =gv.Fpy;
Fpz         =gv.Fpz;
Fpmx        =gv.Fpmx;
Fpmy        =gv.Fpmy;
Fpmz        =gv.Fpmz; 
nstep       =gv.nstep; 

% Code for getting parameters.
g_r              =[0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];     % 惯性系与body frame之间变换矩阵
g_prec           =diag([1 1 1 1]); %Initial configuration
eta_prec         =zeros(6,1);      %initial speed 

num_piece        = 2;  
num_disc         = 20;          % This should be big enough.
N                = num_piece;

GIM         = zeros(6*num_piece);        % mass matrix
GCM1        = zeros(6*num_piece);        % Coriolis matrix 1
GCM2        = zeros(6*num_piece);        % Coriolis matrix 2
Tau         = zeros(6*num_piece,1);      % Internal force and actuation force L(F_a-F_i)
GM          = zeros(6*num_piece,6);      % Gravitational matrix
ECL         = zeros(6*num_piece,1);      % External Concentrated Load

Xi          =z(1:6*num_piece,:);  
Xidot       =z(6*num_piece+1:12*num_piece,:);

for m = 1:num_piece 
    
    for n = 1:num_piece 
        
        M_mn   = zeros(6);
        C1_mn  = zeros(6);
        C2_mn  = zeros(6);
%       Tau_m  = zeros(6,1);
        G_m    = zeros(6,6);
%       ECL_m  = zeros(6,1);
        xi_array = [];
        theta_array = [];
        xidot_array = [];
        
        for iii = 1:N           % xi and theta in every piece
            xin         = Xi(6*(iii-1)+1:6*(iii-1)+6,:);     % xi and theta in centain section, produced and saved           
            xidotn      = Xidot(6*(iii-1)+1:6*(iii-1)+6,:);           
            kn          = xin(1:3);        
            thetan      = sqrt(kn'*kn);
           
            xi_array    = [xi_array,xin];
            theta_array = [theta_array,thetan];
            xidot_array = [xidot_array,xidotn];  
        end
        
        for i = max(m,n):N
            
            L_i_1 = X(i);                %  index of X     X=linspace(0,L,num_piece);
            
            X1    = L_i_1;                 %  L_i_1=L_{i-1};
                           
            for ii = 1:num_disc
            
                X1 = X1 + dX;              %  The actual position of disc.
                
                invAdjg_array  = [];
                intdAdjg_array = [];
                
                for jj = 1:i

                    invAdjg        = piecewise_invAdjoint(X1,theta_array(:,jj),xi_array(:,jj));   %%  Ad^{-1}_g
                    
                    intdAdjg       = piecewise_ADJ(X1,theta_array(:,jj),xi_array(:,jj));          %%  T_g
                    
%                   g_prec         = g_prec*piecewise_expmap(X1,theta_array(:,jj),xi_array(:,jj));
                    
%                   eta_prec       = invAdjg*(eta_prec+intdAdjg*xidot_array(:,jj));

                    invAdjg_array  = [invAdjg_array,invAdjg];

                    intdAdjg_array = [intdAdjg_array,intdAdjg];

                end    
                
                % generate Jacobian matrix
          
                J_i = [];                                 %  equation 20
             
                    for j = i:-1:1                        %  non-zero items            
                        
                        item = intdAdjg_array(:,i-j+1);   %  prod

                        for k = 1:j                      
                            
                            item = invAdjg_array(k)*item;   
                        end
               
                        J_i   = [J_i,item];
                    end
            
                    for j=1:N-i                        %  zero items       

                        item = zeros(6);
                       
                        J_i = [J_i,item];
                        
                    end 
                    
                 Sm = J_i(m);                           %  Compute S. where S is the function of the abscissi X
                 
                 Sn = J_i(n);
                                
                 M_mn  = M_mn + Sm' * M* Sn * dX;                                                  % equation 31                  
                 C1_mn = C1_mn + Sm'* matrix_coadj(eta_prec+intdAdjg*xidot_array(:,jj))*M* Sn * dX;% equation 32
                 C2_mn = C2_mn + Sm'* M * matrix_adj(intdAdjg*xidot_array(:,jj))*Sn * dX;          % equation 33
                 G_m   = G_m + Sm' * M * matrix_Adjoint(g_prec^-1)*dX;                             % equation 35
                
            end
            
            %   Actuation and internal load
                
                if t<=tact                                      %    turn                        
                  Fan         =[Famx(i);Famy(i);Famz(i);Fax(i);Fay(i);Faz(i)]*t/tact;
                end
                
                if (t>tact && t<=trel+tact)                     % release / maintenance  
                  Fan   =[Famx(i);Famy(i);Famz(i);Fax(i);Fay(i);Faz(i)];
                end
                
                  Fin   =Eps*(xin-xi_star)+Ipsi*xidotn;

                  Tau_m = L*(Fan-Fin);                          %    equation 29

                %equation 29 任意一个section \tau_n等于从该段开始驱动力的总和与该段的内力，再乘该段长度
                % External Concentrated load
                Fpn     =[Fpmx(i);Fpmy(i);Fpmz(i);Fpx(i);Fpy(i);Fpz(i)];
                ECL_m   = Sm' * Fpn;      % equation 30

        end       
                       
        GIM(6*m-5:6*m, 6*n-5:6*n)  = M_mn;   %  Fill in the integrated block.
        GCM1(6*m-5:6*m, 6*n-5:6*n) = C1_mn;
        GCM2(6*m-5:6*m, 6*n-5:6*n) = C2_mn;
        Tau(6*m-5:6*m,1)           = Tau_m;      
        GM(6*m-5:6*m,1:6)          = G_m;    %   equation 26       
        ECL(6*m-5:6*m,1)           = ECL_m;
        
    end
end

z1       =Xidot;
z2       =GIM^-1*(Tau+GM*matrix_Adjoint(g_r^-1)*Gra+ECL-(GCM1+GCM2)*Xidot);
dz       =[z1;z2];

 gv.GIM  = GIM;   % mass matrix
 gv.GM   = GM;    % Gravitational matrix
 gv.GCM1 = GCM1;  % Coriolis matrix 1
 gv.GCM2 =GCM2;   % Coriolis matrix 2
 gv.Tau  =Tau;    % Internal force and actuation force L(F_a-F_i)
 gv.ECL =ECL;     % External Concentrated Load
 
end