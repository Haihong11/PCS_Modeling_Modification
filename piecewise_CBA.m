function dz = piecewise_CBA(t,z)
global gv

t

% global variable
L           =gv.L;
Eps         =gv.Eps;
Ipsi        =gv.Ipsi;
M           =gv.M;       %  M    =   ro_arm*diag([I J J A A A]); % mass matrix of the disc
xi_star     =gv.xi_star;
Gra         =gv.Gra;
dX          =gv.dX;      %  dX   = L/(num_disc-1);
X           =gv.X;       %  X   = linspace(0,L,num_piece+1);        % [m]
tact        =gv.tact;    %  loading time
trel        =gv.trel;    %  relaxing time
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
g           =gv.g;
eta         =gv.eta;
nstep       =gv.nstep;

% Code for getting parameters.

num_piece        = 2;  
num_disc         = 20;          % This should be big enough.
N                = num_piece;

GIM         = zeros(6*num_piece);        % mass matrix
GCM1        = zeros(6*num_piece);        % Coriolis matrix 1
GCM2        = zeros(6*num_piece);        % Coriolis matrix 2
Tau         = zeros(6*num_piece,1);      % Internal force and actuation force L(F_a-F_i)
GM          = zeros(6*num_piece,6);      % Gravitational matrix
ECL         = zeros(6*num_piece,1);      % External Concentrated Load



Xi          =z(1:6*num_piece);
Xidot       =z(6*num_piece+1:12*num_piece);
g_r         =[0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];    % cantilever
g_prec      =diag([1 1 1 1]); %初始位形
eta_prec    =zeros(6,1);
g_prec_list = [g_prec];
eta_prec_list = [eta_prec];

% Record current joint states.
xi_array = [];
theta_array = [];
xidot_array = [];
for iii = 1:N
    
    % xi and theta in every piece

    xin         = Xi(6*(iii-1)+1:6*(iii-1)+6,:);     % xi and theta in centain section, produced and saved           
    xidotn      = Xidot(6*(iii-1)+1:6*(iii-1)+6,:);           
    kn          = xin(1:3);
    thetan      = sqrt(kn'*kn);

    xi_array    = [xi_array,xin];
    xidot_array = [xidot_array,xidotn];
    theta_array = [theta_array,thetan];
    
    X1 = X(iii);
    for ii=1:num_disc
        X1 = X1 + dX;
        invAdjg    =piecewise_invAdjoint(X1,thetan,xin);

        intdAdjg   =piecewise_ADJ(X1,thetan,xin);

    %-----所有disc的位置组成g矩阵-----------%一个disc对应是4×4的矩阵，这样一个g就会有4×disc×piece列（位置）
    % g有4×n行，因为每一个时刻就会有一条构型曲线。（时间）
    %每一个disc在不同时刻的位置和速度
        g(4*(nstep-1)+1:4*(nstep-1)+4,4*(iii-1)*num_disc+4*(ii-1)+1:4*(iii-1)*num_disc+4*(ii-1)+4)... %%
            =g_r*g_prec*piecewise_expmap(X1,thetan,xin);
        eta(6*(nstep-1)+1:6*(nstep-1)+6,(iii-1)*num_disc+ii)...
            =invAdjg*(eta_prec+intdAdjg*xidotn);
    end
    invAdjgn_last   =piecewise_invAdjoint(X(iii+1),thetan,xin);
    intdAdjgn_last  =piecewise_ADJ(X(iii+1),thetan,xin);
    g_prec          =g_prec*piecewise_expmap(X(iii+1),thetan,xin);%%---position
    ADxin           =intdAdjgn_last*xidotn;
    eta_prec        =invAdjgn_last*(eta_prec+ADxin);  %%--speed
    
    g_prec_list = [g_prec_list, g_prec]; % 4x4(N+1).
    eta_prec_list = [eta_prec_list, eta_prec]; % 6x(N+1).
end

% CBA.
for m = 1:num_piece

    for n = 1:num_piece
    
        M_mn   = zeros(6);
        C1_mn  = zeros(6);
        C2_mn  = zeros(6);
        Tau_m  = zeros(6,1);
        G_m    = zeros(6,6);
        ECL_m  = zeros(6,1);
        
        for i = max(m,n):N
            
            L_i_1 = X(i);                  %  index of X     X=linspace(0,L,num_piece);
            
            X1    = L_i_1;                 %  L_i_1=L_{i-1};
            M_mn_i = zeros(6);
            C1_mn_i = zeros(6);
            C2_mn_i = zeros(6);
            G_m_i = zeros(6);
            for ii = 1:num_disc
            
                X1 = X1 + dX;              %  The actual position of disc.
                
                invAdjg_array  = [];
                intdAdjg_array = [];
                Adjg_array     = [];
                
                for jj = 1:i               %  Express S and dotS in the equation
                    
                    invAdjg        = piecewise_invAdjoint(X1,theta_array(:,jj),xi_array(:,jj));   %%  Ad^{-1}_g 
                    
                    intdAdjg       = piecewise_ADJ(X1,theta_array(:,jj),xi_array(:,jj));          %%  T_g
                    
                    Adjg           = piecewise_Adjoint(X1,theta_array(:,jj),xi_array(:,jj));      %%  Ad_g
                       
                    invAdjg_array  = [invAdjg_array,invAdjg];
                    
                    intdAdjg_array = [intdAdjg_array,intdAdjg];
                    
                    Adjg_array     = [Adjg_array,Adjg];
                    
                end
%                                     
%                 g(4*(nstep-1)+1:4*(nstep-1)+4,4*(jj-1)*num_disc+4*(ii-1)+1:4*(jj-1)*num_disc+4*(ii-1)+4)...
%                         =g_r*g_prec*piecewise_expmap(X1,theta_array(:,jj),xi_array(:,jj));
%                 
%                 eta(6*(nstep-1)+1:6*(nstep-1)+6,(jj-1)*num_disc+ii)...
%                         =invAdjg*(eta_prec+intdAdjg*xidot_array(:,jj));
%                
%                 
% %                   eta_prec       = invAdjg*(eta_prec+intdAdjg*xidot_array(:,jj));

                % generate Jacobian matrix
          
                J_i = [];                              %  Equation 20
                
                for j = i:-1:1                         %  non-zero items in certain piece           

                    item = intdAdjg_array(:,i-j+1);   

                    for k = 1:j                       

                        item = invAdjg_array(k)*item;  % prod

                    end

                    J_i   = [J_i,item];                % From left to right
                end

                for j=1:N-i                            %  zero items in certain piece      

                    item = zeros(6);                   

                    J_i = [J_i,item];                                  

                end
                    
                % generate the derivative of Jacobian matrix
          
                dotJ_i = [];                                       %  Equation 20
                
                for j = i:-1:1                                 %  non-zero items in certain piece

                    item = Adjg_array(:,i-j+1).* matrix_adj(eta_prec);  %  integral item

                    for k = 1:j                      

                        item = invAdjg_array(k)*item;          % prod

                    end

                    dotJ_i   = [dotJ_i,item]; 

                end 

                for j=1:N-i                                    %  zero items in certain piece      

                    item     = zeros(6);

                    dotJ_i   = [dotJ_i,item];

                end 
                    
                 Sm = J_i(m);                                      %  Compute S. where S is the function of the abscissi X
                 
                 Sn = J_i(n);
                 
                 dotSn = dotJ_i(n);
                         
                 M_mn_i  = M_mn_i + Sm' * M * Sn * dX;                                                              % equation 31 
                 
                 C1_mn_i = C1_mn_i + Sm'* matrix_coadj(eta_prec_list(:,i)+intdAdjg_array(:,jj)).*xidot_array(:,jj)*M* Sn * dX;  % equation 32
                                         
%                C2_mn = C2_mn + Sm'* M * matrix_adj(intdAdjg*xidot_array(:,jj))*Sn * dX;        

                 C2_mn_i = C2_mn_i + Sm'* M * dotSn * dX;                                                              % equation 33
                 
                 G_m_i   = G_m_i + Sm' * M * matrix_Adjoint(g_prec_list(:, 4*i-3:4*i)^-1)*dX;                                            % equation 35
                 
            end
            
            M_mn = M_mn + M_mn_i;
            C1_mn = C1_mn + C1_mn_i;
            C2_mn = C2_mn + C2_mn_i;
            G_m = G_m + G_m_i;
            
            %    recursive factors
                 
%              invAdjg_last   =piecewise_invAdjoint(X1,theta_array(:,jj),xi_array(:,jj));
%              intdAdjg_last  =piecewise_ADJ(X1,theta_array(:,jj),xi_array(:,jj));
%              g_prec          =g_prec*piecewise_expmap(X1,theta_array(:,jj),xi_array(:,jj));%%---position
%              ADxin           =intdAdjg_last*xidot_array(:,jj);
%              eta_prec        =invAdjg_last*(eta_prec+ADxin);                               %%--speed
                 
            %   Actuation and internal load
                
                if t<=tact         %    turn                        
                  Fan         =[Famx(i);Famy(i);Famz(i);Fax(i);Fay(i);Faz(i)]*t/tact;
                end
                
                if (t>tact && t<=trel+tact)                            % release / maintenance  
                  Fan   =[Famx(i);Famy(i);Famz(i);Fax(i);Fay(i);Faz(i)];
                end
                
                  Fin     = Eps*(xin-xi_star)+Ipsi*xidotn;

                  Tau_m   = Tau_m + L*(Fan-Fin);                          %    Equation 29

                % Equation 29 任意一个section \tau_n等于从该段开始驱动力的总和与该段的内力，再乘该段长度
                % External Concentrated load
                  Fpn     = [Fpmx(i);Fpmy(i);Fpmz(i);Fpx(i);Fpy(i);Fpz(i)];
                  
                  ECL_m     = ECL_m + Sm' *  Fpn;                          %    Equation 30
                  
        end       
        
        GIM(6*m-5:6*m, 6*n-5:6*n)  = M_mn;   %  Fill in the integrated block.
        GCM1(6*m-5:6*m, 6*n-5:6*n) = C1_mn;
        GCM2(6*m-5:6*m, 6*n-5:6*n) = C2_mn;
        Tau(6*m-5:6*m,1)           = Tau_m;  %  only related to the row    
        GM(6*m-5:6*m,1:6)          = G_m;    %  equation 26       
        ECL(6*m-5:6*m,1)           = ECL_m;
           
    end
end

dotz1       =Xidot;
dotz2       =GIM^-1 * (Tau+GM*matrix_Adjoint(g_r^-1)*Gra+ECL-(GCM1+GCM2)*Xidot);
dz       =[dotz1;dotz2];

gv.g     = g;
gv.eta   = eta;
gv.nstep = nstep+1;
gv.GIM   = GIM;   % mass matrix
gv.GM    = GM;    % Gravitational matrix
gv.GCM1  = GCM1;  % Coriolis matrix 1
gv.GCM2  = GCM2;  % Coriolis matrix 2
gv.Tau   = Tau;   % Internal force and actuation force L(F_a-F_i)
gv.ECL   = ECL;   % External Concentrated Load


end