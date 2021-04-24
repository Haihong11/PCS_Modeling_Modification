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
dX          =gv.dX;      %  dX   =   L/(num_disc-1);
X           =gv.X;       %   X   =   linspace(0,L,num_piece+1);        % [m]
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

Xi          =z(1:6*num_piece);           % strain twist 

Xidot       =z(6*num_piece+1:12*num_piece); % derivative of strain
g_r         =[0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];    % cantilever
g_prec      =diag([1 1 1 1]);                 %initial configuration
eta_prec    =zeros(6,1);                      % initial speed

g_prec_list = [g_prec];

eta_prec_list = [eta_prec];

% Record current joint states.

xi_array = [];
theta_array = [];
xidot_array = [];

% Pre-Compute xi and theta for each piece.
for iii = 1:N
      
    xin         = Xi(6*(iii-1)+1:6*(iii-1)+6,:);     % xi and theta in centain section, produced and saved           
    xidotn      = Xidot(6*(iii-1)+1:6*(iii-1)+6,:);           
    kn          = xin(1:3); 
    
    thetan      = sqrt(kn'*kn);
    xi_array    = [xi_array,xin];
    xidot_array = [xidot_array,xidotn];
    theta_array = [theta_array,thetan];

    X1 = X(iii);
    
    for ii=1:num_disc
       
        X1         = X1 + dX;
       
        invAdjg    = piecewise_invAdjoint(X1,thetan,xin);
        
        intdAdjg   = piecewise_ADJ(X1,thetan,xin);
        
        g(4*(nstep-1)+1:4*(nstep-1)+4,4*(iii-1)*num_disc+4*(ii-1)+1:4*(iii-1)*num_disc+4*(ii-1)+4)... %%
                   = g_r*g_prec*piecewise_expmap(X1,thetan,xin);
        
        eta(6*(nstep-1)+1:6*(nstep-1)+6,(iii-1)*num_disc+ii)...
                   = invAdjg*(eta_prec+intdAdjg*xidotn);
        
    end
    
    invAdjgn_last   =piecewise_invAdjoint(X(iii+1),thetan,xin);
    intdAdjgn_last  =piecewise_ADJ(X(iii+1),thetan,xin);
    g_prec          =g_prec*piecewise_expmap(X(iii+1),thetan,xin);
    ADxin           =intdAdjgn_last*xidotn;
    eta_prec        =invAdjgn_last*(eta_prec+ADxin);  
    g_prec_list     = [g_prec_list, g_prec];          % 4x4(N+1)
    eta_prec_list   = [eta_prec_list, eta_prec];      % 6x(N+1)
end

% Compute Jacobian before actual using it.
J_i_list = [];
J_pre = zeros(6,6*num_piece);
for i = 1:N
            
    L_i_1 = X(i);                  %  index of X     X=linspace(0,L,num_piece);
    X1    = L_i_1;                 %  L_i_1=L_{i-1};

    invAdjg_array  = [];
    intdAdjg_array = [];
    Adjg_array     = [];

    for ii = 1:num_disc
        X1 = X1 + dX;              %  The actual position of disc.

        %  Express S and dotS in the equation

        invAdjg        = piecewise_invAdjoint(X1,theta_array(:,iii),xi_array(:,iii));   %%  Ad^{-1}_g(X).

        intdAdjg       = piecewise_ADJ(X1,theta_array(:,iii),xi_array(:,iii));          %%  T_g(X).

%               Adjg           = piecewise_Adjoint(X1,theta_array(:,iii),xi_array(:,iii));

        invAdjg_array  = [invAdjg_array,invAdjg];

        intdAdjg_array = [intdAdjg_array,intdAdjg];

%               Adjg_array     = [Adjg_array,Adjg];

        % generate Jacobian matrix
        J_i = J_pre;            %  Equation 20

        J_i(:, 6*i-5:6*i) = intdAdjg;   % T_g(X).

        J_i = invAdjg * J_i;
        
        J_i_list = [J_i_list; J_i]; % 6*piece*disc x 6*piece.
        
        if ii == num_disc
            J_pre = J_i;
        end
    end
end


% CBA.   
for m = 1:num_piece 

    for n = 1:num_piece 
    
        M_mn   = zeros(6);
        C1_mn  = zeros(6);
        C2_mn  = zeros(6);
%       Tau_m  = zeros(6,1);
        G_m    = zeros(6,6);
%       ECL_m  = zeros(6,1);
        
        for i = max(m,n):N
            
            L_i_1 = X(i);                  %  index of X     X=linspace(0,L,num_piece);
            X1    = L_i_1;                 %  L_i_1=L_{i-1};
            
            M_mn_i = zeros(6);
            C1_mn_i = zeros(6);
            C2_mn_i = zeros(6);
            G_m_i = zeros(6);
            
            for ii = 1:num_disc
                
                X1 = X1 + dX;              %  The actual position of disc.
                
                intdAdjg       = piecewise_ADJ(X1,theta_array(:,iii),xi_array(:,iii));          %%  T_g(X).

                disc_index = num_disc*(i-1)+ii;
                J_i = J_i_list(6*disc_index-5:6*disc_index,:);
                
                Sm = J_i(:,6*m-5:6*m);                                      %  Compute S. where S is the function of the abscissi X

                Sn = J_i(:,6*n-5:6*n);
                %                dotSn = dotJ_i(n);

                M_mn_i  = M_mn_i + Sm' * M * Sn * dX;                                                                % equation 31 

                C1_mn_i = C1_mn_i + Sm'* matrix_coadj(eta_prec_list(:,i)+intdAdjg.*xidot_array(:,i))*M* Sn * dX;  % equation 32

                C2_mn = C2_mn + Sm'* M * matrix_adj(intdAdjg.*xidot_array(:,i))*Sn * dX;        

    %                C2_mn_i = C2_mn_i + Sm'* M *dotSn * dX;                                                              % equation 33

                G_m_i   = G_m_i + Sm' * M * matrix_Adjoint(g_prec_list(:, 4*i-3:4*i)^-1)*dX;                                            % equation 35
                
            end
            
            M_mn = M_mn + M_mn_i;
            C1_mn = C1_mn + C1_mn_i;
            C2_mn = C2_mn + C2_mn_i;
            G_m = G_m + G_m_i;
        end
        
    	%   Actuation and internal load

        if t<=tact         %    turn                        
          Fan         =[Famx(i);Famy(i);Famz(i);Fax(i);Fay(i);Faz(i)]*t/tact;
        end

        if (t>tact && t<=trel+tact)                                   % release / maintenance  
          Fan   =[Famx(i);Famy(i);Famz(i);Fax(i);Fay(i);Faz(i)];
        end

        Fin     = Eps*(xin-xi_star)+Ipsi*xidotn;

        Tau_m   = L*(Fan-Fin);                                      %    Equation 29

        % Equation 29 任意一个section \tau_n等于从该段开始驱动力的总和与该段的内力，再乘该段长度

        Fpn     = [Fpmx(i);Fpmy(i);Fpmz(i);Fpx(i);Fpy(i);Fpz(i)];   % External Concentrated load

        ECL_m     = Sm' *  Fpn;                                     % Equation 30

        GIM(6*m-5:6*m, 6*n-5:6*n)  = M_mn;   %  Fill in the integrated block.
        GCM1(6*m-5:6*m, 6*n-5:6*n) = C1_mn;
        GCM2(6*m-5:6*m, 6*n-5:6*n) = C2_mn;
        Tau(6*m-5:6*m,1)           = Tau_m;  %  only related to the row    
        GM(6*m-5:6*m,1:6)          = G_m;    %  equation 26       
        ECL(6*m-5:6*m,1)           = ECL_m;        
    end
end

dotz1       = Xidot;
dotz2       = GIM^-1*(Tau+GM*matrix_Adjoint(g_r^-1)*Gra+ECL-(GCM1+GCM2)*Xidot);
dz          = [dotz1;dotz2];

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