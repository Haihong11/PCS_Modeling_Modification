function dz = piecewise_CBA(t,z)
global gv

t

%--------------global variable----------------%

L           = gv.L;
Eps         = gv.Eps;
Ipsi        = gv.Ipsi;
M           = gv.M;       %  M    =   ro_arm*diag([I J J A A A]); % mass matrix of the disc
xi_star     = gv.xi_star;
Gra         = gv.Gra;
dX          = gv.dX;      %  dX   =   L/num_piece/num_disc;
X           = gv.X;       % [m] X =   linspace(0,L,num_piece+1); end point
tact        = gv.tact;    %  loading time
trel        = gv.trel;    %  relaxing time
Fax         = gv.Fax;
Fay         = gv.Fay;
Faz         = gv.Faz;
Famx        = gv.Famx;
Famy        = gv.Famy;
Famz        = gv.Famz;
Fpx         = gv.Fpx;
Fpy         = gv.Fpy;
Fpz         = gv.Fpz;
Fpmx        = gv.Fpmx;
Fpmy        = gv.Fpmy;
Fpmz        = gv.Fpmz;
g           = gv.g;
eta         = gv.eta;
nstep       = gv.nstep; 
num_piece   = gv.num_piece;
num_disc    = gv.num_disc;

% Code for getting parameters. 
                  
N           = num_piece;
GIM         = zeros(6*num_piece);        % Generalized inertia matrix
GCM1        = zeros(6*num_piece);        % Generalized Coriolis matrix 1
GCM2        = zeros(6*num_piece);        % Generalized Coriolis matrix 2
Tau         = zeros(6*num_piece,1);      % Internal force and actuation force L(F_a-F_i)
GM          = zeros(6*num_piece,6);      % Gravitational matrix
ECL         = zeros(6*num_piece,1);      % External Concentrated Load

% Current system state

Xi          = z(1:6*num_piece);          % strain twist 
Xidot       = z(6*num_piece+1:12*num_piece);            % derivative of strain
g_r         = [0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];    % the map between the spatial frame and the base frame

% g (w.r.t initial frame) and eta (w.r.t body frame) of the beginning of each piece

g_pre       = diag([1 1 1 1]);                          % initial configuration  
eta_pre     = zeros(6,1);                               % initial velocity

g_pre_list   = [];                                       

eta_pre_list = []; 

% Record current joint states.

xi_array    = [];
theta_array = [];
xidot_array = [];

% Pre-Compute xi and theta for each piece. Pre-Compute g and eta for each disc %%

for i = 1:N
      
    xin         = Xi(6*(i-1)+1:6*(i-1)+6,:);     % xi and theta in centain section           
    xidotn      = Xidot(6*(i-1)+1:6*(i-1)+6,:);           
    kn          = xin(1:3); 
    
    thetan      = sqrt(kn'*kn);
    xi_array    = [xi_array,xin];                    % produced and saved
    xidot_array = [xidot_array,xidotn];
    theta_array = [theta_array,thetan];

    X1 = X(i);                                     % the beginning of abscissa of each piece
    
    for ii=1:num_disc
        
        X1         = X1 + dX;
        
        invAdjg    = piecewise_invAdjoint(X1,thetan,xin); % Ad^{-1}_g
        
        intdAdjg   = piecewise_ADJ(X1,thetan,xin);        % T_g
        
        g(4*(nstep-1)+1:4*(nstep-1)+4,4*(i-1)*num_disc+4*(ii-1)+1:4*(i-1)*num_disc+4*(ii-1)+4)... 
                   = g_r*g_pre*piecewise_expmap(X1,thetan,xin);  %     
        
        eta(6*(nstep-1)+1:6*(nstep-1)+6,(i-1)*num_disc+ii)... 
                   = invAdjg*(eta_pre+intdAdjg*xidotn);          %
        
    end
       
    invAdjg_next   = piecewise_invAdjoint(X(i+1),thetan,xin);
    intdAdjg_next  = piecewise_ADJ(X(i+1),thetan,xin);
    g_pre          = g_pre*piecewise_expmap(X(i+1),thetan,xin);
    eta_pre        = invAdjg_next*(eta_pre+intdAdjg_next*xidotn);  
    g_pre_list     = [g_pre_list, g_pre];                     % 4x4(N+1)?  4x4N
    eta_pre_list   = [eta_pre_list, eta_pre];                 %  6x(N+1)?  6xN
      
end

% Compute Jacobian before actual using it.

J_i_list = [];                     % Save Jacobian matrix of all discs

J_pre    = zeros(6,6*num_piece);   % initial Jacobian matrix

for i = 1:N                     
            
    L_i_1 = X(i);                  % index of X     X=linspace(0,L,num_piece+1);
    X1    = L_i_1;                 % L_i_1=L_{i-1};
    invAdjg_array  = [];
    intdAdjg_array = [];

    for ii = 1:num_disc
        
        X1             = X1 + dX;                                                   %  The actual position of disc.

        invAdjg        = piecewise_invAdjoint(X1,theta_array(:,i),xi_array(:,i));   %  Ad^{-1}_g(X).

        intdAdjg       = piecewise_ADJ(X1,theta_array(:,i),xi_array(:,i));          %  T_g(X).

        invAdjg_array  = [invAdjg_array,invAdjg];

        intdAdjg_array = [intdAdjg_array,intdAdjg];

        % generate Jacobian matrix
        
        J_i = J_pre;                    %  Equation 20
        
        J_i(:, 6*i-5:6*i) = intdAdjg;   %  T_g(X).
        
        J_i = invAdjg * J_i;
        
        J_i_list = [J_i_list; J_i];     %  6*piece*disc x 6*piece.
        
        if ii == num_disc                
          
           J_pre = J_i;                 
           
        end
    end
   
end

% Compute the derivative of Jacobian before actual using it.

% dotJ_i_list = [];                     % Save derivative of Jacobian matrix of all discs
% 
% dotJ_pre    = zeros(6,6*num_piece);   % initial Jacobian matrix
% 
% for i = 1:N                           
%             
%     L_i_1 = X(i);                  % index of X     X=linspace(0,L,num_piece+1);
%     X1    = L_i_1;                 % L_i_1=L_{i-1};
%     Adjg_array      = [];
%     invAdjg_array   = [];
%     adj_eta_array   = [];
%     
%     for ii = 1:num_disc
%     
%         X1             = X1 + dX;                                                   %  The actual position of disc.
% 
%         Adjg           = piecewise_Adjoint(X1,theta_array(:,i),xi_array(:,i));      %  Ad^_g(X).
%         
%         adj_eta        = matrix_adj(eta_pre_list(:,i)+intdAdjg.*xidot_array(:,i));
%         
%         invAdjg        = piecewise_invAdjoint(X1,theta_array(:,i),xi_array(:,i));   %  Ad^{-1}_g(X).  
%         
%         Adjg_array     = [Adjg_array,Adjg];
%         
%         adj_eta_array  = [adj_eta_array,adj_eta];
%         
%         invAdjg_array  = [invAdjg_array,invAdjg];
%         
%         % generate derivative of Jacobian matrix
%         
%         dotJ_i = dotJ_pre;                    %  Equation 20
%         
%         dotJ_i(:, 6*i-5:6*i) = intdAdjg;      %  integral
%         
%         dotJ_i = invAdjg * dotJ_i;
%         
%         dotJ_i_list = [dotJ_i_list; dotJ_i];     %  6*piece*disc x 6*piece.
%         
%         if ii == num_disc                
%           
%            dotJ_pre = dotJ_i;                 
%            
%         end
%     end   
% end

% Composite Body Algorithm (CBA).

for m = 1:num_piece 
    
    for n = 1:num_piece 
    
        M_mn   = zeros(6);
        C1_mn  = zeros(6);
        C2_mn  = zeros(6);
        G_m    = zeros(6,6);
        
        for i = max(m,n):N
            
            L_i_1   = X(i);               %  index of X     X=linspace(0,L,num_piece);
            X1      = L_i_1;              %  L_i_1=L_{i-1};
            M_mn_i  = zeros(6);
            C1_mn_i = zeros(6);
            C2_mn_i = zeros(6);
            G_m_i   = zeros(6);
            
            % The integral for the ith piece
            
            for ii = 1:num_disc
                
                X1         = X1 + dX;     %  The actual position of disc.
                
                intdAdjg   = piecewise_ADJ(X1,theta_array(:,i),xi_array(:,i));  %  T_g(X).
                
                 g_i       = piecewise_expmap(X1,theta_array(:,i),xi_array(:,i));
                 
                disc_index = num_disc*(i-1)+ii;
                
                J_i        = J_i_list(6*disc_index-5:6*disc_index,:);
                
                Sm         = J_i(:,6*m-5:6*m);            %  Compute S. where S is the function of the abscissi X
                
                Sn         = J_i(:,6*n-5:6*n);
                
%               dotSn      = dotJ_i(:,6*n-5:6*n);
                              
                M_mn_i     = M_mn_i + Sm' * M * Sn * dX;                                                                % equation 31 
                
                C1_mn_i    = C1_mn_i + Sm'* matrix_coadj(eta_pre_list(:,i)+intdAdjg.*xidot_array(:,i))*M* Sn * dX;      % equation 32
                
                C2_mn_i    = C2_mn_i + Sm'* M * matrix_adj(intdAdjg.*xidot_array(:,i))*Sn * dX;        
                
    %           C2_mn_i    = C2_mn_i + Sm'* M *dotSn * dX;                                                               % equation 33
                
                G_m_i      = G_m_i + Sm' * M * matrix_Adjoint((g_pre_list(:, 4*i-3:4*i)*g_i)^-1)*dX;                     % g is disc                                         % equation 35
                
            end
            
            M_mn  = M_mn  + M_mn_i;
            C1_mn = C1_mn + C1_mn_i;
            C2_mn = C2_mn + C2_mn_i;
            G_m   = G_m  + G_m_i;
        
        end
          
    	%  Actuation and Internal Load

        if t<=tact                                                    % turn                        
          Fan         =[Famx(i);Famy(i);Famz(i);Fax(i);Fay(i);Faz(i)]*t/tact;
        end

        if (t>tact && t<=trel+tact)                                % release / maintenance  
          Fan   =[Famx(i);Famy(i);Famz(i);Fax(i);Fay(i);Faz(i)];
        end

        Fin     = Eps*(xin-xi_star)+Ipsi*xidotn;
        
        % Equation 29  \tau_n = L*(Fan-Fin)
        
        Tau_m   = L*(Fan-Fin);                                      % Equation 29  control quantity      

        Fpn     = [Fpmx(i);Fpmy(i);Fpmz(i);Fpx(i);Fpy(i);Fpz(i)];   % External Concentrated load

        ECL_m     = Sm' *  Fpn;                                     % Equation 30
        
        %  Fill in the integrated block. 
        
        GIM(6*m-5:6*m, 6*n-5:6*n)  = M_mn;   
        GCM1(6*m-5:6*m, 6*n-5:6*n) = C1_mn;
        GCM2(6*m-5:6*m, 6*n-5:6*n) = C2_mn;
        Tau(6*m-5:6*m,1)           = Tau_m;     
        GM(6*m-5:6*m,1:6)          = G_m;                           %  equation 26       
        ECL(6*m-5:6*m,1)           = ECL_m;
                      
     end       
    
end

%-----State equation------%

dotz1       = Xidot;
dotz2       = GIM^-1*(Tau+GM*matrix_Adjoint(g_r^-1)*Gra+ECL-(GCM1+GCM2)*Xidot);
dz          = [dotz1;dotz2];

gv.g     = g;
gv.eta   = eta;
gv.nstep = nstep+1;
gv.GIM   = GIM;    % Generalized inertia matrix
gv.GM    = GM;     % Gravitational matrix
gv.GCM1  = GCM1;   % Generalized Coriolis matrix 1
gv.GCM2  = GCM2;   % Generalized Coriolis matrix 2
gv.Tau   = Tau;    % Internal force and actuation force L(F_a-F_i)
gv.ECL   = ECL;    % External Concentrated Load

end