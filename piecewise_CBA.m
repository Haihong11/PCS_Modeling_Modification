function dz = piecewise_CBA(t,z)
global gv

t

% global variable
L           =gv.L;
Eps         =gv.Eps;
Ipsi        =gv.Ipsi;
M           =gv.M;
xi_star     =gv.xi_star;
Gra         =gv.Gra;
dX          =gv.dX; 
X           =gv.X;
tact        =gv.tact;  %加载时间
trel        =gv.trel;  %松弛时间
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

% Code for getting parameters.
g_r              =[0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];     % 惯性系与body frame之间变换矩阵
g_prec           =diag([1 1 1 1]); %Initial configuration
eta_prec         =zeros(6,1);      %initial speed 

num_piece = 2;
num_disc  = 20; % This should be big enough.
N         = num_piece;

GIM         = zeros(6*num_piece);        % mass matrix
GCM1        = zeros(6*num_piece);        % Coriolis matrix 1
GCM2        = zeros(6*num_piece);        % Coriolis matrix 2
Tau         =zeros(6*num_npie,1);        % Internal force and actuation force L(F_a-F_i)
GM          =zeros(6*num_npie,6);        % Gravitational matrix
ECL         =zeros(6*num_npie,1);        % External Concentrated Load

Xi          =z(1:6*num_npie,:);  
Xidot       =z(6*num_npie+1:12*num_npie,:);

for jj= 1:num_npie
    xin    = Xi(6*(jj-1)+1:6*(jj-1)+6,:);
    xidotn = Xidot(6*(jj-1)+1:6*(jj-1)+6,:);
    kn     = xin(1:3);
    thetan = sqrt(kn'*kn);
    
for m = 1:num_piece 
    for n = 1:num_piece 
        
        M_mn   = zeros(6);
        C1_mn  = zeros(6);
        C2_mn  = zeros(6);
        Tau_m  = zeros(6,1);
        G_m    = zeros(6,6);
        ECL_m  = zeros(6,1);
        
        for i = max(m,n):N
            
            X = Li;
            
            for ii = 1:num_disc
                
                invAdjg_ii    = piecewise_invAdjoint(X(ii),thetan,xin);
                intdAdjg_ii   = piecewise_ADJ(X(ii),thetan,xin); 
%S=Ad^{-1}_g*T_g    where intdAdjg_ii is T_g。这里如何区分Sm和Sn的表达式，同一个section不同的disc对应的S只与X(ii)有关
                
                Sm = invAdjg_ii*intdAdjg_ii;  % Compute S. where S is the function of the abscissi X
                Sn = invAdjg_ii*intdAdjg_ii;
                
                X = X + dx;                           % The actual position of disc.

                M_mn  = M_mn + Sm' * M* Sn * dX;       % equation 31  
                C1_mn = C1_mn + Sm'* matrix_coadj(eta_prec+intdAdjg_ii*xidotn)*M* Sn * dX;% equation 32
                C2_mn = C2_mn + Sm'* M * matrix_adj(intdAdjg_ii*xidotn)*Sn * dX;          % equation 33
                G_m   = G_n + Sm' * M * invAdjg_ii*dX;     % equation 35
                
            end
            
                % Actuation and internal load
                
                if t<=tact             %turn                        
                     Fan         =[Famx(jj);Famy(jj);Famz(jj);Fax(jj);Fay(jj);Faz(jj)]*t/tact;
                end
                if (t>tact && t<=trel+tact)                     % release / maintenance  
                Fan  =[Famx(jj);Famy(jj);Famz(jj);Fax(jj);Fay(jj);Faz(jj)];
                end
                Fin   =Eps*(xin-xi_star)+Ipsi*xidotn;
  
                Tau_m = L*(Fan-Fin);   %equation 29 
                %equation 29 任意一个section \tau_n等于从该段开始驱动力的总和与该段的内力，再乘该段长度
                % External Concentrated load
                Fpn   =[Fpmx(jj);Fpmy(jj);Fpmz(jj);Fpx(jj);Fpy(jj);Fpz(jj)];
                ECL_m   = Sm' * Fpn;      % equation 30 
        end
        
        GIM(6*m-5:6*m, 6*n-5:6*n)  = M_mn;   % Fill in the integrated block.
        GCM1(6*m-5:6*m, 6*n-5:6*n) = C1_mn;
        GCM2(6*m-5:6*m, 6*n-5:6*n) = C2_mn;
        Tau(6*n-5:6*n,1)           = Tau_m;
        GM(6*n-5:6*n, 6)           = G_m* matrix_Adjoint(g_prec^-1)*Gra;  %equation 26
        ECL(6*n-5:6*n,1)           = ECL_m;
    end
end

%-----------------recursive factors---------------------%

invAdjg1_last    =piecewise_invAdjoint(X(num_disc),thetan,xin);

if num_npie ~= 1
    g_prec          =g_prec*piecewise_expmap(X(num_disc),thetan,xin);
    eta_prec        =invAdjg1_last*(eta_prec+intdAdjg1_last*xidot1);
end


z1       =Xidot;
z2       =GIM^-1*(Tau+GM*matrix_Adjoint(g_r^-1)*Gra+ECL-(GCM1+GCM2)*Xidot);
dz       =[z1;z2];

end

end

gv.GIM  = GIM;   % mass matrix
gv.GM   = GM;    % Gravitational matrix
gv.GCM1 = GCM1;  % Coriolis matrix 1
gv.GCM2 =GCM2;   % Coriolis matrix 2
gv.Tau  =Tau;    % Internal force and actuation force L(F_a-F_i)
gv.ECL =ECL;     % External Concentrated Load