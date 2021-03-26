function dz = piecewise_derivatives(t,z)
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
num_disc    =gv.num_disc;
num_piece   =gv.num_piece;
tact        =gv.tact;
trel        =gv.trel;
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

% 初始化动力学方程系数矩阵参数
% disp

GIM      =zeros(6*num_piece,6*num_piece);  %%广义的惯性矩阵 generalized inertia matrix（GIM）
GCM1     =zeros(6*num_piece,6*num_piece);  %% 广义科氏力矩阵1 generalized Coriolis matrix (GCM1)
GCM2     =zeros(6*num_piece,6*num_piece);  %% 广义科氏力矩阵2
Tau      =zeros(6*num_piece,1);       % 内力和驱动力 L(F_a-F_i)
GM       =zeros(6*num_piece,6);       % 重力矩阵 Gravitational matrix
ECL      =zeros(6*num_piece,1);       % 外部集中力 （External Concentrated Load）

% 初始化运动学参数

g_r              =[0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];     % 惯性系与body frame之间变换矩阵
Jaco_prec        =diag([1 1 1 1 1 1 zeros(1, 6*(num_piece-1))]);  % 初始化雅可比矩阵
g_prec           =diag([1 1 1 1]); %初始位型
eta_prec         =zeros(6,1);     %初始化速度 
adetan_prec      =zeros(6*num_piece,6*num_piece); 

%-------------------------------------------------------------------------
%计算每个piece的系数矩阵

% 计算第一段的质量矩阵和科氏力矩阵

Xi               =z(1:6*num_piece,:);  
Xidot            =z(6*num_piece+1:12*num_piece,:);
xi1              =Xi(1:6,:);  %%第一个关节位置 所有列代表时间
xidot1           =Xidot(1:6,:);
k1               =xi1(1:3);   %旋量xi中的角应变
theta1           =sqrt(k1'*k1);  %公式 theta^2=(k1)^2+(k2)^2+(k3)^3

MasX             =zeros(6,6*num_disc); 
LMasX            =zeros(6,6*num_disc); 
LRMasX           =zeros(6,6*num_disc); %质量
LRCo1X           =zeros(6,6*num_disc);  %科氏力
Mas_prec         =zeros(6,6); %初始化disc质量
LMas_prec        =zeros(6,6);
LRMas_prec       =zeros(6,6);
LRCo1_prec       =zeros(6,6);
for ii=1:num_disc
    coAdjg1_here                     =piecewise_coAdjoint(X(ii),theta1,xi1);
    invAdjg1_here                    =piecewise_invAdjoint(X(ii),theta1,xi1);
    intdAdjg1_here                   =piecewise_ADJ(X(ii),theta1,xi1);    
    
    Mas_here                         =coAdjg1_here*M*invAdjg1_here;  %一个disc的mass matrix
    trapz                            =dX*(Mas_prec+Mas_here)/2;  
    MasX(:,6*(ii-1)+1:6*num_disc)    =MasX(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]); 
    Mas_prec                         =Mas_here;
  
    LRMas_here                       =intdAdjg1_here'*Mas_here*intdAdjg1_here; % S^T_n*M*S_n展开 其中S_n=invAdjg*intdAdjg
    trapz                            =dX*(LRMas_prec+LRMas_here)/2; 
    LRMasX(:,6*(ii-1)+1:6*num_disc)  =LRMasX(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]); %上三角求和
    LRMas_prec                       =LRMas_here;
    
    LMas_here                        =intdAdjg1_here'*Mas_here;  %求重力矩阵要用到。实际上它等于S^T_nM {Ad^{-1}_g}论文公式35
    trapz                            =dX*(LMas_prec+LMas_here)/2;
    LMasX(:,6*(ii-1)+1:6*num_disc)   =LMasX(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]);
    LMas_prec                        =LMas_here;

    % coriolises 1   论文公式32   科氏力的系数矩阵有两项，研究第一段的时候没有考虑coriolises 2
    LRCo1_here                       =intdAdjg1_here'*matrix_coadj(eta_prec+intdAdjg1_here*xidot1)*Mas_here*intdAdjg1_here;
    trapz                            =dX*(LRCo1_prec+LRCo1_here)/2;
    LRCo1X(:,6*(ii-1)+1:6*num_disc)  =LRCo1X(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]);
    LRCo1_prec                       =LRCo1_here;
end

LRMasX           =LRMasX-repmat(LRMasX(:,1:6),[1,num_disc]);  %%
LRMas            =LRMasX(:,6*(num_disc-1)+1:6*num_disc);

LMasX            =LMasX-repmat(LMasX(:,1:6),[1,num_disc]);  %求和
LMas             =LMasX(:,6*(num_disc-1)+1:6*num_disc); 
                                           
LRCo1X           =LRCo1X-repmat(LRCo1X(:,1:6),[1,num_disc]);
LRCo1            =LRCo1X(:,6*(num_disc-1)+1:6*num_disc);

% Actuation load internal load and tip load of the first
% piece这里是驱动力和内力，可以不用考虑

if t<=tact                                      % turn
    Fa1         =[Famx(1);Famy(1);Famz(1);Fax(1);Fay(1);Faz(1)]*t/tact;
end
if (t>tact && t<=trel+tact)                     % release / hold 
    Fa1         =[Famx(1);Famy(1);Famz(1);Fax(1);Fay(1);Faz(1)];
end
Fi1             =Eps*(xi1-xi_star)+Ipsi*xidot1;   %%计算内力
Fp1             =[Fpmx(1);Fpmy(1);Fpmz(1);Fpx(1);Fpy(1);Fpz(1)];
    
% dynamic coefficients update
invAdjg1_last    =piecewise_invAdjoint(X(num_disc),theta1,xi1); %这里实际上是个6×6的矩阵，在子程序中有其计算公式
invAdjg1R_last   =blkdiag(invAdjg1_last(1:3,1:3),invAdjg1_last(4:6,4:6));%不清楚这里为什么要写成对角矩阵块
intdAdjg1_last   =piecewise_ADJ(X(num_disc),theta1,xi1);  %这里也是个6×6的矩阵，在子程序中有其计算公式
MasB             =blkdiag(LRMas,zeros(6*(num_piece-1),6*(num_piece-1)));
GIM              =GIM+Jaco_prec'*MasB*Jaco_prec;  %在上面已经乘过一次S_n了，S_n是雅可比矩阵中的项，每看懂这里为什么要这样写。
Co1B             =blkdiag(LRCo1,zeros(6*(num_piece-1),6*(num_piece-1))); %同上
GCM1             =GCM1+Jaco_prec'*Co1B*Jaco_prec;
GraB             =[LMas;zeros(6*(num_piece-1),6)];
GM               =GM+Jaco_prec'*GraB*matrix_Adjoint(g_prec^-1);
Tau              =Tau+[L*(Fa1-Fi1);zeros(6*(num_piece-1),1)]; 
TipB             =[(invAdjg1_last*intdAdjg1_last)'*(invAdjg1R_last*Fp1);zeros(6*(num_piece-1),1)];
ECL              =ECL+Jaco_prec'*TipB;

%----------------------------------------------------------------------
% recursive factors
if num_piece ~= 1
    Jaco_prec       =blkdiag(invAdjg1_last*intdAdjg1_last,zeros(6*(num_piece-1),6*(num_piece-1)))*Jaco_prec+...
                     blkdiag(zeros(6,6),diag([1 1 1 1 1 1]),zeros(6*(num_piece-2),6*(num_piece-2)));   %%%如果微段不等于1，这个雅可比更新为什么是这样的
    g_prec          =g_prec*piecewise_expmap(X(num_disc),theta1,xi1);
    eta_prec        =invAdjg1_last*(eta_prec+intdAdjg1_last*xidot1);
end

%--------------------------------------------------------------------------
% masses, coriolises 1, coriolises 2 from the second piece onwards
for jj=2:num_piece                   
    xin            =Xi(6*(jj-1)+1:6*(jj-1)+6,:);
    xidotn         =Xidot(6*(jj-1)+1:6*(jj-1)+6,:);
    kn             =xin(1:3);
    thetan         =sqrt(kn'*kn);
    
    MasX            =zeros(6,6*num_disc);
    LMasX           =zeros(6,6*num_disc);  %左
    RMasX           =zeros(6,6*num_disc);  %右
    LRMasX          =zeros(6,6*num_disc); 
    Co1X            =zeros(6,6*num_disc);
    LCo1X           =zeros(6,6*num_disc);
    RCo1X           =zeros(6,6*num_disc);
    LRCo1X          =zeros(6,6*num_disc);
    Co2X            =zeros(6,6*num_disc);
    LCo2X           =zeros(6,6*num_disc);
    Mas_prec        =zeros(6,6);
    LMas_prec       =zeros(6,6);
    RMas_prec       =zeros(6,6);
    LRMas_prec      =zeros(6,6);
    Co1_prec        =zeros(6,6);
    LCo1_prec       =zeros(6,6);
    RCo1_prec       =zeros(6,6);
    LRCo1_prec      =zeros(6,6);
    Co2_prec        =zeros(6,6);
    LCo2_prec       =zeros(6,6);
    for ii=1:num_disc
        coAdjgn_here                     =piecewise_coAdjoint(X(ii),thetan,xin);
        invAdjgn_here                    =piecewise_invAdjoint(X(ii),thetan,xin);
        intdAdjgn_here                   =piecewise_ADJ(X(ii),thetan,xin);
        
        % masses
        Mas_here                         =coAdjgn_here*M*invAdjgn_here;
        trapz                            =dX*(Mas_prec+Mas_here)/2;
        MasX(:,6*(ii-1)+1:6*num_disc)    =MasX(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]);
        Mas_prec                         =Mas_here;
        
        LMas_here                        =intdAdjgn_here'*Mas_here; 
        trapz                            =dX*(LMas_prec+LMas_here)/2;
        LMasX(:,6*(ii-1)+1:6*num_disc)   =LMasX(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]);
        LMas_prec                        =LMas_here;  
        
        RMas_here                        =Mas_here*intdAdjgn_here;  %%与LMas_here连接    这边应该是一个disc的左截面质量
        trapz                            =dX*(RMas_prec+RMas_here)/2;
        RMasX(:,6*(ii-1)+1:6*num_disc)   =RMasX(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]);
        RMas_prec                        =RMas_here;
                
        LRMas_here                       =intdAdjgn_here'*Mas_here*intdAdjgn_here;  %%左右截面结合
        trapz                            =dX*(LRMas_prec+LRMas_here)/2;
        LRMasX(:,6*(ii-1)+1:6*num_disc)  =LRMasX(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]);
        LRMas_prec                       =LRMas_here;
        
        % coriolises 1
        Co1_here                         =matrix_coadj(eta_prec+intdAdjgn_here*xidotn)*Mas_here;
        trapz                            =dX*(Co1_prec+Co1_here)/2;
        Co1X(:,6*(ii-1)+1:6*num_disc)    =Co1X(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]);
        Co1_prec                         =Co1_here;
        LCo1_here                        =intdAdjgn_here'*Co1_here;
        trapz                            =dX*(LCo1_prec+LCo1_here)/2;
        LCo1X(:,6*(ii-1)+1:6*num_disc)   =LCo1X(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]);
        LCo1_prec                        =LCo1_here;
        RCo1_here                        =Co1_here*intdAdjgn_here;
        trapz                            =dX*(RCo1_prec+RCo1_here)/2;
        RCo1X(:,6*(ii-1)+1:6*num_disc)   =RCo1X(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]);
        RCo1_prec                        =RCo1_here;
        LRCo1_here                       =intdAdjgn_here'*Co1_here*intdAdjgn_here;  %%第一科氏力矩阵（第二段）
        trapz                            =dX*(LRCo1_prec+LRCo1_here)/2;
        LRCo1X(:,6*(ii-1)+1:6*num_disc)  =LRCo1X(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]);
        LRCo1_prec                       =LRCo1_here;
        
        % coriolises 2 
        Co2_here                         =Mas_here*matrix_adj(intdAdjgn_here*xidotn);
        trapz                            =dX*(Co2_prec+Co2_here)/2;
        Co2X(:,6*(ii-1)+1:6*num_disc)    =Co2X(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]);
        Co2_prec                         =Co2_here;
        LCo2_here                        =intdAdjgn_here'*Co2_here;
        trapz                            =dX*(LCo2_prec+LCo2_here)/2;
        LCo2X(:,6*(ii-1)+1:6*num_disc)   =LCo2X(:,6*(ii-1)+1:6*num_disc)+repmat(trapz,[1,num_disc-(ii-1)]);
        LCo2_prec                        =LCo2_here;
    end
    MasX            =MasX-repmat(MasX(:,1:6),[1,num_disc]);
    LMasX           =LMasX-repmat(LMasX(:,1:6),[1,num_disc]);
    RMasX           =RMasX-repmat(RMasX(:,1:6),[1,num_disc]);
    LRMasX          =LRMasX-repmat(LRMasX(:,1:6),[1,num_disc]);
    
    Co1X            =Co1X-repmat(Co1X(:,1:6),[1,num_disc]);
    LCo1X           =LCo1X-repmat(LCo1X(:,1:6),[1,num_disc]);
    RCo1X           =RCo1X-repmat(RCo1X(:,1:6),[1,num_disc]);
    LRCo1X          =LRCo1X-repmat(LRCo1X(:,1:6),[1,num_disc]);
    Co2X            =Co2X-repmat(Co2X(:,1:6),[1,num_disc]);
    LCo2X           =LCo2X-repmat(Co2X(:,1:6),[1,num_disc]);
    Mas             =MasX(:,6*(num_disc-1)+1:6*num_disc);
    LMas            =LMasX(:,6*(num_disc-1)+1:6*num_disc);
    RMas            =RMasX(:,6*(num_disc-1)+1:6*num_disc);
    LRMas           =LRMasX(:,6*(num_disc-1)+1:6*num_disc);
    Co1             =Co1X(:,6*(num_disc-1)+1:6*num_disc);
    LCo1            =LCo1X(:,6*(num_disc-1)+1:6*num_disc);
    RCo1            =RCo1X(:,6*(num_disc-1)+1:6*num_disc);
    LRCo1           =LRCo1X(:,6*(num_disc-1)+1:6*num_disc);
    Co2             =Co2X(:,6*(num_disc-1)+1:6*num_disc);
    LCo2            =LCo2X(:,6*(num_disc-1)+1:6*num_disc);
    
    % Actuation and internal load
    if t<=tact                                     
        Fan         =[Famx(jj);Famy(jj);Famz(jj);Fax(jj);Fay(jj);Faz(jj)]*t/tact;
    end
    if (t>tact && t<=trel+tact)                     % release / maintenance  
        Fan         =[Famx(jj);Famy(jj);Famz(jj);Fax(jj);Fay(jj);Faz(jj)];
    end
    Fin             =Eps*(xin-xi_star)+Ipsi*xidotn;
    Fpn             =[Fpmx(jj);Fpmy(jj);Fpmz(jj);Fpx(jj);Fpy(jj);Fpz(jj)];
    
    % Actuation and internal load subsequent
    if jj~=num_piece
        if t<=tact                                      % tack
            Fan_suc     =[Famx(jj+1);Famy(jj+1);Famz(jj+1);Fax(jj+1);Fay(jj+1);Faz(jj+1)]*t/tact;
        end
        if (t>tact && t<=trel+tact)                     % release / maintenance
            Fan_suc     =[Famx(jj+1);Famy(jj+1);Famz(jj+1);Fax(jj+1);Fay(jj+1);Faz(jj+1)];
        end
    else
        Fan_suc     =[0;0;0;0;0;0];
    end
    
    % 动力学系数更新
    
    invAdjgn_last   =piecewise_invAdjoint(X(num_disc),thetan,xin);
    invAdjgnR_last  =blkdiag(invAdjgn_last(1:3,1:3),invAdjgn_last(4:6,4:6));
    invAdjgprec     =matrix_Adjoint(g_prec^-1);
    invAdjgprecR    =blkdiag(invAdjgprec(1:3,1:3),invAdjgprec(4:6,4:6));
    intdAdjgn_last  =piecewise_ADJ(X(num_disc),thetan,xin);
    MasB            =blkdiag([repmat(Mas,[jj-1 jj-1]) repmat(RMas,[jj-1 1]); repmat(LMas,[1 jj-1]) LRMas],...
                     zeros(6*(num_piece-jj),6*(num_piece-jj)));
    GIM             =GIM+Jaco_prec'*MasB*Jaco_prec;
    Co1B            =blkdiag([repmat(Co1,[jj-1 jj-1]) repmat(RCo1,[jj-1 1]); repmat(LCo1,[1 jj-1]) LRCo1],...
                     zeros(6*(num_piece-jj),6*(num_piece-jj)));
    GCM1            =GCM1+Jaco_prec'*Co1B*Jaco_prec;
    Co2B            =blkdiag([repmat(Co2,[jj-1 jj-1]) zeros(6*(jj-1),6); repmat(LCo2,[1 jj-1]) zeros(6,6)],...
                     zeros(6*(num_piece-jj),6*(num_piece-jj)))+MasB*adetan_prec;
    GCM2            =GCM2+Jaco_prec'*Co2B*Jaco_prec;
    GraB            =[repmat(Mas,[jj-1 1]);LMas;zeros(6*(num_piece-jj),6)];
    GM              =GM+Jaco_prec'*GraB*matrix_Adjoint(g_prec^-1);
    ForB            =[repmat(invAdjgn_last'*(Fan-Fan_suc),[jj-1 1]);(invAdjgn_last*intdAdjgn_last)'*(Fan-Fan_suc);...
                      zeros(6*(num_piece-jj),1)];                                               
    Tau             =Tau+Jaco_prec'*ForB-[zeros(6*(jj-1),1);L*Fin;zeros(6*(num_piece-jj),1)];

    TipB            =[repmat(invAdjgn_last'*(invAdjgnR_last*invAdjgprecR*Fpn),[jj-1 1]);...
                     (invAdjgn_last*intdAdjgn_last)'*(invAdjgnR_last*invAdjgprecR*Fpn);zeros(6*(num_piece-jj),1)];
    ECL             =ECL +Jaco_prec'*TipB;
    
    %----------------------------------------------------------------------
    % recursive factors
    prec2prec       =[];
    for ii=1:jj-1
        prec2prec   =blkdiag(prec2prec,invAdjgn_last);
    end
    Jaco_prec       =blkdiag(prec2prec,invAdjgn_last*intdAdjgn_last,zeros(6*(num_piece-jj),6*(num_piece-jj)))*Jaco_prec+...
                     blkdiag(zeros(6*jj,6*jj),repmat(eye(6),[jj~=num_piece jj~=num_piece]),zeros(6*(num_piece-jj-1),6*(num_piece-jj-1)));
    g_prec          =g_prec*piecewise_expmap(X(num_disc),thetan,xin);
    ADxin           =intdAdjgn_last*xidotn;
    eta_prec        =invAdjgn_last*(eta_prec+ADxin);
    prec2prec       =[];
    prec2prec_inv   =[];
    for ii=1:num_piece
        prec2prec     =blkdiag(prec2prec,invAdjgn_last);
        prec2prec_inv =blkdiag(prec2prec_inv,invAdjgn_last^-1);
        for zz=1:num_piece
            if (ii+zz == jj)
                adetan_prec(6*(ii-1)+1:6*(ii-1)+6,6*(zz-1)+1:6*(zz-1)+6)  =matrix_adj(ADxin);
            end
        end
    end
    adetan_prec     =prec2prec*adetan_prec*prec2prec_inv;
end

% calculate the derivatives

z1       =Xidot;
z2       =GIM^-1*(Tau+GM*matrix_Adjoint(g_r^-1)*Gra+ECL-(GCM1+GCM2)*Xidot);
dz       =[z1;z2];

gv.GIM  = GIM;
gv.GM   = GM;
gv.GCM1 = GCM1;
gv.GCM2 =GCM2;
gv.Tau  =Tau;
gv.ECL =ECL;