close all
format long
clc
clear
warning off
global gv
tic

%-------------------------------------------------------------------------
disp('Pre-processing')

%---------Get geometrical parameters----------%

E        = 110e3;                          % [Pa]   Young modulus
mu       = 5e3;                            % [Pa*s] Shear viscosity modulus
Poi      = 0;                              % [-]    Poisson
G        = E/(2*(1+Poi));                  % [Pa]   Shearing modulus
R        = 10e-3;                          % [m]    Radius of cross-section
L        = 250e-3;                         % [m]    Length of soft robot
num_piece= 2;                              % number of piece
num_disc = 20;                             % This should be big enough.                        
X        = linspace(0,L,num_piece+1);      % End point;
A        = pi*R^2;                         % [m^2] Cross-section area
J        = pi*R^4/4;                       % [m^4] The second moment of the area
I        = pi*R^4/2;                       % [m^4] The polar meomet of the area
ro_arm   = 2000;                           % [kg/m^3]
Gra      = [0;0;0;-9.81;0;0];              % [m/s^2]gravity acceleration twist（relative to spatial frame）
Eps      = diag([G*I E*J E*J E*A G*A G*A]);% stifness matrix
Ipsi     = mu*diag([I 3*J 3*J 3*A A A]);   % viscosity matrix
M        = ro_arm*diag([I J J A A A]);     % screw inertia matrix

%-----------The undeformed strain twist-----%

xi_star    =[0;0;0;1;0;0];

% ---------Actuation force (body frame)------%

tact        =0.5;                      % [s] torque time in dir z o y 
trel        =0.5;                      % [s] relaxation time 
Fax         =0*[0 0 0 0 0 0];              % [N] contraction load
Fay         =0*[0 0 0 0 0 0];              % [N] lateral y load 0.1
Faz         =0*[0 0 0 0 0 0];              % [N] lateral z load 0.01
Famx        =0*[0 0 0 0 0 0];              % [Nm] torsion torque 0.001
Famy        =0*[0 0 0 0 0 0];              % [Nm] bending torque
Famz        =0*[0 0 0 0 0 0];              % [Nm] bending torque 0.005

% -------External tip load (base (X=0) frame)-------%

Fpx         =0*[0 0 0 0 0 0];              % [N] contraction load
Fpy         =0.01*[0 5 0 0 0 0];           % [N] lateral y load
Fpz         =0.01*[0 0 0 0 0 0];           % [N] lateral z load
Fpmx        =0*[0 0 0 0 0 0];              % [Nm] torsion torque
Fpmy        =0*[0 0 0 0 0 0];              % [Nm] bending torque
Fpmz        =0*[0 0 0 0 0 0];              % [Nm] bending torque

%--------------Numerical setting---------------%

time        =1;                       % [s]
nsol        =time*10^2+2;                % Number of solution      
tspan       =linspace(0,time,nsol);     
dX          =L/num_piece/num_disc;       % Length of disc

% ------------Velocity and Position-------------%

g           =zeros(4*nsol,4*num_disc*num_piece);
eta         =zeros(6*nsol,num_disc*num_piece);
nstep       =1;

%--------------global variable---------------%

gv.ro_arm      = ro_arm;
gv.Gra         = Gra;
gv.L           = L;
gv.X           = X;
gv.R           = R;
gv.xi_star     = xi_star;
gv.A           = A;
gv.Eps         = Eps;
gv.Ipsi        = Ipsi;
gv.M           = M;
gv.nsol        = nsol;
gv.num_disc    = num_disc;
gv.num_piece   = num_piece;
gv.dX          = dX;
gv.time        = time;
gv.tspan       = tspan;
gv.tact        = tact;
gv.trel        = trel;
gv.Fax         = Fax;
gv.Fay         = Fay;
gv.Faz         = Faz;
gv.Famx        = Famx;
gv.Famy        = Famy;
gv.Famz        = Famz;
gv.Fpx         = Fpx;
gv.Fpy         = Fpy;
gv.Fpz         = Fpz;
gv.Fpmx        = Fpmx;
gv.Fpmy        = Fpmy;
gv.Fpmz        = Fpmz;
gv.g           = g;
gv.eta         = eta;
gv.nstep       = nstep;
gv.num_piece   = num_piece;
gv.num_disc    = num_disc;

disp('Time-advancing')

options        = odeset('RelTol',1e-3,'AbsTol',1e-3);

%------------------Initial conditions---------------%

xi_0           = [0;0;0;1;0;0];  % given variable of the joint, column vector
xidot_0        = [0;0;0;0;0;0];
          
ini_cond       = [repmat(xi_0',[1,num_piece]) repmat(xidot_0',[1,num_piece])];  % whether transpose or not

%------------------Integrate------------------%

[t,z]          = ode45(@piecewise_CBA,tspan,ini_cond,options);

toc

disp('Post-processing')

% mkdir('data');     % Create new directory 
% 
% % %----------Define 6 strain components----------%
% 
% s         =zeros(nsol,num_piece);
% v         =zeros(nsol,num_piece);
% k         =zeros(nsol,num_piece);
% q         =zeros(nsol,num_piece);
% p         =zeros(nsol,num_piece);
% r         =zeros(nsol,num_piece);
% 
% for ii=1:num_piece    
%     s(:,ii)   =z(:,6*(ii-1)+1); 
%     v(:,ii)   =z(:,6*(ii-1)+2);
%     k(:,ii)   =z(:,6*(ii-1)+3);
%     q(:,ii)   =z(:,6*(ii-1)+4);
%     p(:,ii)   =z(:,6*(ii-1)+5);
%     r(:,ii)   =z(:,6*(ii-1)+6);    
% end
% xidot    =z(:,6*num_piece+1:12*num_piece);

% %-----------Visualization-----------%
% for ii=1:num_piece
%  
%     figure
%     plot(t,s(:,ii))
%     grid on
%     title(strcat('torsion of piece',num2str(ii)))
%     xlabel('t [s]')
%     ylabel('s [1/m]')
%     print('-djpeg')
% 
%     figure
%     plot(t,v(:,ii))
%     grid on
%     title(strcat('curvature on y of piece',num2str(ii)))
%     xlabel('t [s]')
%     ylabel('v [1/m]')
%     print('-djpeg')
% 
%     figure
%     plot(t,k(:,ii))
%     grid on
%     title(strcat('curvature on z of piece',num2str(ii)))
%     xlabel('t [s]')
%     ylabel('k [1/m]')
%     print('-djpeg')
% 
%     figure
%     plot(t,q(:,ii))
%     grid on
%     title(strcat('longitudinal strain of piece',num2str(ii)))
%     xlabel('t [s]')
%     ylabel('q [1]')
%     print('-djpeg')
% 
%     figure
%     plot(t,p(:,ii))
%     grid on
%     title(strcat('tras y strain of piece',num2str(ii)))
%     xlabel('t [s]')
%     ylabel('p [1]')
%     print('-djpeg')
%  
%     figure
%     plot(t,r(:,ii))
%     grid on
%     title(strcat('tras z strain of piece',num2str(ii)))
%     xlabel('t [s]')
%     ylabel('r [1]')
%     print('-djpeg')
%     
% end

%----------------Video-------------------%
 
% save('data\t-z','t','z')
% save('data\strain rate','t','xidot');
% save('data\configuration','t','g');
% save('data\velocity','t','eta');
% save('data\torsion','t','s');
% save('data\curvature on y','t','v');
% save('data\curvature on z ','t','k');
% save('data\longitudinal strain','t','q');
% save('data\tras y strain','t','p');
% save('data\tras z strain','t','r');
% 
% video         = VideoWriter(strcat('data\Dynamics'));
% FrameRate     = 10^2;        % FPS
% open(video)
% 
% %--------screen resolution-----------% 
% scrsz         = get(0,'ScreenSize');  % 1  1  1536  864
% % scrsz(3)--1280； scrsz(4)---800
% figure('Position',[scrsz(3)/12 2*scrsz(4)/48 11*scrsz(3)/6 9*scrsz(4)/10])
% 
% %[left, bottom, width, height]
% 
% angle          = linspace(0,2*pi,180);
% 
% for ii=1:nsol       % for every moment
%  % clf  
%     g1         = g(4*(ii-1)+1:4*(ii-1)+4,:);
%     
%  %--------Camera position, shooting point and perspective-----------%
%     
%     set(gca,'CameraPosition',[0 0 -L],'CameraTarget',[0 0 0],'CameraUpVector',[1 0 0])
%     axis equal
%     grid on
%     hold on
%     xlabel('X [m]')
%     ylabel('Y [m]')
%     zlabel('Z [m]')
%     title(strcat('t= ',num2str(t(ii))))
%     
%     % cantilever
%     axis ([-num_piece*L num_piece*L 0 1.5*num_piece*L -L L])  
%     
%     % drawing the soft mamipulator
%     
%     robot    = [zeros(1,180) 0;R*sin(angle) 0;R*cos(angle) 0;ones(1,180) 1];
%     
%     for zz = 1:num_piece
%         
%         for jj = 1:num_disc
%         
%             robot1  = g1(:,4*num_disc*(zz-1)+4*(jj-1)+1:4*num_disc*(zz-1)+4*(jj-1)+4)*robot;
%             plot3(robot1(1,:),robot1(2,:),robot1(3,:),'Color',[1-mod(zz,2),0,mod(zz,2)])
%         
%         end
%         drawnow
%     end
%     drawnow 
%     
%     F   = getframe(gcf);  % get image 
%     
%     writeVideo(video,F);
% end
% 
% close(video);