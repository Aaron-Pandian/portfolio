%system parameters
g=32.2;
theta_1=2.4;
u_1=399*1.68780986;%knots to ft/s
W=636636;
m=W/g;
I_xx=1.82e7;
I_zz=4.97e7;
I_xz=9.7e5;
S=5500;
b=195.7;
rho=1.2673e-3;
q=0.5*rho*u_1^2;

C_yb=-0.9;
C_yp=0.0;
C_yr=0.0;
C_yda=0.0;
C_ydr=0.12;
C_lb=-0.16;
C_lp=-0.34;
C_lr=0.13;
C_lda=0.013;
C_ldr=0.008;
C_nb=0.16;
C_nTb=0.0;
C_np=-0.026;
C_nr=-0.28;
C_nTr=0.0;
C_nda=0.0018;
C_ndr=-0.1;





%B747 high cruise stability derivatives
Y_b=u_1*q*S/(m*u_1)*C_yb;
Y_p=q*S*b/(2*m*u_1)*C_yp;
Y_r=q*S*b/(2*m*u_1)*C_yr;
L_b=u_1*q*S*b/(I_xx*u_1)*C_lb;
L_p=q*S*b*b/(2*I_xx*u_1)*C_lp;
L_r=q*S*b*b/(2*I_xx*u_1)*C_lr;
N_b=u_1*q*S*b/(I_zz*u_1)*C_nb;
N_Tb=u_1*q*S*b/(I_zz*u_1)*C_nTb;
N_p=q*S*b/(I_zz*u_1)*C_np;
N_r=q*S*b*b/(2*I_zz*u_1)*C_nr;
N_Tr=q*S*b*b/(2*I_zz*u_1)*C_nTr;
Y_da=q*S/(m*u_1)*C_yda;
Y_dr=u_1*q*S/(m*u_1)*C_ydr;
L_da=q*S*b/I_xx*C_lda;
L_dr=q*S*b/I_xx*C_ldr;
N_da=q*S*b/I_zz*C_nda;
N_dr=q*S*b/I_zz*C_ndr;



c2=cosd(-theta_1)^2;
s2=sind(-theta_1)^2;
sc=sind(-theta_1)*cosd(-theta_1);
I_zzb=I_zz;
I_xxb=I_xx;
I_xzb=I_xz;

%Rotate moments of inertia from body frame to trimmed frame
I_zz=c2*I_zzb+s2*I_xxb+2*I_xzb*sc;
I_xx=c2*I_xxb+s2*I_zzb-2*I_xzb*sc;
I_xz=(I_zzb-I_xxb)*sc+I_xzb*(c2-s2);


%Mass Matrix
M = [  u_1      0               0       0   0 ;
        0       1         -(I_xz/I_xx)  0   0 ;
        0   -(I_xz/I_zz)        1       0   0 ;
        0       0               0       1   0 ;
        0       0               0       0   1];

%Stiffness Matrix
R = [  Y_b      Y_p        Y_r-u_1    g*cosd(theta_1)   0;
       L_b      L_p         L_r              0          0;
     N_b+N_Tb   N_p       N_r+N_Tr           0          0;
        0        1      tand(theta_1)        0          0;
        0        0      secd(theta_1)        0          0];

%Control Sensitivity Matrix
F = [ Y_da  Y_dr;
      L_da  L_dr;
      N_da  N_dr;
      0      0  ;
      0      0  ];

%Compute Matrices for Linearized System
A=M\R;
B=M\F;
C=eye(5);
sys=ss(A,B,C,0);

%Step Response of Linear System
[ystep,tstep]=step(sys,300);
%Impulse Response of Linear System
[yimpulse,timpulse]=impulse(sys,300);

figure(1)
set(gcf,'Color','w')
subplot(5,2,1)
plot(tstep,ystep(:,1,1))
title('\beta Aileron Step Response','FontSize',16)
ylabel('$\Delta \beta$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,3)
plot(tstep,ystep(:,2,1))
title('p Aileron Step Response','FontSize',16)
ylabel('$\Delta p$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,5)
plot(tstep,ystep(:,3,1))
title('r Aileron Step Response','FontSize',16)
ylabel('$\Delta r$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,7)
plot(tstep,ystep(:,4,1))
title('\phi Aileron Step Response','FontSize',16)
ylabel('$\Delta \phi$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,9)
plot(tstep,ystep(:,5,1))
title('\psi Aileron Step Response','FontSize',16)
ylabel('$\Delta \psi$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,2)
plot(tstep,ystep(:,1,2))
title('\beta Rudder Step Response','FontSize',16)
ylabel('$\Delta \beta$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,4)
plot(tstep,ystep(:,2,2))
title('p Rudder Step Response','FontSize',16)
ylabel('$\Delta p$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,6)
plot(tstep,ystep(:,3,2))
title('r Rudder Step Response','FontSize',16)
ylabel('$\Delta r$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,8)
plot(tstep,ystep(:,4,2))
title('\phi Rudder Step Response','FontSize',16)
ylabel('$\Delta \phi$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,10)
plot(tstep,ystep(:,5,2))
title('\psi Rudder Step Response','FontSize',16)
ylabel('$\Delta \psi$','interpreter','latex')
xlabel('Time (s)')
grid on

figure(2)
set(gcf,'Color','w')
subplot(5,2,1)
plot(timpulse,yimpulse(:,1,1))
title('\beta Aileron Impulse Response','FontSize',16)
ylabel('$\Delta \beta$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,3)
plot(timpulse,yimpulse(:,2,1))
title('p Aileron Impulse Response','FontSize',16)
ylabel('$\Delta p$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,5)
plot(timpulse,yimpulse(:,3,1))
title('r Aileron Impulse Response','FontSize',16)
ylabel('$\Delta r$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,7)
plot(timpulse,yimpulse(:,4,1))
title('\phi Aileron Impulse Response','FontSize',16)
ylabel('$\Delta \phi$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,9)
plot(timpulse,yimpulse(:,5,1))
title('\psi Aileron Impulse Response','FontSize',16)
ylabel('$\Delta \psi$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,2)
plot(timpulse,yimpulse(:,1,2))
title('\beta Rudder Impulse Response','FontSize',16)
ylabel('$\Delta \beta$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,4)
plot(timpulse,yimpulse(:,2,2))
title('p Rudder Impulse Response','FontSize',16)
ylabel('$\Delta p$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,6)
plot(timpulse,yimpulse(:,3,2))
title('r Rudder Impulse Response','FontSize',16)
ylabel('$\Delta r$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,8)
plot(timpulse,yimpulse(:,4,2))
title('\phi Rudder Impulse Response','FontSize',16)
ylabel('$\Delta \phi$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,10)
plot(timpulse,yimpulse(:,5,2))
title('\psi Rudder Impulse Response','FontSize',16)
ylabel('$\Delta \psi$','interpreter','latex')
xlabel('Time (s)')
grid on

%% Part 2 
%Mass Matrix
M2 = [  u_1      0               0       0   0 ;
        0       1               0     0      0 ;
        0   -(I_xz/I_zz)        1       0   0 ;
        0       0               0       1   0 ;
        0       0               0       0   1];

%Stiffness Matrix
R2 = [  Y_b      Y_p        Y_r-u_1    g*cosd(theta_1)   0;
       0      L_p         0              0          0;
     N_b+N_Tb   N_p       N_r+N_Tr           0          0;
        0        1      tand(theta_1)        0          0;
        0        0      secd(theta_1)        0          0];

%Control Sensitivity Matrix
F2 = [ Y_da  0;
      L_da  0;
      N_da  N_dr;
      0      0  ;
      0      0  ];

%Compute Matrices for Linearized System
A=M2\R2;
B=M2\F2;
C=eye(5);
sys=ss(A,B,C,0);

%Step Response of Linear System
[ystep,tstep]=step(sys,300);
%Impulse Response of Linear System
[yimpulse,timpulse]=impulse(sys,300);

figure(3);
set(gcf,'Color','w')
subplot(5,2,1)
plot(tstep,ystep(:,2,1))
title('p Aileron Step Response','FontSize',16)
ylabel('$\Delta p$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,2)
plot(timpulse,yimpulse(:,2,1))
title('p Aileron Impulse Response','FontSize',16)
ylabel('$\Delta p$','interpreter','latex')
xlabel('Time (s)')
subplot(5,2,3)
plot(tstep,ystep(:,2,2))
title('p Rudder Step Response','FontSize',16)
ylabel('$\Delta p$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,4)
plot(timpulse,yimpulse(:,2,2))
title('p Rudder Impulse Response','FontSize',16)
ylabel('$\Delta p$','interpreter','latex')
xlabel('Time (s)')
grid on

%% Part 4 

%Mass Matrix
M4 = [  1      0               0       0   0 ;
        0       1               0     0      0 ;
        0   0        1       0   0 ;
        0       0               0       1   0 ;
        0       0               0       0   1];

%Stiffness Matrix
R4 = [  0      0        -u_1    g   0;
       L_b      L_p         L_r              0          0;
     N_b   N_p       N_r           0          0;
        0        1      0        0          0;
        0        0      0        0          0];

%Control Sensitivity Matrix
F1 = [ Y_da  0;
      L_da  0;
      N_da  0;
      0      0  ;
      0      0  ];

F4 = [ Y_da  Y_dr;
      L_da  L_dr;
      N_da  N_dr;
      0      0  ;
      0      0  ];

%Compute Matrices for Linearized System
A=M4\R4;
B=M4\F4;
C=eye(5);
sys=ss(A,B,C,0);

%Step Response of Linear System
[ystep,tstep]=step(sys,300);
%Impulse Response of Linear System
[yimpulse,timpulse]=impulse(sys,300);

figure(4);
set(gcf,'Color','w')
subplot(5,2,1)
plot(tstep,ystep(:,4,1))
title('\phi Aileron Step Response','FontSize',16)
ylabel('$\Delta \phi$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,2)
plot(tstep,ystep(:,4,2))
title('\phi Rudder Step Response','FontSize',16)
ylabel('$\Delta \phi$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,3)
plot(timpulse,yimpulse(:,4,1))
title('\phi Aileron Impulse Response','FontSize',16)
ylabel('$\Delta \phi$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,4)
plot(timpulse,yimpulse(:,4,2))
title('\phi Rudder Impulse Response','FontSize',16)
ylabel('$\Delta \phi$','interpreter','latex')
xlabel('Time (s)')
grid on

%% Part 6

%Mass Matrix
M6 = [  1      0;
        0       1];

%Stiffness Matrix
R6 = [  Y_b/u_1              (Y_r-u_1)/u_1;
     N_b        N_r];

%Control Sensitivity Matrix
F6 = [ Y_da/u_1  0;
      0  N_dr];

%Compute Matrices for Linearized System
A6=M6\R6;
B6=M6\F6;
C6=eye(2);
sys=ss(A6,B6,C6,0);

%Step Response of Linear System
[ystep,tstep]=step(sys,300);
%Impulse Response of Linear System
[yimpulse,timpulse]=impulse(sys,300);

figure(5)
set(gcf,'Color','w')
subplot(5,2,1)
plot(tstep,ystep(:,1,1))
title('\beta Aileron Step Response','FontSize',16)
ylabel('$\Delta \beta$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,2)
plot(tstep,ystep(:,2,1))
title('r Aileron Step Response','FontSize',16)
ylabel('$\Delta r$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,3)
plot(tstep,ystep(:,1,2))
title('\beta Rudder Step Response','FontSize',16)
ylabel('$\Delta \beta$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,4)
plot(tstep,ystep(:,2,2))
title('r Rudder Step Response','FontSize',16)
ylabel('$\Delta r$','interpreter','latex')
xlabel('Time (s)')
grid on

figure(6)
set(gcf,'Color','w')
subplot(5,2,1)
plot(timpulse,yimpulse(:,1,1))
title('\beta Aileron Impulse Response','FontSize',16)
ylabel('$\Delta \beta$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,2)
plot(timpulse,yimpulse(:,2,1))
title('r Aileron Impulse Response','FontSize',16)
ylabel('$\Delta r$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,3)
plot(timpulse,yimpulse(:,1,2))
title('\beta Rudder Impulse Response','FontSize',16)
ylabel('$\Delta \beta$','interpreter','latex')
xlabel('Time (s)')
grid on
subplot(5,2,4)
plot(timpulse,yimpulse(:,2,2))
title('r Rudder Impulse Response','FontSize',16)
ylabel('$\Delta r$','interpreter','latex')
xlabel('Time (s)')
grid on

% Computation 
e = eig(A6);
