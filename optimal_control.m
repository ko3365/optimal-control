clc;clear;close;

% System Initialization
M = 2; m = 1; l = 4; g = 9.8;
A = [0 1 0; (M+m)*g/(M*l) 0 0; -m*g/M 0 0];
B = [0; -1/(M*l);1/M];
L = [0;1/(m*l);0];
C = [0 0 1];
%xuo = -inv([A B;C 0])*[L;0];
%thetao = xuo(1)*1
%uo = xuo(4)*1

%LQR Control Gain Calculation
Q = C'*C;
[Ko_bar,X,~] = lqr(A,B,Q,1);
Ko = -1*Ko_bar;
K1 = [Ko -1]*inv([A B; C 0])*[L;0];

%Kalman Gain Calculation
Ak = [A L; zeros(1,4)];
Bk = [0;0;0;1];
Ck = [C 0];
V = 1;
W = Bk*Bk';
[F_bar,Y,~] = lqr(Ak',Ck',W,V);
F = -F_bar';

% Basic Design
K = [Ko K1];
Bcl = [B;0];
x0 = [0;0;0;1;0;0;0;0];
model_basic = ss([Ak Bcl*K; -F*Ck Ak+Bcl*K+F*Ck],[],[],[]);
[~,t_simul,x_simul] = initial(model_basic,x0,20);
u = zeros(length(x_simul),1);
for j=1:length(x_simul)
    u(j) = K*x_simul(j,5:8)';
end
figure(1)
plot(t_simul,[x_simul(:,3),20*x_simul(:,4),20*x_simul(:,8),180*x_simul(:,1)/pi,u]);
grid on
legend('$\dot{x_c}$','20$d$','20$\hat{d}$','$\theta$','$u$','interpreter','latex');

% Design with Additional Tuning Parameter
alpha = 0.16663;
beta = 0.08;

Q = C'*C;
R = alpha^2;
W = Bk*Bk';
V = beta^2;

[ko,~,~] = lqr(A,B,Q,R);
ko_design = -ko;
k1_design = [ko_design -1]*inv([A B; C 0])*[L;0];
K_design = [ko_design, k1_design];
[f_bar,~,~] = lqr(Ak',Ck',W,V);
F_design = -f_bar';

A_d = [A B*K_design; -F_design*C Ak+Bcl*K_design+F_design*Ck];
B_d = [L zeros(3,1); zeros(4,1) -beta*F_design];
C_d = [C zeros(1,4); zeros(1,3) alpha*K_design];
D_d = [0 0 ; 0 0];

model_design = ss(A_d,B_d,C_d,D_d);
t = linspace(0,20,1000);
u1 = ones(1,length(t));
u2 = 0.01*sin(20*pi*t)/beta;
x0 = zeros(7,1);
[yy,tt,xx] = lsim(model_design,[u1' u2'], t, x0);

xc = zeros(1,length(tt));
dt = 0.02;
for i=2:length(tt)
    xc(i) = xc(i-1)+xx(i,3)*dt;
end
figure(2)
plot(tt,[xx(:,3),8*u1',8*xx(:,7),180*xx(:,1)/pi,yy(:,2)/alpha,xc']);
xlim([0 10])
grid on
legend('$\dot{x_c}$','8$d$','8$\hat{d}$','$\theta$','$u$','$x_c$','interpreter','latex');

figure(3)
plot(tt,yy(:,2)/alpha);
xlim([5 15])
ylim([-1.3,-0.7])
grid on

