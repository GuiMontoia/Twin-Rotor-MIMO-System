%% Modelo completo
% clear all;
% Constantes
psi0=40;
phi0=45;
I1=6.8e-2;
I2=2e-2;
a1=0.0135;
b1=0.0924;
a2=0.02;
b2=0.09;
Mg=0.32;
B1psi=6e-3;
B2psi=1e-3;
B1phi=1e-1;
B2phi=1e-2;
Kgy=0.05;
k1=1.1;
k2=0.8;
T11=1.1;
T10=1;
T21=1;
T20=1;
Tp=2;
T0=3.5;
kc=-0.2;
%%
% Variables 
x=sym('x', [7, 1]);
x_diff=sym('x_diff', [7, 1]);
u=sym('u',[2,1]);
assume(u,'real')
% assume(u<2.5 & u>-2.5)

% Condiciones iniciales
x_diff0=[0 0 0 0 0 0 0]';
x(1:4)=[deg2rad(psi0) 0 deg2rad(phi0) 0]';
% Ecuaciones para calcular condiciones iniciales
x_diff(1)=x(2);
x_diff(2) = 1/I1*(a1*x(5)^2*sign(x(5))+b1*x(5)-Mg*sin(x(1))-B1psi*x(2)... 
    -B2psi*sign(x(2))-Kgy*(a1*x(5)^2*sign(x(5))+b1*x(5))*x(4)*cos(x(1)));
x_diff(3)=x(4);
x_diff(4) = 1/I2*(a2*x(6)^2*sign(x(6))+b2*x(6)-B1phi*x(4)... 
    -B2phi*sign(x(4))-x(7));
x_diff(5)=1/T11*(k1*u(1)-T10*x(5));
x_diff(6) = 1/T21*(k2*u(2)-T20*x(6));
% x_diff(7)=1/Tp*(kc*x(5)*(1-T0*T10/T11)+kc*T0*k1*u(1)/T11-x(7));
x_diff(7) = 1/Tp* (-x(7) +kc *T0*x_diff(5)+kc*x(5));
% Condiciones iniciales
[x5_ini, x6_ini, x7_ini, u1_ini, u2_ini] = ...
    solve(x_diff==x_diff0, [x(5) x(6) x(7) u(1) u(2)])
%%
% Redefinimos las variables de estado que antes habian 
% tomado el valor de las condiciones iniciales
x=sym('x', [7, 1]);
x_diff=sym('x_diff', [7, 1]);
u=sym('u',[2,1]);
assume(u,'real')
%assume(x < 0 | x>0) % Sale una delta de Dirac
% Salidas
y=sym('y',[2,1]);
y(1)=x(1);
y(2)=x(3);

x_diff(1)=x(2);
x_diff(2) = 1/I1*(a1*x(5)^2+b1*x(5)-Mg*sin(x(1))-B1psi*x(2)... 
    -B2psi*0-Kgy*(a1*x(5)^2+b1*x(5))*x(4)*cos(x(1)));
x_diff(3)=x(4);
x_diff(4) = 1/I2*(a2*x(6)^2+b2*x(6)-B1phi*x(4)... 
    -B2phi*0-x(7));
x_diff(5) = 1/T11*(k1*u(1)-T10*x(5));
x_diff(6) = 1/T21*(k2*u(2)-T20*x(6));
% x_diff(7)=1/Tp*(kc*x(5)*(1-T0*T10/T11)+kc*T0*k1*u(1)/T11-x(7));
x_diff(7) = 1/Tp* (-x(7) +kc *T0*x_diff(5)+kc*x(5));
% Matrices linealizadas
A=jacobian(x_diff,x);
B=jacobian(x_diff,u);
C=jacobian(y,x);
D=jacobian(y,u);
% Las variables de estado evaluadas
x1=deg2rad(psi0);
x2=0;
x3=deg2rad(phi0);
x4=0;
u1=u1_ini;
u2=u2_ini;
x5=x5_ini;
x6=x6_ini;
x7=x7_ini;
% x1=0; x2=0; x3=0; x4=0; u1=0; u2=0; x5=0; x6=0; x7=0;

% Sustituinos en las derivadas (linealizar)
A_e=double(eval(A))
B_e=double(eval(B))
C_e=double(eval(C))
D_e=double(eval(D))
% Modelo de espacio de estados
states = {'pitch' 'diff_pitch' 'yaw' 'diff_yaw' ...
    'tau1' 'tau2' 'MR'};
inputs = {'u1' 'u2'};
outputs = {'pitch' 'yaw'};

sys_mimo = ss(A_e,B_e,C_e,D_e,'statename',states,...
'inputname',inputs,...
'outputname',outputs);
%%
% Funciones de tranferencia (no sale lo que deberia)
G=tf(sys_mimo);
zpk(minreal(G))

%%

incremento_u_11 = 0.1;
incremento_u_22 = 0.1;
u_11 = double(u1_ini);
u_22 = double(u2_ini);
cond_ini=double([ x1 x2 x3 x4 x5 x6 x7]');