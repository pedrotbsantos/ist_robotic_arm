% Daniel Fortunato (81498), Teodoro Dias (81723), Pedro Santos (84162)

function [J, dejp, dejr] = jacobianDetermination(i, j, w)

syms th0 th1 th2 th3 th4 th5

a1=0; a2=0; a3=40; a4=-120; a5=0; a6=0; a7=20; a8=0; a9=0; a10=0; 
alpha1=0; alpha2=0; alpha3=-sym(pi)/2; alpha4=0; alpha5=0; alpha6=-sym(pi)/2; alpha7=0; alpha8=0; alpha9=sym(pi)/2; alpha10=-sym(pi)/2;
d1=99; d2=0; d3=0; d4=0; d5=0; d6=-40; d7=0; d8=195; d9=0; d10=10; 
theta1=th0; theta2=-pi/2; theta3=th1; theta4=th2; theta5=-pi/2; theta6=0; theta7=th3; theta8=pi/2; theta9=th4; theta10=th5;

% Matrizes de rotação e translação
% Matrizes T01, T12, T45, T56 e T78 são auxiliares para que haja apenas uma rotação em torno de um
% eixo entre frames e para que não haja falta de alguns parâmetros ou para
% que não haja alphas e thetas que não podem ser medidos 

%T12, T45 and T78 hardcoded to avoid errors representing 0 in matrices 
A1=Transform(a1,alpha1,d1,theta1);
T12= [0,1,0,0;-1,0,0,0,0;0,0,1,0;0,0,0,1];
T23=Transform(a3,alpha3,d3,theta3);
T34=Transform(a4,alpha4,d4,theta4);
T45= [0,1,0,0;-1,0,0,0;0,0,1,0;0,0,0,1;];
T56=Transform(a6,alpha6,d6,theta6);
T67=Transform(a7,alpha7,d7,theta7);
T78=[0,-1,0,0;1,0,0,0;0,0,1,d8;0,0,0,1];
T89=Transform(a9,alpha9,d9,theta9);
T910=Transform(a10,alpha10,d10,theta10);

%numbers represented by a maximum of 3 digits
digits(3);

%determine joint transformation matrices
T02=vpa(A1*T12*T23);
T03=vpa(A1*T12*T23*T34);
T04=vpa(A1*T12*T23*T34*T45*T56*T67);
T05=vpa(A1*T12*T23*T34*T45*T56*T67*T78*T89);
T06=vpa(A1*T12*T23*T34*T45*T56*T67*T78*T89*T910)

%R matrices
z0 = A1(1:3,3); 
z1 = T02(1:3,3);
z2 = T03(1:3,3);
z3 = T04(1:3,3);
z4 = T05(1:3,3);
z5 = T06(1:3,3);

%derivatives of position in order of theta
Jp11 = diff(T06(1,4),th0);
Jp12 = diff(T06(1,4),th1);
Jp13 = diff(T06(1,4),th2);
Jp14 = diff(T06(1,4),th3);
Jp15 = diff(T06(1,4),th4);
Jp16 = diff(T06(1,4),th5);

Jp21 = diff(T06(2,4),th0);
Jp22 = diff(T06(2,4),th1);
Jp23 = diff(T06(2,4),th2);
Jp24 = diff(T06(2,4),th3);
Jp25 = diff(T06(2,4),th4);
Jp26 = diff(T06(2,4),th5);

Jp31 = diff(T06(3,4),th0);
Jp32 = diff(T06(3,4),th1);
Jp33 = diff(T06(3,4),th2);
Jp34 = diff(T06(3,4),th3);
Jp35 = diff(T06(3,4),th4);
Jp36 = diff(T06(3,4),th5);

%Linear Velocity Jacobian
Jp = simplify([Jp11, Jp12, Jp13, Jp14, Jp15, Jp16;Jp21, Jp22, Jp23, Jp24, Jp25, Jp26;Jp31, Jp32, Jp33, Jp34, Jp35, Jp36]);

%Angular Velocity Jacobian
Jr = simplify([z0,z1,z2,z3,z4,z5]);

%Jacobian
J = simplify([Jp;Jr])

%Jp = ([ 40*cos(th0) - 120*cos(th0)*cos(th1) - 10*cos(th3)*sin(th0)*sin(th4) + 10*cos(th1 + th2)*cos(th0)*cos(th4) + 10*cos(th0)*sin(th3)*sin(th4)*sin(th1+th2) + 155*cos(th0)*cos(th1+th2) + 20*cos(th0)*sin(th1+th2) , -5.0*sin(th0)*(31*sin(th1+th2) -4*cos(th1+th2)  + 2*cos(th4)*sin(th1+th2) - 2*sin(th4)*sin(th3)*cos(th1+th2)), 5*sin(th0)*(2*cos(th1+th2)*(sin(th3)*sin(th4)+2)-(2*cos(th4)+31)*sin(th1+th2)),     10*sin(th4)*(sin(th0)*cos(th3)*sin(th1+th2)-cos(th0)*sin(th3)),     10*(sin(th0)*(sin(th3)*cos(th4)*sin(th1+th2)-sin(th4)*cos(th1+th2))+cos(th0)*cos(th3)*cos(th4)), 0;...
  %  40*sin(th0) - 120*cos(th1)*sin(th0) + 20*sin(th0)*sin(th1+th2) + 155*sin(th0)*cos(th1+th2) + 10*cos(th0)*cos(th3)*sin(th4) + 10*sin(th0)*sin(th3)*sin(th4)*sin(th1+th2) + 10*cos(th1 + th2)*cos(th4)*sin(th0),  5*cos(th0)*(31*sin(th1+th2) -4*cos(th1+th2) - 24*sin(th1)   + 2*cos(th4)*sin(th1+th2) -2*sin(th4)*sin(th3)*cos(th1+th2)),  5*cos(th0)*(31*sin(th1+th2) - 4*cos(th1+th2) + 2*cos(th4)*sin(th1+th2) -2*sin(th4)*sin(th3)*cos(th1+th2)), -10*sin(th4)*(sin(th0)*sin(th3) + cos(th3)*cos(th0)*sin(th1+th2)), 10*cos(th4)*(cos(th3)*sin(th0) - sin(th3)*cos(th0)*sin(th1+th2)) + 10*sin(th4)*cos(th0)*cos(th1+th2), 0; ...
   % 0, 120*cos(th1) - 20*sin(th1 + th2) - 155*cos(th1 + th2) - 10*sin(th4)*sin(th3)*sin(th1+th2) - 10*cos(th1 + th2)*cos(th4), (-10*sin(th4)*sin(th3)-20)*sin(th1+th2)+(-155-10*cos(th4))*cos(th1+th2), 10*cos(th1 + th2)*cos(th3)*sin(th4), 10*sin(th1 + th2)*sin(th4) + 10*cos(th1 + th2)*cos(th4)*sin(th3), 0])
   
dejp = simplify(det([J(1:3,i),J(1:3,j),J(1:3,w)]))
dejr = simplify(det([J(4:6,i), J(4:6,j), J(4:6,w)]))

end
