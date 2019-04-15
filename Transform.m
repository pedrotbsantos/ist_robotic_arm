
% Daniel Fortunato (81498), Teodoro Dias (81723), Pedro Santos (84162)

function [Matrix] = Transform(a, alpha, d, theta)

r11=cos(theta);
r12=-sin(theta);
r13=0;
r14=a;
r21=sin(theta)*cos(alpha);
r22=cos(theta)*cos(alpha);
r23=-sin(alpha);
r24=-sin(alpha)*d;
r31=sin(theta)*sin(alpha);
r32=cos(theta)*sin(alpha);
r33=cos(alpha);
r34=cos(alpha)*d;
r41=0;
r42=0;
r43=0;
r44=1;

Matrix = [r11, r12, r13, r14; r21, r22, r23, r24; r31, r32, r33, r34; r41, r42, r43, r44]; 
