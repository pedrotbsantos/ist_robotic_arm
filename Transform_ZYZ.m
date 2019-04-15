
% Daniel Fortunato (81498), Teodoro Dias (81723), Pedro Santos (84162)

function [Matrix_ap_bt_gm] = Transform_ZYZ (x_, y_,z_,ap,bt,gm)

r11=cos(ap)*cos(bt)*cos(gm)-sin(ap)*sin(gm);
r12=-cos(gm)*sin(ap)-cos(ap)*cos(bt)*sin(gm);
r13=cos(ap)*sin(bt);
r14=x_;
r21=sin(ap)*cos(bt)*cos(gm)+cos(ap)*sin(gm);
r22=-sin(ap)*cos(bt)*sin(gm)+cos(ap)*cos(gm);
r23=sin(ap)*sin(bt);
r24=y_;
r31=-sin(bt)*cos(gm);
r32=sin(bt)*sin(gm);
r33=cos(bt);
r34=z_;


Matrix_ap_bt_gm=[r11, r12, r13,r14; r21, r22, r23, r24; r31, r32, r33, r34; 0,0,0,1];

end 