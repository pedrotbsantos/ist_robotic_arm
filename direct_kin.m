% Cinemática direta 

% Daniel Fortunato (81498), Teodoro Dias (81723), Pedro Santos (84162)

function [x, y, z, oa, ob, oc] = direct_kin(th0, th1, th2, th3, th4, th5)


sym(pi)

a1=0; a2=0; a3=40; a4=-120; a5=0; a6=0; a7=20; a8=0; a9=0; a10=0; 
alpha1=0; alpha2=0; alpha3=-pi/2; alpha4=0; alpha5=0; alpha6=-pi/2; alpha7=0; alpha8=0; alpha9=pi/2; alpha10=-pi/2;
d1=99; d2=0; d3=0; d4=0; d5=0; d6=-40; d7=0; d8=195; d9=0; d10=10; 
theta1=th0; theta2=-pi/2; theta3=th1; theta4=th2; theta5=-pi/2; theta6=0; theta7=th3; theta8=pi/2; theta9=th4; theta10=th5;

% Matrizes de rotação e translação
% Matrizes T12, T45, T56 e T78 são auxiliares para que haja apenas uma rotação em torno de um
% eixo entre frames e para que não haja falta de alguns parâmetros ou para
% que não haja alphas e thetas que não podem ser medidos 

T01=Transform(a1,alpha1,d1,theta1);
T12=Transform(a2,alpha2,d2,theta2);
T23=Transform(a3,alpha3,d3,theta3);
T34=Transform(a4,alpha4,d4,theta4);
T45=Transform(a5,alpha5,d5,theta5);
T56=Transform(a6,alpha6,d6,theta6);
T67=Transform(a7,alpha7,d7,theta7);
T78=Transform(a8,alpha8,d8,theta8);
T89=Transform(a9,alpha9,d9,theta9);
T910=Transform(a10,alpha10,d10,theta10);


Tfinal=T01*T12*T23*T34*T45*T56*T67*T78*T89*T910;


%Tfinal dá uma matriz 4x4 onde a três primeiras linhas são uma matriz de
%rotação e um vetor de translação em relação ao frame inicial 

%Temos agora de descobrir a orientação da última joint (end-effector) pelos ângulos de
%Euler 

%Ângulos de Euler estão a ser calculados usando a metodologia Z-Y-Z 

%Há duas formas de calcular os ângulos de Euler. Para isso temos a opção de
%termos Beta positivo ou Beta negativo. Neste caso foram calculados das
%duas formas só por uma opção de escolha 

betapos=atan2(sqrt(Tfinal(1,3)^2 +Tfinal(2,3)^2), Tfinal(3,3));
betaneg=atan2(-sqrt(Tfinal(1,3)^2 +Tfinal(2,3)^2), Tfinal(3,3));
beta = [betapos, betaneg]; 

alph=[0,0];
gama=[0,0];


for i=1:2
    if (sin(beta(i))>10^(-10))
       alph(i)=atan2(Tfinal(2,3), Tfinal(1,3));
       gama(i)=atan2(Tfinal(3,2),-Tfinal(3,1));
    elseif(sin(beta(i))<-10^(-10))
       alph(i)=atan2(-Tfinal(2,3),-Tfinal(1,3));
       gama(i)=atan2(-Tfinal(3,2),Tfinal(3,1));
    elseif((abs(sin(beta(i)))<10^(-10) && cos(beta(i))==1))
       alph(i)=0;
       gama(i)=atan2(Tfinal(2,1),Tfinal(1,1));
    else
       alph(i)=0;
       gama(i)=-atan2(-Tfinal(2,1),-Tfinal(1,1));
    end
end


%DÚVIDA NO FACTO DE HAVER DOIS VALORES DE BETA QUE VAI LEVAR A TER DOIS
%VALORES DE ALPHA E DE GAMA 
    
x = Tfinal(1,4); %posição x
y = Tfinal(2,4); %poisição y
z = Tfinal(3,4); %posição z

oa = round(alph,10); %orientação alpha
ob = round(beta,10); %orientação beta
oc = round(gama,10); %orientação gama 



pos_final=[x, y, z] %posição final a ser mostrada na janela de comandos 
orient_final=[oa, ob, oc] %orientção final a ser mostrada na janela de comandos (mostra as duas 
%opções de ângulos seguidas - alpha1 alpha2 beta1 beta2 gama1 gama2)

end

