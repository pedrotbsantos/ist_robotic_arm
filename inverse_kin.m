% Cinemática Inversa 

% Daniel Fortunato (81498), Teodoro Dias (81723), Pedro Santos (84162)

function [th0, th1, th2, th3, th4, th5] = inverse_kin (x, y, z, oap, obt, ogm)


sol_matrix=zeros(6,8); %inicialização da matriz que mostra todas as soluções possíveis 
%para a posição que é dada nos parâmetros de entrada


% parâmetros das matrizes de transformação que vai ser utilizados para calcular 
% a matriz 
a1=0; a2=0; a3=40; a4=-120; a5=0; a6=0; a7=20; a8=0; a9=0; a10=0;
alpha1=0; alpha2=0; alpha3=-pi/2; alpha4=0; alpha5=0; alpha6=-pi/2; alpha7=0; alpha8=0; alpha9=pi/2; alpha10=-pi/2;
d1=99; d2=0; d3=0; d4=0; d5=0; d6=-40; d7=0; d8=195; d9=0; d10=10; 


%posição (x, y, z)
x6 = x; 
y6 = y; 
z6 = z; 

%orietação (alpha, beta, gama) - é escolhido o primeiro resultado da orientação 
alpha = oap; 
beta = obt; 
gama = ogm; 

%Método Geométrico

%Passar do último frame para o penúltimo a partir dos ângulos de Euler

% -10 é a translação do ponto 6 para o ponto 5
T65 = Transform_ZYZ(0,0,-d10,0,0,0); 

%Transformação do frame de referência para o último frame 
T06=Transform_ZYZ(x6, y6, z6, alpha,beta,gama);

%Transform_ZYZ é uma função que tem a posição e a orientção de um ponto como
%parâmetros de entrada e nos dá a matriz de transformação a partir do 
%método Z-Y-Z



T05 = T06*T65; %Transformação do frame 6 para o frame 5

x5 = T05(1,4); %posição x do ponto 5
y5 = T05(2,4); %posição y do ponto 5
z5 = T05(3,4); %posição z do ponto 5

ro5 = sqrt(y5^2 + x5^2); % Distância do frame de referência ao frame 5 

%O primeiro ângulo(th0) é responsável por colocar o braço na posição (x5,y5, z5)
if y5 == 0 && x5 == 0
    th0 = Inf; %Singularity
else
    th0 = [atan2(x5,-y5), atan2(x5,-y5)+pi]; % Ângulo do frame de referência (explicado no relatório)
    
    %código para meter as soluções de theta0 na matriz de soluções 
    for i=1:4
        sol_matrix(1,i) = th0(1);
        sol_matrix(1,2*i) = th0(2);
        sol_matrix(1,5) = th0(2);
        sol_matrix(1,7) = th0(2);
    end
   
end


syms t1 t2;

%resolução do sistema de equações para calcular theta1 e theta2 (explicado
% no relatório todo o processo criativo por detrás das equações resultantes para o
%cálculo dos ângulos pela forma geométrica 

S1 = solve([99 + 120*sin(t1) - 20*sin(3*(pi/2)-t2-t1) + 155*sin(t1+t2-pi)-z5 == 0, -40 + 120*cos(t1) + 20*cos(3*(pi/2)-t2-t1) + 155*cos(t1+t2-pi)  == -ro5 ],[t1, t2]);
aux1=double(S1.t1);
aux2=double(S1.t2);
       
S2 = solve([99 + 120*sin(t1) - 20*sin(3*(pi/2)-t2-t1) + 155*sin(t1+t2-pi)-z5 == 0, -40 + 120*cos(t1) + 20*cos(3*(pi/2)-t2-t1) + 155*cos(t1+t2-pi) == ro5 ],[t1, t2]);
aux3=double(S2.t1);
aux4=double(S2.t2);


%substituição pela palavra "NaN" (not a number) para todos os theta1 e
%theta2 que são soluções imaginárias (para os casos em que o braço não
%consegue chegar ao ponto)

if imag(aux1(1))~=0 
    aux1(1)=NaN;
end

if imag(aux2(1))~=0 
    aux2(1)=NaN;
end

if imag(aux3(1))~=0 
    aux3(1)=NaN;
end

if imag(aux4(1))~=0 
    aux4(1)=NaN;
end

if imag(aux1(2))~=0 
    aux1(2)=NaN;
end

if imag(aux2(2))~=0 
    aux2(2)=NaN;
end

if imag(aux3(2))~=0 
    aux3(2)=NaN;
end

if imag(aux4(2))~=0 
    aux4(2)=NaN;
end

th1=[aux1,aux3]; %ângulo theta1 
th2=[aux2,aux4]; %ângulo theta2 


% Parâmetros das matrizes de transformação da cinamática direta para os
% frames das joints 0 a 3 
theta1=th0; theta2=-pi/2; theta3=th1; theta4=th2; theta5=-pi/2; theta6=0;  


%ciclos for para percorrer todos os theta0, theta1 e theta2 calculados pelo método
%geométrico
for i=1:2
    for j=1:2
           
        %Matrizes de transformação desde a frame de referência até à frame da joint 3    
        T01=Transform(a1,alpha1,d1,theta1(i));
        T12=Transform(a2,alpha2,d2,theta2);
        T23=Transform(a3,alpha3,d3,theta3(j,i));
        T34=Transform(a4,alpha4,d4,theta4(j,i));
        T45=Transform(a5,alpha5,d5,theta5);
        T56=Transform(a6,alpha6,d6,theta6);

        %Matriz de tranformação entre a frame de referência e a frame da
        %joint 3 (calculada através da cinemática direta e dos ângulos theta0,
        %theta1 e theta2 dados pelo método geométrico
        T03=T01*T12*T23*T34*T45*T56;

        inv_T03=[(T03(1:3,1:3))',-(T03(1:3,1:3))'*T03(1:3,4); 0,0,0,1]; %inversão da matriz 
        %calculada acima

        
        %Método Algébrico 
        
        left_matrix1=inv_T03*T06; %multiplicação da matriz inversa da transformação do frame de 
        %referência para o frame da joint 3 pela matriz de tranformação do
        %frame de referência para o último frame que é dada pela posição e pela orientação 
        %do último frame em relação ao frame de referência 
        
        th4(j,i,:)=[acos(left_matrix1(3,3));-acos(left_matrix1(3,3))]; %theta4 calculado pelo 
        %método algébrico (explicado no relatório o porquê desta igualdade)
        
        
        %ciclo para percorrer todos os th4 
        for k=1:2
            t4_s = round(double(sin(th4(j,i,k))),3);
            
            %se theta4 for zero temos infinitas combinações entre theta3 e
            %theta5, então se isso acontecer forçamos o theta3 a zero 
            if t4_s == 0
                th3(j,i,k) = 0;
                th5(j,i,k)=atan2((-left_matrix1(2,2)),(-left_matrix1(1,2)));
            
            %cálculo de theta3 e theta5 para o caso de th4 não ser zero
            %explicação para as igualidades abaixo no relatório
            else
                th3(j,i,k)=atan2(sign(t4_s)*left_matrix1(1,3),sign(t4_s)*(-left_matrix1(2,3)));
                th5(j,i,k)=atan2(sign(t4_s)*(-left_matrix1(3,2)),sign(t4_s)*left_matrix1(3,1));
            end  
            
        end
    end
end


%código para as soluções de theta1, theta2, theta3, theta4 e theta5 na matriz de soluções(little bit 
%hard coded)

i=1;
sol_matrix(2,i) = aux1(1);
sol_matrix(3,i) = aux2(1);
sol_matrix(2,i+1) = aux1(2);
sol_matrix(3,i+1) = aux2(2);
sol_matrix(2,i+2) = aux1(1);
sol_matrix(3,i+2) = aux2(1);
sol_matrix(2,i+3) = aux1(2);
sol_matrix(3,i+3) = aux2(2);

i=5;
sol_matrix(2,i) = aux3(1);
sol_matrix(3,i) = aux4(1);
sol_matrix(2,i+1) = aux3(2);
sol_matrix(3,i+1) = aux4(2);
sol_matrix(2,i+2) = aux3(1);
sol_matrix(3,i+2) = aux4(1);
sol_matrix(2,i+3) = aux3(2);
sol_matrix(3,i+3) = aux4(2);

%th3

sol_matrix(4,1) = th3(1,1,1);
sol_matrix(4,3) = th3(1,1,2);
sol_matrix(4,2) = th3(2,1,1);
sol_matrix(4,4) = th3(2,1,2);

sol_matrix(4,5) = th3(1,2,1);
sol_matrix(4,7) = th3(1,2,2);
sol_matrix(4,6) = th3(2,2,1);
sol_matrix(4,8) = th3(2,2,2);

%th4

sol_matrix(5,1) = th4(1,1,1);
sol_matrix(5,3) = th4(1,1,2);
sol_matrix(5,2) = th4(2,1,1);
sol_matrix(5,4) = th4(2,1,2);

sol_matrix(5,5) = th4(1,2,1);
sol_matrix(5,7) = th4(1,2,2);
sol_matrix(5,6) = th4(2,2,1);
sol_matrix(5,8) = th4(2,2,2);

%th5

sol_matrix(6,1) = th5(1,1,1);
sol_matrix(6,3) = th5(1,1,2);
sol_matrix(6,2) = th5(2,1,1);
sol_matrix(6,4) = th5(2,1,2);

sol_matrix(6,5) = th5(1,2,1);
sol_matrix(6,7) = th5(1,2,2);
sol_matrix(6,6) = th5(2,2,1);
sol_matrix(6,8) = th5(2,2,2);


sol_matrix



% Testar novamente a direta para todas as soluções da matriz  

    for k=1:8


        tht0=sol_matrix(1,k);
        tht1=sol_matrix(2,k);
        tht2=sol_matrix(3,k);
        tht3=sol_matrix(4,k);
        tht4=sol_matrix(5,k);
        tht5=sol_matrix(6,k);

        direct_kin(tht0, tht1, tht2, tht3, tht4, tht5);

    end



end


