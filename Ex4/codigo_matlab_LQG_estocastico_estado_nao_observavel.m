%EXERCÍCIO COMPUTACIONAL - LQG estocástico com estados não
%observáveis
clc, clear all, close all

Gu=tf([450], conv([1 5],[1 7 90])); % Sistema Continuo
Ge = tf([3 2000],([1 150 2000])); % Ruído é filtrado por
Ts=1 ;% periodo de amostragem
Gud=c2d(Gu,Ts); % converte sistema contínuo --> discreto (c2d = continuos to discrete)
Ged=c2d(Ge,Ts);
ss(Gud);%converte o sistema para o espaço de estados
ss(Ged);%converte o sistema para o espaço de estados

A = [1.981  -0.7561  0.4317; 2  0  0; 0  0.5  0];
B = [0.25; 0; 0];
C = [0 0.1306  0.05385];
D = sqrt(0.04);
%covariância do ruído do processo
Q =eye(3);
%covariância do ruído do sensor
R=eye(3)*0.04;

%para o controle
f=1;
F=f;
S = Q;%Matriz S inicial
xk=[30;22;16];%estado inicial

%para o filtro
P = A*Q*A';%covariância inicial do erro

norma = [];%armazena a norma do vetor de estados
norma = [norma; norm(xk)];
xk1=[];xk2=[];xk3=[];%armazena individualmente cada elemento do vetor de estados 
xk1=[xk1; xk(1,1)];xk2=[xk2; xk(2,1)];xk3=[xk3; xk(3,1)];
custo=[];%custo a cada iteração
custo_estage=0;%normrnd(0,0.04)
custo = [custo; custo_estage];
errocov=[];

for i=1:299   
S = A'*S*A + D'*D - (A'*S*B + D'*F)*inv(B'*S*B + F'*F)*(B'*S*A + F'*D);
M = inv(B'*S*B + F'*F)*(B'*S*A + F'*D)
noise_measure = [normrnd(0,0.04);normrnd(0,0.04);normrnd(0,0.04)];%ruído na medição do sensor

%ganho de kalman
H=eye(3);
KG = P*H'/(H*P*H'+R);%ganho de kalman P/(P + R)
%observações com ruído
yk = C*xk + D*noise_measure;
P = (eye(3) - KG*H)*P;
%covariância do erro de estimação do estado
errocov=[errocov; H*P*H'];
%controle, portanto calculamos o estado seguinte por:
uk = -M*xk;
%estado atual
xk = A*xk + B*uk + KG*(yk - C*xk);
P = A*P*A' + Q -(H*P*H'+R);
xk1=[xk1; xk(1,1)];xk2=[xk2; xk(2,1)];xk3=[xk3; xk(3,1)];
norma = [norma;norm(xk) ];

custo_estage = custo_estage + ((norm(yk + f*uk)).^2)+ trace(C*P*C');%soma das duas primeiras parcelas do custo
custo =[custo; custo_estage];
end

custo_estage = custo_estage + trace(P*Q);%Terceira parcela do custo
%custo médio
disp('Custo médio')
custo_estage/i

%Cálculo dos auto valores de (A-BM)
disp('Auto-valores da matriz quadrada A-BM')
eig(A-B*M)

figure(1)
plot(norma)
grid on
title('Norma do vetor de estados Xk ao longo do tempo K')
ylabel('Norma ||Xk||')
xlabel('tempo K (s)')

figure(2)
plot(xk1,'k')
grid on
title('Elemento 1 do vetor de estados')
xlabel('tempo K (s)')

figure(3)
plot(xk2,'r')
grid on
title('Elemento 2 do vetor de estados')
xlabel('tempo K (s)')

figure(4)
plot(xk3,'g')
grid on
title('Elemento 3 do vetor de estados')
xlabel('tempo K (s)')

figure(5)
plot(custo,'b')
grid on
title('Evolução do custo')
xlabel('tempo K (s)')

 figure(6)
 plot3(xk1,xk2,xk3), grid on
 xlabel('Xk1')
 ylabel('Xk2')
 zlabel('Xk3')
