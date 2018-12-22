% Exercicio 4 - LQG
clc, close all, clear all

%% Sistema Escolhido n�o-param�trico
Ts = 1; % Tempo a cada amostra
Gu = tf(540, conv([1 6],[1 5 72]));     % F.T. continua
Gud = c2d(Gu,Ts);                          % F.T. discreta
sys = ss(Gud);       % Espa�o de Estados 

Co = ctrb(sys.A,sys.B) 
rank(Co) % testa controlabilidade
Ob = obsv(sys.A,sys.C)
rank(Ob) % testa observabilidade

%% LQ Deterministico (wk=0) (x_k+1 = sys.a * x_k + sys.b * u_k)

% i)(Finito com N=3)

N = 3;  
Q = eye(3); %Q = I

% c�lculo dos ganhos de controle �timo M(j) e S(N)
S{N} = Q;
A = sys.a; B=sys.b; F = 1000; D = 1;

for j=N-1:-1:1
    M{j} = ((B'*S{j+1}*A + F'*D)/(B'*S{j+1}*B + F'*F));
    S{j} = A'*S{j+1}*A + D'*D - (A'*S{j+1}*B + D'*F)*(B'*S{j+1}*A + F'*D);
end

% c�lculo do filtro de Kalman
F = sys.a; G = sys.b; H = sys.C; D =  sys.D;
P{1} = F*Q*F'; % covariancia inicial 
L=[0 0 1; 0 1 1; 1 1 1];
R=1;
x{1} = sqrt(5)*L*randn(3,1); % gera condi��o inicial
y{1} = H*x{1}; % gera observa��o inicial
xa = [1 2 3]'; % valor do estado a priori [k|k-1]
xe{1} = xa;

for k=1:N-1
    % Measurement update
    K{k} = P{k}*H'/(H*P{k}*H'+R); % Ganhos de kalman
    P{k+1} = (eye(3) - K{k}*H)*P{k}; % Proxima Cov 
    
    % evolu��o do sistema verdadeiro:
    u{k} = -M{k}*xa; % controle �timo que minimiza o custo
    x{k+1} = F*x{k} + G*u{k}; 
    y{k+1} = H*x{k+1};

    %valor estimado:
    xa = F*xe{k} + G*u{k};    % x[k+1|k]
    xe{k+1}= xa + K{k}*(y{k+1}-H*xa);  % x[k+1|k+1]

    % custo no est�gio
    c(k) = x{k}'*Q*x{k} + u{k}'*R*u{k};
end
% custo final
c(N) = x{N}'*S{N}*x{N};

%% ii)(Infinito)

A = sys.a; B=sys.b; D = 1; F =  1000;
% teste de estabilizavel e detectavel
%(A,B) � estabilizavel pois (A,B) � control�vel
D_ = (eye(3)-F*inv(F'*F)*F')*D;
A_ = A-(B*inv(F'*F)*F'*D);
a=A(1,:);
Ob2 = obsv(D_,a)
rank(Ob2) % testa observabilidade para verificar detecbilidade
eig(A_) % verifica o p�lo n�o observ�vel para analise
aa = [-400*eye(3) - a; D_]% verifica a detectabilidade no 
rank(aa)                  % no p�lo fora do circulo unitario
[X,L,G] = dare(A,B,Q,R) % Valor do Ricatti
Mx=inv(B'*X*B+F'*F)*(B'*X*A+F'*D); 
x_inf = [1 2 3]';
u_inf=-Mx*x_inf %controle �timo infinito
V_inf=x_inf'*X*x_inf %custo no infinito

