% Exercicio 4 - LQG
clc, close all, clear all

%% Sistema Escolhido não-paramétrico
Ts = 1; % Tempo a cada amostra
Gu = tf(540, conv([1 6],[1 5 72]));     % F.T. continua
Gud = c2d(Gu,Ts);                          % F.T. discreta
sys = ss(Gud);       % Espaço de Estados 

Co = ctrb(sys.A,sys.B) 
rank(Co) % testa controlabilidade
Ob = obsv(sys.A,sys.C)
rank(Ob) % testa observabilidade

%% Caso Ruido I)

% ruido de medida y_k = C*x_k + sqrt(var_y)* e_k  com e_k ~ N(0,1)
var_y = 0.04;

% ruido do sistema x_k+1 = sys.a * x_k + sys.b * u_k + F*w_k 
F = [1 0 0;0 1 0;0 0 1];
Cor = ctrb(sys.A,F); rank(Cor) % testa controlabilidade através do ruído

% desempenho E[ sum_1^(N-1) [x_k'*Qc*x_k +u_k*Rc*u_k] + x_N'*Sf*x_N]
% problema inicia-se no instante k=1 com x_1 ~ N( m_1, P_1)
L = [1 0 0; 1 1 0;1 1 1];
P_1 = 500*L*L';
m_1 = [0 0 0]';

% i)

N = 3;  
R = var_y;
Q = eye(3); %Q = I

% cálculo dos ganhos de controle ótimo M(j) e S(N)
S{N} = Q;
A = sys.a; B=sys.b; D = 1; F=1; C=sys.c;

for j=N-1:-1:1
    M{j} = ((B'*S{j+1}*A + F'*D)/(B'*S{j+1}*B + F'*F));
    S{j} = A'*S{j+1}*A + D'*D - (A'*S{j+1}*B + D'*F)*(B'*S{j+1}*A + F'*D);
end

% cálculo do ganho do filtro de Kalman
P{1} = P_1;
for k=1:N-1
    % Measurement update
    K{k} = A*P{k}*C'/(C*P{k}*C'+R);
    P{k+1} = A*P{k}*A'+ F*Q*F' - K{k}*(C*P{k}*C'+R)*K{k}';
end

% Controle
x{1} = sqrt(500)*L*randn(3,1); % gera condição inicial
y{1} = C*x{1} + sqrt(var_y)*randn; % gera observação inicial
xa = m_1; % valor do estado a priori [k|k-1]
xe{1} = xa;
     
for k=1:N-1
    % evolução do sistema verdadeiro:
    u{k} = -M{k}*xa;
    x{k+1} = A*x{k} + B*u{k} + F* randn(3,1);
    y{k+1} = C*x{k+1} + sqrt(var_y)*randn;

    %valor estimado:
    xa = A*xe{k} + B*u{k};    % x[k+1|k]
    xe{k+1}= xa + K{k}*(y{k+1}-C*xa);  % x[k+1|k+1]

    % custo no estágio
    c(k) = x{k}'*Q*x{k} + u{k}'*R*u{k};
end
% custo final
c(N) = x{N}'*S{N}*x{N};
custo = cumsum(c);
plot(c,'*')

%% ii)(Infinito)

A = sys.a; B=sys.b; D = 1; F = 1;
% teste de estabilizavel e detectavel
%(A,B) é estabilizavel pois (A,B) é controlável
D_ = (eye(3)-F*inv(F'*F)*F')*D;
A_ = A-(B*inv(F'*F)*F'*D);
a=A_(1,:);
Ob2 = obsv(D_,a)
rank(Ob2) % testa observabilidade para verificar detecbilidade
eig(A_) % verifica o pólo não observável para analise
aa = [-3.9*eye(3) - a; D_]% verifica a detectabilidade no 
rank(aa)                  % no pólo fora do circulo unitario
[X,L,G] = dare(A,B,Q,R) % Valor do Ricatti
Mx=inv(B'*X*B+F'*F)*(B'*X*A+F'*D); 
x_inf = [1 2 3]';
u_inf = -Mx*x_inf %controle ótimo infinito
V_inf = x_inf'*X*x_inf %custo no infinito

%% descontado(Descontado)

p = 0.0001;
p2 = p^(1/2);

A = p2*A; B = p2*B; F = p2*F; D= p2*D;
[Xp,L,G] = dare(A,B,Q,R) % Valor do Ricatti
Mxp=inv(B'*Xp*B+F'*F)*(B'*Xp*A+F'*D); 
x_infp = [1 2 3]';
u_infp = -Mxp*x_infp %controle ótimo infinito
V_infp = x_infp'*Xp*x_infp %custo no infinito
