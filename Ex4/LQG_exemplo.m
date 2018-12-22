clear all
sys=drss(5) % produz um sistema a tempo discreto e estável de 5a. ordem
lambda = max(abs(eig(sys.a)))
sys.a =sys.a*1.5/lambda  % transforma o sistema escolhido num sistem instável
sys.c = [1 0 0 0 0]; % observação somente da 1a. variável de estado

Co = ctrb(sys.A,sys.B); rank(Co) % testa controlabilidade
Ob = obsv(sys.A,sys.C); rank(Ob) % testa observabilidade

% ruido de medida y_k = C*x_k + sqrt(var_y)* e_k  com e_k ~ N(0,1)
var_y = 2;

% ruido do sistema x_k+1 = sys.a * x_k + sys.b * u_k + F*w_k com dim(w_k)=3
F = [1 0 0
    0 0 0
    0 1 0
    0 0 0
    0 0 1];
Cor = ctrb(sys.A,F); rank(Cor) % testa controlabilidade através do ruído

% desempenho E[ sum_1^(N-1) [x_k'*Qc*x_k +u_k*Rc*u_k] + x_N'*Sf*x_N]
% problema inicia-se no instante k=1 com x_1 ~ N( m_1, P_1)

Qc=diag([1 2 4 4 2]);
Obc = obsv(sys.A,sqrt(Qc)); rank(Obc) % testa observabilidade através do custo

display('Pausa para verificação'), keyboard

N=20;
Sf=1e+3*diag([5 5 2 2 1]);
Rc=1e-3;
%c. i.
L = [1 0 0 0 0; 1 1 0 0 0; 1 1 1 0 0; 1 1 1 1 0; 1 1 1 1 1];
P_1 = 500*L*L';
m_1 = [0 0 0 0 0]';

A = sys.a; B=sys.b;
% cálculo dos ganhos de controle ótimo M(j)
S{N} = Sf;
for j=N-1:-1:1
    Mn = A'*S{j+1}*B/(B'*S{j+1}*B+Rc);
    M{j} = Mn';
    S{j} = A'*S{j+1}*A + Qc - M{j}'*(B'*S{j+1}*B+Rc)*M{j};
end

C = sys.c; R = var_y; Q = eye(3);
% cálculo do ganho do filtro de Kalman
P{1} = P_1;
for k=1:N-1
    % Measurement update
    K{k} = A*P{k}*C'/(C*P{k}*C'+R);
    P{k+1} = A*P{k}*A'+ F*Q*F' - K{k}*(C*P{k}*C'+R)*K{k}';
end

%%%%%%%%%%%%%%%%%%%%%% Simulação %%%%%%%%%%%%%%%%%%%
for m=1:8 % n.o de simulações
    x{1} = sqrt(500)*L*randn(5,1); % gera condição inicial
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
        c(k) = x{k}'*Qc*x{k} + u{k}'*Rc*u{k};
    end
    % custo final
    c(N) = x{N}'*Sf*x{N};
    custo(m,:)= cumsum(c);
end
 plot(custo')