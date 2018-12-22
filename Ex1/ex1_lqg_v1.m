%% IA 856 Controle Linear Quadratico Gaussiano Ex.1
clear all,close all,clc

% Criacao das variaveis:

syms d, syms f, syms q, syms x0, syms x1, syms x2
syms a  % alpha
syms b  % beta

% gerando S(j):

s2 = q; % S(2)
A  = (a*s2*b + d*f); % A,B,C apenas para organizar S(1)
B  = (b*s2*b + f*f)^(-1);
C  = (b*s2*a+f*d);
s1 = a*s2*a + d*d -A*B*C; %S(1)
D  = (a*s1*b + d*f); % A,B,C apenas para organizar S(1)
E  = (b*s1*b + f*f)^(-1);
F  = (b*s1*a+f*d);
s0 = a*s1*a + d*d -D*E*F; %S(0)

% gerando M(j):

M1 = ((b*s2*b +f*f)^(-1))*(b*s2*a + d*d);
M0 = ((b*s1*b +f*f)^(-1))*(b*s1*a + d*d);

% gerando u(j):

u1 = -M1*x1
u0 = -M0*x0

% gerando Vj(x)

v2 = x2*s2*x2
v1 = x1*s1*x1
v0 = x0*s0*x0

%calculo do J

J=((d*x0+f*u0)^2)+((d*x1+f*u1)^2) + v2

