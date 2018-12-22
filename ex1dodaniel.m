%% Exemplo numérico para o Teorema da verificação
clear all
close all
clc
%onde, na formulação da dinâmica do sistema: a = alfa e b = beta
%atribuindo valores
d=1; f=2; q=1; x0=3; a=5; b=7;
%Cálculo das matrizes S(j) da equação de Riccati
S2 = q;
S1 = a*S2*a + d*d - (a*S2*b + d*f)*((b*S2*b + f*f).^(-1))*(b*S2*a+f*d);
S0 = a*S1*a + d*d - (a*S1*b + d*f)*((b*S1*b + f*f).^(-1))*(b*S1*a+f*d);
%O controle ótimo é dado por uj(x) = - M(j)*x
M1 = ((b*S2*b +f*f).^(-1))*(b*S2*a + f*d);
M0 = ((b*S1*b +f*f).^(-1))*(b*S1*a + f*d);
%portanto o controle ótimo é dado por:
u0 = -M0*x0;
x1 = a*x0 + b*u0;
u1 = -M1*x1;
x2 = a*x1 + b*u1;
%A fim de verificar a otimalidade do custo, temos que mostrar que
%V0(x0)=J_2(u*), Por Bellman temos:
V2 = x2*S2*x2;
V1 = x1*S1*x1;
V0 = x0*S0*x0;%V0(x0)
%O custo ótimo é cálculado por J_2(u*)
J = ((d*x0+f*u0)^2) + ((d*x1+f*u1)^2) + V2;%J_2(u*)