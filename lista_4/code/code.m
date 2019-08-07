clc
clear all
close all

%%
%ALUNO: GABRIEL DA SILVA CHAVES
%EXERCÃ?CIO: 
linewidth = 5;
%%
%ESPECIFICAÃÃO DO FILTRO DIGITAL

A_p1 = 1;   %Atenuacao da banda passante 1 (dB)
A_p2 = 2;   %Atenuacao da banda passante 2 (dB)
A_r  = 40;  %Atenuacao da banda de rejeicao (dB)

W_p1 = 850;  %Frequencia passante 1 (rad/s)
W_r1 = 980;  %Frequencia de rejeicao 1 (rad/s)
W_r2 = 1020; %Frequencia de rejeicao 2 (rad/s)
W_p2 = 1150; %Frequencia passante 2 (rad/s)
W_s  = 1e4;  %Frequencia de amostragem (rad/s)

%%
%CONVERSÃO PARA FREQUÃNCIAS DIGITAIS

w_r1 = 2*pi*W_r1/W_s; %Frequencia de rejeicao 1 (rad/Sa)
w_p1 = 2*pi*W_p1/W_s; %Frequencia passante 1 (rad/Sa)
w_p2 = 2*pi*W_p2/W_s; %Frequencia passante 2 (rad/Sa)
w_r2 = 2*pi*W_r2/W_s; %Frequencia de rejeicao 2 (rad/Sa)

%%
%CONVERSÃO PARA FREQUÃNCIAS ANALÃGICAS

T = 2*pi/W_s; %Periodo de amostragem (s)

w_ar1 = 2/T*tan(w_r1/2); %Frequencia de rejeicao 1 (rad/s)
w_ap1 = 2/T*tan(w_p1/2); %Frequencia passante 1 (rad/s)
w_ap2 = 2/T*tan(w_p2/2); %Frequencia passante 2 (rad/s)
w_ar2 = 2/T*tan(w_r2/2); %Frequencia de rejeicao 2 (rad/s)

%%
%SIMETRIA

if     w_ar1*w_ar2 < w_ap1*w_ap2
    w_ap2 = w_ar1*w_ar2/w_ap1;   %Diminui a frequencia passante
elseif w_ar1*w_ar2 > w_ap1*w_ap2
    w_ap1 = w_ar1*w_ar2/w_ap2;   %Aumenta a frequencia passante
end

%%
%CHEBYSHEV

W_0  = sqrt(w_ar1*w_ar2);                     %Frequencia central, geometricamente
a    = 1;                                     %Fator de normalizacao
Bw   = w_ap2 - w_ap1;                         %Banda passante do filtro
W_pn = 1/a;                                   %W_p normalizado
W_rn = (1/a)*(w_ap2 - w_ap1)/(w_ar2 - w_ar1); %W_r normalizado

A_p = min(A_p1,A_p2); %Maior atenuacao das bandas;

e = sqrt(10^(0.1*A_p) - 1);                                  %Calculando o epsilon
n = ceil(acosh(sqrt((10^(0.1*A_r) - 1)/(e^2)))/acosh(W_rn)); %Menor ordem do filtro

i     = 0 : 2*n - 1;

x_1 = (2*i + 1)*pi/(2*n);
x_2 = (1/n)*asinh(1/e);

s_n   = sin(x_1)*sinh(x_2) + 1j*cos(x_1)*cosh(x_2); %Todos os zeros do polinÃŽmio A(s')A(-s')
s_aux = find(real(s_n)<0);                          %Indices dos zeros com parte real negativa
p   = s_n(s_aux);                                   %Zeros com parte real negativa

if mod(n,2) == 0
    H_0 = 10^(-0.05*A_p)*prod(-p); %Ganho de H(s') para n par
else
    H_0 = prod(-p);                %Ganho de H(s') para n impar
end

[Bn, An] = zp2tf([],p,H_0); %Coeficientes de H(s')

[Bd, Ad] = lp2bs(Bn,An,W_0,Bw); %Desnormalizacao de filtro passa-baixas para rejeita-faixa

N = 2^13;

[H_an,W1] = freqs(Bn,An,N); %H(s')
[H_ad,W2] = freqs(Bd,Ad,N); %H(s)

[B,A]    = bilinear(Bd,Ad,1/T); %Coeficientes de H(z)
[H_d,W3] = freqz(B,A,N);        %H(z)

figure('name','Resposta em magnitude - Filtro passa-baixas analogico.','position',[0 0 800 800])
plot(W1,20*log10(abs(H_an)),'LineWidth', linewidth)
xlabel('Frequência (rad/s)','FontSize',20)
ylabel('Magnitude (dB)','FontSize',20)
axis([min(W1) max(W1) -50 max(abs(H_an))+2])
set(gca,'Fontsize',15)

figure('name','Resposta em magnitude - Filtro rejeita-faixa analÃ³gico.','position',[0 0 800 800])
plot(W2,20*log10(abs(H_ad)),'LineWidth', linewidth)
xlabel('Frequência (rad/s)','FontSize',20)
ylabel('Magnitude (dB)','FontSize',20)
axis([w_ap1-50 w_ap2+50 -50 max(abs(H_ad))+2])
set(gca,'Fontsize',15)

figure('name','Resposta em magnitude - Filtro rejeita-faixa digital.','position',[0 0 800 800])
hold on
plot(W3*W_s/(2*pi),20*log10(abs(H_d)),'LineWidth', linewidth)
line([850, 850],[-100, 100],'LineStyle','--','LineWidth',2.0,'Color','r')
line([1150, 1150],[-100, 100],'LineStyle','--','LineWidth',2.0,'Color','r')
line([980, 980],[-100, 100],'LineStyle','--','LineWidth',2.0,'Color','r')
line([1020, 1020],[-100, 100],'LineStyle','--','LineWidth',2.0,'Color','r')
xlabel('Frequência (Hz)','FontSize',20)
ylabel('Magnitude (dB)','FontSize',20)
axis([W_p1-50 W_p2+50 -50 max(abs(H_d))+2])
set(gca,'Fontsize',15)