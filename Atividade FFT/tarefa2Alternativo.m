clc
clear all
close all


%ETAPA 1 - CRIAÇÃO DO SINAL ANÁLOGICO E DEFINIÇÃO DO NÚMERO DE AMOSTRAS E FREQUÊNCIA DE AMOSTRAGEM
na = 1 : 1 : 100; %100 AMOSTRAS em nosso sinal analógico

m  = 16;     %Quantidade de amostras por ciclo definida pelo operador
fa = m * 60; %Frequência de amostragem
ta = 1/fa;   %Período de amostragem
A1  = (220/sqrt(3))*sqrt(2);  %Amplitude da 1 Harmonica [V]
Fase1 = deg2rad(45); %Fase da 1 Harmonica
A3  = A1/3; %Amplitude da 3 Harmonica [V]
Fase3 = deg2rad(15); %Fase da 3 Harmonica
t = na*ta; %Vetor de tempo para cada Amostra

SinalAn = A1 * cos (2*pi*60*t + Fase1) + A3 * cos (2*pi*3*60*t + Fase3); %vetor formado pela sinal senoidal analógico
SinalAn = SinalAn + 2.5; %Soma Metade do fundo de escala para conversor unipolar (0 a 5 V)

q = 5/(2^12); %Bits por volts
SinalDig = round(SinalAn/q); %converte a tensao na representacao do AD unipolar

%Plotagem do Sinal Analógico e Sinal Digital, utilizando na como o eixo X
figure(1); subplot(211); plot(na, SinalAn);    grid; ylabel('Tensao AD');
           subplot(212); stairs(na, SinalDig); grid; ylabel('Amostras AD');

%ETAPA 2 - Transformada Discreta de Fourier

for i=1:1:m %Coeficientes 1 Harmonica
  Coef_b_cos1(i) = (sqrt(2)/m)*cos((2*pi*(i-1)*1)/m);
  Coef_b_sin1(i) = (sqrt(2)/m)*sin((2*pi*(i-1)*1)/m);
end

for i=1:1:m %Coeficientes 3 Harmonica
  Coef_b_cos3(i) = (sqrt(2)/m)*cos((2*pi*(i-1)*3)/m);
  Coef_b_sin3(i) = (sqrt(2)/m)*sin((2*pi*(i-1)*3)/m);
end

%Código de plotagem dos gráficos retirado do exemplo do professor
figure(2); 
subplot(221); stem(Coef_b_cos1); grid; ylabel('Coefs. cosseno'); 
title('Coefs. da 1a. harm.');
subplot(222); stem(Coef_b_cos3); grid;
title('Coefs. da 3a. harm.');
subplot(223); stem(Coef_b_sin1); grid; ylabel('Coefs. seno');
subplot(224); stem(Coef_b_sin3); grid;

x = zeros(1,m);
for z = 1:1:100
    %Zerando variaveis para evitar problemas na recursão
    y_r1 = 0;
    y_im1 = 0;
    y_r3 = 0;
    y_im3 = 0;
    for j = 2:1:m
      x(m-(j-2)) = x(m-(j-1));
      %janela deslizante
    end
    x(1) = SinalDig(z);
    for k = 1:1:m
      %Cálculo de Y para 1 Harmonica
      y_r1 =  x(k)*Coef_b_cos1(k)+y_r1; 
      y_im1 = x(k)*Coef_b_sin1(k)+y_im1;
      %Cálculo de Y para 3 Harmonica
      y_r3 =  x(k)*Coef_b_cos3(k)+y_r3;
      y_im3 = x(k)*Coef_b_sin3(k)+y_im3;
    end
    yy_complexo1 = complex(y_r1,y_im1);
    yy1(z) = abs(yy_complexo1);
    y_fase1(z) = angle(yy_complexo1);
    yy_complexo3 = complex(y_r3,y_im3);
    yy3(z) = abs(yy_complexo3);
    y_fase3(z) = angle(yy_complexo3);
    
end

% ETAPA 3 - APRESENTAÇÃO DOS RESULTADOS 
for i = 1:length(yy1)-1
    % Verificando se o elemento i é igual ao próximo elemento i+1
    if yy1(i) == yy1(i+1)
        tempoz = i; %amostra em que o regime permanente é estabelecido
        Tempoo = tempoz*ta;
        disp("Regime Permanente estabelecido na amostra"); disp(i);
        disp("Ou seja o tempo de atraso digital é (em segundos):"); disp(Tempoo);
        break;
    end
end


%Para uma melhor padronização dos resultados, utilizamos o layout dos gráficos que o professor definiu em seu modelo.
figure(3); subplot(211);
plot(na, SinalDig-2048, na, yy1); grid; ylabel('Unidades do AD');
legend('Sinal digital sem offset', 'Modulo do fasor (eficaz)');
subplot(212);
plot(na, (SinalDig-2048)*q, na, yy1*q); grid; ylabel('Volts');
legend('Sinal original', 'Modulo do fasor (eficaz)');

figure(4);
plot(na, yy1*q*sqrt(2), na, yy3*q*sqrt(2)); grid; 
title('Modulos da 1a e 3a hamonica');
ylabel('Volts de pico');

disp('Módulo do fasor da 1a. harm.  na amostra 60');  disp(yy1(60)*q*sqrt(2));
disp('Módulo do fasor da 3a. harm.   na amostra 60');  disp(yy3(60)*q*sqrt(2));