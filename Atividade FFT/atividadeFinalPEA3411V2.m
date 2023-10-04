clc
clear
close all

%ETAPA 1 - CRIAÇÃO DO SINAL ANÁLOGICO E DEFINIÇÃO DO NÚMERO DE AMOSTRAS E FREQUÊNCIA DE AMOSTRAGEM:

f0 = 60; % Frequência fundamental do sinal (60 Hz)
tfinal = 0.1; % Duração do sinal (0.1 s)
m = 2000; % Número de amostras por ciclo da fundamental do sinal
fa = m*f0; % Frequência de amostragem para se criar o sinal analógico
t = 0:1/fa:tfinal; % Vetor de tempos para produzir as amostras desses sinais de tempo continuo
V1 = (220/sqrt(3))*sqrt(2); % Amplitude da fundamental do sinal de valor 220/sqrt(3) VRMS
Fase1 = deg2rad(45); %Fase da fundamental (45°)
V3 = V1/3; % Amplitude da terceira harmônica do sinal
Fase3 = deg2rad(15); %Fase da 3 Harmonica (15°)
fundamental = V1*cos(2*pi*f0*t+Fase1);
terceira = V3*cos(2*pi*3*f0*t+Fase3);
sinal = fundamental + terceira;

% Simula a entrada desses dois sinais secundarios em entradas analogicas de um IED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1)Condicionamento: diminuir a intensidade dos sinais que vem do TP 

sinalIn = sinal * (6/220);     %TP interno de 220:6V
% Ajustando a tensão com um divisor resistivo de 1/5
sinalIn = sinalIn * (1 / 5);

% Adicionando um offset aos canais para eles possuirem tensao analogica UNIPOLAR
% Admitindo eletronica com tensoes entre 0 a 5V, usaremos um offset de 2,5 V.
sinalIn = sinalIn + 2.5;

% 2) Protecao dos sinais condicionados, para limitar os seus valores entre as necessidades
%   dos filtros e da eletronica. Supondo uma eletronica para sinais entre 0 e 5 V

% Simula o circuito de grampeamento
sinalInLim = min(sinalIn, 5.5);     % Valor maximo permitido eh 5 + 0.5 V
sinalInLim = max(sinalInLim,-0.5);  % Valor minimo permitido eh 0.0 - 0.5 V

% 3) Filtragem: simula a filtragem passa baixa para fazer o anti-aliasing
%   Considerando uma frequencia de amostragem (que sera feita
%   posteriormente) com fa = 16*60 = 960 Hz (ou seja, amostragem de 16
%   amostras por ciclo)
%   Considerando um  conversor AD (que sera simulado posteriormente) com 12 bits n = 12

% Especificações do filtro
frequencia_banda_passagem = 60; % Frequência limite da banda de passagem em Hz (deseja-se atenuar todos as harmônicas fora a fundamental)
Amax = 3;            % Atenuação maxima admissivel da banda de passagem em dB (valor escolhido em 
                                 % 3 dB por conta da definição de frequência de corte
frequencia_banda_rejeicao = 180; % Frequência limite da banda de rejeição em Hz
Amin = 40;           % Atenuação minima admissivel da banda de rejeicao em dB (40 dB foi escolhido por ser uma atenuação que faz com o sinal deixe de ser relevante)
                                 % Esse valor não precisa ser maior que 1/2^n

% Calculo das frequencias angulares
Wp = 2*pi*frequencia_banda_passagem;   % Frequencia da banda de passagem em rad/s
Ws = 2*pi*frequencia_banda_rejeicao;   % Frequencia da banda de rejeicao em rad/s

% Projeto do filtro butterworth com funcao de transferencia e componentes ideais
% Ordem do filtro segundo as especificacoes
[nFiltro, Wn] = buttord(Wp, Ws, Amax, Amin, 's');

% Projeto dos polinomios da funcao de transferencia do filtro
[num, den] = butter(nFiltro, Wn, 'low', 's');

% Montagem da funcao de transferencia
filtroPB = tf(num,den);

%Simulando o funcionamento dos filtros anti-aliasing passa baixa no sinal 
sinalFil = lsim(filtroPB, sinalInLim, t);

figure(1);
plot(t, sinal);
grid
title('Sinal original')

figure(2);
plot(t,sinalInLim, t, sinalFil);
grid
legend('Sinal condicionado e limitado', 'Sinal filtrado')
title('Condicionamento, protecao e filtragem do sinal')


% ETAPA 2 - Transformada Discreta de Fourier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Etapa de simulacao da amostragem e digitalização
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%A) Simular o processo de amostragem (sample & hold)
m = 16; % Deseja-se ter 16 amostras por ciclo da onda fundamental
faFiltro = m*60; %freq de amostragem desejada [Hz]
taFiltro = 1/faFiltro;
% Fator de decimação - nesse caso que estamos emulando o comportamento do tempo
% continuo, vamos pegar amostras igualmente espacadas nos vetores dos sinais,
% a cada 'decimation_factor' numero de amostras
decimation_factor = round(fa/faFiltro);
% Inicialização do vetor de dados amostrados - prepara o processo de amostragem
Tsignal_smp = [];   %vetor para armazenar as amostras do tempo (auxiliar)
signal_smp = [];   %vetor para armazenar as amostras da tensao
% Loop para subamostragem ou decimacao
for i = 1:decimation_factor:length(sinalFil)    %varre o vetor de tempo continuo, saltando de decimation_factor
    signal_smp = [signal_smp, sinalFil(i)];   %pega uma amostra de tensao
    Tsignal_smp = [Tsignal_smp, t(i)];        %tambem amostro para saber o tempo de cada amostra
end

% Resultados
% Observacao do sinal no tempo contino e sinal amostrado por circuito de sample & hold
figure(3)
stairs(Tsignal_smp, signal_smp); hold on; grid;
stairs(t, sinalFil);
legend('Sinal amostrado', 'sinal continuo');
title('Processo de amostragem do S&H')

%B) Simulacao da digitalizacao com conversor AD
n=12;        %12 bits de resolucao - 2^12 simbolos diferentes
q=5/(2^n); % Quanta do conversor em Volts / simbolo

signal_dig = round(signal_smp/q);

figure(4)
stem(Tsignal_smp, signal_dig); hold on; grid;
title('Amostras digitalizadas com AD');

% Transformada discreta de Fourier

for i=1:1:m %Coeficientes 1 Harmonica
  Coef_b_cos1(i) = (sqrt(2)/m)*cos((2*pi*(i-1)*1)/m);
  Coef_b_sin1(i) = (sqrt(2)/m)*sin((2*pi*(i-1)*1)/m);
end

for i=1:1:m %Coeficientes 3 Harmonica
  Coef_b_cos3(i) = (sqrt(2)/m)*cos((2*pi*(i-1)*3)/m);
  Coef_b_sin3(i) = (sqrt(2)/m)*sin((2*pi*(i-1)*3)/m);
end

for i=1:1:m %Coeficientes 5 Harmonica
  Coef_b_cos5(i) = (sqrt(2)/m)*cos((2*pi*(i-1)*5)/m);
  Coef_b_sin5(i) = (sqrt(2)/m)*sin((2*pi*(i-1)*5)/m);
end

for i=1:1:m %Coeficientes 7 Harmonica
  Coef_b_cos7(i) = (sqrt(2)/m)*cos((2*pi*(i-1)*7)/m);
  Coef_b_sin7(i) = (sqrt(2)/m)*sin((2*pi*(i-1)*7)/m);
end

%Código de plotagem dos gráficos retirado do exemplo do professor
figure(5); 
subplot(221); stem(Coef_b_cos1); grid; ylabel('Coefs. cosseno'); 
title('Coefs. da 1a. harm.');
subplot(222); stem(Coef_b_cos3); grid;
title('Coefs. da 3a. harm.');

subplot(223); stem(Coef_b_sin1); grid; ylabel('Coefs. seno');
subplot(224); stem(Coef_b_sin3); grid;

x = zeros(1,m);
for z = 1:1:length(signal_dig)
    %Zerando variaveis para evitar problemas na recursão
    y_r1 = 0;
    y_im1 = 0;
    y_r3 = 0;
    y_im3 = 0;
    for j = 2:1:m
      x(m-(j-2)) = x(m-(j-1));
      %janela deslizante
    end
    x(1) = signal_dig(z);
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

% APRESENTAÇÃO DOS RESULTADOS: 
for i = 1:length(yy1)-1
    % Verificando se o elemento i é igual ao próximo elemento i+1
    if yy1(i) == yy1(i+1)
        tempoz = i; %amostra em que o regime permanente é estabelecido
        Tempoo = tempoz*taFiltro;
        disp("Regime Permanente estabelecido na amostra"); disp(i);
        disp("Ou seja o tempo de atraso digital é (em segundos):"); disp(Tempoo);
        break;
    end
end

na = 1:1:length(signal_dig);

% Para uma melhor padronização dos resultados, utilizamos o layout dos gráficos que o professor definiu em seu modelo.
figure(6); subplot(211);
plot(na, signal_dig-2048, na, yy1); grid; ylabel('Unidades do AD');
legend('Sinal digital sem offset', 'Modulo do fasor (eficaz)');
subplot(212);
plot(na, (signal_dig-2048)*q, na, yy1*q); grid; ylabel('Volts');
legend('Sinal original', 'Modulo do fasor (eficaz)');

figure(7);
plot(na, yy1*q*sqrt(2), na, yy3*q*sqrt(2)); grid; 
title('Modulos da 1a e 3a hamonica');
legend('Fundamental','3ª harmônica');
ylabel('Volts de pico');

disp('Módulo do fasor da 1a. harm.  na amostra 60');  disp(yy1(60)*q*sqrt(2));
disp('Módulo do fasor da 3a. harm.   na amostra 60');  disp(yy3(60)*q*sqrt(2));
