%%% Parametros "globais"
sampling_period = 4e-6;  % Periodo de amostragem de 4 s

tempo_deboucing = input ("Escolha um tempo de debounce (maior que 10 microssegundos!)"); %tempo de debouncing ajustavel

duration = ceil(tempo_deboucing/sampling_period)
%O parâmetro duration corresponde ao numero minimo de bits que o codigo passara a reconhecer como uma mudança real no nivel do sinal

t = (0:numel(DigData)-1) * sampling_period; %vetor auxiliar de tempos para a plotagem dos gráficos


%% Testagem para os dados no arquivos Data2.m

printf("Debouncing dos arquivos de dados 2: ")

%Leitura do arquivo:
Data2

%Chamada da funcao que faz o debounce e guarda os eventos
[debounced_data, event_log] = debounce_and_register_events(DigData, duration, sampling_period);

%Plot do grafico com o sinal sem debounce e com debounce

figure
plot(t, DigData, 'b*', t, debounced_data, 'ro');
xlabel('Tempo (s)');
ylabel('Nivel');
legend('Sinal original', 'Sinal debounced');
title('Estudo dos sinais dos dados "Dados2"');

% Apresenta os eventos:
tamanhos = size(event_log);
for i = 1:tamanhos(1) %itera as linhas
  for j = 1:(tamanhos(2)-1)
    printf("O tempo de ocorrência do evento %d foi de: %f segundo\n", i, event_log(i,j));
    if event_log(i,j+1) == 0
      printf("O nivel do sinal foi de 1 para 0 \n\n\n");
    else
      printf("O nivel do sinal foi de 0 para 1 \n\n\n");
    endif
  endfor
endfor

%% Testagem para os dados no arquivos Data3.m

printf("Debouncing dos arquivos de dados 3: ")

%Leitura do arquivo:
Data3

%Chamada da funcao que faz o debounce e guarda os eventos
[debounced_data, event_log] = debounce_and_register_events(DigData, duration, sampling_period);

%Plot do grafico com o sinal sem debounce e com debounce

figure
plot(t, DigData, 'b*', t, debounced_data, 'ro');
xlabel('Tempo (s)');
ylabel('Nivel');
legend('Sinal original', 'Sinal debounced');
title('Estudo dos sinais dos dados "Dados3"');

% Apresenta os eventos:
tamanhos = size(event_log);
for i = 1:tamanhos(1) %itera as linhas
  for j = 1:(tamanhos(2)-1)
    printf("O tempo de ocorrência do evento %d foi de: %f segundo\n", i, event_log(i,j));
    if event_log(i,j+1) == 0
      printf("O nivel do sinal foi de 1 para 0 \n\n\n");
    else
      printf("O nivel do sinal foi de 0 para 1 \n\n\n");
    endif
  endfor
endfor

figure(1)
figure(2)
