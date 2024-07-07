%% Limpa a memoria e carrega os dados
clear % Limpa a memoria
clc % Limpa a janela de comando
tic
% Carrega os dados da planta
load('u.mat'); % Carrega os dados de entrada
load('y.mat'); % Carrega os dados de saida
load('u_.mat'); % Carrega dados da curva estática
load('y_.mat'); % Carrega dados da curva estática
t = (1:length(u))'; % Vetor de amostras
%% Inicializa as variaveis do Algoritmo GK
c = 12; % Numero de clusters
rou = ones(c,1);
N = length(u); % Tamanho dos dados
m = 1.2; % Poderamento exponencial
epsilon = 1e-2; % Tolerancia para termino do algoritmo
Z = [u';y']; % Matriz de dados Z
n = length(Z(:,1)); % Dimensao dos vetores de dados
l = 1; % primeira iteracao
% % % % % Norma Euclidiana
% % % % % Inicializa uma nova matriz de particao aleatoria
U{l} = rand(c,N); % Inicializa a matriz de particao fuzzy
for j = 1:N % Normaliza a matriz de particao
    U{l}(:,j) = U{l}(:,j)/(sum(U{l}(:,j)));
end
%save('U.mat','U'); % Salva a matriz de particao
% % % % %
% % % % % Carrega uma mesma matriz de particao anterior
%load('U.mat'); % Carrega uma mesma matriz de particao anterior
% % % % %
%% Realiza o algoritmo GK
prox_it = 1;
while prox_it
    l = l + 1 % Incrementa a iteracao
    V{l,1} = zeros(n,c); % Inicializa a matriz de prototipos atual
    for i = 1:c % Calcula a matriz de prototipos atual
        num = zeros(n,1); den = 0; % Para calcular o somatorio
        for k = 1:N % Efetua os somatorios
             num = num + U{l-1}(i,k)^m*Z(:,k);
             den = den + U{l-1}(i,k)^m;
        end
        V{l}(:,i) = num/den; % Calcula o i-esimo vetor de prototipos atual
    end
    for i = 1:c % Calcula a matriz de covariância fuzzy
        num = zeros(n,n); den = 0; % Para calcular o somatorio
        for k = 1:N % Efetua os somatorios
            num = num + U{l-1}(i,k)^m*(Z(:,k)-V{l}(:,i))*(Z(:,k)-V{l}(:,i))';
            den = den + U{l-1}(i,k)^m;
        end
        F{i} = num/den; % Calcula a i-esima matriz de covariancia fuzzy
    end
    D = zeros(c,N); % Inicializa a Matriz de distancias
    for i = 1:c % Calcula as distancias dos dados aos prototipos
        for k = 1:N
            D(i,k) = sqrt((Z(:,k)-V{l}(:,i))'*(rou(i)*det(F{i}))^(1/n)*inv(...
                F{i})*(Z(:,k)-V{l}(:,i)));
        end
    end
    U{l,1} = zeros(c,N); % Inicializa a matriz de particao atual
    for k = 1:N
        % Armazena os indices de elementos nulos em D(:,k)
        indcn = find(D(:,k)==0);
        if isempty(indcn) % Se nao houver elementos nulos
            for i = 1:c
                for j = 1:c
                    U{l}(i,k) = U{l}(i,k) + (D(i,k)/D(j,k))^(2/(m-1));
                end
                U{l}(i,k) = 1/U{l}(i,k);
            end
        else % Se houver elementos nulos em D(:,k)
            % Divide igualmente a pertinencia nos clusters em que D(i,k)=0,
            % e deixa pertinencia 0 nos demais clusters
            for i = indcn
                U{l}(i,k) = 1/(length(indcn));
            end
        end
    end
    norma = norm(U{l}-U{l-1})
    %if max(max(abs(U{l}-U{l-1}))) < epsilon
    if norma < epsilon
        prox_it = 0;
    end
end
tempo = toc;
figure
plot(V{end}(1,:),V{end}(2,:),'b*'); % Plota os prototipos
hold on
plot(u_,y_,'k');
centros_GK = V{end};
save('centros_GK.mat','centros_GK');