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
max_u = 5; min_u = 0; % Valor maximo e minimo da entrada para normalizacao
max_y = 2; min_y = -2; % Valor maximo e minimo da saida para normalizacao
u_n = (u-min_u)/(max_u-min_u); % normalizando
y_n = (y-min_y)/(max_y-min_y); % normalizando
%% Inicializa as variaveis do Algoritmo FCM
c = 12; % Numero de clusters
N = length(u); % Tamanho dos dados
m = 1.2; % Poderamento exponencial
w = 1.2;
a = 0.5; b = 0.5;
epsilon = 1e-2; % Tolerancia para termino do algoritmo
Z = [u';y']; % Matriz de dados Z
n = length(Z(:,1)); % Dimensao dos vetores de dados
l = 1; % primeira iteracao
% % % % % Norma Euclidiana
A = eye(n); %Matriz de norma induzida (Identidade)
% % % % % 
% % % % Norma diagonal
% dp = zeros(n,1); % Calcula os desvios padrão dos dados
% for i = 1:n
%     dp(i) = std(Z(i,:));
% end
% A = diag((1./dp).^2); % Matriz de norma induzida (Diagonal)
% A(2,2) = 1e-6; % Apelacao
% % % %
% % % % % Norma de Mahalonobis
%Zb = mean(Z,2); % Media dos dados zi, i = 1:N
%R = zeros(n,n); % Inicializa a matriz de covariancia
% for k = 1:N
%     R = R + (Z(:,k)-Zb)*(Z(:,k)-Zb)';
% end
% R = R/N;
% A = inv(R); % Matriz de norma induzida (Mahalonobis)
% % % % %
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
% Calcula os primeiros centros
V{l,1} = zeros(n,c); % Inicializa a matriz de prototipos atual
for i = 1:c % Calcula a matriz de prototipos atual
    num = zeros(n,1); den = 0; % Para calcular o somatorio
    for k = 1:N % Efetua os somatorios
        num = num + U{l}(i,k)^m*Z(:,k);
        den = den + U{l}(i,k)^m;
    end
    V{l}(:,i) = num/den; % Calcula o i-esimo vetor de prototipos atual
end
%% Realiza o algoritmo FCM
prox_it = 1;
while prox_it
    l = l + 1 % Incrementa a iteracao
    D = zeros(c,N); % Inicializa a Matriz de distancias
    for i = 1:c % Calcula as diatancias dos dados aos prototipos
        for k = 1:N
            D(i,k) = (Z(:,k)-V{l-1}(:,i))'*A*(Z(:,k)-V{l-1}(:,i));
        end
    end
    U{l,1} = zeros(c,N); % Inicializa a matriz de particao atual
    for k = 1:N
        % Armazena os indices de elementos nulos em D(:,k)
        indcn = find(D(:,k)==0);
        if isempty(indcn) % Se nao houver elementos nulos
            for i = 1:c
                for j = 1:c
                    U{l}(i,k) = U{l}(i,k) + (D(i,k)/D(j,k))^(1/(m-1));
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
    V{l,1} = zeros(n,c); % Inicializa a matriz de prototipos atual
    for i = 1:c % Calcula a matriz de prototipos atual
        num = zeros(n,1); den = 0; % Para calcular o somatorio
        for k = 1:N % Efetua os somatorios
            num = num + U{l}(i,k)^m*Z(:,k);
            den = den + U{l}(i,k)^m;
        end
        V{l}(:,i) = num/den; % Calcula o i-esimo vetor de prototipos atual
    end
    norma = norm(V{l}-V{l-1})
    %if max(max(abs(U{l}-U{l-1}))) < epsilon
    if norma < epsilon
        prox_it = 0;
    end
end

%% Realiza o algoritmo PFGG
UU{1} = U{end}; VV{1} = V{end};
prox_it = 1;
l = 1;
for i = 1:c
    num = zeros(n,n);
    den1 = 0;
    den2 = 0;
    for k = 1:N
        num = num + UU{1}(i,k)^m*(Z(:,k)-VV{1}(:,i))*(Z(:,k)-VV{1}(:,i))';
        den1 = den1 + UU{1}(i,k)^m;
        den2 = den2 + 1/N*UU{1}(i,k);
    end
    beta(i,1) = sqrt(det(num/den1))/den2;
end
while prox_it
    l = l + 1 % Incrementa a iteracao
    for i = 1:c
        num = zeros(n,n);
        den1 = 0;
        p(i,1) = 0;
        SF{i,1} = zeros(n,n);
        for k = 1:N
            num = num + UU{l-1}(i,k)^m*(Z(:,k)-VV{l-1}(:,i))*(Z(:,k)-VV{l-1}(:,i))';
            den1 = den1 + UU{l-1}(i,k)^m;
            p(i) = p(i) + 1/N*UU{l-1}(i,k);
        end
        SF{i} = num/den1;
    end
    D = zeros(c,N); % Inicializa a Matriz de distancias]
    Dx = D;
    for i = 1:c % Calcula as diatancias dos dados aos prototipos
        for k = 1:N
            D(i,k) = sqrt(det(SF{i}))/p(i)*exp(0.5*(Z(:,k)-VV{l-1}(:,i))'*inv(SF{i})*(Z(:,k)-VV{l-1}(:,i)));
            Dx(i,k) = (Z(:,k)-VV{l-1}(:,i))'*(Z(:,k)-VV{l-1}(:,i));
        end
    end
    UU{l,1} = zeros(c,N); % Inicializa a matriz de particao atual
    for k = 1:N
        % Armazena os indices de elementos nulos em D(:,k)
        indcn = find(D(:,k)==0);
        if isempty(indcn) % Se nao houver elementos nulos
            for i = 1:c
                if isinf(D(i,k))
                    UU{l}(i,k) = 0;
                else
                    for j = 1:c
                        UU{l}(i,k) = UU{l}(i,k) + (D(i,k)/D(j,k))^(1/(m-1));
                    end
                    UU{l}(i,k) = 1/UU{l}(i,k);
                end
            end
        else % Se houver elementos nulos em D(:,k)
            % Divide igualmente a pertinencia nos clusters em que D(i,k)=0,
            % e deixa pertinencia 0 nos demais clusters
            for i = indcn
                UU{l}(i,k) = 1/(length(indcn));
            end
        end
    end
    T{l} = zeros(c,N);
    for i = 1:c
        for k = 1:N
            T{l}(i,k) = 1/(1+(b/beta(i)*D(i,k))^(1/(w-1)));
        end
    end
    VV{l,1} = zeros(n,c); % Inicializa a matriz de prototipos atual
    for i = 1:c % Calcula a matriz de prototipos atual
        num = zeros(n,1); den = 0; % Para calcular o somatorio
        for k = 1:N % Efetua os somatorios
            num = num + (a*UU{l}(i,k)^m+b*T{l}(i,k)^w)*Z(:,k);
            den = den + a*UU{l}(i,k)^m+b*T{l}(i,k)^w;
        end
        VV{l}(:,i) = num/den; % Calcula o i-esimo vetor de prototipos atual
    end
    norma = norm(VV{l}-VV{l-1})
    %if max(max(abs(U{l}-U{l-1}))) < epsilon
    if norma < epsilon
        prox_it = 0;
    end
end
tempo = toc
%% Plota os resultados
figure(1)
plot(VV{end}(1,:),VV{end}(2,:),'b*'); % Plota os prototipos
hold on
plot(u_,y_,'k');
centros_PFGG = VV{end};
save('centros_PFGG.mat','centros_PFGG');