%% Limpa a memoria e carrega os dados
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
% Amostras em regime permanente
max_u = 5; min_u = 0; % Valor maximo e minimo da entrada para normalizacao
max_y = 2; min_y = -2; % Valor maximo e minimo da saida para normalizacao
u_n = (u-min_u)/(max_u-min_u); % normalizando
y_n = (y-min_y)/(max_y-min_y); % normalizando
%% Inicializa as variaveis do Algoritmo Subtrativo
x = [u';y']; % Vetor de dados
x_n = [u_n';y_n']; % Vetor de dados normalizados
ra = 0.07; % Raio de vizinhança para calculo do potencial
alpha = 4/ra^2; % Parametro do calculo do potencial dependente do raio
rb = 1.5*ra; % Raio de vizinhança para subtração dos potenciais
beta = 4/rb^2; % Parametro do calculo subtrativo dos potenciais
epsilon_sup = 0.5; % Parametro que representa o percentual minimo do 
% primeiro potencial a ser sempre aceito como novo centro de cluster
epsilon_inf = 0.14; % Parametro que representa o percentual maximo do 
% primeiro potencial em que um novo centro e sempre rejeitado
N = length(x_n(1,:)); % Quantidade de dados
%% Inicia as iteracoes do algoritmo subtrativo
prox_it = 1; % Habilita a proxima iteracao
l = 1; % primeira iteracao
while prox_it
    l
    if l == 1 % Se for a primeira iteracao, calcular todos os potenciais
        % iniciais pela formula exponencial
        P{l} = zeros(N,1);
        for i = 1:N
            i
            for j = 1:N
                P{l}(i) = P{l}(i) + exp(-alpha*norm(x_n(:,i)-x_n(:,j))^2);
            end
        end
        % load('P.mat'); % carrega os potenciais
        c = 1; % Inicializa a quantidade de clusters
        indc_centros(c,1) = find(P{l}==max(P{l})); % Encontra o indice do
        % centro do primeiro cluster
    else % Se não, calcular os potenciais por subtracao e verificar se o
        % dado de potencial maximo sao aceitos como novos centros
        P{l} = zeros(N,1); %Inicializa os potenciais atuais
        for i = 1:N % Calcula os novos potenciais de forma subtrativa
            P{l}(i) = P{l-1}(i) - P{l-1}(indc_centros(c))*exp(-beta*...
                norm(x_n(:,i)-x_n(:,indc_centros(c)))^2);
        end
        sit_laco_def = 0;
        while ~sit_laco_def % Enquanto nao for definido se teremos novo
            % centro ou se finalizaremos com os clusters ja existentes
            indc_cand = find(P{l}==max(P{l})); % Candidato a novo centro
            if P{l}(indc_cand(1)) > epsilon_sup*P{1}(indc_centros(1)) % Se o
                % potencial for maior que o limite superior, adicionar novo
                % centro de cluster
                c = c + 1; % atualiza a quantidade de clusters
                indc_centros(c,1) = indc_cand(1); % identifica o ind do novo
                % centro
                sit_laco_def = 1; % Sai do laco while no final, pois encon-
                %trou-se o centro
            elseif P{l}(indc_cand(1)) < epsilon_inf*P{1}(indc_centros(1))
                % Se for menor que o limite inferior, rejeitar centro
                sit_laco_def = 1; % Sai do laco while no final, pois defi-
                % niu-se que nao havera mais centros
                prox_it = 0; % finaliza o algoritmo no final do laco extern
            else % Se nao, verifica a distancia para o centro mais proximo
                % e analisa as condicoes para ser aceito ou rejeitado
                dmin = 2; % inicializa dmin
                for i = indc_centros % Busca a dist para o centr mais prxm
                    if norm(x_n(:,indc_cand(1))-x_n(:,i)) < dmin
                        dmin = norm(x_n(:,indc_cand(1))-x_n(:,i));
                    end
                    % Encontra a distancia minima no final do laco for
                end
                if (dmin/ra + P{l}(indc_cand(1))/P{1}(indc_centros(1))) >= 1
                % Se essa condicao for satisfeita, aceita-se o centro
                c = c + 1; % atualiza a quantidade de clusters
                indc_centros(c,1) = indc_cand(1); % identifica o ind do novo
                % centro
                sit_laco_def = 1; % Sai do laco while no final, pois encon-
                %trou-se o centro
                else 
                    % Se não, rejeita o centro e testa o proximo de maior
                    % potencial, sit_laco_def permanece em 0.
                    P{l}(indc_cand(1)) = 0; % Zera o potencial do centro 
                    % rejeitado
                end
            end
        end
    end  
    l = l + 1;
end
tempo = toc
% Determinando as pertinências dos dados em relação aos clusters
% encontrados
U = zeros(c,N); %Inicializando a matriz de pertinencias
for k = 1:N
    for i = 1:c
        % Calcula a pertinencia do dado k em relacao a todos os clusters
        % usando a funcao exponencial
        U(i,k) = exp(-alpha*norm(x_n(:,k)-x_n(:,indc_centros(i)))^2);
    end
    % Normaliza para que a soma das pertinencias seja 1 
    U(:,k) = U(:,k)/sum(U(:,k));
end
figure
plot(x(1,indc_centros),x(2,indc_centros),'b*'); % Plota os centros
hold on
plot(u_,y_,'k'); % Plota a curva estatica
centros_Sub = x(:,indc_centros);
save('centros_Sub.mat','centros_Sub');