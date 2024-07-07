%% Limpa a memoria e carrega os dados
clear % Limpa a memoria
clc % Limpa a janela de comando
close all
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
% u_n = u;
% y_n = y;
%% Inicializa as variaveis do Algoritmo ETS
k = 1; % Amostra atual
R = 1; % Numero de regras inicial
xn = [u_n';y_n']; % Vetor de dados normalizados
n = length(xn(:,1)); % Tamanho de dados de entrada
xcn{1} = [xn(:,1)]; % Primeiro centro de cluster
x = [u';y']; % Vetor de dados original
% Passa os valores de centros para a escala original
% xc{1}(1,:) = xcn{1}(1,:)*(max(u)-min(u))+min(u);
% xc{1}(2,:) = xcn{1}(2,:)*(max(temp)-min(temp))+min(temp);
D(1) = 1; % Potencial do primeiro dado
Dc{1} = D; % Potencial do primeiro centro
N = length(xn(1,:)); % Tamando dos dados
r = 0.25; % raio de influencia de um cluster
alpha = 4/r^2; % Parametro dependente de r
c(1) = 0; % Valor inicial para calculo do potencial
beta(:,1) = zeros(n,1); % Valor inicial para calculo do potencial
Nt = 0;
sigma(:,1) = ones(n,1);
sigma2 = sigma.^2;
raio = r^2/sqrt(8)*ones(n,1);
S = 1;
rho = 0.5;
%% Realizacao do algoritmo
for k = 2:N
    k
    % Calculo dos potenciais
    vv(k,1) = norm(xn(:,k))^2; % Calcula vv(k) atual
    c(k,1) = c(k-1) + norm(xn(:,k-1))^2; % calcula sigma(k) recursi-
    % vamente
    beta(:,k) = beta(:,k-1) + xn(:,k-1); % Calculo recursivo de beta(k)
    v(k,1) = xn(:,k)'*beta(:,k); % Calcula v(k) atual
    D(k,1) = (k-1)/((k-1)*(vv(k)+1)+c(k)-2*v(k)); % Calculo do 
    % potencial do dado atual
%     plot(k,P(k),'bo')
    for l = 1:R % Atualiza os potenciais dos centros ja existentes com
        % base no novo dado
        Dc{k,1}(l,1) = (k-1)/(k-1+(k-2)*(1/Dc{k-1}(l)-1)+...
           norm(xn(:,k)-xcn{k-1}(:,l))^2);
        % Dc{k,1}(l,1) = (k-1)/(k-1+(k-2)*(1/Dc{k-1}(l)-1)+...
        %     norm(xn(:,k)-xn(:,k-1))^2);
    end
    xcn{k,1} = xcn{k-1}; % Deixa inicialmente os centros no instante atual
    % iguais aos do instante anterior (sera ou nao modificado)
    % Avalia as condicoes para evolucao da estrutura
    if D(k) > max(Dc{k}) || D(k) < min(Dc{k}) % Se o potencial do dado atual for maior que todos
        % os potenciais de centros
        R = R + 1; % Aumenta o numero de regras
        xcn{k,1}(:,R) = xn(:,k); % Centro da nova regra
        Dc{k,1}(R,1) = D(k); % Potencial do novo centro
        sigma(:,R) = ones(n,1);
        sigma2(:,R) = sigma(:,R).^2;
        raio(:,R) = r^2/sqrt(8)*ones(n,1);
        S(R,1) = 0;
        for i = 1:(R-1)
            condicaoB = 1;
            for j = 1:n
                mu(j,i) = exp(-(xcn{k}(j,i)-xn(j,k))^2/(2*raio(j,i)^2));
                if mu(j,i) < exp(-1)
                    condicaoB = 0;
                end
            end
            if condicaoB == 1
                indc_rem = i;
                break
            end
        end
        if condicaoB == 1
            xcn{k,1}(:,indc_rem) = []; % Centro da nova regra
            Dc{k,1}(indc_rem) = []; % Potencial do novo centro
            sigma(:,indc_rem) = [];
            sigma2(:,indc_rem) = [];
            raio(:,indc_rem) = [];
            mu(:,indc_rem) = [];
            S(indc_rem) = [];
            R = R-1;
        end
    end
    mu = ones(R,1);
    mu_max = 0;
    for i = 1:R
        for j = 1:n
            mu(i) = mu(i)*exp(-(xn(j,k)-xcn{k}(j,i))^2/(2*raio(j,i)^2));
        end
        if mu(i) > mu_max
            mu_max = mu(i);
            indc_prox = i;
        end
    end
    S(indc_prox) = S(indc_prox) + 1;
    for i = 1:n
        sigma2(i,indc_prox) = (S(indc_prox)-1)/S(indc_prox)*sigma2(i,indc_prox) + ...
            1/S(indc_prox)*(xn(i,k)-xcn{k}(i,indc_prox))^2;
        sigma(i,indc_prox) = sqrt(sigma2(i,indc_prox));
    end
    raio(:,indc_prox) = rho*raio(:,indc_prox) + (1-rho)*sigma(:,indc_prox);
end
tempo = toc
xc(1,:) = xcn{end}(1,:)*(max_u-min_u)+min_u;
xc(2,:) = xcn{end}(2,:)*(max_y-min_y)+min_y;
% xc(1,:) = xcn{end}(1,:);
% xc(2,:) = xcn{end}(2,:);
figure
plot(xc(1,:),xc(2,:),'b*')
hold on 
plot(u_,y_,'k'); % Plota a curva estatica
centros_ETS_plus = xc;
save('centros_ETS_plus.mat','centros_ETS_plus');