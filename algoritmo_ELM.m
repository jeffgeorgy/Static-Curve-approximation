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
xcn = [xn(:,1)]; % Primeiro centro de cluster
x = [u';y']; % Vetor de dados original
% Passa os valores de centros para a escala original
% xc{1}(1,:) = xcn{1}(1,:)*(max(u)-min(u))+min(u);
% xc{1}(2,:) = xcn{1}(2,:)*(max(temp)-min(temp))+min(temp);
N = length(xn(1,:)); % Tamando dos dados
r = 0.1; % raio de influencia de um cluster
alpha(:,1) = xn(:,1);
beta(1,1) = norm(xn(:,1))^2;
c(1,1) = 1;
sigma(1,1) = 0;
%% Realizacao do algoritmo
for k = 2:N
    k
    dist = zeros(R,1);
    s1 = zeros(R,1);
    for i = 1:R
        dist(i) = norm(xn(:,k)-xcn(:,i));
        if dist(i) < (max(sigma(i),r)+r)
            s1(i) = 1;
        end
    end
    if ~norm(s1)
        R = R + 1;
        xcn(:,R) = xn(:,k);
        sigma(R,1) = 0;
        c(R,1) = 1;
        alpha(:,R) = xn(:,k);
        beta(R,1) = norm(xn(:,k))^2;
    else
        dmin = 1e6;
        p = 0;
        mean = 0;
        variance = 0;
        for i = 1:R
            if s1(i)==1 && dist(i)<dmin
                dmin = dist(i);
                p = i;
            end
        end
        beta(p,1) = beta(p,1) + norm(xn(:,k))^2;
        alpha(:,p) = alpha(:,p) + xn(:,k);
        mean = (c(p)*xcn(:,p) + xn(:,k))/(c(p)+1);
        variance = (beta(p)+c(p)*norm(mean)^2-2*mean'*alpha(:,p))/(c(p)+1);
        xcn(:,p) = mean;
        sigma(p,1) = variance;
        c(p) = c(p) + 1;
        dist = zeros(R,1);
        s2 = zeros(R,1);
        for j = 1:R
            dist(j) = norm(xcn(:,p)-xcn(:,j));
            if dist(j) < (max(sigma(p),r)+max(sigma(j),r)) && j~=p
                s2(j) = 1;
            end
        end
        if norm(s2)
            dmin = 1e6;
            q = 0;
            for i = 1:R
                if s2(i)==1 && dist(i)<dmin
                    dmin = dist(i);
                    q = i;
                end
            end
            beta(R+1,1) = beta(p) + beta(q);
            alpha(:,R+1) = alpha(:,p) + alpha(:,q);
            c(R+1,1) = c(p) + c(q);
            xcn(:,R+1) = (c(p)*xcn(:,p)+c(q)*xcn(:,q))/c(R+1);
            sigma(R+1,1) = (beta(R+1)+(c(R+1)-1)*norm(xcn(:,R+1))^2-2*xcn(:,R+1)'*alpha(:,R+1))/c(R+1);
            beta(min(p,q)) = [];
            beta(max(p,q)) = [];
            alpha(:,min(p,q)) = [];
            alpha(:,max(p,q)) = [];
            c(min(p,q)) = [];
            c(max(p,q)) = [];
            xcn(:,min(p,q)) = [];
            xcn(:,max(p,q)) = [];
            sigma(min(p,q)) = [];
            sigma(max(p,q)) = [];
            R = R - 1;
        end
    end
end
tempo = toc
xc(1,:) = xcn(1,:)*(max_u-min_u)+min_u;
xc(2,:) = xcn(2,:)*(max_y-min_y)+min_y;
% xc(1,:) = xcn{end}(1,:);
% xc(2,:) = xcn{end}(2,:);
figure
plot(xc(1,:),xc(2,:),'b*')
hold on 
plot(u_,y_,'k'); % Plota a curva estatica
centros_ELM = xc;
save('centros_ELM.mat','centros_ELM');