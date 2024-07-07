%% Limpa a memoria e carrega os dados
clear % Limpa a memoria
clc % Limpa a janela de comando
close all
% Carrega os dados da planta
load('u.mat'); % Carrega os dados de entrada
load('y.mat'); % Carrega os dados de saida
load('u_.mat'); % Carrega dados da curva estática
load('y_.mat'); % Carrega dados da curva estática
load('centros_CM.mat');
load('centros_nrFCM.mat');
load('centros_GK.mat');
load('centros_PFGG.mat');
load('centros_Sub.mat');
load('centros_ETS_plus.mat');
load('centros_JG.mat');
load('centros_ELM.mat');
t = (1:length(u))'; % Vetor de amostras
[ord_centros_CM,ind_CM] = sort(centros_CM(1,:));
[ord_centros_nrFCM,ind_nrFCM] = sort(centros_nrFCM(1,:));
[ord_centros_GK,ind_GK] = sort(centros_GK(1,:));
[ord_centros_PFGG,ind_PFGG] = sort(centros_PFGG(1,:));
[ord_centros_Sub,ind_Sub] = sort(centros_Sub(1,:));
[ord_centros_JG,ind_JG] = sort(centros_JG(1,:));
[ord_centros_ELM,ind_ELM] = sort(centros_ELM(1,:));
[ord_centros_ETS_plus,ind_ETS_plus] = sort(centros_ETS_plus(1,:));
ord_centros_CM(2,:) = centros_CM(2,ind_CM);
ord_centros_nrFCM(2,:) = centros_nrFCM(2,ind_nrFCM);
ord_centros_GK(2,:) = centros_GK(2,ind_GK);
ord_centros_PFGG(2,:) = centros_PFGG(2,ind_PFGG);
ord_centros_Sub(2,:) = centros_Sub(2,ind_Sub);
ord_centros_JG(2,:) = centros_JG(2,ind_JG);
ord_centros_ELM(2,:) = centros_ELM(2,ind_ELM);
ord_centros_ETS_plus(2,:) = centros_ETS_plus(2,ind_ETS_plus);

for i = 1:length(ord_centros_CM(1,:))
    if ord_centros_CM(1,i) <= 1.5
        saida_CM(1,i) = 2*tanh(2*ord_centros_CM(1,i));
        erro_CM(1,i) = saida_CM(1,i) - ord_centros_CM(2,i);
    else
        saida_CM(1,i) = -2*(exp(ord_centros_CM(1,i))-1)/(exp(ord_centros_CM(1,i))+1);
        erro_CM(1,i) = saida_CM(1,i) - ord_centros_CM(2,i);
    end
    if ord_centros_nrFCM(1,i) <= 1.5
        saida_nrFCM(1,i) = 2*tanh(2*ord_centros_nrFCM(1,i));
        erro_nrFCM(1,i) = saida_nrFCM(1,i) - ord_centros_nrFCM(2,i);
    else
        saida_nrFCM(1,i) = -2*(exp(ord_centros_nrFCM(1,i))-1)/(exp(ord_centros_nrFCM(1,i))+1);
        erro_nrFCM(1,i) = saida_nrFCM(1,i) - ord_centros_nrFCM(2,i);
    end
    if ord_centros_GK(1,i) <= 1.5
        saida_GK(1,i) = 2*tanh(2*ord_centros_GK(1,i));
        erro_GK(1,i) = saida_GK(1,i) - ord_centros_GK(2,i);
    else
        saida_GK(1,i) = -2*(exp(ord_centros_GK(1,i))-1)/(exp(ord_centros_GK(1,i))+1);
        erro_GK(1,i) =  saida_GK(1,i)  - ord_centros_GK(2,i);
    end
    if ord_centros_PFGG(1,i) <= 1.5
        saida_PFGG(1,i) = 2*tanh(2*ord_centros_PFGG(1,i));
        erro_PFGG(1,i) = saida_PFGG(1,i) - ord_centros_PFGG(2,i);
    else
        saida_PFGG(1,i) = -2*(exp(ord_centros_PFGG(1,i))-1)/(exp(ord_centros_PFGG(1,i))+1);
        erro_PFGG(1,i) = saida_PFGG(1,i) - ord_centros_PFGG(2,i);
    end
    if ord_centros_Sub(1,i) <= 1.5
        saida_Sub(1,i) = 2*tanh(2*ord_centros_Sub(1,i));
        erro_Sub(1,i) = saida_Sub(1,i) - ord_centros_Sub(2,i);
    else
        saida_Sub(1,i) = -2*(exp(ord_centros_Sub(1,i))-1)/(exp(ord_centros_Sub(1,i))+1);
        erro_Sub(1,i) =  saida_Sub(1,i)  - ord_centros_Sub(2,i);
    end
    if ord_centros_JG(1,i) <= 1.5
        saida_JG(1,i) = 2*tanh(2*ord_centros_JG(1,i));
        erro_JG(1,i) = saida_JG(1,i) - ord_centros_JG(2,i);
    else
        saida_JG(1,i) = -2*(exp(ord_centros_JG(1,i))-1)/(exp(ord_centros_JG(1,i))+1);
        erro_JG(1,i) =  saida_JG(1,i)  - ord_centros_JG(2,i);
    end
    if ord_centros_ELM(1,i) <= 1.5
        saida_ELM(1,i) = 2*tanh(2*ord_centros_ELM(1,i));
        erro_ELM(1,i) = saida_ELM(1,i) - ord_centros_ELM(2,i);
    else
        saida_ELM(1,i) = -2*(exp(ord_centros_ELM(1,i))-1)/(exp(ord_centros_ELM(1,i))+1);
        erro_ELM(1,i) = saida_ELM(1,i) - ord_centros_ELM(2,i);
    end
    if ord_centros_ETS_plus(1,i) <= 1.5
        saida_ETS_plus(1,i) = 2*tanh(2*ord_centros_ETS_plus(1,i));
        erro_ETS_plus(1,i) = saida_ETS_plus(1,i) - ord_centros_ETS_plus(2,i);
    else
        saida_ETS_plus(1,i) = -2*(exp(ord_centros_ETS_plus(1,i))-1)/(exp(ord_centros_ETS_plus(1,i))+1);
        erro_ETS_plus(1,i) =  saida_ETS_plus(1,i)  - ord_centros_ETS_plus(2,i);
    end
end

figure
plot(centros_CM(1,:),centros_CM(2,:),'b*','LineWidth',1.5)
hold on
plot(u_,y_,'k','LineWidth',1);
legend('Clusters centers','Static curve')
xlabel('u')
ylabel('y')
figure
plot(centros_nrFCM(1,:),centros_nrFCM(2,:),'b*','LineWidth',1.5)
hold on
plot(u_,y_,'k','LineWidth',1);
legend('Clusters centers','Static curve')
xlabel('u')
ylabel('y')
figure
plot(centros_GK(1,:),centros_GK(2,:),'b*','LineWidth',1.5)
hold on
plot(u_,y_,'k','LineWidth',1);
legend('Clusters centers','Static curve')
xlabel('u')
ylabel('y')
figure
plot(centros_PFGG(1,:),centros_PFGG(2,:),'b*','LineWidth',1.5)
hold on
plot(u_,y_,'k','LineWidth',1);
legend('Clusters centers','Static curve')
xlabel('u')
ylabel('y')
figure
plot(centros_Sub(1,:),centros_Sub(2,:),'b*','LineWidth',1.5)
hold on
plot(u_,y_,'k','LineWidth',1);
legend('Clusters centers','Static curve')
xlabel('u')
ylabel('y')
figure
plot(centros_JG(1,:),centros_JG(2,:),'b*','LineWidth',1.5)
hold on
plot(u_,y_,'k','LineWidth',1);
legend('Clusters centers','Static curve')
xlabel('u')
ylabel('y')
figure
plot(centros_ETS_plus(1,:),centros_ETS_plus(2,:),'b*','LineWidth',1.5)
hold on
plot(u_,y_,'k','LineWidth',1);
legend('Clusters centers','Static curve')
xlabel('u')
ylabel('y')
figure
plot(centros_ELM(1,:),centros_ELM(2,:),'b*','LineWidth',1.5)
hold on
plot(u_,y_,'k','LineWidth',1);
legend('Clusters centers','Static curve')
xlabel('u')
ylabel('y')
%legend('$C$-\textit{means}','GK','Subtrativo','Algoritmo proposto','ETS+',"Curva est\'atica",'Interpreter','latex')
% xlabel('u','Interpreter','latex')
% ylabel('y','Interpreter','latex')

% figure
% plot(ord_centros_CM(1,:),erro_CM,'b-.')
% hold on
% plot(ord_centros_GK(1,:),erro_GK,'b--')
% plot(ord_centros_Sub(1,:),erro_Sub,'b')
% plot(ord_centros_JG(1,:),erro_JG,'g')
% plot(ord_centros_ETS_plus(1,:),erro_ETS_plus,'g--')
% legend('C-means','GK','Subtractive','Proposed algorithm','ETS+')

% figure
% plot(ord_centros_Sub(1,:),erro_Sub,'b')
% hold on
% plot(ord_centros_JG(1,:),erro_JG,'g')
% legend('Subtractive','Proposed algorithm')
% xlabel('Input','FontSize',12)
% ylabel('Error','FontSize',12)

RMSE_CM = sqrt(mean(erro_CM.^2))
RMSE_nrFCM = sqrt(mean(erro_nrFCM.^2))
RMSE_GK = sqrt(mean(erro_GK.^2))
RMSE_PFGG = sqrt(mean(erro_PFGG.^2))
RMSE_Sub = sqrt(mean(erro_Sub.^2))
RMSE_JG = sqrt(mean(erro_JG.^2))
RMSE_ELM = sqrt(mean(erro_ELM.^2))
RMSE_ETS_plus = sqrt(mean(erro_ETS_plus.^2))

MAE_CM = mean(abs(erro_CM))
MAE_nrFCM = mean(abs(erro_nrFCM))
MAE_GK = mean(abs(erro_GK))
MAE_PFGG = mean(abs(erro_PFGG))
MAE_Sub = mean(abs(erro_Sub))
MAE_JG = mean(abs(erro_JG))
MAE_ELM = mean(abs(erro_ELM))
MAE_ETS_plus = mean(abs(erro_ETS_plus))

VAF_CM = 100*(1-var(erro_CM)/var(saida_CM))
VAF_nrFCM = 100*(1-var(erro_nrFCM)/var(saida_nrFCM))
VAF_GK = 100*(1-var(erro_GK)/var(saida_GK))
VAF_PFGG = 100*(1-var(erro_PFGG)/var(saida_PFGG))
VAF_Sub = 100*(1-var(erro_Sub)/var(saida_Sub))
VAF_JG = 100*(1-var(erro_JG)/var(saida_JG))
VAF_ELM = 100*(1-var(erro_ELM)/var(saida_ELM))
VAF_ETS_plus = 100*(1-var(erro_ETS_plus)/var(saida_ETS_plus))

