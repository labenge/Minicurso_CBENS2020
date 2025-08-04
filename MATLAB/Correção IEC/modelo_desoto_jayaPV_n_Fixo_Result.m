%clc; clear; close all;

%% Qual é o valor ideal para RUNS, pop, maxFes considerando a confiabilidade 
% da resposta e seu tempo de execução?
%%

warning('off')
%load 'IEC.mat';
%load 'Curvas_Corrig_STC_SANDIA_alfa_Beta_regress_alfa_beta_excel.mat';
% load 'Desoto_demo.mat';

%IVCurves = Result;
%% Criterios de ajuste


IVCurves = Result(find([Result.E]>800)); % Só fica as curvas acima de 800 W/m2
IVCurves = IVCurves(randi([1 length(IVCurves)],300)); % seleciona algumas curvas para diminuir o conjunto
ncurves = length(IVCurves);     % Numero de curvas 

%idealidade = 1.124650878157024; % Idealidade calculada (SANDIA 1:750 curvas)

%idealidade = 1.117198815096043;  % Idealidade calculada (SANDIA Todas as curvas)

%idealidade = 1.122152200633934;  % Idealidade calculada (SANDIA Regressão simples)

idealidade = 1.123754955255396 % Calculado regressão idealidade 

% idealidade = 1.110633153; % Idealidade calculada (SANDIA PVSyst)

% Limites para Corrente I0

% minimo = 3.3003e-09;
% %maximo = 5.0e-06;
% maximo = 1.3053e-06;

minimo = -50*10^-6;
maximo =  50*10^-6;
I0i = 1.436e-24;
I0Ti = 0.1225;
 
% I0 = I0 . e^{I0T . T)
% var =      [Iph, I0   , I0T,           n,     nT, Rs0, RsT,    Rp0,  RpT,    a]; % Rs = Rs0 + T. RsT 
limite_inf = [7,       0,   0,  idealidade,      0,   0,   0,      0,    0,    0]; % limite inferior
limite_sup = [8,  maximo,   0,  idealidade,      0,  .4,   0,   9000,    0,    0]; % limite superior

RUNS = 30;
pop = 20;
maxFes = 5000; % idealmente acima de 20.000
graphic = false;

% valores próximos a Isc
val = randi([1 14],5,1);
% valores próximos a Imp
val = [val;randi([15 41],20,1)];
% valores próximos a Voc
val = [val;randi([42 56],5,1)];
% aí se seleciona os valores, ainda tem dois problemas: os numeros se
% repetem e não estão em ordem... faz sentido algo nesse sentido?
Vmed = IVCurves.V;
Imed = IVCurves.I;
Suns = IVCurves(1).E*ones(size(Imed));
Tmed = IVCurves(1).Tc*ones(size(Imed));
for k = 2:ncurves
    Vmed = [Vmed; IVCurves(k).V];
    Imed = [Imed; IVCurves(k).I];
    Suns = [Suns; IVCurves(k).E *ones(size(IVCurves(k).I))];
    Tmed = [Tmed; IVCurves(k).Tc*ones(size(IVCurves(k).I))];
end

    %IVCurves(k).FF = (([IVCurves(k).Imp])*([IVCurves(k).Vmp]))/(([IVCurves(k).Isc])*([IVCurves(k).Voc]));
    tic;
    [Iph, I0, I0T, n, nT, Rs0, RsT, Rp0, RpT, a, RMSE] =  Jaya_PV_Mod(Vmed, Imed, Suns/1000, Specs.Ns, Tmed, limite_inf, limite_sup, RUNS, pop, maxFes);
    toc
    fprintf('\n Matriz = %i, RUNS = %i, pop = %i e maxFes = %1.0e',size(Imed,1),RUNS,pop,maxFes);
    
    TEXTO = ["Iph", "I0", "I0T", "n", "nT", "Rs0", "RsT", "Rp0", "RpT", "a"];
    
    fprintf('\n VAriáveis com o Limite Inferior Atingido');
    TEXTO(find(~([Iph, I0, I0T, n, nT, Rs0, RsT, Rp0, RpT, a]-limite_inf)))
    
    fprintf('\n VAriáveis com o Limite Superior Atingido');
    TEXTO(find(~([Iph, I0, I0T, n, nT, Rs0, RsT, Rp0, RpT, a]-limite_sup)))

    
    fprintf('\n Iph(A) = %.3e',Iph);
    fprintf('\n I0(uA) = %.3e . exp (%.3e . T)',I0,I0T);
    fprintf('\n n = %.3e + %.3e . T',n,nT);
    fprintf('\n Rs = %.3e + %.3e . T',Rs0,RsT);
    fprintf('\n Rp = %.3e + %.3e . T',Rp0,RpT);
    fprintf('\n alpha = %.3e',a);
    fprintf('\n RMSE = %.5e \n',RMSE);

    %%%%%%%%% CRIA O VETOR DE CORRENTE ESTIMADA
    k = 1.3806503e-23;    % Boltzmann [J/K]
    q = 1.60217646e-19;   % Electron charge [C]
    Iest = (Suns/1000).*(Iph * (1 + a*(Tmed - 25))) - I0 .* exp((Tmed+273.15) * I0T).*(exp((Vmed+Imed.*(Rs0 + (Tmed+273.15) * RsT))./(n.*(Specs.Ns*k*(Tmed+273.15)/q)))-1) - ...
    (Vmed+Imed.*(Rs0 + (Tmed+273.15) * RsT))./(Rp0 + (Tmed+273.15) * RpT);
    MSE = (Imed - Iest).^2;
    
    %% PLOTA a corrente medida pela estimada
    figure
    plot(Imed,Iest,'Marker','+','LineStyle','none')
    xlabel('I medida')
    ylabel('I estimada')
    
    figure
    plot(Imed,MSE,'Marker','+','LineStyle','none')
    xlabel('I medida')
    ylabel('MSE')
    
    figure
    plot(Tmed,(Rs0 + (Tmed+273.15) * RsT),'Marker','+','LineStyle','none')
    xlabel('Temperature (°C)')
    ylabel('Rs (Ohms)')
    
    figure
    plot(Tmed,(Rp0 + (Tmed+273.15) * RpT),'Marker','+','LineStyle','none')
    xlabel('Temperature (°C)')
    ylabel('Rp (Ohms)')
    
    figure
    plot(Vmed,Imed)
    hold on
    plot(Vmed,Iest,'Marker','s','LineStyle','none')
    hold off
    xlabel('Tensão')
    ylabel('Corrente')
    
%     figure
%     plot([Vmed(1000:1200);Vmed(11000:11200)],[Imed(1000:1200);Imed(11000:11200)])
%     hold on
%     plot([Vmed(1000:1200);Vmed(11000:11200)],[Iest(1000:1200);Iest(11000:11200)],'Marker','s','LineStyle','none')
%     hold off
%     xlabel('Tensão')
%     ylabel('Corrente')


%     if graphic
%       figure(k)
%      
%       plot(IVCurves(k).V, IVCurves(k).I,'o');
%       hold on
%       kb = 1.38066e-23; 
%       q = 1.60218e-19;
%       Vt = Specs.Ns*kb*(IVCurves(k).Tc+273.15)/q;
%       fmodelo = @(V,I) (I - Iph + I0*(exp((V+I*Rs)/(n*Vt))-1) + (V+I*Rs)/Rp);
%       ezplot(fmodelo,[0 IVCurves(k).Voc 0 (IVCurves(k).Isc*1.03)])
%       title(['Curva ', num2str(k)])
%       xlim([0, IVCurves(k).Voc+1])
%       ylim([0, IVCurves(k).Isc+0.3])
%       drawnow
%     end
% 
% 
% [lin, col] = size(Model);
% %% Graficos em função da temperatura
% 
% Qtd_curvas_plot = 400; % Escolhe quantidade de curvas para plotagem
% 
% if (Qtd_curvas_plot > ncurves)
%     Qtd_curvas_plot = ncurves;
% else
%     Qtd_curvas_plot = Qtd_curvas_plot;
% end
% 
% figure; 
% plot([IVCurves(1:Qtd_curvas_plot).Tc], [Model(1:Qtd_curvas_plot).Rs],'*')
% title(['Rs vs T'])
% xlabel('Temperature (°C)')
% ylabel('Rs (Ohms)')
% 
% figure; 
% plot([IVCurves(1:Qtd_curvas_plot).Tc], [Model(1:Qtd_curvas_plot).Rp],'*')
% title(['Rp vs T'])
% xlabel('Temperature (°C)')
% ylabel('Rp (Ohms)')
% 
% figure; 
% plot([IVCurves(1:Qtd_curvas_plot).Tc], [Model(1:Qtd_curvas_plot).Iph],'*')
% title(['Iph vs T'])
% xlabel('Temperature (°C)')
% ylabel('Iph (A)')
% 
% figure; 
% plot([IVCurves(1:Qtd_curvas_plot).Tc], [Model(1:Qtd_curvas_plot).I0]*10^6,'*')
% title(['I0 vs T'])
% xlabel('Temperature (°C)')
% ylabel('I0 (uA)')
% 
% figure; 
% plot([IVCurves(1:Qtd_curvas_plot).Tc], [Model(1:Qtd_curvas_plot).n],'*')
% title(['n vs T'])
% xlabel('Temperature (°C)')
% ylabel('n')
% 
% %% Subplot Temperatura
% figure;
% subplot(2,2,1)
% plot([IVCurves(1:Qtd_curvas_plot).Tc], [Model(1:Qtd_curvas_plot).Rs],'r.')
% title(['Rs vs T'])
% xlabel('Temperature (°C)')
% ylabel('Rs (Ohms)')
% 
% subplot(2,2,2)
% plot([IVCurves(1:Qtd_curvas_plot).Tc], [Model(1:Qtd_curvas_plot).Rp],'r.')
% title(['Rp vs T'])
% xlabel('Temperature (°C)')
% ylabel('Rp (Ohms)')
% 
% subplot(2,2,3) 
% plot([IVCurves(1:Qtd_curvas_plot).Tc], [Model(1:Qtd_curvas_plot).Iph],'r.')
% title(['Iph vs T'])
% xlabel('Temperature (°C)')
% ylabel('Iph (A)')
% 
% subplot(2,2,4)
% plot([IVCurves(1:Qtd_curvas_plot).Tc], [Model(1:Qtd_curvas_plot).I0]*10^6,'r.')
% title(['I0 vs T'])
% xlabel('Temperature (°C)')
% ylabel('I0 (uA)')
% 
% %% Subplot Irradiancia
% figure;
% subplot(2,2,1)
% plot([IVCurves(1:Qtd_curvas_plot).Ee], [Model(1:Qtd_curvas_plot).Rs],'b.')
% title(['Rs vs E'])
% xlabel('irradiancia (W/m^2)')
% ylabel('Rs (Ohms)')
% 
% subplot(2,2,2)
% plot([IVCurves(1:Qtd_curvas_plot).Ee], [Model(1:Qtd_curvas_plot).Rp],'b.')
% title(['Rp vs E'])
% xlabel('irradiancia (W/m^2)')
% ylabel('Rp (Ohms)')
% 
% subplot(2,2,3) 
% plot([IVCurves(1:Qtd_curvas_plot).Ee], [Model(1:Qtd_curvas_plot).Iph],'b.')
% title(['Iph vs E'])
% xlabel('irradiancia (W/m^2)')
% ylabel('Iph (A)')
% 
% subplot(2,2,4)
% plot([IVCurves(1:Qtd_curvas_plot).Ee], [Model(1:Qtd_curvas_plot).I0]*10^6,'b.')
% title(['I0 vs E'])
% xlabel('irradiancia (W/m^2)')
% ylabel('I0 (uA)')

%% Graficos em função da irradiancia
% figure; 
% plot([IVCurves(1:Qtd_curvas_plot).Ee], [Model(1:Qtd_curvas_plot).Rs],'*')
% title(['Rs vs E'])
% xlabel('irradiancia (W/m^2)')
% ylabel('Rs (Ohms)')
% 
% figure; 
% plot([IVCurves(1:Qtd_curvas_plot).Ee], [Model(1:Qtd_curvas_plot).Rp],'*')
% title(['Rp vs E'])
% xlabel('irradiancia (W/m^2)')
% ylabel('Rp (Ohms)')
% 
% figure; 
% plot([IVCurves(1:Qtd_curvas_plot).Ee], [Model(1:Qtd_curvas_plot).Iph],'*')
% title(['Iph vs E'])
% xlabel('irradiancia (W/m^2)')
% ylabel('Iph (A)')
% 
% figure; 
% plot([IVCurves(1:Qtd_curvas_plot).Ee], [Model(1:Qtd_curvas_plot).I0]*10^6,'*')
% title(['I0 vs E'])
% xlabel('irradiancia (W/m^2)')
% ylabel('I0 (uA)')
% 
% figure; 
% plot([IVCurves(1:Qtd_curvas_plot).Ee], [Model(1:Qtd_curvas_plot).n],'*')
% title(['n vs E'])
% xlabel('irradiancia (W/m^2)')
% ylabel('n')