clear all; close all; %clc;
%%
warning('off')
load 'Desoto_demo.mat';

curvas = IVCurves;

%% Escolha do Modelo
% 1 - Fixo;              n Calculado [Sandia]
% 2 - Temperatura(1);    n = n0 + (T * nT) Modelo [PVSyst]
% 3 - Temperatura(2);    n = n0 * (T/T0) [bai, 2014]
Model_n = 1;

% 1 - Fixo;            Rs Fixo [Sandia]
% 2 - Temperatura;     Rs Varia com Temperatura  Rs = Rs0+(T*RsT)
% 3 - Temp e Irrad;    Rs(T,E) Rs = Rs0*(T/T0)*[1 - b*ln(E/E0) [tossa, 2014]       
Model_Rs = 3;

% 1 - Temperatura;     Rp = Rp0 + (RpT * T)   
% 2 - Irradiância;     Rp = Rp0 * (E0/E)     [Sandia]
Model_Rp = 2;


%% Criterios de ajuste

% Seleciona curvas no intervalo de temperatura (28 a 65 ºC) de 1 em 1 grau
IntervT = [2872 2875 2895 2881 1493 1495 1496 1497 2287 5 8 10 2296 3571 30 34 3270 ...
    1535 730 78 758 66 2211 1113 1437 3117 1727 508 1799 607 1331 438 1142 ...
    1929 858 1300 3351 1867];
% Seleciona curvas no intervalo de Irradiancia (400W/m^2 a 1357W/m^2) de 20 em 20W/m^2 
IntervEe = [3525 2874 2731 749 747 3571 1485 2790 2286 1480 3543 1476 731 3556 724 ...
    3423 21 3532 2231 2415 38 2545 695 64 2902 1424 1546 776 1399 1381 2531 2957 ...
    2105 878 1836 368 3106 2629 2750 3293 3231 2624 3201 2680 3254 2679 2686];
InterCurvas = [IntervT IntervEe];

%IVCurves = IVCurves(find([IVCurves.E]>400)); % Só fica as curvas acima de 800 W/m2
%IVCurves = IVCurves(randi([1 length(IVCurves)],100,1)); % seleciona algumas curvas para diminuir o conjunto
IVCurves = [IVCurves([1 10 30 50 100 150 380])];
%IVCurves = [IVCurves([InterCurvas])];
%IVCurves = [IVCurves(1:end)];
ncurves = length(IVCurves);     % Numero de curvas 

idealidade = 1.123754955255396; % Calculado regressão idealidade 

%% Definição do limites do Modelo
inf_Iph = 7;  inf_I0 = 1e-9; inf_Eg0 = 0.8;  inf_a = 0;   inf_b = 0;
sup_Iph = 8;  sup_I0 = 1e-6;  sup_Eg0 = 1.15; sup_a = 0.1; sup_b = 0;

if (Model_n == 1)
    inf_n = idealidade; sup_n = idealidade;
    inf_nT = 0; sup_nT = 0;
end
if (Model_n == 2)
    inf_n = 1; sup_n = 2;
    inf_nT = 0; sup_nT = 0.1;
end
if (Model_n == 3)
    inf_n = 1; sup_n = 2;
    inf_nT = 0; sup_nT = 0.1;
end

if (Model_Rs == 1)
    inf_Rs0 = 0.2; sup_Rs0 = 0.3;
    inf_RsT = 0; sup_RsT = 0;
end
if (Model_Rs == 2)
    inf_Rs0 = 0.2; sup_Rs0 = 0.3;
    inf_RsT = 0; sup_RsT = 0.1;
end
if (Model_Rs == 3)
    inf_Rs0 = 0.2; sup_Rs0 = 0.3;
    inf_RsT = 0; sup_RsT = 0;
    inf_b = -1; sup_b = 0;
end

if (Model_Rp == 1)
    inf_Rp0 = 10; sup_Rp0 = 900;
    inf_RpT = -3; sup_RpT = 3;
end
if (Model_Rp == 2)
    inf_Rp0 = 0; sup_Rp0 = 900;
    inf_RpT = 0; sup_RpT = 0;
end
% var =      [Iph,     I0,      Eg0,     n,     nT,     Rs0,     RsT,     Rp0,     RpT,     a,     b    ]; 
limite_inf = [inf_Iph, inf_I0,  inf_Eg0, inf_n, inf_nT, inf_Rs0, inf_RsT, inf_Rp0, inf_RpT, inf_a, inf_b]; % limite inferior
limite_sup = [sup_Iph, sup_I0,  sup_Eg0, sup_n, sup_nT, sup_Rs0, sup_RsT, sup_Rp0, sup_RpT, sup_a, sup_b]; % limite superior

RUNS = 30;
pop = 20;
maxFes = 20000; % idealmente acima de 20.000
graphic = false;

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

    tic;
    [Iph, I0, Eg0, n, nT, Rs0, RsT, Rp0, RpT, a, b, RMSE, xmelhores] =  Jaya_PV_Mod(Vmed, Imed, Suns, Specs.Ns, Tmed, limite_inf, limite_sup, RUNS, pop, maxFes, Model_n, Model_Rs, Model_Rp);
    toc
    fprintf('\n Matriz = %i, RUNS = %i, pop = %i e maxFes = %1.0e',size(Imed,1),RUNS,pop,maxFes);
    
    TEXTO = ["Iph", "I0", "Eg0", "n", "nT", "Rs0", "RsT", "Rp0", "RpT", "a", "b"];
    
    fprintf('\n VAriáveis com o Limite Inferior Atingido');
    TEXTO(find(~([Iph, I0, Eg0, n, nT, Rs0, RsT, Rp0, RpT, a, b]-limite_inf)))
    
    fprintf('\n VAriáveis com o Limite Superior Atingido');
    TEXTO(find(~([Iph, I0, Eg0, n, nT, Rs0, RsT, Rp0, RpT, a, b]-limite_sup)))

        %%%%%%%%% CRIA O VETOR DE CORRENTE ESTIMADA
    k = 1.3806503e-23;    % Boltzmann [J/K]
    q = 1.60217646e-19;   % Electron charge [C]
    keV = k * 6.24150934e18; % Converção J/K to eV/K
    dEgdT = 0.0002677; % dependência da temperatura para energia do bangap
    
    Iest = (Suns/1000).*(Iph * (1 + a*(Tmed - 25))) - I0 .* ((Tmed+273.15)/298.15).^3 .* exp( Eg0 /keV .* (1/298.15 - (1-dEgdT*(Tmed-25))./(Tmed+273.15))).*(exp((Vmed+Imed.*(Rs0 + (Tmed+273.15) * RsT))./(n.*(Specs.Ns*k*(Tmed+273.15)/q)))-1) - (Vmed+Imed.*(Rs0 + (Tmed+273.15) * RsT))./(Rp0 + (Tmed+273.15) * RpT);
    MSE = (Imed - Iest).^2;
    
    fprintf('\n Iph(E,T)(A) = %.3e',Iph);
    fprintf('\n I0(T) (A) = %.3e',I0);
    fprintf('\n Eg0(T) = %0.3f',Eg0);
    fprintf('\n alpha = %.3e',a);
    if Model_Rs ==3
        fprintf('\n beta = %.3e',b);
    end
    if Model_n == 1
        fprintf('\n n = %.3e ',n);
    end
    if Model_n == 2
        fprintf('\n n(T) = %.3e + %.3e . T',n,nT);
    end
    if Model_n == 3
        fprintf('\n n = %.3e .(T/T0)',n);
    end
    if Model_Rs == 1
        fprintf('\n Rs = %.3e ',Rs0);
    end
    if Model_Rs == 2
        fprintf('\n Rs(T) = %.3e + %.3e . T',Rs0,RsT);
    end
    if Model_Rs == 3
        fprintf('\n Rs(E,T) = %.3e *(T/T0)*[1 - b*ln(E/E0)',Rs0);
    end
    if Model_Rp == 1
    fprintf('\n Rp(T) = %.3e + %.3e . T',Rp0,RpT);
    end
    if Model_Rp == 2
        fprintf('\n Rp(E) = %.3e * (E0/E)',Rp0);
    end   
    fprintf('\n RMSE = %.5e \n',RMSE);


    
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
    plot(Suns,(Rp0.*(1000./Suns)),'Marker','+','LineStyle','none')
    xlabel('Temperature (°C)')
    ylabel('Rp (Ohms)')
    
    figure
    plot(Vmed,Imed)
    hold on
    plot(Vmed,Iest,'Marker','s','LineStyle','none')
    hold off
    xlabel('Tensão')
    ylabel('Corrente')
    
    


