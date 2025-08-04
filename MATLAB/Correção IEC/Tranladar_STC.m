clc; clear; close all;
%load 'IVCurves_Ponto_Quente_reduzido.mat';
load('Desoto_demo.mat')
NomeArquivoSaida = 'Curvas_Corrig_STC_SANDIA_alfa_Beta_regress.mat';
%% %%%%%%%%%------------------------------------------------------------
%
% Salvar graficos com o formato .emf
NomeArquivoImagem = 'STC_SANDIA.emf';
%
%% %%%%%%%%-------------------------------------------------------------
%nCurves = 750;
IVCurves = IVCurves([1:754]);

%% Seleciona Curvas numa faixa de Irradiancia:
Esel_min = 400;
Esel_max = 1500;
id = NaN; k = 1;
for i = 1:length(IVCurves)
    if ((IVCurves(i).Ee >= Esel_min) & (IVCurves(i).Ee <= Esel_max))
        id(k) = i;
        k = k + 1;
    end
end
IVCurves = IVCurves(id);


%for i = 1:nCurves
%IVCurves(i).Tc = 25.0000;
%end

% Specs.aIsc = 0.054/100;			%Valor do datasheet do modulo
% Specs.bVoc = -0.343/100;		%Valor do datasheet do modulo

%Specs.bVoc = -0.42/100; % [1/K] beta supostamente corrigido

%% Modulo Sandia
%Voc0 = 21.52;    % Voc medido no stc Datasheet
%  Specs.aIsc_rel = 0.054/100;  % [relativo 1/K]
%  Specs.bVoc_rel = -0.343/100;  % [relativo 1/K]

    %  Alfa e Beta Regressão Linear
% Specs.aIsc_rel = 0.06815/100;  % [relativo 1/K]
% Specs.bVoc_rel = -0.35139/100;  % [relativo 1/K]

    %  Alfa e Beta Sandia
% Specs.aIsc_rel =  0.06006351/100;  % [relativo 1/K]
% Specs.bVoc_rel = -0.37139/100;  % [relativo 1/K]

    %  Alfa e Beta Regressão 
% Specs.aIsc_rel =  0.060199999999998/100;  % [relativo 1/K]%0.004620952000000;
% Specs.bVoc_rel = -0.348504186290205/100;;  % [relativo 1/K]%-0.074291649329415;

% Media Alfa e Beta Excel através de Voc e Isc alfa e beta médio%
Specs.aIsc_rel =   0.06020000/100;
Specs.bVoc_rel = -0.35468047/100;

    %  Alfa e Beta Média Gaussiana 5%
% Specs.aIsc_rel =  0.06112360/100;  % [relativo 1/K]
% Specs.bVoc_rel = -0.34794506/100;  % [relativo 1/K]

%    %  Alfa e Beta Sandia PVSyst
% Specs.aIsc_rel =  0.06006351/100;  % [relativo 1/K]
% Specs.bVoc_rel = -0,37914313/100;  % [relativo 1/K]


nCurves = length([IVCurves.Isc]);
for i = 1:nCurves  
plot([IVCurves(i).V],[IVCurves(i).I], '-')
hold on
xlim([0 max([IVCurves(:).Voc])+1])
ylim([0 max([IVCurves(:).Isc])*1.03])
end
title('dados de entrada');

% %Coeficiente calculados através de Alfa e Beta por Regressão Linear:
% a  = 0.0504404;		% Valor Calculado pelo IEC_60891_Proc_2_parameretros_Rs_a
% Rs = 0.2761017;		% Valor Calculado pelo IEC_60891_Proc_2_parameretros_Rs_a
% kappa = 0.0038721;  % Valor Calculado pelo IEC_60891_Parametros_kappa

%Coeficiente calculados através de Alfa e Beta por Media Gaussiana:
% a  = 0.0287813;		% Valor Calculado pelo IEC_60891_Proc_2_parameretros_Rs_a
% Rs = 0.2887619;		% Valor Calculado pelo IEC_60891_Proc_2_parameretros_Rs_a
% kappa = 0.0038059;  % Valor Calculado pelo IEC_60891_Parametros_kappa

% a  = 0.0512100;		% Valor Calculado pelo IEC_60891_Proc_2_parameretros_Rs_a
% Rs = 0.2751221;		% Valor Calculado pelo IEC_60891_Proc_2_parameretros_Rs_a
% kappa = 0.0048828;  % Valor Calculado pelo IEC_60891_Parametros_kappa

% Coeficientes calculados com alfa e beta obtidos através de regressão
% a  = 0.0778828;		% Valor Calculado pelo IEC_60891_Proc_2_parameretros_Rs_a
% Rs = 0.2668732;		% Valor Calculado pelo IEC_60891_Proc_2_parameretros_Rs_a
% kappa = 0.0036689;  % Valor Calculado pelo IEC_60891_Parametros_kappa

 %Media Alfa e Beta Excel através de Voc e Isc alfa e beta médio%
a  = 0.0647002;
Rs = 0.2779663;
kappa = 0.0021924;
 
 
figure(2)
for i=1:nCurves
    %Result(i) = transladaCurvas(IVCurves(i), IVCurves(i).Tc, 1000, Specs.aIsc_rel, Specs.bVoc_rel, Rs, a, kappa, false);
    Result(i) = transladaCurvas(IVCurves(i), 25, 1000, Specs.aIsc_rel, Specs.bVoc_rel, Rs, a, kappa, false);
    plot([Result(i).V],[Result(i).I], '.-')
    hold on
    xlim([0 max([Result(:).Voc])+1])
    ylim([0 max([Result(:).Isc])*1.03])
end
idi = NaN; ki = 1;
for i = 1:length(Result)
    if (Result(i).Isc > 1)
        idi(ki) = i;
        ki = ki + 1;
    end
end
IVCurves = Result(idi);
%IVCurves = Result;
title('Todas as curvas corrigidas para stc')

% curvas selecionadas
IVtemp1 = IVCurves(:);
% u1 = [1;0;0;1;0;0;1;0;1;0;1;0;0;1] == 1;
% IVtemp1 = IVtemp1(u1);
fig = figure;
%set(gcf, 'Position', get(0, 'Screensize'));
% linestyle = {'-','o','*','x','+', 'v','d','^','s','>','<'};
estiloLinha = {'--','-','-','-','--', '-','d','^','s','>','<'};
for i = 1:length(IVtemp1)
    k = rem(i,length(estiloLinha));
    p = plot([IVtemp1(i).V],[IVtemp1(i).I]);
    p.LineWidth = 1;
    hold on
    xlim([0 max([Result(:).Voc])+1])
    ylim([0 max([Result(:).Isc])*1.03])
   
end
xlabel('Tensão (V)', 'FontSize',12,'FontWeight','bold');
ylabel('Corrente (A)', 'FontSize',12,'FontWeight','bold');



title('Todas as curvas corrigidas para stc') 

% V = [curvaReferencia.V]; I=[curvaReferencia.I];
% [~,id] = max(V.*I);
% p = plot(V(id), I(id), 'ok');
% p.MarkerSize = 8;
% p = plot(V(end), I(end), '*k');
% p.MarkerSize = 5;

% leg1 = legend([pontoQuente,'PMP','Voc']);
% leg1.Location = 'southwest';
                 
magnifyOnFigure(fig, 'initialPositionSecondaryAxes', [330 68 90 110],...
                     'initialPositionMagnifier',[476 46 16 40], ...
                     'secondaryAxesFaceColor', [1 1 1]);
                 
magnifyOnFigure(fig, 'initialPositionSecondaryAxes', [180 240 90 110],...
                     'initialPositionMagnifier',[395 330 30 40], ...
                     'secondaryAxesFaceColor', [1 1 1]);                  

                 

                 % 
% IVtemp2 = IVCurves([15:28]);
% u2 = u1;
% IVtemp2 = IVtemp2(u2);
% figure
% for i = 1:length(IVtemp2)
%     plot([IVtemp2(i).V],[IVtemp2(i).I], '.-')
%     hold on
%     xlim([0 max([Result(:).Voc])+1])
%     ylim([0 max([Result(:).Isc])*1.03])
% end
% 
% IVtemp3 = IVCurves([29:43]);
% u3 = [0;1;0;0;1;0;0;1;0;1;0;1;0;0;1] == 1;
% IVtemp3 = IVtemp3(u3);
% figure
% for i = 1:length(IVtemp3)
%     plot([IVtemp3(i).V],[IVtemp3(i).I], '.-')
%     hold on
%     xlim([0 max([Result(:).Voc])+1])
%     ylim([0 max([Result(:).Isc])*1.03])
% end
% 
% IVflir = IVCurves([44:74]);
% u4 = [1;0;0;0;0;0;0;1;0;0;0;0;1;0;0;0;0;0;0;0;1;0;0;0;0;1;0;0;1;0;0] == 1;
% IVflir = IVflir(u4);
% figure
% for i = 1:length(IVflir)
%     plot([IVflir(i).V],[IVflir(i).I], '.-')
%     hold on
%     xlim([0 max([Result(:).Voc])+1])
%     ylim([0 max([Result(:).Isc])*1.03])
% end

% for i = 1: nCurves
%     Result(i).V = wrev(Result(i).V);
%     Result(i).I = wrev(Result(i).I);
% end
Result = IVtemp1;
save(NomeArquivoSaida,'Result', 'Specs')
% saveas(fig,NomeArquivoImagem);
disp('Arquivo salvo')

%%
function newCurves = transladaCurvas(IVCurves, Tobj, Eobj, alfa, beta, Rs, a, kappa, plotResult)
% IVCurves estrutura com as curvas a serem transladadas
% temperatura objetivo em ºC
% Irradiancia objetivo em W/m^2
% Rs resistencia série
newCurves = IVCurves;
for i = 1:length([IVCurves(:).Tc])
    I0 = [IVCurves(i).I];
    V0 = [IVCurves(i).V];
    T0 = IVCurves(i).Tc; % temperatura atual
    E0 = IVCurves(i).E;  % irradiancia atual
    [Itransl] = [I0] * (1 + alfa*(Tobj - T0)) * Eobj / E0;
    [Vtransl] = [V0] + max(V0)*(beta*(Tobj - T0) + a*log(Eobj/E0)) - Rs*(Itransl - I0) - kappa*Itransl*(Tobj - T0);
    
    % novo ponto de maxima potencia(aproximado)
    % para maior precisao, implementar uma reconstrucao da curva num dominio
    % continuo
    [~, mpp_id] = max(Vtransl .* Itransl);
    newCurves(i).Imp = Itransl(mpp_id);
    newCurves(i).Vmp = Vtransl(mpp_id);
    %---
%     obtem novos valores para Voc e Isc por extrapolacao linear
%     [Vtemp, id] = unique(Vtransl);
%     Itemp = Itransl(id);
%     [Itemp, id] = unique(Itemp);
%     Vtemp = Vtemp(id);
%     Isc = interp1(Vtemp, Itemp, 0, 'linear', 'extrap');
%     Voc = interp1(Itemp, Vtemp, 0, 'linear', 'extrap');
    
    uI = Itransl >= 0.99*max(Itransl);
    x = [ones(length(Vtransl(uI)),1), Vtransl(uI)];
    y = Itransl(uI);
    theta = x\y;
    Isc = theta(1);
    
    uV = Vtransl >= 0.99*max(Vtransl);
    x = [ones(length(Itransl(uV)),1), Itransl(uV)];
    y = Vtransl(uV);
    theta = x\y;
    Voc = theta(1);
    
    newCurves(i).I = [Isc; Itransl; 0];
    newCurves(i).V = [0; Vtransl; Voc];
    newCurves(i).Isc = Isc;
    newCurves(i).Voc = Voc;
    %----
    newCurves(i).Tc = Tobj;
    newCurves(i).E = Eobj;
    
end
if plotResult
    figure
    hold on;
    plot([newCurves(:).V], [newCurves(:).I], '*')
    xlim([0 max([newCurves(1).V])+1])
    ylim([0 max([newCurves(1).I])*1.03])
end
end

