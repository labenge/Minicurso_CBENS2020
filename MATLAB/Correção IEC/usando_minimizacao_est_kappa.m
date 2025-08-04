clearvars
clc; clear; close all
% Voc0 = 21.53; % Voc medido no stc
% aIsc = 0.054/100;  % [1/K]
% bVoc = -0.343/100; % [1/K]
%bVoc = -0.42/100; % [1/K] beta supostamente corrigido

%% Modulo Sandia
%Voc0 = 21.52;    % Voc medido no stc Datasheet
load('Desoto_demo.mat')

% aIsc = 0.054/100;  % [relativo 1/K]
% bVoc = -0.343/100;  % [relativo 1/K]

% aIsc = 0.06815/100;
% bVoc = -0.35139/100;

% Media Alfa e Beta Gaussiana 5%
% aIsc =  0.06112360/100;
% bVoc =-0.34794506/100; 

% Alfa e Beta Regressão Robusta (991 a 1010 e T<65º)Curvas (1:754) 1 dia
% bVoc = - 0.348504186290205/100;    %-0.074291649329415;
% aIsc = 0.060199999999998/100;      %0.004620952000000;

% Media Alfa e Beta Excel através de Voc e Isc alfa e beta médio%
aIsc =   0,06020000/100;
bVoc =-0.35468047/100; 

Ns = 36;
Np = 1;
%load('IVCurves_IEC_reduzido.mat')
%curvas = IVCurves([2922, 147, 806]); % CURVAS com IRRADIANCIAS constante e DIFERENTES TEMPERATURAS

%% Seleciona curvas com Irradiancias constantes Próximo de 1000W/m²
k = 1;
Tsel = 25;



Esel = 999.5; %Bom resultado
Tol = 0.5;

% Esel = 999.5; %Bom resultado
% Tol = 0.5;
id = NaN;
for i = 1:length([IVCurves(:).Tc])
    if((IVCurves(i).Ee >= Esel-Tol) & (IVCurves(i).Ee <= Esel+Tol) & (IVCurves(i).Tc >= Tsel))
        id(k) = i;
        k = k + 1;    
    end 
end
if(isnan(id))
    disp('nenhuma curva encontrada');
else
    curvas = IVCurves([id]);
end
%% identificar menor Temperatura
[~, idmenor] = min([curvas(:).Tc]);
for i = 1:length(curvas)
    plot([curvas(i).V],[curvas(i).I], '-')
    xlim([0 curvas(1).Voc+1])
    ylim([0 curvas(idmenor).Isc*1.03])
    hold on
end
Tobj = curvas(idmenor).Tc;
Eobj = curvas(idmenor).E;


%% Determinacao de Kappa
Rs = 0; a = 0;
kappa_objfun = @(kappa) max_P_diff_for_Rs(curvas, Tobj, Eobj, aIsc, bVoc, Rs, a, kappa);
options = optimset('TolX', 1e-6);
kappa = fminsearch(kappa_objfun, 0, options);
fprintf('kappa = %.7f\n', kappa);
erro_kappa = kappa_objfun(kappa);
fprintf('erro_kappa  = %.7f\n', erro_kappa);
%% Resultado
IVtransl = transladaCurvas(curvas, Tobj, Eobj, aIsc, bVoc, 0, 0, kappa, 1);
%%
function IVCurves = transladaCurvas(IVCurves, Tobj, Eobj, alfa, beta, Rs, a, kappa, plotResult)
% IVCurves estrutura com as curvas a serem transladadas
% temperatura objetivo em ºC
% Irradiancia objetivo em W/m^2
% Rs resistencia série
for i = 1:length([IVCurves(:).Tc])
    I0 = [IVCurves(i).I];
    V0 = [IVCurves(i).V];
    T0 = IVCurves(i).Tc; % temperatura atual
    E0 = IVCurves(i).E;  % irradiancia atual
    [Itransl] = [I0] * (1 + alfa*(Tobj - T0)) * Eobj / E0;
    [Vtransl] = [V0] + max(V0)*(beta*(Tobj - T0) + a*log(Eobj/E0)) - Rs*(Itransl - I0) - kappa*Itransl*(Tobj - T0);
    IVCurves(i).Itransl = Itransl;
    IVCurves(i).Vtransl = Vtransl;
end
if plotResult
    figure(2)
    for i = 1:length(IVCurves)
        plot([IVCurves(i).Vtransl],[IVCurves(i).Itransl], '-')
        xlim([0 IVCurves(1).Voc+1])
        ylim([0 IVCurves(1).Isc*1.03])
        hold on
    end
    
end
end

%% correcao para Ponto de Máxima potência
function val = max_P_diff_for_Rs(IVCurves, Tobj, Eobj, alfa, beta, Rs, a,kappa)
    IVCurves = transladaCurvas(IVCurves, Tobj, Eobj, alfa, beta, Rs, a, kappa, 0);
    for i = 1:length([IVCurves(:).Tc])
        Pmp(i) = max([IVCurves(i).Vtransl] .* [IVCurves(i).Itransl]);
    end
    val = max(abs((Pmp - median(Pmp))./median(Pmp)));
end

