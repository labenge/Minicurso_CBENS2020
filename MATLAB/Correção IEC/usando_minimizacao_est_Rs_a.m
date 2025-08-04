clc; clear; close all

%% Modulo Sandia
%Voc0 = 21.53; % Voc medido no stc



%%
Ns = 36;
Np = 1;

load('Desoto_demo.mat')

% aIsc = 0.054/100;  % [relativo 1/K]
% bVoc = -0.343/100;  % [relativo 1/K]

aIsc = 0.06815/100;
bVoc = -0.35139/100;


[~,idx] = sort([IVCurves(:).Tc]);
IVCurves = IVCurves(idx);
%curvas = IVCurves([2406, 2781, 146]);  % CURVAS com TEMPERATURAS constante e DIFERENTES IRRADIANCIAS
%curvas = IVCurves([2900, 2296, 1514]);
%% Seleciona curvas com Temperaturas constantes Proximo a Tsel

% Tsel = 49.7;
% Esel = 800;
% Tol = 0.5;
Tsel =  40;%+4
Esel = 700;
Tol = 1;

id = NaN; k = 1;
for i = 1:length(IVCurves)
    if((IVCurves(i).Tc >= Tsel-Tol) & (IVCurves(i).Tc <= Tsel+Tol) & (IVCurves(i).Ee >= Esel))
        id(k) = i;  
        k = k + 1;
    end
    IVCurves(i).legenda = [num2str(i), '-    ' num2str(IVCurves(i).Tc,2), 'ºC   ', num2str(IVCurves(i).E,4) 'W/m^2 ',];
end
if(isnan(id))
    disp('nenhuma curva encontrada');
else
     curvas = IVCurves(id);
     
end
%% identificar maior irradiancia
[~,idx] = sort([curvas(:).E]);
curvas = curvas(idx);
[~, idmaior] = max([curvas(:).E]);
figure(1)
for i = 1:length(curvas)
    plot([curvas(i).V],[curvas(i).I], '-')
    xlim([0 curvas(1).Voc+1])
    ylim([0 curvas(idmaior).Isc*1.03])
    hold on
end
leng1 = legend({curvas(:).legenda});
leng1.Location = 'northeastoutside';
numColumns = 1;
if length(curvas)> 45 
    numColumns = ceil(length(curvas)/45)
end
leng1.NumColumns = numColumns;
title('Sem correção - dados brutos');

Tobj = curvas(idmaior).Tc;
Eobj = curvas(idmaior).E;
%% Translacionar curvas
IVtransl = transladaCurvas(curvas, Tobj, Eobj, aIsc, bVoc, 0, 0, 0, 1);
%% Determinacao de a
a_objfun = @(a) max_Voc_diff_for_a(IVtransl, Tobj, Eobj, bVoc, a);
options = optimset('TolX', 1e-6);
a = fminsearch(a_objfun, 0, options);
fprintf('a  = %.7f\n', a);
erro_a = a_objfun(a);
fprintf('erro_a  = %.7f\n', erro_a);
%% Determinacao de Rs
Rs_objfun = @(Rss) max_P_diff_for_Rs(IVtransl, Tobj, Eobj, aIsc, bVoc, Rss, a);
options = optimset('TolX', 1e-6);
Rs = fminsearch(Rs_objfun, 1, options);
fprintf('Rs = %.7f\n', Rs);
erroRS = Rs_objfun(Rs);
fprintf('erro_RS  = %.7f\n', erroRS);
%%
IVtransl = transladaCurvas(curvas, Tobj, Eobj, aIsc, bVoc, Rs, a, 0, 1);

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
   figure
    for i = 1:length([IVCurves(:).Tc])
        plot([IVCurves(i).Vtransl], [IVCurves(i).Itransl], '-')
        hold on
        xlim([0 max([IVCurves(1).Vtransl])+1])
        ylim([0 max([IVCurves(1).Itransl])*1.03])
    end
end
end

%% correcao para tensao de circuito aberto
function val = max_Voc_diff_for_a(IVCurves, Tobj, Eobj, beta, a)
IVtransl =  IVCurves;
for i = 1:length([IVCurves(:).Tc])
    I0 = IVCurves(i).I;
    V0 = IVCurves(i).V;
    T0 = IVCurves(i).Tc; % temperatura atual
    E0 = IVCurves(i).E;  % irradiancia atual
    [Vtransl] = [V0] + max(V0)*(beta*(Tobj - T0) + a*log(Eobj/E0)); %%- Rs*(Itransl - I0); %%- kappa*Itransl*(Tobj - T0);
    IVtransl(i).V = Vtransl;
end
for i = 1:length([IVCurves(:).Tc])
    
    Voc(i) = max([IVtransl(i).V]);
end
val = max(abs((Voc - median(Voc))./median(Voc)));
end

%% correcao para Ponto de Máxima potência
function val = max_P_diff_for_Rs(IVCurves, Tobj, Eobj, alfa, beta, Rs, a)
    IVCurves = transladaCurvas(IVCurves, Tobj, Eobj, alfa, beta, Rs, a, 0, 0);
    for i = 1:length([IVCurves(:).Tc])
        Pmp(i) = max([IVCurves(i).Vtransl] .* [IVCurves(i).Itransl]);
    end
    val = max(abs((Pmp - median(Pmp))./median(Pmp)));
end

