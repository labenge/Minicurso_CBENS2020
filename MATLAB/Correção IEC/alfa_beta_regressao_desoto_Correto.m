clearvars
clc; clear; close all
%load('IEC.mat')
load 'Desoto_demo.mat'
%% Voc0 = 21.53; Usar Voc medido nas condicoes stc ou Voc obtido atraves da
%% regressao???
%% Seleciona curvas com irradiancia próximas a 1000W/m^2
curvas = IVCurves;
IVCurves = IVCurves(1:754);
k = 1;
for i = 1:length([IVCurves(:).Tc])
    %if(IVCurves(i).Ee<=1001.5 & IVCurves(i).Ee>=998.5 & IVCurves(i).Tc<=65)
    if(IVCurves(i).Ee<=1010 & IVCurves(i).Ee>=990 & IVCurves(i).Tc<=85)
        id(k) = i;
        k = k + 1;    
    end 
end
curvas2 = IVCurves(id);
T0 = 25; %temperatura stc
E0 = 1000;
E = [IVCurves(id).E];
Temp = [IVCurves(id).Tc];
Voc = [IVCurves(id).Voc];
Isc = [IVCurves(id).Isc];
k = 1.38066E-23; % J/K
q = 1.60218E-19; % c
Vth = [Specs.Ns*k*(Temp + 273.15)/q];
n = 1.2; % for cSi cell

plot(Temp - T0, Voc - n*Vth.*log(E/E0),'*')
hold on
title('Voc vs Temp')
%% Regressao para obter beta e Voc0
x = [ones(length(Temp),1), [Temp.' - T0]];
y = [Voc - n*Vth.*log(E/E0)].';

coeficiente = x\y;
Voc_regress = coeficiente(1);
beta = coeficiente(2);
beta_relativo = coeficiente(2)/Voc_regress*100;
fprintf('Resultados da regressao:\n');
fprintf('Voc_stc:       %.5f   V\n', Voc_regress);
fprintf('beta:          %.5f   V/ºC\n', beta);
fprintf('beta relativo: %.5f   %%/ºC\n\n',beta_relativo) 


 disp('regressao robusta')
    X =  [Temp.' - T0];
    Y = Voc - n*Vth.*log(E/E0);
    theta2 = robustfit(X,Y.');
    Voc0_robust = theta2(1)
    beta_robust = theta2(2)
    beta_rel = beta_robust/Voc0_robust*100
temp_regress = [Temp.' - T0];
Voc_calc  = Voc0_robust + beta_robust*temp_regress;
plot(temp_regress, Voc_calc, 'r-')

% Plot do resultado da regressao
% temp_regress = [min(Temp)-1:0.1:max(Temp)+1] - T0;
temp_regress = [Temp.' - T0];
Voc_calc1  = Voc_regress + beta*temp_regress;
plot(temp_regress, Voc_calc1, 'c-')


%% Regressao para obter alfa e Isc0
figure(2)
plot(Temp - T0, Isc,'*')
hold on
title('Isc vs Temp')
xi = [ones(length(Temp),1), [Temp.' - T0]];
yi = [Isc.*(E0./E)];

coeficiente = xi\yi.';
Isc_regress = coeficiente(1);
alfa = coeficiente(2);
alfa_relativo = coeficiente(2)/Isc_regress*100;
fprintf('Resultados da regressao:\n');
fprintf('Isc_stc:       %.5f   A\n', Isc_regress);
fprintf('alfa:          %.5f   A/ºC\n', alfa);
fprintf('alfa relativo: %.5f   %%/ºC\n',alfa_relativo) 

% Plot do resultado da regressao
% temp_regress = [min(Temp)-1:0.1:max(Temp)+1] - T0;
temp_regress = [Temp.' - T0];
Isc_calc  = Isc_regress + alfa*temp_regress;
plot(temp_regress, Isc_calc, 'r-')

disp('regressao robusta')
    Xi =  [Temp.' - T0];
    Yi = [Isc.*(E0./E)];
    theta3 = robustfit(Xi,Yi.');
    Isc_robust = theta3(1)
    alfa_robust = theta3(2)
    alfa_rel = alfa_robust/Isc_robust*100
    temp_regress = [Temp.' - T0];
    Isc_calc1  = Isc_regress + alfa_robust*temp_regress;
    plot(temp_regress, Isc_calc1, 'c-')