function [Iph, I0, I0T, n, nT, Rs0, RsT, Rp0, RpT, a, RMSE] = Jaya_PV_Mod(Vmed, Imed, Suns, Ns, T, lim_inf, lim_sup, runs, pop_size, max_fes)
%%
% Descricao:
%   Jaya_PV_Original computa os cinco parametros para o modelo de diodo unico
%   utilizando o algoritmo de Jaya em sua versao original.
%
% Entradas:
%   Vmed - vetor de tensao medida [V]
%   Imed - vetor de corrente medida [A]
%   Ns   - numero de celulas em serie
%   T    - temperatura do modulo  [°C]
%   lim_inf - vetor dos limites inferiores dos parametros a
%      serem estimados. Formato: [Iph_min, I0_min, n_min, Rs_min, Rp_min];
%   lim_sup - vetor dos limites superiores dos parametros a
%      serem estimados. Formato: [Iph_max, I0_max, n_max, Rs_max, Rp_max];
%   runs - quantidade de rodadas
%   pop_size - tamanho da populacao de cada rodada
%   max_fes - quantidade maxima de avaliacao da funcao objetivo
%
% Saidas:
%   Iph - corrente fotogerada [A]
%   I0 - corrente de saturacao reversa [A]
%   n - indice de idealidade do diodo
%   Rs - resistencia serie [Ohms]
%   Rp - resistencia paralelo [Ohms]
%   RMSE - erro quadradico medio
%% Constantes
T = T + 273.15;       % Tempeture [K]
% showEvolution = false; % deseja visualizar os graficos?
%%
% Default values for runs, pop_size and max_fes.
if isnan(runs)
    runs = 30;
end
if isnan(pop_size)
    pop_size = 20;
end
if isnan(max_fes)
    max_fes = 50000;
end
%%
maxGen = floor(max_fes/pop_size); % quantidade maxima de geracoes
var = length(lim_inf); % var é a quantidade de variaves.
%var = [Iph, I0, n, Rs, Rp];
mini = lim_inf;       % limite inferior
maxi = lim_sup;       % limite superior
for runsCounter = 1:runs
    % inicializacao da populacao
    x = mini.*ones(pop_size,var) + (maxi-mini).*ones(pop_size,var).*rand(pop_size,var);
    f = myobj(x, Vmed, Imed, T, pop_size, Ns, Suns); % fitness inicial
    for genCounter = 1:maxGen
        % Cria nova populacao
        [~, id] = sort(f); %    id_best = id(1);                            id_worst = id(end);
        xnew = x + rand(pop_size,var).*(x(id(1), :) - abs(x)) - rand(pop_size,var).*(x(id(end), :) - abs(x));
        
        % Mantem parametros dentro dos limites
        for i = 1:var
            xnew(xnew(:,i) < mini(i), i) = mini(i);
            xnew(xnew(:,i) > maxi(i), i) = maxi(i);
        end
        
        % Mantem somente individuos melhores
        fnew = myobj(xnew, Vmed, Imed, T, pop_size, Ns, Suns);
        u = fnew < f;
        % A populacao nova da iteracao corrente é a antiga da proxima geracao
        f(u) = fnew(u);
        x(u, 1:var) = xnew(u, 1:var);
        
        % Verificar se erro quadratico estabilizou
        
    end
    [fmin, id] = min(f);
    xmelhores(runsCounter,:) =  [x(id,:), fmin]; % Vetor que armazena o melhor de cada rodada
    %% Plotar graficos a cada fim de geração
    %    disp('Parametros obtidos no ciclo i')
    %    Iph = x(id,1)
    %    I0 = x(id,2)
    %    n = x(id,3)
    %    Rs = x(id,4)
    %    Rp = x(id,5)
    %    figure(runs)
    %    plot(Vmed, Imed,'o');
    %    hold on
    %    fmodelo = @(V,I) (I - Iph*ones(size(I)) + I0*(exp((V+I*Rs)/(n*Vt))-1) + (V+I*Rs)/Rp);
    % fimplicit(fmodelo,[0 Vmed(end) 0 Imed(1)]);
    %    ezplot(fmodelo,[0 20 0 7.8]); % Mudar
    %temp = input('digite algo');
    %    xlim([0 Vmed(end)+1])
    %    ylim([0, Imed(1) + 0.2])
    %%
    %     [val,ind] = min(fopt);   %
    %     Fes(runs) = pop_size*ind;     % number of function evaluations
    %     best(runs) = val;
end
% valores medios dos parametros
%disp('Parametros medios')
%Iph_med = mean(xmelhores(:,1))
%I0_med = mean(xmelhores(:,2))
%n_med = mean(xmelhores(:,3))
%Rs_med = mean(xmelhores(:,4))
%Rp_med = mean(xmelhores(:,5))
%f_med = mean(xmelhores(:,6))

%% Melhor resultado
[~, id] = min(xmelhores(:,11));
Iph = xmelhores(id,1);
I0  = xmelhores(id,2);
I0T = xmelhores(id,3);
n   = xmelhores(id,4);
nT  = xmelhores(id,5);
Rs0 = xmelhores(id,6);
RsT = xmelhores(id,7);
Rp0 = xmelhores(id,8);
RpT = xmelhores(id,9);
a = xmelhores(id,10);
RMSE = sqrt(xmelhores(id,11));
%disp('\n Melhor resultado')
%disp('[Iph, I0, n, Rs, Rp, fitness]');
%disp(xmelhores(id,:))
%fprintf('\n Iph(A) = %f',xmelhores(id,1));
%fprintf('\n I0(uA) = %f',xmelhores(id,2)*10^6);
%fprintf('\n n = %f',xmelhores(id,3));
%fprintf('\n Rs(ohms) = %f',xmelhores(id,4));
%fprintf('\n Rph(ohms) = %f',xmelhores(id,5));
%fprintf('\n fmin = %f',xmelhores(id,6));



%% Mostra resultados em termos de desvio padrao, media, minimo e max fitness
%bbest = min(best);
%mbest = mean(best);
%wbest = max(best);
%stdbest = std(best);
%mFes = mean(Fes);
%fprintf('\n\n best = %f',bbest);
%fprintf('\n mean = %f',mbest);
%fprintf('\n worst = %f',wbest);
%fprintf('\n std. dev. = %f',stdbest);
%fprintf('\n mean Fes = %f',mFes);
%fprintf('\n');
end
%%
function [f] = myobj(x, Vmed, Imed, T, pop_size, Ns, Suns)
% sqrt foi removido para melhorar performace
MSE = zeros(1,pop_size); % pre-alocacao de memoria
l = length(Imed);
for i = 1:pop_size
    a = x(i,10);
    Iph = Suns.*(x(i,1) * (1 + a*(T - 298.15))); %x(i,1);
    I0 = x(i,2) .* exp( T * x(i,3));
    n  = x(i,4) + T * x(i,5);
    Rs = x(i,6) + T * x(i,7);
    Rs(Rs<0)=0;
    Rp = x(i,8) + T * x(i,9);
    Rp(Rp<0)=0;
    k = 1.3806503e-23;    % Boltzmann [J/K]
    q = 1.60217646e-19;   % Electron charge [C]
    Vt = Ns*k*T/q;
    %%%%%%%%%%%% Sem coeficientes de temperatura a, b e Suns
    MSE(i) = sum((Imed - Iph + I0.*(exp((Vmed+Imed.*Rs)./(n.*Vt))-1) + ...
    (Vmed+Imed.*Rs)./Rp).^2)/l;
%     DEN(i) = sum((Imed - mean(Imed)).^2)/l;
    %%%%%%%%%%%% COM
%     MSE(i) = sum((Imed - Iph + I0*(exp((Vmed+Imed*Rs)./(n.*Vt))-1) + 
%         (Vmed+Imed*Rs)/Rp).^2)/l;
end
f = MSE; % MSE./DEN; %
end