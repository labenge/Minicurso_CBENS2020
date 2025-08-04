load('Desoto_demo.mat')
Vmed = IVCurves.V;
Imed = IVCurves.I;
Voc  = IVCurves(1).Voc*ones(size(Imed));
Suns = IVCurves(1).E*ones(size(Imed));
Tmed = IVCurves(1).Tc*ones(size(Imed));
for k = 2:ncurves
    Vmed = [Vmed; IVCurves(k).V];
    Imed = [Imed; IVCurves(k).I];
    Voc  = [Voc; IVCurves(k).Voc*ones(size(IVCurves(k).I))];
    Suns = [Suns; IVCurves(k).E *ones(size(IVCurves(k).I))];
    Tmed = [Tmed; IVCurves(k).Tc*ones(size(IVCurves(k).I))];
end
Itransl = Imed .* (1 + a*(25 - Tmed)) ./ Suns * 1000;
Vtransl = Vmed + Voc.*(b*(25 - Tmed) + aa*log(1./Suns)) - Rs*(Itransl - Imed) - kappa*Itransl.*(25 - Tmed);
%%%%%%%%%%%%%%%%%%% TENHO DÚVIDAS QUANTO A ESSE VOC VARIANDO... 
%%%%%%%%%%%%%%%%%%% NAO TERIA QUE SER UM VALOR UNICO?
k = 1.3806503e-23;    % Boltzmann [J/K]
q = 1.60217646e-19;   % Electron charge [C]
Vt = Specs.Ns*k*(Tmed+273.15)/q;
%%%%%%%%%%%% Sem coeficientes de temperatura a, b e Suns
RMSE = (sum((Itransl - Iph + I0.*(exp((Vtransl+Itransl.*Rs)./(n.*Vt))-1) + ...
    (Vtransl+Itransl.*Rs)./Rp).^2)/length(Imed))^.5;

    fprintf('\n Iph(A) = %.3e',Iph);
    fprintf('\n I0(uA) = %.3e . exp (%.3e . T)',I0,I0T);
    fprintf('\n n = %.3e + %.3e . T',n,nT);
    fprintf('\n Rs = %.3e + %.3e . T',Rs0,RsT);
    fprintf('\n Rp = %.3e + %.3e . T',Rp0,RpT);
    fprintf('\n alpha = %.3e',a);
    fprintf('\n beta = %.3e',b);
    fprintf('\n kappa = %.3e',kappa);
    fprintf('\n a = %.3e',aa);
    fprintf('\n RMSE = %.5e \n',RMSE);