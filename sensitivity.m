function [df, dBAT] = sensitivity(param)

% This sensitivty function includes an assumed fixed T1prime, so we fix the
% f in T1prime to a sensible value.
fixedF = 50/6000; 

T1prime = 1./((1./param.T1t)+(fixedF./param.lamda));
M = 2*param.M0B*param.alpha*T1prime.*exp(-param.BAT./param.T1b);

%%%% Initialise %%%%
df = zeros(size(param.t));
dBAT = zeros(size(param.t));

%%%% for t between deltaT and tau plus deltaT %%%%
tRegion = param.t > param.BAT & param.t <= (param.tau+param.BAT);
t = param.t(tRegion);
% if sum(size(param.BAT)>1) > 1 % Check whether multiple values are being used
if any(size(param.BAT)>1)
    BAT = param.BAT(tRegion);
    f = param.f(tRegion);
    M_temp = M(tRegion);
else
    BAT = param.BAT;
    f = param.f;
    M_temp = M;
end
T1prime_temp = T1prime;

df(tRegion) = ...
    M_temp.*(1 - exp((BAT-t)./T1prime_temp));

dBAT(tRegion) = ...
    M_temp.*f.*((-1./param.T1b) - exp((BAT-t)./T1prime_temp).*((1./T1prime_temp)-(1./param.T1b)));

%%%% for t greater than tau plus deltaT %%%%
tRegion = param.t > (param.tau+param.BAT);
t = param.t(tRegion);
if any(size(param.tau)>1)
    tau = param.tau(tRegion);
else
    tau = param.tau;
end
% if sum(size(param.BAT)>1) > 1
if any(size(param.BAT)>1)
    BAT = param.BAT(tRegion);
    f = param.f(tRegion);
    M_temp = M(tRegion);
else
    BAT = param.BAT;
    f = param.f;
    M_temp = M;
end
T1prime_temp = T1prime;

df(tRegion) = ...
    M_temp.*exp((-t+tau+BAT)./T1prime_temp).*(1-exp(-tau./T1prime_temp));

dBAT(tRegion) = ...
    M_temp.*f.*(1-exp(-tau./T1prime_temp)).*exp((BAT+tau-t)./T1prime_temp).*((1./T1prime_temp)-(1./param.T1b));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end