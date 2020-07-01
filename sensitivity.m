%% (P)CASL Buxton model CBF an ATT sensitivity functions
%
% [df, dBAT] = sensitivity(param)
%
% in:
%      param - struct of Buxton model parameters
%
% out:
%      df   - CBF sensitivity function
%      dBAT - ATT sensitivity function (BAT = Bolus Arrival Time)
%
% param:
%      .tau    - (τ) labelling duration
%      .t      - time of sampling from start of labelling (tau + PLD)
%      .f      - cerebral blood flow
%      .BAT    - (Δt) bolus arrival time
%      .T1b    - T1 of blood
%      .T1t    - T1 of tissue
%      .alpha  - (α) labelling efficiency (ratio)
%      .lambda - (λ) equilibrium tissue/blood partition coefficient of water
%                (ratio)
%      .M0B    - equilibrium magnetization of arterial blood (imaging units)
%
% units:
%      - Any unit of time can be used, but it should be consistent. If the
%        unit of time is, say, seconds, f should be in units of s^-1.
%
% dimensions:
%      - All parameters can be scalar or matrices. Any parameters that are
%        matrices must be of the same size.
%      - t should always be the same size as any non-scalar parameter, even
%        if it only has one value, because indices are derived from param.t.
%
% References: 
%             Woods et al., A General Framework for Optimizing Arterial
%             Spin Labeling MRI Experiments, MRM 2019
%
% Written by Joseph G. Woods, FMRIB, Oxford, 2016

function [df, dBAT] = sensitivity(param)

fixedF = 50/6000; % outflow in T1' is fixed to 50 mL/100g/min

% The apparent relaxation rate (includes outflow)
T1prime = 1 ./ ( (1./param.T1t) + (fixedF./param.lambda) );

% Signal scaling common to all parts of the sensitivity curves
c = 2 * param.M0B .* T1prime .* param.alpha .* exp( -param.BAT ./ param.T1b );

% If c is scalar, make it robust to the use of indices below
if numel(c) > 1; c = @(ind) c(ind);
else;            c = @(ind) c     ; end

% Initialise the sensitivity arrays
df   = zeros(size(param.t));
dBAT = zeros(size(param.t));

% Δt < t ≤ τ + Δt
tInd       = param.t > param.BAT & param.t <= (param.tau+param.BAT);
a          = exp( -(param.t-param.BAT) ./ T1prime );
qss        = 1 - a;
b          = param.f .* ( -1./param.T1b - a.*(1./T1prime-1./param.T1b) );
if numel(b) > 1; b = b(tInd); end
df(tInd)   = c(tInd) .* qss(tInd);
dBAT(tInd) = c(tInd) .* b;

% τ + Δt < t
tInd       = param.t > (param.tau+param.BAT);
a          = exp( -(param.t-param.tau-param.BAT) ./ T1prime );
b          = param.f .* (1./T1prime-1./param.T1b);
qss        = 1 - exp( -param.tau ./ T1prime );
if numel(b)   > 1; b   = b(tInd)  ; end
if numel(qss) > 1; qss = qss(tInd); end
df(tInd)   = c(tInd) .* a(tInd) .* qss;
dBAT(tInd) = c(tInd) .* a(tInd) .* b .* qss;

end