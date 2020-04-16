function [tau, PLD_effective] = bolus_fit(param, num_blocks, maxTime, PLD, stepSize)
% Cost function for variable TI for te-pCASL

if nargin < 5
    stepSize = 0.001; % 1 ms
end

options = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
lb = zeros(1,num_blocks); lb(1) = maxTime;
ub = maxTime*ones(1,num_blocks);
tau_0 = linspace(maxTime,0,num_blocks);

tau = fmincon(@te_cost,tau_0,[],[],[],[],lb,ub,@constraint,options);

tau = round(tau ./ stepSize) .* stepSize; % Round to nearest stepSize

PLD_effective = cumsum(tau(end:-1:1)) - tau(end:-1:1) + PLD;
PLD_effective = PLD_effective(end:-1:1);

    function cost = te_cost(tau)
        
        PLD_eff = cumsum(tau(end:-1:1)) - tau(end:-1:1) + PLD;
        PLD_eff = PLD_eff(end:-1:1);
        
        total_Mz = param.T1b*exp(-PLD_eff./param.T1b).*(1 - exp(-tau./param.T1b));
        cost = 1E6 * sum(abs(diff(total_Mz)));
    end

    function [c,ceq] = constraint(tau)
        
        % Compute % nonlinear inequalities at x.
        c = diff(tau); % The label durations must be decreasing
        
        % Compute nonlinear equalities at x.
        ceq = 0;
        
    end

end
