function [bestPLD,bestTau,bestminVariance] = OED_PCASL_SeqmultiLD_LOptimal(varargin)

time = tic;

[param, BATDist, distWeight, scanTime, A, allSlice, nPLD, lims, slicedt] = parseInputs(varargin,nargin);
allSliceL = length(allSlice);

% Using what we know about the optimal design to limit the distribution we
% search over.
tauTry = round( ( lims.tauLB:lims.tauStep:lims.tauUB ), 5)';
tauTryL = length(tauTry);

TTStep = 10; % We are only going to calculate the FIM for this many tau at a time, to save memory

distL = length(BATDist);

%% Now minimise the variance over the distribution

bestminVariance = 1e99;
disp(['nPLD = ' ns(nPLD)])

PLD = (rand(nPLD,1) * (lims.PLDUB-lims.PLDLB)) + lims.PLDLB; % Random array between lims.PLDLB and lims.PLDUB
PLD = sort(PLD);
PLD = round(PLD ./ lims.PLDStep) .* lims.PLDStep;

tau = rand(nPLD,1) + 0.8; % Random array between 0.8 and 1.8
tau = round(tau ./ lims.tauStep) .* lims.tauStep;

continueFlag = true;
countPLDLoop = 0;
while continueFlag
    countPLDLoop = countPLDLoop + 1;
    oldPLD = PLD;
    oldTau = tau;
    
    for ii = randperm(nPLD)
     
        if ii == 1
            PLDTry = round((lims.PLDLB:lims.PLDStep:PLD(2)), 5)';
        elseif ii == nPLD
            PLDTry = round((PLD(ii-1):lims.PLDStep:lims.PLDUB), 5)';
        else
            PLDTry = round((PLD(ii-1):lims.PLDStep:PLD(ii+1)), 5)';
        end
        PLDTryL = length(PLDTry);  
     
        costMean = zeros(PLDTryL, tauTryL);
        for TT = 1:TTStep:tauTryL
            disp(['PLD loop = ' ns(countPLDLoop) ' | PLD opt = ' ns(ii) '/' ns(nPLD) ' | LD opt = ' ns(TT) '/' ns(tauTryL)])
            
            tauTryInd = TT:min( (TT+TTStep-1),  tauTryL); % Make sure not to exceed dimensions
            tauTryIndL = length(tauTryInd);
            
            % Adjust dimensions
            param.PLD = repmat(PLD, 1, PLDTryL, min(TTStep,tauTryIndL));
            param.tau = repmat(tau, 1, PLDTryL, min(TTStep,tauTryIndL));
            distWeightCurr = permute( repmat( distWeight', PLDTryL, 1, tauTryIndL ), [1,3,2] );
            
            otherInd = [1:ii-1,ii+1:nPLD];

            param.PLD(otherInd,:,:) = repmat(PLD(otherInd),1,PLDTryL,tauTryIndL);
            param.PLD(ii,:,:) = repmat(PLDTry',1,1,tauTryIndL);

            param.tau(otherInd,:,:) = repmat(tau(otherInd),1,PLDTryL,tauTryIndL);
            param.tau(ii,:,:) = permute(repmat(tauTry(tauTryInd)',1,1,PLDTryL),[1,3,2]);

            variance = zeros(PLDTryL, min(TTStep,tauTryIndL), distL, allSliceL);
            for kk = 1:allSliceL
                slice = allSlice(kk);

                param.t = param.tau + param.PLD + ((slice-1)*slicedt);
                [variance(:,:,:,kk)] = Hessian_LOptimal_analytical(param,A,scanTime,slice,slicedt);
            end
            
            variance = variance * 6000 * 6000; % Change into (ml/100g/min)^2
            
            % Take mean of generalised variance across slices
            varianceMean = mean(variance,4);
            
            if any(varianceMean(:)==0)
                warning('Some variances are zero. Setting to inf...')
                varianceMean(varianceMean==0) = inf;
            end
            
            % Take mean of generalised variance across the BAT distribution
            cost = distWeightCurr .* varianceMean;
            cost(distWeightCurr==0) = 0;              % To correct for 0*nan in distWeight .* CovCurr
            cost(isnan(cost)) = inf;                  % To correct for 0*inf in distWeight .* CovCurr
            costMean(:, tauTryInd) = mean( cost ,3 ); % Weighted mean
            
        end
        
        % Find the Tau and PLD that leads to the minimum generalised variance
        [minVariance,jj] = min(costMean(:));
        [row,col] = ind2sub(size(costMean),jj);
        
        % Save the optimal PLD and LD
        PLD(ii) = PLDTry(row);
        tau(ii) = tauTry(col);
        
        disp(['Tau = ' array2string(tau,'; ')])
        disp(['PLD = ' array2string(PLD,'; ')])
        
    end % End PLD loop

    % If the PLDs and LDs stay constant, then exit
    if sum(PLD~=oldPLD)==0 && sum(tau~=oldTau)==0
        continueFlag = false;
    end
    
end

if (bestminVariance-minVariance)/bestminVariance > 1e-12
    
    bestTau = tau;
    bestPLD = PLD;
    
    disp(['Best Tau = ' array2string(bestTau,'; ')])
    disp(['Best PLD = ' array2string(bestPLD,'; ')])
    
    bestminVariance = minVariance;
    
    param.tau = bestTau;
    param.PLD = bestPLD;
    param.t = param.tau + param.PLD;
    disp(['numAv = ' ns(TRWeightingOrNAveFloor(param,scanTime,1,1,slicedt))])
end

toc(time)

%%
%%%% Function to set up inputs %%%%
    function [param, BATDist, distWeight, scanTime, A, allSlice, nPLD, lims, slicedt] = parseInputs(varargins,nargins)
        param = varargins{1};
        BATDist = varargins{2}(:);
        param.BAT = BATDist;
        if nargins>2; distWeight = varargins{3}(:); else; distWeight = ones(size(BATDist),1); end
        if nargins>3; scanTime = varargins{4}; else; scanTime = 300; end
        if nargins>4; A = varargins{5};        else; A = [1,1;1,1]; end
        if nargins>5; allSlice = varargins{6}; else; allSlice = 1; end
        if nargins>6; nPLD = varargins{7};     else; nPLD = 6; end
        if nargins>7
            lims = varargins{8};
            %if(~isfield(lims,'maxIter'));lims.maxIter=50;end
            if(~isfield(lims,'tauStep'));lims.tauStep=0.025;end
            if(~isfield(lims,'PLDStep'));lims.PLDStep=0.025;end
            if(~isfield(lims,'tauLB'));lims.tauLB=0.1;end
            if(~isfield(lims,'tauUB'));lims.tauUB=6;end
            if(~isfield(lims,'PLDLB'));lims.PLDLB=0;end
            if(~isfield(lims,'PLDUB'));lims.PLDUB=BATDist(end)+lims.PLDStep;end % If 3s is longest BAT, we will see the peak signal with 3s PLD
        else
            lims = struct('PLDStep',0.025,'PLDLB',0,'PLDUB',BATDist(end)+0.025,...
                'tauStep',0.025,'tauLB',0.1,'tauUB',6); %'maxIter',50
        end
        if nargins>8; slicedt = varargins{9};
        else
            slicedt = 0.053125;
            warning(['OED_CASL_LOptimal_PLDTAU_acrossSlices.m: No slicedt specified. Using default slicedt = ' ns(slicedt)])
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%

end
