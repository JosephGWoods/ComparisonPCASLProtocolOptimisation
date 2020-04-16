function [bestPLD,bestTau,bestminVariance] = OED_PCASL_Hadfixed_LOptimal(varargin)

time = tic;

[param, BATDist, distWeight, scanTime, A, allSlice, nPLD, lims, slicedt] = parseInputs(varargin,nargin);
allSliceL = length(allSlice);

% Using what we know about the optimal design to limit the distribution we
% search over.
tauTry = round((lims.tauLB:lims.tauStep:lims.tauUB),5)';
tauTryL = length(tauTry);
PLDTry = round((lims.PLDLB:lims.PLDStep:lims.PLDUB),5)';
PLDTryL = length(PLDTry);

distL = length(BATDist);

%% Now minimise the variance over the distribution

bestminVariance = 1e99;
for numPLD = nPLD
    disp(['numPLD = ' ns(numPLD)])
    
    % Set the initial block length and PLD
    tau = -1; % Currently set up for fixed block duration
    PLD = -1;
    
    distWeightCurr = permute( repmat( distWeight', PLDTryL, 1, tauTryL ), [1,3,2] );
    clear distWeight
    
    continueFlag = true;
    while continueFlag
        oldTau = tau;
        oldPLD = PLD; % Each time through all of the samples, save them
        
        param.tau = permute(repmat(tauTry, 1, nPLD, PLDTryL), [2,3,1]);
        
        % For time-encoded data, there are effective PLDs
        PLD_effectiveShift = cumsum(param.tau(end:-1:1,:,:)) - param.tau(end:-1:1,:,:); % Count up backwards
        PLD_effectiveShift = PLD_effectiveShift(end:-1:1,:,:); % Order so the blocks are in order
        PLDTry_curr = PLD_effectiveShift + repmat(PLDTry(:)', nPLD, 1, tauTryL);
        clear PLD_effectiveShift
        
        for ii = 1%1:numPLD
            
            variance = zeros(PLDTryL, tauTryL, distL, allSliceL);
            for kk = 1:allSliceL
                slice = allSlice(kk);
                
                param.PLD = PLDTry_curr + ((slice-1)*slicedt);
                param.t = param.tau + param.PLD;
                [variance(:,:,:,kk)] = Hessian_LOptimal_analytical(param,A,scanTime,slice,slicedt);
            end % End slice loop
            
            % For comparison between nPLD, divide by the number of
            % images will be used for decoded: the result of signal
            % averaging on the variance.
            variance = variance / ((nPLD + 1)/2);
            
            variance = variance * 6000 * 6000; % Change into (ml/100g/min)^2
            
            % Take mean of generalised variance across the slices
            varianceMean = mean(variance,4);
            
            if any(varianceMean(:)==0)
                warning('Some variances are zero. Setting to inf...')
                varianceMean(varianceMean==0) = inf;
            end
            
            % Take mean of generalised variance across the BAT distribution
            cost = distWeightCurr .* varianceMean;
            cost(distWeightCurr==0) = 0; % To correct for 0*nan in distWeight .* CovCurr
            cost(isnan(cost)) = inf;     % To correct for 0*inf in distWeight .* CovCurr
            costMean = mean( cost , 3 ); % Weighted mean
            
            % Find the Tau and PLD that leads to the minimum generalised variance
            [minVariance,jj] = min(costMean(:));
            [row,col] = ind2sub(size(costMean),jj);
        
            % Save the optimal PLD and LD
            PLD = PLDTry(row);
            tau = tauTry(col);
            
            disp(['Tau = ' array2string(tau,'; ')])
            disp(['PLD = ' array2string(PLD,'; ')])
            
        end % End PLD loop
        
        if sum(PLD~=oldPLD)==0 && (tau~=oldTau)==0
            continueFlag = false;
        end
        
    end % End while
    
    if (bestminVariance-minVariance)/bestminVariance > 1e-12
        
        bestTau = tau;
        bestPLD = PLD;

        disp(['Best Tau = ' array2string(bestTau,'; ')])
        disp(['Best PLD = ' array2string(bestPLD,'; ')])

        bestminVariance = minVariance;
        
        param.tau = repmat(bestTau, numPLD, 1);
        PLD_effectiveShift = cumsum(param.tau(end:-1:1)) - param.tau(end:-1:1); % Count up backwards
        PLD_effectiveShift = PLD_effectiveShift(end:-1:1,:,:);
        param.PLD = PLD_effectiveShift + bestPLD;
        param.t = param.tau + param.PLD;
        
        disp(['numAv = ' ns(TRWeightingOrNAveFloor(param,scanTime,1,1,slicedt))])
    end
    
    toc(time)
end


%%
%%%% Function to set up inputs %%%%

    function [param, BATDist, distWeight, scanTime, A, allSlice, nPLD, lims, slicedt] = parseInputs(varargins,nargins)
        param = varargins{1};
        BATDist = varargins{2}(:);
        param.BAT = BATDist;
        if nargins>2; distWeight = varargins{3}(:); else; distWeight = ones(size(BATDist)); end
        if nargins>3; scanTime = varargins{4}; else; scanTime = 300; end
        if nargins>4; A = varargins{5};        else; A = [1,1;1,1]; end
        if nargins>5; allSlice = varargins{6}; else; allSlice = 1; end
        if nargins>6; nPLD = varargins{7}; else; nPLD = 7; end
        if nargins>7
            lims = varargins{8};
            %if(~isfield(lims,'maxIter'));lims.maxIter=50;end
            if(~isfield(lims,'tauStep'));lims.tauStep=0.025;end
            if(~isfield(lims,'PLDStep'));lims.PLDStep=0.025;end
            if(~isfield(lims,'tauLB'));lims.tauLB=0.1;end
            if(~isfield(lims,'tauUB'));lims.tauUB=4/nPLD(end);end
            if(~isfield(lims,'PLDLB'));lims.PLDLB=0;end
            if(~isfield(lims,'PLDUB'));lims.PLDUB=0.5;end % If 3s is longest BAT, we will see the peak signal with 3s PLD
        else
            lims = struct('PLDStep',0.025,'PLDLB',0,'PLDUB',BATDist(end)+0.025,...
                'tauStep',0.025,'tauLB',0.1,'tauUB',4/nPLD(end));
        end
        if nargins>8; slicedt = varargins{9};
        else
            slicedt = 0.053125;
            warning(['OED_CASL_PLD_NLLS_TRW_WDist_Floor_Unnormalised_acrossSlices.m: No slicedt specified. Using default slicedt = ' ns(slicedt)])
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%

end
