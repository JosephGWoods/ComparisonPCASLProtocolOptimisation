function [bestPLD,bestTau,bestminVariance] = OED_PCASL_Hadfreelunch_fixed_LOptimal(varargin)

time = tic;

% nPLD is the same as param.num_enc for time-encoded protocols
[param, BATDist, distWeight, scanTime, A, allSlice, nPLD, lims, slicedt] = parseInputs(varargin,nargin);
allSliceL = length(allSlice);

% Using what we know about the optimal design to limit the distribution we
% search over.
tauTry = round((lims.tauLB:lims.tauStep:lims.tauUB),5)'; % The duration of the 2nd encoded LD
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
    
    param.tau = permute(repmat(tauTry, 1, numPLD, PLDTryL), [2,3,1]); % Dimensions nPLD x PLDTryL x tauTryL
    param.tau(1,:,:) = 1.8; % White Paper tag duration (Should make this a user defined parameter)
    param.PLD = repmat(PLDTry(:)', 1, 1, tauTryL);
    [~,param.PLD] = generateTETimings(param.tau , param.PLD, numPLD, 'standard');
    
    distWeight = permute( repmat( distWeight', PLDTryL, 1, tauTryL ), [1,3,2] );
    
    continueFlag = true;
    while continueFlag
        oldTau = tau;
        oldPLD = PLD; % Each time through all of the samples, save them

        for ii = 1 % Don't need to loop through timepoints for time-encoded
            
            variance = zeros(PLDTryL, tauTryL, distL, allSliceL);
            for kk = 1:allSliceL
                slice = allSlice(kk);
                
                param.t = param.tau + param.PLD + ((slice-1)*slicedt);
                [variance(:,:,:,kk)] = Hessian_LOptimal_analytical(param,A,scanTime,slice,slicedt);
            end % End slice loop
            
            % For comparison between nPLD, divide by the number of
            % images will be used for decoded: the result of signal
            % averaging on the variance.
            variance = variance / ((numPLD + 1)/2);
            
            variance = variance * 6000 * 6000; % Change into (ml/100g/min)^2
            
            % Take mean of generalised variance across the slices
            varianceMean = mean(variance,4);
            
            if any(varianceMean(:)==0)
                warning('Some variances are zero. Setting to inf...')
                varianceMean(varianceMean==0) = inf;
            end
            
            % Take mean of generalised variance across the BAT distribution
            cost = distWeight .* varianceMean;  % Weight by the ATT distribution
            cost(distWeight==0) = 0;            % To correct for 0*nan in distWeight .* varianceMean
            cost(isnan(cost)) = inf;            % To correct for 0*inf in distWeight .* varianceMean
            costMean = mean( cost , 3 );        % Weighted mean
            
            % Find the Tau and PLD that leads to the minimum generalised variance
            [minVariance,jj] = min(costMean(:));
            [row,col] = ind2sub(size(costMean),jj);
        
            % Save the optimal PLD and LD
            PLD = PLDTry(row);
            tau = tauTry(col);
            
            disp(['Tau = ' array2string(tau,'; ')])
            disp(['PLD = ' array2string(PLD,'; ')])
            
        end % End PLD loop
        
        % If the PLDs and LDs stay constant, then exit
        if sum(PLD~=oldPLD)==0 && (tau~=oldTau)==0
            continueFlag = false;
        end
        
    end % End while
    
    %if (bestminEffLS-minEffLS)/bestminEffLS > 1e-12, update
    if (bestminVariance-minVariance)/bestminVariance > 1e-12
        
        bestTau = param.tau(:,row,col);
        bestPLD = PLD;

        disp(['Best Tau = ' array2string(bestTau,'; ')])
        disp(['Best PLD = ' array2string(bestPLD,'; ')])

        bestminVariance = minVariance;
        
        param.tau = bestTau;
        [~,param.PLD] = generateTETimings(param.tau , bestPLD, numPLD, 'standard');
        param.t = param.tau + param.PLD;
        [numAv,TotalTR] = TRWeightingOrNAveFloor(param,scanTime,1,1,slicedt);
        disp(['numAv = ' ns(numAv)])
        disp(['Scan time = ' ns(numAv*TotalTR)])
        
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
        if nargins>4; A = varargins{5};        else; A = [1,0;0,0]; end
        if nargins>5; allSlice = varargins{6}; else; allSlice = 1; end
        if nargins>6; nPLD = varargins{7}; else; nPLD = 7; end
        if nargins>7
            lims = varargins{8};
            %if(~isfield(lims,'maxIter'));lims.maxIter=50;end
            if(~isfield(lims,'tauStep'));lims.tauStep=0.025;end
            if(~isfield(lims,'PLDStep'));lims.PLDStep=0.025;end
            if(~isfield(lims,'tauLB'));lims.tauLB=0.1;end
            if(~isfield(lims,'tauUB'));lims.tauUB=1.8;end
            if(~isfield(lims,'PLDLB'));lims.PLDLB=0;end
            if(~isfield(lims,'PLDUB'));lims.PLDUB=1;end
        else
            lims = struct('PLDStep',0.025,'PLDLB',0,'PLDUB',1,...
                'tauStep',0.025,'tauLB',0.1,'tauUB',1.8);
        end
        if nargins>8; slicedt = varargins{9};
        else
            slicedt = 0.053125;
            warning(['OED_CASL_PLD_NLLS_TRW_WDist_Floor_Unnormalised_acrossSlices.m: No slicedt specified. Using default slicedt = ' ns(slicedt)])
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%

end
