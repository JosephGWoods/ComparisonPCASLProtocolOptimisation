function [bestPLD,bestTau,bestminVariance] = OED_PCASL_Hybridvariable_LOptimal(varargin)

time = tic;

[param, BATDist, distWeight, scanTime, A, allSlice, nPLD, lims, slicedt] = parseInputs(varargin,nargin);
allSliceL = length(allSlice);

% Using what we know about the optimal design to limit the distribution we
% search over.
PLDTry = round((lims.PLDLB:lims.PLDStep:lims.PLDUB),5)';
PLDTryL = length(PLDTry);

distL = length(BATDist);

TTStep = 10; % We are only going to calculate the FIM for this many tau at a time, to save memory

%% Now minimise the variance over the distribution

bestminVariance = 1e99;
numPLD = nPLD;
disp(['numPLD = ' ns(numPLD)])

% Set the initial block length and PLD
PLD = repmat(lims.PLDLB, param.multiPLD, 1);

tau = rand(numPLD,1) * lims.tauUB; % Spread the LDs between 0 s and the LD UB
tau(tau<lims.tauLB) = lims.tauLB; % Correct short labels
tau = round(tau ./ lims.tauStep) .* lims.tauStep; % Round to nearest stepSize
tau = sort(tau, 'descend');
tau = repmat(tau, param.multiPLD, 1);

continueFlag = true;
countLoop = 0;
while continueFlag
    countLoop = countLoop + 1;
    
    oldTau = tau;
    oldPLD = PLD; % Each time through all of the samples, save them
    
    for ii = 1:param.multiPLD
        
        ind = (ii-1)*param.num_enc+1:(ii-1)*param.num_enc+param.num_enc;
        
        for ll = randperm(numPLD)

            if ll==1
                tauTry = round((tau(ind(2)):lims.tauStep:lims.tauUB),5)';
            elseif ll == numPLD
                tauTry = round((lims.tauLB:lims.tauStep:tau(ind(numPLD-1))),5)';
            else
                tauTry = round((tau(ind(ll+1)):lims.tauStep:tau(ind(ll-1))),5)';
            end
            tauTryL = length(tauTry);
            
            costMean = zeros(PLDTryL, tauTryL);
            for TT = 1:TTStep:tauTryL
                disp(['Loop = ' ns(countLoop) ' | multi-PLD = ' ns(ii) '/' ns(param.multiPLD) ' | block = ' ns(ll) '/' ns(numPLD) ' | LD = ' ns(TT) '/' ns(tauTryL)])
                
                tauTryInd = TT:min( (TT+TTStep-1),  tauTryL); % Make sure not to exceed dimensions
                tauTryIndL = length(tauTryInd);
                
                distWeightCurr = permute( repmat( distWeight', PLDTryL, 1, min(TTStep,tauTryIndL) ), [1,3,2] );
                
                param.tau         = repmat(tau, 1, PLDTryL, min(TTStep,tauTryIndL));
                param.tau(ind(ll),:,:) = repmat(tauTry(tauTryInd)' , PLDTryL, 1);
                
                % Update the effective PLDs
                PLDTry_curr = zeros(size(param.tau));
                for pp = [1:(ii-1),(ii+1):param.multiPLD]
                    indCurr = (pp-1)*numPLD+1 : (pp-1)*numPLD + numPLD;
                    PLDTry_curr(indCurr,:,:) = cumsum(param.tau(indCurr(end:-1:1),:,:)) - param.tau(indCurr(end:-1:1),:,:); % Count up backwards
                    PLDTry_curr(indCurr,:,:) = PLDTry_curr(indCurr(end:-1:1),:,:); % Order so the blocks are in order
                    PLDTry_curr(indCurr,:,:) = PLDTry_curr(indCurr,:,:) + PLD(pp);
                end
                PLDTry_curr(ind,:,:) = cumsum(param.tau(ind(end:-1:1),:,:)) - param.tau(ind(end:-1:1),:,:); % Count up backwards
                PLDTry_curr(ind,:,:) = PLDTry_curr(ind(end:-1:1),:,:); % Order so the blocks are in order
                PLDTry_curr(ind,:,:) = PLDTry_curr(ind,:,:) + repmat(PLDTry(:)', numPLD, 1, min(TTStep,tauTryIndL));
                
                variance = zeros(PLDTryL, min(TTStep,tauTryIndL), distL, allSliceL);
                for kk = 1:allSliceL
                    slice = allSlice(kk);
                       
                    param.PLD = PLDTry_curr + ((slice-1)*slicedt);
                    param.t = param.tau + param.PLD;
                    [variance(:,:,:,kk)] = Hessian_LOptimal_analytical(param,A,scanTime,slice,slicedt);
                    
                end % End slice loop
                
                variance = variance / ((numPLD + 1)/2);
                
                variance = variance * 6000 * 6000; % Change into (ml/100g/min)^2
                
                % Take mean of generalised variance across the slices
                varianceMean = mean(variance,4);
                
                if any(varianceMean(:)==0)
                    warning('Some variances are zero. Setting to inf...')
                    varianceMean(varianceMean==0) = inf;
                end
                
                % Take mean of generalised variance across the BAT distribution
                cost = distWeightCurr .* varianceMean;
                cost(distWeightCurr==0) = 0;            % To correct for 0*nan in distWeight .* CovCurr
                cost(isnan(cost)) = inf;                  % To correct for 0*inf in distWeight .* CovCurr
                costMean(:, tauTryInd) = mean( cost , 3 );    % Weighted mean
                
            end % end loop through TT
            
            % Find the Tau and PLD that leads to the minimum generalised variance
            [minVariance,jj] = min(costMean(:));
            [row,col] = ind2sub(size(costMean),jj);
            
            % Save the optimal PLD and LD
            PLD(ii) = PLDTry(row);
            tau(ind(ll)) = tauTry(col);
            
            disp(['Tau = ' array2string(tau,'; ')])
            disp(['PLD = ' array2string(PLD,'; ')])
            
        end % End PLD loop
        
    end
    
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
    
    PLD_effectiveShift = zeros(numPLD*param.multiPLD,1);
    PLDTry_curr = zeros(numPLD*param.multiPLD,1);
    for pp = [1:(ii-1),(ii+1):param.multiPLD]
        indCurr = (pp-1)*numPLD+1 : (pp-1)*numPLD + numPLD;
        PLD_effectiveShift(indCurr) = cumsum(bestTau(indCurr(end:-1:1))) - bestTau(indCurr(end:-1:1)); % Count up backwards
        PLD_effectiveShift(indCurr) = PLD_effectiveShift(indCurr(end:-1:1)); % Order so the blocks are in order
        PLDTry_curr(indCurr) = PLD_effectiveShift(indCurr) + bestPLD(pp);
    end
    
    param.PLD = PLDTry_curr;
    param.tau = bestTau;
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
