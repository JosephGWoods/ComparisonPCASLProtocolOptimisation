function [TRWeight,TotalTR] = TRWeightingOrNAveFloor(param,scanTime,tDim,slice,slicedt)

% TRWeight is the number of averages
if nargin < 5
    slicedt = 0.053125;
    warning(['TRWeightingOrNAveFloor.m: No slicedt specified. Using default slicedt = ' ns(slicedt)])
    if nargin < 4
        slice = 1;
        warning(['TRWeightingOrNAveFloor.m: No slice specified. Using default slice = ' ns(slice)])
        if nargin < 3
        	tDim = 1;
        	warning(['TRWeightingOrNAveFloor.m: No tDim specified. Using default tDim = ' ns(tDim)])
        end
    end
end

tReadout = 0.58293 + 0.05472; % TGSE 20 slices, 47 EPI factor, with bandwidth = 2298 Hz + presat

% I round the TR since there are occaisionally computational rounding
% errors. I have used 5 decimal places to all dense BAT sampling, but I am
% very unlikely to go finer than 0.001 density.

switch param.filename

    case 'var_te_pCASL'
        
        [~,TEInd] = max(param.t(1:param.num_enc,1,1,1,1),[],1);
        TR = squeeze(param.t(TEInd,:,:,:,:)) + tReadout - ((slice-1)*slicedt);
        TotalTR = round(TR.*(param.num_enc+1),5); % Multiply by the number of images that have to be acquired

    case 'var_te_pCASL_nPLD'
        
        [~,TEInd] = max(param.t(1:param.num_enc,1,1,1,1),[],1);
        TR = 0;
        for ii = 1:param.multiPLD
            ind = TEInd + (ii-1)*param.num_enc;
            TR = TR + squeeze(param.t(ind,:,:,:,:)) + tReadout - ((slice-1)*slicedt);
        end
        TotalTR = round(TR.*(param.num_enc+1),5); % Multiply by the number of images that have to be acquired
        
    case 'var_multi_pCASL'
        
        TR = param.t + tReadout - ((slice-1)*slicedt);
        TotalTR = round((2*(sum(TR,tDim))),5);
        
end

TRWeight = floor(round(squeeze(scanTime./TotalTR),5));

