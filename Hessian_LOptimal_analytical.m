function objective = Hessian_LOptimal_analytical(param,A,scanTime,slice,slicedt)
% This function calculates an approximation of the Hessian, finds the
% inverse (which is proportional to the covariance matrix) then multiplies
% it by the user specified L-optimal matrix.
%
% param is the struct containing the variables,
%
% Here A is a weighting for the trace
% Written by Joseph G. Woods, FMRIB, Oxford, July 2017

if nargin < 5
    slicedt = 0.053125;
    warning(['HessianWithTRWLOptimalFloorTrace.m: No slicedt specified. Using default slicedt = ' ns(slicedt)])
    if nargin < 4
        slice = 1;
        warning(['HessianWithTRWLOptimalFloorTrace.m: No slice specified. Using default slice = ' ns(slice)])
        if nargin < 3
            scanTime = 300;
            warning(['HessianWithTRWLOptimalFloorTrace.m: No scanTime specified. Using default scanTime = ' ns(scanTime)])
            if nargin < 2
                A = [1,0;0,0];
            end
        end
    end
end

tDim = 1;
param = HessianParseInputs(param);

% Multiply TR by weights to get the correct weighted TotalTR
TRWeight = TRWeightingOrNAveFloor(param,scanTime,tDim,slice,slicedt);
TRWeight = TRWeight./(param.noise^2);

[df, dBAT] = sensitivity(param);
%%%% Form det(Hessian) %%%%
 
% Form det(Hessian)
% Hessian dimensions can be: 4 x PLD x Tau x BAT x f
H(1,1,:,:,:,:) = TRWeight .* squeeze(sum(df   .* df  , tDim));
H(1,2,:,:,:,:) = TRWeight .* squeeze(sum(df   .* dBAT, tDim));
H(2,1,:,:,:,:) = TRWeight .* squeeze(sum(dBAT .* df  , tDim));
H(2,2,:,:,:,:) = TRWeight .* squeeze(sum(dBAT .* dBAT, tDim));

% Use arrayfun to apply inv to data efficiently
Hsize   = size(H);
tmp     = reshape(H,2,2,prod(Hsize(3:end)));
inverse = arrayfun(@(i) A.*jw_inv(tmp(:,:,i)),1:size(tmp,3),'UniformOutput',false);
tmp     = cell2mat(inverse);
cov     = reshape(tmp,2,2,prod(Hsize(3:end)));

cov = abs(cov);
cov(isnan(cov)) = inf; % To correct for inf*0 errors in A.*inverse

r = arrayfun(@(i) (trace(cov(:,:,i))) ,1:size(cov,3),'UniformOutput',true);

if length(Hsize(3:end)) == 1
    objective = reshape(r,1,Hsize(3));
else
    objective = reshape(r,Hsize(3:end));
end

end
