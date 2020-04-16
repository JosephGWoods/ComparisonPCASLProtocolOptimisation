function [param,A_BAT,A_f,nDimOnes] = HessianParseInputs(param,tDim)

%%%% Function to set up inputs %%%%
lT = size(param.t,tDim);
param.t = squeeze(param.t);
param.tau = squeeze(param.tau);
if length(param.tau)==1
    param.tau = param.tau*ones(size(param.t));
end
param.f = param.f(:);
A_f = param.f;
param.BAT = param.BAT(:);
A_BAT = param.BAT;

% Resize: tau needs to be the same size as t
lF = length(A_f);
lBAT = length(A_BAT);
tSize = size(param.t);
nDimOnes = ones(1,sum(tSize>1)); % Doesn't include last dimension if it is 1 (i.e. nx1 array)
if isempty(nDimOnes); nDimOnes = 1;
elseif tSize(1)==1 && tSize(2)>1; nDimOnes = [1,1];
end
param.t = repmat(param.t,[nDimOnes,lBAT,lF]);
param.tau = repmat(param.tau,[nDimOnes,lBAT,lF]); 
% Make f and BAT the same dimensions as t
% shiftdim(A,2) rotates the dims to the left by 2
param.f = shiftdim(repmat(param.f',[lBAT,1,tSize(logical(nDimOnes))]),2);
param.BAT = shiftdim(repmat(param.BAT,[1,lF,tSize(logical(nDimOnes))]),2);

if tSize(1)==1
    warning('Only one TI. Adjusting arrays accordingly')
    param.f = param.f(:)';
    param.BAT = param.BAT(:)';
end

% % Calculate the number of averages
% tReadout = 1;
% if strcmp(param.filename,'var_te_pCASL')
%     TR = param.t(end,:,:) + tReadout;   % TimeEnc specific TR. Extra dimensions relate to tau and PLD parameter space
%     N_ave = 300./((nBlocks+1)*TR);
% else
%     TR = param.t+tReadout;
%     N_ave = 300./(2*(sum(TR,tDim)));
% end
% N_ave = squeeze(N_ave);% No longer multiply N_ave by each part of the
%                             % Hessian and take the determinant. The
%                             % rounding errors could cause small
%                             % determinants (~10^-12) to become negative.
%                             % Now multiply after taking the determinant.
end