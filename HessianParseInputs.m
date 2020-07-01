function param = HessianParseInputs(param)

%%%% Function to set up inputs %%%%
param.t   = squeeze(param.t);
param.tau = squeeze(param.tau);
if length(param.tau)==1
    param.tau = param.tau*ones(size(param.t));
end
param.f   = param.f(:);
param.BAT = param.BAT(:);

% Resize: tau needs to be the same size as t
lF    = length(param.f);
lBAT  = length(param.BAT);
tSize = size(param.t);
nDimOnes = ones(1,sum(tSize>1)); % Doesn't include last dimension if it is 1 (i.e. nx1 array)
if     isempty(nDimOnes);         nDimOnes = 1;
elseif tSize(1)==1 && tSize(2)>1; nDimOnes = [1,1]; end
param.t   = repmat(param.t  ,[nDimOnes,lBAT,lF]);
param.tau = repmat(param.tau,[nDimOnes,lBAT,lF]); 

% Make f and BAT the same dimensions as t
% shiftdim(A,2) rotates the dims to the left by 2
param.f   = shiftdim(repmat(param.f' ,[lBAT,1,tSize(logical(nDimOnes))]),2);
param.BAT = shiftdim(repmat(param.BAT,[1,lF,tSize(logical(nDimOnes))])  ,2);

if tSize(1)==1
    warning('Only one TI. Adjusting arrays accordingly')
    param.f   = param.f(:)';
    param.BAT = param.BAT(:)';
end

end