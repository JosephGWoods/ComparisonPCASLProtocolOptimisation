function [LD_out,PLD_out] = generateTETimings(LD , PLD, nBlocks, order)
%
% GenerateTETimings outputs the corresponding label durations and PLDs in
% natural order: last block with the shortest PLD is first, etc.
%
% Written by Joseph G. Woods, FMRIB, Oxford, 2018
%
% Inputs:
%         LD:       (n-dimensional array)
%            The label durations of the TE blocks in chronological order
%            (the order they would be played out within the TR).
%            1st dimension should the the LDs from a single experiment.
%            Other dimensions can be used for other experiments to
%            calculate timings in one go to avoid looping.
%
%         PLD:      (n-dimensional array)
%            The final PLD(s) after the last block. If length(PLD) > 1,
%            then we assume that multiple PLDs are to be used. These are
%            ordered sequentially. 1st dimension is reserved for hybrid
%            style sequential PLDs.
%
%         nBlocks:  (scalar)
%            The number of time-encoded LDs. For a Hadamard 8x7 matrix,
%            this would be 7.
%
%         order:    (string)
%            The order that the output effective LDs and PLDs are returned.
%            'standard': Standard ordering after decoding (longest to shortest PLD)
%            'natural' : Reverse ordered (shortest to longest PLD with TR)

if exist('order','var') || isempty(order)
    order = 'standard';
end

% Save the original array dimensions
siz_LD  = size(LD);
siz_PLD = size(PLD);

% Reshape for ease of calculation
LD  = reshape(LD, siz_LD(1), []);
PLD = reshape(PLD, siz_PLD(1), []);

% Make sure the array sizes match
if numel(PLD)>1 && any(siz_LD(2:end) ~= siz_PLD(2:end))
    error('2nd and higher dimensions of LDs and PLDs do not match');
end

% Repmat LDs along 1st dimension to match the sequential PLDs
multiPLD = siz_PLD(1);
if siz_LD(1) < multiPLD*nBlocks
    LD = repmat(LD, multiPLD, 1);
    siz_LD(1) = multiPLD*nBlocks;
end

% Loop through any sequential PLDs and calculate timings
LD_out = zeros(size(LD));
PLD_out = zeros(size(LD));
for nn = 1:multiPLD
    
    % Find the LDs for the current sequential PLD
    ind = ((nn-1)*nBlocks)+1:nn*nBlocks;
    
    % Count up backwards to calculate the effective PLDs shifts due to
    % following time-encoded LDs.
    PLD_effectiveShift = cumsum(LD(ind(end:-1:1), :), 1) - LD(ind(end:-1:1), :);
    
    % Add the final PLD
    PLD_out(ind, :) = PLD_effectiveShift + repmat(PLD(nn,:), nBlocks, 1);
    
    % Reorder appropriately
    if strcmp(order, 'standard')
        LD_out(ind, :)  = LD(ind, :);
        PLD_out(ind, :) = PLD_out(ind(end:-1:1), :);
    elseif strcmp(order, 'natural')
        LD_out(ind, :)  = LD(ind(end:-1:1), :);
    end
    
end

% Reshape arrays to the original dimensions
LD_out  = reshape(LD_out, siz_LD);
PLD_out = reshape(PLD_out, siz_LD);

