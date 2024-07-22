function [T_matrix] = calculate_TP(numStates,microstates,options)
% Authors:
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Z�rich, Psychologisches Institut, Methoden der
% Plastizit�tsforschung.
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% September 2017.

% Copyright (C) 2017 Andreas Pedroni, andreas.pedroni@uzh.ch.
% Adapted: Ruchella Kock, Leiden University, 15/11/2023
arguments
    numStates;
    microstates;
    options.parameter = 'tp';
end

states = my_RLE(microstates);
states = states(states ~= 0);
seqLen = length(states);

% get the durations of the sequences before transition
if strcmp(options.parameter, 'duration')
    index = [1,1+find(diff(microstates))];
    durations = diff(index);
elseif  strcmp(options.parameter, 'duration_after')
    index = [find(diff(microstates)),length(microstates)];
    durations = diff(index);
end



if contains(options.parameter, 'duration')
    % prepare output matrix
    tr = cell(numStates);
    % concat the average durations 
    for count = 1:seqLen-1
        tr{states(count),states(count+1)} = [tr{states(count),states(count+1)},durations(count)];
    end
    % calculate the median transition duration
    T_matrix = cellfun(@(x) median(x, 'omitnan'),tr);
elseif strcmp(options.parameter, 'tp')
    % prepare output matrix
    tr = zeros(numStates);
    % count up the transitions from the state path
    for count = 1:seqLen-1
        tr(states(count),states(count+1)) = tr(states(count),states(count+1)) + 1;
    end
    trRowSum = sum(tr,2);
    % if we don't have any values then report zeros instead of NaNs.
    trRowSum(trRowSum == 0) = -inf;
    % normalize to give frequency estimate.
    T_matrix = tr./repmat(trRowSum,1,numStates);
end

end
function [d,c]=my_RLE(x)
    %% RLE
    % This function performs Run Length Encoding to a strem of data x.
    % [d,c]=rl_enc(x) returns the element values in d and their number of
    % apperance in c. All number formats are accepted for the elements of x.
    % This function is built by Abdulrahman Ikram Siddiq in Oct-1st-2011 5:15pm.

    if nargin~=1
        error('A single 1-D stream must be used as an input')
    end

    ind=1;
    d(ind)=x(1);
    c(ind)=1;

    for i=2 :length(x)
        if x(i-1)==x(i)
            c(ind)=c(ind)+1;
        else ind=ind+1;
            d(ind)=x(i);
            c(ind)=1;
        end
    end

end