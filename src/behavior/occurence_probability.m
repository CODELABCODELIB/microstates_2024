function [f,occurence_state,occurence_triad] = occurence_probability(triad,microstates,state)

selected_microstates = microstates(triad(1):triad(3));
% find the first index of each microstate sequence
first_index = [1,1+find(diff(selected_microstates))]; 

% how many times does the microstate appear 
occurence_state = sum(selected_microstates(first_index) == state);
% all the potential microstates occuring during the triad
occurence_triad = length(selected_microstates(first_index));

f = occurence_state/occurence_triad;
end