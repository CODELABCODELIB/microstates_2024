function [f,duration_state,duration_triad] = duration_probability(triad,microstates,state)  
% full duration of the triad
duration_triad = triad(3)-triad(1);
% duration of the multiple microstates appearing
duration_state = sum(microstates(triad(1):triad(3)) == state);
f = duration_state/duration_triad;
end