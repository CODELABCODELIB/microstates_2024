function f = background_probability(triad,microstates,state)
f = microstates(triad(3)) == state;
end