function ADJ = find_adj_matrix(bins, neigh)
% Creates an adjacency matrix for a JID.
%    [ADJ = find_adj_matrix(bins, neigh);
%
%    Positional parameters:
% 
%      bins             number of bins in the JID.
%      neigh            size of the neighborhood in bins.
%    
%    Output:  
% 
%      ADJ              adjacency matrix of size (bins^2, bins^2)
%
% Enea Ceolini, Leiden University, 26/05/2021

a = 0:bins * bins-1;

ADJ = zeros(bins * bins, bins * bins);

for k = 0:bins * bins - 1
    c = (abs(int32(a / bins) - int32(k / bins)) <= neigh);
    d = (abs(int32(mod(a, bins)) - int32(mod(k, bins))) <= neigh);
    ADJ(:, k + 1) = c & d;
end