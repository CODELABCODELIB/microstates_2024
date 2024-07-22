function [xi,bin_edges] = create_jid_grid(BINS,MIN_H,MAX_H)
% Calcuates the Joint Interval Distribution (JID) for the list of taps in
% input. 
%    jid = taps2JID(taps, ...);
%    Positional parameters:
% 
%      taps           An array of tap times in ms. 
%    
%    Optional input parameters:  
% 
%      'Bins'           Number of bins (per side) in the JID. Basically the 
%                       length of the side of the JID matrix. 
%      'MIN_H'           The minimum delta(t) value to consider in the JID 
%                       space expressed in log10(ms) space. Default 1.5 ~ 30 ms 
%                       10 ^ 1.5 = 31.6.
%      'MAX_H'           The maximum delta(t) value to consider in the JID 
%                       space expressed in log10 space. Default 5 ~  100 s 
%                       10 ^ 5 = 100000.
%
%   Returns: log transformed grid of size BINS, bin_edges of JID 
%
%
% Enea Ceolini, Leiden University, 26/05/2021
% output bins values in linear space

arguments
BINS = 50;
MIN_H = 1.5;
MAX_H = 5.0;
end

gridx = linspace(MIN_H, MAX_H, BINS);
bin_edges = 10 .^ gridx;

[x1, x2] = meshgrid(gridx, gridx);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
end


