function [DCalign, ind_max]=alignDC(DC)

% INPUT
% DC: matrix whose columns are the 52 week unaligned filtered surges
% we want to align all surges, such their the global maximum occur at time
% instant 26

% OUTPUT
% DCaligned: matrix whose columns are the 52 week aligned filtered surges 
% ind_max: time indices of the global maximum of each surge. We will use
% this data to collect statistics on the time variability of the surge peaks

%maxDC=max(DC);  % obtains the global maximum of each column of DC. Then align 

ind_max = zeros(size(DC, 2), 1);  % Initialize the index array for global maxima
for i = 1:size(DC, 2)
    [~, ind_max(i)] = max(DC(:, i));  % Find the index of the global maximum for each surge
end

% Align the surges based on the global maximum indices
DCalign = zeros(size(DC));  % Initialize the aligned matrix
for i = 1:size(DC, 2)
    shift = 26 - ind_max(i);  % Calculate the shift needed to align the maximum
    DCalign(:, i) = circshift(DC(:, i), shift);  % Shift the surge
end