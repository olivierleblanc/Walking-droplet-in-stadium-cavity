function [diam] = diam_from_speed (GdivGF, speed)

close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Olivier Leblanc
% Date : 01/11/2020
%
% Function :
% Computes the droplet diameter from its mean speed and the memory.
%
% Inputs :
% - GdivGF : ratio between the bath vertical acceleration and the Faraday threshold.
% - speed : the droplet speed, in mm/s
%
% Outputs :
% - diam : the diameter, mm
%
%
% Examples :
%
% [diam] = diam_from_speed (0.8, 12)
%
% Options :
%   /
%
% Parameters :
folder = 'C:\Users\leblanco.OASIS\Documents\IngeCivilMaster2\Memoire\IEEE_Paper\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% G_div_GF
load([folder, 'G_div_GF_Bush_Molacek.mat'])
% speed_vs_acc curves 
load([folder, 'speed_vs_acc_data_Bush_Molacek.mat'])
% diameters related to curves
load([folder, 'related_diameters_Bush_Molacek.mat'])

%% Find 2 nearest G/GF
dist = abs(GdivGF-G_div_GF);
[val1, ind1] = min(dist);
dist(ind1) = max(dist);
[val2, ind2] = min(dist);
w1 = val2./(val1+val2);
w2 = 1-w1;

%% Interpolated speeds and find 2 nearest
speeds = w1.*speed_vs_acc(:,ind1)+w2.*speed_vs_acc(:,ind2);

dist = abs(speed-speeds);
[val3, ind3] = min(dist);
dist(ind3) = max(dist);
[val4, ind4] = min(dist);
w3 = val4./(val3+val4);
w4 = 1-w3;

diam = w3.*diameters(ind3)+w4.*diameters(ind4);

end

