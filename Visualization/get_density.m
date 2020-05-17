%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Olivier Leblanc
% Date : 24/04/2020
%
% Function :
% This script aims to draw a 2D histogram of the droplet positions
% computed by 'Particle_Tracking.m' on the image of the cavity.
%
% Inputs :
% - 'drop_position.mat'
%
% outputs :
% /.
%
% Options :
plot_trajectory =  0;
plot_pos_histogram = 1;
plot_R_histogram = 1;
plot_L_histogram = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc;

directory = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Experiments\02-05-20\';
filename = 'drop_position.mat';
videoname = 'chaotic.MTS';
load([directory, filename]);

obj = VideoReader([directory, videoname]);
I1 = read(obj, 1);

% Droplet properties
D = 0.887e-3; % Diameter
R0 = D/2;
rho = 0.95e3; % Silicon oil density
m = rho*4*pi*R0^3/3;

% Compute center of trajectory
val = mean(drop_position);
xC = val(1);
yC = val(2);

dt = 1/25; % Time between two frames.
v = (drop_position(2:end, :) - drop_position(1:end-1, :)).*pixel2mm./dt; % Droplet speed
p = m.*v; % Momentum
R = sqrt(sum( (drop_position-[xC, yC])'.^2 )).*pixel2mm; % Droplet distance to the center of its trajectory
r = (drop_position-[xC, yC]);
L = r(1:end-1,1).*p(:,2) - r(1:end-1,2).*p(:,1); % Angular momentum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (plot_trajectory)
    f1 = figure(1); hold on;
    set(f1, 'position',[0 200 500 500]);
    imshow(I1); axis on; hold on;
    title('Reference image');
    
    figure(1); hold on;
    scatter(drop_position(1,1),drop_position(1,2), 'b');
    n = length(drop_position(:,1));
    xx=[drop_position(:,1) drop_position(:,1)];
    yy=[drop_position(:,2) drop_position(:,2)];
    zz=zeros(size(xx));
    hs=surf(xx,yy,zz,[1:n; 1:n]','EdgeColor','interp') %// color binded to index values
    plot(xC, yC, '.g', 'MarkerSize', 10.0);
    colormap('jet')
    drawnow;
    legend('Initial position', 'Trajectory', 'Center of the trajectory');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( plot_pos_histogram)
    f2 = figure(2); hold on;
    set(f2, 'position',[600 200 500 500]);
    imshow(I1); colorbar; hold on;
    hist3(drop_position, 'Nbins', [45,45], 'CDataMode','auto','FaceColor','interp');
    axis on;
    hold on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3 = figure(3); hold on;
set(f3, 'position',[600 200 500 500]);
hv = histogram(v(:,1).^2+v(:,2).^2, 50);
hv.FaceColor = 'b';
% h.EdgeColor = 'None';
xlabel('v [m/s]', 'FontSize', 14.0);
ylabel('# of occurences', 'FontSize', 14.0);
legend('Histogram of the droplet speed', 'FontSize', 16.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( plot_R_histogram)
    f4 = figure(4); hold on;
    set(f4, 'position',[600 200 500 500]);
    h = histogram(R, 400);
    h.FaceColor = 'b';
    h.EdgeColor = 'None';
    xlabel('R [m]', 'FontSize', 14.0);
    ylabel('# of occurences', 'FontSize', 14.0);
    legend('Histogram of R', 'FontSize', 16.0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( plot_L_histogram)
    f5 = figure(5); hold on;
    set(f5, 'position',[600 200 500 500]);
    h = histogram(L, 50);
    h.FaceColor = 'r';
    % h.EdgeColor = 'None';
    xlabel('L [kg.m^2/s]', 'FontSize', 14.0);
    ylabel('# of occurences', 'FontSize', 14.0);
    legend('Histogram of L', 'FontSize', 16.0);
end


