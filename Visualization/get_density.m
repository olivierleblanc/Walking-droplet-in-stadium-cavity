%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Olivier Leblanc
% Date : 06/11/2020
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
plot_trajectory =  1;
colortype = 'index'; % 'index' or 'speed'
plot_pos_histogram = 1;
plot_speed_histogram = 0;
plot_R_histogram = 0;
plot_L_histogram = 0;

define_resolution = 1;
analyze_video = 0; % 0 for analyzing a directory of images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
nbins = 50;
% Constants
ref_length = 0.095;
h = 6.62607015*1e-34; % Planck's constant
h_ = h/(2*pi);

close all;
clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% directory = 'C:\Users\leblanco.OASIS\Documents\IngeCivilMaster2\Memoire\Liege_Data\J3_Corrected_Positions\';
% filename = 'High_mem_A_position_corrected.mat';
% text = 'C:\Users\leblanco.OASIS\Documents\IngeCivilMaster2\Memoire\Liege_Data\J3_30-10\G_4_00\';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directory = 'C:\Users\leblanco.OASIS\Documents\IngeCivilMaster2\Memoire\Liege_Data\J7_20_11\';
filename = '3_G_4_50_D_positions.mat';
text = 'C:\Users\leblanco.OASIS\Documents\IngeCivilMaster2\Memoire\Liege_Data\J7_20_11\3_G_4_50_D\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([directory, filename]);
% drop_position = drop_position_corrected;
imnames = dir([text, '*.jpg']);

if (analyze_video)
    videoname = '00062_f80_1910_M_D105.MTS';
    obj = VideoReader([directory, videoname]);
    I1 = read(obj, 1);
else
    I1 = imread([text, imnames(1).name]);
end
    rgbImage = cat(3, I1,I1,I1);

% Droplet properties
D = 0.887e-3; % Diameter
R0 = D/2;
rho = 0.95e3; % Silicon oil density
m = rho*4*pi*R0^3/3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (define_resolution)
    figure(1); hold on;
    imshow(I1);
    axis on;
    title('\color{red} please click on the left side of the gray base', 'FontSize', 20);
    [xl, yl] = ginput(1);
    title('\color{red} please click on the right side of the gray base', 'FontSize', 20);
    [xr, yr] = ginput(1);
    pixel2mm = ref_length/(sqrt((xr-xl)^2+(yr-yl)^2));
    close(figure(1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute center of trajectory
val = mean(drop_position);
xC = val(1);
yC = val(2);

% dt = 1/25; % Time between two frames.
% v = (drop_position(2:end, :) - drop_position(1:end-1, :)).*pixel2mm./dt; % Droplet speed
% norm_v = sqrt(v(:,1).^2+v(:,2).^2);
% ind = find(norm_v>0.2);
% norm_v(ind, : ) =[];
% drop_position(ind, :)=[];
% v(ind, :)=[];
% p = m.*v; % Momentum
% R = sqrt(sum( (drop_position-[xC, yC])'.^2 )).*pixel2mm; % Droplet distance to the center of its trajectory
% r = (drop_position-[xC, yC]);
% L = r(1:end-1,1).*p(:,2) - r(1:end-1,2).*p(:,1); % Angular momentum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (plot_trajectory)
    f1 = figure(1); hold on;
    set(f1, 'position',[0 200 500 500]);
    imshow(rgbImage); axis on; hold on;
%     axis([200, 850, 150,500]);
%     title('Reference image', 'interpreter', 'latex');
    
    scatter(drop_position(1,1),drop_position(1,2), 'b');
    n = length(drop_position(:,1));
    xx=[drop_position(:,1) drop_position(:,1)];
    yy=[drop_position(:,2) drop_position(:,2)];
    zz=zeros(size(xx));
    if (colortype == 'index')
        hs=surf(xx,yy,zz,[1:n; 1:n]','EdgeColor','interp') %// color binded to index values
    else
        norm_v2 = [ norm_v(1); norm_v];
        hs=surf(xx,yy,zz,[norm_v2 norm_v2],'EdgeColor','interp'); %// color binded to droplet's speed
    end
    plot(xC, yC, '.k', 'MarkerSize', 10.0);  
    colormap('jet')
    drawnow;
    legend('Initial position', 'Trajectory', 'Center of the trajectory', 'interpreter', 'latex', 'Location', 'SouthWest');
    latexfig(gca, pixel2mm);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( plot_pos_histogram)
    f2 = figure(2); hold on;
    set(f2, 'position',[600 200 500 500]);
    imshow(rgbImage); axis on; hold on;
    
    drop_pos = repmat(drop_position,20,1);
    hist3(drop_pos, 'Nbins', [nbins,nbins], 'CDataMode','auto', 'FaceColor','interp', 'EdgeColor', 'None');
    colorbar; hold on; 
%     axis([200, 850, 280,600]);
    hold on;
    latexfig(gca, pixel2mm);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (plot_speed_histogram)
    f3 = figure(3); hold on;
    set(f3, 'position',[600 200 550 500]);
    set(gca,'TickLabelInterpreter','Latex', 'FontSize', 12.0);
    norm_v = sqrt(v(:,1).^2+v(:,2).^2);
    hv = histogram(norm_v, 50);
    hv.FaceColor = 'b';
    % h.EdgeColor = 'None';
    xlabel('v [m/s]', 'FontSize', 16.0, 'interpreter', 'latex');
    ylabel('# of occurrences', 'FontSize', 16.0, 'interpreter', 'latex');
    % legend('Histogram of the droplet speed', 'FontSize', 14.0, 'interpreter', 'latex');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( plot_R_histogram)
    f4 = figure(4); hold on;
    set(f4, 'position',[600 200 500 500]);
    h = histogram(R, 400);
    h.FaceColor = 'b';
    h.EdgeColor = 'None';
    xlabel('R [m]', 'FontSize', 14.0, 'interpreter', 'latex');
    ylabel('\# of occurrences', 'FontSize', 14.0, 'interpreter', 'latex');
    legend('Histogram of R', 'FontSize', 16.0, 'interpreter', 'latex', 'FontSize', 12.0);
    set(gca,'TickLabelInterpreter','Latex');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( plot_L_histogram)
    f5 = figure(5); hold on;
    set(f5, 'position',[600 200 500 500]);
    h = histogram(L, 50);
    h.FaceColor = 'r';
    % h.EdgeColor = 'None';
    xlabel('L [kg.m^2/s]', 'FontSize', 14.0, 'interpreter', 'latex');
    ylabel('# of occurrences', 'FontSize', 14.0, 'interpreter', 'latex');
    legend('Histogram of L', 'FontSize', 16.0, 'interpreter', 'latex');
    set(gca,'TickLabelInterpreter','Latex');
end


%% Work on empirical wavefunction
[X,Y] = meshgrid(1:nbins,1:nbins);

[h, b] = hist3(drop_position, 'Nbins', [nbins,nbins], 'CDataMode','auto', 'FaceColor','interp');
h = rot90(h,1);
pdf = h./sum(sum(h));
phi = sqrt(pdf);

% % Momentum operator
% [px,py] = gradient(phi);
% px = -h_.*px;
% py = -h_.*py;
% p = sqrt(px.^2 + py.^2);
% f11 = figure(11); hold on;
% % set(f6, 'position',[200 100 1000 700]);
% % surf(X,Y,p);
% quiver(X,Y,px,py); colorbar; axis on;

% Position expectation
span = [b{1}(end)-b{1}(1), b{2}(end)-b{2}(1)]./nbins;
E_x = sum(sum(X.*pdf));
E_y = sum(sum(Y.*pdf));
figure(1); hold on;
plot(E_x.*span(1)+b{1}(1), E_y.*span(2)+b{2}(1), '.g', 'MarkerSize', 12.0); 

% % Fourier transform of wavefunction
% FT_phi = fftshift(fft2(phi));
% 
% % Kinetic energy term of the hamiltonian
% cin = -(h_^2/2).*del2(phi, pixel2mm); 
% 
% f6 = figure(6); hold on;
% set(f6, 'position',[200 100 1000 700]);
% subplot(121); hold on;
% surf(pdf); colorbar; axis on;
% set(gca,'TickLabelInterpreter','Latex');
% subplot(122); hold on;
% surf(cin); colorbar; axis on;
% set(gca,'TickLabelInterpreter','Latex');

% f7 = figure(7); hold on;
% set(f7, 'position',[200 100 1000 700]);
% subplot(121); hold on;
% surf(real(FT_phi)); colorbar; axis on;
% set(gca,'TickLabelInterpreter','Latex');
% subplot(122); hold on;
% surf(imag(FT_phi)); colorbar; axis on;
% set(gca,'TickLabelInterpreter','Latex');

% f8 = figure(8); hold on;
% set(f8, 'position',[200 100 1000 700]);
% subplot(121); hold on;
% for i = 1:50
%     plot(abs(FT_phi(:,i)));
% end
% subplot(122); hold on;
% for i = 1:50
%     plot(abs(FT_phi(i,:)));
% end

[V,~] = eig(cin*phi);
E = eig(cin*phi);
figure(9); hold on;
scatter(real(E), imag(E)); grid on;
figure(10); hold on;
plot(E(find(imag(E)==0))); grid on;
