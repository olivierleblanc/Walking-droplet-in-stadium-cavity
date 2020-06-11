%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Olivier Leblanc
% Date : 16/05/2020
%
% Function :
% This script aims to analyze images or a video of the vibrating bath to
% estimate the free surface profile of the liquid, in this case silicon oil.
%
% Inputs :
% - image or video directory
% - ref_num :  Number of the reference image
% - cur_num :  Number of the current image
% - ref_length : length of a reference object to compute the resolution
%
% Outputs :
% /.
%
% Options :
% analyze_video  % 0 for analyzing a directory of images
define_resolution = 0;
define_ROI = 0;
lighting_effects = 0;  % Enhance surface visualization with ligthning effect
true_height = 0;    % 1 : gives true reconstructed height. 0 : gives height difference from equilibrium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze_video = 0; % 0 for analyzing a directory of images
% folder = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Visualization\Images\';
% directory = '10-04-20_circle';   % Directory name
% filename = 'circle';
% format = '.jpeg';
% ref_num = 10;   % Number of the reference image
% cur_num = 196;
% ref_length = 0.02;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze_video = 0; % 0 for analyzing a directory of images
% folder = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Visualization\Images\';
% dir = '09-04-20_App';   % Directory name
% filename = 'app';
% ref_num = 0;   % Number of the reference image
% cur_num = 502;
% ref_length = 0.03;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze_video = 0; % 0 for analyzing a directory of images
% folder = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Visualization\Images\';
% dir = '10-04-20_drop';   % Directory name
% filename = 'drop';
% ref_num = 10;   % Number of the reference image
% cur_num = 1261;
% ref_length = 0.02;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze_video = 0; % 0 for analyzing a directory of images
% folder = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Visualization\Images\';
% dir = '10-04-20_faraday';   % Directory name
% filename = 'faraday';
% ref_num = 10;   % Number of the reference image
% % cur_num = 290;
% cur_num = 800;
% ref_length = 0.02;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
window_size = 16;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analyze_video = 0; % 0 for analyzing a directory of images
folder = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Visualization\Images\';
directory = '10-04-20_circle';   % Directory name
filename = 'circle';
format = '.jpeg';
ref_num = 10;   % Number of the reference image
cur_num = 196;
ref_length = 0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Physical system values
lambda = 4.72e-3;   % Typical wavelength
h0 = 6e-3; % Water level
n = 1.4; % Silicon oil refractive index
alpha = 1-1/n;  % See paper
hp = 6.5e-3;  % Liquid depth
H = 0.62;  % Distance pattern-camera
h_star = 1/(1/(alpha*hp)-1/H);

if (analyze_video)
    obj = VideoReader([folder, directory, filename, format]);
    % % Automatically find the number of image files in the folder.
    % /!\ When the video is a MTS file, it gives 50 FPS instead of the true 25 FPS.
    number_of_images = floor(obj.Duration*obj.FrameRate/2);
else
    number_of_images = numel(dir([folder, directory,'\*', format]))-1;    % Numbers of images in the folder or to treat
end

% Load reference image
if (analyze_video)
    Iref = read(obj, ref_num);
else
    Iref = imread([folder, directory,'\',filename,'_',sprintf('%06d',ref_num), format]);
end
Iref_gray = rgb2gray(Iref);
% Finds the orientation and resolution of the images.
[Im_res_y, Im_res_x, ~] = size(Iref_gray);
% Load current image
if (analyze_video)
    I2 = read(obj, cur_num(1));
else
    I2 = imread([folder, directory,'\',filename,'_',sprintf('%06d',cur_num(1)), format]);
end
I2_gray = rgb2gray(I2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (define_resolution)
    figure(1); hold on;
    imshow(I2_gray);
    axis on;
    title('\color{red} please click on the left side of the gray base', 'FontSize', 20);
    [xl, yl] = ginput(1);
    title('\color{red} please click on the right side of the gray base', 'FontSize', 20);
    [xr, yr] = ginput(1);
    pixel2mm = ref_length/(sqrt((xr-xl)^2+(yr-yl)^2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (define_ROI)
    figure(1); hold on;
    imshow(I2_gray);
    axis on;
    title('\color{red} Now select your ROI', 'FontSize', 20);
    rect = round(getrect());
    x_ul = rect(1);
    y_ul = rect(2);
    wx = rect(3);
    hy = rect(4);
    N = wx*hy;
    
    ROI_indx = x_ul:x_ul+wx;
    ROI_indy = y_ul:y_ul+hy;
    
    close(figure(1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Iref_ROI = Iref(ROI_indy, ROI_indx);
I2_ROI = I2(ROI_indy, ROI_indx);

dx = zeros(wx, hy);
dy = zeros(wx, hy);

%% Plot

% Load current image
if (analyze_video)
    I2 = read(obj, cur_num);
else
    I2 = imread([folder, directory,'\',filename,'_',sprintf('%06d',cur_num),'.jpeg']);
end
I2_gray = rgb2gray(I2);

% First correct for an eventual shift between reference and current image
Iref_gray2 = double(Iref_gray);
m_ref = mean(mean(Iref_gray2));
Iref_window = Iref_gray2-m_ref;
fourier_Iref_corr = fftshift(fft2(Iref_window));
I2_gray2 = double(I2_gray);
m_2 = mean(mean(I2_gray2));
I2_window = I2_gray2-m_2;
fourier_I2_corr = fftshift(fft2(I2_window));
fourier_prod = fourier_Iref_corr.*conj(fourier_I2_corr);
corr = fftshift(real(ifft2(ifftshift(fourier_prod))));
[val y] = max(corr);
[val2 x] = max(val);
y = y(x);
% [X,Y] = meshgrid(1:Im_res_x, 1:Im_res_y);
% figure(10); hold on;
% si = surf(X,Y, corr); colorbar;
% si.EdgeColor = 'none';
ROI_indy2 = ROI_indy-(y-Im_res_y/2)+1;
ROI_indx2 = ROI_indx-(x-Im_res_x/2)+1;
I2_ROI = I2(ROI_indy2, ROI_indx2);
x-Im_res_x/2
y-Im_res_y/2
Iref_gray = double(Iref_gray);
I2_gray = double(I2_gray);

f1 = figure(1); hold on;
set(f1, 'position',[0 200 500 500]);
imshow(Iref_ROI);
axis on;
% title('Reference image');
% xlabel('X [m]');
% ylabel('Y [m]');
ax = gca;
xTick = get(ax,'XTick');
set(ax,'XTick',downsample(xTick,2));
yTick = get(ax,'YTick');
set(ax,'YTick',downsample(yTick,2));
ax.XTickLabel = round(ax.XTick.*pixel2mm.*1e4)./1e4;
ax.YTickLabel = round(ax.YTick.*pixel2mm.*1e4)./1e4;
drawnow;

f2 = figure(2); hold on;
set(f2, 'position',[400 200 500 500]);
imshow(I2_ROI);
axis on;
% title('Modified image');
ax = gca;
xTick = get(ax,'XTick');
set(ax,'XTick',downsample(xTick,2));
yTick = get(ax,'YTick');
set(ax,'YTick',downsample(yTick,2));
ax.XTickLabel = round(ax.XTick.*pixel2mm.*1e4)./1e4;
ax.YTickLabel = round(ax.YTick.*pixel2mm.*1e4)./1e4;
drawnow;

f20 = figure(20); hold on;
set(f20, 'position',[400 200 500 500]);
axis on;
title('Modified image');
imshow(I2_ROI-Iref_ROI);
axis on;
drawnow;

% f11 = figure(11); hold on; colorbar;
% [X,Y] = meshgrid(1:window_size,1:window_size);
% f12 = figure(12); hold on; colorbar;
% s = handle( surf( NaN(3) ) ) ;     %// create an empty "surface" object
% set( s , 'XData',X , 'YData',Y , 'ZData',X );
% sc = scatter3( 1,1,1, 'r' ) ;
% set( sc , 'XData', 1 , 'YData',1, 'ZData', 1 );

% For each pixel of the reference image ROI
for i = 1:wx  %wx
    for j = 1:hy  %hy
        % for i = 90:100  %wx
        %     for j = 100:110  %hy
        Iref_window = Iref_gray(ROI_indy(j)-window_size/2:ROI_indy(j)+window_size/2-1 ,...
            ROI_indx(i)-window_size/2:ROI_indx(i)+window_size/2-1);
        m_ref = mean(mean(Iref_window));
        std_ref = sqrt( sum(sum( (Iref_window-m_ref).^2 ))./(N) );
        Iref_window = Iref_window-m_ref;
        fourier_Iref_corr = fftshift(fft2(Iref_window));
        
        I2_window = I2_gray(ROI_indy2(j)-window_size/2:ROI_indy2(j)+window_size/2-1 ,...
            ROI_indx2(i)-window_size/2:ROI_indx2(i)+window_size/2-1);
        m_2 = mean(mean(I2_window));
        std_2 = sqrt( sum(sum( (I2_window-m_2).^2 ))./(N) );
        I2_window = I2_window-m_2;
        fourier_I2_corr = fftshift(fft2(I2_window));
        % Pas oublier que la corrélation c'est pas exactement une
        % convolution donc dans Fourier il faut le complexe conj d'un des deux
        fourier_prod = fourier_Iref_corr.*conj(fourier_I2_corr);
        corr = fftshift(real(ifft2(ifftshift(fourier_prod))));
        [val y] = max(corr);
        [val2 x] = max(val);
        y = y(x);
        
        %         %         surf(X,Y, corr); colorbar;
        %         s.ZData = corr ;
        %         title(['pixel (', num2str(i),',', num2str(j),')']);
        %         legend(['max at (', num2str(x), ',', num2str(y), ')']);
        %         %         scatter3(x,y,val2+10,'r');
        %         set( sc , 'XData', x , 'YData',y, 'ZData', val2+10 );
        %         drawnow;
        
        dx(i,j) = x-window_size/2-1;
        dy(i,j) = y-window_size/2-1;
        depl(i,j) = sqrt(dx(i,j)^2+dy(i,j)^2);
    end
end

% Rescale displacement in true distance [m]
dx = dx.*pixel2mm;
dy = dy.*pixel2mm;
% Inverse gradient
dhdx = -dx./h_star;
dhdy = -dy./h_star;
% Remove mean to avoid integrating a constant
dhdx = dhdx-mean(mean(dhdx));
dhdy = dhdy-mean(mean(dhdy));

if (true_height)
    h = intgrad2(dhdy,dhdx,pixel2mm,pixel2mm,h0);
else
    h = intgrad2(dhdy,dhdx,pixel2mm,pixel2mm,0);
end

% Naive fourier filtering
ff = fftshift(fft2(h));
ff(abs(ff)<0.005) = 0;
h = real(ifft2(ifftshift(ff)));

%% Plots
% Meshgrid creation
[X,Y] = meshgrid(1:hy,1:wx);

f3 = figure(3); hold on;
set(f3, 'position',[400 100 500 500]);
s2 = surf(X, Y, dx);
s2.EdgeColor = 'none';
colorbar;
title('dx');
ax = gca;
xTick = get(ax,'XTick');
set(ax,'XTick',downsample(xTick,2));
yTick = get(ax,'YTick');
set(ax,'YTick',downsample(yTick,2));
ax.XTickLabel = ax.XTick.*pixel2mm;
ax.YTickLabel = ax.YTick.*pixel2mm;

f4 = figure(4); hold on;
set(f4, 'position',[800 100 500 500]);
s3 = surf(X, Y, dy);
s3.EdgeColor = 'none';
colorbar;
title('dy');

% f5 = figure(5); hold on;
% set(f5, 'position',[800 100 500 500]);
% s4 = surf(X,Y, abs(ff));
% s4.EdgeColor = 'none';
% colorbar;
% title('FFT2');

f7 = figure(7); hold on;
set(f7, 'position',[800 200 500 500]);
s3 = surf(X,Y,h);
s3.EdgeColor = 'none';
% title(['Surface profile of ', num2str(cur_num), '-th image']);
grid on;
if (lighting_effects)
    light               % create a light
    lighting gouraud    % preferred method for lighting curved surfaces
end
L = max(wx, hy);
axis([0 L 0 L]);
direction = [0 0 1];
rotate(s3, direction,-90);
colorbar;
% xlabel('x');
% ylabel('y');
ax = gca;
xTick = get(ax,'XTick');
set(ax,'XTick',downsample(xTick,2));
yTick = get(ax,'YTick');
set(ax,'YTick',downsample(yTick,2));
ax.XTickLabel = round(ax.XTick.*pixel2mm.*1e4)./1e4;
ax.YTickLabel = round(ax.YTick.*pixel2mm.*1e4)./1e4;
