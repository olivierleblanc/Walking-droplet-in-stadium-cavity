close all;
clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = '13-04-20_DOD';
% im_name = 'dod';
% image_format = '.jpeg';
% ref_length = 0.02;
% first_image = 1093;
% last_image = 1936;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = '09-04-20_App';
% im_name = 'app';
% image_format = '.jpeg';
% resolution = 0.02;   % 1 pixel = a micrometer
% first_image = 120;
% last_image = 1500;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Options
find_pixel2mm = 0;

show_tracked_trajectory = 1;

%% Parameters
subsample = 1;  % No if 1. Put 10 for test.
Image_dist = 5; % Initial spacing (in images numbers) between the actual reference image and the analyzed image
Num_droplets = 1; % Number of droplets to track

% Informations varying with the analyzed video sequence :
folder_name = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Visualization\Images\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory name
directory = '09-04-20_App';
im_name = 'app';
image_format = '.jpeg';
resolution = 0.02;   % 1 pixel = a micrometer
first_image = 140;
last_image = 1500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

text = [folder_name, directory, '\', im_name, '_'];
% Automatically find the number of image files in the folder.
number_of_images = numel(dir([folder_name, directory,'\*', image_format]))-1;    % Numbers of images in the folder or to treat

% ------------------------------------------------------------------------------
if (find_pixel2mm)
    % Let the user click on some reference point to estimate the image's pixel size
    I1 = imread([text, sprintf('%06d',first_image), image_format]);
    [Im_res_y, Im_res_x, ~] = size(I1); % Finds the orientation and resolution of the images.
    figure(1); hold on;
    imshow(I1);
    axis on;
    title('\color{red} please click on the left side of the gray base', 'FontSize', 20);
    [xl, yl] = ginput(1);
    title('\color{red} please click on the right side of the gray base', 'FontSize', 20);
    [xr, yr] = ginput(1);
    pixel2mm = ref_length/(xr-xl);
    
    corr_size = round(1e-3/pixel2mm); % One droplet is more or less 1mm diameter
    corr_size = corr_size + (mod(corr_size, 2)==1);
    half_corr_size = round(corr_size/2);
    ROI_size = 100;
    half_ROI_size = round(ROI_size/2);
    shift_x = round((Im_res_x - corr_size)/2);
    shift_y = round((Im_res_y - corr_size)/2);
    
    
    close(figure(1));
end
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Let the user select the droplet of interest.
% Define centers of initial regions of interest (ROI)
figure(10); hold on;
I1 = imread([folder_name, directory, '\', im_name, '_', sprintf('%06d',first_image), image_format]);
[Im_res_y, Im_res_x, ~] = size(I1); % Finds the orientation and resolution of the images.
imshow(I1);
colorbar;
axis on;
title('\color{red} please click on the droplet to define its initial coordinates', 'FontSize', 20);
[xi, yi] = ginput(1);
ROI_centers = [round(yi); round(xi)];
close(figure(10));
% ------------------------------------------------------------------------------

% Strucure containing y ROI values on the first row, and x ROI values on
% the second row, repeated Num_droplets times.
ROI = repmat([[ROI_centers(1,1)-half_ROI_size+1 : ROI_centers(1,1)+half_ROI_size] ;...
    [ROI_centers(2,1)-half_ROI_size+1 : ROI_centers(2,1)+half_ROI_size]],1,Num_droplets);
corrzone = repmat([[ROI_centers(1,1)-half_corr_size+1 : ROI_centers(1,1)+half_corr_size-1] ;...
    [ROI_centers(2,1)-half_corr_size+1 : ROI_centers(2,1)+half_corr_size-1]],1,Num_droplets);
%% Code

% Initialization
drop_position = zeros(Num_droplets*last_image-first_image,2);
image_ind = zeros(last_image-first_image);

i=first_image;
%     while (i+Image_dist < last_image)
while (i+Image_dist < first_image+30)
    i
    if (i+Image_dist>number_of_images)
        break;
    end
    % Refresh the reference image for next image to analyze
    Iref = rgb2gray(imread([text, sprintf('%06d',i), image_format]));
    Iref_corr = Iref([corrzone(1,:)],[corrzone(2,:)]);
    Iref_ROI = Iref([ROI(1,:)],[ROI(2,:)]);
    % Load the new image
    I2 = rgb2gray(imread([text, sprintf('%06d',i+Image_dist), image_format]));
    I2_ROI = I2([ROI(1,:)],[ROI(2,:)]);
    
    %%%%%%%%%%%%%%%%% Correlation in Fourier domain %%%%%%%%%%%%%%%%%%%%%%%
    Iref_corr_double = double(Iref_corr);
    m_ref = mean(mean(Iref_corr_double));
    Iref_corr_double = Iref_corr_double-1.2.*m_ref;
    fourier_Iref_corr = fftshift(fft2(Iref_corr_double, ROI_size, ROI_size));
    I2_ROI_double = double(I2_ROI);
    m_2 = mean(mean(I2_ROI_double));
    I2_ROI_double = I2_ROI_double-1.2.*m_2;
    fourier_I2_corr = fftshift(fft2(I2_ROI_double));
    fourier_prod = fourier_Iref_corr.*conj(fourier_I2_corr);
    corr = real(ifft2(ifftshift(fourier_prod)));
    % Normalization
    mean_corr = mean(mean(corr));
    std_corr = sqrt( sum(sum((corr-mean_corr).^2))/(corr_size^2) );
    corr = (corr-mean_corr)./std_corr;
    corr = rot90(corr,3);
    
    %%%%%%%%%%%%%%%%%%%%%%% Image difference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I4 = I2_ROI-Iref_ROI;
    % Filtering
    % I4 = double(I4);
    mean_I4 = mean(mean(I4));
    std_I4 = sqrt( sum(sum( (I4-mean_I4).^2 ))/(ROI_size^2) );
    % I4 = (I4-mean_I4)./std(I4);
    I4 (I4 < mean_I4+15.*std_I4) = 0;
    % I4 = uint8(255.*I4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Combining %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [val y] = max(abs(corr.*double(I4)));
    [val2 x] = max(val);
    y = y(x);
    
    x_true = x + ROI_centers(2,1) - half_ROI_size;
    y_true = y + ROI_centers(1,1) - half_ROI_size;
    
    % Saving current position
    drop_position(i+Image_dist,:) = [x_true,y_true];
    image_ind(i+Image_dist) = i+Image_dist;
    
    % Updates
    i = i+Image_dist;
    ROI_centers = [y_true; x_true];
    ROI = repmat([[ROI_centers(1,1)-half_ROI_size+1 : ROI_centers(1,1)+half_ROI_size] ;...
        [ROI_centers(2,1)-half_ROI_size+1 : ROI_centers(2,1)+half_ROI_size]],1,Num_droplets);
    corrzone = repmat([[ROI_centers(1,1)-half_corr_size+1 : ROI_centers(1,1)+half_corr_size-1] ;...
        [ROI_centers(2,1)-half_corr_size+1 : ROI_centers(2,1)+half_corr_size-1]],1,Num_droplets);
end

drop_position(find(drop_position(:,1)==0), :) = [];
image_ind(find(image_ind==0)) = [];

% Meshgrid creation
[X,Y] = meshgrid(1:ROI_size,1:ROI_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure(1); hold on;
set(f1, 'position',[0 200 500 500]);
imshow(Iref); axis on;
title('Reference image');

f10 = figure(10); hold on;
set(f10, 'position',[0 200 500 500]);
imshow(I2); axis on; hold on;
title('Current image');
scatter(x_true, y_true, 'r');

f2 = figure(2); hold on;
set(f2, 'position',[0 100 500 500]);
title('Reference image');
imshow(Iref_ROI); axis on; hold on;
plot([half_ROI_size-half_corr_size, half_ROI_size+half_corr_size],...
    [half_ROI_size+half_corr_size, half_ROI_size+half_corr_size],'b', 'LineWidth', 2.0);   % Upper_side
plot([half_ROI_size-half_corr_size, half_ROI_size+half_corr_size],...
    [half_ROI_size-half_corr_size, half_ROI_size-half_corr_size],'b', 'LineWidth', 2.0);   % Lower_side
plot([half_ROI_size-half_corr_size, half_ROI_size-half_corr_size],...
    [half_ROI_size-half_corr_size, half_ROI_size+half_corr_size],'b', 'LineWidth', 2.0);   % Left_side
plot([half_ROI_size+half_corr_size, half_ROI_size+half_corr_size],...
    [half_ROI_size-half_corr_size, half_ROI_size+half_corr_size],'b', 'LineWidth', 2.0);   % Right_side

% f12 = figure(12); hold on;
% set(f12, 'position',[0 100 500 500]);
% title('Reference image');
% imshow(Iref_corr); axis on;

f4 = figure(4); hold on;
set(f4, 'position',[600 100 500 500]);
imshow(I4); axis on;
title('Difference');

f3 = figure(3); hold on;
set(f3, 'position',[300 100 500 500]);
imshow(I2_ROI); axis on;
title('Current image'); hold on;
scatter(x, y, 'r');

f5 = figure(5); hold on;
set(f5, 'position',[700 300 500 500]);
axis on; colorbar;
s5 = surf(X,Y, corr);
s5.EdgeColor = 'None';
direction = [0 0 1];
% rotate(s5, direction,-180);
colorbar;
xlabel('x');
ylabel('y');
set(gca,'Ydir','reverse');

f6 = figure(6); hold on;
set(f6, 'position',[1100 50 500 500]);
axis on; colorbar;
s6 = surf(X,Y, corr.*double(I4));
s6.EdgeColor = 'None';
hold on;
scatter3(x, y, val2+10, 'r'); hold on;
set(gca,'Ydir','reverse');

if (show_tracked_trajectory)
    figure(7); hold on;
    imshow(I1); hold on;
    axis on;
    set(gca,'Ydir','reverse');
    scatter(drop_position(1,1),drop_position(1,2), 'b');
    n = length(drop_position(:,1));
    % put "+ 50" only to shift tracked trajectory and compare with mean of images
    p = plot(drop_position(:,1)  ,drop_position(:,2) , 'r', 'LineWidth', 1.5);
    % modified jet-colormap
    cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
    drawnow;
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    legend('Initial position', 'Trajectory');
end


