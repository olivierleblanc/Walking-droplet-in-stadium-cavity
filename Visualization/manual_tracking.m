%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Olivier Leblanc
% Date : 17/05/2020
%
% Function :
% This script aims to let the user manually estimating the walking droplet position at any time.
%
% Inputs :
% - image or video directory
%
% Outputs :
% /.
%
% Options :
analyze_video = 1; % 0 for analyzing a directory of images
estim_diameter = 0; % If the use wishes to estimate the droplet diameter
show_tracked_trajectory = 1;
show_true_trajectory = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc;

%% Parameters
Image_dist = 2; % Initial spacing (in images numbers) between the actual reference image and the analyzed image
Num_droplets = 1; % Number of droplets to track

% Informations varying with the analyzed video sequence :
folder_name = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Visualization\Video_Recordings\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory name
directory = '';
filename = '18-05-20_eigen';
image_format = '.MTS';
ref_length = 0.09;
first_image = 3;
last_image = 600;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = VideoReader([folder_name, directory, filename, image_format]);
% Automatically find the number of image files in the folder.
if (analyze_video)
    % /!\ When the video is a MTS file, it gives 50 FPS instead of the true 25 FPS.
    number_of_images = floor(obj.Duration*obj.FrameRate/2);
else
    number_of_images = numel(dir([folder_name, directory,'\*', image_format]))-1;    % Numbers of images in the folder or to treat
end

if (last_image > number_of_images)
    last_image = number_of_images;
end

% ------------------------------------------------------------------------------
% Let the user select the droplet of interest.
% Define centers of initial regions of interest (ROI)
figure(10); hold on;
if (analyze_video)
    I1 = read(obj, first_image-2);
else
    I1 = imread([folder_name, directory, '\', filename, '_', sprintf('%06d',first_image), image_format]);
end
[Im_res_y, Im_res_x, ~] = size(I1); % Finds the orientation and resolution of the images.
imshow(I1);
colorbar;
axis on;
title('\color{red} please click on the droplet to define its initial coordinates', 'FontSize', 20);
[xi, yi] = ginput(1);
ROI_centers = [round(yi); round(xi)];
% Strucure containing y ROI values on the first row, and x ROI values on
% the second row, repeated Num_droplets times.
ROI = repmat([[ROI_centers(1,1)-half_ROI_size+1 : ROI_centers(1,1)+half_ROI_size] ;...
    [ROI_centers(2,1)-half_ROI_size+1 : ROI_centers(2,1)+half_ROI_size]],1,Num_droplets);
close(figure(10));
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
if (estim_diameter)
    % Let the user click on two sides of the droplet to estimate its diameter
    if (analyze_video)
        I1 = read(obj, first_image-1);
    else
        I1 = imread([text, sprintf('%06d',first_image), image_format]);
    end
    I1_corr = I1([corrzone(1,:)],[corrzone(2,:)]);
    I1_ROI = I1([ROI(1,:)],[ROI(2,:)]);
    figure(1); hold on;
    imshow(I1_corr);
    axis on;
    title('\color{red} please click on one side of the droplet', 'FontSize', 20);
    [xi,yi] = getline
    diameter = sqrt( (xi(1)-xi(2))^2+(yi(1)-yi(2))^2 )*pixel2mm;
    close(figure(1));
end
% ------------------------------------------------------------------------------

%% Code

% Initialization
drop_position = zeros(Num_droplets*(last_image),2);
image_ind = zeros(last_image-first_image,1);

figure(1); hold on;
colorbar;
axis on;

i=first_image;
while (i+Image_dist < last_image)
    i
    if (i+Image_dist>number_of_images)
        break;
    end
    % Refresh the reference image for next image to analyze
    if (analyze_video)
        Iref = rgb2gray(read(obj, i));
    else
        Iref = rgb2gray(imread([text, sprintf('%06d',i), image_format]));
    end
    Iref_ROI = Iref([ROI(1,:)],[ROI(2,:)]);
    
    imshow(Iref_ROI);
    title(['\color{red} Click on the droplet. Image :', num2str(i)], 'FontSize', 20);
    [x, y] = ginput(1);
    x = round(x);
    y = round(y);
    
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
end

drop_position(find(drop_position(:,1)==0), :) = [];
image_ind(find(image_ind==0)) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------
if (show_true_trajectory)
    tic;
    % Shows the mean image between the first 'nb_for_mean' images
    space_for_mean = 20;
    nb_for_mean = ((last_image-first_image)/(space_for_mean+1));
    if (analyze_video)
        Iref = double(read(obj, first_image));
        Ipre = double(read(obj, first_image));
    else
        Iref = double(imread([text, sprintf('%06d', first_image), image_format]));
        Ipre = double(imread([text, sprintf('%06d', first_image), image_format]));
    end
    for n = 1:nb_for_mean
        %         Inew = double(imread([text, sprintf('%06d', first_image+space_for_mean*n), image_format]));
        Inew = double(read(obj, first_image+space_for_mean*n));
        I_diff = Inew-Ipre;
        Ipre = Inew;
        % Mettre un facteur qui dépende de l'intensité moyenne!!
        %         weights = I_diff./(100/nb_for_mean);
        weights = I_diff./(1000/nb_for_mean);
        Iref = Iref + Inew.*(1 + weights);
    end
    Iref = uint8(round(Iref./nb_for_mean));
    toc;
end
% ----------------------------------------------------------------------

f1 = figure(1); hold on;
set(f1, 'position',[0 200 500 500]);
imshow(Iref); axis on;
title('Reference image');

f10 = figure(10); hold on;
set(f10, 'position',[0 200 500 500]);
imshow(I2); axis on; hold on;
title('Current image');
scatter(x_true, y_true, 'r');

f4 = figure(4); hold on;
set(f4, 'position',[600 100 500 500]);
imshow(I4); axis on;
title('Difference');

f3 = figure(3); hold on;
set(f3, 'position',[300 100 500 500]);
imshow(I2_ROI); axis on;
title('Current image'); hold on;
scatter(x, y, 'r');

% ----------------------------------------------------------------------
if (show_tracked_trajectory)
    figure(1); hold on;
    %     set(gca,'Ydir','reverse');
    %     scatter(drop_position(1,1),drop_position(1,2), 'b');
    scatter(round(xi), round(yi), 'b');
    n = length(drop_position(:,1));
    % put "+ 50" only to shift tracked trajectory and compare with mean of images
    p = plot(drop_position(:,1)  ,drop_position(:,2) , 'r', 'LineWidth', 1.0); hold on;
    % modified jet-colormap
    cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
    drawnow;
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd);
    legend('Initial position', 'Trajectory');
end
% ----------------------------------------------------------------------
