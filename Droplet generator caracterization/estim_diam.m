%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Olivier Leblanc
% Date : 24/04/2020
% 
% Function : 
% This script aims to estimate the droplet diameter inside the
% cavity by letting the user draw a line along the droplet diameter
% multiple times.
%
% Inputs : 
% - the chosen images for the diameter estimation
% - a reference distance to compute the actual resolution of the camera
%
% outputs : 
% - the average of the estimated diameters for each image.
%
% Options :
% /
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Informations varying with the analyzed video sequence :
folder_name = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Experiments\24-04-20_Me5_9';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory name
directory = '';
im_name = 'drop';
image_format = '.jpg';
ref_length = 0.09;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Num_droplets = 1;
text = [folder_name, directory, '\', im_name, '_'];
% Automatically find the number of image files in the folder.
number_of_images = numel(dir([folder_name, directory,'\*', image_format]));    % Numbers of images in the folder or to treat
% ------------------------------------------------------------------------------
% Let the user click on some reference point to estimate the image's pixel size
I1 = imread([text, sprintf('%06d',1), image_format]);
[Im_res_y, Im_res_x, ~] = size(I1); % Finds the orientation and resolution of the images.
figure(1); hold on;
imshow(I1);
axis on;
title('\color{red} please click on the left side of the gray base', 'FontSize', 20);
[xl, yl] = ginput(1);
title('\color{red} please click on the right side of the gray base', 'FontSize', 20);
[xr, yr] = ginput(1);
pixel2mm = ref_length/(xr-xl);

corr_size = round(4e-3/pixel2mm); % One droplet is more or less 1mm diameter
corr_size = corr_size + (mod(corr_size, 2)==1);
half_corr_size = round(corr_size/2);
ROI_size = 100;
half_ROI_size = round(ROI_size/2);
shift_x = round((Im_res_x - corr_size)/2);
shift_y = round((Im_res_y - corr_size)/2);

close(figure(1));
% ------------------------------------------------------------------------------

diams = zeros(1, number_of_images);

for i = 1:number_of_images
    % ------------------------------------------------------------------------------
    % Let the user select the droplet of interest.
    % Define centers of initial regions of interest (ROI)
    figure(10); hold on;
    I1 = imread([folder_name, directory, '\', im_name, '_', sprintf('%06d',i), image_format]);
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
    
    % ------------------------------------------------------------------------------
    % Let the user click on two sides of the droplet to estimate its diameter
    I1 = imread([text, sprintf('%06d', i), image_format]);
    I1_corr = I1([corrzone(1,:)],[corrzone(2,:)]);
    I1_ROI = I1([ROI(1,:)],[ROI(2,:)]);
    f1 = figure(1); hold on;
    set(f1, 'position',[500 500 800 800]);
    imshow(I1_corr);
    axis on;
    title('\color{red} please click on one side of the droplet', 'FontSize', 20);
    [xi,yi] = getline;
    diams(i) = sqrt( (xi(1)-xi(2))^2+(yi(1)-yi(2))^2 )*pixel2mm;
    meandiams = mean(diams);
    close(figure(1));
    % ------------------------------------------------------------------------------
    
end

disp([num2str(number_of_images), ' images have been analyzed.']);
disp(['The average over the droplets diameters estimations is : ', num2str(meandiams), ' [m]']);

