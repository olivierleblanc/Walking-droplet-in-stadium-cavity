% Characterization of the droplet on-demand generator

close all;
clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = 'nozzle12_width5_height40';
% im_name = 'fall';
% image_format = '.jpeg';
% image_ind = [237, 314, 669, 1006, 1183];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = 'nozzle12_width5_height40_v2';
% im_name = 'fall';
% image_format = '.jpeg';
% image_ind = [195, 545, 1029, 1132, 1457];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = 'test';
% im_name = 'fall';
% image_format = '.jpeg';
% image_ind = [2, 7, 10, 13, 16, 19];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = 'test';
% im_name = 'fall';
% image_format = '.jpeg';
% image_ind = [4, 8, 12, 16];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Options
find_pixel2mm = 1;
estim_diam = 1;

%% Parameters

% window_size = [20,20];
width = 300; % For the zoom
ref_length = 15.5; % [mm]
% ref_length = 18.0; % [mm]

folder_name = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Droplet_Generator\Caracterization\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory name
directory = 'nozzle12_width5_height40_v2';
im_name = 'fall';
image_format = '.jpeg';
image_ind = [195, 545, 1029, 1132, 1457];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Automatically find the number of files in the folder.
% Numbers of images in the folder or to treat.
number_of_images = numel(dir([folder_name, directory,'\*', image_format]))-1;    % Numbers of images in the folder or to treat

text = [folder_name, directory, '\', im_name, '_'];

%% Code

if (find_pixel2mm)
% ------------------------------------------------------------------------------
% Let the user click on the hexagon of the nozzle to estimate its pixel size
I1 = imread([folder_name, directory, '\', im_name, '_', sprintf('%06d',0), image_format]);
[Im_res_y, Im_res_x, ~] = size(I1); % Finds the orientation and resolution of the images.
figure(1); hold on;
imshow(I1);
axis on;
% title('\color{red} please click on the hexagon of the nozzle', 'FontSize', 20);
% [xi, yi] = ginput(1);
% xi = round(xi)
% yi = round(yi)
% close(figure(1));
% figure(1); hold on;
% imshow(I1(yi-100/2:yi+100/2, xi-width/2:xi+width/2, :));
% axis on;
title('\color{red} please click on the left side hexagon of the nozzle', 'FontSize', 20);
[xl, yl] = ginput(1);
title('\color{red} please click on the right side hexagon of the nozzle', 'FontSize', 20);
[xr, yr] = ginput(1);
pixel2mm = ref_length/(xr-xl);
close(figure(1));
% ------------------------------------------------------------------------------
end

if (estim_diam)
    I1 = imread([folder_name, directory, '\', im_name, '_', sprintf('%06d',0), image_format]);
    droplet_diameters = zeros(1, length(image_ind));
    for i = 1:length(image_ind)
        fig1 = figure(1); hold on;
        I2 = imread([folder_name, directory, '\', im_name, '_', sprintf('%06d',image_ind(i)), image_format]);
        imshow(5.*rgb2gray(I1-I2));
        axis on;
        title('\color{red} Please evaluate the diameter of the droplet in pixels then close the figure', 'FontSize', 20);
        
        while(ishghandle(fig1))
            pause(1);
        end
        disp(['Analyzing ', num2str(image_ind(i)), '-th image']);
        prompt = 'What is the droplet diameter in pixels? ';
        read = input(prompt);
        if (~isempty(read))
        pixel_length = read;
        droplet_diameters(i) = pixel_length * pixel2mm;
        end
%         disp(['The diameter of this droplet is ', num2str(droplet_diameter(i)), 'mm']);
    end
    
    droplet_diameters(droplet_diameters==0) = [];
    mean_diam = mean(droplet_diameters)
end

% i=0;
% while (i < number_of_images)
%     I2 = imread([folder_name, directory, '\', im_name, '_', sprintf('%06d',i), image_format]);
%     Idiff = rgb2gray(I1-I2);
%     % Apply mask
%     Idiff = Idiff(yi:end, xi-3*width: xi+3*width);
%     % Recovery of cm reference image to obtain pixel resolution
%     i
%     
%     % Finding of droplet(s) in the image
%     var = sum(sum((Idiff-mean(mean(Idiff))).^2))/(Im_res_x*Im_res_y);
%     
%     top_left_corner = [1,1];
%     top_left_corner2 = top_left_corner + [xi-3*width, yi];
%     while(sum(top_left_corner < size(Idiff))==2)
%         if (mean(mean(Idiff(top_left_corner:top_left_corner+window_size)>20*var)))           
%             disp('yes');
%             figure(2); hold on;
%             imshow(I2);
%             hold on;
%             plot([top_left_corner2(1) top_left_corner2(1)], [top_left_corner2(2) top_left_corner2(2)+window_size(2)],'g', 'LineWidth', 3.0);   % Left-side
%             plot([top_left_corner2(1)+window_size(1) top_left_corner2(1)+window_size(1)], [top_left_corner2(2) top_left_corner2(2)+window_size(2)],'g', 'LineWidth', 3.0);  % right-side
%             plot([top_left_corner2(1) top_left_corner2(1)+window_size(1)], [top_left_corner2(2) top_left_corner2(2)],'g', 'LineWidth', 3.0);  % upper-side
%             plot([top_left_corner2(1) top_left_corner2(1)+window_size(1)], [top_left_corner2(2)+window_size(2) top_left_corner2(2)+window_size(2)],'g', 'LineWidth', 3.0);  % lower-side
%             legend('Found droplet');
%             axis on;
%             title(['Found droplet for the ', num2str(i), '-th image']);
%         end
%         
%         top_left_corner = top_left_corner + window_size;
%     end
%     
%     % Estimation of diameter
%     
%     i=i+1;
% end






