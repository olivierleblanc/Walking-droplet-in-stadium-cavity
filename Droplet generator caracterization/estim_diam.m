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
define_resolution = 1;
show_charac = 0; % Show the characterization result for different droplets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% analyze_video = 0;  % 0 for analyzing a directory of images
% folder = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Experiments\24-04-20_Me5_9';
% filename = 'drop';
% format = '.jpg';
% ref_length = 0.09;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% analyze_video = 1;  % 0 for analyzing a directory of images
% folder = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Experiments\24-04-20_Me6_1';
% filename = 'P1060751';
% format = '.MOV';
% ref_length = 0.085;
% list = 3000:50:4000;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory name
analyze_video = 0;  % 0 for analyzing a directory of images
folder = 'C:\Users\leblanco.OASIS\Documents\IngeCivilMaster2\Memoire\Liege_Data\J2_23-10\2_G_4_00\';
filename = 'P1060751';
format = '.jpg';
ref_length = 0.095;
list = 100:50:1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Num_droplets = 1;
text = [folder, filename, '_'];
if (analyze_video)
    obj = VideoReader([folder, filename, format]);
    number_of_images = length(list);
else
    names = dir([folder, '*', format]);
    number_of_images = numel(dir([folder, '*', format]))-1;    % Numbers of images in the folder or to treat
end

f1 = figure(1); hold on;
f1 = figure(1); hold on;

% ------------------------------------------------------------------------------
if (define_resolution)
    % Let the user click on some reference point to estimate the image's pixel size
    if (analyze_video)
        I1 = read(obj, 1);
    else
        I1 = imread([folder, names(1).name]);
    end
    [Im_res_y, Im_res_x, ~] = size(I1); % Finds the orientation and resolution of the images.
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
end
% ------------------------------------------------------------------------------

diams = zeros(1, number_of_images);

for i = 1:number_of_images
    j = list(i);
    % ------------------------------------------------------------------------------
    % Let the user select the droplet of interest.
    % Define centers of initial regions of interest (ROI)
    if (analyze_video)
        I1 = read(obj, j);
    else
        I1 = imread([folder, names(i).name]);
    end
    figure(1);
    imshow(I1);
%     colorbar;
    axis on;
    title('\color{red} please click on the droplet to define its initial coordinates', 'FontSize', 20);
    [xi, yi] = ginput(1);
    ROI_centers = [round(yi); round(xi)];
    % ------------------------------------------------------------------------------
    
    % Strucure containing y ROI values on the first row, and x ROI values on
    % the second row, repeated Num_droplets times.
    ROI = repmat([[ROI_centers(1,1)-half_ROI_size+1 : ROI_centers(1,1)+half_ROI_size] ;...
        [ROI_centers(2,1)-half_ROI_size+1 : ROI_centers(2,1)+half_ROI_size]],1,Num_droplets);
    corrzone = repmat([[ROI_centers(1,1)-half_corr_size+1 : ROI_centers(1,1)+half_corr_size-1] ;...
        [ROI_centers(2,1)-half_corr_size+1 : ROI_centers(2,1)+half_corr_size-1]],1,Num_droplets);
    
    % ------------------------------------------------------------------------------
    % Let the user click on two sides of the droplet to estimate its diameter
    if (analyze_video)
        I1 = read(obj, j);
    else
        I1 = imread([folder, names(i).name]);
    end 
    I1_corr = I1([corrzone(1,:)],[corrzone(2,:)]);
    % Interpolate the image by 'factor'
    im = double(I1_corr);
    factor = 20;
    imup = interpft(interpft(im, factor*size(im,1), 1), factor*size(im,2), 2);
    rescaled = uint8(interp1([min(imup(:)), max(imup(:))], [min(im(:)), max(im(:))], imup));
    figure(2);
    imshow(rescaled);
    axis on;
    title('\color{red} please click on one side of the droplet', 'FontSize', 20);
    [xi,yi] = getline;
    diams(i) = sqrt( (xi(1)-xi(2))^2+(yi(1)-yi(2))^2 )*pixel2mm/factor;
    meandiams = mean(diams);
    vardiams = mean((diams-meandiams).^2);
    sigdiams = sqrt(vardiams);
    % ------------------------------------------------------------------------------
end
close(figure(1));
close(figure(2));

disp([num2str(number_of_images), ' images have been analyzed.']);
disp(['The average over the droplets diameters estimations is : ', num2str(meandiams), ' [m]']);

figure(1); hold on; 
plot(1e3.*diams, 'r'); 
plot(1e3.*meandiams.*ones(1,length(diams)), 'b');
xlabel('Index [-]', 'interpreter', 'latex', 'FontSize', 12.0);
ylabel('Diameter [mm]', 'interpreter', 'latex', 'FontSize', 12.0);
legend('diameter estimations', 'mean', 'interpreter', 'latex', 'Location', 'NorthWest', 'FontSize', 14.0)
set(gca,'TickLabelInterpreter','Latex', 'FontSize', 14.0);

if (show_charac)
    load('C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Droplet_Generator\Characterization_In_Bath\droplets.mat');
    figure(2); hold on;
    plot(droplets, 'g');
    plot(mean(droplets).*ones(1,length(droplets)), 'm');
    xlabel('Droplet number [-]', 'interpreter', 'latex', 'FontSize', 12.0);
    ylabel('Diameter [mm]', 'interpreter', 'latex', 'FontSize', 12.0);
    legend('diameter estimations', 'mean', 'interpreter', 'latex', 'Location', 'NorthEast', 'FontSize', 14.0)
    set(gca,'TickLabelInterpreter','Latex', 'FontSize', 14.0);
end



