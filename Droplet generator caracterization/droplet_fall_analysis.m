% Characterization of the droplet on-demand generator

close all;
clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = 'ajd';
% im_name = 'ajd';
% image_format = '.jpeg';
% numdrop = 20; % Number of created droplets
% first_image = 0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Options
find_pixel2mm = 0;
estim_diam = 0;
show_images = 0;
show_found_droplets = 0; %1-estim_diam

give_boxplot = 1;

%% Parameters

% window_size = [20,20];
% width = 300; % For the zoom
subsample = 4; % Step between the analyzed images
ref_length = 60; % [mm] : the reference length is the width of the gray base of the DOD
folder_name = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Droplet_Generator\Caracterization\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory name
directory = 'nozzle3_pulse5000_voltage50';
im_name = 'drop';
image_format = '.jpeg';
numdrop = 20; % Number of created droplets
first_image = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Automatically find the number of files in the folder.
% Numbers of images in the folder or to treat.
number_of_images = numel(dir([folder_name, directory,'\*', image_format]))-1;    % Numbers of images in the folder or to treat

text = [folder_name, directory, '\', im_name, '_'];

%% Code

if (find_pixel2mm)
    % ------------------------------------------------------------------------------
    % Let the user click on the base of the DOD to estimate the image's pixel size
    I1 = imread([folder_name, directory, '\', im_name, '_', sprintf('%06d',first_image), image_format]);
    [Im_res_y, Im_res_x, ~] = size(I1); % Finds the orientation and resolution of the images.
    figure(1); hold on;
    imshow(I1);
    axis on;
    title('\color{red} please click on the left side of the gray base', 'FontSize', 20);
    [xl, yl] = ginput(1);
    title('\color{red} please click on the right side of the gray base', 'FontSize', 20);
    [xr, yr] = ginput(1);
    pixel2mm = ref_length/(xr-xl);
    
    title('\color{red} Now select your ROI', 'FontSize', 20);
    rect = round(getrect());
    y_ul = rect(1);
    x_ul = rect(2);
    wy = rect(3);
    lx = rect(4);
    
    ROI_indx = x_ul:x_ul+lx;
    ROI_indy = y_ul:y_ul+wy;
    
    close(figure(1));
    % ------------------------------------------------------------------------------
end

if (estim_diam)
    % Preallocation
    image_with_drop = zeros(2*numdrop, 1);
    estim_diameter = zeros(2*numdrop, 1);
    
    % Use morphological operators to remove salt noise
    square_22 = [1 1; 1 1];
    square_33 = [1 1 1; 1 1 1; 1 1 1];
    square_44 = [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1];
    rect_42 = [1 1; 1 1; 1 1; 1 1];
    line2 = [1 1 1];
    
    kernel = square_22;
    
    % Process images to estimate droplets' diameters
    I1 = imread([folder_name, directory, '\', im_name, '_', sprintf('%06d',first_image), image_format]);
    I1 = I1(ROI_indx, ROI_indy, :);
    
    fig1 = figure(1); hold on;
    set(fig1, 'position',[0 200 500 500]);
    
    i = 1;
    for ind=first_image:subsample:number_of_images
        % for ind=2:subsample:150
%             for ind=324:subsample:324
        
        I2 = imread([folder_name, directory, '\', im_name, '_', sprintf('%06d',ind), image_format]);
        I2 = I2(ROI_indx, ROI_indy, :);
        I3 = rgb2gray(abs(I1-I2));
        
        % Convert the difference image in binary image
        line_I3 = double(reshape(I3, numel(I3), 1));
        mu = mean(line_I3);
        std_I3 = std(line_I3);
        log = logical(double(I3>3*std_I3).*double(I3>3));
        I4 = I3;
        I4(log) = 255;
        I4(I4~=255) = 0;
        
        % Remove salt and pepper noise
        I5 = remove_salt_and_pepper_noise(I4, kernel);
        I5 = remove_salt_and_pepper_noise(I5, rect_42);
        I5 = I5./255;
        
        % We found a droplet!
        if (sum(sum(I5))>1000)
            image_with_drop(i) = ind;
            % Compute widths
            horiz_sum = sum(I5')';
            % Remove outliers
            horiz_sum(horiz_sum==0) = [];
            mu_horiz_sum = mean(horiz_sum);
            std_horiz_sum = std(horiz_sum);
            horiz_sum(horiz_sum>mu_horiz_sum +std_horiz_sum/2) = [];
            horiz_sum(horiz_sum<mu_horiz_sum -std_horiz_sum/2) = [];
            % Estimate the diameter
            pixel_length = mean(horiz_sum);
            estim_diameter(i) = pixel_length * pixel2mm;
            i=i+1;
        end
        
        %                    figure(2); hold on;
        %                    hist(double(line_I3), 100);
        
        if (show_images)
            % Plots
            subplot(1,4,1); hold on;
            title('Original');
            imshow(I2);
            axis on;
            subplot(1,4,2); hold on;
            title('Difference with ref.');
            imshow(I3.*5);
            axis on;
            subplot(1,4,3); hold on;
            title('Binary image');
            imshow(I4);
            axis on;
            subplot(1,4,4); hold on;
            title('Wo salt and pepper');
            imshow(I5.*255);
            axis on;
            set(fig1, 'Name',['Analyzing ', num2str(ind), '-th image']);
            drawnow;
        end
    end
    
    axis on;
    image_with_drop(image_with_drop==0) = [];
    estim_diameter(estim_diameter==0) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (show_found_droplets)
    fig3 = figure(3); hold on;
    set(fig3, 'position',[800 200 500 500]);
    for i=1:length(image_with_drop)
        I2 = imread([folder_name, directory, '\', im_name, '_', sprintf('%06d',image_with_drop(i)), image_format]);
        I2 = I2(ROI_indx, ROI_indy, :);
        I3 = rgb2gray(abs(I1-I2));
        
        % Convert the difference image in binary image
        line_I3 = double(reshape(I3, numel(I3), 1));
        mu = mean(line_I3);
        std_I3 = std(line_I3);
        log = logical(double(I3>3*std_I3).*double(I3>3));
        I4 = I3;
        I4(log) = 255;
        I4(I4~=255) = 0;
        
        % Remove salt and pepper noise
        I5 = remove_salt_and_pepper_noise(I4, kernel);
        I5 = remove_salt_and_pepper_noise(I5, rect_42);
        I5 = I5./255;
        
        % Plots
        subplot(1,4,1); hold on;
        title('Original');
        imshow(I2);
        axis on;
        subplot(1,4,2); hold on;
        title('Difference with ref.');
        imshow(I3.*5);
        axis on;
        subplot(1,4,3); hold on;
        title('Binary image');
        imshow(I4);
        axis on;
        subplot(1,4,4); hold on;
        title('Wo salt and pepper');
        imshow(I5.*255);
        axis on;
        set(fig3, 'Name',['Analyzing ', num2str(image_with_drop(i)), '-th image']);
        drawnow;
        pause(1.5);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (give_boxplot)
    myFiles = dir([folder_name, '\*', '.mat']);
    number_of_mat_files = numel(dir([folder_name, '\*', '.mat']))-1;
    figure(2); hold on;
    arr = 0;
    x = 0;
    
    for k = 1:number_of_mat_files
        baseFileName = myFiles(k).name;
        disp([num2str(k), ' : ', baseFileName]);
        fullFileName = fullfile(folder_name, baseFileName);
        arr2 = cell2mat(struct2cell(load(fullFileName)));
        arr = [arr; arr2];
        x = [x; k*ones(length(arr2),1)];
    end
    arr(1)=[];
    x(1) = [];
    boxplot(arr, x); hold on;
    xlabel('Test index');
    ylabel('Droplet diameter [mm]');
    ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dilation = dilate(image, kernel)
dilation = imcomplement(imfilter(imcomplement(image),kernel,'conv'));
%     dilation(dilation>0) = 1;
end

function erosion = erode(image, kernel)
erosion = imfilter(image,kernel,'conv');
%     erosion(erosion>0) = 1;
end

function closed = closing(image, kernel)
closed = dilate( erode(image, kernel), kernel );
end

function opened = opening(image, kernel)
opened = erode( dilate(image, kernel), kernel );
end

function cleaned = remove_salt_and_pepper_noise(image, kernel)
cleaned = opening( closing(image, kernel), kernel);
end


