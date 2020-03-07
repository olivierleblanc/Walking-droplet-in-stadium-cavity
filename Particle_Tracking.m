% Particle tracking
% Author : Olivier Leblanc
% Date : 19-02-2020


% Explains the principle.


close all;
clc;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = 'Tracking2';
% im_name = 'track2';
% image_format = '.jpeg';
% resolution = 60.0e-6;   % 1 pixel = a micrometer
% first_image = 100;
% last_image = 1500;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = '25-11-19_droplet1';
% im_name = 'droplet1';
% image_format = '.jpeg';
% resolution = 150.0e-6;   % 1 pixel = a micrometer
% first_image = 434;
% last_image = 1500;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = '26-11-19_circle';
% im_name = 'circle';
% image_format = '.jpeg';
% resolution = 24.3e-6;   % 1 pixel = a micrometer
% first_image = 1;
% last_image = 1500;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Note : Pour le test ci-dessous, on remarque que le tracking jump de la
% grosse à la petite goutte autour de "index_test"=169!

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = '03-12-19-DeuxGouttes';
% im_name = 'deux';
% image_format = '.jpeg';
% resolution = 3.7037e-05;   % 1 pixel = a micrometer
% first_image = 1;
% last_image = 1500;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot options

show_true_trajectory = 1;
show_tracked_trajectory = 1;
show_verify_position = 1;
show_final_ROI = 0;

%% Options
subsample = 10;  % No if 1. Put 10 for test.

%% Useful values

vmax_drop = 0.01;   % 1cm/s
displ_max = 2e-4;   % Maximum displacement of a droplet between two frames

%% Parameters
Num_droplets = 1; % Number of droplets to track
Intensity_Threshold = 40;   % Found by trials

% Informations varying with the analyzed video sequence :

folder_name = 'C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Visualization\Images\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory name
directory = '21-02-20_stadium_random';
im_name = 'random';
image_format = '.jpeg';
resolution = 1.3889e-04;   % 1 pixel = a micrometer
first_image = 1;
last_image = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

text = [folder_name, directory, '\', im_name, '_'];

% ------------------------------------------------------------------------------
% Let the user select the droplet of interest.
% Define centers of initial regions of interest (ROI)
figure(10); hold on;
I1 = imread([folder_name, directory, '\', im_name, '_', sprintf('%06d',first_image), image_format]);
I2 = imread([folder_name, directory, '\', im_name, '_', sprintf('%06d',first_image+4), image_format]);
I1=rgb2gray(I1-I2);
[Im_res_y, Im_res_x, ~] = size(I1); % Finds the orientation and resolution of the images.
imshow(5.*I1);
colorbar;
axis on;
title('\color{red} please click on the droplet to define its initial coordinates', 'FontSize', 20);
[xi, yi] = ginput(1);
ROI_centers = [round(yi); round(xi)];
close(figure(10));
% ------------------------------------------------------------------------------

Image_dist = 7; % Initial spacing (in images numbers) between the actual reference image and the analyzed image

ROI_size = round(Image_dist*displ_max/resolution);
half_ROI_size = round(ROI_size/2);
shift_x = round((Im_res_x - ROI_size)/2);
shift_y = round((Im_res_y - ROI_size)/2);

% Automatically find the number of files in the folder.
% Numbers of images in the folder or to treat.
number_of_images = numel(dir([folder_name, directory,'\*', image_format]))-1;    % Numbers of images in the folder or to treat

% Strucure containing y ROI values on the first row, and x ROI values on
% the second row, repeated Num_droplets times.
ROI = repmat([[ROI_centers(1,1)-half_ROI_size+1 : ROI_centers(1,1)+half_ROI_size-1] ;...
    [ROI_centers(2,1)-half_ROI_size+1 : ROI_centers(2,1)+half_ROI_size-1]],1,Num_droplets);

%% Code

if (show_tracked_trajectory)

    % Initialization
    drop_position = zeros(Num_droplets*last_image,2);
    image_ind = zeros(last_image);
    
    tic;
    
    i=first_image;
    while (i+Image_dist < last_image) 
        i
        droplet_found = 0;
        % Refresh the reference image for next image to analyze
        Iref = imread([text, sprintf('%06d',i), image_format]);
        
        while (~droplet_found)
            
            if (i+Image_dist>number_of_images)
                break;
            end
            % Load the new image
            I2 = imread([text, sprintf('%06d',i+Image_dist), image_format]); 
            % absolute value is not a good idea because otherwise the absence is as much important as the presence.
            I3 = I2-Iref;    % Compute the difference
            I3 = rgb2gray(I3);  % Convert into grayscale
            
            % Finds droplets' new positions
            zone = I3([ROI(1,:)],[ROI(2,:)]);
            maxima = max(max(zone));
            
            % Verifies that the maximum value corresponds to a droplet,
            % otherwise, tries with the next frame.
            if (maxima(1)>=Intensity_Threshold)
                Image_dist
                %             disp('yes');
                [y,x] = find(I3([ROI(1,:)],[ROI(2,:)])==maxima(1));
                y = y(1) + shift_y;
                x = x(1) + shift_x;
                drop_position(i+Image_dist,:) = [x,y];
                image_ind(i+Image_dist) = i+Image_dist;
                %             Image_dist
                droplet_found = 1;
                i = i+Image_dist;   % First image where all droplets are found becomes new reference image
                Image_dist = subsample;
                ROI_size = round(Image_dist*displ_max/resolution); % Width and height of ROI depending on time distance between the two images
                half_ROI_size = round(ROI_size/2);
                % Updates ROI_centers
                ROI_centers = [y; x];
                shift_y = ROI_centers(1,:)-half_ROI_size;
                shift_x = ROI_centers(2,:)-half_ROI_size;
                ROI = repmat([[ROI_centers(1,1)-half_ROI_size+1 : ROI_centers(1,1)+half_ROI_size-1] ;...
                    [ROI_centers(2,1)-half_ROI_size+1 : ROI_centers(2,1)+half_ROI_size-1]],1,Num_droplets);
                
            else
                %             disp('try again')
                Image_dist = Image_dist+subsample;
                % Increases ROI size
                ROI_size = round(Image_dist*displ_max/resolution); % Width and height of ROI depending on time distance between the two images
                half_ROI_size = round(ROI_size/2);
                shift_y = ROI_centers(1,:)-half_ROI_size;
                shift_x = ROI_centers(2,:)-half_ROI_size;
                ROI = repmat([[ROI_centers(1,1)-half_ROI_size+1 : ROI_centers(1,1)+half_ROI_size-1] ;...
                    [ROI_centers(2,1)-half_ROI_size+1 : ROI_centers(2,1)+half_ROI_size-1]],1,Num_droplets);               
            end
            
            % Cut the ROI if goes out the image
            vec1 = [find(ROI(1,:)<=0)];
            ROI(:,vec1) =[];
            vec2 = [find(ROI(2,:)<=0)];
            ROI(:,vec2)=[];
            shift_x=shift_x+length(vec1)+length(vec2);  % Ajout un shift pour tenir compte des elements supprimes
            shift_y=shift_y+length(vec1)+length(vec2);
        end
        
    end
    
    toc;
    
    drop_position(find(drop_position(:,1)==0), :) = [];
    image_ind(find(image_ind==0)) = [];
end

if (show_true_trajectory)
    tic;  
    % Shows the mean image between the first 'nb_for_mean' images
    space_for_mean = 20;
    nb_for_mean = (last_image/(space_for_mean+1));
    
    I1 = double(imread([text, sprintf('%06d', first_image), image_format]));
    Ipre = double(imread([text, sprintf('%06d', first_image), image_format]));
    for n = 1:nb_for_mean
        Inew = double(imread([text, sprintf('%06d', first_image+space_for_mean*n), image_format]));
        I_diff = Inew-Ipre;
        Ipre = Inew;
        % Mettre un facteur qui dépende de l'intensité moyenne!!
        %         weights = I_diff./(100/nb_for_mean);
        weights = I_diff./(1000/nb_for_mean);
        I1 = I1 + Inew.*(1 + weights);
    end
    I1 = uint8(round(I1./nb_for_mean));
    toc; 
else
    I1 = uint8(imread([text, sprintf('%06d', first_image), image_format]));
end

%% Visualization

if (show_tracked_trajectory)
    figure(1); hold on;
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

if (show_final_ROI)
    % Shows the ROI at last frame
    fig2 = figure(2); hold on;
    set(fig2, 'position',[200 200 500 500]);
    imshow(I2); hold on;
    plot([ROI(2,1), ROI(2,1)], [ROI(1,1), ROI(1,end)],'b', 'LineWidth', 2.0);   % Left-side
    plot([ROI(2,end), ROI(2,end)], [ROI(1,1), ROI(1,end)],'b', 'LineWidth', 2.0);  % right-side
    plot([ROI(2,1), ROI(2,end)], [ROI(1,1), ROI(1,1)],'b', 'LineWidth', 2.0);  % upper-side
    plot([ROI(2,1), ROI(2,end)], [ROI(1,end), ROI(1,end)],'b', 'LineWidth', 2.0);  % lower-side
    legend('Region of interest (ROI)');
    axis on;  
end

% Shows the considered ROI at the 'index_test' frame number, allows to
% manually verify the tracking.
if (show_verify_position)
    n_iter = 10;
    
    fig3 = figure(3); hold on;
    set(fig3, 'position',[600 200 500 500]);
    axis on;
    
    iter = 1;
    while iter < n_iter
        index_test = round( (first_image + round(length(image_ind)*iter/n_iter)));
        I3 = imread([text, sprintf('%06d', image_ind(index_test)), image_format]);
        imshow(I3); hold on;
        title(['Found position for the ', num2str(image_ind(index_test)), '-th image']);
        scatter(drop_position(index_test,1),drop_position(index_test,2), 'b');
        iter = iter+1;
        pause(2);
    end
    
%     fig4 = figure(4); hold on;
%     set(fig4, 'position',[400 50 500 500]);
%     axis on;
%     index_test = 2;
%     I3 = imread([text, sprintf('%06d', image_ind(index_test)), image_format]);
%     imshow(I3); hold on;
%     title(['Found position for the ', num2str(image_ind(index_test)), '-th image']);
%     scatter(drop_position(index_test,1),drop_position(index_test,2), 'b');
end
