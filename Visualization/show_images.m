close all;
clc;

%% Parameters
% Informations varying with the analyzed video sequence :
folder_name = 'C:\Users\leblanco.OASIS\Documents\IngeCivilMaster2\Memoire\Liege_Data\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory name
directory = '2_G_4_00';
image_format = '.jpg';
ref_length = 0.09;
first_image = 1;
last_image = 473;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obj = VideoReader([folder_name, directory, filename, image_format]);
% text = [folder_name, directory, '\', filename, '_'];
% Automatically find the number of image files in the folder.
number_of_images = numel(dir([folder_name, directory,'\*', image_format]))-1;    % Numbers of images in the folder or to treat
text = [folder_name, directory, '\'];

names = dir([folder_name, directory,'\*', image_format]);

if (last_image > number_of_images)
    last_image = number_of_images;
end

figure(1); hold on;

for i = first_image:last_image
    I1 = imread([text, names(i).name]);
    imshow(I1(100:end-100, 700:1700));
    axis on; hold on;
    title(['Image number ', num2str(i)], 'FontSize', 20);
    drawnow;
end
