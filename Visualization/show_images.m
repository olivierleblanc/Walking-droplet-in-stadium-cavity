%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Olivier Leblanc
% Date : 19/11/2020
%
% Script :
% Allows to show the images contained in a directory one by one.
%
% Inputs :
% /.
% Outputs :
% /.
% Examples :
% /.
%
% Options :
%
% Parameters :
time_delay = 0.2; % Between two frames
folder_name = 'C:\Users\leblanco.OASIS\Documents\IngeCivilMaster2\Memoire\Liege_Data\J6_18_11\';
directory = '2_G_4_00';
image_format = '.jpg';
first_image = 1;
last_image = 1e5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc;

text = [folder_name, directory, '\'];
names = dir([text,'*', image_format]);
number_of_images = numel(names)-1;  % Numbers of images in the folder or to treat
if (last_image > number_of_images)
    last_image = number_of_images;
end

names = dir([folder_name, directory,'\*', image_format]);

if (last_image > number_of_images)
    last_image = number_of_images;
end

figure(1); hold on;

for i = first_image:last_image
    I1 = imread([text, names(i).name]);
    imshow(I1);
    axis on; hold on;
    title(['Image number ', num2str(i)], 'FontSize', 20);
    drawnow;
    pause(time_delay);
end
