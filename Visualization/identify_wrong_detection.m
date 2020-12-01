%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Olivier Leblanc
% Date : 01/11/2020
%
% Script :
% Let the user identify a wrong position detection and replace it by the
% mean between its two neighbours in time. 
% The position is graphically selected by drawing a polygon around it.
%
%
% Examples :
%
close all;
clc;
%
% Options :
modify_pos = 1; % Put 1 for changing "drop_position" variable. 0 otherwise.
% Parameters :
directory = 'C:\Users\leblanco.OASIS\Documents\IngeCivilMaster2\Memoire\Liege_Data\J7_20_11\';
im_folder = '2_G_4_30_B';
% Note :
% The path names here below may be slightly modified depending on the names
% of the position variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos = load([directory, im_folder, '_positions.mat']);
text = [directory, im_folder];
imnames = dir([text, '\*.jpg']);
I = imread([text, '\', imnames(1).name]);


fig1 = figure(1); hold on;
set(fig1, 'position',[0 0 800 800]);
imshow(I); hold on;
X = pos.drop_position(:,1);
Y = pos.drop_position(:,2);
p = plot(X , Y , 'r', 'LineWidth', 1.0); 
n = length(X);
cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
drawnow;
set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
axis on;
axis([0,2050,0,2050]);

title(['\color{red} Draw polygon around point'], 'FontSize', 20);
cutpolygon = drawpolygon();
pos_polygon = cutpolygon.Position;
[in, on] = inpolygon(X, Y, pos_polygon(:,1), pos_polygon(:,2));
ind = find(in)

close(fig1);

X(ind) = (X(ind-1)+X(ind+1))/2;
Y(ind) = (Y(ind-1)+Y(ind+1))/2;

figure(1); 
imshow(I); hold on;
p = plot(X , Y , 'r', 'LineWidth', 1.0); 
n = length(X);
cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
drawnow;
set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
axis on;
axis([0,2050,0,2050]);

if (modify_pos)
   drop_position = [X,Y];
end

