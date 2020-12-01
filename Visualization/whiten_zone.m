function [] = whiten_zone (path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Olivier Leblanc
% Date : 01/11/2020
%
% Function :
% Allows the user to whiten some artifacts in an image.
% The artifact zone is selected by drawing a polygon.
%
% Inputs :
% - path : path to the directory containing all images to be modified
%
% Outputs :
% /.
%
% Options :
test = 1;   % If 1, applies on the first image of the directory only. 
            % If 0, applies on all images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    image_format = '.jpg';
    number_of_images = numel(dir([path,'*', image_format]))-1;  % Numbers of images in the folder or to treat
    names = dir([path,'*', image_format]);
    
    
    I = imread([path, names(1).name]);
    [H, W, ~] = size(I);
    matx = repmat(1:W,H,1);
    maty = repmat(transpose(1:H),1,W);
    
    
    figure(13); imshow(I); hold on;
    title(['\color{red} Draw the polygon around the artifact'], 'FontSize', 20);
    cutpolygon = drawpolygon();
    pos_polygon = cutpolygon.Position;
    [in, on] = inpolygon(matx, maty, pos_polygon(:,1), pos_polygon(:,2));
    
    if(test)
       number_of_images = 1; 
    end
    
    for i = 1 : number_of_images
        i
        filename = [path, names(i).name];
        % Read images
        I = imread(filename);
                
        I(in) = 255; 
        imwrite(I,filename);
    end
    
    close(figure(13));

end

