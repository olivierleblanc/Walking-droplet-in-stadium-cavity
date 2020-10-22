% Particle tracking
% Author : Olivier Leblanc
% Date : 19-02-2020


% Explains the principle.


close all;
clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = 'Tracking2';
% image_format = '.jpeg';
% resolution = 60.0e-6;   % 1 pixel = a micrometer
% first_image = 100;
% last_image = 1500;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = '25-11-19_droplet1';
% image_format = '.jpeg';
% resolution = 150.0e-6;   % 1 pixel = a micrometer
% first_image = 434;
% last_image = 1500;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Directory name
% directory = '26-11-19_circle';
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
% image_format = '.jpeg';
% resolution = 3.7037e-05;   % 1 pixel = a micrometer
% first_image = 1;
% last_image = 1500;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot options

show_true_trajectory = 0;
show_tracked_trajectory = 1;
show_verify_position = 0;
show_final_ROI = 0;

%% Options
subsample = 1;  % No if 1. Put 10 for test.
define_resolution = 0;
define_ROI = 0;

%% Parameters
Num_droplets = 1; % Number of droplets to track
Intensity_Threshold = 40;   % Found by trials
halfwidth = 80; % For the ROI size

% Informations varying with the analyzed video sequence :
folder_name = 'C:\Users\leblanco.OASIS\Documents\IngeCivilMaster2\Memoire\Liege_Data\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directory = '2_G_4_00';
image_format = '.jpg';
image_color = 0;
ref_length = 0.095;
first_image = 1;
last_image = 200;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

text = [folder_name, directory, '\'];
names = dir([folder_name, directory,'\*', image_format]);
number_of_images = numel(dir([text,'*', image_format]))-1;  % Numbers of images in the folder or to treat
if (last_image > number_of_images)
    last_image = number_of_images;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (define_resolution)
    figure(1); hold on;
    imshow(imread([text, names(first_image).name]));
    axis on;
    title('\color{red} please click on the left side of the gray base', 'FontSize', 20);
    [xl, yl] = ginput(1);
    title('\color{red} please click on the right side of the gray base', 'FontSize', 20);
    [xr, yr] = ginput(1);
    pixel2mm = ref_length/(sqrt((xr-xl)^2+(yr-yl)^2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (define_ROI)
    % Define centers of initial regions of interest (ROI)
    figure(10); hold on;
    I1 = imread([text, names(first_image).name]);
    I2 = imread([text, names(first_image+1).name]);
    if (image_color)
        I1=rgb2gray(I1-I2);
    else
        I1 = I1-I2;
    end
    [Im_res_y, Im_res_x, ~] = size(I1); % Finds the orientation and resolution of the images.
    imshow(5.*I1);
    colorbar;
    axis on;
    title('\color{red} Select the ROI', 'FontSize', 20);
    rect = round(getrect());
    x_ul = rect(1);
    y_ul = rect(2);
    wx = rect(3);
    hy = rect(4);
    N = wx*hy;

    ROI_indx = x_ul:x_ul+wx;
    ROI_indy = y_ul:y_ul+hy;

    imshow(5.*I1(ROI_indy,ROI_indx));
    axis on;
    title('\color{red} please click on the droplet to define its initial coordinates', 'FontSize', 20);
    [xi, yi] = ginput(1);
    ROI_centers_init = [round(y_ul+yi); round(x_ul+xi)];
    close(figure(10));
end

ROI_centers = ROI_centers_init;
ROIx = max(1,ROI_centers_init(2)-halfwidth) : min(2040,ROI_centers_init(2)+halfwidth);
ROIy = max(1,ROI_centers_init(1)-halfwidth) : min(2048,ROI_centers_init(1)+halfwidth);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use morphological operators to remove salt noise
square_22 = [1 1; 1 1];
square_33 = [1 1 1; 1 1 1; 1 1 1];
square_44 = [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1];
rect_42 = [1 1; 1 1; 1 1; 1 1];
line3 = [1 1 1];

kernel = square_22;

matx1 = repmat(1:Im_res_x,Im_res_y,1);
maty1 = repmat(transpose(1:Im_res_y),1,Im_res_x);

% fig10 = figure(10); hold on;
% set(fig10, 'position',[1000 0 500 500]);
% colorbar; axis on; hold on;
% 
% fig11 = figure(11); hold on;
% set(fig11, 'position',[500 0 500 500]);
% colorbar; axis on; hold on;
% 
% fig12 = figure(12); hold on;
% set(fig12, 'position',[0 0 500 500]);
% colorbar; axis on; hold on;

fig13 = figure(13); hold on;
set(fig13, 'position',[500 1000 1000 1000]);
% legend('Brightest pixel in the image difference','Gravity center of the mask','Most correlated point','Darkest point inside the ROI');

imref = imread([text, names(21).name]);
imref = imref(180:380, 1010:1220);
H = length(imref(:,1));
W = length(imref(1,:));
[X,Y] = meshgrid(1:H, 1:W);
% figure(13); hold on; imshow(imref); axis on;
imref2 = double(imref);
mref = mean(mean(imref2));
imref_window = imref2-mref;
fourier_imref_corr = fftshift(fft2(imref_window));
sigma = 0.001;

drop_position = zeros(Num_droplets*last_image,2);

for i = first_image : first_image+300
    I1 = imread([text, names(i).name]);
    I2 = imread([text, names(i+1).name]);
    if (image_color)
        I1=rgb2gray(I1-I2);
    else
        I1 = I1-I2;
    end

    % Compute the maximum difference value
    [val y] = max(I1(ROIy, ROIx));
    [val2 x] = max(val);
    y = y(x);
    
    % Compute the mass center of the mask
    mask = ones( Im_res_y, Im_res_x );
    mask(5.*I1<30)=0;
    mask = remove_salt_and_pepper_noise(mask, kernel);
    mask = remove_salt_and_pepper_noise(mask, rect_42);
    mask = remove_salt_and_pepper_noise(mask, rect_42');
    mask([1:ROIy(1)-1, ROIy(end)+1:end], [1:ROIx(1)-1, ROIx(end)+1:end]) = 0;
    matx = matx1.*mask;
    maty = maty1.*mask;
    matx(~mask)=[];
    maty(~mask)=[];
    xm = round(mean(mean(matx)));
    ym = round(mean(mean(maty)));
    indy = max(1,ym-halfwidth):min(2048,ym+halfwidth);
    indx = max(1,xm-halfwidth):min(2040,xm+halfwidth);
    
    % Compute the correlation with imref
    I2b = double(I2(indy,indx));
    m2 = mean(mean(I2b));
    I2_window = I2b-m2;
    fourier_I2_corr = fftshift(fft2(I2_window,H,W));
    fourier_prod = fourier_I2_corr.*conj(fourier_imref_corr);
    corr = fftshift(real(ifft2(ifftshift(fourier_prod))));
    corr = rot90(corr);
    [val xcorr] = max(abs(corr));
    [val2 ycorr] = max(val);
    xcorr = xcorr(ycorr);
    
    % Compute the darkest pixel in the ROI
%     [val yb] = min(double(I2(indy,indx)).*exp((matx1(indy,indx)-x).^2+(maty1(indy,indx)-y).^2));
    [val yb] = min(I2(indy,indx));
    [val2 xb] = min(val);
    yb = yb(xb);
   
    %%%%%%%%%%%%% Figures %%%%%%%%%%%%%%
%     figure(14); hold on;
%     s = surf(X,Y,corr);
%     s.EdgeColor = 'none';
%     scatter3(ycorr, xcorr, val2, 'k');
%     
%     figure(10);
%     imshow(5.*I1); hold on;
%     scatter(x+ROIx(1),y+ROIy(1),'b');
%     scatter(xm,ym,'m');
%     scatter(x+ROIx(1)-H+xcorr, y+ROIy(1)-ycorr,'g');
%     scatter(indx(1)+xb,indy(1)+yb,'r');
%     plot([indx(1), indx(1)], [indy(1), indy(end)],'b', 'LineWidth', 2.0);   % Left-side
% 	plot([indx(end), indx(end)], [indy(1), indy(end)],'b', 'LineWidth', 2.0);  % right-side
% 	plot([indx(1), indx(end)], [indy(1), indy(1)],'b', 'LineWidth', 2.0);  % upper-side
% 	plot([indx(1), indx(end)], [indy(end), indy(end)],'b', 'LineWidth', 2.0);  % lower-side
%     plot([ROIx(1), ROIx(1)], [ROIy(1), ROIy(end)],'r', 'LineWidth', 2.0);   % Left-side
% 	plot([ROIx(end), ROIx(end)], [ROIy(1), ROIy(end)],'r', 'LineWidth', 2.0);  % right-side
% 	plot([ROIx(1), ROIx(end)], [ROIy(1), ROIy(1)],'r', 'LineWidth', 2.0);  % upper-side
% 	plot([ROIx(1), ROIx(end)], [ROIy(end), ROIy(end)],'r', 'LineWidth', 2.0);  % lower-side
%     hold on; axis on;
%     
%     figure(11); imshow(mask); scatter(xm,ym,'m'); axis on;
%     
%     figure(12); 
%     imshow(I2); hold on;
%     scatter(x+ROIx(1),y+ROIy(1),'b');
%     scatter(xm,ym,'m');
%     scatter(x+ROIx(1)-(xcorr-round(H/2))+40, ycorr+y+ROIy(1)-round(W/2)+30,'g');
%     scatter(indx(1)+xb,indy(1)+yb,'r');
%     plot([indx(1), indx(1)], [indy(1), indy(end)],'b', 'LineWidth', 2.0);   % Left-side
% 	plot([indx(end), indx(end)], [indy(1), indy(end)],'b', 'LineWidth', 2.0);  % right-side
% 	plot([indx(1), indx(end)], [indy(1), indy(1)],'b', 'LineWidth', 2.0);  % upper-side
% 	plot([indx(1), indx(end)], [indy(end), indy(end)],'b', 'LineWidth', 2.0);  % lower-side
%     plot([ROIx(1), ROIx(1)], [ROIy(1), ROIy(end)],'r', 'LineWidth', 2.0);   % Left-side
% 	plot([ROIx(end), ROIx(end)], [ROIy(1), ROIy(end)],'r', 'LineWidth', 2.0);  % right-side
% 	plot([ROIx(1), ROIx(end)], [ROIy(1), ROIy(1)],'r', 'LineWidth', 2.0);  % upper-side
% 	plot([ROIx(1), ROIx(end)], [ROIy(end), ROIy(end)],'r', 'LineWidth', 2.0);  % lower-side
%     axis on;
%     title(['N°',num2str(i)]);
    
%     figure(16); hold on;
    imshow(I2(indy,indx)); hold on;
%     scatter(x-indx(1)+ROIx(1),y-indy(1)+ROIy(1),'b');
    scatter(xm-ROIx(1),ym-ROIy(1),'m');
%     scatter(xcorr, ycorr,'g');
    scatter(xb,yb,'r');
    axis on;
    drawnow;
    
    drop_position(i,:) = [xb+ROIx(1),yb+ROIy(1)];
    
    ROIx = max(1,xb+indx(1)-halfwidth) : min(2040,xb+indx(1)+halfwidth);
    ROIy = max(1,yb+indy(1)-halfwidth) : min(2048,yb+indy(1)+halfwidth);
end

drop_position(find(drop_position(:,1)==0), :) = [];




% Image_dist = 1; % Initial spacing (in images numbers) between the actual reference image and the analyzed image
% 
% % ROI_size = round(Image_dist*displ_max/pixel2mm);
% % half_ROI_size = round(ROI_size/2);
% % shift_x = round((Im_res_x - ROI_size)/2);
% % shift_y = round((Im_res_y - ROI_size)/2);
% 
% % Strucure containing y ROI values on the first row, and x ROI values on
% % the second row, repeated Num_droplets times.
% ROI = repmat([[ROI_centers(1,1)-half_ROI_size+1 : ROI_centers(1,1)+half_ROI_size-1] ;...
%     [ROI_centers(2,1)-half_ROI_size+1 : ROI_centers(2,1)+half_ROI_size-1]],1,Num_droplets);
% 
% %% Code
% 
% if (show_tracked_trajectory)
% 
%     % Initialization
%     drop_position = zeros(Num_droplets*last_image,2);
%     image_ind = zeros(last_image);
%     
%     tic;
%     
%     i=first_image;
% %     while (i+Image_dist < last_image) 
%     while (i+Image_dist < 20) 
%         i
%         droplet_found = 0;
%         % Refresh the reference image for next image to analyze
% %         Iref = imread([text, sprintf('%06d',i), image_format]);
%         Iref = imread([text, names(i).name]);
%         
%         while (~droplet_found)
%             
%             if (i+Image_dist>number_of_images)
%                 break;
%             end
%             % Load the new image
%             I2 = imread([text, names(i+Image_dist).name]);
%             I3 = I2-Iref;    % Compute the difference
%             figure(12); hold on;
%             imshow(I3);
%             if (image_color)
%                 I3 = rgb2gray(I3);  % Convert into grayscale
%             end
%             
%             % Finds droplets' new positions
% %             zone = I3([ROI(1,:)],[ROI(2,:)]);
% %             maxima = max(max(zone));
%             mask = ones( length(I3(:,1)), length(I3(1,:)) );
%             mask(5.*I1<80)=0;
%             [val y] = max(I1);
%             [val2 x] = max(val);
%             y = y(x);
%             
%             % Verifies that the maximum value corresponds to a droplet,
%             % otherwise, tries with the next frame.
%             if (maxima(1)>=Intensity_Threshold)
%                 Image_dist
%                 %             disp('yes');
%                 [y,x] = find(I3([ROI(1,:)],[ROI(2,:)])==maxima(1));
%                 y = y(1) + shift_y;
%                 x = x(1) + shift_x;
%                 drop_position(i+Image_dist,:) = [x,y];
%                 image_ind(i+Image_dist) = i+Image_dist;
%                 %             Image_dist
%                 droplet_found = 1;
%                 i = i+Image_dist;   % First image where all droplets are found becomes new reference image
%                 Image_dist = subsample;
%                 ROI_size = round(Image_dist*displ_max/pixel2mm); % Width and height of ROI depending on time distance between the two images
%                 half_ROI_size = round(ROI_size/2);
%                 % Updates ROI_centers
%                 ROI_centers = [y; x];
%                 shift_y = ROI_centers(1,:)-half_ROI_size;
%                 shift_x = ROI_centers(2,:)-half_ROI_size;
%                 ROI = repmat([[ROI_centers(1,1)-half_ROI_size+1 : ROI_centers(1,1)+half_ROI_size-1] ;...
%                     [ROI_centers(2,1)-half_ROI_size+1 : ROI_centers(2,1)+half_ROI_size-1]],1,Num_droplets);
%                 
%             else
%                 %             disp('try again')
%                 Image_dist = Image_dist+subsample;
%                 % Increases ROI size
%                 ROI_size = round(Image_dist*displ_max/pixel2mm); % Width and height of ROI depending on time distance between the two images
%                 half_ROI_size = round(ROI_size/2);
%                 shift_y = ROI_centers(1,:)-half_ROI_size;
%                 shift_x = ROI_centers(2,:)-half_ROI_size;
%                 ROI = repmat([[ROI_centers(1,1)-half_ROI_size+1 : ROI_centers(1,1)+half_ROI_size-1] ;...
%                     [ROI_centers(2,1)-half_ROI_size+1 : ROI_centers(2,1)+half_ROI_size-1]],1,Num_droplets);               
%             end
%             
%             % Cut the ROI if goes out the image
%             vec1 = [find(ROI(1,:)<=0)];
%             ROI(:,vec1) =[];
%             vec2 = [find(ROI(2,:)<=0)];
%             ROI(:,vec2)=[];
%             shift_x=shift_x+length(vec1)+length(vec2);  % Ajout un shift pour tenir compte des elements supprimes
%             shift_y=shift_y+length(vec1)+length(vec2);
%         end
%         
%     end
%     
%     toc;
%     
%     drop_position(find(drop_position(:,1)==0), :) = [];
%     image_ind(find(image_ind==0)) = [];
% end
% 
% %% Visualization
% 
% if (show_true_trajectory)
%     tic;  
%     % Shows the mean image between the first 'nb_for_mean' images
%     space_for_mean = 20;
%     nb_for_mean = (last_image/(space_for_mean+1));
%     
%     I1 = double(imread([text, names(first_image).name]));
%     Ipre = double(imread([text, names(first_image).name]));
%     for n = 1:nb_for_mean
%         Inew = double(imread([text, names(first_image+space_for_mean*n).name]));
%         I_diff = Inew-Ipre;
%         Ipre = Inew;
%         % Mettre un facteur qui dépende de l'intensité moyenne!!
%         %         weights = I_diff./(100/nb_for_mean);
%         weights = I_diff./(1000/nb_for_mean);
%         I1 = I1 + Inew.*(1 + weights);
%     end
%     I1 = uint8(round(I1./nb_for_mean));
%     toc; 
% else
%     I1 = uint8(imread([text, names(first_image).name]));
% end

if (show_tracked_trajectory)
    figure(1); hold on;
    imshow(I2); hold on;
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
        I3 = imread([text, names(image_ind(index_test)).name]);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dilation = dilate(image, kernel)
    dilation = imcomplement(imfilter(imcomplement(image),kernel,'conv'));
    dilation(dilation>0) = 1;
    dilation(dilation<0) = 0;
end

function erosion = erode(image, kernel)
    erosion = imfilter(image,kernel,'conv');
    erosion(erosion>0) = 1;
    erosion(erosion<0) = 0;
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

