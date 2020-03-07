%% Script description
% Correlation based surface visualization
% Author : Olivier Leblanc
% Date : 19-11-2019

clear all;
close all;
clc;

% Feature1 : Proximity in gray_level
% Feature2 : Proximity in image's position
% Feature3 : Linear correlation


%% Good views

%%%%%%% 13-11-19 - droplet

% % Files parameters
% dir = '13-11-19_Droplet';   % Directory name
% filename = 'droplet';
% ref_num = 0;   % Number of the reference image
% cur_num = 800;
% Resolution = 185e-6;  % Real size of a pixel
%
% % Processing parameters
% num_tests = 5;
% neighbour_size = 10;     % For correlation
% x_size = 300;            % For ROI
% y_size = 300;
% ROI_center = [480, 540];


%%%%%%% 13-11-19 - airgap

% dir = '13-11-19_AirGap';   % Directory name
% filename = 'airgap';
% ref_num = 0;   % Number of the reference image
% cur_num = 316;
% Resolution = 185e-6;  % Real size of a pixel
%
% % Processing parameters
% num_tests = 5;
% neighbour_size = 10;     % For correlation
% x_size = 300;            % For ROI
% y_size = 300;
% ROI_center = [695, 790];


%%%%%%% 13-11-19 - faraday

% % Files parameters
% dir = '13-11-19_Faraday';   % Directory name
% filename = 'faraday';
% ref_num = 0;   % Number of the reference image
% cur_num = 476;
% Resolution = 185e-6;  % Real size of a pixel
%
% % Processing parameters
% num_tests = 5;
% neighbour_size = 10;     % For correlation
% x_size = 300;            % For ROI
% y_size = 300;
% ROI_center = [520, 1200];


%%%%%%% 25-11-19 - zoom

% % Files parameters
% dir = '25-11-19_Zoom';   % Directory name
% filename = 'zoom';
% ref_num = 0;   % Number of the reference image
% cur_num = 1022;
% Resolution = 93e-6;  % Real size of a pixel
% 
% % Processing parameters
% num_tests = 5;
% neighbour_size = 10;     % For correlation
% x_size = 300;            % For ROI
% y_size = 300;
% ROI_center = [960-100, 540+100];

%% Options

lighting_effects = 0;  % Enhance surface visualization with ligthning effect
true_height = 0;    % 1 : gives true reconstructed height. 0 : gives height difference from equilibrium
color_processing = 1;   % Compute distance between the positions of the point with RGB
testweight = 0;
use_correl = 0; % increases the computational time a lot
fourier_filtering = 1;  % Spatial filtering in Fourier domain.
hard_fourier_thresholding = 1;  % Needs 'fourier_filtering' activated.

%% Parameters

dir = '13-11-19_AirGap';   % Directory name
filename = 'airgap';
ref_num = 0;   % Number of the reference image
cur_num = 316;
Resolution = 185e-6;  % Real size of a pixel

% Processing parameters
num_tests = 5;
neighbour_size = 10;     % For correlation
x_size = 300;            % For ROI
y_size = 300;
ROI_center = [695, 790];
percentile_to_cut = 0.99;   % Percentile under which the fourier values are put to 0

% Physical system values
lambda = 4.75e-3;   % Typical wavelength
h0 = 5e-3; % Water level
n = 1.4; % Silicon oil refractive index
alpha = 1-1/n;  % See paper
hp = 6e-3;  % Liquid depth
H = 1;  % Distance pattern-camera
h_star = 2/(1/(alpha*hp)-1/H);


%% Test of different weights
if (testweight)
    
    x_size = 100;
    y_size = 100;
    
    % Meshgrid creation
    [X,Y] = meshgrid(1:x_size,1:y_size);
    
    % Spatial filtering in Fourier domain.
    half_ring_width = 3;
    radius_fourier = round(lambda/Resolution/1080*x_size);
    image_center = [x_size/2+1, y_size/2+1];
    
    if (fourier_filtering)
        circlePixels = (Y - image_center(2)).^2 ...
            + (X - image_center(1)).^2 <= (radius_fourier-half_ring_width).^2;
        circlePixels2 = (Y - image_center(2)).^2 ...
            + (X - image_center(1)).^2 >= (2*radius_fourier).^2;
        circlePixels = circlePixels + circlePixels2;
    end
    
    dx = zeros(y_size, x_size);
    dy = zeros(y_size, x_size);
    
    % Load reference image
    Iref = imread(['C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Visualization\Images\',dir,'\',filename,'_00',sprintf('%04d',ref_num),'.jpeg']);
    Iref_gray = rgb2gray(Iref);
    
    if (color_processing)
        Iref_ROI = Iref(ROI_center(1)-y_size/2+1:ROI_center(1)+y_size/2 , ROI_center(2)-x_size/2+1 : ROI_center(2)+x_size/2, :);
    else
        Iref_ROI = Iref(ROI_center(1)-y_size/2+1:ROI_center(1)+y_size/2 , ROI_center(2)-x_size/2+1 : ROI_center(2)+x_size/2);
    end
    
    % Load current image
    I2 = imread(['C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Visualization\Images\',dir,'\',filename,'_00',sprintf('%04d',cur_num),'.jpeg']);
    I2_gray = rgb2gray(I2);
    
    if (color_processing)
        I2_ROI = I2(ROI_center(1)-y_size/2+1:ROI_center(1)+y_size/2 , ROI_center(2)-x_size/2+1 : ROI_center(2)+x_size/2, :);
    else
        I2_ROI = I2(ROI_center(1)-y_size/2+1:ROI_center(1)+y_size/2 , ROI_center(2)-x_size/2+1 : ROI_center(2)+x_size/2);
    end
    
    for test = 1:num_tests
        w1 = 10*test/num_tests^2;   % Weight of proximity in gray level
        w2 = 1-10*test/num_tests^2;   % Weight of proximity in image's position
        w3 = 0.33;  % Weight for correlation
        
        tic;
        
        for i = 1:y_size    % y-axis
            for j = 1:x_size    % x-axis
                
                % Squeeze cuts the useless 3rd dimension of the result
                % Here, ft1 takes each difference in color separately
                if (color_processing)
                    ft1 = double(squeeze(mean(permute(Iref_ROI(i,j,:)-I2_ROI, [3, 1, 2]))));
                else
                    ft1 = double(Iref_ROI(i,j)-I2_ROI);
                end
                
                ft2 = repmat((i-(1:y_size))',1,x_size).^2+repmat((j-(1:x_size)),y_size,1).^2;
                % Normalize features
                mean_ft1 = mean(mean(ft1));
                mean_ft2 = mean(mean(ft2));
                ft1 = (ft1-mean_ft1)./(sqrt(sum(sum((ft1-mean_ft1).^2))./(y_size*x_size)));
                ft2 = (ft2-mean_ft2)./(sqrt(sum(sum((ft2-mean_ft2).^2))./(y_size*x_size)));
                
                cost = w1.*ft1+w2.*ft2;
                
                if (use_correl)
                    
                    fourier_Iref_corr = fftshift(fft2(...
                        Iref_gray(ROI_center(1)+i-x_size/2:ROI_center(1)+i+x_size/2-1 ,...
                        ROI_center(2)+j-y_size/2:ROI_center(2)+j+y_size/2-1)));
                    
                    fourier_I2_corr = fftshift(fft2(...
                        I2_gray(ROI_center(1)+j-x_size/2:ROI_center(1)+j+x_size/2-1 ,...
                        ROI_center(2)+i-y_size/2:ROI_center(2)+i+y_size/2-1)));
                    % Pas oublier que la corrélation c'est pas exactement une
                    % convolution donc dans Fourier il faut le complexe conj d'un des deux
                    fourier_prod = fourier_Iref_corr*conj(fourier_I2_corr);
                    fourier_prod(interrogation_window_size-1, interrogation_window_size-1) = 0;
                    corr_mat = real(ifft2(fftshift(fourier_prod)));
                    corr_mat = corr_mat./max(max(corr_mat));
                    %                 ft3 = correlation(Iref, I2, neighbour_size, i, j, ROI_center, x_size, y_size);
                    ft3 = 1./corr_mat;
                    cost = w1.*ft1+w2.*ft2+w3.*ft3;
                end
                
                [val y] = min(cost);
                [val2 x] = min(val);
                y = y(x);
                dx(j, i) = j-x;
                dy(j, i) = i-y;
            end
        end
        
        toc;    % Displays time needed to computed those two 4D matrices.
        
        % Spatial filtering in Fourier domain.
        inner_radius = 3;
        radius_fourier = round(lambda/Resolution/1080*x_size);
        image_center = [x_size/2+1, y_size/2+1];
        
        fourier_dx = fftshift(fft2(dx));
        fourier_dy = fftshift(fft2(dy));
        
        if (fourier_filtering)
            circlePixels = (Y - image_center(2)).^2 ...
                + (X - image_center(1)).^2 <= (inner_radius).^2;
            circlePixels2 = (Y - image_center(2)).^2 ...
                + (X - image_center(1)).^2 >= (2*radius_fourier).^2;
            circlePixels = circlePixels + circlePixels2 ;
            filtered_fourier_dx = fourier_dx;
            filtered_fourier_dy = fourier_dy;
            filtered_fourier_dx(find(circlePixels==1)) = 0;
            filtered_fourier_dy(find(circlePixels==1)) = 0;
            if (hard_fourier_thresholding)
                val_dx = sort(reshape(real(filtered_fourier_dx), [x_size*y_size, 1]));
                val_dy = sort(reshape(real(filtered_fourier_dy), [x_size*y_size, 1]));
                hard_fourier_threshold_dx = val_dx(round(percentile_to_cut*x_size*y_size));
                hard_fourier_threshold_dy = val_dy(round(percentile_to_cut*x_size*y_size));
                
                % Hard thresholding
                filtered_fourier_dx(filtered_fourier_dx<hard_fourier_threshold_dx)=0;
                filtered_fourier_dy(filtered_fourier_dy<hard_fourier_threshold_dy)=0;
            end
            
            dx = real(ifft2(ifftshift(filtered_fourier_dx)));
            dy = real(ifft2(ifftshift(filtered_fourier_dy)));
        else
            dx = real(ifft2(ifftshift(fourier_dx)));
            dy = real(ifft2(ifftshift(fourier_dy)));
        end
        
        % Rescale displacement in true distance [m]
        dx = dx.*Resolution;
        dy = dy.*Resolution;
        
        % Correct for unmatched corners
        depl(1:round(x_size/20),1:round(y_size/20))=0;
        depl(1:round(x_size/20),y_size-round(y_size/20):y_size)=0;
        depl(x_size-round(x_size/20):x_size, 1:round(y_size/20))=0;
        depl(x_size-round(x_size/20):x_size, y_size-round(y_size/20):y_size)=0;
        
        dx(1:round(x_size/20),1:round(y_size/20))=0;
        dx(1:round(x_size/20),y_size-round(y_size/20):y_size)=0;
        dx(x_size-round(x_size/20):x_size, 1:round(y_size/20))=0;
        dx(x_size-round(x_size/20):x_size, y_size-round(y_size/20):y_size)=0;
        
        dy(1:round(x_size/20),1:round(y_size/20))=0;
        dy(1:round(x_size/20),y_size-round(y_size/20):y_size)=0;
        dy(x_size-round(x_size/20):x_size, 1:round(y_size/20))=0;
        dy(x_size-round(x_size/20):x_size, y_size-round(y_size/20):y_size)=0;
        
        % Inverse gradient
        dhdx = -dx./h_star;
        dhdy = -dy./h_star;
        if (true_height)
            h = intgrad2(dhdy,dhdx,Resolution,Resolution,h0);
        else
            h = intgrad2(dhdy,dhdx,Resolution,Resolution,0);
        end
        
        figure(50+test); hold on;
        s = surf(X, Y, h);
        if (lighting_effects)
            s.EdgeColor = 'none';
            light               % create a light
            lighting gouraud    % preferred method for lighting curved surfaces
        end
        direction = [0 0 1];
        rotate(s, direction,-90);
        colorbar;
        
    end
    
    f1 = figure(1); hold on;
    imshow(Iref_ROI);
    axis on;
    title('Reference image');
    
    f2 = figure(2); hold on;
    imshow(I2_ROI);
    axis on;
    title('Modified image');
    
    
else
    w1 = 0.8;   % Weight of proximity in gray level
    w2 = 0.2;   % Weight of proximity in image's position
    w3 = 0.05;  % Weight for correlation
    
    % Load reference image
    Iref = imread(['C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Visualization\Images\',dir,'\',filename,'_00',sprintf('%04d',ref_num),'.jpeg']);
    Iref_gray = rgb2gray(Iref);
    
    if (color_processing)
        Iref_ROI = Iref(ROI_center(1)-y_size/2+1:ROI_center(1)+y_size/2 , ROI_center(2)-x_size/2+1 : ROI_center(2)+x_size/2, :);
    else
        Iref_ROI = Iref(ROI_center(1)-y_size/2+1:ROI_center(1)+y_size/2 , ROI_center(2)-x_size/2+1 : ROI_center(2)+x_size/2);
    end
    
    % Load current image
    I2 = imread(['C:\Users\Leblanc\Documents\IngeCivilMaster2\Memoire\Visualization\Images\',dir,'\',filename,'_00',sprintf('%04d',cur_num),'.jpeg']);
    I2_gray = rgb2gray(I2);
    
    if (color_processing)
        I2_ROI = I2(ROI_center(1)-y_size/2+1:ROI_center(1)+y_size/2 , ROI_center(2)-x_size/2+1 : ROI_center(2)+x_size/2, :);
    else
        I2_ROI = I2(ROI_center(1)-y_size/2+1:ROI_center(1)+y_size/2 , ROI_center(2)-x_size/2+1 : ROI_center(2)+x_size/2);
    end
    
    tic;
    
    depl = zeros(x_size, y_size);
    dx = zeros(x_size, y_size);
    dy = zeros(x_size, y_size);
    
    % For each pixel of the reference image
    for i = 1:y_size    % y-axis
        for j = 1:x_size    % x-axis
            
            % Squeeze cuts the useless 3rd dimension of the result
            % Here, ft1 takes each difference in color separately
            if (color_processing)
                ft1 = double(squeeze(mean(permute(Iref_ROI(i,j,:)-I2_ROI, [3, 1, 2]))));
            else
                ft1 = double(Iref_ROI(i,j)-I2_ROI);
            end
            ft2 = repmat((i-(1:y_size))',1,x_size).^2+repmat((j-(1:x_size)),y_size,1).^2;
            
            % Normalize features
            mean_ft1 = mean(mean(ft1));
            mean_ft2 = mean(mean(ft2));
            ft1 = (ft1-mean_ft1)./(sqrt(sum(sum((ft1-mean_ft1).^2))./(y_size*x_size)));
            ft2 = (ft2-mean_ft2)./(sqrt(sum(sum((ft2-mean_ft2).^2))./(y_size*x_size)));
            
            cost = w1.*ft1+w2.*ft2;
            
            % Apply linear cross-correlation in Fourier domain.
            if (use_correl)
                
                fourier_Iref_corr = fftshift(fft2(...
                    Iref_gray(ROI_center(1)+i-x_size/2:ROI_center(1)+i+x_size/2-1 ,...
                    ROI_center(2)+j-y_size/2:ROI_center(2)+j+y_size/2-1)));
                
                fourier_I2_corr = fftshift(fft2(...
                    I2_gray(ROI_center(1)+j-x_size/2:ROI_center(1)+j+x_size/2-1 ,...
                    ROI_center(2)+i-y_size/2:ROI_center(2)+i+y_size/2-1)));
                % Pas oublier que la corrélation c'est pas exactement une
                % convolution donc dans Fourier il faut le complexe conj d'un des deux
                fourier_prod = fourier_Iref_corr*conj(fourier_I2_corr);
                fourier_prod(x_size/2-1, y_size/2-1) = 0;
                ft3 = real(ifft2(fftshift(fourier_prod)));
                %                 ft3 = correlation(Iref, I2, neighbour_size, i, j, ROI_center, x_size, y_size);
                % Normalization
                mean_ft3 = mean(mean(ft3));
                ft3 = (ft3-mean_ft3)./(sqrt(sum(sum((ft3-mean_ft3).^2))./(y_size*x_size)));
                ft3 = 1-ft3;
                cost = w1.*ft1+w2.*ft2+w3.*ft3;
            end
            
            [val y] = min(cost);
            [val2 x] = min(val);
            y = y(x);
            dx(j, i) = j-x;
            dy(j, i) = i-y;
            depl(j, i) = sqrt((i-y)^2+(j-x)^2);
        end
    end
    
    toc;
    
    % Meshgrid creation
    [X,Y] = meshgrid(1:x_size,1:y_size);
    
    % Spatial filtering in Fourier domain.
    inner_radius = 3;
    radius_fourier = round(lambda/Resolution/1080*x_size);
    image_center = [x_size/2+1, y_size/2+1];
    
    fourier = fftshift(fft2(depl));
    fourier_dx = fftshift(fft2(dx));
    fourier_dy = fftshift(fft2(dy));
    
    if (fourier_filtering)
        circlePixels = (Y - image_center(2)).^2 ...
            + (X - image_center(1)).^2 <= (inner_radius).^2;
        circlePixels2 = (Y - image_center(2)).^2 ...
            + (X - image_center(1)).^2 >= (2*radius_fourier).^2;
        circlePixels = circlePixels + circlePixels2 ;
        filtered_fourier = fourier;
        filtered_fourier(find(circlePixels==1)) = 0;
        filtered_fourier_dx = fourier_dx;
        filtered_fourier_dy = fourier_dy;
        filtered_fourier_dx(find(circlePixels==1)) = 0;
        filtered_fourier_dy(find(circlePixels==1)) = 0;
        
        if (hard_fourier_thresholding)
            val_dx = sort(reshape(real(filtered_fourier_dx), [x_size*y_size, 1]));
            val_dy = sort(reshape(real(filtered_fourier_dy), [x_size*y_size, 1]));
            hard_fourier_threshold_dx = val_dx(round(percentile_to_cut*x_size*y_size));
            hard_fourier_threshold_dy = val_dy(round(percentile_to_cut*x_size*y_size));
            
            % Hard thresholding
            filtered_fourier_dx(filtered_fourier_dx<hard_fourier_threshold_dx)=0;
            filtered_fourier_dy(filtered_fourier_dy<hard_fourier_threshold_dy)=0;
            filtered_fourier(filtered_fourier<sqrt(hard_fourier_threshold_dx^2+hard_fourier_threshold_dy^2))=0;
        end
        
        depl = real(ifft2(ifftshift(filtered_fourier)));
        dx = real(ifft2(ifftshift(filtered_fourier_dx)));
        dy = real(ifft2(ifftshift(filtered_fourier_dy)));
    else
        depl = real(ifft2(ifftshift(fourier)));
        dx = real(ifft2(ifftshift(fourier_dx)));
        dy = real(ifft2(ifftshift(fourier_dy)));
    end
    
    
    % Rescale displacement in true distance [m]
    depl = depl.*Resolution;
    dx = dx.*Resolution;
    dy = dy.*Resolution;
    
    % Correct for unmatched corners
    depl(1:round(x_size/20),1:round(y_size/20))=0;
    depl(1:round(x_size/20),y_size-round(y_size/20):y_size)=0;
    depl(x_size-round(x_size/20):x_size, 1:round(y_size/20))=0;
    depl(x_size-round(x_size/20):x_size, y_size-round(y_size/20):y_size)=0;
    
    dx(1:round(x_size/20),1:round(y_size/20))=0;
    dx(1:round(x_size/20),y_size-round(y_size/20):y_size)=0;
    dx(x_size-round(x_size/20):x_size, 1:round(y_size/20))=0;
    dx(x_size-round(x_size/20):x_size, y_size-round(y_size/20):y_size)=0;
    
    dy(1:round(x_size/20),1:round(y_size/20))=0;
    dy(1:round(x_size/20),y_size-round(y_size/20):y_size)=0;
    dy(x_size-round(x_size/20):x_size, 1:round(y_size/20))=0;
    dy(x_size-round(x_size/20):x_size, y_size-round(y_size/20):y_size)=0;
    
    % Inverse gradient
    dhdx = -dx./h_star;
    dhdy = -dy./h_star;
    if (true_height)
        h = intgrad2(dhdy,dhdx,Resolution,Resolution,h0);
    else
        h = intgrad2(dhdy,dhdx,Resolution,Resolution,0);
    end
    
    
    %% Plot
    
    f1 = figure(1); hold on;
    imshow(Iref_ROI);
    axis on;
    title('Reference image');
    
    f2 = figure(2); hold on;
    imshow(I2_ROI);
    axis on;
    title('Modified image');
    
    figure(5); hold on;
    s = surf(X, Y, depl);
    s.EdgeColor = 'none';
    direction = [0 0 1];
    rotate(s, direction,-90);
    colorbar;
    title('Surface profile of the displacement norm');
    
    figure(6); hold on;
    s2 = surf(X, Y, abs(fourier));
    s2.EdgeColor = 'none';
    colorbar;
    title('FFT of displacement norm after filtering')
    
    f8 = figure(8); hold on;
    s2 = surf(X, Y, abs(fourier_dx));
    s2.EdgeColor = 'none';
    colorbar;
    title('dx fft without Fourier filtering');
    
    f9 = figure(9); hold on;
    s2 = surf(X, Y, abs(fourier_dy));
    s2.EdgeColor = 'none';
    colorbar;
    title('dy fft without Fourier filtering');
    
    if (fourier_filtering)
        f10 = figure(10); hold on;
        s2 = surf(X, Y, abs(filtered_fourier_dx));
        s2.EdgeColor = 'none';
        colorbar;
        title('dx fft after Fourier filtering');
        
        f11 = figure(11); hold on;
        s2 = surf(X, Y, abs(filtered_fourier_dy));
        s2.EdgeColor = 'none';
        colorbar;
        title('dy fft after Fourier filtering');
    end
    
    f7 = figure(7); hold on;
    s3 = surf(X, Y, h);
    s3.EdgeColor = 'none';
    if (lighting_effects)
        light               % create a light
        lighting gouraud    % preferred method for lighting curved surfaces
    end
    direction = [0 0 1];
    rotate(s3, direction,-90);
    colorbar;
    title('Surface profile');
    xlabel('x');
    ylabel('y');
    
    if (use_correl)
        [yi, ya] = meshgrid(1:length(fourier_I2_corr(1,:)),1:length(fourier_I2_corr(:,1)));
        f25 = figure(25); hold on;
        s4 = surf(yi, ya, abs(fourier_prod));
        s4.EdgeColor = 'none';
        colorbar;
        title('correlation');
    end
    
    
    % Places figures
    movegui(f1,[0 700-y_size]);
    movegui(f2,[2*x_size 700-y_size]);
    movegui(f7,[1000 400]);
    movegui(f8,[0 100]);
    movegui(f9,[650 100]);
    movegui(f10,[100 150]);
    movegui(f11,[750 150]);
end



function correl_mat = correlation(I1, I2, interrogation_window_size, m, n, ROI_center, x_size, y_size)
% Compute the correlation matrix between I1(m,n) and I2
correl_mat = zeros(y_size, x_size);
for i = 1 : x_size
    for j = 1 : y_size
        I1_ROI = I1(ROI_center(1)+m-interrogation_window_size/2:ROI_center(1)+m+interrogation_window_size ,...
            ROI_center(2)+n-interrogation_window_size/2:ROI_center(2)+n+interrogation_window_size);
        I2_ROI = I2(ROI_center(1)+j-interrogation_window_size/2:ROI_center(1)+j+interrogation_window_size ,...
            ROI_center(2)+i-interrogation_window_size/2:ROI_center(2)+i+interrogation_window_size);
        I1_mean = mean(mean(I1_ROI));
        I2_mean = mean(mean(I2_ROI));
        correl_mat(j,i) = sum(sum((I1_ROI-I1_mean).*(I2_ROI-I2_mean)))...
            ./sqrt(sum(sum((I1_ROI-I1_mean).^2)).*sum(sum((I2_ROI-I2_mean).^2)));
    end
end
end




