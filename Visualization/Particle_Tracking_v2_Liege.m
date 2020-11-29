% Particle tracking
% Author : Olivier Leblanc
% Date : 19-11-2020

% This routine is dedicated to the computing of walking droplet's position
% from images captured on the experimental setup in Microfluidics Lab at ULiege.


% Explains the principle.

close all;
clc;

%% Plot options

show_each_ROI = 1;
show_each_track = 1;
show_each_mask = 1;

show_true_trajectory = 0;
show_tracked_trajectory = 1;
show_verify_position = 0;
show_final_ROI = 0;

%% Options
define_ROI = 1;
im_find_circles = 0;

%% Parameters
Num_droplets = 1; % Number of droplets to track
halfwidth = 90; % For the ROI size

% Informations varying with the analyzed video sequence :
folder_name = 'C:\Users\leblanco.OASIS\Documents\IngeCivilMaster2\Memoire\Liege_Data\J5_17_11\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directory = 'G_4_40';
image_format = '.jpg';
image_color = 0;
first_image = 1;
last_image = 1e5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

width = 2*halfwidth;
text = [folder_name, directory, '\'];
names = dir([folder_name, directory,'\*', image_format]);
number_of_images = numel(dir([text,'*', image_format]))-1;  % Numbers of images in the folder or to treat
if (last_image > number_of_images)
    last_image = number_of_images;
end
time_spent = zeros(1,number_of_images);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (define_ROI)
    % Define centers of initial regions of interest (ROI)
    figure(10); hold on;
    I1 = imread([text, names(first_image).name]);
    I2 = imread([text, names(first_image+1).name]);
    imshow(I1);
    axis on;
    title('\color{red} Select the ROI of the stadium', 'FontSize', 20);
    rect = round(getrect());
    x_ul = rect(1);
    y_ul = rect(2);
    wx = rect(3);
    hy = rect(4);
    N = wx*hy;

    xmin = x_ul;
    xmax = x_ul+wx;
    ymin = y_ul;
    ymax = y_ul+hy;
    I1 = I1(ymin:ymax, xmin:xmax);
    [Im_res_y, Im_res_x, ~] = size(I1); % Finds the orientation and resolution of the images.
    imshow(I1);
    axis on;
    title('\color{red} please click on the droplet to define its initial coordinates', 'FontSize', 20);
    [xi, yi] = ginput(1);
    ROI_centers_init = [round(yi); round(xi)];
    close(figure(10));
end

ROIx = max(1,ROI_centers_init(2)-halfwidth) : min(Im_res_x,ROI_centers_init(2)+halfwidth);
ROIy = max(1,ROI_centers_init(1)-halfwidth) : min(Im_res_y,ROI_centers_init(1)+halfwidth);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (show_tracked_trajectory)
    % Use morphological operators to remove salt noise
    square_22 = [1 1; 1 1];
    square_33 = [1 1 1; 1 1 1; 1 1 1];
    square_44 = [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1];
    rect_42 = [1 1; 1 1; 1 1; 1 1];
    line3 = [1 1 1];

    kernel = square_22;

    matx1 = repmat(1:Im_res_x,Im_res_y,1);
    maty1 = repmat(transpose(1:Im_res_y),1,Im_res_x);


    if (show_each_track + show_each_ROI + show_each_mask)
        fig16 = figure(16); hold on;
        set(fig16, 'position',[0 50 1000 600]);
        ax0 = gca;
        set(ax0,'visible','off');
        % Top-left
        axes1 = axes('Parent',fig16,...
            'Position',[0.05 0.55 0.35 0.35]);
        hold(axes1,'on');
        set(axes1,'TickLabelInterpreter','Latex','FontSize',16); 
        box(axes1,'on');
        axis(axes1,'ij');
        xlim(axes1,[0 wx]);
        ylim(axes1,[0 hy]);
        set(axes1,'TickLabelInterpreter','Latex','FontSize',16); 
        % Bottom-left
        axes2 = axes('Parent',fig16,...
            'Position',[0.05 0.1 0.35 0.35]);
        hold(axes2,'on');
        set(axes2,'TickLabelInterpreter','Latex','FontSize',16); 
        axis(axes2,'ij');
        set(axes2,'TickLabelInterpreter','Latex','FontSize',16); 
        xlim(axes2,[0 2*halfwidth]);
        ylim(axes2,[0 2*halfwidth]);
        % Right
        axes3 = axes('Parent',fig16,'Position',[0.23 0.1 0.8 0.8]);
        hold(axes3,'on');
        axis(axes3,'ij');
        set(axes3,'TickLabelInterpreter','Latex','FontSize',16); 
        xlim(axes3,[0 wx]);
        ylim(axes3,[0 hy]);
        % Legend 
        p1=plot(xm,ym,'^m','Parent',axes3, 'MarkerSize',5,'MarkerFaceColor','m','MarkerEdgeColor','k','LineWidth',0.6);
        p2=plot(indx(1)+xb,indy(1)+yb,'^','Parent',axes3, 'MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',0.6);
        p3=plot([indx(1),indx(1)], [indy(1),indy(end)],'Parent',axes3,'color','m', 'LineWidth', 2.0);   
        p4=plot([ROIx(1), ROIx(1)], [ROIy(1), ROIy(end)],'Parent',axes3,'color','b', 'LineWidth', 2.0);  
        p5=plot([ROIx(1)-width, ROIx(1)-width], [ROIy(1)-width, ROIy(end)+width],'Parent',axes1,'color','g', 'LineWidth', 2.0);   
        [h,icons]=legend(axes3, [p1,p2,p4,p5,p3], {'$\mathbf{x}_{m}$','$\mathbf{x}_{pos}$','$\mathcal{R}_{1}$','$\mathcal{R}_{1,ext}$','$\mathcal{R}_{2}$'}, 'Interpreter','Latex','FontSize', 18.0);
        icons = findobj(icons,'Type','line');
        set(icons,'MarkerSize',12);
        set(h,'Position',[0.74 0.665 0.116 0.217],...
            'Interpreter','latex','AutoUpdate','off');
    end
    % legend('Brightest pixel in the image difference','Gravity center of the mask','Most correlated point','Darkest point inside the ROI');

%     imref = imread([text, names(21).name]);
%     imref = imref(180:380, 1010:1220);
%     H = length(imref(:,1));
%     W = length(imref(1,:));
%     [X,Y] = meshgrid(1:H, 1:W);
%     % figure(13); hold on; imshow(imref); axis on;
%     imref2 = double(imref);
%     mref = mean(mean(imref2));
%     imref_window = imref2-mref;
%     fourier_imref_corr = fftshift(fft2(imref_window));
%     sigma = 0.001;

    drop_position = zeros(Num_droplets*last_image,2);

    I1 = imread([text, names(first_image).name]);
    I1 = I1(ymin:ymax,xmin:xmax);
    % ------------- Starting loop --------------
    for i = first_image : first_image+5 %last_image
        i
        tic;
        
        % Read images
        I2 = imread([text, names(i+1).name]);
        I2 = I2(ymin:ymax,xmin:xmax);
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

        % Restrict the mask around the previous location
        mask([1:max(1,ROIy(1)-width), min(xmax, ROIy(end)+width):end], :) = 0;
        mask(:, [1:max(1, ROIx(1)-width), min(ymax, ROIx(end)+width):end]) = 0;
        matx = matx1.*mask;
        maty = maty1.*mask;
        
        xm = round(sum(sum(matx))/sum(sum(matx~=0)));
        ym = round(sum(sum(maty))/sum(sum(maty~=0)));
        indy = max(1,ym-halfwidth):min(Im_res_y,ym+halfwidth);
        indx = max(1,xm-halfwidth):min(Im_res_x,xm+halfwidth);

%     %     % Compute the correlation with imref
%         I2b = double(I2(indy,indx));
%         m2 = mean(mean(I2b));
%         I2_window = I2b-m2;
%         fourier_I2_corr = fftshift(fft2(I2_window,H,W));
%         fourier_prod = fourier_I2_corr.*conj(fourier_imref_corr);
%         corr = fftshift(real(ifft2(ifftshift(fourier_prod))));
%         corr = rot90(corr);
%         [val xcorr] = max(abs(corr));
%         [val2 ycorr] = max(val);
%         xcorr = xcorr(ycorr);

        % Compute the droplet's position, either with imfindcircles or darkest
        [val yb] = min(I2(indy,indx));
        [val2 xb] = min(val);
        yb = yb(xb);        
        if (im_find_circles)
            [centers, radii, metric] = imfindcircles(round(I2(indy,indx)./4), [2, 5], 'Method', 'TwoStage', 'ObjectPolarity','dark');
            if (length(radii)>0)
                xb = round(centers(1,1));
                yb = round(centers(1,2));
            end
        end

        if (show_each_mask)
            imshow(mask, 'Parent', axes1);
            plot(xm,ym, '^b','Parent',axes1, 'MarkerSize',8,'MarkerFaceColor','m','MarkerEdgeColor','k','LineWidth',1);
            axis(axes1, 'on'); 
            plot([ROIx(1)-width, ROIx(1)-width], [ROIy(1)-width, ROIy(end)+width],'Parent',axes1,'color','g', 'LineWidth', 2.0);   % Left-side
            plot([ROIx(end)+width, ROIx(end)+width], [ROIy(1)-width, ROIy(end)+width],'Parent',axes1,'color','g', 'LineWidth', 2.0);  % right-side
            plot([ROIx(1)-width, ROIx(end)+width], [ROIy(1)-width, ROIy(1)-width],'Parent',axes1,'color','g', 'LineWidth', 2.0);  % upper-side
            plot([ROIx(1)-width, ROIx(end)+width], [ROIy(end)+width, ROIy(end)+width],'Parent',axes1,'color','g', 'LineWidth', 2.0);  % lower-side
        end

        if (show_each_track)
            imshow(I2, 'Parent', axes3);
            p1=plot(xm,ym,'^m','Parent',axes3, 'MarkerSize',5,'MarkerFaceColor','m','MarkerEdgeColor','k','LineWidth',0.6);
%             scatter(x+ROIx(1),y+ROIy(1),'Parent',axes3,'b');
            plot(indx(1)+xb,indy(1)+yb,'^','Parent',axes3, 'MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',0.6);
            plot([indx(1),indx(1)], [indy(1),indy(end)],'Parent',axes3,'color','m', 'LineWidth', 2.0);   % Left-side
            plot([indx(end), indx(end)], [indy(1), indy(end)],'Parent',axes3,'color','m', 'LineWidth', 2.0);  % right-side
            plot([indx(1), indx(end)], [indy(1), indy(1)],'Parent',axes3,'color','m', 'LineWidth', 2.0);  % upper-side
            plot([indx(1), indx(end)], [indy(end), indy(end)],'Parent',axes3,'color','m', 'LineWidth', 2.0);  % lower-side
            plot([ROIx(1), ROIx(1)], [ROIy(1), ROIy(end)],'Parent',axes3,'color','b', 'LineWidth', 2.0);   % Left-side
            plot([ROIx(end), ROIx(end)], [ROIy(1), ROIy(end)],'Parent',axes3,'color','b', 'LineWidth', 2.0);  % right-side
            plot([ROIx(1), ROIx(end)], [ROIy(1), ROIy(1)],'Parent',axes3,'color','b', 'LineWidth', 2.0);  % upper-side
            plot([ROIx(1), ROIx(end)], [ROIy(end), ROIy(end)],'Parent',axes3,'color','b', 'LineWidth', 2.0);  % lower-side
            axis(axes3, 'on');
            %             title(ax0, ['N°',num2str(i)]);
          end  

        if (show_each_ROI)
              imshow(I2(indy, indx), 'Parent', axes2);
              plot([1 1], [1 width],'Parent',axes2,'color','m', 'LineWidth', 2.0);   % Left-side
              plot([width width], [1, width],'Parent',axes2,'color','m', 'LineWidth', 2.0);  % right-side
              plot([1 width], [1 1],'Parent',axes2,'color','m', 'LineWidth', 2.0);  % upper-side
              plot([1 width], [width width],'Parent',axes2,'color','m', 'LineWidth', 2.0);  % lower-side
            if (im_find_circles)
                if (length(radii)>=2)
                    viscircles(centers(1,:), radii(1),'EdgeColor','g');
                    viscircles(centers(2,:), radii(2),'EdgeColor','m');
                end
            end
            plot(xm-indx(1),ym-indy(1), '^b','Parent',axes2, 'MarkerSize',10,'MarkerFaceColor','m','MarkerEdgeColor','k','LineWidth',1.2);
            plot(xb,yb,'^r','Parent',axes2, 'MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',1.2);
%             plot(x-indx(1)+ROIx(1),y-indy(1)+ROIy(1),'^m','Parent',axes2, 'MarkerSize',10,'MarkerFaceColor','m','MarkerEdgeColor','k');
%             scatter(xcorr, ycorr,'g');
            axis(axes2, 'on');
        end

        drawnow;

        drop_position(i,:) = [xb+indx(1),yb+indy(1)];

        ROIx = max(1,xb+indx(1)-halfwidth) : min(Im_res_x,xb+indx(1)+halfwidth);
        ROIy = max(1,yb+indy(1)-halfwidth) : min(Im_res_y,yb+indy(1)+halfwidth);
        
        I1 = I2;
        
        time_spent(i) = toc;
    end

    drop_position(find(drop_position(:,1)==0), :) = [];
    drop_position(:,1) = drop_position(:,1) + xmin;
    drop_position(:,2) = drop_position(:,2) + ymin;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (show_tracked_trajectory)
    figure(1); hold on;
    imshow(I2); hold on;
    axis on;
    set(gca,'Ydir','reverse');
    scatter(drop_position(1,1)-xmin,drop_position(1,2)-ymin, 'b');
    n = length(drop_position(:,1));
    % put "+ 50" only to shift tracked trajectory and compare with mean of images
    p = plot(drop_position(:,1)-xmin  ,drop_position(:,2)-ymin , 'r', 'LineWidth', 1.0);
    % modified jet-colormap
    cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
    drawnow;
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
%     legend('Initial position', 'Trajectory');
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

if (show_true_trajectory)
    tic;  
    % Shows the mean image between the first 'nb_for_mean' images
    space_for_mean = 20;
    nb_for_mean = (last_image/(space_for_mean+1));
    
    I1 = double(imread([text, names(first_image).name]));
    Ipre = double(imread([text, names(first_image).name]));
    for n = 1:nb_for_mean
        Inew = double(imread([text, names(first_image+space_for_mean*n).name]));
        I_diff = Inew-Ipre;
        Ipre = Inew;
        % Mettre un facteur qui dépende de l'intensité moyenne!!
        %         weights = I_diff./(100/nb_for_mean);
        weights = I_diff./(100/nb_for_mean);
        I1 = I1 + Inew.*(1 + weights);
    end
    I1 = uint8(round(I1./nb_for_mean));
    toc; 
    figure(11); hold on;
    imshow(I1); hold on;
    axis on;
else
    I1 = uint8(imread([text, names(first_image).name]));
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

