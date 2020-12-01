function [resolution, diameter] = draw_stadium(directory, imfolder, M, is_cell, stadindex )

    % This function aims to make a clean plot of the droplet trajectory
    % inside the stadium. The stadium is cut out by whitening the zones
    % outside of the geometry. The stadium is also rotated to horizontality.
    % 
    % Inputs :
    %       - dictory : path of the day of experiment
    %       - imfolder : name of the experiment folder
    %       - M : distance to Faraday theshold   
    %       - is_cell : boolean, 1 if cell array, 0 if simple array.
    %       - stadindex : index of the stadium used (crf paper), 1,...,6
    % Outputs :
    %       - resolution : the pixel resolution of the images
    %       - diameter : the droplet diameter computed from its mean speed.
    %
    
    close all; clc;
        
    %% Options
    add_sym = 0;
    jetmap = 0;
    %% Parameters
    padL = 200;
    ref_lengths = [95, 57, 38, 28.5, 19, 14.25, 76].*1e-3;
    
    %% Reference length for the resolution
    ref_length = ref_lengths(stadindex);
    
    %% User manipulations
    
    if (~is_cell)
        matname = [directory, '\', imfolder, 'drop_positions.mat'];
        drop_position = load(matname);
        drop_position = drop_position.drop_position;
        X = drop_position(:,1);
        Y = drop_position(:,2);
    else
        matname = [directory, '\', imfolder, '_pos_cell.mat'];
        pos_cell = load(matname);
        pos_cell = pos_cell.pos_cell;
    end
    
    path = [directory, '\', imfolder, '\'];
    imnames = dir([path, '*.jpg']);
    [path, imnames(1).name]
    im = imread([path, imnames(1).name]);
    im = padarray(im, [padL, padL], 0);

    % Find stadium orientation
    figure(14); imshow(im); axis on; hold on;
    % Correct the padding
    if (~is_cell)
        X = X + padL;
        Y = Y + padL;
        plot(X, Y, 'r');
    else
        for j = 1 : length(pos_cell) 
            pos_cell{j} = pos_cell{j}+repmat([padL, padL], length(pos_cell{j}(:,1)),1);
            Xcor = pos_cell{j}(:,1);
            Ycor = pos_cell{j}(:,2);
            p = plot(Xcor, Ycor, 'r');
        end
    end
    title('\color{red} Create line with the principal axis from left to right', 'FontSize', 20);
    principal_axis = getline;
    
    % Rotate the stadium
    dx = principal_axis(2,1)-principal_axis(1,1);
    dy = principal_axis(2,2)-principal_axis(1,2);
    alpha = asin(dy/(sqrt(dx^2+dy^2)));
    imrot = imrotate( im , alpha*180/pi, 'crop' );
    imshow(imrot); axis on; hold on;
    [H, W, ~] = size(im);
    halfW = round(W/2);
    halfH = round(H/2);
    
    if (~is_cell)
        pos = [halfW; halfH] + [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[X-halfW, Y-halfH]';
        X=pos(1,:);
        Y=pos(2,:);
        plot(X, Y, 'r');
    else
        for j = 1 : length(pos_cell) 
            postemp = pos_cell{j};
            pos_cell{j} = [halfW; halfH]...
                + [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]...
                *[postemp(:,1)-halfW, postemp(:,2)-halfH]';           
            pos = pos_cell{j};
            Xcor = pos(1,:);
            Ycor = pos(2,:);
            p = plot(Xcor, Ycor, 'r');
        end
    end
    
    % Zoom on stadium
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
    imrot = imrot(ymin:ymax, xmin:xmax);
    close(figure(14));
    
    if (~is_cell)
        X=X-xmin;
        Y=Y-ymin;
    else
        for j = 1 : length(pos_cell) 
            pos_cell{j} = pos_cell{j}-[xmin;ymin];
        end
    end
    
    [H, W, ~] = size(imrot);
    N = H*W;
    halfW = round(W/2);
    halfH = round(H/2);
  
    % Select the upper and lower sidelines
    figure(14);
    imshow(imrot);  hold on;
    if (~is_cell)
        plot(X, Y, 'r');
    else
        for j = 1 : length(pos_cell) 
            Xcor = pos_cell{j}(1,:);
            Ycor = pos_cell{j}(2,:);
            p = plot(Xcor, Ycor, 'r');
        end
    end
    axis on; hold on;
    % Define the four points of the stadium  
    title('\color{red} Draw the rectangle separating the semicircles', 'FontSize', 20);
    poly = drawpolygon();
    corners4 = poly.Position;

    %% Stadium properties      
    % Correct for orientation again
    dxu = corners4(2,1)-corners4(1,1);
    dyu = corners4(2,2)-corners4(1,2);
    alphau = acos(dxu/(sqrt(dxu^2+dyu^2)));
    dxl = corners4(3,1)-corners4(4,1);
    dyl = corners4(3,2)-corners4(4,2);
    alphal = acos(dxl/(sqrt(dxl^2+dyl^2)));
    alpha2 = (alphau+alphal)/2;
    corners4=corners4';
    corners4 = [halfW; halfH] + [cos(alpha2) sin(alpha2); -sin(alpha2) cos(alpha2)]*...
        [corners4(1,:)-halfW; corners4(2,:)-halfH];
    corners4=corners4';
    imrot = imrotate( imrot , alpha2.*180/pi, 'crop' );
    imshow(imrot); axis on; hold on;
    if (~is_cell)
        pos = [halfW; halfH] + [cos(alpha2) sin(alpha2); -sin(alpha2) cos(alpha2)]*[X-halfW; Y-halfH];
        X=pos(1,:);
        Y=pos(2,:);
        plot(X, Y, 'r');
    else
        for j = 1 : length(pos_cell) 
            postemp = pos_cell{j};
            pos_cell{j} = [halfW; halfH]...
                + [cos(alpha2) sin(alpha2); -sin(alpha2) cos(alpha2)]...
                *[postemp(:,1)-halfW, postemp(:,2)-halfH]';           
            pos = pos_cell{j};
            Xcor = pos(1,:);
            Ycor = pos(2,:);
            p = plot(Xcor, Ycor, 'r');
        end
    end
    
    % Circle's radii
    radius_l = sqrt((corners4(1,1)-corners4(4,1))^2 + (corners4(1,2)-corners4(4,2))^2)/2;
    radius_r = sqrt((corners4(2,1)-corners4(3,1))^2 + (corners4(2,2)-corners4(3,2))^2)/2;
    % Circles' centers
    cl = [(corners4(1,1)+corners4(4,1))/2, (corners4(1,2)+corners4(4,2))/2];
    cr = [(corners4(2,1)+corners4(3,1))/2, (corners4(2,2)+corners4(3,2))/2];
    % Angle of orientation of the stadium
    alpha_l = acos((corners4(1,2)-corners4(4,2))/(2*radius_l));
    alpha_r = acos((corners4(2,2)-corners4(3,2))/(2*radius_r));
    % Angles of orientation for the semicircles around their centers
    thetas_l = linspace(pi/2-alpha_l, pi+alpha_l, 100);
    thetas_r = linspace(-pi/2-alpha_r, pi/2-alpha_r,100);

    % Whiten the surroundings of the stadium
    matx = repmat(1:W,H,1);
    maty = repmat(transpose(1:H),1,W);
    matx = reshape(matx, N, 1);
    maty = reshape(maty, N, 1);
    mask = inpolygon(matx,maty, corners4(:,1),corners4(:,2));
    mask = mask | ((matx-cl(1)).^2 + (maty-cl(2)).^2 < radius_l.^2);
    mask = mask | ((matx-cr(1)).^2 + (maty-cr(2)).^2 < radius_r.^2);
    mask = reshape(mask, H, W);
    if (length(size(im)==2)) % If not RGB
        imrot = reshape(imrot, N,1);
    else
        imrot = reshape(imrot, N, 3);
    end
    imrot(~mask, :) = 255;
    if (length(size(im)==2)) % If not RGB
        imrot = reshape(imrot, H, W, 1);
    else
        imrot = reshape(imrot, H, W, 3);
    end
    close(figure(14));
    
    %% Show final trajectory
    figure(14); hold on;
    imshow(imrot); hold on;
    % The two lines
    plot(corners4(1:2,1), corners4(1:2,2), 'k', 'LineWidth', 1.5); 
    plot(corners4(3:4,1), corners4(3:4,2), 'k', 'LineWidth', 1.5);
    % The two semicircles
    plot([cl(1)+radius_l*cos(thetas_r)], [cl(2)+radius_l*sin(thetas_r)], 'k', 'LineWidth', 1.5);
    plot([cr(1)-radius_r*cos(thetas_r)], [cr(2)-radius_r*sin(thetas_r)], 'k', 'LineWidth', 1.5);
    % The trajectory
    if (~is_cell)
        p = plot(X, Y, 'r');
        if (jetmap)
            n = length(X);
            cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
            drawnow;
            set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd);
        end
    else
        for j = 1 : length(pos_cell) 
            Xcor = pos_cell{j}(1,:);
            Ycor = pos_cell{j}(2,:);
            p = plot(Xcor, Ycor, 'r');
        end
    end
    axis([corners4(1,1)-1.1*radius_l, corners4(2,1)+1.1*radius_r,...
        corners4(1,2)-10, corners4(3,2)+10]);
    axis on; hold on;
    set(gca,'TickLabelInterpreter','Latex', 'FontSize', 14.0);
    
    %% Normalized trajectory
    fig15 = figure(15); hold on;
    set(fig15, 'position',[200 200 800 400]);
    
    % The two lines
    plot([-0.5, 0.5], [0.5, 0.5], 'k', 'LineWidth', 1.5); 
    plot([-0.5, 0.5], [-0.5, -0.5], 'k', 'LineWidth', 1.5);
    % The two semicircles
    r=0.5;
    thetas = linspace(-pi/2, pi/2, 100);
    plot([0.5+r.*cos(thetas)], [r.*sin(thetas)], 'k', 'LineWidth', 1.5);
    plot([-0.5-r.*cos(thetas)], [-r.*sin(thetas)], 'k', 'LineWidth', 1.5);
    % The trajectory
    c = (cl+cr)./2;
    W = 2.*(cl(1)-cr(1));
    H = corners4(1,2)-corners4(4,2);
    
    if (~is_cell)
        X=X-c(1);
        Y=Y-c(2);
        Xnorm = X.*2./W;
        Ynorm = Y./H;
        p = plot(Xnorm, Ynorm, 'r');
        if (jetmap)
            drawnow;
            set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
        end
        if (add_sym)
            p = plot(-Xnorm, -Ynorm, 'r');
            if (jetmap)
                drawnow;
                set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd);
            end
            Xnorm = [Xnorm, -Xnorm];
            Ynorm = [Ynorm, -Ynorm];
        end
    else
        for j = 1 : length(pos_cell) 
            pos_cell{j} = pos_cell{j}-c';
            Xnorm = pos_cell{j}(1,:).*2./W;
            Ynorm = pos_cell{j}(2,:)./H;
            p = plot(Xnorm, Ynorm, 'r');
            if (add_sym)
                p = plot(-Xnorm, -Ynorm, 'r');
                pos_cell{j} = [pos_cell{j}, -pos_cell{j}];
            end
        end
    end
    hold on;
    axis off;
    
    
    %% 2D histogram
    fig16 = figure(16); hold on;
    set(fig16, 'position',[600 200 800 400]);
    % Estimate a continuous pdf from the discrete data
    W = 120;
    H = 60;
    N = W*H;
    xi = linspace(-1,1,W);
    yi = linspace(-0.5,0.5,H);
    [xxi,yyi] = meshgrid(xi,yi);
    pdfxy = zeros(size(xxi));
    
    % Whiten the surroundings of the stadium
    matx = reshape(xxi, N, 1);
    maty = reshape(yyi, N, 1);
    mask = inpolygon(matx,maty, [-0.5,0.5,0.5,-0.5],[0.5,0.5,-0.5,-0.5]);
    mask = mask | ((matx+0.5).^2 + maty.^2 < (0.5).^2);
    mask = mask | ((matx-0.5).^2 + maty.^2 < (0.5).^2);
    mask = reshape(mask, size(xxi));

    if (is_cell)   
        Xnorm = []; Ynorm = [];
        for j = 1 : length(pos_cell) 
           pos_cell{j} = pos_cell{j}-c';
           Xnorm = [Xnorm, pos_cell{j}(1,:).*2./W];
           Ynorm = [Ynorm, pos_cell{j}(2,:)./H];
        end
    end
    
    for i = 1:length(Xnorm)
        [~,xmin] = min((Xnorm(i)-xxi(1,:)).^2);
        [~,ymin] = min((Ynorm(i)-yyi(:,1)).^2);
        pdfxy(ymin, xmin) = pdfxy(ymin, xmin)+1;
    end
    pdfxy = pdfxy./length(Xnorm);
    pdfxy(~mask)=nan;
    mesh(xxi,yyi,pdfxy, 'FaceColor', 'flat');
%         c = hot;
%         c = getPyPlot_cMap('Greys');
        c = getPyPlot_cMap('viridis');
%         c = flipud(c); % Eventually flip the colormap
    colormap(c);
        
    r=0.49;
    lwidth = 4;
    % The two lines
    plot([-0.5, 0.5], [r, r], 'k', 'LineWidth', lwidth); 
    plot([-0.5, 0.5], [-r, -r], 'k', 'LineWidth', lwidth);
    % The two semicircles
    thetas = linspace(-pi/2, pi/2, 100);
    plot([0.5+r.*cos(thetas)], [r.*sin(thetas)], 'k', 'LineWidth', lwidth);
    plot([-0.5-r.*cos(thetas)], [-r.*sin(thetas)], 'k', 'LineWidth', lwidth);
%     plot(xxi(~mask), yyi(~mask), '.w', 'MarkerSize',15);
    axis([-1.01, 1.01, -0.51, 0.51]);
    colorbar;
    axis off;
    
    %% Some characteristics
    % Compute center of trajectory
    xC = mean(Xnorm);
    yC = mean(Ynorm);
    
    dist_ul_ll = sqrt((corners4(1,1)-corners4(4,1)).^2 + (corners4(1,2)-corners4(4,2)).^2) ;
    resolution = ref_length./dist_ul_ll;
    
    if (~is_cell)
%         dt = 1/2; % Time between two frames.
        dt = 1/25;
        vx = (X(2:end) - X(1:end-1)).*resolution./dt; 
        vy = (Y(2:end) - Y(1:end-1)).*resolution./dt; 
        norm_v = sqrt(vx.^2+vy.^2).*1e3;
%         ind = find(norm_v>0.2);
%         norm_v(ind, : ) =[];
%         drop_position(ind, :)=[];
%         v(ind, :)=[];

%         p = m.*v; % Momentum
%         R = sqrt(sum( (drop_position-[xC, yC])'.^2 )).*pixel2mm; % Droplet distance to the center of its trajectory
%         r = (drop_position-[xC, yC]);
%         L = r(1:end-1,1).*p(:,2) - r(1:end-1,2).*p(:,1); % Angular momentum
        
        f3 = figure(3); hold on;
        set(f3, 'position',[600 200 550 500]);
        set(gca,'TickLabelInterpreter','Latex', 'FontSize', 12.0);
        hv = histogram(norm_v, 50);
        hv.FaceColor = 'b';
        % h.EdgeColor = 'None';
        xlabel('v [mm/s]', 'FontSize', 16.0, 'interpreter', 'latex');
        ylabel('# of occurrences', 'FontSize', 16.0, 'interpreter', 'latex');
        % legend('Histogram of the droplet speed', 'FontSize', 14.0, 'interpreter', 'latex');
        
        speed = mean(norm_v) % [mm/s] 
        G_div_GF = (M-1)/M;
        diameter = diam_from_speed (G_div_GF, speed);
    end
    
end

