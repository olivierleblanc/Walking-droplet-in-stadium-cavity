function [positions, thetas, sins] = classorbits( xi, yi, ri, si, numbounds, colorplot )

close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Olivier Leblanc
% Date : 24/05/2020
%
% Function :
% Allows to simulate the trajectories of an inelastic ball inside a
% circular, elliptical or stadium geometry by choosing the initial
% conditions. Plots the result both in the space and in the phase space
%
% Inputs :
% - (xi, yi) : 2D coordinates of the initial position
% - (ri, si) : 2D coordinates of the initial direction of motion
% - numbounds : number of bounds on the borders
% - colorplot : color of the line describing the trajectories
%
% Outputs :
% /.
%
%
% Examples
%
% %%%% The simple horizontal line %%%% 
% [positions, thetas, sins] = classorbit( 0, 0, 1, 0, 2, 'b' );
%
% %%%% The diamond %%%% 
% [positions, thetas, sins] = classorbit( 0, -1, 2, 1, 16, 'b' );
%
% %%%% The rectangle %%%% 
% [positions, thetas, sins] = classorbit( 1+sqrt(2)/2, 0, 0, 1, 10, 'b' );
%
% %%%% Few iterations in chaotic case %%%% 
% [positions, thetas, sins] = classorbit( 0.5, 0, 0.1, 1, 100, 'b' );
%
% %%%% Lots of iterations in chaotic case %%%% 
% [positions, thetas, sins] = classorbit( 0, 0, 1, 1, 5000, 'b' );
%
% [positions, thetas, sins] = classorbit( 0, 0.99, 1, 0, 16, 'b' );
% 
%
% Options
cavity = 'stadium'; % 'ellipse' or 'stadium'
yax = 'sin'; % Choose between 'sin' and 'cos'
colorscatter = 'r';
scattersize = 10;
dev = 1;  % Distance between the two semicircles for a stadium.  (0 for circle)
excplus = 0; % Difference from one for the 'a' term of ellipse eq.: (x/a)^2+y^2=1

show_subplots = 0; % Shows two examples for the report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Code

p = xi;
q = yi;
r=ri;
s=si;

% axes('position',[0 0 1 1]);
fig1 = figure(1); hold on;
set(gca,'TickLabelInterpreter','Latex', 'FontSize', 16.0);
set(fig1, 'Visible','off');
xlabel('x', 'interpreter', 'latex', 'FontSize', 18.0);
ylabel('y', 'interpreter', 'latex', 'FontSize', 18.0);
title('Trajectories in the geometry', 'interpreter', 'latex', 'FontSize', 16.0);

% Plot the cavity
if (strcmp(cavity, 'stadium'))
    set(fig1, 'position',[50 400 300+dev*300 300]); hold on;
    axis([-1-dev 1+dev -1 1]);
    plot([-dev dev],[1,1],'k');
    plot([-dev dev],[-1,-1],'k');
    y=linspace(-1,1,100);
    xx=-sqrt(1.-y.^2)-dev;
    plot(xx,y,'k');
    xx=sqrt(1.-y.^2)+dev;
    plot(xx,y,'k');
elseif (strcmp(cavity, 'ellipse'))
    set(fig1, 'position',[50 400 300+excplus*300 300]); hold on;
    axis([-1-excplus 1+excplus -1 1]);
    y=linspace(-1,1,100);
    xx=sqrt((1+excplus)^2*(1.-y.^2));
    plot(xx,y,'k');
    xx=-sqrt((1+excplus)^2*(1.-y.^2));
    plot(xx,y,'k');
end
hold on;

% Phase-space plot
fig2 = figure(2); hold on;
axis([-0.5 3.5 -1.2 1.2]);
grid on;
set(gca,'TickLabelInterpreter','Latex', 'FontSize', 18.0);
set(fig2, 'position',[800 200 500 500]);
set(fig2, 'Visible','off');
xlabel('$\theta$', 'interpreter', 'latex', 'FontSize', 20.0);
if (strcmp(yax, 'sin'))
    ylabel('$\sin (i)$', 'interpreter', 'latex', 'FontSize', 20.0);
elseif (strcmp(yax, 'cos'))
    ylabel('$ \cos (i)$', 'interpreter', 'latex', 'FontSize', 20.0);
end
title('Phase space plot', 'interpreter', 'latex', 'FontSize', 18.0);
hold on;
sz = 25;
c = linspace(1,10,numbounds);

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (numbounds < 75000)
    % Preallocation
    positions = zeros(2, numbounds+1);
    thetas = zeros(1, numbounds);
    sins = zeros(1, numbounds);

    % Initialization
    positions(:,1) = [p;q];
    
    set(0, 'CurrentFigure', fig1);
    for iter=1:numbounds
        % Find data for the next bound
        if (strcmp(cavity, 'stadium'))
            [po qo ro so sinus]=nextref(p,q,r,s);
        elseif (strcmp(cavity, 'ellipse'))
            [po qo ro so sinus]=nextellipse(p,q,r,s);
        end
            
        % Save current data
        thetas(iter) = atan2(qo,po);
        if (strcmp(yax, 'sin'))
            sins(iter) = sinus;
        elseif (strcmp(yax, 'cos'))
            sins(iter) = cos(asin(sinus));
        end
        positions(: , iter+1) = [po; qo];
        plot([p po], [q qo], colorplot, 'LineWidth', 0.05);
        
        % Update current condition
        p = po;
        q = qo;
        r=ro;
        s=so;
    end
    
    set(0, 'CurrentFigure', fig2);
    % scatter(thetas, sins,sz, c,'filled');
    scatter(thetas, sins, scattersize, colorscatter,'filled');
    
else
    
    %         for iter=1:numbounds
    %             % Find data for the next bound
    %             [po qo ro so sinus]=nextref(p,q,r,s);
    %             theta = atan2(qo,po);;
    %
    %             set(0, 'CurrentFigure', fig1);
    %             plot([p po], [q qo], colorplot, 'LineWidth', 1.0);
    %
    %             % scatter(thetas, sins,sz, c,'filled');
    %             set(0, 'CurrentFigure', fig2);
    %             scatter(theta, sinus, 2, 'r','filled');
    %
    %             % Update current condition
    %             p = po;
    %             q = qo;
    %             r=ro;
    %             s=so;
    %         end
    
    % --------- Faster to have two loops to avoid switching --------------
    % --------- between the two figures at each iteration ----------------
    % Additional note : applying scatter at each iteration is very very slow.
    
    set(0, 'CurrentFigure', fig1);
    for iter=1:numbounds
        % Find data for the next bound
        [po qo ro so sinus]=nextref(p,q,r,s);
        plot([p po], [q qo], colorplot, 'LineWidth', 1.0);
        
        % Update current condition
        p = po;
        q = qo;
        r=ro;
        s=so;
    end
    
    set(0, 'CurrentFigure', fig2);
    for iter=1:numbounds
        % Find data for the next bound
        [po qo ro so sinus]=nextref(p,q,r,s);
        theta = atan2(qo,po);
        scatter(theta, sinus, 2, 'r', 'filled');
        
        % Update current condition
        p = po;
        q = qo;
        r=ro;
        s=so;
    end
       
    % To return smt
    positions = 0;
    thetas = 0;
    sins = 0;
    
end

set(fig1, 'Visible','on');
set(fig2, 'Visible','on');
toc;

if (show_subplots)
%%%%%%%%%%%%%%%%%%%%%% Subplots %%%%%%%%%%%%%%%%%%%%%%%
fig3 = figure(3); hold on;
set(gca,'TickLabelInterpreter','Latex', 'FontSize', 12.0);

% stadiums
subplot(2,2,1); hold on;
title('Trajectories in the geometry', 'interpreter', 'latex');
    plot([-dev dev],[1,1],'k');
    plot([-dev dev],[-1,-1],'k');
    y=linspace(-1,1,100);
    xx=-sqrt(1.-y.^2)-dev;
    plot(xx,y,'k');
    xx=sqrt(1.-y.^2)+dev;
    plot(xx,y,'k');
    ylabel('y', 'interpreter', 'latex');
    
subplot(2,2,3); hold on;
    plot([-dev dev],[1,1],'k');
    plot([-dev dev],[-1,-1],'k');
    y=linspace(-1,1,100);
    xx=-sqrt(1.-y.^2)-dev;
    plot(xx,y,'k');
    xx=sqrt(1.-y.^2)+dev;
    plot(xx,y,'k');
xlabel('x', 'interpreter', 'latex');
ylabel('y', 'interpreter', 'latex');

% First orbit
subplot(2,2,1);
thetas = zeros(1, 2);
sins = zeros(1, 2);
for iter=1:2
    % Find data for the next bound
    [po qo ro so sinus]=nextref(p,q,r,s); 
    % Save current data
    thetas(iter) = atan2(qo,po);
    sins(iter) = sinus;
    plot([p po], [q qo], colorplot, 'LineWidth', 1.0);
    % Update current condition
    p = po;
    q = qo;
    r=ro;
    s=so;
end

% First phase space
subplot(2,2,2); grid on; axis([-0.5 3.5 -1.2 1.2]); hold on;
scatter(thetas, sins, scattersize, colorscatter,'filled');
title('Phase space plot', 'interpreter', 'latex');
ylabel('$\sin (i)$', 'interpreter', 'latex');

% Second orbits
p = 0; q = -1; r = 2; s = 1;
subplot(2,2,3);
thetas = zeros(1, 4);
sins = zeros(1, 4);
for iter=1:4
    % Find data for the next bound
    [po qo ro so sinus]=nextref(p,q,r,s); 
    % Save current data
    thetas(iter) = atan2(qo,po);
    sins(iter) = sinus;
    plot([p po], [q qo], colorplot, 'LineWidth', 1.0);
    % Update current condition
    p = po;
    q = qo;
    r=ro;
    s=so;
end

% Second phase space
subplot(2,2,4); grid on; axis([-0.5 3.5 -1.2 1.2]); hold on;
scatter(thetas, sins, scattersize, colorscatter,'filled');
xlabel('$\theta$', 'interpreter', 'latex');
ylabel('$\sin (i)$', 'interpreter', 'latex');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function taken from "Qstad library" and modified to computed sin_i and
% handle any spacing between the two semicircles.
    function [po qo ro so sin_i]=nextref(p,q,r,s)
        if abs(s)>.000001 % Condition car on veut pas diviser par 0
            up=p+r/s*(sign(s)-q); %finds the x point point on y=+/-1 on the ray
            if up>=-dev&&up<=dev %and reflects back if hitting the straight edge
                po=up;
                qo=sign(s);
                ro=r;
                so=-s;
                %----------------------------------------
                % Computation of sin_i
                uvec = [r s 0];
                vvec = [0 -qo 0];
                crossprod2 = cross(uvec,vvec);
                i2 = (atan2(crossprod2(3),-dot(uvec,vvec)));
                sin_i = sin(i2);
                %----------------------------------------
                
            else %we are hitting a curved portion
                sp=sign(up); %sets left or right circle
                frac=r/s;
                a2 = frac^2+1;
                b2 = 2*frac*(p-dev*sp-q*frac);
                cst = (p-dev*sp-frac*q)^2-1;
                qo = (-b2+sign(s)*sqrt(b2^2-4*a2*cst))/(2*a2);
                po = p+frac*(qo-q);
                
                %----------------------------------------
                % Computation of sin_i
                uvec = [r s 0];
                vvec = [-po+dev*sign(po) -qo 0];
                crossprod2 = cross(uvec,vvec);
                i2 = (atan2(crossprod2(3),-dot(uvec,vvec)));
                sin_i = sin(i2);
                %----------------------------------------
                
                P = [cos(i2) sin_i; -sin_i cos(i2)];
                val = P*[-po+dev*sign(po); -qo];
                ro = val(1);
                so = val(2);
                
            end
        else
            qo = q;
            po = sign(r)*sqrt(1-qo^2)+dev*sign(r);
            
            
            %----------------------------------------
            % Computation of sin_i
            uvec = [r s 0];
            vvec = [-po+dev*sign(po) -qo 0];
            crossprod2 = cross(uvec,vvec);
            i2 = (atan2(crossprod2(3),-dot(uvec,vvec)));
            sin_i = sin(i2);
            %----------------------------------------
            
            P = [cos(i2) sin_i; -sin_i cos(i2)];
            val = P*[-po+dev*sign(po); -qo];
            ro = val(1);
            so = val(2);
        end
        
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [po qo ro so sin_i]=nextellipse(p,q,r,s)
        if abs(s)>.000001 % Condition car on veut pas diviser par 0
                frac=r/s;
                a2 = frac^2+(1+excplus)^2;
                b2 = 2*frac*(p-q*frac);
                cst = (p-frac*q)^2-(1+excplus)^2;
                qo = (-b2+sign(s)*sqrt(b2^2-4*a2*cst))/(2*a2);
                po = p+frac*(qo-q);
        else
            qo = q;
            po = sign(r)*sqrt((1+excplus)^2*(1-qo^2));
        end  
            
            %----------------------------------------
            % Computation of sin_i
            uvec = [r s 0];
            vvec = [-po/(1+excplus)^2 -qo 0];
            crossprod2 = cross(uvec,vvec);
            i2 = (atan2(crossprod2(3),-dot(uvec,vvec)));
            sin_i = sin(i2);
            %----------------------------------------
            
            P = [cos(i2) sin_i; -sin_i cos(i2)];
            val = P*[-po/(1+excplus)^2; -qo];
            ro = val(1);
            so = val(2);
    end

end
