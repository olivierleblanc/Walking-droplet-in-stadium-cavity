function [positions, thetas, sins] = classorbits( xi, yi, ri, si, numbounds, colorplot )

close all;
clc;

%% Examples

% [positions, thetas, sins] = classorbit( 0.5, 0, 0.1, 1, 100, 'b' );

% [positions, thetas, sins] = classorbit( 0, 0, 1, 1, 5000, 'b' );

% [positions, thetas, sins] = classorbit( 0, 0.99, 1, 0, 16, 'b' );

%% Options

cavity = 'stadium'; % 'ellipse' or 'stadium'
yax = 'sin'; % Choose between 'sin' and 'cos'
colorscatter = 'g';
scattersize = 20;
dev = 1;  % Distance between the two semicircles for a stadium.  (0 for circle)
excplus = 0; % Difference from one for the 'a' term of ellipse eq.: (x/a)^2+y^2=1

%% Code

p = xi;
q = yi;
r=ri;
s=si;

% axes('position',[0 0 1 1]);
fig1 = figure(1); hold on;
set(fig1, 'Visible','off');
xlabel('x');
ylabel('y');
title('Trajectories in the geometry');

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
set(fig2, 'position',[800 200 500 500]);
set(fig2, 'Visible','off');
xlabel('\theta');
if (strcmp(yax, 'sin'))
    ylabel('sin(i)');
elseif (strcmp(yax, 'cos'))
    ylabel('cos(i)');
end
title('Phase space plot');
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
        plot([p po], [q qo], colorplot, 'LineWidth', 1.0);
        
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
