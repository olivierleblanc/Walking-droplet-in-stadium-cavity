% Intensity variation wrt Fresnel coeffs
% Author : Olivier Leblanc
% Date : 9-11-2019

close all;
clc;

n1 = 1;
n2 = 1.4;

theta_B = atan(1.4)*180/pi; % Brewster angle [°]

dtheta_max = 0.7578*pi/180;  % Maximum deviation angle

%% Light intensity of the parallel component at rest
theta1 = linspace(0,90,900)*pi/180;    % Incidence angle from vertical line
theta2 = asin(sin(theta1)./n2);    % Transmission angle from vertical line

E1r_par = (n1.*cos(theta2)-n2.*cos(theta1))./(n1.*cos(theta2)+n2.*cos(theta1));
I1r_par = E1r_par.^2;

figure(4); hold on;
plot(theta1.*180/pi, I1r_par, 'b', 'LineWidth', 1.2);
xlabel('$ \theta_1 [^\circ]$', 'interpreter', 'latex', 'FontSize', 14.0);
ylabel('$I_0/I_i$', 'interpreter', 'latex', 'FontSize', 14.0);
title('\textbf{Intensity reflected with the parallel component}', 'interpreter', 'latex');
set(gca,'TickLabelInterpreter','Latex', 'FontSize', 14.0);

%% Test for fixed maximum slope and varying incidence angle

theta1 = linspace(0,90,900)*pi/180;    % Incidence angle from vertical line
theta1_true = theta1 + dtheta_max;

% Coefficients for big slope
theta2 = -dtheta_max + asin(sin(theta1_true)./n2);    % Transmission angle from vertical line

theta2_true = theta2 + dtheta_max;
E1r_par = (n1.*cos(theta2_true)-n2.*cos(theta1_true))./(n1.*cos(theta2_true)+n2.*cos(theta1_true));
E1r_per = (n1.*cos(theta1_true)-n2.*cos(theta2_true))./(n1.*cos(theta1_true)+n2.*cos(theta2_true));

% Coefficients for flat surface
theta2_flat = asin(sin(theta1)./n2);
E1r_par_flat = (n1.*cos(theta2_flat)-n2.*cos(theta1))./(n1.*cos(theta2_flat)+n2.*cos(theta1));
E1r_per_flat = (n1.*cos(theta1)-n2.*cos(theta2_flat))./(n1.*cos(theta1)+n2.*cos(theta2_flat));

I_var_par = (E1r_par./E1r_par_flat).^2;
I_var_per = (E1r_per./E1r_per_flat).^2;



% Plot the result
figure(1); hold on;
plot(theta1*180/pi, I_var_per, 'r', 'LineWidth', 1.5);
xlabel('\theta_1 [°]', 'FontSize', 14.0);
ylabel('I_{\Delta}/I_0', 'FontSize', 14.0);
title('Maximum ratio for perpendicular component');
legend(['\Delta_{\theta} = ', num2str(dtheta_max*180/pi), '°'])
set(gca,'TickLabelInterpreter','Latex', 'FontSize', 14.0);

f2 = figure(2); hold on;
plot(theta1*180/pi, I_var_par, 'b', 'LineWidth', 1.5);
plot([theta_B theta_B+2*dtheta_max*180/pi],[5 5], 'r', 'LineWidth', 1.5);
xlabel('$\theta_1 [^\circ]$','Interpreter', 'Latex', 'FontSize', 14.0);
ylabel('$I_{\Delta}/I_0$', 'Interpreter', 'Latex', 'FontSize', 14.0);
set(gca,'TickLabelInterpreter','Latex', 'FontSize', 12.0);
title('Maximum ratio for parallel component');
legend(['\Delta_{\theta} = ', num2str(dtheta_max*180/pi), '°'], 'range', 'FontSize', 12.0);

hf2 = axes('parent',f2,'position',[0.2 0.5 0.3 0.3]); % normalized units are used for position
                                                     % parent - handle to main figure
                                                     % position - [x,y,width,height] of inset
plot(theta1*180/pi, I_var_par, 'b', 'LineWidth', 1.5); hold on;
plot([theta_B theta_B+2*dtheta_max*180/pi],[5 5], 'r', 'LineWidth', 1.5);
axis([54 58 0 70])
set(gca,'TickLabelInterpreter','Latex', 'FontSize', 12.0);




%% Test for fixed incidence angle and varying slope

theta1 = 55.22*pi/180;    % Incidence angle from vertical line

dtheta = linspace(-dtheta_max, dtheta_max, 1000);
% dtheta = linspace(-90, 90, 1000)*pi/180;

theta1_true = theta1 + dtheta;

% Coefficients for all slopes
theta2 = -dtheta + asin(sin(theta1_true)./n2);
theta2_true = theta2 + dtheta;
E1r_par = (n1.*cos(theta2_true)-n2.*cos(theta1_true))./(n1.*cos(theta2_true)+n2.*cos(theta1_true));
E1r_per = (n1.*cos(theta1_true)-n2.*cos(theta2_true))./(n1.*cos(theta1_true)+n2.*cos(theta2_true));

% Coefficients for flat surface
theta2_flat = asin(sin(theta1)./n2);
E1r_par_flat = (n1.*cos(theta2_flat)-n2.*cos(theta1))./(n1.*cos(theta2_flat)+n2.*cos(theta1));
E1r_per_flat = (n1.*cos(theta1)-n2.*cos(theta2_flat))./(n1.*cos(theta1)+n2.*cos(theta2_flat));

I_var_par = (E1r_par./E1r_par_flat).^2;
I_var_per = (E1r_per./E1r_per_flat).^2;



% Plot the result
figure(3); hold on;
plot(dtheta*180/pi, I_var_par, 'b', 'LineWidth', 1.5);
plot(dtheta*180/pi, I_var_per, 'r', 'LineWidth', 1.5);
xlabel('$\Delta_{\theta} [^\circ]$', 'interpreter', 'latex', 'FontSize', 12.0);
ylabel('$I_{\Delta}/I_0$', 'Interpreter', 'Latex', 'FontSize', 14.0);
legend('parallel', 'perpendicular', 'FontSize', 16.0, 'Location', 'NorthWest');
title(['Variation of intensity wrt \Delta_{\theta} for \theta_1 = ', num2str(theta1*180/pi), '°']);
set(gca,'TickLabelInterpreter','Latex', 'FontSize', 14.0);


% %% Correlation based visualization
% 
% 
% alpha = 1-n1/n2;
% H = 1.3;          % Distance between camera and dotted pattern
% dz_dx = 0.0132277;  % Maximum slope
% h0 = 4.75e-3;   % Liquid heigt
% hg = 1e-3;  % Glass height
% ng = 1.51; % Glass refractive index
% 
% 
% ha = linspace(0,2e-2,100);
% 
% hp = h0 + n2./ng.*hg + n2./n1.*ha;
% 
% dx = dz_dx.*(alpha.*hp.*H)./(H-alpha.*hp);
% 
% figure(4); hold on;
% plot(ha, dx, 'b', 'LineWidth', 1.5);
% xlabel('ha');
% ylabel('dx');
% % legend('');




















