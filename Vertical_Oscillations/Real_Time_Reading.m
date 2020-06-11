%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Olivier Leblanc
% Date : 17/05/2020
%
% Function :
% Reads and plots the values received from analog input A0
% from the Arduino Uno in real-time. "AnalogReadSerial" code in the
% examples of Arduino IDE must run and the serial monitor must be closed.
%
% Inputs :
maxVal = 400;           % number of values in plot 
n_mult = 1;             % number of times we plot maxVal measurements
baudrate = 250000;      % baudrate imposed by the Arduino Uno code
%
% Outputs :
% /.
%
% Options :
newreading = 1;         % If 0 : only show the FFT of already saved data
True_meas = 1;          % If 0, use artificial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc;

%% Reading of new values

if (newreading)
    % Some initializations
    count = 1;
    Acc_X = 1000*ones(1,maxVal);
    Acc_Y = zeros(1,maxVal);
    Acc_Z = -1000*ones(1,maxVal);
    up = zeros(1,maxVal); % Enveloppe of Acc_Z
    amplitude_estimate = zeros(1,maxVal);
    fs = 1/(1470*10^(-6));  % Sampling frequency --> Obtained from Arduino Uno CODE
    tmax = maxVal/fs;
    f = (-round(maxVal/2)+1:round(maxVal/2)).*fs./(maxVal);
    dist_to_80 = (f-80).^2;
    [~, ind80] = min(dist_to_80);
    n_harm = floor(maxVal/(2*(ind80-round(maxVal/2)))) - 1;
    ind_harm = zeros(1,n_harm);
    for i=1:n_harm
        dist = (f-(i+1)*80).^2;
        [~, ind_harm(i)] = min(dist);
    end
    ind_harm = ind_harm + 1; % Correction
    
    % Figure configurations
    fig2 = figure(2); hold on;
    set(fig2, 'position',[800 200 500 500]);
    %     xlim([0 fs/2]);
    xlabel('f [Hz]');
    ylabel('Acceleration spectrum [mg]');
    
    fig1 = figure(1);
    set(fig1, 'position',[200 200 500 500]);
    xlim([0 tmax]);
    ylim([-8000 8000]);
    xlabel('t [s]');
    ylabel('Acceleration [mg]');
    
    % Creation of lines and its properties
    x=1:maxVal;
    t = (1:maxVal)./fs;
    if (~True_meas)
        %        mat = [500;1000;1500].*(sin(2*pi*80.*t)+sin(2*pi*160.*t)+sin(2*pi*240.*t));
        mat = [500;1000;1500].*(sin(2*pi*80.*t).*(1+0.1.*(sin(2*pi*30.*t))));
    end
    lhx = line(t,Acc_X);
    lhy = line(t,Acc_Y);
    lhz = line(t,Acc_Z);
    line_amp = line(t, up);
    line_mean = line(t, amplitude_estimate);
    set(lhx,'LineWidth',1.5);
    set(lhy,'LineWidth',1.5);
    set(lhz,'LineWidth',1.5);
    set(line_amp,'LineWidth',1.0);
    set(line_mean,'LineWidth',1.2);
    lhx.Color = 'b';
    lhy.Color = 'g';
    lhz.Color = 'r';
    line_amp.Color = 'k';
    line_mean.Color = 'c';
    set(lhx, 'DisplayName', 'Acc-X');
    set(lhy, 'DisplayName', 'Acc-Y');
    set(lhz, 'DisplayName', 'Acc-Z');
    set(line_amp, 'DisplayName', 'Enveloppe of Acc-Z');
    set(line_mean, 'DisplayName', 'Mean amplitude of Acc-Z');
    legend();
    
    if (True_meas)
        % Creation and opening of the serial port
        instrreset;
        s=serial('COM3', 'BaudRate', baudrate);
        fopen(s);
    end
    
    %     tic
    times = zeros(1,maxVal);
    count = 1;
    mult = 1;
    while (mult < n_mult+1)
        % tic;
        
        if(True_meas)
            axis = str2num(fscanf(s));
        else axis = mat(:,count);
        end
        if (~isempty(axis) && length(axis)==3)
            Acc_X(count) = axis(1);
            Acc_Y(count) = axis(2);
            Acc_Z(count) = axis(3);
            %             times(count) = toc;
            count=count+1;
            
            % Plot
            set(lhx, 'ydata', Acc_X);
            set(lhy, 'ydata', Acc_Y);
            set(lhz, 'ydata', Acc_Z);
            
            % Obtention signal envelope to compute the mean amplitude of the oscillations
            ind = count-20:count;
            ind(ind<1) = ind(ind<1)+maxVal; % Wrap ind values inside [1, maxVal]
            %                 [up(ind), ~] = envelope(Acc_Z(ind), 5, 'peak'); % maybe add 'peak'.. 10 is the number of samples between maxima
            [up(ind), ~] = max(Acc_Z(ind));
            set(line_amp, 'ydata', up);
            if (mod(count, 20)==0)
                amplitude_estimate(ind) = mean(up(ind))-1000; % [mg]
                set(line_mean, 'ydata', amplitude_estimate+1000);
            end
            %                 plot(t, amplitude_estimate.*ones(1,length(t)), 'k', 'LineWidth', 2.0);
            
            if (count==maxVal) % Reset
                % Spectrum analysis
                Acc_X_fourier = abs(fftshift(fft(Acc_X)));
                Acc_Y_fourier = abs(fftshift(fft(Acc_Y)));
                Acc_Z_fourier = abs(fftshift(fft(Acc_Z)));
                
                % Cut DC value
                Acc_X_fourier(round(maxVal/2)+1)=0;
                Acc_Y_fourier(round(maxVal/2)+1)=0;
                Acc_Z_fourier(round(maxVal/2)+1)=0;
                figure(2);
                plot(f,Acc_X_fourier(1:maxVal), 'm', 'LineWidth', 2.0);
                plot(f,Acc_Y_fourier(1:maxVal), 'g', 'LineWidth', 2.0);
                plot(f,Acc_Z_fourier(1:maxVal), 'r', 'LineWidth', 2.0)
                legend('FFT of Acc_X','FFT of Acc_Y', 'FFT of Acc_Z');
                
                % Taux de distorsion harmonique sur les maxVal dernières mesures
                harm_distorsion = (sqrt(sum(Acc_Z_fourier(ind_harm).^2))/Acc_Z_fourier(ind80));
                
                count=1;
                mult = mult+1;
            end
        end
        drawnow;
    end
    
    if (True_meas)
        fclose(s);
    end
end

%% Spectrum analysis
Acc_X_fourier = abs(fftshift(fft(Acc_X)));
Acc_Y_fourier = abs(fftshift(fft(Acc_Y)));
Acc_Z_fourier = abs(fftshift(fft(Acc_Z)));

% Cut DC value
Acc_X_fourier(round(maxVal/2)+1)=0;
Acc_Y_fourier(round(maxVal/2)+1)=0;
Acc_Z_fourier(round(maxVal/2)+1)=0;

x = (-round(maxVal/2)+1:round(maxVal/2)-1).*fs./(maxVal);

fig2 = figure(2); hold on;
set(fig2, 'position',[800 200 500 500]);
xlim([0 fs/2]);
xlabel('f [Hz]');
ylabel('Acceleration spectrum [mg]')
plot(x,Acc_X_fourier(1:maxVal-1), 'm', 'LineWidth', 2.0);
plot(x,Acc_Y_fourier(1:maxVal-1), 'g', 'LineWidth', 2.0);
plot(x,Acc_Z_fourier(1:maxVal-1), 'r', 'LineWidth', 2.0)
legend('FFT of Acc_X','FFT of Acc_Y', 'FFT of Acc_Z');

