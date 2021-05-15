%% Computes the rank of an all-pole noisy signal mapped into a Hankel matrix with various vector length.
%
%
%HISTORY:
% 2021-05-04: Lucas Abdalah
%

%% Initialization
close all;
clearvars; clc;

%% General Setup
Fs = 1e3;       % Sampling Frequency (Hz)
Fn = Fs/2;      % Nyquist Frequency (Hz)
T = 1/Fs;       % Sampling period
% N = 1e2;        % Length of signal
N = 1e2;        % Length of signal
t = (0:N-1)*T;  % Time vector

%% Noise values
noise = [310, 50];

%% Poles setup
itPoles = 2:1:8;    % Number of poles
spacings = 100;     % Space between poles
wAll = logspace(-4,-1,spacings);

%% Info storing
hankelSample = zeros(length(itPoles),length(wAll)); % Store the number of samples to build a Hankel matrice with rank = number of poles
LHsAll = zeros(length(itPoles),length(wAll));       % Store rank LHs Hankel

%% Main algorithm
for ii = 1:1:length(itPoles)
    
    Nexp = itPoles(ii); % Number of exponential
    fprintf(['.............. Rank => ', num2str(Nexp),'..............\n']) 
    
    for jj = 1:1:length(wAll)
        
        %% Poles Setup
        wn = wAll(jj);      % Normalized frequency distance  \Delta \omega
        
        %% Sum of Complex Exponentials
        f = [linspace(0, wn*(Fn),Nexp)'];   % Linear Frequency
        comp = exp(1j*2*pi*f.*t);           % Poles
        s=sum(comp);  % Signal
        
        %% Noisy scenario with SNR = 310dB;
        y = awgn(s,noise(1),'measured');

        %% Study - Hankel
        HsRank = zeros(N,1);    % Store the rank for a signal with jj samples
        Hs = Hankelize(y);      % Hankelize the whole signal
        for kk=1:1:N % Computing Hankel matrix rank for a signal with jj samples 
            HsRank(kk,1) = rank(Hankelize(y(1:kk)));  % Rank computing with s(1:jj)
        end
       
        %% Saving Hankel sample and current LHs rank
        if sum(find(HsRank == Nexp,1)) == 0 % Verify if the distance is enough to find the rank = number of poles
            hankelSample(ii,jj) = 150; % Replace when N isn't enough to find the rank = number of poles <======= here
        else
            hankelSample(ii,jj) = find(HsRank == Nexp,1);
        end

    end

    str = ['L = ', num2str(Nexp)]; % Save the legend
    leg{ii} = str;

end

%% Plot results - SNR = 310
h = figure;
subplot(2,1,1)
hold on
myColor = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};

myLineStyle = {'--','-',':','-.'};
myMarker = {'o', 'square', '^', 'x'};
for ii = 1:1:length(itPoles)
    semilogx(wAll,hankelSample(ii,:),...
        'Color', myColor{ii},...        
        'LineStyle', myLineStyle{rem(ii,4)+1},...
        'LineWidth',1.0,...
        'Marker',myMarker{rem(ii,4)+1},...
        'MarkerFaceColor', myColor{ii},...
        'MarkerSize', 3);
end
hold off
xlabel('\Delta\omega (\times \pi radians/samples)');
ylabel('N (Sample size)');
ylim([0 80])
% legend(leg,'Location', 'bestoutside') % legend(leg,'Location', 'Best')
grid

%% Main algorithm 2
for ii = 1:1:length(itPoles)
    
    Nexp = itPoles(ii); % Number of exponential
    fprintf(['.............. Rank => ', num2str(Nexp),'..............\n']) 
    
    for jj = 1:1:length(wAll)
        
        %% Poles Setup
        wn = wAll(jj);      % Normalized frequency distance  \Delta \omega
        
        %% Sum of Complex 
        f = [linspace(0, wn*(Fn),Nexp)'];   % Linear Frequency
        comp = exp(1j*2*pi*f.*t);           % Poles
        s=sum(comp);  % Signal
        
        %% Noisy scenario with SNR = 50dB;
        y = awgn(s,noise(2),'measured');

        %% Study - Hankel
        HsRank = zeros(N,1);    % Store the rank for a signal with jj samples
        Hs = Hankelize(y);      % Hankelize the whole signal
        for kk=1:1:N % Computing Hankel matrix rank for a signal with jj samples 
            HsRank(kk,1) = rank(Hankelize(y(1:kk)));  % Rank computing with s(1:jj)
        end
       
        %% Saving Hankel sample and current LHs rank
        if sum(find(HsRank == Nexp,1)) == 0 % Verify if the distance is enough to find the rank = number of poles
            hankelSample(ii,jj) = 150; % Replace when N isn't enough to find the rank = number of poles <======= here
        else
            hankelSample(ii,jj) = find(HsRank == Nexp,1);
        end

    end

    str = ['L = ', num2str(Nexp)]; % Save the legend
    leg{ii} = str;

end

%% Plot results - SNR = 50
subplot(2,1,2)
hold on
myColor = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};

myLineStyle = {'--','-',':','-.'};
myMarker = {'o', 'square', '^', 'x'};
for ii = 1:1:length(itPoles)
    semilogx(wAll,hankelSample(ii,:),...
        'Color', myColor{ii},...        
        'LineStyle', myLineStyle{rem(ii,4)+1},...
        'LineWidth',1.0,...
        'Marker',myMarker{rem(ii,4)+1},...
        'MarkerFaceColor', myColor{ii},...
        'MarkerSize', 3);
end
hold off
xlabel('\Delta\omega (\times \pi radians/samples)');
ylabel('N (Sample size)');
ylim([0 20]);
yticks([0:5:20]);
legend(leg,'Location', 'northeastoutside','Orientation', 'Vertical');
% lgd = legend(leg,'Location', 'northeastoutside','Orientation', 'Vertical','Position', [0.95 0. 0.0 1], 'Units','normalized');
lgd = legend(leg,'Location', 'northeastoutside','Orientation', 'Horizontal','Position', [0.5 0.47 0.0 1], 'Units','normalized');
legend boxoff
grid


%% Hankel samples (rank = poles) in a noisy scenario - Export Result Figure
savefigPath = 'C:\Users\lukin\Documents\GitHub\Hankel-Matrices\matlab\fig';
outNamePDF = 'noisy_proximity-PoleDistance_vs_Samples4';
datetimestamp = datestr(now, 'yyyy-mm-dd');
outfilename = strjoin({[savefigPath '\' outNamePDF '-' datetimestamp]});
savefig_tight(h,outfilename);