%% For fixed \Delta\omega = 0.05, we build a signal with L Construir um sinal de L polos.. Percorrer as matrizes de Hankel e calcular os N samples necessários para obter R=L.
%
%
%HISTORY:
% 2021-05-15: Lucas Abdalah
%

%% Preparing Matlab ambient
close all; clearvars; clc;
%% Stopwacth
tStart = tic;

%% General Setup
Fs = 1e3;       % Sampling Frequency (Hz)
Fn = Fs/2;      % Nyquist Frequency (Hz)
T = 1/Fs;       % Sampling period
% N = 1e2;        % Length of signal
N = 1e2;        % Length of signal
t = (0:N-1)*T;  % Time vector

%% Poles setup
L = 2:1:8;  % Number of poles
wn = 5e-2;   % For a fixed omega instead -> wAll = logspace(-4,-1,spacings);
%% Noise values
SNRnoise = 10:5:350;
%% Monte Carlo Scenario
Nrun = 100;
yRMC = zeros(Nrun,length(L),length(SNRnoise));

%% Monte Carlo Run
for rmc = 1:1:Nrun
    fprintf('... Monte Carlo Run -> %d \n', rmc);
    % fprintf('||  L | SNR || \n');

    %% Info storing
    hankelSample = zeros(length(L),length(SNRnoise)); % Store the number of samples to build a Hankel matrice with rank = number of poles
    LHsAll = zeros(N,length(L),length(SNRnoise));       % Store rank LHs Hankel

    %% Main algorithm
    for ii = 1:1:length(L)
        
        % fprintf('..............\n');
        Nexp = L(ii); % Number of exponential
        
        for jj = 1:1:length(SNRnoise)

            % fprintf('|| %2d | %d ||\n',Nexp, SNRnoise(jj))
            
            %% Sum of Complex Exponentials for a fixed wn
            f = [linspace(0, wn*(Fn),Nexp)'];   % Linear Frequency
            comp = exp(1j*2*pi*f.*t);           % Poles
            s=sum(comp);  % Signal
            
            %% Noisy scenario with noise in dB SNRnoise(jj);
            y = awgn(s,SNRnoise(jj),'measured');

            %% Study - Hankel
            HsRank = zeros(N,1);    % Store the rank for a signal with jj samples
            % Hs = Hankelize(y);    % Hankelize the whole signal
            for kk=1:1:N % Computing Hankel matrix rank for a signal with jj samples 
                HsRank(kk,1) = rank(Hankelize(y(1:kk))); % Rank computing with s(1:jj)
            end
            
            LHsAll(:,ii,jj) = HsRank;

            %% Saving Hankel sample and current LHs rank
            if sum(find(HsRank == Nexp,1)) == 0 % Verify if the distance is enough to find the rank = number of poles    
                hankelSample(ii,jj) = 5*N; % When N isn't enough to find the rank = number of poles <======= here
                % disp('==> debug code:01')
            else
                hankelSample(ii,jj) = find(HsRank == Nexp,1);
                % disp('==> debug code:02')
            end
            
        end
        
        %% Set the legend
        if rmc == 1
            str = ['L = ', num2str(Nexp)]; % Save the legend
            leg{ii} = str;
        end
    
    end

    yRMC(rmc,:,:) = hankelSample;
    
end

%% Monte Carlo - elapsed time 
tEnd = toc(tStart);
fprintf('|%d rmc, %d poles, %d SNRs| \n',Nrun, length(L),length(SNRnoise))
fprintf('---> elapsed time: %.2f s <---\n',tEnd)

%% Mean for the Monte Carlo Samples
yRMC = (mean(yRMC,1));
yRMC = squeeze(yRMC);

%% For all SNR 
hFigure = figure;
myColor = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560],...
[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};

myLineStyle = {'--','-',':','-.'};
myMarker = {'o', 'square', '^', 'x'};
for ii = 1:1:length(L)
    plot(SNRnoise,yRMC(ii,:),...
    'Color', myColor{ii},...        
    'LineStyle', myLineStyle{rem(ii,4)+1},...
    'LineWidth',1.0,...
    'Marker',myMarker{rem(ii,4)+1},...
    'MarkerFaceColor', myColor{ii},...
    'MarkerSize', 3);
    hold on
end
hold off
xlabel('SNR (dB)');
ylabel('Sample size, N');
ylim([0 0.8*N])
legend(leg,'Location', 'best') % legend(leg,'Location', 'Best')
legend boxoff
grid

%% Hankel samples (rank = poles) in a noisy scenario - Export Result Figure
savefigPath = 'C:\Users\lukin\Documents\GitHub\Hankel-Matrices\matlab\fig';
outNamePDF = ['SNR_proximity-MonteCarlo-with-',num2str(Nrun),'runs-'];
datetimestamp = datestr(now, 'yyyy-mm-dd');
outfilename = strjoin({[savefigPath '\' outNamePDF '-' datetimestamp]});
% savefig_tight(h,outfilename);

%% Save the Monte Carlo Run 
save([outNamePDF,datetimestamp],'N','L','SNRnoise','yRMC','leg');