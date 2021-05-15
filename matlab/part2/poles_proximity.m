%% Computes the rank of an all-pole signal mapped into a Hankel matrix with various vector length.
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

%% Poles setup
itPoles = 2:1:8;    % Number of poles 
spacings = 100;     % Space between poles
wAll = logspace(-4,-1,spacings);

%% Info storing
hankelSample = zeros(length(itPoles),length(wAll)); % Store the number of samples to build a Hankel matrice with rank = number of poles
LHsAll = zeros(length(itPoles),length(wAll));       % Store rank LHs Hankel
fail2Compute = 1; % Store the number of fails to compute rank = Nexp

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

        %% Study - Hankel
        HsRank = zeros(N,1);    % Store the rank for a signal with jj samples
        Hs = Hankelize(s);      % Hankelize the whole signal
        for kk=1:1:N % Computing Hankel matrix rank for a signal with jj samples 
            HsRank(kk,1) = rank(Hankelize(s(1:kk)));  % Rank computing with s(1:jj)
        end
       
        %% Saving Hankel sample and current LHs rank
        if sum(find(HsRank == Nexp,1)) == 0 % Verify if the distance is enough to find the rank = number of poles
            fail2Compute=fail2Compute+1;  
            % fprintf(['we are in! => Fail no. ', num2str(fail2Compute),'\n']);
            hankelSample(ii,jj) = 150; % Replace when N isn't enough to find the rank = number of poles <======= here
        else
            hankelSample(ii,jj) = find(HsRank == Nexp,1);
            % fprintf(['Rank =>', num2str(Nexp),'\t','Sample =>', num2str(hankelSample(ii,jj)),'\n']) 
        end

    end

    str = ['L = ', num2str(Nexp)]; % Save the legend
    leg{ii} = str;

end

%% Plot results
h = figure;
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
    % if ii > 5
    % line([min(wAll) max(wAll)], (2*itPoles(ii)-1)*[1 1],...
    %     'Color', [0.1 0.1 0.1],...
    %     'LineStyle','--',...
    %     'HandleVisibility','off');
    % end
end
hold off
xlabel('\Delta\omega (\times \pi radians/samples)');
ylabel('N (Sample size)');
ylim([0 80])
legend(leg) % legend(leg,'Location', 'Best')
legend boxoff
grid

%% Hankel samples (rank = poles) - Export Result Figure
savefigPath = 'C:\Users\lukin\Documents\GitHub\Hankel-Matrices\matlab\fig';
outNamePDF = 'poles_proximity-PoleDistance_vs_Samples';
datetimestamp = datestr(now, 'yyyy-mm-dd');
outfilename = strjoin({[savefigPath '\' outNamePDF '-' datetimestamp]});
savefig_tight(h,outfilename);