%% Plot the results of SNR_proximity.m
%
%
%HISTORY:
% 2021-05-15: Lucas Abdalah
%

%% Preparing Matlab ambient
clearvars; clc;
%% Hankel samples (rank = poles) in a noisy scenario - Import Result .mat
Nrun = 10;
savefigPath = 'C:\Users\lukin\Documents\GitHub\Hankel-Matrices\matlab\fig';
outNamePDF = ['SNR_proximity-MonteCarlo-with-',num2str(Nrun),'runs-'];
datetimestamp = datestr(now, 'yyyy-mm-dd');

%% Load the Monte Carlo Run 
load([outNamePDF,datetimestamp]);

%% For all SNR 
h = figure;
myColor = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560],...
[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};

myLineStyle = {'--','-',':','-.'};
myMarker = {'o', 'square', '^', 'x'};
for ii = 1:1:length(L)
    semilogx(SNRnoise,yRMC(ii,:),...
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

outfilename = strjoin({[savefigPath '\' outNamePDF '-' datetimestamp]});
% savefig_tight(h,outfilename);
