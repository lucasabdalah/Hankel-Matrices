% Computes eigendecomposition of a simple 2-pole Hankel matrix and plots its condition
% number vs sample size and pole separation.
%
%
% HISTORY:
%
% 2021/04/29: - created by Vicente Zarzoso, UCA, CNRS, I3S, France.

%
close all; clearvars; clc;

% simulation parameters
num_dw = 5;                         % number of angular separation values
num_N = 10;                         % number of sample size values
dw = logspace(-4, 0, num_dw);       % pole angular separation
N = round(logspace(1, 3, num_N));   % sample size

cond_num = zeros(num_dw, num_N);

% loop over parameters
for i = 1:num_dw
    dw_current = dw(i);
    str = ['\Delta\omega = ', num2str(dw_current)];
    leg{i} = str;      % for legend
    disp(['===>>> Computing ', str])
   
    for j = 1:num_N
        
        % generate Vandermonde matrix of current parameters
        N_current = N(j);
        V = [ones(N_current, 1), exp(1j*(0:N_current-1)'*dw_current)];
        
        s = svd(V*V.');    % SVD of associated Hankel matrix
        
        cond_num(i, j) = abs(s(1)/s(2));
    end
    
end

% plot results - vz 
% h = figure;
% for i = 1:num_dw
%     semilogy(N, cond_num(i, :)); hold on
% end
% grid; axis tight;
% % set(gca, 'FontSize', 14)
% xlabel('sample size, N');
% ylabel('condition number');
% legend(leg)

%% Plot results - la
h = figure;
myColor = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};
myLineStyle = {'--','-',':','-.'};
myMarker = {'o', 'square', '^', 'x'};

for ii = 1:num_dw
    semilogy(N, cond_num(ii, :),...
        'Color', myColor{ii},...        
        'LineStyle', myLineStyle{rem(ii,4)+1},...
        'LineWidth',1.0,...
        'Marker',myMarker{rem(ii,4)+1},...
        'MarkerFaceColor', myColor{ii},...
        'MarkerSize', 3);
    hold on;
end
hold off
xlabel('N (Sample size)');
ylabel('\sigma (Condition number)');
legend(leg) % legend(leg,'Location', 'Best')
grid; axis tight;
ylim([0.5 1e8])


%% conditionNumber - Export Result Figure
savefigPath = 'C:\Users\lukin\Documents\GitHub\Hankel-Matrices\matlab\fig';
outNamePDF = 'eigen_hankel-conditionNumber2';
datetimestamp = datestr(now, 'yyyy-mm-dd');
outfilename = strjoin({[savefigPath '\' outNamePDF '-' datetimestamp]});
savefig_tight(h,outfilename);