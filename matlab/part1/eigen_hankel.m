% Computes eigendecomposition of a simple 2-pole Hankel matrix and plots its condition
% number vs sample size and pole separation.
%
%
% HISTORY:
%
% 2021/04/29: - created by Vicente Zarzoso, UCA, CNRS, I3S, France.
% 2021/05/15: - modification by Lucas Abdalah, UFC, GTEL, Brazil. 
% Description: To Insert Graphics and verify the relation between Cond_num vs. H matrix

%
close all; clearvars; clc;

% simulation parameters
num_dw = 7;                         % number of angular separation values
num_N = 100;                         % number of sample size values
dw = logspace(-8.0, -7.0, num_dw);       % pole angular separation
N = round(logspace(1, 3, num_N));   % sample size

cond_num = zeros(num_dw, num_N);
Hrank = zeros(num_dw, num_N);
rank_cond = zeros(num_dw, 2);
% loop over parameters
for i = 1:num_dw

    dw_current = dw(i);
    str = ['\Delta\omega = ', num2str(dw_current, '%.2e')];
    leg{i} = str;      % for legend
    disp(['===>>> Computing ', str])
   
    for j = 1:num_N
        % generate Vandermonde matrix of current parameters
        N_current = N(j);
        V = [ones(N_current, 1), exp(1j*(0:N_current-1)'*dw_current)];
        H = V*V.';
        s = svd(H);    % SVD of associated Hankel matrix
        cond_num(i, j) = abs(s(1)/s(2));
        Hrank(i,j) = rank(H);
        
        if j == 1
            % disp(['===> if'])
            if Hrank(i,j) == 2;
                rank_cond(i,:) = [N_current;cond_num(i,j)];       
            end
        else
            % disp(['===> else'])
            if Hrank(i,j) ~= 2
                if Hrank(i,j-1) ~= 2
                    rank_cond(i,:) = [N_current;cond_num(i,j)];
                end
            end
        end
        

    end
    
end

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
%% Scatter the point where the rank changes
scatter(rank_cond(:,1),rank_cond(:,2),...
        'Marker', 'd', ...
        'MarkerFaceColor', [0,0,0],...
        'MarkerEdgeColor',[0,0,0],...
        'SizeData', 60);
xlabel('Sample size, N');
ylabel('Condition number, \sigma');
xlim([min(N)-1 max(N)+1])
ylim([min(cond_num(:))-1 max(cond_num(:))+1])

xticks([0:10:100])

leg{num_dw+1} = 'R = L'; 
legend(leg, 'Location', 'Northeastoutside') % legend(leg,'Location', 'Best')
legend boxoff
grid;
hold off

%% conditionNumber - Export Result Figure
savefigPath = 'C:\Users\lukin\Documents\GitHub\Hankel-Matrices\matlab\fig';
outNamePDF = 'eigen_hankel-conditionNumber3';
datetimestamp = datestr(now, 'yyyy-mm-dd');
outfilename = strjoin({[savefigPath '\' outNamePDF '-' datetimestamp]});
% savefig_tight(h,outfilename);