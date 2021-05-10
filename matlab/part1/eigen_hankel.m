% Computes eigendecomposition of a simple 2-pole Hankel matrix and plots its condition
% number vs sample size and pole separation.
%
%
% HISTORY:
%
% 2021/04/29: - created by Vicente Zarzoso, UCA, CNRS, I3S, France.

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

% plot results
h = figure;
for i = 1:num_dw
    semilogy(N, cond_num(i, :)); hold on
end
grid; axis tight;
% set(gca, 'FontSize', 14)
xlabel('sample size, N');
ylabel('condition number');
legend(leg)

%% conditionNumber - Export Result Figure
% savefigPath = 'C:\Users\lukin\Documents\GitHub\Hankel-Matrices\matlab\fig';
% outNamePDF = 'eigen_hankel-conditionNumber';
% datetimestamp = datestr(now, 'yyyy-mm-dd');
% outfilename = strjoin({[savefigPath '\' outNamePDF '-' datetimestamp]});
% savefig_tight(h,outfilename);