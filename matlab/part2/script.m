%% Computes the rank of an all-pole signal mapped into a Hankel matrix with various vector length.
%
%
%HISTORY:
% 2021-04-12: Lucas Abdalah
%

%% Initialization
close all;
clearvars; clc;

%% Add mandatory paths - 1st run
% addpath C:\Users\lukin\Documents\MATLAB\internship\Tensor_ECGs\myfunctions\TSP\semi_synthetic\rankHs

%% Directories path to save results
% path_fig = 'C:\Users\lukin\Documents\GitHub\Hankel-Matrices\matlab\fig';

%% General Setup
Fs = 1e3;       % Sampling Frequency (Hz)
Fn = Fs/2;      % Nyquist Frequency (Hz)
T = 1/Fs;       % Sampling period
N = 1e2;        % Length of signal
t = (0:N-1)*T;  % Time vector

%% Poles Setup
Nexp = 3; % Number of exponential
wn = 1e-4; % Normalized frequency distance  \Delta \omega
f = [linspace(0, wn*(Fn),Nexp)']; % Linear Frequency
rankLimiar = 2*Nexp - 1; % Literature limiar

%% Sum of Complex Exponentials
comp =  exp(1j*2*pi*f.*t); % Poles
% addConjugate = false; 
% if addConjugate == true
%   comp = comp + exp(-1j*2*pi.*f*t); % Complex conjugate
% end
s=sum(comp);  % Signal

%% Frequency Domain
n = 2^nextpow2(N);
sf = fft(s,n);
f = Fs*(0:(n/2))/n;
om = linspace(0,1,n/2+1);
P = abs(sf/n);

%% Study - Hankel
Hs = Hankelize(s);
LHs = rank(Hs);
HsRank = zeros(N,1);
for ii=1:1:N
  HsRank(ii,1) = rank(Hankelize(s(1:ii)));  % Rank computing with s(1:ii)
end

%% Figures
h(1) = figure;
subplot(2,1,1)
plot(om,P(1:n/2+1), 'HandleVisibility','off')
title([num2str(Nexp),' poles - \Delta\omega=' num2str(wn/(Nexp-1))]);
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('|S(f)|')
xlim([0 max(om)])
grid minor

subplot(2,1,2)
plot(1:N,HsRank,...
    'HandleVisibility','off');
hold on
line([rankLimiar rankLimiar],[0 LHs+1],'Color','red','LineStyle','--'); % Limiar delimiter
hold off
xlabel('Samples in Hs matrix');
ylabel('Hankel matrix rank');
yticks([1:floor(LHs/2):LHs])
xlim([0 N/4])
ylim([0 LHs+1])
legend([num2str(rankLimiar) ' samples'], 'Location', 'Best');
grid on