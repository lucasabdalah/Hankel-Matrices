close all; clc; 
eigen_hankel
h1 = figure; 
N_screeplot = 20;
semilogy(s(1:N_screeplot));
% hold on;
% stem(s);
xticks(0:2:length(s(1:N_screeplot)));
grid minor; 
axis tight;