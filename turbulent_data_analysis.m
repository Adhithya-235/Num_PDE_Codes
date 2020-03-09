%========================================================================%

% Playing around with turbulent data

%========================================================================%

clear
close all
clc

%% LOAD DATA

hwd1 = load('Hotwire_data1.mat');
hwd2 = load('Hotwire-data2.mat');
hwd3 = load('Hotwire-data3.mat');

%% PLOTTING HISTOGRAM

figure(1)
h = histogram(hwd3.u, 101);
xlabel('u');
ylabel('count');
grid on
axis square
axis tight

%% PLOTTING PDF

 p = histcounts(hwd3.u, 101, 'normalization', 'pdf');
 u_vals = h.BinEdges + (h.BinWidth/2);
 u_vals = u_vals(1:end-1);
 figure(2)
 plot(u_vals, p, 'k-', 'linewidth', 1.2)
 xlabel('u')
 ylabel('probability density')
 grid on
 axis square
 axis tight
 
 %% CALCULATE MOMENTS
 
 % --> FIRST MOMENT
 
 u_mean = trapz(u_vals, u_vals.*p);
 u_mean_dir = mean(hwd3.u);
 
 % --> SECOND CENTRAL MOMENT
 % --> The central moments are approximately the same as the moments here 
 %     because the mean is nearly zero. 
 
 u_var = trapz(u_vals, p.*(u_vals.^2));
 u_var_dir = var(hwd3.u);
 
 % --> THIRD CENTRAL MOMENT, SKEWNESS
 
 u_skew = (trapz(u_vals, p.*(u_vals.^3)))/(u_var^(3/2));
 u_skew_dir = skewness(hwd3.u);
 
 % --> FOURTH CENTRAL MOMENT, KURTOSIS
 
 u_kurt = (trapz(u_vals, p.*(u_vals.^4)))/(u_var^2);
 u_kurt_dir = kurtosis(hwd3.u);
 
 %% AUTOCORRELATION FUNCTION
 
 u_t = hwd3.u;
 Nlags = 100;
 
 % --> CREATE LAG MATRIX
 
 U_SHIFT = zeros(length(u_t),Nlags);
 U_SHIFT(:,1) = u_t;
 for i = 2:Nlags
     
     U_SHIFT(:,i) = circshift(U_SHIFT(:,i-1),-1);
     
 end
 
 % --> CALCULATE AUTOCORRELATION FUNCTION
 
 rho_T = zeros(Nlags,1);
 
 for i = 1:Nlags
     
     rho_T(i) = (mean(u_t.*U_SHIFT(:,i)))/(std(u_t)*std(U_SHIFT(:,i)));
     
 end
 
 % --> PLOT AUTOCORRELATION FUNCTION
 
 figure(3)
 plot(rho_T,'k-','linewidth',1.2)
 
 

