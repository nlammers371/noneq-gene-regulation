clear
close all

FigPath = '../../fig/validation_v2/';
mkdir(FigPath);

% specify core simulation parameters
t_sim = 1e5; % # time points
tau_c = 50; % time scale of extrinsic fluctuations (in seconds)
dT = .1; % resolution of simulation (must be a lot smaller than rates)
t_vec = linspace(0,t_sim*dT,t_sim);
mu_c = 1; % average activator concentration
sigma_c = .5; % variation in concentration
n_sim = 100;

%%%%%%%%%%%%%%%%%
% Generate correlated noise vector computationally
%%%%%%%%%%%%%%%%%

% kernel to introduce temporal correlation
corr_kernel = exp(-(cumsum(ones(1,5*tau_c/dT)*dT)-dT)/tau_c);


% draw uncorrelated random gaussian vec
gauss_vec = normrnd(mu_c,sigma_c,1,t_sim+4*numel(corr_kernel));
% generate correlated noise vec 
corr_gauss_vec = conv(gauss_vec,corr_kernel,'full');
norm_vec = conv(ones(size(gauss_vec)),corr_kernel,'full');
corr_gauss_vec = corr_gauss_vec(2*numel(corr_kernel):end-3*numel(corr_kernel))./...
    norm_vec(2*numel(corr_kernel):end-3*numel(corr_kernel));

% scale up noise to desired level
corr_gauss_diffs = corr_gauss_vec-mean(corr_gauss_vec);
signs = sign(corr_gauss_diffs);
corr_gauss_diffs_new = sqrt((corr_gauss_diffs.^2).*sigma_c^2/var(corr_gauss_diffs)).*signs;
corr_gauss_vec_new = corr_gauss_diffs_new + mu_c;
corr_gauss_vec_new(corr_gauss_vec_new<0) = 0;


%%%%%%%%%%%%%%%%%%%
% Check that vector has expected autocorrelation 
a = autocorr(corr_gauss_vec_new,numel(corr_kernel)-1);

corr_fig = figure;
hold on
plot(t_vec(1:numel(corr_kernel)),corr_kernel);
plot(t_vec(1:numel(corr_kernel)),a,'--');

legend('simulation','theoretical expectation')
xlabel('lag (seconds)')
ylabel('autocorrelation')
set(gca,'FontSize',14)
grid on
saveas(corr_fig,[FigPath 'correlated_noise_validation.png'])




%%%%%%%%%%%%%%%%%%%%
% Examine variance of sum of vecs samped at different frequencies
gauss_array = normrnd(0,sigma_c,t_sim+4*numel(corr_kernel),n_sim);
% generate correlated noise vec 
corr_gauss_array_raw = conv2(gauss_array,corr_kernel','full');
% normalize
norm_vec = conv(ones(size(gauss_array(:,1))),corr_kernel,'full');
norm_mat = repmat(norm_vec,1,n_sim);
corr_gauss_array_raw = corr_gauss_array_raw(2*numel(corr_kernel):end-3*numel(corr_kernel),:)./...
    norm_mat(2*numel(corr_kernel):end-3*numel(corr_kernel),:);

% scale up noise to desired level    
signs = sign(corr_gauss_array_raw);
corr_gauss_diffs_new = sqrt((corr_gauss_array_raw.^2).*sigma_c^2./var(corr_gauss_array_raw)).*signs;
corr_gauss_array = corr_gauss_diffs_new + mu_c;
corr_gauss_array(corr_gauss_array<0) = 0;

%%

n_vec = round(logspace(0,3,100));
snr_vec = NaN(size(n_vec));
for n = 1:numel(n_vec)
    cr_ds = corr_gauss_array(1:n_vec(n):t_sim,:);
    snr_vec(n) = var(sum(cr_ds))/size(cr_ds,1);
end
term_check = NaN(size(snr_vec));
pd_var = repelem(sigma_c^2,numel(snr_vec));
for n = 1:numel(n_vec)
    term_check(n) = sum(exp(-dT*(n_vec(n):n_vec(n):t_sim) / tau_c));
    pd_var(n) = pd_var(n)+ 2*sigma_c^2 *term_check(n);
end

%%
t_c = 1;
t2 = linspace(0,10,1e5);
dT2 = median(diff(t2));
test = exp(-t2/t_c);
sum_vec = NaN(size(snr_vec));
for n = 1:numel(n_vec)
    sum_vec(n) = sum(test(n_vec(n):n_vec(n):end));
end
    



