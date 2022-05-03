% script to validate analytic expression for variance of general 4 state
% Markov chain using stochastic simulations
clear
close all
addpath('../utilities')
% define parameters
c1 = 10;
c0 = 1;
r1 = 1;
r2 = 1;
r3 = 1;
k1 = 1;
k3 = 1;
k4 = 1;
r4 = 1;
k2 = (r1*r2*r3*r4) / (k1*k4*k3);

rate_vec = [k1,k2,k3,k4,r1,r2,r3,r4];
rate_vec_high = rate_vec .* [1 c1 1 1 1 1 1  c1];
rate_vec_low = rate_vec .* [1 c0 1 1 1 1 1  c0];

R = @(k,r) [-k(1)-r(4)   r(1)             0               k(4); 
             k(1)        -r(1)-k(2)      r(2)              0
              0           k(2)         -r(2)-k(3)         r(3)
             r(4)          0               k(3)        -r(3)-k(4) ];                
R_high = R(rate_vec_high(1:4),rate_vec_high(5:8));
R_low = R(rate_vec_low(1:4),rate_vec_low(5:8));
% calculate mean and variance for high and low c mRNA production
rHigh = fourStateProduction(rate_vec,c1,1);
sigmaHigh = sqrt(fourStateVariance(rate_vec,c1,1));
rLow = fourStateProduction(rate_vec,c0,1);
sigmaLow = sqrt(fourStateVariance(rate_vec,c0,1));



% simulation parameters
T = 500;
n_sim = 1000;
t_grid = 1:T;
% simulate trajectories and track movement of log likelihood ratio over
% time
production_array_high = NaN(T,n_sim);
production_array_low = NaN(T,n_sim);

state_array_high = NaN(T,n_sim);
state_array_low = NaN(T,n_sim);

logL_array_high = NaN(T,n_sim);
logL_array_low = NaN(T,n_sim);

production = [0 0 1 0];
state_options = 1:4;
% iterate
for n = 1:n_sim
    if mod(n,10) == 0
        disp(n)
    end
    state_array_high(1,n) = randsample(state_options,1);
    state_array_low(1,n) = randsample(state_options,1);
    
    production_array_high(1,n) = 0;
    production_array_low(1,n) = 0;
    
    logL_array_high(1,n) = 0;
    logL_array_low(1,n) = 0;
    
    total_time = 0;    
    for t = t_grid(2:end)
        t_fin = 1;
        t_inc_low = 0;
        pd_low = 0;
        t_inc_high = 0;
        pd_high = 0;
        cs_low = state_array_low(t-1,n);
        cs_high = state_array_high(t-1,n);
        % simulate low chain first
        while t_inc_low < t_fin
            % draw jump time
            dt = exprnd(1/-R_low(cs_low,cs_low),1);
            t_inc_low = t_inc_low + dt;
            pd_low = pd_low + min(dt,t_fin-t_inc_low+dt)*production(cs_low);      
            if t_inc_low < t_fin              
                options = state_options(state_options~=cs_low);
                cs_low = randsample(options,1,true,R_low(options,cs_low));                                                                                              
            end            
        end
        production_array_low(t,n) = pd_low;
        pd_low_tot = sum(production_array_low(1:t,n))/t;
        state_array_low(t,n) = cs_low;
        logL_array_low(t,n) = log(sigmaHigh/sigmaLow) - .5*t*(((rHigh-pd_low_tot)/sigmaHigh)^2 - ((rLow-pd_low_tot)/sigmaLow)^2);
        % high chain
        while t_inc_high < t_fin
            % draw jump time
            dt = exprnd(1/-R_high(cs_high,cs_high),1);
            t_inc_high = t_inc_high + dt;
            pd_high = pd_high + min(dt,t_fin-t_inc_high+dt)*production(cs_high);      
            if t_inc_high < t_fin                                                                                                                       
                options = state_options(state_options~=cs_high);
                cs_high = randsample(options,1,true,R_high(options,cs_high)); 
            end            
        end
        production_array_high(t,n) = pd_high;
        pd_high_tot = sum(production_array_high(1:t,n))/t;
        state_array_high(t,n) = cs_high;
        logL_array_high(t,n) = log(sigmaHigh/sigmaLow) - .5*t*(((rHigh-pd_high_tot)/sigmaHigh)^2 - ((rLow-pd_high_tot)/sigmaLow)^2);
    end       
end

cum_pd_high = cumsum(production_array_high);
cum_pd_low = cumsum(production_array_low);

snr_sim_high = nanmean(cum_pd_high,2) ./ nanstd(cum_pd_high,[],2);
snr_sim_low = nanmean(cum_pd_low,2) ./ nanstd(cum_pd_low,[],2);

sigma_vec_high = sigmaHigh*sqrt(t_grid);
sigma_vec_low = sigmaLow*sqrt(t_grid);

r_vec_high = rHigh*t_grid;
r_vec_low = rLow*t_grid;

pSuccess_vec = probSuccess(rLow,rHigh,sigmaLow,sigmaHigh,t_grid);

r_critical_vec = (sigma_vec_high.*r_vec_low+ sigma_vec_low.*r_vec_high) ./ (sigma_vec_high+sigma_vec_low);

% flag cases in which simulated trajectories cross critical r
pd_high_flags = cum_pd_high < repmat(r_critical_vec',1,n_sim);
pd_low_flags = cum_pd_low >= repmat(r_critical_vec',1,n_sim);

pSuccessSim_vec = 1-nanmean([pd_high_flags pd_low_flags],2);

% calcualte simulated and exptected logL divergence
logL_vec_sim = nanmean([-logL_array_low logL_array_high],2);
logL_vec_analytic = .25*(((rLow-rHigh)/sigmaLow)^2+((rLow-rHigh)/sigmaHigh)^2)*t_grid;

%%
figPath = '../../fig/information_figs/';
mkdir(figPath)
% make figures
cm = jet(128);
snr_fig = figure;
hold on
p1 = plot(snr_sim_high,'Color',cm(120,:),'LineWidth',1.5);
plot(r_vec_high./sigma_vec_high,'--','Color',cm(120,:),'LineWidth',1.5);

p2 = plot(snr_sim_low,'Color',cm(35,:),'LineWidth',1.5);
plot(r_vec_low./sigma_vec_low,'--','Color',cm(35,:),'LineWidth',1.5);

grid on
xlabel('time (au)')
ylabel('SNR')
legend([p1 p2],'high region','low region')
saveas(snr_fig,[figPath 'snr_validation.png'])


pSuccess_fig = figure;
hold on
plot(pSuccessSim_vec,'LineWidth',1.5);
plot(pSuccess_vec,'--','LineWidth',1.5);

grid on
xlabel('time (au)')
ylabel('success probability')
legend('simulation','analytic solution')
saveas(pSuccess_fig,[figPath 'pSuccess_validation.png'])

%%
mv_path = [figPath 'pd_frames/'];
mkdir(mv_path);
cmb = brewermap(9,'Set2');
close all
hbins = linspace(0,1);

for i = 5:15:T
    pd_fig = figure('Visible','off');
    hold on
    plot(cum_pd_low(1:i,:),'Color',[cmb(3,:) .1])
    plot(cum_pd_high(1:i,:),'Color',[cmb(2,:) .1])    
%     grid on
    xlabel('time')
    ylabel('accumulated mRNA')    
    xlim([0 T])
    ylim([0 250])
    legend('high region','low region','Location','northwest')
    set(gca,'Fontsize',14)
    saveas(pd_fig,[mv_path 'cm_pd_frame' sprintf('%03d',i) '.tif'])
    
    pSuccess_fig = figure('Visible','off');   
    plot(pSuccessSim_vec(1:i),'LineWidth',1.5);
%     grid on
    xlabel('time')
    ylabel('success probability')     
    set(gca,'Fontsize',14)
    xlim([0 T])
    ylim([.5 1])
    set(gca,'Fontsize',14)
    saveas(pSuccess_fig,[mv_path 'pSuccess_frame' sprintf('%03d',i) '.tif'])
    
    hist_fig = figure('Visible','off');   
    hold on
    h0 = histogram(cum_pd_low(i,:)/i,hbins,'EdgeAlpha',0,'FaceColor',cmb(3,:),'Normalization','probability');
    h1 = histogram(cum_pd_high(i,:)/i,hbins,'EdgeAlpha',0,'FaceColor',cmb(2,:),'Normalization','probability');
    xlabel('accumulated mRNA (normalized)')
    ylabel('share')
    xlim([0,1])
    ylim([0,.2])
    set(gca,'Fontsize',14)
    saveas(hist_fig,[mv_path 'hist_frame' sprintf('%03d',i) '.tif'])
end