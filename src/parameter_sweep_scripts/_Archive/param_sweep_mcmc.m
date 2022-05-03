% script to use MCMC sampling to explore performance landscape of
% proofreading networks
clear
close all
% specify paths
addpath('../utilities')
outPath = '../../out/parameter_search/';
mkdir(outPath)

% specify permitted model architectures (determined by location of specific
% edges)
m_space_full = [repelem(1:4,4)' repmat((1:4)',4,1)]; 
m_space_ne = m_space_full(m_space_full(:,1)~=m_space_full(:,2),:); % opposing specific edges may not connect same two states
m_space = m_space_ne(~ismember(m_space_ne(:,1)-m_space_ne(:,2),[1,-3]),:); % cannot have both edges leaving a state be specific
% calculate expected cycle flux from expenditure of 1 ATP per cycle
k_b = 1.381e-23;
T = 273;
j_per_atp = 1e-19;
% cycle_constant = j_per_atp / (k_b*T); % why is this so big?
cycle_constant = 1000;
% specify sampling hyperparameters
log_sigma = 1.05;
model_spec_vec = [5]; % 2 corresponds to Hopfield/Ninio architecture
n_mcmc = 1; % number of independent runs to undertake
optimality_metric = 'critical rate';%'fidelity'; % 'fidelity' or 'logL divergence'
op_index = 1;
temperature = 30; % dictates acceptance rate
if strcmpi(optimality_metric,'logL divergence')
    op_index = 2;
    temperature = 10; % dictates acceptance rate
elseif strcmpi(optimality_metric,'critical rate')
    op_index = 3;
    temperature = 10; % dictates acceptance rate
end

prior_weights = ones(1,8)/8; % controls initial size of each edg
rate_bounds = [1e-8, 1e3]; % constrain transition rate magnitude
max_iters = 1e5;
delta_mu_in = 2; % ratio of high-to-low concentrations (assumed fixed)
sampling_vec = 1:8; % designate edges to sample (holds r4 fixed at 1)
excluded_vec = ~ismember(1:8,sampling_vec);
r_options = sampling_vec(sampling_vec>4);
k_options = sampling_vec(sampling_vec<=4);
% iterate through specified model specs
simulation_results = struct;
for mdi = 1:numel(model_spec_vec)
    spec_edges = m_space(model_spec_vec(mdi),:);
    % initialize arrays
    rate_array = NaN(max_iters,8,n_mcmc);
    % metric array 
    metric_array = NaN(max_iters,4,n_mcmc);    
    % record relevant hyperparameters
    simulation_results(mdi).log_simga = log_sigma;
    simulation_results(mdi).model_spec = spec_edges;
    simulation_results(mdi).temperature = temperature;
    simulation_results(mdi).optimality_metric = optimality_metric;
    % iterate through specified number of independent MCMC runs
    for mci = 1:n_mcmc
        % initialize rates
        pass = false;
        while ~pass
            rate_vec_init = rand(1,8)*rate_bounds(end);                
            rate_vec_init(excluded_vec) = 1;
            ratio_init = prod(rate_vec_init(1:4))/prod(rate_vec_init(5:8));
            flux_factor = cycle_constant/ratio_init;
            candidates = find(rate_vec_init(5:8)/flux_factor >= rate_bounds(1) & ~ismember(5:8,excluded_vec));
            if ~isempty(candidates)
                rate_vec_init(4+candidates(1)) = rate_vec_init(4+candidates(1)) / flux_factor;
                pass = true;
            end
        end               
        rate_array(1,:,mci) = rate_vec_init;        
        % calculate performance_metrics
        metric_array(1,:,mci) = calculatePerformanceMetrics(delta_mu_in,spec_edges,rate_vec_init(1:4),rate_vec_init(5:8));        
        % perform MCMC sampling
        for n = 2:max_iters            
            rates_orig = rate_array(n-1,:,mdi);
            % draw jump proposals for each rate
            prop_factors = lognrnd(log(1),log(log_sigma),1,4);
            r_pairs = randsample(r_options,numel(k_options),true);
            prop_k_vec = rates_orig(1:4).*prop_factors;
            prop_k_vec(prop_k_vec<rate_bounds(1)) = rate_bounds(1);
            prop_k_vec(prop_k_vec>rate_bounds(2)) = rate_bounds(2);
            % attempt MH move for each param independly
            rates_current = rates_orig;            
            metrics_current = metric_array(n-1,:,mdi);
            iter = 1;
            for i = k_options
                rates_prop = rates_current;
                rates_prop(i) = prop_k_vec(i);
                r_new = rates_current(r_pairs(iter))*prop_factors(i);
                if r_new < rate_bounds(1) || r_new > rate_bounds(2)
                    continue
                end
                rates_prop(r_pairs(iter)) = r_new;                             
                metrics_prop = calculatePerformanceMetrics(delta_mu_in,spec_edges,rates_prop(1:4),rates_prop(5:8));
                % make decision
%                 prop_ratio = exp(-(metrics_current(op_index)-metrics_prop(op_index))/temperature);
                prop_ratio = (metrics_prop(op_index)/metrics_current(op_index))^temperature;
                if prop_ratio > rand()
                    rates_current = rates_prop;                    
                    metrics_current = metrics_prop;
                end                              
            end    
            % record
            rate_array(n,:,mdi) = rates_current;
            metric_array(n,:,mdi) = metrics_current;  
        end
    end
end

if op_index == 1
    plot(exp(metric_array(:,op_index)))
else
    plot(metric_array(:,op_index))
end

%%
cmb = brewermap(128,'Spectral');
figure(1);
colormap(cmb)
scatter(exp(metric_array(:,1))/2,metric_array(:,3),20,1:size(metric_array,1),'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',0)
xlabel('cooperativity')
ylabel('decision rate')
grid on
set(gca,'Fontsize',12)