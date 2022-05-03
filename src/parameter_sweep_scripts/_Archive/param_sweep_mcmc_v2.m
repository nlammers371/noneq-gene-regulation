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

% specify sampling hyperparameters
log_sigma = 1.05;
model_spec_vec = [5]; % 5 corresponds to Hopfield/Ninio architecture
n_mcmc = 1; % number of independent runs to undertake
optimality_metric = 'fidelity';%'fidelity'; % 'fidelity' or 'logL divergence'
op_index = 1;
temperature = 30; % dictates acceptance rate
if strcmpi(optimality_metric,'logL divergence')
    op_index = 2;
    temperature = 10; % dictates acceptance rate
elseif strcmpi(optimality_metric,'critical rate')
    op_index = 3;
    temperature = 10; % dictates acceptance rate
end
% constrained to equilibrium?
equilibrium_flag = 1;
% sampling params and hyperparams
prior_weights = ones(1,8)/8; % controls initial size of each edg
rate_bounds = [1e-4, 1e4]; % constrain transition rate magnitude
max_iters = 1e5;
delta_mu_in = 2; % ratio of high-to-low concentrations (assumed fixed)
% initialize beta distribution to sample rate value multipliers
proposal_dist = makedist('Beta',100,100);
% iterate through specified model specs
simulation_results = struct;
for mdi = 1:numel(model_spec_vec)
    % determine specific edges
    spec_edges = m_space(model_spec_vec(mdi),:);
    % initialize arrays
    rate_array = NaN(max_iters,8,n_mcmc);
    % metric array 
    metric_array = NaN(max_iters,4,n_mcmc);    
    % record relevant hyperparameters
    simulation_results(mdi).log_simga = log_sigma;
    simulation_results(mdi).model_spec = spec_edges;    
    simulation_results(mdi).optimality_metric = optimality_metric;
    % iterate through specified number of independent MCMC runs
    for mci = 1:n_mcmc
        % initialize rates 
        rate_vec_init = rand(1,8)*rate_bounds(end);   
        if equilibrium_flag
            rate_vec_init(end) = prod(rate_vec_init(1:4)) / prod(rate_vec_init(5:7));  
        end                   
        rate_array(1,:,mci) = rate_vec_init;        
        % calculate performance_metrics
        metric_array(1,:,mci) = calculatePerformanceMetrics(delta_mu_in,spec_edges,rate_vec_init(1:4),rate_vec_init(5:8));        
        % perform MCMC sampling
        for n = 2:max_iters            
            rates_orig = rate_array(n-1,:,mdi);
            % draw jump proposals for each rate
            prop_factors = 2*random(proposal_dist,1,8);
            prop_rates = rates_orig.*prop_factors;
            if equilibrium_flag
                prop_rates(end) = prod(prop_rates(1:4)) / prod(prop_rates(5:7));
            end       
            % impose size constraints
%             prop_rates(prop_rates<rate_bounds(1)) = rates_orig(prop_rates<rate_bounds(1));
%             prop_rates(prop_rates>rate_bounds(2)) = rates_orig(prop_rates>rate_bounds(2));
            % attempt MH move for all rates simultaneously
            metrics_prop = calculatePerformanceMetrics(delta_mu_in,spec_edges,...
                prop_rates(1:4),prop_rates(5:8)); 
            % initialize temp vars
            metrics_current = metric_array(n-1,:,mdi);
            rates_current = rates_orig;
            % calculate relative performance
            prop_ratio = (metrics_prop(op_index)/metrics_current(op_index))^temperature;
            if prop_ratio > rand()
                rates_current = prop_rates;                    
                metrics_current = metrics_prop;
            end   
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