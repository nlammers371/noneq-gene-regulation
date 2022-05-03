% script to validate analytic expression for rate of log likelihood
% divergence for (a) system with full "microscopic" transitions and (b)
% system where only accumulate mRNA is available to drive a decision

clear
close all
addpath(genpath('../utilities'))

% make figure path
FigurePath = '../../fig/validation/';
mkdir(FigurePath)
% make path to save data
DataPath = '../../out/validation/';
mkdir(DataPath)
% %%%%%%%%%%%%%%%%%%%%%%% define core parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

T = 1e4;
n_conditions = 1;%10;
n_reps = 100; % number of replicates per parameter set
n_sim = n_conditions^2;
c0 = 1;
c1_vec = c0*logspace(0,1,n_conditions+2); % range of "high" concentrations
c1_vec = 1.1;%c1_vec(2:end-1);

c_true = c1_vec; % For now make the simplifying assumption that lower concentration is the "true" one
tau_cycle = 1; % enforce consistent cycle kinetics
koff_vec = 1e1;%1./linspace(0,tau_cycle,n_conditions+2); % unbinding range
koff_vec = 1e-4;%koff_vec(2:end-1);
kon_vec = 1e-4;%1./(tau_cycle-1./koff_vec);

% set target error rate to 1%
minError = .01; 
K_target = log((1-minError)/minError);

% initialize parallel workers
myCluster = parcluster('local');
NumWorkers = myCluster.NumWorkers;
% p = gcp('nocreate');
% if isempty(p)
%   parpool(ceil(NumWorkers/3));
% end

% %%%%%%%%%%%%%%% generate lists of parameters for simulation %%%%%%%%%%%%%
sim_results = struct;

iter = 1;

for tc = 1:n_conditions
    c1 = c1_vec(tc);
    for k1 = 1:n_conditions
        kon = kon_vec(k1);       
        koff = koff_vec(k1);
        % update lists
        sim_results(iter).c1 = c1;            
        sim_results(iter).koff = koff;
        sim_results(iter).kon = kon;        
        iter = iter + 1;   
    end
end

%% Calculate theoretical expectations for each condition
for n = 1:n_sim
    
    % extract simulation parametes
    kon = sim_results(n).kon;
    koff = sim_results(n).koff;
    c1 = sim_results(n).c1;
    
    % calculate analytic predictions
    [sim_results(n).r0, sim_results(n).r1, sim_results(n).rt, sim_results(n).v0, ...
      sim_results(n).v1, sim_results(n).vt, sim_results(n).V_mRNA, sim_results(n).D_mRNA,...
      sim_results(n).V_micro, sim_results(n).D_micro, sim_results(n).tau_micro, sim_results(n).tau_mRNA]...
                          = calculateDriftStats2State(kon,koff,c0,c1,c_true,K_target);
end                      
                  
% calculate simulation durations and resampling resolutions based upon 
% estimated decision times
tau_array = vertcat([sim_results.tau_mRNA],[sim_results.tau_micro]);
resamp_increment = 10.^round(log10(min(tau_array,[],1)/10));
tau_max = ceil(5*max(tau_array,[],1));

for n = 1:n_sim
  sim_results(n).T = tau_max(n);
  sim_results(n).dT = resamp_increment(n);
end

%% %%%%%%%%%%%%% Conduct simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_iter = 1;
state_options = 0:1;

% simulate trajectories and track movement of log likelihood ratio over
% time
% rng(234)
% iterate
for n = 1:n_sim    
  
    % extract simulation parametes
    kon = sim_results(n).kon;
    koff = sim_results(n).koff;
    c1 = sim_results(n).c1;
 
    % for convenience
    tau_vec = [1/(c_true*kon) 1/koff];
    
    % generate grid for resampling
    T = sim_results(n).T;
    dT = sim_results(n).dT;
    time_grid = 0:dT:T;
    sim_results(n).time_grid = time_grid;
    
    % initialize vectors to store results                      
    sim_results(n).total_mRNA = NaN(length(time_grid),n_reps);    
    sim_results(n).logL_micro = NaN(length(time_grid),n_reps);
    sim_results(n).logL_micro2 = NaN(length(time_grid),n_reps);
    sim_results(n).logL_mRNA = NaN(length(time_grid),n_reps);
    
    % extract analytic quantities
    r0 = sim_results(n).r0;
    r1 = sim_results(n).r1;
    rt = sim_results(n).rt;
    v0 = sim_results(n).v0;
    v1 = sim_results(n).v1;
    
    for r = 1:n_reps
    
        % initialize temporary structures        
        jump_vec = [0];
        total_mRNA = [0];
        logL_micro = [0];
        logL_micro2 = [0];
        logL_mRNA = [0];
        total_time = 0;
        % initialize binding state using SS vector    
        state_vec = [randsample(state_options,1,true)];%,[1-rt, rt])];

        while total_time < T
          
            current_state = state_vec(end);

            % draw jump time and toggle state
            dt = exprnd(tau_vec(current_state+1),1);
            next_state = state_options(state_options~=current_state);

            total_time = total_time + dt;

            jump_vec(end+1) = dt;
            state_vec(end+1) = next_state; 
            total_mRNA(end+1) = total_mRNA(end) + dt*current_state;
            
            % update micro logL 
            if next_state == 1 
              k0 = (kon*c0*koff)/(kon*c0+koff);
              k1 = (kon*c1*koff)/(kon*c1+koff);
              dt2 = sum(jump_vec(end-1:end));
              logL_micro(end+1) = logL_micro(end) + log(c1/c0)-kon*dt*(c1-c0);
              logL_micro2(end+1) = logL_micro2(end) + log(c1*(kon*c0+koff)/(c0*(kon*c1+koff)))-dt2*(kon*koff*c1/(c1*kon+koff)-kon*koff*c0/(c0*kon+koff));
%               logL_mRNA(end+1) = logL_mRNA(end) + log(c1*(kon*c0+koff)/(c0*(kon*c1+koff))) + ...
%                dt2.*(-0.5*(r1-r_mean)^2*k1^2 + 0.5*(r0-r_mean)^2*k0^2);
            else
              logL_micro(end+1) = logL_micro(end);
              logL_micro2(end+1) = logL_micro2(end);
%               logL_mRNA(end+1) = logL_mRNA(end);
            end

            % update the mRNA logL 
            r_curr = dt*current_state;
            r_mean = total_mRNA(end)/total_time;
            logL_mRNA(end+1) = log(sqrt(v0/v1)) + 0.5*total_time.*(-(r1-r_mean)^2/v1 + (r0-r_mean)^2/v0); 
             
            
        end    
        
        % downsample and store quantities of interest
        cumulative_time = cumsum(jump_vec);
        sim_results(n).total_mRNA(:,r) = interp1(cumulative_time,total_mRNA,time_grid);
        sim_results(n).logL_mRNA(:,r) = interp1(cumulative_time,logL_mRNA,time_grid);
        sim_results(n).logL_micro(:,r) = interp1(cumulative_time,logL_micro,time_grid);       
        sim_results(n).logL_micro2(:,r) = interp1(cumulative_time,logL_micro2,time_grid);       
    end    
end


%% Check simulations against theoretical predictions

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% use linear fits to estimate simulated logL drift rates
%%%%%%%%%%%%%%%%%%%%%%%%

predicted_drift_vec_micro = [sim_results.V_micro];
sim_drift_vec_micro = NaN(1,n_sim);
predicted_drift_vec_mRNA = [sim_results.V_mRNA];
sim_drift_vec_mRNA = NaN(1,n_sim);
parfor n = 1:n_sim
  time_array = repmat(sim_results(n).time_grid',1,n_reps);
  logL_mRNA = sim_results(n).logL_mRNA;
  logL_micro = sim_results(n).logL_micro;
  
  sim_results(n).logL_micro_mean = nanmean(sim_results(n).logL_micro,2);
  sim_results(n).logL_mRNA_mean = nanmean(sim_results(n).logL_mRNA,2);
  % take second half only to avoid initial transients
  start_i = ceil(length(logL_mRNA)/2);
  
  % perform linear fits
  micro_fit = fitlm(reshape(time_array(start_i:end,:),[],1),reshape(logL_micro(start_i:end,:),[],1));
  mRNA_fit = fitlm(reshape(time_array(start_i:end,:),[],1),reshape(logL_mRNA(start_i:end,:),[],1));
  
  % store
  sim_drift_vec_micro(n) = micro_fit.Coefficients.Estimate(2);
  sim_drift_vec_mRNA(n) = mRNA_fit.Coefficients.Estimate(2);
  
  sim_results(n).sim_logL_drift_micro = sim_drift_vec_micro(n);
  sim_results(n).sim_logL_drift_mRNA = sim_drift_vec_mRNA(n);
end
  
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% find K first passage times for simualted data
%%%%%%%%%%%%%%%%%%%%%%%%

predicted_tau_vec_micro = [sim_results.tau_micro];
sim_tau_vec_micro = NaN(1,n_sim);
predicted_tau_vec_mRNA = [sim_results.tau_mRNA];
sim_tau_vec_mRNA = NaN(1,n_sim);

for n = 1:n_sim
  time_grid = sim_results(n).time_grid;
  logL_mRNA = sim_results(n).logL_mRNA;
  logL_micro = sim_results(n).logL_micro;
    
  % find hit times
  sim_results(n).micro_times = NaN(1,size(logL_mRNA,2));
  sim_results(n).mRNA_times = NaN(1,size(logL_mRNA,2));
  for m = 1:size(logL_mRNA,2)
      micro_ind = find(logL_micro(:,m)<=-K_target,1);
      mRNA_ind = find(logL_mRNA(:,m)<=-K_target,1);
      
      if ~isempty(micro_ind)
        sim_results(n).micro_times(m) = time_grid(micro_ind);
      end
      if ~isempty(mRNA_ind)
        sim_results(n).mRNA_times(m) = time_grid(mRNA_ind);
      end
  end

  sim_tau_vec_micro(n) = nanmean(sim_results(n).micro_times);
  sim_tau_vec_mRNA(n) = nanmean(sim_results(n).mRNA_times);
  
  
end

sim_results = rmfield(sim_results,'logL_mRNA');
sim_results = rmfield(sim_results,'logL_micro');
% save([DataPath 'decision_rate_2state.mat'],'sim_results')
%% %%%%%%%%%%%%%%% make validation figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
micro_drift_fig = figure;
cmap = brewermap(9,'Set2');

scatter(predicted_drift_vec_micro,sim_drift_vec_micro,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')

ylabel('mRNA logL drift rate (simulation)')
xlabel('mRNA logL drift rate (theoretical prediction)')
grid on

set(gca,'Fontsize',14)
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';

StandardFigurePBoC([],gca);
micro_drift_fig.InvertHardcopy = 'off';
saveas(micro_drift_fig,[FigurePath 'micro_logL_validation.png'])
saveas(micro_drift_fig,[FigurePath 'micro_logL_validation.pdf'])


mRNA_drift_fig = figure;

scatter(predicted_drift_vec_mRNA,sim_drift_vec_mRNA,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')

ylabel('mRNA logL drift rate (simulation)')
xlabel('mRNA logL drift rate (theoretical prediction)')
grid on

set(gca,'Fontsize',14)
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';

StandardFigurePBoC([],gca);
mRNA_drift_fig.InvertHardcopy = 'off';
saveas(mRNA_drift_fig,[FigurePath 'mRNA_logL_validation.png'])
saveas(mRNA_drift_fig,[FigurePath 'mRNA_logL_validation.pdf'])

% Now deicision times

% make validation figures
micro_time_fig = figure;

scatter(predicted_tau_vec_micro/60,sim_tau_vec_micro/60,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')

ylabel('simulated micro decision time (minutes)')
xlabel('predicted micro decision time (minutes)')
grid on

set(gca,'Fontsize',14)
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';

StandardFigurePBoC([],gca);
micro_time_fig.InvertHardcopy = 'off';
saveas(micro_time_fig,[FigurePath 'micro_decision_time_validation.png'])
saveas(micro_time_fig,[FigurePath 'micro_decision_time_validation.pdf'])


mRNA_time_fig = figure;

scatter(predicted_tau_vec_mRNA/60,sim_tau_vec_mRNA/60,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')

ylabel('simulated mRNA decision time (minutes)')
xlabel('predicted mRNA decision time (minutes)')
grid on

set(gca,'Fontsize',14)
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';

StandardFigurePBoC([],gca);
mRNA_time_fig.InvertHardcopy = 'off';
saveas(mRNA_time_fig,[FigurePath 'mRNA_decision_time_validation.png'])
saveas(mRNA_time_fig,[FigurePath 'mRNA_decision_time_validation.pdf'])

