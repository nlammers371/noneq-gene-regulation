% Plot results of sharpness parameter sweeps
clear 
close all
addpath(genpath('../utilities/'))

% %%%%%%%%%%%%%%%%  set relative read and write paths %%%%%%%%%%%%%%%%%%%%
DataPath = ['../../out/bivariate_parameter_sweeps_n4_OR' filesep];
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\';
% DropboxFolder = 'S:\Nick\Dropbox\Nonequilibrium\Nick\manuscript\';
FigPath = [DropboxFolder 'fourStateSweeps' filesep];
mkdir(FigPath);

% %%%%%%%%%%%%%%%%  get metric names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,metric_names] = calculateMetricsMultiState([]);
flux_index = find(strcmp(metric_names,'Flux'));
decision_rate_index = find(strcmp(metric_names,'DecisionRateNorm'));
decision_time_index = find(strcmp(metric_names,'DecisionTimeNorm'));
sharpness_index = find(strcmp(metric_names,'Sharpness'));
cycle_time_index = find(strcmp(metric_names,'CycleTime'));
precision_index = find(strcmp(metric_names,'Precision'));

% %%%%%%%%%%%%%%%%  Set plot parameters and constants %%%%%%%%%%%%%%%%%%%%

n_plot = 3e3; % number of points to plot
markerAlpha = 0.5; % marker transparency
markerSize = 40; % marker size

% %%%%%%%%%%%%%%%%  Load Sweep results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metric_one_name = metric_names{sharpness_index};
metric_two_name = metric_names{precision_index};

% generate load names
load_name_neq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
              '_eq0.mat'];   
load_name_eq = ['param_sweep_results_' metric_one_name '_' metric_two_name ...
              '_eq1.mat'];               

% load data
load([DataPath load_name_neq]);
load([DataPath load_name_eq]);

% process neq data
metric_array_neq = vertcat(sim_struct_neq.metric_array);   
sharpness_vec_neq = metric_array_neq(:,sharpness_index);
info_vec_neq = metric_array_neq(:,decision_rate_index);
plot_filter_neq = sharpness_vec_neq>=0;% & precision_vec >=0;

% process eq data
metric_array_eq = vertcat(sim_struct_eq.metric_array);   
sharpness_vec_eq = metric_array_eq(:,sharpness_index);
info_vec_eq = metric_array_eq(:,decision_rate_index);
plot_filter_eq = sharpness_vec_eq>=0;% & precision_vec >=0;

%  %%%%%%%%%%%%%%%%  Find optimal eq and non-eq architectures %%%%%%%%%%%%%
[~,best_eq_index] = nanmax(info_vec_eq.*plot_filter_eq);
[~,best_neq_index] = nanmax(info_vec_neq.*plot_filter_neq);

rate_array_eq = vertcat(sim_struct_eq.rate_array);
rate_array_neq = vertcat(sim_struct_neq.rate_array);

% get cycle times
tau_opt_eq = metric_array_eq(best_eq_index,cycle_time_index);
tau_opt_neq = metric_array_neq(best_neq_index,cycle_time_index);

% extract rates
opt_eq_rates = rate_array_eq(best_eq_index,:)*tau_opt_eq/60;
opt_neq_rates = rate_array_neq(best_neq_index,:)*tau_opt_neq/60;

%  %%%%%%%%%%%%%%%%%%%%%  Simulate network dynamics %%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('../'))
rmpath(genpath('../utilities/metricFunctions/'));
addpath(['../utilities/metricFunctions/n4_OR/']);
% generate path to save metric functions 
savePath = '../utilities/metricFunctions/n4_OR/';
mkdir(savePath);

% define some basic parameters. For convenience we will start in 6 state
% frame and then truncate
activeStates = [3 4];
nStates = 6; % this is a dummy variable. We will truncate to 4 states later on
sensingEdges = [2,1 ; 3,4];

% define symbolic variables
RSymFull = sym('k%d%d', [nStates nStates],'positive');

% zero out forbidden transitions
[fromRef,toRef] = meshgrid(1:nStates,1:nStates);
diffRef = abs(fromRef-toRef);
toRefHalved = toRef<=nStates/2;
fromRefHalved = fromRef<=nStates/2;
permittedConnections = (diffRef==1 & toRefHalved==fromRefHalved) | diffRef==nStates/2;

% permute these connections to follow a more intuitive labeling scheme
indexCurr = 1:nStates;
indexAdjusted = circshift(indexCurr,floor(nStates/4));
indexAdjusted = [indexAdjusted(1:nStates/2) fliplr(indexAdjusted(nStates/2+1:end))];
[~,si] = sort(indexAdjusted);
permittedConnections = permittedConnections(si,si);

% update transition matrix to reflect network topology
RSymRaw = RSymFull.*permittedConnections;
  
for i = 1:nStates
    for j = 1:nStates
        if i~=j 
            string = char(RSymRaw(i,j));
            if length(string)>1
                syms(string,'positive')
            end
        end
    end
end
        
%incorporate sensing edges   
syms c positive
sensingIndices = sub2ind(size(RSymRaw),sensingEdges(:,1),sensingEdges(:,2));
RSym = RSymRaw;
RSym(sensingIndices) = c*RSym(sensingIndices);

% Now truncate
RSym = RSym(1:4,1:4);
nStates = 4;
% add diagonal terms
RSym(eye(size(RSym))==1) = -sum(RSym);

% turn this into a function
RSymFun = matlabFunction(RSym);


rate_array = (vertcat(opt_neq_rates, opt_eq_rates));
% n_sim = 10;
n_reps = 50;
t_sim = 300*60;
time_grid = 1:10:t_sim;
state_options = 1:4;
pd_states = [3 4];

% initialize arrays
sim_struct = struct;
id_cell = {'non-eq optimum','eq optimum'};
c_val = 1;

% initialize waitbar
% h = waitbar(0,'Simulating transcription network dynamics...');
simPath = [FigPath 'logL_sim_struct.mat'];
if ~exist(simPath)

    rng(123);
    % iterate through time points  
    for n = 1:2%:n_sim
    %   waitbar(n/n_sim,h);          
      tic
      % calculate steady state stats  
      sim_struct(n).inputMat = [c_val rate_array(n,:)];
      sim_struct(n).id = id_cell{n};
      inputCell = mat2cell(sim_struct(n).inputMat,size(sim_struct(n).inputMat,1),ones(1,size(sim_struct(n).inputMat,2)));

      R_matrix = RSymFun(inputCell{:});    
      [V,D] = eig(R_matrix);
      [~,mi] = max(real(diag(D)));
      ss = V(:,mi)/sum(V(:,mi));
      jump_weight_array = R_matrix;
      jump_weight_array(eye(4)==1) = 0;

      % extract state dwell vector
      dwell_time_vec = -diag(R_matrix)'.^-1;       

      % initilize array to store cumulative mRNA
      total_mRNA_array = NaN(length(time_grid),n_reps);

      parfor rep = 1:n_reps

        % initialize vector to store results
        jump_time_vec = [];
        state_vec = [randsample(state_options,1,true,ss)];    
        total_time = 0;   
        tic
    %     wb = waitbar(0,'simulating network dynamics...');
        while total_time < t_sim     
    %       waitbar(total_time/t_sim,wb);
          % extract state info
          state_curr = state_vec(end);

          % calculate jump time means
          tau = dwell_time_vec(state_curr);    

          % randomly draw jump times
          dt = exprnd(tau);
          total_time = total_time + dt;

          % randomly draw destination state            
          next_state = randsample(state_options,1,true,jump_weight_array(:,state_curr));      

          % assigne states
          state_vec(end+1) = next_state;
          jump_time_vec(end+1) = dt;    
        end  
        delete(wb)
        % downsample and save
        cumulative_time = cumsum(jump_time_vec);
        cumulative_mRNA = cumsum(jump_time_vec.*double(ismember(state_vec(1:end-1),pd_states)));
        total_mRNA_array(:,rep) = interp1(cumulative_time,cumulative_mRNA,time_grid);
      end  

      % record average values
      sim_struct(n).r_mean = mean(total_mRNA_array(end,:))/time_grid(end);
      sim_struct(n).r_var = var(total_mRNA_array(end,:))/time_grid(end);
      sim_struct(n).mRNA_array = total_mRNA_array;

      toc
    end


% tau_vec = [60/tau_opt_neq 60/tau_opt_eq];

    for n = 1:2
        inputMat0 = [c_val rate_array(n,:)];  
        inputMat1 = [c_val*1.1 rate_array(n,:)];  
        inputCell1 = mat2cell(inputMat1,size(inputMat1,1),ones(1,size(inputMat1,2)));
        inputCell0 = mat2cell(inputMat0,size(inputMat0,1),ones(1,size(inputMat0,2)));

        % calculate production rate         
        r1 = productionRateFunction(inputCell1{:});
        r0 = productionRateFunction(inputCell0{:});

        % calculate variance    
        v1 = intrinsicVarianceFunction(inputCell1{:});%*tau_vec(n);
        v0 = intrinsicVarianceFunction(inputCell0{:});%*tau_vec(n);

        sim_struct(n).m1_vec = (r1*time_grid');
        sim_struct(n).v1_vec = v1*time_grid';

        sim_struct(n).m0_vec = (r0*time_grid');
        sim_struct(n).v0_vec = v0*time_grid';

        % calculate (un-normalized) log likelihoods
        sim_struct(n).logL0 = -0.5*(sim_struct(n).m0_vec-sim_struct(n).mRNA_array).^2./sim_struct(n).v0_vec;
        sim_struct(n).logL1 = -0.5*(sim_struct(n).m1_vec-sim_struct(n).mRNA_array).^2./sim_struct(n).v1_vec;
        sim_struct(n).logL_delta = sim_struct(n).logL0-sim_struct(n).logL1;
    end    

    disp('saving...')
    save(simPath,'sim_struct')
else
    load(simPath);
end

%% Make plots
close all

plot_indices_eq = [2 7 12];
plot_indices_neq = [1 2 11];

minError = .1; 
K_target = log((1-minError)/minError);

% find absorption times
logLDeltaPlotNeq = sim_struct(1).logL_delta;
absorption_times_neq = NaN(1,n_reps);
absorption_signs_neq = NaN(1,n_reps);
for n = 1:n_reps
    at = find(logLDeltaPlotNeq(:,n)>K_target | logLDeltaPlotNeq(:,n)<-K_target,1);
    if ~isempty(at)
        absorption_times_neq(n) = at;
        logLDeltaPlotNeq(at+1:end,n) = NaN;
        if logLDeltaPlotNeq(at,n)<-K_target
            absorption_signs_neq(n) = -K_target;
        else
            absorption_signs_neq(n) = K_target;
        end
    end
end

logLDeltaPlotEq = sim_struct(2).logL_delta;
absorption_times_eq = NaN(1,n_reps);
absorption_signs_eq = NaN(1,n_reps);
for n = 1:n_reps
    at = find(logLDeltaPlotEq(:,n)>K_target | logLDeltaPlotEq(:,n)<-K_target,1);
    if ~isempty(at)
        absorption_times_eq(n) = at;
        logLDeltaPlotEq(at+1:end,n) = NaN;
        if logLDeltaPlotEq(at,n)<-K_target
            absorption_signs_eq(n) = -K_target;
        else
            absorption_signs_eq(n) = K_target;
        end
            
    end
end

for i = 2
    logL_fig = figure;
    cmap = brewermap(8,'Set2');
    hold on

    % plot the decision thresholds
    plot(time_grid/60,repelem(1,length(time_grid)),'-','Color','k','LineWidth',1);
    plot(time_grid/60,repelem(exp(K_target),length(time_grid)),'--','Color','k','LineWidth',1.5);
    plot(time_grid/60,repelem(exp(-K_target),length(time_grid)),'--','Color','k','LineWidth',1.5);

    % plot trends
    if i == 2
        p1 = plot(time_grid/60,exp(logLDeltaPlotNeq(:,plot_indices_neq)),'-','Color',[cmap(2,:) 1],'LineWidth',3);
    end
    p2 = plot(time_grid/60,exp(logLDeltaPlotEq(:,plot_indices_eq)),'-','Color',[cmap(3,:) 1],'LineWidth',3);

    if i == 2        
        scatter(absorption_times_neq(plot_indices_neq)/6,exp(absorption_signs_neq(plot_indices_neq))...
                                  ,'MarkerEdgeColor','k','MarkerFaceColor','w')
    end
    scatter(absorption_times_eq(plot_indices_eq)/6,exp(absorption_signs_eq(plot_indices_eq))...
                                  ,'MarkerEdgeColor','k','MarkerFaceColor','w')                            

    if i == 2
        legend([p1(1) p2(2)],'non-equilibrium','equilibrium','Location','southwest')
    else
        legend([p2(2)],'equilibrium','Location','southwest')
    end
    ylim([minError/5 5/minError])
    % ylim([-1.1*K_target 1.1*K_target])
    xlim([0 240])
    set(gca,'Yscale','log')
    grid on
    yl = ylabel('$P_0/P_1$');
    % yl = ylabel('$log(\mathcal{L})$');
    yl.Interpreter = 'Latex';
    xlabel('time (minutes)')
    box on
    set(gca,'Fontsize',14)
    ax = gca;
    ax.YColor = 'black';
    ax.XColor = 'black';
    % set(gca,'YTick',[10 100 1000])
    set(gca,'Color',[228,221,209]/255) 
    logL_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    if i == 2      
        saveas(logL_fig,[FigPath 'logL_plot_eq_vs_neq.png'])
        saveas(logL_fig,[FigPath 'logL_plot_eq_vs_neq.pdf'])
    else
        saveas(logL_fig,[FigPath 'logL_plot_eq.png'])
        saveas(logL_fig,[FigPath 'logL_plot_eq.pdf'])
    end
end