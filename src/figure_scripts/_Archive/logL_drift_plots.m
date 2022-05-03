% script to generate illustrative log likelihood plots
clear
close all
addpath(genpath('../utilities'))

% make figure path
FigurePath = '../../fig/logL_drift_plots/';
mkdir(FigurePath)

%% %%%%%%%%%%%%%%%%%%%%%%% define core parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

T = 1e4;
n_conditions = 3;
n_reps = 10; % number of replicates per parameter set
c0 = 1;
c1_vec = c0*[1.05 1.1 1.15]; % range of "high" concentrations

c_true_vec = c1_vec; % For now make the simplifying assumption that lower concentration is the "true" one
koff = 2;
kon = 2;

% set target error rate to 1%
minError = .01; 
K_target = log((1-minError)/minError);

%% Calculate theoretical expectations for each condition
sim_results = struct;
for n = 1:length(c1_vec)
    
    c1 = c1_vec(n);
    c_true = c_true_vec(n);
    sim_results(n).c1 = c1;
    sim_results(n).c0 = c0;
    sim_results(n).kon = kon;
    sim_results(n).koff = koff;
    sim_results(n).ct = c_true_vec(n);
    % calculate analytic predictions
    [sim_results(n).r0, sim_results(n).r1, sim_results(n).rt, sim_results(n).v0, ...
      sim_results(n).v1, sim_results(n).vt, sim_results(n).V_mRNA, sim_results(n).D_mRNA,...
      sim_results(n).V_micro, sim_results(n).D_micro, sim_results(n).tau_micro, sim_results(n).tau_mRNA]...
                          = calculateDriftStats2State(kon,koff,c0,c1,c_true,K_target);
end                      
                  
% calculate simulation durations and resampling resolutions based upon 
% estimated decision times
T_vec = [sim_results.tau_micro];
resamp_increment = 10.^round(log10(T_vec/500));
tau_ceil = ceil(2*T_vec);

for n = 1:length(c1_vec)
  sim_results(n).T = max(tau_ceil);
  sim_results(n).dT = min(resamp_increment);
end

%% %%%%%%%%%%%%%%%%%%%%%%% Conduct simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%

state_options = 0:1;

% simulate trajectories and track movement of log likelihood ratio over
% time
rng(245)
% iterate
for n = 1:length(c1_vec)    
  
    % extract simulation parametes
    kon = sim_results(n).kon;
    koff = sim_results(n).koff;
    c1 = sim_results(n).c1;
    c_true = sim_results(n).ct;
    
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
              logL_micro(end+1) = logL_micro(end) + log(c1/c0)-kon*dt*(c1-c0);             
            else
              logL_micro(end+1) = logL_micro(end);
            end
        end    
        
        % downsample and store quantities of interest
        cumulative_time = cumsum(jump_vec);
        sim_results(n).total_mRNA(:,r) = interp1(cumulative_time,total_mRNA,time_grid);        
        sim_results(n).logL_micro(:,r) = interp1(cumulative_time,logL_micro,time_grid);               
    end    
end

%% %%%%%%%%%%%%%%%%%% Make movie of logL diffusion
% define colors
blue = [190 201 224]/255;
red = [246 141 100]/255;
green = [203 220 170]/255;
gray = [0.7020    0.7020    0.7020];
cmap1 = [brighten(green,-.2) ; blue ;red];


% path
frame_path = [FigurePath 'movie_frames' filesep];
mkdir(frame_path);

% movie parameters
sim_index = 2;
plot_res = 30;
plot_times = 0:30:60*60;
example_indices = [2 4 6];
master_iter = 1;

for e = 1:length(example_indices)
  
    trace_index = example_indices(e);
    prev_indices = example_indices(1:e-1);

    time_vec = sim_results(sim_index).time_grid;    
   
    last_time_index_raw = find(sim_results(sim_index).logL_micro(:,trace_index)>=K_target,1);
    if isempty(last_time_index_raw)
      last_time_index_raw = length(time_vec);
    end
    last_time_index = find(plot_times>=time_vec(last_time_index_raw),1);
    
    for p = 1:last_time_index
        close all
        currTime = plot_times(p);
        if p < last_time_index   
          last_index = max([1 find(time_vec>=currTime,1)]);
        else
          last_index = last_time_index_raw;
        end
          
        temp_fig = figure('Visible','off');
        hold on  
        
        % plot previous trends
        if ~isempty(prev_indices)
            plot(time_vec/60, sim_results(sim_index).logL_micro(:,prev_indices),'Color',[cmap1(sim_index,:) .5],'LineWidth',1); 
        end
        
        % plot log likelihood trends
        plot(time_vec(1:last_index)/60, sim_results(sim_index).logL_micro(1:last_index,trace_index),'Color',cmap1(sim_index,:),'LineWidth',1.5);   
        if p < last_time_index                          
          scatter(time_vec(last_index)/60, sim_results(sim_index).logL_micro(last_index,trace_index),75,'o','MarkerEdgeColor','k')
          scatter(time_vec(last_index)/60, sim_results(sim_index).logL_micro(last_index,trace_index),10,'o','MarkerEdgeColor','k','MarkerFaceColor','k')
        else
          scatter(time_vec(last_index)/60, sim_results(sim_index).logL_micro(last_index,trace_index),75,'o','MarkerEdgeColor','r')
          scatter(time_vec(last_index)/60, sim_results(sim_index).logL_micro(last_index,trace_index),10,'o','MarkerEdgeColor','r','MarkerFaceColor','r')
        end

        % plot the decision thresholds
        plot(time_grid/60,repelem(0,length(time_grid)),'-','Color','k','LineWidth',1);
        plot(time_grid/60,repelem(K_target,length(time_grid)),'--','Color','k','LineWidth',1.5);
        plot(time_grid/60,repelem(-K_target,length(time_grid)),'--','Color','k','LineWidth',1.5);

        ylim([-1.2*K_target 1.2*K_target])
        xlim([0 60])
        yl = ylabel('$log(\mathcal{L})$');
        yl.Interpreter = 'Latex';
        xlabel('time (minutes)')
        box on
        set(gca,'Fontsize',14)
        ax = gca;
        ax.YColor = 'black';
        ax.XColor = 'black';
        set(gca,'YTick',-5:2.5:5)
        grid on
        
        set(gca,'Color',[228,221,209]/255) 
        temp_fig.InvertHardcopy = 'off';
        set(gcf,'color','w');
        
        saveas(temp_fig,[frame_path 'logL_frame' sprintf('%03d',master_iter) '.tif'])
        master_iter = master_iter + 1; 
    end
end

%% 
close all

logL_fig = figure;
hold on
for i = 1:3
  % plot log liklihood trends
  p1 = plot(time_grid/60, sim_results(sim_index).logL_micro,'Color',[cmap1(sim_index,:) .7],'LineWidth',1);   
  if i > 1
    p2 = plot(time_grid/60, nanmean(sim_results(sim_index).logL_micro,2),'Color',brighten(cmap1(sim_index,:),-.6),'LineWidth',2.5);
  end
  if i > 2
    p3 = plot(time_grid/60, time_grid*sim_results(sim_index).V_micro,'-.','Color',brighten(cmap1(sim_index,:),-.8),'LineWidth',2.5);
  end
  
  % plot the decision thresholds
  plot(time_grid/60,repelem(0,length(time_grid)),'-','Color','k','LineWidth',1);
  plot(time_grid/60,repelem(K_target,length(time_grid)),'--','Color','k','LineWidth',1.5);
  plot(time_grid/60,repelem(-K_target,length(time_grid)),'--','Color','k','LineWidth',1.5);
  
  if i > 2
    legend([p1(1) p2 p3],'individual time series','simulation average','theoretical prediction','Color','w','Location','southeast')
  end
  
  ylim([-1.2*K_target 1.2*K_target])
  xlim([0 60])
  yl = ylabel('$log(\mathcal{L})$');
  yl.Interpreter = 'Latex';
  xlabel('time (minutes)')
  box on
  set(gca,'Fontsize',14)
  ax = gca;
  ax.YColor = 'black';
  ax.XColor = 'black';
  set(gca,'YTick',-5:2.5:5)
  grid on
  % set(gca,'Xscale','log')
  set(gca,'Color',[228,221,209]/255) 
  logL_fig.InvertHardcopy = 'off';
  set(gcf,'color','w');   

  saveas(logL_fig,[FigurePath 'spec_logL_traces' num2str(i) '.png'])
  saveas(logL_fig,[FigurePath 'spec_logL_traces' num2str(i) '.pdf'])
end
%% %%%%%%%%%%%%%%%%%% Plot of different c's
close all


logL_fig = figure;
hold on

% plot log liklihood trends
lgd_str = {};
for n = 1:length(c1_vec)
  plot(time_grid/60, sim_results(n).logL_micro,'Color',[cmap1(n,:) .7],'LineWidth',1);    
  lgd_str{n} = ['$\frac{c_1}{c_0} = ' num2str(c1_vec(n)) '$'];
end
for n = 1:length(c1_vec)
  p(n) = plot(time_grid/60, nanmean(sim_results(n).logL_micro,2),'Color',brighten(cmap1(n,:),-.7),'LineWidth',2.5);
end
% plot the decision thresholds
plot(time_grid/60,repelem(0,length(time_grid)),'-','Color','k','LineWidth',1);
plot(time_grid/60,repelem(K_target,length(time_grid)),'--','Color','k','LineWidth',1.5);
plot(time_grid/60,repelem(-K_target,length(time_grid)),'--','Color','k','LineWidth',1.5);

lgd = legend(p,lgd_str{:},'Location','southeast');
lgd.Interpreter = 'Latex';

ylim([-1.2*K_target 1.2*K_target])
xlim([0 60])
yl = ylabel('$log(\mathcal{L})$');
yl.Interpreter = 'Latex';
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14)
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
set(gca,'YTick',-5:2.5:5)
grid on
% set(gca,'Xscale','log')
set(gca,'Color',[228,221,209]/255) 
logL_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(logL_fig,[FigurePath 'all_logL_traces.png'])
saveas(logL_fig,[FigurePath 'all_logL_traces.pdf'])

%% make little exponential figure
exp_fig = figure;
hold on
exp_vec = exprnd(1/(c1_vec(1)*kon),1,1e4);
h = histogram(exp_vec,'Normalization','probability','FaceColor',brighten(blue,-.4),'EdgeAlpha',1);
plot(h.BinEdges,exp(-c1_vec(1)*kon.*h.BinEdges)/sum(exp(-c1_vec(1)*kon.*h.BinEdges)),'Color','k','LineWidth',1.5)

xlim([0 10])
ylabel('probability');
xl = xlabel('$\tau_{off}$ (seconds)');
xl.Interpreter = 'Latex';

box on
set(gca,'Fontsize',14)
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
% set(gca,'YTick',-5:2.5:5)
% grid on
% set(gca,'Xscale','log')
set(gca,'Color',[228,221,209]/255) 
exp_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(exp_fig,[FigurePath 'exp_dist.png'])
saveas(exp_fig,[FigurePath 'exp_dist.pdf'])



