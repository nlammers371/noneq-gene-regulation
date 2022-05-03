% script to generate illustrative log likelihood plots
clear
close all
addpath(genpath('../utilities'))

% make figure path
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\Figure 1\';

% path for profile movie
ProfileMoviePath = [DropboxFolder 'ProfileMovie' filesep];
mkdir(ProfileMoviePath)

% path for decision point trends
DecisionPointMoviePath = [DropboxFolder 'DecisionPointMovie' filesep];
mkdir(DecisionPointMoviePath)

%% %%%%%%%%%%%%%%%%%%%%%%% define core parameters %%%%%%%%%%%%%%%%%%%%%%%%%
rng(245); % initialize random number generator for consistency
T = 1e3; % simulation time
n_conditions = 3;
n_reps = 10; % number of replicates per parameter set

c0 = 1;
c1_vec = c0*[1.05 1.1 1.15]; % range of "high" concentrations
c_range = unique([linspace(0,2) c1_vec]);

c_true_vec = c1_vec; % For now make the simplifying assumption that lower concentration is the "true" one
koff = 2;
kon = 2;

% set target error rate to 1%
minError = .01; 
K_target = log((1-minError)/minError);

% %%%%%%%%%%%% "simulate" full profile and point trends %%%%%%%%%%%%%%%%%%
close all
profile_sigma = 2*kon.*c_range*koff ./ (kon.*c_range+koff).^3;
profile_rate = kon.*c_range ./ (kon.*c_range+koff);
profile_array = NaN(length(c_range),n_reps,T);
% for t = 1:T
%     for n = 1:n_reps
%         if t == 1
%             profile_array(:,n,t) = normrnd(profile_rate',profile_sigma');
%         else
%             profile_array(:,n,t) = profile_array(:,n,t-1)+ normrnd(profile_rate',profile_sigma'./sqrt(t));
%         end
%     end
% end


state_options = [0 1];
time_grid = 1:T;
wb = waitbar(0,'simulating 2 state trajectories...');
for c = 1:length(c_range)
    waitbar(c/T,wb)
    kon_curr = c_range(c)*kon+1e-10;
    tau_vec = [1/kon_curr 1/koff];
    for r = 1:n_reps

        % initialize temporary structures                
        total_mRNA = [0];        
        total_time = 0;
        jump_vec = [0];
        % initialize binding state using SS vector    
        state_vec = [randsample(state_options,1,true,[koff*tau_vec(1) kon_curr*tau_vec(2)])];%,[1-rt, rt])];

        while total_time < T

            current_state = state_vec(end);

            % draw jump time and toggle state
            dt = exprnd(tau_vec(current_state+1),1);
            next_state = state_options(state_options~=current_state);

            total_time = total_time + dt;

            jump_vec(end+1) = dt;
            state_vec(end+1) = next_state; 
            total_mRNA(end+1) = total_mRNA(end) + dt*current_state;         
        end    

        % downsample and store quantities of interest
        cumulative_time = cumsum(jump_vec);
        profile_array(c,r,:) = interp1(cumulative_time,total_mRNA,time_grid);                            
    end     
end
delete(wb);

%%
close all
n = 3;
[~, c1_index] = min(abs(c_range-c1_vec(n)));
[~, c0_index] = min(abs(c_range-c0));

% make profile and point plot movie frames
wb = waitbar(0,'printing simulated mRNA frames...');
for t = 1:10:T
    close all
    
    waitbar(t/T,wb)
    
    profile_frame = figure('Position',[200 200 512 256],'Visible','off');
    cmap = brewermap([],'Paired');
    cmap1 = brewermap([],'Set2');
    hold on
    plot(c_range,profile_array(:,:,t)/t,'Color',cmap1(8,:),'LineWidth',1)  
    scatter(repelem(c_range(c1_index),n_reps),profile_array(c1_index,:,t)/t,'MarkerEdgeColor',cmap(9,:),'LineWidth',2)  
    scatter(repelem(c_range(c0_index),n_reps),profile_array(c0_index,:,t)/t,'MarkerEdgeColor',cmap(10,:),'LineWidth',2)  
    set(gca,'FontSize',14)
%     set(gca, 'xtick', [],'ytick', [])
    xlabel('activator concentration (C)')
    ylabel('normalized protein (M/t)')
    set(gca,'Color',[228,221,209]/255) 

    ax = gca;
    % ax.XAxis.MinorTickValues = 10.^(-5:.1:5);
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    profile_frame.InvertHardcopy = 'off';
    set(gcf,'color','w');
    ylim([0 .8])
%     text(1.75,.05,[sprintf('%03d',round(t/60)) ' min'])
    saveas(profile_frame,[ProfileMoviePath 'profile_frame' sprintf('%04d',t) '.tif'])
    if t == T
        saveas(profile_frame,[DropboxFolder 'profile_frame' sprintf('%04d',t) '.pdf'])
    end
    
    
    point_frame = figure('Position',[100 100 512 256],'Visible','off');    
    
    hold on
    t_vec = 1:t;
    c1_activity_vec = permute(profile_array(c1_index,:,1:t),[3 2 1])./(1:t)';
    c0_activity_vec = permute(profile_array(c0_index,:,1:t),[3 2 1])./(1:t)';
    p1 = plot(t_vec/60,c1_activity_vec,'-','Color',cmap(9,:),'LineWidth',1);
    p0 = plot(t_vec/60,c0_activity_vec,'-','Color',cmap(10,:),'LineWidth',1);  
      
    legend([p1(1) p0(1)],'c1','c0','Location','northeast')
    set(gca,'FontSize',14)
%     set(gca, 'xtick', [],'ytick', [])
    xlabel('time')
    ylabel('normalized protein (M/t)')
    set(gca,'Color',[228,221,209]/255) 
%     ylim([0 1.05*profile_rate(c1_index)*T])
    ylim([.25 .75])
    xlim([0 T/60])
    ax = gca;
    % ax.XAxis.MinorTickValues = 10.^(-5:.1:5);
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    point_frame.InvertHardcopy = 'off';    
    set(gcf,'color','w');   
%     saveas(point_frame,[DecisionPointMoviePath 'point_frame' sprintf('%04d',t) '.tif'])
    if t == T
        saveas(point_frame,[DropboxFolder 'point_frame' sprintf('%04d',t) '.pdf'])
    end
        
end
delete(wb);
%% Calculate theoretical expectations for each condition
n = 3;

c1 = c1_vec(n);
c_true = c_true_vec(n);

[~, c1_index] = min(abs(c_range-c1_vec(n)));
[~, c0_index] = min(abs(c_range-c0));

% calculate analytic predictions
[r0, r1, rt, v0, ...
  v1, vt, V_mRNA, D_mRNA,...
  V_micro, D_micro, tau_micro, tau_mRNA]...
                      = calculateDriftStats2State(kon,koff,c0,c1,c_true,K_target);
                 

% calculate logL for simulated trajectories 

time_grid = 1:T;

% extract sumulated mRNA trends
c0_mRNA = permute(profile_array(c0_index,:,:),[3 2 1]);

m1_vec = (r1*time_grid');
sigma1_vec = v1*time_grid';

m0_vec = (r0*time_grid');
sigma0_vec = v0*time_grid';

% calculate (un-normalized) log likelihoods
logL0 = -0.5*(m0_vec-c0_mRNA).^2./sigma0_vec;
logL1 = -0.5*(m1_vec-c0_mRNA).^2./sigma1_vec;
logLDelta = logL0 - logL1;
logLDelta_orig = logLDelta;
for n = 1:n_reps
    at = find(logLDelta(:,n)>K_target | logLDelta(:,n)<-K_target,1);
    if ~isempty(at)
        absorption_times(n) = at;
        logLDelta(at+1:end,n) = NaN;
    end
end


% Make time series plots of decision process
close all
plot_indices = 4;
plot_time_vec = [50 100 500 900];
r_vec = linspace(0,1,1e3);

for p = 1:4%1:length(plot_time_vec)
    plot_time = plot_time_vec(p);
    
    logL_fig = figure('Position',[100 100 256 256]);
    cmap = brewermap([],'Paired');
    hold on

    % plot the decision thresholds
    plot(time_grid,repelem(1,length(time_grid)),'-','Color','k','LineWidth',1);
%     plot(time_grid/60,repelem(exp(K_target),length(time_grid)),'--','Color','k','LineWidth',1.5);
%     plot(time_grid/60,repelem(exp(-K_target),length(time_grid)),'--','Color','k','LineWidth',1.5);

    % plot trends
%     logLDeltaPlot = logLDelta(:,plot_indices);
    plot(time_grid(1:plot_time),exp(logLDelta_orig(1:plot_time,plot_indices)),'-','Color','k','LineWidth',3)
    scatter(plot_time,exp(logLDelta_orig(plot_time,plot_indices)),'MarkerEdgeColor','k','MarkerFaceColor','w')
%     scatter(absorption_times(plot_indices)/60,repelem(K_target,length(plot_indices)),'MarkerEdgeColor','k','MarkerFaceColor','w')
    % plot(time_grid/60,time_grid*V_mRNA,'-','Color',cmap(8,:),'LineWidth',2)

%     ylim([-1.1*K_target 1.1*K_target])
    
    ylim([10^-4 10^4])
    yl = ylabel('$P_0/P_1$');
    yl.Interpreter = 'Latex';
    xlabel('time')
    box on
    set(gca,'Fontsize',14)
    ax = gca;
    ax.YColor = 'black';
    ax.XColor = 'black';
    set(gca,'YTick',[1e-4 1e-2 1e0 1e2 1e4])
    set(gca,'Yscale','log')
    grid on
    set(gca,'Color',[228,221,209]/255) 
    xlim([0 T])
    logL_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    
    saveas(logL_fig,[DropboxFolder 'logL_plot_t' sprintf('%04d',plot_time) '.png'])
    saveas(logL_fig,[DropboxFolder 'logL_plot_t' sprintf('%04d',plot_time) '.pdf'])            
  
  
    %%%%%% Probability plot
    
    prob_fig = figure('Position',[100 200 256 256]);    
    % cmap1 = brewermap([],'set2');
    hold on    
    % c0_mRNA(plot_time,plot_indices)
    p0_vec = normpdf(r_vec,r0,sqrt(v0/plot_time));%gaussmf(r_vec,[sqrt(v0/plot_time) r0]);%exp(-0.5*(r0-r_vec).^2./v0*plot_time) .* sqrt(plot_time) ./ sqrt(2*pi*v0);    
    p0_vec = p0_vec/sum(p0_vec);
    p1_vec = normpdf(r_vec,r1,sqrt(v1/plot_time));
    p1_vec = p1_vec/sum(p1_vec);

    % plot trends
    f0 = fill([r_vec fliplr(r_vec)],[p0_vec,zeros(size(p0_vec))],cmap(10,:),'LineWidth',2,'FaceAlpha',0.5);%,'EdgeColor',cmap(2,:))
    f1 = fill([r_vec fliplr(r_vec)],[p1_vec,zeros(size(p1_vec))],cmap(9,:),'LineWidth',2,'FaceAlpha',0.5);%,'EdgeColor',cmap(8,:))
    % scatter(absorption_times(plot_indices)/60,repelem(K_target,length(plot_indices)),'MarkerEdgeColor','k','MarkerFaceColor','w')

    r_sim = c0_mRNA(plot_time,plot_indices)/plot_time;
    [rt,rti] = min(abs(r_vec-r_sim));

    plot(repmat(r_sim,100,1),linspace(0,0.04),'--','Color','k','LineWidth',2)

    plot(r_vec(1:rti),repelem(p0_vec(rti),rti),'--','Color','k','LineWidth',2)
    scatter(r_sim,p0_vec(rti),'MarkerEdgeColor','k','MarkerFaceColor',cmap(3,:))

    plot(r_vec(1:rti),repelem(p1_vec(rti),rti),'--','Color','k','LineWidth',2)
    scatter(r_sim,p1_vec(rti),'MarkerEdgeColor','k','MarkerFaceColor',cmap(5,:))

%     text(.66, .003, ['t = ' num2str(plot_time)],'Fontsize',14)    
    legend([f0 f1],'m0 (right)','m1 (wrong)')
    ylim([0 .04])
    xlim([.35 .65])
    grid on
    xlabel('normalized protein (M/t)')
    ylabel('probability');
    % xlabel('time (minutes)')
    box on
    set(gca,'Fontsize',14)
    ax = gca;
    ax.YColor = 'black';
    ax.XColor = 'black';
    % set(gca,'YTick',-5:2.5:5)
    set(gca,'Color',[228,221,209]/255) 
    prob_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    
    saveas(prob_fig,[DropboxFolder 'prob_dist_t' sprintf('%04d',plot_time) '.png'])
    saveas(prob_fig,[DropboxFolder 'prob_dist_t' sprintf('%04d',plot_time) '.pdf'])
    
    
%     first_fig = figure('Position',[100 200 512 256]);    
%     % cmap1 = brewermap([],'set2');
%     hold on    
% 
%     r_sim = c0_mRNA(plot_time,plot_indices)/plot_time;
%     [rt,rti] = min(abs(r_vec-r_sim));
% 
%     plot(repmat(r_sim,100,1),linspace(0,0.04),'--','Color','k','LineWidth',2)
% 
% %     plot(r_vec(1:rti),repelem(p0_vec(rti),rti),'--','Color','k','LineWidth',2)
% %     scatter(r_sim,p0_vec(rti),'MarkerEdgeColor','k','MarkerFaceColor',cmap(3,:))
% % 
% %     plot(r_vec(1:rti),repelem(p1_vec(rti),rti),'--','Color','k','LineWidth',2)
% %     scatter(r_sim,p1_vec(rti),'MarkerEdgeColor','k','MarkerFaceColor',cmap(5,:))
% 
% 
% %     legend([f0 f1],'m0 (right)','m1 (wrong)')
%     text(.66, .003, ['t = ' num2str(plot_time)],'Fontsize',14)
%     ylim([0 .04])
%     xlim([.25 .75])
%     grid on
%     xlabel('normalized protein (M/t)')
%     ylabel('probability');
%     % xlabel('time (minutes)')
%     
%     box on
%     set(gca,'Fontsize',14)
%     ax = gca;
%     ax.YColor = 'black';
%     ax.XColor = 'black';
%     % set(gca,'YTick',-5:2.5:5)
%     
%     set(gca,'Color',[228,221,209]/255) 
%     first_fig.InvertHardcopy = 'off';
%     set(gcf,'color','w');
%     
%     saveas(first_fig,[DropboxFolder 'rsim_dist_t' sprintf('%04d',plot_time) '.png'])
%     saveas(first_fig,[DropboxFolder 'rsim_dist_t' sprintf('%04d',plot_time) '.pdf'])
end

%% Make additional logL plots

% make plot
close all
plot_indices = [4];

logL_fig = figure('Position',[100 200 512 256]);
hold on

% plot the decision thresholds
plot(time_grid/60,repelem(1,length(time_grid)),'-','Color','k','LineWidth',1);
plot(time_grid/60,repelem(exp(K_target),length(time_grid)),'--','Color','k','LineWidth',1.5);
plot(time_grid/60,repelem(exp(-K_target),length(time_grid)),'--','Color','k','LineWidth',1.5);
        
% plot trends
logLDeltaPlot = logLDelta_orig(:,plot_indices);
plot(time_grid/60,exp(logLDeltaPlot),'-','Color','k','LineWidth',3)
scatter(absorption_times(plot_indices)/60,repelem(exp(K_target),length(plot_indices)),'MarkerEdgeColor','k','MarkerFaceColor','w')
% plot(time_grid/60,time_grid*V_mRNA,'-','Color',cmap(8,:),'LineWidth',2)

ylim([10^-4 10^4])
yl = ylabel('$P_0/P_1$');
yl.Interpreter = 'Latex';
xlabel('time')
box on
set(gca,'Fontsize',14)
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
set(gca,'YTick',[1e-4 1e-2 1e0 1e2 1e4])
set(gca,'Yscale','log')
grid on
set(gca,'Color',[228,221,209]/255) 
xlim([0 T/60])
logL_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(logL_fig,[DropboxFolder 'logL_plot_with_K.png'])
saveas(logL_fig,[DropboxFolder 'logL_plot_with_K.pdf'])


other_indices = ~ismember(1:n_reps,plot_indices);

logL_fig = figure('Position',[100 200 512 256]);
hold on

% plot the decision thresholds
plot(time_grid/60,repelem(1,length(time_grid)),'-','Color','k','LineWidth',1);
plot(time_grid/60,repelem(exp(K_target),length(time_grid)),'--','Color','k','LineWidth',1.5);
plot(time_grid/60,repelem(exp(-K_target),length(time_grid)),'--','Color','k','LineWidth',1.5);
        
% plot trends
logLDeltaPlot = logLDelta_orig(:,plot_indices);
plot(time_grid/60,exp(logLDeltaPlot),'-','Color',[0 0 0 1],'LineWidth',3)
scatter(absorption_times(plot_indices)/60,repelem(exp(K_target),length(plot_indices)),'MarkerEdgeColor','k','MarkerFaceColor','w')

logLDeltaPlotAlt = logLDelta_orig(:,other_indices);
plot(time_grid/60,exp(logLDeltaPlotAlt),'-','Color',[0 0 0 .2],'LineWidth',1)
% scatter(absorption_times(other_indices)/60,repelem(exp(K_target),sum(other_indices)),'MarkerEdgeColor','k','MarkerFaceColor','w')
% plot(time_grid/60,time_grid*V_mRNA,'-','Color',cmap(8,:),'LineWidth',2)

 ylim([10^-4 10^4])
yl = ylabel('$P_0/P_1$');
yl.Interpreter = 'Latex';
xlabel('time')
box on
set(gca,'Fontsize',14)
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
set(gca,'YTick',[1e-4 1e-2 1e0 1e2 1e4])
set(gca,'Yscale','log')
grid on
set(gca,'Color',[228,221,209]/255) 
xlim([0 T/60])
logL_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(logL_fig,[DropboxFolder 'logL_plot_with_K_others.png'])
saveas(logL_fig,[DropboxFolder 'logL_plot_with_K_others.pdf'])


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