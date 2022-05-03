% script to generate illustrative log likelihood plots
clear
close all
addpath(genpath('../utilities'))

% make figure path
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\SPRT\';

% path for profile movie
ProfileMoviePath = [DropboxFolder 'ProfileMovie' filesep];
mkdir(ProfileMoviePath)

%% %%%%%%%%%%%%%%%%%%%%%%% define core parameters %%%%%%%%%%%%%%%%%%%%%%%%%
rng(245); % initialize random number generator for consistency
T = 1e3; % simulation time
n_conditions = 3;
n_reps = 100; % number of replicates per parameter set

c0 = 1;
c1_vec = c0*[1.05 1.1 1.15]; % range of "high" concentrations
c_range = c1_vec;%unique([linspace(0,2) c1_vec]);

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

state_options = [0 1];
time_grid = 1:T;
wb = waitbar(0,'simulating 2 state trajectories...');
for c = 1:length(c_range)
    waitbar(c/length(c_range),wb)
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
for t = unique([1:10:T T])
    close all
    
    waitbar(t/T,wb)
    
%     profile_frame = figure('Position',[200 200 512 256]);%,'Visible','off');
%     cmap = brewermap([],'Paired');
%     cmap1 = brewermap([],'Set2');
%     hold on
%     plot(c_range,profile_array(:,:,t)/t,'Color',cmap1(8,:),'LineWidth',1)  
%     scatter(repelem(c_range(c1_index),n_reps),profile_array(c1_index,:,t)/t,'MarkerEdgeColor',cmap(9,:),'LineWidth',2)  
%     scatter(repelem(c_range(c0_index),n_reps),profile_array(c0_index,:,t)/t,'MarkerEdgeColor',cmap(10,:),'LineWidth',2)  
%     set(gca,'FontSize',14)
% %     set(gca, 'xtick', [],'ytick', [])
%     xlabel('activator concentration (C)')
%     ylabel('normalized protein (M/t)')
%     set(gca,'Color',[228,221,209]/255) 
%     set(gca,'xscale','log')
%     ax = gca;
%     % ax.XAxis.MinorTickValues = 10.^(-5:.1:5);
%     ax.YAxis(1).Color = 'k';
%     ax.XAxis(1).Color = 'k';
%     profile_frame.InvertHardcopy = 'off';
%     set(gcf,'color','w');
%     ylim([0 .8])
% %     text(1.75,.05,[sprintf('%03d',round(t/60)) ' min'])
%     saveas(profile_frame,[ProfileMoviePath 'profile_frame' sprintf('%04d',t) '.tif'])
%     if t == T
%         saveas(profile_frame,[DropboxFolder 'profile_frame' sprintf('%04d',t) '.pdf'])
%     end
%     
    
    point_frame = figure('Position',[100 100 512 256],'Visible','off'); 
    cmap = brewermap([],'Paired');
    cmap1 = brewermap([],'Set2');
    y_max = ceil(nanmax(profile_array(c1_index,:,end))/100)*100;
    hold on
    t_vec = 1:t;
    c1_activity_vec = permute(profile_array(c1_index,:,1:t),[3 2 1]);
    c0_activity_vec = permute(profile_array(c0_index,:,1:t),[3 2 1]);
    p1 = plot(t_vec/60,c1_activity_vec,'-','Color',[cmap(9,:) .5],'LineWidth',1);
    p0 = plot(t_vec/60,c0_activity_vec,'-','Color',[cmap(10,:) .5],'LineWidth',1);  
      
    legend([p1(1) p0(1)],'c1','c0','Location','southeast')
    set(gca,'FontSize',14)
%     set(gca, 'xtick', [],'ytick', [])
    xlabel('time')
    ylabel('accumulated mRNA')
    set(gca,'Color',[228,221,209]/255) 
    ylim([0 y_max])
    xlim([0 T/60])
    ax = gca;
    
    % ax.XAxis.MinorTickValues = 10.^(-5:.1:5);
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    point_frame.InvertHardcopy = 'off';    
    set(gcf,'color','w');   
    saveas(point_frame,[ProfileMoviePath 'point_frame' sprintf('%04d',t) '.tif'])
    if t == T
        saveas(point_frame,[DropboxFolder 'point_frame' sprintf('%04d',t) '.pdf'])
    end
        
end
delete(wb);
