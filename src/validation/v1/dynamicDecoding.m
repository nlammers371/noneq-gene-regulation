clear
close all
addpath('../utilities')
% script to validate assymptotic variance model for 2 state promoter
% and make c <--> m map figures
figPath = '../../fig/intro_figs/';
mkdir(figPath);
% define system parameters
kon = 1/120; % s^-1
koff = 1/30; % s^-1
c_vec = linspace(1,20,100);
% define simulation parameters
n_sim = 100;
T = 60*50;
Tres = 1;
% define m->c mapping function
c_fun = @(x) (x(1)*x(4)) / (x(3)*(x(2)-x(1)));
% initialize arrays
state_array = NaN(numel(c_vec),n_sim,T);
mRNA_array = zeros(numel(c_vec),n_sim,T);
c_guess_array = NaN(numel(c_vec),n_sim,T);
% check predicted mean reponse
c_grid_full = linspace(c_vec(1),c_vec(end),100);
mRNA_ref_profile = c_grid_full.*kon ./ (c_grid_full.*kon + koff);


% Perform simulations
r_vec = [0 1];
options = [1 2];
% iterate
for c = 1:numel(c_vec)
    K = [kon*c_vec(c) koff];
    tic
    for n = 1:n_sim     
        state_array(c,n,1) = ceil(2*rand());
        mRNA_array(c,n,1) = r_vec(state_array(c,n,1))*Tres;
        c_guess_array(c,n,1) = c_fun([mRNA_array(c,n,1),1,K(1)/c_vec(c),K(2)]);
        cs = state_array(c,n,1);
        for t = 2:Tres:T
            dt = exprnd(1/K(cs));
            less_than = (dt < Tres)*1;      
            s_vec = [cs options(options~=cs)];
            cs = s_vec(less_than+1);
            state_array(c,n,t) = cs;
            mRNA_array(c,n,t) = mRNA_array(c,n,t-1) + r_vec(cs)*Tres;                        
            c_guess_array(c,n,t) = c_fun([mRNA_array(c,n,t),T,K(1)/c_vec(c),K(2)]);
            % try alternative approach
            [~,mi] = min(abs(mRNA_array(c,n,t)-mRNA_ref_profile*t*Tres));
            c_guess_array(c,n,t) = c_grid_full(mi);
        end        
    end
    toc
end
%% Map to bcd levels

c_guess_array1 = NaN(numel(c_vec),n_sim,T);
for t = 1000%1:T
    mRNA_slice = mRNA_array(:,:,t);
    dist_array = pdist2(mRNA_slice(:),mRNA_ref_profile'*t*Tres);
    index_vec = 1:numel(mRNA_slice);
    [~,mi_vec] = min(dist_array,[],2);
    c_slice = c_guess_array1(:,:,t);
    c_slice(index_vec) = c_grid_full(mi_vec);
    c_guess_array1(:,:,t) = c_slice;
end


%% Make plots
close all
m_bins = linspace(1,max(mRNA_array(:)),numel(c_vec));
c_grid_array = repmat(c_vec',1,n_sim);
plot_time = 40;
mRNA_vec = reshape(mRNA_array(:,:,plot_time*60),[],1);
N = histcounts2(mRNA_vec,c_grid_array(:),m_bins',c_vec');
N_norm = N / sum(N(:));
N_sm = imgaussfilt(N_norm,1);
N_sm(N_sm==0) = NaN;
% make figure of input/output distribution
in_out_fig = figure;
cmap = flipud(brewermap(128,'Spectral'));
colormap(cmap);
p = pcolor(N_sm);
% p.EdgeAlpha = 0;
h = colorbar;
ylabel('hb mRNA (au)')
xlabel('[Bcd] (au)')
caxis([0,prctile(N_norm(:),99.5)])
grid on
saveas(in_out_fig,[figPath 'Bcd_in_hb_out_t' num2str(plot_time) '.tif'])

%%
slice_index = 50;
slice_times = (5:3:50)*60;
bcd_slice = N_sm(50,:);

bcd_profile_array = NaN(numel(slice_times),numel(c_vec)-1);
for t = 6%1:numel(slice_times)
    plane = mRNA_array(:,:,slice_times(t));
    N = imgaussfilt(histcounts2(plane(:),c_grid_array(:),m_bins',c_vec'),1);
    slice = N(slice_index,:);
    bcd_profile_array(t,:) = slice / sum(slice);
end
%%
bcd_inf_fig = figure;
bar(c_vec(1:end-1) + (c_vec(2)-c_vec(1))/2,bcd_slice / nansum(bcd_slice),1,...
    'FaceColor',cmap(115,:),'FaceAlpha',.5,'EdgeAlpha',0)
xlim([0 10])
grid on
xlabel('[Bcd]')
ylabel('share')

