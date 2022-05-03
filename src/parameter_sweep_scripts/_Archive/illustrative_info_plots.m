% this script generates illustrative plots for information metrric
clear 
close all

FigPath = '../../fig/intrinsic_noise_plots/';
mkdir(FigPath);

% basic sim parameters
n_examples = 5; 
mid_ex = ceil(n_examples/2);
n_range = [0 logspace(0,1,n_examples-1)];%0:2:2*n_examples;
c_vec = linspace(0,1,101);
snr_range = fliplr(logspace(-1.5,-.2,n_examples));
noise_vec = 0.5 .* snr_range;

% set plot colors
pboc = [228 221 210]/256;

% function to generate sigmoids
sigmoid_fun = @(c,n) c.^n ./ (mean(c).^n + c.^n);

% initialize array
sigmoid_array = NaN(length(c_vec),n_examples);
s_array = NaN(n_examples,2);
n_array = NaN(n_examples,2);
% generate sigmoids
for i = 1:n_examples
  sigmoid_array(:,i) = sigmoid_fun(c_vec,n_range(i));
  dfs = diff(sigmoid_array(:,i));
  % calculate information rates
  s_array(i,1) = mean([dfs(49) dfs(50)]);  
  if i == mid_ex
    s_array(:,2) = mean([dfs(49) dfs(50)]);
  end
end
n_array(:,1) = 1./noise_vec(mid_ex);
n_array(:,2) = 1./noise_vec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate sharpness modulation figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate example porfiles
sfig1 = figure;
hold on
cmap1 = brewermap(n_examples,'Blues');
cmap2 = brewermap(n_examples,'Reds');
cmap3 = brewermap(9,'Set2');
% generate noise fill ranges
for i = 1:n_examples
  ub = sigmoid_array(:,i) + noise_vec(mid_ex);
  lb = sigmoid_array(:,i) - noise_vec(mid_ex);
  fill([c_vec fliplr(c_vec)],[ub' fliplr(lb')],cmap1(i,:),'EdgeAlpha',0,'FaceAlpha',0.3)
end
% plot respons curves
for i = 1:n_examples
  plot(c_vec,sigmoid_array(:,i),'Color',cmap1(i,:),'LineWidth',1.5)
end
xlabel('activator concentration (c)')
ylabel('transcription rate (r)')
set(gca,'FontSize',14)
set(gca,'Color',pboc)
box on
ylim([-.05 1.05])
sfig1.InvertHardcopy = 'off';
saveas(sfig1,[FigPath 'sensitivity_mod_profiles.png'])
saveas(sfig1,[FigPath 'sensitivity_mod_profiles.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate noise modulation figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate example profiles
nfig1 = figure;
hold on
% generate noise fill ranges
for i = 1:n_examples
  ub = sigmoid_array(:,mid_ex) + noise_vec(i);
  lb = sigmoid_array(:,mid_ex) - noise_vec(i);
  fill([c_vec fliplr(c_vec)],[ub' fliplr(lb')],cmap2(i,:),'EdgeAlpha',0,'FaceAlpha',0.5)
end

% plot respons curves
plot(c_vec,sigmoid_array(:,mid_ex),'Color','k','LineWidth',1)

xlabel('activator concentration (c)')
ylabel('transcription rate (r)')
set(gca,'FontSize',14)
set(gca,'Color',pboc)
ylim([-.05 1.05])
box on
nfig1.InvertHardcopy = 'off';
saveas(nfig1,[FigPath 'noise_mod_profiles.png'])
saveas(nfig1,[FigPath 'noise_mod_profiles.pdf'])

% plot info for each
% generate info heatmap
ms = 1.1*max(max(diff(sigmoid_array)));
mp = 1.1*max(1./noise_vec);
ms_vec = linspace(0,ms);
mp_vec = linspace(0,mp)';
info_array = ms_vec .* mp_vec;

n_levels = 10;
levels = linspace(0,max(info_array(:)),n_levels+1);


info_fig = figure;
cmap4 = flipud(brewermap(n_levels,'Spectral'));
colormap(cmap4)
hold on
% make contour plot
contour(repmat(ms_vec,length(mp_vec),1),repmat(mp_vec,1,length(ms_vec)),info_array,n_levels,'Linewidth',2,'ShowText','Off')

for i = 1:n_examples
  scatter(s_array(i,1),n_array(i,1),100,'s','MarkerFaceColor',cmap1(i,:),'MarkerEdgeColor','k');
  scatter(s_array(i,2),n_array(i,2),75,'MarkerFaceColor',cmap2(i,:),'MarkerEdgeColor','k');  
end

xlabel('sensitivity (dr/dc)')
ylabel('precision (\sigma_{int}^{-1})')
set(gca,'FontSize',14)
set(gca,'Color',pboc)
box on
info_fig.InvertHardcopy = 'off';
saveas(info_fig,[FigPath 'info_contours.png'])
saveas(info_fig,[FigPath 'info_contours.pdf'])
