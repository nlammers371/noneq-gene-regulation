clear
close all

DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\Nonequilibrium\Nick\manuscript\Figure 1\';

%%%%%%%%%%%%%%%%%
% Sharpness panel
%%%%%%%%%%%%%%%%%

n_vec = [1 3 6 12];
c_vec = linspace(0,1,1e3);
response_array = NaN(length(c_vec),length(n_vec));

for n = 1:length(n_vec)
    response_array(:,n) = c_vec.^n_vec(n) ./ (0.5^n_vec(n) + c_vec.^n_vec(n));
end

sharpness_fig = figure;
hold on
cmap = brewermap(8,'greens');

for i = 1:length(n_vec)
    plot(c_vec,response_array(:,i),'Color',cmap(i+1,:),'LineWidth',3)
end    

set(gca,'FontSize',14)
set(gca, 'xtick', [],'ytick', [])
xlabel('input (C)')
ylabel('output (M)')
set(gca,'Color',[228,221,209]/255) 

ax = gca;
% ax.XAxis.MinorTickValues = 10.^(-5:.1:5);
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
sharpness_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
ylim([0 1.05])

saveas(sharpness_fig,[DropboxFolder 'sharpness_panel.pdf'])
saveas(sharpness_fig,[DropboxFolder 'sharpness_panel.png'])

%%%%%%%%%%%%%%%%%
%% noise panel
%%%%%%%%%%%%%%%%%

noise_vec = fliplr(linspace(.0,.25,4));
n = 3;
response_vec = c_vec.^n ./ (0.5^n + c_vec.^n);


noise_fig = figure;
hold on
cmap = brewermap(8,'reds');

for i = 1:length(noise_vec)
    noisy_response = response_vec + normrnd(0,noise_vec(i),1,length(response_vec));
    noisy_response(noisy_response<0) = 0;
    plot(c_vec,noisy_response,'Color',cmap(i+1,:),'LineWidth',2)
end   

set(gca,'FontSize',14)
set(gca, 'xtick', [],'ytick', [])
xlabel('input (C)')
ylabel('output (M)')
set(gca,'Color',[228,221,209]/255) 
ylim([0 1.05])
ax = gca;
% ax.XAxis.MinorTickValues = 10.^(-5:.1:5);
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
noise_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(noise_fig,[DropboxFolder 'noise_panel.pdf'])
saveas(noise_fig,[DropboxFolder 'noise_panel.png'])

%%%%%%%%%%%%%%%%%
%% specificity panel
%%%%%%%%%%%%%%%%%

response_array = NaN(length(c_vec),length(n_vec));
n = 3;
response_vec = c_vec.^n ./ (0.5^n + c_vec.^n);


spec_fig = figure;
hold on
cmap = brewermap(9,'blues');

plot(c_vec,response_vec,'Color',cmap(8,:),'LineWidth',3)
plot(c_vec,response_vec*.3,'Color',cmap(2,:),'LineWidth',3)

legend('right activator (c_r)','wrong activator (c_w)','Location','northwest')
set(gca,'FontSize',14)
set(gca, 'xtick', [],'ytick', [])
xlabel('input (C)')
ylabel('output (M)')
set(gca,'Color',[228,221,209]/255) 
ylim([0 1.05])
ax = gca;
% ax.XAxis.MinorTickValues = 10.^(-5:.1:5);
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
spec_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(spec_fig,[DropboxFolder 'specificity_panel.pdf'])
saveas(spec_fig,[DropboxFolder 'specificity_panel.png'])