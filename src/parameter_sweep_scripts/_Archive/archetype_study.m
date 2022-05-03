clear
close all
load('C:\Users\nlamm\projects\noneq-transcription\out\bivariate_parameter_sweeps_v3\param_sweep_results_Sharpness_Precision.mat');

%% extract useful arrays
simulation_results_sv = simulation_results;
metric_array = simulation_results_sv.metric_array;
rate_array = simulation_results_sv.rate_array;

% find minimal variance network
[max_p,mi_p] = max(metric_array(:,4));
[max_s,mi_s] = max(metric_array(:,1));
precise_rates = rate_array(mi_p,:);
sharp_rates = rate_array(mi_s,:);
c_val_precise = fourStateHalfMaxGeneral(precise_rates);
c_val_sharp = fourStateHalfMaxGeneral(sharp_rates);

% define rate matrix 
R_precise =     [-precise_rates(3)-precise_rates(end)       c_val_precise*precise_rates(2)               0                            precise_rates(5); 
             precise_rates(end)            -c_val_precise*precise_rates(2)-precise_rates(7)     precise_rates(1)                           0
              0                              precise_rates(7)     -precise_rates(1)-c_val_precise*precise_rates(6)          precise_rates(4)
             precise_rates(3)                          0                    c_val_precise*precise_rates(6)        -precise_rates(4)-precise_rates(5) ];
         
         
R_sharp =     [-sharp_rates(3)-sharp_rates(end)       c_val_sharp*sharp_rates(2)               0                            sharp_rates(5); 
             sharp_rates(end)            -c_val_sharp*sharp_rates(2)-sharp_rates(7)     sharp_rates(1)                           0
              0                              sharp_rates(7)     -sharp_rates(1)-c_val_sharp*sharp_rates(6)          sharp_rates(4)
             sharp_rates(3)                          0                    c_val_sharp*sharp_rates(6)        -sharp_rates(4)-sharp_rates(5) ];         
         
         
% calculate predicted production rates and variance
p_rate_predicted = fourStateProductionGeneral(precise_rates,c_val_precise);
p_precision_predicted = 1/sqrt(fourStateVarianceGeneral(precise_rates,c_val_precise))
p_sharp_predicted = fourStateSharpnessGeneral(precise_rates,c_val_precise)*c_val_precise

s_rate_predicted = fourStateProductionGeneral(sharp_rates,c_val_sharp);
s_precision_predicted = 1/sqrt(fourStateVarianceGeneral(sharp_rates,c_val_sharp))
s_sharp_predicted = fourStateSharpnessGeneral(sharp_rates,c_val_sharp)*c_val_sharp

% calculate full SS vec
[V,D] = eig(R_precise);
[~,mi] = max(real(diag(D)));
ss_precise = V(:,mi)/sum(V(:,mi))      


[V,D] = eig(R_sharp);
[~,mi] = max(real(diag(D)));
ss_sharp = V(:,mi)/sum(V(:,mi)) 

%%

load('C:\Users\nlamm\projects\noneq-transcription\out\bivariate_parameter_sweeps_v3\param_sweep_results_Flux_Information.mat');

simulation_results_info = simulation_results;
metric_array = simulation_results_info.metric_array;
rate_array = simulation_results_info.rate_array;

% find minimal variance network
[max_i,mi_i] = max(metric_array(:,5));
info_rates = rate_array(mi_i,:);
c_val_info = 1;
% c_val_info = fourStateHalfMaxGeneral(info_rates);

% define rate matrix 
R_info =     [-info_rates(3)-info_rates(end)       c_val_info*info_rates(2)               0                            info_rates(5); 
             info_rates(end)            -c_val_info*info_rates(2)-info_rates(7)     info_rates(1)                           0
              0                              info_rates(7)     -info_rates(1)-c_val_info*info_rates(6)          info_rates(4)
             info_rates(3)                          0                    c_val_info*info_rates(6)        -info_rates(4)-info_rates(5) ];
         
         
% calculate predicted production rates and variance
i_rate_predicted = fourStateProductionGeneral(info_rates,c_val_info)
i_precision_predicted = 1/sqrt(fourStateVarianceGeneral(info_rates,c_val_info))
i_sharp_predicted = fourStateSharpnessGeneral(info_rates,c_val_info)*c_val_info


% calculate full SS vec
[V,D] = eig(R_info);
[~,mi] = max(real(diag(D)));
ss_info = V(:,mi)/sum(V(:,mi))      

