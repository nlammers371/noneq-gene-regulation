clear 
close all

project = 'mHMMeve2_weka_inf_2018_05_07';
% File Paths
writePath = ['../../out/input_output/' project '/'];

Bcd_path = 'E:\Jonathan\Dropbox\ZeldaLizHernan\Data\';
mkdir(writePath);
% load trace-related data
% load([dataPath '/analysis_data_mRNA.mat'])
bcdFileList = dir([Bcd_path '*-BcdGFP-HisRFP']);
%%% Now generate Bcd gradient data
% ap vector
ap_vec_src = 1:1000;
ap_vec_out = 100:900;
ap_ft = ismember(ap_vec_src,ap_vec_out);
% data at native resolution
time_full = [];
ap_full = [];
bcd_full_norm = [];
id_full = [];
particle_id_full = [];
% data at interp res
time_interp = [];
ap_interp = [];
id_interp = [];
for i = 1:numel(bcdFileList)
    f_path = [Bcd_path bcdFileList(i).name];    
    if isfolder(f_path)
        load([f_path '/CompiledNuclei.mat'])
        load([f_path '/' bcdFileList(i).name '_lin.mat']);           
        RawTF = AllTracesVector;                
        % Generate AP vector
        ap_temp = [];
        time_temp = [];
        bcd_temp = []; 
        p_id_temp = [];
        schnitz_ref = [CompiledNuclei.schnitz];
        for j = 1:numel(schnitzcells)
            cn_ind = find(schnitz_ref==j);
            if ~isempty(cn_ind)
                apv = schnitzcells(j).APpos;
                bcdv = CompiledNuclei(cn_ind).FluoMax';
                s_frames = schnitzcells(j).frames;
                cn_frames = CompiledNuclei(cn_ind).Frames;
                s_frames(apv>max(ap_vec_out)/1000|apv<min(ap_vec_out)/1000) = NaN;
                cn_frames(isnan(bcdv)) = NaN;
                % time vec
                tv = ElapsedTime(cn_frames(ismember(cn_frames,s_frames)&cn_frames>=nc14)) - ElapsedTime(nc14);
                if numel(tv) < 2
                    continue
                end          
                apv = 1000*apv(ismember(s_frames,cn_frames)&s_frames>=nc14);
                bcdv = bcdv(ismember(cn_frames,s_frames)&cn_frames>=nc14);                
                % record
                ap_temp = [ap_temp apv];
                bcd_temp = [bcd_temp bcdv];                
                time_temp = [time_temp tv];
                p_id_temp = [p_id_temp i + repelem(j / 1e4,numel(apv))];
            end
        end            
        bcd_full_norm = [bcd_full_norm bcd_temp];
        id_full = [id_full repelem(i,numel(bcd_temp))];
        time_full = [time_full time_temp];
        ap_full = [ap_full ap_temp];
        particle_id_full = [particle_id_full p_id_temp];
    end    
end

bcd_table = array2table([particle_id_full' round(time_full',3) round(ap_full',2) bcd_full_norm'],'VariableNames',...
    {'nucleus_id', 'time', 'ap', 'bcd'});
writetable(bcd_table,'../dat/bcd_table.csv')