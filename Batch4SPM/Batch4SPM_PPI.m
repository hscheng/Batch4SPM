function Batch4SPM_PPI()
%% function to run PPI in batch mode
% ppi analysis follow 4 steps:
% 1. first level (multiple sessions should be merged)
% 2. gen voi and extract timecource
% 3. gen ppi struct
% 4. gen ppi glm model
% 
% written by hongshengcheng.math@gmail.com
% created date:20141201

    clear,clc;
fprintf('Welcome to our lab:  <a href="http://www.qiujlab.com/">Qiu Jiang Lab</a>\nhttp://www.qiujlab.com/\n');
fprintf('<a href="http://www.creativitybrain.com/">Center for Creativity and Brain</a>\nhttp://www.creativitybrain.com/\n');
disp('Any question, Please contact me at hongshengcheng.math@gmail.com');
disp('****************************************************************');

    h = figure;
    set(h,'color','w','numbertitle','off','menubar','none', 'resize','off',...
        'name','PPI','position',[400 400 300 150]);%[left, bottom, width, height]

    Button1 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.7 0.7 0.2],'fontsize',10,'fontweight','bold',...
        'String','generate voi',...
        'TooltipString','generate voi mat and img files,then extract timecourse');
    Button2 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.4 0.7 0.2],'fontsize',10,'fontweight','bold',...
        'String','gen ppi struct',...
        'TooltipString','batch to gen ppi struct');
    Button3 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.1 0.7 0.2],'fontsize',10,'fontweight','bold',...
        'String','ppi glm stat',...
        'TooltipString','perform glm');

    set(Button1,'Callback',@batch_gen_voi_callback);
    set(Button2,'Callback',@batch_gen_ppi_struct_callback);
    set(Button3,'Callback',@batch_ppi_glm_callback);
end
% http://cn.mathworks.com/help/matlab/creating_guis/write-callbacks-using-the-programmatic-workflow.html

function batch_gen_voi_callback(hObject,callbackdata)
    close all;
    spm_mat_list = ScanImg2Cell('Select 1st level spm mat','mat');
    sub_num = length(spm_mat_list);
    
    % voi parametres
    prompt = {'voi name e.g. DLPFC','center coords e.g. [12 24 53]','radius of voi e.g.6','F Contrast Index'};
    dlg_title = 'voi paras';
    num_lines = [8 30;8 30;1 30;1 30];
    def = {'','','6','15'};
    paras = inputdlg(prompt,dlg_title,num_lines,def);

    voi_names_str = paras{1};
    voi_num = size(voi_names_str,1);
    voi_names = cell(voi_num,1);
    for kk = 1:voi_num
        str = voi_names_str(kk,:);
        str(isspace(str)) = [];
        voi_names{kk} = str;
    end
    
    voi_coords = str2num(paras{2});
    voi_radius = str2num(paras{3});
    fcon_idx = str2num(paras{4});

%% save data info for check
    GenVOI_DataInfo.spm_mat_list = spm_mat_list;
    GenVOI_DataInfo.voi_coord_list = voi_coords;
    GenVOI_DataInfo.voi_name_list = voi_names;
    GenVOI_DataInfo.voi_radius = voi_radius;
    GenVOI_DataInfo.fcon_idx = fcon_idx;
    time_mark = datestr(clock,30);
    save(['GenVOI_DataInfo_',time_mark],'GenVOI_DataInfo');
 
    % batch to gen voi ang gen ppi struct
    for ii  = 1:sub_num
        sub_spm_mat = spm_mat_list{ii};      
        if matlabpool('size') >1
            parfor jj = 1:voi_num
                gen_voi(sub_spm_mat,fcon_idx,voi_names{jj},voi_coords(jj,:),voi_radius);
            end
        else
            for jj = 1:voi_num
                gen_voi(sub_spm_mat,fcon_idx,voi_names{jj},voi_coords(jj,:),voi_radius);
            end
        end
    end

end

function batch_gen_ppi_struct_callback(hObject,callbackdata)
    close all;
    sub_dir_list = MergeDir;
    
%     voi_mat_list = ScanImg2Cell('Select VOI Mat Files','mat');
    para = inputdlg('input VOI mat name','Enter',[10 30],{''});
    ppi_temp = para{1};
    mat_num = size(ppi_temp,1);
    mat_list = cell(length(sub_dir_list),mat_num);
    for kk = 1:mat_num
        str = ppi_temp(kk,:);
        str(isspace(str)) = [];
        mat_list(:,kk) = strcat(sub_dir_list,'\',str,'.mat');
    end
    
    voi_mat_list = reshape(mat_list',numel(mat_list),1);
    
    
    con_weight = gen_conweight();
%     ppi_prefix = input('ppi prefix : ','s');
    voi_num = length(voi_mat_list);
    con_num = length(con_weight);
    % batch to gen voi ang gen ppi struct
    for ii = 1:con_num
        con_name = con_weight(1,ii).name;
        con_mat = con_weight(1,ii).mat;
        
        if matlabpool('size') >1
            parfor jj  = 1:voi_num
                voi_mat_path = voi_mat_list{jj};
                [sub_dir,voi_mat_name] = fileparts(voi_mat_path);
                spm_mat_path = fullfile(sub_dir,'SPM.mat');
                ppi_name = [voi_mat_name,'_',con_name];
                gen_ppi_struct(spm_mat_path,voi_mat_path,ppi_name,con_mat);
            end
        else
            for jj  = 1:voi_num
                voi_mat_path = voi_mat_list{jj};
                [sub_dir,voi_mat_name] = fileparts(voi_mat_path);
                spm_mat_path = fullfile(sub_dir,'SPM.mat');
                ppi_name = [voi_mat_name,'_',con_name];
                gen_ppi_struct(spm_mat_path,voi_mat_path,ppi_name,con_mat);
            end
        end
    end

end

function batch_ppi_glm_callback(hObject,callbackdata)
%     dir_1st = uigetdir(pwd,'Select 1st Level Data Dir');
%     dir_data = fileparts(dir_1st);
%     dir_output = input('ppi output dir name e.g. Stat_PPI :','s');
%     if isempty(dir_output)
%         dir_output = 'Stat_PPI';
%     end
%     dir_ppi = fullfile(dir_data,dir_output);
%     mkdir(dir_ppi);
    close all;
    dir_ppi = uigetdir(pwd,'Select PPI Output Dir');
    
    model_spec_list = ScanImg2Cell('Select step1 model spec mat','mat');
    
    para = inputdlg('input ppi name','Enter',[10 30],{''});
    ppi_temp = para{1};
    ppi_mat_num = size(ppi_temp,1);
    ppi_mat_list = cell(ppi_mat_num,1);
    for kk = 1:ppi_mat_num
        str = ppi_temp(kk,:);
        str(isspace(str)) = [];
        ppi_mat_list{kk} = str;
    end

    consess = gen_consess;
    
%     spm('defaults','fmri');
%     spm_jobman('initcfg');
    sub_num = length(model_spec_list);

    time_mark = datestr(clock,30);
    log_file = fullfile(dir_ppi,['log_',time_mark,'.txt']);

for jj = 1:ppi_mat_num
    ppi_mat_name = ppi_mat_list{jj};
    if matlabpool('size') >1
        parfor  ii = 1:sub_num
            sub_model_spec = model_spec_list{ii};
            sub_dir = fileparts(sub_model_spec);
            [~,sub_id] = fileparts(sub_dir);
            ppi_mat = fullfile(sub_dir,[ppi_mat_name,'.mat']); %PPI_sub01_ppi_2
            ppi_data = importdata(ppi_mat);
            sub_ppi_dir = fullfile(dir_ppi,sub_id,ppi_mat_name);
            if exist(sub_ppi_dir,'dir')
                if ~exist(fullfile(sub_ppi_dir,'spmT_0001.img'),'file')
                    delete(fullfile(sub_ppi_dir,'SPM.mat'));
                    ppi_glm_model(sub_model_spec,sub_ppi_dir,ppi_data,consess);
                end
            else
                mkdir(sub_ppi_dir);
                ppi_glm_model(sub_model_spec,sub_ppi_dir,ppi_data,consess);
            end
            log_str = sprintf('%s work done\n',sub_ppi_dir);
            WriteRunLog(log_file,log_str);
        end
    else
        for  ii = 1:sub_num
            sub_model_spec = model_spec_list{ii};
            sub_dir = fileparts(sub_model_spec);
            [~,sub_id] = fileparts(sub_dir);
            ppi_mat = fullfile(sub_dir,[ppi_mat_name,'.mat']); %PPI_sub01_ppi_2
            ppi_data = importdata(ppi_mat);
            sub_ppi_dir = fullfile(dir_ppi,sub_id,ppi_mat_name);
            if exist(sub_ppi_dir,'dir')
                if ~exist(fullfile(sub_ppi_dir,'spmT_0001.img'),'file')
                    delete(fullfile(sub_ppi_dir,'SPM.mat'));
                    ppi_glm_model(sub_model_spec,sub_ppi_dir,ppi_data,consess);
                end
            else
                mkdir(sub_ppi_dir);
                ppi_glm_model(sub_model_spec,sub_ppi_dir,ppi_data,consess);
            end
            log_str = sprintf('%s work done\n',sub_ppi_dir);
            WriteRunLog(log_file,log_str);
        end
    end % matlabpool
end % ppi_mat_num
end

function voi_mat_path = gen_voi(spm_mat_path,fcon_idx,voi_name,sphere_center,sphere_radius)
%% function to define singal voi
% voi_path = Sub_VOI(spm_mat_path,con_idx,voi_prefix,sphere_center,sphere_radius)
% voi_name e.g. 'V2'
% con_idx e.g. 3,means the 3rd contrast
% sphere_center e.g. [15 -78 -9]
% sphere_radius e.g. 6

    sub_dir = fileparts(spm_mat_path);
    voi_mat_path = fullfile(sub_dir,['VOI_',voi_name,'_1.mat']);
    
    matlabbatch{1}.spm.util.voi.spmmat = {spm_mat_path};
    matlabbatch{1}.spm.util.voi.adjust = fcon_idx; % index of f contrast,default is 1
    matlabbatch{1}.spm.util.voi.session = 1;
    matlabbatch{1}.spm.util.voi.name = voi_name; % 'V2';
    
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = sphere_center; % [15 -78 -9];
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = sphere_radius; % 6;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
    matlabbatch{1}.spm.util.voi.expression = 'i1';
       
    % save matlabbatch
    voi_mat_list = dir(fullfile(sub_dir,['batch_voi_',voi_name,'*.mat']));
    if ~isempty(voi_mat_list)
        mkdir(sub_dir,'old_backup')
        time_mark = datestr(clock,30);
        voi_mat_file = voi_mat_list.name;
        [~,mat_name] = fileparts(voi_mat_file);
        
        source_dir = fullfile(sub_dir,voi_mat_file);
        target_dir = fullfile(sub_dir,'old_backup',[mat_name,'_backup_',time_mark,'.mat']);
        movefile(source_dir,target_dir);
    end
    save(fullfile(sub_dir,['batch_voi_',voi_name,'.mat']),'matlabbatch');
    % run batch
%     spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
end

function ppi_path = gen_ppi_struct(spm_mat_path,voi_mat,ppi_name,con_weight)
%% function to generate ppi structure
    clear matlabbatch

    sub_dir = fileparts(spm_mat_path);
    ppi_path = fullfile(sub_dir,['PPI_',ppi_name,'.mat']);
    
    matlabbatch{1}.spm.stats.ppi.spmmat = {spm_mat_path};
    matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {voi_mat};
    matlabbatch{1}.spm.stats.ppi.type.ppi.u = con_weight; % [2 1 -1; 3 1 1];
    matlabbatch{1}.spm.stats.ppi.name = ppi_name; % 'V2x(Att-NoAtt)';
    matlabbatch{1}.spm.stats.ppi.disp = 0;

    % save matlabbatch
    ppi_mat_list = dir(fullfile(sub_dir,['batch_ppi_',ppi_name,'.mat']));
    if ~isempty(ppi_mat_list)
        mkdir(sub_dir,'old_backup')
        time_mark = datestr(clock,30);
        ppi_mat_file = ppi_mat_list.name;
        [~,mat_name] = fileparts(ppi_mat_file);
        
        source_dir = fullfile(sub_dir,ppi_mat_file);
        target_dir = fullfile(sub_dir,'old_backup',[mat_name,'_backup_',time_mark,'.mat']);
        movefile(source_dir,target_dir);
    end
    save(fullfile(sub_dir,['batch_ppi_',ppi_name,'.mat']),'matlabbatch');
    % run batch
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    
end

function ppi_glm_model(model_spec_path,ppi_output_dir,PPI,consess)
%% funcction to rumm ppi glm model
% model_spec_path,1st level model spec mat path
% ppi_output_dir fullpath of ppi output
% PPI,ppi data struct
% consess,which contain the contrast info

[~,sub_id] = fileparts(ppi_output_dir);

regress_data(1).name = 'PPI_ppi';regress_data(1).val = PPI.ppi;
regress_data(2).name = 'PPI_Y';regress_data(2).val = PPI.Y;
regress_data(3).name = 'PPI_P';regress_data(3).val = PPI.P;

%% load old matlabbatch in 1st level
matlabbatch = importdata(model_spec_path);%matlabbatch in work space

old_sess = matlabbatch{1}.spm.stats.fmri_spec.sess;
if ismember('regress',fieldnames(old_sess))
    old_reg_num = length(old_sess.regress);
    if old_reg_num>0
        regress_data(4:3+old_reg_num) = old_sess.regress;
    end
end
    
%% modify old model spec batch
matlabbatch{1}.spm.stats.fmri_spec.dir = {ppi_output_dir}; % change the output dir
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct([]); % delete old condtions,onsets durations
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = regress_data; % regress ppi data
% matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''}; % no multiple regress data

%% model estimation
% get the SPM.mat Subject Dir
spm_mat_path = fullfile(ppi_output_dir,'SPM.mat');

matlabbatch{2}.spm.stats.fmri_est.spmmat = {spm_mat_path};
matlabbatch{2}.spm.stats.fmri_est.method.Classical=1;

%% contrast
matlabbatch{3}.spm.stats.con.spmmat = {spm_mat_path};
matlabbatch{3}.spm.stats.con.consess = consess;
matlabbatch{3}.spm.stats.con.delete = 0;

%% save and run batch
ppi_stat_mat_list = dir(fullfile(ppi_output_dir,['batch_ppi_glm_',sub_id,'.mat']));
if ~isempty(ppi_stat_mat_list)
    mkdir(ppi_output_dir,'old_backup')
    time_mark = datestr(clock,30);
    ppi_stat_mat_file = ppi_stat_mat_list.name;
    [~,mat_name] = fileparts(ppi_stat_mat_file);

    source_dir = fullfile(ppi_output_dir,ppi_stat_mat_file);
    target_dir = fullfile(ppi_output_dir,'old_backup',[mat_name,'_backup_',time_mark,'.mat']);
    movefile(source_dir,target_dir);
end
save(fullfile(ppi_output_dir,['batch_ppi_glm_',sub_id,'.mat']),'matlabbatch');
% run batch
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

clear all;
end

function consess = gen_consess()
%% function to define multiple f contrast and t contrast
    % define input dlg
    prompt = {'Enter the name of contrast:',...
        'Enter the contrast matrix:'};
    dlg_title = 'Input parameters';
    numlines = 1;
    para_default = {'ppi','1 0 0 0'};% e.g.kron([1 0 0 0]',ones(90,1))

    % add more regeressors
    Mark2Run = 1;
    Runtime = 1;
    while (Mark2Run == 1)
        choice_add = questdlg('Add Contrast?', 'Questdlg','T_Contrast','F_Contrast','No','No');
        switch choice_add
            case 'T_Contrast' 
                para = inputdlg(prompt,dlg_title,numlines,para_default);
                consess{1,Runtime}.tcon.name = para{1};
                consess{1,Runtime}.tcon.convec = str2num(para{2});
                Runtime = Runtime+1;
            case 'F_Contrast'
                para = inputdlg(prompt,dlg_title,numlines,para_default);
                consess{1,Runtime}.fcon.name = para{1};
                consess{1,Runtime}.fcon.convec = {eval(para{2})};
                Runtime = Runtime+1;
            case 'No'
                Mark2Run = 0;
        end
    end
end

function con_weight = gen_conweight()
    prompt = {'con name','contrast weight mat'};
    dlg_title = 'Input';
    num_lines = 1;
    def = {'con1','[1 1 -1;2 1 1]'};
    
    % add more contrast
    Mark2Run = 1;
    Runtime = 1;
    while (Mark2Run == 1)
        choice_add = questdlg('Add Contrast?', 'Questdlg','Yes','No','No');
        switch choice_add
            case 'Yes' 
                paras = inputdlg(prompt,dlg_title,num_lines,def);
                con_weight(1,Runtime).name = paras{1};
                con_weight(1,Runtime).mat = eval(paras{2});
                Runtime = Runtime+1;
            case 'No'
                Mark2Run = 0;
        end
    end
end

function WriteRunLog(logfilename,loginfo)
%% function to write log file
% logfilename:log.txt
% loginfo:string which you want to write in log file

    if(exist(logfilename,'file')==0) 
       fid = fopen(logfilename,'w+'); % create
    else %
       fid = fopen(logfilename,'a'); % add and write
    end
%    fwrite(fid,loginfo); 
    fprintf(fid,'%s\n',loginfo);
    fclose(fid); 
end