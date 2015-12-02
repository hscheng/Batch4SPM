function GenVOI(varargin)
% GenVOI('sphere') 
% 
% written by hongshengcheng.math@gmail.com
% created date 20141201


    if numel(varargin) == 1
        voi_type = varargin{1};
    else
        voi_type = questdlg('voi type','Quest','sphere','thresh','sphere');
    end
    
    switch voi_type
        case 'sphere'
            gen_voi_sphere;
        case 'thresh'
            gen_voi_thresh;
    end
end
function voi_mat_path = gen_voi_sphere(spm_mat_path,fcon_idx,voi_name,sphere_center,sphere_radius)
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
    
%     matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {''};
%     matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = con_idx; % 3;
%     matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none'; % 'FWE';
%     matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.05;
%     matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
%     matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = sphere_center; % [15 -78 -9];
%     matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = sphere_radius; % 6;
%     matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.local.spm = 1; %fixed
%     matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
    
    % save matlabbatch
    voi_mat_list = dir(fullfile(sub_dir,['voi_',voi_name,'*.mat']));
    if ~isempty(voi_mat_list)
        mkdir(sub_dir,'old_backup')
        time_mark = datestr(clock,30);
        voi_mat_file = voi_mat_list.name;
        [~,mat_name] = fileparts(voi_mat_file);
        
        source_dir = fullfile(sub_dir,voi_mat_file);
        target_dir = fullfile(sub_dir,'old_backup',[mat_name,'_backup_',time_mark,'.mat']);
        movefile(source_dir,target_dir);
    end
    save(fullfile(sub_dir,['voi_',voi_name,'.mat']),'matlabbatch');
    % run batch
%     spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
end

function voi_mat_path = gen_voi_thresh(spm_mat_path,fcon_idx,voi_name,con_idx)
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
    
    matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {''};
    matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = con_idx; % 3;
    matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none'; % 'FWE';
    matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.05;
    matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
    
    % save matlabbatch
    voi_mat_list = dir(fullfile(sub_dir,['voi_',voi_name,'*.mat']));
    if ~isempty(voi_mat_list)
        mkdir(sub_dir,'old_backup')
        time_mark = datestr(clock,30);
        voi_mat_file = voi_mat_list.name;
        [~,mat_name] = fileparts(voi_mat_file);
        
        source_dir = fullfile(sub_dir,voi_mat_file);
        target_dir = fullfile(sub_dir,'old_backup',[mat_name,'_backup_',time_mark,'.mat']);
        movefile(source_dir,target_dir);
    end
    save(fullfile(sub_dir,['voi_',voi_name,'.mat']),'matlabbatch');
    % run batch
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
end