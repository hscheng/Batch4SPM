function Modify_Contrast(varargin)
%% function to modify the modelspec and contrast mat and run batch
% Modify_Contrast(modelspec_loc,run_del_idx)
% e.g. Modify_Contrast('E:\Stat_1stLevel\sub50\Stat_step1_ModelSpec_sub50.mat',[2 3]);
% 
% Note:
% 1.this function is only support for the Batch_Preproc4SPM
% 2.if you need delete some run in your data(e.g. headmotion supass 3.)
% 
% written by Hongsheng-Cheng | hongshengcheng.math@gmail.com
% 
% create date:20140712

if numel(varargin) == 2
    batch_mat_loc = varargin{1};
    run_del_idx = varargin{2};
    [mat_path,mat_name,ext] = fileparts(batch_mat_loc);
    mat_fullname = [mat_name,ext];
else
    % read the modelspec mat file
    [mat_fullname,mat_path] = uigetfile('*.mat','Select Step1 mat file');
    batch_mat_loc = fullfile(mat_path,mat_fullname);
    run_del_idx = input('run index to delete e.g.[5 6] : ');

end

run_num_del = length(run_del_idx);
disp('The mat you want to modify is:');
fprintf('\n%s',batch_mat_loc);

% ************************* Backup Old Files ******************************
mkdir(mat_path,'oldbackup');
BackupDir = fullfile(mat_path,'oldbackup');

% check *.img files
img_path = fullfile(mat_path,'*.img');
if ~isempty(dir(img_path))
    movefile(img_path,BackupDir);
end
% check *.hdr files
hdr_path = fullfile(mat_path,'*.hdr');
if ~isempty(dir(hdr_path))
    movefile(hdr_path,BackupDir);
end
% check *.nii files
nii_path = fullfile(mat_path,'*.nii');
if ~isempty(dir(nii_path))
    movefile(nii_path,BackupDir);
end
% check SPM.mat files
spmmat_path = fullfile(mat_path,'SPM.mat');
if ~isempty(dir(spmmat_path))
    movefile(spmmat_path,BackupDir);
end
% *************************************************************************

load(batch_mat_loc);

% get the postfix of the mat file
postfix_temp = regexp(mat_fullname,'sub','split');

run_num = size(matlabbatch{1,1}.spm.stats.fmri_spec.sess,2);
if (run_num_del >= run_num)||(max(run_del_idx)>run_num)
    error('Something wrong in your data');
else
    matlabbatch{1,1}.spm.stats.fmri_spec.sess(:,run_del_idx) = [];
end


% modify contrst
tcon_num = size(matlabbatch{1,3}.spm.stats.con.consess,2);
for ii =1:tcon_num
    contrast_temp = matlabbatch{1,1}.spm.stats.con.consess{1,ii}.tcon.convec;
    contrast_len = length(contrast_temp)/run_num;
    mark_del = (run_num-run_num_del)*contrast_len+1;
    contrast_temp(mark_del:end) =[];
    matlabbatch{1,1}.spm.stats.con.consess{1,ii}.tcon.convec = contrast_temp;
    clear contrast_temp;

end
    save(fullfile(mat_path,[mat_name,'_modified.mat']),'matlabbatch');
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    clear matlabbatch
    clear;clc;
    disp('work done');
end