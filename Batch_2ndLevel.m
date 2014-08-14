function Batch_2ndLevel_v2_20140722()
% function to run 2nd level with batch
% written by Hongsheng-Cheng | hscheng.ok@gmail.com
% created date: 20140721
clear;clc;

disp('Any question, Please contact me at hscheng.ok@gmail.com');
fprintf('Welcome to our lab:  <a href="http://www.qiujlab.com/">Qiu Jiang Lab</a>\nhttp://www.qiujlab.com/\n');

% make sure your data img have the same demension as mask
[mask_filename,mask_filepath] = uigetfile('*.img;*.nii','Select Explict Mask');
explit_mask_loc = fullfile(mask_filepath,mask_filename);

GrpNum = input('How many groups in your data [1 or 2] : ');
switch GrpNum
    case 1
        stat_choose = input('Stat Choose: 1-1sample-ttest;2-multireg : ');
        if stat_choose == 1
            GrpName1 = input('Input Group Name : ','s');
            scan_grp1 = ScanImg2Cell;
            output_dir = fullfile(pwd,['1Sample_ttest_',GrpName1]);
            mkdir(output_dir)
            stat_1sample(output_dir,scan_grp1,explit_mask_loc,GrpName1)
        else
            % read behavior data
            [labels,datamat] = ReadXlsData();
            [cov,consess] = gen_cov_convec(labels,datamat);
            
            grp_img_name = input('Input group img name : ','s');
            ImgScans{1,1} = ScanImg2Cell;
            output_dir = fullfile(pwd,['Multireg_',grp_img_name,'_',labels{1}]);
            mkdir(output_dir);
            OutputDirs{1,1} = output_dir;
            Mark_Runtime = 2;
            % read multiple scans of img
            Mark2Run = 1;
            while (Mark2Run == 1) && (Mark_Runtime > 1)
                choice_add = questdlg('Would you like to add more file?', 'Questdlg','Yes','No','No');
                switch choice_add
                    case 'Yes' 
                        grp_img_name = input('Input group img name : ','s');
                        ImgScans{Mark_Runtime,1} = ScanImg2Cell;
                        output_dir = fullfile(pwd,['Multireg_',grp_img_name,'_',labels{1}]);
                        mkdir(output_dir);
                        OutputDirs{Mark_Runtime,1} = output_dir;
                        Mark_Runtime = Mark_Runtime+1;
                    case 'No'
                    Mark2Run = 0;
                end
            end
            
            grp_num = size(ImgScans,1);
            for kk = 1:grp_num
                stat_multireg(OutputDirs{kk,1},ImgScans{kk,1},explit_mask_loc,cov,consess)
            end
        end
   
    case 2
        % scan img file list in 2 groups
        GrpName1 = input('Input Group A Name : ','s');
        scan_grp1 = ScanImg2Cell;
        GrpName2 = input('Input Group B Name : ','s');
        scan_grp2 = ScanImg2Cell;
        
        output_dir_grp1 = fullfile(pwd,['1Sample_ttest_',GrpName1]);
        mkdir(output_dir_grp1)
        stat_1sample(output_dir_grp1,scan_grp1,explit_mask_loc,GrpName1)
        
        output_dir_grp2 = fullfile(pwd,['1Sample_ttest_',GrpName2]);
        mkdir(output_dir_grp2)
        stat_1sample(output_dir_grp2,scan_grp2,explit_mask_loc,GrpName2)
        
        output_dir_merge = fullfile(pwd,'1Sample_ttest_Merge');
        mkdir(output_dir_merge)
        scan_merge = [scan_grp1;scan_grp2];
        stat_1sample(output_dir_merge,scan_merge,explit_mask_loc,GrpName2)
        
        output_dir_2sample = fullfile(pwd,['2Sample_ttest_',GrpName1,'-',GrpName2]);
        mkdir(output_dir_2sample)
        stat_2sample(output_dir_2sample,scan_grp1,scan_grp2,explit_mask_loc,GrpName1,GrpName2)
    otherwise
        error('Wrong Input!');
end

end

function stat_multireg(output_dir,scan_grp1,explit_mask_loc,cov,consess)

matlabbatch{1, 1}.spm.stats.factorial_design.dir = {output_dir};
matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.scans = scan_grp1;
matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.mcov = cov;
% matlabbatch{1,1}.spm.stats.factorial_design.cov = cov;
matlabbatch{1, 1}.spm.stats.factorial_design.masking.em = {explit_mask_loc};

save(fullfile(output_dir,'step1_factorial_design.mat'),'matlabbatch');
spm_mat_list = dir(fullfile(output_dir,'SPM*.mat'));
if ~isempty(spm_mat_list)
    delete(fullfile(output_dir,spm_mat_list.name));
end

spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
clear matlabbatch;

% *************************** model estimation ****************************
% get the SPM.mat Subject Dir
spm_mat_list = dir(fullfile(output_dir,'SPM*.mat'));
spm_mat_path = fullfile(output_dir,spm_mat_list.name);

matlabbatch{1}.spm.stats.fmri_est.spmmat = {spm_mat_path};
matlabbatch{1}.spm.stats.fmri_est.method.Classical=1;

save(fullfile(output_dir,'step2_modelest.mat'),'matlabbatch');

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
clear matlabbatch;

% *************************** contrast manager ****************************

matlabbatch{1}.spm.stats.con.spmmat = {spm_mat_path};
matlabbatch{1, 1}.spm.stats.con.consess = consess;
matlabbatch{1}.spm.stats.con.delete = 0;

save(fullfile(output_dir,'step3_conrast.mat'),'matlabbatch');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
clear matlabbatch;

end

function stat_1sample(output_dir,scan_grp1,explit_mask_loc,grpname)

% ******************* factorial design specification **********************
matlabbatch{1, 1}.spm.stats.factorial_design.dir = {output_dir};
matlabbatch{1, 1}.spm.stats.factorial_design.des.t1.scans = scan_grp1;
% explicit mask
matlabbatch{1, 1}.spm.stats.factorial_design.masking.em = {explit_mask_loc};

save(fullfile(output_dir,'step1_factorial_design.mat'),'matlabbatch');
spm_mat_list = dir(fullfile(output_dir,'SPM*.mat'));
if ~isempty(spm_mat_list)
    delete(fullfile(output_dir,spm_mat_list.name));
end
spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
clear matlabbatch;
% *************************** model estimation ****************************
% get the SPM.mat Subject Dir
spm_mat_list = dir(fullfile(output_dir,'SPM*.mat'));
spm_mat_path = fullfile(output_dir,spm_mat_list.name);

matlabbatch{1}.spm.stats.fmri_est.spmmat = {spm_mat_path};
matlabbatch{1}.spm.stats.fmri_est.method.Classical=1;

save(fullfile(output_dir,'step2_modelest.mat'),'matlabbatch');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
clear matlabbatch;

% *************************** contrast manager ****************************

matlabbatch{1}.spm.stats.con.spmmat = {spm_mat_path};
matlabbatch{1, 1}.spm.stats.con.consess{1, 1}.tcon.name = grpname;
matlabbatch{1, 1}.spm.stats.con.consess{1, 1}.tcon.convec = 1;
matlabbatch{1}.spm.stats.con.delete = 0;

save(fullfile(output_dir,'step3_conrast.mat'),'matlabbatch');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
clear matlabbatch;
end

function stat_2sample(output_dir,scan_grp1,scan_grp2,explit_mask_loc,grpname1,grpname2)

% ******************* factorial design specification **********************
matlabbatch{1, 1}.spm.stats.factorial_design.dir = {output_dir};
matlabbatch{1, 1}.spm.stats.factorial_design.des.t2.scans1 = scan_grp1;
matlabbatch{1, 1}.spm.stats.factorial_design.des.t2.scans2 = scan_grp2;
% explicit mask
matlabbatch{1, 1}.spm.stats.factorial_design.masking.em = {explit_mask_loc};

save(fullfile(output_dir,'step1_factorial_design.mat'),'matlabbatch');
spm_mat_list = dir(fullfile(output_dir,'SPM*.mat'));
if ~isempty(spm_mat_list)
    delete(fullfile(output_dir,spm_mat_list.name));
end

spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
clear matlabbatch;
% *************************** model estimation ****************************
% get the SPM.mat Subject Dir
spm_mat_list = dir(fullfile(output_dir,'SPM*.mat'));
spm_mat_path = fullfile(output_dir,spm_mat_list.name);

matlabbatch{1}.spm.stats.fmri_est.spmmat = {spm_mat_path};
matlabbatch{1}.spm.stats.fmri_est.method.Classical=1;

save(fullfile(output_dir,'step2_modelest.mat'),'matlabbatch');

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
clear matlabbatch;

% *************************** contrast manager ****************************

matlabbatch{1}.spm.stats.con.spmmat = {spm_mat_path};
matlabbatch{1, 1}.spm.stats.con.consess{1, 1}.tcon.name = [grpname1,'-',grpname2];
matlabbatch{1, 1}.spm.stats.con.consess{1, 1}.tcon.convec = [1 -1];
matlabbatch{1, 1}.spm.stats.con.consess{1, 2}.tcon.name = [grpname2,'-',grpname1];
matlabbatch{1, 1}.spm.stats.con.consess{1, 2}.tcon.convec = [-1 1];
matlabbatch{1}.spm.stats.con.delete = 0;

save(fullfile(output_dir,'step3_contrast.mat'),'matlabbatch');

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
clear matlabbatch;
end

function img_list_cell = ScanImg2Cell()
% scan img files and reture the cell format list
    img_list = spm_select(Inf,'image','Select Img Files','',pwd,'.*',1);
    img_num = size(img_list,1);
    img_list_cell = cell(img_num,1);
    for ii = 1:img_num
        img_list_cell{ii} = [img_list(ii,:),',1'];
    end
end

function [Labels,DataMat] = ReadXlsData()
    prompt = {'Enter matrix :','Enter Labels :'};
    dlg_title = 'Input';
    num_lines = 1;
    Paras = inputdlg(prompt,dlg_title,num_lines);

    DataRaw = char(Paras{1,1});
    DataStr = regexp(DataRaw,' ','split');
    RowNum = size(DataStr,2);
    LabelTemp = char(Paras{2,1});
    Labels = regexp(LabelTemp,'\t','split');
    ColNum = size(Labels,2);

    DataMat = zeros(RowNum,ColNum);
    for ii = 1:RowNum
        DataTempStr = regexp(DataStr{ii},'\t','split');
        for jj = 1:ColNum
            DataMat(ii,jj) = cellfun(@str2num,DataTempStr(jj));
        end
    end
end

function [cov,consess] = gen_cov_convec(Labels,DataMat)
%     [Labels,DataMat] = ReadXlsData();
    var_num = length(Labels);
    
    % gen cov struct
    for jj = 1:var_num;
        cov(1,jj).c = DataMat(:,jj);
        cov(1,jj).cname = Labels{jj};
        cov(1,jj).iCFI = 1;
        cov(1,jj).iCC = 1;
        
        % initial a null convec
        % the first column is mean
        convec_null = zeros(1,var_num+1);
        consess{1,jj}.tcon.name = Labels{jj};
        convec_null(jj+1) = 1;
        consess{1,jj}.tcon.convec = convec_null;
        
    end
end
