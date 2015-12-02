function Batch4SPM_2ndLevel()
%% function to run 2nd level
% support t-test,multireg,and 2+ group corr interaction
% written by hongshengcheng.math@gmail.com
%
% 20141125,add paired glm regression(voxel matched)
% 20141111,add anova,fix error in corrinteract,and more poweful
% 20141105,moduled functions,which will be more powerfull
% 20141026,add new function to calcu paired img list correlation
% 20141012,fix bug in function ScanImg2Cell,add a new function PairedCorr
% 20141002,multiple regession support multi-behavior and multi-imgs;also fix the ReadXlsData function
% 20140819,add paired ttest
% 20140817,only save 1 batch mat file in res dir
% created date: 20140721

clear all;clc;
fprintf('Welcome to our lab:  <a href="http://www.qiujlab.com/">Qiu Jiang Lab</a>\nhttp://www.qiujlab.com/\n');
disp('Any question, Please contact me at hongshengcheng.math@gmail.com');
disp('****************************************************************');


% make sure your data img have the same demension as mask
[mask_filename,mask_filepath] = uigetfile('*.img;*.nii','Select Explict Mask');
if mask_filename ~= 0
    explicit_mask = fullfile(mask_filepath,mask_filename);
else
    explicit_mask = '';
end
grp_mode = questdlg('Group Type?', 'Questdlg','MultiGroups','2Groups','PairedGrps','MultiGroups');

LuckyNum = input('Enter your Lucky Number : ');
tic;
if LuckyNum ~= 1111
    switch grp_mode
        case 'MultiGroups'
            BatchStat_MultiGrps(explicit_mask);
        case '2Groups'
            BatchStat_2Sample(explicit_mask);
        case 'PairedGrps'
            BatchStat_Paired2Sample(explicit_mask);
    end
else
    switch grp_mode
        case 'MultiGroups'
            BatchStat_MultiGrps_Sup(explicit_mask);
        case '2Groups'
            BatchStat_2SampleStat_Sup(explicit_mask);
        case 'PairedGrps'
            BatchStat_Paired2Sample_Sup(explicit_mask);
    end
    clear all;
end
toc;
end

function BatchStat_MultiGrps(explicit_mask)
%% funtion for multiple groups (also fit for 1 sample/group)
% 1.1 sample t-test
% 2.multiple regression
% 3,multiple groups corrlation interaction

% choose stat type
stat_choose = questdlg('Stat Choose', 'Questdlg',...
    'MultiGrps_1Sample_ttest','MultiGrps_Regression','MultiGrps_ANOVA',...
    'MultiGrps_1Sample_ttest');

switch stat_choose
    case 'MultiGrps_1Sample_ttest'
        % read multiple brain img list
        [ImgListNames,ImgListScans] = ReadMultiImgList();
        ImgListNum = size(ImgListScans,1);
        
        for ii = 1:ImgListNum
            output_dir = fullfile(pwd,['1Sample_ttest_',ImgListNames{ii,1}]);
            mkdir(output_dir)
            stat_1sample_ttest(output_dir,ImgListScans{ii,1},ImgListNames{ii,1},explicit_mask)
        end
    case 'MultiGrps_Regression'
        % read multiple brain img list
        [ImgListNames,ImgListScans] = ReadMultiImgList();
        ImgListNum = size(ImgListScans,1);
        
        % read behavior data
        [labels_iv,data_iv,labels_cv,data_cv] = ReadMultiBehData();
        iv_num = length(labels_iv);

        % each img grp and each behavior in iv(independent variable)
        for ii = 1:ImgListNum
            for jj = 1:iv_num
                [cov,consess] = gen_cov_convec(labels_iv(jj),data_iv(:,jj),labels_cv,data_cv);
                OutputDir = fullfile(pwd,['Multireg_',ImgListNames{ii,1},'_',labels_iv{jj}]);
                mkdir(OutputDir);
                stat_multireg(OutputDir,ImgListScans{ii,1},explicit_mask,cov,consess)
            end
        end
    case 'MultiGrps_ANOVA'
        [GrpNames,GrpScans] = ReadMultiImgList();
            
        % read behavior data
        [labels_iv,data_iv,labels_cv,data_cv] = ReadMultiBehData();
        iv_num = length(labels_iv);

        for jj = 1:iv_num
            cov = gen_cov_convec(labels_iv(jj),data_iv(:,jj),labels_cv,data_cv);
            OutputDir = fullfile(pwd,['OnewayANOVA_1Grp_',labels_iv{jj}]);
            mkdir(OutputDir);
            stat_1way_anova(OutputDir,GrpScans,GrpNames,explicit_mask,cov)
        end
end
end

function BatchStat_2Sample(explicit_mask)
%% funtion for 2 sample stat(not paired)
% 1.2 sample t-test
% 2.2 sample correlation interaction

name_grp1 = input('Input Group A Name : ','s');
name_grp2 = input('Input Group B Name : ','s');
GrpAIndex = input('Group A index e.g. 1:30 : ');

% read multiple img list
[ImgListNames,ImgListScans] = ReadMultiImgList();
ImgListNum = size(ImgListScans,1);

stat_choose = questdlg('Stat Choose', 'Questdlg','2Sample_ttest','CorrelationInteraction','2Sample_ttest');
switch stat_choose
    case '2Sample_ttest'
        all_test_mark = questdlg('Also do 1 sample ttest? ','Quest','Yes','No','No');
        for ii = 1:ImgListNum
            GrpImgList = ImgListScans{ii,1};
            scan_grp1 = GrpImgList(GrpAIndex,:);
            GrpBIndex = GrpAIndex(end)+1:length(GrpImgList);
            scan_grp2 = GrpImgList(GrpBIndex,:);

            if strcmp(all_test_mark,'Yes') % one sample ttest
                output_dir_grp1 = fullfile(pwd,['1Sample_ttest_',name_grp1,'_',ImgListNames{ii,1}]);
                mkdir(output_dir_grp1)
                stat_1sample_ttest(output_dir_grp1,scan_grp1,name_grp1,explicit_mask)

                output_dir_grp2 = fullfile(pwd,['1Sample_ttest_',name_grp2,'_',ImgListNames{ii,1}]);
                mkdir(output_dir_grp2)
                stat_1sample_ttest(output_dir_grp2,scan_grp2,name_grp2,explicit_mask)

                output_dir_merge = fullfile(pwd,['1Sample_ttest_Merge','_',ImgListNames{ii,1}]);
                mkdir(output_dir_merge)               
                scan_merge = GrpImgList; % scan_merge = [scan_grp1;scan_grp2];
                stat_1sample_ttest(output_dir_merge,scan_merge,name_grp2,explicit_mask)
            end

            output_dir_2sample = fullfile(pwd,['2Sample_ttest_',name_grp1,'-',name_grp2,'_',ImgListNames{ii,1}]);
            mkdir(output_dir_2sample)
            stat_2sample_ttest(output_dir_2sample,scan_grp1,scan_grp2,name_grp1,name_grp2,explicit_mask)
        end

    case  'CorrelationInteraction'
        % read behavior data
        [labels_iv,data_iv,labels_cv,data_cv] = ReadMultiBehData();
        iv_num = length(labels_iv);

        % each img grp and each behavior in iv(independent variable)
        for ii = 1:ImgListNum
            GrpImgList = ImgListScans{ii,1};
            scan_grp1 = GrpImgList(GrpAIndex,:);
            GrpBIndex = GrpAIndex(end)+1:length(GrpImgList);
            scan_grp2 = GrpImgList(GrpBIndex,:);

            for jj = 1:iv_num
                cov = gen_cov_convec(labels_iv(jj),data_iv(:,jj),labels_cv,data_cv);
                OutputDir = fullfile(pwd,['CorrInteract_',ImgListNames{ii,1},'_',labels_iv{jj}]);
                mkdir(OutputDir);
                stat_2sample_corr_interact(OutputDir,scan_grp1,scan_grp2,name_grp1,name_grp2,explicit_mask,cov)
            end
        end
end
end

function BatchStat_Paired2Sample(explicit_mask)
%% funtion for paired 2 sample stat
% 1.paired t-test
% 2.paired img correlation (voxel matched)
% 3.paired imgs general linear model regression(voxel matched)

    % scan img file list in 2 groups
    % make sure the 2 Img List have the same order
    name_grp1 = input('Input Group A Name : ','s');
    scan_grp1 = ScanImg2Cell('Select Group 1 Img List');
    name_grp2 = input('Input Group B Name : ','s');
    scan_grp2 = ScanImg2Cell('Select Group 2 Img List');
    
    stat_choose = questdlg('Stat Choose', 'Questdlg','Paired_ttest','Paired_Corr','Paired_GlmReg','Paired_ttest');

    switch stat_choose
        case 'Paired_ttest'
            OutputDir = fullfile(pwd,['Paired_ttest_',name_grp1,'-',name_grp2]);
            mkdir(OutputDir)
            stat_paired_ttest(OutputDir,scan_grp1,scan_grp2,name_grp1,name_grp2,explicit_mask) 
        case 'Paired_Corr'
            OutputDir = fullfile(pwd,['Paired_Corr_',name_grp1,'_',name_grp2]);
            mkdir(OutputDir);
            block_mark = input('block comput for very large mat? 1-Yes;0-No : ');
            stat_paired_corr(OutputDir,scan_grp1,scan_grp2,name_grp1,name_grp2,explicit_mask,block_mark)
        case 'Paired_GlmReg'
            tic;
            OutputDir = fullfile(pwd,['Paired_GlmReg_',name_grp1,'_',name_grp2]);
            mkdir(OutputDir);
            block_mark = input('block comput for very large mat? 1-Yes;0-No : ');
            stat_paired_glmreg(OutputDir,scan_grp1,scan_grp2,name_grp1,name_grp2,explicit_mask,block_mark);
            toc;
    end
end

function BatchStat_MultiGrps_Sup(explicit_mask)
%% funtion for multiple groups (also fit for 1 sample/group)
% 1.1 sample t-test
% 2.multiple regression
% 3.multiple groups corrlation interaction(F test)

% read multiple para_path,xls
% a para_path.xls should be in format:
% 1st column is the group name,2nd,3rd...is the parametre path
% 1st row is the name of parametres
% make sure all the paras have the same order

[GrpNames,ParaNames,ParaCell] = ReadMultiPathInfo;
% when there is only 1 grp in para_path.xls,then paracell is 1*n cell
% each cell is sub_num*1 cell
ParaNum = length(ParaNames);

stat_choose = questdlg('Stat Choose', 'Questdlg',...
    'MultiGrps_1Sample_ttest','MultiGrps_Regression','MultiGrps_ANOVA',...
    'MultiGrps_1Sample_ttest');
switch stat_choose
    case 'MultiGrps_1Sample_ttest'       
        for ii = 1:ParaNum
            output_dir = fullfile(pwd,['1Sample_ttest_',ParaNames{ii}]);
            mkdir(output_dir)
            stat_1sample_ttest(output_dir,ParaCell{:,ii},ParaNames{ii},explicit_mask);
        end
    case 'MultiGrps_Regression'
        % read behavior data
        [labels_iv,data_iv,labels_cv,data_cv] = ReadMultiBehData();
        iv_num = length(labels_iv);

        for ii = 1:ParaNum
            GrpScans = ParaCell(:,ii);
            for jj = 1:iv_num
                [cov,consess] = gen_cov_convec(labels_iv(jj),data_iv(:,jj),labels_cv,data_cv);
                OutputDir = fullfile(pwd,['Multireg_',ParaNames{ii},'_',labels_iv{jj}]);
                mkdir(OutputDir);
                stat_multireg(OutputDir,GrpScans{1},explicit_mask,cov,consess);
            end
        end
    case 'MultiGrps_ANOVA'
        % read behavior data
        [labels_iv,data_iv,labels_cv,data_cv] = ReadMultiBehData();
        iv_num = length(labels_iv);

        for ii = 1:ParaNum
            GrpScans = ParaCell(:,ii);
            for jj = 1:iv_num
                cov = gen_cov_convec(labels_iv(jj),data_iv(:,jj),labels_cv,data_cv);
                OutputDir = fullfile(pwd,['OnewayANOVA_',ParaNames{ii},'_',labels_iv{jj}]);
                mkdir(OutputDir);
                stat_1way_anova(OutputDir,GrpScans{1},GrpNames,explicit_mask,cov);
            end
        end
end
end

function BatchStat_2SampleStat_Sup(explicit_mask)
%% funtion for 2 sample stat(not paired)
% 1.2 sample t-test
% 2.2 sample correlation interaction

% read multiple para_path,xls
[GrpNames,ParaNames,ParaCell] = ReadMultiPathInfo;
% when there is only 1 grp in para_path.xls,then paracell is 1*n cell
% each cell is sub_num*1 cell
ParaNum = length(ParaNames);

name_grp1 = GrpNames{1};
name_grp2 = GrpNames{2};

stat_choose = questdlg('Stat Choose', 'Questdlg','2Sample_ttest','CorrelationInteraction','2Sample_ttest');
switch stat_choose
    case '2Sample_ttest'
        all_test_mark = questdlg('Also do 1 sample ttest? ','Quest','Yes','No','No');
        for ii = 1:ParaNum
            GrpScans = ParaCell(:,ii);
            GrpScan1 = GrpScans{1};
            GrpScan2 = GrpScans{2};
            if strcmp(all_test_mark,'Yes') % one sample ttest
                output_dir_grp1 = fullfile(pwd,['1Sample_ttest_',name_grp1,'_',ParaNames{ii}]);
                mkdir(output_dir_grp1)
                stat_1sample_ttest(output_dir_grp1,GrpScan1,name_grp1,explicit_mask)

                output_dir_grp2 = fullfile(pwd,['1Sample_ttest_',name_grp2,'_',ParaNames{ii}]);
                mkdir(output_dir_grp2)
                stat_1sample_ttest(output_dir_grp2,GrpScan2,name_grp2,explicit_mask)

%                 output_dir_merge = fullfile(pwd,['1Sample_ttest_Merge','_',ParaNames{ii}]);
%                 mkdir(output_dir_merge)               
%                 % scan_merge = GrpScans; 
%                 scan_merge = [GrpScan1;GrpScan2];
%                 stat_1sample_ttest(output_dir_merge,scan_merge,name_grp1,explicit_mask)
            end

            output_dir_2sample = fullfile(pwd,['2Sample_ttest_',name_grp1,'_',name_grp2,'_',ParaNames{ii}]);
            mkdir(output_dir_2sample)
            stat_2sample_ttest(output_dir_2sample,GrpScan1,GrpScan2,name_grp1,name_grp2,explicit_mask)
        end

    case  'CorrelationInteraction'
        % read behavior data
        [labels_iv,data_iv,labels_cv,data_cv] = ReadMultiBehData();
        iv_num = length(labels_iv);

        % each img grp and each behavior in iv(independent variable)
        for ii = 1:ParaNum
            GrpScans = ParaCell(:,ii);
            GrpScan1 = GrpScans{1};
            GrpScan2 = GrpScans{2};

            for jj = 1:iv_num
                cov = gen_cov_convec(labels_iv(jj),data_iv(:,jj),labels_cv,data_cv);
                output_dir = fullfile(pwd,['CorrInteract_',ParaNames{ii},'_',labels_iv{jj}]);
                mkdir(output_dir);
                stat_2sample_corr_interact(output_dir,GrpScan1,GrpScan2,name_grp1,name_grp2,explicit_mask,cov)
            end
        end
end
end

function BatchStat_Paired2Sample_Sup(explicit_mask)
%% funtion for paired 2 sample stat
% 1.paired t-test
% 2.paired img correlation (voxel matched)
% 3.paired imgs general linear model regression(voxel matched)

% read multiple para_path,xls
[GrpNames,ParaNames,ParaCell] = ReadMultiPathInfo;
% when there is only 1 grp in para_path.xls,then paracell is 1*n cell
% each cell is sub_num*1 cell
ParaNum = length(ParaNames);

name_grp1 = GrpNames{1};
name_grp2 = GrpNames{2};

stat_choose = questdlg('Stat Choose', 'Questdlg','Paired_t-test','Paired_Corr','Paired_GlmReg','Paired_t-test');
switch stat_choose
    case 'Paired_t-test'
        for ii = 1:ParaNum
            GrpScans = ParaCell(:,ii);
            OutputDir = fullfile(pwd,['Paired_ttest_',ParaNames{ii},'_',name_grp1,'-',name_grp2]);
            mkdir(OutputDir)
            stat_paired_ttest(OutputDir,GrpScans{1},GrpScans{2},name_grp1,name_grp2,explicit_mask) 
        end
    case 'Paired_Corr'
        block_mark = input('block comput for very large mat? 1-Yes;0-No : ');
        for ii = 1:ParaNum
            GrpScans = ParaCell(:,ii);
            OutputDir = ['Paired_Corr_',ParaNames{ii},'_',name_grp1,'_',name_grp2];
            mkdir(OutputDir);
            stat_paired_corr(OutputDir,GrpScans{1},GrpScans{2},name_grp1,name_grp2,explicit_mask,block_mark)
        end
    case 'Paired_GlmReg'
        block_mark = input('block comput for very large mat? 1-Yes;0-No : ');
        for ii = 1:ParaNum
            GrpScans = ParaCell(:,ii);
            tic;
            OutputDir = fullfile(pwd,['Paired_GlmReg_',ParaNames{ii},'_',name_grp1,'_',name_grp2]);
            mkdir(OutputDir);
            stat_paired_glmreg(OutputDir,GrpScans{1},GrpScans{2},name_grp1,name_grp2,explicit_mask,block_mark);
            toc;
        end
end
end

function stat_multireg(output_dir,scan_grp,explicit_mask,cov,consess)
%% function to run multiple regression
    matlabbatch{1,1}.spm.stats.factorial_design.dir = {output_dir};
    matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.scans = scan_grp;
    matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.mcov = cov;
    % matlabbatch{1,1}.spm.stats.factorial_design.cov = cov;
    matlabbatch{1,1}.spm.stats.factorial_design.masking.em = {explicit_mask};

    spm_mat_list = dir(fullfile(output_dir,'SPM.mat'));
    if ~isempty(spm_mat_list)
        delete(fullfile(output_dir,spm_mat_list.name));
    end

    % model estimation
    % get the SPM.mat Subject Dir
    spm_mat_path = fullfile(output_dir,'SPM.mat');

    matlabbatch{1,2}.spm.stats.fmri_est.spmmat = {spm_mat_path};
    matlabbatch{1,2}.spm.stats.fmri_est.method.Classical=1;

    %  contrast manager
    matlabbatch{1,3}.spm.stats.con.spmmat = {spm_mat_path};
    matlabbatch{1,3}.spm.stats.con.consess = consess;
    matlabbatch{1,3}.spm.stats.con.delete = 0;

    % save and run batch
    save(fullfile(output_dir,'batch_multi_glm_reg.mat'),'matlabbatch');
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    clear all;%matlabbatch;
end

function stat_1sample_ttest(output_dir,scan_grp,name_grp,explicit_mask)
%% function to run One sample ttest

    % factorial design specification
    matlabbatch{1,1}.spm.stats.factorial_design.dir = {output_dir};
    matlabbatch{1,1}.spm.stats.factorial_design.des.t1.scans = scan_grp;
    % explicit mask
    matlabbatch{1,1}.spm.stats.factorial_design.masking.em = {explicit_mask};

    spm_mat_list = dir(fullfile(output_dir,'SPM.mat'));
    if ~isempty(spm_mat_list)
        delete(fullfile(output_dir,spm_mat_list.name));
    end

    % Model estimation
    % get the SPM.mat Subject Dir
    spm_mat_path = fullfile(output_dir,'SPM.mat');

    matlabbatch{1,2}.spm.stats.fmri_est.spmmat = {spm_mat_path};
    matlabbatch{1,2}.spm.stats.fmri_est.method.Classical=1;

    % Contrast manager
    matlabbatch{1,3}.spm.stats.con.spmmat = {spm_mat_path};
    matlabbatch{1,3}.spm.stats.con.consess{1,1}.tcon.name = name_grp;
    matlabbatch{1,3}.spm.stats.con.consess{1,1}.tcon.convec = 1;
    matlabbatch{1,3}.spm.stats.con.delete = 0;

    % save and run batch
    save(fullfile(output_dir,'batch_1sample_ttest.mat'),'matlabbatch');
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    clear all;%matlabbatch;
end

function stat_2sample_ttest(output_dir,scan_grp1,scan_grp2,name_grp1,name_grp2,explicit_mask)
%% function for 2 sample ttest

    % factorial design specification
    matlabbatch{1,1}.spm.stats.factorial_design.dir = {output_dir};
    matlabbatch{1,1}.spm.stats.factorial_design.des.t2.scans1 = scan_grp1;
    matlabbatch{1,1}.spm.stats.factorial_design.des.t2.scans2 = scan_grp2;
    % explicit mask
    matlabbatch{1,1}.spm.stats.factorial_design.masking.em = {explicit_mask};

    spm_mat_list = dir(fullfile(output_dir,'SPM.mat'));
    if ~isempty(spm_mat_list)
        delete(fullfile(output_dir,spm_mat_list.name));
    end

    % Model estimation
    % get the SPM.mat Subject Dir
    spm_mat_path = fullfile(output_dir,'SPM.mat');

    matlabbatch{1,2}.spm.stats.fmri_est.spmmat = {spm_mat_path};
    matlabbatch{1,2}.spm.stats.fmri_est.method.Classical=1;

    % Contrast manager
    matlabbatch{1,3}.spm.stats.con.spmmat = {spm_mat_path};
    matlabbatch{1,3}.spm.stats.con.consess{1,1}.tcon.name = [name_grp1,'-',name_grp2];
    matlabbatch{1,3}.spm.stats.con.consess{1,1}.tcon.convec = [1 -1];
    matlabbatch{1,3}.spm.stats.con.consess{1,2}.tcon.name = [name_grp2,'-',name_grp1];
    matlabbatch{1,3}.spm.stats.con.consess{1,2}.tcon.convec = [-1 1];
    matlabbatch{1,3}.spm.stats.con.delete = 0;

    % save and run batch
    save(fullfile(output_dir,'batch_2sample_ttest.mat'),'matlabbatch');
    
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    clear all;%matlabbatch;
end

function stat_2sample_corr_interact(output_dir,scan_grp1,scan_grp2,name_grp1,name_grp2,explicit_mask,cov)
%% function to run 2 sample correlation interacion
% when covariate have total brain volumn the matlabbatch should be fixed

    matlabbatch{1,1}.spm.stats.factorial_design.dir = {output_dir};
    matlabbatch{1,1}.spm.stats.factorial_design.des.anova.icell(1,1).scans = scan_grp1;
    matlabbatch{1,1}.spm.stats.factorial_design.des.anova.icell(1,2).scans = scan_grp2;
    matlabbatch{1,1}.spm.stats.factorial_design.cov = cov;
    matlabbatch{1,1}.spm.stats.factorial_design.masking.em = {explicit_mask};
    % 20141111,error fixed
    for jj = 1:size(cov,2)
        %in corr interact analysis,TBV should not interact with factor 1
        if strcmp(cov(1,jj).cname,'TBV')~=1;%TBV for total brain volunme
            matlabbatch{1,1}.spm.stats.factorial_design.cov(1,jj).iCFI = 2;%interact with factor 1
        end
    end
    spm_mat_list = dir(fullfile(output_dir,'SPM.mat'));
    if ~isempty(spm_mat_list)
        delete(fullfile(output_dir,spm_mat_list.name));
    end

    % model estimation
    % get the SPM.mat Subject Dir
    spm_mat_path = fullfile(output_dir,'SPM.mat');

    matlabbatch{1,2}.spm.stats.fmri_est.spmmat = {spm_mat_path};
    matlabbatch{1,2}.spm.stats.fmri_est.method.Classical=1;

    %  contrast manager
    matlabbatch{1,3}.spm.stats.con.spmmat = {spm_mat_path};
    matlabbatch{1,3}.spm.stats.con.consess{1,1}.tcon.name = [name_grp1,'-',name_grp2];
    matlabbatch{1,3}.spm.stats.con.consess{1,1}.tcon.convec = [0 0 1 -1];
    matlabbatch{1,3}.spm.stats.con.consess{1,2}.tcon.name = [name_grp2,'-',name_grp1];
    matlabbatch{1,3}.spm.stats.con.consess{1,2}.tcon.convec = [0 0 -1 1];
    matlabbatch{1,3}.spm.stats.con.delete = 0;

    % save and run batch
    save(fullfile(output_dir,'batch_2sample_corr_interact.mat'),'matlabbatch');

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    clear all;%matlabbatch;

end

function stat_1way_anova(output_dir,grp_scans,grp_names,explicit_mask,cov)
%% function for multiple grps anova
% factorial design specification
% when covariate have total brain volumn the matlabbatch should be fixed

    grp_num = length(grp_scans);%grp_scans is cell in cell format
    
    matlabbatch{1,1}.spm.stats.factorial_design.dir = {output_dir};
    
    for ii = 1:grp_num
        matlabbatch{1,1}.spm.stats.factorial_design.des.anova.icell(1,ii).scans = grp_scans{ii};
    end
    
    matlabbatch{1,1}.spm.stats.factorial_design.cov = cov;
    % 20141111,error fixed
    for jj = 1:size(cov,2)
        %in corr interact analysis,TBV should not interact with factor 1
        if strcmp(cov(1,jj).cname,'TBV')~=1;%TBV for total brain volunme
            matlabbatch{1,1}.spm.stats.factorial_design.cov(1,jj).iCFI = 2;%interact with factor 1
        end
    end
    
    % explicit mask
    matlabbatch{1,1}.spm.stats.factorial_design.masking.em = {explicit_mask};
    
    spm_mat_list = dir(fullfile(output_dir,'SPM.mat'));
    if ~isempty(spm_mat_list)
        delete(fullfile(output_dir,spm_mat_list.name));
    end

    % Model estimation
    % get the SPM.mat Subject Dir
    spm_mat_path = fullfile(output_dir,'SPM.mat');

    matlabbatch{1,2}.spm.stats.fmri_est.spmmat = {spm_mat_path};
    matlabbatch{1,2}.spm.stats.fmri_est.method.Classical=1;

    % Contrast manager
    matlabbatch{1,3}.spm.stats.con.spmmat = {spm_mat_path};
    
    for ii = 1:length(grp_scans)
        matlabbatch{1,3}.spm.stats.con.consess{1,ii}.fcon.name = grp_names{ii};
        diag_mat = -ones(grp_num,1)/(grp_num-1);%contrast parameter sums equal to zero
        diag_mat(ii) = 1;
        con_mat = [zeros(grp_num) diag(diag_mat)];
        matlabbatch{1,3}.spm.stats.con.consess{1,ii}.fcon.convec = {con_mat};
    end
    matlabbatch{1,3}.spm.stats.con.delete = 0;

    % save and run batch
    save(fullfile(output_dir,'batch_anova_multigrps.mat'),'matlabbatch');

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    clear all;%matlabbatch;
end

function stat_paired_ttest(output_dir,scan_grp1,scan_grp2,name_grp1,name_grp2,explicit_mask)
%% funtion for paired ttest

    % factorial design specification
    matlabbatch{1,1}.spm.stats.factorial_design.dir = {output_dir};
    for ii = 1:length(scan_grp1)
        matlabbatch{1,1}.spm.stats.factorial_design.des.pt.pair(1,ii).scans{1,1} = scan_grp1{ii};
        matlabbatch{1,1}.spm.stats.factorial_design.des.pt.pair(1,ii).scans{2,1} = scan_grp2{ii};
    end

    % explicit mask
    matlabbatch{1,1}.spm.stats.factorial_design.masking.em = {explicit_mask};

    spm_mat_list = dir(fullfile(output_dir,'SPM.mat'));
    if ~isempty(spm_mat_list)
        delete(fullfile(output_dir,spm_mat_list.name));
    end

    % Model estimation
    % get the SPM.mat Subject Dir
    spm_mat_path = fullfile(output_dir,'SPM.mat');

    matlabbatch{1,2}.spm.stats.fmri_est.spmmat = {spm_mat_path};
    matlabbatch{1,2}.spm.stats.fmri_est.method.Classical=1;

    % Contrast manager
    matlabbatch{1,3}.spm.stats.con.spmmat = {spm_mat_path};
    matlabbatch{1,3}.spm.stats.con.consess{1,1}.tcon.name = [name_grp1,'-',name_grp2];
    matlabbatch{1,3}.spm.stats.con.consess{1,1}.tcon.convec = [1 -1];
    matlabbatch{1,3}.spm.stats.con.consess{1,2}.tcon.name = [name_grp2,'-',name_grp1];
    matlabbatch{1,3}.spm.stats.con.consess{1,2}.tcon.convec = [-1 1];
    matlabbatch{1,3}.spm.stats.con.delete = 0;

    % save and run batch
    save(fullfile(output_dir,'batch_paired_ttest.mat'),'matlabbatch');
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    clear all;%matlabbatch;
end

function stat_paired_corr(output_dir,scan_grp1,scan_grp2,name_grp1,name_grp2,explicit_mask,block_mark)
%% function to perform paired imgs (voxel matched) correlation
% function to corr 1 group 2 kind img.e.g.reho alff
% make sure 2 groups have same img number

    DataStruct = struct();
    DataStruct = setfield(DataStruct,name_grp1,scan_grp1);
    DataStruct = setfield(DataStruct,name_grp2,scan_grp2);
    
    [Data_Grp1,MaskIdx] = ReadMultiImg(scan_grp1,explicit_mask);
    Data_Grp2 = ReadMultiImg(scan_grp2,explicit_mask);
    DataStruct.Data_Grp1 = Data_Grp1;
    DataStruct.Data_Grp2 = Data_Grp2;
    save(fullfile(output_dir,'DataStruct.mat'), 'DataStruct');
    
    % check data dimesion
    % PROD is only supported for floating point input.
    if prod(double(size(Data_Grp1)==size(Data_Grp2))) == 0
        error('2 ImgLists dimesion not same');
    end

    if isempty(explicit_mask)
        ArrayLen = size(Data_Grp1,2);
    else
        ArrayLen = length(MaskIdx);
    end
    

%     block_mark = input('Block Comput for Very large matrices? 1-Yes;0-No : ');
    
    if block_mark == 1
        % set block size for block-comput
        BlockSize = 5000;
        BlockNum = ceil(ArrayLen/BlockSize);

        % pre allocate memory
        r_merge = zeros(1,ArrayLen);
        p_merge = zeros(1,ArrayLen);

        for jj = 1:BlockNum
            disp([num2str(BlockNum-jj),' Blocks Left']);
            % get the array index
            if jj < BlockNum
                ArrayIdx = BlockSize*(jj-1)+1:BlockSize*jj;
            else
                ArrayIdx = BlockSize*(jj-1)+1:ArrayLen;
            end
            Temp1 = Data_Grp1(:,ArrayIdx);
            Temp2 = Data_Grp2(:,ArrayIdx);

            % comput general linear model regression paras
            % arrayfun faster than cellfun
%             [r,p] = cellfun(@corr,num2cell(Temp1,1),num2cell(Temp2,1));%SubNum*BlockSize
            [r,p] = arrayfun(@(k) corr(Temp1(:,k),Temp2(:,k)),1:length(ArrayIdx));
            r_merge(ArrayIdx) = r;
            p_merge(ArrayIdx) = p;
        end
    else
%         [r_merge,p_merge] = cellfun(@corr,num2cell(GrpData1,1),num2cell(GrpData2,1));
        [r_merge,p_merge] = arrayfun(@(k) corr(Data_Grp1(:,k),Data_Grp2(:,k)),1:ArrayLen);
    end
    
    V = spm_vol(scan_grp1{1});
    Y = spm_read_vols(V);
    % reshape the merged result mat for img-written
    if isempty(explicit_mask)
        r_last = reshape(r_merge,size(Y));r_last(isnan(r_last)) = 0;
        p_last = reshape(p_merge,size(Y));p_last(isnan(p_last)) = 0;
    else
        r_last = zeros(size(Y));r_last(MaskIdx) = r_merge;
        p_last = zeros(size(Y));p_last(MaskIdx) = p_merge;
    end
    % write res imgs

    V.fname = 'res_R_value.nii';spm_write_vol(V,r_last);
    V.fname = 'res_P_value.nii';spm_write_vol(V,p_last);
    p_last = 1/p_last;p_last(isinf(p_last))=0;
    V.fname = 'res_P_recip_val.nii';spm_write_vol(V,p_last);
    movefile('res*.nii',output_dir);
    clear all;

end

function stat_paired_glmreg(output_dir,scan_grp1,scan_grp2,name_grp1,name_grp2,explicit_mask,block_mark)
%% function to perform paired imgs(voxel matched) glm regression
% make sure 2 groups' imgs are matched

    DataStruct = struct();
    DataStruct = setfield(DataStruct,name_grp1,scan_grp1);
    DataStruct = setfield(DataStruct,name_grp2,scan_grp2);

    %GrpData is ImgNum*voxel_num
    [Data_Grp1,MaskIdx] = ReadMultiImg(scan_grp1,explicit_mask);
    Data_Grp2 = ReadMultiImg(scan_grp2,explicit_mask);
    DataStruct.Data_Grp1 = Data_Grp1;
    DataStruct.Data_Grp2 = Data_Grp2;
    save(fullfile(output_dir,'DataStruct.mat'), 'DataStruct');
    
	% PROD is only supported for floating point input.
    if prod(double(size(Data_Grp1)==size(Data_Grp2))) == 0
        error('2 ImgLists dimesion not same');
    end

    if isempty(explicit_mask)
        ArrayLen = size(Data_Grp1,2);
    else
        ArrayLen = length(MaskIdx);
    end
    
%     block_mark = input('Block Comput for Very large matrices? 1-Yes;0-No : ');
    
    if block_mark == 1
        % set block size according to computer memory
        BlockSize = 5000;
        BlockNum = ceil(ArrayLen/BlockSize);
        % pre allocate memory
        b_merge = zeros(1,ArrayLen);
        p_merge = zeros(1,ArrayLen);
        t_merge = zeros(1,ArrayLen);
%         para_merge = cell(1,ArrayLen);
        for jj = 1:BlockNum
            disp([num2str(BlockNum-jj),' Blocks Left']);
            % get the array index
            if jj < BlockNum
                ArrayIdx = BlockSize*(jj-1)+1:BlockSize*jj;
            else
                ArrayIdx = BlockSize*(jj-1)+1:ArrayLen;
            end
            Temp1 = Data_Grp1(:,ArrayIdx);
            Temp2 = Data_Grp2(:,ArrayIdx);

            % comput general linear model regression paras
%             [b,p,t] = cellfun(@vm_glmfit,num2cell(Temp1,1),num2cell(Temp2,1));%SubNum*BlockSize
            [b,p,t] = arrayfun(@(i) vm_glmfit(Temp1(:,i),Temp2(:,i)),1:length(ArrayIdx));
            b_merge(ArrayIdx) = b;
            p_merge(ArrayIdx) = p;
            t_merge(ArrayIdx) = t;
%             para_temp = arrayfun(@(i) wvm_glmfit(Temp1(:,i),Temp2(:,i)),1:length(ArrayIdx),'UniformOutput',0);
%             para_merge(ArrayIdx) = para_temp;

        end
%         MergeParasMat = cell2mat(para_merge');
%         b_merge = MergeParasMat(:,1);
%         p_merge = MergeParasMat(:,2);
%         t_merge = MergeParasMat(:,3);
    else % do not use block comput
%         MergeParas = arrayfun(@(i) wvm_glmfit(Data_Grp1(:,i),Data_Grp2(:,i)),[1:ArrayLen],'UniformOutput',0);
%         MergeParasMat = cell2mat(MergeParas');
%         b_merge = MergeParasMat(:,1);
%         p_merge = MergeParasMat(:,2);
%         t_merge = MergeParasMat(:,3);
        [b_merge,p_merge,t_merge] = arrayfun(@(k) vm_glmfit(Data_Grp1(:,k),Data_Grp2(:,k)),1:ArrayLen);
    end
    V = spm_vol(scan_grp1{1});
    Y = spm_read_vols(V);
    % reshape the merged result mat for img-written
    if isempty(explicit_mask)
        b_last = reshape(b_merge,size(Y));b_last(isnan(b_last)) = 0;
        p_last = reshape(p_merge,size(Y));p_last(isnan(p_last)) = 0;
        t_last = reshape(t_merge,size(Y));t_last(isnan(t_last)) = 0;
    else
        b_last = zeros(size(Y));b_last(MaskIdx) = b_merge;
        p_last = zeros(size(Y));p_last(MaskIdx) = p_merge;
        t_last = zeros(size(Y));t_last(MaskIdx) = t_merge;
    end
    % write res imgs

    V.fname = 'res_BetaRes.nii';spm_write_vol(V,b_last);
    V.fname = 'res_P_value.nii';spm_write_vol(V,p_last);
    V.fname = 'res_T_value.nii';spm_write_vol(V,t_last);
    movefile('res*.nii',output_dir);
    clear all;
    
end

% function m = wvm_glmfit(x,y)
%     [a1,a2,a3] = glmfit(x,y);
%     b = a1(2);
%     p = a3.p(2);
%     t = a3.t(2);
%     m = [b,p,t];
% end
function [b,p,t] = vm_glmfit(x,y)
    [O1,O2,O3] = glmfit(x,y);
    b = O1(2);
    p = O3.p(2);
    t = O3.t(2);
end


% function Matmerge = ReadMultiImg(ImgList)
%     ImgNum = length(ImgList);
%     for ii = 1:ImgNum
%         Vtemp = spm_vol(ImgList{ii});
%         Vmerge(ii,1) = Vtemp;
%     end
%     Matmerge = spm_read_vols(Vmerge);
% end
function [Y_merge,MaskIdx] = ReadMultiImg(ImgList,MaskImg)
%% function to read multiple images into datamat

    V = spm_vol(ImgList{1});
    Y = spm_read_vols(V);
    
    ImgNum = length(ImgList);
    if isempty(MaskImg)
        Y_merge = zeros(ImgNum,numel(Y));
        for ii = 1:ImgNum
            V_temp = spm_vol(ImgList{ii});
            Y_temp = spm_read_vols(V_temp);
            Y_merge(ii,:) = reshape(Y_temp,1,numel(Y));
        end
        MaskIdx = [];
    else
        
        Vm = spm_vol(MaskImg);
        Ym = spm_read_vols(Vm);
        
        % check dim matched or not
        if prod(double(size(Y)==size(Ym))) == 0
            error('Mask dimesion different from the Imgs');
        end
        
        MaskIdx = find(Ym~=0)'; % MaskIdx is 1*N format
        MaskLen = length(MaskIdx);
        Y_merge = zeros(ImgNum,MaskLen);
        
        for ii = 1:ImgNum
            V_temp = spm_vol(ImgList{ii});
            Y_temp = spm_read_vols(V_temp);
            Y_merge(ii,:) = Y_temp(MaskIdx);
        end
        
    end

end

function img_list_cell = ScanImg2Cell(varargin)
% scan img files and reture the cell format list
% ScanImg2Cell(dlg_title)
% written by hongshengcheng.math@gmail.com

    if numel(varargin) == 1
        dlg_title = varargin{1};
    else
        dlg_title = 'Select Img Files';
    end
    img_list = spm_select(Inf,'image',dlg_title,'',pwd,'.*',1);
    img_num = size(img_list,1);
    img_list_cell = cell(img_num,1);
    for ii = 1:img_num
        % modified 20141012
        % fix a bug which will generate extra space
        str_temp = img_list(ii,:);
        str_ok = str_temp(~isspace(str_temp));
        img_list_cell{ii} = str_ok;
    end
end

function [labels_iv,data_iv,labels_cv,data_cv] = ReadMultiBehData()
%% function to read multiple behavior and brain img data
% IV:independent variables
% CV:covariable variables

    [labels_iv,data_iv] = ReadXlsData('Input the Independent Variables');
    
    cv_mark = questdlg('Add covariates?', 'Questdlg','Yes','No','No');

    if strcmp(cv_mark,'Yes') == 1
        [labels_cv,data_cv] = ReadXlsData('Input the Covariates');
    else
        labels_cv = [];
        data_cv = [];
    end
end

function [ImgListNames,ImgListScans] = ReadMultiImgList()
%% function to read multiple brain img list

    img_list_name = input('Input group img name : ','s');
    ImgListScans{1,1} = ScanImg2Cell;
    ImgListNames{1,1} = img_list_name;
    
    % read multiple scans of img
    Mark2Run = 1;
    Mark_Runtime = 2;
    while (Mark2Run == 1)
        choice_add = questdlg('Add more files?', 'Questdlg','Yes','No','No');
        switch choice_add
            case 'Yes' 
                img_list_name = input('Input group img name : ','s');
                ImgListScans{Mark_Runtime,1} = ScanImg2Cell;
                ImgListNames{Mark_Runtime,1} = img_list_name;
                Mark_Runtime = Mark_Runtime+1;
            case 'No'
            Mark2Run = 0;
        end
    end
end

function [Labels,DataMat] = ReadXlsData(varargin)
%% function to copy data direct from xls file
% make sure the end of the data should be cleaned

    if isempty(varargin)
        dlg_title = 'Input Data from Excel';
    else
        dlg_title = varargin{1};
    end
    
    prompt = {'Enter Labels :','Enter Matrix :'};
    Paras = inputdlg(prompt,dlg_title,[1 100;10 100]);

    LabelRaw = Paras{1,1};
    space_mark = strfind(LabelRaw,' ');
    if isempty(space_mark)
        Labels = regexp(LabelRaw,'\t','split');
    else
        Labels = regexp(LabelRaw,' ','split');
    end
    Labels = Labels'; % n*1 cell mat
    
    DataRaw = Paras{2,1};
    DataMat = str2num(DataRaw);
%     if read raw data into one line, script follow will be fine
%     DataStr = regexp(DataRaw,' ','split');
%     RowNum = size(DataStr,2);
%     ColNum = size(Labels,1);
%     DataMat = zeros(RowNum,ColNum);
%     for ii = 1:RowNum
%         DataTempStr = regexp(DataStr{ii},'\t','split');
%         for jj = 1:ColNum
%             DataMat(ii,jj) = cellfun(@str2num,DataTempStr(jj));
%         end
%     end
end

function [grp_names,para_names,para_cell] = ReadMultiPathInfo(varargin)
%% function to read multiple path info
% Input
% Select Path.xls file
% Output
% grp_names,para_names,para_cell
% a para_path.xls should be in format:
% 1st column is the group name,2nd,3rd...is the parametre path
% 1st row is the name of parametres
% make sure all the paras have the same order

    if isempty(varargin)
        [filename,filepath] = uigetfile('*.xls;*.xlsx','Select Path file');
        para_path_loc = fullfile(filepath,filename);
    else
        para_path_loc = varargin{1};
    end

    [num,txt,raw] = xlsread(para_path_loc);

    para_names = raw(1,2:end);%the 1st column is the grp mark;1*n format cell
    para_num = length(para_names);

    para_paths = raw(2:end,2:end);
    grp_info = raw(2:end,1);

    grp_table = tabulate(grp_info);%get the grp index infomation

    grp_num = size(grp_table,1);
    grp_names = grp_table(:,1);
    grp_idx = cell(grp_num,1);
    para_cell = cell(grp_num,para_num);

    for ii = 1:grp_num
        idx_temp= find(strcmp(grp_info,grp_names(ii))==1);
        grp_idx{ii} = idx_temp;
        for jj = 1:para_num
            para_cell{ii,jj} = para_paths(idx_temp,jj);
        end
    end
end

function [cov,consess] = gen_cov_convec(label_iv,data_iv,label_cv,data_cv)

    if isempty(label_cv)
        labels = label_iv;
        data = data_iv;
    else
        labels = [label_iv;label_cv];
        data = [data_iv data_cv];
    end
    var_num = length(labels);
    
    % gen cov struct
    for jj = 1:var_num;
        cov(1,jj).c = data(:,jj);
        cov(1,jj).cname = labels{jj};
        cov(1,jj).iCFI = 1;
        cov(1,jj).iCC = 1;       
    end
    consess{1,1}.tcon.name = label_iv{1};
%     consess{1,1}.tcon.convec = [1 zeros(1,var_num-1)];
    % the first column is mean
    consess{1,1}.tcon.convec = [0 1 zeros(1,var_num-1)];
end
% function [cov,consess] = gen_cov_convec(Labels,DataMat)
% %     [Labels,DataMat] = ReadXlsData();
%     var_num = length(Labels);
%     
%     % gen cov struct
%     for jj = 1:var_num;
%         cov(1,jj).c = DataMat(:,jj);
%         cov(1,jj).cname = Labels{jj};
%         cov(1,jj).iCFI = 1;
%         cov(1,jj).iCC = 1;
%         
%         % initial a null convec
%         % the first column is mean
%         convec_null = zeros(1,var_num+1);
%         consess{1,jj}.tcon.name = Labels{jj};
%         convec_null(jj+1) = 1;
%         consess{1,jj}.tcon.convec = convec_null;
%         
%     end
% end