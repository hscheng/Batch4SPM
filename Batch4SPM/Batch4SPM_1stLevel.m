function Batch4SPM_1stLevel()
%% function to run 1st level
% Batch4SPM_1stLevel(DataDir,OutputDir)
% DataDir,preproced smoothed data dir
% OutputDir,store each subjet's 1st level result files
% 
% % data should be organized as follows:
% if each subject have MULTIPLE runs,then:
%     wkdir\sub*\run*\*.dcm, e.g. wkdir\sub01\run1..  wkdir\sub01\run2.
% if each subject have SINGLE run,then:
%     wkdir\sub*\*., e.g. wkdir\sub01\..;wkdir\sub02\..;wkdir\sub03\...
% 
% 20150611,modify a bug when read the preprocessed data file
% 20141215,modify error to load onset
% 20141205,fit for single run,which could use dparsf preproc res files
% 20150518,fit both for block design and rapid event related design
% 
% written by hongshengcheng.math@gmail.com
% created date:20141201

clc;
fprintf('Welcome to our lab:  <a href="http://www.qiujlab.com/">Qiu Jiang Lab</a>\nhttp://www.qiujlab.com/\n');
fprintf('<a href="http://www.creativitybrain.com/">Center for Creativity and Brain</a>\nhttp://www.creativitybrain.com/\n');
disp('Any question, Please contact me at hongshengcheng.math@gmail.com');
disp('****************************************************************');

%     spm('defaults','fmri');
%     spm_jobman('initcfg');
    

    DataDir = uigetdir(pwd,'Select Data Dir');
    [WorkDir,~,~] = fileparts(DataDir);
    [SubDataPathAll,SubID] = MergeDir(DataDir);
    SubNum = size(SubDataPathAll,1);
    
    % set the output dir as you need
%     OutputDir = input('Name of OutputDir (e.g.Stat_1stLevel) : ','s');
%     if isempty(OutputDir)
%         OutputDir = 'Stat_1stLevel';
%     end
%     FirstLevelResDir = fullfile(ParentDir,OutputDir);% 20140920'Stat_1stLevel'
%     mkdir(FirstLevelResDir);

    FirstLevelResDir = uigetdir(pwd,'select 1st level output dir');
    cd(FirstLevelResDir)
    
    % firstleve setup
    firstlevel_setup = Setup_FirstLevel(SubID);
    time_mark = datestr(clock,30);
%     save(fullfile(FirstLevelResDir,['firstlevel_setup_',time_mark,'.mat']),'firstlevel_setup');
    save(fullfile(FirstLevelResDir,['firstlevel_setup_',time_mark,'.mat']),'firstlevel_setup');
    
    OnesetMode = firstlevel_setup.OnesetMode;
    
    run_num = firstlevel_setup.run_num;
%     run_mark = questdlg('how many runs in your data','quest','multiple','single','multiple');
%         
%     if (strcmp(run_mark,'multiple') == 1)&&(OnesetMode == 1)

    if (run_num>1)&&(OnesetMode == 1)
    %% multiple runs and all subjects share the same onset
        firstlevel_setup.cond =  firstlevel_setup.onsets;%20141130
        if matlabpool('size')>1
            parfor ii = 1:SubNum
                stat_1st_multiruns(WorkDir,FirstLevelResDir,SubID{ii},SubDataPathAll{ii},firstlevel_setup);
            end
        else
            for ii = 1:SubNum
                stat_1st_multiruns(WorkDir,FirstLevelResDir,SubID{ii},SubDataPathAll{ii},firstlevel_setup);
            end
        end
    elseif (run_num>1)&&(OnesetMode ~= 1)
    %% multiple runs and NOT all subjects share the same onset
        % SubID = fieldnames(cond_all_subs);
        % In case some subject have contrast while not have brain imgs
        cond_all_subs = firstlevel_setup.onsets;
        for ii = 1:SubNum
            cond_sub = getfield(cond_all_subs,SubID{ii});
            firstlevel_setup.cond = cond_sub;
            stat_1st_multiruns(WorkDir,FirstLevelResDir,SubID{ii},SubDataPathAll{ii},firstlevel_setup);
        end
     elseif (run_num==1)&&(OnesetMode == 1)
     %% each subject have single run and all subjects share the same onset
        firstlevel_setup.cond =  firstlevel_setup.onsets;%20141130
        if matlabpool('size')>1
            parfor ii = 1:SubNum
                stat_1st_multiruns(WorkDir,FirstLevelResDir,SubID{ii},SubDataPathAll{ii},firstlevel_setup);
            end
        else
            for ii = 1:SubNum
                stat_1st_singlerun(WorkDir,FirstLevelResDir,SubID{ii},SubDataPathAll{ii},firstlevel_setup);
            end
        end
    elseif (run_num==1)&&(OnesetMode ~= 1)
    %% each subject have single run and NOT all subjects share the same onset
            % SubID = fieldnames(cond_all_subs);
            % In case some subject have contrast while not have brain imgs
            cond_all_subs = firstlevel_setup.onsets;
            for ii = 1:SubNum
                cond_sub = getfield(cond_all_subs,SubID{ii});
                firstlevel_setup.cond = cond_sub;
                stat_1st_singlerun(WorkDir,FirstLevelResDir,SubID{ii},SubDataPathAll{ii},firstlevel_setup);
            end
    end
    cd(FirstLevelResDir)

end

function firstlevel_setup = Setup_FirstLevel(SubID)
%% function to set up first level parameters

% ******************* setup model specification parameters ****************
%     % add 20150518,ahong
%     design_type = questdlg('Design Type','Quest','RapidEvent','Block','RapidEvent');
%     firstlevel_setup.design_type = design_type;
    design_unit = questdlg('Units for design','SelectionMode','scans','seconds','scans');
    para_dlg=inputdlg({'Units for design',...
        'Interscan interval(TR)',...
        'Microtime resolution(Slice Num)',...
        'Microtime onset(Middle Slice Order)',...
        'Duration(how many TR,Length = Cond_Num or 1)',...
        'Condition Num'},...
        'Please input the parameter', ...
        [1 50;1 50;1 50;1 50;1 50;1 50],{design_unit,'2','32','16','0','4'}) ;

    firstlevel_setup.units = design_unit;
    firstlevel_setup.RT = str2num(para_dlg{2});
    firstlevel_setup.fmri_t = str2num(para_dlg{3});
    firstlevel_setup.fmri_t0 = str2num(para_dlg{4});

    % firstlevel_setup.duration = str2num(para_dlg{5});
    % 20141014,fix bug for diff condition may have diff duration
    duration_temp = str2num(para_dlg{5});
    cond_num = str2num(para_dlg{6});
    % check duration set
    if length(duration_temp) == cond_num
        duration = duration_temp;
        display('duration is OK');
    elseif length(duration_temp) == 1% all conditon have same duration
        duration = repmat(duration_temp,1,cond_num);
    else
        error('duration wrong,length should equal to cond_num or 1');
    end
    firstlevel_setup.duration = duration;
        
%**************************** read onset.xls ******************************
% generate cond struct which contain the onsets and duration set up

    OnesetMark = questdlg('All subjects share the same onsets?','Ask','Yes','No','Yes');
    % if you choose No, it means that not all subject have same onsets
    % and the run.xls file should be in format
    % SubID Condition1 SubID Conditon_2 ...
    % First Raw should all be strings
    % default Yes--all subjects shares the same onsets
    
    % read the run xls file
    [run_filename,run_filepath] = uigetfile('*.xls;*.xlsx','Select Onset.xls file');
    run_fullpath = fullfile(run_filepath,run_filename);
    
    [~,sheets] = xlsfinfo(run_fullpath);
    %get run_num
    [~,run_num]=size(sheets);
    firstlevel_setup.run_num = run_num;%20141215
    
    if strcmp(OnesetMark,'Yes') == 1
        firstlevel_setup.OnesetMode = 1;
        if ispc == 1 %os = 'windows'; 
            % read the run xls file
            cond_sub = gen_cond(run_fullpath,duration);
        else % for os is linux
            cond_sub = gen_cond_linux(duration);
        end 
        
        firstlevel_setup.onsets = cond_sub;

    else
        firstlevel_setup.OnesetMode = 0;
        cond_all_subs.mark = 1; 
        % if not, there will be an error report cond_all_subs is not exist
        SubNum = length(SubID);
        for ii = 1:SubNum
            msg_head = sprintf('Reading %s Onset Info...',SubID{ii});
            disp(msg_head);
            WriteRunLog('log_OnsetCheck.txt',msg_head);
            
            cond_one_sub = gen_cond_advanced(run_fullpath,SubID{ii},duration);          
            cond_all_subs = setfield(cond_all_subs,SubID{ii},cond_one_sub); 
        end
        firstlevel_setup.onsets = cond_all_subs;
    end
        
    % run merge mode
    merge_mark = questdlg('Concatenate runs ?','Choose','NO','Yes','NO');
    switch merge_mark
        case 'NO'
            run_merge = 0;
        case 'Yes'
            run_merge = 1;
            % block regressors for multiple runs or sessions
            reg_choose = questdlg('add session regressors?','Quest','Yes','No','No');
            firstlevel_setup.reg_choose = reg_choose;
    end
    firstlevel_setup.run_merge = run_merge;

%************************** read contrast.txt ****************************
% generarate consess struct which contain the contrast information
contrast_mode = questdlg('contrast kind','Quest','t_con_only','f_con_only','both_con','t_con_only');
switch contrast_mode
    case 't_con_only'
        [con_filename,con_filepath] = uigetfile('*.txt','Select T Contrast.txt file');
        con_txt_loc = fullfile(con_filepath,con_filename);
        consess_1run = gen_consess_tcon(con_txt_loc);
    case 'f_con_only'
        consess_1run = gen_consess_fcon(); %cond_num = size(cond_sub,2);eye(cond_num)
    case 'both_con'
        [con_filename,con_filepath] = uigetfile('*.txt','Select T Contrast.txt file');
        con_txt_loc = fullfile(con_filepath,con_filename);
        consess_t = gen_consess_tcon(con_txt_loc);
        consess_f = gen_consess_fcon();
        consess_1run = [consess_t consess_f];
end
firstlevel_setup.consess = consess_1run;

end

function stat_1st_multiruns(WorkDir,OutputDir,SubID,SubDataPath,firstlevel_setup)
%% function for one subject 1st level

    % mkdir to store the resluts
    SubOutputDir = fullfile(OutputDir,SubID);
    mkdir(SubOutputDir);
    cd(SubOutputDir);
    
    SubRpDir = fullfile(WorkDir,'RealignParameter',SubID);
    
    cond_sub = firstlevel_setup.cond;
    run_merge = firstlevel_setup.run_merge;%1 merge;0 not;
    cond_num = size(cond_sub,2);

    runlist_smooth = dir(fullfile(SubDataPath,'run*'));
    run_num = length(runlist_smooth);
    run_len = zeros(1,run_num);
    scans_merge = [];
    
%% model specification

    % output dir
    matlabbatch{1}.spm.stats.fmri_spec.dir={SubOutputDir};%cell format

    %timing
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = firstlevel_setup.units; % 'scans'
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = firstlevel_setup.RT; % 2
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = firstlevel_setup.fmri_t; % 32
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = firstlevel_setup.fmri_t0; % 16
%     % when concatenate runs, there should have a scaling process
%     % which will set all the runs comparetable
%     if run_merge == 1
%         matlabbatch{1,1}.spm.stats.fmri_spec.global = 'Scaling';
%     end
    

    % specify images that need transformation
    for ii = 1:run_num % run num
        run_id = runlist_smooth(ii,1).name;
        % 20141215 modify for error load onset
        run_mark = str2num(run_id(end));%run1,run2,...run9...
        % modified at 20150518,by hscheng
        % scanslist = dir(fullfile(SubDataPath,run_id,'swra*.nii'));
        scanslist = dir(fullfile(SubDataPath,run_id,'*.nii'));
        % add,20150813,search for img file ,if nii file is not found
        if isempty(scanlist)
            scanlist = dir(fullfile(SubDataPath,run_id,'*.img'));
        end
        % read each runs' img files
        trials_run = length(scanslist);
        run_len(ii) = trials_run;
        scans_run = cell(trials_run,1);
        for k = 1:trials_run
            scans_run{k,1} = fullfile(SubDataPath,run_id,scanslist(k).name);
        end
        
        % locate the rp txt file in each run
        rp_filename = fullfile(SubRpDir,['rp_a',SubID,'_',run_id,'*.txt']);
        rp_run_loc = dir(rp_filename);
        rp_run_path = fullfile(SubRpDir,rp_run_loc.name);
        
        if run_merge == 1 % merge the runs img files in one scan and rp files
            scans_merge = [scans_merge;scans_run];
            
            % merge rp files
            rp_temp = importdata(rp_run_path);

            % modify the condition set
            if ii == 1
                cond_merge = cond_sub(run_mark,:);
                trials_num = trials_run;
                rp_merge = rp_temp;   
            else
                for jj =1:cond_num
                    cond_merge(1,jj).onset = [cond_merge(1,jj).onset; cond_sub(run_mark,jj).onset + trials_num];
                end
            trials_num = trials_num + trials_run;
            rp_merge = [rp_merge;rp_temp];
            end
        else
            matlabbatch{1}.spm.stats.fmri_spec.sess(1,ii).scans = scans_run;
            matlabbatch{1}.spm.stats.fmri_spec.sess(1,ii).cond = cond_sub(run_mark,:);
            matlabbatch{1}.spm.stats.fmri_spec.sess(1,ii).multi_reg = {rp_run_path};
        end
    end

    if run_merge == 1 % block regression only for merger runs condition
        rp_merge_txt = fullfile(SubOutputDir,['rp_merge_',SubID,'.txt']);
        save(rp_merge_txt,'rp_merge','-ascii')

        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans_merge;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = cond_merge;
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {rp_merge_txt};
        
        % add for PPI,20141129
        % reg_choose = questdlg('add session regressors?','Quest','Yes','No','No');
        reg_choose = firstlevel_setup.reg_choose;
        if strcmp(reg_choose,'Yes') == 1
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress = gen_sessionreg(run_num,run_len);
        end
    end
    
%% model estimation

    % get the SPM.mat Subject Dir
    spm_mat = fullfile(SubOutputDir,'SPM.mat');
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {spm_mat};
    matlabbatch{2}.spm.stats.fmri_est.method.Classical=1;

%% contrast
    consess = firstlevel_setup.consess;
    
    if run_merge == 0 % when not merge multiple runs
        con_num = size(consess,2);
        for jj =1:con_num
            con_mark = fieldnames(consess{1,jj});
            switch con_mark
                case 'tcon'
                    contrast_1run = consess{1,jj}.tcon.convec;
                    contrast_allrun = repmat(contrast_1run,1,run_num);
                    consess{1,jj}.tcon.convec = contrast_allrun;
                case 'fcon'
                    contrast_1run = consess{1,jj}.fcon.convec{1};
                    contrast_allrun = repmat(contrast_1run,1,run_num);
                    consess{1,jj}.fcon.convec{1} = contrast_allrun;
            end
            clear contrast_1run;
        end
    end
    
    matlabbatch{3}.spm.stats.con.spmmat = {spm_mat};
    matlabbatch{3}.spm.stats.con.consess = consess;
    matlabbatch{3}.spm.stats.con.delete = 0;
    
%% save and run matlabbatch
    batch_mat_list = dir(fullfile(SubOutputDir,'batch_1stlevel_*.mat'));
    if ~isempty(batch_mat_list)
        mkdir(SubOutputDir,'old_backup')
        time_mark = datestr(clock,30);
        batch_mat_file = batch_mat_list.name;
        [~,mat_name] = fileparts(batch_mat_file);
        
        source_dir = fullfile(SubOutputDir,batch_mat_file);
        target_dir = fullfile(SubOutputDir,'old_backup',[mat_name,'_backup_',time_mark,'.mat']);
        movefile(source_dir,target_dir);
    end
    save(fullfile(SubOutputDir,['batch_1stlevel_',SubID,'.mat']),'matlabbatch');
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    
    clear all;
    disp('One Subject works finished')
end

function stat_1st_singlerun(WorkDir,OutputDir,SubID,SubDataPath,firstlevel_setup)
%% function for one subject 1st level

    % mkdir to store the resluts
    SubOutputDir = fullfile(OutputDir,SubID);
    mkdir(SubOutputDir);
    cd(SubOutputDir);
        
    SubRpDir = fullfile(WorkDir,'RealignParameter',SubID);
    
    cond_sub = firstlevel_setup.cond;
    
%% model specification

    % output dir
    matlabbatch{1}.spm.stats.fmri_spec.dir={SubOutputDir};%cell format
    %timing
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = firstlevel_setup.units; % 'scans'
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = firstlevel_setup.RT; % 2
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = firstlevel_setup.fmri_t; % 32
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = firstlevel_setup.fmri_t0; % 16
%     % when concatenate runs, there should have a scaling process
%     % which will set all the runs comparetable
%     if run_merge == 1
%         matlabbatch{1,1}.spm.stats.fmri_spec.global = 'Scaling';
%     end
    

    % specify images that need transformation
    % modified at 20150518,by hscheng
    % scanslist = dir(fullfile(SubDataPath,'swra*.nii'));
    scanslist = dir(fullfile(SubDataPath,'*.nii'));
    % add,20150813,search for img file ,if nii file is not found
    if isempty(scanlist)
        scanlist = dir(fullfile(SubDataPath,'*.img'));
    end
    % read each runs' img files
    trials_run = length(scanslist);
    scans_run = cell(trials_run,1);
    for k = 1:trials_run
        scans_run{k,1} = fullfile(SubDataPath,scanslist(k).name);
    end

    % locate the rp txt file in each run
    rp_filename = fullfile(SubRpDir,'rp_*.txt');
    rp_run_loc = dir(rp_filename);
    rp_run_path = fullfile(SubRpDir,rp_run_loc.name);
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans_run;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = cond_sub;
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {rp_run_path};

%% model estimation

    % get the SPM.mat Subject Dir
    spm_mat = fullfile(SubOutputDir,'SPM.mat');
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {spm_mat};
    matlabbatch{2}.spm.stats.fmri_est.method.Classical=1;

%% contrast
    consess = firstlevel_setup.consess;
    matlabbatch{3}.spm.stats.con.spmmat = {spm_mat};
    matlabbatch{3}.spm.stats.con.consess = consess;
    matlabbatch{3}.spm.stats.con.delete = 0;
    
%% save and run matlabbatch
    batch_mat_list = dir(fullfile(SubOutputDir,'batch_1stlevel_*.mat'));
    if ~isempty(batch_mat_list)
        mkdir(SubOutputDir,'old_backup')
        time_mark = datestr(clock,30);
        batch_mat_file = batch_mat_list.name;
        [~,mat_name] = fileparts(batch_mat_file);
        
        source_dir = fullfile(SubOutputDir,batch_mat_file);
        target_dir = fullfile(SubOutputDir,'old_backup',[mat_name,'_backup_',time_mark,'.mat']);
        movefile(source_dir,target_dir);
    end
    save(fullfile(SubOutputDir,['batch_1stlevel_',SubID,'.mat']),'matlabbatch');
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    
    clear all;
    disp('One Subject works finished')
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

function cond = gen_cond(run_xls_loc,duration)
%% generate cond struct which contains the onset and duration
% only for all subject shares the same onsets and duration,
% make sure duration is condition_num*1 format
% 
% Format of the run.xls file shoul be
% 1 condition 1 column,first row is the name of the condition
% just like
% cond1 cond2 cond3 cond4
%   2     5      6     8
%   3     7      9     11
%   12    15    16     20
% so the clumn num is the conditon num
% 1 sheet for 1 run, do not left blank sheet in run.xls

    [~,sheets] = xlsfinfo(run_xls_loc);

    %get run_num
    [~,run_num]=size(sheets);

    % sess.cond which contain your conditions and onset setup

    for nn=1:run_num
        [onset,cond_name,~] = xlsread(run_xls_loc,sheets{1,nn});%load each run's onset

        [~,col_num] = size(cond_name);
        cond_num = col_num;

        %disp the information of each run
        msg = sprintf('In %s ,your data have %d conditions.', sheets{1,nn}, cond_num);
        disp(msg)

        for mm = 1:cond_num
        % this part is modified 20140718
            cond(nn,mm).name = cond_name{1,mm};

            %set oneset time point
            onset_list = onset(:,mm);%onset n*1
            onset_last = onset_list(isfinite(onset_list));
            if ~isempty(onset_last)
                cond(nn,mm).onset = onset_last;
            else
                cond(nn,mm).onset = 1;
                msg_check = sprintf('In run %d condition %d,blank onset exists!!',nn,mm);
                disp(msg_check);
                WriteRunLog('log_OnsetCheck.txt',msg_check);
            end

            %set duration
%             cond(nn,mm).duration = duration;
            % 20141014 modified, for different condition may have different
            % duration
            cond(nn,mm).duration = duration(mm);

            cond(nn,mm).tmod = 0;
            cond(nn,mm).pmod = struct([]);
        end
    end
%     save cond cond  
end

function cond = gen_cond_advanced(run_xls_loc,SubID,duration)
%% gen cond if not all subjects have the same onsets
% the run.xls format should like below :
% SubID     Cond_1	SubID       Cond_2 ...
% sub_24901 2       sub_24901	3
% sub_24901	4       sub_24901	5
% sub_24901	6       sub_24901	7
% sub_24906	18      sub_24906	19
% sub_24906	20      sub_24906	21
% sub_24906	22      sub_24906	23
% sub_24906	24      sub_24906	25
% ...
% 1 sheet 1 run,do not leave blank sheet in run.xls
% make sure duration is condition_num*1 format
[~,sheets] = xlsfinfo(run_xls_loc);

%get run_num
[~,run_num]=size(sheets);

% sess.cond which contain your conditions and onset setup
% [content, index] = vlookup(raw, 'sub01', 2, 1);
for nn = 1:run_num
    [onset,cond_name,raw_all] = xlsread(run_xls_loc,sheets{1,nn});%load each run's onset

    [~,col_num]=size(cond_name);
        cond_num = col_num/2;
%         %disp the information of each run
%         msg = sprintf('In %s ,your data have %d conditions.', sheets{1,nn}, cond_num);
%         disp(msg)

        for mm = 1:cond_num
            cond_mat = raw_all(2:end,(mm-1)*2+1:mm*2);
            [content, index] = vlookup(cond_mat,SubID, 2, 1);
            cond(nn,mm).name = cond_name{1,mm*2};

            %set oneset time point
            cond_onset = onset(:,(mm-1)*2+1);%onset n*1
            if ~isempty(index)
                cond(nn,mm).onset = cond_onset(index);
            else
                cond(nn,mm).onset = 1;
                msg_check = sprintf('In run %d condition %d,blank onset exists!!',nn,mm);
                disp(msg_check);
                WriteRunLog('log_OnsetCheck.txt',msg_check);
            end

            %set duration
%             cond(nn,mm).duration = duration;
            % 20141014,different condition may have different duration
            cond(nn,mm).duration = duration(mm);
            cond(nn,mm).tmod=0;
            cond(nn,mm).pmod=struct([]);
        end
 end

end

function cond = gen_cond_linux(duration)
%% generate cond struct which contains the onset and duration
% date:20140820
% Format of the onset.txt should be[under linux cannot read *.xls file]
% run1_cond1 run1_cond2 run2_cond1 run2_cond2 ...
%      1         3          4          6
%      5         6          8          9
%      10        12         15         20
% make sure each run have same condition number
% make sure duration is condition_num*1 format

    [name,path] = uigetfile('*.txt','Select Onset.txt file');
    onset_txt_loc = fullfile(path,name);
    M = importdata(onset_txt_loc);
    onset = M.data;
    cond_name = M.textdata;
    run_num = input('RUN Number : ');
    cond_num = size(onset,2)/run_num;

    % sess.cond which contain your conditions and onset setup

    for nn = 1:run_num
        for mm = 1:cond_num
        % this part is modified 20140718
            cond(nn,mm).name = cond_name{1,cond_num*(nn-1)+mm};

            %set oneset time point
            onset_list = onset(:,cond_num*(nn-1)+mm);%onset n*1
            onset_last = onset_list(isfinite(onset_list));
            if ~isempty(onset_last)
                cond(nn,mm).onset = onset_last;
            else
                cond(nn,mm).onset = 1;
                msg_check = sprintf('In run %d condition %d,blank onset exists!!',nn,mm);
                disp(msg_check);
                WriteRunLog('log_OnsetCheck.txt',msg_check);
            end
            %set duration
%             cond(nn,mm).duration = duration;
            % 20141014,different condition may have different duration
            cond(nn,mm).duration = duration(mm);
            cond(nn,mm).tmod = 0;
            cond(nn,mm).pmod = struct([]);
        end
    end
%     save cond cond  
end

function consess_t = gen_consess_tcon(con_txt_loc)
%% function to gen consess,which contains the contrast infomation
% Format of the contrast.txt file shoul be like
% WF-WZ 0 1 0 -1 0 0 0 0 0 0
% WF-JF -1 1 0 0 0 0 0 0 0 0
% WZ-JZ 0 0 -1 1 0 0 0 0 0 0
%
% 6 zeros in each contrast conditon[6 headmotion parameters]
% first string of each row is the name of contrast
% Note:Only support the T Contrast

    txt_info = importdata(con_txt_loc);

    tcon_num = size(txt_info.data,1);
    tcon_name = txt_info.textdata;
    % disp contrast infomation
    msg = sprintf('You have %d constracts.', tcon_num);
    disp(msg)
    convec_1run = txt_info.data;
%     this part move to the contrast part
%     in case some subject have less or more runs
% %  ********************* merge_choose ***************************
%     merge_choose=questdlg('Concatenate runs ?',...
%         'Make a choice',...
%         'NO','Yes','NO');
%     switch merge_choose
%         case 'NO'
%             run_merge = 0;
%             % make sure contrast txt file add 6 zeros in each contrast
%             convec_all = repmat(txt_info.data,1,run_num);
% 
%         case 'Yes'
%             run_merge = 1;
%             convec_all = txt_info.data;
%     end

    % gen consess struct
    %make sure that consess is 1*n cell 'struct'
    for kk = 1:tcon_num
        consess_t{1,kk}.tcon.name = tcon_name{kk,1};
        consess_t{1,kk}.tcon.convec = convec_1run(kk,:);
        consess_t{1,kk}.tcon.sessrep = 'none';
    end
end

function consess_f = gen_consess_fcon()
%% function to define multiple f contrast
    % define input dlg
    prompt = {'Enter the name of F-contrast:',...
        'Enter the contrast matrix:'};
    dlg_title = 'Input parameters';
    numlines = 1;
    para_default = {'fcon','eye(3)'};% e.g.kron([1 0 0 0]',ones(90,1))
    para = inputdlg(prompt,dlg_title,numlines,para_default);

    consess_f{1,1}.fcon.name = para{1};
    consess_f{1,1}.fcon.convec = {eval(para{2})};
    % add more regeressors
    Mark2Run = 1;
    Runtime = 2;
    while (Mark2Run == 1)
        choice_add = questdlg('Add more f contrst?', 'Questdlg','Yes','No','No');
        switch choice_add
            case 'Yes' 
                para = inputdlg(prompt,dlg_title,numlines,para_default);
                consess_f{1,Runtime}.fcon.name = para{1};
                consess_f{1,Runtime}.fcon.convec = {eval(para{2})};
                Runtime = Runtime+1;
            case 'No'
            Mark2Run = 0;
        end
    end
end

function consess = gen_consess()
%% function to define multiple f contrast and t contrast
    % define input dlg
    prompt = {'Enter the name of contrast:',...
        'Enter the contrast matrix:'};
    dlg_title = 'Input parameters';
    numlines = 1;
    para_default = {'tcon','1 1 0 0'};% e.g.kron([1 0 0 0]',ones(90,1))

    % add more regeressors
    Mark2Run = 1;
    Runtime = 1;
    while (Mark2Run == 1)
        choice_add = questdlg('Add more contrast?', 'Questdlg','T_Contrast','F_Contrast','No','No');
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

function regress = gen_sessionreg(session_num,session_len)
%% function to ass multiple regressors for 1st level analysis
% regress = firstlevel_autoaddreg(session_num,session_len)
% session_num e.g.4
% session_len e.g. [90 90 90 90]
% written by hongshengcheng.math@gmail.com
if length(session_len)~= session_num
    error('para error')
end

session_idx = cumsum(session_len);

for ii = 1: session_num

    regress(ii).name = ['Session',int2str(ii)];
    val_init = zeros(session_idx(end),1);
    if ii == 1;
        val_init(1:session_idx(ii)) = 1;
    else
        val_init(session_idx(ii-1)+1:session_idx(ii)) = 1;
    end
    regress(ii).val = val_init;
end
% regress(2).name = 'Session 2';
% regress(2).val = kron([0 1 0 0]',ones(90,1));
% regress(3).name = 'Session 3';
% regress(3).val = kron([0 0 1 0]',ones(90,1));


end

function [content, index] = vlookup(m, e, column, lookcolumn)
%% VLOOKUP the function as vlookup in Excel
% http://www.mathworks.com/matlabcentral/fileexchange/29233-vlookup-similar-to-ms-excel-function
% this function is needed when generate cond advaneced function
% 
%   [content, index] = vlookup(m, e, column, lookcolumn) look for 'e' in 
%   the 'lookcolumn'-th column of 'm', and return the coresponding
%   'column'-th element of 'm' in the same row.
%
%   the 'm' could be a numeric matrix of a cell matrix.
% 
%   lookcolumn is 1 by default if omitted.
% 
% Example:
% 
%     m = {1, 'a', [2 3];
%     2, 'b', 'cd'
%     3, 'a', true;};
%      [content, index] = vlookup(m, 'a', 3, 2) then
%     content = {[2 3], 1};
%     index = [1;3]

% Copyright: zhang@zhiqiang.org, 2010
% author: http://zhiqiang.org/blog/tag/matlab

if isempty(m) || isempty(e), return; end
if nargin <= 3, lookcolumn = 1; end

isechar = ischar(e);
assert(isechar || isnumeric(e), 'the second parameter must be a string or numeric');

if iscell(m)
    content = {}; index = [];
    if isechar
        index = find(strcmp(e, m(:, lookcolumn)));
        content = m(index, column);
    else
        for i = 1:size(m, 1)
            if isnumeric(m{i, lookcolumn}) && m{i, lookcolumn} == e
                content = [content; m(i, column)]; %#ok<*AGROW>
                index = [index; i];
            end
        end
    end
else
    assert(~isechar, 'When the first para is a matrix, the second para must be numeric');
    
    index = find(m(:, lookcolumn) == e);
    content = m(index, column);
end
end
