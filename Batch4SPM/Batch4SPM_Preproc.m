function Batch4SPM_Preproc()
%% function for preprocess of task fmri data
% start at the nii 3D img file,FunImg
% output 4 dirs respond to 4 steps
% FunImg_Step1_Slice;FunImg_Step2_Realign;% FunImg_Step3_Normalize;% FunImg_Step4_Smooth
% RealignParameter
% 
% data should be organized as follows:
%     wkdir\sub*\run*\*.dcm, e.g. wkdir\sub01\run1..  wkdir\sub01\run2.
% if each subject has single run,try DPASFA/dpabi
% 
% written by hongshengcheng.math@gmail.com
% 20150411,add dcm2nii
% created date:20141201

fprintf('Welcome to our lab:  <a href="http://www.qiujlab.com/">Qiu Jiang Lab</a>\nhttp://www.qiujlab.com/\n');
disp('Any question, Please contact me at hongshengcheng.math@gmail.com');
disp('****************************************************************');

%     spm('defaults','fmri');
%     spm_jobman('initcfg');
    
    start_file_type = questdlg('start file type','quest','dcm','nii','dcm');
    DataDir = uigetdir(pwd,'Select Data Dir');
    cd(DataDir)
    [ParentDir,~,~] = fileparts(DataDir);
    if strcmp(start_file_type,'dcm')==1
        Batch4SPM_dcm2nii(DataDir);
        FunImgDir = fullfile(ParentDir,'FunImg');
        clc 
        continue_mark = questdlg('Image format convert done','quest',...
            'continue for preproc','do nothing','continue for preproc');
        if strcmp(continue_mark,'continue for preproc')==1
            [SubDataPathList,SubID] = MergeDir(FunImgDir);
        end
    else
        [SubDataPathList,SubID] = MergeDir(DataDir);
    end
    
    
    SubNum = size(SubDataPathList,1);
    
    % set up preproc parameters
    preproc_setup = Setup_Preproc();   
    time_mark = datestr(clock,30); % add time mark to mat file
    save(['preproc_setup_',time_mark,'.mat'],'preproc_setup');
    
    if matlabpool('size')>1
        parfor ii = 1:SubNum
            Sub_Preproc(ParentDir,SubDataPathList{ii},SubID{ii},preproc_setup);    
        end
    else
        spm fmri;
        for ii = 1:SubNum
            Sub_Preproc(ParentDir,SubDataPathList{ii},SubID{ii},preproc_setup);
        end
    end
    
    cd(ParentDir)
    run_modules = preproc_setup.run_modules;
    % if start after reaglin then check normalization
    if min(run_modules) < 4
        % gen Normalize check imgs use the mean normalized imgs
        Dir_NormalizeCheck = fullfile(ParentDir,'RealignParameter','NormalizeCheck');
        NormalizeCheckImg(Dir_NormalizeCheck);
    end
    
end

function preproc_setup = Setup_Preproc()
%% function to setup preproc parameters
    para_dlg = inputdlg({'TR',...
        'Slice Order,e.g.[2:2:32 1:2:31]',...
        'Referance Slice [Time mid slice]',...
        'Smooth FWHM.e.g.[8 8 8]',...
        'Register to Mean [1-Yes;0-No]'},...
        'Please input the parameter of Preproc', ...
        [1 50;1 50;1 50;1 50;1 50],{'2',num2str([2:2:32 1:2:31]),'32',num2str([8 8 8]),'1'}) ;

    preproc_setup.tr = str2num(para_dlg{1});

    slice_order = str2num(para_dlg{2});
    preproc_setup.so = slice_order;

    slcie_num = size(slice_order,2);
    preproc_setup.nslices = slcie_num;

    preproc_setup.refslice = str2num(para_dlg{3}); 

    preproc_setup.fwhm = str2num(para_dlg{4});

    preproc_setup.bb = [-90 -126 -72;90 90 108];
    
    % PET Imgs set rtm 1 i.e. register to mean
    % MRI typically be 0
    % MRI images are typically registered to the first image. 
    % The more accurate way would be to use a two pass  procedure
    % but this probably wouldn't improve the results so much
    % and would take twice as long
    preproc_setup.rgm = str2num(para_dlg{5});
    % select run modules
    [ run_modules,~] =listdlg('ListString',{'Slice Timing','Realigh_Estimate&Reslice','Normalise_Estimate&Write','Smooth'},...
    'Name','Run Modules','SelectionMode','Multiple',...
    'ListSize',[200 100],'InitialValue',1:4);
    preproc_setup.run_modules = run_modules;

end

%% ********************* One Subject Preprocess ********************
function Sub_Preproc(WorkDir,SubDataPath,SubID,preproc_setup)
% batch for preprocessing fMRI data in SPM8

% ****************** Make Dirs to store Results Imgs **********************

Dir_Slice = fullfile(WorkDir,'FunImg_Step1_Slice',SubID);
mkdir(Dir_Slice);

Dir_Realign = fullfile(WorkDir,'FunImg_Step2_Realign',SubID);
mkdir(Dir_Realign);

Dir_Normalize = fullfile(WorkDir,'FunImg_Step3_Normalize',SubID);
mkdir(Dir_Normalize);

Dir_Smooth = fullfile(WorkDir,'FunImg_Step4_Smooth',SubID);
mkdir(Dir_Smooth);

Dir_RealignParameter = fullfile(WorkDir,'RealignParameter',SubID);
mkdir(Dir_RealignParameter);

Dir_rp_Report = fullfile(WorkDir,'RealignParameter','rp_Report');
mkdir(Dir_rp_Report);

Dir_HeadMotionPlot = fullfile(WorkDir,'RealignParameter','HeadMotionPlot');
mkdir(Dir_HeadMotionPlot);

Dir_NormalizeCheck = fullfile(WorkDir,'RealignParameter','NormalizeCheck');
mkdir(Dir_NormalizeCheck);

mkdir(Dir_rp_Report,'3mm_cutoff');
mkdir(Dir_rp_Report,'2mm_cutoff');
mkdir(Dir_rp_Report,'1mm_cutoff');
run_modules = preproc_setup.run_modules;
% run module
if ~isempty(find(run_modules == 1, 1))
    SubPreproc_SliceTiming(WorkDir,SubDataPath,SubID,preproc_setup)
end
if ~isempty(find(run_modules == 2, 1))
    Sub_Realign(WorkDir,SubDataPath,SubID,preproc_setup)
end
if ~isempty(find(run_modules == 3, 1))
    Sub_Normalise(WorkDir,SubDataPath,SubID,preproc_setup)
end
if ~isempty(find(run_modules == 4, 1))
    Sub_Smooth(WorkDir,SubDataPath,SubID,preproc_setup)
end
end
function SubPreproc_SliceTiming(WorkDir,SubDataPath,SubID,preproc_setup)
clear matlabbatch
cd(SubDataPath)

% list all the runs
% make sure you don't have any files named 'run*' under the same directory 
% dir setting
Dir_Slice = fullfile(WorkDir,'FunImg_Step1_Slice',SubID);

runlist_rawnii = dir(fullfile(SubDataPath,'run*'));
run_num = length(runlist_rawnii);
sub_scans = cell(run_num,1);
for ii = 1:run_num
    
    % run_id = ['run',int2str(ii)];
    % if subdata only have run2 run4 run5..,there will be an error
    % modified date:20140720
    
    run_id = runlist_rawnii(ii,1).name;
    scanslist = dir(fullfile(SubDataPath,run_id,'*.nii'));
    run_tr_num = length(scanslist);
    run_scans = cell(run_tr_num,1);
    for kk = 1:run_tr_num
        run_scans{kk,1} = [fullfile(SubDataPath,run_id,scanslist(kk).name),',1'];
    end
    sub_scans{ii,1} = run_scans;

end

% fill in fields of structure matlabbatch
matlabbatch{1}.spm.temporal.st.scans = sub_scans;
% number of slices
matlabbatch{1}.spm.temporal.st.nslices = preproc_setup.nslices; % default 32
% TR in seconds
matlabbatch{1}.spm.temporal.st.tr = preproc_setup.tr; % default 2
% TE = TR - (TR/nslices)
tr = preproc_setup.tr;
nslices = matlabbatch{1}.spm.temporal.st.nslices;
matlabbatch{1}.spm.temporal.st.ta = tr - (tr/nslices);

% slice order very important. depend on EPI sequence used
% here it's interleaved
% ask JC what order you should use
% slice order
matlabbatch{1}.spm.temporal.st.so = preproc_setup.so; % [2:2:32 1:2:31]

% reference slice
% usally the mid-time slice
matlabbatch{1}.spm.temporal.st.refslice = preproc_setup.refslice; % 32

% save batch file for review
slicetiming_mat_list = dir(fullfile(SubDataPath,'Step1_Slicetiming*.mat'));
if ~isempty(slicetiming_mat_list)
    mkdir(SubDataPath,'old_backup')
    time_mark = datestr(clock,30); %add time mark to the backup files
    slicetiming_mat_file = slicetiming_mat_list.name;
    [~,mat_name] = fileparts(slicetiming_mat_file);
    source_dir = fullfile(SubDataPath,slicetiming_mat_file);
    target_dir = fullfile(SubDataPath,'old_backup',[mat_name,'_backup_',time_mark,'.mat']);
    movefile(source_dir,target_dir);
end
save(fullfile(SubDataPath,['Step1_Slicetiming_',SubID,'.mat']),'matlabbatch');

% run batch
% spm('defaults','fmri');
% spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
clear matlabbatch
% move file to New Folders
for ii = 1:run_num
    run_id = runlist_rawnii(ii,1).name;
    sliced_file = fullfile(SubDataPath,run_id,'a*.nii');
    slice_dir = fullfile(Dir_Slice,run_id);
    mkdir(slice_dir);
    
    movefile(sliced_file,slice_dir)
end
end
function Sub_Realign(WorkDir,SubDataPath,SubID,preproc_setup)
clear matlabbatch
cd(SubDataPath);

% dir setting
Dir_Slice = fullfile(WorkDir,'FunImg_Step1_Slice',SubID);
Dir_Realign = fullfile(WorkDir,'FunImg_Step2_Realign',SubID);
Dir_RealignParameter = fullfile(WorkDir,'RealignParameter',SubID);
Dir_rp_Report = fullfile(WorkDir,'RealignParameter','rp_Report');
Dir_HeadMotionPlot = fullfile(WorkDir,'RealignParameter','HeadMotionPlot');

% gen the sliced img list in each run
runlist_slice = dir(fullfile(Dir_Slice,'run*'));
run_num = length(runlist_slice);
% loop the runs
for ii = 1:run_num

    run_id = runlist_slice(ii,1).name;
    % temporal corrected files has prefix 'a'
    % sub id is sub*, so the prefix should be 'asub'
    scanslist = dir(fullfile(Dir_Slice,run_id,'a*.nii'));
    scan_num = length(scanslist);
    scans = cell(scan_num,1);
    for kk = 1:scan_num
        scans{kk,1} = [fullfile(Dir_Slice,run_id,scanslist(kk).name),',1'];
    end
        
    matlabbatch{1}.spm.spatial.realign.estwrite.data{1,ii} = scans;
end

% PET Imgs set rtm 1 i.e. register to mean
% MRI typically be 0
% MRI images are typically registered to the first image. 
% The more accurate way would be to use a two pass  procedure
% but this probably wouldn't improve the results so much
% and would take twice as long
% preproc_setup.rgm
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = preproc_setup.rgm;% 1

% save batch file
step2_mat_list = dir(fullfile(SubDataPath,'Step2_Realign*.mat'));
if ~isempty(step2_mat_list)
    mkdir(SubDataPath,'old_backup')
    time_mark=datestr(clock,30);
    step2_mat_file = step2_mat_list.name;
    [~,mat_name] = fileparts(step2_mat_file);
    source_dir = fullfile(SubDataPath,step2_mat_file);
    dest_dir = fullfile(SubDataPath,'old_backup',[mat_name,'_backup_',time_mark,'.mat']);
    movefile(source_dir,dest_dir);
end

save(fullfile(SubDataPath,['Step2_Realign_',SubID,'.mat']),'matlabbatch');
% run batch
% spm('defaults','fmri');
% spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
clear matlabbatch;
% move file to New Folders
mean_file = fullfile(Dir_Slice,runlist_slice(1,1).name,'meana*.*');
movefile(mean_file,Dir_RealignParameter)

mean_img_temp = dir(fullfile(Dir_RealignParameter,'meana*.nii'));
mean_img_path = fullfile(Dir_RealignParameter,mean_img_temp.name);
    
for ii = 1:run_num
    run_id = runlist_slice(ii,1).name;
    % Calculate FD
    run_rp_list = dir(fullfile(Dir_Slice,run_id,'rp*.txt'));
    run_rp_path = fullfile(Dir_Slice,run_id,run_rp_list.name);
    [FD_VanDijk,FD_Jenkinson,FD_Power,HeadMotion] = FD_Calculate(run_rp_path,mean_img_path);
            
    file_id = [SubID,'_',run_id];
    save(fullfile(Dir_RealignParameter,['FD_VanDijk_',file_id,'.txt']),'FD_VanDijk','-ascii');
    save(fullfile(Dir_RealignParameter,['FD_Jenkinson_',file_id,'.txt']),'FD_Jenkinson','-ascii');
    save(fullfile(Dir_RealignParameter,['FD_Power_',file_id,'.txt']),'FD_Power','-ascii');
    save(fullfile(Dir_RealignParameter,['HeadMotion_',file_id,'.txt']),'HeadMotion');
    
% **************************** rp_Report [Start] **************************
    rp_data_temp = importdata(run_rp_path);
    rp_data_temp(:,4:6) = rp_data_temp(:,4:6)*180/pi; 
    rp_data_temp = abs(rp_data_temp);
%     [x,y] = find(rp_data_temp >= cutoff);
    % 3mm cutoff
    [x3,~] = find(rp_data_temp >= 3);
    if ~isempty(x3)
        trial_mark3 = unique(x3);
%         save(fullfile(Dir_rp_Report,'3mm_cutoff',[file_id,'_3mm.txt']),'trial_mark3','-ascii');
        txt_filename_3mm = fullfile(Dir_rp_Report,'3mm_cutoff',[file_id,'_3mm.txt']);
        dlmwrite(txt_filename_3mm,trial_mark3,'precision', '%.0f','delimiter','\t','newline','pc')
    end
    % 2mm_cutoff
    [x2,~] = find(rp_data_temp >= 2);
    if ~isempty(x2)
        trial_mark2 = unique(x2);
%         save(fullfile(Dir_rp_Report,'2mm_cutoff',[file_id,'_2mm.txt']),'trial_mark2','-ascii');
        txt_filename_2mm = fullfile(Dir_rp_Report,'2mm_cutoff',[file_id,'_3mm.txt']);
        dlmwrite(txt_filename_2mm,trial_mark2,'precision', '%.0f','delimiter','\t','newline','pc')
    end
    % 2mm_cutoff
    [x1,~] = find(rp_data_temp >= 1);
    if ~isempty(x1)
        trial_mark1 = unique(x1);
%         save(fullfile(Dir_rp_Report,'1mm_cutoff',[file_id,'_1mm.txt']),'trial_mark1','-ascii');
        txt_filename_1mm = fullfile(Dir_rp_Report,'1mm_cutoff',[file_id,'_3mm.txt']);
        dlmwrite(txt_filename_1mm,trial_mark1,'precision', '%.0f','delimiter','\t','newline','pc')
    end
 % **************************** rp_Report [End] ***************************   
 
    realigned_file = fullfile(Dir_Slice,run_id,'ra*.nii');
    rp_file = fullfile(Dir_Slice,run_id,'rp*.txt');
    
    run_realign_dir = fullfile(Dir_Realign,run_id);
    mkdir(run_realign_dir);
    
    movefile(realigned_file,run_realign_dir)
    movefile(rp_file,Dir_RealignParameter)
    
end

if matlabpool('size') <= 1
    % the ps file will be generated by spm,so the spm_graph should be open
    % if not , there will be error occur
    % convert the ps file to pdf 
    % which contain the headmotion and normalization check
    ps_file = dir(fullfile(SubDataPath,'*.ps'));
    movefile(fullfile(SubDataPath,ps_file.name),fullfile(Dir_HeadMotionPlot,[SubID,'_headmotion_plot.ps']))
end
end
function Sub_Normalise(WorkDir,SubDataPath,SubID,preproc_setup)
clear matlabbatch
cd(SubDataPath);
%  dir setting
Dir_Realign = fullfile(WorkDir,'FunImg_Step2_Realign',SubID);
Dir_Normalize = fullfile(WorkDir,'FunImg_Step3_Normalize',SubID);
Dir_RealignParameter = fullfile(WorkDir,'RealignParameter',SubID);
Dir_NormalizeCheck = fullfile(WorkDir,'RealignParameter','NormalizeCheck');

% gen realigned img file list in each run
runlist_realign = dir(fullfile(Dir_Realign,'run*'));
run_num = length(runlist_realign);
scans = [];
% specify images that need transformation
for ii = 1:run_num
    run_id = runlist_realign(ii,1).name;
    scanslist = dir(fullfile(Dir_Realign,run_id,'ra*.nii'));
    run_tr_num = length(scanslist);
    run_scans = cell(run_tr_num,1);
    for kk = 1:run_tr_num
        run_scans{kk} = [fullfile(Dir_Realign,run_id,scanslist(kk).name),',1'];
    end
    scans = [scans;run_scans];

end

% specify the transformation you want to apply
mean_img_temp = dir(fullfile(Dir_RealignParameter,'mean*.nii'));
mean_img_path = fullfile(Dir_RealignParameter,mean_img_temp.name);

matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {mean_img_path};

matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = scans;

% find the EPT Template location
EPI_Template_loc= [which('EPI.nii'),',1'];
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template = {EPI_Template_loc};
% set the bounding box
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = preproc_setup.bb;% [-90 -126 -72;90 90 108];
% set the voxel size
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [3,3,3];
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';

% save batch
normalize_mat_list = dir(fullfile(SubDataPath,'Step3_Normalise*.mat'));
if ~isempty(normalize_mat_list)
    mkdir(SubDataPath,'old_backup')
    time_mark = datestr(clock,30);
    normalize_mat_file = normalize_mat_list.name;
    [~,mat_name] = fileparts(normalize_mat_file);
    
    source_dir = fullfile(SubDataPath,normalize_mat_file);
    target_dir = fullfile(SubDataPath,'old_backup',[mat_name,'_backup_',time_mark,'.mat']);
    movefile(source_dir,target_dir);
end

save(fullfile(SubDataPath,['Step3_Normalise_',SubID,'.mat']),'matlabbatch');
% run batch
% spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch)

clear matlabbatch;
% move file to New Folders

for ii = 1:run_num
    run_id = runlist_realign(ii,1).name;
    normalize_file = fullfile(Dir_Realign,run_id,'wra*.nii');
    normalize_dir = fullfile(Dir_Normalize,run_id);
    mkdir(normalize_dir);
    movefile(normalize_file,normalize_dir)
end

if matlabpool('size') <= 1 % ParallelMode == 1
    % the ps file will be generated by spm,so the spm_graph should be open
    % if not , there will be error occur
    % convert the ps file to pdf 
    % which contain the headmotion and normalization check
    ps_file = dir(fullfile(SubDataPath,'*.ps'));
    movefile(fullfile(SubDataPath,ps_file.name),fullfile(Dir_NormalizeCheck,[SubID,'_normalize_check.ps']))
end

%************* Comput the mean normalized imgs in each run ****************
cd(Dir_Normalize);
% gen normalized img file list in each run
runlist_normalize = dir(fullfile(Dir_Normalize,'run*'));

for ii = 1:run_num%run num
  
    run_id = runlist_normalize(ii,1).name;
    scanslist = dir(fullfile(Dir_Normalize,run_id,'wra*.nii'));
    scan_num = length(scanslist);
    for kk = 1:scan_num
        run_scans = fullfile(Dir_Normalize,run_id,scanslist(kk).name);
        % this part is modified at 20140719
        % comput the mean normalized imgs in each run
        V_wmean = spm_vol(run_scans);
        [Y_wmean_temp,~] = spm_read_vols(V_wmean);
        if kk ==1
            Y_wmean_total = Y_wmean_temp;
        else
            Y_wmean_total = Y_wmean_total + Y_wmean_temp;
        end      
    end
    
    % write the wmean*.nii imgs, which is the mean normalized image
    Y_wmean = Y_wmean_total/scan_num;
    V_wmean.fname = [SubID,'_',run_id,'_wmean.nii'];
    spm_write_vol(V_wmean,Y_wmean);
    movefile(fullfile(Dir_Normalize,V_wmean.fname),Dir_NormalizeCheck);
end
%**************************************************************************

end
function Sub_Smooth(WorkDir,SubDataPath,SubID,preproc_setup)
clear matlabbatch list
cd(SubDataPath);

% dir setting
Dir_Normalize = fullfile(WorkDir,'FunImg_Step3_Normalize',SubID);
Dir_Smooth = fullfile(WorkDir,'FunImg_Step4_Smooth',SubID);

% gen normalized img file list in each run
runlist_normalize = dir(fullfile(Dir_Normalize,'run*'));
run_num = length(runlist_normalize);

scans = [];
for ii = 1:run_num

    run_id = runlist_normalize(ii,1).name;
    scanslist = dir(fullfile(Dir_Normalize,run_id,'wra*.nii'));
    
    run_tr_num = length(scanslist);
    run_scans = cell(run_tr_num,1);
    for kk = 1:run_tr_num
        run_scans{kk,1} = [fullfile(Dir_Normalize,run_id,scanslist(kk).name),',1'];     
    end
    scans = [scans;run_scans];
end

matlabbatch{1}.spm.spatial.smooth.data = scans;
% specify smoothing kernel
matlabbatch{1}.spm.spatial.smooth.fwhm = preproc_setup.fwhm; % [8 8 8]

% save batch
smooth_mat_list = dir(fullfile(SubDataPath,'Step4_Smooth*.mat'));
if ~isempty(smooth_mat_list)
    mkdir(SubDataPath,'old_backup')
    time_mark = datestr(clock,30);
    smooth_mat_file = smooth_mat_list.name;
    [~,mat_name,~] = fileparts(smooth_mat_file);
    
    source_dir = fullfile(SubDataPath,smooth_mat_file);
    target_dir = fullfile(SubDataPath,'old_backup',[mat_name,'_backup_',time_mark,'.mat']);
    movefile(source_dir,target_dir);
end

save(fullfile(SubDataPath,['Step4_Smooth_',SubID,'.mat']),'matlabbatch');

% run batch
spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch)
clear matlabbatch;
% move file to New Folders

for ii = 1:run_num
    run_id = runlist_normalize(ii,1).name;
    smooth_file = fullfile(Dir_Normalize,run_id,'swra*.nii');
    smooth_dir = fullfile(Dir_Smooth,run_id);
    mkdir(smooth_dir);
    movefile(smooth_file,smooth_dir)
end

% % generate the imgs to check normalization
% cd(Dir_NormalizeCheck)
% wmean_img_list = dir('*wmean.nii');
% wmean_img_num = length(wmean_img_list);
% 
% for kk = 1:wmean_img_num
%     wmean_img_path = fullfile(pwd,wmean_img_list(kk).name);
%     spm_image('init',wmean_img_path);
%     [pathstr, name, ~] = fileparts(wmean_img_path);
% %         export_fig(fullfile(pathstr, name),'-tif','-nocrop')% cannot be excuted in parallel
%     f=getframe(gcf);
%     imwrite(f.cdata,fullfile(pathstr, [name,'.jpg']));
% end
clear

end

%% ***************** HeadMotion Parameters Comput ******************
function [FD_VanDijk,FD_Jenkinson,FD_Power,HeadMotion] = FD_Calculate(rp_txt_file,mean_img_file)
%% function to calculate framewise displacement (FD) 
    RP = load(rp_txt_file);
    MaxRP = max(abs(RP));
    MaxRP(4:6) = MaxRP(4:6)*180/pi;

    MeanRP = mean(abs(RP));
    MeanRP(4:6) = MeanRP(4:6)*180/pi;

    %Calculate FD Van Dijk (Van Dijk, K.R., Sabuncu, M.R., Buckner, R.L., 2012. The influence of head motion on intrinsic functional connectivity MRI. Neuroimage 59, 431-438.)
    RPRMS = sqrt(sum(RP(:,1:3).^2,2));
    MeanRMS = mean(RPRMS);

    FD_VanDijk = abs(diff(RPRMS));
    FD_VanDijk = [0;FD_VanDijk];

    MeanFD_VanDijk = mean(FD_VanDijk);

    %Calculate FD Power (Power, J.D., Barnes, K.A., Snyder, A.Z., Schlaggar, B.L., Petersen, S.E., 2012. Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. Neuroimage 59, 2142-2154.) 
    RPDiff=diff(RP);
    RPDiff=[zeros(1,6);RPDiff];
    RPDiffSphere=RPDiff;
    RPDiffSphere(:,4:6)=RPDiffSphere(:,4:6)*50;
    FD_Power=sum(abs(RPDiffSphere),2);

    MeanFD_Power = mean(FD_Power);

    NumberFD_Power_05 = length(find(FD_Power>0.5));
    PercentFD_Power_05 = length(find(FD_Power>0.5)) / length(FD_Power);
    NumberFD_Power_02 = length(find(FD_Power>0.2));
    PercentFD_Power_02 = length(find(FD_Power>0.2)) / length(FD_Power);

    %Calculate FD Jenkinson (FSL's relative RMS) (Jenkinson, M., Bannister, P., Brady, M., Smith, S., 2002. Improved optimization for the robust and accurate linear registration and motion correction of brain images. Neuroimage 17, 825-841. Jenkinson, M. 1999. Measuring transformation error by RMS deviation. Internal Technical Report TR99MJ1, FMRIB Centre, University of Oxford. Available at www.fmrib.ox.ac.uk/analysis/techrep for downloading.)
 
    FD_Jenkinson = y_FD_Jenkinson(rp_txt_file,mean_img_file);
    MeanFD_Jenkinson = mean(FD_Jenkinson);


    HeadMotion = [MaxRP,MeanRP,MeanRMS,MeanFD_VanDijk,MeanFD_Power,NumberFD_Power_05,PercentFD_Power_05,NumberFD_Power_02,PercentFD_Power_02,MeanFD_Jenkinson];

end
function [rel_rms, abs_rms] = y_FD_Jenkinson(RealignmentParameterFile,ReferenceImage)
%% function [rel_rms, abs_rms] = y_FD_Jenkinson(RealignmentParameterFile,ReferenceImage)
% Calculate FD Jenkinson (relative RMS) and absolute RMS based on SPM's realignment parameters
% Reference: Jenkinson, M., Bannister, P., Brady, M., Smith, S., 2002. Improved optimization for the robust and accurate linear registration and motion correction of brain images. Neuroimage 17, 825-841.
%            Jenkinson, M. 1999. Measuring transformation error by RMS deviation. Internal Technical Report TR99MJ1, FMRIB Centre, University of Oxford. Available at www.fmrib.ox.ac.uk/analysis/techrep for downloading.
% Input:
% 	RealignmentParameterFile  -   The realignment parameter file for a given participant generated by SPM. E.g., rp***.txt
%   ReferenceImage            -   The reference image for realignment (usually the first time point (one-pass) or the mean image after an initial motion correction (two-pass))
% Output:
%	rel_rms      -   relative RMS (FD Jenkinson)
%	abs_rms      -   absolute RMS
%-----------------------------------------------------------
% Written by YAN Chao-Gan 120930.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com


rmax = 80.0; %The default radius (as in FSL) of a sphere represents the brain

RP=load(RealignmentParameterFile);
nTimePoint=size(RP,1);
sinq1=sin(RP(:,4));
sinq2=sin(RP(:,5));
sinq3=sin(RP(:,6));
cosq1=cos(RP(:,4));
cosq2=cos(RP(:,5));
cosq3=cos(RP(:,6));

[RefData, RefHead] = rest_ReadNiftiImage(ReferenceImage,1);
center = RefHead.mat*([0.5*(RefHead.dim(1));0.5*(RefHead.dim(2));0.5*(RefHead.dim(3));1]);
center = center(1:3); %Get the coordinate for the center

abs_rms = zeros(nTimePoint,1);
for t=1:nTimePoint

    M1=[1       0        0     0;...
        0    cosq1(t)  sinq1(t)  0;...
        0    -sinq1(t) cosq1(t)  0;...
        0       0        0     1;];
    
    M2=[cosq2(t)  0    sinq2(t)     0;...
        0        1       0        0;...
        -sinq2(t) 0    cosq2(t)     0;...
        0       0        0        1;];
    
    M3=[cosq3(t)   sinq3(t)   0     0;...
        -sinq3(t)  cosq3(t)   0     0;...
        0           0       1     0;...
        0           0       0     1;];
    
    MT=[1    0     0     RP(t,1);...
        0    1     0     RP(t,2);...
        0    0     1     RP(t,3);...
        0    0     0     1;];
    
    M_RigidBodyTransform=MT*M1*M2*M3;
    
    MA1=eye(4);
    MA2=(M_RigidBodyTransform);
    
%     M = MA1*inv(MA2) - eye(4);
    M = MA1/MA2 - eye(4);
    
    A = M(1:3,1:3);
    
    T = M(1:3,4);
    
    abs_rms(t) = sqrt(rmax*rmax/5*trace(A'*A) + (T+A*center)'*(T+A*center));
end


rel_rms = zeros(nTimePoint-1,1);
for t=2:nTimePoint
    M1=[1       0        0     0;...
        0    cosq1(t)  sinq1(t)  0;...
        0    -sinq1(t) cosq1(t)  0;...
        0       0        0     1;];
    
    M2=[cosq2(t)  0    sinq2(t)     0;...
        0        1       0        0;...
        -sinq2(t) 0    cosq2(t)     0;...
        0       0        0        1;];
    
    M3=[cosq3(t)   sinq3(t)   0     0;...
        -sinq3(t)  cosq3(t)   0     0;...
        0           0       1     0;...
        0           0       0     1;];
    
    MT=[1    0     0     RP(t,1);...
        0    1     0     RP(t,2);...
        0    0     1     RP(t,3);...
        0    0     0     1;];
    
    M_RigidBodyTransform=MT*M1*M2*M3;
    
    
    M1=[1       0        0     0;...
        0    cosq1(t-1)  sinq1(t-1)  0;...
        0    -sinq1(t-1) cosq1(t-1)  0;...
        0       0        0     1;];
    
    M2=[cosq2(t-1)  0    sinq2(t-1)     0;...
        0        1       0        0;...
        -sinq2(t-1) 0    cosq2(t-1)     0;...
        0       0        0        1;];
    
    M3=[cosq3(t-1)   sinq3(t-1)   0     0;...
        -sinq3(t-1)  cosq3(t-1)   0     0;...
        0           0       1     0;...
        0           0       0     1;];
    
    MT=[1    0     0     RP(t-1,1);...
        0    1     0     RP(t-1,2);...
        0    0     1     RP(t-1,3);...
        0    0     0     1;];
    
    M_RigidBodyTransform_1=MT*M1*M2*M3;
    
    MA1=(M_RigidBodyTransform_1);
    MA2=(M_RigidBodyTransform);
    
%     M = MA1*inv(MA2) - eye(4);
    M = MA1/MA2 - eye(4);
    
    A = M(1:3,1:3);
    
    T = M(1:3,4);
    
    rel_rms(t-1) = sqrt(rmax*rmax/5*trace(A'*A) + (T+A*center)'*(T+A*center));
    
end

rel_rms=[0;rel_rms]; %The FD_Jenkinson at time point t means the movement from time point t-1 to time point t. (Put the FD_Jenkinson for the first time point to "0".)
end
function NormalizeCheckImg(Dir_NormalizeMeanImg)
%% generate the imgs to check normalization use NormalizeMeanImg
    cd(Dir_NormalizeMeanImg)
    wmean_img_list = dir('*wmean.nii');
    wmean_img_num = length(wmean_img_list);

    for kk = 1:wmean_img_num
        wmean_img_path = fullfile(pwd,wmean_img_list(kk).name);
        spm_image('init',wmean_img_path);
        [pathstr, name, ~] = fileparts(wmean_img_path);
    %         export_fig(fullfile(pathstr, name),'-tif','-nocrop')% cannot be excuted in parallel
        f=getframe(gcf);
        imwrite(f.cdata,fullfile(pathstr, [name,'.jpg']));
    end
end
