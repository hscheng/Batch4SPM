function ExtractPerSigChange(varargin)
% function ExtractPerSigChange(dur,ROI_List,SubMatList)

% cond_num = 4; % how many conditions in your experiment

if numel(varargin)<3
    % duration in seconds of event to estimate for
    % dur = 4;% maybe you want to set it to 0,ref function event_signal
    dur = input('duration(in seconds : ');
    % Select roi mat files
    ROI_List = ScanImg2Cell('Slecet ROI Img files'); % setting arrary for roi file name.
    % select the spm mat list in txt format
    SubMatList = ScanImg2Cell('Slecet SPM mat files','mat');
else
    dur = varargin{1};
    ROI_List = varargin{2};
    SubMatList = varargin{3};
end

ROI_Num = length(ROI_List);
SubNum = length(SubMatList);

% Initial mat to store result
Beta_All = cell(SubNum,1);
PerSigChange_All = cell(SubNum,1);

for roi_idx = 1:ROI_Num   
    roi_mat = nii2mat(ROI_List{roi_idx});
    for sub_idx =1:SubNum; % subjects' number for loop
        spm_mat = SubMatList{sub_idx}; % select subject's design
        [sub_persigchange,sub_beta] = ComputPerSigChange(spm_mat,roi_mat,dur);
        Beta_All{sub_idx} = sub_beta;
%         sub_persigchange = ComputPerSigChange(spm_mat,roi_mat,dur);
        PerSigChange_All{sub_idx} = sub_persigchange;
     end

% save beta values & %signal change
format short;
[~,roi_name,~] = fileparts(roi_mat);
% dlmwrite(['beta_',roi_name,'.txt'],Beta_All,'delimiter','\t','precision',6);
% dlmwrite(['Per_Signal_Change_',roi_name,'.txt'],PerSigChange_All,'delimiter','\t','precision',6);
% 
% % xlswrite is time-consuming, you can delete the following line.
% xlswrite(['beta_',roi_name,'.xls'],Beta_All,1); 
% xlswrite(['Per_Signal_Change_',roi_name,'.xls'],PerSigChange_All,1); 

% save(['beta_',roi_name,'.mat'],'Beta_All');
save(['PerSignalChange_',roi_name,'.mat'],'PerSigChange_All');
end

msgbox('All Work Done !!','ALL DONE');
clc;
end

% function PerSigChange = ComputPerSigChange(spm_mat,roi_mat,dur)
function [PerSigChange,Beta,marsS] = ComputPerSigChange(spm_mat,roi_mat,dur)
%% http://marsbar.sourceforge.net/faq.html#how-is-the-percent-signal-change-calculated

% dur = 2;
% spm_mat = 'E:\TNT_RawData\DataProc_Task\Stat_1stLevel\sub02\SPM.mat';
% roi_mat = 'ROI_-27_-27_-6_6mm.mat';

    % Make marsbar design object
    D = mardo(spm_mat);
    % Mkae marsbar ROI object
    R = maroi(roi_mat);
    % Fetch data into marsbar data object, 
    Y = get_marsy(R,D,'mean');
    % Get contrasts from original design
    xCon = get_contrasts(D);
    % Estimate design on ROI data
    E = estimate(D,Y);
    % Put contrasts from original design back into design object
    E = set_contrasts(E,xCon);
    % Get desing betas
    Beta = betas(E); 
    % Get stats and stuff for all contrasts into statistics structure
    marsS = compute_contrasts(E, 1:length(xCon)); 
    
%     % return beta for all events in design
%     for conditions = 1:cond_num;
%         Beta(conditions) = Beta(conditions,1);
%     end

%% extract percent signal change from design
    %get definitions of all events in model
    [e_specs, e_names] = event_specs(E);
    n_events = size(e_specs,2);
    % dur = 0;
    % return percent signal esimate for all events in design
    pct_ev = zeros(1,n_events);
    for e_s = 1:n_events
        pct_ev(e_s) = event_signal(E,e_specs(:,e_s),dur);
    end
%     fprintf('n_events is %d\n',n_events);
    PerSigChange = pct_ev;
    
%% extract all the FIR timecourses from design
%     % Bin size in seconds for FIR
%     bin_size = tr(E);
%     % Length of FIR in seconds
%     fir_length = 24;
%     % Number of FIR time bins to cover length of FIR
%     bin_no = fir_length / bin_size;
%     % Options - here 'single' FIR model, return estimated
%     % http://marsbar.sourceforge.net/apidocs/marsbar/@mardo/event_fitted_fir.html
%     opts = struct('single', 1, 'percent', 1);
%     % Return time courses for all events in fir_tc matrix
%     for e_s = 1:n_events
%       fir_tc(:, e_s) = event_fitted_fir(E, e_specs(:,e_s), bin_size,bin_no, opts);
%     end
%     FIR_TimeCourse = fir_tc;
    
% If your events have the same name across sessions
% and you want to average across the events with the same name:
%     % Get compound event types structure
%     ets = event_types_named(E);
%     n_event_types = length(ets);
%     for e_t = 1:n_event_types
%        fir_tc(:, e_t) = event_fitted_fir(E, ets(e_t).e_spec, bin_size,bin_no, opts);
%     end
end

function TimeCourse = GetRawTC(spm_mat,roi_mat)
    rois = maroi('load_cell', roi_mat); % make maroi ROI objects
    des = mardo(spm_mat);  % make mardo design object
    mY = get_marsy(rois{:}, des, 'mean'); % extract data into marsy data object
    TimeCourse  = summary_data(mY);  % get summary time course(s)
end

function mat_file = nii2mat(img_file)
% nii_file should be binary img file
% http://sourceforge.net/p/marsbar/mailman/message/31299330/

    % img_file = 'my_image.nii';
    [pathstr,filename] = fileparts(img_file);
    output = maroi_image(struct('vol', spm_vol(img_file), 'binarize',0,'func', 'img'));
    % saveroi(output, 'my_image_roi.mat')
    mat_file = fullfile(pathstr,[filename,'.mat']);
    saveroi(output,mat_file);
end