%% function for conditionize the onset
% for subject in multiple run and multiple conditions
% 
% use when not all subjects share the same onset
% if all subjects share same onset,so you just need one subject info
% 
% written by hongshengcheng.math@gmail.com

clear;clc;
disp('running...please wait ...');

run_idx = 1; % which column is the run mark
sub_id_idx = 2; % the colunm index of sub id
tr_onset_idx = 3; % the colunm index of onset tr
cond_idx = 12; %12th column is the condition index


%% DO NOT MODIFY CODE BELOW
% get the raw onset xls
% which should have colums to identify run_mark,sub_id,cond_mark,onset
[xlsname,xlspath] = uigetfile('*.xls;*.xlsx','select data');
[~,~,raw] = xlsread(fullfile(xlspath,xlsname),1);

% get rid of the head information
raw_txt = raw(2:end,:);

% here,1st column is the run index,like run1,run2,run3,...
run_col = raw_txt(:,run_idx);
run_list = unique(run_col);
run_num = length(run_list);

% here,12th column is the condition index,like neut_t,neut_nt,nega_t,...
cond_col = raw_txt(:,cond_idx); %12
cond_list = unique(cond_col);
cond_num = length(cond_list);

for ii = 1:run_num
    data_cond = cell(cond_num,1);
    data_len = zeros(cond_num,1);
    for jj = 1:cond_num
        data_idx = strcmp(run_col,run_list(ii)).*strcmp(cond_col,cond_list(jj));
        % here the 2nd,3rd col is subid and tr
        % change the index as you need
        data_temp = raw_txt(logical(data_idx),[sub_id_idx tr_onset_idx]);
        % add head information
        data = [{'SubID',cond_list{jj}};data_temp];
        % if all subject have the same onset,then
        % data = [cond_list(jj);data_temp];
        data_cond{jj} = data;
        data_len(jj) = size(data,1);
    end
    % resize the onset
    data_merge = cell(max(data_len(:)),cond_num*2);
    for jj = 1:cond_num
        data_merge(1:data_len(jj),(jj-1)*2+1:(jj-1)*2+2) = data_cond{jj};
    end
    % write into xls
    xlswrite(fullfile(xlspath,['Onset_',xlsname]),data_merge,['run',int2str(ii)]);
end
%%
clear all;
clc;
disp('All work done')