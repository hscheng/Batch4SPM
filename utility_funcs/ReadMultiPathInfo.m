function [grp_names,para_names,para_cell] = ReadMultiPathInfo()
%% function to read multiple path info
% Input
% Select Path.xls file
% Output
% grp_names,para_names,para_cell

    clear;clc;

    [name,path] = uigetfile('*.xls;*.xlsx','Select Path file');
    filepath = fullfile(path,name);

    [num,txt,raw] = xlsread(filepath);

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