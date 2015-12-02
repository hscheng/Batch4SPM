function [grp_names,para_names,para_cell] = ReadMultiPathInfo(varargin)
%% function to read multiple path info
% [grp_names,para_names,para_cell] = ReadMultiPathInfo('parapath.xls')
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