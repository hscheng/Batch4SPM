function [SubDir_fullpath,SubFolders] = MergeDir(varargin)
%% function to merge subdir full path and get the subfolder names list
% [SubDir_fullpath,SubFolders] = MergeDir(parent_dir,'part')
% mode  'part',default,user need to choose subfolders manually
%       'all',auto select all subfolders

% written by hongshengcheng.math@gmail.com
% 20150411,add a mark to select part of subfolders or all
% 20141201,add listdlg to select dir

    switch numel(varargin)
        case 0
            parent_dir = uigetdir(pwd,'select data dir');
            mode = 'part';
        case 1
            parent_dir = varargin{1};
            mode = 'part';
        case 2
            parent_dir = varargin{1};
            mode = varargin{2};
    end
    
    
    TempList = dir(parent_dir); % get the dir struct
    ISub = [TempList(:).isdir]; % filter dir
    folders_list = {TempList(ISub).name}'; % get all dir name
    folders_list(1:2) = []; %remove . and ..
    
    % merge parent_dir and subfolders to fullpath
    if ispc % OS is windows
        fullpath_list = strcat(parent_dir,'\',folders_list);
    else % linux
        fullpath_list = strcat(parent_dir,'/',folders_list);
    end
    
    if strcmp(mode,'part')==1 %if we just need part of subfolders
        dir_idx = listdlg('PromptString','Select Data Dirs',...
                    'ListSize',[300 300],...
                    'SelectionMode','multiple',...
                    'ListString',fullpath_list,...
                    'InitialValue',[1:length(fullpath_list)]);
        SubDir_fullpath = fullpath_list(dir_idx);
        SubFolders = folders_list(dir_idx);
    else % all subfolders
        SubDir_fullpath = fullpath_list;
        SubFolders = folders_list;
    end
end