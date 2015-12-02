function img_list_cell = ScanImg2Cell(varargin)
%% function to scan img files and reture the n*1 cell format list
% ScanImg2Cell(dlg_title)
% written by hongshengcheng.math@gmail.com
% INPUT
% Select Img List
% OUTPUT
% N*1 Format cell img list
% 20141130,suit for multikind files

    if numel(varargin) == 1
        dlg_title = varargin{1};
        filetype = 'image';
    elseif numel(varargin) == 2
        dlg_title = varargin{1};
        filetype = varargin{2}; %'any','mat','image'
    else
        dlg_title = 'Select Img Files';
        filetype = 'image';
    end
    img_list = spm_select(Inf,filetype,dlg_title,'',pwd,'.*',1);
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