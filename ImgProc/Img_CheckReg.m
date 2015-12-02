function Img_CheckReg()
% function to check image registration use spm_check_registration
% check results will be saved in jpg format

% [filename,pathname] = uigetfile('*.txt','Select the file list in txt format');
% FileList = importdata(fullfile(pathname,filename));
FileList = ScanImg2Cell();
ResDir = 'CheckRegRes';
mkdir(ResDir)
for ii = 1:length(FileList)
    FilePathTemp = FileList{ii};
    [pathstr,FileNameTemp,~] = fileparts(FilePathTemp);
    [~,sub_dir] = fileparts(pathstr);
    spm_check_registration(FilePathTemp);
	% save result in jpg format
    f=getframe(gcf);
    imwrite(f.cdata,fullfile(ResDir, [sub_dir,FileNameTemp,'.jpg']));
    close all;

%     mkdir(FileNameTemp);
%     copyfile(FilePathTemp,FileNameTemp)
end

msgbox('All work Done!',':)');
end
function img_list_cell = ScanImg2Cell(varargin)
% scan img files and reture the cell format list
% 
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