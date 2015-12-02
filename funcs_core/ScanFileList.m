function PathinCell = ScanFileList(varargin)
%% function to copy data direct from xls file
% make sure the end of the data should be cleaned
% written by hongshengcheng.math@gmail.com
if isempty(varargin)
    dlg_title = 'Select Files';
else
    dlg_title = varargin{1};
end

choice_img = questdlg('Are Files in one folder ?',dlg_title,'Yes','No','No');

    switch choice_img
        case 'Yes'
            FileFilter = input('File Filter (e.g.*.mat): ','s');
            if isempty(FileFilter)
                FileFilter = '*.*';
            end
            % modified date 20140828
            % Raw Imgs refer to the valued template in each parameters
            [RawImgFileName,RawImgFilePath] = uigetfile(FileFilter,dlg_title,'MultiSelect','on');
            RawImgFullPath = strcat(RawImgFilePath,RawImgFileName);
            if iscell(RawImgFullPath)
                PathinCell = RawImgFullPath';
            else
                PathinCell = {RawImgFullPath};
            end
        case 'No'
            Paras = inputdlg({'Enter fullpath list :'},'Input path list',[10 100]);
            PathListRaw = Paras{1,1};
            
            FileNum = size(PathListRaw,1);
            PathinCell = cell(FileNum,1);
            
            for ii = 1:FileNum
                % modified 20141012
                % fix a bug which will generate extra space
                str_temp = PathListRaw(ii,:);
                str_ok = str_temp(~isspace(str_temp));
                PathinCell{ii} = str_ok;
            end
    end