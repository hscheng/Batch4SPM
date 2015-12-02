function [Labels,DataMat] = ReadXlsData(varargin)
%% function to copy data direct from xls file
% make sure the end of the data should be cleaned

    if isempty(varargin)
        dlg_title = 'Input Data from Excel';
    else
        dlg_title = varargin{1};
    end
    
    prompt = {'Enter Labels :','Enter Matrix :'};
    Paras = inputdlg(prompt,dlg_title,[1 100;10 100]);

    LabelRaw = Paras{1,1};
    space_mark = strfind(LabelRaw,' ');
    if isempty(space_mark)
        Labels = regexp(LabelRaw,'\t','split');
    else
        Labels = regexp(LabelRaw,' ','split');
    end
    Labels = Labels'; % n*1 cell mat
    
    DataRaw = Paras{2,1};
    DataMat = str2num(DataRaw);
%     if read raw data into one line, script follow will be fine
%     DataStr = regexp(DataRaw,' ','split');
%     RowNum = size(DataStr,2);
%     ColNum = size(Labels,1);
%     DataMat = zeros(RowNum,ColNum);
%     for ii = 1:RowNum
%         DataTempStr = regexp(DataStr{ii},'\t','split');
%         for jj = 1:ColNum
%             DataMat(ii,jj) = cellfun(@str2num,DataTempStr(jj));
%         end
%     end
end