% n*m format to (m*n)*1 format data
% first column is mark

[xlsname,xlspath] = uigetfile('*.xls;*.xlsx','select data');
[~,~,raw] = xlsread(fullfile(xlspath,xlsname),1);

num = raw(2:end,2:end);
[row_num,col_num] = size(num);

new_data = cell(numel(num),2);

raw_mark = raw(2:end,1);
sub_id = raw(1,2:end);

for ii = 1:col_num
    new_mark = strcat(sub_id{ii},'_',raw_mark);
    new_data((ii-1)*row_num+1:ii*row_num,1) = new_mark;
    new_data((ii-1)*row_num+1:ii*row_num,2) = num(:,ii);
end
xlswrite('output.xls',new_data);