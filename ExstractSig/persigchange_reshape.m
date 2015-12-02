% load the percent signal change result mat first

cond_num = 4;
run_num = 6;

%% DO NOT Modify code below
col_num = cond_num*run_num;

sub_num = size(PerSigChange_All,1);

MatData = zeros(sub_num,col_num);

for ii = 1:sub_num
    % you should load the data mat first
    data_temp = PerSigChange_All{ii,1};
    idx = length(data_temp);
    MatData(ii,1:idx)= data_temp;
end