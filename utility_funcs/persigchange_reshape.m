sub_num = 30;
run_num = 6;
cond_num =4;
MatData = zeros(sub_num,run_num*cond_num);
for ii = 1:sub_num
    data_temp = PerSigChange_All{ii,1};
    idx = length(data_temp);
    MatData(ii,1:idx)= data_temp;
end