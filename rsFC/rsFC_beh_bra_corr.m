function rsFC_beh_bra_corr(zfc_mat)
%% funtion to correlation the roi2roi FC and behavior

[beh_label,beh_data] = ReadXlsData('Input Behavior Data');
beh_num = length(beh_label);
% rsFC_BehData.beh_label = beh_label;
% rsFC_BehData.beh_data = beh_data;
save('rsFC_BehData','beh_label','beh_data')

corr_chos = input('1.corr; 2. partial_corr,Please enter the NO.: ');

if corr_chos == 1
    [beh_fc_r, beh_fc_p] = corr(beh_data,zfc_mat);
else
    [cov_label,cov_data] = ReadXlsData('Input Behavior Data');
    [beh_fc_r, beh_fc_p] = partialcorr(beh_data,zfc_mat,cov_data);
%     rsFC_BehData.cov_label = cov_label;
%     rsFC_BehData.cov_data = cov_data;
    save('rsFC_BehData','cov_label','cov_data','-append');
end

% save rsFC_BehData rsFC_BehData

    rsFC_Res = struct();

    for ii = 1:beh_num;beh_fc_r_3d(:,:,ii)=squareform(beh_fc_r(ii,:));end
    for ii = 1:beh_num;beh_fc_p_3d(:,:,ii)=squareform(beh_fc_p(ii,:));end
    rsFC_Res.beh_fc_r = beh_fc_r;
    rsFC_Res.beh_fc_p = beh_fc_p;
    rsFC_Res.beh_fc_r_3d = beh_fc_r_3d;
    rsFC_Res.beh_fc_p_3d = beh_fc_p_3d;
%     save('rsFC_Res','beh_fc_r','beh_fc_p','beh_fc_r_3d','beh_fc_p_3d');
    
%deviide the corr matrix into n*1
for beh_idx = 1:beh_num
    %find the sig p
    temp_p = beh_fc_p(beh_idx,:);
    
    %find the sig r
    temp_r = beh_fc_r(beh_idx,:);
    
    res_r = squareform(temp_r);
    res_p = squareform(temp_p);
    res_p = res_p + tril(ones(size(res_p)));
    
    [x,y] = find(res_p<0.05);
    sig_XY =[x,y];
    sig_p = res_p(res_p<0.05);
    sig_r = res_r(res_p<0.05);
    
    if ~isempty(sig_p)
        res_all = [sig_XY,sig_r,sig_p];
        rsFC_Res = setfield(rsFC_Res,beh_label{beh_idx},res_all);
        str_temp = strcat(beh_label{beh_idx},'$',num2str(res_all));
        cell_temp = cellstr(str_temp);
        cell2txt({cell_temp},'result')
    end

end
   save rsFC_Res  rsFC_Res
   corr_chos = input('would you want to correct the correlation? Y / N ','s');
   if strcmp(corr_chos,'y')||strcmp(corr_chos,'Y')==1
       roi_num = size(res_r,1);
       adjust_corr(beh_fc_p,roi_num);
   else
   end
end