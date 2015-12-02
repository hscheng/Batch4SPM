function rsFC_roi2roi()
% function rsFC_roi2roi(roi_num,sub_num,DataFullPath)
% written by hongshengcheng.math@gmail.com
%
% 20141124,fix some error
% 20140611,merge the functions in one script
% 20131210,add a time mark to the result files
% 20131209,modified a bug that can save 0 in a corr mat
%
% this function can get the triu matrix of the functional connectivity of ROI2ROI
% there are 2 functions, one for ttest between 2 groups; another is to corr the behaviour data with the FC 

% before runing this script,you may should prepare three file:1.raw roi2roi fcanalysis results;
% 2.the behaviour data;3.covrariate data(if you want to do a partial corr between behaviour and fc)
% all the 3 files should be stroed in mat or txt format
    clear;clc;

    
%     [path_filename,path_filepath] = uigetfile({'*.txt;*.mat'},'Please choose the fMRI data corr path file( in mat/txt format)');
%     zFC_FileList = importdata(fullfile(path_filepath,path_filename));
    zFC_FileList = ScanImg2Cell('Slect zFC txt Files','mat');
    file_num = length(zFC_FileList);
    
    % merge the correlation marix into 2d and 3d mat
    for ii = 1:file_num
        data_temp = importdata(zFC_FileList{ii});
        if ii == 1
            roi_num = size(data_temp,1);
            % initial the mat to store res
            rsFC_Data_3D = zeros(roi_num,roi_num,file_num);
            rsFC_Data_2D = zeros(file_num,roi_num*(roi_num-1)/2);
        end
        % add eps
        data_temp = data_temp+eps;

        %get the uptri matrix which diag not include, then rotate the matrix
        data_triu = triu(data_temp,1)';

        %line the data like[1,2] [1,3]...[1,n],[2,3],[2,4]...[2,n]...
        mat_2d_temp = data_triu(data_triu~=0)'; 

        %merge a 3d matrix
        rsFC_Data_3D(:,:,ii) = data_temp-eps;
        rsFC_Data_2D(ii,:) = mat_2d_temp-eps;
    end

    % add a time mark to the res dir
    dir_res = ['rsFC_stat_',datestr(clock,30)];
    mkdir(dir_res);
    cd(dir_res)

    
    % save data information
%     rsFC_Data.rsFC_Data_3D = rsFC_Data_3D;
%     rsFC_Data.rsFC_Data_2D = rsFC_Data_2D;
%     rsFC_Data.filelist = zFC_FileList;
%     save rsFC_Data rsFC_Data
    save('rsFC_Data','rsFC_Data_3D','rsFC_Data_2D','zFC_FileList');
    
    % further stat 
    % Construct a questdlg with three options
    choice = questdlg('What would you like to do?', ...
        'rsFC_roi2roi', ...
        '2-group T-test','behavior&brain corr','No thank you','No thank you');
    % Handle response
    switch choice
        case '2-group T-test'
            group_ttest(rsFC_Data_2D);
        case 'behavior&brain corr'
            beh_bra_corr(rsFC_Data_2D);
        case 'No thank you'
            close
    end

end

function beh_bra_corr(zfc_mat)
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

function adjust_corr(corr_beh_bra_p,roi_num)
%% function to adjust the correlation between FC and behavior
   ad_chos=input('choose the type of correlation?1.fdr 2. fdr-2-step 3. bonferroni Enter the num : ');
   ad_p=input('set the p value(e.g.0.05,0.01): ');
   switch ad_chos
       case 1
           disp('fdr correction  ref:Benjamini & Hochberg (1995) and the Benjamini & Yekutieli (2001) procedure');

           met_chos=input('what kind of method do you want to use?  1.pdep 2.dep Enter the num : ');
           if met_chos==1
               ad_method='pdep';
           else
               ad_method='dep';
           end
           [h, crit_p, adj_p]=fdr_bh(corr_beh_bra_p,ad_p,ad_method,'yes')
       case 2
           disp('fdr 2-step correction ref: "two-stage" Benjamini, Krieger, & Yekutieli (2006)');
           [h, crit_p]=fdr_bky(corr_beh_bra_p,ad_p,'yes');
       case 3
           disp('bonferroni corrction ref :Bonferroni-Holm (1979)');
           [crit_p, h]=bonf_holm(corr_beh_bra_p,ad_p)
       otherwise
           %nothing to do
   end
   index = h==1
   res_h=zeros(roi_num);
   k=1;
   j=roi_num-2;
   for i=1:roi_num-1
        res_h(i,i+1:roi_num)=h(k:k+j);
        k=k+j+1;
        j=j-1;
   end
    save('res_h','res_h');
end

function ttest_res = group_ttest(zfc_2d_mat)
%% [func]ttest between 2 groups roi2roi functional connectivtiy
    %last modified date: 20140611
   
    % define two groups by enter the index
    grp1_idx = input('input the group 1 data index, eg,1:12  :');

    % devide the data into 2 groups
    grp1_data = zfc_2d_mat(grp1_idx,:);
    grp2_data = zfc_2d_mat(grp1_idx(end)+1:end,:);
    
    paired_mark = input('Paired Data?1-Yes,2-No :' );
    if paired_mark == 1
        [h, p] = ttest(grp1_data,grp2_data);
    else
        [h, p] = ttest2(grp1_data,grp2_data);
    end

    % change the shape to roi_num*roi_num matrix
    res_h = squareform(h);
    res_p = squareform(p);
        
    %h_p=cat(3,res_h,res_p); 
    [x,y] = find(res_h==1);
    sig_h = [x,y];
    sig_p = res_p((res_h==1));
    ttest_res = [sig_h,sig_p];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %check the four outputs for your results
    %save('h_p','h_p');%you just need check these 2 results
    save('res_p','res_p');
    save('res_h','res_h');
    save('ttest_res','ttest_res');
end