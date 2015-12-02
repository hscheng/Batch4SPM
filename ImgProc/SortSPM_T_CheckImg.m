% function SortSPM_T_CheckImg()

clear;clc;
ResDir = 'SortedImgs';
mkdir(ResDir);
spmT_dir = uigetdir(pwd,'spm T imgs');
spmT_idx = input('spmT Index e.g.[1:14] : ');
spmT_tag = reshape(sprintf('%04d',spmT_idx),4,[]).';
spmT_num = length(spmT_idx);
for ii = 1:spmT_num
    spmT_tag_temp = spmT_tag(ii,:);
    spmT_dir_sort = fullfile(ResDir,['spmT_',spmT_tag_temp]);
    mkdir(spmT_dir_sort);
    spmT_img = fullfile(spmT_dir,['*spmT_',spmT_tag_temp,'.jpg']);
    movefile(spmT_img,spmT_dir_sort);
end

clear;clc;
disp('All Work Done')