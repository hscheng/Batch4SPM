function ExtractSig_Voxels()
%% function to extract voxels values besed on mask
% Written by hongshengcheng.math@gmail.com
% 
% INPUT
% Select Imgs to extract and mask
% OUTPUT
% a folder named Voxel_Signal_+time which stored the result files
% 
% 20141109,modify a bug to read img file
% modified date 20140401

% select img list
ImgList = ScanImg2Cell('Select Imgs');
ImgNum = size(ImgList,1);

% choose the mask file,mask must be choosed
MaskLoc = spm_select(1,'image','Select Mask Img','',pwd,'.*',1);

% read the mask
Vmask = spm_vol(MaskLoc);
[Ymask,XYZmask] = spm_read_vols(Vmask);
[Xn_mask,Yn_mask,Zn_mask] = size(Ymask);

% read one img to check the dimension match or not
V = spm_vol(ImgList{1});
TR_num = size(V,1);
fprintf('Your img file have %d TRs\n',TR_num);
Y = spm_read_vols(V);
[xn,yn,zn] = size(Y);

% check the dimension match or not
if xn==Xn_mask && yn==Yn_mask && zn/TR_num==Zn_mask
    Index = find(Ymask~=0);
    if ~isempty(Index)
        [I,J,K] = ind2sub(size(Ymask),Index);
        XYZ_idx = [I J K]';
        XYZ_idx(4,:) = 1;
    end
else
    errordlg('Mask File and Img Files dimension not match','File Error')
end

Voxel_Signal = cell(ImgNum,1);

% add a wait bar
hwait = waitbar(0,'Please Wait ...');
step = ImgNum/100;

%read each img file
for ii = 1:ImgNum
    
    Vin = spm_vol(ImgList{ii});
    
    fprintf('Extracting time series for %s\n',ImgList{ii});
    
    SigTemp = spm_get_data(Vin, XYZ_idx);
    
    Voxel_Signal{ii,1} = SigTemp;
    
    disp('Done!')
    
    % wait bar
    PerStr = fix(ii/step);
    caption_str = ['Running ...',num2str(PerStr),'%'];
    waitbar(ii/ImgNum,hwait,caption_str);
    pause(0.05);
    
end
close(hwait);% close the wait bar
Voxel_Signal_All_Subs = cell2mat(Voxel_Signal(:,1));

% ExtractInfo struct to store some infomation for check
ExtractInfo.dataset = ImgList;
ExtractInfo.sub_num = ImgNum;
ExtractInfo.mask = MaskLoc;

% add a time mark to the res_mat
disp('Saving the result, Please wait ....')
time_mark = datestr(clock,30);
ResDir = ['ExtractSig_Voxels_',time_mark];
mkdir(ResDir);
cd(ResDir);
save(['ExtractInfo','_',time_mark],'ExtractInfo');
save(['Voxel_Signal_',time_mark],'Voxel_Signal_All_Subs');

% generate a txt file which contain the signal
txt_filename = ['Voxel_Signal_',time_mark,'.txt'];
dlmwrite(txt_filename,Voxel_Signal_All_Subs,'delimiter', '\t','newline','pc','precision',6);

clear;clc;
msgbox('All work done! Please check your res files','congradulation!')

end