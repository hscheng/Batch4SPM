function ExtractSig_Template()
% The raw name of this function is gretna_mean_tc, which is used to calculate the mean regional time course in ROIs.
% Raw Authors:Jinhui WANG, NKLCNL, BNU, BeiJing, 2011/10/23, Jinhui.Wang.1982@gmail.com
% ==========================================================================
% INPUT
% Selce Img List and Template
% OUTPUT
% folder named with the prefix of 'ExtractSig_Template'
% 
% modified by hongshengcheng.math@gmail.com
% batch script for extract template signal
%
% 20141109,modify select file into cell format list
% modified date 20140421

% select img list
ImgList = ScanImg2Cell('Select Imgs');
ImgNum = size(ImgList,1);

% choose the template
TemplateLoc = spm_select(1,'image','Select Template Img','',pwd,'.*',1);
% read template
Vtem = spm_vol(TemplateLoc);
[Ytem, XYZtem] = spm_read_vols(Vtem);
Ytem(isnan(Ytem)) = 0;

% set the roi index
ROI_Index = input('Please Enter the ROI Index (e.g. 1:1024 ): ');
% MNI_Coord store the ROI coordinates
MNI_Coord = cell(length(ROI_Index),1);
ROINum = length(ROI_Index);

for jj = 1:ROINum
    Region = ROI_Index(jj);
    
    Index = find(Ytem == Region);
    
    if ~isempty(Index)
        [I,J,K] = ind2sub(size(Ytem),Index);
        XYZ = [I J K]';
        XYZ(4,:) = 1;
        MNI_Coord{jj,1} = XYZ;
    else
        error (['There are no voxels in ROI' blanks(1) num2str(ROI_Index(jj)) ', please specify ROIs again']);
    end
%     Index=find(Ytem==Region);
%     XYZ=XYZtem(:,Index);
%     XYZ(4,:) = 1;
%     MNI_coord{jj,1} = XYZ;
end

% Mean_ROI to store ROI signal
Mean_ROI = cell(ImgNum,1);

% http://blog.sina.com.cn/s/blog_4d633dc70100nwcw.html
hwait = waitbar(0,'Please Wait ...');
step = ImgNum/100;

for ii = 1:ImgNum
       
    Vin = spm_vol(ImgList{ii});
    
    Mean_sig = zeros(size(Vin,1),length(ROI_Index));
    
    for kk = 1:length(ROI_Index)
        VY = spm_get_data(Vin,MNI_Coord{kk,1});
        Mean_sig(:,kk) = mean(VY,2);
    end
    
    Mean_ROI{ii,1} = Mean_sig;
    
    % wait bar
    PerStr = fix(ii/step);
    str = ['Running ...',num2str(PerStr),'% Completed'];
    waitbar(ii/ImgNum,hwait,str);
    pause(0.05);
end
close(hwait);

Mean_Signal=cell2mat(Mean_ROI(:,1));
Mean_Signal(isnan(Mean_Signal))=0;

% ExtractInfo struct to store some infomation for check
ExtractInfo.dataset = ImgList;
ExtractInfo.img_num = ImgNum;
ExtractInfo.template = TemplateLoc;
ExtractInfo.roi_index = ROI_Index;

disp('Saving the result, Please wait ....')

% add a time mark to the res_mat
time_mark = datestr(clock,30);
ResDir = ['ExtractSig_Template_',time_mark];
mkdir(ResDir);
cd(ResDir);

save(['ExtractInfo','_',time_mark,'.mat'],'ExtractInfo');
save (['Mean_Signal','_',time_mark,'.mat'],'Mean_Signal');

%generate a txt file which contain the signal
txt_filename=['Mean_Signal_',int2str(length(ROI_Index)),'_ROIs_',time_mark,'.txt'];
dlmwrite(txt_filename,Mean_Signal,'delimiter', '\t','newline','pc','precision',6);

clear;clc
msgbox('All work done! Please check res files','congradulation!')
end