% function to check the SPM T brain map 
% written by Hongsheng-Cheng | hongshengcheng.math@gmail.com
% v3:20141102 fit for AhongToolkit
% modified date: 20140721
% created date:20140716

clear;clc;


prompt = {'Threshold','Name of OutputDir','Prefix of OutPut'};
dlg_title = 'Input';
num_lines = 1;
def = {'2.58','SPM_T_Check_Results','CR_'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
% parameters can be changed
threshold = str2num(answer{1});
dir_output = answer{2};
prefix_outimg = answer{3};

% make a folder to store results
if exist(fullfile(pwd,dir_output),'dir')~=7
    mkdir(dir_output);
end

% [list_filename,list_filepath] = uigetfile('*.txt','Select the file list');
% file_list = importdata(fullfile(list_filepath,list_filename));

ImgList = ScanImg2Cell;
file_num = size(ImgList,1);

[Mask_filename,Mask_filepath] = uigetfile('*.nii;*.img','Select Brain Mask');
Mask_fullpath = fullfile(Mask_filepath,Mask_filename);
V_mask = spm_vol(Mask_fullpath);
[Y_mask,XYZ_mask] = spm_read_vols(V_mask);

% P_005 = 2.58;
% P_001 = 3.09;

% Check the name be unique or not
% [pathstr1,ImgName1] = fileparts(ImgList{1});
% [pathstr2,ImgName2] = fileparts(ImgList{2});
% if strcmp(ImgName1,ImgName2) ~= 1
%     mark_name = 1;% it means can use the img filename directly
% elseif strcmp(pathstr1,pathstr2) ~= 1
%     mark_name =2;
% else error('please check your image file');
% end
    
Voxel_Num_Report = zeros(file_num,2);
for ii = 1:file_num   
%     if mark_name == 1
%             [~,ImgName,~] = fileparts(ImgList{ii});
%             output_filename = [prefix_outimg,ImgName,'.jpg'];
%         else
%             [pathstr,~,~] = fileparts(ImgList{ii});
%             [~,subdir,~] = fileparts(pathstr);
%             output_filename = [prefix_outimg,subdir,'.jpg'];
%     end
% %%
    [pathstr,ImgName,~] = fileparts(ImgList{ii});
    [~,subdir,~] = fileparts(pathstr);
    output_filename = [prefix_outimg,subdir,'_',ImgName,'.jpg'];

    xjview(ImgList{ii});
    f = getframe(gcf);
    imwrite(f.cdata,fullfile(dir_output, output_filename));
    close all;

%%
    V = spm_vol(ImgList{ii});
    [Y,XYZ] = spm_read_vols(V);
    Y_Masked = Y.*Y_mask;
    Posi_Voxel_Num = size(find(Y_Masked > threshold),1);
    Nega_Voxel_Num = size(find(Y_Masked < -threshold),1);
    Voxel_Num_Report(ii,:) = [Posi_Voxel_Num Nega_Voxel_Num];
    
end
save Voxel_Num_Report Voxel_Num_Report
% save('Voxel_Num_Report.txt','Voxel_Num_Report','-ascii'); %20140721
dlmwrite('Voxel_Num_Report.txt',Voxel_Num_Report,'precision', '%.0f','delimiter','\t','newline','pc')
clear;clc;
disp('All work done!');