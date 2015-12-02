function Img_4Dto3D()
%% this function can transform 4d nii files to 3d nii files
% Note:for SPM8 use only
%
% spm_file_part
% Convert a 4D volume file into a series of 3D volume files
% FUNCTION Vo = spm_file_split(V, odir)
% V           - filename or spm_vol struct
% odir        - output directory [default: same as input]
% 
% Vo          - spm_vol struct array of output files
% 
% written by hongshengcheng.math@gmail.com
% 20150523,add data organization type


ImgList_4D = ScanImg2Cell('Select 4D Img Files');
ResDir = fullfile(pwd,'Converted_3D_files');
mkdir(ResDir)

% Check the name be unique or not
data_org = questdlg('data organization','quest',...
    'all 4D files in one folder','single 4D file in each subject',...
    'all 4D files in one folder');

if strcmp(data_org,'all 4D img files in one folder')
    mark_name = 1;% it means can use the img filename directly
else
    mark_name =2; % use the subject dir name
end
    
for ii = 1:length(ImgList_4D)
    if mark_name == 1
        [~,ImgName,~] = fileparts(ImgList_4D{ii});
        OutputDir = fullfile(ResDir,ImgName);
    else
        [pathstr,~,~] = fileparts(ImgList_4D{ii});
        [~,subdir,~] = fileparts(pathstr);
        OutputDir = fullfile(ResDir,subdir);
    end
    
    mkdir(OutputDir);
    %default outputdir is same as input
    spm_file_split(ImgList_4D{ii,1},OutputDir);

end
clear;clc;
msgbox('All Work Done!','Note')
end