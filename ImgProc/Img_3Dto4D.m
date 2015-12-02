function Img_3Dto4D()
% batch to convert 3d img/hdr pair files to 4d nii file
% also can convert hdr/img to nii
% written by Hongsheng-Cheng | hongshengcheng.math@gmail.com
% created date:20140725

Raw3D_dir = uigetdir(pwd,'Select 3D Img File Dir');
[SubDir_fullpath,SubFolders] = MergeDir(Raw3D_dir);

DestDir = uigetdir(pwd,'Select Dest Dir to store Merged file');

SubNum = size(SubFolders,1);

for ii = 1:SubNum
    SubDir = SubDir_fullpath{ii,1};
    cd(SubDir);
    img_list = dir(fullfile(SubDir,'*.img'));
    img_num = size(img_list);
    ImgList_3D = cell(img_num,1);

    for jj = 1:img_num
        ImgList_3D{jj,1} = [fullfile(SubDir,img_list(jj,1).name),',1'];
    end
    FileName_4D = [SubFolders{ii,1},'_4d.nii'];
    Write4DNii(ImgList_3D,FileName_4D)
    
    SubDirNew = fullfile(DestDir,SubFolders{ii,1});
    mkdir(SubDirNew);
    movefile(FileName_4D,fullfile(SubDirNew,FileName_4D));
end
end

function Write4DNii(ImgList,OutputFileName)        
    matlabbatch{1,1}.spm.util.cat.vols = ImgList; %cell n*1;
    matlabbatch{1,1}.spm.util.cat.name = OutputFileName;%output 4d img filename
    matlabbatch{1,1}.spm.util.cat.dtype = 4;%16 bit
%     spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    clear matlabbatch;
end

function [SubDir_fullpath,SubFolders] = MergeDir(parent_dir)
%% function to merge subdir full path and get the subfolder names list
    TempList=dir(parent_dir); %get the dir struct

    ISub=[TempList(:).isdir];

    SubFolders={TempList(ISub).name}';
    SubFolders(1:2)=[]; %remove . and ..

    SubDir_fullpath=strcat(parent_dir,'\',SubFolders);
end