function Img_dcm2nii()
%% funtion to convert dcm file to nii format
% only for data dirs without subfolders
% written by hongshengcheng.math@gmail.com
% 20141004: fix some bugs

    RawDcmDir = uigetdir(pwd,'Select Raw DICOM Img Dir');
    % make FunImg dir to store transfered imgs
    [WorkDir,~,~] = fileparts(RawDcmDir);
    FunImgPath = fullfile(WorkDir,'FunImg');
    mkdir(FunImgPath);
    
    
    delete_img_num = input('Remove First Img Number : ');
    
    % get the location of the dcm2nii.exe
    dcm2nii_loc = which('dcm2nii.exe');
    [pathstr, name, ~] = fileparts(dcm2nii_loc);
    dcm2nii_path = fullfile(pathstr,name);
    
    cd(RawDcmDir)
    [FullSubDataPath,SubID] = MergeDir(RawDcmDir);

    SubNum = size(FullSubDataPath,1);

% ********************** Matlabpool Settings ******************************
    % check whether matlabpool is runing
    matlabpool_state = matlabpool('size');

    % if matlabpool is active, then close it
    if matlabpool_state > 0
        matlabpool close;
    end

    % worker number set
    % 0 for the max; default is 1
    worker_num = input('How many workers to use?(Default:1) : ');

    if isempty(worker_num)|| (worker_num==1)
        for ii = 1:SubNum
            Sub_Dcm2Nii(FullSubDataPath{ii},FunImgPath,delete_img_num,dcm2nii_path)
        end
    else
        if worker_num > 1
            matlabpool('open',worker_num);
        else
            matlabpool local;
        end
        parfor ii = 1:SubNum
            Sub_Dcm2Nii(FullSubDataPath{ii},FunImgPath,delete_img_num,dcm2nii_path)
        end
    end
    
    % check whether matlabpool is runing
    matlabpool_state = matlabpool('size');

    % if matlabpool is active, then close it
    if matlabpool_state > 0
        matlabpool close;
    end
    cd(RawDcmDir)
    clc;
    disp('All work Done!');
end

function Sub_Dcm2Nii(SourceDir,OutputDir,ImgNum2Del,dcm2nii_exe_loc)
%% function to conver dicom img to nii format
    cd(SourceDir);
    [~,SubDirTemp,~] = fileparts(SourceDir);

    % mkdir to store Img files
    ImgDir4D = fullfile(OutputDir,SubDirTemp);
    mkdir(ImgDir4D);

    % dcm 2 nii
    SourceFiles = fullfile(SourceDir,'*.dcm');

    % get the operation system info and merge dcm2nii cmd
    if ispc 
        dcm2nii_cmd = [dcm2nii_exe_loc,' -a y -d y -n y -v y -o ',ImgDir4D,' ',SourceFiles];
    else
        dcm2nii_cmd = ['dcm2nii -a y -d y -n y -g n -v y -o ',ImgDir4D,' ',SourceFiles];
    end 

    % run the merged cmd
    dos(dcm2nii_cmd);

    % locate the 4d file
    niifilelist = dir(fullfile(ImgDir4D,'*.nii'));
    niifilepath = fullfile(ImgDir4D,niifilelist.name);

    % delete first imgs
    if(ImgNum2Del > 0)
        V = spm_vol(niifilepath);
        V(1:ImgNum2Del) = [];
        Y = spm_read_vols(V);
%         Y(:,:,:,1:ImgNum2Del) = [];
        delete(niifilepath);
        rest_Write4DNIfTI(Y,V(1),V(1).fname);

    end
end

function [SubDir_fullpath,SubFolders] = MergeDir(parent_dir)
%% function to merge subdir full path and get the subfolder names list
    TempList=dir(parent_dir); %get the dir struct

    ISub = [TempList(:).isdir];

    SubFolders = {TempList(ISub).name}';
    SubFolders(1:2) = []; %remove . and ..
    
    % get the operation system info
    if ispc % true if the os is windows
        SubDir_fullpath = strcat(parent_dir,'\',SubFolders);
    else 
        SubDir_fullpath = strcat(parent_dir,'/',SubFolders);
    end 
end

function [flag] = rest_Write4DNIfTI(Data,Head,imageOUT)
% Write 4D NIfTI file Based on SPM's nifti
% %------------------------------------------------------------------------
% Write data (Data) with a specified header (Head) into a image file with format 
% of Nifti 1.1. The data (Data) should be 3D matrix, the header (Head) should 
% be a structure the same as SPM5. If the filename (imageOUT) is with 
% extra name as '.img', then it will generate two files (header and
% data seperately), or else, '.nii', it will generate single file with
% header 
% and data together.
%
% Usage: [flag] = y_Write4DNIfTI(Data,Head,imageOUT)
%
% Input:
% 1. Data -  Data of 4D matrix to write
% 2. Head - a structure containing image volume information, the structure
%    is the same with a structure have read
%    The elements in the structure are:
%       Head.fname - the filename of the image. If the filename is not set, 
%                    just use the parameter.
%       Head.dt    - A 1x2 array.  First element is datatype (see spm_type).
%                 The second is 1 or 0 depending on the endian-ness.
%       Head.mat   - a 4x4 affine transformation matrix mapping from
%                 voxel coordinates to real world coordinates.
%       Head.pinfo - plane info for each plane of the volume.
%              Head.pinfo(1,:) - scale for each plane
%              Head.pinfo(2,:) - offset for each plane
%                 The true voxel intensities of the jth image are given
%                 by: val*Head.pinfo(1,j) + Head.pinfo(2,j)
%              Head.pinfo(3,:) - offset into image (in bytes).
%                 If the size of pinfo is 3x1, then the volume is assumed
%                 to be contiguous and each plane has the same scalefactor
%                 and offset.
%              The scale and intercept will be changed according to the
%              data to write
% 3. imageOUT - the path and filename of image file to output [path\*.img or *.nii]
% Output:
% 1. flag - a flag for all done, 1: successful, 0: fail
% ------------------------------------------------------------------------
% Written by YAN Chao-Gan 120301. Based on SPM's nifti.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com


dat = file_array;
dat.fname = imageOUT;
dat.dim   = size(Data);
if isfield(Head,'dt')
    dat.dtype  = Head.dt(1);
else % If data type is defined by the nifti command
    dat.dtype  = Head.dat.dtype;
end

dat.offset  = ceil(348/8)*8;

NIfTIObject = nifti;
NIfTIObject.dat=dat;
NIfTIObject.mat=Head.mat;
NIfTIObject.mat0 = Head.mat;
NIfTIObject.descrip = Head.descrip;

create(NIfTIObject);
dat(:,:,:,:)=Data;
end
