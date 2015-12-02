function Batch4SPM_dcm2nii(varargin)
%% function to convert multiple *.dcm files to *.nii files
% Batch4SPM_dcm2nii(RawDcmDir,ImgNum2Del,subdir_mark)
% RawDcmDir     raw dcm file data dir
% ImgNum2Del	for rest or task data,default,0
% subdir_mark   'multiple runs'-have subfolders
%               'single run',without subfolders,default,'sub'
% 
% make sure dcm2nii.exe is in Matlab path
% make sure the default seting is to 4D nii file
% 
% written by hongshengcheng.math@gmail.com
% 20150411,fix subdir_mark
% 20150408,fix movefile rename bug
% created date:20141201


    switch numel(varargin)
        case 0
            RawDcmDir = uigetdir(pwd,'Select Data Dir');
            ImgNum2Del = input('Remove First Img Number : ');
            subdir_mark = questdlg('run number in one subject?','quest','multiple runs','single run','multiple runs');
        case 1
            RawDcmDir = varargin{1};
            ImgNum2Del = input('Remove First Img Number : ');
            subdir_mark = questdlg('run number in one subject?','quest','multiple runs','single run','multiple runs');
        case 2
            RawDcmDir = varargin{1};
            ImgNum2Del = varargin{2};
            subdir_mark = questdlg('run number in one subject?','quest','multiple runs','single run','multiple runs');
        case 3
            RawDcmDir = varargin{1};
            ImgNum2Del = varargin{2};
            subdir_mark = varargin{3};
    end
    
    % get the location of the dcm2nii.exe
    dcm2nii_loc = which('dcm2nii.exe');
    [pathstr, name, ~] = fileparts(dcm2nii_loc);
    dcm2nii_path = fullfile(pathstr,name);
    
    % make FunImg dir to store transfered imgs
    [WorkDir,~,~] = fileparts(RawDcmDir);
    FunImgDir = fullfile(WorkDir,'FunImg');
    mkdir(FunImgDir);
    
    SubDataPathList = MergeDir(RawDcmDir);

    SubNum = size(SubDataPathList,1);

    if matlabpool('size')>1
        parfor ii = 1:SubNum
            SubDcmDir= SubDataPathList{ii};
            if strcmp(subdir_mark,'multiple runs') == 1
                Sub_dcm2nii_sub(SubDcmDir,FunImgDir,ImgNum2Del,dcm2nii_path);
            else
                Sub_dcm2nii_nosub(SubDcmDir,FunImgDir,ImgNum2Del,dcm2nii_path);
            end
        end
    else
        for ii = 1:SubNum
            SubDcmDir= SubDataPathList{ii};
            if strcmp(subdir_mark,'multiple runs') == 1
                Sub_dcm2nii_sub(SubDcmDir,FunImgDir,ImgNum2Del,dcm2nii_path);
            else
                Sub_dcm2nii_nosub(SubDcmDir,FunImgDir,ImgNum2Del,dcm2nii_path);
            end
        end
    end
    cd(RawDcmDir)
end

function Sub_dcm2nii_sub(SubDcmDir,FunImgDir,ImgNum2Del,dcm2nii_exe_loc)
%% function to convert one subject dcm2nii,multiple sessions/runs
% Sub_dcm2nii_sub(SubDcmDir,FunImgDir,ImgNum2Del,dcm2nii_exe_loc)
% SubDcmDir         subjects raw dcm img file path
% FunImgDir         path to store 3d nii img file
% ImgNum2Del        deleter the first img number
% dcm2nii_exe_loc   location of dcm2nii.exe
% 
% each subject data have subfolders
% dcm files will be convert to 4D nii files first
% then transform 4D file to 3d files
% then delete first imgs if you want


    cd(SubDcmDir)

    [~,sub_id,~] = fileparts(SubDcmDir);
    
    run_path_list = MergeDir(SubDcmDir,'all');

    run_num = size(run_path_list,1);

    for jj = 1:run_num
        run_dir = run_path_list{jj};
        cd(run_dir);
        [~,run_id,~] = fileparts(run_dir);

        % mkdir to store Img files
        ImgDir4D = fullfile(FunImgDir,sub_id,run_id,'Img4D');
        if exist(ImgDir4D,'dir')==7
            delete(ImgDir4D);
        end

        mkdir(ImgDir4D);

        ImgDir3D = fullfile(FunImgDir,sub_id,run_id);

        % dcm 2 nii
        SourceFiles = fullfile(run_dir,'*.dcm');
        
        % get the operation system info
        if ispc
            dcm2nii_cmd = [dcm2nii_exe_loc,' -a y -d y -n y -v y -o ',ImgDir4D,' ',SourceFiles];
        else 
            dcm2nii_cmd = ['dcm2nii -a y -d y -n y -g n -v y -o ',ImgDir4D,' ',SourceFiles];
        end 
        
        % dcm2nii -a y -d y -n y -v y -o D:\Temp\Raw\sub01\ D:\Temp\Raw\sub01\sub01_run6\*.dcm

        % run the merged cmd
        dos(dcm2nii_cmd);

        % locate the 4d file
        niifilelist = dir(fullfile(ImgDir4D,'*.nii'));
%         nii_oldpath = fullfile(ImgDir4D,niifilelist.name);%20150408

        % rename the 4d nii files
        nii_newname = [sub_id,'_',run_id,'.nii'];
        nii_newpath = fullfile(ImgDir4D,nii_newname);
%         movefile(nii_oldpath,nii_newpath) %20150408
        cd(ImgDir4D);movefile(niifilelist.name,nii_newname);

        %convert 4d nii files to 3d nii files
        spm_file_split(nii_newpath,ImgDir3D);
        disp('Convert 4D Files to 3D Files ....')

        % delete first 5 imgs
        % you can change following settings as you need
        cd(ImgDir3D);
        if(ImgNum2Del > 0)
            img_list = dir('sub*run*.nii');
            delete(img_list(1:ImgNum2Del).name)
        end
%             delete('sub*run*00001.nii')
%             delete('sub*run*00002.nii')
%             delete('sub*run*00003.nii')
%             delete('sub*run*00004.nii')
%             delete('sub*run*00005.nii')
    end
end

function Sub_dcm2nii_nosub(SubDcmDir,FunImgDir,ImgNum2Del,dcm2nii_exe_loc)
%% function to convert one subject dcm2nii,singal session/run
% each subject data DO NOT have subfolders
%
% Sub_dcm2nii_sub(SubDcmDir,FunImgDir,ImgNum2Del,dcm2nii_exe_loc)
% SubDcmDir         subjects raw dcm img file path
% FunImgDir         path to store 3d nii img file
% ImgNum2Del        deleter the first img number
% dcm2nii_exe_loc   location of dcm2nii.exe
% 
    cd(SubDcmDir);
    [~,sub_id,~] = fileparts(SubDcmDir);

    % mkdir to store Img files
    ImgDir3D = fullfile(FunImgDir,sub_id);
    mkdir(ImgDir3D);

    % dcm 2 nii
    SourceFiles = fullfile(SubDcmDir,'*.dcm');

    % get the operation system info and merge dcm2nii cmd
    if ispc 
        dcm2nii_cmd = [dcm2nii_exe_loc,' -a y -d y -n y -v y -o ',ImgDir3D,' ',SourceFiles];
    else
        dcm2nii_cmd = ['dcm2nii -a y -d y -n y -g n -v y -o ',ImgDir3D,' ',SourceFiles];
    end 

    % run the merged cmd
    dos(dcm2nii_cmd);

    % locate the 4d file
    niifilelist = dir(fullfile(ImgDir3D,'*.nii'));
    nii_oldpath = fullfile(ImgDir3D,niifilelist.name);

    % rename the 4d nii files
    nii_newname = [sub_id,'.nii'];
    ImgDir4D = fullfile(ImgDir3D,'Img_4D');
    mkdir(ImgDir4D);
    nii_newpath = fullfile(ImgDir4D,nii_newname);
    movefile(nii_oldpath,nii_newpath);

    %convert 4d nii files to 3d nii files
    spm_file_split(nii_newpath,ImgDir3D);
    disp('Convert 4D Files to 3D Files ....')

    % delete first 5 imgs
    % you can change following settings as you need
    cd(ImgDir3D);
    if(ImgNum2Del > 0)
        img_list = dir('sub*.nii');
        delete(img_list(1:ImgNum2Del).name)
    end

%     % delete first imgs
%     if(ImgNum2Del > 0)
%         V = spm_vol(niifilepath);
%         V(1:ImgNum2Del) = [];
%         Y = spm_read_vols(V);
%         delete(niifilepath);
%         rest_Write4DNIfTI(Y,V(1),V(1).fname);
%     end

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