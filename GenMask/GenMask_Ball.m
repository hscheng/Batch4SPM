function GenMask_Ball()
% Modified by hscheng [hongshengcheng.math@gmail.com]
% v2,20141107,make the script more powerful
% Modified date 20140611
% Created date 20140424
% function to create multiple Ball Mask and merge a template
% 3*3*3 mm
% before run this script,please make sure the REST toolbox in your path
%
% INPUT
% MNI Corrdnate list(txt file),which format should be
% 3 4 5
% 6 8 9
% ...
% 
% OUTPUT
% Merged_Template.nii

    clear;clc;
%     [MNI_name,MNI_path] = uigetfile('*.txt','Select your MNI Coordinate txt file');
%     MNI_List=importdata(fullfile(MNI_path,MNI_name));

    MNI_Coord_List = ReadMultiCoords();
    
    % Check the Input
    if size(MNI_Coord_List,2)~=3
        error('MNI Coordinateformat is wrong ,Please check');
    end
    
    ROI_Num = size(MNI_Coord_List,1);

    ROI_Radious = input('Please Enter the ROI Radious(mm) e.g. 6 : ');
       
    voxel_size = input('Enter the Voxel Size of Mask (e.g. 3): ');
    template_dir = which('Template_1mm.nii');
    
    if voxel_size == 1
        img_dir = template_dir;
    else
        template_path = fileparts(template_dir);
        img_dir = fullfile(template_path,['Template_',num2str(voxel_size),'mm.nii']);
        if ~exist(img_dir)
            V_out = y_Reslice(template_dir,img_dir,repmat(voxel_size,1,3),0,'ImageItself');
        end
    end
    Vm = spm_vol(img_dir);
    Ym = spm_read_vols(Vm);
    
    % merge the ROIs or NOT
    merge_choose = questdlg('Merge ROIs into Template?','Choose','NO','Yes','NO');
    switch merge_choose
        case 'NO'
            merge_mark = 0;
        case 'Yes'
            merge_mark = 1;
    end
    % [IMG Vm] = rest_ReadNiftiImage(img_dir);  
    for roi_idx = 1:ROI_Num
        MNIcoord = MNI_Coord_List(roi_idx,:);
        ROI_Temp = BallSelect(MNIcoord,ROI_Radious,Ym,Vm);
    %     ImgName=['Ball_Mask_',num2str(roi_idx),'.nii'];
    %     BallSelect(MNIcoord,ROI_Radious,img_dir,ImgName)
        if merge_mark == 1
            if roi_idx == 1
                Y_Merge = ROI_Temp;
            else
                Y_Merge = Y_Merge+ROI_Temp*roi_idx;
            end
        else       
            Vm.fname = ['ROI_',num2str(MNIcoord(1)),'_',num2str(MNIcoord(2)),'_',num2str(MNIcoord(3)),'_',num2str(ROI_Radious),'mm.nii'];
            spm_write_vol(Vm,ROI_Temp); 
        end    
    end
    
    if merge_mark == 1
        Vm.fname=['Merged_Template_',num2str(ROI_Num),'rois_',num2str(ROI_Radious),'mm.nii'];
        % Vm.dt=[16 0]; %[2 0] means binary image
        % Vm.mat(1,:)=[-3,0,0,93];
        % Vm.pinfo=[1;0;352;];
        spm_write_vol(Vm,Y_Merge);       
    end
    clear;
    disp('All work done, please check result files')
end

function MAT_out = BallSelect(MNI_coord,radius,IMG,HDR)
%% function BallSelect(MNI_coord,radius)
% MNI_coord could be an 1*3 or a 3*1 matrix;
% The unit of radius is mm;
% Note that the image loading in should be a normalized image;
% and the voxelsize should all equal!

% img_dir = spm_select(1,'image','MNI Ball mask image');
% [IMG HDR] = rest_ReadNiftiImage(img_dir);   
% Actually, We think that the image comes from spm8, so its matrix must be
% like this: [-a 0 0 -ori(1); 0 a 0 ori(2); 0 0 a ori(3)]; So this
% function is only suit for spm8's result!
DIMs = size(IMG);
MAT_out = zeros(DIMs);
clear IMG;
matrix_hdr = HDR.mat;
inv_mat = matrix_hdr([1:3],[1:3]);
ori_mat = matrix_hdr([1:3],4);
if inv_mat(1,1)<0
    inv_mat(1,:) = inv_mat(1,:).*-1;
    ori_mat(1) = ori_mat(1).*-1;
end
v_size = abs(inv_mat(1,1));
max_v = round(radius/v_size);
MNI_size = size(MNI_coord);
if MNI_size(1)<MNI_size(2);
    mni_coord = MNI_coord';
else
    mni_coord = MNI_coord;
end
phy_coord = round(inv(inv_mat)*(mni_coord-ori_mat));
if phy_coord(1)<0||phy_coord(1)>DIMs(1)||phy_coord(2)<0||phy_coord(2)>DIMs(2)...
        ||phy_coord(3)<0||phy_coord(3)>DIMs(3)
    error('Not a correct input MNI coordinate!');
end

if (phy_coord(1)-max_v>=1)&&(phy_coord(1)+max_v<=DIMs(1)) %% Which means that the ball is complete in the brain in x dims
    st_x = phy_coord(1)-max_v;
    ed_x = phy_coord(1)+max_v;
elseif (phy_coord(1)-max_v<1)&&(phy_coord(1)+max_v<=DIMs(1)) %% Which means that the ball is in the left side of the x dims, in this condition, the ball would lose the left side in x dims!
    st_x = 1;
    ed_x = phy_coord(1)+max_v;
elseif (phy_coord(1)-max_v>=1)&&(phy_coord(1)+max_v>DIMs(1)) %% Which means that the ball is in the right side of the x dims, in this condition, the ball would lose the right side in x dims!
    st_x = phy_coord(1)-max_v;
    ed_x = DIMs(1);
else      %% Which means that the ball is so large that contained the whole brain in x dims!
    st_x = 1;
    ed_x = DIMs(1);
end

if (phy_coord(2)-max_v>=1)&&(phy_coord(1)+max_v<=DIMs(2)) %% Which means that the ball is complete in the brain in y dims
    st_y = phy_coord(2)-max_v;
    ed_y = phy_coord(2)+max_v;
elseif (phy_coord(2)-max_v<1)&&(phy_coord(2)+max_v<=DIMs(2)) %% Which means that the ball is in the left side of the y dims, in this condition, the ball would lose the left side in y dims!
    st_y = 1;
    ed_y = phy_coord(2)+max_v;
elseif (phy_coord(2)-max_v>=1)&&(phy_coord(2)+max_v>DIMs(2)) %% Which means that the ball is in the right side of the y dims, in this condition, the ball would lose the right side in y dims!
    st_y = phy_coord(2)-max_v;
    ed_y = DIMs(2);
else      %% Which means that the ball is so large that contained the whole brain in y dims!
    st_y = 1;
    ed_y = DIMs(2);
end

if (phy_coord(3)-max_v>=1)&&(phy_coord(3)+max_v<=DIMs(3)) %% Which means that the ball is complete in the brain in z dims
    st_z = phy_coord(3)-max_v;
    ed_z = phy_coord(3)+max_v;
elseif (phy_coord(3)-max_v<1)&&(phy_coord(3)+max_v<=DIMs(3)) %% Which means that the ball is in the left side of the z dims, in this condition, the ball would lose the left side in z dims!
    st_z = 1;
    ed_z = phy_coord(2)+max_v;
elseif (phy_coord(3)-max_v>=1)&&(phy_coord(3)+max_v>DIMs(3)) %% Which means that the ball is in the right side of the z dims, in this condition, the ball would lose the right side in z dims!
    st_z = phy_coord(2)-max_v;
    ed_z = DIMs(3);
else      %% Which means that the ball is so large that contained the whole brain in z dims!
    st_z = 1;
    ed_z = DIMs(3);
end
    

% Next part is to make a ball
for i = st_x:ed_x
    for j = st_y:ed_y
        for k = st_z:ed_z
            if norm(([i, j, k]' - phy_coord).*v_size)<=radius
                MAT_out(i,j,k) = 1;
            end
        end
    end
end

% HDR.mat(1,:) = HDR.mat(1,:).*-1;
if HDR.mat(1,1)<0         %% Here we filp x dims for its matrix is negative!
    MAT_out = flipdim(MAT_out,1);
end
if HDR.mat(2,2)<0         %% Here we filp y dims for its matrix is negative!
    MAT_out = flipdim(MAT_out,2);
end
if HDR.mat(3,3)<0         %% Here we filp z dims for its matrix is negative!
    MAT_out = flipdim(MAT_out,3);
end
% rest_WriteNiftiImage(MAT_out,HDR,'Ball_mask.nii');
% rest_WriteNiftiImage(MAT_out,HDR,outputname)
end

function CoordList = ReadMultiCoords(varargin)
%% function to copy data direct from xls file
% make sure the end of the data should be cleaned

    if isempty(varargin)
        dlg_title = 'Enter Multiple Coordinates';
    else
        dlg_title = varargin{1};
    end
    Paras = inputdlg('Enter Coordinates(e.g. 6,7,8)',dlg_title,[10 30]);
    RawInput = Paras{1,1};
    CoordList = str2num(RawInput);
end