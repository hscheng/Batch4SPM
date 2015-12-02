function GenMask_Template()
% written by H.S.Cheng hongshengcheng.math@gmail.com
% 20140828: output modify
% date: 20131206
%
% function gen_Mask_from_template
% A Template have multi indexs, and if you want some of the indexed ROI to
% make a Mask---This script will be help
%
% INPUT
% Template
%
% OUTPUT
% ROI_Mask_index
[filename,filepath] = uigetfile({'*.nii;*.img'},'Please select template');
Vtemp = spm_vol(fullfile(filepath,filename));
[Ytemp,XYZtemp] = spm_read_vols(Vtemp);
[dx,dy,dz] = size(Ytemp);

ROI_Index = input('Please enter roi index(e.g.[1:3 6:10]): ');
roi_num = length(ROI_Index);
Mask_mode = questdlg('What kind of mask?', 'Questdlg','No Merge','Merge','No Merge');
switch Mask_mode
    case 'No Merge'
        for ii = 1:roi_num
            Ymask = zeros(dx,dy,dz);
            index = find(Ytemp==ROI_Index(ii));
            Ymask(index) = 1;
            Vmask = Vtemp;
            Vmask.fname = ['Mask_ROI_',num2str(ROI_Index(ii)),'.nii'];
            spm_write_vol(Vmask,Ymask);
        end
    case 'Merge'
        for ii = 1:roi_num
            Ymask = zeros(dx,dy,dz);
            if ii == 1
                index = find(Ytemp==ROI_Index(ii));
            else
                index_temp = find(Ytemp==ROI_Index(ii));
                index = [index;index_temp];
            end
            Ymask(index) = 1;
            Vmask = Vtemp;
            Vmask.fname = 'Merged_Mask.nii';
            spm_write_vol(Vmask,Ymask);
        end
end
disp('All work Done!')
end