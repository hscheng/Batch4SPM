function Img_Reslice()
    % hld - interpolation method. 0: Nearest Neighbour. 1: Trilinear.
%     NewVoxSize = [3 3 3];
%     hld = 1;


    prompt = {'New Voxel Size','interpolation method. 0: Nearest Neighbour. 1: Trilinear'};
    dlg_title = 'Input';
    num_lines = 1;
    def = {num2str([3 3 3]),'0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    % parameters can be changed
    NewVoxSize = str2num(answer{1});
    hld = str2num(answer{2});
    
    ImgList = ScanFileList('Selece Imgs to Reslice');

    
    % Check the name be unique or not
    if length(ImgList) == 1
        mark_name = 1;
    else
        [pathstr1,ImgName1] = fileparts(ImgList{1});
        [pathstr2,ImgName2] = fileparts(ImgList{2});
        if strcmp(ImgName1,ImgName2) ~= 1
            mark_name = 1;% it means can use the img filename directly
        elseif strcmp(pathstr1,pathstr2) ~= 1
            mark_name = 2;
        else error('please check your image file');
        end
    end
    
    for ii = 1:length(ImgList)
        if mark_name == 1
            [~,ImgName,~] = fileparts(ImgList{ii});
            output_filename = ['Resliced_',ImgName,'.nii'];
        else
            pathstr = fileparts(ImgList{ii});
            [~,subdir,~] = fileparts(pathstr);
            output_filename = ['Resliced_',subdir,'.nii'];
        end
        
        y_Reslice(ImgList{ii},output_filename,NewVoxSize,hld,'ImageItself');
    end
    
    mkdir('ReslicedRes');
    movefile('Resliced*.nii','ReslicedRes');

    clear;
    clc;
    msgbox('All Work Done!',':)');
end
