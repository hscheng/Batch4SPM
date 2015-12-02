function Img_ArthCalcu(varargin)
%% function to calculate paired imgs
% written by hongshengcheng.math@gmail.com
% v3,20141116,support read para_path.xls
% v2:20141101 add inputdlg
% created date 20141029
% if numel(varargin)<1


% ParaNumMark = input('1-SinglePara;2-MultiParas(use para_path.xls) : ');
%if use para_path.xls,make sure the order is matched
ParaNumMark = questdlg('ParaNum','Mode','SinglePara','MultiParas','SinglePara');
switch ParaNumMark
    case 'SinglePara' % if only one parametre to perform paired ttest
        prompt = {'Prefix of OutputImgs','Name of OutputDir','Comput Way:+-*/'};
        dlg_title = 'Input';
        num_lines = 1;
        def = {'ArthCalcu','CalcuRes',''};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        % parameters can be changed
        prefix_outimg = answer{1};
        dir_output = answer{2};
        computway = answer{3}; % can be '+','-','*','/'

        mkdir(dir_output);

        % make sure 2 img list is well paired OR list2 choose single img
        ImgList1 = ScanImg2Cell('Select ImgList 1');
        ImgList2 = ScanImg2Cell('Select ImgList 2');
        ComputImgList(ImgList1,ImgList2,computway,prefix_outimg,dir_output)
    case 'MultiParas'
        [GrpNames,ParaNames,ParaCell] = ReadMultiPathInfo;
        computway = input('computway(+-*/) : ','s');
        prefix_outimg = input('prefix of output img file : ','s');
        for ii = 1:length(ParaNames)
            GrpScans = ParaCell(:,ii);
            OutputDir = fullfile(pwd,['CalcuRes_',ParaNames{ii}]);
            mkdir(OutputDir)
            ComputImgList(GrpScans{1},GrpScans{2},computway,prefix_outimg,OutputDir) 
        end
end
end

function ComputImgList(ImgList1,ImgList2,computway,OutPrefix,OutputDir)
    % Save data
    Data.computway = computway;
    Data.ImgList1 = ImgList1;
    Data.ImgList2 = ImgList2;
    save(fullfile(OutputDir,'Data.mat'),'Data')
    
    % Check the name be unique or not
    [pathstr1,ImgName1] = fileparts(ImgList1{1});
    [pathstr2,ImgName2] = fileparts(ImgList1{2});
    if strcmp(ImgName1,ImgName2) ~= 1
        mark_name = 1;% it means can use the img filename directly
    elseif strcmp(pathstr1,pathstr2) ~= 1
        mark_name =2;
    else error('please check your image file');
    end
    
    % ImgList2's length should equal to 1 OR length of ImgList1      
    if length(ImgList2) == 1
        ImgList2 = repmat(ImgList2,length(ImgList1),1);
    elseif  length(ImgList1) ~= length(ImgList2);
        error('Make sure two lists have same length OR list2 have only 1 img file');
    end
    
    for ii = 1:length(ImgList1)
        if mark_name == 1
            [~,ImgName,~] = fileparts(ImgList1{ii});
            OutputImgName = [OutPrefix,ImgName,'.nii'];
        else
            [pathstr,~,~] = fileparts(ImgList1{ii});
            [~,subdir,~] = fileparts(pathstr);
            OutputImgName = [OutPrefix,subdir,'.nii'];
        end

        ComputImg(ImgList1{ii},ImgList2{ii},OutputImgName,computway)%minus images
    end
        
        movefile([OutPrefix,'*.nii'],OutputDir);
        if strcmp(computway,'-')
            movefile('Ratio*.nii',OutputDir);%only use for the computway is minus
        end

    clc;
    disp('All Work Done!');
end

function ComputImg(InputImg1,InputImg2,OutputImgName,computway)
%% function to calculate 2 imgs
% Img1      fullpath of img file location
% OutputImg fullpath of output img file name
% computway can be '+','-','*','/'

V1 = spm_vol(InputImg1);
Y1 = spm_read_vols(V1);

V2 = spm_vol(InputImg2);
Y2 = spm_read_vols(V2);
switch computway
    case '+'
        Y_out = Y1 + Y2;
    case '-'
        Y_out = Y1 - Y2;
        
        %comput change ratio
        Y_ratio = Y_out./Y1;
        Y_ratio(Y1==0) = 0;
        % write img
        V_out = V1;
        V_out.fname = ['Ratio_',OutputImgName];
        spm_write_vol(V_out,Y_ratio);
        
    case '*'
        Y_out = Y1 .* Y2;
    case '/'
        Y_out = Y1 ./ Y2;
    otherwise
        error('parameter computway setting wrong !');
end
V_out = V1;
V_out.fname = OutputImgName;
Voutput = spm_write_vol(V_out,Y_out);
end

function img_list_cell = ScanImg2Cell(varargin)
% scan img files and reture the cell format list
% 
    if numel(varargin) == 1
        dlg_title = varargin{1};
    else
        dlg_title = 'Select Img Files';
    end
    img_list = spm_select(Inf,'image',dlg_title,'',pwd,'.*',1);
    img_num = size(img_list,1);
    img_list_cell = cell(img_num,1);
    for ii = 1:img_num
        % modified 20141012
        % fix a bug which will generate extra space
        str_temp = img_list(ii,:);
        str_ok = str_temp(~isspace(str_temp));
        img_list_cell{ii} = str_ok;
    end
end
