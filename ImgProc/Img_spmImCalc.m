function Img_spmImCalc()
%% function to run image calcu in batch mode
% note:only for paired imgs' calculation
% written by hongshengcheng.math@gmail.com
% created date : 20141103

    prompt = {'Prefix of OutputImgs','Name of OutputDir','Expression(e.g.i1+i2)'};
    dlg_title = 'Input';
    num_lines = 1;
    def = {'ArthCalcu','CalcuRes',''};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    % parameters can be changed
    prefix_outimg = answer{1};
    out_dir = answer{2};
    expression = answer{3}; 
    
    mkdir(out_dir);
    
    % make sure 2 img list is well paired OR list2 choose single img
    ImgList1 = ScanImg2Cell('Select ImgList 1');
    ImgList2 = ScanImg2Cell('Select ImgList 2');
    
    % Save data
    Data.computway = expression;
    Data.ImgList1 = ImgList1;
    Data.ImgList2 = ImgList2;
    save(fullfile(out_dir,'Data.mat'),'Data')
    
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
            out_name = [prefix_outimg,ImgName,'.nii'];
        else
            [pathstr,~,~] = fileparts(ImgList1{ii});
            [~,subdir,~] = fileparts(pathstr);
            out_name = [prefix_outimg,subdir,'.nii'];
        end

        matlabbatch{1,ii}.spm.util.imcalc.input = [ImgList1(ii);ImgList2(ii)];
        matlabbatch{1,ii}.spm.util.imcalc.output = out_name;
        matlabbatch{1,ii}.spm.util.imcalc.outdir = {out_dir};
        matlabbatch{1,ii}.spm.util.imcalc.expression = expression;
    end
    % save and run batch
    save(fullfile(out_dir,'ImgCalcuBatch.mat'),'matlabbatch');
%     spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    clear matlabbatch;
%     movefile([prefix_outimg,'*.nii'],out_dir);
    clc;
    disp('All Work Done!');
end
