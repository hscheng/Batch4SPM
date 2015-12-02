function Batch4SPM_1stLevel_addCon()
%% function to add contrast for 1st level analysis
% if data have multiple sessions,6 zeros add
% written by hongshengcheng.math@gmail.com
% 20141219,modify an bug in con_mark
% created date:20141201
% 
%     spm('defaults','fmri');
%     spm_jobman('initcfg');
    
    spm_mat_list = ScanImg2Cell('Slect SPM Mat Files','mat');
    consess = gen_consess;
    con_num = size(consess,2);
    
    for ii = 1:length(spm_mat_list)
        spm_mat_path = spm_mat_list{ii};
        sub_dir = fileparts(spm_mat_path);
        [~,sub_id] = fileparts(sub_dir);
        SPM = importdata(spm_mat_path);
        session_num = length(SPM.nscan);
%         cond_num = length(SPM.Sess(1,1).U);
%         con_vec_len = size(SPM.xCon(1,15).c,1);
                
        for jj =1:con_num
            con_mark = fieldnames(consess{1,jj});
            % con_mark is cell format,20141219
            switch con_mark{1}
                case 'tcon'
                    contrast_1run = consess{1,jj}.tcon.convec;
                    contrast_allrun = repmat(contrast_1run,1,session_num);
                    consess{1,jj}.tcon.convec = contrast_allrun;
                case 'fcon'
                    contrast_1run = consess{1,jj}.fcon.convec{1};
                    contrast_allrun = repmat(contrast_1run,1,session_num);
                    consess{1,jj}.fcon.convec{1} = contrast_allrun;
            end
            clear contrast_1run;
        end

    matlabbatch{1}.spm.stats.con.spmmat = {spm_mat_path};
    matlabbatch{1}.spm.stats.con.consess = consess;
    
    % save and run matlabbatch
    time_mark = datestr(clock,30);
    save(fullfile(sub_dir,['add_contrast_',sub_id,time_mark,'.mat']),'matlabbatch');
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    end
    
    msgbox('All work finished!','tip');
end

function consess = gen_consess()
%% function to define multiple f contrast and t contrast
    % define input dlg
    prompt = {'Enter the name of contrast:',...
        'Enter the contrast matrix:'};
    dlg_title = 'Input parameters';
    numlines = 1;
    para_default = {'tcon','1 0 0 0'};% e.g.kron([1 0 0 0]',ones(90,1))

    % add more regeressors
    Mark2Run = 1;
    Runtime = 1;
    while (Mark2Run == 1)
        choice_add = questdlg('Add more contrast?', 'Questdlg','T_Contrast','F_Contrast','No','No');
        switch choice_add
            case 'T_Contrast' 
                para = inputdlg(prompt,dlg_title,numlines,para_default);
                consess{1,Runtime}.tcon.name = para{1};
                consess{1,Runtime}.tcon.convec = str2num(para{2});
                Runtime = Runtime+1;
            case 'F_Contrast'
                para = inputdlg(prompt,dlg_title,numlines,para_default);
                consess{1,Runtime}.fcon.name = para{1};
                consess{1,Runtime}.fcon.convec = {eval(para{2})};
                Runtime = Runtime+1;
            case 'No'
                Mark2Run = 0;
        end
    end
end

function consess_t = gen_consess_tcon(con_txt_loc)
%% function to gen consess,which contains the contrast infomation
% Format of the contrast.txt file shoul be like
% WF-WZ 0 1 0 -1 0 0 0 0 0 0
% WF-JF -1 1 0 0 0 0 0 0 0 0
% WZ-JZ 0 0 -1 1 0 0 0 0 0 0
%
% 6 zeros in each contrast conditon[6 headmotion parameters]
% first string of each row is the name of contrast
% Note:Only support the T Contrast

    txt_info = importdata(con_txt_loc);

    tcon_num = size(txt_info.data,1);
    tcon_name = txt_info.textdata;
    % disp contrast infomation
    msg = sprintf('You have %d constracts.', tcon_num);
    disp(msg)
    convec_1run = txt_info.data;
    
    consess_t = cell(1,k);
    for kk = 1:tcon_num
        consess_t{1,kk}.tcon.name = tcon_name{kk,1};
        consess_t{1,kk}.tcon.convec = convec_1run(kk,:);
        consess_t{1,kk}.tcon.sessrep = 'none';
    end
end