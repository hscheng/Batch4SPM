function Batch4SPM_BMS()
%% function to rum bayesian model selection 

% written by hongshengcheng.math@gmail.com
% created date:20141203

    clear,clc;
    h = figure;
    set(h,'color','w','numbertitle','off','menubar','none', 'resize','off',...
        'name','BMS','position',[400 400 300 150]);%[left, bottom, width, height]

    Button1 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.7 0.7 0.2],'fontsize',10,'fontweight','bold',...
        'String','generate voi',...
        'TooltipString','generate voi mat and img files,then extract timecourse');
    Button2 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.4 0.7 0.2],'fontsize',10,'fontweight','bold',...
        'String','gen ppi struct',...
        'TooltipString','batch to gen ppi struct');
    Button3 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.1 0.7 0.2],'fontsize',10,'fontweight','bold',...
        'String','ppi glm stat',...
        'TooltipString','perform glm');

    set(Button1,'Callback',@batch_gen_voi_callback);
    set(Button2,'Callback',@batch_gen_ppi_struct_callback);
    set(Button3,'Callback',@batch_ppi_glm_callback);
end
% http://cn.mathworks.com/help/matlab/creating_guis/write-callbacks-using-the-programmatic-workflow.html

function batch_gen_voi_callback(hObject,callbackdata)
    

end

function gen_dcm_model_space()

sub_dir_list = MergeDir();
sub_num = length(sub_dir_list);

para = inputdlg('DCM model names (each sub)','Enter',[10 30],{''});
str_temp = para{1};
dcm_mod_num = size(str_temp,1);
dcm_mod_list = cell(dcm_mod_num,1);
for kk = 1:dcm_mod_num
    str = str_temp(kk,:);
    str(isspace(str)) = [];
    dcm_mod_list{kk} = str;
end
    
for ii = 1:sub_num
    for jj = 1:dcm_mod_num
        subj(ii).sess(1).model(jj).fname = ...
            fullfile(sub_dir_list{ii},dcm_mod_list{jj});
    end
end
time_mark = datestr(clock,30);
save(['dcm_model_space_',time_mark],'subj');
end