function Utility_GenMask()
    clear,clc;
    h = figure;
    set(h,'color','w','numbertitle','off','menubar','none', 'resize','off',...
        'name','Create Mask','position',[400 400 300 200]);%[left, bottom, width, height]

    Button1 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.15 0.7 0.15],'fontsize',10,'fontweight','bold',...
        'String','BatchThreshImg2Mask',...
        'TooltipString','set threshold and gen masks');
    Button2 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.35 0.7 0.15],'fontsize',10,'fontweight','bold',...
        'String','MultiImgOverLapMask',...
        'TooltipString','multiple imgs overlap base on a threshold');
    Button3 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.55 0.7 0.15],'fontsize',10,'fontweight','bold',...
        'String','Batch2CreateBallMask',...
        'TooltipString','create ball mask in MNI space');
    Button4 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.75 0.7 0.15],'fontsize',10,'fontweight','bold',...
        'String','GenMaskFromTemplate',...
        'TooltipString','gen mask based on template idx');

    set(Button1,'Callback',@GenMask_Threshold_Callback);
    set(Button2,'Callback',@GenMask_MultiImgOverLap_Callback);
    set(Button3,'Callback',@GenMask_Ball_Callback);
    set(Button4,'Callback',@GenMask_Template_Callback);
end
% http://cn.mathworks.com/help/matlab/creating_guis/write-callbacks-using-the-programmatic-workflow.html
function GenMask_Threshold_Callback(hObject,callbackdata)
    GenMask_Threshold;
end

function GenMask_MultiImgOverLap_Callback(hObject,callbackdata)
    GenMask_MultiImgOverLap;
end

function GenMask_Ball_Callback(hObject,callbackdata)
    GenMask_Ball;
end

function GenMask_Template_Callback(hObject,callbackdata)
    GenMask_Template;
end

% function Batch_CreateMask()
%    myvar = 5;
%    figure
%    uicontrol('Style','pushbutton',...
%              'Callback',{@pushbutton_callback,myvar});
% end
% function pushbutton_callback(hObject,callbackdata,x)
%    display(x);
% end