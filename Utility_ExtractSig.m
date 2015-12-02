function Utility_ExtractSig()
    clear,clc;
    h = figure;
    set(h,'color','w','numbertitle','off','menubar','none', 'resize','off',...
        'name','ExtractSig','position',[450 450 250 200]);%[left, bottom, width, height]

    Button1 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.75 0.7 0.15],'fontsize',10,'fontweight','bold',...
        'String','rex',...
        'TooltipString','for multi-masks');
    Button2 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.55 0.7 0.15],'fontsize',10,'fontweight','bold',...
        'String','ExtractSig_Template',...
        'TooltipString','extract signal besed on template index');
    Button3 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.35 0.7 0.15],'fontsize',10,'fontweight','bold',...
        'String','ExtractSig_Voxels',...
        'TooltipString','extract signal in mask');
    Button4 = uicontrol('Parent',h,'Style','pushbutton','Units','normalized',...
        'Position',[0.15 0.15 0.7 0.15],'fontsize',10,'fontweight','bold',...
        'String','ExtractPerSigChange',...
        'TooltipString','extract percent signal change');

    set(Button1,'Callback',@Rex_Callback);
    set(Button2,'Callback',@ExtractSig_Template_Callback);
    set(Button3,'Callback',@ExtractSig_Voxels_Callback);
    set(Button4,'Callback',@ExtractPerSigChange_Callback);
end
% http://cn.mathworks.com/help/matlab/creating_guis/write-callbacks-using-the-programmatic-workflow.html
function Rex_Callback(hObject,callbackdata)
    rex;
end

function ExtractSig_Template_Callback(hObject,callbackdata)
    ExtractSig_Template;
end

function ExtractSig_Voxels_Callback(hObject,callbackdata)
    ExtractSig_Voxels;
end

function ExtractPerSigChange_Callback(hObject,callbackdata)
    ExtractPerSigChange;
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