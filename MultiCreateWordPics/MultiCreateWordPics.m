% written by hscheng [hongshengcheng.math@gmail.com]
%
% 20150405,modify a bug when make results dir
% modified date 20140506
% create date 20140430
% function to create multiple words for an experiment
% 
% depend toolbox: export_fig
% http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig
% 
% before run this script,you should prepare a excel file
% which 1st column is the Words
% and 2nd column is the img file names
%
% output figure format is bmp image files
% (which will be saved in the Results folder in the working dir)
% and you can change the format as you need

clear;
clc;


fprintf('Welcome to our lab:  <a href="http://www.qiujlab.com/">Qiu Jiang Lab</a>\nhttp://www.qiujlab.com/\n');
disp('****************************************************************')

% make a folder to store results figures
if exist(fullfile(pwd,'Results'),'dir')~=7
    mkdir Results
end

% figure parameters input
prompt = {'Font size[30+]','Font Color','Background Color','FontName'};
dlg_title = 'Input the Parameter';
num_lines = [1 50];
def = {'50','white','black','ËÎÌו'};
fig_paras = inputdlg(prompt,dlg_title,num_lines,def);

% get the parameters
FontSize=str2num(fig_paras{1}); %#ok<ST2NM>
FontColor=fig_paras{2};
BackgroundColor=fig_paras{3};
FontName=fig_paras{4};



[Word_filename,Word_filepath]=uigetfile('*.xls;*.xlsx','Please choose your word list file');
Word_fileloc = fullfile(Word_filepath,Word_filename);
%     [name,path] = uigetfile('*.txt','Select Onset.txt file');
%     onset_txt_loc = fullfile(path,name);
%     M = importdata(onset_txt_loc);
%     onset = M.data;
%     cond_name = M.textdata;

[num,txt,WordList]=xlsread(Word_fileloc,1);
% WordList = importdata(Word_fileloc);

WordNum=size(WordList,1);

if size(WordList,2)==1
    WordNames=WordList;
    FileNames=WordList;
else
    WordNames=WordList(:,1);
    FileNames=WordList(:,2);
end


for ii=1:WordNum
    
    % http://stackoverflow.com/questions/4638077/display-text-on-screen-matlab
    
    % To keep the background as in your figure, use the command set.
    % set(gcf, 'InvertHardCopy', 'off');
    
    word_len = length(WordNames{ii});
    figure('Color',BackgroundColor, 'Menu','none','Position',[0,0,FontSize*(word_len+2),FontSize*2]);

%     figure('Color',BackgroundColor, 'Menu','none','Position',[0,0,FontSize*4,FontSize*2]);

    text(0.5-1/32, 0.5-1/16,WordNames{ii}, 'FontSize',FontSize, 'FontName',FontName,...
        'Color',FontColor,'FontWeight','bold', ...
        'HorizontalAlignment','Center', 'VerticalAlignment','Middle')

    axis off
% ************************************************************************%
%     figure('Color','k', 'Menu','none','Position',[0,0,200,100])
% 
%     text(0.5, 0.5, WordNames{ii}, 'FontSize',50, 'FontName','ËÎÌו',...
%         'Color','white','FontWeight','bold', ...
%         'HorizontalAlignment','Center', 'VerticalAlignment','Middle')
% 
%     axis off
% ************************************************************************%

    % https://sites.google.com/site/oliverwoodford/software/export_fig
%     export_fig(fullfile(pwd,FileNames{ii}),'-bmp','-nocrop') %20141004
    f=getframe(gcf);
    imwrite(f.cdata,fullfile(pwd,'Results',[FileNames{ii},'.bmp']));
    
    close all %close all the figures

end
% close all %close all the figures
% movefile('*.bmp','Results');%20141004
% clear;
clc;
clear all;
disp('All work done,please check your files')
msgbox('All work done,please check your files','Congratulations!')