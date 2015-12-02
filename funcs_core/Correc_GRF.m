function Correc_GRF()
% this script is run based on rest toolbox


[mask_name,mask_path] = uigetfile('*.img;*.nii','Select Mask File');
mask_file = fullfile(mask_path,mask_name);

prompt = {'two tail','voxel level p','cluster level p','output dir name'};
dlg_title = 'Input';
num_lines = 1;
def = {'1','0.01','0.05','GRF_ThreshedImgRes'};
para = inputdlg(prompt,dlg_title,num_lines,def);
para{5} = mask_file;

StatImgList = ScanImg2Cell('Select Stat Imgs');
tic;
% if matlabpool_mode == 1
if matlabpool('size')>1
    parfor ii = 1:length(StatImgList)
        GRF_Thresh_Img(StatImgList{ii},para)
    end
else
    for ii = 1:length(StatImgList)
        GRF_Thresh_Img(StatImgList{ii},para)
    end
end
toc;

if matlabpool('size') > 0
    matlabpool close;
end
clear all;
end
function GRF_Thresh_Img(StatImg,para)
is_two_tag = str2num(para{1});%1-two tail;0-one tail
voxel_p = str2num(para{2});
cluster_p = str2num(para{3});
output_dir = para{4};
mask_file = para{5};
% is_two_tag = 1;%1-two tail;0-one tail
% voxel_p = 0.01;
% cluster_p = 0.05;
% output_dir = uigetdir(pwd,'Select output dir');
% output_dir = 'GRF_ThreshedImgRes';
mkdir(output_dir);

[Outdata,VoxelSize,Header] = rest_readfile(StatImg,1);

if isfield(Header,'descrip')
    headinfo = Header.descrip; 
    testDf2 = 0;
    if ~isempty(strfind(headinfo,'{T_['))% dong 100331 begin
        testFlag='T';
        Tstart=strfind(headinfo,'{T_[')+length('{T_[');
        Tend=strfind(headinfo,']}')-1;
        testDf = str2num(headinfo(Tstart:Tend));
    elseif ~isempty(strfind(headinfo,'{F_['))
        testFlag='F';
        Tstart=strfind(headinfo,'{F_[')+length('{F_[');
        Tend=strfind(headinfo,']}')-1;
        F_Df = str2num(headinfo(Tstart:Tend));
        testDf=F_Df(1,1);
        testDf2=F_Df(1,2);
    elseif ~isempty(strfind(headinfo,'{R_['))
        testFlag='R';
        Tstart=strfind(headinfo,'{R_[')+length('{R_[');
        Tend=strfind(headinfo,']}')-1;
        testDf = str2num(headinfo(Tstart:Tend));
    elseif ~isempty(strfind(headinfo,'{Z_['))
        testFlag='Z';
        Tstart=strfind(headinfo,'{Z_[')+length('{Z_[');
        Tend=strfind(headinfo,']}')-1;
        testDf = str2num(headinfo(Tstart:Tend));
    end
end
[pathstr,filename,ext] = fileparts(StatImg);
[pathstr2,subdir] = fileparts(pathstr);
output_name = fullfile(output_dir,[subdir,'_',filename]);
rest_GRF_Threshold(StatImg,...   %rest_GRF_Threshold(Data,...
    voxel_p,...
    is_two_tag,...
    cluster_p,...
    output_name,...
    mask_file,...
    testFlag,...
    testDf,...
    testDf2,...
    VoxelSize,...
    Header); %rest_GRF_Threshold by YAN Chao-Gan
end