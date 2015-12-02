% written by Hongsheng Cheng |hscheng.ok@gmail.com
% modified date: 20140624
% created date£º20140401

function ExstractFullDirPath()

wkdir = pwd;

DirMode = input('Dir Mode;1 all sub folders; 2 choosed sub folders : ');
% if you choose the mode 2,you need prepare a txt sub dir list
switch DirMode
    case 1
        DirAll = gen_subdir(wkdir);
        save DirAll DirAll
        
    case 2
        [dirname,dirpath] = uigetfile('*.txt;*.mat','Select dir list file');
        DirAll = importdata(fullfile(dirpath,dirname));
        
end


MergeMode = input('Merge Mode;1 full file path; 2 full dir path : ');
% when choose merge full file path ,the dir shoul be only parent dir

switch MergeMode
    case 1
        filetype = input('file type [e.g. *.img] : ','s');
        FullFilePath = merge_fullfilepath(DirAll,filetype);
        cd (wkdir)
        save FullFilePath FullFilePath
        
    case 2
        FullDirPath=merge_fulldir(DirAll);
        cd (wkdir)
        save FullDirPath FullDirPath
        
end
end
function FullFilePath=merge_fullfilepath(DirList,filetype)

    DirNum=size(DirList,1);
    for ii = 1:DirNum

        TempDir=DirList{ii};
        cd(TempDir)

%         TempList=dir('*.img'); %get the dir struct
        TempList=dir(filetype); %get the dir struct

        ISub=[TempList(:).isdir];

        IsFile=logical(1-ISub);

        NameFiles={TempList(IsFile).name}';

    %     NameFolds={TempDir(ISub).name}';
    %     NameFolds(1:2)=[]; %remove . and ..
    %     FileNum=sum(IsFile);

        FullFilePathTemp=strcat(DirList{ii},NameFiles);

        if ii == 1
            FullFilePath = FullFilePathTemp;
        else
            FullFilePath=[FullFilePath;FullFilePathTemp];
        end
    end

end

function FullDir=merge_fulldir(DirList)
    DirNum=size(DirList,1);
    for ii = 1:DirNum

        TempDir=DirList{ii};
        cd(TempDir)

        TempList=dir(TempDir); %get the dir struct

        ISub=[TempList(:).isdir];

        NameFolds={TempList(ISub).name}';
        NameFolds(1:2)=[]; %remove . and ..

        FullDirTemp=strcat(DirList{ii},'\',NameFolds);

        if ii == 1
            FullDir = FullDirTemp;
        else
            FullDir=[FullDir;FullDirTemp];
        end
    end

end

function FullSubDir=gen_subdir(wkdir)
Dir=uigetdir(wkdir,'Select your dir');

TempList=dir(Dir); %get the dir struct

ISub=[TempList(:).isdir];

NameFolds={TempList(ISub).name}';
NameFolds(1:2)=[]; %remove . and ..
FullSubDir=strcat(Dir,'\',NameFolds);
end

function DirecSubFile()
TempDir = uigetdir(pwd,'Select Data Dir');

filetype ='*.nii';
file_mark = fullfile(TempDir,filetype);
TempList=dir(file_mark); %get the dir struct

ISub=[TempList(:).isdir];

IsFile=logical(1-ISub);

NameFiles={TempList(IsFile).name}';

FullFilePath = strcat(TempDir,NameFiles);
save FullFilePath FullFilePath
end