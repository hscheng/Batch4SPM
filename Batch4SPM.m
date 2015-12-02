function varargout = Batch4SPM(varargin)
% AHONGTOOLKIT MATLAB code for AHONGTOOLKIT.fig
%      AHONGTOOLKIT, by itself, creates a new AHONGTOOLKIT or raises the existing
%      singleton*.
%
%      H = AHONGTOOLKIT returns the handle to a new AHONGTOOLKIT or the handle to
%      the existing singleton*.
%
%      AHONGTOOLKIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AHONGTOOLKIT.M with the given input arguments.
%
%      AHONGTOOLKIT('Property','Value',...) creates a new AHONGTOOLKIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AHONGTOOLKIT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AHONGTOOLKIT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AHONGTOOLKIT

% Last Modified by GUIDE v2.5 02-Dec-2014 11:02:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Batch4SPM_OpeningFcn, ...
                   'gui_OutputFcn',  @Batch4SPM_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before AHONGTOOLKIT is made visible.
function Batch4SPM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AHONGTOOLKIT (see VARARGIN)

% Choose default command line output for AHONGTOOLKIT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
clc;
fprintf('Welcome: %s, \nVersion: %s\n', getenv('USERNAME'),'v0.3d_20150611');
fprintf('Welcome to our lab:  <a href="http://www.qiujlab.com/">Qiu Jiang Lab</a>\nhttp://www.qiujlab.com/\n');
disp('Any question, Please contact me at hongshengcheng.math@gmail.com');
disp('****************************************************************');
% disp('Initializing spm........');
% spm('defaults','fmri');
% spm_jobman('initcfg');
% disp('ready for work');

% UIWAIT makes AHONGTOOLKIT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Batch4SPM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% Batch for SPM:Preproc,1stLevel,2ndLevel
% --- Executes on button press in Batch4SPM_Preproc.
function Batch4SPM_Preproc_Callback(hObject, eventdata, handles)
% hObject    handle to Batch4SPM_Preproc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
design_type = questdlg('Design Type','Quest','RapidEvent','Block','RapidEvent');
if strcmp(design_type,'RapidEvent')
    disp('Rapid Event Related Design')
    Batch4SPM_Preproc;
else
    disp('Block Design')
    Batch4SPM_Preproc_Block;
end
clc;disp('All Work Done');

% --- Executes on button press in Batch4SPM_1stLevel.
function Batch4SPM_1stLevel_Callback(hObject, eventdata, handles)
% hObject    handle to Batch4SPM_1stLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Batch4SPM_1stLevel;
clc;msgbox('First Level Work Done!',':)')

% --- Executes on button press in Batch4SPM_2ndLevel.
function Batch4SPM_2ndLevel_Callback(hObject, eventdata, handles)
% hObject    handle to Batch4SPM_2ndLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Batch4SPM_2ndLevel;
clc;disp('All Work Done');

% --- Executes on button press in Batch4SPM_DCM.
function Batch4SPM_DCM_Callback(hObject, eventdata, handles)
% hObject    handle to Batch4SPM_DCM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Batch4SPM_DCM;
clc;

% --- Executes on button press in Batch4SPM_PPI.
function Batch4SPM_PPI_Callback(hObject, eventdata, handles)
% hObject    handle to Batch4SPM_PPI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Batch4SPM_PPI;
clc;

%% Exstract Signal
% % --- Executes on button press in rex.
% function rex_Callback(hObject, eventdata, handles)
% % hObject    handle to rex (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% clear;clc;
% rex;
% clc;disp('All Work Done');
% 
% % --- Executes on button press in ExtractPerSigChange.
% function ExtractPerSigChange_Callback(hObject, eventdata, handles)
% % hObject    handle to ExtractPerSigChange (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% clear;clc;
% ExtractPerSigChange;
% clc;disp('All Work Done');
% 
% % --- Executes on button press in ExtractSig_Voxels.
% function ExtractSig_Voxels_Callback(hObject, eventdata, handles)
% % hObject    handle to ExtractSig_Voxels (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% clear;clc;
% ExtractSig_Voxels;
% clc;disp('All Work Done');
% 
% % --- Executes on button press in ExtractSig_Template.
% function ExtractSig_Template_Callback(hObject, eventdata, handles)
% % hObject    handle to ExtractSig_Template (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% clear;clc;
% ExtractSig_Template;
% clc;disp('All Work Done');

%% Brain Areas Report
% --- Executes on button press in BA_Rep_Thresh.
function BA_Rep_Thresh_Callback(hObject, eventdata, handles)
% hObject    handle to BA_Rep_Thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
BA_Rep_Thresh;
clc;disp('All Work Done');

% --- Executes on button press in BA_Rep_TemplateCluster.
function BA_Rep_TemplateCluster_Callback(hObject, eventdata, handles)
% hObject    handle to BA_Rep_TemplateCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
BA_Rep_TemplateCluster;
clc;disp('All Work Done');

%% Img Viewer
% --- Executes on button press in Img_CheckReg.
function Img_CheckReg_Callback(hObject, eventdata, handles)
% hObject    handle to Img_CheckReg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Img_CheckReg;
clc;disp('All Work Done');

% --- Executes on button press in Img_CheckSpmT.
function Img_CheckSpmT_Callback(hObject, eventdata, handles)
% hObject    handle to Img_CheckSpmT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Img_CheckSpmT;
clc;disp('All Work Done');

% --- Executes on button press in DPABI_VIEW.
function DPABI_VIEW_Callback(hObject, eventdata, handles)
% hObject    handle to DPABI_VIEW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
DPABI_VIEW;

% --- Executes on button press in rest_sliceviewer.
function rest_sliceviewer_Callback(hObject, eventdata, handles)
% hObject    handle to rest_sliceviewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
rest_sliceviewer;

%% Utilities
% --- Executes on button press in DPABI_Utilities.
function DPABI_Utilities_Callback(hObject, eventdata, handles)
% hObject    handle to DPABI_Utilities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
DPABI_Utilities;

% --- Executes on button press in Rest_Utilities.
function Rest_Utilities_Callback(hObject, eventdata, handles)
% hObject    handle to Rest_Utilities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
rest_Utilities_gui;

% --- Executes on button press in Batch_genMask.
function Batch_genMask_Callback(hObject, eventdata, handles)
% hObject    handle to Batch_genMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Utility_GenMask;
clc;disp('All Work Done');

% --- Executes on button press in Utility_ExtractSig.
function Utility_ExtractSig_Callback(hObject, eventdata, handles)
% hObject    handle to Utility_ExtractSig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Utility_ExtractSig;
clc;disp('All Work Done');

%% Batch exe
% --- Executes on button press in BulkRenameUtility.
function BulkRenameUtility_Callback(hObject, eventdata, handles)
% hObject    handle to BulkRenameUtility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
exepath = which('BulkRenameUtility.exe');
winopen(exepath);
clear;clc;

% --- Executes on button press in BatchDelFile.
function BatchDelFile_Callback(hObject, eventdata, handles)
% hObject    handle to BatchDelFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
exepath = which('批量delete.exe');
winopen(exepath);
clear;clc;

% --- Executes on button press in BatchCpFile_1Folder.
function BatchCpFile_1Folder_Callback(hObject, eventdata, handles)
% hObject    handle to BatchCpFile_1Folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
exepath = which('全文件复制.exe');
winopen(exepath);
clear;clc;

% --- Executes on button press in BatchCpFilewithDir.
function BatchCpFilewithDir_Callback(hObject, eventdata, handles)
% hObject    handle to BatchCpFilewithDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exepath = which('保文件夹复制.exe');
winopen(exepath);
clear;clc;

% --- Executes on button press in BatchGenFilePath.
function BatchGenFilePath_Callback(hObject, eventdata, handles)
% hObject    handle to BatchGenFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
exepath = which('批量导出文件路径.exe');
winopen(exepath);
clear;clc;

%% Img process
% --- Executes on button press in Img_ArthCalcu.
function Img_ArthCalcu_Callback(hObject, eventdata, handles)
% hObject    handle to Img_ArthCalcu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Img_ArthCalcu;
clc;disp('All Work Done');

% --- Executes on button press in Img_spmImCalc.
function Img_spmImCalc_Callback(hObject, eventdata, handles)
% hObject    handle to Img_spmImCalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Img_spmImCalc;
clc;disp('All Work Done');

% --- Executes on button press in Img_dcm2nii.
function Img_dcm2nii_Callback(hObject, eventdata, handles)
% hObject    handle to Img_dcm2nii (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Img_dcm2nii;
clc;disp('All Work Done');

% --- Executes on button press in Img_Reslice.
function Img_Reslice_Callback(hObject, eventdata, handles)
% hObject    handle to Img_Reslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Img_Reslice;
clc;disp('All Work Done');

% --- Executes on button press in Img_Convert.
function Img_Convert_Callback(hObject, eventdata, handles)
% hObject    handle to Img_Convert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;

choice = questdlg('Choose Mode','Img Convert',...
    'Convert dcm to nii','Convert 3D to 4D','Convert 4D to 3D',...
    'Convert dcm to nii') ;
switch choice
    case 'Convert dcm to nii'
        Batch4SPM_dcm2nii; % Img_dcm2nii;
    case 'Convert 3D to 4D'
        Img_3Dto4D;
    case 'Convert 4D to 3D'
        Img_4Dto3D;
    otherwise
        error('wrong input')
end

clc;disp('All Work Done');

% --- Executes on button press in Correc_GRF.
function Correc_GRF_Callback(hObject, eventdata, handles)
% hObject    handle to Correc_GRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
Correc_GRF;
clc;disp('All Work Done');

% --- Executes on button press in rsFC_roi2roi.
function rsFC_roi2roi_Callback(hObject, eventdata, handles)
% hObject    handle to rsFC_roi2roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear;clc;
rsFC_roi2roi;
clc;disp('All Work Done');
