function GenMask_Threshold()
%% function to create masked by threshold

    ImgList = ScanImg2Cell('Select Imgs to Make Mask');
    threshold = input('Input Threshold : ');
    ResDir = 'ThresholdedMask';
    mkdir(ResDir);
    for ii = 1:length(ImgList)
        [PathStr,ImgName] = fileparts(ImgList{ii});
        [~,SubDir] = fileparts(PathStr);
        outputname = ['tMask_',SubDir,'_',ImgName,'.nii'];
        GenMaskThresh(ImgList{ii},threshold,outputname);
    end
    movefile('tMask_*.nii',ResDir);
    
end

function GenMaskThresh(ImgFile,Threshold,OutputName)
%% function to merge multiple binary mask imgs into templates
% written by hongshengcheng.math@gmail.com
% INPUT
% ImgFile Threshold
% OUTPUT
% Mask.nii

    V =  spm_vol(ImgFile);
    Y = spm_read_vols(V);
    Mask = zeros(size(Y));
    if Threshold == 0
        Mask(Y~=0) = 1;
    elseif Threshold > 0
        Index = sort([find(Y >= Threshold);find(Y <= -Threshold)]);
        Mask(Index) = 1;
    elseif Threshold < 0
        Mask(Y<0) = 1;
    end
    
    V.fname = OutputName;
    spm_write_vol(V,Mask);
end

