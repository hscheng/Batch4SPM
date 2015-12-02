function GenMask_MultiImgOverLap()
% http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0023175
% only suppot thresh surpass

    ImgList = ScanImg2Cell('Select Img Files');

    thresh = input('set the threshold : ');
    
    % set default thresh value
    if isempty(thresh)
        thresh = 0.1;
    end
    save ImgListInfo
    
    ImgNum = length(ImgList);

    for ii = 1:ImgNum
        V = spm_vol(ImgList{ii});
        Y = spm_read_vols(V);
        if ii == 1
            Ymask = Y > thresh;
        else
            Ymask = Ymask.*(Y > thresh);
        end
    end
    V.fname = ['Mask_overlap_thesh',num2str(thresh),'.nii'];
    spm_write_vol(V,Ymask);
end