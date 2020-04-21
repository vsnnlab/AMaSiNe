clear all; clc;
load('Step1_Outline_result.mat')
warning('off')
STEP_0_Parameters;
img_name=Img_filename_list;

if ~isempty(img_IDs_reBoundary)
    for manual_ii=1:size(img_IDs_reBoundary,2)
        
        img_structure=imread(img_name{img_IDs_reBoundary(manual_ii), ...
            Color_Channel_Structure});
        try
            img_structure=rgb2gray(img_structure);
        end
        img_structure=imadjust(img_structure,stretchlim(img_structure,0.01),[0 1]);
        
        
        unsharp_filter_size=round([400, 400]/xy_pix);
        if sum(mod(unsharp_filter_size,2))==0
            unsharp_filter_size=unsharp_filter_size+1;
        end
        
        img_structure=img_structure+...
            (img_structure-imgaussfilt((img_structure),75/xy_pix,'FilterSize',unsharp_filter_size,...
            'FilterDomain','frequency'))*5;
        
        thresh_raw_img=graythresh(img_structure);
        thresh_raw_img=thresh_raw_img*threshold_scale;
        
        [img_structure, BW_mask ] = manual_bubble_removal(img_structure);
        img_structure=(img_structure.*uint8(BW_mask));
        
        img_structure = imbinarize(img_structure, thresh_raw_img);
        img_structure = 255*uint8(img_structure);
        
        [ BWoutline, BWobject_pad] = SliceBoundaryDetection(img_structure,  xy_pix);
        img_structure(~BWobject_pad)=0; %% zerofy regions outside the slice of interest
        [window_bnd] = CutOutBlank(BWobject_pad);   %% cut out the regions outside slice to save memory resources
        
        if ~isempty(window_bnd)
            BW_boundary_coord_cut=find(BWoutline(window_bnd(1):window_bnd(2),window_bnd(3):window_bnd(4)));
            img_info(img_IDs_reBoundary(manual_ii)).slice_window=window_bnd;
            img_info(img_IDs_reBoundary(manual_ii)).bnd_pix_ind=BW_boundary_coord_cut;
        else
            disp(strcat('Image No.',num2str(img_IDs_reBoundary(manual_ii)),' is excluded'));
            
            img_info(img_IDs_reBoundary(manual_ii)).slice_window=[];
            img_info(img_IDs_reBoundary(manual_ii)).bnd_pix_ind=[];
        end
        
        figure; imshowpair(img_structure,  BWobject_pad,'montage')
    end
end

save('Step1_Outline_result','img_info','-v7.3');
