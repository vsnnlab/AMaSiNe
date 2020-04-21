%%%%% Step 1 - Slice Boundary Detection
%%%%% Input - Your Images
%%%%% Output - 1. Rectangular Window Vertices of Each Slice
%%%            2. Point Coordinates of the Outermost Contour of Each Slice

clear all; close all; clc;
STEP_0_Parameters;
toolbox_chk
warning('off')
img_name=Img_filename_list(img_format);

h_progress = waitbar(0,'Slice Boundary Detection');

manual_list=[];

for img_ID=1:size(img_name, 1)
    
    img_structure=imread(img_name{img_ID, Color_Channel_Structure});
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
    thresh_raw_img=thresh_raw_img*1.5;
    
    img_structure = imbinarize(img_structure, thresh_raw_img);
    img_structure = 255*uint8(img_structure);
    
    if size(img_structure,1)*xy_pix<4000 && size(img_structure,2)*xy_pix<4000   %%% Neglect Images that are too small
        img_info(img_ID).slice_window=[];
        img_info(img_ID).bnd_pix_ind=[];
    else
        %%%%%%% 1. Rough CONTRAST ADJUSTMENT of RAW Images and Manual Removal of Bubbles %%%
        I_contrast=img_structure;
        %%%%%%% 2. SLICE BOUNDARY DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [ BWoutline, BWobject_pad] = SliceBoundaryDetection( I_contrast,  xy_pix);
        I_contrast(~BWobject_pad)=0; %% zerofy regions outside the slice of interest
        [window_bnd] = CutOutBlank(BWobject_pad);   %% cut out the regions outside slice to save memory resources
        
        if ~isempty(window_bnd) && (window_bnd(2)-window_bnd(1))*(window_bnd(4)-window_bnd(3))*xy_pix^2 > 2000^2
            BW_boundary_coord_cut=find(BWoutline(window_bnd(1):window_bnd(2),window_bnd(3):window_bnd(4)));
            img_info(img_ID).slice_window=window_bnd;
            img_info(img_ID).bnd_pix_ind=BW_boundary_coord_cut;
        else
            disp(strcat('Image No.',num2str(img_ID),' is to be processed manually'));
            manual_list=[manual_list img_ID];
        end
    end
    waitbar(img_ID / length(img_name),h_progress,...
        strcat('Slice Boundary Detection: ',{num2str(img_ID)},{' of '},{num2str(length(img_name))},{' images done'}))
    
end

close(h_progress)
if ~isempty(manual_list)
    for manual_ii=1:size(manual_list,2)
        
        img_structure=imread(img_name{manual_list(manual_ii), Color_Channel_Structure});
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
            'FilterDomain','frequency'))*30;
        
        thresh_raw_img=graythresh(img_structure);
        thresh_raw_img=thresh_raw_img*0.7;
        
        [img_structure, BW_mask ] = manual_bubble_removal(img_structure);
        img_structure=(img_structure.*uint8(BW_mask));
        
        img_structure = imbinarize(img_structure, thresh_raw_img);
        img_structure = 255*uint8(img_structure);
        
        [ BWoutline, BWobject_pad] = SliceBoundaryDetection(img_structure,  xy_pix);
        img_structure(~BWobject_pad)=0; %% zerofy regions outside the slice of interest
        [window_bnd] = CutOutBlank(BWobject_pad);   %% cut out the regions outside slice to save memory resources
        
        if ~isempty(window_bnd)
            BW_boundary_coord_cut=find(BWoutline(window_bnd(1):window_bnd(2),window_bnd(3):window_bnd(4)));
            img_info(manual_list(manual_ii)).slice_window=window_bnd;
            img_info(manual_list(manual_ii)).bnd_pix_ind=BW_boundary_coord_cut;
        else
            disp(strcat('Image No.',num2str(manual_list(manual_ii)),' is excluded'));
        end
        
    end
end

save('Step1_Outline_result','img_info','-v7.3');

