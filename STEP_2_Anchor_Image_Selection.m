warning('off')
addpath(genpath(pwd))
img_name=Img_filename_list;
load('Step1_Outline_result.mat')
STEP_0_Parameters;

ref_resc=1/0.5;
downscaled_xy_pix=ref_atlas_vox_res*ref_resc;


for img_ID=1:size(img_info,2)
    if ~isempty(img_info(img_ID).slice_window)
        %%%% load raw image %%%%
        img_act=imread(img_name{img_ID, Color_Channel_Structure});
        
        try
            img_act=rgb2gray(img_act);
        end
        
        %%%% Leave the slice part only %%%%
        img_act_pad=zeros(size(img_act));
        
        img_act=img_act(img_info(img_ID).slice_window(1):img_info(img_ID).slice_window(2),...
            img_info(img_ID).slice_window(3):img_info(img_ID).slice_window(4));
        
        img_act_pad=img_act_pad(img_info(img_ID).slice_window(1):img_info(img_ID).slice_window(2),...
            img_info(img_ID).slice_window(3):img_info(img_ID).slice_window(4));
        img_act_pad(img_info(img_ID).bnd_pix_ind)=1;
        
        
        img_act_pad=(imfill(img_act_pad));
        
        img_act_pad=uint8(logical(img_act_pad));
        img_act=img_act.*(img_act_pad);
    
        
        %%%% Contrast adjustment and Edge-boosting by unsharpening %%%%
        I_adj=imadjust(img_act,stretchlim(img_act,0),[0 1]);

        I_adj= I_adj+( I_adj-imgaussfilt( I_adj,0.5*100,'FilterSize',[3 3]*101,...
            'FilterDomain','frequency'))*30; 
        I_adj=imresize(I_adj, (xy_pix)/(downscaled_xy_pix) );
        figure; imshow(I_adj)
       
        title(strcat('Img #.',num2str(img_ID)))
    end
end
