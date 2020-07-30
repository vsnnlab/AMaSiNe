profile on
close all; clear all
addpath(genpath(pwd));

warning('off')
load('ANO_roi_edge');
load('Step_4_Angle_Search_Result');
load('Step_4to5_Cell_Detection_Result.mat') %% Load Cell detection From STEP4to5
STEP_0_Parameters;

img_name=Img_filename_list;
mkdir Image_Analysed
mkdir Image_Analysed_ROI_absent

elastixParams = {
    strcat(pwd,'\Core_Functions\elastix_functions\warping_parameters_Affine.txt'),...
    strcat(pwd,'\Core_Functions\elastix_functions\warping_parameters_BSpline.txt')}; %%% elastix

%% Prepartion : indices for both the anchor and non-anchor imgs
img_idx=anc_img_IDs;
ap_found=max_APpos_stage_final;
img_AP=[];
if length(img_idx) == 1
    img_AP = ap_found;
else
    for img_ID=1:length(img_idx)-1
        img_AP=[img_AP, linspace(ap_found(img_ID),ap_found(img_ID+1),...
            img_idx(img_ID+1)-img_idx(img_ID)+1)];
    end
end

img_AP=round(unique(img_AP));
if Slice_AP_orPA==1
    img_idx=min(img_idx):max(img_idx);
else
    img_idx=max(img_idx):-1:min(img_idx);
end

%% LOAD ATLAS
if strcmp(Structure_stain,'DAPI') || strcmp(Structure_stain,'Nissl')
    [VOL, metaVOL] = nrrdread('ara_nissl_25_2017.nrrd');
    VOL=uint8(uint16(rot90(permute(VOL,[3 1 2]),3))/(2^8));
elseif strcmp(Structure_stain,'AutoF')
    load('AutofluoresenceAtlas.mat');
    VOL=rot90(permute(VOL,[3 1 2]),3);
else
    error('Staining method not recognized');
end

yaw_found=yaw_stage5_max;
pitch_found=pitch_stage5_max;

tform_yaw=transform_matrix_0822( yaw_found,[0 1 0]);
tform_pitch=transform_matrix_0822(pitch_found,[1 0 0] );
tform_combined=mtimes(tform_yaw, tform_pitch);
tf_atlas= affine3d(tform_combined);
VOL_rot=imwarp(VOL,tf_atlas,'cubic');
ANO_rot=imwarp(ANO_roi_edge,tf_atlas,'nearest');


%% MAIN PART
parfor_progress(length(img_AP));
errorneous=false(1,length(img_AP));
img_essence(1:length(img_AP))=struct;

parfor img_ID=1:length(img_AP)
    img_ID
    if ~isempty(img_info(img_idx(img_ID)).slice_window)
        
        img_essence(img_ID).img_AP_pos= img_AP(img_ID);
        
        current_ap=img_AP(img_ID)-size(VOL,3)/2;
        current_ap=round(size(VOL_rot,3)/2+current_ap*...
            cosd(pitch_found)*cosd(yaw_found));
        
        %%%%%%%%%%%%%%%%%%%%%%%% Atlas Slice Prepartion for Transformation %%%%%%%%%%%%%%%%%%%%
        tform_general_resc_factor=1;
        img_ref=uint8(squeeze(VOL_rot(:,:,current_ap)));
        img_ref=imadjust(img_ref,stretchlim(img_ref,0),[0 1]);
        img_ref=imresize(img_ref,tform_general_resc_factor);
        img_ref=padarray(img_ref,round([3000 3000]/...
            (ref_atlas_vox_res/tform_general_resc_factor)));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Get Annotated Slice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        img_ANO=uint8(squeeze(ANO_rot(:,:,current_ap)));
        img_ANO=padarray(img_ANO,round([3000 3000]/(ref_atlas_vox_res)));
        img_ANO_downscale=img_ANO;
        img_ANO=imresize(img_ANO,(ref_atlas_vox_res/xy_pix),...
            'method','nearest','Antialiasing',false);
        
        %%%%%%%%%%%%%%%%%%%%%%%% Slice Preparation for Transformation %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        img_act=imread(img_name{img_idx(img_ID), Color_Channel_Structure});
        try
            img_act=rgb2gray(img_act);
        end
        
        img_act_pad=zeros(size(img_act));
        img_act=img_act(img_info(img_idx(img_ID)).slice_window(1):...
            img_info(img_idx(img_ID)).slice_window(2),...
            img_info(img_idx(img_ID)).slice_window(3):...
            img_info(img_idx(img_ID)).slice_window(4));
        img_act_pad=img_act_pad(img_info(img_idx(img_ID)).slice_window(1):...
            img_info(img_idx(img_ID)).slice_window(2),...
            img_info(img_idx(img_ID)).slice_window(3):...
            img_info(img_idx(img_ID)).slice_window(4));
        img_act_pad(img_info(img_idx(img_ID)).bnd_pix_ind)=1;
        img_act_pad=(imfill(img_act_pad));
        img_act_pad=uint8(logical(img_act_pad));
        img_act=img_act.*(img_act_pad);
        img_act=padarray(img_act,round([3000 3000]/(xy_pix)));
        img_act=imadjust(img_act,stretchlim(img_act,0.00),[0 1]);
        img_act=img_act+( img_act-imgaussfilt(img_act,0.5*201,...
            'FilterSize',[3 3]*603,'FilterDomain','frequency'))*5;
        img_act=imresize(img_act,xy_pix/ref_atlas_vox_res*...
            tform_general_resc_factor);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Find Transformation Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elastixDir_currImg = fullfile(pwd,...
            strcat('\img_warping_log\img_ID_',num2str(img_ID)));  %%% elastix
        mkdir(elastixDir_currImg);
        
        try
            [matched_img,transform_params]=elastix(img_act, img_ref,elastixDir_currImg,elastixParams);
            matched_img=uint8(matched_img);
        end
        
        if isempty(matched_img) || isempty(transform_params)
            errorneous(img_ID)=true;
            continue;
        end
        
        img_essence(img_ID).transform_params_downscaled=transform_params;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%% Get DownScaled Img for 3D Rec. %%%%%%%%%%%%%%%%%%%%%%%
        
        img_act_downscale = matched_img;
        %%% Match size of images
        whos_bigger=size(img_ANO_downscale)-size(img_act_downscale);
        if whos_bigger(1)>=0; x_dim=size(img_ANO_downscale,1);
        else; x_dim=size(img_act_downscale,1); end
        if whos_bigger(2)>=0; y_dim=size(img_ANO_downscale,2);
        else; y_dim=size(img_act_downscale,2); end
        
        img_act_sizeMatch=[x_dim,y_dim]-size(img_act_downscale);
        img_act_sizeMatch(img_act_sizeMatch<0)=0;
        img_ANO_sizeMatch=[x_dim,y_dim]-size(img_ANO_downscale);
        img_ANO_sizeMatch(img_ANO_sizeMatch<0)=0;
        img_act_downscale=padarray(img_act_downscale,...
            img_act_sizeMatch,'post');
        
        img_essence(img_ID).transformed_img_downscaled=...
            img_act_downscale;
        
        %%%%%%%%%%%% Warp Images in Their Original Pix Resolution %%%%%%%%%%%%%%%%%%%%%%%
        
        %%% modify transform parameters for warping IMAGES %%%
        mov_scale_factor=(ref_atlas_vox_res/xy_pix);
        
        transform_params_img_ori_scale=...
            img_essence(img_ID).transform_params_downscaled
        transform_params_img_ori_scale.TransformParameters{1, 1}.Size=...
            fliplr(size(img_ANO)); %% elastix vs matlab...
        transform_params_img_ori_scale.TransformParameters{1, 1}.TransformParameters(5:6)  = ...
            transform_params_img_ori_scale.TransformParameters{1, 1}.TransformParameters(5:6)*mov_scale_factor;
        transform_params_img_ori_scale.TransformParameters{1, 1}.CenterOfRotationPoint  = ...
            round((fliplr(size(img_ref))/2)*mov_scale_factor);
        transform_params_img_ori_scale.TransformParameters{1, 2}.Size=...
            fliplr(size(img_ANO)); %% elastix vs matlab...
        transform_params_img_ori_scale.TransformParameters{1, 2}.GridSpacing  = ...
            transform_params_img_ori_scale.TransformParameters{1, 2}.GridSpacing*mov_scale_factor;
        transform_params_img_ori_scale.TransformParameters{1, 2}.GridOrigin  = ...
            transform_params_img_ori_scale.TransformParameters{1, 2}.GridOrigin*mov_scale_factor;
        transform_params_img_ori_scale.TransformParameters{1, 2}.TransformParameters  = ...
            transform_params_img_ori_scale.TransformParameters{1, 2}.TransformParameters*mov_scale_factor;
        
        %%% Recalculate transform parameters for warping CELL COORDINATES %%%
        transform_params_cellPos_downscaled = ...
            invertElastixTransform(elastixDir_currImg);
        %%%NOTE: I know the line above seems redundant and costs you
        %%%extra few minutes for recalculation, but transforming the
        %%%cell coordinates using "transform_params_ori_scale" (used for
        %%%image warping) kept giving me weird results. It also doesn't
        %%%work the other way round (i.e.img warping with
        %%%"transform_params_cellPos_downscaled")
        
        img_essence(img_ID).transform_params_img_ori_scale= ...
           transform_params_img_ori_scale;
        img_essence(img_ID).transform_params_cellPos_downscaled= ...
            transform_params_cellPos_downscaled;
        
        
        %%% Load/Warp Images and Cell Pos in Original Resolution %%%
        for color_ch_ID=1:length(Color_Channel_Interest)
            cells_in_orig_img = ...
                cell_detection_rs(img_ID).Color_Cells(color_ch_ID).cell_locations;
            
            img_Color=(imread(img_name{img_idx(img_ID),...
                Color_Channel_Interest(color_ch_ID)}));
            try; img_Color=rgb2gray(img_Color); end
            img_Color=img_Color(img_info(img_idx(img_ID)).slice_window(1):...
                img_info(img_idx(img_ID)).slice_window(2),...
                img_info(img_idx(img_ID)).slice_window(3):...
                img_info(img_idx(img_ID)).slice_window(4));
            img_Color=img_Color.*img_act_pad;
            img_Color=padarray(img_Color,round([3000 3000]/(xy_pix)));
            
            
            [img_warped_no_scale,~]=transformix(img_Color,...
                transform_params_img_ori_scale);
            img_warped_no_scale=uint8(img_warped_no_scale);
            
            [cell_detected_all,~]=transformix(cells_in_orig_img/mov_scale_factor, ...
                transform_params_cellPos_downscaled);
            cell_detected_all=round(cell_detected_all.OutputPoint*...
                mov_scale_factor); %%% Looks a bit weird, but works....
            
            %%% Match size of annotation and warped images
            whos_bigger=size(img_ANO)-size(img_warped_no_scale);
            if whos_bigger(1)>=0; x_dim=size(img_ANO,1);
            else; x_dim=size(img_warped_no_scale,1); end
            if whos_bigger(2)>=0; y_dim=size(img_ANO,2);
            else; y_dim=size(img_warped_no_scale,2); end
            
            img_warped_sizeMatch=[x_dim,y_dim]-size(img_warped_no_scale);
            img_warped_sizeMatch(img_warped_sizeMatch<0)=0;
            img_ANO_sizeMatch=[x_dim,y_dim]-size(img_ANO);
            img_ANO_sizeMatch(img_ANO_sizeMatch<0)=0;
            img_ANO=padarray(img_ANO,img_ANO_sizeMatch,'post');
            img_warped_no_scale=padarray(img_warped_no_scale,...
                img_warped_sizeMatch,'post');
            
            %%% nope nope nope don't use the below...
            %                 [injection_site_pts] = injection_volume_pix( img_ch_thumbnail );
            %                 img_essence(img_ID).injection.Color(color_ch_ID).valid_pts=injection_site_pts;
            %                 img_essence(img_ID).injection.Color(color_ch_ID).image_thumbnail=img_ch_thumbnail;
            %                 img_essence(img_ID).Color_image_downscaled(color_ch_ID).img=img_warped_no_scale;
            
            
            %%%% for process monitoring (IMG) %%%%
            monitoring_img=img_warped_no_scale;
            monitoring_img_raw=monitoring_img;
            
            ROI_map=cat(3,...
                false(size(img_ANO)),logical(img_ANO),logical(img_ANO));
            
            monitoring_img=imadjust(monitoring_img,...
                stretchlim(monitoring_img,0),[0 1]);
            monitoring_img=cat(3,...
                monitoring_img,monitoring_img,monitoring_img);
            monitoring_img(ROI_map)=255;
            
            %%%% for process monitoring (CELL POS) %%%%
            if ~isempty(cell_detected_all)
                out_bnd_alpha = ref_boundarypad_0809_step5(...
                    img_warped_no_scale, xy_pix );
                out_bnd=inShape(out_bnd_alpha,...
                    cell_detected_all(:,2),cell_detected_all(:,1));
                cell_detected_all=cell_detected_all(~out_bnd,:);
            end
            
            if ~isempty(cell_detected_all)
                cell_detected_all_pos_ind=sub2ind(size(img_warped_no_scale),...
                    cell_detected_all(:,2),cell_detected_all(:,1));
            else
                cell_detected_all_pos_ind=[];
            end
            
            img_essence(img_ID).Color_Cells(color_ch_ID).cell_locations=...
                cell_detected_all;
            
            %%%% for process monitoring (All TOGETHER) %%%%
            image_analyzed_ROI=figure;   %%% with ROI boundaries overlaid
            imshow(monitoring_img,[]); hold on
            title(strcat({'Image Name : '}, ...
                img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)},...
                {'   - Ch No. '},num2str(color_ch_ID)),'Interpreter', 'none');
            if ~isempty(cell_detected_all)
                scatter(cell_detected_all(:,1),...
                    cell_detected_all(:,2),9,'r','filled')
            end
            
            image_analyzed=figure; %%% w/o ROI boundaries overlaid
            imshow(monitoring_img_raw,[]); hold on
            title(strcat({'Image Name : '}, ...
                img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)},...
                {'   - Ch No. '},num2str(color_ch_ID)),'Interpreter', 'none');
            if ~isempty(cell_detected_all)
                scatter(cell_detected_all(:,1),...
                    cell_detected_all(:,2),9,'r','filled')
            end
            
            try
                save_name=strcat('/Image_Analysed/',...
                    img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)},'.fig');
                saveas( image_analyzed_ROI ,[pwd  save_name]);
                close(image_analyzed_ROI)
                
                save_name=strcat('/Image_Analysed_ROI_absent/',...
                    img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)},'.fig');
                saveas( image_analyzed,[pwd  save_name]);
                close(image_analyzed)
            end
        end
        
    end
    parfor_progress;
end

delete parfor_progress.txt
save('Step_5_Cell_Detection_Result','img_essence','-v7.3');

profile viewer