profile on
close all
clear all

warning('off')
STEP_0_Parameters;
load('ANO_roi_edge');
load('Step_4_Angle_Search_Result');
img_name=Img_filename_list;
mkdir Image_Analysed
mkdir Image_Analysed_ROI_absent

%% Prepartion : indices for both the anchor and non-anchor imgs
img_idx=anc_img_IDs;
ap_found=max_APpos_stage_final;
img_AP=[];
if length(img_idx) == 1
    img_AP = ap_found;
else
    for img_ID=1:length(img_idx)-1
        img_AP=[img_AP, linspace(ap_found(img_ID),ap_found(img_ID+1),img_idx(img_ID+1)-img_idx(img_ID)+1)];
    end
end

img_AP=round(unique(img_AP));
if Slice_AP_orPA==1
    img_idx=min(img_idx):max(img_idx);
else
    img_idx=max(img_idx):-1:min(img_idx);
end
xy_pix_resc_factor = 1;
%% MAIN PART
for img_ID=1:length(img_AP)
    disp(['Cell Detection: ', num2str(img_ID), '/', num2str(length(img_AP))]);
    img_ID
    for color_ch_ID=1:length(Color_Channel_Interest)
        if ~isempty(img_info(img_idx(img_ID)).slice_window)
            cell_detection_rs(img_ID).img_AP_pos= img_AP(img_ID);
            
            img_Color=(imread(img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)}));
            try
                img_Color=rgb2gray(img_Color);
            end                        
            img_Color_pad=zeros(size(img_Color));      
            img_Color=img_Color(img_info(img_idx(img_ID)).slice_window(1):img_info(img_idx(img_ID)).slice_window(2),...
                img_info(img_idx(img_ID)).slice_window(3):img_info(img_idx(img_ID)).slice_window(4));
            
            img_Color_pad=img_Color_pad(img_info(img_idx(img_ID)).slice_window(1):img_info(img_idx(img_ID)).slice_window(2),...
                img_info(img_idx(img_ID)).slice_window(3):img_info(img_idx(img_ID)).slice_window(4));
            img_Color_pad(img_info(img_idx(img_ID)).bnd_pix_ind)=1;
            img_Color_pad=(imfill(img_Color_pad));
            img_Color_pad=uint8(logical(img_Color_pad));
        
            img_Color=img_Color.*img_Color_pad;
            img_Color=padarray(img_Color,round([3000 3000]/(xy_pix)));
            
            atlas_exp_scale_fac=(25/xy_pix)*xy_pix_resc_factor;
            
            img_Color=imresize(img_Color,xy_pix_resc_factor);
            
            %%% detect cells across the whole slice image
            cell_detected_all= SomaDetection0827 (img_Color, xy_pix, soma_radius,cell_det_thresh);
            
            if ~isempty(cell_detected_all)
                out_bnd_alpha = ref_boundarypad_0809_step5( img_Color, xy_pix/xy_pix_resc_factor );
                out_bnd=inShape(out_bnd_alpha,cell_detected_all(:,2),cell_detected_all(:,1));
                cell_detected_all=cell_detected_all(~out_bnd,:);
            end
            
            cell_detected_all= round(cell_detected_all);
            
            if ~isempty(cell_detected_all)
                cell_detected_all_pos_ind=sub2ind(size(img_Color),...
                    cell_detected_all(:,2),cell_detected_all(:,1));
            else
                cell_detected_all_pos_ind=[];
            end
            
            cell_detection_rs(img_ID).Color_Cells(color_ch_ID).cell_locations=cell_detected_all;
            
            figure; imshow(img_Color);
            hold on;
            plot(cell_detected_all(:,1), cell_detected_all(:,2), 'ro');
            
            image_analyzed=figure;
            imshow(img_Color,[]); hold on
            title(strcat({'Image Name : '}, img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)},...
                {'   - Ch No. '},num2str(color_ch_ID)),'Interpreter', 'none');
            if ~isempty(cell_detected_all)
                scatter(cell_detected_all(:,1),cell_detected_all(:,2),9,'r','filled')
            end
                        
            save_name=strcat('/Image_Analysed_ROI_absent/',img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)},'.fig');
            saveas( image_analyzed,[pwd  save_name]);
            close(image_analyzed)
            
% % %             [B,RB]=imwarp(A,tform);
% % %             [xdataT,ydataT]=transformPointsForward(tform,xdata,ydata);
% % %             [xdataI,ydataI]=worldToIntrinsic(RB,xdataT,ydataT);
        end
    end
end
save('Step_4to5_Cell_Detection_Result','cell_detection_rs','-v7.3');