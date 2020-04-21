clear all; close all
addpath(genpath(pwd));
load('Step_4_Angle_Search_Result.mat')
load('Step_6_Cell_Detection_Result.mat')
load ANO.mat;
load Step5_ANO_info.mat;
load ROI_region_pixels.mat;
STEP_0_Parameters;
xy_pix=xy_pix/xy_pix_resc_factor;

%% Load Atlases
[VOL, metaVOL] = nrrdread('ara_nissl_25_2017.nrrd');
VOL=double(rot90(permute(VOL,[3 1 2]),3));
VOL=padarray(VOL,round([3000 3000]/ref_atlas_vox_res));

yaw_found=yaw_stage5_max; pitch_found=pitch_stage5_max;
tform_yaw=transform_matrix_0822( yaw_found,[0 1 0]);
tform_pitch=transform_matrix_0822(pitch_found,[1 0 0] );
tform_combined=mtimes(tform_yaw, tform_pitch);
tf_atlas= affine3d(tform_combined);
VOL_rot=imwarp(VOL,tf_atlas,'nearest');

%% Cell-Sorting
cell_pos=[];
Color_ch=struct('cell_pos',cell_pos);
Color_ch=repmat(Color_ch,[length(Color_Channel_Interest),1]);
cell_locations=struct('Color_ch',Color_ch);

for img_ID=1:size(img_essence,2)
    for color_ch_id=1:length(Color_Channel_Interest)          
        try
            cell_loc=img_essence(img_ID).Color_Cells(color_ch_id).cell_locations*(xy_pix)/ref_atlas_vox_res;
            cell_loc=[cell_loc  img_essence(img_ID).img_AP_pos*ones(size(cell_loc(:,1)))];
            cell_locations.Color_ch(color_ch_id).cell_pos=...
                vertcat(cell_locations.Color_ch(color_ch_id).cell_pos, cell_loc);
        end
    end
end
            
%% Rotate data points 

transformation_x=inv([cosd(-yaw_found) 0 sind(-yaw_found);0 1 0; -sind(-yaw_found) 0 cosd(-yaw_found)]);
transformation_y=inv([1 0 0; 0 cosd(-pitch_found) -sind(-pitch_found); 0 sind(-pitch_found) cosd(-pitch_found)]);

transform_mat=mtimes(transformation_x,transformation_y);
rot_center=[size(VOL_rot,2) size(VOL_rot,1) size(VOL,3)]/2;

for color_ch_id=1:length(Color_Channel_Interest)
    if ~isempty(cell_locations.Color_ch(color_ch_id).cell_pos)
        for cell_id=1:size(cell_locations.Color_ch(color_ch_id).cell_pos,1)
            cell_locations.Color_ch(color_ch_id).cell_pos(cell_id,:)=...
                mtimes(transform_mat, (cell_locations.Color_ch(color_ch_id).cell_pos(cell_id,:) ...
                -rot_center)')'+rot_center;
        end
    end
end

%% Scale and shift data points 

shift_calibration=(size(VOL_rot)-size(VOL))/2;
shift_calibration= [shift_calibration(2), -shift_calibration(1), 0];


for color_ch_id=1:length(Color_Channel_Interest)
    if ~isempty( cell_locations.Color_ch(color_ch_id).cell_pos)
        for cell_id=1:size(cell_locations.Color_ch(color_ch_id).cell_pos,1)
            cell_locations.Color_ch(color_ch_id).cell_pos(cell_id,2)=...
                -cell_locations.Color_ch(color_ch_id).cell_pos(cell_id,2) ;
            
            cell_locations.Color_ch(color_ch_id).cell_pos(cell_id,3)=...
                -cell_locations.Color_ch(color_ch_id).cell_pos(cell_id,3) ;
            
            cell_locations.Color_ch(color_ch_id).cell_pos(cell_id,:)=...
                25*(cell_locations.Color_ch(color_ch_id).cell_pos(cell_id,:)...
                -[348,  -116,  -214]-shift_calibration);
        end
    end
end

%% Annotate Cells into ROIs 
for color_ch_id=1:length(Color_Channel_Interest)
    cells_color_ch=cell_locations.Color_ch(color_ch_id).cell_pos;
    total_no_cells=size(cells_color_ch,1);
    proc_blocks=1:1000:total_no_cells;
    proc_blocks=[proc_blocks total_no_cells];
    
    cell_ROI_id=nan(total_no_cells,1);
    D_troubleshoot=[];
    
    for block_ii=2:length(proc_blocks)
         if block_ii~=length(proc_blocks)
            cell_pos_tempo=cells_color_ch(proc_blocks(block_ii-1):proc_blocks(block_ii)-1,:);
            [D_tempo_troubleshoot,ano_pix_idx] = pdist2(ROI_shell_coord(:,1:3),cell_pos_tempo,'euclidean','Smallest',1);
            cell_ROI_id(proc_blocks(block_ii-1):proc_blocks(block_ii)-1)=ROI_shell_coord(ano_pix_idx,4);
            
            D_troubleshoot=[D_troubleshoot D_tempo_troubleshoot];
        else
            cell_pos_tempo=cells_color_ch(proc_blocks(block_ii-1):proc_blocks(block_ii),:);
            [D_tempo_troubleshoot,ano_pix_idx] = pdist2(ROI_shell_coord(:,1:3),cell_pos_tempo,'euclidean','Smallest',1);
            cell_ROI_id(proc_blocks(block_ii-1):proc_blocks(block_ii))=ROI_shell_coord(ano_pix_idx,4);
            
            D_troubleshoot=[D_troubleshoot D_tempo_troubleshoot];
        end

    end
    
    cell_locations.Color_ch(color_ch_id).cell_pos=[cell_locations.Color_ch(color_ch_id).cell_pos cell_ROI_id];
      for ii=1:size(region_ID_list,2)     
        
        cells_ROI(ii).ROI_name=region_name_list(ii).name;
        cells_idx_current_ROI=find(ismember(cell_ROI_id,region_ID_list(ii).list));
        cells_ROI_tempo=cell_locations.Color_ch(color_ch_id).cell_pos(cells_idx_current_ROI,1:3);
        
        if Slice_AP_orPA==1
            cells_ROI_tempo(:,1)=-cells_ROI_tempo(:,1);
        end
        cells_ROI_tempo_L_idx=find(cells_ROI_tempo(:,1)>0);
        cells_ROI_tempo_R_idx=find(cells_ROI_tempo(:,1)<=0);
        cells_ROI(ii).Color_ch(color_ch_id).RL(1).cells_pos=cells_ROI_tempo(cells_ROI_tempo_L_idx,:);
        cells_ROI(ii).Color_ch(color_ch_id).RL(2).cells_pos=cells_ROI_tempo(cells_ROI_tempo_R_idx,:);
      end  
end

save('Step_7_3D_Reconstructed_cells_Result','cells_ROI','-v7.3');
