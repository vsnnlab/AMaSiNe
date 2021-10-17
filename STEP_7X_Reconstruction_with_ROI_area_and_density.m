%% Written by Woochul Choi at 21/10/01
%%% Aim of this script is to
%%% 1) Calculate the area of each ROIs & calculate the density! (#/mm^2) in a single slice
%%% 2) If multiple slice contributes to a specific ROI, then calculate the density (#/mm^2) in each slice and later sum them up.
clear all; close all
addpath(genpath(pwd));
STEP_0_Parameters;
STEP4_rs = matfile('Step_4_Angle_Search_Result.mat');
load('Step_5_Cell_Detection_Result.mat')
load ANO.mat;
load Step5_ANO_info.mat;
load ROI_region_pixels.mat;
STEP_0_Parameters;
%% Load Atlases
[VOL, metaVOL] = nrrdread('ara_nissl_25_2017.nrrd');
VOL=double(rot90(permute(VOL,[3 1 2]),3));
VOL=padarray(VOL,round([3000 3000]/ref_atlas_vox_res));

yaw_found=STEP4_rs.yaw_stage5_max; pitch_found=STEP4_rs.pitch_stage5_max;
tform_yaw=transform_matrix_0822( yaw_found,[0 1 0]);
tform_pitch=transform_matrix_0822(pitch_found,[1 0 0] );
tform_combined=mtimes(tform_yaw, tform_pitch);
tf_atlas= affine3d(tform_combined);
VOL_rot=imwarp(VOL,tf_atlas,'nearest');
ANO_rot=imwarp(ANO,tf_atlas,'nearest');

%%% Rotational Matrix
transformation_x=inv([cosd(-yaw_found) 0 sind(-yaw_found);0 1 0; -sind(-yaw_found) 0 cosd(-yaw_found)]);
transformation_y=inv([1 0 0; 0 cosd(-pitch_found) -sind(-pitch_found); 0 sind(-pitch_found) cosd(-pitch_found)]);
transform_mat=mtimes(transformation_x,transformation_y);
rot_center=[size(VOL_rot,2) size(VOL_rot,1) size(VOL,3)]/2;
%%% Shift Matrix
shift_calibration=(size(VOL_rot)-size(VOL))/2;
shift_calibration= [shift_calibration(2), -shift_calibration(1), 0];

EachSlice_Results = struct;
EachSlice_Results(1:size(img_essence,2))=struct;
for img_ID=1:size(img_essence,2)
    img_ID
    img_AP_pos= img_essence(img_ID).img_AP_pos;
    temp_Left_ROI_Area = nan(1, numel(region_ID_list));
    temp_Right_ROI_Area = nan(1, numel(region_ID_list));
    temp_Detected_Cell_Left_ROI = nan(1, numel(region_ID_list));
    temp_Detected_Cell_Right_ROI = nan(1, numel(region_ID_list));
    
    EachSlice_Results(img_ID).img_AP_pos = img_AP_pos;
    EachSlice_Results(img_ID).ROI_names = region_name_list;
    EachSlice_Results(img_ID).Left_ROI_Area = temp_Left_ROI_Area;
    EachSlice_Results(img_ID).Right_ROI_Area = temp_Right_ROI_Area;
    
    EachSlice_Results(img_ID).Left_ROI_NumCell = temp_Detected_Cell_Left_ROI;
    EachSlice_Results(img_ID).Right_ROI_NumCell = temp_Detected_Cell_Right_ROI;
    %% Counting Ch1 Neurons
    cell_loc_ch1=img_essence(img_ID).Color_Cells(1).cell_locations*(xy_pix)/ref_atlas_vox_res;
    cell_loc_ch1_3D=[cell_loc_ch1  img_AP_pos*ones(size(cell_loc_ch1(:,1)))];
    rotated_cell_loc_ch1 = mtimes(transform_mat, (cell_loc_ch1_3D-rot_center)')'+rot_center;
    shifted_cell_loc_ch1 = rotated_cell_loc_ch1;
    shifted_cell_loc_ch1(:,2) = -rotated_cell_loc_ch1(:,2);
    shifted_cell_loc_ch1(:,3) = -rotated_cell_loc_ch1(:,3);
    shifted_cell_loc_ch1 = 25*(shifted_cell_loc_ch1-[348,-116,-214]-shift_calibration);
    
    if Slice_AP_orPA==-1
        compensated_cell_loc_ch1 = shifted_cell_loc_ch1;
    else
        compensated_cell_loc_ch1 = -shifted_cell_loc_ch1;
    end
    
    Left_Cell_Index = compensated_cell_loc_ch1(:,1)>0;
    Right_Cell_Index = compensated_cell_loc_ch1(:,1)<=0;
    
    [D_tempo_troubleshoot,ano_pix_idx] = pdist2(ROI_shell_coord(:,1:3),shifted_cell_loc_ch1,'euclidean','Smallest',1);
    cell_ROI_ID = ROI_shell_coord(ano_pix_idx,4);
    
    temp_Detected_Cell_Left_ROI = zeros(1, numel(region_ID_list));
    temp_Detected_Cell_Right_ROI = zeros(1, numel(region_ID_list));
    for roi_iter=1:size(region_ID_list,2)
        cells_idx_current_ROI=find(ismember(cell_ROI_ID,region_ID_list(roi_iter).list));
        boolIndex = false(size(cell_ROI_ID,1), 1);
        boolIndex(cells_idx_current_ROI) = 1;
        
        num_L = sum(boolIndex & Left_Cell_Index);
        num_R = sum(boolIndex & Right_Cell_Index);
        
        temp_Detected_Cell_Left_ROI(1, roi_iter) = num_L;
        temp_Detected_Cell_Right_ROI(1, roi_iter) = num_R;
    end
    
    EachSlice_Results(img_ID).Left_ROI_NumCell = temp_Detected_Cell_Left_ROI;
    EachSlice_Results(img_ID).Right_ROI_NumCell = temp_Detected_Cell_Right_ROI;
    %% Calculating ROI Area
    current_ap=img_AP_pos-size(VOL,3)/2;
    current_ap=round(size(VOL_rot,3)/2+current_ap*...
        cosd(pitch_found)*cosd(yaw_found));
    
    img_ANO=(squeeze(ANO_rot(:,:,current_ap)));
    img_ANO=padarray(img_ANO,round([3000 3000]/(ref_atlas_vox_res))); %% Annotated atlas image;
    
    %%% Define Left and Right from whole area
    wholearea_ROI= find(~isnan(img_ANO));
    [row,col] = ind2sub(size(img_ANO), wholearea_ROI);
    pos_3D=[col, row, img_AP_pos*ones(numel(col), 1)];
    rotated_pos = mtimes(transform_mat, (pos_3D-rot_center)')'+rot_center;
    shifted_pos = rotated_pos;
    shifted_pos(:,2) = -rotated_pos(:,2);
    shifted_pos(:,3) = -rotated_pos(:,3);
    shifted_pos = 25*(shifted_pos-[348,-116,-214]-shift_calibration);
    
    Left_Hemi_Index = shifted_pos(:,1)>0;
    Right_Hemi_Index = shifted_pos(:,1)<=0;
    
    L_ROI_Area = nan(1,size(region_ID_list,2));
    R_ROI_Area = nan(1,size(region_ID_list,2));
    for roi_iter=1:size(region_ID_list,2)
        position_idx_current_ROI=find(ismember(img_ANO,region_ID_list(roi_iter).list));
        boolIndex = false(numel(img_ANO), 1);
        boolIndex(position_idx_current_ROI) = 1;
        
        L_ROI_Area(roi_iter) = sum(Left_Hemi_Index & boolIndex);
        R_ROI_Area(roi_iter) = sum(Right_Hemi_Index & boolIndex);
    end
    
    L_ROI_Area = L_ROI_Area * (ref_atlas_vox_res/1000)^2;
    R_ROI_Area = R_ROI_Area * (ref_atlas_vox_res/1000)^2;
    
    EachSlice_Results(img_ID).Left_ROI_Area = L_ROI_Area;
    EachSlice_Results(img_ID).Right_ROI_Area = R_ROI_Area;
    
    EachSlice_Results(img_ID).Left_ROI_Density = temp_Detected_Cell_Left_ROI./ L_ROI_Area;
    EachSlice_Results(img_ID).Right_ROI_Density = temp_Detected_Cell_Right_ROI./ R_ROI_Area;
end
%% Whole Result Summary
WholeSlice_Area_Density_RS = struct;
WholeSlice_Area_Density_RS.ROI_names = region_name_list;

L_ROI_Area_Whole = zeros(1, numel(region_name_list));
R_ROI_Area_Whole = zeros(1, numel(region_name_list));
L_NumCell_Whole = zeros(1, numel(region_name_list));
R_NumCell_Whole = zeros(1, numel(region_name_list));
for img_ID = 1:size(img_essence,2)
    L_ROI_Area_Whole = L_ROI_Area_Whole + EachSlice_Results(img_ID).Left_ROI_Area;
    R_ROI_Area_Whole = R_ROI_Area_Whole + EachSlice_Results(img_ID).Right_ROI_Area;
    L_NumCell_Whole = L_NumCell_Whole + EachSlice_Results(img_ID).Left_ROI_NumCell;
    R_NumCell_Whole = R_NumCell_Whole + EachSlice_Results(img_ID).Right_ROI_NumCell;
end

WholeSlice_Area_Density_RS.Left_ROI_Area = L_ROI_Area_Whole;
WholeSlice_Area_Density_RS.Right_ROI_Area = R_ROI_Area_Whole;
WholeSlice_Area_Density_RS.Left_ROI_NumCell = L_NumCell_Whole;
WholeSlice_Area_Density_RS.Right_ROI_NumCell = R_NumCell_Whole;
WholeSlice_Area_Density_RS.Left_ROI_Density = L_NumCell_Whole ./ L_ROI_Area_Whole;
WholeSlice_Area_Density_RS.Right_ROI_Density = R_NumCell_Whole ./ R_ROI_Area_Whole;

save('STEP7X_ROI_Area_and_Density_result.mat', 'WholeSlice_Area_Density_RS', 'EachSlice_Results');
%% Write Excel
disp('********All Matlab File is saved********');
disp('*****Now Excel file will be written*****');
disp('');

maindir = cd;
% Your_File_Name = input('Enter your filename: ', 's');
% if isempty(Your_File_Name)
    Your_File_Name = 'STEP7X_ROI_Area_and_Density_Result.xlsx';
% else
%     Your_File_Name = [Your_File_Name, '.xlsx'];
% end
disp(['Excel file will be saved as ', Your_File_Name]);
excel_name = Your_File_Name;
excel_export_filename = [maindir, '\', excel_name];

%%% Title
writematrix('ROI area and density from whole brain', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', 'A1');

%%% ROI Names
writematrix('ROI Name', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', 'A2');
ROI_Name_Cell = squeeze(struct2cell(region_name_list));
Range = ['A4:A', num2str(numel(region_name_list)+3)];
writecell(ROI_Name_Cell, excel_export_filename, 'Sheet', 'WholeBrain', 'Range', Range);

%%% Detected Cell Number
writematrix('Detected Cell #', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', 'C2');
writematrix('Left hemi', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', 'C3');
writematrix('Right hemi', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', 'D3');
Range = ['C4:C', num2str(numel(region_name_list)+3)];
writematrix(WholeSlice_Area_Density_RS.Left_ROI_NumCell', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', Range);
Range = ['D4:D', num2str(numel(region_name_list)+3)];
writematrix(WholeSlice_Area_Density_RS.Right_ROI_NumCell', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', Range);

%%% ROI Area (mm^2)
writematrix('ROI Area Sum (mm^2 * # Slice)', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', 'F2');
writematrix('Left hemi', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', 'F3');
writematrix('Right hemi', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', 'G3');
Range = ['F4:F', num2str(numel(region_name_list)+3)];
writematrix(WholeSlice_Area_Density_RS.Left_ROI_Area', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', Range);
Range = ['G4:G', num2str(numel(region_name_list)+3)];
writematrix(WholeSlice_Area_Density_RS.Right_ROI_Area', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', Range);

%%% ROI Density (# / mm^2)
writematrix('ROI Density (# / (mm^2 * # Slice))', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', 'I2');
writematrix('Left hemi', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', 'I3');
writematrix('Right hemi', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', 'J3');
Range = ['I4:I', num2str(numel(region_name_list)+3)];
writematrix(WholeSlice_Area_Density_RS.Left_ROI_Density', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', Range);
Range = ['J4:J', num2str(numel(region_name_list)+3)];
writematrix(WholeSlice_Area_Density_RS.Right_ROI_Density', excel_export_filename, 'Sheet', 'WholeBrain', 'Range', Range);

disp('*****Now Whole Brain Result is done*****');

%%% Each slice as a separate sheet
disp('*****Now Writing Each Slice Data*****');
for img_ID = 1:size(img_essence,2)
    img_ID
    SheetName = ['Slice #', num2str(img_ID)];
    %%% Title
    writematrix(['Slice #', num2str(img_ID), 'ROI area and density'], excel_export_filename, 'Sheet', SheetName, 'Range', 'A1');
    
    %%% ROI Names
    writematrix('ROI Name', excel_export_filename, 'Sheet', SheetName, 'Range', 'A2');
    ROI_Name_Cell = squeeze(struct2cell(region_name_list));
    Range = ['A4:A', num2str(numel(region_name_list)+3)];
    writecell(ROI_Name_Cell, excel_export_filename, 'Sheet', SheetName, 'Range', Range);
    
    %%% Detected Cell Number
    writematrix('Detected Cell #', excel_export_filename, 'Sheet', SheetName, 'Range', 'C2');
    writematrix('Left hemi', excel_export_filename, 'Sheet', SheetName, 'Range', 'C3');
    writematrix('Right hemi', excel_export_filename, 'Sheet', SheetName, 'Range', 'D3');
    Range = ['C4:C', num2str(numel(region_name_list)+3)];
    writematrix(EachSlice_Results(img_ID).Left_ROI_NumCell', excel_export_filename, 'Sheet', SheetName, 'Range', Range);
    Range = ['D4:D', num2str(numel(region_name_list)+3)];
    writematrix(EachSlice_Results(img_ID).Right_ROI_NumCell', excel_export_filename, 'Sheet', SheetName, 'Range', Range);
    
    %%% ROI Area (mm^2)
    writematrix('ROI Area Sum (mm^2)', excel_export_filename, 'Sheet', SheetName, 'Range', 'F2');
    writematrix('Left hemi', excel_export_filename, 'Sheet', SheetName, 'Range', 'F3');
    writematrix('Right hemi', excel_export_filename, 'Sheet', SheetName, 'Range', 'G3');
    Range = ['F4:F', num2str(numel(region_name_list)+3)];
    writematrix(EachSlice_Results(img_ID).Left_ROI_Area', excel_export_filename, 'Sheet', SheetName, 'Range', Range);
    Range = ['G4:G', num2str(numel(region_name_list)+3)];
    writematrix(EachSlice_Results(img_ID).Right_ROI_Area', excel_export_filename, 'Sheet', SheetName, 'Range', Range);
    
    %%% ROI Density (# / mm^2)
    writematrix('ROI Density (#/mm^2)', excel_export_filename, 'Sheet', SheetName, 'Range', 'I2');
    writematrix('Left hemi', excel_export_filename, 'Sheet', SheetName, 'Range', 'I3');
    writematrix('Right hemi', excel_export_filename, 'Sheet', SheetName, 'Range', 'J3');
    Range = ['I4:I', num2str(numel(region_name_list)+3)];
    writematrix(EachSlice_Results(img_ID).Left_ROI_Density', excel_export_filename, 'Sheet', SheetName, 'Range', Range);
    Range = ['J4:J', num2str(numel(region_name_list)+3)];
    writematrix(EachSlice_Results(img_ID).Right_ROI_Density', excel_export_filename, 'Sheet', SheetName, 'Range', Range);
end

disp('*****All Excel File is saved*****');