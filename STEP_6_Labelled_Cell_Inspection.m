function varargout = STEP_6_Labelled_Cell_Inspection(varargin)
% STEP_6_LABELLED_CELL_INSPECTION MATLAB code for STEP_6_Labelled_Cell_Inspection.fig
%      STEP_6_LABELLED_CELL_INSPECTION, by itself, creates a new STEP_6_LABELLED_CELL_INSPECTION or raises the existing
%      singleton*.
%
%      H = STEP_6_LABELLED_CELL_INSPECTION returns the handle to a new STEP_6_LABELLED_CELL_INSPECTION or the handle to
%      the existing singleton*.
%
%      STEP_6_LABELLED_CELL_INSPECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STEP_6_LABELLED_CELL_INSPECTION.M with the given input arguments.
%
%      STEP_6_LABELLED_CELL_INSPECTION('Property','Value',...) creates a new STEP_6_LABELLED_CELL_INSPECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before STEP_6_Labelled_Cell_Inspection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to STEP_6_Labelled_Cell_Inspection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help STEP_6_Labelled_Cell_Inspection

% Last Modified by GUIDE v2.5 11-Feb-2018 19:09:03

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @STEP_6_Labelled_Cell_Inspection_OpeningFcn, ...
    'gui_OutputFcn',  @STEP_6_Labelled_Cell_Inspection_OutputFcn, ...
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

% --- Executes just before STEP_6_Labelled_Cell_Inspection is made visible.
function STEP_6_Labelled_Cell_Inspection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to STEP_6_Labelled_Cell_Inspection (see VARARGIN)

% Choose default command line output for STEP_6_Labelled_Cell_Inspection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using STEP_6_Labelled_Cell_Inspection.
STEP_0_Parameters
handles.Color_channel_menu.String = Name_Channels(Color_Channel_Interest);

if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end


% UIWAIT makes STEP_6_Labelled_Cell_Inspection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = STEP_6_Labelled_Cell_Inspection_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in Update.
function Update_Callback(hObject, eventdata, handles)
% hObject    handle to Update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

current_idx=evalin('base','current_idx');
save_flag=evalin('base','save_flag');

if current_idx==save_flag
    save_yn='Yes';
else
    save_yn = questdlg('Do you wish to continue without saving your work on current image?','Continue?','Yes','No','No');
end

switch save_yn
    case 'Yes'
        
        axes(handles.Image_Current);
        cla;
        
        contents = get(handles.list_imgs,'String');
        img2update = contents{get(handles.list_imgs,'Value')};
        img2update_tif=img2update;
        img2update = strcat(img2update,'.fig');
        
        od=open(img2update);
        img_current=getimage(od);
        
        fig_title=get(gca);
        fig_title=fig_title.Title.String;
        
        axesObjs = get(od, 'Children');  %axes handles
        dataObjs = get(axesObjs, 'Children');
        objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
        
        try
            xdata_scatter = get(dataObjs, 'XData');  xdata_scatter=xdata_scatter{1};
            ydata_scatter = get(dataObjs, 'YData');  ydata_scatter=ydata_scatter{1};
        end
        
        close(od);
        
        hold on
        
        imshow(img_current, 'Parent', handles.Image_Current);
        
        try
            scatter(handles.Image_Current,xdata_scatter,ydata_scatter,25,'r','filled');
        end
        
        title(fig_title,'Interpreter', 'none')
        
        h_imcontrast=imcontrast(handles.Image_Current);
        
        list_of_imgs=evalin('base','list_of_imgs');
        current_idx=find(strcmp(list_of_imgs,img2update_tif));
        assignin('base','current_idx',current_idx);
        
    case 'No'
end

% --- Executes on selection change in list_imgs.
function list_imgs_Callback(hObject, eventdata, handles)
% hObject    handle to list_imgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_imgs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_imgs


% --- Executes during object creation, after setting all properties.
function list_imgs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_imgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'empty, click the INITIALIZE button!'});


% --- Executes on button press in Next_Image.
function Next_Image_Callback(hObject, eventdata, handles)
% hObject    handle to Next_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

current_idx=evalin('base','current_idx');
save_flag=evalin('base','save_flag');

if current_idx==save_flag
    save_yn='Yes';
else
    save_yn = questdlg('Do you wish to continue without saving your work on current image?','Continue?','Yes','No','No');
end

switch save_yn
    case 'Yes'
        list_of_imgs=evalin('base','list_of_imgs');
        
        if current_idx<length(list_of_imgs)
            
            current_idx=current_idx+1;
            axes(handles.Image_Current);
            cla;
            img2update = strcat(list_of_imgs{current_idx},'.fig');
            
            od=open(img2update);
            img_current=getimage(od);
            
            fig_title=get(gca);
            fig_title=fig_title.Title.String;
            
            axesObjs = get(od, 'Children');  %axes handles
            dataObjs = get(axesObjs, 'Children');
            objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
            xdata_scatter = get(dataObjs, 'XData');
            ydata_scatter = get(dataObjs, 'YData');
            
            close(od);
            
            imshow(img_current, 'Parent', handles.Image_Current);
            title(fig_title,'Interpreter', 'none')
            hold on
            
            try
                xdata_scatter=xdata_scatter{1};
                ydata_scatter=ydata_scatter{1};
                scatter(handles.Image_Current,xdata_scatter,ydata_scatter,25,'r','filled');
            end
            
            h_imcontrast=imcontrast(handles.Image_Current);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            assignin('base','current_idx',current_idx);
            
        else
            warndlg('There is no image next to this one','Last Image');
        end
        
    case 'No'
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --- Executes on button press in Prev_Image.
function Prev_Image_Callback(hObject, eventdata, handles)
% hObject    handle to Prev_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current_idx=evalin('base','current_idx');
save_flag=evalin('base','save_flag');

if current_idx==save_flag
    save_yn='Yes';
else
    save_yn = questdlg('Do you wish to continue without saving your work on current image?','Continue?','Yes','No','No');
end

switch save_yn
    case 'Yes'
        
        list_of_imgs=evalin('base','list_of_imgs');
        
        if current_idx>1
            
            current_idx=current_idx-1;
            
            axes(handles.Image_Current);
            cla;
            img2update = strcat(list_of_imgs{current_idx},'.fig');
            
            od=open(img2update);
            img_current=getimage(od);
            
            fig_title=get(gca);
            fig_title=fig_title.Title.String;
            
            axesObjs = get(od, 'Children');  %axes handles
            dataObjs = get(axesObjs, 'Children');
            objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
            xdata_scatter = get(dataObjs, 'XData');  xdata_scatter=xdata_scatter{1};
            ydata_scatter = get(dataObjs, 'YData');  ydata_scatter=ydata_scatter{1};
            
            close(od);
            
            hold on
            
            imshow(img_current, 'Parent', handles.Image_Current);
            scatter(handles.Image_Current,xdata_scatter,ydata_scatter,25,'r','filled');
            title(fig_title,'Interpreter', 'none')
            
            h_imcontrast=imcontrast(handles.Image_Current);
            
            assignin('base','current_idx',current_idx);
        else
            warndlg('There is no image previous to this one','First Image');
        end
        
    case 'No'
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in Add_cell.
function Add_cell_Callback(hObject, eventdata, handles)
% hObject    handle to Add_cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

brush off
exit_tf=false;

while true
    if exit_tf==3
        break;
    end
    [x,y,exit_tf] = ginput(1);
    if exit_tf~=3
        scatter(handles.Image_Current,x,y,'r','filled');
    end
end
set( gcf, 'pointer','arrow' )



% --- Executes on button press in Remove_Cell.
function Remove_Cell_Callback(hObject, eventdata, handles)
% hObject    handle to Remove_Cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h= brush;
set(h,'Color',[0 1 1],'Enable','on');

% --- Executes on button press in SAVE_tempo.tm
function SAVE_tempo_Callback(hObject, eventdata, handles)
% hObject    handle to SAVE_tempo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Get Data from current handle %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Save flag %%%
STEP_0_Parameters;
current_idx=evalin('base','current_idx');
ch_interest_idx=evalin('base','ch_interest_idx');
ch_interest_idx=find(Color_Channel_Interest==ch_interest_idx);
assignin('base','save_flag',current_idx)
%%%
brush off
od=(handles.Image_Current);
img_current=getimage(od);
axesObjs = get(od, 'Children');  %axes handles
xdata_scatter_cat=[]; ydata_scatter_cat=[];

for axesObjs_ii=1:length(axesObjs)-1
    xdata_scatter = axesObjs(axesObjs_ii).XData;
    ydata_scatter = axesObjs(axesObjs_ii).YData;
    
    xdata_scatter_cat=[xdata_scatter_cat xdata_scatter];
    ydata_scatter_cat=[ydata_scatter_cat ydata_scatter];
    
end

cell_pos=[ xdata_scatter_cat;  ydata_scatter_cat]';
nan_idx=isnan(ydata_scatter_cat);
cell_pos(nan_idx,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% save_tempo Cell Position %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STEP_0_Parameters;
cd Image_Analysed_ROI_absent

img_essence=evalin('base','img_essence');

% if Slice_AP_orPA==1
current_idx_mod=current_idx;
% elseif Slice_AP_orPA==-1
%     current_idx_mod=length(img_essence)-current_idx+1;
% end

img_essence(current_idx_mod).Color_Cells(ch_interest_idx).cell_locations=cell_pos;

assignin('base','img_essence',img_essence)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% save_tempo Figure %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_of_imgs=evalin('base','list_of_imgs');

fig_name_save=strcat(list_of_imgs{current_idx},'.fig');  %% name when saving

image_reassessed=figure;
imshow(img_current); hold on;

if ~isempty(cell_pos)
    scatter(cell_pos(:,1),cell_pos(:,2),16,'r','filled');
end

saveas(image_reassessed,fig_name_save)

close(image_reassessed);



cd ..
cd Image_Analysed

fig_name_roi=strcat(list_of_imgs{current_idx},'.fig');
        
od=open(fig_name_roi);
img_current_roi=getimage(od);
close(od);

image_reassessed_roi=figure;
imshow(img_current_roi); hold on;

if ~isempty(cell_pos)
    scatter(cell_pos(:,1),cell_pos(:,2),16,'r','filled');
end

saveas(image_reassessed_roi,fig_name_roi)

close(image_reassessed_roi);

cd ..
cd Image_Analysed_ROI_absent
        





% --- Executes on button press in Initialize.
function Initialize_Callback(hObject, eventdata, handles)
% hObject    handle to Initialize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

patient=warndlg('Give me a moment.. Loading Data','Be Patient');

STEP_0_Parameters;
load('Step_4_Angle_Search_Result.mat');
load('Step_5_Cell_Detection_Result.mat');

assignin('base','img_essence',img_essence)

parent_folder=pwd;

%%%%% which color channel to analyze %%%%%
STEP_0_Parameters;
color_ch_contents = get(handles.Color_channel_menu,'String')
color_ch_interest = color_ch_contents{get(handles.Color_channel_menu,'Value')}
ch_interest_idx=strfind(Name_Channels, color_ch_interest)
ch_interest_idx = find(not(cellfun('isempty', ch_interest_idx)))
assignin('base','ch_interest_idx',ch_interest_idx)

%%%%%%%%%%%%%%%%
cd Image_Analysed_ROI_absent

img_idx=anc_img_IDs;
ap_found=max_APpos_stage_final;
img_AP=[];
for img_ID=1:length(img_idx)-1
    img_AP=[img_AP, linspace(ap_found(img_ID),ap_found(img_ID+1),img_idx(img_ID+1)-img_idx(img_ID)+1)];
end

img_AP=round(unique(img_AP));
if Slice_AP_orPA==1
    img_idx=min(img_idx):max(img_idx);
else
    img_idx=max(img_idx):-1:min(img_idx);
end
%%%%%%%%%%%%%

list_of_imgs={img_name{img_idx, ch_interest_idx}};
handles.list_imgs.String = list_of_imgs;

close(patient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Load First Image after Initiation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axes(handles.Image_Current);
cla;

img2update = list_of_imgs{1};
img2update = strcat(img2update,'.fig');

od=open(img2update);
img_current=getimage(od);

fig_title=get(gca);
fig_title=fig_title.Title.String;

axesObjs = get(od, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children');
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
xdata_scatter = get(dataObjs, 'XData');
ydata_scatter = get(dataObjs, 'YData');

close(od);

imshow(img_current, 'Parent', handles.Image_Current);
title(fig_title,'Interpreter', 'none')
hold on
try
    xdata_scatter=xdata_scatter{1};
    ydata_scatter=ydata_scatter{1};
    
    scatter(handles.Image_Current,xdata_scatter,ydata_scatter,25,'r','filled');
end
current_idx=1;

assignin('base','parent_folder',parent_folder);
assignin('base','list_of_imgs',list_of_imgs);
assignin('base','current_idx',current_idx);
assignin('base','save_flag',0);

h_imcontrast=imcontrast(handles.Image_Current);



% --------------------------------------------------------------------
function brush_cell_remove_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to brush_cell_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h= brush;
set(h,'Color',[0 1 1],'Enable','on');


% --- Executes on button press in save_all.
function save_all_Callback(hObject, eventdata, handles)
% hObject    handle to save_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

save_yn = questdlg('Are you sure to save the current data and overwrite previously saved one?','Save?','Yes','No','No');

switch save_yn
    case 'Yes'
        patient=warndlg('Give me a moment.. Saving Data','Saving...');
        
        img_essence=evalin('base','img_essence');
        parent_folder=evalin('base','parent_folder');
        
        cd(parent_folder)
        save('Step_5_Cell_Detection_Result','img_essence','-v7.3');
        
        cd Image_Analysed_ROI_absent
        
        close(patient)
    case 'No'

end
