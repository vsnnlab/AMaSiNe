function varargout = STEP_8_Data_Visualisation(varargin)
% STEP_8_DATA_VISUALISATION MATLAB code for STEP_8_Data_Visualisation.fig
%      STEP_8_DATA_VISUALISATION, by itself, creates a new STEP_8_DATA_VISUALISATION or raises the existing
%      singleton*.
%
%      H = STEP_8_DATA_VISUALISATION returns the handle to a new STEP_8_DATA_VISUALISATION or the handle to
%      the existing singleton*.
%
%      STEP_8_DATA_VISUALISATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STEP_8_DATA_VISUALISATION.M with the given input arguments.
%
%      STEP_8_DATA_VISUALISATION('Property','Value',...) creates a new STEP_8_DATA_VISUALISATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before STEP_8_Data_Visualisation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to STEP_8_Data_Visualisation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help STEP_8_Data_Visualisation

% Last Modified by GUIDE v2.5 11-Feb-2018 19:09:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @STEP_8_Data_Visualisation_OpeningFcn, ...
    'gui_OutputFcn',  @STEP_8_Data_Visualisation_OutputFcn, ...
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


% --- Executes just before STEP_8_Data_Visualisation is made visible.
function STEP_8_Data_Visualisation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to STEP_8_Data_Visualisation (see VARARGIN)

% Choose default command line output for STEP_8_Data_Visualisation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

main_folder_dir=pwd;
subfolders=dir(main_folder_dir);
addpath(genpath(main_folder_dir));
cd(main_folder_dir);

load Step5_ANO_info.mat;
load shp_1.mat;
assignin('base','shp_1',shp_1);

Brain_Files=dir(fullfile(pwd,'*.mat'));
Brain_Files={Brain_Files.name};

set(handles.Brain_Pool_list, 'String', Brain_Files);
set(handles.Brain_Pool_list,'Max',10,'Min',1);
set(handles.Brain_Selected_List,'Max',10,'Min',1);

region_name_display={};

for ii=1:length(region_name_list)
    name_str_tempo=strcat({blanks(3*generation(ii))},{region_name_list(ii).name});
    region_name_display{ii}=name_str_tempo{1};
end

assignin('base','region_name_display',region_name_display)

set(handles.ROI_Pool_List, 'String', region_name_display);
set(handles.ROI_Pool_List,'Max',inf,'Min',1);
set(handles.ROI_Selected_List,'Max',inf,'Min',1);


% UIWAIT makes STEP_8_Data_Visualisation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = STEP_8_Data_Visualisation_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Brain_Pool_list.
function Brain_Pool_list_Callback(hObject, eventdata, handles)
% hObject    handle to Brain_Pool_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Brain_Pool_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Brain_Pool_list
% to string

% --- Executes during object creation, after setting all properties.
function Brain_Pool_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Brain_Pool_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ROI_Pool_List.
function ROI_Pool_List_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_Pool_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ROI_Pool_List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROI_Pool_List


% --- Executes during object creation, after setting all properties.
function ROI_Pool_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROI_Pool_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ROI_Selected_List.
function ROI_Selected_List_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_Selected_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ROI_Selected_List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROI_Selected_List


% --- Executes during object creation, after setting all properties.
function ROI_Selected_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROI_Selected_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Remove_ROI.
function Remove_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to Remove_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents_existing = cellstr(get(handles.ROI_Selected_List,'String'));

contents_removal = {contents_existing{get(handles.ROI_Selected_List,'Value')}};

contents_removal_names= intersect(contents_removal,contents_existing, 'stable');

Index = (contains(contents_existing,contents_removal_names));
contents_existing=contents_existing(~Index);

set(handles.ROI_Selected_List, 'String', contents_existing);



% --- Executes on selection change in Brain_Selected_List.
function Brain_Selected_List_Callback(hObject, eventdata, handles)
% hObject    handle to Brain_Selected_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Brain_Selected_List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Brain_Selected_List


% --- Executes during object creation, after setting all properties.
function Brain_Selected_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Brain_Selected_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Add_ROI.
function Add_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to Add_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents_existing = cellstr(get(handles.ROI_Selected_List,'String'));

contents_update = cellstr(get(handles.ROI_Pool_List,'String'));
contents_update = {contents_update{get(handles.ROI_Pool_List,'Value')}};

[~, contents_novel_idx]= setdiff(contents_update,contents_existing, 'stable');

for ii=1:length(contents_novel_idx)
    contents_existing{length(contents_existing)+1}=contents_update{contents_novel_idx(ii)};
end
contents_existing=contents_existing';

Index = (contains(contents_existing,'Currently Empty'));
contents_existing=contents_existing(~Index);

set(handles.ROI_Selected_List, 'String', contents_existing);

% --- Executes on button press in Add_brain.
function Add_brain_Callback(hObject, eventdata, handles)
% hObject    handle to Add_brain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents_existing = cellstr(get(handles.Brain_Selected_List,'String'));

contents_update = cellstr(get(handles.Brain_Pool_list,'String'));
contents_update = {contents_update{get(handles.Brain_Pool_list,'Value')}};

[~, contents_novel_idx]= setdiff(contents_update,contents_existing, 'stable');

for ii=1:length(contents_novel_idx)
    contents_existing{length(contents_existing)+1}=contents_update{contents_novel_idx(ii)};
end
contents_existing=contents_existing';

Index = (contains(contents_existing,'Currently Empty'));
contents_existing=contents_existing(~Index);

set(handles.Brain_Selected_List, 'String', contents_existing);



% --- Executes on button press in Remove_brain.
function Remove_brain_Callback(hObject, eventdata, handles)
% hObject    handle to Remove_brain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents_existing = cellstr(get(handles.Brain_Selected_List,'String'));

contents_removal = {contents_existing{get(handles.Brain_Selected_List,'Value')}};

contents_removal_names= intersect(contents_removal,contents_existing, 'stable');

Index = (contains(contents_existing,contents_removal_names));
contents_existing=contents_existing(~Index);

set(handles.Brain_Selected_List, 'String', contents_existing);


% --- Executes on button press in Visualize.
function Visualize_Callback(hObject, eventdata, handles)
% hObject    handle to Visualize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_folder_dir=pwd;
STEP_0_Parameters;
subfolders=dir(main_folder_dir);
addpath(genpath(main_folder_dir));
cd(main_folder_dir);

load Step5_ANO_info.mat;

region_name_display=evalin('base','region_name_display');

figure; hold on;

shp_1=evalin('base','shp_1');
plot(shp_1.RL(1).shp,'EdgeColor','none','FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.05);
        %%%%%% YH %%%%%%
        %%%% If you want to see your ROI only --> comment out line 302
        %%%% ('plot(shp_1.....')

brains_interest = cellstr(get(handles.Brain_Selected_List,'String'));
ROI_interest = cellstr(get(handles.ROI_Selected_List,'String'));
ROI_idx=find(endsWith(region_name_display,ROI_interest)); %%%%%

ROI_color_hex=region_color(ROI_idx);

ROI_color_list=nan(length(ROI_color_hex),3);

for ROI_color_idx=1:length(ROI_color_hex)
    ROI_hex_tempo=region_color{ROI_color_idx};
    r=hex2dec(ROI_hex_tempo(1:2)); g=hex2dec(ROI_hex_tempo(3:4)); b=hex2dec(ROI_hex_tempo(5:6));
    ROI_rgb_tempo=[r g b]/255;
    ROI_color_list(ROI_color_idx,:)=ROI_rgb_tempo;
end

Mouse_Color={'r','g','b','c','m','y','k'};

ROI_names=region_name_display(ROI_idx);
ROI_names=ROI_names';

assignin('base','brains_interest',brains_interest);
assignin('base','ROI_names',ROI_names);

for ii=1:length(brains_interest)
    brain_name=brains_interest{ii};
    load(brain_name);
 
    for jj=1:length(ROI_idx)
        
        ROI_file_name=strcat('shp_',num2str(ROI_idx(jj)));
        
        if ii==1
            if ~strcmp(ROI_file_name,'shp_1')
                load(strcat(ROI_file_name,'.mat'));
            end
            shp_tempo=eval(ROI_file_name);
            
            for rl=1:length(shp_tempo.RL)
                plot(shp_tempo.RL(rl).shp,'FaceColor',ROI_color_list(jj,:),'FaceAlpha',0.2,'EdgeColor','none');
            end
        end
        
        
        for color_ch_idx=1:length(Color_Channel_Interest)
            scatter3(cells_ROI(ROI_idx(jj)).Color_ch(color_ch_idx).RL(1).cells_pos(:,1),...
                cells_ROI(ROI_idx(jj)).Color_ch(color_ch_idx).RL(1).cells_pos(:,2), ...
                cells_ROI(ROI_idx(jj)).Color_ch(color_ch_idx).RL(1).cells_pos(:,3),...
                16,Mouse_Color{color_ch_idx},'filled');
            %
            scatter3(cells_ROI(ROI_idx(jj)).Color_ch(color_ch_idx).RL(2).cells_pos(:,1),...
                cells_ROI(ROI_idx(jj)).Color_ch(color_ch_idx).RL(2).cells_pos(:,2), ...
                cells_ROI(ROI_idx(jj)).Color_ch(color_ch_idx).RL(2).cells_pos(:,3),...
                16,Mouse_Color{color_ch_idx},'filled');
            
            no_cells{jj,2+3*(ii-1)*length(Color_Channel_Interest)+(color_ch_idx-1)*3}=size(cells_ROI(ROI_idx(jj)).Color_ch(color_ch_idx).RL(1).cells_pos,1);
            no_cells{jj,3+3*(ii-1)*length(Color_Channel_Interest)+(color_ch_idx-1)*3}=size(cells_ROI(ROI_idx(jj)).Color_ch(color_ch_idx).RL(2).cells_pos,1);
        end
               
    end
      
    camlight right
    assignin('base','no_cells',no_cells);

    
    
end


% --- Executes on button press in save_xls.
function save_xls_Callback(hObject, eventdata, handles)
% hObject    handle to save_xls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
STEP_0_Parameters;

file_name=get(handles.xls_name,'String');
file_name=file_name{1};

brains_interest=evalin('base','brains_interest');
no_cells=evalin('base','no_cells');
ROI_names=evalin('base','ROI_names');
no_channels=length(Color_Channel_Interest);


for ii=1:length(brains_interest)
    Mouse_name{2+3*(ii-1)*no_channels}=brains_interest{ii};
    for jj=1:no_channels   %%% temporary for beta (multi_channel)

        Ch_name{2+3*(ii-1)*(no_channels)+3*(jj-1)}=Name_Channels{Color_Channel_Interest(jj)};
        RL{3+3*(ii-1)*(no_channels)+3*(jj-1)}='Right hemi (cells)';
        RL{2+3*(ii-1)*(no_channels)+3*(jj-1)}='Left hemi (cells)';
    end
end

Mouse_name{end+length(RL)-length(Mouse_name)}=[];
Ch_name{end+1}=[];

column_name_tempo_1=vertcat(Mouse_name,Ch_name,RL);
column_name_tempo_2={'ROI name'; []; []};

column_name=horzcat(column_name_tempo_2,column_name_tempo_1);
main_info=horzcat(ROI_names,no_cells);
to_xls=vertcat(column_name,main_info);

% % % % save('STEP_8_Cells_in_ROI.mat', 'to_xls'); 
xlswrite(file_name,to_xls);

% --- Executes on button press in inj_vol_check.
function inj_vol_check_Callback(hObject, eventdata, handles)
% hObject    handle to inj_vol_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of inj_vol_check



function xls_name_Callback(hObject, eventdata, handles)
% hObject    handle to xls_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of xls_name as text
%        str2double(get(hObject,'String')) returns contents of xls_name as a double


% --- Executes during object creation, after setting all properties.
function xls_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xls_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
