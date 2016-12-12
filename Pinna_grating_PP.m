function varargout = Pinna_grating_PP(varargin)
% PINNA_GRATING_PP MATLAB code for Pinna_grating_PP.fig
%      PINNA_GRATING_PP, by itself, creates a new PINNA_GRATING_PP or raises the existing
%      singleton*.
%
%      H = PINNA_GRATING_PP returns the handle to a new PINNA_GRATING_PP or the handle to
%      the existing singleton*.
%
%      PINNA_GRATING_PP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PINNA_GRATING_PP.M with the given input arguments.
%
%      PINNA_GRATING_PP('Property','Value',...) creates a new PINNA_GRATING_PP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Pinna_grating_PP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Pinna_grating_PP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Pinna_grating_PP

% Last Modified by GUIDE v2.5 12-Dec-2016 12:35:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Pinna_grating_PP_OpeningFcn, ...
                   'gui_OutputFcn',  @Pinna_grating_PP_OutputFcn, ...
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


% --- Executes just before Pinna_grating_PP is made visible.
function Pinna_grating_PP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Pinna_grating_PP (see VARARGIN)

% Choose default command line output for Pinna_grating_PP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Pinna_grating_PP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Pinna_grating_PP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in select_directory.
function select_directory_Callback(hObject, eventdata, handles)
% hObject    handle to select_directory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ResultDirName = uigetdir();
set(handles.file_path,'String',ResultDirName);


function file_path_Callback(hObject, eventdata, handles)
% hObject    handle to file_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of file_path as text
%        str2double(get(hObject,'String')) returns contents of file_path as a double


% --- Executes during object creation, after setting all properties.
function file_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function one_trials_Callback(hObject, eventdata, handles)
% hObject    handle to one_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of one_trials as text
%        str2double(get(hObject,'String')) returns contents of one_trials as a double


% --- Executes during object creation, after setting all properties.
function one_trials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to one_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ResultDir = get(handles.file_path,'String');
one_trials = str2double(get(handles.one_trials,'String'));
duration = str2double(get(handles.duration,'String'));
match_time = str2double(get(handles.match_time,'String'));
angle_pattern = str2num(get(handles.angles_pattern,'String'));
move_speed_i = str2num(get(handles.move_speed,'String'));
angle_speed_i = str2num(get(handles.angle_speed,'String'));
is_binary_mask = get(handles.is_binary_mask,'Value');
mask_diameter= str2num(get(handles.mask_diameter,'String'));
mask_xpos = str2num(get(handles.mask_xpos,'String'));
mask_ypos = str2num(get(handles.mask_ypos,'String'));
use_eyelink = get(handles.use_eyelink,'Value');
calib_file = get(handles.gamma_path,'String');
fix_radius = str2num(get(handles.fix_radius,'String'));
Pinna_grating_main(angle_pattern,move_speed_i,angle_speed_i,...
	ResultDir,one_trials,duration,match_time,is_binary_mask,...
	mask_diameter,mask_xpos,mask_ypos,use_eyelink,fix_radius,...
	calib_file);

% --- Executes on button press in analysis.
function analysis_Callback(hObject, eventdata, handles)
% hObject    handle to analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ResultDir = get(handles.file_path,'String');
CommonResultDir = get(handles.common_path,'String');
cwd=pwd;
ana_angle_tuning(ResultDir,CommonResultDir);
cd(cwd);

% --- Executes on button press in demo.
function demo_Callback(hObject, eventdata, handles)
% hObject    handle to demo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
is_with_line  = get(handles.is_with_line,'Value');
benchmark  = get(handles.benchmark_demo,'Value');
Pinna_grating_withlineO(is_with_line,benchmark);

% --- Executes on button press in is_with_line.
function is_with_line_Callback(hObject, eventdata, handles)
% hObject    handle to is_with_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of is_with_line



function duration_Callback(hObject, eventdata, handles)
% hObject    handle to duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of duration as text
%        str2double(get(hObject,'String')) returns contents of duration as a double


% --- Executes during object creation, after setting all properties.
function duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function match_time_Callback(hObject, eventdata, handles)
% hObject    handle to match_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of match_time as text
%        str2double(get(hObject,'String')) returns contents of match_time as a double


% --- Executes during object creation, after setting all properties.
function match_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to match_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select_common_path.
function select_common_path_Callback(hObject, eventdata, handles)
% hObject    handle to select_common_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CommonDirName = uigetdir();
set(handles.common_path,'String',CommonDirName);

function common_path_Callback(hObject, eventdata, handles)
% hObject    handle to common_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of common_path as text
%        str2double(get(hObject,'String')) returns contents of common_path as a double


% --- Executes during object creation, after setting all properties.
function common_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to common_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in IsOnlyExp.
function IsOnlyExp_Callback(hObject, eventdata, handles)
% hObject    handle to IsOnlyExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IsOnlyExp


% --- Executes on button press in IsOnlyCon.
function IsOnlyCon_Callback(hObject, eventdata, handles)
% hObject    handle to IsOnlyCon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IsOnlyCon


% --- Executes on button press in IsOnlyCCW.
function IsOnlyCCW_Callback(hObject, eventdata, handles)
% hObject    handle to IsOnlyCCW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IsOnlyCCW


% --- Executes on button press in IsOnlyCW.
function IsOnlyCW_Callback(hObject, eventdata, handles)
% hObject    handle to IsOnlyCW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IsOnlyCW



function move_speed_Callback(hObject, eventdata, handles)
% hObject    handle to move_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of move_speed as text
%        str2double(get(hObject,'String')) returns contents of move_speed as a double


% --- Executes during object creation, after setting all properties.
function move_speed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to move_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function angle_speed_Callback(hObject, eventdata, handles)
% hObject    handle to angle_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_speed as text
%        str2double(get(hObject,'String')) returns contents of angle_speed as a double


% --- Executes during object creation, after setting all properties.
function angle_speed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function angles_pattern_Callback(hObject, eventdata, handles)
% hObject    handle to angles_pattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angles_pattern as text
%        str2double(get(hObject,'String')) returns contents of angles_pattern as a double


% --- Executes during object creation, after setting all properties.
function angles_pattern_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angles_pattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in is_binary_mask.
function is_binary_mask_Callback(hObject, eventdata, handles)
% hObject    handle to is_binary_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_binary_mask



function mask_diameter_Callback(hObject, eventdata, handles)
% hObject    handle to mask_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mask_diameter as text
%        str2double(get(hObject,'String')) returns contents of mask_diameter as a double


% --- Executes during object creation, after setting all properties.
function mask_diameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mask_xpos_Callback(hObject, eventdata, handles)
% hObject    handle to mask_xpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mask_xpos as text
%        str2double(get(hObject,'String')) returns contents of mask_xpos as a double


% --- Executes during object creation, after setting all properties.
function mask_xpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_xpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mask_ypos_Callback(hObject, eventdata, handles)
% hObject    handle to mask_ypos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mask_ypos as text
%        str2double(get(hObject,'String')) returns contents of mask_ypos as a double


% --- Executes during object creation, after setting all properties.
function mask_ypos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_ypos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select_gamma.
function select_gamma_Callback(hObject, eventdata, handles)
% hObject    handle to select_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f,p] = uigetfile();
set(handles.gamma_path,'String',[p f]);


function gamma_path_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma_path as text
%        str2double(get(hObject,'String')) returns contents of gamma_path as a double


% --- Executes during object creation, after setting all properties.
function gamma_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in proceduraldemo.
function proceduraldemo_Callback(hObject, eventdata, handles)
is_with_line  = get(handles.is_with_line,'Value');
benchmark  = get(handles.benchmark_demo,'Value');
Pinna_grating_withlineP(is_with_line,benchmark);

% --- Executes on button press in benchmark_demo.
function benchmark_demo_Callback(hObject, eventdata, handles)



% --- Executes on button press in use_eyelink.
function use_eyelink_Callback(hObject, eventdata, handles)
% hObject    handle to use_eyelink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_eyelink



function fix_radius_Callback(hObject, eventdata, handles)
% hObject    handle to fix_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fix_radius as text
%        str2double(get(hObject,'String')) returns contents of fix_radius as a double


% --- Executes during object creation, after setting all properties.
function fix_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fix_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
