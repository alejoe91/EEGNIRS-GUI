function varargout = CLASSIFICATION(varargin)
% CLASSIFICATION MATLAB code for CLASSIFICATION.fig
%      CLASSIFICATION, by itself, creates a new CLASSIFICATION or raises the existing
%      singleton*.
%
%      H = CLASSIFICATION returns the handle to a new CLASSIFICATION or the handle to
%      the existing singleton*.
%
%      CLASSIFICATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLASSIFICATION.M with the given input arguments.
%
%      CLASSIFICATION('Property','Value',...) creates a new CLASSIFICATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CLASSIFICATION_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CLASSIFICATION_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CLASSIFICATION

% Last Modified by GUIDE v2.5 30-Nov-2014 20:15:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CLASSIFICATION_OpeningFcn, ...
    'gui_OutputFcn',  @CLASSIFICATION_OutputFcn, ...
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


% --- Executes just before CLASSIFICATION is made visible.
function CLASSIFICATION_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CLASSIFICATION (see VARARGIN)

% Choose default command line output for CLASSIFICATION
handles.output = hObject;

if ~isempty(varargin)
    
    mainGuiInput = find(strcmp(varargin, 'Classification'));
    
    if ~(isempty(mainGuiInput)) || (length(varargin) < mainGuiInput) || (~ishandle(varargin{mainGuiInput-1}))
        % Remember the handle, and adjust our position
        handles.main = varargin{mainGuiInput-1};
        
        mainHandlesdata = guidata(handles.main);
        
        handles.datacells = mainHandlesdata.datacells;
        
        handles.preprocessedfile = mainHandlesdata.preprocessedfile;
        
        %eeg-nirs directory
        handles.dirEEGNIRS = mainHandlesdata.dir;
        
        %set handle variable to the current type of EEG (gnd,cmnref,slf) data selected
        handles.xxeeg = handles.datacells{1,2}{1,1};
        handles.eeglabels = handles.datacells{1,1}{1,1}{4,1};
        
        
        %if data passed by main gui, eeg already referenced to gnd
        set(handles.ref_check,'Value',1);
    end
    
end

%% Create Classification folder (if it doesn't exist already)

if ~exist([handles.dirEEGNIRS '\CLASSIFICATION'],'dir')
    mkdir(handles.dirEEGNIRS,'CLASSIFICATION');
end

handles.dir = [handles.dirEEGNIRS '\CLASSIFICATION'];


% declare features and trials handles
handles.features = cell(1,2);
handles.trials = [];
handles.string_trial = [];
handles.current_string_trial = [];
handles.from_trial = 0;
handles.to_trial = 0;

%set default reference
set(handles.ref_list,'Value',1);

%declare current reference variable (initialized at GND)
handles.current_ref = get(handles.ref_list,'Value');

%store all the different referenced signals, once computed
handles.xxeegref = [];
handles.xxeegsf = [];

%selected eeg signal (for features and view)
handles.xxeeg_current = handles.xxeeg;

%selected spatial configuration
handles.spatial_config = 0;

% set default tw, overlap and number of time windows
set(handles.tw_menu,'Value',2);
set(handles.over_menu,'Value',1);
set(handles.n_tw,'Value',1);
handles.tw = 1;
handles.overlap = 0;
handles.noftw = 1;

%set select view default strings

set(handles.feature_selection,'Value',1);
set(handles.signal_selection,'Value',1);

%time signals
set(handles.select_view,'Value',1);
handles.prep_select_view = get(handles.select_view,'String');

%features select view
str = cell(1,3);
str{1} = 'NONE';
str{2} = 'TOTAL POWER';
str{3} = 'DELTA';
str{4} = 'THETA';
str{5}= 'ALFA';
str{6} = 'BETA';
str{7} = 'GAMMA';
str{8} = 'GAMMA LOW';
str{9} = 'GAMMA HIGH';
str{10} = 'HbO AVERAGE';
str{11} = 'HbO SLOPE';
str{12} = 'HbR AVERAGE';
str{13} = 'HbR SLOPE';

handles.features_select_view = str;
handles.eegfeatures_select_view = str(1:9);
handles.nirsfeatures_select_view = [str{1} str(10:end)];

%set laplacian filter menu invisible
set(handles.choosefilt,'Visible','off');
set(handles.viewfilter,'Visible','off');
set(handles.newbutton,'Visible','off');

%set trial by trial button to 0 and relative components invisible
set(handles.trial_by_trial,'Value',0);
set(handles.selecttrialspanel,'Visible','off');
set(handles.trials_menu,'Visible','off');
set(handles.include_trial,'Visible','off');

%set view panel visibility off
set(handles.view_panel,'Visible','off');

%declare global variable to hide plots
handles.hideplot = [];
%1: one axis, 2:two axis
handles.typeplot = 0;

%save default view_panel dimension and position
handles.viewpanelpos = get(handles.view_panel,'Position');

%preprocesseddata
handles.preprocessdata = 1;
%featurecomputed
handles.featurecomputed = 0;


%CREATE Spatial filters folder and default filters .txt

if ~exist([handles.dir '\SpatialFilter Configuration'],'dir')
    mkdir(handles.dir,'SpatialFilter Configuration');
    handles.configdir = [handles.dir '\SpatialFilter Configuration'];
    
    %write default config (1,2,3)
    
    % Configuration_1
    
    fid = fopen([handles.configdir '\SpatialFilter_1.txt'],'wt');
    
    if fid~=-1
        
        fprintf(fid, '%d', 4);
        fprintf(fid,'\n');
        fprintf(fid, '%d', 9);
        fprintf(fid,'\n');
        
        fprintf(fid,'%s','Cz Fc1 Fc2 Cp1 Cp2');
        fprintf(fid,'\n');
        fprintf(fid,'%s','C3 Fc1 Fc5 Cp1 Cp5');
        fprintf(fid,'\n');
        fprintf(fid,'%s','C4 Fc2 Fc6 Cp3 Cp6');
        fprintf(fid,'\n');
        fprintf(fid,'%s','Fc1 F3 Fz C3 Cz');
        fprintf(fid,'\n');
        fprintf(fid,'%s','Fc2 F4 Fz C4 Cz');
        fprintf(fid,'\n');
        fprintf(fid,'%s','Cp1 C3 Cz P3 Pz');
        fprintf(fid,'\n');
        fprintf(fid,'%s','Cp2 C4 Cz P3 Pz');
        fprintf(fid,'\n');
        fprintf(fid,'%s','Fz Fc1 Fc2 F3 F4');
        fprintf(fid,'\n');
        fprintf(fid,'%s','Pz Cp1 Cp2 P3 P4');
        fprintf(fid,'\n');
        
        fprintf(fid, '%d', 3);
        fprintf(fid,'\n');
        fprintf(fid, '%d', 4);
        fprintf(fid,'\n');
        
        fprintf(fid,'%s','F3 Fc5 Fc1 Fz');
        fprintf(fid,'\n');
        fprintf(fid,'%s','F4 Fc6 Fc2 Fz');
        fprintf(fid,'\n');
        fprintf(fid,'%s','P3 Cp5 Cp1 Cz');
        fprintf(fid,'\n');
        fprintf(fid,'%s','P4 Cp6 Cp2 Pz');
        fprintf(fid,'\n');
        
        fprintf(fid, '%d', 2);
        fprintf(fid,'\n');
        fprintf(fid, '%d', 4);
        fprintf(fid,'\n');
        
        fprintf(fid,'%s','Fc5 F3 C3');
        fprintf(fid,'\n');
        fprintf(fid,'%s','Fc6 F4 C4');
        fprintf(fid,'\n');
        fprintf(fid,'%s','Cp5 C3 P3');
        fprintf(fid,'\n');
        fprintf(fid,'%s','Cp6 C4 P4');
        
        str = { '1' };        
        
    end
    
else
    
    handles.configdir = [handles.dir '\SpatialFilter Configuration'];
    
    count = 0;
    str = [];
    
    while (exist([handles.configdir '\SpatialFilter_' num2str(count + 1) '.txt'],'file'))
        count = count + 1;
        str{count} = num2str(count);
    end
end

%set default configuration ID
set(handles.choosefilt,'String',str);
set(handles.choosefilt,'Value',1);

% Update handles structure
guidata(hObject, handles);




% --- Outputs from this function are returned to the command line.
function varargout = CLASSIFICATION_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in ref_list.
function ref_list_Callback(hObject, eventdata, handles)
% hObject    handle to ref_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ref_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ref_list

if get(hObject,'Value') == 3
    set(handles.choosefilt,'Visible','on');
    set(handles.viewfilter,'Visible','on');
    set(handles.newbutton,'Visible','on');
else
    set(handles.choosefilt,'Visible','off');
    set(handles.viewfilter,'Visible','off');
    set(handles.newbutton,'Visible','off');
end

set(handles.ref_check,'Value',0);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ref_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in go_ref_button.
function go_ref_button_Callback(hObject, eventdata, handles)
% hObject    handle to go_ref_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.current_ref ~= get(handles.ref_list,'Value')
    
    %xxeeg = handles.xxeeg;
    
    switch get(handles.ref_list,'Value')
        
        case 1
            
            handles.xxeeg_current = handles.xxeeg;
            
            handles.current_ref = get(handles.ref_list,'Value');
            set(handles.ref_check,'Value',1);
            
        case 2
            
            %if the common reference signal has not been computed yet
            if isempty(handles.xxeegref)
                
                sub_matrix = mean(handles.xxeeg,2);
                for i = 1:size(handles.xxeeg,2)
                    handles.xxeegref(:,i) = handles.xxeeg(:,i) - sub_matrix;
                end
                
            end
            
            handles.xxeeg_current = handles.xxeegref;
            
            handles.current_ref = get(handles.ref_list,'Value');
            set(handles.ref_check,'Value',1);
            
            
        case 3
            
            %if the spatial filtered signal with the selected configuration
            %has not been computed yet
            if isempty(handles.xxeegref) || handles.spatial_config ~= get(handles.choosefilt,'Value');
                
                filtertxt = [handles.configdir '\SpatialFilter_' num2str(get(handles.choosefilt,'Value')) '.txt'];
                [xxeegsf, filteredchan] = laplacian_sf (handles.xxeeg, handles.eeglabels, filtertxt);
                handles.xxeegsf = xxeegsf;
                
            end
            
            handles.xxeeg_current = handles.xxeegsf;
            
            handles.current_ref = get(handles.ref_list,'Value');
            set(handles.ref_check,'Value',1);
            
    end
    
    guidata(hObject, handles);
    
end


% --- Executes on button press in ref_check.
function ref_check_Callback(hObject, eventdata, handles)
% hObject    handle to ref_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ref_check


% --- Executes on selection change in choosefilt.
function choosefilt_Callback(hObject, eventdata, handles)
% hObject    handle to choosefilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns choosefilt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from choosefilt




% --- Executes during object creation, after setting all properties.
function choosefilt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to choosefilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in tw_menu.
function tw_menu_Callback(hObject, eventdata, handles)
% hObject    handle to tw_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tw_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tw_menu

val = get(hObject,'Value');

switch val
    case 1
        handles.tw = 0.5;
    case 2
        handles.tw = 1;
    case 3
        handles.tw = 1.5;
    case 4
        handles.tw = 2;
end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function tw_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tw_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in over_menu.
function over_menu_Callback(hObject, eventdata, handles)
% hObject    handle to over_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(hObject,'Value');

switch val
    case 1
        handles.overlap = 0;
    case 2
        handles.overlap = 25;
    case 3
        handles.overlap = 50;
    case 4
        handles.overlap = 75;
end

guidata(hObject,handles);

% Hints: contents = cellstr(get(hObject,'String')) returns over_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from over_menu


% --- Executes during object creation, after setting all properties.
function over_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to over_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in n_tw.
function n_tw_Callback(hObject, eventdata, handles)
% hObject    handle to n_tw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(hObject,'Value');

switch val
    case 1
        handles.noftw = 1;
    case 2
        handles.noftw = 2;
    case 3
        handles.noftw = 3;
    case 4
        handles.noftw = 4;
    case 5
        handles.noftw = 5;
end

guidata(hObject,handles);

% Hints: contents = cellstr(get(hObject,'String')) returns n_tw contents as cell array
%        contents{get(hObject,'Value')} returns selected item from n_tw


% --- Executes during object creation, after setting all properties.
function n_tw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_tw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in select_view.
function select_view_Callback(hObject, eventdata, handles)
% hObject    handle to select_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns select_view contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_view

val = get(hObject,'Value');
str = get(hObject,'String');

val_sig = get(handles.signal_selection,'Value');
val_feat = get(handles.feature_selection,'Value');

%remove events
set(handles.eventsbutton,'Value',0);

switch val_sig
    
    case 1
        
        switch str{val}
            
            case 'NONE'
                if handles.typeplot==1
                    plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                    set(plot_hide,'Visible','off');
                    legend('off');
                end
                if ~isempty(handles.hideplot) && handles.typeplot==2
                    axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                    set(axes_hide,'Visible','off');
                    legend('off');
                end
                if handles.typeplot==3
                    plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                    plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                    set(plot_hide1,'Visible','off');
                    set(plot_hide2,'Visible','off');
                    legend(handles.sub1,'off');
                    legend(handles.sub2,'off');
                end
                
                set(handles.channel1,'String', 'Channel1');
                set(handles.channel2,'String', 'Channel2');
                set(handles.channel1,'Value', 1);
                set(handles.channel2,'Value', 1);
                handles.typeplot=1;
                
            case 'EEG'
                if handles.preprocessdata
                    % update channel1 and channel2 popupmenu
                    eegchanstr = [];
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    eegchanstr{1} = 'NONE';
                    for i = 1:length(eeglabels)
                        
                        eegchanstr{i+1} = eeglabels{i};
                        
                    end
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel2,'String', 'Channel2');
                    set(handles.channel2,'Value', 1);
                    set(handles.channel2,'Visible', 'off');
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                end
                
            case 'HbO'
                if handles.preprocessdata
                    % update channel1 and channel2 popupmenu
                    nirschanstr = [];
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    nirschanstr{1} = 'NONE';
                    for i = 1:length(nirslabels)
                        
                        nirschanstr{i+1} = nirslabels{i};
                        
                    end
                    set(handles.channel1,'String', nirschanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel2,'String', 'Channel2');
                    set(handles.channel2,'Value', 1);
                    set(handles.channel2,'Visible', 'off');
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                end
                
            case 'HbR'
                if handles.preprocessdata
                    % update channel1 and channel2 popupmenu
                    nirschanstr = [];
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    nirschanstr{1} = 'NONE';
                    for i = 1:length(nirslabels)
                        
                        nirschanstr{i+1} = nirslabels{i};
                        
                    end
                    set(handles.channel1,'String', nirschanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel2,'String', 'Channel2');
                    set(handles.channel2,'Value', 1);
                    set(handles.channel2,'Visible', 'off');
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                end
                
            case 'HbT'
                if handles.preprocessdata
                    % update channel1 and channel2 popupmenu
                    nirschanstr = [];
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    nirschanstr{1} = 'NONE';
                    for i = 1:length(nirslabels)
                        
                        nirschanstr{i+1} = nirslabels{i};
                        
                    end
                    set(handles.channel1,'String', nirschanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel2,'String', 'Channel2');
                    set(handles.channel2,'Value', 1);
                    set(handles.channel2,'Visible', 'off');
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                end
                
            case 'EEG + HbO'
                if handles.preprocessdata
                    % update channel1 and channel2 popupmenu
                    eegchanstr = [];
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    eegchanstr{1} = 'NONE';
                    for i = 1:length(eeglabels)
                        
                        eegchanstr{i+1} = eeglabels{i};
                        
                    end
                    % update channel1 and channel2 popupmenu
                    nirschanstr = [];
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    nirschanstr{1} = 'NONE';
                    for i = 1:length(nirslabels)
                        
                        nirschanstr{i+1} = nirslabels{i};
                        
                    end
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel2,'String', nirschanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel2,'Value', 1);
                    set(handles.channel1,'Visible', 'on');
                    set(handles.channel2,'Visible', 'on');
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                end
                
            case 'EEG + HbR'
                if handles.preprocessdata
                    % update channel1 and channel2 popupmenu
                    eegchanstr = [];
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    eegchanstr{1} = 'NONE';
                    for i = 1:length(eeglabels)
                        
                        eegchanstr{i+1} = eeglabels{i};
                        
                    end
                    % update channel1 and channel2 popupmenu
                    nirschanstr = [];
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    nirschanstr{1} = 'NONE';
                    for i = 1:length(nirslabels)
                        
                        nirschanstr{i+1} = nirslabels{i};
                        
                    end
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel2,'String', nirschanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel2,'Value', 1);
                    set(handles.channel1,'Visible', 'on');
                    set(handles.channel2,'Visible', 'on');
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                end
                
            case 'HbO + HbR'
                if handles.preprocessdata
                    % update channel1 and channel2 popupmenu
                    nirschanstr = [];
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    nirschanstr{1} = 'NONE';
                    for i = 1:length(nirslabels)
                        
                        nirschanstr{i+1} = nirslabels{i};
                        
                    end
                    set(handles.channel1,'String', nirschanstr);
                    set(handles.channel2,'String', nirschanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel2,'Value', 1);
                    set(handles.channel1,'Visible', 'on');
                    set(handles.channel2,'Visible', 'on');
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                end
                
            case 'EEG + EEG'
                if handles.preprocessdata
                    % update channel1 and channel2 popupmenu
                    eegchanstr = [];
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    eegchanstr{1} = 'NONE';
                    for i = 1:length(eeglabels)
                        
                        eegchanstr{i+1} = eeglabels{i};
                        
                    end
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel2,'String', eegchanstr);
                    set(handles.channel2,'Value', 1);
                    set(handles.channel2,'Visible', 'on');
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                end
                
            case 'HbO + HbO'
                if handles.preprocessdata
                    % update channel1 and channel2 popupmenu
                    nirschanstr = [];
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    nirschanstr{1} = 'NONE';
                    for i = 1:length(nirslabels)
                        
                        nirschanstr{i+1} = nirslabels{i};
                        
                    end
                    set(handles.channel1,'String', nirschanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel2,'String', nirschanstr);
                    set(handles.channel2,'Value', 1);
                    set(handles.channel2,'Visible', 'on');
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                end
                
            case 'HbR + HbR'
                if handles.preprocessdata
                    % update channel1 and channel2 popupmenu
                    nirschanstr = [];
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    nirschanstr{1} = 'NONE';
                    for i = 1:length(nirslabels)
                        
                        nirschanstr{i+1} = nirslabels{i};
                        
                    end
                    set(handles.channel1,'String', nirschanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel2,'String', nirschanstr);
                    set(handles.channel2,'Value', 1);
                    set(handles.channel2,'Visible', 'on');
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                end
        end
        
    case 2
        
        switch str{val}
            
            case 'NONE'
                
                if handles.typeplot==1
                    plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                    set(plot_hide,'Visible','off');
                    legend('off');
                end
                if ~isempty(handles.hideplot) && handles.typeplot==2
                    axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                    set(axes_hide,'Visible','off');
                    legend('off');
                end
                if handles.typeplot==3
                    plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                    plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                    set(plot_hide1,'Visible','off');
                    set(plot_hide2,'Visible','off');
                    legend(handles.sub1,'off');
                    legend(handles.sub2,'off');
                end
                
                set(handles.channel1,'String', 'Channel1');
                set(handles.channel2,'Visible', 'off');
                set(handles.channel1,'Value', 1);
                handles.typeplot=1;
                handles.typeplot=1;
                
            case 'TOTAL POWER'
                if handles.featurecomputed
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                    % update channel1 and channel2 popupmenu
                    eegchanstr = [];
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    eegchanstr{1} = 'NONE';
                    for i = 1:length(eeglabels)
                        
                        eegchanstr{i+1} = eeglabels{i};
                        
                    end
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel2,'Visible', 'off');
                    set(handles.channel1,'Value', 1);
                    handles.typeplot=1;
                    
                end
                
            case 'DELTA'
                if handles.featurecomputed
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                    % update channel1 and channel2 popupmenu
                    eegchanstr = [];
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    eegchanstr{1} = 'NONE';
                    for i = 1:length(eeglabels)
                        
                        eegchanstr{i+1} = eeglabels{i};
                        
                    end
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel2,'Visible', 'off');
                    set(handles.channel1,'Value', 1);
                    handles.typeplot=1;
                    
                end
                
            case 'THETA'
                if handles.featurecomputed
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                    % update channel1 and channel2 popupmenu
                    eegchanstr = [];
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    eegchanstr{1} = 'NONE';
                    for i = 1:length(eeglabels)
                        
                        eegchanstr{i+1} = eeglabels{i};
                        
                    end
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel2,'Visible', 'off');
                    set(handles.channel1,'Value', 1);
                    handles.typeplot=1;
                    
                end
                
            case 'ALFA'
                if handles.featurecomputed
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                    % update channel1 and channel2 popupmenu
                    eegchanstr = [];
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    eegchanstr{1} = 'NONE';
                    for i = 1:length(eeglabels)
                        
                        eegchanstr{i+1} = eeglabels{i};
                        
                    end
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel2,'Visible', 'off');
                    set(handles.channel1,'Value', 1);
                    handles.typeplot=1;
                    
                end
                
            case 'BETA'
                if handles.featurecomputed
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                    % update channel1 and channel2 popupmenu
                    eegchanstr = [];
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    eegchanstr{1} = 'NONE';
                    for i = 1:length(eeglabels)
                        
                        eegchanstr{i+1} = eeglabels{i};
                        
                    end
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel2,'Visible', 'off');
                    set(handles.channel1,'Value', 1);
                    handles.typeplot=1;
                    
                end
                
            case 'GAMMA'
                if handles.featurecomputed
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                    % update channel1 and channel2 popupmenu
                    eegchanstr = [];
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    eegchanstr{1} = 'NONE';
                    for i = 1:length(eeglabels)
                        
                        eegchanstr{i+1} = eeglabels{i};
                        
                    end
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel2,'Visible', 'off');
                    set(handles.channel1,'Value', 1);
                    handles.typeplot=1;
                    
                end
                
            case 'GAMMA LOW'
                if handles.featurecomputed
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                    % update channel1 and channel2 popupmenu
                    eegchanstr = [];
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    eegchanstr{1} = 'NONE';
                    for i = 1:length(eeglabels)
                        
                        eegchanstr{i+1} = eeglabels{i};
                        
                    end
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel2,'Visible', 'off');
                    set(handles.channel1,'Value', 1);
                    handles.typeplot=1;
                    
                end
                
            case 'GAMMA HIGH'
                if handles.featurecomputed
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                    % update channel1 and channel2 popupmenu
                    eegchanstr = [];
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    eegchanstr{1} = 'NONE';
                    for i = 1:length(eeglabels)
                        
                        eegchanstr{i+1} = eeglabels{i};
                        
                    end
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel1,'Value', 1);
                    set(handles.channel1,'String', eegchanstr);
                    set(handles.channel2,'Visible', 'off');
                    set(handles.channel1,'Value', 1);
                    handles.typeplot=1;
                    
                end
                
            case 'HbO AVERAGE'
                if handles.featurecomputed
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                    % update channel1 and channel2 popupmenu
                    nirschanstr = [];
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    nirschanstr{1} = 'NONE';
                    for i = 1:length(nirslabels)
                        
                        nirschanstr{i+1} = nirslabels{i};
                        
                    end
                    set(handles.channel1,'String', nirschanstr);
                    set(handles.channel2,'Visible', 'off');
                    set(handles.channel1,'Value', 1);
                    handles.typeplot=1;
                    
                end
                
            case 'HbO SLOPE'
                if handles.featurecomputed
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                    % update channel1 and channel2 popupmenu
                    nirschanstr = [];
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    nirschanstr{1} = 'NONE';
                    for i = 1:length(nirslabels)
                        
                        nirschanstr{i+1} = nirslabels{i};
                        
                    end
                    set(handles.channel1,'String', nirschanstr);
                    set(handles.channel2,'Visible', 'off');
                    set(handles.channel1,'Value', 1);
                    handles.typeplot=1;
                    
                end
                
            case 'HbR AVERAGE'
                if handles.featurecomputed
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                    % update channel1 and channel2 popupmenu
                    nirschanstr = [];
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    nirschanstr{1} = 'NONE';
                    for i = 1:length(nirslabels)
                        
                        nirschanstr{i+1} = nirslabels{i};
                        
                    end
                    set(handles.channel1,'String', nirschanstr);
                    set(handles.channel2,'Visible', 'off');
                    set(handles.channel1,'Value', 1);
                    handles.typeplot=1;
                    
                end
                
            case 'HbR SLOPE'
                if handles.featurecomputed
                    
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend(handles.sub1,'off');
                        legend(handles.sub2,'off');
                    end
                    
                    % update channel1 and channel2 popupmenu
                    nirschanstr = [];
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    nirschanstr{1} = 'NONE';
                    for i = 1:length(nirslabels)
                        
                        nirschanstr{i+1} = nirslabels{i};
                        
                    end
                    set(handles.channel1,'String', nirschanstr);
                    set(handles.channel2,'Visible', 'off');
                    set(handles.channel1,'Value', 1);
                    handles.typeplot=1;
                    
                end
                
        end
        
        
        
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function select_view_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in channel1.
function channel1_Callback(hObject, eventdata, handles)
% hObject    handle to channel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channel1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channel1

val = get(handles.select_view,'Value');
str = get(handles.select_view,'String');

val_sig = get(handles.signal_selection,'Value');
val_feat = get(handles.feature_selection,'Value');

%remove events
%set(handles.eventsbutton,'Value',0);

if (~get(handles.trial_by_trial,'Value'))

switch val_sig
    
    case 1
        
        switch str{val}
            case 'NONE'
                if handles.typeplot==1
                    plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                    set(plot_hide,'Visible','off');
                    legend('off');
                end
                if ~isempty(handles.hideplot) && handles.typeplot==2
                    axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                    set(axes_hide,'Visible','off');
                    legend('off');
                end
                if handles.typeplot==3
                    plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                    plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                    set(plot_hide1,'Visible','off');
                    set(plot_hide2,'Visible','off');
                    legend(handles.sub1,'off');
                    legend(handles.sub2,'off');
                end
            case 'EEG'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                tteeg = handles.datacells{1,2}{1,2};
                xxeeg = handles.datacells{1,2}{1,1};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,tteeg,xxeeg(:,val1 - 1),'Color',[0 0.5 0]), xlabel('time [s]'), ylabel('EEG [uV]'), title(['EEG: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                
                guidata(hObject, handles);
            case 'HbO'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                ttnirs = handles.datacells{1,2}{1,5};
                xxnirs = handles.datacells{1,2}{1,3};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ttnirs,xxnirs(:,val1 - 1),'r','linewidth',2), xlabel('time [s]'), ylabel('HbO [mM/L]'),title(['HbO: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
            case 'HbR'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                ttnirs = handles.datacells{1,2}{1,5};
                xxnirsR = handles.datacells{1,2}{1,4};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ttnirs,xxnirsR(:,val1 - 1),'b','linewidth',2), xlabel('time [s]'), ylabel('HbR [mM/L]'),title(['HbR: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
            case 'HbT'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                ttnirs = handles.datacells{1,2}{1,5};
                xxnirs = handles.datacells{1,2}{1,3};
                xxnirsR = handles.datacells{1,2}{1,4};
                xxnirsT = xxnirs - xxnirsR;
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ttnirs,xxnirsT(:,val1 - 1),'g','linewidth',2), xlabel('time [s]'), ylabel('HbT [mM/L]'),title(['HbT: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                end
                guidata(hObject, handles);
            case 'EEG + HbO'
                if get(handles.channel2,'Value')~=1
                    val1 = get(handles.channel1,'Value');
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    tteeg = handles.datacells{1,2}{1,2};
                    xxeeg = handles.datacells{1,2}{1,1};
                    
                    val2 = get(handles.channel2,'Value');
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    ttnirs = handles.datacells{1,2}{1,5};
                    xxnirs = handles.datacells{1,2}{1,3};
                    set(handles.viewplot,'Visible','on');
                    
                    if val1 ~= 1 && val2 ~= 1
                        axes(handles.viewplot);
                        [hAx,hLine1,hLine2] = plotyy(handles.viewplot,tteeg,xxeeg(:,val1 - 1),ttnirs,xxnirs(:,val2 - 1),'plot');
                        set(hLine1,'color',[0 0.5 0]);
                        set(hLine2,'color','red');
                        set(hAx,{'ycolor'},{[0 0.5 0];'r'}) ;
                        ylabel(hAx(1),'EEG [uV]') % left y-axis
                        ylabel(hAx(2),'HbO [mM/L]') % right y-axis
                        legend('EEG','HbO');
                        xlabel('time [s]'), title(['EEG: ' eeglabels{val1-1} ' HbO: ' nirslabels{val2-1}]);
                        grid(handles.viewplot,'on');
                        handles.hideplot = hAx;
                        handles.typeplot=2;
                    end
                    guidata(hObject, handles);
                end
            case 'EEG + HbR'
                if get(handles.channel2,'Value')~=1
                    val1 = get(handles.channel1,'Value');
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    tteeg = handles.datacells{1,2}{1,2};
                    xxeeg = handles.datacells{1,2}{1,1};
                    
                    val2 = get(handles.channel2,'Value');
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    ttnirs = handles.datacells{1,2}{1,5};
                    xxnirsR = handles.datacells{1,2}{1,4};
                    set(handles.viewplot,'Visible','on');
                    
                    if val1 ~= 1 && val2 ~= 1
                        axes(handles.viewplot);
                        [hAx,hLine1,hLine2] = plotyy(handles.viewplot,tteeg,xxeeg(:,val1 - 1),ttnirs,xxnirsR(:,val2 - 1),'plot');
                        set(hLine1,'color',[0 0.5 0]);
                        set(hLine2,'color','blue');
                        set(hAx,{'ycolor'},{[0 0.5 0];'b'}) ;
                        ylabel(hAx(1),'EEG [uV]') % left y-axis
                        ylabel(hAx(2),'HbR [mM/L]') % right y-axis
                        legend('EEG','HbR');
                        xlabel('time [s]'), title(['EEG: ' eeglabels{val1-1} ' HbR: ' nirslabels{val2-1}]);
                        grid(handles.viewplot,'on');
                        
                        handles.hideplot = hAx;
                        handles.typeplot=2;
                    end
                    
                    guidata(hObject, handles);
                end
            case 'HbO + HbR'
                if get(handles.channel2,'Value')~=1
                    val1 = get(handles.channel1,'Value');
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    ttnirs = handles.datacells{1,2}{1,5};
                    xxnirs = handles.datacells{1,2}{1,3};
                    
                    val2 = get(handles.channel2,'Value');
                    xxnirsR = handles.datacells{1,2}{1,4};
                    set(handles.viewplot,'Visible','off');
                    
                    if val1 ~= 1 && val2 ~= 1
                        
                        
                        set(handles.sub1,'Visible','on');
                        axes(handles.sub1);
                        plot(ttnirs,xxnirsR(:,val1 - 1),'b');
                        hold on
                        plot(ttnirs,xxnirs(:,val1 - 1),'r');
                        hold off
                        
                        legend('HbO','HbR');
                        ylabel('Hb [mM/L]'), title(['HbO: ' nirslabels{val1-1} ' HbR: ' nirslabels{val1-1}]);
                        grid(handles.sub1,'on');
                        
                        set(handles.sub2,'Visible','on');
                        axes(handles.sub2)
                        plot(ttnirs,xxnirs(:,val2 - 1),'r');
                        hold on
                        plot(ttnirs,xxnirsR(:,val2 - 1),'b');
                        hold off
                        legend('HbO','HbR');
                        xlabel('time [s]'), ylabel('Hb [mM/L]'), title(['HbO: ' nirslabels{val2-1} ' HbR: ' nirslabels{val2-1}]);
                        grid(handles.sub2,'on');
                        
                        handles.typeplot = 3;
                    end
                    
                    guidata(hObject, handles);
                end
                
            case 'EEG + EEG'
                if get(handles.channel2,'Value')~=1
                    val1 = get(handles.channel1,'Value');
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    tteeg = handles.datacells{1,2}{1,2};
                    xxeeg = handles.datacells{1,2}{1,1};
                    
                    val2 = get(handles.channel2,'Value');
                    set(handles.viewplot,'Visible','on');
                    
                    if val1 ~= 1 && val2 ~= 1
                        axes(handles.viewplot);
                        plot(handles.viewplot,tteeg,xxeeg(:,val1 - 1),'Color',[0 0.5 0])
                        hold on
                        plot(handles.viewplot,tteeg,xxeeg(:,val2 - 1),'Color',[0 0.2 0]);
                        hold off
                        
                        legend('EEG-1','EEG-2');
                        xlabel('time [s]'), ylabel('EEG [uV]'), title(['EEG-1: ' eeglabels{val1-1} ' EEG-2: ' eeglabels{val2-1}]);
                        grid(handles.viewplot,'on');
                        
                        handles.typeplot=1;
                    end
                    guidata(hObject, handles);
                end
                
            case 'HbO + HbO'
                if get(handles.channel2,'Value')~=1
                    val1 = get(handles.channel1,'Value');
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    ttnirs = handles.datacells{1,2}{1,5};
                    xxnirs = handles.datacells{1,2}{1,3};
                    
                    val2 = get(handles.channel2,'Value');
                    set(handles.viewplot,'Visible','on');
                    
                    if val1 ~= 1 && val2 ~= 1
                        axes(handles.viewplot);
                        plot(handles.viewplot,ttnirs,xxnirs(:,val1 - 1),'r');
                        hold on
                        plot(ttnirs,xxnirs(:,val2 - 1),'Color',[0.5,0,0.5]);
                        hold off
                        
                        legend('HbO-1','HbO-2');
                        xlabel('time [s]'), ylabel('HbO [mM/L]'), title(['HbO-1: ' nirslabels{val1-1} ' HbO-2: ' nirslabels{val2-1}]);
                        grid(handles.viewplot,'on');
                        
                        handles.typeplot=1;
                    end
                    guidata(hObject, handles);
                end
                
            case 'HbR + HbR'
                if get(handles.channel2,'Value')~=1
                    val1 = get(handles.channel1,'Value');
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    ttnirs = handles.datacells{1,2}{1,5};
                    xxnirsR = handles.datacells{1,2}{1,4};
                    
                    val2 = get(handles.channel2,'Value');
                    set(handles.viewplot,'Visible','on');
                    
                    if val1 ~= 1 && val2 ~= 1
                        axes(handles.viewplot);
                        plot(handles.viewplot,ttnirs,xxnirsR(:,val1 - 1),'b');
                        hold on
                        plot(ttnirs,xxnirsR(:,val2 - 1),'Color',[0,0.5,0.5]);
                        hold off
                        
                        legend('HbR-1','HbR-2');
                        xlabel('time [s]'), ylabel('HbR [mM/L]'), title(['HbR-1: ' nirslabels{val1-1} ' HbR-2: ' nirslabels{val2-1}]);
                        grid(handles.viewplot,'on');
                        
                        handles.typeplot=1;
                    end
                    guidata(hObject, handles);
                end
                
                guidata(hObject, handles);
                
        end
        
    case 2
        
        switch str{val}
            case 'NONE'
                if handles.typeplot==1
                    plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                    set(plot_hide,'Visible','off');
                    legend('off');
                end
                if ~isempty(handles.hideplot) && handles.typeplot==2
                    axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                    set(axes_hide,'Visible','off');
                    legend('off');
                end
                if handles.typeplot==3
                    plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                    plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                    set(plot_hide1,'Visible','off');
                    set(plot_hide2,'Visible','off');
                    legend(handles.sub1,'off');
                    legend(handles.sub2,'off');
                end
            case 'TOTAL POWER'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                ts = handles.features{1,2}{1,1};
                P = handles.features{1,2}{1,3};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts,P(:,val1 - 1,1),'Color',[0.6 0.6 0.6],'linewidth',2), xlabel('time [s]'), ylabel('Total Power [mW]^2/Hz'), title(['Feature: TOTAL POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
            case 'DELTA'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                ts = handles.features{1,2}{1,1};
                P = handles.features{1,2}{1,3};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts,P(:,val1 - 1,2),'Color',[0.5 0.5 0.8],'linewidth',2), xlabel('time [s]'), ylabel('Delta Power [mW]^2/Hz'), title(['Feature: Delta POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'THETA'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                ts = handles.features{1,2}{1,1};
                P = handles.features{1,2}{1,3};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts,P(:,val1 - 1,3),'Color',[0.5 0.8 0.8],'linewidth',2), xlabel('time [s]'), ylabel('Theta Power [mW]^2/Hz'), title(['Feature: Theta POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'ALFA'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                ts = handles.features{1,2}{1,1};
                P = handles.features{1,2}{1,3};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts,P(:,val1 - 1,4),'Color',[0.2 0.8 0.2],'linewidth',2), xlabel('time [s]'), ylabel('Alfa Power [mW]^2/Hz'), title(['Feature: Alfa POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'BETA'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                ts = handles.features{1,2}{1,1};
                P = handles.features{1,2}{1,3};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts,P(:,val1 - 1,5),'Color',[0.8 0.8 0.2],'linewidth',2), xlabel('time [s]'), ylabel('Beta Power [mW]^2/Hz'), title(['Feature: Beta POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'GAMMA'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                ts = handles.features{1,2}{1,1};
                P = handles.features{1,2}{1,3};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts,P(:,val1 - 1,6),'Color',[0.8 0.2 0.2],'linewidth',2), xlabel('time [s]'), ylabel('Gamma Power [mW]^2/Hz'), title(['Feature: Gamma POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'GAMMA LOW'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                ts = handles.features{1,2}{1,1};
                P = handles.features{1,2}{1,3};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts,P(:,val1 - 1,1),'Color',[0.8 0.2 0.4],'linewidth',2), xlabel('time [s]'), ylabel('Gamma-Low Power [mW]^2/Hz'), title(['Feature: Low-Gamma POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'GAMMA HIGH'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                ts = handles.features{1,2}{1,1};
                P = handles.features{1,2}{1,3};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts,P(:,val1 - 1,1),'Color',[0.8 0.2 0.8],'linewidth',2), xlabel('time [s]'), ylabel('Gamma-High Power [mW]^2/Hz'), title(['Feature: High-Gamma POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'HbO AVERAGE'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                ts = handles.features{1,2}{1,2};
                H = handles.features{1,2}{1,5};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts,H(:,val1 - 1),'b','linewidth',2), xlabel('time [s]'), ylabel('HbO Average [mM/L*s]'),title(['HbO Average: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'HbO SLOPE'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                numchan = length(nirslabels);
                ts = handles.features{1,2}{1,2};
                H = handles.features{1,2}{1,5};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts,H(:,numchan + val1),'b','linewidth',2), xlabel('time [s]'), ylabel('HbO Slope [mM/L*s]'),title(['HbO Slope: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'HbR AVERAGE'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                numchan = length(nirslabels);
                ts = handles.features{1,2}{1,2};
                H = handles.features{1,2}{1,5};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts,H(:,2*numchan + val1),'b','linewidth',2), xlabel('time [s]'), ylabel('HbR Average [mM/L*s]'),title(['HbR Average: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'HbR SLOPE'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                numchan = length(nirslabels);
                ts = handles.features{1,2}{1,2};
                H = handles.features{1,2}{1,5};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts,H(:,3*numchan + val1),'b','linewidth',2), xlabel('time [s]'), ylabel('HbR Slope [mM/L*s]'),title(['HbR Slope: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
        end
        
end

else
    
    %NEED TO COMPUTE TRIALS AND
    %BLOCKS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                trials = get(handles.trials_menu,'String');  
                value = get(handles.trials_menu,'Value'); 
                selected_trial = str2num(trials{value});
                tbeg = handles.trials(selected_trial,1);
                tfin = handles.trials(selected_trial,2);
                
    
switch val_sig
    
    case 1
        
        switch str{val}
            case 'NONE'
                if handles.typeplot==1
                    plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                    set(plot_hide,'Visible','off');
                    legend('off');
                end
                if ~isempty(handles.hideplot) && handles.typeplot==2
                    axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                    set(axes_hide,'Visible','off');
                    legend('off');
                end
                if handles.typeplot==3
                    plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                    plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                    set(plot_hide1,'Visible','off');
                    set(plot_hide2,'Visible','off');
                    legend(handles.sub1,'off');
                    legend(handles.sub2,'off');
                end
            case 'EEG'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                tteeg = handles.datacells{1,2}{1,2};
                xxeeg = handles.datacells{1,2}{1,1};
                set(handles.viewplot,'Visible','on');
                
                inter = find(tteeg >= tbeg & tteeg <= tfin);
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,tteeg(inter),xxeeg(inter,val1 - 1),'Color',[0 0.5 0]), xlabel('time [s]'), ylabel('EEG [uV]'), title(['EEG: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                
                guidata(hObject, handles);
            case 'HbO'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                ttnirs = handles.datacells{1,2}{1,5};
                xxnirs = handles.datacells{1,2}{1,3};
                set(handles.viewplot,'Visible','on');
                
                inter = find(ttnirs >= tbeg & ttnirs <= tfin);
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ttnirs(inter),xxnirs(inter,val1 - 1),'r','linewidth',2), xlabel('time [s]'), ylabel('HbO [mM/L]'),title(['HbO: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
            case 'HbR'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                ttnirs = handles.datacells{1,2}{1,5};
                xxnirsR = handles.datacells{1,2}{1,4};
                set(handles.viewplot,'Visible','on');
                
                inter = find(ttnirs >= tbeg & ttnirs <= tfin);
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ttnirs(inter),xxnirsR(inter,val1 - 1),'b','linewidth',2), xlabel('time [s]'), ylabel('HbR [mM/L]'),title(['HbR: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
            case 'HbT'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                ttnirs = handles.datacells{1,2}{1,5};
                xxnirs = handles.datacells{1,2}{1,3};
                xxnirsR = handles.datacells{1,2}{1,4};
                xxnirsT = xxnirs - xxnirsR;
                set(handles.viewplot,'Visible','on');
                
                inter = find(ttnirs >= tbeg & ttnirs <= tfin);
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ttnirs(inter),xxnirsT(inter,val1 - 1),'g','linewidth',2), xlabel('time [s]'), ylabel('HbT [mM/L]'),title(['HbT: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                end
                guidata(hObject, handles);
            case 'EEG + HbO'
                if get(handles.channel2,'Value')~=1
                    val1 = get(handles.channel1,'Value');
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    tteeg = handles.datacells{1,2}{1,2};
                    xxeeg = handles.datacells{1,2}{1,1};
                    
                    intereeg = find(tteeg >= tbeg & tteeg <= tfin);
                    
                    val2 = get(handles.channel2,'Value');
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    ttnirs = handles.datacells{1,2}{1,5};
                    xxnirs = handles.datacells{1,2}{1,3};
                    set(handles.viewplot,'Visible','on');
                    
                    internirs = find(ttnirs >= tbeg & ttnirs <= tfin);
                    
                    if val1 ~= 1 && val2 ~= 1
                        axes(handles.viewplot);
                        [hAx,hLine1,hLine2] = plotyy(handles.viewplot,tteeg(intereeg),xxeeg(intereeg,val1 - 1),ttnirs(internirs),xxnirs(internirs,val2 - 1),'plot');
                        set(hLine1,'color',[0 0.5 0]);
                        set(hLine2,'color','red');
                        set(hAx,{'ycolor'},{[0 0.5 0];'r'}) ;
                        ylabel(hAx(1),'EEG [uV]') % left y-axis
                        ylabel(hAx(2),'HbO [mM/L]') % right y-axis
                        legend('EEG','HbO');
                        xlabel('time [s]'), title(['EEG: ' eeglabels{val1-1} ' HbO: ' nirslabels{val2-1}]);
                        grid(handles.viewplot,'on');
                        handles.hideplot = hAx;
                        handles.typeplot=2;
                    end
                    guidata(hObject, handles);
                end
            case 'EEG + HbR'
                if get(handles.channel2,'Value')~=1
                    val1 = get(handles.channel1,'Value');
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    tteeg = handles.datacells{1,2}{1,2};
                    xxeeg = handles.datacells{1,2}{1,1};
                    
                    intereeg = find(tteeg >= tbeg & tteeg <= tfin);
                    
                    val2 = get(handles.channel2,'Value');
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    ttnirs = handles.datacells{1,2}{1,5};
                    xxnirsR = handles.datacells{1,2}{1,4};
                    set(handles.viewplot,'Visible','on');
                    
                    internirs = find(ttnirs >= tbeg & ttnirs <= tfin);
                    
                    if val1 ~= 1 && val2 ~= 1
                        axes(handles.viewplot);
                        [hAx,hLine1,hLine2] = plotyy(handles.viewplot,tteeg(intereeg),xxeeg(intereeg,val1 - 1),ttnirs(internirs),xxnirsR(internirs,val2 - 1),'plot');
                        set(hLine1,'color',[0 0.5 0]);
                        set(hLine2,'color','blue');
                        set(hAx,{'ycolor'},{[0 0.5 0];'b'}) ;
                        ylabel(hAx(1),'EEG [uV]') % left y-axis
                        ylabel(hAx(2),'HbR [mM/L]') % right y-axis
                        legend('EEG','HbR');
                        xlabel('time [s]'), title(['EEG: ' eeglabels{val1-1} ' HbR: ' nirslabels{val2-1}]);
                        grid(handles.viewplot,'on');
                        
                        handles.hideplot = hAx;
                        handles.typeplot=2;
                    end
                    
                    guidata(hObject, handles);
                end
            case 'HbO + HbR'
                if get(handles.channel2,'Value')~=1
                    val1 = get(handles.channel1,'Value');
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    ttnirs = handles.datacells{1,2}{1,5};
                    xxnirs = handles.datacells{1,2}{1,3};
                    
                    val2 = get(handles.channel2,'Value');
                    xxnirsR = handles.datacells{1,2}{1,4};
                    set(handles.viewplot,'Visible','off');
                    
                    inter = find(ttnirs >= tbeg & ttnirs <= tfin);
                    
                    if val1 ~= 1 && val2 ~= 1
                        
                        
                        set(handles.sub1,'Visible','on');
                        axes(handles.sub1);
                        plot(ttnirs(inter),xxnirsR(inter,val1 - 1),'b');
                        hold on
                        plot(ttnirs(inter),xxnirs(inter,val1 - 1),'r');
                        hold off
                        
                        legend('HbO','HbR');
                        ylabel('Hb [mM/L]'), title(['HbO: ' nirslabels{val1-1} ' HbR: ' nirslabels{val1-1}]);
                        grid(handles.sub1,'on');
                        
                        set(handles.sub2,'Visible','on');
                        axes(handles.sub2)
                        plot(ttnirs(inter),xxnirs(inter,val2 - 1),'r');
                        hold on
                        plot(ttnirs(inter),xxnirsR(inter,val2 - 1),'b');
                        hold off
                        legend('HbO','HbR');
                        xlabel('time [s]'), ylabel('Hb [mM/L]'), title(['HbO: ' nirslabels{val2-1} ' HbR: ' nirslabels{val2-1}]);
                        grid(handles.sub2,'on');
                        
                        handles.typeplot = 3;
                    end
                    
                    guidata(hObject, handles);
                end
                
            case 'EEG + EEG'
                if get(handles.channel2,'Value')~=1
                    val1 = get(handles.channel1,'Value');
                    eeglabels = handles.datacells{1,1}{1,1}{4,1};
                    tteeg = handles.datacells{1,2}{1,2};
                    xxeeg = handles.datacells{1,2}{1,1};
                    
                    val2 = get(handles.channel2,'Value');
                    set(handles.viewplot,'Visible','on');
                    
                    inter = find(tteeg >= tbeg & tteeg <= tfin);
                    
                    if val1 ~= 1 && val2 ~= 1
                        axes(handles.viewplot);
                        plot(handles.viewplot,tteeg(inter),xxeeg(inter,val1 - 1),'Color',[0 0.5 0])
                        hold on
                        plot(handles.viewplot,tteeg(inter),xxeeg(inter,val2 - 1),'Color',[0 0.2 0]);
                        hold off
                        
                        legend('EEG-1','EEG-2');
                        xlabel('time [s]'), ylabel('EEG [uV]'), title(['EEG-1: ' eeglabels{val1-1} ' EEG-2: ' eeglabels{val2-1}]);
                        grid(handles.viewplot,'on');
                        
                        handles.typeplot=1;
                    end
                    guidata(hObject, handles);
                end
                
            case 'HbO + HbO'
                if get(handles.channel2,'Value')~=1
                    val1 = get(handles.channel1,'Value');
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    ttnirs = handles.datacells{1,2}{1,5};
                    xxnirs = handles.datacells{1,2}{1,3};
                    
                    val2 = get(handles.channel2,'Value');
                    set(handles.viewplot,'Visible','on');
                    
                    inter = find(ttnirs >= tbeg & ttnirs <= tfin);
                    
                    if val1 ~= 1 && val2 ~= 1
                        axes(handles.viewplot);
                        plot(handles.viewplot,ttnirs(inter),xxnirs(inter,val1 - 1),'r');
                        hold on
                        plot(ttnirs(inter),xxnirs(inter,val2 - 1),'Color',[0.5,0,0.5]);
                        hold off
                        
                        legend('HbO-1','HbO-2');
                        xlabel('time [s]'), ylabel('HbO [mM/L]'), title(['HbO-1: ' nirslabels{val1-1} ' HbO-2: ' nirslabels{val2-1}]);
                        grid(handles.viewplot,'on');
                        
                        handles.typeplot=1;
                    end
                    guidata(hObject, handles);
                end
                
            case 'HbR + HbR'
                if get(handles.channel2,'Value')~=1
                    val1 = get(handles.channel1,'Value');
                    nirslabels = handles.datacells{1,1}{2,1}{4,1};
                    ttnirs = handles.datacells{1,2}{1,5};
                    xxnirsR = handles.datacells{1,2}{1,4};
                    
                    val2 = get(handles.channel2,'Value');
                    set(handles.viewplot,'Visible','on');
                    
                    inter = find(ttnirs >= tbeg & ttnirs <= tfin);
                    
                    if val1 ~= 1 && val2 ~= 1
                        axes(handles.viewplot);
                        plot(handles.viewplot,ttnirs(inter),xxnirsR(inter,val1 - 1),'b');
                        hold on
                        plot(ttnirs(inter),xxnirsR(inter,val2 - 1),'Color',[0,0.5,0.5]);
                        hold off
                        
                        legend('HbR-1','HbR-2');
                        xlabel('time [s]'), ylabel('HbR [mM/L]'), title(['HbR-1: ' nirslabels{val1-1} ' HbR-2: ' nirslabels{val2-1}]);
                        grid(handles.viewplot,'on');
                        
                        handles.typeplot=1;
                    end
                    guidata(hObject, handles);
                end
                
                guidata(hObject, handles);
                
        end
        
    case 2
        ts_nirs = handles.features{1,2}{1,2};
        ts_eeg = handles.features{1,2}{1,1};
        P = handles.features{1,2}{1,3};
        H = handles.features{1,2}{1,5};
                        
        intereeg = find(ts_eeg >= tbeg & ts_eeg <= tfin);              
        internirs = find(ts_nirs >= tbeg & ts_nirs <= tfin);
                
        switch str{val}
            case 'NONE'
                if handles.typeplot==1
                    plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                    set(plot_hide,'Visible','off');
                    legend('off');
                end
                if ~isempty(handles.hideplot) && handles.typeplot==2
                    axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                    set(axes_hide,'Visible','off');
                    legend('off');
                end
                if handles.typeplot==3
                    plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                    plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                    set(plot_hide1,'Visible','off');
                    set(plot_hide2,'Visible','off');
                    legend(handles.sub1,'off');
                    legend(handles.sub2,'off');
                end
            case 'TOTAL POWER'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                set(handles.viewplot,'Visible','on');
                
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts_eeg(intereeg),P(intereeg,val1 - 1,1),'Color',[0.6 0.6 0.6],'linewidth',2), xlabel('time [s]'), ylabel('Total Power [mW]^2/Hz'), title(['Feature: TOTAL POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
            case 'DELTA'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts_eeg(intereeg),P(intereeg,val1 - 1,2),'Color',[0.5 0.5 0.8],'linewidth',2), xlabel('time [s]'), ylabel('Delta Power [mW]^2/Hz'), title(['Feature: Delta POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'THETA'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts_eeg(intereeg),P(intereeg,val1 - 1,3),'Color',[0.5 0.8 0.8],'linewidth',2), xlabel('time [s]'), ylabel('Theta Power [mW]^2/Hz'), title(['Feature: Theta POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'ALFA'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts_eeg(intereeg),P(intereeg,val1 - 1,4),'Color',[0.2 0.8 0.2],'linewidth',2), xlabel('time [s]'), ylabel('Alfa Power [mW]^2/Hz'), title(['Feature: Alfa POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'BETA'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts_eeg(intereeg),P(intereeg,val1 - 1,5),'Color',[0.8 0.8 0.2],'linewidth',2), xlabel('time [s]'), ylabel('Beta Power [mW]^2/Hz'), title(['Feature: Beta POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'GAMMA'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts_eeg(intereeg),P(intereeg,val1 - 1,6),'Color',[0.8 0.2 0.2],'linewidth',2), xlabel('time [s]'), ylabel('Gamma Power [mW]^2/Hz'), title(['Feature: Gamma POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'GAMMA LOW'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts_eeg(intereeg),P(intereeg,val1 - 1,1),'Color',[0.8 0.2 0.4],'linewidth',2), xlabel('time [s]'), ylabel('Gamma-Low Power [mW]^2/Hz'), title(['Feature: Low-Gamma POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'GAMMA HIGH'
                val1 = get(handles.channel1,'Value');
                eeglabels = handles.datacells{1,1}{1,1}{4,1};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts_eeg(intereeg),P(intereeg,val1 - 1,1),'Color',[0.8 0.2 0.8],'linewidth',2), xlabel('time [s]'), ylabel('Gamma-High Power [mW]^2/Hz'), title(['Feature: High-Gamma POWER Chan: ' eeglabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'HbO AVERAGE'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};                
                
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts_nirs(internirs),H(internirs,val1 - 1),'b','linewidth',2), xlabel('time [s]'), ylabel('HbO Average [mM/L*s]'),title(['HbO Average: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'HbO SLOPE'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                numchan = length(nirslabels);
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts_nirs(internirs),H(internirs,numchan + val1),'b','linewidth',2), xlabel('time [s]'), ylabel('HbO Slope [mM/L*s]'),title(['HbO Slope: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'HbR AVERAGE'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                numchan = length(nirslabels);
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts_nirs(internirs),H(internirs,2*numchan + val1),'b','linewidth',2), xlabel('time [s]'), ylabel('HbR Average [mM/L*s]'),title(['HbR Average: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
                
            case 'HbR SLOPE'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                numchan = length(nirslabels);
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ts_nirs(internirs),H(internirs,3*numchan + val1),'b','linewidth',2), xlabel('time [s]'), ylabel('HbR Slope [mM/L*s]'),title(['HbR Slope: ' nirslabels{val1-1}]);
                    grid(handles.viewplot,'on');
                    handles.typeplot=1;
                else
                    if handles.typeplot==1
                        plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                        set(plot_hide,'Visible','off');
                        legend('off');
                    end
                    if ~isempty(handles.hideplot) && handles.typeplot==2
                        axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                        set(axes_hide,'Visible','off');
                        legend('off');
                    end
                    if handles.typeplot==3
                        plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                        plot_hide2 = [handles.sub2;get(handles.sub1,'Children')];
                        set(plot_hide1,'Visible','off');
                        set(plot_hide2,'Visible','off');
                        legend('off');
                    end
                end
                guidata(hObject, handles);
        end
        
end

%restore events
eventsbutton_Callback(hObject, eventdata, handles);
    
end

% --- Executes during object creation, after setting all properties.
function channel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in channel2.
function channel2_Callback(hObject, eventdata, handles)
% hObject    handle to channel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channel2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channel2

channel1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function channel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in signal_selection.
function signal_selection_Callback(hObject, eventdata, handles)
% hObject    handle to signal_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns signal_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from signal_selection

%the callback will be called also by feature selection callback -> can't
%use hObject
val = get(handles.signal_selection,'Value');

switch (val)
    
    case 1
        
        if handles.typeplot==1
            plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
            set(plot_hide,'Visible','off');
            legend('off');
        end
        if ~isempty(handles.hideplot) && handles.typeplot==2
            axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
            set(axes_hide,'Visible','off');
            legend('off');
        end
        if handles.typeplot==3
            plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
            plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
            set(plot_hide1,'Visible','off');
            set(plot_hide2,'Visible','off');
            legend(handles.sub1,'off');
            legend(handles.sub2,'off');
        end
        set(handles.select_view,'Visible','on');
        set(handles.channel2,'Visible','off');
        set(handles.select_view,'String',handles.prep_select_view);
        set(handles.select_view,'Value',1);
        set(handles.channel1,'Value',1);
        
    case 2
        
        if handles.featurecomputed
            
            if handles.typeplot==1
                plot_hide = [handles.viewplot;get(handles.viewplot,'Children')];
                set(plot_hide,'Visible','off');
                legend('off');
            end
            if ~isempty(handles.hideplot) && handles.typeplot==2
                axes_hide = [handles.hideplot(1);handles.hideplot(2);get(handles.hideplot(1),'Children');get(handles.hideplot(2),'Children')];
                set(axes_hide,'Visible','off');
                legend('off');
            end
            if handles.typeplot==3
                plot_hide1 = [handles.sub1;get(handles.sub1,'Children')];
                plot_hide2 = [handles.sub2;get(handles.sub2,'Children')];
                set(plot_hide1,'Visible','off');
                set(plot_hide2,'Visible','off');
                legend(handles.sub1,'off');
                legend(handles.sub2,'off');
            end
            
            % val_feat is 1 for eeg_nirs, 2 for eeg only and 3 for nirs only
            
            val_feat = handles.features{1,1}{1,4};
            
            switch val_feat
                
                case 1
                    set(handles.select_view,'Visible','on');
                    set(handles.channel2,'Visible','off');
                    set(handles.select_view,'String',handles.features_select_view);
                    set(handles.select_view,'Value',1);
                    set(handles.channel1,'Value',1);
                case 2
                    set(handles.select_view,'Visible','on');
                    set(handles.channel2,'Visible','off');
                    set(handles.select_view,'String',handles.eegfeatures_select_view);
                    set(handles.select_view,'Value',1);
                    set(handles.channel1,'Value',1);
                case 3
                    set(handles.select_view,'Visible','on');
                    set(handles.channel2,'Visible','off');
                    set(handles.select_view,'String',handles.nirsfeatures_select_view);
                    set(handles.select_view,'Value',1);
                    set(handles.channel1,'Value',1);
            end
            
        else
            msgbox('Compute features before!!!');
        end
end

guidata(hObject,handles);

% --- Executes on button press in eventsbutton.
function eventsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to eventsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eventsbutton

if (~get(handles.trial_by_trial,'Value'))

if get(hObject,'Value')    
    
    if handles.typeplot == 1
        
        evt = handles.datacells{1,2}{1,6};
        
        axes(handles.viewplot);
        
        hold on
        hx = graph2d.constantline(evt(:,1), 'LineStyle','-', 'Color',[0 0 0], 'linewidth',2);
        changedependvar(hx,'x');
        hold off
        
    elseif    handles.typeplot == 2
        
        evt = handles.datacells{1,2}{1,6};
        
        axes(handles.hideplot(1));
        
        hold on
        hx = graph2d.constantline(evt(:,1), 'LineStyle','-', 'Color',[0 0 0], 'linewidth',2);
        changedependvar(hx,'x');
        hold off
        
        axes(handles.hideplot(2));
        
        hold on
        hx = graph2d.constantline(evt(:,1), 'LineStyle','-', 'Color',[0 0 0], 'linewidth',2);
        changedependvar(hx,'x');
        hold off
        
        
    elseif    handles.typeplot == 3
        
        evt = handles.datacells{1,2}{1,6};
        
        axes(handles.sub1);
        
        hold on
        hx = graph2d.constantline(evt(:,1), 'LineStyle','-', 'Color',[0 0 0], 'linewidth',2);
        changedependvar(hx,'x');
        hold off
        
        axes(handles.sub2);
        
        hold on
        hx = graph2d.constantline(evt(:,1), 'LineStyle','-', 'Color',[0 0 0], 'linewidth',2);
        changedependvar(hx,'x');
        hold off
        
    end
    
else
    
    channel1_Callback(hObject, eventdata, handles)
    
end

else
    %find trial interval
                trials = get(handles.trials_menu,'String');  
                value = get(handles.trials_menu,'Value'); 
                selected_trial = str2num(trials{value});
                tbeg = handles.trials(selected_trial,1);
                tfin = handles.trials(selected_trial,2);
                evt = handles.datacells{1,2}{1,6};
                
                inter = find(evt(:,1)>=tbeg & evt(:,1)<=tfin);
                
  if get(hObject,'Value')    
    
    if handles.typeplot == 1     
        
        axes(handles.viewplot);
        
        hold on
        hx = graph2d.constantline(evt(inter,1), 'LineStyle','-', 'Color',[0 0 0], 'linewidth',2);
        changedependvar(hx,'x');
        hold off
        
    elseif    handles.typeplot == 2
        
        axes(handles.hideplot(1));
        
        hold on
        hx = graph2d.constantline(evt(inter,1), 'LineStyle','-', 'Color',[0 0 0], 'linewidth',2);
        changedependvar(hx,'x');
        hold off
        
        axes(handles.hideplot(2));
        
        hold on
        hx = graph2d.constantline(evt(inter,1), 'LineStyle','-', 'Color',[0 0 0], 'linewidth',2);
        changedependvar(hx,'x');
        hold off
        
        
    elseif    handles.typeplot == 3
        
        axes(handles.sub1);
        
        hold on
        hx = graph2d.constantline(evt(inter,1), 'LineStyle','-', 'Color',[0 0 0], 'linewidth',2);
        changedependvar(hx,'x');
        hold off
        
        axes(handles.sub2);
        
        hold on
        hx = graph2d.constantline(evt(inter,1), 'LineStyle','-', 'Color',[0 0 0], 'linewidth',2);
        changedependvar(hx,'x');
        hold off
        
    end
    
else
    
    channel1_Callback(hObject, eventdata, handles)
    
  end             
                
              
end

% --- Executes during object creation, after setting all properties.
function signal_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to signal_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in view_button.
function view_button_Callback(hObject, eventdata, handles)
% hObject    handle to view_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of view_button

if get(hObject, 'Value')
    set(handles.view_panel,'Visible','on');
    set(handles.sub1,'Visible','off');
    set(handles.sub2,'Visible','off');
else
    set(handles.view_panel,'Visible','off');
end


% --- Executes on button press in expand_butt.
function expand_butt_Callback(hObject, eventdata, handles)
% hObject    handle to expand_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of expand_butt
if get(hObject,'Value')
    set(handles.processingpanel,'Visible','off');
    set(handles.featurespanel,'Visible','off');
    set(handles.view_button,'Visible','off');
    set(handles.filenamespanel,'Visible','off');
    set(handles.featurefilepanel,'Visible','off');
    set(handles.view_panel,'Position',[0 0 1 1]);
else
    set(handles.view_panel,'Position',handles.viewpanelpos);
    set(handles.processingpanel,'Visible','on');
    set(handles.featurespanel,'Visible','on');    
    set(handles.filenamespanel,'Visible','on');
    set(handles.featurefilepanel,'Visible','on');
    set(handles.view_button,'Visible','on');
end

guidata(hObject,handles);


% --- Executes on button press in features_butt.
function features_butt_Callback(hObject, eventdata, handles)
% hObject    handle to features_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%waitbar
    h = waitbar(0,'Please wait...');

% time window in sec
tw = handles.tw;

% overlap in %
overlap = handles.overlap;

noftw = handles.noftw;

% LOAD variables from datacells.mat

% eeg info

Fs = handles.datacells{1,1}{1,1}{2,1};
neegchan = handles.datacells{1,1}{1,1}{3,1};
eeglabels = handles.datacells{1,1}{1,1}{4,1};

Fsnirs = handles.datacells{1,1}{2,1}{2,1};
nnirschan = handles.datacells{1,1}{2,1}{3,1};
nirslabels = handles.datacells{1,1}{2,1}{4,1};

xxeeg = handles.xxeeg_current;
tteeg = handles.datacells{1,2}{1,2};
xxnirs = handles.datacells{1,2}{1,3};
xxnirsR = handles.datacells{1,2}{1,4};
ttnirs = handles.datacells{1,2}{1,5};
evt = handles.datacells{1,2}{1,6};


%round events at time window resolution
acc = tw;
evt(:,1) = round(evt(:,1)/acc)*acc;

val_feat = get(handles.feature_selection,'Value');

switch val_feat
    
    case 1
        
        % EEG POWER
        % calc power time course
        nweegpower = floor(tw * Fs);
        noverlap = floor( (overlap / 100) * nweegpower);
        
        freqs = [];
        numfreqbands = 8;
        
        for kchan = 1:neegchan,
            fprintf('Calculating eeg power, chan %d/%d\n',kchan,neegchan);
            %freqs = [];
            [SS,Freq,Tim,Pow] = spectrogram(xxeeg(:,kchan),nweegpower,noverlap,freqs,Fs);
            %Pow = Pow/max(max(Pow));
            Pow = 10.*log(Pow/min(min(Pow)));
            
            if kchan==1
                xxeegbp = zeros ( length(Pow), neegchan, numfreqbands );
            end
            
            % loop over frequency bands
            for ifreqs = 1:numfreqbands
                if 1==ifreqs, freqlolimit = 0; freqhilimit = 80;
                elseif 2==ifreqs, freqlolimit = 0; freqhilimit = 4;
                elseif 3==ifreqs, freqlolimit = 4; freqhilimit = 8;
                elseif 4==ifreqs, freqlolimit = 8; freqhilimit = 12;
                elseif 5==ifreqs, freqlolimit = 12; freqhilimit = 30;
                elseif 6==ifreqs, freqlolimit = 30; freqhilimit = 80;
                elseif 7==ifreqs, freqlolimit = 30; freqhilimit = 50;
                elseif 8==ifreqs, freqlolimit = 50; freqhilimit = 80;
                end
                
                %EEG power  ( time x chans x freqbands )
                indxs = find( freqlolimit<Freq & Freq<=freqhilimit );
                totalbandpower = sum( Pow(indxs,:), 1 );
                xxeegbp(:,kchan,ifreqs) = totalbandpower;
                %xxeegbp(:,kchan,ifreqs) = smooth(totalbandpower,5);
            end
            
            
        waitbar(kchan/neegchan*0.25);
            
        end
        
        % Create TIME SAMPLES
        
        if overlap~=0
            acc = tw * overlap/100;
        else
            acc = tw;
        end
        
        % BLOCKS and TRIALS DVISION
        
        % Discard first event (sync)
        
        %DIVIDE IN BLOCKS
        
        evt = evt(2:end,:);
        
        begin_block_ind = find(evt(:,2)==10);
        begin_block = evt(begin_block_ind,1);
        nblocks = length(begin_block);
        
        %blocks contains the beginning and end of each block
        blocks = zeros(nblocks-1,2);
        block_length = evt(end,1) - begin_block(end);
        
        for i=1:nblocks
            if i~=nblocks
                blocks(i,:) = begin_block(i:i+1);
            else
                blocks(i,:) = [begin_block(i), begin_block(i)+block_length];
            end
        end
        
        %DIVIDE EACH BLOCK IN 10 TRIALS
        
        %each trial contains: beginning and end, belonging block, motor task
        %performed
        
        begin_trial_ind = find(evt(:,2)==1 & evt(:,1)<=blocks(end,2));
        begin_trial = evt(begin_trial_ind,1);
        
        ntrials = length(begin_trial);
        trials = zeros(ntrials,4);
        trial_length = begin_trial(2) - begin_trial(1);
        
        for i = 1:ntrials
            if i~=ntrials
                trials(i,1:2) = begin_trial(i:i+1);
                trials(i,3) = find(blocks(:,1)<trials(i,1) & trials(i,1)<blocks(:,2));
            else
                trials(i,1:2) = [begin_trial(i), begin_trial(i)+trial_length];
                trials(i,3) = find(blocks(:,1)<trials(i,1) & trials(i,1)<blocks(:,2));
            end
        end
        
        %ADD MOTOR TASK LABEL TO EACH TRIAL
        
        %for the different motor tasks (event 3->8)
        for j = 3:8
            timeoftask = evt(evt(:,2)==j,1);
            
            for l = 1:length(timeoftask)
                trials(trials(:,1)<timeoftask(l) & timeoftask(l)<trials(:,2),4) = j;
            end
            
        end
        
        blocks(end,2) = trials(end,2);
        
        %create trial menu string
        string_trial = cell(1,length(trials));
        
        for i=1:length(trials)
            string_trial{i} = i; 
        end
        
        
        %compute the initial resting state average to be subtracted
        ttmean = find (ttnirs < 10);
        restnirs = mean (xxnirs(ttmean, :));
        restnirsR = mean (xxnirsR(ttmean, :));
        
        waitbar(0.5);
        
        % DIFFERENT organization of the features depending on the number of time
        % windows selected
        
        switch (noftw)
            
            case 1
                
                ts1 = linspace(tw/2, tteeg(end) - tw/2, size(xxeegbp,1));
                ts = ts1;
                
                % Calculate NIRS features
                
                % 1 TW
                % average & slope
                
                xxnirsav = zeros (length(ts1),nnirschan);
                xxnirsRav = zeros (length(ts1),nnirschan);
                xxnirssl = zeros (length(ts1),nnirschan);
                xxnirsRsl = zeros (length(ts1),nnirschan);
                
                
                for ichan = 1:nnirschan,
                    for tss = 1 : length(ts1)
                        % Average
                        
                        % substitute the initial average with the average
                        % of the "ready"
                        
                        indx = find (ttnirs >= ts1(tss)-tw/2 & ttnirs < ts1(tss)+tw/2);
                        
                        xxnirsav(tss,ichan) = mean (xxnirs(indx,ichan)) - restnirs(ichan);
                        xxnirsRav(tss,ichan) = mean (xxnirsR(indx,ichan)) - restnirsR(ichan);
                        
                        % Slope
                        xxnirssl(tss,ichan) = xxnirs(indx(end),ichan) - xxnirs(indx(1),ichan);
                        xxnirsRsl(tss,ichan) = xxnirsR(indx(end),ichan) - xxnirsR(indx(1),ichan);
                        
                        
                        %coefficient of the linear regression curve
                        %         P = polyfit(ttnirs(indx),xxnirs(indx),1);
                        %         PR = polyfit(ttnirs(indx),xxnirs(indx),1);
                        %
                        %         xxnirssl(tss,ichan) = P(1);
                        %         xxnirsRsl(tss,ichan) = P(2);
                        
                                                                       
                    end
                    waitbar(0.5 + ichan/nnirschan*0.25);
                end
                
                
                % CREATE FEATURES MATRIXES
                
                % 	CREATE LABELS array from events
                
                labels = zeros (length(ts1),1);
                
                for i = 1:length(trials)
                    
                    ind_beg = find(ts1==trials(i,2)-6);
                    ind_end = find(ts1==trials(i,2));
                    
                    labels(ind_beg:ind_end) = trials(i,4);
                    
                end
                
                % UNWRAP EEGBANDPOWER
                
                
                % 1TW: eegbp = ts x neegchan x nfreq --> all in one line
                
                Xeeg = zeros(length(ts1),neegchan*numfreqbands);
                
                
                for ii = 1:neegchan
                    
                    Xeeg ( : , (ii-1)*numfreqbands + 1 : ii*numfreqbands) = xxeegbp ( : , ii , :);
                    
                end
                
                %Standardization
                
                % for r = 1:size(Xeeg1,1)
                % Xeeg1(r,:) = (Xeeg1(r,:) - mean(Xeeg1))/std(Xeeg1);
                % end
                
                % UNWRAP NIRS AVERAGE AND SLOPE
                
                Xnirs = zeros(length(ts1) , nnirschan*4); % 4 = 2 features x Hbo and Hb chans
                
                Xnirs ( : , 1 : nnirschan) = xxnirsav ( : , :);
                Xnirs ( : , nnirschan + 1 : nnirschan*2) = xxnirssl ( : , :);
                Xnirs ( : , nnirschan*2 + 1 : nnirschan*3) = xxnirsRav ( : , :);
                Xnirs ( : , nnirschan*3 + 1 : end) = xxnirsRsl ( : , :);
                
                
                
            case 2
                
                %compute ts1 to display eeg powers
                ts1 = linspace(tw/2, tteeg(end) - tw/2, size(xxeegbp,1));
                
                ts = linspace(tw/2, tteeg(end) - tw/2, floor(size(xxeegbp,1)/2));
                ts = round(ts/acc)*acc;
                
                % 2 TW
                % average & slope
                
                xxnirsav = zeros(length(ts),nnirschan);
                xxnirsRav = zeros(length(ts),nnirschan);
                xxnirssl = zeros(length(ts),nnirschan);
                xxnirsRsl = zeros(length(ts),nnirschan);
                
                
                for ichan = 1:nnirschan,
                    for tss = 1 : length(ts)
                        % Average
                        indx = find (ttnirs >= ts(tss)-tw & ttnirs < ts(tss)+tw);
                        
                        xxnirsav(tss,ichan) = mean (xxnirs(indx,ichan)) - restnirs(ichan);
                        xxnirsRav(tss,ichan) = mean (xxnirsR(indx,ichan)) - restnirsR(ichan);
                        
                        % Slope
                        xxnirssl(tss,ichan) = xxnirs(indx(end),ichan) - xxnirs(indx(1),ichan);
                        xxnirsRsl(tss,ichan) = xxnirsR(indx(end),ichan) - xxnirsR(indx(1),ichan);
                        
                        %         %coefficient of the linear regression curve
                        %         P = polyfit(ttnirs(indx),xxnirs(indx),1);
                        %         PR = polyfit(ttnirs(indx),xxnirs(indx),1);
                        %
                        %         xxnirssl(tss,ichan) = P(1);
                        %         xxnirsRsl(tss,ichan) = P(2);
                        %
                        %         nofstw2(tss) = length(indx);
                        
                    end
                    waitbar(0.5 + ichan/nnirschan*0.25);
                end
                
                % CREATE FEATURES MATRIXES
                
                % 	CREATE LABELS array from events
                
                labels = zeros (length(ts),1);
                
                
                for i = 1:length(trials)
                    
                    ind_beg = find(ts>=trials(i,2)-6.5 & ts<=trials(i,2)-6);
                    ind_end = find(ts>=trials(i,2)-0.5 & ts<=trials(i,2));
                    
                    labels(ind_beg:ind_end) = trials(i,4);
                    
                end
                
                % UNWRAP EEGBANDPOWER
                
                Xeeg = zeros(length(ts),neegchan*numfreqbands*2);
                
                for jj = 1 :size(Xeeg,1)
                    for ii = 1:neegchan
                        
                        Xeeg ( jj , (ii-1)*numfreqbands + 1 : ii*numfreqbands) = xxeegbp ( 2*(jj - 1) + 1 , ii , :);
                        Xeeg ( jj , (ii-1)*numfreqbands + 1 + neegchan*numfreqbands  : ii*numfreqbands + neegchan*numfreqbands) = xxeegbp ( 2*(jj - 1) + 2 , ii , :);
                        
                    end
                end
                
                %Standardization
                
                % for r = 1:size(Xeeg2,1)
                % Xeeg2(r,:) = (Xeeg2(r,:) - mean(Xeeg2))/std(Xeeg2);
                % end
                
                % UNWRAP NIRS AVERAGE AND SLOPE
                
                Xnirs = zeros(length(ts) , nnirschan*4); % 4 = 2 features x Hbo and Hb chans
                
                Xnirs ( : , 1 : nnirschan) = xxnirsav ( : , :);
                Xnirs ( : , nnirschan + 1 : nnirschan*2) = xxnirssl ( : , :);
                Xnirs ( : , nnirschan*2 + 1 : nnirschan*3) = xxnirsRav ( : , :);
                Xnirs ( : , nnirschan*3 + 1 : end) = xxnirsRsl ( : , :);
                
                
                
            case 3
                
                %compute ts1 to display eeg powers
                ts1 = linspace(tw/2, tteeg(end) - tw/2, size(xxeegbp,1));
                
                ts = linspace(tw/2, tteeg(end) - tw/2, floor(size(xxeegbp,1)/3));
                ts = round(ts/acc)*acc;
                
                % 3 TW
                % average & slope
                
                xxnirsav = zeros(length(ts),nnirschan);
                xxnirsRav = zeros(length(ts),nnirschan);
                xxnirssl = zeros(length(ts),nnirschan);
                xxnirsRsl = zeros(length(ts),nnirschan);
                
                
                for ichan = 1:nnirschan,
                    for tss = 1 : length(ts)
                        % Average
                        indx = find (ttnirs >= ts(tss)- tw*1.5 & ttnirs < ts(tss)+tw*1.5);
                        
                        xxnirsav(tss,ichan) = mean (xxnirs(indx,ichan)) - restnirs(ichan);
                        xxnirsRav(tss,ichan) = mean (xxnirsR(indx,ichan)) - restnirsR(ichan);
                        
                        %Slope
                        xxnirssl(tss,ichan) = xxnirs(indx(end),ichan) - xxnirs(indx(1),ichan);
                        xxnirsRsl(tss,ichan) = xxnirsR(indx(end),ichan) - xxnirsR(indx(1),ichan);
                        
                        %         %coefficient of the linear regression curve
                        %         P = polyfit(ttnirs(indx),xxnirs(indx),1);
                        %         PR = polyfit(ttnirs(indx),xxnirs(indx),1);
                        %
                        %         xxnirssl(tss,ichan) = P(1);
                        %         xxnirsRsl(tss,ichan) = P(2);
                        %
                        %         nofstw3(tss) = length(indx);
                        
                    end
                    waitbar(0.5 + ichan/nnirschan*0.25);
                end
                
                % CREATE FEATURES MATRIXES
                
                % 	CREATE LABELS array from events
                
                labels = zeros (length(ts),1);
                
                for i = 1:length(trials)
                    
                    ind_beg = find(ts>=trials(i,2)-6.5 & ts<=trials(i,2)-5.5);
                    ind_end = find(ts>=trials(i,2)-0.5 & ts<=trials(i,2)+0.5);
                    
                    labels(ind_beg:ind_end) = trials(i,4);
                    
                end
                
                % UNWRAP EEGBANDPOWER
                
                Xeeg = zeros(length(ts),neegchan*numfreqbands*3);
                
                for jj = 1 :size(Xeeg,1)
                    for ii = 1:neegchan
                        
                        Xeeg ( jj , (ii-1)*numfreqbands + 1 : ii*numfreqbands) = xxeegbp ( 3*(jj - 1) + 1 , ii , :);
                        Xeeg ( jj , (ii-1)*numfreqbands + 1 + neegchan*numfreqbands  : ii*numfreqbands + neegchan*numfreqbands) = xxeegbp ( 3*(jj - 1) + 2 , ii , :);
                        Xeeg ( jj , (ii-1)*numfreqbands + 1 + 2*neegchan*numfreqbands  : ii*numfreqbands + 2*neegchan*numfreqbands) = xxeegbp ( 3*(jj - 1) + 3 , ii , :);
                    end
                end
                
                %Standardization
                
                % for r = 1:size(Xeeg3,1)
                % Xeeg3(r,:) = (Xeeg3(r,:) - mean(Xeeg3))/std(Xeeg3);
                % end
                
                % UNWRAP NIRS AVERAGE AND SLOPE
                
                Xnirs = zeros(length(ts) , nnirschan*4); % 4 = 2 features x Hbo and Hb chans
                
                Xnirs ( : , 1 : nnirschan) = xxnirsav ( : , :);
                Xnirs ( : , nnirschan + 1 : nnirschan*2) = xxnirssl ( : , :);
                Xnirs ( : , nnirschan*2 + 1 : nnirschan*3) = xxnirsRav ( : , :);
                Xnirs ( : , nnirschan*3 + 1 : end) = xxnirsRsl ( : , :);
                
                
            case 4
                
                %compute ts1 to display eeg powers
                ts1 = linspace(tw/2, tteeg(end) - tw/2, size(xxeegbp,1));
                
                ts = linspace(tw/2, tteeg(end) - tw/2, floor(size(xxeegbp,1)/4));
                ts = round(ts/acc)*acc;
                
                % 3 TW
                % average & slope
                
                xxnirsav = zeros(length(ts),nnirschan);
                xxnirsRav = zeros(length(ts),nnirschan);
                xxnirssl = zeros(length(ts),nnirschan);
                xxnirsRsl = zeros(length(ts),nnirschan);
                
                
                for ichan = 1:nnirschan,
                    for tss = 1 : length(ts)
                        % Average
                        indx = find (ttnirs >= ts(tss)- tw*2 & ttnirs < ts(tss)+tw*2);
                        
                        xxnirsav(tss,ichan) = mean (xxnirs(indx,ichan)) - restnirs(ichan);
                        xxnirsRav(tss,ichan) = mean (xxnirsR(indx,ichan)) - restnirsR(ichan);
                        
                        %Slope
                        xxnirssl(tss,ichan) = xxnirs(indx(end),ichan) - xxnirs(indx(1),ichan);
                        xxnirsRsl(tss,ichan) = xxnirsR(indx(end),ichan) - xxnirsR(indx(1),ichan);
                        
                        %         %coefficient of the linear regression curve
                        %         P = polyfit(ttnirs(indx),xxnirs(indx),1);
                        %         PR = polyfit(ttnirs(indx),xxnirs(indx),1);
                        %
                        %         xxnirssl(tss,ichan) = P(1);
                        %         xxnirsRsl(tss,ichan) = P(2);
                        %
                        %         nofstw3(tss) = length(indx);
                        
                    end
                    waitbar(0.5 + ichan/nnirschan*0.25);
                end
                
                % CREATE FEATURES MATRIXES
                
                % 	CREATE LABELS array from events
                
                labels = zeros (length(ts),1);
                
                for i = 1:length(trials)
                    
                    ind_beg = find(ts>=trials(i,2)-6.5 & ts<=trials(i,2)-5.5);
                    ind_end = find(ts>=trials(i,2)-0.5 & ts<=trials(i,2)+0.5);
                    
                    labels(ind_beg:ind_end) = trials(i,4);
                    
                end
                
                % UNWRAP EEGBANDPOWER
                
                Xeeg = zeros(length(ts),neegchan*numfreqbands*3);
                
                for jj = 1 :size(Xeeg,1)
                    for ii = 1:neegchan
                        
                        Xeeg ( jj , (ii-1)*numfreqbands + 1 : ii*numfreqbands) = xxeegbp ( 4*(jj - 1) + 1 , ii , :);
                        Xeeg ( jj , (ii-1)*numfreqbands + 1 + neegchan*numfreqbands  : ii*numfreqbands + neegchan*numfreqbands) = xxeegbp ( 4*(jj - 1) + 2 , ii , :);
                        Xeeg ( jj , (ii-1)*numfreqbands + 1 + 2*neegchan*numfreqbands  : ii*numfreqbands + 2*neegchan*numfreqbands) = xxeegbp ( 4*(jj - 1) + 3 , ii , :);
                        Xeeg ( jj , (ii-1)*numfreqbands + 1 + 3*neegchan*numfreqbands  : ii*numfreqbands + 3*neegchan*numfreqbands) = xxeegbp ( 4*(jj - 1) + 4 , ii , :);

                    end
                end
                
                %Standardization
                
                % for r = 1:size(Xeeg3,1)
                % Xeeg3(r,:) = (Xeeg3(r,:) - mean(Xeeg3))/std(Xeeg3);
                % end
                
                % UNWRAP NIRS AVERAGE AND SLOPE
                
                Xnirs = zeros(length(ts) , nnirschan*4); % 4 = 2 features x Hbo and Hb chans
                
                Xnirs ( : , 1 : nnirschan) = xxnirsav ( : , :);
                Xnirs ( : , nnirschan + 1 : nnirschan*2) = xxnirssl ( : , :);
                Xnirs ( : , nnirschan*2 + 1 : nnirschan*3) = xxnirsRav ( : , :);
                Xnirs ( : , nnirschan*3 + 1 : end) = xxnirsRsl ( : , :);
                
        end
        
        
        
    case 2
        
    case 3
        
end

waitbar(0.9);

infocell = cell(1,4);

infocell{1} = tw;
infocell{2} = overlap;
infocell{3} = noftw;
infocell{4} = val_feat;

featuredata = cell(1,8);

featuredata{1} = ts1;
featuredata{2} = ts;
featuredata{3} = xxeegbp;
featuredata{4} = Xeeg;
featuredata{5} = Xnirs;
featuredata{6} = labels;

handles.features{1} = infocell;
handles.features{2} = featuredata;

handles.trials = trials;
handles.string_trial = string_trial;
handles.current_string_trial = string_trial;
handles.from_trial = 1;
handles.to_trial = length(trials);

handles.featurecomputed = 1;

% update features file path text and save check box
    set(handles.feature_filepaths,'String', ['Features: ' 'No name' ]);
    set(handles.feature_save_check,'Value',0);
    set(handles.feature_save_check,'Visible','on');
    
waitbar(1);

msgbox('Features computed CORRECTLY!');

close(h);

guidata(hObject,handles);



% --- Executes on button press in create_arff_butt.
function create_arff_butt_Callback(hObject, eventdata, handles)
% hObject    handle to create_arff_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in viewfilter.
function viewfilter_Callback(hObject, eventdata, handles)
% hObject    handle to viewfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

open([handles.configdir '\SpatialFilter_' num2str(get(handles.choosefilt,'Value')) '.txt']);


% --- Executes on button press in newbutton.
function newbutton_Callback(hObject, eventdata, handles)
% hObject    handle to newbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = get(handles.choosefilt,'String');
last = str2num(str{end});

fid = fopen([handles.configdir '\SpatialFilter_' num2str(last + 1) '.txt'],'w');

fclose(fid);

open([handles.configdir '\SpatialFilter_' num2str(last + 1) '.txt']);

str{last+1} = num2str(last+1);

set(handles.choosefilt,'String',str);

guidata(hObject,handles);


% --- Executes on selection change in feature_selection.
function feature_selection_Callback(hObject, eventdata, handles)
% hObject    handle to feature_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns feature_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from feature_selection

signal_selection_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function feature_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to feature_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in filter_check.
function filter_check_Callback(hObject, eventdata, handles)
% hObject    handle to filter_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filter_check


% --- Executes on button press in feature_save_check.
function feature_save_check_Callback(hObject, eventdata, handles)
% hObject    handle to feature_save_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of feature_save_check


% --- Executes on button press in save_features.
function save_features_Callback(hObject, eventdata, handles)
% hObject    handle to save_features (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.featurecomputed
    data = handles.features;
    %check if Pre Processed Data folder exists
    if ~exist([handles.dir '\Features'],'dir')
        mkdir(handles.dir,'Features');
    end
    
    [file, path] = uiputfile('*.mat','Save features',[handles.dir '\Features']);
    if file ~=0
        save([path file],'data');
        
        % update file path text
        set(handles.feature_filepaths,'String', ['Features: ' file ]);
        
        set(handles.feature_save_check,'Value',1);
        set(handles.feature_save_check,'Visible','on');
        
        msgbox('SUCCESSFULLY saved features!');
    else
        
        % update file path text
        set(handles.filepaths,'String', ['Features: ' 'No name' ]);
        
        set(handles.save_check,'Value',0);
        set(handles.save_check,'Visible','on');
        msgbox('UNSUCCESSFULLY saved features!');
    end
else
    msgbox('WARNING! Compute features before saving them...');
end


% --- Executes on selection change in trials_menu.
function trials_menu_Callback(hObject, eventdata, handles)
% hObject    handle to trials_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trials_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trials_menu

channel1_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function trials_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trials_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in trial_by_trial.
function trial_by_trial_Callback(hObject, eventdata, handles)
% hObject    handle to trial_by_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trial_by_trial

if get(hObject,'Value')
    set(handles.selecttrialspanel,'Visible','on');
    set(handles.trials_menu,'String',handles.string_trial);
    set(handles.trials_menu,'Visible','on');
    set(handles.include_trial,'Visible','on');
else
    set(handles.selecttrialspanel,'Visible','off');
    set(handles.trials_menu,'Visible','off');
    set(handles.include_trial,'Visible','off');
end

guidata(hObject,handles);    
    
% --- Executes on button press in include_trial.
function include_trial_Callback(hObject, eventdata, handles)
% hObject    handle to include_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of include_trial



function from_trial_Callback(hObject, eventdata, handles)
% hObject    handle to from_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of from_trial as text
%        str2double(get(hObject,'String')) returns contents of from_trial as a double

handles.from_trial = str2num(get(hObject,'String'));
handles.current_string_trial = handles.string_trial(handles.from_trial:handles.to_trial);
set(handles.trials_menu,'String',handles.current_string_trial);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function from_trial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to from_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function to_trial_Callback(hObject, eventdata, handles)
% hObject    handle to to_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of to_trial as text
%        str2double(get(hObject,'String')) returns contents of to_trial as a double

handles.to_trial = str2num(get(hObject,'String'));
handles.current_string_trial = handles.string_trial(handles.from_trial:handles.to_trial);
set(handles.trials_menu,'String',handles.current_string_trial);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function to_trial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to to_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_filtered.
function save_filtered_Callback(hObject, eventdata, handles)
% hObject    handle to save_filtered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %check if Pre Processed Data folder exists
    if ~exist([handles.dir '\Filtered Data'],'dir')
        mkdir(handles.dir,'Filtered Data');
    end
    
    switch handles.current_ref
        
        case 1
            save_data = 0;
        case 2
            save_data = 1;
            suffix = '_CAR';
        case 3
            save_data = 1;
            suffix = '_SLFILT';
    end
    
    if save_data
    
    [file, path] = uiputfile('*.mat','Save filtered data',[handles.dir '\Filtered Data\' handles.preprocessedfile suffix]);
    if file ~=0
        
        data = handles.xxeeg_current;
        save([path file],'data');
        
        % update file path text
        set(handles.filetext,'String', ['Filtered file: ' file ]);
        
        set(handles.filter_check,'Value',1);
        set(handles.filter_check,'Visible','on');
                
        msgbox('SUCCESSFULLY saved features!');
    else
        
        % update file path text
        set(handles.filepaths,'String', ['Features: ' 'No name' ]);
        
        set(handles.save_check,'Value',0);
        set(handles.save_check,'Visible','on');
        msgbox('UNSUCCESSFULLY saved features!');
    end
    else
        msgbox('GND reference is already saved!');
    end
