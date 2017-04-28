function varargout = EEGNIRS_GUI(varargin)
% EEGNIRS_GUI MATLAB code for EEGNIRS_GUI.fig
%      EEGNIRS_GUI, by itself, creates a new EEGNIRS_GUI or raises the existing
%      singleton*.
%
%      H = EEGNIRS_GUI returns the handle to a new EEGNIRS_GUI or the handle to
%      the existing singleton*.
%
%      EEGNIRS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EEGNIRS_GUI.M with the given input arguments.
%
%      EEGNIRS_GUI('Property','Value',...) creates a new EEGNIRS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EEGNIRS_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EEGNIRS_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EEGNIRS_GUI

% Last Modified by GUIDE v2.5 07-Jun-2015 16:46:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @EEGNIRS_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @EEGNIRS_GUI_OutputFcn, ...
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


% --- Executes just before EEGNIRS_GUI is made visible.
function EEGNIRS_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EEGNIRS_GUI (see VARARGIN)

% Choose default command line output for EEGNIRS_GUI
handles.output = hObject;

% ADD function path
addpath([pwd '\functions']);

% default GUI values: CHANGE HERE DEFAULT GUI VALUES
set(handles.slider_low_eeg,'Value',0.5);
set(handles.lowfreqeegval,'String',num2str(get(handles.slider_low_eeg,'Value')));

set(handles.slider_high_eeg,'Value',80);
set(handles.hifreqeegval,'String',num2str(get(handles.slider_high_eeg,'Value')));

set(handles.slider_low_nirs,'Value',0.01);
set(handles.lowfreqnirsval,'String',num2str(get(handles.slider_low_nirs,'Value')));

set(handles.slider_high_nirs,'Value',0.5);
set(handles.hifreqnirsval,'String',num2str(get(handles.slider_high_nirs,'Value')));

set(handles.signal_selection,'Value',4);

set(handles.notcheeg,'Value',1);
set(handles.notchnirs,'Value',2);

%load default filter values
handles.eegfiltLval = get(handles.slider_low_eeg,'Value');
handles.eegfiltHval = get(handles.slider_high_eeg,'Value');
handles.nirsfiltLval = get(handles.slider_low_nirs,'Value');
handles.nirsfiltHval = get(handles.slider_high_nirs,'Value');
handles.filteegN = 60;
handles.filtnirsN = 0;

%declare datacells
handles.datacells = cell( 1,  2 );

%set default plot state
%none -> no plot.
set(handles.select_view,'Value',1);
handles.prep_select_view = get(handles.select_view,'String');

% Select View STRINGS
% raw NIRS select view
str = cell(1,3);
str{1} = 'NONE';
str{2} = 'WL1';
str{3} = 'WL2';
handles.raw_select_view = str;
% spectra select view
str = cell(1,3);
str{1} = 'NONE';
str{2} = 'HbO';
str{3} = 'HbR';
handles.nirs_power_select_view = str;

% set subplot axes invisible
set(handles.sub1,'Visible','off');
set(handles.sub2,'Visible','off');

% save default viewpanel dimension and position
handles.viewpanelpos = get(handles.viewpanel,'Position');

% declare global variable to hide plots
handles.hideplot = [];
% 1: one axis, 2:two axis, 3: 2 subplots
handles.typeplot = 0;

%set default configuration
str = get(handles.nirsID,'String');
handles.configID = str2double(str{3});

%set save check box visibility off
set(handles.save_check,'Visible','off');

%set load raw data flag to 0
handles.loadrawdata = 0;
%set preprocess data flag to 0
handles.preprocessdata = 0;


% handles variables to plot and filter eeg
handles.xxeegRaw = [];
handles.xxeegf = [];
handles.xxeeg = [];

handles.BP = 0;

% set eeg filter panel invisible
set(handles.eegfilter,'Visible','off');


%% CREATE EEGNIRS_LAB directory

current = pwd;

if ~isempty(strfind(current,'Documents'))
    k = strfind(current,'Documents');
    current = [current(1:k-1) 'Documents'];
    if ~exist([current '\EEGNIRS_LAB'],'dir')
        mkdir(current,'EEGNIRS_LAB');
    end
elseif ~isempty(strfind(current,'Documenti'))
    k = strfind(current,'Documenti');
    current = [current(1:k-1) 'Documenti'];
    if ~exist([current '\EEGNIRS_LAB'],'dir')
        mkdir([current '\EEGNIRS_LAB'],'EEGNIRS_LAB');
    end
    %%add more for different languages and OS
end

%update eegnirs directory path
handles.dir = [current '\EEGNIRS_LAB'];

%% CREATE, if non existing, configuration folder and write there default config
% files
if ~exist([handles.dir '\NIRS Configuration'],'dir')
    mkdir(handles.dir,'NIRS Configuration');
    handles.configdir = [handles.dir '\NIRS Configuration'];
    
    %write default config (1,2,3)
    
    % Configuration_1
    
    fid = fopen([handles.configdir '\Configuration_1.txt'],'wt');
    
    if fid~=-1
        
        fprintf(fid, '%s', 'Configuration 1');
        fprintf(fid,'\n');
        fprintf(fid, '%s\t%s\t%s', 'S','D','Label');
        fprintf(fid,'\n');
        
        src_det = [1    1;2 2;3 3;4 4;5 5;6	6;7	7;8	8;9	9;10	10;11	11;12	12;17	17; ...
            18	18;19	19;20	20;21   21;22   22;23   23;24   24];
        
        lab = cell(20,1);
        
        lab{1} = 'Fp1';
        lab{2} = 'F3';
        lab{3} = 'F7';
        lab{4} = 'T7';
        lab{5} = 'C3';
        lab{6} = 'P7';
        lab{7} = 'O1';
        lab{8} = 'O2';
        lab{9} = 'P3';
        lab{10} = 'Pz';
        lab{11} = 'P4';
        lab{12} = 'SS';
        lab{13} = 'FP2';
        lab{14} = 'Fz';
        lab{15} = 'Cz';
        lab{16}	= 'F4';
        lab{17} = 'F8';
        lab{18} = 'T8';
        lab{19} = 'C4';
        lab{20} = 'T6';
        
        format = '%d\t%d\t%s';
        
        for i = 1:length(lab)
            
            fprintf(fid,format, src_det(i,1),src_det(i,2),lab{i});
            fprintf(fid,'\n');
            
        end
        fclose(fid);
    end
    
    % Configuration_2 --> SUPER channels
    
    fid = fopen([handles.configdir '\Configuration_2.txt'],'wt');
    
    if fid~=-1
        
        fprintf(fid, '%s', 'Configuration 2');
        fprintf(fid,'\n');
        fprintf(fid, '%s\t%s\t%s', 'S','D','Label');
        fprintf(fid,'\n');
        
        src_det = [1    1;2 2;3 3;4 4;5 5;6	6;7	7;8	8;9	9;10	10;11	11;12	12;13	13;14	14;15	15;16	16;17	17; ...
            18	18;19	19;20	20;3	23;4	24;11	21;12	22];
        
        lab = cell(24,1);
        
        lab{1} = 'Fp1';
        lab{2} = 'F3';
        lab{3} = 'F7';
        lab{4} = 'T3';
        lab{5} = 'C3';
        lab{6} = 'T5';
        lab{7} = 'P3';
        lab{8} = 'O1';
        lab{9} = 'Fp2';
        lab{10} = 'F4';
        lab{11} = 'F8';
        lab{12} = 'T4';
        lab{13} = 'C4';
        lab{14} = 'T6';
        lab{15} = 'P4';
        lab{16} = '02';
        lab{17} = 'Fz';
        lab{18} = 'Cz';
        lab{19} = 'Pz';
        lab{20}	= 'SSuper';
        lab{21} = 'F7Super';
        lab{22} = 'T3Super';
        lab{23} = 'F8Super';
        lab{24} = 'T4Super';
        
        format = '%d\t%d\t%s';
        
        for i = 1:length(lab)
            
            fprintf(fid,format, src_det(i,1),src_det(i,2),lab{i});
            fprintf(fid,'\n');
            
        end
        fclose(fid);
        
    end
    
    % Configuration_3 --> ALE motor tasks
    
    fid = fopen([handles.configdir '\Configuration_3.txt'],'wt');
    
    if fid~=-1
        
        fprintf(fid, '%s', 'Configuration 3');
        fprintf(fid,'\n');
        fprintf(fid, '%s\t%s\t%s', 'S','D','Label');
        fprintf(fid,'\n');
        
        src_det = [1    1;1 2;2 1;1 3;3 2;2	3;3	3;2	4;4	3;3	5;4	4;4	5;5 4;4	6;6	5;5	6;6	6; ...
            7    7;7 8;8 7;7 9;9 8;8	9;9	9;8	10;10	9;9	11;10	10;10	11;11	10;10	12;12	11;11	12;12	12];
        
        lab = cell(34,1);
        
        lab{1} = 'Fc1A';
        lab{2} = 'Fc3A';
        lab{3} = 'Fc1M';
        lab{4} = 'Fc3M';
        lab{5} = 'Fc3L';
        lab{6} = 'C1A';
        lab{7} = 'C3A';
        lab{8} = 'C1M';
        lab{9} = 'C3M';
        lab{10} = 'C3L';
        lab{11} = 'Cp1A';
        lab{12} = 'Cp3A';
        lab{13} = 'Cp1M';
        lab{14} = 'Cp3M';
        lab{15} = 'Cp3L';
        lab{16} = 'Cp1P';
        lab{17} = 'Cp3P';
        lab{18} = 'Fc2A';
        lab{19} = 'Fc4A';
        lab{20}	= 'Fc2M';
        lab{21} = 'Fc4M';
        lab{22} = 'Fc4L';
        lab{23} = 'C2A';
        lab{24} = 'C4A';
        lab{25} = 'C2M';
        lab{26} = 'C4M';
        lab{27} = 'C4L';
        lab{28} = 'Cp2A';
        lab{29} = 'Cp4A';
        lab{30} = 'Cp2M';
        lab{31} = 'Cp4M';
        lab{32} = 'Cp4L';
        lab{33} = 'Cp2P';
        lab{34} = 'Cp4P';
        
        format = '%d\t%d\t%s';
        
        for i = 1:length(lab)
            
            fprintf(fid,format, src_det(i,1),src_det(i,2),lab{i});
            fprintf(fid,'\n');
            
        end
        fclose(fid);
        
    end
    
    strnirs = { '1' ; '2' ; '3'};
    
else
    % count how many configuration files are there in NIRS configuration
    % folder
    handles.configdir = [handles.dir '\NIRS Configuration'];
    
    count = 0;
    strnirs = [];
    
    while (exist([handles.configdir '\Configuration_' num2str(count + 1) '.txt'],'file'))
        count = count + 1;
        strnirs{count} = num2str(count);
    end
    
end

handles.configdir = [handles.dir '\NIRS Configuration'];


%% CREATE, if non existing, eeg channels folder and write there default eeg config
% files

if ~exist([handles.dir '\EEG Configuration'],'dir')
    mkdir(handles.dir,'EEG Configuration');
    handles.eegconfigdir = [handles.dir '\EEG Configuration'];
    
    %write default config (1,2)
    
    % Configuration_1
    
    fid = fopen([handles.eegconfigdir '\Configuration_1.txt'],'wt');
    
    if fid~=-1
        
        fprintf(fid, '%s', 'Motor Cortex Configuration');
        fprintf(fid,'\n');
        fprintf(fid, '%s', 'EEG Channel Label');
        fprintf(fid,'\n');
        
        lab = cell(17,1);
        
        lab{1} = 'F3';
        lab{2} = 'Fz';
        lab{3} = 'F4';
        lab{4} = 'Fc5';
        lab{5} = 'Fc1';
        lab{6} = 'Fc2';
        lab{7} = 'Fc6';
        lab{8} = 'C3';
        lab{9} = 'Cz';
        lab{10} = 'C4';
        lab{11} = 'Cp5';
        lab{12} = 'Cp1';
        lab{13} = 'Cp2';
        lab{14} = 'Cp6';
        lab{15} = 'P3';
        lab{16}	= 'Pz';
        lab{17} = 'P4';
        
        format = '%s';
        
        for i = 1:length(lab)
            fprintf(fid,format,lab{i});
            fprintf(fid,'\n');
        end
        fclose(fid);
    end
    
    % Configuration_2
    
    fid = fopen([handles.eegconfigdir '\Configuration_2.txt'],'wt');
    
    if fid~=-1
        
        fprintf(fid, '%s', 'Standard Configuration');
        fprintf(fid,'\n');
        fprintf(fid, '%s', 'EEG Channel Label');
        fprintf(fid,'\n');
        
        lab = cell(19,1);
        
        lab{1} = 'Fp1';
        lab{2} = 'Fp2';
        lab{3} = 'F4';
        lab{4} = 'C4';
        lab{5} = 'P4';
        lab{6} = 'O2';
        lab{7} = 'F8';
        lab{8} = 'T4';
        lab{9} = 'T6';
        lab{10} = 'Fz';
        lab{11} = 'Cz';
        lab{12} = 'Pz';
        lab{13} = 'F3';
        lab{14} = 'C3';
        lab{15} = 'P3';
        lab{16}	= 'O1';
        lab{17} = 'F7';
        lab{18} = 'T3';
        lab{19} = 'T5';
        
        format = '%s';
        
        for i = 1:length(lab)
            fprintf(fid,format,lab{i});
            fprintf(fid,'\n');
        end
        fclose(fid);
        
    end
    
    streeg = { '1' ; '2' };
    
else
    % count how many configuration files are there in NIRS configuration
    % folder
    handles.eegconfigdir = [handles.dir '\EEG Configuration'];
    
    count = 0;
    streeg = [];
    
    while (exist([handles.eegconfigdir '\Configuration_' num2str(count + 1) '.txt'],'file'))
        count = count + 1;
        streeg{count} = num2str(count);
    end
    
end

%set default NIRS configuration IDs
set(handles.nirsID,'String',strnirs);
set(handles.nirsID,'Value',3);

%set default EEG configuration IDs
set(handles.eegID,'String',streeg);
set(handles.eegID,'Value',1);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EEGNIRS_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EEGNIRS_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
varargout{1} = handles.datacells;


% --- Executes on button press in load_raw.
function load_raw_Callback(hObject, eventdata, handles)
% hObject    handle to load_raw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileeeg, patheeg] = uigetfile('*.edf','Select EEG raw data file (.edf)','C:\Users\alessio\Google Drive\UHshared\NIRS_EEG\Data\Ale_raw_data\');

if fileeeg ~=0
    
    fullpatheeg = [patheeg fileeeg];
    handles.fileeeg = fullpatheeg;
    [filenirs, pathnirs] = uigetfile...
        ('*.wl1','Select fNIRS raw data file (.wl1)','C:\Users\alessio\Google Drive\UHshared\NIRS_EEG\Data\Ale_raw_data\');
    
    if filenirs ~=0
        
        % find substring withouth extension for nirs processing
        fullpathnirs = [pathnirs filenirs];
        l = length(fullpathnirs);
        fullpathnirs = fullpathnirs(1:l-4);
        handles.filenirs = fullpathnirs;
        
        
        
        %flag to allow preprocess of the data
        handles.loadrawdata = 1;
        handles.preprocessdata = 0;
        
        % clean axes pointer
        handles.hideplot = [];
        
        % update file path text
        set(handles.filepaths,'String', ['EEG: ' fileeeg  char(10) 'fNIRS: ' filenirs]);
        
        %hide save check box
        set(handles.save_check,'Value',0);
        set(handles.save_check,'Visible','off');
        
        % Update handles structure
        guidata(hObject, handles);
        
        msgbox('Raw data loaded CORRECTLY!');
        
    else
        msgbox('Raw data NOT loaded CORRECTLY!');
    end
    
else
    msgbox('Raw data NOT loaded CORRECTLY!');
    
end

% --- Executes on button press in load_preproc.
function load_preproc_Callback(hObject, eventdata, handles)
% hObject    handle to load_preproc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~exist([handles.dir '\Pre Processed Data'],'dir')
    mkdir(handles.dir,'Pre Processed Data');
end

[filepreproc, path] = uigetfile('*.mat','Select pre processed data MATLAB file',[handles.dir '\Pre Processed Data']);

if filepreproc ~= 0
    
    tmp = load([path filepreproc]);
    handles.datacells = tmp.data;
    
    handles.preprocessedfile = filepreproc;
    
    %flag to allow preprocess of the data
    handles.preprocessdata = 1;
    
    handles.xxeeg = handles.datacells{1,2}{1,1};
    handles.xxeegRaw = handles.datacells{1,2}{1,14};
    
    % clean axes pointer
    handles.hideplot = [];
    
    % update file path text
    set(handles.filepaths,'String', ['Pre processed: ' filepreproc ]);
    %save check box
    set(handles.save_check,'Visible','on');
    set(handles.save_check,'Value',1);
    
    msgbox('Pre processed data loaded CORRECTLY!');
    
else
    msgbox('Pre processed data NOT loaded CORRECTLY!');
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in preprocess_button.
function preprocess_button_Callback(hObject, eventdata, handles)
% hObject    handle to preprocess_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.loadrawdata
    
    %     if handles.preprocessdata~=1
    
    %waitbar
    h = waitbar(0,'Please wait...');
    
    % GET AND PRE PROCESS DATA
    
    % EEG
    
    %clip 4 seconds at beginning and end
    dtbeg = 4;
    dtend = 4;
    % EEG  Bandpas filter and notch 
    filterlofreq = handles.eegfiltLval;
    filterhifreq = handles.eegfiltHval; %
    filternotch = handles.filteegN;
    file1 = handles.fileeeg;
    [eegchanlabels, eegchanlocs,xxeeg, xxeegf, evteeg, Fs] = getpreprocess_eegALE( file1, dtbeg, dtend, filterlofreq, filterhifreq, filternotch );
        
    % extract EEG array based on labels
    configeeg = [handles.eegconfigdir '\Configuration_' num2str(get(handles.eegID,'Value')) '.txt'];    
    [eeg_chan_array, eegchanlabels] = extract_eeg_array(configeeg, eegchanlabels);
    
    xxeegf = xxeegf(:,eeg_chan_array);
    xxeeg = xxeeg(:,eeg_chan_array);
    
    numchaneeg = length(eeg_chan_array);
    
    waitbar(0.25);
    
    
    % Compute eeg power spectra
    
    Feeg = [];
    nweegpower = Fs;
    
    for i=1:numchaneeg
        [SS,Feeg,Tim,P] = spectrogram(xxeegf(:,i),nweegpower,0,Feeg,Fs);
        if i == 1
            Peeg = zeros(length(Feeg),numchaneeg);
        end
        Peeg(:,i) = 10.*log10(mean(P,2));
        waitbar(0.25 + 0.25/numchaneeg*i);
    end
    
    
    waitbar(0.5);
    
    % NIRS
    
    % Bandpass filter and notch 
    filterlofreq = handles.nirsfiltLval;
    filterhifreq = handles.nirsfiltHval; %
    filternotch = handles.filtnirsN;
    % GET sample rate
    display(handles.filenirs);
    Fsnirs = str2double( getfieldval( [handles.filenirs '_config.txt'], 'SamplingRate', '=') );
    %  EVENTS filename
    filenirsevt = [ handles.filenirs  '.evt' ];
    % LABELS and SRC_DET PAIR from Configuration FILE
    filename = [handles.configdir '\Configuration_' num2str(handles.configID) '.txt'];
    [src_det_pair, nirschanlabels] = extract_nirs_labels(filename);
    numchan = length(nirschanlabels);
    numsrc = str2double( getfieldval( [handles.filenirs '_config.txt'], 'source_N', '=') );
    numdet = str2double( getfieldval( [handles.filenirs '_config.txt'], 'detector_N', '=') );
    
    % LOAD DATA
    mydat = load( [ handles.filenirs '.wl1' ]  );
    xx = fnnirs_reduce_data_ALE( src_det_pair, nirschanlabels, numsrc, numdet, mydat );
    xx_wl1 = xx;
    mydat = load( [ handles.filenirs '.wl2' ]  );
    xx = fnnirs_reduce_data_ALE( src_det_pair, nirschanlabels, numsrc, numdet, mydat );
    xx_wl2 = xx;
    
    
    % convert to HbO, HbR
    [xxo, xxr] = mybeerlambert2( xx_wl1, xx_wl2 );
    
    %compute nirs power spectra
    
    Fn = linspace(0,Fsnirs/2,100);
    res = length(Fn);
    
    Pxxo = zeros(res,numchan);
    Pxxr = zeros(res,numchan);
    
    if (numchan ~= size(xxo,2))
       error('Inconsistent INFO: NIRS number of channels in config file is wrong') 
       close (h);
    end
    
    for i=1:numchan
        [So,Fn,To,Po] = spectrogram(xxo(:,i),hamming(512),256,Fn,Fsnirs);
        [Sr,Fn,Tr,Pr] = spectrogram(xxr(:,i),hamming(512),256,Fn,Fsnirs);
        Pxxo(:,i) = 10.*log10(mean(Po,2));
        Pxxr(:,i) = 10.*log10(mean(Pr,2));
        waitbar(0.5 + 0.25/numchan*i);
    end
    
    
    % preprocess
    xxnirsf = mybutterandiirnotchfiltersALE( xxo, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs, 'nirs' );
    xxnirsfR = mybutterandiirnotchfiltersALE( xxr, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs, 'nirs' );
    
    waitbar(0.75);
    
    %% SYNC
    
    % EEG event sample numbers
    
    Vevthresh = 1000;
    %finds the samples above threshold 
    nevteeg = remove_contiguous(  find( Vevthresh<evteeg )  );
    tteeg = (1:length(evteeg))./Fs;
    tevteeg = nevteeg./Fs;
    
    %Use first event to sync
    nbegeeg = nevteeg(1);
    
    % NIRS event sample numbers
    myevtdat = load( filenirsevt );
    nevtnirs = zeros(size(myevtdat,1),2);
    nevtnirs(:,1) = myevtdat(:,1);
    % Save the type of trigger assiociated to an event
    for i = 1:size(myevtdat,1)
        ind = find(myevtdat(i,:)==1);
        %if 2 events occur at the same time, it means the trigger was 255
        if ~isempty(ind) && length(ind)==1
            nevtnirs(i,2) = ind - 1;
        elseif ~isempty(ind) && length(ind)>1
            nevtnirs(i,2) = 10;
        end
    end
    
    
    ttnirs = (1:length(xxnirsf))./Fsnirs;
    tevtnirs = [nevtnirs(:,1)./Fsnirs nevtnirs(:,2)];
    nbegnirs = nevtnirs(1,1);
    
    % clip eeg data at first event
    xxeegfclip = xxeegf( nbegeeg:end, : );
    xxeegclip = xxeeg( nbegeeg:end, : );
    tteegclip = tteeg( nbegeeg:end );
    tteegclip = tteegclip - tteegclip(1);
    
    % clip nirs data at first event
    xxnirsfclip = xxnirsf(nbegnirs:end, :);
    xxnirsfRclip = xxnirsfR(nbegnirs:end, :);
    
    xx_wl1 = xx_wl1(nbegnirs:end, :);
    xx_wl2 = xx_wl2(nbegnirs:end, :);
    
    xxo = xxo(nbegnirs:end, :);
    xxr = xxr(nbegnirs:end, :);
    
    ttnirsclip = ttnirs( nbegnirs:end );
    ttnirsclip = ttnirsclip - ttnirsclip(1);
    
    % clip the tails --> ttnirsclip and tteegclip have same lenght
    if ttnirsclip(end)>tteegclip(end),
        % since nirs is longer, clip end of nirs
        nnn = min(find( tteegclip(end)<ttnirsclip,1)) - 1;
        if ~isempty(nnn)
            ttnirsclip = ttnirsclip(1:nnn);
            xxnirsfclip = xxnirsfclip(1:nnn, :);
            xxnirsfRclip = xxnirsfRclip(1:nnn, :);
            xx_wl1 = xx_wl1(1:nnn, :);
            xx_wl2 = xx_wl2(1:nnn, :);
            xxo = xxo(1:nnn, :);
            xxr = xxr(1:nnn, :);
        end
    else
        % the other way around
        nnn = min(find( ttnirsclip(end)<tteegclip,1) ) - 1;
        if ~isempty(nnn)
            tteegclip = tteegclip(1:nnn);
            xxeegfclip = xxeegfclip(1:nnn, :);
            xxeegclip = xxeegclip(1:nnn, :);
        end
    end
    
    %cut at the last floor second
    lastsec = floor(ttnirsclip(end));
    %cut EEG
    nnn = min( find( lastsec < tteegclip,1) ) - 1;
    if ~isempty(nnn)
        tteegclip = tteegclip(1:nnn);
        xxeegfclip = xxeegfclip(1:nnn,:);
        xxeegclip = xxeegclip(1:nnn, :);
    end
    %cut NIRS
    mmm = min( find( lastsec < ttnirsclip,1) ) - 1;
    if ~isempty(mmm)
        ttnirsclip = ttnirsclip(1:mmm);
        xxnirsfclip = xxnirsfclip(1:mmm,:);
        xxnirsfRclip = xxnirsfRclip(1:mmm, :);
        xx_wl1 = xx_wl1(1:mmm, :);
        xx_wl2 = xx_wl2(1:mmm, :);
        xxo = xxo(1:mmm, :);
        xxr = xxr(1:mmm, :);
    end
    
    waitbar(0.90);
    
    %% assign info data
    
    infocells = cell(2,1);
    infocellseeg = cell(4,1);
    infocellsnirs = cell(4,1);
    
    infocellseeg{1} = handles.fileeeg;
    infocellseeg{2} = Fs;
    infocellseeg{3} = numchaneeg;
    infocellseeg{4} = eegchanlabels;
    
    infocellsnirs{1} = handles.filenirs;
    infocellsnirs{2} = Fsnirs;
    infocellsnirs{3} = numchan;
    infocellsnirs{4} = nirschanlabels;
    
    infocells{1} = infocellseeg;
    infocells{2} = infocellsnirs;
    
    synceddata = cell(1,16);
    synceddata{1} = xxeegfclip;
    synceddata{2} = tteegclip;
    synceddata{3} = xxnirsfclip;
    synceddata{4} = xxnirsfRclip;
    synceddata{5} = ttnirsclip;
    synceddata{6} = [tevtnirs(:,1) - tevtnirs(1,1) , tevtnirs(:,2)];
    synceddata{7} = xx_wl1;
    synceddata{8} = xx_wl2;
    synceddata{9} = Pxxo;
    synceddata{10} = Pxxr;
    synceddata{11} = Fn;
    synceddata{12} = Peeg;
    synceddata{13} = Feeg;
    synceddata{14} = xxeegclip;
    synceddata{15} = xxo;
    synceddata{16} = xxr;
    
    waitbar(1);
    
    %% assign return data
    % info
    handles.datacells{1} = infocells;
    % synced data
    handles.datacells{2} = synceddata;
    
    handles.preprocessdata = 1;
    
    handles.xxeeg = xxeegfclip;
    handles.xxeegRaw = xxeegclip;
    
    % update file path text and save check box
    set(handles.filepaths,'String', ['Pre processed: ' 'No name' ]);
    set(handles.save_check,'Value',0);
    set(handles.save_check,'Visible','on');
    
    handles.axes_hide = [];
    
    guidata(hObject, handles);
    
    close(h);
    
    msgbox('Raw data pre processed CORRECTLY!');
    %     else
    %         msgbox('WARNING! Data have already been pre processed');
    %
    %     end
else
    msgbox('WARNING! Load raw data for pre processing!');
end

% --- Executes on selection change in notchnirs.
function notchnirs_Callback(hObject, eventdata, handles)
% hObject    handle to notchnirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns notchnirs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from notchnirs

val = get(hObject,'Value');
str = get(hObject,'String');
switch str{val}
    case 'YES (0.1 Hz)'
        handles.filtnirsN = 0.1;
    case 'NO'
        handles.filtnirsN = 0;
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function notchnirs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to notchnirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in notcheeg.
function notcheeg_Callback(hObject, eventdata, handles)
% hObject    handle to notcheeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns notcheeg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from notcheeg

val = get(hObject,'Value');
str = get(hObject,'String');
switch str{val}
    case 'YES (60 Hz)'
        handles.filteegN = 60;
    case 'YES (50 Hz)'
        handles.filteegN = 50;
    case 'NO'
        handles.filteegN = 0;
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function notcheeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to notcheeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function slider_low_eeg_Callback(hObject, eventdata, handles)
% hObject    handle to slider_low_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.eegfiltLval = get(handles.slider_low_eeg,'Value');
set(handles.lowfreqeegval,'String',num2str(get(handles.slider_low_eeg,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_low_eeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_low_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider_high_eeg_Callback(hObject, eventdata, handles)
% hObject    handle to slider_high_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.eegfiltHval = get(handles.slider_high_eeg,'Value');
set(handles.hifreqeegval,'String',num2str(get(handles.slider_high_eeg,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_high_eeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_high_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider_low_nirs_Callback(hObject, eventdata, handles)
% hObject    handle to slider_low_nirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.nirsfiltLval = get(handles.slider_low_nirs,'Value');
set(handles.lowfreqnirsval,'String',num2str(get(handles.slider_low_nirs,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_low_nirs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_low_nirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider_high_nirs_Callback(hObject, eventdata, handles)
% hObject    handle to slider_high_nirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.nirsfiltHval = get(handles.slider_high_nirs,'Value');
set(handles.hifreqnirsval,'String',num2str(get(handles.slider_high_nirs,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_high_nirs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_high_nirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function filepaths_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filepaths (see GCBO)
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
set(handles.filter_button,'Value',0);
set(handles.band_power_button,'Value',0);
handles.BP = 0;

if val_sig == 4 && (val==2 || val == 5 || val == 6 || val == 8)
    set(handles.eegfilter,'Visible','on');
else
    set(handles.eegfilter,'Visible','off');
end

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
                set(handles.channel2,'Visible', 'off');
                set(handles.channel1,'Value', 1);
                handles.typeplot=1;
                
            case 'WL1'
                if handles.preprocessdata
                    
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
                
            case 'WL2'
                if handles.preprocessdata
                    
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
                
            case 'HbO'
                if handles.preprocessdata
                    
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
                
            case 'HbR'
                if handles.preprocessdata
                    
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
        
    case 3
        
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
        
    case 4
        
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
            case 'WL1'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                ttnirs = handles.datacells{1,2}{1,5};
                xxnirs = handles.datacells{1,2}{1,7};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ttnirs,xxnirs(:,val1 - 1),'r','linewidth',2), xlabel('time [s]'), ylabel('Wavelength 1 [mW]'),title(['RAW Wavelength 1: ' nirslabels{val1-1}]);
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
            case 'WL2'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                ttnirs = handles.datacells{1,2}{1,5};
                xxnirsR = handles.datacells{1,2}{1,8};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,ttnirs,xxnirsR(:,val1 - 1),'b','linewidth',2), xlabel('time [s]'), ylabel('Wavelength 2 [mW]'),title(['RAW Wavelength 1: ' nirslabels{val1-1}]);
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
            case 'HbO'
                val1 = get(handles.channel1,'Value');
                nirslabels = handles.datacells{1,1}{2,1}{4,1};
                F = handles.datacells{1,2}{1,11};
                Pxxo = handles.datacells{1,2}{1,9};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,F,Pxxo(:,val1 - 1),'r','linewidth',2), xlabel('frequency [Hz]'), ylabel('HbO [mM/L]^2/Hz'),title(['SPECTRUM HbO: ' nirslabels{val1-1}]);
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
                F = handles.datacells{1,2}{1,11};
                Pxxr = handles.datacells{1,2}{1,10};
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    plot(handles.viewplot,F,Pxxr(:,val1 - 1),'b','linewidth',2), xlabel('frequency [Hz]'), ylabel('HbR [mM/L]^2/Hz'),title(['Spectrum HbR: ' nirslabels{val1-1}]);
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
        
    case 3
        
        val1 = get(handles.channel1,'Value');
        eeglabels = handles.datacells{1,1}{1,1}{4,1};
        Feeg = handles.datacells{1,2}{1,13};
        Peeg = handles.datacells{1,2}{1,12};
        set(handles.viewplot,'Visible','on');
        
        if val1 ~= 1
            axes(handles.viewplot);
            plot(handles.viewplot,Feeg,Peeg(:,val1 - 1),'Color',[0 0.5 0]), xlabel('frequency [Hz]'), ylabel('EEG [uV^2/Hz]'), title(['EEG: ' eeglabels{val1-1}]);
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
        
        
    case 4
        
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
                xxeeg = handles.xxeeg;
                set(handles.viewplot,'Visible','on');
                
                if val1 ~= 1
                    axes(handles.viewplot);
                    if handles.BP
                        plot(handles.viewplot,tteeg,xxeeg(:,val1 - 1),'Color',[0 0.5 0],'linewidth',2), xlabel('time [s]'), ylabel('EEG [uV]'), title(['EEG: ' eeglabels{val1-1}]);
                        grid(handles.viewplot,'on');
                    else
                        plot(handles.viewplot,tteeg,xxeeg(:,val1 - 1),'Color',[0 0.5 0],'linewidth',1), xlabel('time [s]'), ylabel('EEG [uV]'), title(['EEG: ' eeglabels{val1-1}]);
                        grid(handles.viewplot,'on');
                    end
                    
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
                    xxeeg = handles.xxeeg;
                    
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
                        if handles.BP
                            set(hLine1,'linewidth',2);
                        end
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
                    xxeeg = handles.xxeeg;
                    
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
                        if handles.BP
                            set(hLine1,'linewidth',2);
                        end
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
                    xxeeg = handles.xxeeg;
                    
                    val2 = get(handles.channel2,'Value');
                    set(handles.viewplot,'Visible','on');
                    
                    if val1 ~= 1 && val2 ~= 1
                        axes(handles.viewplot);
                        if handles.BP
                            plot(handles.viewplot,tteeg,xxeeg(:,val1 - 1),'Color',[0 0.5 0],'Linewidth',2)
                        else
                            plot(handles.viewplot,tteeg,xxeeg(:,val1 - 1),'Color',[0 0.5 0])
                        end
                        hold on
                        if handles.BP
                            plot(handles.viewplot,tteeg,xxeeg(:,val2 - 1),'Color',[0 0.2 0],'Linewidth',2)
                        else
                            plot(handles.viewplot,tteeg,xxeeg(:,val2 - 1),'Color',[0 0.2 0])
                        end
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

val = get(handles.select_view,'Value');
str = get(handles.select_view,'String');

%remove events
set(handles.eventsbutton,'Value',0);

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
    case 'EEG + HbO'
        if get(handles.channel1,'Value')~=1
            val1 = get(handles.channel1,'Value');
            eeglabels = handles.datacells{1,1}{1,1}{4,1};
            tteeg = handles.datacells{1,2}{1,2};
            xxeeg = handles.xxeeg;
            
            val2 = get(hObject,'Value');
            nirslabels = handles.datacells{1,1}{2,1}{4,1};
            ttnirs = handles.datacells{1,2}{1,5};
            xxnirs = handles.datacells{1,2}{1,3};
            set(handles.viewplot,'Visible','on');
            
            if val1 ~= 1 && val2 ~=1
                axes(handles.viewplot);
                [hAx,hLine1,hLine2] = plotyy(handles.viewplot,tteeg,xxeeg(:,val1 - 1),ttnirs,xxnirs(:,val2 - 1),'plot');
                set(hLine1,'color',[0 0.5 0]);
                set(hLine2,'color','red');
                set(hAx,{'ycolor'},{[0 0.5 0];'r'}) ;
                ylabel(hAx(1),'EEG [uV]') % left y-axis
                ylabel(hAx(2),'HbO [mM/L]') % right y-axis
                if handles.BP
                    set(hLine1,'linewidth',2);
                end
                legend('EEG','HbO');
                xlabel('time [s]'), title(['EEG: ' eeglabels{val1-1} ' HbO: ' nirslabels{val2-1}]), grid on;
                
                handles.hideplot = hAx;
                handles.typeplot=2;
            end
            guidata(hObject, handles);
        end
    case 'EEG + HbR'
        if get(handles.channel1,'Value')~=1
            val1 = get(handles.channel1,'Value');
            eeglabels = handles.datacells{1,1}{1,1}{4,1};
            tteeg = handles.datacells{1,2}{1,2};
            xxeeg = handles.xxeeg;
            
            val2 = get(hObject,'Value');
            nirslabels = handles.datacells{1,1}{2,1}{4,1};
            ttnirs = handles.datacells{1,2}{1,5};
            xxnirsR = handles.datacells{1,2}{1,4};
            set(handles.viewplot,'Visible','on');
            
            if val1 ~= 1 && val2 ~=1
                axes(handles.viewplot);
                [hAx,hLine1,hLine2] = plotyy(handles.viewplot,tteeg,xxeeg(:,val1 - 1),ttnirs,xxnirsR(:,val2 - 1),'plot');
                set(hLine1,'color',[0 0.5 0]);
                set(hLine2,'color','blue');
                set(hAx,{'ycolor'},{[0 0.5 0];'b'}) ;
                ylabel(hAx(1),'EEG [uV]') % left y-axis
                ylabel(hAx(2),'HbR [mM/L]') % right y-axis
                if handles.BP
                    set(hLine1,'linewidth',2);
                end
                legend('EEG','HbR');
                xlabel('time [s]'), title(['EEG: ' eeglabels{val1-1} ' HbR: ' nirslabels{val2-1}]), grid on;
                
                handles.hideplot = hAx;
                handles.typeplot=2;
            end
            guidata(hObject, handles);
        end
    case 'HbO + HbR'
        if get(handles.channel1,'Value')~=1
            val1 = get(handles.channel1,'Value');
            nirslabels = handles.datacells{1,1}{2,1}{4,1};
            ttnirs = handles.datacells{1,2}{1,5};
            xxnirs = handles.datacells{1,2}{1,3};
            
            val2 = get(hObject,'Value');
            xxnirsR = handles.datacells{1,2}{1,4};
            set(handles.viewplot,'Visible','off');
            
            if val1 ~= 1 && val2 ~=1
                set(handles.sub1,'Visible','on');
                axes(handles.sub1)
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
        if get(handles.channel1,'Value')~=1
            val1 = get(handles.channel1,'Value');
            eeglabels = handles.datacells{1,1}{1,1}{4,1};
            tteeg = handles.datacells{1,2}{1,2};
            xxeeg = handles.xxeeg;
            
            val2 = get(hObject,'Value');
            set(handles.viewplot,'Visible','on');
            
            if val1 ~= 1 && val2 ~=1
                axes(handles.viewplot);
                if handles.BP
                    plot(handles.viewplot,tteeg,xxeeg(:,val1 - 1),'Color',[0 0.5 0],'linewidth',2)
                else
                    plot(handles.viewplot,tteeg,xxeeg(:,val1 - 1),'Color',[0 0.5 0])
                end
                hold on
                if handles.BP
                    plot(handles.viewplot,tteeg,xxeeg(:,val2 - 1),'Color',[0 0.2 0],'linewidth',2)
                else
                    plot(handles.viewplot,tteeg,xxeeg(:,val2 - 1),'Color',[0 025 0])
                end
                hold off
                
                legend('EEG-1','EEG-2');
                xlabel('time [s]'), ylabel('EEG [uV]'), title(['EEG-1: ' eeglabels{val1-1} ' EEG-2: ' eeglabels{val2-1}]), grid on;
                
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
            
            val2 = get(hObject,'Value');
            set(handles.viewplot,'Visible','on');
            
            if val1 ~= 1 && val2 ~=1
                axes(handles.viewplot);
                plot(handles.viewplot,ttnirs,xxnirs(:,val1 - 1),'r');
                hold on
                plot(ttnirs,xxnirs(:,val2 - 1),'Color',[0.5,0,0.5]);
                hold off
                
                legend('HbO-1','HbO-2');
                xlabel('time [s]'), ylabel('HbO [mM/L]'), title(['HbO-1: ' nirslabels{val1-1} ' HbO-2: ' nirslabels{val2-1}]), grid on;
                
                handles.typeplot=1;
            end
            guidata(hObject, handles);
        end
    case 'HbR + HbR'
        if get(handles.channel1,'Value')~=1
            val1 = get(handles.channel1,'Value');
            nirslabels = handles.datacells{1,1}{2,1}{4,1};
            ttnirs = handles.datacells{1,2}{1,5};
            xxnirsR = handles.datacells{1,2}{1,4};
            
            val2 = get(hObject,'Value');
            set(handles.viewplot,'Visible','on');
            
            if val1 ~= 1 && val2 ~=1
                axes(handles.viewplot);
                plot(handles.viewplot,ttnirs,xxnirsR(:,val1 - 1),'b');
                hold on
                plot(ttnirs,xxnirsR(:,val2 - 1),'Color',[0,0.5,0.5]);
                hold off
                
                legend('HbR-1','HbR-2');
                xlabel('time [s]'), ylabel('HbR [mM/L]'), title(['HbR-1: ' nirslabels{val1-1} ' HbR-2: ' nirslabels{val2-1}]), grid on;
                
                handles.typeplot=1;
            end
            guidata(hObject, handles);
        end
        guidata(hObject, handles);
        
end

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


% --- Executes on button press in save_preproc.
function save_preproc_Callback(hObject, eventdata, handles)
% hObject    handle to save_preproc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.preprocessdata
    data = handles.datacells;
    %check if Pre Processed Data folder exists
    if ~exist([handles.dir '\Pre Processed Data'],'dir')
        mkdir(handles.dir,'Pre Processed Data');
    end
    
    [file, path] = uiputfile('*.mat','Save pre processed data',[handles.dir '\Pre Processed Data\'...
        '_eegfilt-' num2str(handles.eegfiltLval) '-' num2str(handles.eegfiltHval) '_nirsfilt-'...
        num2str(handles.nirsfiltLval) '-' num2str(handles.nirsfiltHval)]);
    
    handles.preprocessedfile = file;
    
    if file ~=0
        save([path file],'data');
        
        % update file path text
        set(handles.filepaths,'String', ['Pre processed: ' file ]);
        
        set(handles.save_check,'Value',1);
        set(handles.save_check,'Visible','on');
        
        msgbox('SUCCESSFULLY saved preprocessed data!');
    else
        
        % update file path text
        set(handles.filepaths,'String', ['Pre processed: ' 'No name' ]);
        
        set(handles.save_check,'Value',0);
        set(handles.save_check,'Visible','on');
        msgbox('UNSUCCESSFULLY saved preprocessed data!');
    end
else
    msgbox('WARNING! Preprocess data before saving them...');
end

guidata(hObject,handles);


% --- Executes on button press in save_check.
function save_check_Callback(hObject, eventdata, handles)
% hObject    handle to save_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_check

% --- Executes on button press in class_button.
function class_button_Callback(hObject, eventdata, handles)
% hObject    handle to class_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.preprocessdata ~= 0
    CLASSIFICATION (handles.figure1,'Classification');
end


% --- Executes on button press in corr_button.
function corr_button_Callback(hObject, eventdata, handles)
% hObject    handle to corr_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.preprocessdata ~= 0
    CORRELATION;
end

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in nirsID.
function nirsID_Callback(hObject, eventdata, handles)
% hObject    handle to nirsID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nirsID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nirsID

val = get(hObject,'Value');
str = get(hObject,'String');
handles.configID = str2double(str{val});

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function nirsID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nirsID (see GCBO)
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


% --- Executes on button press in view_config_button.
function view_config_button_Callback(hObject, eventdata, handles)
% hObject    handle to view_config_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

open([handles.configdir '\Configuration_' num2str(handles.configID) '.txt']);


% --- Executes on button press in add_newconfig_button.
function add_newconfig_button_Callback(hObject, eventdata, handles)
% hObject    handle to add_newconfig_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = get(handles.nirsID,'String');
last = str2num(str{end});

fid = fopen([handles.configdir '\Configuration_' num2str(last + 1) '.txt'],'w');

fprintf(fid, '%s', ['Configuration ' num2str(last + 1)]);
fprintf(fid,'\n');
fprintf(fid, '%s\t%s\t%s', 'S','D','Label');
fprintf(fid,'\n');

fclose(fid);

open([handles.configdir '\Configuration_' num2str(last + 1) '.txt']);

str{last+1} = num2str(last+1);

set(handles.nirsID,'String',str);

guidata(hObject,handles);




% --- Executes on button press in expand_butt.
function expand_butt_Callback(hObject, eventdata, handles)
% hObject    handle to expand_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of expand_butt
if get(hObject,'Value')
    set(handles.loadpanel,'Visible','off');
    set(handles.preprocesspanel,'Visible','off');
    set(handles.filenamespanel,'Visible','off');
    set(handles.analyzepanel,'Visible','off');
    set(handles.viewpanel,'Position',[0 0 1 1]);
else
    set(handles.viewpanel,'Position',handles.viewpanelpos);
    set(handles.loadpanel,'Visible','on');
    set(handles.preprocesspanel,'Visible','on');
    set(handles.filenamespanel,'Visible','on');
    set(handles.analyzepanel,'Visible','on');
end

guidata(hObject,handles);


% --- Executes when viewpanel is resized.
function viewpanel_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to viewpanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% CHANGE this callback to automatically deal with multiple events

% --- Executes on button press in eventsbutton.
function eventsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to eventsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eventsbutton

if get(hObject,'Value')
    
    if get(handles.signal_selection,'Value') == 1 || get(handles.signal_selection,'Value') == 4
        
        if handles.typeplot == 1
            
            evt = handles.datacells{1,2}{1,6};
            
            ev0 = evt(evt(:,2) == 10,1);
            ev1 = evt(evt(:,2) == 1,1);
            ev3 = evt(evt(:,2) == 3,1);
            ev4 = evt(evt(:,2) == 4,1);
            ev7 = evt(evt(:,2) == 7,1);
            ev8 = evt(evt(:,2) == 8,1);
            
            axes(handles.viewplot);
            
            hold on
            y1=get(gca,'ylim');
            plot([ev0 ev0], y1, 'LineStyle','-', 'Color',[0 0 0], 'linewidth',4) 
            if length(ev1) ~= 0
                plot([ev1 ev1], y1, 'LineStyle','-', 'Color',[0 0 0], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev3 ev3], y1, 'LineStyle','-', 'Color',[1 0 0], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev4 ev4], y1, 'LineStyle','-', 'Color',[0 0 1], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev7 ev7], y1, 'LineStyle','-', 'Color',[1 0 0.8], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev8 ev8], y1, 'LineStyle','-', 'Color',[0 0.8 1], 'linewidth',4)
            end
            hold off
            
        elseif    handles.typeplot == 2
            
            evt = handles.datacells{1,2}{1,6};
            
            ev0 = evt(evt(:,2) == 10,1);
            ev1 = evt(evt(:,2) == 1,1);
            ev3 = evt(evt(:,2) == 3,1);
            ev4 = evt(evt(:,2) == 4,1);
            ev7 = evt(evt(:,2) == 7,1);
            ev8 = evt(evt(:,2) == 8,1);
            
            axes(handles.hideplot(1));
            
            hold on
            y1=get(gca,'ylim');
           if length(ev1) ~= 0
                plot([ev1 ev1], y1, 'LineStyle','-', 'Color',[0 0 0], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev3 ev3], y1, 'LineStyle','-', 'Color',[1 0 0], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev4 ev4], y1, 'LineStyle','-', 'Color',[0 0 1], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev7 ev7], y1, 'LineStyle','-', 'Color',[1 0 0.8], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev8 ev8], y1, 'LineStyle','-', 'Color',[0 0.8 1], 'linewidth',4)
            end
            hold off
            hold off
            
            
        elseif    handles.typeplot == 3
            
            evt = handles.datacells{1,2}{1,6};
            
            ev0 = evt(evt(:,2) == 10,1);
            ev1 = evt(evt(:,2) == 1,1);
            ev3 = evt(evt(:,2) == 3,1);
            ev4 = evt(evt(:,2) == 4,1);
            ev7 = evt(evt(:,2) == 7,1);
            ev8 = evt(evt(:,2) == 8,1);
            
            axes(handles.sub1);
            
            hold on
            y1=get(gca,'ylim');
            if length(ev1) ~= 0
                plot([ev1 ev1], y1, 'LineStyle','-', 'Color',[0 0 0], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev3 ev3], y1, 'LineStyle','-', 'Color',[1 0 0], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev4 ev4], y1, 'LineStyle','-', 'Color',[0 0 1], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev7 ev7], y1, 'LineStyle','-', 'Color',[1 0 0.8], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev8 ev8], y1, 'LineStyle','-', 'Color',[0 0.8 1], 'linewidth',4)
            end
            hold off
            hold off
            
            axes(handles.sub2);
            
            hold on
            y1=get(gca,'ylim');
            if length(ev1) ~= 0
                plot([ev1 ev1], y1, 'LineStyle','-', 'Color',[0 0 0], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev3 ev3], y1, 'LineStyle','-', 'Color',[1 0 0], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev4 ev4], y1, 'LineStyle','-', 'Color',[0 0 1], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev7 ev7], y1, 'LineStyle','-', 'Color',[1 0 0.8], 'linewidth',4)
            end
            if length(ev1) ~= 0
                plot([ev8 ev8], y1, 'LineStyle','-', 'Color',[0 0.8 1], 'linewidth',4)
            end
            hold off
            hold off
            
        end
        
    else
        set(hObject,'Value',0);
    end
    
else
    
    channel1_Callback(hObject, eventdata, handles)
    
end



% --- Executes on selection change in signal_selection.
function signal_selection_Callback(hObject, eventdata, handles)
% hObject    handle to signal_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns signal_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from signal_selection

val = get(hObject,'Value');

set(handles.eegfilter,'Visible','off');
set(handles.filter_button,'Value',0);
set(handles.band_power_button,'Value',0);
handles.BP = 0;

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
        set(handles.select_view,'String',handles.raw_select_view);
        set(handles.select_view,'Value',1);
        set(handles.channel1,'Value',1);
    case 2
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
        set(handles.select_view,'String',handles.nirs_power_select_view);
        set(handles.select_view,'Value',1);
        set(handles.channel1,'Value',1);
    case 3
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
        set(handles.select_view,'Visible','off');
        set(handles.channel2,'Visible','off');
        set(handles.select_view,'Value',1);
        set(handles.channel1,'Value',1);
    case 4
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
        set(handles.channel2,'Visible','on');
        set(handles.select_view,'String',handles.prep_select_view);
        set(handles.select_view,'Value',1);
        set(handles.channel1,'Value',1);
        set(handles.channel2,'Value',1);
end

guidata(hObject,handles);



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


% --- Executes on selection change in eegID.
function eegID_Callback(hObject, eventdata, handles)
% hObject    handle to eegID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns eegID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from eegID


% --- Executes during object creation, after setting all properties.
function eegID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eegID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in view_eegconfig_button.
function view_eegconfig_button_Callback(hObject, eventdata, handles)
% hObject    handle to view_eegconfig_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

open([handles.eegconfigdir '\Configuration_' num2str(get(handles.eegID,'Value')) '.txt']);

% --- Executes on button press in add_eegconfig_button.
function add_eegconfig_button_Callback(hObject, eventdata, handles)
% hObject    handle to add_eegconfig_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = get(handles.eegID,'String');
last = str2num(str{end});

fid = fopen([handles.eegconfigdir '\Configuration_' num2str(last + 1) '.txt'],'w');

fprintf(fid, '%s', ' <Name> ') ;
fprintf(fid,'\n');
fprintf(fid, '%s', 'EEG Channel Label');
fprintf(fid,'\n');

fclose(fid);

open([handles.eegconfigdir '\Configuration_' num2str(last + 1) '.txt']);

str{last+1} = num2str(last+1);

set(handles.eegID,'String',str);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function lowfreqeegval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowfreqeegval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in filter_band.
function filter_band_Callback(hObject, eventdata, handles)
% hObject    handle to filter_band (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filter_band contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filter_band
filter_button_Callback(handles.filter_button,eventdata,handles);


% --- Executes during object creation, after setting all properties.
function filter_band_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_band (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in filter_button.
function filter_button_Callback(hObject, eventdata, handles)
% hObject    handle to filter_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filter_button

if get(hObject,'Value')
    
    band = get(handles.filter_band,'Value');
    Fs = handles.datacells{1,1}{1,1}{2,1};
    Fn = Fs/2;
    
    switch band
        
        case 1
            [b,a] = butter(4, [0.5/Fn 4/Fn],'bandpass');
        case 2
            [b,a] = butter(4, [4/Fn 8/Fn],'bandpass');
        case 3
            [b,a] = butter(4, [8/Fn 15/Fn],'bandpass');
        case 4
            [b,a] = butter(4, [8/Fn 12/Fn],'bandpass');
        case 5
            [b,a] = butter(4, [15/Fn 30/Fn],'bandpass');
        case 6
            [b,a] = butter(4, [0.5/Fn 4/Fn],'bandpass');
            
    end
    
    handles.xxeegf = filtfilt(b,a,handles.xxeegRaw);
    handles.xxeeg = handles.xxeegf;
    
else
    handles.xxeeg = handles.datacells{1,2}{1,1};
    set(handles.band_power_button,'Value',0);
end

band_power_button_Callback(handles.band_power_button,eventdata,handles);

guidata(hObject,handles);




% --- Executes on button press in band_power_button.
function band_power_button_Callback(hObject, eventdata, handles)
% hObject    handle to band_power_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of band_power_button

if get(hObject,'Value') && get(handles.filter_button,'Value')
    
    smooth = 500;
    smoothNum = 1/smooth*ones(1,smooth);
    
    bp = filtfilt(smoothNum,1,abs(handles.xxeegf));
    handles.xxeeg = bp;
    %for increasing linewidth of the plot
    handles.BP = 1;
    
else
    set(hObject,'Value',0);
    handles.xxeeg = handles.xxeegf;
    handles.BP = 0;
end

channel1_Callback(handles.channel1,eventdata, handles);
guidata(hObject,handles);



function hifreqnirsval_Callback(hObject, eventdata, handles)
% hObject    handle to hifreqnirsval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hifreqnirsval as text
%        str2double(get(hObject,'String')) returns contents of hifreqnirsval as a double
set(handles.slider_high_nirs,'Value',str2double(get(hObject,'String')));
guidata(hObject,handles);


function lowfreqnirsval_Callback(hObject, eventdata, handles)
% hObject    handle to lowfreqnirsval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowfreqnirsval as text
%        str2double(get(hObject,'String')) returns contents of lowfreqnirsval as a double
set(handles.slider_low_nirs,'Value',str2double(get(hObject,'String')));
guidata(hObject,handles);


function lowfreqeegval_Callback(hObject, eventdata, handles)
% hObject    handle to lowfreqeegval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowfreqeegval as text
%        str2double(get(hObject,'String')) returns contents of lowfreqeegval as a double
set(handles.slider_low_eeg,'Value',str2double(get(hObject,'String')));
guidata(hObject,handles);


function hifreqeegval_Callback(hObject, eventdata, handles)
% hObject    handle to hifreqeegval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hifreqeegval as text
%        str2double(get(hObject,'String')) returns contents of hifreqeegval as a double
set(handles.slider_high_eeg,'Value',str2double(get(hObject,'String')));
guidata(hObject,handles);
