clear;
format compact;
fprintf('\n');
%% compute and store
dtbegeeg = 4; %sec
dtendeeg = 4;
filterlofreqeeg = .5;%Hz
filterhifreqeeg = 80;
filternotcheeg = 60;
numsampleswindoweegpower = 200; % num EEG samples (0.5 s if sampled at 400Hz)
%numsampleswindoweegpower = 125; % num EEG samples (0.5 tw)
filterlofreqnirs = 0.01;
filterhifreqnirs = 0.5;
filternotchnirs = 0;

%time window (sec)
tw = 1;
%overlap (%)
overlap = 0;

infostring = [ '_eegfilt' num2str(filterlofreqeeg) '-' num2str(filterhifreqeeg) '_nirsfilt' num2str(filterlofreqnirs)...
    '-' num2str(filterhifreqnirs) '_tw' num2str(tw) '_overlap' num2str(overlap)];
flag_viewnirs = 0;
flag_vieweeg = 0;



%ALE: useless to assign experimentid..comment only for maintaining the list
%visible-->assign before for cycle
%experimentid = [703,704,803,804,1003,1004,1102,1103,1112,1113,1122,1123,1132,1133,1134,1135];

for experimentid = [1135],
    
    datacells = geteegnirsdataALE( experimentid,dtbegeeg,dtendeeg,filterlofreqeeg,filterhifreqeeg,filternotcheeg,...
        numsampleswindoweegpower,filterlofreqnirs,filterhifreqnirs , filternotchnirs, flag_viewnirs,flag_vieweeg, ...
        tw , overlap...
        );
    
    % save data
    matlabdatafolder = 'C:\Users\alessio\Documents\Ale\UH\Thesis\NIRS_EEG\matlab\matlabdata';
    savestring = ['save '  matlabdatafolder '\nirseeg_preprocessed_' num2str(experimentid) infostring '.mat ' ' datacells'];
    fprintf('%s\n', savestring);
    
    eval( savestring );
end

return;
%%
%{
    DATA STRUCTURE
Each file is one experiment.
The file is a mat file. When the mat file is loaded, it reads a cell
array datacells which is 1 by 22. THe first dimension is always 1, higher
dimensions not used. The second dimension stores data about the experiment.
The data is as follows:

%%   data
    datacells{count,1} = infocells;

    % synced data
    datacells{count,2} = tteegclip;
    datacells{count,3} = xxeegfclip;
    datacells{count,4} = ttnirsclip;
    datacells{count,5} = xxnirsfclip;
    datacells{count,6} = xxnirsfRclip;
    % eeg resampled to nirs time points
    datacells{count,7} = xxeegfclipds;
    % eeg band power
    datacells{count,8} = xxeegfclipbandpower;
    % eeg band power resampled to nirs time points
    datacells{count,9} = xxeegfclipbandpowerds;
    % event times - data clipped at first event, thus first event time is 0
    datacells{count,10} = eventtimes;

    %%  info data
    infocells = cell(3,1);
    infocells{1} = infocellseeg;
    infocells{2} = infocellsnirs;
    infocellseeg = cell(8,1);
    infocellsnirs = cell(8,1);

    infocellsnirs{1} = fpath;
    infocellsnirs{2} = fpath2;
    infocellsnirs{3} = fname;
    infocellseeg{4} = fpatheeg;
    infocellseeg{5} = fnameeeg;
    infocellsnirs{6} = badchans;
    infocellseeg{7} = Fs;
    infocellsnirs{7} = Fsnirs;
    
%}
