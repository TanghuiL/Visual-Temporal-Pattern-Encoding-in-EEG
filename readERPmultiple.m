function [wave1,wave2,tmw] = readERPmultiple
% Reads two sets of ERPs from .dat format to matrix 
% Input:
% erps.dat     - ERP wave in dat format{channels x timepoints}
% stt          - epoch starting time point
% endt         - epoch ending time point
% ft           - sampling rate
%
% Output:
% wavePool         - ERP waveform (trial x channel x timepoints)
% tmw          - epoch window [ms] of ERP waveform (e.g. -100:2:500 ms)
% -------------------------------------------------------------------
% copyright (c) Huizhen Tang, e-mail: eagtang2007@gmail.com, Nov-7-2017

%% Import single subject's average ERP of condition 1/stimulus type 1
[fname1,pname1] = uigetfile(...
    { '*.dat','average ERPs of condition 1/stimulus type 1 (*.dat)';'*.*','All Files' }, ...
    'Select mat File(s)', 'Multiselect','on');
% Abort if the user hit 'Cancel'
if isequal(fname1,0)||isequal(pname1,0),
    disp('Aborted.');
    return;
end

%Import single subject's average ERP of condition 2/stimulus type 2
[fname2,pname2] = uigetfile(...
    { '*.dat','average ERPs of condition 1/stimulus type 1 (*.dat)';'*.*','All Files' }, ...
    'Select mat File(s)', 'Multiselect','on');
% Abort if the user hit 'Cancel'
if isequal(fname2,0)||isequal(pname2,0),
    disp('Aborted.');
    return;
end
%% Specify epoch windows and sampling rate
prompt = {'Epoch starts at (e.g. -100 ms)', 'Epoch ends at (e.g. 500 ms)'...
    'EEG reacording sampling rate (e.g. 500 Hz)'};
dlg_title = 'Epoch definition';
num_lines = 1;
defAns = {'-100','800','500'};
answer = inputdlg(prompt,dlg_title,num_lines,defAns);
% Abort if the user clicks 'Cancel'.
if isempty(answer), disp('Aborted.');
    return;
end
[stt status] = str2num(answer{1});
if ~status  % Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end
[endt status] = str2num(answer{2});
if ~status  % Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end
[ft status] = str2num(answer{3});
if ~status  % Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end

%% Pool the averaged ERPs of two conditions/two stimulus types in a single set
% Read the text files
wave1 = []; wave2 = [];
for fn1 = 1:length(fname1)
    ffile = fullfile(pname1,fname1{fn1});
    fprintf(1,'Importing %s\n',ffile);
    [fid,msg] = fopen(ffile, 'rt'); % get file indentifier. Here 'r' sets the permission, which is the type of access of the file, read, write, append, or update. attching a 't' specifies whether to open files in binary or text mode.
    if fid == -1,
        fprintf(1,'Error opening dat file "%s":', ffile);
        disp(msg);
        return;
    end
    r_dt = textscan(fid, '%s'); % put identifier "fid" into a cell array.
    fclose(fid); % unload fid identifier
    data = r_dt{1};  val = [];  j = 0;
    for i = 1:length(data),
        v = str2num(data{i}); % convert string to number. Non number string would therefore return as miss/empty value
        if ~isempty(v)
            j = j + 1;
            val(j) = v; % only save the non-empty v value to val vector, skip it when it is empty
        end
    end
    tmw = stt:(1000/ft):(endt+1);
    tm = length(tmw); % tm is the total number of time points
    elect = length(val)/tm; % elect is the total number of electrodes
    wave1(fn1,:,:) = reshape(val,[elect,tm]);
    fprintf(1,'Finished importing %s\n',ffile);
end

for fn2 = 1:length(fname2)
    ffile = fullfile(pname2,fname2{fn2});
    fprintf(1,'Importing %s\n',ffile);
    [fid,msg] = fopen(ffile, 'rt'); % get file indentifier. Here 'r' sets the permission, which is the type of access of the file, read, write, append, or update. attching a 't' specifies whether to open files in binary or text mode.
    if fid == -1,
        fprintf(1,'Error opening dat file "%s":', ffile);
        disp(msg);
        return;
    end
    r_dt = textscan(fid, '%s'); % put identifier "fid" into a cell array.
    fclose(fid); % unload fid identifier
    data = r_dt{1};  val = [];  j = 0;
    for i = 1:length(data),
        v = str2num(data{i}); % convert string to number. Non number string would therefore return as miss/empty value
        if ~isempty(v)
            j = j + 1;
            val(j) = v; % only save the non-empty v value to val vector, skip it when it is empty
        end
    end
    tmw = stt:(1000/ft):(endt+1);
    tm = length(tmw); % tm is the total number of time points
    elect = length(val)/tm; % elect is the total number of electrodes
    wave2(fn2,:,:) = reshape(val,[elect,tm]);
    fprintf(1,'Finished importing %s\n',ffile);
end
end
