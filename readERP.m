function [wave,tmw] = readERP
% Reads ERP waveform from .dat format to matrix with
% Input:
% erps.dat     - ERP wave in dat format{channels x timepoints}
% stt          - epoch starting time point
% endt         - epoch ending time point
% ft           - sampling rate
%
% Output:
% wave         - ERP waveform (subject x channel x timepoints)         
% tmw          - epoch window [ms] of ERP waveform (e.g. -100:2:500 ms)
% -------------------------------------------------------------------
% copyright (c) Huizhen Tang, e-mail: eagtang2007@gmail.com, Nov-7-2017

%% Load and read individual ERP waveform of a specific condition
[filename,pathname] = uigetfile({ '*.dat*', 'ERP waveform'; ...
    '*.*','All Files' }, 'Select .dat file(s)', ...
    'Multiselect','off');
% Abort if the user hit 'Cancel'
if isequal(filename,0)||isequal(pathname,0),
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

%% read the text file of the selected ERP waveform 
wave = [];
ffile = fullfile(pathname,filename);
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
wave(:,:) = reshape(val,[elect,tm]);
fprintf(1,'Finished importing %s\n',ffile);
end