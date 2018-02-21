% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % Matching trigger codes extracted from EXPO xml file % % % % % % %
% % % % % % % % % % And then modify event file for EEG data % % % % % % % %
% % % % % % % % % % VisA study (for Active Condition) % % %  % % % % % % %
% % % % % % % % % % % Albert Einstein College of Medicine % % % % % % % % %
% % % % % % Last updated on 07/06/2016 by Huizhen Tang (Joann) % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear

%% Selecting event file exported from NeuroScan % % % % % % % % %
[filename1,pathname1] = uigetfile(...
    { '*.ev2','Event file exported from NeuroScan';'*.*','All Files' }, ...
    'Select .ev2 file(s)', ...
    'Multiselect','off');
% Abort if the user hit 'Cancel'
if isequal(filename1,0)||isequal(pathname1,0),
    disp('Aborted.');
    return;
end
% cd(pathname1)
disp(filename1)

%% Slecting mat file with the triiger codes extracted from EXPO xml file % % % % % % % %
[filename2,pathname2] = uigetfile(...
    { '*.mat','Trigger codes extracted from EXPO (*.mat)';'*.*','All Files' }, ...
    'Select .mat file(s)', ...
    'Multiselect','off');
% Abort if the user hit 'Cancel'
if isequal(filename2,0)||isequal(pathname2,0),
    disp('Aborted.');
    return;
end
load([pathname2 filename2])
disp(filename2)

%% Reads data from an open text file indentified by the file
%Reading event file exported from NeuroScan .ev2 file % % % % % %
ffile = fullfile(pathname1,filename1);
fprintf(1,'Processing %s\n',ffile);
[fid,msg] = fopen(ffile, 'rt'); %%% get file indentifier. Here 'r' sets the permission, which is the type of access of the file, read, write, append, or update. attching a 't' specifies whether to open files in binary or text mode.
if fid == -1,
    fprintf(1,'Error opening ev2 file "%s":', ffile);
    disp(msg);
    return;
end
% identifier "fid" into a cell array.
triggers = textscan(fid, '%s'); %%% Before reading a file with textscan, you must open the file with the fopen function. fopen supplies the fid input required by textscan. When you are finished reading from the file, close the file by calling fclose(fid).
fclose(fid);

%% reshape data dimension
star = 6 + 1; %%% 6 is the number of colums (this line is to take out the headers
data = triggers{1}(star:end); %%% taking the data from begining to the end
for j=1:length(data),
    dd(j) = str2num(data{j});
end
num_row = length(dd)/6;
for i=1:num_row,
    D = dd(((i-1)*6 + 1):i*6);
    etrgg(i,:) = D;
end
%% matching trigger codes from NeuroScan and EXPO % % % % % % % % % % % %
etime = etrgg(:,6); %% get etime from etrgg for the first time
intvl_ini = etime(2) - etime(1);
intvl_end = etime(end) - etime(end-1);
if intvl_ini > 0.08 && intvl_ini < 0.12
    etrgg(1,:) = [];
else
    disp('First trigger is missing')
end

if intvl_end > 0.08 && intvl_end < 0.12
    etrgg((end-1):end,:) = [];
    
elseif intvl_end > 1.70 && intvl_end < 1.80
    etrgg(end,:) = [];
    disp('Final testing trigger is missing')
elseif intvl_end > 2
    disp('Final two testing triggers and at least one stimulus trigger are missing')
else
    disp('Final two testing triggers are missing')
end
etime = etrgg(:,6);  %% get etime from etrgg for the second time after deleting the trigger check codes
etrgg(1,2) = trgg(1);
k = 1;
for ii = 2:length(etime)
    k = k + 1;
    interval = etime(ii) - etime(ii-1);
    if interval > 0.650 && interval < 0.750
        k = k;
        etrgg(ii,2) = trgg(k);
    elseif interval > 1.250 && interval < 1.500
        k = k + 1;
        etrgg(ii,2) = trgg(k);
    elseif interval > 1.95 && interval < 2.25
        k = k + 2;
        etrgg(ii,2) = trgg(k);
    elseif interval > 2.6 && interval < 3.0
        k = k + 3;
        etrgg(ii,2) = trgg(k);
    elseif interval > 3.25 && interval < 3.75
        k = k + 4;
        etrgg(ii,2) = trgg(k);
    else
    end
end

%% Construct output filename
[~,name,~] = fileparts(filename1);
outfname = fullfile(pathname1,[name '.evt']);

if exist(outfname) ~= 0,
    fprintf('File creation failed: file %s already exists\n', outfname);
    msgbox(['File ' outfname ' already exists'],'Cannot create file','error');
    return;
end

fprintf(1,'Creating %s\n', [name '.evt']);

fid=fopen(outfname, 'w');
fprintf(fid, '# \t Type \t Response \t Acc \t RT \t Offset \n');
fclose(fid);

dlmwrite(outfname, etrgg, 'delimiter', '\t', 'precision', '%6.3f', '-append');

%% code ends here

