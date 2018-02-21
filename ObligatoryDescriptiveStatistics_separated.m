% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % For ploting ERP waveforms % % % % % % % % % % %
% % % % % % % % % % % % % % VisA study 1 % % %  % % % % % %  % % % % %
% % % % % % % % % % % % Dr Sussman's lab % % % % % % % % % % % % % % %
% % % % % % % % % Albert Einstein College of Medicine % % % % % % % % %
% % % % Last updated on 9/13/2017 by Huizhen Tang (Joann) % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% clear workspace
clear
%% defining electrodes
chls = {'FPz', 'Fz', 'Cz', 'Pz', 'Oz', 'FP1', 'FP2', 'F7', 'F8', 'PO9', ...
    'PO10', 'FC5', 'FC6', 'FC1', 'FC2', 'T7', 'T8', 'C3', 'C4', 'CP5', ...
    'CP6', 'CP1', 'CP2', 'P7', 'P8', 'P3', 'P4', 'O1', 'O2', 'LM', 'RM', 'EOG' };
%% specify the number of waveforms
prompt = {'Epoch starts at (e.g. -100 ms)','Epoch ends at (e.g. 800 ms)',...
    'Sampling rate (e.g. 500 Hz)'};
dlg_title = 'Waveform properties';
num_lines = 1;
defAns = {'-100','800','500'};
answer = inputdlg(prompt,dlg_title,num_lines,defAns);%%% If the user clicks the Cancel button to close an input dialog box,
% % % Abort if the user clicks 'Cancel'.
if isempty(answer), disp('Aborted.');
    return;
end
[stt status] = str2num(answer{1});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end
[endt status] = str2num(answer{2});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end
[ft status] = str2num(answer{3});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end
%% Get ERP waveforms
wave1 = []; wave2 = [];
% select the 1st group of waveforms
[fname1,pathname1] = uigetfile(...
    { '*std-O.dat*','Select ERP waveforms';'*std-O.dat','All dat Files' }, ...
    'Select the first group of waveforms','Multiselect','on');
% Abort if the user hit 'Cancel'
if isequal(fname1,0)||isequal(pathname1,0),disp('Aborted.'); return; end
% select the 2nd group of waveforms
[fname2,pathname2] = uigetfile(...
    { '*std-X3.dat*','Select ERP waveform';'*std-X3.dat','All dat Files' }, ...
    'Select the second group of waveforms','Multiselect','on');
% Abort if the user hit 'Cancel'
if isequal(fname2,0)||isequal(pathname2,0),disp('Aborted.'); return; end
%% Abort if the total numbers of waveforms selected in the two groups are not equal.
length(fname1)
length(fname2)
if ~isequal(length(fname1),length(fname2)),disp('Aborted.Unequal total number of waveforms');
    return;
end
%% Calculate mean
for sjn = 1:length(fname1)
    %% Load and Read the 1st waveform
    ffile = fullfile(pathname1,fname1{sjn});
    fprintf(1,'Processing %s\n',fname1{sjn});
    [fid,msg] = fopen(ffile, 'rt'); %%% get file indentifier. Here 'r' sets the permission, which is the type of access of the file, read, write, append, or update. attching a 't' specifies whether to open files in binary or text mode.
    if fid == -1,
        fprintf(1,'Error opening dat file "%s":', ffile);
        disp(msg);
        return;
    end
    % identifier "fid" into a cell array.
    data = textscan(fid, '%s'); %%% Before reading a file with textscan, you must open the file with the fopen function. fopen supplies the fid input required by textscan. When you are finished reading from the file, close the file by calling fclose(fid).
    fclose(fid);
    %% reshape vector format data to a matrix
    k = 0; wave_dt = []; data1 = data{1}(100:end);
    for j = 1:length(data1)
        dt = str2num(data1{j});
        if ~isempty(dt), k=k+1; wave_dt(k) = dt; end
    end
    n_row = length(wave_dt)/length(chls); %length(chls) is the total number of channels.
    for n = 1:n_row,
        D = wave_dt(((n-1)*length(chls)+1):n*length(chls));
        wave1(sjn,n,:) = D;
    end
    %% Load and Read the 2nd waveform
    for sjm = 1:length(fname2)
        if fname2{sjm}(12:17) == fname1{sjn}(12:17)
            ffile = fullfile(pathname2,fname2{sjm});
            fprintf(1,'Processing %s\n',fname2{sjm});
            [fid,msg] = fopen(ffile, 'rt'); %%% get file indentifier. Here 'r' sets the permission, which is the type of access of the file, read, write, append, or update. attching a 't' specifies whether to open files in binary or text mode.
            if fid == -1,
                fprintf(1,'Error opening dat file "%s":', ffile);
                disp(msg);
                return;
            end
            % identifier "fid" into a cell array.
            data = textscan(fid, '%s'); %%% Before reading a file with textscan, you must open the file with the fopen function. fopen supplies the fid input required by textscan. When you are finished reading from the file, close the file by calling fclose(fid).
            fclose(fid);
            %% reshape vector format data to a matrix
            k = 0; wave_dt = []; data1 = data{1}(100:end);
            for j = 1:length(data1)
                dt = str2num(data1{j});
                if ~isempty(dt), k=k+1; wave_dt(k) = dt; end
            end
            n_row = length(wave_dt)/length(chls); %length(chls) is the total number of channels.
            for n = 1:n_row,
                D = wave_dt(((n-1)*length(chls)+1):n*length(chls));
                wave2(sjn,n,:) = D;
            end
        end
    end
end
wave_dif = wave1-wave2;
%% descriptive statistics
stopper = 1;
while stopper > 0
% specify the parameters for computing ERP mean value
prompt = {'Component label','Peak latency (e.g. 146 ms) of wave1',...
    'Peak latency (e.g. 146 ms) of wave2','Latency interval (e.g. 40 ms)',...
    'Electrode','Condition','Orientation change'};
dlg_title = 'ERP mean value';
num_lines = 1;
defAns = {'p1','146','146','40','5','PiR','10'};
answer = inputdlg(prompt,dlg_title,num_lines,defAns);%%% If the user clicks the Cancel button to close an input dialog box,
% % % Abort if the user clicks 'Cancel'.
if isempty(answer), disp('Aborted.');
    return;
end
epName = answer{1};
[pkl1 status] = str2num(answer{2});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end
[pkl2 status] = str2num(answer{3});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end
[intvl status] = str2num(answer{4});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end
[chl status] = str2num(answer{5});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end
cond = answer{6};
ori = answer{7};
%% Compute the specified ERP mean value
pkl1 = (pkl1 + 100)/(1000/ft); pkl2 = (pkl2 + 100)/(1000/ft);
erp1 = squeeze(mean(wave1(:,pkl1-intvl/((1000/ft)*2):pkl1+intvl/((1000/ft)*2),chl),2));%calculate mean
erp1sd = squeeze(std(wave1(:,pkl1-intvl/((1000/ft)*2):pkl1+intvl/((1000/ft)*2),chl),1,2));%calculate standard deviation
erp2 = squeeze(mean(wave2(:,pkl2-intvl/((1000/ft)*2):pkl2+intvl/((1000/ft)*2),chl),2));%mean
erp2sd = squeeze(std(wave2(:,pkl2-intvl/((1000/ft)*2):pkl2+intvl/((1000/ft)*2),chl),1,2));%standard deviation
erpDif = erp1 - erp2;
erpDifsd = squeeze(std(wave_dif(:,pkl1-intvl/((1000/ft)*2):pkl1+intvl/((1000/ft)*2),chl),1,2));%standard deviation
erp = cat(2,erp1,erp2,erp1sd,erp2sd,erpDif,erpDifsd);
gmerp = mean(erp,1);
[h,p,ci,stats] = ttest(erpDif);
erp = cat(1,erp,gmerp);
%% Write the values to a spreadsheet 
col_header = {'Task','Orientation change','Electrode','Std-O','Std-X3','Std-O(sd)',...
    'Std-X3(sd)','Difference','Difference(sd)'};
row_header(1:length(fname1),1) = {cond};
row_header(1:length(fname1),2) = {ori};
row_header(length(fname1)+1,3) = {'group average'};
row_header(length(fname1)+2:length(fname1)+5,3) = {'h ( t test significant if 1)','p value','t value','df'};
row_header(length(fname1)+2:length(fname2)+5,4) = {h,p,stats.tstat,stats.df};
if chl == 5
    row_header(1:length(fname1),3) = {'Oz'};
xlswrite('stdO-stdX3-Oz',col_header,epName,'A1');
xlswrite('stdO-stdX3-Oz',row_header,epName,'A2');
xlswrite('stdO-stdX3-Oz',erp,epName,'D2');
elseif chl == 10
    row_header(1:length(fname1),3) = {'PO9'};
xlswrite('stdO-stdX3-PO9',col_header,epName,'A1');
xlswrite('stdO-stdX3-PO9',row_header,epName,'A2');
xlswrite('stdO-stdX3-PO9',erp,epName,'D2'); %Dont change the order of writing. Data should be the last to write. Otherwise it would be overwritten by row_header
elseif chl == 11
    row_header(1:length(fname1),3) = {'PO10'};
xlswrite('stdO-stdX3-PO10',col_header,epName,'A1');
xlswrite('stdO-stdX3-PO10',row_header,epName,'A2');
xlswrite('stdO-stdX3-PO10',erp,epName,'D2'); %Dont change the order of writing. Data should be the last to write. Otherwise it would be overwritten by row_header
end
%% Decide whether continue to compute the mean for a new ERP component
prompt = {'Continue to calcuate a new ERP? (1= YES; 0=NO)'};
dlg_title = 'new ERP';
num_lines = 1;
defAns = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,defAns);
%Abort if the user clicks 'cancel'
if isempty(answer),disp('Aborted');return;end
[stopper status] = str2num(answer{1});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end
end
