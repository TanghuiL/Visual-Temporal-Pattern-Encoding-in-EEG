% ------------------------------------------------------------------------
% Compute Global Map DISSimilarity (GMD/DISS) from epoched EEG waves
% For multiple pair-comparions for VisA study.
%
% gmd are computed to describe the difference between two different time
% points within one condition or the same time point between two different
% conditions
% ------------------------------------------------------------------------
% Compute Global Map DISSimilarity (GMD/DISS) from epoched EEG waves
% Input (by calling the computeGFP function):
% type         - comparison between two time points or between conditions
%                 type = 1, between two time points t1 and t2
%                 type = 2, between two conditions c1 and c2
% wave         - waveform data in 2D (channels x timepoints) array
% tmw          - epoch window [ms] of ERP waveform (e.g. -100:2:500 ms)
% gfp          - global field power of each time point
% nch          - selected electrodes for gfp calculation(ie.EOG excluded)

% Output:
% gmd          - global map dissimilarity index
%
% gmd can be computed to describe the difference between two different time
% points within one condition or the same time point between two different
% conditions
%
% copyright (c) Huizhen Tang, e-mail: eagtang2007@gmail.com, Feb-27-2018
% ------------------------------------------------------------------------

clear all

%% compute gfp and gmd for each ERP
chls = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ...
    21 22 23 24 25 26 27 28 29 30 31]; % chls specifies the electrodes included in the calculation

% call readERP function to convert text (.dat) waveform to 2D array.
[wave1 tmw1 fname1] = readERP; %Call the "readERP" function to convert dt file to 2D array into wave1
[wave2 tmw2 fname2] = readERP; %Call again the function to convert dt file to 2D array into wave2

% calculate gfp for first ERP
nch = 0;
gfp1 = computeGFP(wave1,tmw1,nch);
% calculate gmd for first ERP
gmd1 = []; aveu = []; wave = wave1; gfp = gfp1; tmw = tmw1;
for t = 1:tmw
    aveu(t) = mean(wave(:,t));
end
for t1 = 1:tmw
    for t2 = 1:tmw
        difuv = [];
        for n = 1:length(chls)
            difuv(n) = ((wave(n,t1) - aveu(t1))/gfp(t1) - (wave(n,t2) - aveu(t2))/gfp(t2))^2;
        end
        gmd1(t1,t2) = sqrt(sum(difuv)/n);
    end
end
% calculate gfp for second ERP
gfp2 = computeGFP(wave2,tmw2,nch);
% calculate gmd for second ERP
gmd2 = []; aveu = []; wave = wave2; gfp = gfp2;tmw = tmw2;
for t = 1:tmw
    aveu(t) = mean(wave(:,t));
end
for t1 = 1:tmw
    for t2 = 1:tmw
        difuv = [];
        for n = 1:length(chls)
            difuv(n) = ((wave(n,t1) - aveu(t1))/gfp(t1) - (wave(n,t2) - aveu(t2))/gfp(t2))^2;
        end
        gmd2(t1,t2) = sqrt(sum(difuv)/n);
    end
end

%% paired comparison
gmdPaired = []; aveu1 = []; aveu2 = []; u = []; v = [];
for t = 1:tmw1
    aveu1(t) = mean(wave1(:,t));
    aveu2(t) = mean(wave2(:,t));
end
for t = 1:tmw1
    difuv = [];
    for n = 1:length(chls)
        difuv(n) = ((wave1(n,t) - aveu1(t))/gfp1(t) - (wave2(n,t) - aveu2(t))/gfp2(t))^2;
    end
    gmdPaired(t) = sqrt(sum(difuv)/n);
end

%% Calculate the gmd of successive maps for each ERP
% first ERP
gmd1_suc = [];
for t1 = 1:(size(gmd1,1)-10)
    gmd1_suc(t1) = gmd1(t1,t1+10);    
end
% second ERP
gmd2_suc = [];
for t1 = 1:(size(gmd2,1)-10)
    gmd2_suc(t1) = gmd2(t1,t1+10);    
end

% Normalization of gmd_suc (in reference to prestimulus baseline -100 to 0)
gmd1_suc_nrm = []; 
for t1 = 1:length(gmd1_suc)
    gmd1_suc_nrm(t1) = (gmd1_suc(t1)-mean(gmd1_suc(1:50)))/mean(gmd1_suc(1:50));
end
% second ERP
gmd2_suc_nrm = [];
for t1 = 1:length(gmd2_suc)
    gmd2_suc_nrm(t1) = (gmd2_suc(t1)-mean(gmd2_suc(1:50)))/mean(gmd2_suc(1:50)); 
end


%% Calculate the gmd of a specificy latency band as function of time
%% Specify the latency window for change detection 
prompt = {'Comparison (1 = stdXX; 2=stdOX; 3=edOstdX; 4 =ldXstdX',...
    'latency window lower end (ms)', 'latency window upper end (ms)',...
    'band (1=simple adaptation; 2 = change detection'};
dlg_title = 'latency windows';
num_lines = 1;
defAns = {'4','268','308','1'};
answer = inputdlg(prompt,dlg_title,num_lines,defAns);%%% If the user clicks the Cancel button to close an input dialog box,
% % % Abort if the user clicks 'Cancel'.
if isempty(answer), disp('Aborted.');
    return;
end
[cond status] = str2num(answer{1});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end
[lat_st status] = str2num(answer{2});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end
[lat_end status] = str2num(answer{3});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end

band = answer{4};
if strcmp(band,'1')
    band = 'simpleAdaptation';
elseif strcmp(band,'2')
    band = 'changeDetection';
else 
end 

% latency window convert to current x axis 
lat_st = lat_st/2 + 50;
lat_end = lat_end/2 + 50;
% first ERP
gmd1_band = []; 
for t1 = 1:size(gmd1,1)
    gmd1_band(t1) = mean(gmd1(t1,lat_st:lat_end),2);    
end
% second ERP
gmd2_band = [];
for t1 = 1:size(gmd2,1)
    gmd2_band(t1) = mean(gmd2(t1,lat_st:lat_end),2);    
end

%% ploting
if strcmp('i',fname1(2))
    st = 5;
else st = 4
end
fz = 24;
rgb = [0 0 0;255 0 255;0 255 255;0 255 0; 0 0 255; 200 200 200]; rgb = rgb/255; 
% magenta [255 0 255]; cyan [0 255 255]; green [0 255 0]; blue [0 0 255]
if cond == 1
    lineColor = rgb([1,end],:);
elseif cond == 2
    lineColor = rgb([1,end],:);
elseif cond == 3
    lineColor = rgb([2,1],:);
elseif cond == 4 
    lineColor = rgb([3,1],:);
else
end 


% plot gmd1
fz = 36;
figure(1)
imagesc(gmd1,[0 2])
axis xy
xlim([0 400]); xticks([0 50 100 150 200 250 300 350 400])
xticklabels({'-100','0','100','200','300','400','500','600','700'})
ylim([0 400]); yticks([50 100 150 200 250 300 350 400])
yticklabels({'0','100','200','300','400','500','600','700'})
set(gca,'fontsize',fz-5,'fontname','Times New Roman')
ylabel('Time (ms)','fontsize',fz,'fontname','Times New Roman');
xlabel('Time (ms)','fontsize',fz,'fontname','Times New Roman');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 8];
saveas(fig,[fname1(st:end-4) '-gmd.png']);
saveas(fig,[fname1(st:end-4) '-gmd.fig']);

% plot gmd2
fz = 36;
figure(2)
imagesc(gmd2,[0 2])
axis xy
xlim([0 400]); xticks([0 50 100 150 200 250 300 350 400])
xticklabels({'-100','0','100','200','300','400','500','600','700'})
ylim([0 400]); yticks([50 100 150 200 250 300 350 400])
yticklabels({'0','100','200','300','400','500','600','700'})
set(gca,'fontsize',fz-5,'fontname','Times New Roman')
ylabel('Time (ms)','fontsize',fz,'fontname','Times New Roman');
xlabel('Time (ms)','fontsize',fz,'fontname','Times New Roman');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 8];
saveas(fig,[fname2(st:end-4) '-gmd.png']);
saveas(fig,[fname2(st:end-4) '-gmd.fig']);

fz = 30;
% Paired gmd1 vs. gmd2
figure(3)
plot(1:450,gmdPaired(1:450),'-','color',rgb(1,:),'LineWidth',4)
xlim([0 400]); xticks([0 50 100 150 200 250 300 350 400])
xticklabels({'','0','100','200','300','400','500','600','700'})
% lgd = legend([fname1(st:end-4) '-vs-' fname2(st:end-4)],'Location','NorthOutside');
% lgd.Orientation = 'vertical'; lgd.Box = 'off'; lgd.FontSize = 16;
grid on
set(gca,'fontsize',fz-5,'fontname','Times New Roman')
ylabel('DISS','fontsize',fz,'fontname','Times New Roman');
xlabel('Time (ms)','fontsize',fz,'fontname','Times New Roman');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 4];
saveas(fig,[fname1(st:end-4) '-vs-' fname2(st:end-4) '-gmd.png']);
saveas(fig,[fname1(st:end-4) '-vs-' fname2(st:end-4) '-gmd.fig']);

%% plotting gmd for different components (latency bands)
figure
plot(1:450,gmd1_band(1:450),'-','color',lineColor(1,:),'LineWidth',4)
hold on 
plot(1:450,gmd2_band(1:450),'-','color',lineColor(2,:),'LineWidth',4)
xlim([0 400]); xticks([0 50 100 150 200 250 300 350 400])
xticklabels({'','0','100','200','300','400','500','600','700'})
% lgd = legend([fname1(st:end-4) fname2(st:end-4)],'Location','NorthOutside');
% lgd.Orientation = 'vertical'; lgd.Box = 'off'; lgd.FontSize = 16;
grid on
set(gca,'fontsize',fz-5,'fontname','Times New Roman')
ylabel('DISS','fontsize',fz,'fontname','Times New Roman');
xlabel('Time (ms)','fontsize',fz,'fontname','Times New Roman');
hold off 
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 4];
saveas(fig,[fname1(st:end-4) 'vs' fname2(st:end-4) '-gmd-band.png']);
saveas(fig,[fname1(st:end-4) 'vs' fname2(st:end-4) '-gmd-band.fig']);

% plot gmd1_band 
figure
plot(1:450,gmd1_band(1:450),'-','color',rgb(1,:),'LineWidth',4)
xlim([0 400]); xticks([0 50 100 150 200 250 300 350 400])
xticklabels({'','0','100','200','300','400','500','600','700'})
% lgd = legend([fname1(st:end-4) fname2(st:end-4)],'Location','NorthOutside');
% lgd.Orientation = 'vertical'; lgd.Box = 'off'; lgd.FontSize = 16;
grid on
set(gca,'fontsize',fz-5,'fontname','Times New Roman')
ylabel('DISS','fontsize',fz,'fontname','Times New Roman');
xlabel('Time (ms)','fontsize',fz,'fontname','Times New Roman');
hold off 
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 4];
saveas(fig,[fname1(st:end-4) '-gmd-band-' band '.png']);
saveas(fig,[fname1(st:end-4) '-gmd-band-' band '.fig']);

% plot gmd2_band 
figure
plot(1:450,gmd2_band(1:450),'-','color',rgb(1,:),'LineWidth',4)
xlim([0 400]); xticks([0 50 100 150 200 250 300 350 400])
xticklabels({'','0','100','200','300','400','500','600','700'})
% lgd = legend([fname1(st:end-4) fname2(st:end-4)],'Location','NorthOutside');
% lgd.Orientation = 'vertical'; lgd.Box = 'off'; lgd.FontSize = 16;
grid on
set(gca,'fontsize',fz-5,'fontname','Times New Roman')
ylabel('DISS','fontsize',fz,'fontname','Times New Roman');
xlabel('Time (ms)','fontsize',fz,'fontname','Times New Roman');
hold off 
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 4];
saveas(fig,[fname2(st:end-4) '-gmd-band-' band '.png']);
saveas(fig,[fname2(st:end-4) '-gmd-band-' band '.fig']);


% plot gmd of successive maps 
figure
plot(1:401,gmd1_suc_nrm(1:401),'-','color',lineColor(1,:),'LineWidth',4)
hold on 
plot(1:401,gmd2_suc_nrm(1:401),'-','color',lineColor(2,:),'LineWidth',4)
xlim([0 400]); xticks([0 50 100 150 200 250 300 350 400])
xticklabels({'','0','100','200','300','400','500','600','700'})
% lgd = legend([fname1(st:end-4) fname2(st:end-4)],'Location','NorthOutside');
% lgd.Orientation = 'vertical'; lgd.Box = 'off'; lgd.FontSize = 16;
grid on
set(gca,'fontsize',fz-5,'fontname','Times New Roman')
ylabel('DISS','fontsize',fz,'fontname','Times New Roman');
xlabel('Time (ms)','fontsize',fz,'fontname','Times New Roman');
hold off 
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 4];
saveas(fig,[fname1(st:end-4) 'vs' fname2(st:end-4) '-gmd-suc.png']);
saveas(fig,[fname1(st:end-4) 'vs' fname2(st:end-4) '-gmd-suc.fig']);



%% Code ends here.





%% functions used in the main code
function [wave,tm,filename] = readERP
%% defining channels
chls = {'FPz', 'Fz', 'Cz', 'Pz', 'Oz', 'FP1', 'FP2', 'F7', 'F8', 'PO9', ...
    'PO10', 'FC5', 'FC6', 'FC1', 'FC2', 'T7', 'T8', 'C3', 'C4', 'CP5', ...
    'CP6', 'CP1', 'CP2', 'P7', 'P8', 'P3', 'P4', 'O1', 'O2', 'LM', 'RM', 'EOG' };
%% select waveform
[filename,pathname] = uigetfile(...
    { '*.dat*','Select ERP waveform';'*.dat','All dat Files' }, ...
    'Select the first waveform','Multiselect','off');
% Abort if the user hit 'Cancel'
if isequal(filename,0)||isequal(pathname,0),disp('Aborted.'); return; end
%% Load and Read the selected waveform
ffile = fullfile(pathname,filename);
fprintf(1,'Processing %s\n',filename);
[fid,msg] = fopen(ffile, 'rt'); %%% get file indentifier. Here 'r' sets the permission, which is the type of access of the file, read, write, append, or update. attching a 't' specifies whether to open files in binary or text mode.
if fid == -1,
    fprintf(1,'Error opening dat file "%s":', ffile);
    disp(msg);
    return;
end
% identifier "fid" into a cell array.
str_dt = textscan(fid, '%s'); %%% Before reading a file with textscan, you must open the file with the fopen function. fopen supplies the fid input required by textscan. When you are finished reading from the file, close the file by calling fclose(fid).
fclose(fid);
%% reshape vector format data to a matrix
k = 0; dt = []; data = str_dt{1};
for i = 1:length(data)
    n_dt = str2num(data{i});
    if ~isempty(n_dt)
        k=k+1;
        dt(k) = n_dt;
    end
end
tm = length(dt)/length(chls); % this gives the total number of time points
elect = length(chls);
wave(:,:) = reshape(dt, [elect, tm]);
end


function gfp = computeGFP(wave,tmw,nch)
% ------------------------------------------------------------------------
% Compute Global Field Power (gfp) from epoched EEG waves
% Input (by calling the readERP function):
% wave         - waveform data in 2D (channels x timepoints) array
% tmw          - epoch window [ms] of ERP waveform (e.g. -100:2:500 ms)
% Output:
% gfp          - global field power
%
% copyright (c) Huizhen Tang, e-mail: eagtang2007@gmail.com, Nov-7-2017
% ------------------------------------------------------------------------
gfp = [];
% %% Specify electrodes for gfp computation
if nch == 0;
    chls = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ...
        21 22 23 24 25 26 27 28 29 30 31]; % chls specifies the electrodes included in the calculation
else disp('Invalid input of nch!')
end
wave = wave(chls, :); % exclude the unselected electrode from wave
% computing gfp for each time point and iterate through all the time points
for t = 1:tmw
    aveu = mean(wave(:,t)); % aveu is the mean voltage at time point tmw(t)
    difu = [];
    for n = 1:length(chls)
        difu(n) = (wave(n,t)-aveu)^2;
    end
    gfp(t) = sqrt(sum(difu)/length(chls));
end

end

%% end of functions 