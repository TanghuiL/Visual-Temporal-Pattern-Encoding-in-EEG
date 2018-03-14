% ------------------------------------------------------------------------
% Run permutation testing on Event-related potentials (ERPs)
% copyright (c) Huizhen Tang, e-mail: eagtang2007@gmail.com, Feb-20-2018
% last updated Mar-09-2018 
% ------------------------------------------------------------------------

%% clean workspace
clear all

%% Specify the paired conditions for comparison
prompt = {'Comparisons (1 = stdXX; 2=stdOX; 3=edOstdX; 4 =ldXstdX'};
dlg_title = 'Comparisons';
num_lines = 1;
defAns = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,defAns);%%% If the user clicks the Cancel button to close an input dialog box,
% % % Abort if the user clicks 'Cancel'.
if isempty(answer), disp('Aborted.');
    return;
end
[cond status] = str2num(answer{1});
if ~status  %%%Handle empty value returned for unsuccessful conversion
    msgbox('Invalid Number','Error in Parameter settings','error');
end

if cond == 1, stim = 'stdXX';
    test_tw_st = 50; % starting time window is 0
test_tw_end = 200; % end time window is 400
elseif cond == 2,  stim = 'stdOX';
    test_tw_st = 50; % starting time window is 0
test_tw_end = 200; % end time window is 400
elseif cond == 3,  stim = 'edOstdX';
    test_tw_st = 50; % starting time window is 0
test_tw_end = 250; % end time window is 400
elseif cond == 4,  stim = 'ldXstdX';
    test_tw_st = 50; % starting time window is 0
test_tw_end = 250; % end time window is 400
else
end 

% the time window for testing is different for different comparison. 
test_tw = length(test_tw_st+1:test_tw_end);

%% Import single subject's average ERP of condition 1/stimulus type 1
[wave1 wave2 tm] = readERPmultiple;

% select electrode PO9 for comparison
ch = 10; 
w1 = squeeze(wave1(:,ch,test_tw_st+1:test_tw_end));
w2 = squeeze(wave2(:,ch,test_tw_st+1:test_tw_end));
erpPool = cat(1,w1,w2); % pool the two sets of waveforms into a single set

%% Calculate the actual observed test statistic (paired-sample t-test)for
% each time point
threshold = 0.025;  % for two-sided test, for an alpha = 0.05 test, the
% threshold needs to be 0.025. For one-sided test,
% the threshold should be 0.05
t_observed = []; h_observed = [];
for t = 1:test_tw
    v1 = w1(:,t);
    v2 = w2(:,t);
    [h,p,ci,stats] = ttest(v1,v2,'Alpha',threshold); % run paired-sample t-test
    h_observed(t) = h;
    t_observed(t) = stats.tstat;   % obtaining the t stats from struct stats;
end
%% calculate the sum of absolute t value for each cluster.
% also locate the latency window for each cluster
tSum_observed = []; tw = [];
% t_obs_sig = t_observed .* h_observed;   % Save the t statistics with a
% % probability larger than the threshold to t_obs_sig
ind = find(h_observed); % Save the index of significant t statistics to ind
kk = 0; ii = 0;
for i = 1:(length(ind)-1)
    if (ind(i+1) - ind(i)) == 1 & i < (length(ind)-1)
        ii = ii + 1;
    elseif ind(i+1) - ind(i) > 1 & ii > 0 & i < (length(ind)-1)% if ii = 0, there is no cluster
        kk = kk + 1;
        tw(1,kk) = ind(i-ii); % starting latency
        tw(2,kk) = ii; % latency window
        tSum_observed(kk) = sum(t_observed(ind(i-ii):ind(i)));
        ii = 0;
    elseif ind(i+1) - ind(i) == 1 & i == (length(ind)-1)
        kk = kk + 1;
        tSum_observed(kk) = sum(t_observed(ind(i-ii):ind(i+1)));
        tw(1,kk) = ind(i-ii); % starting latency
        tw(2,kk) = ii+1; % latency window
        ii = 0;
    else
    end
end

%% permutation testing
n_permute  = 1000;% number of iteration
t_PerAll = []; h_PerAll = [];

for np = 1:n_permute
    %% Random partition
    % Randomly draw as many as the length of fname1 into subset 1 and put
    % the rest of into subset 2
    set1Ind = randperm(size(erpPool,1),size(w1,1)); % randomly drawing k=length(fname1)from n=length(pool) numbers as the index for subset 1
    pInd = randperm(size(erpPool,1)); % shuffle the whole set again and save the index to pInd
    % getting the rest of the numbers as the index for subset 2
    for i = 1:length(pInd)
        for j = 1:length(set1Ind)
            if pInd(i) == set1Ind(j)
                pInd(i) = 0;
            end
        end
    end
    set2Ind = pInd(pInd > 0);
    % Create file name list for each of the subset
    erp_set1 = []; erp_set2 = [];
    for n = 1:length(set1Ind)
        erp_set1(n,:) = erpPool(set1Ind(n),:);
    end
    for k = 1:length(set2Ind)
        erp_set2(k,:) = erpPool(set2Ind(k),:);
    end
    
    %% Calculate the test statistic (npaired-sample t-test) for the current run of random partition
    % each time point
    tPer = []; h_Per = [];
    for t = 1:test_tw
        v1 = erp_set1(:,t);
        v2 = erp_set2(:,t);
        [h,p,ci,stats] = ttest(v1,v2,'Alpha',threshold); %setting the significance level as Alpha = 0.05
        h_Per(t) = h;
        t_Per(t) = stats.tstat;
    end
    
    h_PerAll(np,:) = h_Per; % h_PerAll contains all the test results for a total permutation of np times
    t_PerAll(np,:) = t_Per; % t_PerAll contains all the test statistic for a total permutation of np times
    
    disp('Running')
    np
end

%% creating cluster based histogram
kk = 0; ii = 0; tSum = [];
t_PerAll_sig = t_PerAll .* h_PerAll; 
for iper = 1:n_permute
    tem = h_PerAll(iper,:);
    ind_p= find(tem>0);
    if length(ind_p)>3
        for i = 1:(length(ind_p)-1)
            if ind_p(i+1) - ind_p(i) == 1 & i < length(ind_p)-1
                ii = ii + 1;
            elseif ind_p(i+1) - ind_p(i) > 1 & ii > 0  % if ii == 0, there is no cluster
                kk = kk + 1;
                tSum(kk) = sum(t_PerAll(iper,ind_p(i-ii):ind_p(i)));
                t_PerAll_sig(iper,ind_p(i-ii):ind_p(i)) = 0;
                ii = 0;
            elseif ind_p(i+1) - ind_p(i) == 1 & i == length(ind_p)-1
                kk = kk + 1;
                tSum(kk) = sum(t_PerAll(iper,ind_p(i-ii):ind_p(i+1)));
                t_PerAll_sig(iper,ind_p(i-ii):ind_p(i+1)) = 0;
                ii = 0;
            end
        end
    end
end

%% compared the observed t value with the permuated distribution and
% calculate the cluster-based p value
p_cluster = [];line_t = []; nn = 0;
if ~isempty(tw)
    for ns = 1:length(tSum_observed)
        nlarger = length(find(abs(tSum)>abs(tSum_observed(ns))));
        %         p_cluster(ns) = nlarger/length(tSum);
         p_cluster(ns) = nlarger/length(find(h_PerAll));
%         p_cluster(ns) = nlarger/(length(tSum) + length(find(t_PerAll_sig>0)));
        if p_cluster(ns) < 0.05
            nn = nn + 1;
            line_t(1,nn) = (tw(1,ns) + test_tw_st)*2-100;
            line_t(2,nn) = (tw(1,ns) + tw(2,ns) + test_tw_st)*2 - 100;
        end
    end
end


%% plot the grand mean waveforms
time_axis = -100:2:800;
aveErp1 = squeeze(mean(wave1(:,ch,:),1));
aveErp2 = squeeze(mean(wave2(:,ch,:),1));

figure
rgb = [0 0 0;255 0 255;0 255 255;0 255 0; 0 0 255]; rgb = rgb/255; 
% magenta [255 0 255]; cyan [0 255 255]; green [0 255 0]; blue [0 0 255]

if strcmp(stim, 'stdOX') || strcmp(stim,'stdXX')
% Solid & dash
plot(time_axis,aveErp1,'-','Color',rgb(1,:),'LineWidth',4)
hold on
plot(time_axis,aveErp2,':','Color',rgb(1,:),'LineWidth',4)

elseif strcmp(stim, 'edOstdX')
% Megenta & Black
plot(time_axis,aveErp1,'-','color',rgb(2,:),'LineWidth',4)
hold on
plot(time_axis,aveErp2,'-','color',rgb(1,:),'LineWidth',4)

elseif strcmp(stim, 'ldXstdX')
% Cyan & Black
plot(time_axis,aveErp1,'-','color',rgb(3,:),'LineWidth',4)
hold on
plot(time_axis,aveErp2,'-','color',rgb(1,:),'LineWidth',4)

else
end

% default lim and tick handling for axises 
xlim([-100 700]); 
xticks([0 100 200 300 400 500 600 700]);
ylim([-5 10]);
yticks([-5 -2.5 0 2.5 5 7.5 10])

%% overlaying the gray bars indicating significant time windows 
gray = [0.9 0.9 0.9]; white = [1 1 1];
if ~isempty(line_t)
    for nn = 1:length(line_t(1,:))
        % add lines
        y = get(gca,'YLim');
%         h1 = line([line_t(1,nn) line_t(1,nn)],y);
%         h2 = line([line_t(2,nn) line_t(2,nn)],y);
%         % set properties of lines
%         set([h1 h2],'Color',gray,'LineWidth',0.000001)
        % Add a patch
        patch([line_t(1,nn) line_t(2,nn) line_t(2,nn) line_t(1,nn)],[y(1)+0.5 y(1)+0.5 y(2)-0.5 y(2)-0.5],gray,'EdgeColor',white)
    end
    % The order of the 'Children' of the plot determines which one appears on
    % top need to flip it here
    set(gca,'children',flipud(get(gca,'children')))   
end
% setting graph parameters
fz = 24;
% grid on
set(gca,'fontsize',fz-5,'fontname','Times New Roman')
ylabel('µV','fontsize',fz,'fontname','Times New Roman');
xlabel('Time (ms)','fontsize',fz,'fontname','Times New Roman');
hold off

% save graphs as figures and png images 
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 4];
saveas(fig,[stim '.png']);
saveas(fig,[stim '.fig']);

% save the latency window of the significant difference of amplitude 
if cond == 1, save stdXX_latencyWindows line_t;
elseif cond == 2,  save stdOX_latencyWindows line_t;
elseif cond == 3,  save edOstdX_latencyWindows line_t;
elseif cond == 4,  save ldXstdX_latencyWindows line_t;
else
end 

% if cond == 1, save Late-stdXX_latencyWindows line_t;
% elseif cond == 2,  save Late-stdOX_latencyWindows line_t;
% elseif cond == 3,  save Late-edOstdX_latencyWindows line_t;
% elseif cond == 4,  save Late-ldXstdX_latencyWindows line_t;
% else
% end 
%% Code ends here






