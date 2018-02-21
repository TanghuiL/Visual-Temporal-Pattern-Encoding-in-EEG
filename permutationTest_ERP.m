% ------------------------------------------------------------------------
% Run permutation testing on Event-related potentials (ERPs)
% copyright (c) Huizhen Tang, e-mail: eagtang2007@gmail.com, Feb-20-2018
% ------------------------------------------------------------------------

%% clean workspace
clear all

%% Import single subject's average ERP of condition 1/stimulus type 1
[wave1 wave2 tmw] = readERPmultiple;

ch = 10; % select electrode PO9
w1 = squeeze(wave1(:,ch,:));
w2 = squeeze(wave2(:,ch,:));
erpPool = cat(1,w1,w2); % pool the two sets of waveforms into a single set

%% Calculate the actual observed test statistic (paired-sample t-test)for
% each time point
t_Observed = []; h_Observed = [];
for t = 1:length(tmw)
    v1 = w1(:,t);
    v2 = w2(:,t);
    [h,p,ci,stats] = ttest(v1,v2,'Alpha',0.1); % run paired-sample t-test
    h_Observed(t) = h;
    t_Observed(t) = abs(stats.tstat); % obtaining the t stats from struct stats
end
%% calculate the sum of absolute t value for each cluster.
% also locate the latency window for each cluster
tSum_Observed = [];
ind = find(h_Observed);
kk = 0; ii = 0;
tw = [];
for i = 1:(length(ind)-1)
    if (ind(i+1) - ind(i)) == 1 && i+1 < length(ind)
        ii = ii + 1;
    elseif ind(i+1) - ind(i) > 1 && ii > 0 && i < length(ind) % if ii = 0, there is no cluster
        kk = kk + 1;
        tw(1,kk) = ind(i-ii); % starting latency
        tw(2,kk) = ii; % latency window
        tSum_Observed(kk) = sum(t_Observed(ind(i-ii):ind(i)));
        ii = 0;
    elseif ind(i+1) - ind(i) == 1 && i+1 == length(ind)
        kk = kk + 1;
        tw(1,kk) = ind(i-ii); % starting latency
        tw(2,kk) = ii; % latency window
        tSum_Observed(kk) = sum(t_Observed(ind(i-ii):ind(i+1)));
        ii = 0;
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
    for t = 1:length(tmw)
        v1 = erp_set1(:,t);
        v2 = erp_set2(:,t);
        [h,p,ci,stats] = ttest(v1,v2,'Alpha',0.1); %setting the significance level as Alpha = 0.05
        h_Per(t) = h;
        t_Per(t) = abs(stats.tstat);
    end

    h_PerAll(np,:) = h_Per; % h_PerAll contains all the test results for a total permutation of np times
    t_PerAll(np,:) = t_Per; % t_PerAll contains all the test statistic for a total permutation of np times

    disp('Running')
    np
end

%% creating cluster based histogram
kk = 0; ii = 0; tSum = []; 
for iper = 1:n_permute
    tem = h_PerAll(iper,:);
    ind = find(tem>0);
    if length(ind)>3
    for i = 1:(length(ind)-1)
    if ind(i+1) - ind(i) == 1 && i+1 < length(ind)
        ii = ii + 1;
    elseif ind(i+1) - ind(i) > 1 && ii > 0 && i < length(ind)  % if ii == 0, there is no cluster
        kk = kk + 1;
        tSum(kk) = sum(t_PerAll(iper,ind(i-ii):ind(i)));
        ii = 0; 
    elseif ind(i+1) - ind(i) == 1 && i+1 == length(ind)
        kk = kk + 1;
        tSum(kk) = sum(t_PerAll(iper,ind(i-ii):ind(i+1)));
        ii = 0;
    end
    end
    end
end

%% compared the observed t value with the permuated distribution and
% calculate the cluster-based p value 
p_cluster = [];line_t = []; nn = 0;
if ~isempty(tw)
for ns = 1:length(tSum_Observed)
    nlarger = length(find(tSum>tSum_Observed(ns)));
    p_cluster(ns) = nlarger/length(tSum);
    if p_cluster(ns) < 0.05
        nn = nn + 1;
        line_t(1,nn) = 2*tw(1,ns)-100;
        line_t(2,nn) = 2*(tw(1,ns) + tw(2,ns))-100;
    end
end
end

%% plot the grand mean gfp waveforms
time_axis = -100:2:800;
aveErp1 = mean(w1,1);
aveErp2 = mean(w2,1);
y_max = max(max(aveErp1),max(aveErp2));
gray = [0.9 0.9 0.9];

figure
plot(time_axis,aveErp1,'-','Color','k','LineWidth',4)
hold on
plot(time_axis,aveErp2,':','Color','k','LineWidth',4)
if ~isempty(line_t)
for nn = 1:length(line_t(1,:))
% add lines
y = get(gca,'YLim');
h1 = line([line_t(1,nn) line_t(1,nn)],y);
h2 = line([line_t(2,nn) line_t(2,nn)],y);
% set properties of lines 
set([h1 h2],'Color',gray,'LineWidth',0.000001)
% Add a patch
patch([line_t(1,nn) line_t(2,nn) line_t(2,nn) line_t(1,nn)],[y(1) y(1) y(2) y(2)],gray)
% The order of the 'Children' of the plot determines which one appears on
% top
% Need to flip it here

end
set(gca,'children',flipud(get(gca,'children')))
end

    
    
    
    
    
    
    
