% ------------------------------------------------------------------------
% Run permutation testing on global field power (GFP)
% copyright (c) Huizhen Tang, e-mail: eagtang2007@gmail.com, Feb-14-2018
% ------------------------------------------------------------------------

%% clean workspace
clear all

%% Import single subject's average ERP of condition 1/stimulus type 1
[wave1 wave2 tmw] = readERPmultiple;

% calculate GFP using the computeGFP function for the actual observed data
gfp1 = [];gfp2 = [];
for sj1 = 1:size(wave1,1)
    w1 = squeeze(wave1(sj1,:,:));
    gfp1(sj1,:) = computeGFP(w1,tmw,0);
end
for sj2 = 1:size(wave2,1)
    w2 = squeeze(wave2(sj2,:,:));
    gfp2(sj2,:) = computeGFP(w2,tmw,0);
end
gfpPool = cat(1,gfp1,gfp2); % pool the two sets of waveforms into a single set

%% Calculate the actual observed test statistic (paired-sample t-test)for
% each time point
t_Observed = []; h_Observed = [];
for t = 1:length(tmw)
    g1 = gfp1(:,t);
    g2 = gfp2(:,t);
    [h,p,ci,stats] = ttest(g1,g2,'Alpha',0.1); % run paired-sample t-test
    h_Observed(t) = h;
    t_Observed(t) = stats.tstat; % obtaining the t stats from struct stats
end
%% calculate the sum of absolute t value for each cluster.
% also locate the latency window for each cluster
tSum_Observed = [];
ind = find(h_Observed);
kk = 0; ii = 0;
tw = [];
for i = 1:(length(ind)-1)
    if (ind(i+1) - ind(i)) == 1
        ii = ii + 1;
    elseif ind(i+1) - ind(i) > 1
        kk = kk + 1;
        tw(1,kk) = ind(i-ii); % starting latency
        tw(2,kk) = ii; % latency window
        tSum_Observed(kk) = sum(t_Observed((i-ii):i));
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
    set1Ind = randperm(size(gfpPool,1),size(gfp1,1)); % randomly drawing k=length(fname1)from n=length(pool) numbers as the index for subset 1
    pInd = randperm(size(gfpPool,1)); % shuffle the whole set again and save the index to pInd
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
    gfp_set1 = []; gfp_set2 = [];
    for n = 1:length(set1Ind)
        gfp_set1(n,:) = gfpPool(set1Ind(n),:);
    end
    for k = 1:length(set2Ind)
        gfp_set2(k,:) = gfpPool(set2Ind(k),:);
    end

    %% Calculate the test statistic (npaired-sample t-test) for the current run of random partition
    % each time point
    tPer = []; h_Per = [];
    for t = 1:length(tmw)
        g1 = gfp_set1(:,t);
        g2 = gfp_set2(:,t);
        [h,p,ci,stats] = ttest(g1,g2,'Alpha',0.1); %setting the significance level as Alpha = 0.05
        h_Per(t) = h;
        t_Per(t) = stats.tstat;
    end

    h_PerAll(np,:) = h_Per; % h_PerAll contains all the test results for a total permutation of np times
    t_PerAll(np,:) = t_Per; % t_PerAll contains all the test statistic for a total permutation of np times

    disp('Running')
    np
end

%% creating cluster based histogram
ind = find(h_PerAll);
kk = 0; ii = 0; tSum = []; 
for ip = 1:size(t_PerAll,1)
    for i = 1:(length(ind)-1)
    if (ind(i+1) - ind(i)) == 1
        ii = ii + 1;
    elseif ind(i+1) - ind(i) > 1
        kk = kk + 1;
        tSum(kk) = sum(t_PerAll(ip,(i-ii):i));
        ii = 0;
    end
    end
end

% for i = 1:size(tt,1)
%     k = 0;
%     for j = 1:size(tt,2)
%         if tt(i,j)>0 || tt(i,j)<0
%             k = k + 1;
%             ts = ts + tt(i,j);
%         elseif tt(i,j) == 0 && ts > 0
%             tSum(m) = ts;
%             ts = 0;
%             m = m + 1;
%         elseif tt(i,j) == 0 && ts ==0
%         end
%     end
% end

%% compared the observed t value with the permuated distribution and
% calculate the cluster-based p value 
p_cluster = [];line_t = []; nn = 0;
if ~isempty(tw)
for ns = 1:length(tSum_Observed)
    nlarger = length(find(tSum>tSum_Observed(ns)));
    p_cluster(ns) = nlarger/length(tSum);
    if p_cluster(ns) < 0.05
        nn = nn + 1;
        line_t(1,nn) = tw(1,ns);
        line_t(2,nn) = tw(1,ns) + tw(2,ns);
    end
end
end

%% plot the grand mean gfp waveforms
time_axis = -100:2:800;
avegfp1 = mean(gfp1,1);
avegfp2 = mean(gfp2,1);
y_max = max(max(avegfp1),max(avegfp2));
gray = [0.9 0.9 0.9];

figure
plot(time_axis,avegfp1,'-','Color','k','LineWidth',4)
hold on
plot(time_axis,avegfp2,':','Color','k','LineWidth',4)
if ~isempty(line_t)
for nn = 1:length(line_t(1,:))
% add lines
y = get(gca,'YLim');
h1 = line([line_t(1,nn) line_t(1,nn)],y);
h2 = line([line_t(2,nn) line_t(2,nn)],y);
% set properties of lines 
set([h1 h2],'Color',gray,'LineWidth',0.001)
% Add a patch
patch([line_t(1,nn) line_t(2,nn) line_t(2,nn) line_t(1,nn)],[y(1) y(1) y(2) y(2)],gray)
% The order of the 'Children' of the plot determines which one appears on
% top
% Need to flip it here
set(gca,'children',flipud(get(gca,'children')))
end
end

    
    
    
    
    
    
    
