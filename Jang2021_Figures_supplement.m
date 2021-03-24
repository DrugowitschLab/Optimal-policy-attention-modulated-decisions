%% Creat supplemental figures for Jang et al., 2021
% Optimal policy for attention-modulated decisions explains human fixation behavior
% Run this block first, then any other block below to make the figures

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET CODE DIRECTORY HERE TO RUN ALL SCRIPTS
rootdir = '/Users/anthony_jang/Dropbox (MIT)/Research/Ongoing/att_paper/202009-elife_submission/Code';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath(rootdir));
data_dir = fullfile(rootdir,'data');
optpolicy_dir = fullfile(rootdir,'data/optimal_policy_save');
set(0,'defaultAxesFontSize',20);



%% Figure 4 - Figure supplement 1-3
% Figure S1: Effect of item values on attention switch rate and fixation duration across trials for human data, optimal model, and aDDM
% Figure S2: Effect of passed time on switch probability and fixation duration within trials

% Initialize variables
whichTask = 1;
excludeLastFix = 1;   % exclude last fixation duration, since it's cut short by choosing
useBinnedX = 1;   % Create time bins for switch prob data for smoother results
useBonferroni = 1;  % If set to 0, use FDR instead

% Load behavioral data
taskstr = 'aVBDM';
load(fullfile(data_dir,sprintf('datastruct_%s.mat',taskstr)));
% Load params
params = getParams(whichTask);


% Set aDDM parameters
% 1 = Parameters used in original paper (Krajbich et al., 2010)
% 2 = Set parameters to match snr with Bayesian model (compare performance)
aDDM_param_type = 1;
if aDDM_param_type==1
    dt = 0.001;
    d = 0.0002;
    sig2_origModel = 0.02^2;
    aGamma_k = 0.3;
    k = d/dt;
    sig2_k = sig2_origModel/dt;
    decbound = 1;
else
    dt = params.dt;
    k = 1;
    aGamma_k = sqrt(params.aGamma);
    sig2_k = 2*params.sig2;
    % boundary that maximizes reward rate
    decbound = 11;
end

% Get empty datastruct for both models
iter = 3;  % To increase trials, repeat the human trials X times
dstruct_addm = getEmptyDstruct_realData(dstruct_real,iter);
dstruct_opt = getEmptyDstruct_realData(dstruct_real,iter);
numsubs = length(dstruct_real);
% Get empirical distribution of fixation behavior for each difficulty level
fixdist = getEmpiricalFixationDist(dstruct_real);
maxdectime = 8;

% RT bin info
binstep = 0.5;
binedges_rt = 1:binstep:4;
binctr_rt = binedges_rt(1:end-1)+binstep/2;

% Get all unique value difference & value sum
valdiff_unique = unique(getValdiff(dstruct_real))';
valsum_unique = unique(getValsum(dstruct_real))';

% Load optimal policy information
optpolicy = getOptPolicy_aVBDM(params,optpolicy_dir);

% Simulate behavior
every10Prc = floor(prctile(1:numsubs,10:10:100));
for s = 1:numsubs
    if any(s==every10Prc), fprintf('%.0f%% ',find(s==every10Prc)*10); end
    % Generate trials - draw from the prior distribution from the optimal model
    numtrials = length(dstruct_opt(s).trialnum);
    for ti = 1:numtrials
        % Optimal policy
        [dstruct_opt(s).choice(ti),dstruct_opt(s).rt(ti),dstruct_opt(s).tItem(ti,:),~,dstruct_opt(s).fixdur{ti},dstruct_opt(s).fixitem{ti},~]...
            = runTrial_aVBDM(dstruct_opt(s).itemval(ti,:),optpolicy,params);
        
        % aDDM: simulate behavior for this decision boundary
        thisvaldiff = abs(dstruct_addm(s).itemval(ti,1)-dstruct_addm(s).itemval(ti,2));
        if thisvaldiff > 8, thisvaldiff = 8; end
        fixdist_first = fixdist.all_first{fixdist.valdiff==round(thisvaldiff)};
        fixdist_mid = fixdist.all_mid{fixdist.valdiff==round(thisvaldiff)};
        if rand <= 0.5, y = 1;
        else, y=2;
        end
        % sequence of fixations
        yseq_toMaxTime = nan(1,maxdectime/dt);
        ni = 1;
        while sum(~isnan(yseq_toMaxTime)) < maxdectime/dt
            if ni==1, itemfixdur = fixdist_first(randperm(length(fixdist_first),1));
            else
                itemfixdur = fixdist_mid(randperm(length(fixdist_mid),1));
            end
            itemfixN = round(itemfixdur/dt);
            yseq_toMaxTime(ni:ni+itemfixN-1) = y;
            ni = ni+itemfixN;
            y = 3-y;
        end
        yseq_toMaxTime = yseq_toMaxTime(1:maxdectime/dt);
        % Re-write choice
        [dstruct_addm(s).choice(ti),dstruct_addm(s).rt(ti),dstruct_addm(s).fixitem{ti},dstruct_addm(s).fixdur{ti},dstruct_addm(s).tItem(ti,:),~] = run_aDDM(dstruct_addm(s).itemval(ti,:),aGamma_k,dt,k,sig2_k,decbound,yseq_toMaxTime);
    end
end
fprintf('\n');

% Combine all datastructs to plot them efficiently
allDstructs = {dstruct_real,dstruct_opt,dstruct_addm};
titles_dstruct = {'Data','OptimalModel','aDDM'};
colors_dstruct = {'k','r','b'};
% colors_dstruct = {'k','r','m'};
marksize = 30;

% Figure S1
% Plot all three models separately
stats = struct;
stats.switchrate_rt = cell(1,length(length(allDstructs)));
stats.switchrate_vd = cell(1,length(length(allDstructs)));
stats.switchrate_vs = cell(1,length(length(allDstructs)));
stats.midfix_rt = cell(1,length(length(allDstructs)));
stats.midfix_vd = cell(1,length(length(allDstructs)));
stats.midfix_vs = cell(1,length(length(allDstructs)));
for d = 1:length(allDstructs)
    dstruct = allDstructs{d};
    for s = 1:length(dstruct)
        % compile switchrate and mid fixdur for all trials
        totaltrials = length(dstruct(s).trialnum);
        switchrate_sub = nan(totaltrials,1);
        fixdur_mid_sub = nan(totaltrials,1);
        for ti = 1:totaltrials
            numswitch = length(dstruct(s).fixdur{ti})-1;
            switchrate_sub(ti) = numswitch./dstruct(s).rt(ti);
            if numswitch > 1
                fixdur_mid_sub(ti) = mean(dstruct(s).fixdur{ti}(2:end-1));
            end
        end
        % split by RT
        rt_sub = dstruct(s).rt;
        [~,~,binassign] = histcounts(rt_sub,binedges_rt);
        for b = 1:length(binedges_rt)-1
            switchrate.rt(s,b) = nanmean(switchrate_sub(binassign==b));
            fixdur_mid.rt(s,b) = nanmean(fixdur_mid_sub(binassign==b));
        end
        % split by value difference
        valdiff_sub = getValdiff(dstruct(s));
        for vd = valdiff_unique
            switchrate.valdiff(s,vd==valdiff_unique) = nanmean(switchrate_sub(vd==valdiff_sub));
            fixdur_mid.valdiff(s,vd==valdiff_unique) = nanmean(fixdur_mid_sub(vd==valdiff_sub));
        end
        % split by value sum
        valsum_sub = getValsum(dstruct(s));
        for vs = valsum_unique
            switchrate.valsum(s,vs==valsum_unique) = nanmean(switchrate_sub(vs==valsum_sub));
            fixdur_mid.valsum(s,vs==valsum_unique) = nanmean(fixdur_mid_sub(vs==valsum_sub));
        end
    end
    % Plot switch rate
    fh1 = figure('units','normalized','outerposition',[0,0,0.7,0.5]);
    subplot(1,3,1); pbaspect([1,1,1]); hold on;
    i_use = sum(isnan(switchrate.rt),1) < (1/3)*size(switchrate.rt,1);
    [data_mean,data_se] = getMeanAndSE(switchrate.rt);
    errorbar(binctr_rt(i_use),data_mean(i_use),data_se(i_use),'.-','color',colors_dstruct{d},'markersize',marksize);
    set(gca,'xlim',[binedges_rt(1),binedges_rt(end)]);
    ylabel('Switch rate (s^{-1})'); xlabel('RT');
    % Stats
    stats.switchrate_rt{d} = getSlopeStats(binctr_rt,switchrate.rt);
    title(sprintf('p = %.04f',stats.switchrate_rt{d}.p));
    
    subplot(1,3,2); pbaspect([1,1,1]); hold on;
    i_use = sum(isnan(switchrate.valdiff),1) < (1/3)*size(switchrate.valdiff,1);
    [data_mean,data_se] = getMeanAndSE(switchrate.valdiff);
    xaxis = valdiff_unique(i_use);
    errorbar(xaxis,data_mean(i_use),data_se(i_use),'.-','color',colors_dstruct{d},'markersize',marksize);
    set(gca,'xlim',[xaxis(1)-1,xaxis(end)+1],'xtick',xaxis);
    ylabel('Switch rate (s^{-1})'); xlabel({'Abs(value difference)';'(item 1 - item 2)'});
    % Stats
    stats.switchrate_vd{d} = getSlopeStats(valdiff_unique,switchrate.valdiff);
    title(sprintf('p = %.04f',stats.switchrate_vd{d}.p));
    
    subplot(1,3,3); pbaspect([1,1,1]); hold on;
    i_use = sum(isnan(switchrate.valsum),1) < (1/3)*size(switchrate.valsum,1);
    [data_mean,data_se] = getMeanAndSE(switchrate.valsum);
    xaxis = valsum_unique(i_use);
    errorbar(xaxis,data_mean(i_use),data_se(i_use),'.-','color',colors_dstruct{d},'markersize',marksize);
    set(gca,'xlim',[-1,xaxis(end)+1],'xtick',xaxis);
    ylabel('Switch rate (s^{-1})'); xlabel('Value sum');
    % Stats
    stats.switchrate_vs{d} = getSlopeStats(valsum_unique,switchrate.valsum);
    title(sprintf('p = %.04f',stats.switchrate_vs{d}.p));
    
    % Plot fixation duration
    fh2 = figure('units','normalized','outerposition',[0,0,0.7,0.5]);
    subplot(1,3,1); pbaspect([1,1,1]); hold on;
    i_use = sum(isnan(fixdur_mid.rt),1) < (1/3)*size(fixdur_mid.rt,1);
    [data_mean,data_se] = getMeanAndSE(fixdur_mid.rt);
    errorbar(binctr_rt(i_use),data_mean(i_use),data_se(i_use),'.-','color',colors_dstruct{d},'markersize',marksize);
    set(gca,'xlim',[binedges_rt(1),binedges_rt(end)]);
    ylabel('Middle fixation duration (s)'); xlabel('RT (s)');
    % Stats
    stats.midfix_rt{d} = getSlopeStats(binctr_rt,fixdur_mid.rt);
    title(sprintf('p = %.04f',stats.midfix_rt{d}.p));
    
    subplot(1,3,2); pbaspect([1,1,1]); hold on;
    i_use = sum(isnan(fixdur_mid.valdiff),1) < (1/3)*size(fixdur_mid.valdiff,1);
    [data_mean,data_se] = getMeanAndSE(fixdur_mid.valdiff);
    xaxis = valdiff_unique(i_use);
    errorbar(xaxis,data_mean(i_use),data_se(i_use),'.-','color',colors_dstruct{d},'markersize',marksize);
    set(gca,'xlim',[xaxis(1)-1,xaxis(end)+1],'xtick',xaxis);
    ylabel('Middle fixation duration (s)'); xlabel({'Abs(value difference)';'(item 1 - item 2)'});
    % Stats
    stats.midfix_vd{d} = getSlopeStats(valdiff_unique,fixdur_mid.valdiff);
    title(sprintf('p = %.04f',stats.midfix_vd{d}.p));
    
    subplot(1,3,3); pbaspect([1,1,1]); hold on;
    i_use = sum(isnan(fixdur_mid.valsum),1) < (1/3)*size(fixdur_mid.valsum,1);
    [data_mean,data_se] = getMeanAndSE(fixdur_mid.valsum);
    xaxis = valsum_unique(i_use);
    errorbar(xaxis,data_mean(i_use),data_se(i_use),'.-','color',colors_dstruct{d},'markersize',marksize);
    set(gca,'xlim',[-1,xaxis(end)+1],'xtick',xaxis);
    ylabel('Middle fixation duration (s)'); xlabel('Value sum');
    % Stats
    stats.midfix_vs{d} = getSlopeStats(valsum_unique,fixdur_mid.valsum);
    title(sprintf('p = %.04f',stats.midfix_vs{d}.p));
end


% Figure S2
RTcutoff = 5;  % Cutoff of RT for trials
switchBinaryMat = struct;  % binary matrix with 1's at time point when switch ocurred
switchBinaryMat.all = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
switchBinaryMat.rt_low = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
switchBinaryMat.rt_hi = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
switchBinaryMat.valsum_low = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
switchBinaryMat.valsum_hi = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
switchBinaryMat.valdiff_low = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
switchBinaryMat.valdiff_hi = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
switchBinaryMat.logreg = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
% Fixation duration - store the fixation duration, stored at fixation onset
% time
fixdurMat = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
fixdur_all = {};  % Accumulate all middle fixation durations
% Switch number (normalized) for valdiff & valsum
switchnum_valdiff = {}; switchnum_valsum = {};
% Switch rate (#switches/rt) for valdiff & valsum
switchrate_valdiff = {}; switchrate_valsum = {};
% x-axis
xaxis_orig = 0:params.dt:RTcutoff-params.dt;

% Get switch probability and fixation duration across time
for d = 1:length(allDstructs)
    dstruct = allDstructs{d};
    fixdur_all{d} = [];
    % Switch proportion for value diff & value sum
    itemval_all = [];
    for s = 1:length(dstruct)
        itemval_all = cat(1,itemval_all,dstruct(s).itemval);
    end
    valdiff_unique = unique(abs(itemval_all(:,1)-itemval_all(:,2)))';
    valsum_unique = unique(itemval_all(:,1)+itemval_all(:,2))';
    switchnum_valdiff_d = nan(length(dstruct),length(valdiff_unique));
    switchnum_valsum_d = nan(length(dstruct),length(valsum_unique));
    switchrate_valdiff_d = nan(length(dstruct),length(valdiff_unique));
    switchrate_valsum_d = nan(length(dstruct),length(valsum_unique));
    for s = 1:length(dstruct)
        % All trials
        numtrials = length(dstruct(s).fixdur);
        switchBinaryMat_sub = zeros(numtrials,RTcutoff/params.dt);
        fixdurMat_sub = nan(numtrials,RTcutoff/params.dt);
        for t = 1:length(dstruct(s).fixdur)
            fixdur_trial = dstruct(s).fixdur{t};
            % Switch probability
            switchTimeIndex = round(round(cumsum(fixdur_trial),2)/params.dt) + 1;  % +1 to include RT of 0
            switchTimeIndex(switchTimeIndex>RTcutoff/params.dt) = [];  % remove any time points greater than maxdectime
            if excludeLastFix == 1 && ~isempty(switchTimeIndex)
                switchTimeIndex(end) = [];
                fixdur_trial(end) = [];
            end
            if ~isempty(switchTimeIndex)
                switchBinaryMat_sub(t,switchTimeIndex) = 1;
            end
            
            % Fixation duration
            fixonsetTimeIndex = [1,switchTimeIndex(1:end-1)];
            if length(fixonsetTimeIndex) > 1
                fixdur_this = fixdur_trial(1:length(switchTimeIndex));
                fixdur_this(1) = [];  % remove first fixation, just look at middle fixations
                fixonsetTimeIndex(1) = [];
                fixdurMat_sub(t,fixonsetTimeIndex) = fixdur_this;
            end
        end
        switchBinaryMat.all(s,:,d) = mean(switchBinaryMat_sub,1);
        % Split trials by upper third and lower third in terms of RT, valsum, and valdiff
        % 1. Split RT
        i_lowerthird = dstruct(s).rt < prctile(dstruct(s).rt,33);
        i_upperthird = dstruct(s).rt > prctile(dstruct(s).rt,66);
        switchBinaryMat.rt_low(s,:,d) = mean(switchBinaryMat_sub(i_lowerthird,:),1);
        switchBinaryMat.rt_hi(s,:,d) = mean(switchBinaryMat_sub(i_upperthird,:),1);
        % 2. Split valsum
        valsum = sum(dstruct(s).itemval,2);
        i_lowerthird = valsum < prctile(valsum,33);
        i_upperthird = valsum > prctile(valsum,66);
        switchBinaryMat.valsum_low(s,:,d) = mean(switchBinaryMat_sub(i_lowerthird,:),1);
        switchBinaryMat.valsum_hi(s,:,d) = mean(switchBinaryMat_sub(i_upperthird,:),1);
        % 3. Split valdiff
        valdiff = abs(dstruct(s).itemval(:,1)-dstruct(s).itemval(:,2));
        i_lowerthird = valdiff < prctile(valdiff,33);
        i_upperthird = valdiff > prctile(valdiff,66);
        switchBinaryMat.valdiff_low(s,:,d) = mean(switchBinaryMat_sub(i_lowerthird,:),1);
        switchBinaryMat.valdiff_hi(s,:,d) = mean(switchBinaryMat_sub(i_upperthird,:),1);
        
        % Fixation duration
        fixdurMat(s,:,d) = nanmean(fixdurMat_sub,1);
        fixdur_all{d} = cat(1,fixdur_all{d},fixdurMat_sub(~isnan(fixdurMat_sub)));
        
        % 2. Switch proportion vs valdiff
        % Switch rate: #switches / RT
        switchnum = cellfun('length',dstruct(s).fixdur) - 1;
        valdiff_sub = abs(dstruct(s).itemval(:,1)-dstruct(s).itemval(:,2));
        for vd = valdiff_unique
            i_vd = valdiff_sub==vd;
            switchnum_valdiff_d(s,vd == valdiff_unique) = nanmean(switchnum(i_vd)) / mean(switchnum);  % Normalize by mean number of all switches for subject
            rt_without_lastfix = dstruct(s).rt(i_vd) - cellfun(@(v) v(end), dstruct(s).fixdur(i_vd));
            switchrate_valdiff_d(s,vd == valdiff_unique) = nanmean(switchnum(i_vd)./rt_without_lastfix);
        end
        % 3. Switch prop vs valsum
        valsum_sub = dstruct(s).itemval(:,1)+dstruct(s).itemval(:,2);
        for vs = valsum_unique
            i_vs = valsum_sub==vs;
            switchnum_valsum_d(s,vs == valsum_unique) = nanmean(switchnum(i_vs)) / mean(switchnum);  % Normalize by mean number of all switches for subject
            rt_without_lastfix = dstruct(s).rt(i_vs) - cellfun(@(v) v(end), dstruct(s).fixdur(i_vs));
            switchrate_valsum_d(s,vs == valdiff_unique) = nanmean(switchnum(i_vs)./rt_without_lastfix);
        end
    end
    switchnum_valdiff{d} = switchnum_valdiff_d;
    switchnum_valsum{d} = switchnum_valsum_d;
    switchrate_valdiff{d} = switchrate_valdiff_d;
    switchrate_valsum{d} = switchrate_valsum_d;
end

% Create time bins for switch prob data for smoother results
if useBinnedX == 1
    switchBinaryMat_orig = switchBinaryMat;
    binSize = 4;
    i_bin = discretize(xaxis_orig,length(xaxis_orig)/binSize);
    % Get new x-axis, using bin mean
    xaxis_binned = nan(1,length(unique(i_bin)));
    for b = 1:max(i_bin)
        xaxis_binned(b) = mean(xaxis_orig(b==i_bin));
    end
    % Bin data
    for ff = fieldnames(switchBinaryMat)'
        tempMat = nan(size(switchBinaryMat.(char(ff))(:,:,d),1),max(i_bin),size(switchBinaryMat.(char(ff)),3));
        for d = 1:length(allDstructs)
            for b = 1:max(i_bin)
                tempMat(:,b,d) = mean(switchBinaryMat.(char(ff))(:,b==i_bin,d),2);
            end
        end
        switchBinaryMat.(char(ff)) = tempMat;
    end
    xaxis = xaxis_binned;
else
    xaxis = xaxis_orig;
end

% Plot results

% Switch probability
ylimit_switchprob = [0,0.15];
xlimit_switchprob = [0,RTcutoff];
for d = 1:length(allDstructs)
    fh = figure('units','normalized','outerposition',[0,0,1,0.4]);
    % All trials
    subplot(1,4,1); pbaspect([1,1,1]); hold on;
    [data_mean,data_se] = getMeanAndSE(switchBinaryMat.all(:,:,d));
    errorbar(xaxis,data_mean,data_se,'Color',colors_dstruct{d});
    set(gca,'ylim',ylimit_switchprob,'xlim',xlimit_switchprob);
    xlabel('Time (s)'); ylabel('Switch probability');
    title(titles_dstruct{d});
    
    % All trials, predicted by logistic regression
    %     subplot(1,5,2); hold on;
    %     plot(xaxis,mean(switchBinaryMat.logreg(:,:,d),1),'Color',colors_dstruct{d});
    %     title(titles_dstruct{d});
    % Split by RT
    subplot(1,4,2); pbaspect([1,1,1]); hold on;
    [mean_low,se_low] = getMeanAndSE(switchBinaryMat.rt_low(:,:,d));
    % Show where the low-RT trials cut off in time
    lowrt_zero = mean_low==0; lowrt_zero(xaxis<1) = false;
    lowrt_cutoff = xaxis(find(lowrt_zero,1));
    ax_low = errorbar(xaxis,mean_low,se_low,'g');
    [mean_hi,se_hi] = getMeanAndSE(switchBinaryMat.rt_hi(:,:,d));
    ax_hi = errorbar(xaxis,mean_hi,se_hi,'m');
    plot([lowrt_cutoff,lowrt_cutoff],ylimit_switchprob,'k--');
    set(gca,'ylim',ylimit_switchprob,'xlim',xlimit_switchprob);
    xlabel('Time (s)'); ylabel('Switch probability');
    title('Split by RT');
    % Stats
    [~,pvals] = ttest(switchBinaryMat.rt_low(:,:,d),switchBinaryMat.rt_hi(:,:,d));
    if useBonferroni==1
        i_sig = pvals < 0.05/size(switchBinaryMat.all,1);
    else
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals);  % FDR-corrected for multiple comparisons
        i_sig = adj_p < 0.05;
    end
    y_asterisk = 0.14;  % Where to show stars for significance
    plot(xaxis(i_sig),repmat(y_asterisk,1,sum(i_sig)),'k*');
    legend([ax_low,ax_hi],{'Lower 1/3','Upper 1/3'});
    
    % Split by valsum
    subplot(1,4,3); pbaspect([1,1,1]); hold on;
    [mean_low,se_low] = getMeanAndSE(switchBinaryMat.valsum_low(:,:,d));
    ax_low = errorbar(xaxis,mean_low,se_low,'g');
    [mean_hi,se_hi] = getMeanAndSE(switchBinaryMat.valsum_hi(:,:,d));
    ax_hi = errorbar(xaxis,mean_hi,se_hi,'m');
    set(gca,'ylim',ylimit_switchprob,'xlim',xlimit_switchprob);
    xlabel('Time (s)'); ylabel('Switch probability');
    title('Split by value sum');
    % Stats
    [~,pvals] = ttest(switchBinaryMat.valsum_low(:,:,d),switchBinaryMat.valsum_hi(:,:,d));
    if useBonferroni==1
        i_sig = pvals < 0.05/size(switchBinaryMat.all,1);
    else
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals);  % FDR-corrected for multiple comparisons
        i_sig = adj_p < 0.05;
    end
    y_asterisk = 0.14;  % Where to show stars for significance
    plot(xaxis(i_sig),repmat(y_asterisk,1,sum(i_sig)),'k*');
    legend([ax_low,ax_hi],{'Lower 1/3','Upper 1/3'});
    
    % Split by valdiff
    subplot(1,4,4); pbaspect([1,1,1]); hold on;
    [mean_low,se_low] = getMeanAndSE(switchBinaryMat.valdiff_low(:,:,d));
    ax_low = errorbar(xaxis,mean_low,se_low,'g');
    [mean_hi,se_hi] = getMeanAndSE(switchBinaryMat.valdiff_hi(:,:,d));
    ax_hi = errorbar(xaxis,mean_hi,se_hi,'m');
    set(gca,'ylim',ylimit_switchprob,'xlim',xlimit_switchprob);
    xlabel('Time (s)'); ylabel('Switch probability');
    title('Split by value difference');
    % Stats
    [~,pvals] = ttest(switchBinaryMat.valdiff_low(:,:,d),switchBinaryMat.valdiff_hi(:,:,d));
    if useBonferroni==1
        i_sig = pvals < 0.05/size(switchBinaryMat.all,1);
    else
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals);  % FDR-corrected for multiple comparisons
        i_sig = adj_p < 0.05;
    end
    y_asterisk = 0.14;  % Where to show stars for significance
    plot(xaxis(i_sig),repmat(y_asterisk,1,sum(i_sig)),'k*');
    legend([ax_low,ax_hi],{'Lower 1/3','Upper 1/3'});
end


% Fixation duration
minNumSubs = length(dstruct_opt)/3;
xlimit_fixdur = [0,1.5];
ylimit_fixdur = [0.4,1.8];
fixdur_slope_stats = {};
figure('units','normalized','outerposition',[0,0,1,0.5]);
for d = 1:length(allDstructs)
    subplot(1,3,d); pbaspect([1,1,1]); hold on;
    fixdurMat_d = fixdurMat(:,:,d);
    fixdurMat_d(:,sum(isnan(fixdurMat_d),1) > minNumSubs) = nan;
    [this_mean,this_se] = getMeanAndSE(fixdurMat_d);
    errorbar(xaxis_orig,this_mean,this_se,'.','markersize',marksize,'Color',colors_dstruct{d});
    set(gca,'ylim',ylimit_fixdur,'xlim',xlimit_fixdur);
    title(titles_dstruct{d});
    xlabel('Time within trial (s)'); ylabel('Fixation duration (s)');
    % Stats on slope
    b_all = nan(size(fixdurMat_d,1),1);
    for s = 1:size(fixdurMat_d,1)
        b = glmfit(xaxis_orig',fixdurMat_d(s,:)','normal');
        b_all(s) = b(2);
    end
    [h,p,ci,stats] = ttest(b_all);
    stats.p = p;
    fixdur_slope_stats{d} = stats;
end


% Figure 4 - Figure supplement 3A
% Check for single-fixation trials
maxFixnum = 4;
allnumfix = nan(length(dstruct),maxFixnum,3);
for d = 1:length(allDstructs)
    dstruct = allDstructs{d};
    for s = 1:numsubs
        numfix = cellfun('length',dstruct(s).fixitem);
        for fixnum = 1:maxFixnum
            allnumfix(s,fixnum,d) = mean(numfix==fixnum);
        end
    end
end
% Plot just the first fixation
figure; pbaspect([1,1,1]); hold on;
fixnum = 1;
[this_mean,this_se] = getMeanAndSE(squeeze(allnumfix(:,fixnum,:)));
for d = 1:length(allDstructs)
    errorbar(d,this_mean(d),this_se(d),'.','Color',colors_dstruct{d},'MarkerSize',marksize);
end
title(sprintf('Trials w/ : %d fixation(s)',fixnum));
set(gca,'XLim',[0,length(allDstructs)+1],'XTick',1:length(allDstructs),'XTickLabel',titles_dstruct,'XTickLabelRotation',45);
ylabel('Proportion of trials with single fixation');
% Stats
[~,p,~,stats_singleFixModelComp_models] = ttest2(allnumfix(:,1,2),allnumfix(:,1,3));
[~,p,~,stats_singleFixModelComp_human_addm] = ttest2(allnumfix(:,1,3),allnumfix(:,1,1));
[~,p,~,stats_singleFixModelComp_human_opt] = ttest2(allnumfix(:,1,2),allnumfix(:,1,1));

% Figure 4 - Figure supplement 3B
% Compare mean reward to human
rew_all = {};
fh = figure; pbaspect([1,1,1]); hold on;
for d = 1:length(allDstructs)
    perf = getPerfFromDstruct(allDstructs{d},params);
    rew_all{d} = perf.reward;
    [this_mean,this_se] = getMeanAndSE(perf.reward);
    errorbar(d,this_mean,this_se,'.','Color',colors_dstruct{d},'markersize',marksize);
end
set(gca,'xlim',[0,length(allDstructs)+1],'xtick',1:length(allDstructs),'xticklabel',titles_dstruct,'XTickLabelRotation',45);
ylabel('Mean reward');
% Compare mean reward to human data
[~,p_opt,~,stats_opt] = ttest2(rew_all{2},rew_all{1});
[~,p_addm,~,stats_addm] = ttest2(rew_all{3},rew_all{1});

% Figure 4 - Figure supplement 3C,D
%  Show fixation bias predictions w/ aDDM
fixbiasStats = {};
nbin_rt = 3;
nbin_valsum = 5;
d = 3;
plotOptions = struct;
plotOptions.saveFigures = 0;
plotOptions.saveStrAppend = titles_dstruct{d};
useFixationBins = 0;
fixbiasStats{d} = fixbiaseffects(allDstructs{d},nbin_rt,nbin_valsum,useFixationBins,plotOptions);
sgtitle(titles_dstruct{d},'fontsize',20);

% Figure 4 - Figure supplement 3E
% Show first, second, and third fixation durations for both models
maxNumFixShow = 3;
all_fixDur = nan(numsubs,maxNumFixShow,length(allDstructs));
fh = figure('units','normalized','outerposition',[0 0 0.2 1]);
for d = 1:length(allDstructs)
    for s = 1:numsubs
        for fn = 1:maxNumFixShow
            i_fn = cellfun('length',allDstructs{d}(s).fixdur)>=fn;
            all_fixDur(s,fn,d) = mean(cellfun(@(x) x(fn),allDstructs{d}(s).fixdur(i_fn)));
        end
    end
    subplot(length(allDstructs),1,d); pbaspect([1,1,1]); hold on;
    [this_mean,this_se] = getMeanAndSE(all_fixDur(:,:,d));
    errorbar(1:maxNumFixShow,this_mean,this_se,'.','MarkerSize',marksize,'Color',colors_dstruct{d});
    xlabel('Fixation number'); ylabel('Fixation duration (s)');
    set(gca,'XLim',[0,maxNumFixShow+1],'xtick',1:maxNumFixShow);
    title(titles_dstruct{d});
end

% Figure 4 - Figure supplement 3F
% Show RT distributions and CDF
allRT = {};
maxRT_show = 6;
binEdges = 0:0.1:maxRT_show;
for d = 1:length(allDstructs)
    fh = figure('units','normalized','outerposition',[0 0 0.3 0.7]);
    allRT{d} = [];
    for s = 1:numsubs
        allRT{d} = cat(1,allRT{d},allDstructs{d}(s).rt);
    end
    subplot(2,1,1); hold on;
    histogram(allRT{d},binEdges,'facecolor',colors_dstruct{d});
    xlabel('RT'); ylabel('Count');
    set(gca,'XLim',[0,maxRT_show]);
    title(titles_dstruct{d});
    
    subplot(2,1,2); hold on;
    h = cdfplot(allRT{d});
    h.Color = colors_dstruct{d};
    xlabel('RT'); ylabel('Cumulative distribution');
    set(gca,'XLim',[0,maxRT_show]);
end



