function allstats = fixbiaseffects(dstruct,nbin_rt,nbin_valsum,useFixationBins,plotOptions)
% Show effect of RT (divided by number of fixations) and value sum on choice bias

allstats = struct;

% useFixationBins: For RT, use number of fixations rather than time bins


% Don't use single-fixation trials
removeSingleFixationTrials = 1;

minNumTrials = 10;
minNumSubs = 10;

% Initialize variables
% effect of rt (time bins)
fixbiascoeffs_rt = nan(length(dstruct),nbin_rt);
meanRT_rt = nan(length(dstruct),nbin_rt);
meanValsum_RT_rt = nan(length(dstruct),nbin_valsum);  % valsum of each bin to be used in the RT regression model
b_rt = nan(length(dstruct),1);
b_rt_valsum = nan(length(dstruct),1);
meannumfix_rt = nan(length(dstruct),nbin_rt);
% effect of rt (binned by number of fixations)
fixbiascoeffs_fn = nan(length(dstruct),nbin_rt);
meanRT_fn = nan(length(dstruct),nbin_rt);
meanValsum_RT_fn = nan(length(dstruct),nbin_valsum);  % valsum of each bin to be used in the RT regression model
b_fn = nan(length(dstruct),1);
b_fn_valsum = nan(length(dstruct),1);
% effect of valsum
fixbiascoeffs_valsum = nan(length(dstruct),nbin_valsum);
meanValsum = nan(length(dstruct),nbin_valsum);
meanRT_valsum = nan(length(dstruct),nbin_valsum);  % RT of each bin to be used in the valsum regression model
b_valsum = nan(length(dstruct),1);
b_valsumRT = nan(length(dstruct),1);
for s = 1:length(dstruct)
    % 0. rt (time bins)
    [~,~,index_allBins] = histcounts(dstruct(s).rt,[-inf,quantile(dstruct(s).rt,nbin_rt-1),inf]);
    for bin = 1:nbin_rt
        i_bin = index_allBins==bin;
        % make dstruct with just this number of fixations
        if sum(i_bin) >= minNumTrials
            dstruct_rt = struct;
            for ff = fieldnames(dstruct)'
                if ~isempty(dstruct(s).(char(ff))), dstruct_rt.(char(ff)) = dstruct(s).(char(ff))(i_bin,:); end
            end
            bout_rt = getbehavoutput(dstruct_rt);
            fixbiascoeffs_rt(s,bin) = bout_rt.item1chosen_timeadv_beta;
            meanValsum_RT_rt(s,bin) = mean(sum([dstruct_rt.itemval],2),1);
            meanRT_rt(s,bin) = mean([dstruct_rt.rt]);
            meannumfix_rt(s,bin) = mean(cellfun('length',dstruct_rt.fixitem));
        end
    end
    % Compute slope
    i_use_rt = ~isnan(meanRT_rt(s,:)) & ~isnan(fixbiascoeffs_rt(s,:));
    b = glmfit(meanRT_rt(s,i_use_rt)',fixbiascoeffs_rt(s,i_use_rt)');
    b_rt(s) = b(2);
    % regression coeff after accounting for value sum
    b = glmfit([meanRT_rt(s,i_use_rt)',meanValsum_RT_rt(s,i_use_rt)'],fixbiascoeffs_rt(s,i_use_rt)');
    b_rt_valsum(s) = b(3);
    
    % 1. rt (fixation number)
    for fn = 1:nbin_rt
        i_bin = cellfun('length',dstruct(s).fixitem)==fn;
        % make dstruct with just this number of fixations
        if sum(i_bin) >= minNumTrials
            dstruct_fn = struct;
            for ff = fieldnames(dstruct)'
                if ~isempty(dstruct(s).(char(ff))), dstruct_fn.(char(ff)) = dstruct(s).(char(ff))(i_bin,:); end
            end
            bout_fn = getbehavoutput(dstruct_fn);
            fixbiascoeffs_fn(s,fn) = bout_fn.item1chosen_timeadv_beta;
            meanValsum_RT_fn(s,fn) = mean(sum([dstruct_fn.itemval],2),1);
            meanRT_fn(s,fn) = mean([dstruct_fn.rt]);
        end
    end
    % Compute slope
    i_use_fn = ~isnan(meanRT_fn(s,:)) & ~isnan(fixbiascoeffs_fn(s,:));
    if removeSingleFixationTrials==1, i_use_fn(1) = false; end
    b = glmfit(meanRT_fn(s,i_use_fn)',fixbiascoeffs_fn(s,i_use_fn)');
    b_fn(s) = b(2);
    % regression coeff after accounting for value sum
    b = glmfit([meanRT_fn(s,i_use_fn)',meanValsum_RT_fn(s,i_use_fn)'],fixbiascoeffs_fn(s,i_use_fn)');
    b_fn_valsum(s) = b(3);
    
    % 2. value sum
    valsum_all = sum([dstruct(s).itemval],2);
    valsum_bins = discretize(valsum_all,[-inf,quantile(valsum_all,nbin_valsum-1),inf]);
    for b = 1:nbin_valsum
        i_bin = valsum_bins==b;
        dstruct_valsum = struct;
        for ff = fieldnames(dstruct)'
            if ~isempty(dstruct(s).(char(ff))), dstruct_valsum.(char(ff)) = dstruct(s).(char(ff))(i_bin,:); end
        end
        bout_valsum = getbehavoutput(dstruct_valsum);
        fixbiascoeffs_valsum(s,b) = bout_valsum.item1chosen_timeadv_beta;
        meanValsum(s,b) = mean(sum([dstruct_valsum.itemval],2),1);
        meanRT_valsum(s,b) = mean([dstruct_valsum.rt]);
    end
    % Compute slope
    % Regression using both value sum and fixation number
    i_use_valsum = ~isnan(meanValsum(s,:)) & ~isnan(fixbiascoeffs_valsum(s,:));
    if removeSingleFixationTrials==1, i_use_valsum(1) = false; end
    b = glmfit([meanValsum(s,i_use_valsum)',meanRT_valsum(s,i_use_valsum)'],fixbiascoeffs_valsum(s,i_use_valsum)');
    b_valsum(s) = b(2);
    b_valsumRT(s) = b(3);
end


% Stats compare slopes versus zero
allstats.rt = ttest_full(b_rt);
allstats.rt_valsum = ttest_full(b_rt_valsum);
p_rt = allstats.rt.p;

allstats.fixnum = ttest_full(b_fn);
allstats.fixnum_valsum = ttest_full(b_fn_valsum);
p_fn = allstats.fixnum.p;

allstats.valsum = ttest_full(b_valsum);
allstats.valsumRT = ttest_full(b_valsumRT);
p_valsum = allstats.valsum.p;

figAspect = [1,0.8,1];
figMarkerSize = 30;

% Plot
% fh = figure('unit','normalized','outerposition',[1,1,1,1]);
fh = figure;
if useFixationBins==0
    meannumfix = mean(meannumfix_rt,1);
    % 1. effect of rt (fixation number) on choice bias
    i_use = sum(~isnan(fixbiascoeffs_rt),1) >= minNumSubs;
    [mean_y,se_y] = getMeanAndSE(fixbiascoeffs_rt);
    [mean_x,se_x] = getMeanAndSE(meanRT_rt);
    subplot(1,2,1); hold on;
    errorbar(mean_x(i_use),mean_y(i_use),se_y(i_use),'k.','MarkerSize',figMarkerSize);
    errorbar(mean_x(i_use),mean_y(i_use),se_x(i_use),'horizontal','k.','Marker','none');
    xtix = round(mean_x(i_use),2);
    xtix_numfix = round(meannumfix(i_use),2);
    set(gca,'xtick',xtix,'xlim',[xtix(1)-1,xtix(end)+1],'xticklabel',xtix_numfix,'xticklabelrotation',45);
    % Show number of fixations on x axis
%     set(gca,'xtick',xtix,'xlim',[xtix(1)-1,xtix(end)+1],'xticklabel',xtix_numfix,'xticklabelrotation',45);
    title(sprintf('p = %.03f',p_rt));
    ylabel('Choice bias coefficients'); xlabel('RT [s]');
    pbaspect(figAspect);
elseif useFixationBins==1
    % 1. effect of rt (fixation number) on choice bias
    i_use = sum(~isnan(fixbiascoeffs_fn),1) >= minNumSubs;
    if removeSingleFixationTrials==1, i_use(1) = false; end
    [mean_y,se_y] = getMeanAndSE(fixbiascoeffs_fn);
    [mean_x,se_x] = getMeanAndSE(meanRT_fn);
    subplot(1,2,1); hold on;
    errorbar(mean_x(i_use),mean_y(i_use),se_y(i_use),'k.','MarkerSize',figMarkerSize);
    errorbar(mean_x(i_use),mean_y(i_use),se_x(i_use),'horizontal','k.','Marker','none');
    xtix = round(mean_x(i_use),2);
    set(gca,'xtick',xtix,'xlim',[xtix(1)-1,xtix(end)+1],'xticklabelrotation',45);
    title(sprintf('p = %.03f',p_fn));
    ylabel('Choice bias coefficients'); xlabel('RT [s] (Binned by fixations)');
    pbaspect(figAspect);
end
% 2. effect of value sum on choice bias
i_use = sum(~isnan(fixbiascoeffs_valsum),1) >= minNumSubs;
[mean_y,se_y] = getMeanAndSE(fixbiascoeffs_valsum);
[mean_x,se_x] = getMeanAndSE(meanValsum);
subplot(1,2,2); hold on;
errorbar(mean_x(i_use),mean_y(i_use),se_y(i_use),'k.','MarkerSize',figMarkerSize);
errorbar(mean_x(i_use),mean_y(i_use),se_x(i_use),'horizontal','k.','Marker','none');
xtix = round(mean_x(i_use),2);
set(gca,'xtick',xtix,'xlim',[xtix(1)-1,xtix(end)+1],'xticklabelrotation',45);
title(sprintf('p = %.03f',p_valsum));
ylabel('Choice bias coefficients'); xlabel('Value sum');
pbaspect(figAspect);

if plotOptions.saveFigures == 1
    print(fh,fullfile(plotOptions.figsavedir,sprintf('fixbias_predictions%s.svg',plotOptions.saveStrAppend)),'-dsvg','-r200');
end



end

