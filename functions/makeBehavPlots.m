function allstats = makeBehavPlots(bout,plotOptions)
%Make plots for behavioral data similar to Krajbich et al., 2010

allstats = struct;  % return statistical tests done
subCount = length(bout.rt);

plot_psychometrics = 1;
plot_middleFix = 0;
plot_choiceBias = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic psychometric curves
if plot_psychometrics==1
    fh = figure('units','normalized','outerposition',[0 0 1 0.6]);
    
    % Proportion chosen item 1 vs. value difference
    toplot = bout.item1chosen_valdiff_all;
    usebin_i = sum(~isnan(toplot),1)>=size(toplot,1)/2;  % Remove columns where there are NaNs more than half of all data
    subplot(1,4,1); hold on;
    [this_mean,this_se] = getMeanAndSE(toplot);
    errorbar(bout.valdiff(usebin_i),this_mean(usebin_i),this_se(usebin_i),'k.-','markersize',30);
    set(gca,'xlim',[-7,7]);
    xlabel('Delta'); ylabel('P(chose item 1)'); pbaspect([1,1,1]);
    % Stats
    b_all = nan(subCount,1); x = bout.valdiff(usebin_i);
    for s = 1:subCount
        b = regress(toplot(s,usebin_i)',[ones(length(x),1),x']);
        b_all(s) = b(2);
    end
    allstats.prop1_valDiff = ttest_full(b_all);
    
    % RT vs. absolute value difference
    toplot = bout.rt_valdiff_abs;
    usebin_i = sum(~isnan(toplot),1)>=size(toplot,1)/2;
    subplot(1,4,2); hold on;
    [this_mean,this_se] = getMeanAndSE(bout.rt_valdiff_abs);
    errorbar(bout.valdiff_abs(usebin_i),this_mean(usebin_i),this_se(usebin_i),'k.-','markersize',30);
    set(gca,'xlim',[-1,8]);
    % set(gca,'xlim',[bout.valdiff_abs(1)-1,bout.valdiff_abs(end)+1]);
    xlabel('Delta'); ylabel('RT (s)'); pbaspect([1,1,1]);
    % Stats
    b_all = nan(subCount,1); x = bout.rt_valdiff_abs(usebin_i);
    for s = 1:subCount
        b = regress(toplot(s,usebin_i)',[ones(length(x),1),x']);
        b_all(s) = b(2);
    end
    allstats.rt_absValDiff = ttest_full(b_all);
    
    % Switch count vs absolute value difference
    for use_switchrate = [0,1]
        % Use switchcountRT to show switch rate (#switches per time)
        if use_switchrate==1
            toplot = bout.switchcountRT_valdiff_abs;
        else
            toplot = bout.switchcount_valdiff_abs;
        end
        usebin_i = sum(~isnan(toplot),1)>=size(toplot,1)/2;
        subplot(1,4,3+use_switchrate); hold on;
        [this_mean,this_se] = getMeanAndSE(toplot);
        errorbar(bout.valdiff_abs(usebin_i),this_mean(usebin_i),this_se(usebin_i),'k.-','markersize',30);
        set(gca,'xlim',[-1,8]);
        % set(gca,'xlim',[bout.valdiff_abs(1)-1,bout.valdiff_abs(end)+1]);
        xlabel('Delta');
        if use_switchrate==1
            ylabel('Number of switches per 1s');
        else
            ylabel('Number of switches');
        end
        pbaspect([1,1,1]);
    end
    if plotOptions.saveFigures == 1
        print(fh,fullfile(plotOptions.figsavedir,sprintf('psychometric_curves%s.svg',plotOptions.saveStrAppend)),'-dsvg','-r200');
    end
    % Stats
    b_all = nan(subCount,1); x = bout.rt_valdiff_abs(usebin_i);
    for s = 1:subCount
        b = regress(toplot(s,usebin_i)',[ones(length(x),1),x']);
        b_all(s) = b(2);
    end
    allstats.switchCountRT_absValDiff = ttest_full(b_all);
    
end

% Middle fixation properties
if plot_middleFix==1
    fh = figure('units','normalized','outerposition',[0 0 1 0.6]);
    toplot = bout.mid_fixitemval;
    usebin_i = sum(~isnan(toplot),1)>=size(toplot,1)/2;
    subplot(1,3,1); hold on;
    [out_mean,out_se] = getMeanAndSE(bout.mid_fixitemval);
    bar(bout.val(usebin_i),out_mean(usebin_i),'facecolor','w');
    errorbar(bout.val(usebin_i),out_mean(usebin_i),out_se(usebin_i),'.k','marker','none');
    set(gca,'xlim',[bout.val(1)-1,bout.val(end)+1]); pbaspect([1,1,1]);
    ylabel('Middle fixation duration (ms)'); xlabel('Fixated item value');
    
    toplot = bout.mid_fixunfixdiff;
    usebin_i = sum(~isnan(toplot),1)>=size(toplot,1)/2;
    subplot(1,3,2); hold on;
    [out_mean,out_se] = getMeanAndSE(bout.mid_fixunfixdiff);
    bar(bout.valdiff(usebin_i),out_mean(usebin_i),'facecolor','w');
    errorbar(bout.valdiff(usebin_i),out_mean(usebin_i),out_se(usebin_i),'.k','marker','none');
    set(gca,'xlim',[bout.valdiff(1)-1,bout.valdiff(end)+1]); pbaspect([1,1,1]);
    ylabel('Middle fixation duration (ms)'); xlabel('Fixated - unfixated item value');
    
    toplot = bout.mid_absvaldiff;
    usebin_i = sum(~isnan(toplot),1)>=size(toplot,1)/2;
    subplot(1,3,3); hold on;
    [out_mean,out_se] = getMeanAndSE(bout.mid_absvaldiff);
    bar(bout.valdiff_abs(usebin_i),out_mean(usebin_i),'facecolor','w');
    errorbar(bout.valdiff_abs(usebin_i),out_mean(usebin_i),out_se(usebin_i),'.k','marker','none');
    set(gca,'xlim',[bout.valdiff_abs(1)-1,bout.valdiff_abs(end)+1]); pbaspect([1,1,1]);
    ylabel('Middle fixation duration (ms)'); xlabel('Abs value difference');
    if plotOptions.saveFigures == 1
        %     print(fh,fullfile(plotOptions.figsavedir,sprintf('midfix%s.eps',plotOptions.saveStrAppend)),'-depsc','-r200');
        %     print(fh,fullfile(plotOptions.figsavedir,sprintf('midfix%s.tiff',plotOptions.saveStrAppend)),'-dtiff','-r100');
    end
end

% Choice biases
if plot_choiceBias==1
    vdshow_i = abs(bout.valdiff) <= 5;
    plots = [];
    fh = figure('units','normalized','outerposition',[0 0 1 0.6]);
    subplot(1,3,1); hold on;
    [out_mean,out_se] = getMeanAndSE(bout.item1chosen_valdiff_all);
    plots(1) = errorbar(bout.valdiff(vdshow_i),out_mean(vdshow_i),out_se(vdshow_i),'k.-','markersize',30);
    [out_mean,out_se] = getMeanAndSE(bout.item1chosen_valdiff_lastL);
    plots(2) = errorbar(bout.valdiff(vdshow_i),out_mean(vdshow_i),out_se(vdshow_i),'r.-','markersize',30);
    [out_mean,out_se] = getMeanAndSE(bout.item1chosen_valdiff_lastR);
    plots(3) = errorbar(bout.valdiff(vdshow_i),out_mean(vdshow_i),out_se(vdshow_i),'b.-','markersize',30);
    set(gca,'ylim',[0,1]);
    legend(plots,{'All','Last fixated left','Last fixated right'},'location','northwest');
    xlabel('Left - right value'); ylabel('p(left)');
    pbaspect([1,1,1]);
    
    % Proportion chose 1 vs. fixation time difference
    toplot = bout.item1chosen_timeadv_norm;
    usebin_i = sum(~isnan(toplot),1)>=size(toplot,1)/2;
    % usebin_i = true(size(toplot,2),1);
    subplot(1,3,2); hold on;
    [out_mean,out_se] = getMeanAndSE(bout.item1chosen_timeadv_norm);
    bar(bout.binmeans_item1timeadv(usebin_i),out_mean(usebin_i),'facecolor','w');
    errorbar(bout.binmeans_item1timeadv(usebin_i),out_mean(usebin_i),out_se(usebin_i),'k.','marker','none');
    ylimit = get(gca,'ylim');
    set(gca,'ylim',[-max(abs(ylimit)),max(abs(ylimit))]);
    set(gca,'xlim',[bout.binmeans_item1timeadv(1)-0.5,bout.binmeans_item1timeadv(end)+0.5]); pbaspect([1,1,1]);
    xlabel('Left - right fixation time'); ylabel('Normalized p(left)');
    pbaspect([1,1,1]);
    % Stats
    b_all = nan(subCount,1); x = bout.binmeans_item1timeadv(usebin_i);
    for s = 1:subCount
        b = regress(toplot(s,usebin_i)',[ones(length(x),1),x']);
        b_all(s) = b(2);
    end
    allstats.prop1_fixDiff = ttest_full(b_all);
    
    % toplot = bout.firstchosen_firstfixdur_norm;
    toplot = bout.firstchosen_firstfixdur;
    usebin_i = sum(~isnan(toplot),1)>=size(toplot,1)/2;
    subplot(1,3,3); hold on;
    [out_mean,out_se] = getMeanAndSE(bout.firstchosen_firstfixdur_norm);
    bar(bout.binmeans_firstfixdur(usebin_i),out_mean(usebin_i),'facecolor','w');
    errorbar(bout.binmeans_firstfixdur(usebin_i),out_mean(usebin_i),out_se(usebin_i),'k.','marker','none');
    set(gca,'xlim',[bout.binmeans_firstfixdur(1)-0.1,bout.binmeans_firstfixdur(end)+0.1]); pbaspect([1,1,1]);
    xlabel('First fixation duration'); ylabel('Normalized p(first item)');
    pbaspect([1,1,1]);
    if plotOptions.saveFigures == 1
        print(fh,fullfile(plotOptions.figsavedir,sprintf('fixbias%s.svg',plotOptions.saveStrAppend)),'-dsvg','-r200');
    end
end

end

