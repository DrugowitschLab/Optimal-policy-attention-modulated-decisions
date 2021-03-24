function stats = ttest_full(dataArray)
%Perform t-test that returns all stats info needed for manuscript

[~,p,ci,stats] = ttest(dataArray);
stats.p = p;
stats.ci = ci;
stats.cohen_d = nanmean(dataArray)/nanstd(dataArray);

end

