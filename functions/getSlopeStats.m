function stats_out = getSlopeStats(x,inputMat)
%Given a set of values, compute whether the linear slope is significant across
%participants
% stats_out = getSlopeStats(x,inputMat)
% x: x values corresponding to columns in inputMat
% inputMat: each row corresponds to data from single subject
all_b = nan(size(inputMat,1),1);
for s = 1:size(inputMat,1)
    b = regress(inputMat(s,:)',[ones(length(x),1),x']);
    all_b(s) = b(2);
end
% Perform t-test on all slopes
[h,p,ci,stats_out] = ttest(all_b);
stats_out.p = p;

end

