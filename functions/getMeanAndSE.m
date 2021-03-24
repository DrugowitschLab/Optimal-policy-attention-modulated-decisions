function [out_mean,out_se] = getMeanAndSE(dataMat)
%Given an array of data, get the mean and standard error
% If input is matrix, compute mean and se along first dimension

if size(dataMat,1)==1 || size(dataMat,1)==1  % vector data
    out_mean = nanmean(dataMat);
    out_se = nanstd(dataMat)/sqrt(sum(~isnan(dataMat)));
    
else  % matrix data
    out_mean = nanmean(dataMat,1);
    out_se = nanstd(dataMat,[],1)./sqrt(sum(~isnan(dataMat),1));
end

end

