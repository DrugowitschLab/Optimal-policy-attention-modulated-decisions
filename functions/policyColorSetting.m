function policyColorSetting
%Color settings for imagesc for policies in aVBDM

colors_policy = linspecer(4);
colors_policy = colors_policy([2,1,4,3],:);
colormap(colors_policy); caxis([1,4]);
cBar = colorbar;
cBar.Ticks = 1:4;
% cBar.TickLabels = {'Choose 1','Choose 2','Accumulate','Switch'};
cBar.TickLabels = {'1','2','Ac','Sw'};
end

