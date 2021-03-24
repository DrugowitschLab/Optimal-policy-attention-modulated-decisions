function plotHandle = shadedErrorBars(xData,data_mean,data_se,plotcolor,plotlineWidth,plotlinestyle)
% Takes in the raw time series and plots it to the current figure with
% error bars (SE)
%
% xRange: x range of plot
% timeSeries (events*time)
% color: 'r' or [1,0,0]
% lineWidth: width of line (e.g. 2)
% shadedSE: if 1, use the shaded error bars, if 2, show individual bars, if 0, draw error bars as time series
% keyboard
% Calculate mean & standard error
% keyboard


plotHandle = plot(xData,data_mean,'Color',plotcolor,'LineWidth',plotlineWidth,'linestyle',plotlinestyle);
xStep = xData;
ii    = length(xStep):-1:1;
fillX = [xStep xStep(ii)];
fillY = [data_mean+data_se data_mean(ii)-data_se(ii)];
fillY(isnan(fillY)) = 1;
hFill = fill(fillX,fillY,plotcolor);
set(hFill, 'edgecolor', 'none','FaceAlpha',0.3);


% % faceAlpha = 0.3;
% % 
% % plotHandle = plot(xData,data_mean,'Color',plotcolor,'LineWidth',plotlineWidth,'linestyle',plotlinestyle);
% % 
% % %Calculate the error bars
% % uE = data_mean + data_se;
% % lE = data_mean - data_se;
% % 
% % 
% % %Make the patch (the shaded error bar)
% % yP=[lE,fliplr(uE)];
% % xP=[xData,fliplr(xData)];
% % 
% % %remove nans otherwise patch won't work
% % xP(isnan(yP))=[];
% % yP(isnan(yP))=[];
% % 
% % thispatch=patch(xP,yP,1);
% % 
% % 
% % set(thispatch,'facecolor',plotcolor, ...
% %     'edgecolor','none', ...
% %     'facealpha',faceAlpha, ...
% %     'HandleVisibility', 'off', ...
% %     'Tag', 'shadedErrorBar_patch');
% % 
% % 
% % % %Make pretty edges around the patch.
% % % thisedge(1)=plot(xData,lE,'-');
% % % thisedge(2)=plot(xData,uE,'-');
% % % 
% % % set([thisedge], 'color',linecolor, ...
% % %     'HandleVisibility','off', ...
% % %     'Tag', 'shadedErrorBar_edge');




end

