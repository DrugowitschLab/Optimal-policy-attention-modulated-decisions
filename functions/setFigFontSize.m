function setFigFontSize(fontSize,figureHandle)
% Sets the font size of all text in a currently open figure
    set(gca,'FontSize',fontSize);
%     figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'FontSize',fontSize);
end

