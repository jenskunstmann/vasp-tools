function SetFontsInFigure(fontsize)
% changes the fonts size of all elements in the figure

set(findall(gcf,'type','axes'),'fontsize', fontsize)
set(findall(gcf,'type','text'),'fontSize', fontsize) 

%set(findall(gcf,'type','axes'),'LineWidth', 2)