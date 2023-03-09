function PlotToFileColorEPS(gcf,filename,x,y)
%
% Print the current plot at size (x,y) to an eps file
set(gcf,'PaperUnits','inches')% centimeters
set(gcf,'PaperSize',[x,y])
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperPosition',[0,0,x,y]);
set(gcf, 'renderer', 'painters');
print(filename,'-depsc','-r1200')
% print(gcf,filename,'-r600','-dtiff');
end
