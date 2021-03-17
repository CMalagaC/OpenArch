function PlotABlock(block)
xin = block.xin;
yin = block.yin ;
xex = block.xex ;
yex = block.yex ;

hold on
plot(xin,yin,'k','linewidth',1.5);
plot(xex,yex,'k','linewidth',1.5);

for i = 1:length(xin)
    plot([xin(i),xex(i)],[yin(i),yex(i)],'k','linewidth',0.5);
end

plot([xin(1) xex(1)],[yin(1) yex(1)],'k','linewidth',1.5);
plot([xin(end) xex(end)],[yin(end) yex(end)],'k','linewidth',1.5);
hold off
end