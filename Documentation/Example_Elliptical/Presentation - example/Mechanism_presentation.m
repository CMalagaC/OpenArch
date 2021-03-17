%% OpenArch - Script to present limit thrust line and associated collapse mechanism

%  Contact:
%  T. McLean, thomas.o.mclean@outlook.com
%  C. Malaga-Chuquitaype, c.malaga@imperial.com

clear
clc
close all

%% Load in result from analysis
set(0,'DefaultTextInterpreter','latex')

%self weight

load('Circ03')


hinges = Toolkit.HingeLocations(Result.geom,Result.x,Result.y,4,"extrados");

ARCH = Result.geom;

%% get the intrados and extrados of the arch 
intrados = zeros(1,length(ARCH)+1);
extrados = zeros(1,length(ARCH)+1);
thetas_int = zeros(1,length(ARCH)+1);
thetas_ext = zeros(1,length(ARCH)+1);

for i =1:length(ARCH)
    intrados(i) = ARCH(i).r(1);
    extrados(i) = ARCH(i).r(3);
    thetas_int(i) = ARCH(i).theta(1);
    thetas_ext(i) = ARCH(i).theta(3);
end
intrados(end) = ARCH(end).r(2);
extrados(end) = ARCH(end).r(4);
thetas_int(end) = ARCH(end).theta(2);
thetas_ext(end) = ARCH(end).theta(4);

[xin,yin] = Toolkit.pol2car(intrados,thetas_int);
[xex,yex] = Toolkit.pol2car(extrados,thetas_ext);



%% Split into blocks 

block1 = struct;
block1.xin = xin(1:hinges(2));
block1.yin = yin(1:hinges(2));
block1.xex = xex(1:hinges(2));
block1.yex = yex(1:hinges(2));

block2 = struct;
block2.xin = xin(hinges(2):hinges(3));
block2.yin = yin(hinges(2):hinges(3));
block2.xex = xex(hinges(2):hinges(3));
block2.yex = yex(hinges(2):hinges(3));

block3 = struct;
block3.xin = xin(hinges(3):hinges(4));
block3.yin = yin(hinges(3):hinges(4));
block3.xex = xex(hinges(3):hinges(4));
block3.yex = yex(hinges(3):hinges(4));

block4 = struct;
block4.xin = xin(hinges(4):end);
block4.yin = yin(hinges(4):end);
block4.xex = xex(hinges(4):end);
block4.yex = yex(hinges(4):end);


figure

subplot(1,2,1)
title('Limit thrust line')
axis equal tight
axis off
PlotABlock(block1)
PlotABlock(block2)
PlotABlock(block3)
PlotABlock(block4)
hold on
set(gca,'fontsize',20)
plot(Result.x,Result.y,'b','linewidth',2)
for i = 1:length(hinges)
    x = Result.x(hinges(i));
    y = Result.y(hinges(i));
    plot(x,y,'k.','markersize',14);
    plot([0,x],[0,y],'k--','linewidth',1)
end
xlim([-12,12])
ylim([0, 12])
hold off


%% collapse mechanism
subplot(1,2,2)
title('Collapse mechanism under inertial loading')
hold on
grid on
axis equal tight
axis off

%Adjust angles manually
phi = 0.1;
block1 = RotateBlock(block1,-phi,'ex');

block2 = RotateBlock(block2,phi,'in');
block2 = ShiftBlock(block2,'in',block1.xin(end),block1.yin(end));

block3 = RotateBlock(block3,-phi,'ex');
block3 = ShiftBlock(block3,'ex',block2.xex(end),block2.yex(end));

% block4 = RotateBlock(block4,phi,'in');
block4 = ShiftBlock(block4,'in',block3.xin(end),block3.yin(end));


PlotABlock(block1)
hold on
plot([0,block1.xex(1)],[0,block1.yex(1)],'k--','linewidth',1);
plot(block1.xex(1),block1.yex(1),'k.','markersize',14);

PlotABlock(block2)
hold on
plot([0,block2.xin(1)],[0,block2.yin(1)],'k--','linewidth',1);
plot(block2.xin(1),block2.yin(1),'k.','markersize',14);

PlotABlock(block3)
hold on
plot([0,block3.xex(1)],[0,block3.yex(1)],'k--','linewidth',1);
plot(block3.xex(1),block3.yex(1),'k.','markersize',14);

PlotABlock(block4)
hold on
plot([0,block4.xin(1)],[0,block4.yin(1)],'k--','linewidth',1);
plot(block4.xin(1),block4.yin(1),'k.','markersize',14);

set(gca,'fontsize',20)



 xlim([-12,12])
 ylim([0, 12])