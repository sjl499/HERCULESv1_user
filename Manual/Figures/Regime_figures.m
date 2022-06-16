%script to draw illustrative diagram for email to Prof Hubbard 14/9/15
clear all
close all

%define the dimensions of the ellipses
b=[1.5,1.2];
a=[2,2.4];

%define the radius of the circle
R=1.7;

%plotting 
Nplot=10000;
limits=[-3,3];

%label offsets
xoffset=0.1;
yoffset=0.15;

%find the intersection point(s) with inner ellipse
intery=sqrt((1-(R/a(1))^2)/((1/b(1))^2 - (1/a(1))^2));
interx=sqrt(R^2-intery^2);

%plot the ellipses and circles
% xplot1=linspace(-a(1),a(1), Nplot);
% xplot2=linspace(-a(2),a(2), Nplot);
% xplotR=linspace(-R,R,Nplot);
% yplot1=b(1)*sqrt(1-(xplot1/a(1)).^2);
% yplot2=b(2)*sqrt(1-(xplot2/a(2)).^2);
% yplotR=sqrt(R^2-(xplotR).^2);
% figure(1)
% plot([xplot1, fliplr(xplot1)], [-yplot1, yplot1], 'k-', 'linewidth', 2)
% hold on
% plot([xplotR, fliplr(xplotR)], [-yplotR, yplotR], 'k:', 'linewidth', 2)
% TeXString = texlabel('O');
% text(-2.5*xoffset,-yoffset,TeXString,'FontSize',14)
% TeXString = texlabel('a');
% text(-a(1)+xoffset,-yoffset,TeXString,'FontSize',14)
% TeXString = texlabel('b');
% text(-2*xoffset,b(1)-yoffset,TeXString,'FontSize',14)
% plot([0,0], [-1.1*a(1),1.1*a(1)], 'k-', 'linewidth', 1)
% plot([-1.1*a(1),1.1*a(1)],[0,0], 'k-', 'linewidth', 1)
% xlim(limits)
% ylim(limits)
% axis square
% box off
% set(gca,'Visible','off')
% set(gca, 'XTickLabelMode', 'Manual')
% set(gca, 'XTick', [])
% set(gca, 'YTickLabelMode', 'Manual')
% set(gca, 'YTick', [])


%plot the regime
%I (R=1.0) inside the surface with r<b
% theta=pi*0.25;
% plot([0,R*cos(theta)], [0,R*sin(theta)], 'k:', 'linewidth', 2)
% plot([R*cos(theta)], [R*sin(theta)], 'k.', 'linewidth', 2)
% TeXString = texlabel('(r,mu)');
% text(R*cos(theta)+xoffset,R*sin(theta)+yoffset,TeXString,'FontSize',14)
% plot([0,R*cos(theta)], [0, R*sin(theta)], 'k.', 'markersize',12)
% print('Regime_I.eps', '-depsc2')

%II  (R=1.7) Inside the surface with b<r<a
% theta=pi*0.1;
% plot([0,R*cos(theta)], [0,R*sin(theta)], 'k:', 'linewidth', 2)
% plot([R*cos(theta)], [R*sin(theta)], 'k.', 'linewidth', 2)
% TeXString = texlabel('(r,mu)');
% text(R*cos(theta)-6*xoffset,R*sin(theta)+yoffset,TeXString,'FontSize',14)
% plot([0,interx], [0,intery], 'k-', 'linewidth', 1)
% plot([0,interx], [0,-intery], 'k-', 'linewidth', 1)
% TeXString = texlabel('A');
% text(interx+xoffset,intery+yoffset,TeXString,'FontSize',14)
% TeXString = texlabel('B');
% text(interx+xoffset,-intery-yoffset,TeXString,'FontSize',14)
% plot([0,interx,interx, R*cos(theta)], [0,-intery,intery,R*sin(theta)], 'k.', 'markersize',12)
% print('Regime_II.eps', '-depsc2')


%III  (R=1.7) Outside the surface with b<r<a (The Kong et al. 2013 case)
% theta=pi*0.45;
% plot([0,R*cos(theta)], [0,R*sin(theta)], 'k:', 'linewidth', 2)
% plot([R*cos(theta)], [R*sin(theta)], 'k.', 'linewidth', 2)
% TeXString = texlabel('(r,mu)');
% text(R*cos(theta)+xoffset,R*sin(theta)+1*yoffset,TeXString,'FontSize',14)
% plot([0,interx], [0,intery], 'k-', 'linewidth', 1)
% plot([0,interx], [0,-intery], 'k-', 'linewidth', 1)
% TeXString = texlabel('A');
% text(interx+xoffset,intery+yoffset,TeXString,'FontSize',14)
% TeXString = texlabel('B');
% text(interx+xoffset,-intery-yoffset,TeXString,'FontSize',14)
% plot([0,interx,interx, R*cos(theta)], [0,-intery,intery,R*sin(theta)], 'k.', 'markersize',12)
% print('Regime_III.eps', '-depsc2')

%IV  (R=2.1) Outside the surface with r>a
% theta=pi*0.25;
% plot([0,R*cos(theta)], [0,R*sin(theta)], 'k:', 'linewidth', 2)
% plot([R*cos(theta)], [R*sin(theta)], 'k.', 'linewidth', 2)
% TeXString = texlabel('(r,mu)');
% text(R*cos(theta)+xoffset,R*sin(theta)+yoffset,TeXString,'FontSize',14)
% plot([0,R*cos(theta)], [0, R*sin(theta)], 'k.', 'markersize',12)
% print('Regime_IV.eps', '-depsc2')

%MODEL OUTLINE
xplot1=linspace(-a(1),a(1), Nplot);
xplot2=linspace(-0.7*a(1),0.7*a(1), Nplot);
xplot3=linspace(-0.4*a(1),0.4*a(1), Nplot);
yplot1=b(1)*sqrt(1-(xplot1/a(1)).^2);
yplot2=0.7*b(1)*sqrt(1-(xplot2/(0.7*a(1))).^2);
yplot3=0.4*b(1)*sqrt(1-(xplot3/(0.4*a(1))).^2);
figure(1)
plot([xplot1, fliplr(xplot1)], [-yplot1, yplot1], 'k-', 'linewidth', 2)
hold on
plot([xplot2, fliplr(xplot2)], [-yplot2, yplot2], 'k-', 'linewidth', 2)
plot([xplot3, fliplr(xplot3)], [-yplot3, yplot3], 'k-', 'linewidth', 2)
hold on
%plot([0,0], [-1.1*a(1),1.1*a(1)], 'k-', 'linewidth', 1)
%plot([-1.1*a(1),1.1*a(1)],[0,0], 'k-', 'linewidth', 1)
TeXString = texlabel('rho_0');
text(1.60,0,TeXString,'FontSize',14)
TeXString = texlabel('rho_1');
text(1.00,0,TeXString,'FontSize',14)
TeXString = texlabel('rho_2');
text(0.20,0,TeXString,'FontSize',14)
xlim(limits)
ylim(limits)
axis square
box off
set(gca,'Visible','off')
set(gca, 'XTickLabelMode', 'Manual')
set(gca, 'XTick', [])
set(gca, 'YTickLabelMode', 'Manual')
set(gca, 'YTick', [])
print('HERCULES_cartoon.eps', '-depsc2')



