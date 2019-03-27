function mtrdpfax(xx,yy,zz,fr)
% this function is used to plot three axis (X Y Z) on the polarization
% pattern figure and let the original point (0,0,0) be the center of the
% figure.

% plot the original polarization pattern
figure; hp=surf(xx,yy,zz,fr);
if min(min(fr))>=0
    % if polarization is all positive (for example an explosive source), use a specific colormap
    load('cmapmtrdpos2','cmapmtrdpos2');
    colormap(gca,cmapmtrdpos2);
elseif max(max(fr))<=0
    % if polarization is all negtive (for example an implosion source), use a specific colormap
    load('cmapmtrdneg','cmapmtrdneg');
    colormap(gca,cmapmtrdneg);
else
    % have both positive and negtive polarization in different direction
    load('cmapmtrdp2','cmapmtrdp2');% 256 colormap
    colormap(gca,cmapmtrdp2);
end
set(hp,'LineStyle','none');
set(hp,'FaceColor','interp');
axis equal;
xlabel('X');ylabel('Y');zlabel('Z');
% title('P-wave Polarization Pattern');

%ho1=gcf;
%ho2=figure;
%objects=allchild(ho1);
%copyobj(get(ho1,'children'),ho2);

% decide the maximum value of the three axis (keep constent to three axis)
mx=max(max(abs(xx)));
my=max(max(abs(yy)));
mz=max(max(abs(zz)));
alim=max([mx, my, mz]);
alim=1.5*alim; % slightly bigger than the maximum value of the original data

% remove the original axis
axis off; hold on;

% plot the three axis
% with negtive axis
quiver3(-alim,0,0,2*alim,0,0,'Color','k','LineWidth',2,'MaxHeadSize',0.2,'AutoScale','off'); hold on;
quiver3(0,-alim,0,0,2*alim,0,'Color','k','LineWidth',2,'MaxHeadSize',0.2,'AutoScale','off'); hold on;
quiver3(0,0,-alim,0,0,2*alim,'Color','k','LineWidth',2,'MaxHeadSize',0.2,'AutoScale','off'); hold on;
% without negtive axis
%quiver3(0,0,0,alim,0,0,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on;
%quiver3(0,0,0,0,alim,0,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on;
%quiver3(0,0,0,0,0,alim,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on;

% add text label for the three axis
axtxv=1.1*alim; % control the position of the text label, slightly bigger than the axis limit
text(axtxv,0,0,'X','HorizontalAlignment','center');
text(0,axtxv,0,'Y','HorizontalAlignment','center');
text(0,0,axtxv,'Z','HorizontalAlignment','center');
hold off;

end