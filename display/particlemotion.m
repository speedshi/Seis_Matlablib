function h=particlemotion(data1,data2,data3)
% This function is used to animate the particle movement.
% Input:----------------------------------------------------
% data1: the time serial showing the particle motion in x1 direction, vector;
% data2: the time serial showing the particle motion in x2 direction, vector;
% data3: the time serial showing the particle motion in x3 direction, vector;
% data1, data2 and data3 must have the same length.

data1=data1(:);
nt=length(data1); % number of time points

figure;
h=animatedline('MaximumNumPoints',1,'Marker','o','markersize',10,'markeredgecolor','k','markerfacecolor','b');
% set the axis limit
if nargin==1
    xlim([1 nt]); ylim([min(data1) max(data1)]);
    xlabel('X1'); ylabel('Amplitude');
elseif nargin==2
    data2=data2(:); 
    axis equal;
    axis([min(data1),max(data1),min(data2),max(data2)]);
    xlabel('X1');ylabel('X2'); 
elseif nargin==3
    data2=data2(:);
    data3=data3(:);
    axis equal;
    axis([min(data1),max(data1),min(data2),max(data2),min(data3),max(data3)]);
    xlabel('X1');ylabel('X2');zlabel('X3');
end

for ii=1:nt
    if nargin==1
        addpoints(h,ii,data1(ii));
        line(1:ii,data1(1:ii),'color','k','linewidth',0.6);
    elseif nargin==2
        addpoints(h,data1(ii),data2(ii));
        line(data1(1:ii),data2(1:ii),'color','k','linewidth',0.6);
    elseif nargin==3
        addpoints(h,data1(ii),data2(ii),data3(ii));
        line(data1(1:ii),data2(1:ii),data3(1:ii),'color','k','linewidth',0.6);
    end
    drawnow;
end
clearpoints(h);
end