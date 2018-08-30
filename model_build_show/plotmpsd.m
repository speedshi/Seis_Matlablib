function plotmpsd(thk,vp,vs,den,qp,qs,hlay)
% This function is used to display the layered.
% The function will plot P-wave velocity, S-wave velocity, density, QP and QS of
% the layered model.
% 'thk', 'vp', 'vs', 'den', 'qp' and 'qs' must be a vector and have the same length.
% INPUT--------------------------------------------------- 
% thk: thickness of each layer (m), nl*1;
% vp: P-wave velocity of each layer (m/s), nl*1;
% vs: S-wave velocity of each layer (m/s), nl*1;
% den: density of each layer (kg/m^3), nl*1;
% qp: absorbing value of P-wave of each layer, nl*1;
% qs: absorbing value of S-wave of each layer, nl*1;
% hlay: the layer that you want to height light, could be a vector.
% OUTPUT----------------------------------------------
% seven figures display the vp, vs, den, qp and qs.

if nargin<7
    hlay=0;
end

nl=max(size(thk)); % number of layers
thki=zeros(2*nl,1);
vpi=zeros(2*nl,1);
vsi=zeros(2*nl,1);
deni=zeros(2*nl,1);
qpi=zeros(2*nl,1);
qsi=zeros(2*nl,1);
sthk=sum(thk); % total depth of all layers
zc=[0 sthk];  % boundry of depth axis of the model (m)

for i=1:nl
    if (i<nl)
        thki(2*i)=thki(2*i-1)+thk(i);
        thki(2*i+1)=thki(2*i);
    else
        thki(2*i)=thki(2*i-1)+thk(i);
    end
    
    vpi(2*i-1)=vp(i);
    vpi(2*i)=vp(i);
    vsi(2*i-1)=vs(i);
    vsi(2*i)=vs(i);
    deni(2*i-1)=den(i);
    deni(2*i)=den(i);
    qpi(2*i-1)=qp(i);
    qpi(2*i)=qp(i);
    qsi(2*i-1)=qs(i);
    qsi(2*i)=qs(i);
end

figure;
xx1=vpi; xx2=vsi;
hf=plot(xx1,thki,'b',xx2,thki,'r');
set(hf,'linewidth',1.5);
xlabel('Velocity (m/s)'); ylabel('Depth (m)');
legend('Vp','Vs');
xmin=min([min(xx1) min(xx2)]);
xmax=max([max(xx1) max(xx2)]);
dxl=0.15*(xmax-xmin);
xmin=xmin-dxl;xmax=xmax+dxl;
xlim([xmin xmax]); ylim([zc(1) zc(2)]);
axis ij;

figure;
xx=vpi;
plot(xx,thki,'b','linewidth',1.5);
xlabel('Vp (m/s)');ylabel('Depth (m)');
xmin=min(xx);xmax=max(xx);
dxl=0.15*(xmax-xmin);
xmin=xmin-dxl;xmax=xmax+dxl;
xlim([xmin xmax]); ylim([zc(1) zc(2)]);
axis ij;

figure;
xx=vsi;
plot(xx,thki,'b','linewidth',1.5);
xlabel('Vs (m/s)');ylabel('Depth (m)');
xmin=min(xx);xmax=max(xx);
dxl=0.15*(xmax-xmin);
xmin=xmin-dxl;xmax=xmax+dxl;
xlim([xmin xmax]); ylim([zc(1) zc(2)]);
axis ij;

figure;
xx=deni;
plot(xx,thki,'b','linewidth',1.5);
xlabel('Density (kg/m^3)');ylabel('Depth (m)');
xmin=min(xx);xmax=max(xx);
dxl=0.15*(xmax-xmin);
xmin=xmin-dxl;xmax=xmax+dxl;
xlim([xmin xmax]); ylim([zc(1) zc(2)]);
axis ij;

figure;
xx=qpi;
plot(xx,thki,'b','linewidth',1.5);
xlabel('Qp');ylabel('Depth (m)');
xmin=min(xx);xmax=max(xx);
xmin=xmin-100;xmax=xmax+100;
xlim([xmin xmax]); ylim([zc(1) zc(2)]);
axis ij;

figure;
xx=qsi;
plot(xx,thki,'b','linewidth',1.5);
xlabel('Qs');ylabel('Depth (m)');
xmin=min(xx);xmax=max(xx);
xmin=xmin-100;xmax=xmax+100;
xlim([xmin xmax]); ylim([zc(1) zc(2)]);
axis ij;

figure;
subplot(1,3,1); xx=vpi/1000; yy=thki/1000;
plot(xx,yy,'b','linewidth',1.5); hold on;
xlabel('Vp ($$\mathrm{km}/\mathrm{s}$$)','Interpreter','latex');ylabel('Depth ($$\mathrm{km}$$)','Interpreter','latex');
xmin=min(xx);xmax=max(xx);
dxl=0.15*(xmax-xmin);
xmin=xmin-dxl;xmax=xmax+dxl;
% hight light the target layer
if hlay~=0
    for ii=1:length(hlay)
        fill([xmin xmax xmax xmin],[yy(2*hlay(ii)-1) yy(2*hlay(ii)-1) yy(2*hlay(ii)) yy(2*hlay(ii))],'r','linestyle','none','FaceAlpha',0.35); hold on;
    end
end
xlim([xmin xmax]); ylim([zc(1) zc(2)]/1000);
axis ij;
subplot(1,3,2); xx=vsi/1000; yy=thki/1000;
plot(xx,yy,'b','linewidth',1.5); hold on;
xlabel('Vs ($$\mathrm{km}/\mathrm{s}$$)','Interpreter','latex');
%xlabel('Vs ($$km/s$$)','Interpreter','latex');
xmin=min(xx);xmax=max(xx);
dxl=0.15*(xmax-xmin);
xmin=xmin-dxl;xmax=xmax+dxl;
% hight light the target layer
if hlay~=0
    for ii=1:length(hlay)
        fill([xmin xmax xmax xmin],[yy(2*hlay(ii)-1) yy(2*hlay(ii)-1) yy(2*hlay(ii)) yy(2*hlay(ii))],'r','linestyle','none','FaceAlpha',0.35); hold on;
    end
end
xlim([xmin xmax]); ylim([zc(1) zc(2)]/1000);
axis ij;
subplot(1,3,3); xx=deni/1000; yy=thki/1000;
plot(xx,yy,'b','linewidth',1.5); hold on;
xlabel('Density ($$\mathrm{g}/\mathrm{cm}^3$$)','Interpreter','latex');
xmin=min(xx);xmax=max(xx);
dxl=0.15*(xmax-xmin);
xmin=xmin-dxl;xmax=xmax+dxl;
% hight light the target layer
if hlay~=0
    for ii=1:length(hlay)
        fill([xmin xmax xmax xmin],[yy(2*hlay(ii)-1) yy(2*hlay(ii)-1) yy(2*hlay(ii)) yy(2*hlay(ii))],'r','linestyle','none','FaceAlpha',0.35); hold on;
    end
end
xlim([xmin xmax]); ylim([zc(1) zc(2)]/1000);
axis ij;


end