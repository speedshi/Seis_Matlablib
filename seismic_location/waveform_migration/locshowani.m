function locshowani(migloc,tarxr,taryr,tarzr)
% This function is used to show the location results in an animation way.
% Input:----------------------------------------------------
% migloc: event location results, nevt*5, Origin_time-X-Y-Z-Coherency.
% tarxr: coordinate limit (range) in the X direction;
% taryr: coordinate limit (range) in the Y direction;
% tarzr: coordinate limit (range) in the Z direction.

nevt=size(migloc,1); % number of seismic events

%vid=VideoWriter('locres.avi');
%open(vid);

figure;
for ii=1:nevt
    scatter(migloc(ii,3),migloc(ii,2),'filled','k');
    axis equal; xlim(taryr); ylim(tarxr); box on;
    xlabel('East (km)'); ylabel('North (km)');
    drawnow;
    %mov(ii)=getframe(gcf);
    %writeVideo(vid,mov(ii));
    pause(0.01);
end

%close(vid)

end