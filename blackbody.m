clc;
clear;
close all;

%% constants
c=2.998e8;
h=6.626e-34;
k=1.381e-23;

%% blackbody spectrum

t=2500;

wl=transpose(1:1:10000);
em=1;

bb=2*h*c^2./(wl*10^-9).^5*1./(exp(h*c./(k*t*wl*10^-9))-1)*10^-9;
%bbnorm=bb/max(bb);

[rgb,xyzval]=spectrumtocolor(wl,bb,em);

%% plot
plot(wl,bb)
ax=gca;
title('Blackbody Spectrum & Color')
ylabel('Intensity (W SR^{-1} m^{-2} nm^{-1})');
xlabel('Wavelength (nm)');
ax.Color=rgb;

formatpresplot
legend off

lines = findobj(gcf, 'type', 'line'); 
for i=1:length(lines)
    set(lines(i),'Color',[0,0,0])
end
set(gca,'Color',rgb)

