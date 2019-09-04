clc;
clear;
close all;

%% constants
c=2.998e8; %m/s
h=4.136e-15; %eV s
k=8.617e-5; %eV/K
q=1.602e-19; %C
T_c=300;
T_s=5800;
Omega_s=6.8e-5;
Omega_c=pi;

e_step=0.001;
egs=0:e_step:3; %bandgaps in eV

%% SQ, AM0

%incident power
e_tot=0.1:e_step:8;
BBe_s=2*Omega_s/(h^3*c^2)*e_tot.^3./(exp(e_tot/(k*T_s))-1);
pinc=trapz(e_tot,BBe_s)*q; %W/(m^2) - multiply by q for eV to J
pinc3=pinc;

%initialize efficiency 
eff_am0=zeros(1,length(egs));

%run over bandgaps, find max efficiency
for i=1:length(egs)
    
    eg=egs(i);
    e=eg:e_step:8; 
    v=0:e_step:eg; 
    
    BB_c=2*Omega_c/(h^3*c^2)*e.^2./(exp(e/(k*T_c))-1); %1/(m^2*s)
    j0=q*trapz(e,BB_c); %C/(m^2*s) = A/m^2
    
    BB_s=2*Omega_s/(h^3*c^2)*e.^2./(exp(e/(k*T_s))-1);
    jsc=q*trapz(e,BB_s);
    
    j=-(j0*(exp(v/(k*T_c))-1)-jsc); %just need v b/c k is in eV/K
    %j=-(j0*(exp(v/(k*T_c)))-jsc); %no real difference
    
    pmax=max(j.*v); 
    eff_am0(i)=pmax/pinc;
    
end


%% SQ, AM1.5

load('am1p5_raw.mat');
ev=flip(1240./am1p5_raw(:,1));
am1p5=1./(ev.^2).*flip(1240.*am1p5_raw(:,2)); %W/(m^2 eV)
G=griddedInterpolant(ev,am1p5);
pinc=trapz(ev,am1p5);
pinc2=trapz(am1p5_raw(:,1),am1p5_raw(:,2));

%initialize efficiency 
eff_am1p5=zeros(1,length(egs));

%run over bandgaps, find max efficiency
for i=1:length(egs)
    
    eg=egs(i);
    e=eg:e_step:5; 
    v=0:e_step:eg; 
    
    BB_c=2*Omega_c/(h^3*c^2)*e.^2./(exp(e/(k*T_c))-1); %1/(m^2*s)
    j0=q*trapz(e,BB_c); %C/(m^2*s) = A/m^2
    
    rad_s=G(e)./(q*e); %photons/(m^2*s*eV)
    jsc=q*trapz(e,rad_s); 
    
    %j=-(j0*(exp(v/(k*T_c))-1)-jsc); %just need v b/c k is in eV/K
    j=-(j0*(exp(v/(k*T_c)))-jsc); %no real difference
    
    pmax=max(j.*v); 
    eff_am1p5(i)=pmax/pinc;
    
    if eg==1.337
       J=j;
       V=v;
    end
    
end

%% calculate voc etc.
[c_voc, i_voc] = min(abs(transpose(J)));
voc=round(V(i_voc),3,'significant');

[c_jsc, i_jsc] = min(abs(transpose(V)));
jsc=round(J(i_jsc),3,'significant');

pce=round(max(J.*V)/pinc,3,'significant');
ff=round(max(J.*V)/(voc*jsc),3,'significant');

%% calculate contributions to efficiency

eg=1.337;
e=eg:e_step:5;

hw_avg=trapz(e,G(e))./(jsc);

n_abs=jsc*hw_avg/pinc;
n_therm=(eg+3/2*k*T_c)/hw_avg;
n_thermo=voc/(eg+3/2*k*T_c);

n=n_abs*n_therm*n_thermo*ff;

%% plots
% figure(1)
% plot(V,J/10)
% xlabel('Voltage (V)')
% ylabel('Current Density (mA/cm^2)')
% title('I-V Curve')
% ylim([0,40])
% str = {strcat('V_{oc}=',num2str(voc*1000),' mV'),...
%     strcat('J_{sc}=',num2str(jsc/10),' mA/cm^2'),...
%     strcat('FF=',num2str(ff)),...
%     strcat('PCE=',num2str(pce*100),'%')
%     };
% text(0.5,20,str,'FontSize',16);
% formatpresplot
% legend off
% 
% figure(2)
% plot(e_tot,BBe_s*q,ev,am1p5)
% xlabel('Energy (eV)')
% ylabel('Spectral Irradiance (W m^{-2} eV^{-1})')
% legend('AM0','AM1.5')
% title('Illumination')
% xlim([0,5])
% formatpresplot
% 
% figure(3)
% plot(egs,[eff_am0;eff_am1p5])
% xlabel('Bandgap (eV)')
% ylabel('Max Efficiency (%)')
% title('Detailed Balance Efficiency Limit')
% legend('AM0','AM1.5')
% formatpresplot
% 
% max(eff_am1p5)
% 
