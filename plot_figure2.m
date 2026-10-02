clear
load G:/PCbackup/Northeast_US/data/bathymetry/NEUS_etopo1_twosatgrid.mat
lonb=lon; latb=lat;
load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/fitEta_model/prepareIAPERA5/IAPstericslaERA5_1deg_notrendGMSL.mat lon lat t sla
lon_ssh=lon; lat_ssh=lat; t_ssh=t;
% Rossby wave speed
% negative is westward, positive is eastward
load G:/PCbackup/Northeast_US/data/rossbyspeed/rossrad.mat cR lat lon
cR2=cR; lonc=lon; latc=lat;
clear cR lon lat
lonc(lonc>180)=lonc(lonc>180)-360;
ilat=latc>=20&latc<=65; ilon=lonc>=-80&lonc<=-10;
latc=latc(ilat); lonc=lonc(ilon);
cR2=squeeze(nanmean(cR2(ilon,ilat),1)); % m/s

load G:/CUdesktop/Northeast_US/sl_extreme/climatemode_relation/cr0mxcrpos_decsltgMABSABlp8yr_noIBwindGMSLtrend_IAPsteric19502023.mat lon lat mxcr
% load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/fitEta_model/prepareIAPERA5/IAPstericslaERA5_1deg_notrendGMSL.mat sla t
% etaobs=sla;
load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/fitEta_model/IAPstericERA5/fitEtanotrendGMSLlp8yr_IAPERA5windbuoyc_pfromfitEtav3.mat
load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/reg_openocean/openoceaneffect_iapstericsladecnoIBwindGMSL_19502023.mat
sltg1open0=sltg1open; sltg2open0=sltg2open; 
load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/reg_openocean/openoceaneffect_iapstericsladecnoIBwindGMSL_RW_total_19502023.mat
sltg1open_total=sltg1open; sltg2open_total=sltg2open; 

ok=squeeze(~isnan(etamod(:,:,100)));
topdb=interp2(lonb',latb,topdb',lon',lat)';
mxcr(repmat(topdb,[2,1,1])>-100)=nan;

Nlon=length(lon); Nlat=length(lat); Nt=length(t);
%calculate correlation matrix
cr0=nan([Nlon Nlat]);
pval0=nan([Nlon Nlat]);
for m=1:Nlon
    for n=1:Nlat
        if sum(~isnan(squeeze(etaobs(m,n,:))))>0&&sum(~isnan(squeeze(etamod(m,n,:))))>0
            [cr0(m,n),pval0(m,n)]=cal_corr(t,squeeze(etaobs(m,n,:)),squeeze(etamod(m,n,:)));
        end
    end
end
% VARIANCE RATIO MAP
stdratio=squeeze(nanstd(etamod,1,3))./squeeze(nanstd(etaobs,1,3)); %
stdratio(~ok)=nan;
%  
sltgrem1=sltg1o; sltgrem2=sltg2o;
co1=round(corr(sltgrem1,sltg1open0,'rows','pairwise'),2);
co2=round(corr(sltgrem2,sltg2open0,'rows','pairwise'),2);
so1=round((1-nanvar(sltgrem1-sltg1open0)./nanvar(sltgrem1))*100,1);
so2=round((1-nanvar(sltgrem2-sltg2open0)./nanvar(sltgrem2))*100,1);
co1_total=round(corr(sltgrem1,sltg1open_total,'rows','pairwise'),2);
co2_total=round(corr(sltgrem2,sltg2open_total,'rows','pairwise'),2);
so1_total=round((1-nanvar(sltgrem1-sltg1open_total)./nanvar(sltgrem1))*100,1);
so2_total=round((1-nanvar(sltgrem2-sltg2open_total)./nanvar(sltgrem2))*100,1);

% latitude and longitide of Hovmuller diagram
lat_ssh1=47.5;
lat_ssh2=33.5; 
lon_ssh1=[-55 -37];
lon_ssh2=[-77 -57];
% Killworth et al., 1997. The Speed of Observed and Theoretical Long Extratropical Planetary Waves
cR2_latssh1=interp1(latc,cR2(:),lat_ssh1)*2; % 3 times standard theory
cR2_latssh2=interp1(latc,cR2(:),lat_ssh2)*2; % 2 times standard theory
Re=6371e3; %m
dx_ssh1=1*2*pi*Re*cosd(lat_ssh1)/360;
dx_ssh2=1*2*pi*Re*cosd(lat_ssh2)/360;
dt_ssh1=(lon_ssh1(end)-lon_ssh1(1))*dx_ssh1/abs(cR2_latssh1)/(3600*24*365.25); % year
dt_ssh2=(lon_ssh2(end)-lon_ssh2(1))*dx_ssh2/abs(cR2_latssh2)/(3600*24*365.25); % year
% set up the subplot position
xl=0.13; xr=0.96; yb=0.58; yt=0.96;
dxem=0.04; dyem=0.07;
Nc=2; Nr=2;
axposit1=ax_position(xl,xr,yb,yt,dxem,dyem,Nc,Nr);

xl=0.13; xr=0.89; yb=0.05; yt=0.51;
dxem=0.04; dyem=0.04;
Nc=1; Nr=2;
axposit2=ax_position(xl,xr,yb,yt,dxem,dyem,Nc,Nr);

% calcualte trend of eddy number
figure('units','points','position',[0 0 600 1100])

axes('position',axposit1(1,2,:));
plot_corr(lon,lat,cr0,pval0,mxcr)
hold on
plot(lon_ssh1, [lat_ssh1 lat_ssh1],'-g','LineWidth',6)
plot(lon_ssh2, [lat_ssh2 lat_ssh2],'-g','LineWidth',6)
text(lon_ssh1(end),lat_ssh1,[num2str(lat_ssh1),'^\circN'],'color','g','FontSize',15)
text(lon_ssh2(end),lat_ssh2,[num2str(lat_ssh2),'^\circN'],'color','g','FontSize',15)
text(0.0,1.11,'(a) r(h_{Obs},h_{Model})','FontSize',18,'Units','normalized')

axes('position',axposit1(2,2,:));
plot_std(lon,lat,stdratio)
set(gca,'clim',[0,1])
colormap(gca,'jet');
yticklabels([])
text(0.0,1.11,'(b) STD(h_{Model})/STD(h_{Obs})','FontSize',18,'Units','normalized')

axes('position',axposit1(1,1,:));
plot_lontime(lon_ssh,t_ssh,squeeze(sla(:,lat_ssh==lat_ssh2,:)))
colorbar
hold on
plot(lon_ssh2,[2011,2011-dt_ssh2],'-g','LineWidth',3)
xlim(lon_ssh2)
set(gca,'xtick',[-70 -60])
xticklabels({'70^\circW','60^\circW'})
ylabel('Time')
text(0,1.11,['(c) IAP SLAs at ', num2str(lat_ssh2),'^\circN'],'FontSize',18,'Units','normalized')

%
axes('position',axposit1(2,1,:));
plot_lontime(lon_ssh,t_ssh,squeeze(sla(:,lat_ssh==lat_ssh1,:)))
colorbar
hold on
plot(lon_ssh1,[2006,2006-dt_ssh1],'-g','LineWidth',3)
xlim(lon_ssh1)
yticklabels([])
set(gca,'xtick',[-60 -50 -40])
xticklabels({'60^\circW','50^\circW','40^\circW'})
text(1.11,-0.08,'(m)','FontSize',18,'Units','normalized')
text(0,1.11,['(d) IAP SLAs at ', num2str(lat_ssh1),'^\circN'],'FontSize',18,'Units','normalized')

% -------------time series--------------------------- 
axes('position',axposit2(1,2,:));
plot(t,sltgrem2,'k',t,sltg2open0,'b',t,sltg2open_total,'r','LineWidth',3)
xlim([1950,2024])
ylim([-90,90])
xticklabels([])
ylabel('SLAs (mm)','fontsize',18)
set(gca,'fontsize',18)
legend('\eta^R','\eta^{EST}','\eta^{EST}_{Model}',...
        'Location','NorthEast','Orientation','horizontal','NumColumns',4,'fontsize',14)
text(0.05,0.85,['r(\eta^R,\eta^{EST})=',num2str(co2), ';   r(\eta^R,\eta^{EST}_{Model})=',num2str(co2_total)],...
  'FontSize',15,'Units','normalized')
text(0.05,0.17,['s(\eta^{EST})=',num2str(so2),'%', ';  s(\eta^{EST}_{Model})=',num2str(so2_total),'%'],'FontSize',15,'Units','normalized')
text(0.0,1.09,'(e) North USEC','FontSize',18,'Units','normalized')

axes('position',axposit2(1,1,:));
plot(t,sltgrem1,'k',t,sltg1open0,'b',t,sltg1open_total,'r','linewidth',3)
xlim([1950,2024])
ylim([-90,90])
ylabel('SLAs (mm)','fontsize',18)
set(gca,'fontsize',18)
% legend('\eta^R','\eta^{EST}','\eta^{EST}_{Model}',...
%         'Location','NorthEast','Orientation','horizontal','NumColumns',4,'fontsize',18)
text(0.05,0.85,['r(\eta^R,\eta^{EST})=',num2str(co1), ';   r(\eta^R,\eta^{EST}_{Model})=',num2str(co1_total)],...
  'FontSize',15,'Units','normalized')
text(0.05,0.17,['s(\eta^{EST})=',num2str(so1),'%', ';  s(\eta^{EST}_{Model})=',num2str(so1_total),'%'],'FontSize',15,'Units','normalized')
text(0.0,1.09,'(f) South USEC', 'FontSize',18,'Units','normalized')

set(gcf,'PaperUnits','points','PaperPosition',[0 0 600 1100],'Papersize',[600 1100])
export_fig Figure2.jpg -c [nan,nan,nan,nan] -noinvert

function plot_corr(lon,lat,cr,pval,mxcr)
%cr(pval>0.1)=nan;
%cr(abs(cr)<0.3)=nan;
lon2=lon*ones(1,length(lat));
lat2=ones(length(lon),1)*lat';

imagescnan(lon,lat,cr');
set(gca,'clim',[-1,1])
set(gca,'clim',[-1,1])
colormap(gca,bluewhitered(20));colorbar
set(gca,'YDir','normal')
hold on
basemap(lat(1),lat(end)-lat(1)+1.1,lon(1),lon(end)-lon(1)+1.1,1,0.7*[1 1 1])
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)
hold on
scatter(lon2(pval<=0.1),lat2(pval<=0.1),2,'o',...
    'MarkerEdgeColor',0.7*[1 1 1],...
    'MarkerFaceColor',0.7*[1 1 1]);
%contour(lon,lat,squeeze(mxcr(1,:,:))',[0.6,0.6],'g','linewidth',3)
%contour(lon,lat,squeeze(mxcr(2,:,:))',[0.6,0.6],'c','linewidth',3)
set(gca,'ytick',[10 25 40 55],'xtick',[-80 -50 -20])
xticklabels({'80^\circW','50^\circW','20^\circW'})
yticklabels({'10^\circN','25^\circN','40^\circN','55^\circN'})
axis([-95 -5 5 65])
set(gca,'fontsize',18)
end

function plot_std(lon,lat,std)

imagescnan(lon,lat,std');
set(gca,'clim',[0,1])
colorbar
colormap(gca,'jet')
set(gca,'YDir','normal')
hold on
basemap(lat(1),lat(end)-lat(1)+1.1,lon(1),lon(end)-lon(1)+1.1,1,0.7*[1 1 1])
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)
hold on
set(gca,'ytick',[10 25 40 55],'xtick',[-80 -50 -20])
xticklabels({'80^\circW','50^\circW','20^\circW'})
yticklabels({'10^\circN','25^\circN','40^\circN','55^\circN'})
axis([-95 -5 5 65])
set(gca,'fontsize',18)
end


function plot_lontime(lon,t,y)

    contourf(lon,t,y',30,'LineColor', 'none');
    clim([-0.07 0.07])
    %ytick([2000 2010 2020])
    colormap(gca, nclCM(131))
    ylim([1990,2021])  
    ytick([1995 2005 2015])
    set(gca,'fontsize',18)    
end


function [cr0,pval0]=cal_corr(t,y1,y2)

    %y1=rmTrend(y1);
    %y2=rmTrend(y2);
    %ar1=max(armodel(y(~isnan(y1)),1),0);
    %ar2=max(armodel(y(~isnan(y2)),1),0);
    %dof_ar1=min([Nt*(1-ar1)/(1+ar1),Nt]);
    dof1=dof_autocorr(t,y1);
    dof2=dof_autocorr(t,y2);
    dof=min(dof1,dof2);
    %--------------------------
    % lag=0 correlation
    %--------------------------
    [cr0,~]=corr(y1,y2,'rows','pairwise');
    tr=abs(cr0)*sqrt((dof-2)/(1-cr0^2));
    pval0= tcdf(tr,dof,'upper');% pval is from t test

end
