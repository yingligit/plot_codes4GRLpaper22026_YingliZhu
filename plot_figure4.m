clear
% CESM1 IO results
load G:/CUdesktop/Northeast_US/data/climatemode_2024/modeindex_rmseatrendMVLR19502024.mat
it=t>1950&t<2024; t=t(it); amm=modindex(it,14); amo=modindex(it,15); atl3=modindex(it,23); % used for era5
it1=t>1950&t<2014; t1=t(it1);amm1=amm(it1); amo1=amo(it1); atl31=atl3(it1); % used for  
[~,~,yf]=reg_model([ones(length(t1),1),amm1,atl31],amo1);
amo_ammatl31=amo1-yf;

load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/CESM1-NATL/PSL_natlemmean_19502013_rmseatrendMVLR.mat
slp_natl=slp; lon_natl=lon; lat_natl=lat; t_natl=t; lon_natl(lon_natl>180)=lon_natl(lon_natl>180)-360; clear t

load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/fitEta_model/prepareIAPERA5/dlmreg_1deg_all_NAO_AMOnoATL3AMM_ATL3_AMM.mat
lon_era=lon; lat_era=lat;
buoy_era=squeeze(dthermdlmfit(:,:,:,4)+dhalindlmfit(:,:,:,4));
load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/fitEta_model/prepareIAPERA5/dlmreg_slpuvstress_1deg_all_NAO_AMOnoATL3AMM_ATL3_AMM.mat
slp_era=slpdlmfit(:,:,:,4);
utau_era=utaudlmfit(:,:,:,4);
vtau_era=vtaudlmfit(:,:,:,4);

% load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/fitEta_model/prepareIAPERA5/IAPstericslaERA5_1deg_notrendGMSL.mat
% buoy_era=dtherm+dhalin;
% load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/fitEta_model/prepareIAPERA5/ERA5slpuvstress_1deg.mat
% slp_era=slp;
% utau_era=utau;
% vtau_era=vtau;
% lon_era=lon; lat_era=lat;

t=tmon;

it=ceil((t-floor(t))*12)>=1&ceil((t-floor(t))*12)<=12;
idx1=(amm>=nanstd(amm));
% temporal mean
utau_erax=nanmean(utau_era(:,:,it&idx1),3);
vtau_erax=nanmean(vtau_era(:,:,it&idx1),3);
slp_erax=nanmean(slp_era(:,:,it&idx1),3);
buoy_erax=nanmean(buoy_era(:,:,it&idx1),3);

idx11=(amo1>=nanstd(amo1))&(abs(amo_ammatl31)<nanstd(amo_ammatl31));
slp_natlx=nanmean(slp_natl(:,:,idx11),3);

% load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/CAM_exp/AMM/buoyc_AMM.mat
% buoy1=dtherm+dhalin;
% load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/CAM_exp/AMM/uvstress_AMM.mat
% utau1=utau; vtau1=vtau; 
% load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/CAM_exp/ctrl/buoyc_ctrl.mat
% buoy2=dtherm+dhalin;
% load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/CAM_exp/ctrl/uvstress_ctrl.mat
% utau2=utau; vtau2=vtau; 
% load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/CAM_exp/AMM/slp_AMM.mat
% slp_amm=slp;
load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/CAM_exp/AMM/slp_AMM_members.mat
%slp_amm=squeeze(mean(slp(:,:,:,[3,7,8,11,14,20,9,12,13,17,21,23,24]),4));
slp_amm=squeeze(mean(slp(:,:,:,[3,7,8,11,14,20,9,12,17,23]),4));
load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/CAM_exp/ctrl/slp_ctrl.mat
slp_ctrl=slp;
lon_cam=lon; lat_cam=lat;

isea=[1:12];
%buoy_cam=mean(buoy1(:,:,isea)-buoy2(:,:,isea),3);
slp_camx=mean(slp_amm(:,:,isea)-slp_ctrl(:,:,isea),3);
%utau_cam=mean(utau1(:,:,isea)-utau2(:,:,isea),3);
%vtau_cam=mean(vtau1(:,:,isea)-vtau2(:,:,isea),3);

Ns2=500; % multiply it with data, for the vector scale
Ws=10; % for vector scale shown in vector symbol 
Nwint2=5; % wind plot interval

% set up the subplot position
xl=0.06; xr=0.96; yb=0.1; yt=0.95;
dxem=0.04; dyem=0.09;
Nc=2; Nr=2;
axposit1=ax_position(xl,xr,yb,yt,dxem,dyem,Nc,Nr);

figure('units','points','position',[0 0 900 600])

axes('position',axposit1(1,2,:));
plot_buoywind(lon_era,lat_era,buoy_erax,lon_era,lat_era,utau_erax*Ns2,vtau_erax*Ns2,Nwint2)
rectangle('Position',[-78,20,20,20],'EdgeColor','g','LineWidth',3)
rectangle('Position',[-65,40,40,20],'EdgeColor','g','LineWidth',3)
text(-70,30,'R1','color','k','fontsize',24) % the scale value is 1/Ns*Ws
text(-50,50,'R2','color','k','fontsize',24) % the scale value is 1/Ns*Ws
xticklabels([])
hold on
text(1.02,-0.06,'(m/s)','FontSize',20,'Units','normalized')
text(-91,50,'2\times10^{-2}Nm^{-2}','color','m','fontsize',18) % the scale value is 1/Ns*Ws
quiver(-91,55,Ws,0,0,'color','m','linewidth',2,'MaxHeadSize',1)
text(0,1.08,'(a) Buoyancy and winds (BDLM_{ERA5})','FontSize',20,'Units','normalized')

axes('position',axposit1(2,2,:));
plot_slp(lon_era,lat_era,slp_erax)
xticklabels([])
yticklabels([])
text(1.02,-0.06,'(Pa)','FontSize',20,'Units','normalized')
text(0,1.08,'(b) SLP (BDLM_{ERA5})','FontSize',20,'Units','normalized')

axes('position',axposit1(1,1,:));
plot_slp(lon_natl,lat_natl,slp_natlx)
%set(gca,'clim',[-3e-10,3e-10])
text(1.02,-0.06,'(Pa)','FontSize',20,'Units','normalized')
text(0,1.08,'(c) SLP (CESM1-Atlantic)','FontSize',20,'Units','normalized')

axes('position',axposit1(2,1,:));
plot_slp(lon_cam,lat_cam,slp_camx)
yticklabels([])
hold on
text(1.02,-0.06,'(Pa)','FontSize',20,'Units','normalized')
text(0,1.08,'(d) SLP (CAM5)','FontSize',20,'Units','normalized')

set(gcf,'PaperUnits','points','PaperPosition',[0 0 900 600],'PaperSize',[900 600])
export_fig Figure4.jpg -c [nan,nan,nan,nan] -noinvert

function plot_slp(lon,lat,y)
imagescnan(lon,lat,y');
set(gca,'clim',[-160,160])
colormap(gca,nclCM(143));
colorbar
set(gca,'YDir','normal')
hold on
basemap(lat(1),lat(end)-lat(1),lon(1),lon(end)-lon(1),1,0.7*[1 1 1])
% hline = findobj(gcf, 'type', 'line');
% set(hline,'LineWidth',1)
set(gca,'ytick',[10 25 40 55],'xtick',[-80 -50 -20])
xticklabels({'80^\circW','50^\circW','20^\circW'})
yticklabels({'10^\circN','25^\circN','40^\circN','55^\circN'})
axis([-95 -5 5 65])
set(gca,'fontsize',20)
end

function plot_buoywind(lonlp,latlp,slp,lon,lat,taux,tauy,Nwint)
imagescnan(lonlp,latlp,slp');
%set(gca,'clim',[-40,40])
set(gca,'clim',[-1.9e-9,1.9e-9])
colormap(gca,nclCM(332));
colorbar
set(gca,'YDir','normal')
hold on
quiver(lon(1:Nwint:end),lat(1:Nwint:end),...
    taux(1:Nwint:end,1:Nwint:end)',tauy(1:Nwint:end,1:Nwint:end)'...
       ,0,'color','k', 'MaxHeadSize',8,'linewidth',1)
hold on
basemap(latlp(1),latlp(end)-latlp(1)+1.1,lonlp(1),lonlp(end)-lonlp(1)+1.1,1,0.7*[1 1 1])
hold on
set(gca,'ytick',[10 25 40 55],'xtick',[-80 -50 -20])
xticklabels({'80^\circW','50^\circW','20^\circW'})
yticklabels({'10^\circN','25^\circN','40^\circN','55^\circN'})
axis([-95 -5 5 65])
set(gca,'fontsize',20)
end
