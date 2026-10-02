% output of ..NAO2... is from ...all_NAOnoATL3AMM_AMOnoATL3AMM_ATL3_AMM.mat
% output of ..NAO1... is from ...all_NAO_AMOnoATL3AMM_ATL3_AMM.mat
clear
load G:/PCbackup/Northeast_US/data/bathymetry/NEUS_etopo1_twosatgrid.mat
lonb=lon; latb=lat;
load G:/CUdesktop/Northeast_US/sl_extreme/climatemode_relation/cr0mxcrpos_decsltgMABSABlp8yr_noIBwindGMSLtrend_IAPsteric19502023.mat lon lat mxcr
% load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/fitEta_model/prepareIAPERA5/IAPstericslaERA5_1deg_notrendGMSL.mat sla t
% etaobs=sla;
load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/fitEta_model/IAPstericERA5/fitEtanotrendGMSLlp8yr_IAPERA5windbuoyc_allNAOAMO1_pfromfitEtav3.mat
etamod_nao=etamod_wind+etamod_buoyc;
etaobs_nao=etaobs;
load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/fitEta_model/IAPstericERA5/fitEtanotrendGMSLlp8yr_IAPERA5windbuoyc_allATL3AMM1_pfromfitEtav3.mat
etamod_tna=etamod_wind+etamod_buoyc;
etaobs_tna=etaobs;

load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/reg_openocean/openoceaneffect_iapstericsladecnoIBwindGMSL_RW_allNAOAMO1_19502023.mat
sltg1open_nao=sltg1open; sltg2open_nao=sltg2open; 
load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/reg_openocean/openoceaneffect_iapstericsladecnoIBwindGMSL_RW_allATL3AMM1_19502023.mat
sltg1open_tna=sltg1open; sltg2open_tna=sltg2open; 
sltg1open_naotna=sltg1open_nao+sltg1open_tna;
sltg2open_naotna=sltg2open_nao+sltg2open_tna;

% 
% etamod_nao=etamod_nao+etamod_amo;
% etamod_tna=etamod_tna+etamod_atl3;

%ok=squeeze(~isnan(etamod(:,:,100)));
topdb=interp2(lonb',latb,topdb',lon',lat)';
mxcr(repmat(topdb,[2,1,1])>-100)=nan;

Nlon=length(lon); Nlat=length(lat); Nt=length(t);
%calculate correlation matrix
cr0_nao=nan([Nlon Nlat]);
pval0_nao=nan([Nlon Nlat]);
cr0_tna=nan([Nlon Nlat]);
pval0_tna=nan([Nlon Nlat]);
for m=1:Nlon
    for n=1:Nlat
        if sum(~isnan(squeeze(etaobs_nao(m,n,:))))>0&&sum(~isnan(squeeze(etamod_nao(m,n,:))))>0
            [cr0_nao(m,n),pval0_nao(m,n)]=cal_corr(t,squeeze(etaobs_nao(m,n,:)),squeeze(etamod_nao(m,n,:)));
        end        
        if sum(~isnan(squeeze(etaobs_nao(m,n,:))))>0&&sum(~isnan(squeeze(etamod_tna(m,n,:))))>0
            [cr0_tna(m,n),pval0_tna(m,n)]=cal_corr(t,squeeze(etaobs_tna(m,n,:)),squeeze(etamod_tna(m,n,:)));
        end  
    end
end
% % VARIANCE RATIO MAP
% stdratio=squeeze(nanstd(etamod,1,3))./squeeze(nanstd(etaobs,1,3)); %
% stdratio(~ok)=nan;
%  
sltgrem1=sltg1o; sltgrem2=sltg2o;
co1_nao=round(corr(sltgrem1,sltg1open_nao,'rows','pairwise'),2);
co2_nao=round(corr(sltgrem2,sltg2open_nao,'rows','pairwise'),2);
so1_nao=round((1-nanvar(sltgrem1-sltg1open_nao)./nanvar(sltgrem1))*100,1);
so2_nao=round((1-nanvar(sltgrem2-sltg2open_nao)./nanvar(sltgrem2))*100,1);
co1_tna=round(corr(sltgrem1,sltg1open_tna,'rows','pairwise'),2);
co2_tna=round(corr(sltgrem2,sltg2open_tna,'rows','pairwise'),2);
so1_tna=round((1-nanvar(sltgrem1-sltg1open_tna)./nanvar(sltgrem1))*100,1);
so2_tna=round((1-nanvar(sltgrem2-sltg2open_tna)./nanvar(sltgrem2))*100,1);
co1_naotna=round(corr(sltgrem1,sltg1open_naotna,'rows','pairwise'),2);
co2_naotna=round(corr(sltgrem2,sltg2open_naotna,'rows','pairwise'),2);
so1_naotna=round((1-nanvar(sltgrem1-sltg1open_naotna)./nanvar(sltgrem1))*100,1);
so2_naotna=round((1-nanvar(sltgrem2-sltg2open_naotna)./nanvar(sltgrem2))*100,1);

% co1_nao_b1990=round(corr(sltgrem1(t<1990),sltg1open_nao(t<1990),'rows','pairwise'),2);
% co1_nao_a1990=round(corr(sltgrem1(t>1990),sltg1open_nao(t>1990),'rows','pairwise'),2);
% co2_nao_b1990=round(corr(sltgrem2(t<1990),sltg1open_nao(t<1990),'rows','pairwise'),2);
% co2_nao_a1990=round(corr(sltgrem2(t>1990),sltg1open_nao(t>1990),'rows','pairwise'),2);
% co1_tna_b1990=round(corr(sltgrem1(t<1990),sltg1open_tna(t<1990),'rows','pairwise'),2);
% co1_tna_a1990=round(corr(sltgrem1(t>1990),sltg1open_tna(t>1990),'rows','pairwise'),2);
% co2_tna_b1990=round(corr(sltgrem2(t<1990),sltg1open_tna(t<1990),'rows','pairwise'),2);
% co2_tna_a1990=round(corr(sltgrem2(t>1990),sltg1open_tna(t>1990),'rows','pairwise'),2);

% set up the subplot position
xl=0.08; xr=0.96; yb=0.65; yt=0.96;
dxem=0.04; dyem=0.04;
Nc=2; Nr=1;
axposit1=ax_position(xl,xr,yb,yt,dxem,dyem,Nc,Nr);

xl=0.08; xr=0.96; yb=0.05; yt=0.55;
dxem=0.04; dyem=0.04;
Nc=1; Nr=2;
axposit2=ax_position(xl,xr,yb,yt,dxem,dyem,Nc,Nr);

% calcualte trend of eddy number
figure('units','points','position',[0 0 800 700])

axes('position',axposit1(1,1,:));
plot_corr(lon,lat,cr0_nao,pval0_nao,mxcr)
hold on
contour(lon,lat,squeeze(mxcr(1,:,:))',[0.6,0.6],'g','linewidth',2) % SAB sea level correlation
text(0.02,1.07,'(a) r(h_{Obs},h_{NAO+AMO})','FontSize',20,'Units','normalized')

axes('position',axposit1(2,1,:));
plot_corr(lon,lat,cr0_tna,pval0_tna,mxcr)
hold on
contour(lon,lat,squeeze(mxcr(2,:,:))',[0.6,0.6],'g','linewidth',2) % MAB sea level correlation
yticklabels([])
text(0.02,1.07,'(b) r(h_{Obs},h_{AMM+ATL3})','FontSize',20,'Units','normalized')

% -------------time series--------------------------- 
axes('position',axposit2(1,2,:));
plot(t,sltgrem2,'k',t,sltg2open_nao,'b',t,sltg2open_tna,'r','LineWidth',3)
% hold on
% plot(t,sltg2open_naotna,'Color',0.5*[1 1 1],'LineWidth',4)
xlim([1950,2024])
ylim([-90,90])
xticklabels([])
ylabel('SLAs (mm)','fontsize',20)
set(gca,'fontsize',20)
legend('\eta^R','\eta^{EST}_{NAO+AMO}','\eta^{EST}_{AMM+ATL3}',...
        'Location','NorthEast','Orientation','horizontal','NumColumns',4,'fontsize',18)
text(0.01,0.87,['r(\eta^R,\eta^{EST}_{NAO+AMO})=',num2str(co2_nao), ';  r(\eta^R,\eta^{EST}_{AMM+ATL3})=',num2str(co2_tna)],...
    'FontSize',18,'Units','normalized')
text(0.01,0.12,['s(\eta^{EST}_{NAO+AMO})=',num2str(so2_nao),'%', ';  s(\eta^{EST}_{AMM+ATL3})=',num2str(so2_tna),'%'],...
    'FontSize',18,'Units','normalized')
text(0.01,1.09,'(c) North USEC','FontSize',20,'Units','normalized')

axes('position',axposit2(1,1,:));
plot(t,sltgrem1,'k',t,sltg1open_nao,'b',t,sltg1open_tna,'r','linewidth',3)
% hold on
% plot(t,sltg1open_naotna,'Color',0.5*[1 1 1],'LineWidth',4)
xlim([1950,2024])
ylim([-90,90])
ylabel('SLAs (mm)','fontsize',20)
set(gca,'fontsize',20)
legend('\eta^R','\eta^{EST}_{NAO+AMO}','\eta^{EST}_{AMM+ATL3}',...
        'Location','NorthEast','Orientation','horizontal','NumColumns',4,'fontsize',18)
text(0.01,0.87,['r(\eta^R,\eta^{EST}_{NAO+AMO})=',num2str(co1_nao), ';  r(\eta^R,\eta^{EST}_{AMM+ATL3})=',num2str(co1_tna),...
 ],'FontSize',18,'Units','normalized')
text(0.01,0.12,['s(\eta^{EST}_{NAO+AMO})=',num2str(so1_nao),'%', ';  s(\eta^{EST}_{AMM+ATL3})=',num2str(so1_tna),'%',...
    ],'FontSize',18,'Units','normalized')
text(0.01,1.09,'(d) South USEC', 'FontSize',20,'Units','normalized')

set(gcf,'PaperUnits','points','PaperPosition',[0 0 800 700],'Papersize',[800 700])
export_fig Figure3.jpg -c [nan,nan,nan,nan] -noinvert

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
%contour(lon,lat,squeeze(mxcr(2,:,:))',[0.6,0.6],'c','linewidth',3)
set(gca,'ytick',[10 25 40 55],'xtick',[-80 -50 -20])
xticklabels({'80^\circW','50^\circW','20^\circW'})
yticklabels({'10^\circN','25^\circN','40^\circN','55^\circN'})
axis([-95 -5 5 65])
set(gca,'fontsize',20)
end

function [cr0,pval0]=cal_corr(t,y1,y2)

    %y1=rmTrend(y1);
    %y2=rmTrend(y2);
    %ar1=max(armodel(y(~isnan(y1)),1),0);
    %ar2=max(armodel(y(~isnan(y2)),1),0);
    %dof_ar1=min([Nt*(1-ar1)/(1+ar1),Nt]);
    dof1=dof_autocorr(t,y1);
    dof2=dof_autocorr(t,y2);
    dof=max(min(dof1,dof2),6);
    %--------------------------
    % lag=0 correlation
    %--------------------------
    [cr0,~]=corr(y1,y2,'rows','pairwise');
    tr=abs(cr0)*sqrt((dof-2)/(1-cr0^2));
    pval0= tcdf(tr,dof,'upper');% pval is from t test

end
