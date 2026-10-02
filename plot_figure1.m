clear
load G:/PCbackup/Northeast_US/data/PSMSL/rlr_monthly_2024_paperTNA/NEUS_sltg_19502023_rmseatrendlp8yr.mat
load G:/PCbackup/Northeast_US/data/GMSL_v2/Palmer21CU_19502023_rmtrendlp8yr.mat
gmslmon=gmslmon*1000;
load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/IB_v2/NEUS_ERA5paperTNAsltgib_rmseatrend19502023lp8yr.mat ibtg
load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/coastwinds_v2/sltgyinglimodel_localwind19502023rmseatrendlp8yr_wholeoptpara.mat tausinttg
%load G:/CUdesktop/Northeast_US/sealeveldec_dynamics_19502023/reg_openocean/openoceaneffect_iapstericsladecnoIBwindGMSL_19502023.mat

sltg=sltg-repmat(gmslmon,[1,size(sltg,2)]);
ilat1=lattg<35.5; ilat2=lattg>35.5;
sltg1=nanmean(sltg(:,ilat1),2); sltg2=nanmean(sltg(:,ilat2),2);

sltgloc1=nanmean(ibtg(:,ilat1),2)+nanmean(tausinttg(:,ilat1),2); 
sltgloc2=nanmean(ibtg(:,ilat2),2)+nanmean(tausinttg(:,ilat2),2);
%sltgrem1=sltg1open; sltgrem2=sltg2open ;
sltgrem1=sltg1-sltgloc1; sltgrem2=sltg2-sltgloc2 ;

cloc1=round(corr(sltgloc1,sltg1,'rows','pairwise'),2);
cloc2=round(corr(sltgloc2,sltg2,'rows','pairwise'),2);

crem1=round(corr(sltgrem1,sltg1,'rows','pairwise'),2);
crem2=round(corr(sltgrem2,sltg2,'rows','pairwise'),2);

% sloc1=round(nanstd(sltgloc1),1);
% srem1=round(nanstd(sltgrem1),1);
% ssum1=round(nanstd(sltgsum1),1);
% 
% sloc2=round(nanstd(sltgloc2),1);
% srem2=round(nanstd(sltgrem2),1);

sloc1=round((1-nanvar(sltg1-sltgloc1)./nanvar(sltg1))*100,1);
sloc2=round((1-nanvar(sltg2-sltgloc2)./nanvar(sltg2))*100,1);
srem1=round((1-nanvar(sltg1-sltgrem1)./nanvar(sltg1))*100,1);
srem2=round((1-nanvar(sltg2-sltgrem2)./nanvar(sltg2))*100,1);

% set up the subplot position
xl=0.07; xr=0.4; yb=0.1; yt=0.93;
dxem=0.01; dyem=0.04;
Nc=1; Nr=1;
axposit1=ax_position(xl,xr,yb,yt,dxem,dyem,Nc,Nr);

xl=0.48; xr=0.98; yb=0.1; yt=0.93;
dxem=0.01; dyem=0.08;
Nc=1; Nr=2;
axposit=ax_position(xl,xr,yb,yt,dxem,dyem,Nc,Nr);

% calcualte trend of eddy number
figure('units','points','position',[0 0 900 400])

%----------- domain--------------------------------
axes('position',axposit1(1,1,:));
basemap(25,25,-83,32,1,0.7*[1 1 1]);
hold on
plot(lontg,lattg,'o','color','g',...
    'markerfacecolor','g','markersize',10);
plot(-76:-76+30,35.2*ones(31,1),'--k','linewidth',2);
text(-80,32,'South USEC','Color','k','FontWeight','bold','FontSize',22)
text(-73,39.5,'North USEC','Color','k','FontWeight','bold','FontSize',22)
set(gca,'ytick',[30 40 50],'xtick',[-80 -70])
xticklabels({'80^\circW','70^\circW'})
yticklabels({'30^\circN','40^\circN','50^\circN'})
axis([-86 -60 25 51])
set(gca,'fontsize',20)
text(0.01,1.04,'(a) Tide gauge locations','FontSize',20,'Units','normalized')

% -------------time series--------------------------- 
axes('position',axposit(1,2,:));
plot(t,sltg2,'k',t,sltgloc2,'b',t,sltgrem2,'r','linewidth',2)
xlim([1950,2024])
ylim([-80,80])
xticklabels([])
ylabel('SLAs (mm)')
set(gca,'fontsize',20)
legend('\eta^{O}','\eta^L','\eta^R',...
        'Location','NorthEast','Orientation','horizontal','NumColumns',4)
text(0.01,0.91,['r(\eta^O,\eta^L)=',num2str(cloc2),'\rm;  r(\eta^O,\eta^{R})=',num2str(crem2)],...
  'FontSize',20,'Units','normalized')
text(0.01,0.12,['s(\eta^L)=',num2str(sloc2),'%;  s(\eta^R)=',num2str(srem2),'%'],'FontSize',20,'Units','normalized')
text(0.01,1.09,'(b) North USEC','FontSize',20,'Units','normalized')

axes('position',axposit(1,1,:));
plot(t,sltg1,'k',t,sltgloc1,'b',t,sltgrem1,'r','linewidth',2)
xlim([1950,2024])
ylim([-80,80])
ylabel('SLAs (mm)')
set(gca,'fontsize',20)
legend('\eta^{O}','\eta^L','\eta^R',...
        'Location','NorthEast','Orientation','horizontal','NumColumns',4)
text(0.01,0.91,['r(\eta^O,\eta^L)=',num2str(cloc1),'\rm;  r(\eta^O,\eta^{R})=',num2str(crem1)]...
   ,'FontSize',20,'Units','normalized')
text(0.01,0.12,['s(\eta^L)=',num2str(sloc1),'%;  s(\eta^R)=',num2str(srem1), '%'],'FontSize',20,'Units','normalized')
text(0.01,1.09,'(c) South USEC','FontSize',20,'Units','normalized')

set(gcf,'PaperUnits','points','PaperPosition',[0 0 900 400],'Papersize',[900 400])
export_fig Figure1.jpg -c [nan,nan,nan,nan] -noinvert
