clear
close all
clc
loadpath = {'Y:\GTPase4\Combined Results\RhoA\AnalysisResults_withSmoothedProtrusion_AllDepth_v2_201216y\';'Y:\GTPase4\Combined Results\Rac\AnalysisResults_withSmoothedProtrusion_AllDepth_v2_201216y\'};
savefolder = 'Y:\GTPase4\Combined Results\Summary20210509_smoothed\';
mkdir(savefolder)


load([savefolder(1:end-14) '0329_smoothed\' 'RhoARacDistance6Valuerange5FORCE2FORCEr34p56.mat'], '-regexp', '^(?!savefolder$)\w');


%% correct the normalized velocity
temp1 = nanmean(Rac.SPEEDprotrusionLineNorm(:,13));
Rac.SPEEDprotrusionLineNorm = Rac.SPEEDprotrusionLineNorm - temp1;
Rac.SPEEDprotrusionLineNorm_CI = Rac.SPEEDprotrusionLineNorm_CI - temp1;

temp2 = nanmean(Rac.SPEEDretractionLineNorm(:,13));
Rac.SPEEDretractionLineNorm = Rac.SPEEDretractionLineNorm - temp2;
Rac.SPEEDretractionLineNorm_CI = Rac.SPEEDretractionLineNorm_CI - temp2;

%% define plot colors
color.RhoA = [0 0.4470 0.7410];
color.Rac = [0.4940 0.1840 0.5560];
% color.force1 = [0.8500 0.3250 0.0980];
color.force1 = '#FDA50F';
color.force5 = [0.6350 0.0780 0.1840];
color.speed = '#0B6623';
color.xcorr = 'r';

%% selected figure plot LineNorm
selected.data = {nanmean(Rac.FRETprotrusionLineNorm{1});nanmean(RhoA.FRETprotrusionLineNorm{1});nanmean(Rac.FORCEprotrusionLineNorm{1});nanmean(Rac.FORCEprotrusionLineNorm{5});
    nanmean(Rac.FRETretractionLineNorm{1});nanmean(RhoA.FRETretractionLineNorm{1});nanmean(Rac.FORCEretractionLineNorm{1});nanmean(Rac.FORCEretractionLineNorm{5});
    nanmean(Rac.FRETprotrusionMaxVLineNorm{1});nanmean(RhoA.FRETprotrusionMaxVLineNorm{1});nanmean(Rac.FORCEprotrusionMaxVLineNorm{1});nanmean(Rac.FORCEprotrusionMaxVLineNorm{5});
    nanmean(Rac.FRETretractionMaxVLineNorm{1});nanmean(RhoA.FRETretractionMaxVLineNorm{1});nanmean(Rac.FORCEretractionMaxVLineNorm{1});nanmean(Rac.FORCEretractionMaxVLineNorm{5});
    nanmean(Rac.SPEEDprotrusionLineNorm);nanmean(Rac.SPEEDretractionLineNorm);nanmean(Rac.SPEEDprotrusionMaxVLineNorm);nanmean(Rac.SPEEDretractionMaxVLineNorm)};
selected.ci = {Rac.FRETprotrusionLineNorm_CI{1};RhoA.FRETprotrusionLineNorm_CI{1};Rac.FORCEprotrusionLineNorm_CI{1};Rac.FORCEprotrusionLineNorm_CI{5};
    Rac.FRETretractionLineNorm_CI{1};RhoA.FRETretractionLineNorm_CI{1};Rac.FORCEretractionLineNorm_CI{1};Rac.FORCEretractionLineNorm_CI{5};
    Rac.FRETprotrusionMaxVLineNorm_CI{1};RhoA.FRETprotrusionMaxVLineNorm_CI{1};Rac.FORCEprotrusionMaxVLineNorm_CI{1};Rac.FORCEprotrusionMaxVLineNorm_CI{5};
    Rac.FRETretractionMaxVLineNorm_CI{1};RhoA.FRETretractionMaxVLineNorm_CI{1};Rac.FORCEretractionMaxVLineNorm_CI{1};Rac.FORCEretractionMaxVLineNorm_CI{5};
    Rac.SPEEDprotrusionLineNorm_CI;Rac.SPEEDretractionLineNorm_CI;Rac.SPEEDprotrusionMaxVLineNorm_CI;Rac.SPEEDretractionMaxVLineNorm_CI};
selected.title = {'Rac1-GTP (normalized) at protrusion onset';'RhoA-GTP (normalized) at protrusion onset';'TF_e_d_g_e (normalized) at protrusion onset';'TF_m_a_x_ _l_a_y_e_r (normalized) at protrusion onset';
    'Rac1-GTP (normalized) at retraction onset';'RhoA-GTP (normalized) at retraction onset';'TF_e_d_g_e (normalized) at retraction onset';'TF_m_a_x_ _l_a_y_e_r (normalized) at retraction onset';
    'Rac1-GTP (normalized) at protrusion Vmax onset';'RhoA-GTP (normalized) at protrusion Vmax onset';'TF_e_d_g_e (normalized) at protrusion Vmax onset';'TF_m_a_x_ _l_a_y_e_r (normalized) at protrusion Vmax onset';
    'Rac1-GTP (normalized) at retraction Vmax onset';'RhoA-GTP (normalized) at retraction Vmax onset';'TF_e_d_g_e (normalized) at retraction Vmax onset';'TF_m_a_x_ _l_a_y_e_r (normalized) at retraction Vmax onset';
    'membrane velocity (normalized) at protrusion';'membrane velocity (normalized) at retraction';'membrane velocity (normalized) at protrusion V_m_a_x';'membrane velocity (normalized) at retraction V_m_a_x'};
selected.axis = {[-120 120 -0.15 0.15];[-120 120 -0.15 0.15];[-120 120 -0.3 0.3];[-120 120 -0.15 0.15];
    [-120 120 -0.15 0.15];[-120 120 -0.2 0.2];[-120 120 -0.35 0.35];[-120 120 -0.25 0.25];
    [-120 120 -0.2 0.15];[-120 120 -0.15 0.15];[-120 120 -0.5 0.5];[-120 120 -0.15 0.15];
    [-120 120 -0.2 0.25];[-120 120 -0.25 0.25];[-120 120 -0.25 0.25];[-120 120 -0.3 0.3];
    [-120 120 -1.5 1.5];[-120 120 -1.6 1.6];[-120 120 -2.5 2.5];[-120 120 -2 2]};
selected.ylabel = {'Rac1-GTP';'RhoA-GTP';'TF_e_d_g_e';'TF_m_a_x_ _l_a_y_e_r';
    'Rac1-GTP';'RhoA-GTP';'TF_e_d_g_e';'TF_m_a_x_ _l_a_y_e_r';
    'Rac1-GTP';'RhoA-GTP';'TF_e_d_g_e';'TF_m_a_x_ _l_a_y_e_r';
    'Rac1-GTP';'RhoA-GTP';'TF_e_d_g_e';'TF_m_a_x_ _l_a_y_e_r';
    'Velocity';'Velocity';'Velocity';'Velocity'};
selected.color = {color.Rac;color.RhoA;color.force1;color.force5;
    color.Rac;color.RhoA;color.force1;color.force5;
    color.Rac;color.RhoA;color.force1;color.force5;
    color.Rac;color.RhoA;color.force1;color.force5;
    color.speed;color.speed;color.speed;color.speed};

for i = 1 : size(selected.data,1)
    figure(1)
    set(gcf,'color','w')
    plot_ci((-12:12)*10,[selected.data{i};selected.ci{i}(1,:);selected.ci{i}(2,:)],'PatchColor', selected.color{i}, 'PatchAlpha', 0.3, 'MainLineWidth', 8, 'MainLineStyle', '-', 'MainLineColor', selected.color{i},'LineWidth', 0.5, 'LineStyle','none', 'LineColor', selected.color{i});
    axis(selected.axis{i});
    set(gca,'ycolor',selected.color{i})
%     ylabel(selected.ylabel{i}, 'FontSize',24);
    hold on;
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    set(gca,'FontWeight','bold','FontSize',54);
    xlabel('t(s)', 'FontSize',54);
%     title(selected.title{i}, 'FontSize',24);
    figure(1);  set(gcf,'Position',[0 0 2000 1500]); saveas(gcf,[savefolder num2str(i) ' ' selected.title{i} '.fig']);
    figure(1);  set(gcf,'Position',[0 0 2000 1500]); export_fig(gcf,[savefolder num2str(i) ' ' selected.title{i} '.png']);
    close(1)
end





%% selected figure plot cross-correlation to SPEED(normalized)_selected timerange_mean & individual correlation

data.X = {Rac.FRETprotrusionLineNorm{1}(:,7:19);RhoA.FRETprotrusionLineNorm{1}(:,7:19);Rac.FORCEprotrusionLineNorm{1}(:,7:19);Rac.FORCEprotrusionLineNorm{5}(:,7:19);
    Rac.FRETretractionLineNorm{1}(:,7:19);RhoA.FRETretractionLineNorm{1}(:,7:19);Rac.FORCEretractionLineNorm{1}(:,7:19);Rac.FORCEretractionLineNorm{5}(:,7:19);
    Rac.FRETprotrusionMaxVLineNorm{1}(:,7:19);RhoA.FRETprotrusionMaxVLineNorm{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm{5}(:,7:19);
    Rac.FRETretractionMaxVLineNorm{1}(:,7:19);RhoA.FRETretractionMaxVLineNorm{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm{5}(:,7:19)};
data.Y = {Rac.SPEEDprotrusionLineNorm(:,7:19);RhoA.SPEEDprotrusionLineNorm(:,7:19);Rac.SPEEDprotrusionLineNorm(:,7:19);Rac.SPEEDprotrusionLineNorm(:,7:19);
    Rac.SPEEDretractionLineNorm(:,7:19);RhoA.SPEEDretractionLineNorm(:,7:19);Rac.SPEEDretractionLineNorm(:,7:19);Rac.SPEEDretractionLineNorm(:,7:19);
    Rac.SPEEDprotrusionMaxVLineNorm(:,7:19);RhoA.SPEEDprotrusionMaxVLineNorm(:,7:19);Rac.SPEEDprotrusionMaxVLineNorm(:,7:19);Rac.SPEEDprotrusionMaxVLineNorm(:,7:19);
    Rac.SPEEDretractionMaxVLineNorm(:,7:19);RhoA.SPEEDretractionMaxVLineNorm(:,7:19);Rac.SPEEDretractionMaxVLineNorm(:,7:19);Rac.SPEEDretractionMaxVLineNorm(:,7:19)};
data.X_CI = {Rac.FRETprotrusionLineNorm_CI{1}(:,7:19);RhoA.FRETprotrusionLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionLineNorm_CI{5}(:,7:19);
    Rac.FRETretractionLineNorm_CI{1}(:,7:19);RhoA.FRETretractionLineNorm_CI{1}(:,7:19);Rac.FORCEretractionLineNorm_CI{1}(:,7:19);Rac.FORCEretractionLineNorm_CI{5}(:,7:19);
    Rac.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);RhoA.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm_CI{5}(:,7:19);
    Rac.FRETretractionMaxVLineNorm_CI{1}(:,7:19);RhoA.FRETretractionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm_CI{5}(:,7:19)};
data.Y_CI = {Rac.SPEEDprotrusionLineNorm_CI(:,7:19);RhoA.SPEEDprotrusionLineNorm_CI(:,7:19);Rac.SPEEDprotrusionLineNorm_CI(:,7:19);Rac.SPEEDprotrusionLineNorm_CI(:,7:19);
    Rac.SPEEDretractionLineNorm_CI(:,7:19);RhoA.SPEEDretractionLineNorm_CI(:,7:19);Rac.SPEEDretractionLineNorm_CI(:,7:19);Rac.SPEEDretractionLineNorm_CI(:,7:19);
    Rac.SPEEDprotrusionMaxVLineNorm_CI(:,7:19);RhoA.SPEEDprotrusionMaxVLineNorm_CI(:,7:19);Rac.SPEEDprotrusionMaxVLineNorm_CI(:,7:19);Rac.SPEEDprotrusionMaxVLineNorm_CI(:,7:19);
    Rac.SPEEDretractionMaxVLineNorm_CI(:,7:19);RhoA.SPEEDretractionMaxVLineNorm_CI(:,7:19);Rac.SPEEDretractionMaxVLineNorm_CI(:,7:19);Rac.SPEEDretractionMaxVLineNorm_CI(:,7:19)};
data.Xaxis = {[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -0.33 0.33];[-60 60 -0.15 0.15];
    [-60 60 -0.15 0.15];[-60 60 -0.23 0.23];[-60 60 -0.35 0.35];[-60 60 -0.25 0.25];
    [-60 60 -0.23 0.23];[-60 60 -0.15 0.15];[-60 60 -0.55 0.55];[-60 60 -0.15 0.15];
    [-60 60 -0.25 0.25];[-60 60 -0.25 0.25];[-60 60 -0.25 0.25];[-60 60 -0.33 0.33];};
data.Yaxis = {[-60 60 -1.5 1.5];[-60 60 -1.5 1.5];[-60 60 -1.5 1.5];[-60 60 -1.5 1.5];
    [-60 60 -1.5 1.5];[-60 60 -1.5 1.5];[-60 60 -1.5 1.5];[-60 60 -1.5 1.5];
    [-60 60 -2.5 2.5];[-60 60 -2.5 2.5];[-60 60 -2.5 2.5];[-60 60 -2.5 2.5];
    [-60 60 -2 2];[-60 60 -2 2];[-60 60 -2 2];[-60 60 -2 2];};

data.color = {color.Rac;color.RhoA;color.force1;color.force5;
    color.Rac;color.RhoA;color.force1;color.force5;
    color.Rac;color.RhoA;color.force1;color.force5;
    color.Rac;color.RhoA;color.force1;color.force5;};
data.colorY = color.speed;
data.name = {'Rac1-GTP vs velocity at protrusion';'RhoA-GTP vs velocity at protrusion';'TFedge vs velocity at protrusion';'TFmax layer vs velocity at protrusion';
    'Rac1-GTP vs velocity at retraction';'RhoA-GTP vs velocity at retraction';'TFedge vs velocity at retraction';'TFmax layer vs velocity at retraction';
    'Rac1-GTP vs velocity at protrusionVmax';'RhoA-GTP vs velocity at protrusionVmax';'TFedge vs velocity at protrusionVmax';'TFmax layer vs velocity at protrusionVmax';
    'Rac1-GTP vs velocity at retractionVmax';'RhoA-GTP vs velocity at retractionVmax';'TFedge vs velocity at retractionVmax';'TFmax layer vs velocity at retractionVmax';};
data.ylabel = {'Rac1-GTP';'RhoA-GTP';'TFedge';'TFmax layer';
    'Rac1-GTP';'RhoA-GTP';'TFedge';'TFmax layer';
    'Rac1-GTP';'RhoA-GTP';'TFedge';'TFmax layer';
    'Rac1-GTP';'RhoA-GTP';'TFedge';'TFmax layer';};

SubPlotOrder = [1:8:32,3:8:32,2:8:32,4:8:32];

clearvars r_mean r_individual r_individual_CI lag_mean lag_individual tempr templag

for i =1 : size(data.X,1)
    % cross-correlation of nanmean value
    [r_mean{i},lag_mean{i}] = xcorr(normalize(nanmean(data.X{i})),normalize(nanmean(data.Y{i})),'normalized');



end


% plot cross correlation of nanmean value
for i = 1 : size(data.X,1)
    figure(3)
    set(gcf,'color','w')
    sp1 = subplot(2,1,1);
    sp1.Position = [0.1 0.55 0.8 0.4];
    yyaxis left
    plot_ci((-6:6)*10,[nanmean(data.X{i});data.X_CI{i}(1,:);data.X_CI{i}(2,:)],'PatchColor', data.color{i}, 'PatchAlpha', 0.3, 'MainLineWidth', 5, 'MainLineStyle', '-', 'MainLineColor', data.color{i},'LineWidth', 0.5, 'LineStyle','none', 'LineColor', data.color{i});
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    axis(data.Xaxis{i})
    set(gca,'ycolor',data.color{i})
    yyaxis right
    plot_ci((-6:6)*10,[nanmean(data.Y{i});data.Y_CI{i}(1,:);data.Y_CI{i}(2,:)],'PatchColor', data.colorY, 'PatchAlpha', 0.3, 'MainLineWidth', 5, 'MainLineStyle', '-', 'MainLineColor', data.colorY,'LineWidth', 0.5, 'LineStyle','none', 'LineColor', data.colorY);
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    axis(data.Yaxis{i})
%     xlabel('t(s)', 'FontSize',54);
    set(gca,'ycolor',data.colorY)
    set(gca,'FontWeight','bold', 'FontSize',54);
%     yyaxis left
%     ylabel([data.ylabel{i}], 'FontSize',44);
%     yyaxis right
%     ylabel('Velocity', 'FontSize',44);
    set(gca,'xticklabel',[])

    sp2 = subplot(2,1,2);
    sp2.Position = [0.1 0.25 0.8 0.25];
    plot(lag_mean{i}*10,r_mean{i},'.-', 'color',color.xcorr,'LineWidth', 7,'MarkerSize',10)
    set(gca,'ycolor',data.color{i})
%     ylabel('Correlation', 'FontSize',44);
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    set(gca,'FontWeight','bold');
    axis([-60 60 -1 1]);
    set(gca,'FontWeight','bold','FontSize',54);
    xlabel('t(s)', 'FontSize',54);
%     title(data.name{i}, 'FontSize',24);
    figure(3);  set(gcf,'Position',[0 0 2000 1500]); saveas(gcf,[savefolder num2str(i+20) ' ' data.name{i} '_mean(normalized).fig']);
    figure(3);  set(gcf,'Position',[0 0 2000 1500]); export_fig(gcf,[savefolder num2str(i+20) ' ' data.name{i} '_mean(normalized).png']);
    close(3)
end



%% selected figure plot cross-correlation to Rac1-GTP(normalized)_selected timerange_mean & individual correlation

data.X = {Rac.SPEEDprotrusionLineNorm(:,7:19);RhoA.FRETprotrusionLineNorm{1}(:,7:19);Rac.FORCEprotrusionLineNorm{1}(:,7:19);Rac.FORCEprotrusionLineNorm{5}(:,7:19);
    Rac.SPEEDretractionLineNorm(:,7:19);RhoA.FRETretractionLineNorm{1}(:,7:19);Rac.FORCEretractionLineNorm{1}(:,7:19);Rac.FORCEretractionLineNorm{5}(:,7:19);
    Rac.SPEEDprotrusionMaxVLineNorm(:,7:19);RhoA.FRETprotrusionMaxVLineNorm{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm{5}(:,7:19);
    Rac.SPEEDretractionMaxVLineNorm(:,7:19);RhoA.FRETretractionMaxVLineNorm{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm{5}(:,7:19)};
data.Y = {Rac.FRETprotrusionLineNorm{1}(:,7:19);Rac.FRETprotrusionLineNorm{1}(:,7:19);Rac.FRETprotrusionLineNorm{1}(:,7:19);Rac.FRETprotrusionLineNorm{1}(:,7:19);
    Rac.FRETretractionLineNorm{1}(:,7:19);Rac.FRETretractionLineNorm{1}(:,7:19);Rac.FRETretractionLineNorm{1}(:,7:19);Rac.FRETretractionLineNorm{1}(:,7:19);
    Rac.FRETprotrusionMaxVLineNorm{1}(:,7:19);Rac.FRETprotrusionMaxVLineNorm{1}(:,7:19);Rac.FRETprotrusionMaxVLineNorm{1}(:,7:19);Rac.FRETprotrusionMaxVLineNorm{1}(:,7:19);
    Rac.FRETretractionMaxVLineNorm{1}(:,7:19);Rac.FRETretractionMaxVLineNorm{1}(:,7:19);Rac.FRETretractionMaxVLineNorm{1}(:,7:19);Rac.FRETretractionMaxVLineNorm{1}(:,7:19)};
data.X_CI = {Rac.SPEEDprotrusionLineNorm_CI(:,7:19);RhoA.FRETprotrusionLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionLineNorm_CI{5}(:,7:19);
    Rac.SPEEDretractionLineNorm_CI(:,7:19);RhoA.FRETretractionLineNorm_CI{1}(:,7:19);Rac.FORCEretractionLineNorm_CI{1}(:,7:19);Rac.FORCEretractionLineNorm_CI{5}(:,7:19);
    Rac.SPEEDprotrusionMaxVLineNorm_CI(:,7:19);RhoA.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm_CI{5}(:,7:19);
    Rac.SPEEDretractionMaxVLineNorm_CI(:,7:19);RhoA.FRETretractionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm_CI{5}(:,7:19)};
data.Y_CI = {Rac.FRETprotrusionLineNorm_CI{1}(:,7:19);Rac.FRETprotrusionLineNorm_CI{1}(:,7:19);Rac.FRETprotrusionLineNorm_CI{1}(:,7:19);Rac.FRETprotrusionLineNorm_CI{1}(:,7:19);
    Rac.FRETretractionLineNorm_CI{1}(:,7:19);Rac.FRETretractionLineNorm_CI{1}(:,7:19);Rac.FRETretractionLineNorm_CI{1}(:,7:19);Rac.FRETretractionLineNorm_CI{1}(:,7:19);
    Rac.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);Rac.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);Rac.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);Rac.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);
    Rac.FRETretractionMaxVLineNorm_CI{1}(:,7:19);Rac.FRETretractionMaxVLineNorm_CI{1}(:,7:19);Rac.FRETretractionMaxVLineNorm_CI{1}(:,7:19);Rac.FRETretractionMaxVLineNorm_CI{1}(:,7:19)};
data.Xaxis = {[-60 60 -1.53 1.53];[-60 60 -0.15 0.15];[-60 60 -0.33 0.33];[-60 60 -0.15 0.15];
    [-60 60 -1.53 1.53];[-60 60 -0.23 0.23];[-60 60 -0.35 0.35];[-60 60 -0.25 0.25];
    [-60 60 -2.53 2.53];[-60 60 -0.15 0.15];[-60 60 -0.55 0.55];[-60 60 -0.15 0.15];
    [-60 60 -2.3 2.3];[-60 60 -0.25 0.25];[-60 60 -0.25 0.25];[-60 60 -0.33 0.33];};
data.Yaxis = {[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];
    [-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];
    [-60 60 -0.23 0.23];[-60 60 -0.23 0.23];[-60 60 -0.23 0.23];[-60 60 -0.23 0.23];
    [-60 60 -0.25 0.25];[-60 60 -0.25 0.25];[-60 60 -0.25 0.25];[-60 60 -0.25 0.25];};

data.color = {color.speed;color.RhoA;color.force1;color.force5;
    color.speed;color.RhoA;color.force1;color.force5;
    color.speed;color.RhoA;color.force1;color.force5;
    color.speed;color.RhoA;color.force1;color.force5;};
data.colorY = color.Rac;
data.name = {'Velocity vs Rac1-GTP at protrusion';'RhoA-GTP vs Rac1-GTP at protrusion';'TFedge vs Rac1-GTP at protrusion';'TFmax layer vs Rac1-GTP at protrusion';
    'Velocity vs Rac1-GTP at retraction';'RhoA-GTP vs Rac1-GTP at retraction';'TFedge vs Rac1-GTP at retraction';'TFmax layer vs Rac1-GTP at retraction';
    'Velocity vs Rac1-GTP at protrusionVmax';'RhoA-GTP vs Rac1-GTP at protrusionVmax';'TFedge vs Rac1-GTP at protrusionVmax';'TFmax layer vs Rac1-GTP at protrusionVmax';
    'Velocity vs Rac1-GTP at retractionVmax';'RhoA-GTP vs Rac1-GTP at retractionVmax';'TFedge vs Rac1-GTP at retractionVmax';'TFmax layer vs Rac1-GTP at retractionVmax';};

data.ylabel = {'Velocity';'RhoA-GTP';'TFedge';'TFmax layer';
    'Velocity';'RhoA-GTP';'TFedge';'TFmax layer';
    'Velocity';'RhoA-GTP';'TFedge';'TFmax layer';
    'Velocity';'RhoA-GTP';'TFedge';'TFmax layer';};

SubPlotOrder = [1:8:32,3:8:32,2:8:32,4:8:32];

clearvars r_mean r_individual lag_mean lag_individual tempr templag

for i =1 : size(data.X,1)
    % cross-correlation of nanmean value
    [r_mean{i},lag_mean{i}] = xcorr(normalize(nanmean(data.X{i})),normalize(nanmean(data.Y{i})),'normalized');



end


% plot cross correlation of nanmean value
for i = 1 : size(data.X,1)
    figure(3)
    set(gcf,'color','w')
    sp1 = subplot(2,1,1);
    sp1.Position = [0.1 0.55 0.8 0.4];
    yyaxis left
    plot_ci((-6:6)*10,[nanmean(data.X{i});data.X_CI{i}(1,:);data.X_CI{i}(2,:)],'PatchColor', data.color{i}, 'PatchAlpha', 0.3, 'MainLineWidth', 5, 'MainLineStyle', '-', 'MainLineColor', data.color{i},'LineWidth', 0.5, 'LineStyle','none', 'LineColor', data.color{i});
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    axis(data.Xaxis{i})
    set(gca,'ycolor',data.color{i})
    yyaxis right
    plot_ci((-6:6)*10,[nanmean(data.Y{i});data.Y_CI{i}(1,:);data.Y_CI{i}(2,:)],'PatchColor', data.colorY, 'PatchAlpha', 0.3, 'MainLineWidth', 5, 'MainLineStyle', '-', 'MainLineColor', data.colorY,'LineWidth', 0.5, 'LineStyle','none', 'LineColor', data.colorY);
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    axis(data.Yaxis{i})
%     xlabel('t(s)', 'FontSize',54);
    set(gca,'ycolor',data.colorY)
    set(gca,'FontWeight','bold', 'FontSize',50);
%     yyaxis left
%     ylabel([data.ylabel{i}], 'FontSize',44);
%     yyaxis right
%     ylabel('Velocity', 'FontSize',44);


    set(gca,'xticklabel',[])

    sp2 = subplot(2,1,2);
    sp2.Position = [0.1 0.25 0.8 0.25];

    plot(lag_mean{i}*10,r_mean{i},'.-', 'color',color.xcorr,'LineWidth', 7,'MarkerSize',10)
    set(gca,'ycolor',data.color{i})
%     ylabel('Correlation', 'FontSize',44);
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    set(gca,'FontWeight','bold');
    axis([-60 60 -1 1]);
    set(gca,'FontWeight','bold','FontSize',50);
    xlabel('t(s)', 'FontSize',50);
%     title(data.name{i}, 'FontSize',24);
    figure(3);  set(gcf,'Position',[0 0 2000 1500]); saveas(gcf,[savefolder num2str(i+40) ' ' data.name{i} '_mean(normalized).fig']);
    figure(3);  set(gcf,'Position',[0 0 2000 1500]); export_fig(gcf,[savefolder num2str(i+40) ' ' data.name{i} '_mean(normalized).png']);
    close(3)
end





%% selected figure plot cross-correlation to RhoA-GTP(normalized)_selected timerange_mean & individual correlation

data.X = {Rac.FRETprotrusionLineNorm{1}(:,7:19);RhoA.SPEEDprotrusionLineNorm(:,7:19);RhoA.FORCEprotrusionLineNorm{1}(:,7:19);RhoA.FORCEprotrusionLineNorm{5}(:,7:19);
    Rac.FRETretractionLineNorm{1}(:,7:19);RhoA.SPEEDretractionLineNorm(:,7:19);RhoA.FORCEretractionLineNorm{1}(:,7:19);RhoA.FORCEretractionLineNorm{5}(:,7:19);
    Rac.FRETprotrusionMaxVLineNorm{1}(:,7:19);RhoA.SPEEDprotrusionMaxVLineNorm(:,7:19);RhoA.FORCEprotrusionMaxVLineNorm{1}(:,7:19);RhoA.FORCEprotrusionMaxVLineNorm{5}(:,7:19);
    Rac.FRETretractionMaxVLineNorm{1}(:,7:19);RhoA.SPEEDretractionMaxVLineNorm(:,7:19);RhoA.FORCEretractionMaxVLineNorm{1}(:,7:19);RhoA.FORCEretractionMaxVLineNorm{5}(:,7:19)};
data.Y = {RhoA.FRETprotrusionLineNorm{1}(:,7:19);RhoA.FRETprotrusionLineNorm{1}(:,7:19);RhoA.FRETprotrusionLineNorm{1}(:,7:19);RhoA.FRETprotrusionLineNorm{1}(:,7:19);
    RhoA.FRETretractionLineNorm{1}(:,7:19);RhoA.FRETretractionLineNorm{1}(:,7:19);RhoA.FRETretractionLineNorm{1}(:,7:19);RhoA.FRETretractionLineNorm{1}(:,7:19);
    RhoA.FRETprotrusionMaxVLineNorm{1}(:,7:19);RhoA.FRETprotrusionMaxVLineNorm{1}(:,7:19);RhoA.FRETprotrusionMaxVLineNorm{1}(:,7:19);RhoA.FRETprotrusionMaxVLineNorm{1}(:,7:19);
    RhoA.FRETretractionMaxVLineNorm{1}(:,7:19);RhoA.FRETretractionMaxVLineNorm{1}(:,7:19);RhoA.FRETretractionMaxVLineNorm{1}(:,7:19);RhoA.FRETretractionMaxVLineNorm{1}(:,7:19)};
data.X_CI = {Rac.FRETprotrusionLineNorm_CI{1}(:,7:19);RhoA.SPEEDprotrusionLineNorm_CI(:,7:19);RhoA.FORCEprotrusionLineNorm_CI{1}(:,7:19);RhoA.FORCEprotrusionLineNorm_CI{5}(:,7:19);
    Rac.FRETretractionLineNorm_CI{1}(:,7:19);RhoA.SPEEDretractionLineNorm_CI(:,7:19);RhoA.FORCEretractionLineNorm_CI{1}(:,7:19);RhoA.FORCEretractionLineNorm_CI{5}(:,7:19);
    Rac.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);RhoA.SPEEDprotrusionMaxVLineNorm_CI(:,7:19);RhoA.FORCEprotrusionMaxVLineNorm_CI{1}(:,7:19);RhoA.FORCEprotrusionMaxVLineNorm_CI{5}(:,7:19);
    Rac.FRETretractionMaxVLineNorm_CI{1}(:,7:19);RhoA.SPEEDretractionMaxVLineNorm_CI(:,7:19);RhoA.FORCEretractionMaxVLineNorm_CI{1}(:,7:19);RhoA.FORCEretractionMaxVLineNorm_CI{5}(:,7:19)};
data.Y_CI = {RhoA.FRETprotrusionLineNorm_CI{1}(:,7:19);RhoA.FRETprotrusionLineNorm_CI{1}(:,7:19);RhoA.FRETprotrusionLineNorm_CI{1}(:,7:19);RhoA.FRETprotrusionLineNorm_CI{1}(:,7:19);
    RhoA.FRETretractionLineNorm_CI{1}(:,7:19);RhoA.FRETretractionLineNorm_CI{1}(:,7:19);RhoA.FRETretractionLineNorm_CI{1}(:,7:19);RhoA.FRETretractionLineNorm_CI{1}(:,7:19);
    RhoA.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);RhoA.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);RhoA.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);RhoA.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);
    RhoA.FRETretractionMaxVLineNorm_CI{1}(:,7:19);RhoA.FRETretractionMaxVLineNorm_CI{1}(:,7:19);RhoA.FRETretractionMaxVLineNorm_CI{1}(:,7:19);RhoA.FRETretractionMaxVLineNorm_CI{1}(:,7:19)};
data.Xaxis = {[-60 60 -0.15 0.15];[-60 60 -1.53 1.53];[-60 60 -0.33 0.33];[-60 60 -0.15 0.15];
    [-60 60 -0.15 0.15];[-60 60 -1.53 1.53];[-60 60 -0.35 0.35];[-60 60 -0.25 0.25];
    [-60 60 -0.23 0.23];[-60 60 -2.53 2.53];[-60 60 -0.55 0.55];[-60 60 -0.15 0.15];
    [-60 60 -0.25 0.25];[-60 60 -2.3 2.3];[-60 60 -0.25 0.25];[-60 60 -0.33 0.33];};
data.Yaxis = {[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];
    [-60 60 -0.2 0.2];[-60 60 -0.2 0.2];[-60 60 -0.2 0.2];[-60 60 -0.2 0.2];
    [-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];
    [-60 60 -0.25 0.25];[-60 60 -0.25 0.25];[-60 60 -0.25 0.25];[-60 60 -0.25 0.25];};
data.color = {color.Rac;color.speed;color.force1;color.force5;
    color.Rac;color.speed;color.force1;color.force5;
    color.Rac;color.speed;color.force1;color.force5;
    color.Rac;color.speed;color.force1;color.force5;};
data.colorY = color.RhoA;
data.name = {'Rac1-GTP vs RhoA-GTP at protrusion';'Velocity vs RhoA-GTP at protrusion';'TFedge vs RhoA-GTP at protrusion';'TFmax layer vs RhoA-GTP at protrusion';
    'Rac1-GTP vs RhoA-GTP at retraction';'Velocity vs RhoA-GTP at retraction';'TFedge vs RhoA-GTP at retraction';'TFmax layer vs RhoA-GTP at retraction';
    'Rac1-GTP vs RhoA-GTP at protrusionVmax';'Velocity vs RhoA-GTP at protrusionVmax';'TFedge vs RhoA-GTP at protrusionVmax';'TFmax layer vs RhoA-GTP at protrusionVmax';
    'Rac1-GTP vs RhoA-GTP at retractionVmax';'Velocity vs RhoA-GTP at retractionVmax';'TFedge vs RhoA-GTP at retractionVmax';'TFmax layer vs RhoA-GTP at retractionVmax';};
data.ylabel = {'Rac1-GTP';'Velocity';'TFedge';'TFmax layer';
    'Rac1-GTP';'Velocity';'TFedge';'TFmax layer';
    'Rac1-GTP';'Velocity';'TFedge';'TFmax layer';
    'Rac1-GTP';'Velocity';'TFedge';'TFmax layer';};

SubPlotOrder = [1:8:32,3:8:32,2:8:32,4:8:32];

clearvars r_mean r_individual lag_mean lag_individual tempr templag

for i =1 : size(data.X,1)
    % cross-correlation of nanmean value
    [r_mean{i},lag_mean{i}] = xcorr(normalize(nanmean(data.X{i})),normalize(nanmean(data.Y{i})),'normalized');


end



% plot cross correlation of nanmean value
for i = 1 : size(data.X,1)
    figure(3)
    set(gcf,'color','w')
    sp1 = subplot(2,1,1);
    sp1.Position = [0.1 0.55 0.8 0.4];
    yyaxis left
    plot_ci((-6:6)*10,[nanmean(data.X{i});data.X_CI{i}(1,:);data.X_CI{i}(2,:)],'PatchColor', data.color{i}, 'PatchAlpha', 0.3, 'MainLineWidth', 5, 'MainLineStyle', '-', 'MainLineColor', data.color{i},'LineWidth', 0.5, 'LineStyle','none', 'LineColor', data.color{i});
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    axis(data.Xaxis{i})
    set(gca,'ycolor',data.color{i})
    yyaxis right
    plot_ci((-6:6)*10,[nanmean(data.Y{i});data.Y_CI{i}(1,:);data.Y_CI{i}(2,:)],'PatchColor', data.colorY, 'PatchAlpha', 0.3, 'MainLineWidth', 5, 'MainLineStyle', '-', 'MainLineColor', data.colorY,'LineWidth', 0.5, 'LineStyle','none', 'LineColor', data.colorY);
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    axis(data.Yaxis{i})
%     xlabel('t(s)', 'FontSize',54);
    set(gca,'ycolor',data.colorY)
    set(gca,'FontWeight','bold', 'FontSize',50);
%     yyaxis left
%     ylabel([data.ylabel{i}], 'FontSize',44);
%     yyaxis right
%     ylabel('Velocity', 'FontSize',44);


    set(gca,'xticklabel',[])

    sp2 = subplot(2,1,2);
    sp2.Position = [0.1 0.25 0.8 0.25];

    plot(lag_mean{i}*10,r_mean{i},'.-', 'color',color.xcorr,'LineWidth', 7,'MarkerSize',10)
    set(gca,'ycolor',data.color{i})
%     ylabel('Correlation', 'FontSize',44);
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    set(gca,'FontWeight','bold');
    axis([-60 60 -1 1]);
    set(gca,'FontWeight','bold','FontSize',50);
    xlabel('t(s)', 'FontSize',50);
%     title(data.name{i}, 'FontSize',24);
    figure(3);  set(gcf,'Position',[0 0 2000 1500]); saveas(gcf,[savefolder num2str(i+60) ' ' data.name{i} '_mean(normalized).fig']);
    figure(3);  set(gcf,'Position',[0 0 2000 1500]); export_fig(gcf,[savefolder num2str(i+60) ' ' data.name{i} '_mean(normalized).png']);
    close(3)
end




%% selected figure plot cross-correlation to TFedge(normalized)_selected timerange_mean & individual correlation

data.X = {Rac.FRETprotrusionLineNorm{1}(:,7:19);RhoA.FRETprotrusionLineNorm{1}(:,7:19);Rac.SPEEDprotrusionLineNorm(:,7:19);Rac.FORCEprotrusionLineNorm{5}(:,7:19);
    Rac.FRETretractionLineNorm{1}(:,7:19);RhoA.FRETretractionLineNorm{1}(:,7:19);Rac.SPEEDretractionLineNorm(:,7:19);Rac.FORCEretractionLineNorm{5}(:,7:19);
    Rac.FRETprotrusionMaxVLineNorm{1}(:,7:19);RhoA.FRETprotrusionMaxVLineNorm{1}(:,7:19);Rac.SPEEDprotrusionMaxVLineNorm(:,7:19);Rac.FORCEprotrusionMaxVLineNorm{5}(:,7:19);
    Rac.FRETretractionMaxVLineNorm{1}(:,7:19);RhoA.FRETretractionMaxVLineNorm{1}(:,7:19);Rac.SPEEDretractionMaxVLineNorm(:,7:19);Rac.FORCEretractionMaxVLineNorm{5}(:,7:19)};
data.Y = {Rac.FORCEprotrusionLineNorm{1}(:,7:19);RhoA.FORCEprotrusionLineNorm{1}(:,7:19);Rac.FORCEprotrusionLineNorm{1}(:,7:19);Rac.FORCEprotrusionLineNorm{1}(:,7:19);
    Rac.FORCEretractionLineNorm{1}(:,7:19);RhoA.FORCEretractionLineNorm{1}(:,7:19);Rac.FORCEretractionLineNorm{1}(:,7:19);Rac.FORCEretractionLineNorm{1}(:,7:19);
    Rac.FORCEprotrusionMaxVLineNorm{1}(:,7:19);RhoA.FORCEprotrusionMaxVLineNorm{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm{1}(:,7:19);
    Rac.FORCEretractionMaxVLineNorm{1}(:,7:19);RhoA.FORCEretractionMaxVLineNorm{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm{1}(:,7:19)};
data.X_CI = {Rac.FRETprotrusionLineNorm_CI{1}(:,7:19);RhoA.FRETprotrusionLineNorm_CI{1}(:,7:19);Rac.SPEEDprotrusionLineNorm_CI(:,7:19);Rac.FORCEprotrusionLineNorm_CI{5}(:,7:19);
    Rac.FRETretractionLineNorm_CI{1}(:,7:19);RhoA.FRETretractionLineNorm_CI{1}(:,7:19);Rac.SPEEDretractionLineNorm_CI(:,7:19);Rac.FORCEretractionLineNorm_CI{5}(:,7:19);
    Rac.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);RhoA.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);Rac.SPEEDprotrusionMaxVLineNorm_CI(:,7:19);Rac.FORCEprotrusionMaxVLineNorm_CI{5}(:,7:19);
    Rac.FRETretractionMaxVLineNorm_CI{1}(:,7:19);RhoA.FRETretractionMaxVLineNorm_CI{1}(:,7:19);Rac.SPEEDretractionMaxVLineNorm_CI(:,7:19);Rac.FORCEretractionMaxVLineNorm_CI{5}(:,7:19)};
data.Y_CI = {Rac.FORCEprotrusionLineNorm_CI{1}(:,7:19);RhoA.FORCEprotrusionLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionLineNorm_CI{1}(:,7:19);
    Rac.FORCEretractionLineNorm_CI{1}(:,7:19);RhoA.FORCEretractionLineNorm_CI{1}(:,7:19);Rac.FORCEretractionLineNorm_CI{1}(:,7:19);Rac.FORCEretractionLineNorm_CI{1}(:,7:19);
    Rac.FORCEprotrusionMaxVLineNorm_CI{1}(:,7:19);RhoA.FORCEprotrusionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm_CI{1}(:,7:19);
    Rac.FORCEretractionMaxVLineNorm_CI{1}(:,7:19);RhoA.FORCEretractionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm_CI{1}(:,7:19)};
data.Xaxis = {[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -1.53 1.53];[-60 60 -0.15 0.15];
    [-60 60 -0.15 0.15];[-60 60 -0.23 0.23];[-60 60 -1.53 1.53];[-60 60 -0.25 0.25];
    [-60 60 -0.23 0.23];[-60 60 -0.15 0.15];[-60 60 -2.53 2.53];[-60 60 -0.15 0.15];
    [-60 60 -0.25 0.25];[-60 60 -0.25 0.25];[-60 60 -2.3 2.3];[-60 60 -0.33 0.33];};
data.Yaxis = {[-60 60 -0.3 0.3];[-60 60 -0.3 0.3];[-60 60 -0.3 0.3];[-60 60 -0.3 0.3];
    [-60 60 -0.35 0.35];[-60 60 -0.35 0.35];[-60 60 -0.35 0.35];[-60 60 -0.35 0.35];
    [-60 60 -0.5 0.5];[-60 60 -0.5 0.5];[-60 60 -0.5 0.5];[-60 60 -0.5 0.5];
    [-60 60 -0.25 0.25];[-60 60 -0.25 0.25];[-60 60 -0.25 0.25];[-60 60 -0.25 0.25];};
data.color = {color.Rac;color.RhoA;color.speed;color.force5;
    color.Rac;color.RhoA;color.speed;color.force5;
    color.Rac;color.RhoA;color.speed;color.force5;
    color.Rac;color.RhoA;color.speed;color.force5;};
data.colorY = color.force1;
data.name = {'Rac1-GTP vs TFedge at protrusion';'RhoA-GTP vs TFedge at protrusion';'Velocity vs TFedge at protrusion';'TFmax layer vs TFedge at protrusion';
    'Rac1-GTP vs TFedge at retraction';'RhoA-GTP vs TFedge at retraction';'Velocity vs TFedge at retraction';'TFmax layer vs TFedge at retraction';
    'Rac1-GTP vs TFedge at protrusionVmax';'RhoA-GTP vs TFedge at protrusionVmax';'Velocity vs TFedge at protrusionVmax';'TFmax layer vs TFedge at protrusionVmax';
    'Rac1-GTP vs TFedge at retractionVmax';'RhoA-GTP vs TFedge at retractionVmax';'Velocity vs TFedge at retractionVmax';'TFmax layer vs TFedge at retractionVmax';};
data.ylabel = {'Rac1-GTP';'RhoA-GTP';'Velocity';'TFmax layer';
    'Rac1-GTP';'RhoA-GTP';'Velocity';'TFmax layer';
    'Rac1-GTP';'RhoA-GTP';'Velocity';'TFmax layer';
    'Rac1-GTP';'RhoA-GTP';'Velocity';'TFmax layer';};

SubPlotOrder = [1:8:32,3:8:32,2:8:32,4:8:32];

clearvars r_mean r_individual r_individual_CI lag_mean lag_individual tempr templag

for i =1 : size(data.X,1)
    % cross-correlation of nanmean value
    [r_mean{i},lag_mean{i}] = xcorr(normalize(nanmean(data.X{i})),normalize(nanmean(data.Y{i})),'normalized');


end

% plot cross correlation of nanmean value
for i = 1 : size(data.X,1)
    figure(3)
    set(gcf,'color','w')
    sp1 = subplot(2,1,1);
    sp1.Position = [0.1 0.55 0.8 0.4];
    yyaxis left
    plot_ci((-6:6)*10,[nanmean(data.X{i});data.X_CI{i}(1,:);data.X_CI{i}(2,:)],'PatchColor', data.color{i}, 'PatchAlpha', 0.3, 'MainLineWidth', 5, 'MainLineStyle', '-', 'MainLineColor', data.color{i},'LineWidth', 0.5, 'LineStyle','none', 'LineColor', data.color{i});
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    axis(data.Xaxis{i})
    set(gca,'ycolor',data.color{i})
    yyaxis right
    plot_ci((-6:6)*10,[nanmean(data.Y{i});data.Y_CI{i}(1,:);data.Y_CI{i}(2,:)],'PatchColor', data.colorY, 'PatchAlpha', 0.3, 'MainLineWidth', 5, 'MainLineStyle', '-', 'MainLineColor', data.colorY,'LineWidth', 0.5, 'LineStyle','none', 'LineColor', data.colorY);
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    axis(data.Yaxis{i})
%     xlabel('t(s)', 'FontSize',54);
    set(gca,'ycolor',data.colorY)
    set(gca,'FontWeight','bold', 'FontSize',50);
%     yyaxis left
%     ylabel([data.ylabel{i}], 'FontSize',44);
%     yyaxis right
%     ylabel('Velocity', 'FontSize',44);


    set(gca,'xticklabel',[])

    sp2 = subplot(2,1,2);
    sp2.Position = [0.1 0.25 0.8 0.25];

    plot(lag_mean{i}*10,r_mean{i},'.-', 'color',color.xcorr,'LineWidth', 7,'MarkerSize',10)
    set(gca,'ycolor',data.color{i})
%     ylabel('Correlation', 'FontSize',44);
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    set(gca,'FontWeight','bold');
    axis([-60 60 -1 1]);
    set(gca,'FontWeight','bold','FontSize',50);
    xlabel('t(s)', 'FontSize',50);
%     title(data.name{i}, 'FontSize',24);
    figure(3);  set(gcf,'Position',[0 0 2000 1500]); saveas(gcf,[savefolder num2str(i+80) ' ' data.name{i} '_mean(normalized).fig']);
    figure(3);  set(gcf,'Position',[0 0 2000 1500]); export_fig(gcf,[savefolder num2str(i+80) ' ' data.name{i} '_mean(normalized).png']);
    close(3)
end




%% selected figure plot cross-correlation to TF_maxlayer(normalized)_selected timerange_mean & individual correlation

data.X = {Rac.FRETprotrusionLineNorm{1}(:,7:19);RhoA.FRETprotrusionLineNorm{1}(:,7:19);Rac.FORCEprotrusionLineNorm{1}(:,7:19);Rac.SPEEDprotrusionLineNorm(:,7:19);
    Rac.FRETretractionLineNorm{1}(:,7:19);RhoA.FRETretractionLineNorm{1}(:,7:19);Rac.FORCEretractionLineNorm{1}(:,7:19);Rac.SPEEDretractionLineNorm(:,7:19);
    Rac.FRETprotrusionMaxVLineNorm{1}(:,7:19);RhoA.FRETprotrusionMaxVLineNorm{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm{1}(:,7:19);Rac.SPEEDprotrusionMaxVLineNorm(:,7:19);
    Rac.FRETretractionMaxVLineNorm{1}(:,7:19);RhoA.FRETretractionMaxVLineNorm{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm{1}(:,7:19);Rac.SPEEDretractionMaxVLineNorm(:,7:19)};
data.Y = {Rac.FORCEprotrusionLineNorm{5}(:,7:19);RhoA.FORCEprotrusionLineNorm{5}(:,7:19);Rac.FORCEprotrusionLineNorm{5}(:,7:19);Rac.FORCEprotrusionLineNorm{5}(:,7:19);
    Rac.FORCEretractionLineNorm{5}(:,7:19);RhoA.FORCEretractionLineNorm{5}(:,7:19);Rac.FORCEretractionLineNorm{5}(:,7:19);Rac.FORCEretractionLineNorm{5}(:,7:19);
    Rac.FORCEprotrusionMaxVLineNorm{5}(:,7:19);RhoA.FORCEprotrusionMaxVLineNorm{5}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm{5}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm{5}(:,7:19);
    Rac.FORCEretractionMaxVLineNorm{5}(:,7:19);RhoA.FORCEretractionMaxVLineNorm{5}(:,7:19);Rac.FORCEretractionMaxVLineNorm{5}(:,7:19);Rac.FORCEretractionMaxVLineNorm{5}(:,7:19)};
data.X_CI = {Rac.FRETprotrusionLineNorm_CI{1}(:,7:19);RhoA.FRETprotrusionLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionLineNorm_CI{1}(:,7:19);Rac.SPEEDprotrusionLineNorm_CI(:,7:19);
    Rac.FRETretractionLineNorm_CI{1}(:,7:19);RhoA.FRETretractionLineNorm_CI{1}(:,7:19);Rac.FORCEretractionLineNorm_CI{1}(:,7:19);Rac.SPEEDretractionLineNorm_CI(:,7:19);
    Rac.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);RhoA.FRETprotrusionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm_CI{1}(:,7:19);Rac.SPEEDprotrusionMaxVLineNorm_CI(:,7:19);
    Rac.FRETretractionMaxVLineNorm_CI{1}(:,7:19);RhoA.FRETretractionMaxVLineNorm_CI{1}(:,7:19);Rac.FORCEretractionMaxVLineNorm_CI{1}(:,7:19);Rac.SPEEDretractionMaxVLineNorm_CI(:,7:19)};
data.Y_CI = {Rac.FORCEprotrusionLineNorm_CI{5}(:,7:19);RhoA.FORCEprotrusionLineNorm_CI{5}(:,7:19);Rac.FORCEprotrusionLineNorm_CI{5}(:,7:19);Rac.FORCEprotrusionLineNorm_CI{5}(:,7:19);
    Rac.FORCEretractionLineNorm_CI{5}(:,7:19);RhoA.FORCEretractionLineNorm_CI{5}(:,7:19);Rac.FORCEretractionLineNorm_CI{5}(:,7:19);Rac.FORCEretractionLineNorm_CI{5}(:,7:19);
    Rac.FORCEprotrusionMaxVLineNorm_CI{5}(:,7:19);RhoA.FORCEprotrusionMaxVLineNorm_CI{5}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm_CI{5}(:,7:19);Rac.FORCEprotrusionMaxVLineNorm_CI{5}(:,7:19);
    Rac.FORCEretractionMaxVLineNorm_CI{5}(:,7:19);RhoA.FORCEretractionMaxVLineNorm_CI{5}(:,7:19);Rac.FORCEretractionMaxVLineNorm_CI{5}(:,7:19);Rac.FORCEretractionMaxVLineNorm_CI{5}(:,7:19)};
data.Xaxis = {[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -0.25 0.45];[-60 60 -1.5 1.5];
    [-60 60 -0.15 0.15];[-60 60 -0.2 0.2];[-60 60 -0.43 0.33];[-60 60 -1.5 1.5];
    [-60 60 -0.2 0.15];[-60 60 -0.15 0.15];[-60 60 -0.55 0.63];[-60 60 -2.5 2.5];
    [-60 60 -0.2 0.25];[-60 60 -0.25 0.25];[-60 60 -0.44 0.23];[-60 60 -2 2];};
data.Yaxis = {[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];[-60 60 -0.15 0.15];
    [-60 60 -0.25 0.25];[-60 60 -0.25 0.25];[-60 60 -0.4 0.3];[-60 60 -0.25 0.25];
    [-60 60 -0.25 0.25];[-60 60 -0.25 0.25];[-60 60 -0.35 0.25];[-60 60 -0.25 0.25];
    [-60 60 -0.17 0.17];[-60 60 -0.17 0.17];[-60 60 -0.17 0.17];[-60 60 -0.17 0.17];};
data.color = {color.Rac;color.RhoA;color.force1;color.speed;
    color.Rac;color.RhoA;color.force1;color.speed;
    color.Rac;color.RhoA;color.force1;color.speed;
    color.Rac;color.RhoA;color.force1;color.speed;};
data.colorY = color.force5;
data.name = {'Rac1-GTP vs TFmaxlayer at protrusion';'RhoA-GTP vs TFmaxlayer at protrusion';'TFedge vs TFmaxlayer at protrusion';'Velocity vs TFmaxlayer at protrusion';
    'Rac1-GTP vs TFmaxlayer at retraction';'RhoA-GTP vs TFmaxlayer at retraction';'TFedge vs TFmaxlayer at retraction';'Velocity vs TFmaxlayer at retraction';
    'Rac1-GTP vs TFmaxlayer at protrusionVmax';'RhoA-GTP vs TFmaxlayer at protrusionVmax';'TFedge vs TFmaxlayer at protrusionVmax';'Velocity vs TFmaxlayer at protrusionVmax';
    'Rac1-GTP vs TFmaxlayer at retractionVmax';'RhoA-GTP vs TFmaxlayer at retractionVmax';'TFedge vs TFmaxlayer at retractionVmax';'Velocity vs TFmaxlayer at retractionVmax';};
data.ylabel = {'Rac1-GTP';'RhoA-GTP';'TFedge';'Velocity';
    'Rac1-GTP';'RhoA-GTP';'TFedge';'Velocity';
    'Rac1-GTP';'RhoA-GTP';'TFedge';'Velocity';
    'Rac1-GTP';'RhoA-GTP';'TFedge';'Velocity';};

SubPlotOrder = [1:8:32,3:8:32,2:8:32,4:8:32];

clearvars r_mean r_individual r_individual_CI lag_mean lag_individual tempr templag

for i =1 : size(data.X,1)
    % cross-correlation of nanmean value
    [r_mean{i},lag_mean{i}] = xcorr(normalize(nanmean(data.X{i})),normalize(nanmean(data.Y{i})),'normalized');



end


% plot cross correlation of nanmean value
for i = 1 : size(data.X,1)
    figure(3)
    set(gcf,'color','w')
    sp1 = subplot(2,1,1);
    sp1.Position = [0.1 0.55 0.8 0.4];
    yyaxis left
    plot_ci((-6:6)*10,[nanmean(data.X{i});data.X_CI{i}(1,:);data.X_CI{i}(2,:)],'PatchColor', data.color{i}, 'PatchAlpha', 0.3, 'MainLineWidth', 5, 'MainLineStyle', '-', 'MainLineColor', data.color{i},'LineWidth', 0.5, 'LineStyle','none', 'LineColor', data.color{i});
%     hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    axis(data.Xaxis{i})
    set(gca,'ycolor',data.color{i})
    yyaxis right
    plot_ci((-6:6)*10,[nanmean(data.Y{i});data.Y_CI{i}(1,:);data.Y_CI{i}(2,:)],'PatchColor', data.colorY, 'PatchAlpha', 0.3, 'MainLineWidth', 5, 'MainLineStyle', '-', 'MainLineColor', data.colorY,'LineWidth', 0.5, 'LineStyle','none', 'LineColor', data.colorY);
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    axis(data.Yaxis{i})
%     xlabel('t(s)', 'FontSize',54);
    set(gca,'ycolor',data.colorY)
    set(gca,'FontWeight','bold', 'FontSize',50);
%     yyaxis left
%     ylabel([data.ylabel{i}], 'FontSize',44);
%     yyaxis right
%     ylabel('Velocity', 'FontSize',44);


    set(gca,'xticklabel',[])

    sp2 = subplot(2,1,2);
    sp2.Position = [0.1 0.25 0.8 0.25];

    plot(lag_mean{i}*10,r_mean{i},'.-', 'color',color.xcorr,'LineWidth', 7,'MarkerSize',10)
    set(gca,'ycolor',data.color{i})
%     ylabel('Correlation', 'FontSize',44);
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    set(gca,'FontWeight','bold');
    axis([-60 60 -1 1]);
    set(gca,'FontWeight','bold','FontSize',50);
    xlabel('t(s)', 'FontSize',50);
%     title(data.name{i}, 'FontSize',24);
    figure(3);  set(gcf,'Position',[0 0 2000 1500]); saveas(gcf,[savefolder num2str(i+100) ' ' data.name{i} '_mean(normalized).fig']);
    figure(3);  set(gcf,'Position',[0 0 2000 1500]); export_fig(gcf,[savefolder num2str(i+100) ' ' data.name{i} '_mean(normalized).png']);
    close(3)
end



  

%% selected figure plot LineNorm
selected.data = {nanmean(Rac.SPEEDprotrusionLineNorm(:,1:25))./10;nanmean(Rac.FRETprotrusionLineNorm{1}(:,1:25));nanmean(RhoA.FRETprotrusionLineNorm{1}(:,1:25));nanmean(Rac.FORCEprotrusionLineNorm{1}(:,1:25))./2;
    nanmean(Rac.SPEEDretractionLineNorm(:,1:25))./10;nanmean(Rac.FRETretractionLineNorm{1}(:,1:25));nanmean(RhoA.FRETretractionLineNorm{1}(:,1:25));nanmean(Rac.FORCEretractionLineNorm{1}(:,1:25))./2;
    nanmean(Rac.SPEEDprotrusionMaxVLineNorm(:,1:25))./10;nanmean(Rac.FRETprotrusionMaxVLineNorm{1}(:,1:25));nanmean(RhoA.FRETprotrusionMaxVLineNorm{1}(:,1:25));nanmean(Rac.FORCEprotrusionMaxVLineNorm{1}(:,1:25))./2;
    nanmean(Rac.SPEEDretractionMaxVLineNorm(:,1:25))./10;nanmean(Rac.FRETretractionMaxVLineNorm{1}(:,1:25));nanmean(RhoA.FRETretractionMaxVLineNorm{1}(:,1:25));nanmean(Rac.FORCEretractionMaxVLineNorm{1}(:,1:25))./2;};
selected.ci = {Rac.SPEEDprotrusionLineNorm_CI(:,1:25)./10;Rac.FRETprotrusionLineNorm_CI{1}(:,1:25);RhoA.FRETprotrusionLineNorm_CI{1}(:,1:25);Rac.FORCEprotrusionLineNorm_CI{1}(:,1:25)./2;
    Rac.SPEEDretractionLineNorm_CI(:,1:25)./10;Rac.FRETretractionLineNorm_CI{1}(:,1:25);RhoA.FRETretractionLineNorm_CI{1}(:,1:25);Rac.FORCEretractionLineNorm_CI{1}(:,1:25)./2;
    Rac.SPEEDprotrusionMaxVLineNorm_CI(:,1:25)./10;Rac.FRETprotrusionMaxVLineNorm_CI{1}(:,1:25);RhoA.FRETprotrusionMaxVLineNorm_CI{1}(:,1:25);Rac.FORCEprotrusionMaxVLineNorm_CI{1}(:,1:25)./2;
    Rac.SPEEDretractionMaxVLineNorm_CI(:,1:25)./10;Rac.FRETretractionMaxVLineNorm_CI{1}(:,1:25);RhoA.FRETretractionMaxVLineNorm_CI{1}(:,1:25);Rac.FORCEretractionMaxVLineNorm_CI{1}(:,1:25)./2;};
selected.title = {'protrusion onset';'retraction onset';'protrusion Vmax';'retration Vmax'};
selected.axis = {[-120 120 -0.2 0.2];[-120 120 -0.25 0.25];[-120 120 -0.2 0.2];[-120 120 -0.3 0.3]};
selected.ylabel = {'Velocity';'Rac1-GTP';'RhoA-GTP';'Traction force';
    'Velocity';'Rac1-GTP';'RhoA-GTP';'Traction force';
    'Velocity';'Rac1-GTP';'RhoA-GTP';'Traction force';
    'Velocity';'Rac1-GTP';'RhoA-GTP';'Traction force';};
selected.color = {color.speed;color.Rac;color.RhoA;color.force1;
                color.speed;color.Rac;color.RhoA;color.force1;
                color.speed;color.Rac;color.RhoA;color.force1;
                color.speed;color.Rac;color.RhoA;color.force1};

for i = 1 : 4 %size(selected.data,1)
    figure(1)
    set(gcf,'color','w')
%     k = i+16;
%     plot_ci((-12:12)*10,[selected.data{k}./10;selected.ci{k}(1,:)./10;selected.ci{k}(2,:)./10],'PatchColor', selected.color{k}, 'PatchAlpha', 0.3, 'MainLineWidth', 8, 'MainLineStyle', '-', 'MainLineColor', selected.color{k},'LineWidth', 0.5, 'LineStyle','none', 'LineColor', selected.color{k});
    for j = 1 : 4
        k=(i-1)*4+j;
        hold on
        plot_ci((-12:12)*10,[selected.data{k};selected.ci{k}(1,:);selected.ci{k}(2,:)],'PatchColor', selected.color{k}, 'PatchAlpha', 0.3, 'MainLineWidth', 8, 'MainLineStyle', '-', 'MainLineColor', selected.color{k},'LineWidth', 0.5, 'LineStyle','none', 'LineColor', selected.color{k});
        hold off
    end
%     
    axis(selected.axis{i});
    hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
    set(gca,'FontWeight','bold','FontSize',54);
    xlabel('t(s)', 'FontSize',54);
%     title(selected.title{i}, 'FontSize',24);
    figure(1);  set(gcf,'Position',[0 0 2000 1500]); saveas(gcf,[savefolder num2str(i+120) ' All data ' selected.title{i} '.fig']);
    figure(1);  set(gcf,'Position',[0 0 2000 1500]); export_fig(gcf,[savefolder num2str(i+120) ' All data ' selected.title{i} '.png']);
    close(1)
end

% cross-correlation of nanmean value
for k = 1 : 4
    r_all{k} = [];
    lag_all{k} = [];
    for i =1 : 4
        figure(2)
        set(gcf,'color','w')
        yyaxis left
        
        temp = 1 : 4;
        for j = temp(temp~=i)
            [r_temp,lag_temp] = xcorr(normalize(selected.data{(k-1)*4+j}(:,7:19)),normalize(selected.data{(k-1)*4+i}(:,7:19)),'normalized');
            hold on
            plot(lag_temp,r_temp,'-','color', selected.color{(k-1)*4+j},'LineWidth',7);
            hold off            
        end
        xticks(-10:5:10)
        xticklabels({'-100','-50','0','50','100'})
        hold on;plot(zeros(1,20001),-10000:10000,'--k','LineWidth',3);plot(-1000:1000,zeros(1,2001),'--k','LineWidth',3);hold off;
        axis([-12 12 -1 1])
        set(gca,'ycolor','k')
        ylabel('Cross-corr')
        xlabel('t(s)')
        yyaxis right
        yticks([])
        set(gca,'ycolor',selected.color{(k-1)*4+i})
        ylabel(selected.ylabel{(k-1)*4+i},'color',selected.color{(k-1)*4+i});

        set(gca,'FontSize',54,'FontWeight','bold')
        figure(2);  set(gcf,'Position',[0 0 2000 1500]); saveas(gcf,[savefolder num2str(k+120) ' ' selected.title{k} ' cross-corr to ' selected.ylabel{(k-1)*4+i} '.fig']);
        figure(2);  set(gcf,'Position',[0 0 2000 1500]); export_fig(gcf,[savefolder num2str(k+120) ' ' selected.title{k} ' cross-corr to ' selected.ylabel{(k-1)*4+i} '.png']);
        close(2)
        
    end
end
    



