clear
close all
clc

%%
wo_CHT_data = load('C:\Users\yktsai0121\Desktop\AERO604_abaqus\CaseStudy_ConstrHandling\data\test1_wo\log.mat');
w_CHT_data = load('C:\Users\yktsai0121\Desktop\AERO604_abaqus\CaseStudy_ConstrHandling\data\test1_w\log.mat');
w_TAEA_data = load('C:\Users\yktsai0121\Desktop\AERO604_abaqus\CaseStudy_ConstrHandling\data\test1_w_CTAEA\log.mat');
w_DPEA_data = load('C:\Users\yktsai0121\Desktop\AERO604_abaqus\CaseStudy_ConstrHandling\data\test1_w_cDPEA\log.mat');

figure;
scatter(wo_CHT_data.M.fval_ndom(:,1),wo_CHT_data.M.fval_ndom(:,2),'r^',...
    'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.75)
title('CDP')
xlabel('Parameter: Load [MN]')
ylabel('Objective: Mass [kg]')
ylim([500,2500])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% Make svg be more ambiguous
set(gcf, 'Renderer', 'Painters');

figure
scatter(w_TAEA_data.M.fval_ndom(:,1),w_TAEA_data.M.fval_ndom(:,2),50,'s',...
    'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor','none','MarkerFaceAlpha',0.75);
title('C-TAEA')
xlabel('Parameter: Load [MN]')
ylabel('Objective: Mass [kg]')
ylim([500,2500])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% Make svg be more ambiguous
set(gcf, 'Renderer', 'Painters');

figure
scatter(w_DPEA_data.M.fval_ndom(:,1),w_DPEA_data.M.fval_ndom(:,2),50,'d',...
    'MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerEdgeColor','none','MarkerFaceAlpha',0.75);
title('c-DPEA')
xlabel('Parameter: Load [MN]')
ylabel('Objective: Mass [kg]')
ylim([500,2500])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% Make svg be more ambiguous
set(gcf, 'Renderer', 'Painters');

figure
scatter(w_CHT_data.M.fval_ndom(:,1),w_CHT_data.M.fval_ndom(:,2),'bo',...
    'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.75);
title('Bayesian CHT')
xlabel('Parameter: Load [MN]')
ylabel('Objective: Mass [kg]')
ylim([500,2500])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% Make svg be more ambiguous
set(gcf, 'Renderer', 'Painters');
% legend('CDP','Bayesian CHT','Location','northwest')

figure
plot(wo_CHT_data.M.hvi(1,1:5:end),wo_CHT_data.M.hvi(2,1:5:end),'r-^',...
    'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none'); hold on
plot(w_TAEA_data.M.hvi(1,1:5:end),w_TAEA_data.M.hvi(2,1:5:end),'-s',...
     'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor','none')
plot(w_DPEA_data.M.hvi(1,1:5:end),w_DPEA_data.M.hvi(2,1:5:end),'-d',...
    'MarkerFaceColor',[0.4940 0.1840 0.5560],'color',[0.4940 0.1840 0.5560])
plot(w_CHT_data.M.hvi(1,1:5:end),w_CHT_data.M.hvi(2,1:5:end),'b-o',...
    'MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none')
xlabel('Generation')
ylabel('pHVI')
legend('CDP','C-TAEA','c-DPEA','Bayesian CHT','Location','southeast')
