clear
close all
clc

%%
wo_CHT_data = load('C:\Users\yktsai0121\Desktop\AERO604_abaqus\CaseStudy_AeroWing\data\p3ga_results_wo_CHT\log.mat');
w_TAEA_data = load('C:\Users\yktsai0121\Desktop\AERO604_abaqus\CaseStudy_AeroWing\data\p3ga_results_w_CTAEA\log.mat');
w_DPEA_data = load('C:\Users\yktsai0121\Desktop\AERO604_abaqus\CaseStudy_AeroWing\data\p3ga_results_w_cDPEA\log.mat');
w_CHT_data = load('C:\Users\yktsai0121\Desktop\AERO604_abaqus\CaseStudy_AeroWing\data\p3ga_results_w_CHT\log.mat');

% figure;
% scatter3(wo_CHT_data.M.fval_ndom(:,3),wo_CHT_data.M.fval_ndom(:,1),...
%     wo_CHT_data.M.fval_ndom(:,2),'r^','filled'); hold on
% scatter3(w_CHT_data.M.fval_ndom(:,3),w_CHT_data.M.fval_ndom(:,1),...
%     w_CHT_data.M.fval_ndom(:,2),'bo','filled');
% xlabel('Parameter: air speed [m/s]')
% ylabel('Objective 1: Mass [kg]')
% zlabel('Objective 2: -Safety factor')
% legend('CDP','Bayesian CHT')

% % X_wo = [wo_CHT_data.M.fval_ndom(:,3)];
% % Y_wo = [wo_CHT_data.M.fval_ndom(:,1),wo_CHT_data.M.fval_ndom(:,2)];
% % [model_wo, perf_wo] = Kriging_general(X_wo,Y_wo);
% % 
% % X_w = [w_CHT_data.M.fval_ndom(:,3)];
% % Y_w = [w_CHT_data.M.fval_ndom(:,1),w_CHT_data.M.fval_ndom(:,2)];
% % [model_w, perf_w] = Kriging_general(X_w,Y_w);
% 
% % [X,Y] = meshgrid(linspace(40,65,100),linspace(40,120,100));
% % X = X'; Y = Y';
% % for i = 1:100
% %     for j = 1:100
% %         Z_wo(i,j) = predictor([X(i,j),Y(i,j)], model_wo);
% %         Z_w(i,j) = predictor([X(i,j),Y(i,j)], model_w);
% %     end
% % end
% % figure
% % surf(X,Y,Z_wo)
% % figure
% % surf(X,Y,Z_w)
% 
% %%
% [~,I] = sort(wo_CHT_data.M.fval_ndom(:,3));
% Vwo(:,1) = wo_CHT_data.M.fval_ndom(I,1);
% Vwo(:,2) = wo_CHT_data.M.fval_ndom(I,2);
% Vwo(:,3) = wo_CHT_data.M.fval_ndom(I,3);
% 
% [~,I] = sort(w_CHT_data.M.fval_ndom(:,3));
% Vw(:,1) = w_CHT_data.M.fval_ndom(I,1);
% Vw(:,2) = w_CHT_data.M.fval_ndom(I,2);
% Vw(:,3) = w_CHT_data.M.fval_ndom(I,3);
% 
% lb_para = [40,];
% ub_para = [42,];
% med_para = (lb_para+ub_para)/2;
% for i = 1:5
%     I_wo{i} = find(Vwo(:,3)>lb_para(i) & Vwo(:,3)<ub_para(i));
%     I_w{i} = find(Vw(:,3)>lb_para(i) & Vw(:,3)<ub_para(i));
%     figure
%     scatter(Vwo(I_wo{i},1),Vwo(I_wo{i},2),'r^','filled'); hold on
%     scatter(Vw(I_wo{i},1),Vw(I_wo{i},2),'bo','filled');
%     title(strcat('Air Speed = ',num2str(med_para(i)),' [m/s]'))
% end

% %%
% figure
% scatter(wo_CHT_data.M.hvi(1,:),wo_CHT_data.M.hvi(2,:),'r^','filled'); hold on
% scatter(w_CHT_data.M.hvi(1,:),w_CHT_data.M.hvi(2,:),'bo','filled')
% xlabel('Generation')
% ylabel('pHVI')
% legend('CDP','Bayesian CHT','Location','northwest')
% ylim([0.02,0.18])

figure;
scatter3(wo_CHT_data.M.fval_ndom(:,3),wo_CHT_data.M.fval_ndom(:,1),-wo_CHT_data.M.fval_ndom(:,2),'r^',...
    'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.75)
title('CDP')
xlabel('Parameter: air speed [m/s]')
ylabel('Objective 1: Mass [kg]')
zlabel('Objective 2: Safety factor')
xlim([40,70])
ylim([40,90])
zlim([0,20])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(gcf, 'Renderer', 'Painters');

figure
scatter3(w_TAEA_data.M.fval_ndom(:,3),w_TAEA_data.M.fval_ndom(:,1),-w_TAEA_data.M.fval_ndom(:,2),50,'s',...
    'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor','none','MarkerFaceAlpha',0.75);
title('C-TAEA')
xlabel('Parameter: air speed [m/s]')
ylabel('Objective 1: Mass [kg]')
zlabel('Objective 2: Safety factor')
xlim([40,70])
ylim([40,90])
zlim([0,20])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(gcf, 'Renderer', 'Painters');

figure
scatter3(w_DPEA_data.M.fval_ndom(:,3),w_DPEA_data.M.fval_ndom(:,1),-w_DPEA_data.M.fval_ndom(:,2),50,'d',...
    'MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerEdgeColor','none','MarkerFaceAlpha',0.75);
title('c-DPEA')
xlabel('Parameter: air speed [m/s]')
ylabel('Objective 1: Mass [kg]')
zlabel('Objective 2: Safety factor')
xlim([40,70])
ylim([40,90])
zlim([0,20])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(gcf, 'Renderer', 'Painters');

figure
scatter3(w_CHT_data.M.fval_ndom(:,3),w_CHT_data.M.fval_ndom(:,1),-w_CHT_data.M.fval_ndom(:,2),'bo',...
    'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.75);
title('Bayesian CHT')
xlabel('Parameter: air speed [m/s]')
ylabel('Objective 1: Mass [kg]')
zlabel('Objective 2: Safety factor')
xlim([40,70])
ylim([40,90])
zlim([0,20])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(gcf, 'Renderer', 'Painters');
% legend('CDP','Bayesian CHT','Location','northwest')

figure
plot(wo_CHT_data.M.hvi(1,:),wo_CHT_data.M.hvi(2,:),'r-^',...
    'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none'); hold on
plot(w_TAEA_data.M.hvi(1,:),w_TAEA_data.M.hvi(2,:),'-s',...
     'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor','none')
plot(w_DPEA_data.M.hvi(1,:),w_DPEA_data.M.hvi(2,:),'-d',...
    'MarkerFaceColor',[0.4940 0.1840 0.5560],'color',[0.4940 0.1840 0.5560])
plot(w_CHT_data.M.hvi(1,:),w_CHT_data.M.hvi(2,:),'b-o',...
    'MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none')
xlabel('Generation')
ylabel('pHVI')
legend('CDP','C-TAEA','c-DPEA','Bayesian CHT','Location','southeast')



