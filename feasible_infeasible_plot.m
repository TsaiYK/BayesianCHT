function [X_feasible_new,X_infeasible_new] = feasible_infeasible_plot(X_new,feasible,PlotOrNot)
if nargin<3
    PlotOrNot = 'DontPlot';
end
% Classify based on the feasibility
X_feasible_new = X_new(feasible,:);
X_infeasible_new = X_new;
X_infeasible_new(feasible,:) = [];

% Plot the feasible and infeasible points
if strcmp(PlotOrNot,'Plot')
    figure;
    plot(X_feasible_new(:,1),X_feasible_new(:,2),'b.')
    hold on
    plot(X_infeasible_new(:,1),X_infeasible_new(:,2),'r*')
    xlabel('x'); ylabel('\theta')
    legend('feasible','infeasible')
end