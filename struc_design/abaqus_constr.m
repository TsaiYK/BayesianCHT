function [c,ceq] = abaqus_constr(x,h0,L,yielding_stress,Disp_allowed)
% hL = x(1);
% t = x(2);
% l1 = x(3);
% l2 = x(4);
% r1 = x(5);
% r2 = x(6);
% tol = 0.0001; % unit: m

% c(1) = (h0-hL)/2/L*l1-h0/2+r1+tol*10;
% c(2) = (h0-hL)/2/L*l2-h0/2+r2+tol*5;
% c(3) = r1+r2-(l2-l1)+tol;
% c(4) = r1-l1-tol;
% c(5) = r2+l2-L-tol;

filename_output = 'abaqus_results.txt';
% if c<=0 % check the geometry constraints, if it satisfies, run abaqus
    % Delete if the previous one exists
    if exist(filename_output, 'file')==2
        delete(filename_output);
    end
    
    % Delete if the previous one exists
    if exist('DesignVariables.txt', 'file')==2
        delete('DesignVariables.txt');
    end
    fileDV = fopen('DesignVariables.txt','w');
    
    % Call out python code and run FEA
    fprintf(fileDV,'%.4f\n',x);
    fclose(fileDV);
    system('cd C:\Users\yktsai0121\Desktop\AERO604_abaqus\CaseStudy_ConstrHandling');
    disp('Starting FEA analysis in Abaqus...')
    tic
    system('abaqus cae noGUI=abaqus_py_script.py');
    disp('Done!')
    toc
    
    % Read the data
    fileID = fopen(filename_output,'r');
    output_values = fscanf(fileID,'%f %f');
    fclose(fileID);
    max_stress = output_values(1);
    tipDisp = abs(output_values(2));
    mass = output_values(3);
    c(1) = (max_stress-yielding_stress)/1e6;
    c(2) = (tipDisp-Disp_allowed)*1e3;
    fprintf('max stress: %.4f MPa, tipDisp: %.4f mm, mass: %.4f kg\n',...
        max_stress*1e-6,tipDisp*1e3,mass)
% else
%     disp('Geometrically infeasible!')
%     c(6) = NaN;
% end
    ceq = [];
end