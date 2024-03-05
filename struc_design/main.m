nS = 1;
xDesign = [0.25000, 2.500000, 0.25, 0.75, 0.1, 0.05];
h0 = 0.5; L = 1;
filename_output = 'abaqus_results.txt';
density = 2720; % unit: kg/m^3

%Loops for each set of DVs and runs Abaqus code, calculates min eigenvalue
%and weight for each set and populates arrays
for i = 1:nS
    disp(i);
    % Delete if the previous one exists
    if exist(filename_output, 'file')==2
        delete(filename_output);
    end
    % check the geometry constraints
    c(i,:) = geo_constr(xDesign(i,:),h0,L);
    if c<=0
        % Delete if the previous one exists
        if exist('DesignVariables.txt', 'file')==2
            delete('DesignVariables.txt');
        end
        fileDV = fopen('DesignVariables.txt','w');

        % Call out python code and run FEA
        fprintf(fileDV,'%.4f\n',xDesign(i,:));
        fclose(fileDV);
        system('cd C:\Users\yktsai0121\Desktop\AERO604_abaqus\CaseStudy_ConstrHandling');
        disp('Starting FEA analysis in Abaqus...')
        tic
        system('abaqus cae noGUI=abaqus_py_script.py');
        disp('Done!')
        toc
        
        % Read the data
        fileID = fopen(filename_output,'r');
        output_values = fscanf(fileID,'%f');
        fclose(fileID);
        max_stress(i) = output_values(1);
        mass(i) = mass_fun(xDesign(i,:),h0,L,density);
        
    else
        disp('Geometrically infeasible!')
    end
    
end


