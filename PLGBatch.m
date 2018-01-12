projectName = 'latconv';
unitCell = {'bcz'};
resolution = 12;
strut_dia = 300e-3;
ball_dia = 1.5*strut_dia ;
reps_x = 3;
reps_y = 3;
reps_z = 3;
unitSizeX = 2;
unitSizeY = 2;
unitSizeZ = 2;

%% main loop
for unitInc = 1:length(unitCell);
    currentType = unitCell{unitInc};
    for resInc = resolution;
        for strInc = strut_dia;
            for ballInc = ball_dia;
                for repXInc = reps_x;
                    for repYInc = reps_y;
                        for repZInc = reps_z;
                            for usxInc = unitSizeX;
                                for usyInc = unitSizeY;
                                    for uszInc = unitSizeZ;
                                        name = sprintf('%s_%s_r=%3.1f_d=%3.1f_db=%3.1f_rx=%3.1f_ry=%3.1f_rz=%3.1f_usx=%3.1f_usy=%3.1f_usz=%3.1f',...
                                            projectName,currentType,resInc,strInc,ballInc,...
                                            repXInc,repYInc,repZInc,usxInc,usyInc,uszInc);
                                        obj = PLG(currentType,resInc,strInc,...
                                            1,ballInc,resInc,...
                                            usxInc,usyInc,uszInc,...
                                            repXInc,repYInc,repZInc,...
                                            0,0,0);
                                        %saveCustom(obj,[name,'.custom'],[cd,filesep]);
                                        %saveStl(obj,[name,'.stl'],[cd,filesep]);
                                        getProperties(obj,[name,'.xlsx'])
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end