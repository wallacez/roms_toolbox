% create_rivers.m
%
% Create river forcing NetCDF file.
%
% Author: Z. Wallace
% Created 12 September 2018
% Last edit: 9.13.18 --> hacky workaround to ensure river_time is
%                        unlimited. Fix in the near future.

function create_rivers (ncname) 

     % overwrite netcdf file, if it exists
    if exist(ncname,'file') ~= 0
        delete(ncname)
    end

    % define parameters to be used as dimensions
    xi_rho  = 482;
    eta_rho = 770;
    xi_u    = 481;
    eta_v   = 769;
    xi_psi  = 481;
    eta_psi = 769;
    s_rho   = 40;
    river   = 12;
    river_time    = 0;

    ncid = netcdf.create(ncname,'NETCDF4');
    
    % create dimensions
    d1 = netcdf.defDim(ncid,'xi_rho' ,xi_rho );
    d2 = netcdf.defDim(ncid,'eta_rho',eta_rho);
    d3 = netcdf.defDim(ncid,'xi_u'   ,xi_u   );
    d4 = netcdf.defDim(ncid,'eta_v'  ,e  ta_v  );
    d5 = netcdf.defDim(ncid,'xi_psi' ,xi_psi );
    d6 = netcdf.defDim(ncid,'eta_psi',eta_psi);
    d7 = netcdf.defDim(ncid,'s_rho'  ,s_rho  );
    d8 = netcdf.defDim(ncid,'river'  ,river  );
    d9 = netcdf.defDim(ncid,'river_time',river_time);
    
    % define variables
    varid0  = netcdf.defVar(ncid,'river' ,'double',d8);
    varid1  = netcdf.defVar(ncid,'river_Xposition','double',d8);
    varid2  = netcdf.defVar(ncid,'river_Eposition','double',d8);
    varid3  = netcdf.defVar(ncid,'river_direction','double',d8);    
    varid4  = netcdf.defVar(ncid,'river_Vshape','double',[d8,d7]);
    varid5  = netcdf.defVar(ncid,'river_transport','double',[d8,d9]);
    varid6  = netcdf.defVar(ncid,'river_time','double',d9);
    varid7  = netcdf.defVar(ncid,'river_temp','double',[d8,d7,d9]);
    varid8  = netcdf.defVar(ncid,'river_salt','double',[d8,d7,d9]);
    varid9  = netcdf.defVar(ncid,'river_NO3','double',[d8,d7,d9]);
    varid10  = netcdf.defVar(ncid,'river_NH4','double',[d8,d7,d9]);
    varid11  = netcdf.defVar(ncid,'river_Silica','double',[d8,d7,d9]);
    varid12  = netcdf.defVar(ncid,'river_FeD','double',[d8,d7,d9]);
    
    % define variable attributes
    % river
    netcdf.putAtt(ncid,varid0,'long_name','river runoff identification number')
    
    % river_xpos
    netcdf.putAtt(ncid,varid1,'long_name','river XI-position at rho-points')
    
    % river_epos
    netcdf.putAtt(ncid,varid2,'long_name','river ETA-position at rho-points')
    
    % river_direction
    netcdf.putAtt(ncid,varid3,'long_name','river runoff direction');
    netcdf.putAtt(ncid,varid3,'flag_values', '0, 1');
    netcdf.putAtt(ncid,varid3,'flag_meanings','flow across u-face, flow across v-face');
    netcdf.putAtt(ncid,varid3,'LwSrc_True','flag not used');
    
    % river_Vshape
    netcdf.putAtt(ncid,varid4,'long_name','river runoff mass transport vertical profile');
    
    % river_transport
    netcdf.putAtt(ncid,varid5,'long_name','river runoff mass transport');
    netcdf.putAtt(ncid,varid5,'units','meter-cubed second-1');
    
    % river_time
    netcdf.putAtt(ncid,varid6,'long_name','river runoff time')
    netcdf.putAtt(ncid,varid6,'units','days since 1900-01-01 00:00:00')
    netcdf.putAtt(ncid,varid6,'cycle_length',365.25);
    
    % river_temp
    netcdf.putAtt(ncid,varid7,'long_name','river runoff potential temperature');
    netcdf.putAtt(ncid,varid7,'units','Celsius');
    netcdf.putAtt(ncid,varid7,'time','river_time');

    % river_salt
    netcdf.putAtt(ncid,varid8,'long_name','river runoff salinity');
    netcdf.putAtt(ncid,varid8,'time','river_time');
    
    % river_NO3
    netcdf.putAtt(ncid,varid9,'long_name','river runoff nitrate concentration');
    netcdf.putAtt(ncid,varid9,'units','mmol N m-3');
    netcdf.putAtt(ncid,varid9,'time','river_time');
    
    % river_NH4
    netcdf.putAtt(ncid,varid10,'long_name','river runoff ammonia concentration');
    netcdf.putAtt(ncid,varid10,'units','mmol N m-3');
    netcdf.putAtt(ncid,varid10,'time','river_time');
    
    % river_SiOH
    netcdf.putAtt(ncid,varid11,'long_name','river runoff Silica concentration');
    netcdf.putAtt(ncid,varid11,'units','mmol Si m-3');
    netcdf.putAtt(ncid,varid11,'time','river_time');
    
    % river_FeD
    netcdf.putAtt(ncid,varid12,'long_name','river runoff dissolved iron concentration');
    netcdf.putAtt(ncid,varid12,'units','umol Fe m-3');
    netcdf.putAtt(ncid,varid12,'time','river_time');
    
    % define global attributes
    varid = netcdf.getConstant('GLOBAL');
    
    netcdf.putAtt(ncid,varid,'type'    ,'ROMS river forcing file')
    netcdf.putAtt(ncid,varid,'title'   ,'Patagonia model rivers')
    netcdf.putAtt(ncid,varid,'grd_file','roms_grd_rivers.nc')
    netcdf.putAtt(ncid,varid,'date'    , datestr(today))
    
    netcdf.close(ncid);
    
    % set river runoff values
    riv_time = [15.2188, 45.6563, 76.0938, 106.5313, 136.9688, 167.4063,...
                197.8438, 228.2813, 258.7188, 289.1563, 319.5938, 350.0313];
    ncwrite('roms_rivers.nc','river_time',riv_time);
    
    ncid = netcdf.open('roms_rivers.nc','WRITE');
    riv_num  = 1:12;
    riv_xpos = [277,278,279,280, ...
                153,154,...
                353,...
                231,...
                225,...
                197,...
                187,...
                159];  % Grouped by rivers
    riv_epos = [596,596,596,596,...
                284,284,...
                618,...
                501,...
                483,...
                444,...
                367,...
                326];  % Grouped by rivers
    riv_dir  = [1,1,1,1,1,1,1,1,1,0,0,1];
    transpo = [-6000,-6000,-6000,-6000,...
                -5000,-5000,...
                -2000,...
                -148,...
                -1014,...
                51,...
                5,...
                -790]; % Grouped by rivers, in m3/s
    for i=1:river
        for t=1:12
            riv_trns(i,t) = transpo(i);
        end
    end
        % river_Vshape -- uniform distribution over the vertical 
    for Nsrc=1:12
        for k=1:40
            riv_vshp(Nsrc,k) = 1/s_rho;
        end
    end
    
    % Set 3D tracer variables.
    temps = [20,20,20,20,15,15,12,12,12,12,10,10];
    salt = [5,5,5,5,10,10,5,5,5,5,5,5];
    nit = [35,35,35,35,10,10,20,20,20,20,20,20];
    ammn = [15,15,15,15,2,2,2,2,2,2,2,2];
    sioh = [35,35,35,35,4,4,20,20,20,20,20,20];
    fed = [2,2,2,2,2,2,2,2,2,2,2,2];
    
    for i=1:river
        for k=1:s_rho
            for t=1:length(riv_time)
                riv_temp(i,k,t) = temps(i);
                riv_salt(i,k,t) = salt(i);
                riv_nit(i,k,t) = nit(i);
                riv_ammn(i,k,t) = ammn(i);
                riv_sioh(i,k,t) = sioh(i);
                riv_fed(i,k,t) = fed(i);
            end
        end
    end
    
    % Write data to variables
    netcdf.putVar(ncid,varid0,riv_num);
    netcdf.putVar(ncid,varid1,riv_xpos);
    netcdf.putVar(ncid,varid2,riv_epos);
    netcdf.putVar(ncid,varid3,riv_dir);
    netcdf.putVar(ncid,varid4,riv_vshp);
    netcdf.putVar(ncid,varid5,riv_trns);
    %netcdf.putVar(ncid,varid6,riv_time);
    netcdf.putVar(ncid,varid7,riv_temp);
    netcdf.putVar(ncid,varid8,riv_salt);
    netcdf.putVar(ncid,varid9,riv_nit);
    netcdf.putVar(ncid,varid10,riv_ammn);
    netcdf.putVar(ncid,varid11,riv_sioh);
    netcdf.putVar(ncid,varid12,riv_fed);
    
    netcdf.close(ncid);

end

