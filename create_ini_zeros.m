%
%=========================================================================%
% create_ini_zeros.m
%
% Create uniformly zero initial conditions file for Biology Model.
%
% Initial conditions are taken to be zero everywhere.  Good for testing
% purposes.
%
% Author: Z. Wallace
% Created: 17 September 2018
% Last edit: 
%=========================================================================%
%


function create_ini_zeros(fname)
    
    % overwrite netcdf file, if it exists
    if exist(fname,'file') ~= 0
        delete(fname)
    end    
    ncid = netcdf.create(fname,'NETCDF4');  % create and open new ini file
    netcdf.close(ncid);
    
    %
    %----------------------------------------------------------------------
    % Create dimensions
    %----------------------------------------------------------------------
    %
    % grid dimensions
    xi_rho  = 482;
    xi_u    = 481;
    eta_rho = 770;
    eta_v   = 769;
    s_rho   = 40;
    s_w     = 41;
    time    = 1;
        
    %
    %----------------------------------------------------------------------
    % Create variables and add attributes
    %----------------------------------------------------------------------    
    %
    % Ocean time
    nccreate(fname,'ocean_time', ...
            'Dimensions',{'ocean_time' time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ocean_time','long_name','seconds since 1900-01-01');
    ncwriteatt(fname,'ocean_time','units','seconds');
 
    % Free-surface
    nccreate(fname,'zeta', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'zeta','long_name','free-surface');
    ncwriteatt(fname,'zeta','units','meter');

    % 2D u-direction velocity
    nccreate(fname,'ubar', ...
            'Dimensions',{'xi_u' xi_u 'eta_rho' eta_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ubar','long_name','Vertically integrated u-momentum component');
    ncwriteatt(fname,'ubar','units','meter second-1');

    % 2D v-direction velocity
    nccreate(fname,'vbar', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_v' eta_v 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'vbar','long_name','Vertically integrated v-momentum component');
    ncwriteatt(fname,'vbar','units','meter second-1');

    % 3D u-direction velocity
    nccreate(fname,'u', ...
            'Dimensions',{'xi_u' xi_u 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'u','long_name','u-momentum component');
    ncwriteatt(fname,'u','units','meter second-1');

    % 3D v-direction velocity
    nccreate(fname,'v', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_v' eta_v 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'v','long_name','v-momentum component');
    ncwriteatt(fname,'v','units','meter second-1');

    % temperature
    nccreate(fname,'temp', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'temp','long_name','potential temperature');
    ncwriteatt(fname,'temp','units','meter second-1');

    % salinity
    nccreate(fname,'salt', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'salt','long_name','salinity');
    ncwriteatt(fname,'salt','units','PSU');

    % NO3
    nccreate(fname,'NO3', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NO3','long_name','nitrate concentration');
    ncwriteatt(fname,'NO3','units','mmol N m-3');

    % NH4
    nccreate(fname,'NH4', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NH4','long_name','NH4 concentration');
    ncwriteatt(fname,'NH4','units','mmol N m-3');

    % Phytoplankton
    nccreate(fname,'Phytoplankton', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Phytoplankton','long_name','Phytoplankton concentration');
    ncwriteatt(fname,'Phytoplankton','units','mmol N m-3');

    % Diatoms
    nccreate(fname,'Diatoms', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Diatoms','long_name','Diatom concentration');
    ncwriteatt(fname,'Diatoms','units','mmol N m-3');

    % Bacteria
    nccreate(fname,'Bacteria', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Bacteria','long_name','Bacteria concentration');
    ncwriteatt(fname,'Bacteria','units','mmol N m-3');

    % Chlorophyll-a
    nccreate(fname,'Chlorophyll', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Chlorophyll','long_name','Chlorophyll concentration');
    ncwriteatt(fname,'Chlorophyll','units','mg Chl m-3');

    % Diatom Chlorophyll
    nccreate(fname,'ChlD', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ChlD','long_name','ChlD concentration');
    ncwriteatt(fname,'ChlD','units','mg Chl m-3');

    % DON
    nccreate(fname,'DON', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DON','long_name','DON concentration');
    ncwriteatt(fname,'DON','units','mmol N m-3');

    % DOC
    nccreate(fname,'DOC', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DOC','long_name','DOC concentration');
    ncwriteatt(fname,'DOC','units','mmol C m-3');

    % DetN
    nccreate(fname,'DetN', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetN','long_name','DetN concentration');
    ncwriteatt(fname,'DetN','units','mmol N m-3');

    % DetC
    nccreate(fname,'DetC', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetC','long_name','DetC concentration');
    ncwriteatt(fname,'DetC','units','mmol C m-3');

    % Silica
    nccreate(fname,'Silica', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Silica','long_name','Silica concentration');
    ncwriteatt(fname,'Silica','units','mmol Si m-3');

    % LDetN
    nccreate(fname,'LDetN', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetN','long_name','LDetN concentration');
    ncwriteatt(fname,'LDetN','units','mmol N m-3');

    % LDetC
    nccreate(fname,'LDetC', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetC','long_name','LDetC concentration');
    ncwriteatt(fname,'LDetC','units','mmol C m-3');

    % LDetS
    nccreate(fname,'LDetS', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetS','long_name','LDetS concentration');
    ncwriteatt(fname,'LDetS','units','mmol Si m-3');

    % Microzooplankton
    nccreate(fname,'Microzooplankton', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Microzooplankton','long_name','Microzooplankton concentration');
    ncwriteatt(fname,'Microzooplankton','units','mmol N m-3');

    % Mesozooplankton
    nccreate(fname,'Mesozooplankton', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Mesozooplankton','long_name','Mesozooplankton concentration');
    ncwriteatt(fname,'Mesozooplankton','units','mmol N m-3');

    % Dissolved Iron
    nccreate(fname,'FeD', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'FeD','long_name','Dissolved iron concentration');
    ncwriteatt(fname,'FeD','units','umol Fe m-3');

    % Particulate Iron
    nccreate(fname,'FeP', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'FeP','long_name','Particulate iron concentration');
    ncwriteatt(fname,'FeP','units','umol Fe m-3');


    %
    %----------------------------------------------------------------------
    % Read data into physical variables.  Where data are unavailable, 
    % populate with analytical functions.  Only set biological variables to
    % zero (emulating passive tracers), since we don't want to disrupt the
    % circulation.
    %----------------------------------------------------------------------    
    %
    myroot   = '/home/server/pi/homes/zwallace/ROMS/';
    localdir = 'Project_Patagonia/vinz_in_3-21-18';
    fid = fullfile(myroot,localdir,'rst_PHYSIQUE_NEMURO_CHLAseawif_2001Jun_WOA_VinzFeC_ini.nc'); % get name of Vincent's initial conditions file
    CONST_VAL = zeros(xi_rho,eta_rho,s_rho);

    % Time variables
%    scrum_time = ncread(fid,'scrum_time');
%    ocean_time = scrum_time-367*86400;  % subtract a year to account for indexing difference
    ocean_time = 37064.*86400.;
    
    % zeta
    zeta  = ncread(fid,'zeta');

    % ubar
    ubar  = ncread(fid,'ubar');

    % vbar
    vbar  = ncread(fid,'vbar');

    % u
    u  = ncread(fid,'u');

    % v
    v  = ncread(fid,'v');

    % temp
    temp  = ncread(fid,'temp');

    % salt
    salt  = ncread(fid,'salt');

    % NO3
    NO3  = CONST_VAL;

    % NH4
    NH4  = CONST_VAL;

    % Phytoplankton
    Phytoplankton = CONST_VAL;

    % Diatoms
    Diatoms = CONST_VAL;

    % Bacteria
    Bacteria = CONST_VAL;

    % Chlorophyll
    Chlorophyll = CONST_VAL;

    % ChlD
    ChlD = CONST_VAL;

    % DON
    DON = CONST_VAL;

    % DOC
    DOC  = CONST_VAL;

    % DetN
    DetN = CONST_VAL;

    % DetC
    DetC = CONST_VAL;

    % Silica
    Silica  = CONST_VAL;

    % LDetN
    LDetN = CONST_VAL;

    % LDetC
    LDetC = CONST_VAL;

    % LDetS
    LDetS = CONST_VAL;

    % Microzooplankton
    Microzooplankton = CONST_VAL;

    % Mesozooplankton
    Mesozooplankton = CONST_VAL;

    % Dissolved Iron
    FeD = CONST_VAL;

    % Particulate Iron
    FeP = CONST_VAL;


    %
    %----------------------------------------------------------------------
    % Write data to variables.
    %----------------------------------------------------------------------    
    %
    % bry_time
    ncwrite(fname,'ocean_time',ocean_time);

    % zeta
    ncwrite(fname,'zeta',zeta);

    % ubar
    ncwrite(fname,'ubar',ubar);

    % vbar
    ncwrite(fname,'vbar',vbar);

    % u
    ncwrite(fname,'u',u);

    % v
    ncwrite(fname,'v',v);

    % temp
    ncwrite(fname,'temp',temp);

    % salt
    ncwrite(fname,'salt',salt);

    % NO3
    ncwrite(fname,'NO3',NO3);

    % NH4
    ncwrite(fname,'NH4',NH4);

    % Phytoplankton
    ncwrite(fname,'Phytoplankton',Phytoplankton);

    % Diatoms
    ncwrite(fname,'Diatoms',Diatoms);

    % Bacteria
    ncwrite(fname,'Bacteria',Bacteria);

    % Chlorophyll
    ncwrite(fname,'Chlorophyll',Chlorophyll);

    % ChlD
    ncwrite(fname,'ChlD',ChlD);

    % DON
    ncwrite(fname,'DON',DON);

    % DOC
    ncwrite(fname,'DOC',DOC);

    % DetN
    ncwrite(fname,'DetN',DetN);

    % DetC
    ncwrite(fname,'DetC',DetC);

    % Silica
    ncwrite(fname,'Silica',Silica);

    % LDetN
    ncwrite(fname,'LDetN',LDetN);

    % LDetC
    ncwrite(fname,'LDetC',LDetC);

    % LDetS
    ncwrite(fname,'LDetS',LDetS);

    % Microzooplankton
    ncwrite(fname,'Microzooplankton',Microzooplankton);

    % Mesozooplankton
    ncwrite(fname,'Mesozooplankton',Mesozooplankton);

    % Dissolved Iron
    ncwrite(fname,'FeD',FeD);
    
    % Particulate Iron
    ncwrite(fname,'FeP',FeP);


end  % create_ini_zeros
