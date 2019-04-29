%
%=========================================================================%
% create_bry.m
%
% Create boundary file for Biology Model.  For each physical and biological
% variable, impose conditions at N,S,W and E boundaries, as well as the
% time associated with that variable.  
%
% Boundary conditions for the physics are monthly averages from 2000-2007 
% taken from OFES.  
% Boundary conditions for the biology are monthly climatologies taken from
% WOA(09?)
%
% Author: Z. Wallace
% Created: 13 December 2017
% Last edit: 5.19.28 --> Added cycle_length attribute to time vars.
%                 5.20.18 --> Updated DetFe bry time and conditions.
%                 8.31.18 --> Cleaned up fid paths
%                 11.5.18 --> Added 3rd Fe pool, changed FeP to DetFe
%                 4.22.19 --> Added coccolithophore variables
%=========================================================================%
%


function create_bry_vinz(fname)
    
    ncid = netcdf.create(fname,'NETCDF4');  % create and open new bry file
    netcdf.close(ncid)
    
    % get name of Vincent's boundary file
    myroot   = '/home/server/pi/homes/zwallace/ROMS/';
    localdir = 'Project_Patagonia/orig/vinz_in_3-21-18';
    fid = fullfile(myroot,localdir,'roms_bry_NEMURO_CHLAseawif_2000_2007_WOA_VinzFeC_clm.nc.1');
    fid_cocco=fullfile(myroot,'Project_Patagonia/update/roms_ini_coccos.nc');

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
    
    % time dimensions (12 mon/yr * 8 yr)
%    bry_time   = 96;
   bry_time   = 12;
   NO3_time   = 12;
   NH4_time   = 12;
   phy_time   = 12;
   dia_time   = 12;
   bac_time   = 12;
   chla_time  = 12;
   chld_time  = 12;
   DON_time   = 12;
   DOC_time   = 12;
   DetN_time  = 12;
   DetC_time  = 12;
   sil_time   = 12;
   LDetN_time = 12;
   LDetC_time = 12;
   LDetS_time = 12;
   micro_time = 12;
   meso_time  = 12;
   FeD_time   = 12;
   DetFe_time   = 12;
   LDetFe_time = 12;
   coccos_time = 12;
   ChlC_time = 12;
   PIC_time = 12;
   DIC_time = 12;
   ALK_time = 12;

%    NO3_time   = 96;
%    NH4_time   = 96;
%    phy_time   = 96;
%    dia_time   = 96;
%    bac_time   = 96;
%    chla_time  = 96;
%    chld_time  = 96;
%    DON_time   = 96;
%    DOC_time   = 96;
%    DetN_time  = 96;
%    DetC_time  = 96;
%    sil_time   = 96;
%    LDetN_time = 96;
%    LDetC_time = 96;
%    LDetS_time = 96;
%    micro_time = 96;
%    meso_time  = 96;
%    FeD_time   = 96;
%    DetFe_time   = 96;
    
    %
    %----------------------------------------------------------------------
    % Create variables and add attributes
    %----------------------------------------------------------------------    
    %
    % Open boundary conditions time
    nccreate(fname,'bry_time', ...
            'Dimensions',{'bry_time' bry_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'bry_time','cycle_length',365.25);
    ncwriteatt(fname,'bry_time','long_name','open boundary conditions time');
    ncwriteatt(fname,'bry_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'bry_time','field','bry_time, scalar, series');
 
    % Free-surface
    nccreate(fname,'zeta_west', ...
            'Dimensions',{'eta_rho' eta_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'zeta_west','long_name','free-surface western boundary condition');
    ncwriteatt(fname,'zeta_west','units','meter');
    ncwriteatt(fname,'zeta_west','field','zeta_west, scalar, series');
    ncwriteatt(fname,'zeta_west','time','bry_time');

    nccreate(fname,'zeta_east', ...
            'Dimensions',{'eta_rho' eta_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'zeta_east','long_name','free-surface eastern boundary condition');
    ncwriteatt(fname,'zeta_east','units','meter');
    ncwriteatt(fname,'zeta_east','field','zeta_east, scalar, series');
    ncwriteatt(fname,'zeta_east','time','bry_time');

    nccreate(fname,'zeta_north', ...
            'Dimensions',{'xi_rho' xi_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'zeta_north','long_name','free-surface northern boundary condition');
    ncwriteatt(fname,'zeta_north','units','meter');
    ncwriteatt(fname,'zeta_north','field','zeta_north, scalar, series');
    ncwriteatt(fname,'zeta_north','time','bry_time');

    nccreate(fname,'zeta_south', ...
            'Dimensions',{'xi_rho' xi_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'zeta_south','long_name','free-surface southern boundary condition');
    ncwriteatt(fname,'zeta_south','units','meter');
    ncwriteatt(fname,'zeta_south','field','zeta_south, scalar, series');
    ncwriteatt(fname,'zeta_south','time','bry_time');


    % 2D u-direction velocity
    nccreate(fname,'ubar_west', ...
            'Dimensions',{'eta_rho' eta_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ubar_west','long_name','2D u-momentum western boundary condition');
    ncwriteatt(fname,'ubar_west','units','meter second-1');
    ncwriteatt(fname,'ubar_west','field','ubar_west, scalar, series');
    ncwriteatt(fname,'ubar_west','time','bry_time');

    nccreate(fname,'ubar_east', ...
            'Dimensions',{'eta_rho' eta_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ubar_east','long_name','2D u-momentum eastern boundary condition');
    ncwriteatt(fname,'ubar_east','units','meter second-1');
    ncwriteatt(fname,'ubar_east','field','ubar_east, scalar, series');
    ncwriteatt(fname,'ubar_east','time','bry_time');

    nccreate(fname,'ubar_north', ...
            'Dimensions',{'xi_u' xi_u 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ubar_north','long_name','2D u-momentum northern boundary condition');
    ncwriteatt(fname,'ubar_north','units','meter second-1');
    ncwriteatt(fname,'ubar_north','field','ubar_north, scalar, series');
    ncwriteatt(fname,'ubar_north','time','bry_time');

    nccreate(fname,'ubar_south', ...
            'Dimensions',{'xi_u' xi_u 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ubar_south','long_name','2D u-momentum southern boundary condition');
    ncwriteatt(fname,'ubar_south','units','meter second-1');
    ncwriteatt(fname,'ubar_south','field','ubar_south, scalar, series');
    ncwriteatt(fname,'ubar_south','time','bry_time');


    % 2D v-direction velocity
    nccreate(fname,'vbar_west', ...
            'Dimensions',{'eta_v' eta_v 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'vbar_west','long_name','2D v-momentum western boundary condition');
    ncwriteatt(fname,'vbar_west','units','meter second-1');
    ncwriteatt(fname,'vbar_west','field','vbar_west, scalar, series');
    ncwriteatt(fname,'vbar_west','time','bry_time');

    nccreate(fname,'vbar_east', ...
            'Dimensions',{'eta_v' eta_v 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'vbar_east','long_name','2D v-momentum eastern boundary condition');
    ncwriteatt(fname,'vbar_east','units','meter second-1');
    ncwriteatt(fname,'vbar_east','field','vbar_east, scalar, series');
    ncwriteatt(fname,'vbar_east','time','bry_time');

    nccreate(fname,'vbar_north', ...
            'Dimensions',{'xi_rho' xi_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'vbar_north','long_name','2D v-momentum northern boundary condition');
    ncwriteatt(fname,'vbar_north','units','meter second-1');
    ncwriteatt(fname,'vbar_north','field','vbar_north, scalar, series');
    ncwriteatt(fname,'vbar_north','time','bry_time');

    nccreate(fname,'vbar_south', ...
            'Dimensions',{'xi_rho' xi_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'vbar_south','long_name','2D v-momentum southern boundary condition');
    ncwriteatt(fname,'vbar_south','units','meter second-1');
    ncwriteatt(fname,'vbar_south','field','vbar_south, scalar, series');
    ncwriteatt(fname,'vbar_south','time','bry_time');


    % 3D u-direction velocity
    nccreate(fname,'u_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'u_west','long_name','3D u-momentum western boundary condition');
    ncwriteatt(fname,'u_west','units','meter second-1');
    ncwriteatt(fname,'u_west','field','u_west, scalar, series');
    ncwriteatt(fname,'u_west','time','bry_time');

    nccreate(fname,'u_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'u_east','long_name','3D u-momentum eastern boundary condition');
    ncwriteatt(fname,'u_east','units','meter second-1');
    ncwriteatt(fname,'u_east','field','u_east, scalar, series');
    ncwriteatt(fname,'u_east','time','bry_time');

    nccreate(fname,'u_north', ...
            'Dimensions',{'xi_u' xi_u 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'u_north','long_name','3D u-momentum northern boundary condition');
    ncwriteatt(fname,'u_north','units','meter second-1');
    ncwriteatt(fname,'u_north','field','u_north, scalar, series');
    ncwriteatt(fname,'u_north','time','bry_time');

    nccreate(fname,'u_south', ...
            'Dimensions',{'xi_u' xi_u 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'u_south','long_name','3D u-momentum southern boundary condition');
    ncwriteatt(fname,'u_south','units','meter second-1');
    ncwriteatt(fname,'u_south','field','u_south, scalar, series');
    ncwriteatt(fname,'u_south','time','bry_time');


    % 3D v-direction velocity
    nccreate(fname,'v_west', ...
            'Dimensions',{'eta_v' eta_v 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'v_west','long_name','3D v-momentum western boundary condition');
    ncwriteatt(fname,'v_west','units','meter second-1');
    ncwriteatt(fname,'v_west','field','v_west, scalar, series');
    ncwriteatt(fname,'v_west','time','bry_time');

    nccreate(fname,'v_east', ...
            'Dimensions',{'eta_v' eta_v 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'v_east','long_name','3D v-momentum eastern boundary condition');
    ncwriteatt(fname,'v_east','units','meter second-1');
    ncwriteatt(fname,'v_east','field','v_east, scalar, series');
    ncwriteatt(fname,'v_east','time','bry_time');

    nccreate(fname,'v_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'v_north','long_name','3D v-momentum northern boundary condition');
    ncwriteatt(fname,'v_north','units','meter second-1');
    ncwriteatt(fname,'v_north','field','v_north, scalar, series');
    ncwriteatt(fname,'v_north','time','bry_time');

    nccreate(fname,'v_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'v_south','long_name','3D v-momentum southern boundary condition');
    ncwriteatt(fname,'v_south','units','meter second-1');
    ncwriteatt(fname,'v_south','field','v_south, scalar, series');
    ncwriteatt(fname,'v_south','time','bry_time');


    % temperature
    nccreate(fname,'temp_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'temp_west','long_name','Potential temperature western boundary condition');
    ncwriteatt(fname,'temp_west','units','Celsius');
    ncwriteatt(fname,'temp_west','field','temp_west, scalar, series');
    ncwriteatt(fname,'temp_west','time','bry_time');

    nccreate(fname,'temp_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'temp_east','long_name','Potential temperature eastern boundary condition');
    ncwriteatt(fname,'temp_east','units','Celsius');
    ncwriteatt(fname,'temp_east','field','temp_east, scalar, series');
    ncwriteatt(fname,'temp_east','time','bry_time');

    nccreate(fname,'temp_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'temp_north','long_name','Potential temperature northern boundary condition');
    ncwriteatt(fname,'temp_north','units','Celsius');
    ncwriteatt(fname,'temp_north','field','temp_north, scalar, series');
    ncwriteatt(fname,'temp_north','time','bry_time');

    nccreate(fname,'temp_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'temp_south','long_name','Potential temperature southern boundary condition');
    ncwriteatt(fname,'temp_south','units','Celsius');
    ncwriteatt(fname,'temp_south','field','temp_south, scalar, series');
    ncwriteatt(fname,'temp_south','time','bry_time');


    % salinity
    nccreate(fname,'salt_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'salt_west','long_name','Salinity western boundary condition');
    ncwriteatt(fname,'salt_west','units','PSU');
    ncwriteatt(fname,'salt_west','field','salt_west, scalar, series');
    ncwriteatt(fname,'salt_west','time','bry_time');

    nccreate(fname,'salt_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'salt_east','long_name','Salinity eastern boundary condition');
    ncwriteatt(fname,'salt_east','units','PSU');
    ncwriteatt(fname,'salt_east','field','salt_east, scalar, series');
    ncwriteatt(fname,'salt_east','time','bry_time');

    nccreate(fname,'salt_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'salt_north','long_name','Salinity northern boundary condition');
    ncwriteatt(fname,'salt_north','units','PSU');
    ncwriteatt(fname,'salt_north','field','salt_north, scalar, series');
    ncwriteatt(fname,'salt_north','time','bry_time');

    nccreate(fname,'salt_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'bry_time' bry_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'salt_south','long_name','Salinity southern boundary condition');
    ncwriteatt(fname,'salt_south','units','PSU');
    ncwriteatt(fname,'salt_south','field','salt_south, scalar, series');
    ncwriteatt(fname,'salt_south','time','bry_time');


    % NO3
    nccreate(fname,'NO3_time', ...
            'Dimensions',{'NO3_time' NO3_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NO3_time','cycle_length',365.25);
    ncwriteatt(fname,'NO3_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'NO3_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'NO3_time','field','NO3_time, scalar, series');

    nccreate(fname,'NO3_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'NO3_time' NO3_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NO3_west','long_name','NO3 western boundary condition');
    ncwriteatt(fname,'NO3_west','units','mmol N m-3');
    ncwriteatt(fname,'NO3_west','field','NO3_west, scalar, series');
    ncwriteatt(fname,'NO3_west','time','NO3_time');

    nccreate(fname,'NO3_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'NO3_time' NO3_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NO3_east','long_name','NO3 eastern boundary condition');
    ncwriteatt(fname,'NO3_east','units','mmol N m-3');
    ncwriteatt(fname,'NO3_east','field','NO3_east, scalar, series');
    ncwriteatt(fname,'NO3_east','time','NO3_time');

    nccreate(fname,'NO3_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'NO3_time' NO3_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NO3_north','long_name','NO3 northern boundary condition');
    ncwriteatt(fname,'NO3_north','units','mmol N m-3');
    ncwriteatt(fname,'NO3_north','field','NO3_north, scalar, series');
    ncwriteatt(fname,'NO3_north','time','NO3_time');

    nccreate(fname,'NO3_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'NO3_time' NO3_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NO3_south','long_name','NO3 southern boundary condition');
    ncwriteatt(fname,'NO3_south','units','mmol N m-3');
    ncwriteatt(fname,'NO3_south','field','NO3_south, scalar, series');
    ncwriteatt(fname,'NO3_south','time','NO3_time');


    % NH4
    nccreate(fname,'NH4_time', ...
            'Dimensions',{'NH4_time' NH4_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NH4_time','cycle_length',365.25);
    ncwriteatt(fname,'NH4_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'NH4_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'NH4_time','field','NH4_time, scalar, series');

    nccreate(fname,'NH4_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'NH4_time' NH4_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NH4_west','long_name','NH4 western boundary condition');
    ncwriteatt(fname,'NH4_west','units','mmol N m-3');
    ncwriteatt(fname,'NH4_west','field','NH4_west, scalar, series');
    ncwriteatt(fname,'NH4_west','time','NH4_time');

    nccreate(fname,'NH4_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'NH4_time' NH4_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NH4_east','long_name','NH4 eastern boundary condition');
    ncwriteatt(fname,'NH4_east','units','mmol N m-3');
    ncwriteatt(fname,'NH4_east','field','NH4_east, scalar, series');
    ncwriteatt(fname,'NH4_east','time','NH4_time');

    nccreate(fname,'NH4_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'NH4_time' NH4_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NH4_north','long_name','NH4 northern boundary condition');
    ncwriteatt(fname,'NH4_north','units','mmol N m-3');
    ncwriteatt(fname,'NH4_north','field','NH4_north, scalar, series');
    ncwriteatt(fname,'NH4_north','time','NH4_time');

    nccreate(fname,'NH4_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'NH4_time' NH4_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NH4_south','long_name','NH4 southern boundary condition');
    ncwriteatt(fname,'NH4_south','units','mmol N m-3');
    ncwriteatt(fname,'NH4_south','field','NH4_south, scalar, series');
    ncwriteatt(fname,'NH4_south','time','NH4_time');


    % Bacteria
    nccreate(fname,'bac_time', ...
            'Dimensions',{'bac_time' bac_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'bac_time','cycle_length',365.25);
    ncwriteatt(fname,'bac_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'bac_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'bac_time','field','bac_time, scalar, series');

    nccreate(fname,'Bacteria_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'bac_time' bac_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Bacteria_west','long_name','Bacteria western boundary condition');
    ncwriteatt(fname,'Bacteria_west','units','mmol N m-3');
    ncwriteatt(fname,'Bacteria_west','field','Bacteria_west, scalar, series');
    ncwriteatt(fname,'Bacteria_west','time','bac_time');

    nccreate(fname,'Bacteria_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'bac_time' bac_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Bacteria_east','long_name','Bacteria eastern boundary condition');
    ncwriteatt(fname,'Bacteria_east','units','mmol N m-3');
    ncwriteatt(fname,'Bacteria_east','field','Bacteria_east, scalar, series');
    ncwriteatt(fname,'Bacteria_east','time','bac_time');

    nccreate(fname,'Bacteria_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'bac_time' bac_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Bacteria_north','long_name','Bacteria northern boundary condition');
    ncwriteatt(fname,'Bacteria_north','units','mmol N m-3');
    ncwriteatt(fname,'Bacteria_north','field','Bacteria_north, scalar, series');
    ncwriteatt(fname,'Bacteria_north','time','bac_time');

    nccreate(fname,'Bacteria_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'bac_time' bac_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Bacteria_south','long_name','Bacteria southern boundary condition');
    ncwriteatt(fname,'Bacteria_south','units','mmol N m-3');
    ncwriteatt(fname,'Bacteria_south','field','Bacteria_south, scalar, series');
    ncwriteatt(fname,'Bacteria_south','time','bac_time');

    
    % Phytoplankton
    nccreate(fname,'phy_time', ...
            'Dimensions',{'phy_time' phy_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'phy_time','cycle_length',365.25);
    ncwriteatt(fname,'phy_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'phy_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'phy_time','field','phy_time, scalar, series');

    nccreate(fname,'Phytoplankton_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'phy_time' phy_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Phytoplankton_west','long_name','Phytoplankton western boundary condition');
    ncwriteatt(fname,'Phytoplankton_west','units','mmol N m-3');
    ncwriteatt(fname,'Phytoplankton_west','field','Phytoplankton_west, scalar, series');
    ncwriteatt(fname,'Phytoplankton_west','time','phy_time');

    nccreate(fname,'Phytoplankton_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'phy_time' phy_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Phytoplankton_east','long_name','Phytoplankton eastern boundary condition');
    ncwriteatt(fname,'Phytoplankton_east','units','mmol N m-3');
    ncwriteatt(fname,'Phytoplankton_east','field','Phytoplankton_east, scalar, series');
    ncwriteatt(fname,'Phytoplankton_east','time','phy_time');

    nccreate(fname,'Phytoplankton_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'phy_time' phy_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Phytoplankton_north','long_name','Phytoplankton northern boundary condition');
    ncwriteatt(fname,'Phytoplankton_north','units','mmol N m-3');
    ncwriteatt(fname,'Phytoplankton_north','field','Phytoplankton_north, scalar, series');
    ncwriteatt(fname,'Phytoplankton_north','time','phy_time');

    nccreate(fname,'Phytoplankton_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'phy_time' phy_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Phytoplankton_south','long_name','Phytoplankton southern boundary condition');
    ncwriteatt(fname,'Phytoplankton_south','units','mmol N m-3');
    ncwriteatt(fname,'Phytoplankton_south','field','Phytoplankton_south, scalar, series');
    ncwriteatt(fname,'Phytoplankton_south','time','phy_time');

    
    % Diatoms
    nccreate(fname,'dia_time', ...
            'Dimensions',{'dia_time' dia_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'dia_time','cycle_length',365.25);
    ncwriteatt(fname,'dia_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'dia_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'dia_time','field','dia_time, scalar, series');

    nccreate(fname,'Diatoms_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'dia_time' dia_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Diatoms_west','long_name','Diatoms western boundary condition');
    ncwriteatt(fname,'Diatoms_west','units','mmol N m-3');
    ncwriteatt(fname,'Diatoms_west','field','Diatoms_west, scalar, series');
    ncwriteatt(fname,'Diatoms_west','time','dia_time');

    nccreate(fname,'Diatoms_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'dia_time' dia_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Diatoms_east','long_name','Diatoms eastern boundary condition');
    ncwriteatt(fname,'Diatoms_east','units','mmol N m-3');
    ncwriteatt(fname,'Diatoms_east','field','Diatoms_east, scalar, series');
    ncwriteatt(fname,'Diatoms_east','time','dia_time');

    nccreate(fname,'Diatoms_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'dia_time' dia_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Diatoms_north','long_name','Diatoms northern boundary condition');
    ncwriteatt(fname,'Diatoms_north','units','mmol N m-3');
    ncwriteatt(fname,'Diatoms_north','field','Diatoms_north, scalar, series');
    ncwriteatt(fname,'Diatoms_north','time','dia_time');

    nccreate(fname,'Diatoms_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'dia_time' dia_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Diatoms_south','long_name','Diatoms southern boundary condition');
    ncwriteatt(fname,'Diatoms_south','units','mmol N m-3');
    ncwriteatt(fname,'Diatoms_south','field','Diatoms_south, scalar, series');
    ncwriteatt(fname,'Diatoms_south','time','dia_time');

    
    % Chlorophyll-a
    nccreate(fname,'chla_time', ...
            'Dimensions',{'chla_time' chla_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'chla_time','cycle_length',365.25);
    ncwriteatt(fname,'chla_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'chla_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'chla_time','field','chla_time, scalar, series');

    nccreate(fname,'Chlorophyll_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'chla_time' chla_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Chlorophyll_west','long_name','Chlorophyll western boundary condition');
    ncwriteatt(fname,'Chlorophyll_west','units','mg Chl m-3');
    ncwriteatt(fname,'Chlorophyll_west','field','Chlorophyll_west, scalar, series');
    ncwriteatt(fname,'Chlorophyll_west','time','chla_time');

    nccreate(fname,'Chlorophyll_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'chla_time' chla_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Chlorophyll_east','long_name','Chlorophyll eastern boundary condition');
    ncwriteatt(fname,'Chlorophyll_east','units','mg Chl m-3');
    ncwriteatt(fname,'Chlorophyll_east','field','Chlorophyll_east, scalar, series');
    ncwriteatt(fname,'Chlorophyll_east','time','chla_time');

    nccreate(fname,'Chlorophyll_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'chla_time' chla_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Chlorophyll_north','long_name','Chlorophyll northern boundary condition');
    ncwriteatt(fname,'Chlorophyll_north','units','mg Chl m-3');
    ncwriteatt(fname,'Chlorophyll_north','field','Chlorophyll_north, scalar, series');
    ncwriteatt(fname,'Chlorophyll_north','time','chla_time');

    nccreate(fname,'Chlorophyll_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'chla_time' chla_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Chlorophyll_south','long_name','Chlorophyll southern boundary condition');
    ncwriteatt(fname,'Chlorophyll_south','units','mg Chl m-3');
    ncwriteatt(fname,'Chlorophyll_south','field','Chlorophyll_south, scalar, series');
    ncwriteatt(fname,'Chlorophyll_south','time','chla_time');

    
    % Diatom Chlorophyll
    nccreate(fname,'chld_time', ...
            'Dimensions',{'chld_time' chld_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'chld_time','cycle_length',365.25);
    ncwriteatt(fname,'chld_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'chld_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'chld_time','field','chld_time, scalar, series');

    nccreate(fname,'ChlD_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'chld_time' chld_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ChlD_west','long_name','ChlD western boundary condition');
    ncwriteatt(fname,'ChlD_west','units','mg Chl m-3');
    ncwriteatt(fname,'ChlD_west','field','ChlD_west, scalar, series');
    ncwriteatt(fname,'ChlD_west','time','chld_time');

    nccreate(fname,'ChlD_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'chld_time' chld_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ChlD_east','long_name','ChlD eastern boundary condition');
    ncwriteatt(fname,'ChlD_east','units','mg Chl m-3');
    ncwriteatt(fname,'ChlD_east','field','ChlD_east, scalar, series');
    ncwriteatt(fname,'ChlD_east','time','chld_time');

    nccreate(fname,'ChlD_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'chld_time' chld_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ChlD_north','long_name','ChlD northern boundary condition');
    ncwriteatt(fname,'ChlD_north','units','mg Chl m-3');
    ncwriteatt(fname,'ChlD_north','field','ChlD_north, scalar, series');
    ncwriteatt(fname,'ChlD_north','time','chld_time');

    nccreate(fname,'ChlD_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'chld_time' chld_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ChlD_south','long_name','ChlD southern boundary condition');
    ncwriteatt(fname,'ChlD_south','units','mg Chl m-3');
    ncwriteatt(fname,'ChlD_south','field','ChlD_south, scalar, series');
    ncwriteatt(fname,'ChlD_south','time','chld_time');


    % DON
    nccreate(fname,'DON_time', ...
            'Dimensions',{'DON_time' DON_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DON_time','cycle_length',365.25);
    ncwriteatt(fname,'DON_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'DON_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'DON_time','field','DON_time, scalar, series');

    nccreate(fname,'DON_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'DON_time' DON_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DON_west','long_name','DON western boundary condition');
    ncwriteatt(fname,'DON_west','units','mmol N m-3');
    ncwriteatt(fname,'DON_west','field','DON_west, scalar, series');
    ncwriteatt(fname,'DON_west','time','DON_time');

    nccreate(fname,'DON_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'DON_time' DON_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DON_east','long_name','DON eastern boundary condition');
    ncwriteatt(fname,'DON_east','units','mmol N m-3');
    ncwriteatt(fname,'DON_east','field','DON_east, scalar, series');
    ncwriteatt(fname,'DON_east','time','DON_time');

    nccreate(fname,'DON_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'DON_time' DON_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DON_north','long_name','DON northern boundary condition');
    ncwriteatt(fname,'DON_north','units','mmol N m-3');
    ncwriteatt(fname,'DON_north','field','DON_north, scalar, series');
    ncwriteatt(fname,'DON_north','time','DON_time');

    nccreate(fname,'DON_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'DON_time' DON_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DON_south','long_name','DON southern boundary condition');
    ncwriteatt(fname,'DON_south','units','mmol N m-3');
    ncwriteatt(fname,'DON_south','field','DON_south, scalar, series');
    ncwriteatt(fname,'DON_south','time','DON_time');


    % DOC
    nccreate(fname,'DOC_time', ...
            'Dimensions',{'DOC_time' DOC_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DOC_time','cycle_length',365.25);
    ncwriteatt(fname,'DOC_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'DOC_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'DOC_time','field','DOC_time, scalar, series');

    nccreate(fname,'DOC_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'DOC_time' DOC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DOC_west','long_name','DOC western boundary condition');
    ncwriteatt(fname,'DOC_west','units','mmol C m-3');
    ncwriteatt(fname,'DOC_west','field','DOC_west, scalar, series');
    ncwriteatt(fname,'DOC_west','time','DOC_time');

    nccreate(fname,'DOC_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'DOC_time' DOC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DOC_east','long_name','DOC eastern boundary condition');
    ncwriteatt(fname,'DOC_east','units','mmol C m-3');
    ncwriteatt(fname,'DOC_east','field','DOC_east, scalar, series');
    ncwriteatt(fname,'DOC_east','time','DOC_time');

    nccreate(fname,'DOC_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'DOC_time' DOC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DOC_north','long_name','DOC northern boundary condition');
    ncwriteatt(fname,'DOC_north','units','mmol C m-3');
    ncwriteatt(fname,'DOC_north','field','DOC_north, scalar, series');
    ncwriteatt(fname,'DOC_north','time','DOC_time');

    nccreate(fname,'DOC_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'DOC_time' DOC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DOC_south','long_name','DOC southern boundary condition');
    ncwriteatt(fname,'DOC_south','units','mmol C m-3');
    ncwriteatt(fname,'DOC_south','field','DOC_south, scalar, series');
    ncwriteatt(fname,'DOC_south','time','DOC_time');


    % DetN
    nccreate(fname,'DetN_time', ...
            'Dimensions',{'DetN_time' DetN_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetN_time','cycle_length',365.25);
    ncwriteatt(fname,'DetN_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'DetN_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'DetN_time','field','DetN_time, scalar, series');

    nccreate(fname,'DetN_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'DetN_time' DetN_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetN_west','long_name','DetN western boundary condition');
    ncwriteatt(fname,'DetN_west','units','mmol N m-3');
    ncwriteatt(fname,'DetN_west','field','DetN_west, scalar, series');
    ncwriteatt(fname,'DetN_west','time','DetN_time');

    nccreate(fname,'DetN_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'DetN_time' DetN_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetN_east','long_name','DetN eastern boundary condition');
    ncwriteatt(fname,'DetN_east','units','mmol N m-3');
    ncwriteatt(fname,'DetN_east','field','DetN_east, scalar, series');
    ncwriteatt(fname,'DetN_east','time','DetN_time');

    nccreate(fname,'DetN_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'DetN_time' DetN_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetN_north','long_name','DetN northern boundary condition');
    ncwriteatt(fname,'DetN_north','units','mmol N m-3');
    ncwriteatt(fname,'DetN_north','field','DetN_north, scalar, series');
    ncwriteatt(fname,'DetN_north','time','DetN_time');

    nccreate(fname,'DetN_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'DetN_time' DetN_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetN_south','long_name','DetN southern boundary condition');
    ncwriteatt(fname,'DetN_south','units','mmol N m-3');
    ncwriteatt(fname,'DetN_south','field','DetN_south, scalar, series');
    ncwriteatt(fname,'DetN_south','time','DetN_time');


    % DetC
    nccreate(fname,'DetC_time', ...
            'Dimensions',{'DetC_time' DetC_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetC_time','cycle_length',365.25);
    ncwriteatt(fname,'DetC_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'DetC_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'DetC_time','field','DetC_time, scalar, series');

    nccreate(fname,'DetC_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'DetC_time' DetC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetC_west','long_name','DetC western boundary condition');
    ncwriteatt(fname,'DetC_west','units','mmol C m-3');
    ncwriteatt(fname,'DetC_west','field','DetC_west, scalar, series');
    ncwriteatt(fname,'DetC_west','time','DetC_time');

    nccreate(fname,'DetC_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'DetC_time' DetC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetC_east','long_name','DetC eastern boundary condition');
    ncwriteatt(fname,'DetC_east','units','mmol C m-3');
    ncwriteatt(fname,'DetC_east','field','DetC_east, scalar, series');
    ncwriteatt(fname,'DetC_east','time','DetC_time');

    nccreate(fname,'DetC_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'DetC_time' DetC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetC_north','long_name','DetC northern boundary condition');
    ncwriteatt(fname,'DetC_north','units','mmol C m-3');
    ncwriteatt(fname,'DetC_north','field','DetC_north, scalar, series');
    ncwriteatt(fname,'DetC_north','time','DetC_time');

    nccreate(fname,'DetC_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'DetC_time' DetC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetC_south','long_name','DetC southern boundary condition');
    ncwriteatt(fname,'DetC_south','units','mmol C m-3');
    ncwriteatt(fname,'DetC_south','field','DetC_south, scalar, series');
    ncwriteatt(fname,'DetC_south','time','DetC_time');


    % Silica
    nccreate(fname,'sil_time', ...
            'Dimensions',{'sil_time' sil_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'sil_time','cycle_length',365.25);
    ncwriteatt(fname,'sil_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'sil_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'sil_time','field','sil_time, scalar, series');

    nccreate(fname,'Silica_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'sil_time' sil_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Silica_west','long_name','Silica western boundary condition');
    ncwriteatt(fname,'Silica_west','units','mmol Si m-3');
    ncwriteatt(fname,'Silica_west','field','Silica_west, scalar, series');
    ncwriteatt(fname,'Silica_west','time','sil_time');

    nccreate(fname,'Silica_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'sil_time' sil_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Silica_east','long_name','Silica eastern boundary condition');
    ncwriteatt(fname,'Silica_east','units','mmol Si m-3');
    ncwriteatt(fname,'Silica_east','field','Silica_east, scalar, series');
    ncwriteatt(fname,'Silica_east','time','sil_time');

    nccreate(fname,'Silica_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'sil_time' sil_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Silica_north','long_name','Silica northern boundary condition');
    ncwriteatt(fname,'Silica_north','units','mmol Si m-3');
    ncwriteatt(fname,'Silica_north','field','Silica_north, scalar, series');
    ncwriteatt(fname,'Silica_north','time','sil_time');

    nccreate(fname,'Silica_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'sil_time' sil_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Silica_south','long_name','Silica southern boundary condition');
    ncwriteatt(fname,'Silica_south','units','mmol Si m-3');
    ncwriteatt(fname,'Silica_south','field','Silica_south, scalar, series');
    ncwriteatt(fname,'Silica_south','time','sil_time');


    % LDetN
    nccreate(fname,'LDetN_time', ...
            'Dimensions',{'LDetN_time' LDetN_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetN_time','cycle_length',365.25);
    ncwriteatt(fname,'LDetN_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'LDetN_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'LDetN_time','field','LDetN_time, scalar, series');

    nccreate(fname,'LDetN_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'LDetN_time' LDetN_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetN_west','long_name','LDetN western boundary condition');
    ncwriteatt(fname,'LDetN_west','units','mmol N m-3');
    ncwriteatt(fname,'LDetN_west','field','LDetN_west, scalar, series');
    ncwriteatt(fname,'LDetN_west','time','LDetN_time');

    nccreate(fname,'LDetN_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'LDetN_time' LDetN_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetN_east','long_name','LDetN eastern boundary condition');
    ncwriteatt(fname,'LDetN_east','units','mmol N m-3');
    ncwriteatt(fname,'LDetN_east','field','LDetN_east, scalar, series');
    ncwriteatt(fname,'LDetN_east','time','LDetN_time');

    nccreate(fname,'LDetN_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'LDetN_time' LDetN_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetN_north','long_name','LDetN northern boundary condition');
    ncwriteatt(fname,'LDetN_north','units','mmol N m-3');
    ncwriteatt(fname,'LDetN_north','field','LDetN_north, scalar, series');
    ncwriteatt(fname,'LDetN_north','time','LDetN_time');

    nccreate(fname,'LDetN_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'LDetN_time' LDetN_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetN_south','long_name','LDetN southern boundary condition');
    ncwriteatt(fname,'LDetN_south','units','mmol N m-3');
    ncwriteatt(fname,'LDetN_south','field','LDetN_south, scalar, series');
    ncwriteatt(fname,'LDetN_south','time','LDetN_time');


    % LDetC
    nccreate(fname,'LDetC_time', ...
            'Dimensions',{'LDetC_time' LDetC_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetC_time','cycle_length',365.25);
    ncwriteatt(fname,'LDetC_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'LDetC_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'LDetC_time','field','LDetC_time, scalar, series');

    nccreate(fname,'LDetC_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'LDetC_time' LDetC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetC_west','long_name','LDetC western boundary condition');
    ncwriteatt(fname,'LDetC_west','units','mmol C m-3');
    ncwriteatt(fname,'LDetC_west','field','LDetC_west, scalar, series');
    ncwriteatt(fname,'LDetC_west','time','LDetC_time');

    nccreate(fname,'LDetC_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'LDetC_time' LDetC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetC_east','long_name','LDetC eastern boundary condition');
    ncwriteatt(fname,'LDetC_east','units','mmol C m-3');
    ncwriteatt(fname,'LDetC_east','field','LDetC_east, scalar, series');
    ncwriteatt(fname,'LDetC_east','time','LDetC_time');

    nccreate(fname,'LDetC_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'LDetC_time' LDetC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetC_north','long_name','LDetC northern boundary condition');
    ncwriteatt(fname,'LDetC_north','units','mmol C m-3');
    ncwriteatt(fname,'LDetC_north','field','LDetC_north, scalar, series');
    ncwriteatt(fname,'LDetC_north','time','LDetC_time');

    nccreate(fname,'LDetC_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'LDetC_time' LDetC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetC_south','long_name','LDetC southern boundary condition');
    ncwriteatt(fname,'LDetC_south','units','mmol C m-3');
    ncwriteatt(fname,'LDetC_south','field','LDetC_south, scalar, series');
    ncwriteatt(fname,'LDetC_south','time','LDetC_time');


    % LDetS
    nccreate(fname,'LDetS_time', ...
            'Dimensions',{'LDetS_time' LDetS_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetS_time','cycle_length',365.25);
    ncwriteatt(fname,'LDetS_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'LDetS_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'LDetS_time','field','LDetS_time, scalar, series');

    nccreate(fname,'LDetS_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'LDetS_time' LDetS_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetS_west','long_name','LDetS western boundary condition');
    ncwriteatt(fname,'LDetS_west','units','mmol Si m-3');
    ncwriteatt(fname,'LDetS_west','field','LDetS_west, scalar, series');
    ncwriteatt(fname,'LDetS_west','time','LDetS_time');

    nccreate(fname,'LDetS_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'LDetS_time' LDetS_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetS_east','long_name','LDetS eastern boundary condition');
    ncwriteatt(fname,'LDetS_east','units','mmol Si m-3');
    ncwriteatt(fname,'LDetS_east','field','LDetS_east, scalar, series');
    ncwriteatt(fname,'LDetS_east','time','LDetS_time');

    nccreate(fname,'LDetS_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'LDetS_time' LDetS_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetS_north','long_name','LDetS northern boundary condition');
    ncwriteatt(fname,'LDetS_north','units','mmol Si m-3');
    ncwriteatt(fname,'LDetS_north','field','LDetS_north, scalar, series');
    ncwriteatt(fname,'LDetS_north','time','LDetS_time');

    nccreate(fname,'LDetS_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'LDetS_time' LDetS_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetS_south','long_name','LDetS southern boundary condition');
    ncwriteatt(fname,'LDetS_south','units','mmol Si m-3');
    ncwriteatt(fname,'LDetS_south','field','LDetS_south, scalar, series');
    ncwriteatt(fname,'LDetS_south','time','LDetS_time');


    % Microzooplankton
    nccreate(fname,'micro_time', ...
            'Dimensions',{'micro_time' micro_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'micro_time','cycle_length',365.25);
    ncwriteatt(fname,'micro_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'micro_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'micro_time','field','micro_time, scalar, series');

    nccreate(fname,'Microzooplankton_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'micro_time' micro_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Microzooplankton_west','long_name','Microzooplankton western boundary condition');
    ncwriteatt(fname,'Microzooplankton_west','units','mmol N m-3');
    ncwriteatt(fname,'Microzooplankton_west','field','Microzooplankton_west, scalar, series');
    ncwriteatt(fname,'Microzooplankton_west','time','micro_time');

    nccreate(fname,'Microzooplankton_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'micro_time' micro_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Microzooplankton_east','long_name','Microzooplankton eastern boundary condition');
    ncwriteatt(fname,'Microzooplankton_east','units','mmol N m-3');
    ncwriteatt(fname,'Microzooplankton_east','field','Microzooplankton_east, scalar, series');
    ncwriteatt(fname,'Microzooplankton_east','time','micro_time');

    nccreate(fname,'Microzooplankton_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'micro_time' micro_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Microzooplankton_north','long_name','Microzooplankton northern boundary condition');
    ncwriteatt(fname,'Microzooplankton_north','units','mmol N m-3');
    ncwriteatt(fname,'Microzooplankton_north','field','Microzooplankton_north, scalar, series');
    ncwriteatt(fname,'Microzooplankton_north','time','micro_time');

    nccreate(fname,'Microzooplankton_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'micro_time' micro_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Microzooplankton_south','long_name','Microzooplankton southern boundary condition');
    ncwriteatt(fname,'Microzooplankton_south','units','mmol N m-3');
    ncwriteatt(fname,'Microzooplankton_south','field','Microzooplankton_south, scalar, series');
    ncwriteatt(fname,'Microzooplankton_south','time','micro_time');


    % Mesozooplankton
    nccreate(fname,'meso_time', ...
            'Dimensions',{'meso_time' meso_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'meso_time','cycle_length',365.25);
    ncwriteatt(fname,'meso_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'meso_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'meso_time','field','meso_time, scalar, series');

    nccreate(fname,'Mesozooplankton_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'meso_time' meso_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Mesozooplankton_west','long_name','Mesozooplankton western boundary condition');
    ncwriteatt(fname,'Mesozooplankton_west','units','mmol N m-3');
    ncwriteatt(fname,'Mesozooplankton_west','field','Mesozooplankton_west, scalar, series');
    ncwriteatt(fname,'Mesozooplankton_west','time','meso_time');

    nccreate(fname,'Mesozooplankton_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'meso_time' meso_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Mesozooplankton_east','long_name','Mesozooplankton eastern boundary condition');
    ncwriteatt(fname,'Mesozooplankton_east','units','mmol N m-3');
    ncwriteatt(fname,'Mesozooplankton_east','field','Mesozooplankton_east, scalar, series');
    ncwriteatt(fname,'Mesozooplankton_east','time','meso_time');

    nccreate(fname,'Mesozooplankton_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'meso_time' meso_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Mesozooplankton_north','long_name','Mesozooplankton northern boundary condition');
    ncwriteatt(fname,'Mesozooplankton_north','units','mmol N m-3');
    ncwriteatt(fname,'Mesozooplankton_north','field','Mesozooplankton_north, scalar, series');
    ncwriteatt(fname,'Mesozooplankton_north','time','meso_time');

    nccreate(fname,'Mesozooplankton_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'meso_time' meso_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Mesozooplankton_south','long_name','Mesozooplankton southern boundary condition');
    ncwriteatt(fname,'Mesozooplankton_south','units','mmol N m-3');
    ncwriteatt(fname,'Mesozooplankton_south','field','Mesozooplankton_south, scalar, series');
    ncwriteatt(fname,'Mesozooplankton_south','time','meso_time');


    % Dissolved Iron
    nccreate(fname,'FeD_time', ...
            'Dimensions',{'FeD_time' FeD_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'FeD_time','cycle_length',365.25);
    ncwriteatt(fname,'FeD_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'FeD_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'FeD_time','field','FeD_time, scalar, series');

    nccreate(fname,'FeD_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'FeD_time' FeD_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'FeD_west','long_name','FeD western boundary condition');
    ncwriteatt(fname,'FeD_west','units','umol Fe m-3');
    ncwriteatt(fname,'FeD_west','field','FeD_west, scalar, series');
    ncwriteatt(fname,'FeD_west','time','FeD_time');

    nccreate(fname,'FeD_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'FeD_time' FeD_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'FeD_east','long_name','FeD eastern boundary condition');
    ncwriteatt(fname,'FeD_east','units','umol Fe m-3');
    ncwriteatt(fname,'FeD_east','field','FeD_east, scalar, series');
    ncwriteatt(fname,'FeD_east','time','FeD_time');

    nccreate(fname,'FeD_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'FeD_time' FeD_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'FeD_north','long_name','FeD northern boundary condition');
    ncwriteatt(fname,'FeD_north','units','umol Fe m-3');
    ncwriteatt(fname,'FeD_north','field','FeD_north, scalar, series');
    ncwriteatt(fname,'FeD_north','time','FeD_time');

    nccreate(fname,'FeD_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'FeD_time' FeD_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'FeD_south','long_name','FeD southern boundary condition');
    ncwriteatt(fname,'FeD_south','units','umol Fe m-3');
    ncwriteatt(fname,'FeD_south','field','FeD_south, scalar, series');
    ncwriteatt(fname,'FeD_south','time','FeD_time');


    % Detrital iron
    nccreate(fname,'DetFe_time', ...
            'Dimensions',{'DetFe_time' DetFe_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetFe_time','cycle_length',365.25);
    ncwriteatt(fname,'DetFe_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'DetFe_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'DetFe_time','field','DetFe_time, scalar, series');

    nccreate(fname,'DetFe_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'DetFe_time' DetFe_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetFe_west','long_name','DetFe western boundary condition');
    ncwriteatt(fname,'DetFe_west','units','umol Fe m-3');
    ncwriteatt(fname,'DetFe_west','field','DetFe_west, scalar, series');
    ncwriteatt(fname,'DetFe_west','time','DetFe_time');

    nccreate(fname,'DetFe_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'DetFe_time' DetFe_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetFe_east','long_name','DetFe eastern boundary condition');
    ncwriteatt(fname,'DetFe_east','units','umol Fe m-3');
    ncwriteatt(fname,'DetFe_east','field','DetFe_east, scalar, series');
    ncwriteatt(fname,'DetFe_east','time','DetFe_time');

    nccreate(fname,'DetFe_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'DetFe_time' DetFe_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetFe_north','long_name','DetFe northern boundary condition');
    ncwriteatt(fname,'DetFe_north','units','umol Fe m-3');
    ncwriteatt(fname,'DetFe_north','field','DetFe_north, scalar, series');
    ncwriteatt(fname,'DetFe_north','time','DetFe_time');

    nccreate(fname,'DetFe_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'DetFe_time' DetFe_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetFe_south','long_name','DetFe southern boundary condition');
    ncwriteatt(fname,'DetFe_south','units','umol Fe m-3');
    ncwriteatt(fname,'DetFe_south','field','DetFe_south, scalar, series');
    ncwriteatt(fname,'DetFe_south','time','DetFe_time');

    
    % Large detrital iron
    nccreate(fname,'LDetFe_time', ...
            'Dimensions',{'LDetFe_time' LDetFe_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetFe_time','cycle_length',365.25);
    ncwriteatt(fname,'LDetFe_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'LDetFe_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'LDetFe_time','field','LDetFe_time, scalar, series');

    nccreate(fname,'LDetFe_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'LDetFe_time' LDetFe_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetFe_west','long_name','LDetFe western boundary condition');
    ncwriteatt(fname,'LDetFe_west','units','umol Fe m-3');
    ncwriteatt(fname,'LDetFe_west','field','LDetFe_west, scalar, series');
    ncwriteatt(fname,'LDetFe_west','time','LDetFe_time');

    nccreate(fname,'LDetFe_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'LDetFe_time' LDetFe_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetFe_east','long_name','LDetFe eastern boundary condition');
    ncwriteatt(fname,'LDetFe_east','units','umol Fe m-3');
    ncwriteatt(fname,'LDetFe_east','field','LDetFe_east, scalar, series');
    ncwriteatt(fname,'LDetFe_east','time','LDetFe_time');

    nccreate(fname,'LDetFe_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'LDetFe_time' LDetFe_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetFe_north','long_name','LDetFe northern boundary condition');
    ncwriteatt(fname,'LDetFe_north','units','umol Fe m-3');
    ncwriteatt(fname,'LDetFe_north','field','LDetFe_north, scalar, series');
    ncwriteatt(fname,'LDetFe_north','time','LDetFe_time');

    nccreate(fname,'LDetFe_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'LDetFe_time' LDetFe_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetFe_south','long_name','LDetFe southern boundary condition');
    ncwriteatt(fname,'LDetFe_south','units','umol Fe m-3');
    ncwriteatt(fname,'LDetFe_south','field','LDetFe_south, scalar, series');
    ncwriteatt(fname,'LDetFe_south','time','LDetFe_time');


    %Coccolithophores
    nccreate(fname,'coccos_time', ...
            'Dimensions',{'coccos_time' coccos_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'coccos_time','cycle_length',365.25);
    ncwriteatt(fname,'coccos_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'coccos_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'coccos_time','field','coccos_time, scalar, series');

    nccreate(fname,'Coccos_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'coccos_time' coccos_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Coccos_west','long_name','Coccos western boundary condition');
    ncwriteatt(fname,'Coccos_west','units','mmol N m-3');
    ncwriteatt(fname,'Coccos_west','field','Coccos_west, scalar, series');
    ncwriteatt(fname,'Coccos_west','time','coccos_time');

    nccreate(fname,'Coccos_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'coccos_time' coccos_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Coccos_east','long_name','Coccos eastern boundary condition');
    ncwriteatt(fname,'Coccos_east','units','mmol N m-3');
    ncwriteatt(fname,'Coccos_east','field','Coccos_east, scalar, series');
    ncwriteatt(fname,'Coccos_east','time','coccos_time');

    nccreate(fname,'Coccos_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'coccos_time' coccos_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Coccos_north','long_name','Coccos northern boundary condition');
    ncwriteatt(fname,'Coccos_north','units','mmol N m-3');
    ncwriteatt(fname,'Coccos_north','field','Coccos_north, scalar, series');
    ncwriteatt(fname,'Coccos_north','time','coccos_time');

    nccreate(fname,'Coccos_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'coccos_time' coccos_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Coccos_south','long_name','Coccos southern boundary condition');
    ncwriteatt(fname,'Coccos_south','units','mmol N m-3');
    ncwriteatt(fname,'Coccos_south','field','Coccos_south, scalar, series');
    ncwriteatt(fname,'Coccos_south','time','coccos_time');

    
    % Coccolithophore chlorophyll
    nccreate(fname,'ChlC_time', ...
            'Dimensions',{'ChlC_time' ChlC_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ChlC_time','cycle_length',365.25);
    ncwriteatt(fname,'ChlC_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'ChlC_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'ChlC_time','field','ChlC_time, scalar, series');

    nccreate(fname,'ChlC_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'ChlC_time' ChlC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ChlC_west','long_name','ChlC western boundary condition');
    ncwriteatt(fname,'ChlC_west','units','milligram m-3');
    ncwriteatt(fname,'ChlC_west','field','ChlC_west, scalar, series');
    ncwriteatt(fname,'ChlC_west','time','ChlC_time');

    nccreate(fname,'ChlC_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'ChlC_time' ChlC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ChlC_east','long_name','ChlC eastern boundary condition');
    ncwriteatt(fname,'ChlC_east','units','milligram m-3');
    ncwriteatt(fname,'ChlC_east','field','ChlC_east, scalar, series');
    ncwriteatt(fname,'ChlC_east','time','ChlC_time');

    nccreate(fname,'ChlC_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'ChlC_time' ChlC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ChlC_north','long_name','ChlC northern boundary condition');
    ncwriteatt(fname,'ChlC_north','units','milligram m-3');
    ncwriteatt(fname,'ChlC_north','field','ChlC_north, scalar, series');
    ncwriteatt(fname,'ChlC_north','time','ChlC_time');

    nccreate(fname,'ChlC_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'ChlC_time' ChlC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ChlC_south','long_name','ChlC southern boundary condition');
    ncwriteatt(fname,'ChlC_south','units','milligram m-3');
    ncwriteatt(fname,'ChlC_south','field','ChlC_south, scalar, series');
    ncwriteatt(fname,'ChlC_south','time','ChlC_time');
    
    
    % Particulate Inorganic Carbon
    nccreate(fname,'PIC_time', ...
            'Dimensions',{'PIC_time' PIC_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'PIC_time','cycle_length',365.25);
    ncwriteatt(fname,'PIC_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'PIC_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'PIC_time','field','PIC_time, scalar, series');

    nccreate(fname,'PIC_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'PIC_time' PIC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'PIC_west','long_name','PIC western boundary condition');
    ncwriteatt(fname,'PIC_west','units','mmol C m-3');
    ncwriteatt(fname,'PIC_west','field','PIC_west, scalar, series');
    ncwriteatt(fname,'PIC_west','time','PIC_time');

    nccreate(fname,'PIC_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'PIC_time' PIC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'PIC_east','long_name','PIC eastern boundary condition');
    ncwriteatt(fname,'PIC_east','units','mmol C m-3');
    ncwriteatt(fname,'PIC_east','field','PIC_east, scalar, series');
    ncwriteatt(fname,'PIC_east','time','PIC_time');

    nccreate(fname,'PIC_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'PIC_time' PIC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'PIC_north','long_name','PIC northern boundary condition');
    ncwriteatt(fname,'PIC_north','units','mmol C m-3');
    ncwriteatt(fname,'PIC_north','field','PIC_north, scalar, series');
    ncwriteatt(fname,'PIC_north','time','PIC_time');

    nccreate(fname,'PIC_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'PIC_time' PIC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'PIC_south','long_name','PIC southern boundary condition');
    ncwriteatt(fname,'PIC_south','units','mmol C m-3');
    ncwriteatt(fname,'PIC_south','field','PIC_south, scalar, series');
    ncwriteatt(fname,'PIC_south','time','PIC_time');
    
    
    % Dissolved Inorganic Carbon
    nccreate(fname,'DIC_time', ...
            'Dimensions',{'DIC_time' DIC_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DIC_time','cycle_length',365.25);
    ncwriteatt(fname,'DIC_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'DIC_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'DIC_time','field','DIC_time, scalar, series');

    nccreate(fname,'DIC_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'DIC_time' DIC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DIC_west','long_name','DIC western boundary condition');
    ncwriteatt(fname,'DIC_west','units','mmol C m-3');
    ncwriteatt(fname,'DIC_west','field','DIC_west, scalar, series');
    ncwriteatt(fname,'DIC_west','time','DIC_time');

    nccreate(fname,'DIC_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'DIC_time' DIC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DIC_east','long_name','DIC eastern boundary condition');
    ncwriteatt(fname,'DIC_east','units','mmol C m-3');
    ncwriteatt(fname,'DIC_east','field','DIC_east, scalar, series');
    ncwriteatt(fname,'DIC_east','time','DIC_time');

    nccreate(fname,'DIC_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'DIC_time' DIC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DIC_north','long_name','DIC northern boundary condition');
    ncwriteatt(fname,'DIC_north','units','mmol C m-3');
    ncwriteatt(fname,'DIC_north','field','DIC_north, scalar, series');
    ncwriteatt(fname,'DIC_north','time','DIC_time');

    nccreate(fname,'DIC_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'DIC_time' DIC_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DIC_south','long_name','DIC southern boundary condition');
    ncwriteatt(fname,'DIC_south','units','mmol C m-3');
    ncwriteatt(fname,'DIC_south','field','DIC_south, scalar, series');
    ncwriteatt(fname,'DIC_south','time','DIC_time');

    
    % Alkalinity
    nccreate(fname,'ALK_time', ...
            'Dimensions',{'ALK_time' ALK_time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ALK_time','cycle_length',365.25);
    ncwriteatt(fname,'ALK_time','long_name','Open boundary conditions time');
    ncwriteatt(fname,'ALK_time','units','days since 1900-01-01 00:00:00');
    ncwriteatt(fname,'ALK_time','field','ALK_time, scalar, series');

    nccreate(fname,'ALK_west', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'ALK_time' ALK_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ALK_west','long_name','ALK western boundary condition');
    ncwriteatt(fname,'ALK_west','units','milliequivalents meter-3');
    ncwriteatt(fname,'ALK_west','field','ALK_west, scalar, series');
    ncwriteatt(fname,'ALK_west','time','ALK_time');

    nccreate(fname,'ALK_east', ...
            'Dimensions',{'eta_rho' eta_rho 's_rho' s_rho 'ALK_time' ALK_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ALK_east','long_name','ALK eastern boundary condition');
    ncwriteatt(fname,'ALK_east','units','milliequivalents meter-3');
    ncwriteatt(fname,'ALK_east','field','ALK_east, scalar, series');
    ncwriteatt(fname,'ALK_east','time','ALK_time');

    nccreate(fname,'ALK_north', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'ALK_time' ALK_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ALK_north','long_name','ALK northern boundary condition');
    ncwriteatt(fname,'ALK_north','units','milliequivalents meter-3');
    ncwriteatt(fname,'ALK_north','field','ALK_north, scalar, series');
    ncwriteatt(fname,'ALK_north','time','ALK_time');

    nccreate(fname,'ALK_south', ...
            'Dimensions',{'xi_rho' xi_rho 's_rho' s_rho 'ALK_time' ALK_time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ALK_south','long_name','ALK southern boundary condition');
    ncwriteatt(fname,'ALK_south','units','milliequivalents meter-3');
    ncwriteatt(fname,'ALK_south','field','ALK_south, scalar, series');
    ncwriteatt(fname,'ALK_south','time','ALK_time');         
    
    %
    %----------------------------------------------------------------------
    % Read data into variables.  Where data are unavailable, populate with 
    %   analytical functions.
    %----------------------------------------------------------------------    
    %

    % Time variables
    bry_time   = ncread(fid,'bry_time');
    NO3_time   = bry_time;
    NH4_time   = bry_time;
    phy_time   = bry_time;
    dia_time   = bry_time;
    bac_time   = bry_time;
    chla_time  = bry_time;
    chld_time  = bry_time;
    DON_time   = bry_time;
    DOC_time   = bry_time;
    DetN_time  = bry_time;
    DetC_time  = bry_time;
    sil_time   = bry_time;
    LDetN_time = bry_time;
    LDetC_time = bry_time;
    LDetS_time = bry_time;
    micro_time = bry_time;
    meso_time  = bry_time;
    FeD_time   = bry_time;
    DetFe_time   = bry_time;
    LDetFe_time = bry_time;
    coccos_time = bry_time;
    ChlC_time = bry_time;
    PIC_time = bry_time;
    DIC_time = bry_time;
    ALK_time = bry_time;
    

    % zeta
    zeta_west  = ncread(fid,'zeta_west');
    zeta_east  = ncread(fid,'zeta_east');
    zeta_north = ncread(fid,'zeta_north');
    zeta_south = ncread(fid,'zeta_south');

    % ubar
    ubar_west  = ncread(fid,'ubar_west');
    ubar_east  = ncread(fid,'ubar_east');
    ubar_north = ncread(fid,'ubar_north');
    ubar_south = ncread(fid,'ubar_south');

    % vbar
    vbar_west  = ncread(fid,'vbar_west');
    vbar_east  = ncread(fid,'vbar_east');
    vbar_north = ncread(fid,'vbar_north');
    vbar_south = ncread(fid,'vbar_south');

    % u
    u_west  = ncread(fid,'u_west');
    u_east  = ncread(fid,'u_east');
    u_north = ncread(fid,'u_north');
    u_south = ncread(fid,'u_south');

    % v
    v_west  = ncread(fid,'v_west');
    v_east  = ncread(fid,'v_east');
    v_north = ncread(fid,'v_north');
    v_south = ncread(fid,'v_south');

    % temp
    temp_west  = ncread(fid,'temp_west');
    temp_east  = ncread(fid,'temp_east');
    temp_north = ncread(fid,'temp_north');
    temp_south = ncread(fid,'temp_south');

    % salt
    salt_west  = ncread(fid,'salt_west');
    salt_east  = ncread(fid,'salt_east');
    salt_north = ncread(fid,'salt_north');
    salt_south = ncread(fid,'salt_south');

    % NO3
    NO3_west  = ncread(fid,'NO3_west');
    NO3_east  = ncread(fid,'NO3_east');
    NO3_north = ncread(fid,'NO3_north');
    NO3_south = ncread(fid,'NO3_south');

    % NH4
%    NH4_west  = ncread(fid,'NH4_west');
%    NH4_east  = ncread(fid,'NH4_east');
%    NH4_north = ncread(fid,'NH4_north');
%    NH4_south = ncread(fid,'NH4_south');
    NH4_west  = 0.01*NO3_west;
    NH4_east  = 0.01*NO3_east;
    NH4_north = 0.01*NO3_north;
    NH4_south = 0.01*NO3_south;

    % Phytoplankton
    Phytoplankton_west  = ncread(fid,'SPHY_west');
    Phytoplankton_east  = ncread(fid,'SPHY_east');
    Phytoplankton_north = ncread(fid,'SPHY_north');
    Phytoplankton_south = ncread(fid,'SPHY_south');
%    Phytoplankton_west  = Phytoplankton_west*100;
%    Phytoplankton_east  = Phytoplankton_east*100;
%    Phytoplankton_north = Phytoplankton_north*100;
%    Phytoplankton_south = Phytoplankton_south*100;

    % Diatoms
    Diatoms_west  = ncread(fid,'LPHY_west');
    Diatoms_east  = ncread(fid,'LPHY_east');
    Diatoms_north = ncread(fid,'LPHY_north');
    Diatoms_south = ncread(fid,'LPHY_south');
 %   Diatoms_west  = Diatoms_west*100;
%    Diatoms_east  = Diatoms_east*100;
%    Diatoms_north = Diatoms_north*100;
%    Diatoms_south = Diatoms_south*100;


    % Bacteria
    Bacteria_west  = 0.1*Phytoplankton_west;
    Bacteria_east  = 0.1*Phytoplankton_east;
    Bacteria_north = 0.1*Phytoplankton_north;
    Bacteria_south = 0.1*Phytoplankton_south;

    % Chlorophyll
%    Chlorophyll_west  = ncread(fid,'CHLAS_west');
%    Chlorophyll_east  = ncread(fid,'CHLAS_east');
%    Chlorophyll_north = ncread(fid,'CHLAS_north');
%    Chlorophyll_south = ncread(fid,'CHLAS_south');
   Chlorophyll_west  = Phytoplankton_west*1.59;
   Chlorophyll_east  = Phytoplankton_east*1.59;
   Chlorophyll_north = Phytoplankton_north*1.59;
   Chlorophyll_south = Phytoplankton_south*1.59;


    % ChlD
%    ChlD_west  = ncread(fid,'CHLAL_west');
%    ChlD_east  = ncread(fid,'CHLAL_east');
%    ChlD_north = ncread(fid,'CHLAL_north');
%    ChlD_south = ncread(fid,'CHLAL_south');
   ChlD_west  = Diatoms_west*1.59;
   ChlD_east  = Diatoms_east*1.59;
   ChlD_north = Diatoms_north*1.59;
   ChlD_south = Diatoms_south*1.59;

    % DON
%    DON_west  = ncread(fid,'DON_west');
%    DON_east  = ncread(fid,'DON_east');
%    DON_north = ncread(fid,'DON_north');
%    DON_south = ncread(fid,'DON_south');
   DON_west  = NH4_west;
   DON_east  = NH4_east;
   DON_north = NH4_north;
   DON_south = NH4_south;

    % DOC
    DOC_west  = 6.625*DON_west;
    DOC_east  = 6.625*DON_east;
    DOC_north = 6.625*DON_north;
    DOC_south = 6.625*DON_south;

    % DetN
%    DetN_west  = ncread(fid,'PON_west');
%    DetN_east  = ncread(fid,'PON_east');
%    DetN_north = ncread(fid,'PON_north');
%    DetN_south = ncread(fid,'PON_south');
    DetN_west  = 0.1*DON_west;
    DetN_east  = 0.1*DON_east;
    DetN_north = 0.1*DON_north;
    DetN_south = 0.1*DON_south;

    % DetC
    DetC_west  = 6.625*DetN_west;
    DetC_east  = 6.625*DetN_east;
    DetC_north = 6.625*DetN_north;
    DetC_south = 6.625*DetN_south;

    % Silica
    Silica_west  = ncread(fid,'SiOH_west');
    Silica_east  = ncread(fid,'SiOH_east');
    Silica_north = ncread(fid,'SiOH_north');
    Silica_south = ncread(fid,'SiOH_south');

    % LDetN
    LDetN_west  = 0.1*DetN_west;
    LDetN_east  = 0.1*DetN_east;
    LDetN_north = 0.1*DetN_north;
    LDetN_south = 0.1*DetN_south;

    % LDetC
    LDetC_west  = 6.625*LDetN_west;
    LDetC_east  = 6.625*LDetN_east;
    LDetC_north = 6.625*LDetN_north;
    LDetC_south = 6.625*LDetN_south;

    % LDetS
%    LDetS_west  = ncread(fid,'OPAL_west');
%    LDetS_east  = ncread(fid,'OPAL_east');
%    LDetS_north = ncread(fid,'OPAL_north');
%    LDetS_south = ncread(fid,'OPAL_south');
    LDetS_west  = 0.01*Silica_west;
    LDetS_east  = 0.01*Silica_east;
    LDetS_north = 0.01*Silica_north;
    LDetS_south = 0.01*Silica_south;

    % Microzooplankton
%    Microzooplankton_west  = ncread(fid,'SZOO_west');
%    Microzooplankton_east  = ncread(fid,'SZOO_east');
%    Microzooplankton_north = ncread(fid,'SZOO_north');
%    Microzooplankton_south = ncread(fid,'SZOO_south');
   Microzooplankton_west  = Phytoplankton_west*0.01;
   Microzooplankton_east  = Phytoplankton_east*0.01;
   Microzooplankton_north = Phytoplankton_north*0.01;
   Microzooplankton_south = Phytoplankton_south*0.01;

    % Mesozooplankton
    Mesozooplankton_west  = Diatoms_west*0.01;
    Mesozooplankton_east  = Diatoms_east*0.01;
    Mesozooplankton_north = Diatoms_north*0.01;
    Mesozooplankton_south = Diatoms_south*0.01;

    % Dissolved iron
    FeD_west  = ncread(fid,'FeD_west');
    FeD_east  = ncread(fid,'FeD_east');
    FeD_north = ncread(fid,'FeD_north');
    FeD_south = ncread(fid,'FeD_south');

    % Particulate iron
%    FeSp_west  = ncread(fid,'FeSp_west');
%    FeLp_west  = ncread(fid,'FeLp_west');
%    FeSp_east  = ncread(fid,'FeSp_east');
%    FeLp_east  = ncread(fid,'FeLp_east');
%    FeSp_north = ncread(fid,'FeSp_north');
%    FeLp_north = ncread(fid,'FeLp_north');
%    FeSp_south = ncread(fid,'FeSp_south');
%    FeLp_south = ncread(fid,'FeLp_south');

%    DetFe_west  = FeSp_west  + FeLp_west;
%    DetFe_east  = FeSp_east  + FeLp_east;
%    DetFe_north = FeSp_north + FeLp_north;
%    DetFe_south = FeSp_south + FeLp_north;
    DetFe_west  = FeD_west*0.1;
    DetFe_east  = FeD_east*0.1;
    DetFe_north = FeD_north*0.1;
    DetFe_south = FeD_south*0.1;
    
    % LDetFe
    LDetFe_west  = FeD_west*0.1;
    LDetFe_east  = FeD_east*0.1;
    LDetFe_north = FeD_north*0.1;
    LDetFe_south = FeD_south*0.1;
    
    % set slices to read boundary data from
    west_start = [1,1,1,2];
    west_end = [1,Inf,Inf,Inf];
    east_start = [xi_rho,1,1,2];
    east_end = [Inf,Inf,Inf,Inf];
    north_start = [1,eta_rho,1,2];
    north_end = [Inf,Inf,Inf,Inf];
    south_start = [1,1,1,2];
    south_end = [Inf,1,Inf,Inf];
    
    % Read variables and duplicate so there are 12 copies of data in time
    % dimension.  For now, weird logic replicates the arrays to 16, then
    % pairs down to 12.
    % Coccolithophores
    Coccos_west = squeeze(ncread(fid_cocco,'Coccos',west_start,west_end));
    for k=2:5
        Coccos_west = cat(3,Coccos_west,Coccos_west);
    end
    Coccos_west = Coccos_west(:,:,1:12);
    Coccos_east = squeeze(ncread(fid_cocco,'Coccos',east_start,east_end));
    for k=2:5
        Coccos_east = cat(3,Coccos_east,Coccos_east);
    end
    Coccos_east = Coccos_east(:,:,1:12);
    Coccos_north = squeeze(ncread(fid_cocco,'Coccos',north_start,north_end));
    for k=2:5
        Coccos_north = cat(3,Coccos_north,Coccos_north);
    end
    Coccos_north = Coccos_north(:,:,1:12);
    Coccos_south = squeeze(ncread(fid_cocco,'Coccos',south_start,south_end));
    for k=2:5
        Coccos_south = cat(3,Coccos_south,Coccos_south);
    end    
    Coccos_south = Coccos_south(:,:,1:12);
    
    % Coccolithophore Chlorophyll
    ChlC_west = squeeze(ncread(fid_cocco,'ChlC',west_start,west_end));
    for k=2:5
        ChlC_west = cat(3,ChlC_west,ChlC_west);
    end
    ChlC_west = ChlC_west(:,:,1:12);
    ChlC_east = squeeze(ncread(fid_cocco,'ChlC',east_start,east_end));
    for k=2:5
        ChlC_east = cat(3,ChlC_east,ChlC_east);
    end
    ChlC_east = ChlC_east(:,:,1:12);
    ChlC_north = squeeze(ncread(fid_cocco,'ChlC',north_start,north_end));
    for k=2:5
        ChlC_north = cat(3,ChlC_north,ChlC_north);
    end
    ChlC_north = ChlC_north(:,:,1:12);
    ChlC_south = squeeze(ncread(fid_cocco,'ChlC',south_start,south_end));
    for k=2:5
        ChlC_south = cat(3,ChlC_south,ChlC_south);
    end    
    ChlC_south = ChlC_south(:,:,1:12);  
    
    % Particulate Inorganic Carbon
    PIC_west = squeeze(ncread(fid_cocco,'PIC',west_start,west_end));
    for k=2:5
        PIC_west = cat(3,PIC_west,PIC_west);
    end
    PIC_west = PIC_west(:,:,1:12);
    PIC_east = squeeze(ncread(fid_cocco,'PIC',east_start,east_end));
    for k=2:5
        PIC_east = cat(3,PIC_east,PIC_east);
    end
    PIC_east = PIC_east(:,:,1:12);
    PIC_north = squeeze(ncread(fid_cocco,'PIC',north_start,north_end));
    for k=2:5
        PIC_north = cat(3,PIC_north,PIC_north);
    end
    PIC_north = PIC_north(:,:,1:12);
    PIC_south = squeeze(ncread(fid_cocco,'PIC',south_start,south_end));
    for k=2:5
        PIC_south = cat(3,PIC_south,PIC_south);
    end    
    PIC_south = PIC_south(:,:,1:12);
    
    % Dissolved Inorganic Carbon
    DIC_west = squeeze(ncread(fid_cocco,'DIC',west_start,west_end));
    for k=2:5
        DIC_west = cat(3,DIC_west,DIC_west);
    end
    DIC_west = DIC_west(:,:,1:12);
    DIC_east = squeeze(ncread(fid_cocco,'DIC',east_start,east_end));
    for k=2:5
        DIC_east = cat(3,DIC_east,DIC_east);
    end
    DIC_east = DIC_east(:,:,1:12);
    DIC_north = squeeze(ncread(fid_cocco,'DIC',north_start,north_end));
    for k=2:5
        DIC_north = cat(3,DIC_north,DIC_north);
    end
    DIC_north = DIC_north(:,:,1:12);
    DIC_south = squeeze(ncread(fid_cocco,'DIC',south_start,south_end));
    for k=2:5
        DIC_south = cat(3,DIC_south,DIC_south);
    end    
    DIC_south = DIC_south(:,:,1:12); 
    
    % Alkalinity
    ALK_west = squeeze(ncread(fid_cocco,'ALK',west_start,west_end));
    for k=2:5
        ALK_west = cat(3,ALK_west,ALK_west);
    end
    ALK_west = ALK_west(:,:,1:12);
    ALK_east = squeeze(ncread(fid_cocco,'ALK',east_start,east_end));
    for k=2:5
        ALK_east = cat(3,ALK_east,ALK_east);
    end
    ALK_east = ALK_east(:,:,1:12);
    ALK_north = squeeze(ncread(fid_cocco,'ALK',north_start,north_end));
    for k=2:5
        ALK_north = cat(3,ALK_north,ALK_north);
    end
    ALK_north = ALK_north(:,:,1:12);
    ALK_south = squeeze(ncread(fid_cocco,'ALK',south_start,south_end));
    for k=2:5
        ALK_south = cat(3,ALK_south,ALK_south);
    end    
    ALK_south = ALK_south(:,:,1:12);
    
    %
    %----------------------------------------------------------------------
    % Write data to variables.
    %----------------------------------------------------------------------    
    %
    % bry_time
    ncwrite(fname,'bry_time',bry_time);

    % zeta
    ncwrite(fname,'zeta_west',zeta_west);
    ncwrite(fname,'zeta_east',zeta_east);
    ncwrite(fname,'zeta_north',zeta_north);
    ncwrite(fname,'zeta_south',zeta_south);

    % ubar
    ncwrite(fname,'ubar_west',ubar_west);
    ncwrite(fname,'ubar_east',ubar_east);
    ncwrite(fname,'ubar_north',ubar_north);
    ncwrite(fname,'ubar_south',ubar_south);

    % vbar
    ncwrite(fname,'vbar_west',vbar_west);
    ncwrite(fname,'vbar_east',vbar_east);
    ncwrite(fname,'vbar_north',vbar_north);
    ncwrite(fname,'vbar_south',vbar_south);

    % u
    ncwrite(fname,'u_west',u_west);
    ncwrite(fname,'u_east',u_east);
    ncwrite(fname,'u_north',u_north);
    ncwrite(fname,'u_south',u_south);

    % v
    ncwrite(fname,'v_west',v_west);
    ncwrite(fname,'v_east',v_east);
    ncwrite(fname,'v_north',v_north);
    ncwrite(fname,'v_south',v_south);

    % temp
    ncwrite(fname,'temp_west',temp_west);
    ncwrite(fname,'temp_east',temp_east);
    ncwrite(fname,'temp_north',temp_north);
    ncwrite(fname,'temp_south',temp_south);

    % salt
    ncwrite(fname,'salt_west',salt_west);
    ncwrite(fname,'salt_east',salt_east);
    ncwrite(fname,'salt_north',salt_north);
    ncwrite(fname,'salt_south',salt_south);

    % NO3
    ncwrite(fname,'NO3_time',NO3_time);
    ncwrite(fname,'NO3_west',NO3_west);
    ncwrite(fname,'NO3_east',NO3_east);
    ncwrite(fname,'NO3_north',NO3_north);
    ncwrite(fname,'NO3_south',NO3_south);

    % NH4
    ncwrite(fname,'NH4_time',NH4_time);
    ncwrite(fname,'NH4_west',NH4_west);
    ncwrite(fname,'NH4_east',NH4_east);
    ncwrite(fname,'NH4_north',NH4_north);
    ncwrite(fname,'NH4_south',NH4_south);

    % Phytoplankton
    ncwrite(fname,'phy_time',phy_time);
    ncwrite(fname,'Phytoplankton_west',Phytoplankton_west);
    ncwrite(fname,'Phytoplankton_east',Phytoplankton_east);
    ncwrite(fname,'Phytoplankton_north',Phytoplankton_north);
    ncwrite(fname,'Phytoplankton_south',Phytoplankton_south);

    % Diatoms
    ncwrite(fname,'dia_time',dia_time);
    ncwrite(fname,'Diatoms_west',Diatoms_west);
    ncwrite(fname,'Diatoms_east',Diatoms_east);
    ncwrite(fname,'Diatoms_north',Diatoms_north);
    ncwrite(fname,'Diatoms_south',Diatoms_south);

    % Bacteria
    ncwrite(fname,'bac_time',bac_time);
    ncwrite(fname,'Bacteria_west',Bacteria_west);
    ncwrite(fname,'Bacteria_east',Bacteria_east);
    ncwrite(fname,'Bacteria_north',Bacteria_north);
    ncwrite(fname,'Bacteria_south',Bacteria_south);

    % Chlorophyll
    ncwrite(fname,'chla_time',chla_time);
    ncwrite(fname,'Chlorophyll_west',Chlorophyll_west);
    ncwrite(fname,'Chlorophyll_east',Chlorophyll_east);
    ncwrite(fname,'Chlorophyll_north',Chlorophyll_north);
    ncwrite(fname,'Chlorophyll_south',Chlorophyll_south);

    % ChlD
    ncwrite(fname,'chld_time',chld_time);
    ncwrite(fname,'ChlD_west',ChlD_west);
    ncwrite(fname,'ChlD_east',ChlD_east);
    ncwrite(fname,'ChlD_north',ChlD_north);
    ncwrite(fname,'ChlD_south',ChlD_south);

    % DON
    ncwrite(fname,'DON_time',DON_time);
    ncwrite(fname,'DON_west',DON_west);
    ncwrite(fname,'DON_east',DON_east);
    ncwrite(fname,'DON_north',DON_north);
    ncwrite(fname,'DON_south',DON_south);

    % DOC
    ncwrite(fname,'DOC_time',DOC_time);
    ncwrite(fname,'DOC_west',DOC_west);
    ncwrite(fname,'DOC_east',DOC_east);
    ncwrite(fname,'DOC_north',DOC_north);
    ncwrite(fname,'DOC_south',DOC_south);

    % DetN
    ncwrite(fname,'DetN_time',DetN_time);
    ncwrite(fname,'DetN_west',DetN_west);
    ncwrite(fname,'DetN_east',DetN_east);
    ncwrite(fname,'DetN_north',DetN_north);
    ncwrite(fname,'DetN_south',DetN_south);

    % DetC
    ncwrite(fname,'DetC_time',DetC_time);
    ncwrite(fname,'DetC_west',DetC_west);
    ncwrite(fname,'DetC_east',DetC_east);
    ncwrite(fname,'DetC_north',DetC_north);
    ncwrite(fname,'DetC_south',DetC_south);

    % Silica
    ncwrite(fname,'sil_time',sil_time);
    ncwrite(fname,'Silica_west',Silica_west);
    ncwrite(fname,'Silica_east',Silica_east);
    ncwrite(fname,'Silica_north',Silica_north);
    ncwrite(fname,'Silica_south',Silica_south);

    % LDetN
    ncwrite(fname,'LDetN_time',LDetN_time);
    ncwrite(fname,'LDetN_west',LDetN_west);
    ncwrite(fname,'LDetN_east',LDetN_east);
    ncwrite(fname,'LDetN_north',LDetN_north);
    ncwrite(fname,'LDetN_south',LDetN_south);

    % LDetC
    ncwrite(fname,'LDetC_time',LDetC_time);
    ncwrite(fname,'LDetC_west',LDetC_west);
    ncwrite(fname,'LDetC_east',LDetC_east);
    ncwrite(fname,'LDetC_north',LDetC_north);
    ncwrite(fname,'LDetC_south',LDetC_south);

    % LDetS
    ncwrite(fname,'LDetS_time',LDetS_time);
    ncwrite(fname,'LDetS_west',LDetS_west);
    ncwrite(fname,'LDetS_east',LDetS_east);
    ncwrite(fname,'LDetS_north',LDetS_north);
    ncwrite(fname,'LDetS_south',LDetS_south);

    % Microzooplankton
    ncwrite(fname,'micro_time',micro_time);
    ncwrite(fname,'Microzooplankton_west',Microzooplankton_west);
    ncwrite(fname,'Microzooplankton_east',Microzooplankton_east);
    ncwrite(fname,'Microzooplankton_north',Microzooplankton_north);
    ncwrite(fname,'Microzooplankton_south',Microzooplankton_south);

    % Mesozooplankton
    ncwrite(fname,'meso_time',meso_time);
    ncwrite(fname,'Mesozooplankton_west',Mesozooplankton_west);
    ncwrite(fname,'Mesozooplankton_east',Mesozooplankton_east);
    ncwrite(fname,'Mesozooplankton_north',Mesozooplankton_north);
    ncwrite(fname,'Mesozooplankton_south',Mesozooplankton_south);

    % Dissolved iron
    ncwrite(fname,'FeD_time',FeD_time);
    ncwrite(fname,'FeD_west',FeD_west);
    ncwrite(fname,'FeD_east',FeD_east);
    ncwrite(fname,'FeD_north',FeD_north);
    ncwrite(fname,'FeD_south',FeD_south);

   % Detrital iron
   ncwrite(fname,'DetFe_time',DetFe_time);
   ncwrite(fname,'DetFe_west',DetFe_west);
   ncwrite(fname,'DetFe_east',DetFe_east);
   ncwrite(fname,'DetFe_north',DetFe_north);
   ncwrite(fname,'DetFe_south',DetFe_south);

   % Detrital iron
   ncwrite(fname,'LDetFe_time',LDetFe_time);
   ncwrite(fname,'LDetFe_west',LDetFe_west);
   ncwrite(fname,'LDetFe_east',LDetFe_east);
   ncwrite(fname,'LDetFe_north',LDetFe_north);
   ncwrite(fname,'LDetFe_south',LDetFe_south);
   
    % Coccolithophores
   ncwrite(fname,'coccos_time',coccos_time);
   ncwrite(fname,'Coccos_west',Coccos_west);
   ncwrite(fname,'Coccos_east',Coccos_east);
   ncwrite(fname,'Coccos_north',Coccos_north);
   ncwrite(fname,'Coccos_south',Coccos_south);    
    
    % Coccolithophore Chlorophyll
   ncwrite(fname,'ChlC_time',ChlC_time);
   ncwrite(fname,'ChlC_west',ChlC_west);
   ncwrite(fname,'ChlC_east',ChlC_east);
   ncwrite(fname,'ChlC_north',ChlC_north);
   ncwrite(fname,'ChlC_south',ChlC_south);        
    % Particulate Inorganic Carbon
   ncwrite(fname,'PIC_time',PIC_time);
   ncwrite(fname,'PIC_west',PIC_west);
   ncwrite(fname,'PIC_east',PIC_east);
   ncwrite(fname,'PIC_north',PIC_north);
   ncwrite(fname,'PIC_south',PIC_south);        
   
    % Dissolved Inorganic Carbon
   ncwrite(fname,'DIC_time',DIC_time);
   ncwrite(fname,'DIC_west',DIC_west);
   ncwrite(fname,'DIC_east',DIC_east);
   ncwrite(fname,'DIC_north',DIC_north);
   ncwrite(fname,'DIC_south',DIC_south);    
   
    % Alkalinity   
   ncwrite(fname,'ALK_time',ALK_time);
   ncwrite(fname,'ALK_west',ALK_west);
   ncwrite(fname,'ALK_east',ALK_east);
   ncwrite(fname,'ALK_north',ALK_north);
   ncwrite(fname,'ALK_south',ALK_south);        

end  % create_bry
