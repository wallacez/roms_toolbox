% setup_ini.nc
%
% Create an initial conditions file with appropriate metadata.
% The output NetCDF file will properly setup the variables, dimensions, and
% attributes for nc_inq to query and feed into the d_inital_zw script for
% proper stretching of inital conditions.
%
% Author: Z. Wallace
% Created 13 August 2018
% Rev. hist.: 14 August 2018. --> Write data to variables

function setup_ini(ncname)
    
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
    s_w     = 41;
    time    = 0;

    ncid = netcdf.create(ncname,'NETCDF4');
    
    % create dimensions
    d1 = netcdf.defDim(ncid,'xi_rho' ,xi_rho );
    d2 = netcdf.defDim(ncid,'eta_rho',eta_rho);
    d3 = netcdf.defDim(ncid,'xi_u'   ,xi_u   );
    d4 = netcdf.defDim(ncid,'eta_v'  ,eta_v  );
    d5 = netcdf.defDim(ncid,'xi_psi' ,xi_psi );
    d6 = netcdf.defDim(ncid,'eta_psi',eta_psi);
    d7 = netcdf.defDim(ncid,'s_rho'  ,s_rho  );
    d8 = netcdf.defDim(ncid,'s_w'    ,s_w    );
    d9 = netcdf.defDim(ncid,'time'   ,time   );
    
    % define variables
    varid0  = netcdf.defVar(ncid,'ocean_time' ,'double',d9);
    varid1  = netcdf.defVar(ncid,'spherical'  ,'char'  ,[]);
    varid2  = netcdf.defVar(ncid,'Vtransform' ,'NC_INT',[]);
    varid3  = netcdf.defVar(ncid,'Vstretching','NC_INT',[]);
    varid4  = netcdf.defVar(ncid,'theta_s'    ,'double',[]);
    varid5  = netcdf.defVar(ncid,'theta_b'    ,'double',[]);
    varid6  = netcdf.defVar(ncid,'Tcline'     ,'double',[]);
    varid7  = netcdf.defVar(ncid,'hc'         ,'double',[]);
    varid8  = netcdf.defVar(ncid,'s_rho'      ,'double',d7);
    varid9  = netcdf.defVar(ncid,'s_w'        ,'double',d8);
    varid10 = netcdf.defVar(ncid,'Cs_r'       ,'double',d7);
    varid11 = netcdf.defVar(ncid,'Cs_w'       ,'double',d8);
    varid12 = netcdf.defVar(ncid,'h'          ,'double',[d1 d2]);
    varid13 = netcdf.defVar(ncid,'lon_rho'    ,'double',[d1 d2]);
    varid14 = netcdf.defVar(ncid,'lat_rho'    ,'double',[d1 d2]);
    varid15 = netcdf.defVar(ncid,'zeta'       ,'double',[d1 d2 d9]);
    varid16 = netcdf.defVar(ncid,'ubar'       ,'double',[d1 d2 d9]);
    varid17 = netcdf.defVar(ncid,'vbar'       ,'double',[d1 d2 d9]);
    varid18 = netcdf.defVar(ncid,'u'          ,'double',[d1 d2 d7 d9]);
    varid19 = netcdf.defVar(ncid,'v'          ,'double',[d1 d2 d7 d9]);
    varid20 = netcdf.defVar(ncid,'temp'       ,'double',[d1 d2 d7 d9]);
    varid21 = netcdf.defVar(ncid,'salt'       ,'double',[d1 d2 d7 d9]);
    
    % define variable attributes
    % ocean_time
    netcdf.putAtt(ncid,varid0,'long_name','time since initialization')
    netcdf.putAtt(ncid,varid0,'units','seconds since 1900-01-01 00:00:00')
    
    % spherical
    netcdf.putAtt(ncid,varid1,'long_name','grid type logical switch')
    netcdf.putAtt(ncid,varid1,'flag_values','T, F')
    netcdf.putAtt(ncid,varid1,'flag_meanings','spherical, Cartesian')
    
    % Vtransform
    netcdf.putAtt(ncid,varid2,'long_name','vertical s-coord. transformation equation')
    
    % Vtransform
    netcdf.putAtt(ncid,varid3,'long_name','vertical s-coord. stretching function')
    
    % theta_s
    netcdf.putAtt(ncid,varid4,'long_name','S-coordinate surface control parameter')
    
    % theta_b
    netcdf.putAtt(ncid,varid5,'long_name','S-coordinate bottom control parameter')
    
    % Tcline
    netcdf.putAtt(ncid,varid6,'long_name','S-coordinate surface/bottom layer width')
    
    % hc
    netcdf.putAtt(ncid,varid7,'long_name','S-coordinate parameter, critical depth')
    
    % s_rho
    netcdf.putAtt(ncid,varid8,'long_name','S-coordinate at rho-points')
    netcdf.putAtt(ncid,varid8,'positive' ,'up')
    
    % s_w
    netcdf.putAtt(ncid,varid9,'long_name','S-coordinate at W-points')
    netcdf.putAtt(ncid,varid9,'positive' ,'up')
    
    % Cs_r
    netcdf.putAtt(ncid,varid10,'long_name','S-coordinate stretching curves at rho-points')
    netcdf.putAtt(ncid,varid10,'units'    ,'n.d.')
    
    % Cs_w
    netcdf.putAtt(ncid,varid11,'long_name','S-coordinate stretching curves at W-points')
    netcdf.putAtt(ncid,varid11,'units'    ,'n.d.')
    
    % h
    netcdf.putAtt(ncid,varid12,'long_name','bathymetry at rho-points')
    netcdf.putAtt(ncid,varid12,'units'    ,'meter')
    
    % lon_rho
    netcdf.putAtt(ncid,varid13,'long_name','longitude at rho-points')
    netcdf.putAtt(ncid,varid13,'units'    ,'degree_east')
    
    % lat_rho
    netcdf.putAtt(ncid,varid14,'long_name','latitude at rho-points')
    netcdf.putAtt(ncid,varid14,'units'    ,'degree_north')
    
    % zeta
    netcdf.putAtt(ncid,varid15,'long_name','free-surface')
    netcdf.putAtt(ncid,varid15,'units'    ,'meter')
    netcdf.putAtt(ncid,varid15,'time'     ,'ocean_time')
    
    % ubar
    netcdf.putAtt(ncid,varid16,'long_name','vertically integrated u-momentum component')
    netcdf.putAtt(ncid,varid16,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,varid16,'time'     ,'ocean_time')
    
    % vbar
    netcdf.putAtt(ncid,varid17,'long_name','vertically integrated v-momentum component')
    netcdf.putAtt(ncid,varid17,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,varid17,'time'     ,'ocean_time')
    
    % u
    netcdf.putAtt(ncid,varid18,'long_name','u-momentum component')
    netcdf.putAtt(ncid,varid18,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,varid18,'time'     ,'ocean_time')
    
    % v
    netcdf.putAtt(ncid,varid19,'long_name','v-momentum component')
    netcdf.putAtt(ncid,varid19,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,varid19,'time'     ,'ocean_time')
    
    % temperature
    netcdf.putAtt(ncid,varid20,'long_name','potential temperature')
    netcdf.putAtt(ncid,varid20,'units'    ,'Celsius')
    netcdf.putAtt(ncid,varid20,'time'     ,'ocean_time')
    
    % salinity
    netcdf.putAtt(ncid,varid21,'long_name','salinity')
    netcdf.putAtt(ncid,varid21,'time'     ,'ocean_time')
    
    % define global attributes
    varid = netcdf.getConstant('GLOBAL');
    
    netcdf.putAtt(ncid,varid,'type'    ,'ROMS initial file')
    netcdf.putAtt(ncid,varid,'title'   ,'Patagonia model initial conditions')
    netcdf.putAtt(ncid,varid,'grd_file','roms_grd_rivers.nc')
    netcdf.putAtt(ncid,varid,'date'    , datestr(today))
    
    % obtain data for initial conditions files.
    % set files from which to obtain data.
    iname = 'roms_ini_1900.nc'  % file from which to read initial data
    gname = 'roms_grd_rivers.nc' % file from which to read grid data
    
    % set scalar variables
    ot = (datenum('24-Jun-2001')-datenum('01-Jan-1900'))*86400;
    sp = 'T';
    Vt = 2;  % transformation equation
    Vs = 4;  % stretching equation
    ts = 6.;
    tb = 0.;
    tc = 10;
    hc = tc;
    N  = 40;  % for stretching function
    
    % set vertical s-coordinate parameters.
    [s_r, C_r] = stretching(Vs, ts, tb, hc, N, 0, 1);
    [s_w, C_w] = stretching(Vs, ts, tb, hc, N, 1, 1);
    
    % grid data
    h    = ncread(gname, 'h');  % bathymetry
    lo_r = ncread(gname, 'lon_rho'); 
    la_r = ncread(gname, 'lat_rho');
    
    % 2-D variables
    zeta = ncread(iname, 'zeta');
    ubar = ncread(iname, 'ubar');
    vbar = ncread(iname, 'vbar');
    
    % 3-D variables
    u    = ncread(iname, 'u');
    v    = ncread(iname, 'v');
    temp = ncread(iname, 'temp');
    salt = ncread(iname, 'salt'); 
    
    % setup finished, close file.
    netcdf.close(ncid)  

    % write data to variables
    ncwrite(ncname, 'ocean_time' , ot  )
    ncwrite(ncname, 'spherical'  , sp  )
    ncwrite(ncname, 'Vtransform' , Vt  )
    ncwrite(ncname, 'Vstretching', Vs  )
    ncwrite(ncname, 'theta_s'    , ts  )
    ncwrite(ncname, 'theta_b'    , tb  )
    ncwrite(ncname, 'Tcline'     , tc  )
    ncwrite(ncname, 'hc'         , hc  )
    ncwrite(ncname, 's_rho'      , s_r )
    ncwrite(ncname, 'Cs_r'       , C_r )
    ncwrite(ncname, 's_w'        , s_w )
    ncwrite(ncname, 'Cs_w'       , C_w )
    ncwrite(ncname, 'h'          , h   )
    ncwrite(ncname, 'lon_rho'    , lo_r)
    ncwrite(ncname, 'lat_rho'    , la_r)
    ncwrite(ncname, 'zeta'       , zeta)
    ncwrite(ncname, 'ubar'       , ubar)
    ncwrite(ncname, 'vbar'       , vbar)
    ncwrite(ncname, 'u'          , u   )
    ncwrite(ncname, 'v'          , v   )
    ncwrite(ncname, 'temp'       , temp)
    ncwrite(ncname, 'salt'       , salt)
end

