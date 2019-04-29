% setup_bry.nc
%
% Create a boundary conditions file with appropriate metadata.
% The output NetCDF file will properly setup the variables, dimensions, and
% attributes for the proper stretching function specified by user input.
%
% Author: Z. Wallace | 20 August 2018
% Rev. his.: 8.22.18 --> updated calendar attribute to reflect UNIDATA
%                        standard "proleptic_gregorian"

function setup_bry(ncname)
    
    % overwrite netcdf file, if it exists
    if exist(ncname,'file') ~= 0
        delete(ncname)
    end

    % define parameters to be used as dimensions
    xi_rho    = 482;
    eta_rho   = 770;
    xi_u      = 481;
    eta_v     = 769;
    xi_psi    = 481;
    eta_psi   = 769;
    s_rho     = 40;
    s_w       = 41;
    zeta_time = 12;
    v2d_time  = 12;
    v3d_time  = 12;
    temp_time = 12;
    salt_time = 12;

    cyc_len = 365.25;  % Length of year

    ncid = netcdf.create(ncname,'NETCDF4');
    
    % create dimensions
    d1  = netcdf.defDim(ncid,'xi_rho'   ,xi_rho   );
    d2  = netcdf.defDim(ncid,'eta_rho'  ,eta_rho  );
    d3  = netcdf.defDim(ncid,'xi_u'     ,xi_u     );
    d4  = netcdf.defDim(ncid,'eta_v'    ,eta_v    );
    d5  = netcdf.defDim(ncid,'xi_psi'   ,xi_psi   );
    d6  = netcdf.defDim(ncid,'eta_psi'  ,eta_psi  );
    d7  = netcdf.defDim(ncid,'s_rho'    ,s_rho    );
    d8  = netcdf.defDim(ncid,'s_w'      ,s_w      );
    d9  = netcdf.defDim(ncid,'zeta_time',zeta_time);
    d10 = netcdf.defDim(ncid,'v2d_time' ,v2d_time );
    d11 = netcdf.defDim(ncid,'v3d_time' ,v3d_time );
    d12 = netcdf.defDim(ncid,'temp_time',temp_time);
    d13 = netcdf.defDim(ncid,'salt_time',salt_time);
    
    % define variables
    sp_id  = netcdf.defVar(ncid,'spherical'  ,'char'  ,[]);
    Vt_id  = netcdf.defVar(ncid,'Vtransform' ,'NC_INT',[]);
    Vs_id  = netcdf.defVar(ncid,'Vstretching','NC_INT',[]);
    ts_id  = netcdf.defVar(ncid,'theta_s'    ,'double',[]);
    tb_id  = netcdf.defVar(ncid,'theta_b'    ,'double',[]);
    tc_id  = netcdf.defVar(ncid,'Tcline'     ,'double',[]);
    hc_id  = netcdf.defVar(ncid,'hc'         ,'double',[]);
    sr_id  = netcdf.defVar(ncid,'s_rho'      ,'double',d7);
    sw_id  = netcdf.defVar(ncid,'s_w'        ,'double',d8);
    Cr_id  = netcdf.defVar(ncid,'Cs_r'       ,'double',d7);
    Cw_id  = netcdf.defVar(ncid,'Cs_w'       ,'double',d8);
    h_id   = netcdf.defVar(ncid,'h'          ,'double',[d1 d2]);
    lo_id  = netcdf.defVar(ncid,'lon_rho'    ,'double',[d1 d2]);
    la_id  = netcdf.defVar(ncid,'lat_rho'    ,'double',[d1 d2]);

    % time variables
    zt_id  = netcdf.defVar(ncid,'zeta_time'  ,'double',d9);
    vbt_id = netcdf.defVar(ncid,'v2d_time'   ,'double',d10);
    vt_id  = netcdf.defVar(ncid,'v3d_time'   ,'double',d11);
    t_id   = netcdf.defVar(ncid,'temp_time'  ,'double',d12);
    s_id   = netcdf.defVar(ncid,'salt_time'  ,'double',d13);

    % 2D variables
    zw_id  = netcdf.defVar(ncid,'zeta_west'  ,'double',[d2 d9]);
    ze_id  = netcdf.defVar(ncid,'zeta_east'  ,'double',[d2 d9]);
    zs_id  = netcdf.defVar(ncid,'zeta_south' ,'double',[d1 d9]);
    zn_id  = netcdf.defVar(ncid,'zeta_north' ,'double',[d1 d9]);
    ubw_id = netcdf.defVar(ncid,'ubar_west'  ,'double',[d2 d10]);
    ube_id = netcdf.defVar(ncid,'ubar_east'  ,'double',[d2 d10]);
    ubs_id = netcdf.defVar(ncid,'ubar_south' ,'double',[d3 d10]);
    ubn_id = netcdf.defVar(ncid,'ubar_north' ,'double',[d3 d10]);
    vbw_id = netcdf.defVar(ncid,'vbar_west'  ,'double',[d4 d10]);
    vbe_id = netcdf.defVar(ncid,'vbar_east'  ,'double',[d4 d10]);
    vbs_id = netcdf.defVar(ncid,'vbar_south' ,'double',[d1 d10]);
    vbn_id = netcdf.defVar(ncid,'vbar_north' ,'double',[d1 d10]);

    % 3D variables
    uw_id  = netcdf.defVar(ncid,'u_west'  ,'double',[d2 d7 d11]);
    ue_id  = netcdf.defVar(ncid,'u_east'  ,'double',[d2 d7 d11]);
    us_id  = netcdf.defVar(ncid,'u_south' ,'double',[d3 d7 d11]);
    un_id  = netcdf.defVar(ncid,'u_north' ,'double',[d3 d7 d11]);
    vw_id  = netcdf.defVar(ncid,'v_west'  ,'double',[d4 d7 d11]);
    ve_id  = netcdf.defVar(ncid,'v_east'  ,'double',[d4 d7 d11]);
    vs_id  = netcdf.defVar(ncid,'v_south' ,'double',[d1 d7 d11]);
    vn_id  = netcdf.defVar(ncid,'v_north' ,'double',[d1 d7 d11]);
    tw_id  = netcdf.defVar(ncid,'temp_west'  ,'double',[d2 d7 d12]);
    te_id  = netcdf.defVar(ncid,'temp_east'  ,'double',[d2 d7 d12]);
    ts_id  = netcdf.defVar(ncid,'temp_south' ,'double',[d1 d7 d12]);
    tn_id  = netcdf.defVar(ncid,'temp_north' ,'double',[d1 d7 d12]);
    id_sw  = netcdf.defVar(ncid,'salt_west'  ,'double',[d2 d7 d13]);  % avoid conflict with sw_id
    id_se  = netcdf.defVar(ncid,'salt_east'  ,'double',[d2 d7 d13]);
    id_ss  = netcdf.defVar(ncid,'salt_south' ,'double',[d1 d7 d13]);
    id_sn  = netcdf.defVar(ncid,'salt_north' ,'double',[d1 d7 d13]);
    
    
    % define variable attributes    
    % spherical
    netcdf.putAtt(ncid,sp_id,'long_name','grid type logical switch')
    netcdf.putAtt(ncid,sp_id,'flag_values','T, F')
    netcdf.putAtt(ncid,sp_id,'flag_meanings','spherical, Cartesian')
    
    % Vtransform
    netcdf.putAtt(ncid,Vt_id,'long_name','vertical s-coord. transformation equation')
    
    % Vtransform
    netcdf.putAtt(ncid,Vs_id,'long_name','vertical s-coord. stretching function')
    
    % theta_s
    netcdf.putAtt(ncid,ts_id,'long_name','S-coordinate surface control parameter')
    
    % theta_b
    netcdf.putAtt(ncid,tb_id,'long_name','S-coordinate bottom control parameter')
    
    % Tcline
    netcdf.putAtt(ncid,tc_id,'long_name','S-coordinate surface/bottom layer width')
    
    % hc
    netcdf.putAtt(ncid,hc_id,'long_name','S-coordinate parameter, critical depth')
    
    % s_rho
    netcdf.putAtt(ncid,sr_id,'long_name','S-coordinate at rho-points')
    netcdf.putAtt(ncid,sr_id,'positive' ,'up')
    
    % s_w
    netcdf.putAtt(ncid,sw_id,'long_name','S-coordinate at W-points')
    netcdf.putAtt(ncid,sw_id,'positive' ,'up')
    
    % Cs_r
    netcdf.putAtt(ncid,Cr_id,'long_name','S-coordinate stretching curves at rho-points')
    netcdf.putAtt(ncid,Cr_id,'units'    ,'n.d.')
    
    % Cs_w
    netcdf.putAtt(ncid,Cw_id,'long_name','S-coordinate stretching curves at W-points')
    netcdf.putAtt(ncid,Cw_id,'units'    ,'n.d.')
    
    % h
    netcdf.putAtt(ncid,h_id,'long_name','bathymetry at rho-points')
    netcdf.putAtt(ncid,h_id,'units'    ,'meter')
    
    % lon_rho
    netcdf.putAtt(ncid,lo_id,'long_name','longitude at rho-points')
    netcdf.putAtt(ncid,lo_id,'units'    ,'degree_east')
    
    % lat_rho
    netcdf.putAtt(ncid,la_id,'long_name','latitude at rho-points')
    netcdf.putAtt(ncid,la_id,'units'    ,'degree_north')

    % zeta_time
    netcdf.putAtt(ncid,zt_id,'long_name','free-surface time')
    netcdf.putAtt(ncid,zt_id,'units','days since 1900-01-01 00:00:00')
    netcdf.putAtt(ncid,zt_id,'calendar','proleptic_gregorian')
    netcdf.putAtt(ncid,zt_id,'cycle_length',cyc_len)

    % v2d_time
    netcdf.putAtt(ncid,vbt_id,'long_name','2D momentum time')
    netcdf.putAtt(ncid,vbt_id,'units','days since 1900-01-01 00:00:00')
    netcdf.putAtt(ncid,vbt_id,'calendar','proleptic_gregorian')
    netcdf.putAtt(ncid,vbt_id,'cycle_length',cyc_len)

    % v3d_time
    netcdf.putAtt(ncid,vt_id,'long_name','3D momentum time')
    netcdf.putAtt(ncid,vt_id,'units','days since 1900-01-01 00:00:00')
    netcdf.putAtt(ncid,vt_id,'calendar','proleptic_gregorian')
    netcdf.putAtt(ncid,vt_id,'cycle_length',cyc_len)

    % temp_time
    netcdf.putAtt(ncid,t_id,'long_name','potential temperature time')
    netcdf.putAtt(ncid,t_id,'units','days since 1900-01-01 00:00:00')
    netcdf.putAtt(ncid,t_id,'calendar','proleptic_gregorian')
    netcdf.putAtt(ncid,t_id,'cycle_length',cyc_len)

    % salt_time
    netcdf.putAtt(ncid,s_id,'long_name','salinity time')
    netcdf.putAtt(ncid,s_id,'units','days since 1900-01-01 00:00:00')
    netcdf.putAtt(ncid,s_id,'calendar','proleptic_gregorian')
    netcdf.putAtt(ncid,s_id,'cycle_length',cyc_len)
    
    % zeta_west
    netcdf.putAtt(ncid,zw_id,'long_name','free-surface western boundary condition')
    netcdf.putAtt(ncid,zw_id,'units'    ,'meter')
    netcdf.putAtt(ncid,zw_id,'time'     ,'zeta_time')

    % zeta_east
    netcdf.putAtt(ncid,ze_id,'long_name','free-surface eastern boundary condition')
    netcdf.putAtt(ncid,ze_id,'units'    ,'meter')
    netcdf.putAtt(ncid,ze_id,'time'     ,'zeta_time')

    % zeta_south
    netcdf.putAtt(ncid,zs_id,'long_name','free-surface southern boundary condition')
    netcdf.putAtt(ncid,zs_id,'units'    ,'meter')
    netcdf.putAtt(ncid,zs_id,'time'     ,'zeta_time')

    % zeta_north
    netcdf.putAtt(ncid,zn_id,'long_name','free-surface northern boundary condition')
    netcdf.putAtt(ncid,zn_id,'units'    ,'meter')
    netcdf.putAtt(ncid,zn_id,'time'     ,'zeta_time')
    
    % ubar_west
    netcdf.putAtt(ncid,ubw_id,'long_name','2D u-momentum western boundary condition')
    netcdf.putAtt(ncid,ubw_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,ubw_id,'time'     ,'v2d_time')

    % ubar_east
    netcdf.putAtt(ncid,ube_id,'long_name','2D u-momentum eastern boundary condition')
    netcdf.putAtt(ncid,ube_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,ube_id,'time'     ,'v2d_time')

    % ubar_south
    netcdf.putAtt(ncid,ubs_id,'long_name','2D u-momentum southern boundary condition')
    netcdf.putAtt(ncid,ubs_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,ubs_id,'time'     ,'v2d_time')

    % ubar_north
    netcdf.putAtt(ncid,ubn_id,'long_name','2D v-momentum northern boundary condition')
    netcdf.putAtt(ncid,ubn_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,ubn_id,'time'     ,'v2d_time')
    
    % vbar_west
    netcdf.putAtt(ncid,vbw_id,'long_name','2D v-momentum western boundary condition')
    netcdf.putAtt(ncid,vbw_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,vbw_id,'time'     ,'v2d_time')

    % vbar_east
    netcdf.putAtt(ncid,vbe_id,'long_name','2D v-momentum eastern boundary condition')
    netcdf.putAtt(ncid,vbe_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,vbe_id,'time'     ,'v2d_time')

    % vbar_south
    netcdf.putAtt(ncid,vbs_id,'long_name','2D v-momentum southern boundary condition')
    netcdf.putAtt(ncid,vbs_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,vbs_id,'time'     ,'v2d_time')

    % vbar_north
    netcdf.putAtt(ncid,vbn_id,'long_name','2D v-momentum northern boundary condition')
    netcdf.putAtt(ncid,vbn_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,vbn_id,'time'     ,'v2d_time')

    % u_west
    netcdf.putAtt(ncid,uw_id,'long_name','3D u-momentum western boundary condition')
    netcdf.putAtt(ncid,uw_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,uw_id,'time'     ,'v3d_time')

    % u_east
    netcdf.putAtt(ncid,ue_id,'long_name','3D u-momentum eastern boundary condition')
    netcdf.putAtt(ncid,ue_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,ue_id,'time'     ,'v3d_time')

    % u_south
    netcdf.putAtt(ncid,us_id,'long_name','3D u-momentum southern boundary condition')
    netcdf.putAtt(ncid,us_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,us_id,'time'     ,'v3d_time')

    % u_north
    netcdf.putAtt(ncid,un_id,'long_name','3D v-momentum northern boundary condition')
    netcdf.putAtt(ncid,un_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,un_id,'time'     ,'v3d_time')
    
    % v_west
    netcdf.putAtt(ncid,vw_id,'long_name','3D v-momentum western boundary condition')
    netcdf.putAtt(ncid,vw_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,vw_id,'time'     ,'v3d_time')

    % v_east
    netcdf.putAtt(ncid,ve_id,'long_name','3D v-momentum eastern boundary condition')
    netcdf.putAtt(ncid,ve_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,ve_id,'time'     ,'v3d_time')

    % v_south
    netcdf.putAtt(ncid,vs_id,'long_name','3D v-momentum southern boundary condition')
    netcdf.putAtt(ncid,vs_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,vs_id,'time'     ,'v3d_time')

    % v_north
    netcdf.putAtt(ncid,vn_id,'long_name','3D v-momentum northern boundary condition')
    netcdf.putAtt(ncid,vn_id,'units'    ,'meter second-1')
    netcdf.putAtt(ncid,vn_id,'time'     ,'v3d_time')
    
    % temp_west
    netcdf.putAtt(ncid,tw_id,'long_name','potential temperature western boundary condition')
    netcdf.putAtt(ncid,tw_id,'units'    ,'Celsius')
    netcdf.putAtt(ncid,tw_id,'time'     ,'temp_time')

    % temp_east
    netcdf.putAtt(ncid,te_id,'long_name','potential temperature eastern boundary condition')
    netcdf.putAtt(ncid,te_id,'units'    ,'Celsius')
    netcdf.putAtt(ncid,te_id,'time'     ,'temp_time')

    % temp_south
    netcdf.putAtt(ncid,ts_id,'long_name','potential temperature southern boundary condition')
    netcdf.putAtt(ncid,ts_id,'units'    ,'Celsius')
    netcdf.putAtt(ncid,ts_id,'time'     ,'temp_time')

    % temp_north
    netcdf.putAtt(ncid,tn_id,'long_name','potential temperature northern boundary condition')
    netcdf.putAtt(ncid,tn_id,'units'    ,'Celsius')
    netcdf.putAtt(ncid,tn_id,'time'     ,'temp_time')
    
    % salt_west
    netcdf.putAtt(ncid,id_sw,'long_name','salinity western boundary condition')
    netcdf.putAtt(ncid,id_sw,'time'     ,'salt_time')
    
    % salt_east
    netcdf.putAtt(ncid,id_se,'long_name','salinity eastern boundary condition')
    netcdf.putAtt(ncid,id_se,'time'     ,'salt_time')

    % salt_north
    netcdf.putAtt(ncid,id_sn,'long_name','salinity northern boundary condition')
    netcdf.putAtt(ncid,id_sn,'time'     ,'salt_time')

    % salt_south
    netcdf.putAtt(ncid,id_ss,'long_name','salinity southern boundary condition')
    netcdf.putAtt(ncid,id_ss,'time'     ,'salt_time')

    % define global attributes
    varid = netcdf.getConstant('GLOBAL');
    
    netcdf.putAtt(ncid,varid,'type'    ,'ROMS boundary forcing file')
    netcdf.putAtt(ncid,varid,'title'   ,'Patagonia model boundary forcing')
    netcdf.putAtt(ncid,varid,'grd_file','roms_grd_rivers.nc')
    netcdf.putAtt(ncid,varid,'date'    , datestr(today))
    
    % obtain data for boundary conditions files.
    % set files from which to obtain data.
    bname = 'roms_bry_1900.nc'   % file from which to read initial data
    gname = 'roms_grd_rivers.nc' % file from which to read grid data
    
    % set scalar variables
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
        
    % setup finished, close file.
    netcdf.close(ncid)  

    % Read and grid databoundary data and write to appropriate  
    % variables in newly created boundary file.

    % Read grid data
    h    = ncread(gname, 'h');  % bathymetry
    lo_r = ncread(gname, 'lon_rho'); 
    la_r = ncread(gname, 'lat_rho');

    % Read and set time data
    bry_time  = ncread(bname,'bry_time');
    zeta_time = bry_time;
    v2d_time  = bry_time;
    v3d_time  = bry_time;
    temp_time = bry_time;
    salt_time = bry_time;
    
    % Read 2-D variables
    zeta_west  = ncread(bname, 'zeta_west');
    zeta_east  = ncread(bname, 'zeta_east');
    zeta_south = ncread(bname, 'zeta_south');
    zeta_north = ncread(bname, 'zeta_north');
    ubar_west  = ncread(bname, 'ubar_west');
    ubar_east  = ncread(bname, 'ubar_east');
    ubar_south = ncread(bname, 'ubar_south');
    ubar_north = ncread(bname, 'ubar_north');
    vbar_west  = ncread(bname, 'vbar_west');
    vbar_east  = ncread(bname, 'vbar_east');
    vbar_south = ncread(bname, 'vbar_south');
    vbar_north = ncread(bname, 'vbar_north');
    
    % Read 3-D variables
    u_west  = ncread(bname, 'u_west');
    u_east  = ncread(bname, 'u_east');
    u_south = ncread(bname, 'u_south');
    u_north = ncread(bname, 'u_north');
    v_west  = ncread(bname, 'v_west');
    v_east  = ncread(bname, 'v_east');
    v_south = ncread(bname, 'v_south');
    v_north = ncread(bname, 'v_north');
    temp_west  = ncread(bname, 'temp_west');
    temp_east  = ncread(bname, 'temp_east');
    temp_south = ncread(bname, 'temp_south');
    temp_north = ncread(bname, 'temp_north');
    salt_west  = ncread(bname, 'salt_west');
    salt_east  = ncread(bname, 'salt_east');
    salt_south = ncread(bname, 'salt_south');
    salt_north = ncread(bname, 'salt_north');

    % write data to variables
    % scalars
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

    % time
    ncwrite(ncname, 'zeta_time', zeta_time)
    ncwrite(ncname, 'v2d_time' , v2d_time)
    ncwrite(ncname, 'v3d_time' , v3d_time)
    ncwrite(ncname, 'temp_time', temp_time)
    ncwrite(ncname, 'salt_time', salt_time)

    % 2D variables
    ncwrite(ncname, 'zeta_west'  , zeta_west)
    ncwrite(ncname, 'zeta_east'  , zeta_east)
    ncwrite(ncname, 'zeta_south' , zeta_south)
    ncwrite(ncname, 'zeta_north' , zeta_north)
    ncwrite(ncname, 'ubar_west'  , ubar_west)
    ncwrite(ncname, 'ubar_east'  , ubar_east)
    ncwrite(ncname, 'ubar_south' , ubar_south)
    ncwrite(ncname, 'ubar_north' , ubar_north)
    ncwrite(ncname, 'vbar_west'  , vbar_west)
    ncwrite(ncname, 'vbar_east'  , vbar_east)
    ncwrite(ncname, 'vbar_south' , vbar_south)
    ncwrite(ncname, 'vbar_north' , vbar_north)

    % 3D variables
    ncwrite(ncname, 'u_west'     , u_west)
    ncwrite(ncname, 'u_east'     , u_east)
    ncwrite(ncname, 'u_south'    , u_south)
    ncwrite(ncname, 'u_north'    , u_north)
    ncwrite(ncname, 'v_west'     , v_west)
    ncwrite(ncname, 'v_east'     , v_east)
    ncwrite(ncname, 'v_south'    , v_south)
    ncwrite(ncname, 'v_north'    , v_north)
    ncwrite(ncname, 'temp_west'  , temp_west)
    ncwrite(ncname, 'temp_east'  , temp_east)
    ncwrite(ncname, 'temp_south' , temp_south)
    ncwrite(ncname, 'temp_north' , temp_north)
    ncwrite(ncname, 'salt_west'  , salt_west)
    ncwrite(ncname, 'salt_east'  , salt_east)
    ncwrite(ncname, 'salt_south' , salt_south)
    ncwrite(ncname, 'salt_north' , salt_north)

end

