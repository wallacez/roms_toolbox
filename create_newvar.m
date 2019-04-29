%------------
filename='roms_ini_oceantime.nc'
filename1='roms_ini_oceantime_zw.nc'
%-----------

%Read phytoplankton field and convert it to Chlorophyll
ncid = netcdf.open(filename,'NC_NOWRITE');
varid = netcdf.inqVarID(ncid,'SPHY');
phy= netcdf.getVar(ncid,varid); 
bact=0.1*phy;

%Write the chlorophyll field to the initial condition file
ncid1 = netcdf.open(filename1,'NC_WRITE');
    % ##### Add variables
    % ..... Get dimensions
    dimid1=  netcdf.inqDimID(ncid,'time');
    dimid2=  netcdf.inqDimID(ncid,'s_rho');
    dimid3=  netcdf.inqDimID(ncid,'eta_rho');
    dimid4=  netcdf.inqDimID(ncid,'xi_rho');
    % -- Leave define mode and enter data mode to write data
    %   netcdf.endDef(ncid);
    % -- Create variable
    netcdf.reDef(ncid1);
%    varID = netcdf.defVar(ncid1,'Bacteria','single',[dimid4, dimid3, dimid2, dimid1]);
    varID = netcdf.defVar(ncid1,'Bacteria','single',[dimid1, dimid2, dimid3, dimid4]);
    % -- Leave define mode and enter data mode to write data
    netcdf.endDef(ncid1);
    % -- Write data to variable
    netcdf.putVar(ncid1,varID,bact);
    % -- Re-enter define mode
    netcdf.reDef(ncid1);
    % -- Create attributes associated with the variable
    netcdf.putAtt(ncid1,varID,'long_name','Bacteria Concentration');
    netcdf.putAtt(ncid1,varID,'units','micromole_nitrogen meter-3');
    netcdf.putAtt(ncid1,varID,'time','ocean_time');
    netcdf.putAtt(ncid1,varID,'field','Bacteria, scalar, series');
    % #######    close netCDF
netcdf.close(ncid1);
netcdf.close(ncid);

