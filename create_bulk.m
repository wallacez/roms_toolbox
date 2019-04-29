%=========================================================================%
% create_bulk.m
% Build a netcdf bulk forcing file.
%
% use: function [bulk_file] = create_bulk(fname, mat_file)
%
% Input:
%      fname    (string): name of bulk forcing netcdf file to create.
%      mat_file (string): name of .mat file from which to pull the bulk data.
%      cycle_len(double): number of days after which to recycle forcing data.
%
% Output: 
%      ncfile: netcdf bulk forcing file.
%
% Author: Z. Wallace
% Last edit: 4 April 2018
%            23 April 2018 -> generalize reads from forcing data file.
%                          -> Add cycle_length to parameter list.
%=========================================================================%


function bulk_file = create_bulk(fname, mat_file, cycle_len)

         % Create new NetCDF 4 file with variables
         ncid = netcdf.create(fname, 'NetCDF4');
         
         % Create global attributes
         
         
         % Create dimensions
         xi_rho = netcdf.defDim(ncid,'xi_rho',482);
         eta_rho = netcdf.defDim(ncid,'eta_rho',770);
         xi_psi = netcdf.defDim(ncid,'xi_psi',481);
         eta_psi = netcdf.defDim(ncid,'eta_psi',769);
         xi_u = netcdf.defDim(ncid,'xi_u',481);
         eta_u = netcdf.defDim(ncid,'eta_u',770);
         xi_v = netcdf.defDim(ncid,'xi_v',482);
         eta_v = netcdf.defDim(ncid,'eta_v',769);
         bulk_time = netcdf.defDim(ncid,'bulk_time',netcdf.getConstant('NC_UNLIMITED'));

         % Create variables
         dim_arr = [xi_rho,eta_rho,bulk_time];

         blk_id=netcdf.defVar(ncid,'bulk_time','NC_DOUBLE',bulk_time);
         tair_id=netcdf.defVar(ncid,'tair','NC_DOUBLE',dim_arr);
         rhum_id=netcdf.defVar(ncid,'rhum','NC_DOUBLE',dim_arr);
         prate_id=netcdf.defVar(ncid,'prate','NC_DOUBLE',dim_arr);
         wspd_id=netcdf.defVar(ncid,'wspd','NC_DOUBLE',dim_arr);
         radlw_id=netcdf.defVar(ncid,'radlw','NC_DOUBLE',dim_arr);
         radlwin_id=netcdf.defVar(ncid,'radlw_in','NC_DOUBLE',dim_arr);
         radsw_id=netcdf.defVar(ncid,'radsw','NC_DOUBLE',dim_arr);
         sustr_id=netcdf.defVar(ncid,'sustr','NC_DOUBLE',dim_arr);
         svstr_id=netcdf.defVar(ncid,'svstr','NC_DOUBLE',dim_arr);
         uwnd_id=netcdf.defVar(ncid,'uwnd','NC_DOUBLE',dim_arr);
         vwnd_id=netcdf.defVar(ncid,'vwnd','NC_DOUBLE',dim_arr);
         icecov_id=netcdf.defVar(ncid,'icecov','NC_DOUBLE',dim_arr);

         %% Add attributes 
         % Global
         grd_file = '/home/server/pi/homes/zwallace/ROMS/Project_Patagonia/roms_grd.nc';
         ftype = 'ROMS heat flux bulk forcing file -- 2001-2002';

         netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'title','Patagonia Shelf');
         netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'date',date);
         netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'grd_file',grd_file);
         netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'type',ftype);
         
         % Variable-specific
         netcdf.putAtt(ncid,blk_id,'long_name','bulk formulation execution time');
         netcdf.putAtt(ncid,blk_id,'units','days');
         netcdf.putAtt(ncid,blk_id,'cycle_length',cycle_len);

         netcdf.putAtt(ncid,tair_id,'long_name','surface air temperature');
         netcdf.putAtt(ncid,tair_id,'units','Celsius');

         netcdf.putAtt(ncid,rhum_id,'long_name','relative humidity');
         netcdf.putAtt(ncid,rhum_id,'units','fraction');

         netcdf.putAtt(ncid,prate_id,'long_name','precipitation rate');
         netcdf.putAtt(ncid,prate_id,'units','kg m-2 s-1');

         netcdf.putAtt(ncid,wspd_id,'long_name','wind speed at 10 m');
         netcdf.putAtt(ncid,wspd_id,'units','m s-1');

         netcdf.putAtt(ncid,radlw_id,'long_name','net outgoing longwave radiation');
         netcdf.putAtt(ncid,radlw_id,'units','Watts m-2');
         netcdf.putAtt(ncid,radlw_id,'positive','upward flux, cooling water');

         netcdf.putAtt(ncid,radlwin_id,'long_name','downward longwave radiation');
         netcdf.putAtt(ncid,radlwin_id,'units','Watts m-2');
         netcdf.putAtt(ncid,radlwin_id,'positive','downward flux, warming water');

         netcdf.putAtt(ncid,radsw_id,'long_name','shortwave radiation');
         netcdf.putAtt(ncid,radsw_id,'units','Watts m-2');
         netcdf.putAtt(ncid,radsw_id,'positive','downward flux, heating water');

         netcdf.putAtt(ncid,sustr_id,'long_name','surface u-momentum stress');
         netcdf.putAtt(ncid,sustr_id,'units','Newton m-2');

         netcdf.putAtt(ncid,svstr_id,'long_name','surface v-momentum stress');
         netcdf.putAtt(ncid,svstr_id,'units','Newton m-2');

         netcdf.putAtt(ncid,uwnd_id,'long_name','U-direction wind');
         netcdf.putAtt(ncid,uwnd_id,'units','m s-1');

         netcdf.putAtt(ncid,vwnd_id,'long_name','V-direction wind');
         netcdf.putAtt(ncid,vwnd_id,'units','m s-1');

         netcdf.putAtt(ncid,icecov_id,'long_name','ERA-interim ice coverage');
         netcdf.putAtt(ncid,icecov_id,'units','none');
         
         netcdf.close(ncid)

         %% Write data to variables
         blk_data = mat_file;
         load(blk_data);

         ncwrite(fname,'bulk_time',bt_add);
         ncwrite(fname,'tair',tair_add);
         ncwrite(fname,'rhum',rhum_add);
         ncwrite(fname,'prate',prate_add);  % already converted cm d-1 -> kg m-2 s-1
         ncwrite(fname,'wspd',wspd_add);
         ncwrite(fname,'radsw',radsw_add);
         ncwrite(fname,'radlw',radlw_add);
         ncwrite(fname,'radlw_in',radlw_in_add);
         ncwrite(fname,'sustr',sustr_add);
         ncwrite(fname,'svstr',svstr_add);
         ncwrite(fname,'uwnd',uwnd_add);
         ncwrite(fname,'vwnd',vwnd_add);
         ncwrite(fname,'icecov',icecov_add);
end
