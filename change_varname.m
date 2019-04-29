% ======================================================================= % 
% change_varname.m
% Redefine a netcdf variable name.  Assumes dimensions remain the same.
%
% Parameters:
%   fname (string): name of netcdf file within which is the variable to
%                   change.
%    old_var (string):   Old name of variable.
%    new_var (string):   New name of variable.
%    long_name (string): String that goes in the "long_name" field in the
%                        netcdf file.  Since this is often not the same as 
%                        the name as the new variable (new_var), it is best
%                        to manually enter it.
%    units (string):     String that goes in the "units" field in the 
%                        netcdf file.  Since this is often not the same as 
%                        the name as the new variable (new_var), it is best
%                        to manually enter it.
%
% Author: Z. Wallace
% Last edit: 23 August 2017
% ======================================================================= %

function change_varname(fname, old_var, new_var, long_name, units)
    
    % Ensure correct number of command line args.
    if nargin ~= 5
        msg = ['Incorrect number of parameters.  ', ...
               'Usage: change_varname(fname, old_var, new_var,' ....
               ' long_name, units)'];
        error(msg);
    end
    
    % open netcdf file and enter rename mode
    ncid = netcdf.open(fname, 'WRITE');
    netcdf.reDef(ncid);
    
    % read old variable from netcdf file
    varid = netcdf.inqVarID(ncid, old_var);
    
    % rename old variable to new one
    netcdf.renameVar(ncid, varid, new_var);
    
    % rename attributes 
    netcdf.putAtt(ncid,varid,'long_name',long_name);
    netcdf.putAtt(ncid,varid,'units',units);
    netcdf.putAtt(ncid,varid,'time','ocean_time');
    netcdf.putAtt(ncid,varid,'field',strcat(new_var, ', scalar,',...
                                            ' series'));
%    netcdf.defVarFill(ncid,varid,false,1.e37); Doesn't work on netcdf-3
%    files
                                        
    
    % verify variable name has changed (remove semicolon)
    [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,varid);
    
    netcdf.endDef(ncid)
    netcdf.close(ncid)

end