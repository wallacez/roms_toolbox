%=========================================================================%
% interp2rhogrid.m
% interpolate u- or v-grid to rho-grid.  
%
% use: function [varout] = interp2rhogrid(input_grid, init_grid)
%
% Input:
%      input_grid (3-D array): data to regrid.
%      init_grid (string):     current grid of input_grid data (e.g. 'u').
%
% Output: 
%      rho_grid: data from input_grid regridded to rho-grid.
%
% Author: Z. Wallace
% Last edit: 4 April 2018
%=========================================================================%


function grd_out = interp2rhogrid(grd_in, grd_init)

    msg = ' was passed in as the initial grid';
    
    switch grd_init
    
        case('u')
            disp(strcat(grd_init,msg));
            
            [x,y,z] = size(grd_in);
            grd_out = zeros(x+1,y,z);
            
            % Set end members of new grid to end members of initial grid.  
            % This replaces the lost information at the grid boundary that
            % occurs during interpolation.
            grd_out(1,:,:)   = grd_in(1,:,:);
            grd_out(end,:,:) = grd_in(end,:,:);
            
            % interpolate data from u grid to rho grid
            grd_out(2:end-1,:,:) = (grd_in(1:end-1,:,:) + grd_in(2:end,:,:))*0.5;
    
        case('v')
            disp(strcat(grd_init,msg))

            [x,y,z] = size(grd_in);
            grd_out = zeros(x,y+1,z);
            
            % Set end members of new grid to end members of initial grid.  
            % This replaces the lost information at the grid boundary that
            % occurs during interpolation.
            grd_out(:,1,:)   = grd_in(:,1,:);
            grd_out(:,end,:) = grd_in(:,end,:);
            
            % interpolate data from v grid to rho grid
            grd_out(:,2:end-1,:) = (grd_in(:,1:end-1,:) + grd_in(:,2:end,:))*0.5;
        
        otherwise
            msg = "Invalid init_grid parameter.  Please use 'u' or 'v'";
            disp(msg);
            return 
    end
    
end
