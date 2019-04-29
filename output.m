% output.m
%
% 2-D (depth, lon) pcolor visualizations of model runs in time.  Vertical 
% grid transformed from s-grid to z-grid coordinateds using the set_depth 
% routine.
%
% Usage: 
%   output(AH_flag,inp_var,lat[,t])
%
% Parameters:
%   AH_flag: dtype=string, loop through average ('A') or history ('H') file
%   inp_var: dtype=string, tracer variable to observe
%   lat:     dtype=integer, latitudinal cross-section value
%   t:       dtype=integer, timestep to observe (optional)
%
%   AH_flag determines the file type (average or history) from which to
%       pull tracer variable data.
%   var determines which tracer variable to get from the netcdf file.
%   lat is the constant latitude to observe over (i.e., the latitude over 
%       which a cross section is taken).  
%   t (optional) sets the view to a specfic timestep, where 1 < t < Ntimes 
%       where Ntimes is the numeber of timesteps in the file.
%
% Author: Z. Wallace
% Last edit: 5 July 2017


function [] = output(AH_flag,inp_var,lat,varargin)
% Equivalent to function output(AH_flag,var,lat,varargin);

% get grid data from proper 'ocean_xxx' file 
if strcmp(AH_flag, 'A')
    ncid = netcdf.open('ocean_avg.nc','NOWRITE');
    %ncid = netcdf.open('../Project_Fennel/ocean_avg.nc','NOWRITE');
    
    % get lat/lon data
%     dimid = netcdf.inqDimID(ncid,'xi_rho');
%     [dimname, dimlen] = netcdf.inqDim(ncid,dimid);
%     xi_rho_pts = 1:dimlen;
    
    dimid = netcdf.inqDimID(ncid,'eta_rho');
    [dimname, dimlen] = netcdf.inqDim(ncid,dimid);
    eta_rho_pts = 1:dimlen;
    
    dimid = netcdf.inqDimID(ncid,'s_rho');
    [dimname, dimlen] = netcdf.inqDim(ncid,dimid);
    s_rho_pts = 1:dimlen;
    
    % parameters for calculating z grid
    varname = 'h';
    varid = netcdf.inqVarID(ncid,varname);
    var = netcdf.getVar(ncid,varid,'double');
    
elseif strcmp(AH_flag, 'H')
    ncid = netcdf.open('ocean_his.nc','NOWRITE');
    %ncid = netcdf.open('../Project_Fennel/ocean_his.nc','NOWRITE');
    
    % get lat/lon data
%     dimid = netcdf.inqDimID(ncid,'xi_rho');
%     [dimname, dimlen] = netcdf.inqDim(ncid,dimid);
%     xi_rho_pts = 1:dimlen;

    dimid = netcdf.inqDimID(ncid,'eta_rho');
    [dimname, dimlen] = netcdf.inqDim(ncid,dimid);
    eta_rho_pts = 1:dimlen;
    
    dimid = netcdf.inqDimID(ncid,'s_rho');
    [dimname, dimlen] = netcdf.inqDim(ncid,dimid);
    s_rho_pts = 1:dimlen;
    
    % parameters for calculating z grid
    varname = 'h';
    varid = netcdf.inqVarID(ncid,varname);
    var = netcdf.getVar(ncid,varid,'double');    

else 
    msg = 'Flags to use are A or H';
    error(msg);
end 

% Parameters to change from s-grid to z-grid
V_transform = 2;
V_stretching = 4;
theta_s = 3;
theta_b = 0;
hc = 25;
N = 16;
igrid = 1;
h = var;

% calculate z
[z_grid] = set_depth(V_transform,V_stretching,theta_s,theta_b,hc,N,igrid,h);

% calculate number of timesteps per day (dtdays)
sec_per_day = 86400;

varname = 'dt';     % [s/timestep]
varid   = netcdf.inqVarID(ncid,varname);
dt      = netcdf.getVar(ncid,varid,'double');

dtdays = sec_per_day/dt;    % [timesteps/day]

varname = 'ntimes';     % total timesteps
varid   = netcdf.inqVarID(ncid,varname);
Ntimes  = netcdf.getVar(ncid,varid,'double');

varname = 'nAVG';     % # timesteps between time-averaged records
varid   = netcdf.inqVarID(ncid,varname);
nAVG    = netcdf.getVar(ncid,varid,'double');

varname = 'nHIS';     % # timesteps between snapshot records
varid   = netcdf.inqVarID(ncid,varname);
nHIS    = netcdf.getVar(ncid,varid,'double');

% get variable of interest data
varname = inp_var;
varid = netcdf.inqVarID(ncid,varname);
var   = netcdf.getVar(ncid,varid,'double');

%[x,y] = meshgrid(xi_rho_pts, eta_rho_pts);
[y,z] = meshgrid(eta_rho_pts, s_rho_pts);

% determine step size to use for contourf
min_var = min(min(min(min(var))));
max_var = max(max(max(max(var))));
diff = max_var - min_var     
cf_step = diff/100      % uncomment for debugging purposes

% Check number of arguments are valid
if nargin > 4
    msg = 'Maximumn number of inputs exceeded';
    error(msg)
end

if nargin < 3
    msg = 'Too few input arguments';
    error(msg);
end

% Calculation of day and maximum time dimension changes whether one is
% looking at the average output or the history file.
if strcmp(AH_flag,'A')
    max_time = length(var(1,1,1,:));
    day = nAVG/dtdays;
    if nargin == 3
        t = max_time;
        % Loop through time
        for i=1:t
%            figure
%            pcolor(y',squeeze(z_grid(lat,:,:)),squeeze(var(lat,:,:,i)));
            %cf_step = (max(var(lat,:,:,i))-min(var(lat,:,:,i)))/100
            contourf(y',squeeze(z_grid(lat,:,:)),squeeze(var(lat,:,:,i)),...
                     min_var:cf_step:max_var);
            shading flat;colorbar;caxis([min_var,max_var]);
            title(strcat(varname, ' | ',...
                ' Day: ', num2str(i*day)))
            pause(1)
        end
        
    elseif nargin == 4
        t = varargin{1};
        % error checking
        if(t > max_time)
            msg = strcat('t must be <= ',' ',int2str(max_time));
            error(msg)
            return
        end
        % Display specified timestep
        figure
%        pcolor(y',squeeze(z_grid(lat,:,:)),squeeze(var(lat,:,:,t)));
        var(lat,1,1,t) % surface value of var at timestep t
        contourf(y',squeeze(z_grid(lat,:,:)),squeeze(var(lat,:,:,t)),...
                 min_var:cf_step:max_var);
        shading flat;colorbar;caxis([min_var,max_var]);
        title(strcat(varname, ' | ',...
            ' Day: ', num2str(t*day)))
    end
    
elseif strcmp(AH_flag,'H')
    max_time = length(var(1,1,1,:));
    day = nHIS/dtdays;
    if nargin == 3
        t = max_time;
        %plot_step = dtdays/24;     % [#timesteps/hr]
        % Loop through time
        %for i=1:plot_step:t        % display snapshot once per hour
        for i=1:t
            %pcolor(y',squeeze(z_grid(lat,:,:)),squeeze(var(lat,:,:,i)));
            contourf(y',squeeze(z_grid(lat,:,:)),squeeze(var(lat,:,:,i)),...
                min_var:cf_step:max_var);
            shading flat;colorbar;caxis([min_var,max_var]);
            title(strcat(varname, ' | ',...
                ' Day: ', num2str(i*day)))
            pause(1)
        end
        
    elseif nargin == 4
        t = varargin{1};
        % error checking
        if(t > Ntimes)
            msg = strcat('t must be <= ',' ',int2str(Ntimes));
            error(msg)
            return
        end
        day = nHIS/dtdays;
        
        % Display specified timestep
        figure
%        pcolor(y',squeeze(z_grid(lat,:,:)),squeeze(var(lat,:,:,t)));
%        shading interp;colorbar;caxis([min(min(min(min(var)))),max(max(max(max(var))))]);
        contourf(y',squeeze(z_grid(lat,:,:)),squeeze(var(lat,:,:,t)),...
                 min_var:cf_step:max_var);
        colorbar;shading flat;caxis([min_var,max_var]);
        title(strcat(varname, ' | ',...
            ' Day: ', num2str(t*day)))
    end
end
    
netcdf.close(ncid);

end