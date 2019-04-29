% gen_tracer_month_mean.m
%
% Generate monthly means for tracer variable of choice.
% Input: varname(string): Tracer for which to generate monthly means.
%            month_dix(int): Month over which to generate mean tracer value,
%                                     where 1 is January and 12 is December.
% Output: varout(n-d array): Time-averaged values of tracer "varname" over
%                                            month "month".
%              Additionally, pseudocolor plots of varout are generated.
%
% Author: Z. Wallace
% Created: 2.18.19

function varout = gen_tracer_month_mean(varname, month_idx)

    if nargin ~= 2
        msg = 'USAGE:  gen_tracer_month_mean(varname,month)';
        error(msg);
    end
    % Set general file path variables
    f1 = '/home/server/homes/pi/zwallace';
    f2 = 'ROMS/Project_Patagonia/bulk'; 
    grdfile = fullfile(f1,f2,'roms_grd_rivers.nc');
    
    % Read in grid variables for plotting
    lons = ncread(grdfile,'lon_rho');
    lats = ncread(grdfile,'lat_rho');
    h = ncread(grdfile,'h');
    
    month = {'Jan','Feb','Mar','Apr','May', 'Jun',...
                     'Jul','Aug','Sep','Oct','Nov','Dec'};
    ot = datestr(read_var('ocean_time')./86400,3);  % read in time data as months only
    t = find(ot == month{month_idx},1);  % find first index of tracer variable in selected month
    
    % Read in variables.  Combine certain variables if appropriate.
    if strncmpi(varname,'Total C',7)
        % Combine Chlorophyll terms
        varout = read_var('Chlorophyll')+read_var('ChlD');
        % Remove negatives so 'pcolor' will work when plotting the log10
        % transform
        varout(varout<0) = NaN;
    elseif strncmpi(varname,'Total P',7)
        % Combine Phytoplankton terms
        varout = read_var('Phytoplankton')+read_var('Diatoms');
    else
        varout = read_var(varname);
    end
    
    % Get size of variable to determine time dimension, tdim.
    ndims = size(size(varout)); 
    tdim = ndims(2);
    
    % Plot variables
    figure;colormap('jet');
    if strncmpi(varname,'Total C',7)
        % log10-transform the chlorophyll data
        pcolor(lons,lats,log10(mean(varout(:,:,end,t:t+5),tdim)));
        shading flat;caxis([-1 0.7]);%colorbar
    else
        pcolor(lons,lats,mean(varout(:,:,end,t:t+5),tdim));shading flat;%colorbar
        caxis([0 0.5]) % Bacteria
        %caxis([0 1.5]) % FeD
    end
    xticklabels([]);
    yticklabels([]);
    hold on
    contour(lons,lats,h,[200 200],'w','LineWidth',2);  % Plot 200-m isobath
    title(strcat(varname, " | ", month{month_idx}));
    set(gca,'PlotBoxAspectRatio',[0.8046 1 0.8046]);
    
    return
end