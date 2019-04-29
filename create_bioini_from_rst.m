% create_bioini_from_rst.m
%
% Create an initial conditions file for the biology and physics using
% the spun-up physics file.
%
% Input: bio_file(string): name of file containing bio initial
%                                   conditions.
%           rst_file(string): name of restart file containing spun-up 
%                                   physics.  Biological variables will be
%                                   written to this file.  THIS FILE WILL
%                                   BE MODIFIED.
%
% Author: Z. Wallace
% Created: 3.13.19
% Last edit: 4.18.19  --> Add automatic checking of time records in restart
%                                     file

function create_bioini_from_rst(bio_file,rst_file)
    
    if nargin ~= 2
        msg1 = "Usage: create_bioini_from_rst(bio_file,rst_file).  ";
        msg2 = "Use the HELP function for more details";
        msg = strcat(msg1,msg2);
        error(msg)
    end
    
    nc = ncinfo(bio_file);
    
    % Find index of first biological variable (NO3), since we don't want
    % the physical variables.
    cntr = 1;
    while ~strcmp(nc.Variables(cntr).Name,'NO3')
        cntr = cntr + 1;
    end
    bio_indx = cntr;
    
    % pre-allocate cell arrays for speed
    bio_length = length(bio_indx:length(nc.Variables));
    varnames = cell(1,bio_length);
    vardata_ini = cell(1,bio_length);
    vardata = cell(1,bio_length);
    
   % Read names of biological variables
    for i=bio_indx:length(nc.Variables)
        varnames{i-bio_indx+1}=nc.Variables(i).Name;
    end
    disp('Biological Variables:');
    disp(varnames);
    
    for i=1:length(varnames)
        vardata_ini{i}=ncread(bio_file,varnames{i});
    end

    % Set initial time record in new vardata array
    for i=1:length(vardata_ini)
        vardata{i} = vardata_ini{i};
    end
    
    % Read total number of time records saved in restart file.
    % Copy the biological variables this number of times so
    % physical and biological array time dimensions are compatible.
    num_time_records = length(ncread(rst_file,'ocean_time'));
    
    % Make N copies of the initial condition field to fit with the
    % number of time records in the physics file
    for i=1:length(vardata)
        disp(strcat("Processing variable: ", varnames{i}));
        for j=2:num_time_records
        disp(strcat("Creating time record: ",int2str(j)));
        vardata{i}=cat(4,vardata{i},vardata_ini{i});
        end
    end
    
    disp(strcat("Opening ",rst_file," for writing..."));
    
    % Create biological variables in restart file and write data.
    ncid = netcdf.open(rst_file,'WRITE');
    
    netcdf.reDef(ncid);
    
    % Define variable dimensions
    tdim=netcdf.inqDimID(ncid,'ocean_time');
    sdim=netcdf.inqDimID(ncid,'s_rho');
    edim=netcdf.inqDimID(ncid,'eta_rho');
    xdim=netcdf.inqDimID(ncid,'xi_rho');
    dims=[xdim,edim,sdim,tdim];

    % Create variables in file
    varid = zeros(1,length(varnames));
    for i=1:length(varnames)
        display(strcat("Creating variable ",varnames{i}));
        varid(i) = netcdf.defVar(ncid,varnames{i},'NC_DOUBLE',dims);
    end
    
    netcdf.endDef(ncid)
    
    % Write data to variables
    for i=1:length(varnames)
        display(strcat("Writing data to variable ",varnames{i}));
        netcdf.putVar(ncid,varid(i),vardata{i});
    end
    
    netcdf.close(ncid);
    
end