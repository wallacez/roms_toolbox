% read_var.m
%
% Read the specified timeseries of data from a variable.  The default
% behavior (i.e. nfiles is not specified) is to read data from all files in
% the directory.  To specify a specific subseries of data, give the number
% of files to read as nfiles.  The program will read data from files 
% 1:nfiles.
%
% Input: var_name - name of variable to read from files.
%           nfiles         - number of files to process.
%
% Output: varout - variable containing the specified data within the
%                  specified time series.
%
% Author: Z. Wallace | 17 August 2018
% Rev. his: 8.17.18 --> Double quotes to ensure path names are strings,
%                       not chars.
%           8.20.18  --> Refactored code to use "dir" command for generic
%                        filename compatibility.  This removed the need to
%                        use double quotes to ensure string-like behavior of
%                        pathnames.  Switched back to single quotes for
%                        compatibility with Matlab v2013a.
%                    --> Set "nfiles" as an optional argument; default
%                        behavior is to read all averages files in directory.
%           10.19.18 -->Variable name written to stdout.
%           12.12.18 -->specify starting and ending files.
%           2.18.19  --> Support for 1D (scalar) variables

function varout = read_var(varname, varargin)
    
    % Ensure proper number of input arguments.
    if nargin < 1 || nargin > 3 || nargin == 2
        err = 'USAGE: varout = read_var(var_name[, start_file, end_file]);';
        error(err)
    end
    
    flist = dir('*avg*');  % get only averages files.
    
    if nargin == 1    % default read entire time series of data
        nfiles = numel({flist.name});
        msg1 = strcat('Reading variable: ',varname);
        msg2 = strcat('Number of files to read: ',int2str(nfiles));
        disp(msg1)
        disp(msg2)
    elseif nargin == 3
        start_file = varargin{1};
        end_file   = varargin{2};
        msg = strcat("Reading files: ", ...
                                flist(start_file).name, ...
                                " - ", ...
                                flist(end_file).name);
        disp(msg)
    end
    
    % Read in variable data.  Concatenate variables along the time
    % dimension.
    if nargin == 1
        start_file = 1;
        end_file = nfiles;
    else
         % pass, since start and end files are set 
    end
    
    for i = start_file:end_file 
        fname = flist(i).name;
        msg = strcat('Reading file: ',fname);
        disp(msg)
        if i == start_file
            varout = ncread(fname,varname);
            % Get the number of dimensions to determine the appropriate one
            % on which to concatenate (3rd for 2-D variables, 4th for 3-D
            % variables).
            ndims = size(size(varout));
            catdim = ndims(2);  % position of time dimension
            if catdim == 2
                % 1-D (scalar) variables should be concatenated along
                % first, not second dimension.
                catdim = 1;
            end
        else
            varout = cat(catdim,varout,ncread(fname,varname));
        end
    end
    
end
