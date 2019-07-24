% compute_std_seawifs.m
%
% Compute standard deviation for chlorophyll for specified years.
%
% Input: year_start(int): first year to read data from
%        year_end(int): last year to read data from
%
% Output: chl_std (3D-array): Standard deviation of chlorophyll between
%                             years specified (i.e. year_start and
%                             year_end) for spring (SON) and summer (DJF).
%
% Author: Z. Wallace | 24 Jul 2019

function std_chl_season = compute_std_seawifs(year_start, year_end)
    if year_start < 1997
        error('Data can be no earlier than 1997');
    end
    if year_end > 2010
        error('Data can be no later than 2010');
    end
    
    files=dir;
    
    num_years = year_end-year_start+1;  % inclusive
    
    % Initialize chlorphyll arrays
    chl_sep=zeros(480,528,num_years);
    chl_oct=zeros(480,528,num_years);
    chl_nov=zeros(480,528,num_years);
    chl_dec=zeros(480,528,num_years);
    chl_jan=zeros(480,528,num_years);
    chl_feb=zeros(480,528,num_years);
    
    chl_spring=zeros(480,528,num_years);
    chl_summer=zeros(480,528,num_years);
    std_chl_season = zeros(480,528,2);  % to hold seasonal std.

    cntrs=ones(1,6);
    for t=4:152  % start at 4 to ignore files prepended with a '.'
        year = str2double(char(files(t).name(1:4)));
        month = str2double(char(files(t).name(6:7)));
        if year >=year_start && year <=year_end
            switch(month)
                case 09
                    chl_sep(:,:,cntrs(1)) = ncread(files(t).name,'chlor_a');
                    cntrs(1) = cntrs(1) + 1;
                case 10
                    chl_oct(:,:,cntrs(2)) = ncread(files(t).name,'chlor_a');
                    cntrs(2) = cntrs(2) + 1;
                case 11
                    chl_nov(:,:,cntrs(3)) = ncread(files(t).name,'chlor_a');
                    cntrs(3) = cntrs(3) + 1;
                case 12
                    chl_dec(:,:,cntrs(4)) = ncread(files(t).name,'chlor_a');
                    cntrs(4) = cntrs(4) + 1;
                case 01
                    chl_jan(:,:,cntrs(5)) = ncread(files(t).name,'chlor_a');
                    cntrs(5) = cntrs(5) + 1;
                case 02
                    chl_feb(:,:,cntrs(6)) = ncread(files(t).name,'chlor_a');
                    cntrs(6) = cntrs(6) + 1;
            end           
        end
    end
    
    for t=1:num_years
        chl_spring(:,:,t) = (chl_sep(:,:,t)+chl_oct(:,:,t)+chl_nov(:,:,t))/3;
        chl_summer(:,:,t) = (chl_dec(:,:,t)+chl_jan(:,:,t)+chl_feb(:,:,t))/3;
    end
    
    
    % Compute seasonal means
%     chl_spring = cat(3,chl_sep_mean,chl_oct_mean,chl_nov_mean);
%     chl_summer = cat(3,chl_dec_mean,chl_jan_mean,chl_feb_mean);
    
    std_spring = std(chl_spring,0,3);
    std_summer = std(chl_summer,0,3);
    
    std_chl_season(:,:,1) = std_spring;
    std_chl_season(:,:,2) = std_summer;
    
    return
end