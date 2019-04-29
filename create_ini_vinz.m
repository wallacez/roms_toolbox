%
%=========================================================================%
% create_ini_vinz.m
%
% Create initial conditions file for Biology Model.
%
% Initial conditions for the physics are monthly averages from 2000-2007 
% taken from OFES.  
% Initial conditions for the biology are monthly climatologies taken from
% WOA(09?)
%
% Author: Z. Wallace
% Created: 27 March 2018
% Last edit: 5.16.18 --> Changed LZOO, SZOO and ocean_time.
%                 8.31.18 --> changed to absolute path name for fid.
%                 9.17.18 --> resolved function name ambiguity (was create_ini).
%                 11.5.18 --> added 3rd Fe pool, changed Fep to DetFe
%                 12.31.18--> split fid into physics and biology files to
%                             use spun-up physics.
%                         --> removed NaN mask from restart file, since model
%                             requires mask to be 0 (blows up at first time
%                             step with -NaN everywhere if NaN mask present.
%                         --> changed ocean_time to reflect model spinup
%                             date (23 June 2006).
%                         --> updated pathnames to reflect new directory
%                              structure.
%                 2.1.19 --> Allow ocean_time to start model run in April
%                                  instead of June.
%                            --> Percent of total chlorophyll signal going
%                                   to phyto and dia variable with depth.
%                 2.5.19 --> Chlorophyll from spinup April -> June
%=========================================================================%
%


function create_ini_vinz(fname)
    
    ncid = netcdf.create(fname,'NETCDF4');  % create and open new ini file
    netcdf.close(ncid);
    
    %
    %----------------------------------------------------------------------
    % Create dimensions
    %----------------------------------------------------------------------
    %
    % grid dimensions
    xi_rho  = 482;
    xi_u    = 481;
    eta_rho = 770;
    eta_v   = 769;
    s_rho   = 40;
    s_w     = 41;
    time    = 1;
        
    %
    %----------------------------------------------------------------------
    % Create variables and add attributes
    %----------------------------------------------------------------------    
    %
    % Ocean time
    nccreate(fname,'ocean_time', ...
            'Dimensions',{'ocean_time' time}, ... 
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ocean_time','long_name','seconds since 1900-01-01 00:00:00');
    ncwriteatt(fname,'ocean_time','units','seconds');
 
    % Free-surface
    nccreate(fname,'zeta', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'zeta','long_name','free-surface');
    ncwriteatt(fname,'zeta','units','meter');

    % 2D u-direction velocity
    nccreate(fname,'ubar', ...
            'Dimensions',{'xi_u' xi_u 'eta_rho' eta_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ubar','long_name','Vertically integrated u-momentum component');
    ncwriteatt(fname,'ubar','units','meter second-1');

    % 2D v-direction velocity
    nccreate(fname,'vbar', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_v' eta_v 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'vbar','long_name','Vertically integrated v-momentum component');
    ncwriteatt(fname,'vbar','units','meter second-1');

    % 3D u-direction velocity
    nccreate(fname,'u', ...
            'Dimensions',{'xi_u' xi_u 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'u','long_name','u-momentum component');
    ncwriteatt(fname,'u','units','meter second-1');

    % 3D v-direction velocity
    nccreate(fname,'v', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_v' eta_v 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'v','long_name','v-momentum component');
    ncwriteatt(fname,'v','units','meter second-1');

    % temperature
    nccreate(fname,'temp', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'temp','long_name','potential temperature');
    ncwriteatt(fname,'temp','units','meter second-1');

    % salinity
    nccreate(fname,'salt', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'salt','long_name','salinity');
    ncwriteatt(fname,'salt','units','PSU');

    % NO3
    nccreate(fname,'NO3', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NO3','long_name','nitrate concentration');
    ncwriteatt(fname,'NO3','units','mmol N m-3');

    % NH4
    nccreate(fname,'NH4', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'NH4','long_name','NH4 concentration');
    ncwriteatt(fname,'NH4','units','mmol N m-3');

    % Phytoplankton
    nccreate(fname,'Phytoplankton', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Phytoplankton','long_name','Phytoplankton concentration');
    ncwriteatt(fname,'Phytoplankton','units','mmol N m-3');

    % Diatoms
    nccreate(fname,'Diatoms', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Diatoms','long_name','Diatom concentration');
    ncwriteatt(fname,'Diatoms','units','mmol N m-3');

    % Bacteria
    nccreate(fname,'Bacteria', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Bacteria','long_name','Bacteria concentration');
    ncwriteatt(fname,'Bacteria','units','mmol N m-3');

    % Chlorophyll-a
    nccreate(fname,'Chlorophyll', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Chlorophyll','long_name','Chlorophyll concentration');
    ncwriteatt(fname,'Chlorophyll','units','mg Chl m-3');

    % Diatom Chlorophyll
    nccreate(fname,'ChlD', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'ChlD','long_name','ChlD concentration');
    ncwriteatt(fname,'ChlD','units','mg Chl m-3');

    % DON
    nccreate(fname,'DON', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DON','long_name','DON concentration');
    ncwriteatt(fname,'DON','units','mmol N m-3');

    % DOC
    nccreate(fname,'DOC', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DOC','long_name','DOC concentration');
    ncwriteatt(fname,'DOC','units','mmol C m-3');

    % DetN
    nccreate(fname,'DetN', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetN','long_name','DetN concentration');
    ncwriteatt(fname,'DetN','units','mmol N m-3');

    % DetC
    nccreate(fname,'DetC', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetC','long_name','DetC concentration');
    ncwriteatt(fname,'DetC','units','mmol C m-3');

    % Silica
    nccreate(fname,'Silica', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Silica','long_name','Silica concentration');
    ncwriteatt(fname,'Silica','units','mmol Si m-3');

    % LDetN
    nccreate(fname,'LDetN', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetN','long_name','LDetN concentration');
    ncwriteatt(fname,'LDetN','units','mmol N m-3');

    % LDetC
    nccreate(fname,'LDetC', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetC','long_name','LDetC concentration');
    ncwriteatt(fname,'LDetC','units','mmol C m-3');

    % LDetS
    nccreate(fname,'LDetS', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetS','long_name','LDetS concentration');
    ncwriteatt(fname,'LDetS','units','mmol Si m-3');

    % Microzooplankton
    nccreate(fname,'Microzooplankton', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Microzooplankton','long_name','Microzooplankton concentration');
    ncwriteatt(fname,'Microzooplankton','units','mmol N m-3');

    % Mesozooplankton
    nccreate(fname,'Mesozooplankton', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'Mesozooplankton','long_name','Mesozooplankton concentration');
    ncwriteatt(fname,'Mesozooplankton','units','mmol N m-3');

    % Dissolved Iron
    nccreate(fname,'FeD', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'FeD','long_name','Dissolved iron concentration');
    ncwriteatt(fname,'FeD','units','umol Fe m-3');

    % Detrital Iron
    nccreate(fname,'DetFe', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'DetFe','long_name','Detrital iron concentration');
    ncwriteatt(fname,'DetFe','units','umol Fe m-3');
    
    % Large Detrital Iron
    nccreate(fname,'LDetFe', ...
            'Dimensions',{'xi_rho' xi_rho 'eta_rho' eta_rho 's_rho' s_rho 'time' time},...
            'Datatype','double', ...
            'Format','netcdf4');
    ncwriteatt(fname,'LDetFe','long_name','Large detrital iron concentration');
    ncwriteatt(fname,'LDetFe','units','umol Fe m-3');

    %
    %----------------------------------------------------------------------
    % Read data into variables.  Where data are unavailable, populate with 
    %   analytical functions.
    %----------------------------------------------------------------------    
    %
    myroot   = '/home/server/pi/homes/zwallace/ROMS/Project_Patagonia';
    outroot   = '/home/morrison/data2/zach';
    biodir = 'orig/vinz_in_3-21-18';
    %chldir = 'spinup';
    %chldir= 'spinup_MY25_bio_tempArctic4_alphaLutz_scavDOC_graz0.05_thetamExud_dustx2M05_AprStart';
    hdir = 'bulk';

    fid_bio = fullfile(myroot,biodir,'rst_PHYSIQUE_NEMURO_CHLAseawif_2001Jun_WOA_VinzFeC_ini.nc'); % get name of Vincent's initial conditions file
%    fid_phys= '/home/morrison/data2/zach/phys-QCORRBarnier-swradIcecorr-windminuscurr-emmiss-dqdtVar-rhumFrac-MY25/patagonia_rst.nc';
%    fid_phys= '/home/morrison/data2/zach/phys-QCORRBarnier-swradIcecorr-windminuscurr-emmiss-dqdtVar-rhumFrac/patagonia_rst.nc';
    fid_phys= fullfile(outroot,'phys-QCORRBarnier-swradIcecorr-windminuscurr-emmiss-dqdtVar-rhumFrac-MY25/patagonia_avg_0059.nc');
    %fid_chl = fullfile(myroot,chldir,'chlo_grd_3D_April.mat');
    %fid_chl1 = fullfile(outroot,chldir,'patagonia_avg_0002.nc');
    %fid_chl2 = fullfile(outroot,chldir,'patagonia_avg_0003.nc');
    fid_grd = fullfile(myroot,hdir,'roms_grd_rivers.nc');
    
    %idx_end = 10;  % last index in restart file  
    idx_end = 4;  % 15 April 2006 from avg file
    start_2D=[1,1,idx_end];
    count_2D=[Inf,Inf,1];
    start_3D=[1,1,1,idx_end];
    count_3D=[Inf,Inf,Inf,1];

    C2N = 6.625;   % Redfield C:N ratio

    % Time variables
%    scrum_time = ncread(fid,'scrum_time');
%    ocean_time = scrum_time-367*86400;  % subtract a year to account for indexing difference
%    ocean_time = (datenum('24-Jun-2001')-datenum('01-Jan-1900'))*86400.; % Model start: 24 June 0101
    ocean_time = (datenum('23-Jun-2006')-datenum('01-Jan-1900'))*86400.;  % Model start: 23 June 0106
 %    ocean_time = (datenum('15-Apr-2006')-datenum('01-Jan-1900'))*86400.;  % Model start: 15 April 0106
    
     % Read in bathymetery data for variable %chl going to biomass as a
     % function of depth.
     h = ncread(fid_grd,'h');
     depth_le75 = zeros(xi_rho,eta_rho); % preallocate for speed
     for i=1:xi_rho
        for j=1:eta_rho
            if h(i,j)  <= 75
                depth_le75(i,j) = 1;
            else
                depth_le75(i,j) = 0;
            end
        end
     end
    
    % Load in April chlorophyll data
% 	chl_tot = load(fid_chl,'chlo_grd_3D_April');
%    chl_tot = chl_tot.chlo_grd_3D_April;
%    chl_tot = chl_tot.chlo_grd_3D_April;
    % Load in June spin-up data
% 	chl_tot1 = ncread(fid_chl1,'Chlorophyll',[1,1,1,5],[Inf,Inf,Inf,Inf]) + ...
%                       ncread(fid_chl1,'ChlD',[1,1,1,5],[Inf,Inf,Inf,Inf]);
%     chl_tot2 = ncread(fid_chl2,'Chlorophyll',[1,1,1,1],[Inf,Inf,Inf,4]) + ...
%                       ncread(fid_chl2,'ChlD',[1,1,1,1],[Inf,Inf,Inf,4]);
%     chl_tot = cat(4,chl_tot1,chl_tot2);
%     chl_tot = mean(chl_tot,4);
    chl_tot = 1.59*(ncread(fid_bio,'SPHY')+ncread(fid_bio,'LPHY'));
    % Vary %Chl per phyto group with depth
    chls = zeros(xi_rho,eta_rho,s_rho);
    chld = zeros(xi_rho,eta_rho,s_rho);
    for k = 1:s_rho
        for j = 1:eta_rho
            for i = 1:xi_rho
                if depth_le75(i,j) 
                    % 80% of total chl signal to diatoms
                    chls(i,j,k) = chl_tot(i,j,k)*0.2;
                    chld(i,j,k) = chl_tot(i,j,k)*0.8;
                else
                    % 20% to diatoms
                    chls(i,j,k) = chl_tot(i,j,k)*0.8;
                    chld(i,j,k) = chl_tot(i,j,k)*0.2;
                end
            end
        end
    end
    % Read into variables for writing
    Chlorophyll = chls;
    Chlorophyll(isnan(Chlorophyll)) = 0;
    ChlD = chld;
    ChlD(isnan(ChlD)) = 0;
    
    % zeta
    zeta  = ncread(fid_phys,'zeta',start_2D,count_2D);
    zeta(isnan(zeta)) = 0;

    % ubar
    ubar  = ncread(fid_phys,'ubar',start_2D,count_2D);
    ubar(isnan(ubar)) = 0;

    % vbar
    vbar  = ncread(fid_phys,'vbar',start_2D,count_2D);
    vbar(isnan(vbar)) = 0;

    % u
    u  = ncread(fid_phys,'u',start_3D,count_3D);
    u(isnan(u)) = 0;

    % v
    v  = ncread(fid_phys,'v',start_3D,count_3D);
    v(isnan(v)) = 0;

    % temp
    temp  = ncread(fid_phys,'temp',start_3D,count_3D);
    temp(isnan(temp)) = 0;

    % salt
    salt  = ncread(fid_phys,'salt',start_3D,count_3D);
    salt(isnan(salt)) = 0;

    % NO3
    NO3  = ncread(fid_bio,'NO3');
    NO3(isnan(NO3)) = 0;

    % NH4
%    NH4 = ncread(fid,'NH4');
    NH4  = NO3*0.01;
    NH4(isnan(NH4)) = 0;

    % Phytoplankton
    %Phytoplankton  = ncread(fid_bio,'SPHY');
    Phytoplankton = Chlorophyll./1.59;
    Phytoplankton(isnan(Phytoplankton)) = 0;

    % Diatoms
    %Diatoms  = ncread(fid_bio,'LPHY');
    Diatoms  = ChlD./1.59;
    Diatoms(isnan(Diatoms)) = 0;

    % Bacteria
    Bacteria = 0.1*Phytoplankton;
    Bacteria(isnan(Bacteria)) = 0;

    % Chlorophyll
    % Chlorophyll = Phytoplankton*1.59;
    % ChlD
    %ChlD = Diatoms*1.59;

    % DON
%    DON  = ncread(fid,'DON');
    DON = NH4;
    DON(isnan(DON)) = 0;

    % DOC
    DOC  = C2N*DON;

    % DetN
%    DetN  = ncread(fid,'PON');
    DetN = DON;
    DetN(isnan(DetN)) = 0;

    % DetC
    DetC = C2N*DetN;

    % Silica
    Silica  = ncread(fid_bio,'SiOH');
    Silica(isnan(Silica)) = 0;

    % LDetN
    LDetN = DetN;
    LDetN(isnan(LDetN)) = 0;

    % LDetC
    LDetC = C2N*LDetN;

    % LDetS
%    LDetS  = ncread(fid,'OPAL');
    LDetS = Diatoms*0.1;

    % Microzooplankton
%    Microzooplankton  = ncread(fid,'SZOO');
    Microzooplankton= 0.01*Phytoplankton;
    Microzooplankton(isnan(Microzooplankton)) = 0;

    % Mesozooplankton
%    Mesozooplankton  = ncread(fid,'LZOO');
    %Mesozooplankton(xi_rho,eta_rho,s_rho) = CONST_VAL;
    Mesozooplankton= 0.01*Diatoms;
    Mesozooplankton(isnan(Mesozooplankton)) = 0;

    % Dissolved Iron
    FeD = ncread(fid_bio,'FeD');
    FeD(isnan(FeD)) = 0;

    % Particulate Iron
    %   combine small and large particulate from NEMURO
%    FeSp = ncread(fid,'FeSp');
%    FeLp = ncread(fid,'FeLp');
%
%    FeP = FeSp+FeLp;
    DetFe = 0.1 * FeD;
    DetFe(isnan(DetFe)) = 0;
    
    % Large detrital iron
    LDetFe = 0.1 * DetFe;
    LDetFe(isnan(LDetFe)) = 0;


    %
    %----------------------------------------------------------------------
    % Write data to variables.
    %----------------------------------------------------------------------    
    %
    % bry_time
    ncwrite(fname,'ocean_time',ocean_time);

    % zeta
    ncwrite(fname,'zeta',zeta);

    % ubar
    ncwrite(fname,'ubar',ubar);

    % vbar
    ncwrite(fname,'vbar',vbar);

    % u
    ncwrite(fname,'u',u);

    % v
    ncwrite(fname,'v',v);

    % temp
    ncwrite(fname,'temp',temp);

    % salt
    ncwrite(fname,'salt',salt);

    % NO3
    ncwrite(fname,'NO3',NO3);

    % NH4
    ncwrite(fname,'NH4',NH4);

    % Phytoplankton
    ncwrite(fname,'Phytoplankton',Phytoplankton);

    % Diatoms
    ncwrite(fname,'Diatoms',Diatoms);

    % Bacteria
    ncwrite(fname,'Bacteria',Bacteria);

    % Chlorophyll
    ncwrite(fname,'Chlorophyll',Chlorophyll);

    % ChlD
    ncwrite(fname,'ChlD',ChlD);

    % DON
    ncwrite(fname,'DON',DON);

    % DOC
    ncwrite(fname,'DOC',DOC);

    % DetN
    ncwrite(fname,'DetN',DetN);

    % DetC
    ncwrite(fname,'DetC',DetC);

    % Silica
    ncwrite(fname,'Silica',Silica);

    % LDetN
    ncwrite(fname,'LDetN',LDetN);

    % LDetC
    ncwrite(fname,'LDetC',LDetC);

    % LDetS
    ncwrite(fname,'LDetS',LDetS);

    % Microzooplankton
    ncwrite(fname,'Microzooplankton',Microzooplankton);

    % Mesozooplankton
    ncwrite(fname,'Mesozooplankton',Mesozooplankton);

    % Dissolved Iron
    ncwrite(fname,'FeD',FeD);
    
    % Detrital Iron
    ncwrite(fname,'DetFe',DetFe);
    
    % Large Detrital iron
    ncwrite(fname,'LDetFe',LDetFe);


end  % create_ini_vinz
