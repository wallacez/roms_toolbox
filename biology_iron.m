
clear all

%-----------------------------------------------------------------------
% Options
%-----------------------------------------------------------------------
%
EXUD_MINVAL = 1;  % Lower biomass threshold for exudation
EXUD_LESS = 0;  % 11% d-1 max exudation
QUADRATIC_MORT = 1;  % Quadratic mortality for phyto and dia
QUAD_MORT_ZOO = 1;
GRAZ_INI_LOW = 1;  % Grazers set to 1% of producers
PLOTS = 1;          % Flag to turn on plotting
IRON = 1;           % Flag to turn on iron model
FEC_VAR = 0;   % Turn on variable Fe:C
DUST = 1;
K_LOSS = 1;
APRIL_START = 0; % Start light and temp from April instead of June
FEC_VAR_LIGHT = 0;  % Vary Fe:C with parz (Strzepek and Price 2000)
FAST_REMIN = 0;  % Instantaneous remin. of LDetS --> Si(OH)4 at bottom
SIL_REGEN_BAC = 1;
CHL_PAR_ATT_SUM = 1;
ZOOL_GRAZ_NIGHT = 1;
ALPHA_LIM_FE = 0;
LIEBIG = 0;
ONEILL = 0;
COCCOLITHOPHORES = 1;

if K_LOSS ==1 
    K0 = 0;
    K1 = 0;
    K2 = 1;
    K3 = 0;
    if K0 == 1
        % Constant Chl loss
        chl_dmg = 0.01; % [d-1]
    elseif K1 == 1
        % Loss a function of T-limited growth (-k*aj/Vpt)
        chls_dmg = 0.04;  %[d-1]  Alvarez range 0.01-0.5, chose 0.1
        %chld_dmg = chls_dmg*0.6;
        chld_dmg = chls_dmg;
        chlc_dmg = chls_dmg;
    elseif K2 == 1
        % Loss proportional to light and theta (-k*theta*I)
        % Alvarez range [1d-8 - 2d-6]
        chl_dmg = 1d-8*86400*12;  %[mg C mg-1 Chl (W m-2 d)-1]
    elseif K3 == 1
        % Loss proportional to light saturation (-k*I)
        chl_dmg = 1d-8*86400;  % [m2 J-1] scaled to days
    end
end
if DUST == 1
    MW_Fe = 55.85;  % g / mol
    MW_Si = 28.1;     % g / mol
    Fe_solub = 0.02;   % Soluble Fe fraction in dust
    Si_solub = 0.075;  % Soluble Si fraction in dust
    gSi2gFe = 8.8;  % g Si / g Fe
end
if IRON == 1
    IRON_BACTERIA = 1;
    IRON_SCAVENGING = 1;
    if IRON_SCAVENGING == 1
        SCAV_DOC = 1;
    end
    if FEC_VAR == 1
        SA_Fe=0.6;
        %SB_Fe=8;  % Fe:C for small range 2-8 umol Fe / mmol C
        SB_Fe=13;  % Fe:C for small range ~4-20 umol Fe / mmol C for FeD [0,2]nM
        LA_Fe=0.6;
        %LB_Fe=64;   % Fe:C range for large 20 - 100 umol Fe / mmol C for FeD [0,2]nM
        LB_Fe=13;   % Fe:C range for large 2 - 20 umol Fe / mmol C
        CA_Fe = SA_Fe;
        CB_Fe = SB_Fe;  % Cocco Fe:C same as small phy
        
    end
end

if EXUD_MINVAL == 1
     % Minimum biomass required for exudation.
    exud_min = 0.1;  % mmol N m-3
end

if QUADRATIC_MORT == 1
      phy_mort_max = 0.1d0;   % Maximum 10% mortality d-1
      dia_mort_max = 0.1d0;    % Maximum 10% mortality d-1
      Km_phy = 1.d0; 
      Km_dia = 1.d0;
      if COCCOLITHOPHORES == 1
         coc_mort_max = 0.1d0;
        Km_coc = 1.d0;        
      end
end

%-----------------------------------------------------------------------
% Parameters
%-----------------------------------------------------------------------
%
N = 20;             % Number of vertical levels
%h = 10;             % bottom depth in meters

if IRON == 1
    nbt = 20;       % Number of biological variables with iron
else
    nbt = 17;       % Number of biological variables without iron
end
if IRON == 1
    nsink = 7;      % Number of biological variables that sink with iron
    if COCCOLITHOPHORES == 1
        nsink = 8;  % Include PIC
    end
else
    nsink = 5;      % Number of biological variables that sink without iron
end
dayspersecond = 1./86400; 
bioparfrac = 0.43;  % Fraction of visible light (par0) active in photosynthesis (400nm-700nm wave lengths)
%acof = 0.851;       % Eppley 1972 % acof = 0.59   % Norberg 2004
acof = 0.75;      % d-1
acof_dia = 0.8; % d-1  Berelson 2001
bcof = 0.0633;      % Norberg 2004
topt_phy = 18;      % Topt_phy = 19 % Boyd et al. 2013
topt_dia = 10;       %topt_dia = 11;       %  Topt_dia = 12
 w_phy = 15;         %  w_phy = 32    % Boyd et al. 2013
 w_dia = 11;
 topt_coc = 20;
 w_coc = 15;
%Topt_phy = 14;
%Topt_dia = 9;
%Topt_dia = 7;
%Topt_dia = 8;
%Vptmax = 0.35;  % sets Vpt(Topt) of ~0.8 d-1;
%Vptmax = 1;  % sets Vpt(Topt) of ~2.5d-1, from Eppley
%Vptdiamax = 0.85;  % sets Vptdia(Topt) of ~1.5 d-1;
%Vptmax = 0.5;  
%Vptdiamax = 0.95;  
%Vptdiamax = 1.2;  
rho0 = 1025.0d0;
cp = 3985.0d0;
cnp = 6.0d0;        % Carbon:Nitrogen ratio for phytoplankton
cnd = 4.6d0;        % Carbon:Nitrogen ratio for DOM
cnb = 5.38d0;       % Carbon:Nitrogen ratio for bacteria
cndia = 6.625d0;    % Carbon:Nitrogen ratio for diatoms
cnz = 5.02d0;       % Carbon:Nitrogen ratio for zooplankton
cncoc = 5.d0;
% vp = 0.8d0;         % Phytoplankton maximum growth rate
% vpdia = 1.5d0;      % Diatom maximum growth rate
kdia_no3 = 1.0d0;
%kdia_no3 = 3.0d0;    % Jiang
kp_no3 = 0.59d0;    % Half-saturation for phytoplankton NO3 uptake
kp_nh4 = 3.39d-2;
kdia_nh4 = kp_nh4;  % Kdia_NH4 = 0.5d0    % Jiang
kc_no3 = kp_no3;
kc_nh4 = kp_nh4;
kdia_si  = 3.0d0;   
%kdia_si   = kdia_no3;
kb_c = 0.51d0;      % Half saturation constant for DOCC
kb_don = 1.34d0;    % Bacteria half saturation constant for DON uptake
kb_a = 1.0d0;       % Bacteria half saturation constant for NH4 uptake
kz_i = 1.2d0;       % Half saturation constant for ingestion
kl_i = 0.69d0;      % Large zooplankton half saturation constant for ingestion
anitrof = 0.05;
anitrof_l = 0.025;  % euphotic zone nitrification rate (Jiang et al. 2003)
anitrof_d = 0.1;     % dark nitrification rate (Jiang et al. 2003)
snpd = 2.0;         % Silica:Nitrogen ratio
ldetndis = 0.1;
ldetcdis = 0.1;
ldetsdis = 0.1;
psi = 1.48d0;       % Phytoplankton ammonium inhibition parameter
attsw = 0.04d0;     % Light attenuation due to seawater (1/m)
attphyt = 0.02d0;   % Light attenuation by phytoplankton
attphytdia = 0.02d0;    % Light attenuation by diatoms
if COCCOLITHOPHORES == 1
    attphytcoc = 0.02d0;    % Light attenuation by coccolithophores
end
alpha = 1.3d0;    % Initial slope of the P-I curve
%alpha = 2.76d0;    % Schourup-Kristenset et al. 2014
alphadia = 1.d0;   % Initial slope of the P-I curve for diatoms
%alphadia = 2.3d0;   % Schourup-Kristenset et al. 2014
if COCCOLITHOPHORES == 1
    alphac = alpha;  % Initial slope of P-I curve for coccolithophores
end
prefz_b = 1.14d0;   % Zooplankton preference for bacteria rel to phy
prefz_d = 0.90d0;   % Zooplankton preference for detritus rel to phy 
prefz_c = 0.90d0;   % Zoo. preference for cocco relative to phy
prefl_z = 0.94d0;   % Meso preference for micro rel. to diatoms
prefl_p = 0.9d0;    % Meso preference for phy rel. to small z
prefl_c = 0.9d0;   % Meso preference for cocco rel. to diatoms
gmax = 2.06d0;      % Maximum grazing rate
grazl = 0.6d0;      % Large zooplannkton grazing rate
vb_amino = 0.253d0; % Bacteria maximum growth amino DOM
ggen = 0.628d0;     % Bacteria gross growth efficiency of nitrogen DOC
ggec = 0.18d0;      % Bacteria gross growth efficiency of carbon DOC
mud_dom = 3.1d-2;   % Detrital breakdown to DON
mub_dom = 2.54d-2;  % Detrital breakdown by bacteria
mub_a = 4.76d-3;    % Bacteria specific excretion rate
mub_d = 8.99d-3;    % Bacteria mortality rate
mup_d = 0.072d0;    % Phytoplankton specific mortality rate
amud1 = mup_d;
muz_a = 0.12d0;     % Zooplankton specific mort/excretion rate
%muz_a = 0.2d0;
muzg_a = 0.106d0;   % Growth specific excretion rate
%mul = 0.1d0;        % Large zooplankton mortality
mul = 0.2d0;
mulg_a = 0.1d0;     % Large zooplankton growth specific excretion rate
omegal = 0.15d0;    % Large zooplankton fraction to DON
omegaz_d = 0.18d0;  % Detrital fraction of zooplankton mortality
betaz_b = 0.87d0;   % Zooplankton assimilation efficiency of bacteria
betaz_p = 0.61d0;   % Zooplankton assimilation efficiency of phytoplankton
betaz_d = 0.45d0;   % Zooplankton assimilation efficiency of detritus
betal_dia = 0.72d0; % Large zooplankton assimilation efficiency of diatoms
betal_z = 0.44d0;   % Large zooplankton assimilation efficiency of zooplankton
betal_p = 0.72d0;   % Large zooplankton assimilation efficiency of phytoplankton
if COCCOLITHOPHORES == 1
    betaz_c = betaz_d;          % Zoo. assimilation efficiency of cocco    
    betal_c = betal_z;  % Large zoo assimilation efficiency of cocco
    PIC2C = 0.2;  % Inorganic to organic carbon ratio in coccos (mmol C / mmol C)
    dissPIC = 0.1;  % d-1
end
if EXUD_LESS == 1
    gamma1p = 2.4d-2;   % Fraction of total primary production exuded
    gamma2p = 8.6d-2;  % Fraction of growth-specific primary production exuded
else
    gamma1p = 6.6d-2;   % Fraction of total primary production exuded
    gamma2p = 11.5d-2;  % Fraction of growth-specific primary production exuded
end
gamma3p = 0.538d0;  % Fraction of nutrient-limited primary production exuded
if EXUD_LESS == 1
    gamma1dia = 2.4d-2; % Diatom fraction of total primary production exuded
    gamma2dia = 8.6d-2; % Diatom fraction of growth-specific primary production exuded
else
    gamma1dia = 6.6d-2;  % Diatom fraction of total primary production exuded
    gamma2dia = 11.2d-2; % Diatom fraction of growth-specific primary production exuded
end
gamma3dia = 6.6d-2;  % Diatom fraction of nutrient limited primary production exuded
if COCCOLITHOPHORES == 1
    gamma1c = gamma1p;
    gamma2c = gamma2p;
    gamma3c = gamma3p;
end
epsilonz_a = 0.6d0; % Ammonium fraction of zooplankton excretion
epsilonl = 0.51d0;  % Large zooplankton fraction to ammonium
mx_chlcr = 4.15d-2; % Maximum chlorophyll-a to Carbon ratio
mx_chlcrd = 5.5d-2; % Maximum chlorophyll-a to Carbon ratio (diatoms)
if COCCOLITHOPHORES ==1 
    mx_chlcrc = mx_chlcr;
end

%  Set time-stepping according to the number of iterations.
dt = 900;    % seconds
num_days = 365;
dtdays=dt*dayspersecond;  % [d ts-1]
time_loops = num_days*86400/dt;

if IRON == 1
%  Fe:C for small and large zoo and phyto [mmol Fe (mol C)-1]
%  From Lancelot et al. 2001.
    %k1fe = 0.6d0;     % Half saturation constant for nanoflag. DFe uptake
    %k2fe = 1.2d0;       % Half saturationconstant for diatom DFe uptake
    k1fe = 0.3;
    k2fe = 0.6;
    kcfe = 0.1;
    if FEC_VAR == 1 || FEC_VAR_LIGHT == 1
        fec_phy_diag = zeros(time_loops+1,N);
        fec_dia_diag = zeros(time_loops+1,N);
    else
        fec_phy = 0.0025d0;
        %fec_phy = 0.1d0;
        fec_dia = 0.002d0;
        fec_coc = 0.0025d0;
    end
    fec_zoo = 0.0024d0;  % Lancelot et al. 2001
    fec_bac = 0.0091;   % Tortell et al. 1996
    ldetfedis = 0.1;    % LDetFe, from Jiang et al. 2013
    
    % Diagnostic terms for BSi remineralization
    fl11_16_diag = zeros(time_loops+1,N);  % Mesozoo resp.
    fl16_12_diag = zeros(time_loops+1,N);  % BSi diss/remin
    fl12_13_diag = zeros(time_loops+1,N);  % Si(OH)4 uptake (BSi production)
    % Diagonistic term for chlorophyll light attenuation
    parz_diag = zeros(time_loops+1,N);

    if K_LOSS == 1
        fl8_dmg_diag = zeros(time_loops+1,N);
        fl17_dmg_diag = zeros(time_loops+1,N);
    end
    if IRON_SCAVENGING == 1
        %
        % Scavenging terms for a constant removal of iron via scavenging.
        %      rscav:  Base scavenging rate [d-1]
        %      Fscav:  Fraction of scavenging goint to LDetFe
        %      KeqFeL: Equilibrium constant for the formation of FeL
        %      LTot:   Total ligand concentration [umol L m-3].  Serves as upper
        %              bound for [FeD].
        %
        fscav  = 0.9;        % Moore and Braucher, 2008
        KeqFeL = 1.d4;      % Rue and Bruland, 1995
        if SCAV_DOC == 1
            %
            %Total ligand concentration (LTot) a function of of DOC as in PISCES,
            % and scavenging rate a function of particulate carbon.
            % Assumes a background ligand concentration of 0.6 nM.  
            %
              scav_min = 3.d-5;  % [d-1, PISCES]
              Lam_Fe = 5.d-3;  % [d-1 uM-1, PISCES] slope of scavenging rate of Fe
              LMin = 0.6;  % [nM]
        else
            % constant scavenging rate and ligand concentration
            rscav  = 0.01;  % [d-1]
            LTot = 1;
        end
    else
        fscav = 0;
        rscav = 0;
    end
end

% Set diagnostic growth arrays
aj_diag = zeros(time_loops+1,N);
ajd_diag = zeros(time_loops+1,N);

% Specify x,y location for 1-D run
%x = 100; y = 150;
x = 250; y = 400;
%x = 242; y = 503;
%x = 444; y = 660;

% Read in bathymetry
f1 = '/home/server/homes/pi/zwallace/';
fgrd = 'ROMS/Project_Patagonia/bulk';
grdfile = fullfile(f1,fgrd,'roms_grd_rivers.nc');
h = ncread(grdfile,'h');
h = h(x,y);  % bottom depth in meters
latr=ncread(grdfile,'lat_rho');
latr=latr(x,y);
latr

% Read in temperature data
ft = 'ROMS/zw_tools/temp_spinup_MY25.mat';
tempfile = fullfile(f1,ft);
m = matfile(tempfile);
temp = squeeze(m.temp_spinup_MY25(x,y,:,:));
lt = length(temp(1,:));
if APRIL_START == 1
    % April Start
    ttmp = temp;
    temp = zeros(size(ttmp));
    for z = 1:N
        temp(z,1:ceil(lt/6)) = ttmp(z,ceil(5*lt/6):end);
        temp(z,ceil(lt/6)+1:end) = ttmp(z,1:floor(5*lt/6));
    end
end
% Interpolate 5-day temp data to timestep
ns5dy = 60*60*24*5;  % number of seconds in a 5-day period
temp_dt(time_loops+1,N) = zeros;
for z=1:N
    temp_dt(:,z) = interp1((0:ns5dy:lt*ns5dy-1),temp(z,:),0:dt:num_days*86400);
end

%Read in shortwave data
f2 = 'ROMS/Project_Patagonia/bulk';
fradsw = 'roms_blk_ERA_icecov_RHO_grid_3day_1900.nc';
swfile = fullfile(f1,f2,fradsw);
 % Modify forcing array to start from appropriate month
tmp = ncread(swfile,'radsw');
tmp = squeeze(tmp(x,y,:));  % On shelf
ll = length(tmp);
radsw(ll) = zeros;
if APRIL_START == 1
    % April Start
    radsw(1:ceil(2*ll/3)) = tmp(ceil(ll/3):end);
    radsw(ceil(2*ll/3)+1:end) = tmp(1:floor(ll/3));
else
    % July start
    radsw(1:ll/2) = tmp(ll/2+1:end);
    radsw(ll/2+1:end) = tmp(1:ll/2);
end
% Interpolate to hourly shortwave radiation
%srflx = interp1((0:3600:ll*3600-1),radsw,0:dt:num_days*86400)
% Interpolate 3-day shortwave radiation to time step
ns3dy = 60*60*24*3;  % number of seconds in a 3-day period
srflx = interp1((0:ns3dy:ll*ns3dy-1),radsw,0:dt:num_days*86400);

if DUST == 1
    % Read in iron dust data
    f3 = 'ROMS/Project_Patagonia/spinup';
    fdust = 'roms_frcbio_M2005.nc.1';
    dustfile = fullfile(f1,f3,fdust);
    % Modify forcing array to start from appropriate month
    tmp = ncread(dustfile,'dust');
    tmp = squeeze(tmp(x,y,:));
    ld = length(tmp);  
    dust_dep(ld) = zeros;
    if APRIL_START == 1
        % April Start
       dust_dep(1:2*ld/3) = tmp(ld/3+1:end);
       dust_dep(2*ld/3+1:end) = tmp(1:ld/3);
    else
        % July Start
        dust_dep(1:ld/2) = tmp(ld/2+1:end);
        dust_dep(ld/2+1:end) = tmp(1:ld/2);
    end
    % interpolate monthly dust deposition to timestep
    ns1mo = 60*60*24*(num_days/12);  % [seconds per month]
    Fe_dust = interp1((0:ns1mo:ld*ns1mo-1),dust_dep,0:dt:86400*num_days);
end

%
% Compute hour angle at sunset.  If current hour angle is less than the
% sunset hour angle, then turn off grazing since the model is in daytime.
% Current and sunset hour angles must be in radians and must be absolute
% values.
%
% Equations for current hour and declination angles are taken from ROMS
% kernel (ana_srflx.h) and sunset hour angle is compunted from NOAA Global
% Monitoring Division General Solar Position Calculations document:
%
%: https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF)
%
if ZOOL_GRAZ_NIGHT == 1
    % Compute solar declination and zenith angles
    Dangle = deg2rad(23.44*cos((172-(1:num_days))*2*pi/num_days));  % ROMS
    Zangle = deg2rad(90.833);  % global average
    latr = deg2rad(latr);  % latitude
    
    % Compute hour angle at sunset for each day of the year
    Hangle_sunset = acos(cos(Zangle)./(cos(latr)*cos(Dangle)) - ...
                                         tan(latr)*tan(Dangle));
    % Extend daily sunset angles so they are represented at all time steps.
    Hangle_sunset = repelem(Hangle_sunset,time_loops/num_days);
    
    % pre-allocated hour angle array for speed
    hr_frac_arr = zeros(1,time_loops);
    
    % compute fractional hour
    hr_frac = dt/3600:dt/3600:24;  % [d-1]   
    
    % extend fractional hour so daily cycle is represented at all time steps.       
    for i=1:time_loops/num_days:time_loops
        hr_frac_arr(i:i+time_loops/num_days-1) = hr_frac;
    end
    
   % Compute hour angle (radians) at each fractional hour and extend for days of run.
   Hangle = abs(12.0-hr_frac_arr)*pi/12;  % ROMS 
end

z = -h:h/N:-h/N;
hz = ones(1,N)*h/N;

ino3 = 1;        
inh4 = 2;        
idon = 3;        
idetn = 4;       
ibac = 5;        
iphy = 6;        
izoos = 7;       
ichl = 8;       
idetc = 9;       
idoc = 10;       
izool = 11;     
isil = 12;       
idia = 13;      
ildetn = 14;    
ildetc = 15;   
ildets = 16;    
ichld = 17;
if IRON == 1
    ifed = 18;
    idetfe = 19;
    ildetfe = 20;
end
if COCCOLITHOPHORES == 1
    iCocc = 21;
    iChlC = 22;
    iPIC = 23;
    iDIC = 24;
    iALK = 25;
end
%
% Initial Conditions
%
bio = zeros(time_loops,N,nbt);
% Read in initial conditions for appropriate variables
fbio_ini = fullfile(f1,f3,'roms_ini_MY25_graz0.01.nc');
bio(1,:,ino3) = squeeze(ncread(fbio_ini,'NO3',[x,y,1,1],[1,1,N,1]));
bio(1,:,inh4) = bio(1,:,ino3).*0.1;  % f-ratio 0.8
bio(1,:,idon) = 0.15d0;
bio(1,:,idetn) = 1.0d0;
bio(1,:,ibac) = 0.1d0;
bio(1,:,iphy) = 0.125d0;  % must not be equal to exud_min, or get divison by 0.
if GRAZ_INI_LOW == 1
    % Grazers 1% of producers
    bio(1,:,izoos) = bio(1,:,iphy)*.01;
else
    % Grazers 10% of producers
    bio(1,:,izoos) =  bio(1,:,iphy)*.05;
end
bio(1,:,ichl) = bio(1,:,iphy)*1.59;
bio(1,:,idetc) = bio(1,:,idetn)*cnd;

bio(1,:,idoc) = bio(1,:,idon)*cnd;
if GRAZ_INI_LOW == 1
    % Grazers 1% of producers
    bio(1,:,izool) = bio(1,:,iphy)*.01;
else
    % Grazers 10% of producers
    bio(1,:,izool) = bio(1,:,iphy)*.05;
end
bio(1,:,isil) = bio(1,:,ino3);%squeeze(ncread(fbio_ini,'Silica',[x,y,1,1],[1,1,N,1]));
bio(1,:,idia) = 0.125d0;  % must not be equal to exud_min, or get divison by 0.
bio(1,:,ildetn) = 0.02d0;
bio(1,:,ildetc) = bio(1,:,ildetn)*cnd;
bio(1,:,ildets) = bio(1,:,ildetn);
bio(1,:,ichld) = bio(1,:,idia)*1.59;
if IRON == 1
    %bio(1,:,ifed) = 1.5d0;
    bio(1,:,ifed) = squeeze(ncread(fbio_ini,'FeD',[x,y,1,1],[1,1,N,1]));
    bio(1,:,idetfe) = 0.01d0;
    bio(1,:,ildetfe) = 0.01d0;
end 
if COCCOLITHOPHORES == 1
    bio(1,:,iCocc) = bio(1,:,iphy);
    bio(1,:,iChlC) = bio(1,:,ichl);
    bio(1,:,iPIC) = bio(1,:,iCocc)*cncoc*PIC2C;
    bio(1,:,iDIC) = 2.2d3;
    bio(1,:,iALK) = bio(1,:,iDIC);
end
%
%  Set vertical sinking velocity vector in the same order as the
%  identification vector, idsink.
%
 vsink = 1.35d0;
 vsinkld = 50.0d0;
% vsink = 0;
% vsinkld=0;
wbio(1)=vsink;
wbio(2)=vsink;
wbio(3)=vsinkld;
wbio(4)=vsinkld;
%wbio(5)=vsinkld;
wbio(5)=10.d0;  % LDetS
if IRON == 1
    wbio(6) = vsink;
    wbio(7) = vsinkld;
end
if COCCOLITHOPHORES == 1
    wbio(8) = vsinkld;
end
%
%  Set vertical sinking indentification vector.
%
idsink(1)=idetn;
idsink(2)=idetc;
idsink(3)=ildetn;
idsink(4)=ildetc;
idsink(5)=ildets;
if IRON == 1
    idsink(6) = idetfe;
    idsink(7) = ildetfe;
end
if COCCOLITHOPHORES == 1
    idsink(8) = iPIC;
end
%
%  Compute inverse thickness to avoid repeated divisions.
%
hz_inv = zeros(1,N);
hz_inv2 = zeros(1,N);
hz_inv3 = zeros(1,N);
for k=1:N
        hz_inv(k)=1.0./hz(k);
end
for k=1:N-1
        hz_inv2(k)=1.0./(hz(k)+hz(k+1));
end 
for k=2:N-1
        hz_inv3(k)=1.0./(hz(k-1)+hz(k)+hz(k+1));
end

% Start time stepping
for t = 1:time_loops-1
    if DUST == 1
        % Atmospheric iron deposition
        bio(t+1,N,ifed) = bio(t+1,N,ifed) + ...      
                                    (Fe_solub/(MW_Fe*365.25 * ... 
                                    86400/12)*Fe_dust(t) * ...  
                                    dt*hz_inv(N)*1.d9);    
        bio(t+1,N,isil)= bio(t+1,N,isil) + ...
                                 (8.8*Si_solub/(MW_Si*365.25 * ...
                                 86400/12)*Fe_dust(t) * ...
                                 dt*hz_inv(N)*1.d6);        
    end
%
% Start looping in depth
for k=1:N
%
%-----------------------------------------------------------------------
% Compute biology forcing for current time step.
%-----------------------------------------------------------------------
%
no3   = bio(t,k,ino3); 
nh4   = bio(t,k,inh4); 
don   = bio(t,k,idon); 
detn  = bio(t,k,idetn); 
bact  = bio(t,k,ibac); 
phyto = bio(t,k,iphy); 
zoo   = bio(t,k,izoos); 
chlor = bio(t,k,ichl); 
detc  = bio(t,k,idetc); 
doc   = bio(t,k,idoc); 
zool  = bio(t,k,izool); 
sil   = bio(t,k,isil); 
dia   = bio(t,k,idia); 
ldetn = bio(t,k,ildetn); 
ldetc = bio(t,k,ildetc); 
ldets = bio(t,k,ildets); 
chld  = bio(t,k,ichld); 
if IRON == 1
    fed = bio(t,k,ifed); 
    %fed =max(bio(t,k,ifed), 0.3);  % Model doesn't really drop below this value
    detfe = bio(t,k,idetfe); 
    ldetfe = bio(t,k,ildetfe);
    
    if FEC_VAR == 1
        fec_phy = max(0.002,SB_Fe*bio(t,k,ifed)^SA_Fe/1000);
        fec_phy_diag(t,k) = fec_phy;
        fec_dia = max(0.002,LB_Fe*bio(t,k,ifed)^LA_Fe/1000);
        fec_dia_diag(t,k) = fec_dia;
        fec_coc = max(0.002,CB_Fe*bio(t,k,ifed)^CA_Fe/1000);
    end    
end
if COCCOLITHOPHORES == 1
    cocco = bio(t,k,iCocc);
    chlc = bio(t,k,iChlC);
    PIC= bio(t,k,iPIC);
    DIC = bio(t,k,iDIC);
    ALK = bio(t,k,iALK);
end

% 
%-----------------------------------------------------------------------
% Compute specific cellular chlorophyll-a:carbon ratio (theta,
% mg Chla:mg C).
%-----------------------------------------------------------------------
%
% Small Phytoplankton
%
denom = phyto .* cnp .* 12.0d0;
theta = max(chlor./denom, 0.0d0);
%theta = min(theta,mx_chlcr);
%
% Diatoms
%
denomd = dia .* cndia .* 12.0d0;
thetad = max(chld./denomd, 0.0d0);
%thetad = min(thetad,mx_chlcrd);
%
% Coccolithophores
%
%thetac = 0.026;  % Assume Chl:N of 1.59 with a C:N ratio of 5 
if COCCOLITHOPHORES == 1
    denomc = cocco * cncoc *12.d0;
    thetac = max(chlc/denomc,0.d0);
end

%  Eppley 1974 Temperature modification of phytoplankton growth
%  which gives Vpt =3.26 at 21C and 5.09 at 28C
%         Vpt = 0.851*(1.066**temp)
%         Vptdia = 0.851*(1.066**temp)
%
% Formulation for growth rate taken from Boyd et al. 2013.
%
Vpt = acof*exp(bcof*temp_dt(t,k))*(1-((temp_dt(t,k)-topt_phy)/(w_phy/2))^2);
Vptdia = acof_dia*exp(bcof*temp_dt(t,k))*(1-((temp_dt(t,k)-topt_dia)/(w_dia/2))^2);
Vptc = acof*exp(bcof*temp_dt(t,k))*(1-((temp_dt(t,k)-topt_coc)/(w_coc/2))^2);
%Vptc = 1.5;
%
% Formulation for growth rate from Yvette 1/11/2019
%       vpt = vp.*exp(0.0633.*temp)./(1.0+exp(temp.*2.0-9.0));
%       vptdia = vpdia.*exp(0.0633.*temp)./(1.0+exp(temp.*2.0-9.0));
%         Vpt = Vptmax*exp(bcof*temp_dt(t,k))*exp(-(temp_dt(t,k)-Topt_phy)^2/2.5^2);
%         Vptdia = Vptdiamax*exp(bcof*temp_dt(t,k))*exp(-(temp_dt(t,k)-Topt_dia)^2/2.5^2);
%
% Maintain backgroundtit population by setting growth rate to 0.01 [d-1]
% where it would otherwise be zero.
%
Vpt = max(1.d-4,Vpt);
Vptdia = max(1.d-4,Vptdia);
Vptc = max(1.d-4, Vptc);
%
% Prevent growth rates from exceeding prescribed maxima
%
% vpt = min(vpt,vp);
% vptdia = min(vptdia,vpdia);
%
%-----------------------------------------------------------------------
% Calculate the nutrient limitation terms on small 
%  phytoplankton, diatoms and, if necessary, 
% coccolithophores.
%-----------------------------------------------------------------------
%
% Small Phytoplankton
%
ak = kp_no3/kp_nh4;
deno = kp_no3+no3+nh4.*ak;
q1 = no3.*exp(-psi.*nh4) ./ deno;   % nondimensional NO3 limitation factors
q2 =(nh4.*ak) ./ deno;              % nondimensional NH4 limitation factors
if IRON == 1
    q3 = fed / (k1fe + fed);        % nondimensional FeD limitation
    if LIEBIG == 1
        q  = min(q1 + q2, q3);
    elseif ONEILL == 1
        DIN = (no3+nh4);
        num = DIN*fed;
        denom = (kp_no3*fed) + (k1fe*DIN) + (DIN*fed);
        q = num / denom;
    else
        q = q3*(q1+q2);
    end
else    
    q = q1+q2;
end
if( q == 0.0 )
    fvalue = 0.0;
else
    fvalue = q1./(q1+q2);
end
%
% Diatoms
%
akd = kdia_no3./kdia_nh4;
denod = kdia_no3+no3+nh4.*akd;
qd1 = no3.*exp(-psi.*nh4) ./ denod; % nondimensional NO3 limitation factors
qd2 =(nh4.*akd) ./ denod;           % nondimensional NH4 limitation factors
qd3 = sil ./(sil+kdia_si);          % nondimensional SIL limitation factors
if IRON == 1
    qd4 = fed / (k2fe+fed);
    if LIEBIG == 1
        qd = min([qd1+qd2; qd3; qd4]);    % Law of the minimum -- Z.W. 5.17.18
    elseif ONEILL == 1
        num = DIN*fed*sil;
        denom = (kdia_no3*fed*sil) + (kdia_si*fed*DIN) + (k2fe*DIN*sil) + DIN*fed*sil;
        denom = max(denom,1d-4);  % prevent division by 0
        qd = num/denom;
    else
        qd = qd4*min(qd1+qd2,qd3);
    end
else
    qd = min(qd1+qd2, qd3);
end
if( qd == 0.0 )
    fvalued = 0.0;
else
    fvalued = qd1./(qd1+qd2);
end
if IRON == 1
% Si:N varies between 1:1 (high iron) and 3:1 (low iron)
    lpsin = 2*exp(-(fed-0.1)*5)+1;
    %lpsin = 1;
end
if COCCOLITHOPHORES == 1
    %
    % Coccolithophores
    %
    akc = kc_no3/kc_nh4;
    denoc = kc_no3+no3+nh4.*akc;
    q1c = no3.*exp(-psi.*nh4) ./ denoc;   % nondimensional NO3 limitation factors
    q2c =(nh4.*akc) ./ denoc;              % nondimensional NH4 limitation factors
    if IRON == 1
        q3c = fed / (kcfe + fed);        % nondimensional FeD limitation
        if LIEBIG == 1
            qcoc  = min(q1c + q2c, q3c);
        elseif ONEILL == 1        
            num = DIN*fed;
            denom = (kc_no3*fed) + (kcfe*DIN) + (DIN*fed);
            qcoc = num / denom;
        else
            qcoc = q3c*(q1c+q2c);
        end
    else    
        qcoc = q1c+q2c;
    end
    if( qcoc == 0.0 )
        fvaluec = 0.0;
    else
        fvaluec = q1c./(q1c+q2c);
    end
end
%
%  Calculate Surface Photosynthetically Available Radiation (PAR).  The
%  net shortwave radiation is scaled back to Watts/m2 and multiplied by
%  the fraction that is photosynthetically available, PARfrac.
%
    parsur=bioparfrac*srflx(t);
    %parsur=bioparfrac*0.8*srflx(t);  % Account for PUR
    %parsur=0.1*bioparfrac*srflx(t);  % 10% light level.

if (parsur > 1.0d-6)
%
% Reduce diatom growth in waters shallower than 20m.  This is
% to prevent diatoms from congregating near shallow river mouths
% (e.g. the La Plata) and preventing nutrient flow onto the shelf.
% Kd(PAR) estimated from MODIS-Aqua monthly climatologies as well
% as from the manuscript:
%    Wang et al. 2009, Retrieval of diffuse attenuation coefficient
%    in the Chesapeake Bay and turbid ocean regions for satellite
%    ocean color applications. JGR. 114.
% Added -- Z.W. 8 May 2018.
%
if (h <= 20.0d0) 
    kdpar = 1.5d0;
elseif (h >20.d0 && h < 200.d0)
    kdpar = attsw*2;    
else
    kdpar = attsw;
end
%
% Compute biological light attenuation as a function of chlorophyll.
%
if CHL_PAR_ATT_SUM == 1
    % Attenuate light by the sum of chl at all levels above depth z(k)
    chl_tot = 0;
    for vert_level = k:N
        if COCCOLITHOPHORES == 1
        chl_tot = chl_tot - hz(vert_level) * ...
                       (attphyt*bio(t,vert_level,ichl) + ...
                       attphytdia*bio(t,vert_level,ichld) + ...
                       attphytcoc*bio(t,vert_level,iChlC));
        else
        chl_tot = chl_tot - hz(vert_level) * ...
                       (attphyt*bio(t,vert_level,ichl) + ...
                       attphytdia*bio(t,vert_level,ichld));            
        end
    end
    chl_tot = (1/k)*chl_tot;
    parz_diag(t,k) = chl_tot;
    parz = parsur*exp(z(k)*kdpar + chl_tot);
else
    parz = parsur*exp(z(k).*(kdpar +(attphyt.*chlor+attphytdia.*chld)));
    parz_diag(t,k) = (attphyt.*chlor+attphytdia.*chld)*z(k);
end

if (parz < 1.0d-6)
    parz = 0.0;
end

if FEC_VAR_LIGHT == 1
    % Vary Fe:C with parz (high Fe:C at low parz)
    fec_phy = (20+180*parz^-0.4)/10000;
    fec_phy = min(fec_phy,250d-4);
    fec_phy_diag(t,k) = fec_phy;
    fec_dia = (20+180*parz^-0.4)/1000;
    fec_dia = min(fec_dia, 250d-3);
    fec_dia_diag(t,k) = fec_dia;
end

% Modulate photosynthetic efficiency by ambient Fe concentration
if ALPHA_LIM_FE == 1
    alphaN = alpha * q3;
    alphadiaN = alphadia *qd4;
    aj = Vpt.*theta.*alphaN.*parz./sqrt((Vpt.*Vpt) +(theta.*theta.*alphaN.*alphaN.*parz.*parz));
    ajd = Vptdia.*thetad.*alphadiaN.*parz./sqrt((Vptdia.*Vptdia)+(thetad.*thetad.*alphadiaN.*alphadiaN.*parz.*parz));
else
    aj = Vpt.*theta.*alpha.*parz./sqrt(Vpt.*Vpt +theta.*theta.*alpha.*alpha.*parz.*parz);
    ajd = Vptdia.*thetad.*alphadia.*parz./sqrt(Vptdia.*Vptdia+thetad.*thetad.*alphadia.*alphadia.*parz.*parz);
    if COCCOLITHOPHORES == 1
        ajc = Vptc*thetac*alphac*parz / ...
                  sqrt(Vptc*Vptc+thetac*thetac*alphac*alphac*parz*parz);
    end
end
aj_diag(t,k) = aj;
ajd_diag(t,k) = ajd;

% Check for low light.  If light is low, turn on nitrification.
if(parz >= 0.1 * parsur)
    anitrf = 0.;
    %anitrf = anitrof_l;
else
    anitrf = anitrof;
    %anitrf = anitrof_d;
end
% PAR0 is zero, no light due to earth rotation or earth seaonal variation.
else
    parz  = 0.0;
    aj    = 0.0;
    ajd   = 0.0;
    anitrf = anitrof;   %no light, turn on nitrification
    %anitrf = anitrof_d;
end

% Temp, nutrient, and light-limited growth terms
sigma=aj.*q;
sigmad=ajd.*qd;
if COCCOLITHOPHORES == 1
    sigmac=ajc*qcoc;
end

if COCCOLITHOPHORES ==1
    if (parz<=0. || theta<=0 || thetad<=0 || thetac<=0)
        rho_c =0.;
        rhod_c =0.;
        rhoc_c = 0.d0;
    else
        rho_c = mx_chlcr.*sigma./(alpha.*theta.*parz);
        rhod_c = mx_chlcrd.*sigmad./(alphadia.*thetad.*parz);
        rhoc_c = mx_chlcrc*sigmac/(alphac*thetac*parz);
    end
else
    if (parz<=0. || theta<=0 || thetad<=0)
        rho_c =0.;
        rhod_c =0.;
    else
        rho_c = mx_chlcr.*sigma./(alpha.*theta.*parz);
        rhod_c = mx_chlcrd.*sigmad./(alphadia.*thetad.*parz);
%     rho_c = mx_chlcr.*sigma./(.01*4*theta*parz);
%     rhod_c = mx_chlcrd.*sigmad./(.01*4*thetad*parz);
    end    
end

%
%-----------------------------------------------------------------------
% Calculate small and large zooplankton grazing terms.
%-----------------------------------------------------------------------
%   coeff. and formulation from BATS paper.
if COCCOLITHOPHORES ==1
    value=kz_i.*(phyto+prefz_b*bact+prefz_d *detn+prefz_c*cocco) + ...
               phyto*phyto + ...
               prefz_b*bact*bact + ...
               prefz_d*detn*detn + ...
               prefz_c*cocco*cocco;    
else
    value=kz_i.*(phyto+prefz_b*bact+prefz_d *detn) + ...
               phyto*phyto + ...
               prefz_b *bact*bact + ...
               prefz_d .*detn.*detn;
end
if(value > 0.0d0)
    g1=(gmax.*zoo.*phyto.*phyto)./value;
    g2=(gmax.*zoo.*prefz_b .*bact.*bact)./value;
    g3=(gmax.*zoo.*prefz_d .*detn.*detn)./value;
    if COCCOLITHOPHORES == 1
        GS_coc = (gmax*zoo*prefz_c*cocco*cocco)/value;
    end
else
    g1=0.0d0;
    g2=0.0d0;
    g3=0.0d0;
    if COCCOLITHOPHORES == 1
        GS_coc = 0.d0;
    end
end
%  Include mesozooplankton grazing on phy, dia, and microzoo.
%  Z.W. 11.16.18
%
if COCCOLITHOPHORES == 1
    valuel=kl_i.*(dia+prefl_z.*zoo+prefl_p.*phyto+prefl_c*cocco) + ... 
                dia.*dia + ... 
                prefl_z.*zoo.*zoo + ... 
                prefl_p.*phyto.*phyto + ...
                prefl_c*cocco*cocco;
else
    valuel=kl_i.*(dia+prefl_z.*zoo+prefl_p.*phyto) + ... 
                dia.*dia + ... 
                prefl_z.*zoo.*zoo + ... 
                prefl_p.*phyto.*phyto;
end 
if(valuel > 0.d0)
    g4=grazl .* zool .* dia.*dia ./ valuel;
    g5=grazl .* zool .* prefl_z .* zoo.*zoo ./ valuel;
    g6=grazl .* zool .* prefl_p .* phyto.*phyto ./ valuel;
    if COCCOLITHOPHORES == 1
        GL_coc = grazl*zool*prefl_c*cocco*cocco/valuel;
    end
else
    g4 = 0.0d0;
    g5 = 0.0d0;
    g6 = 0.0d0;
    if COCCOLITHOPHORES == 1
        GL_coc = 0.d0;
    end
end

if ZOOL_GRAZ_NIGHT == 1
   % Determine time of day to see if mesozooplankton are feeding.
   if Hangle(t) < Hangle_sunset(t)  % Day time
       g4 = 0.d0;
       g5 = 0.d0;
       g6 = 0.d0;
       if COCCOLITHOPHORES == 1
           GL_coc = 0.d0;
       end
   end
end

%
%-----------------------------------------------------------------------
% Calculate the bacteria terms.
%-----------------------------------------------------------------------
%
docn = don.*cnd;    %DOC_nitro
docc = max(1.0d-12, doc-docn);  %DOC_carbo
%docc = doc-docn;  %DOC_carbo
u_docc  = docc./(kb_c+docc);    % UDOC
u_don   = vb_amino.*bact.*don./(kb_don+don);    % UDON
u_nitro = u_don.*(1.0d0-(cnd .*ggen ./cnb ));   % UNH4_nitro
u_nh4 =(vb_amino.*bact-u_don).*(nh4./(kb_a+nh4)).* u_docc;
u_carbo =(u_docc.*u_nitro+u_nh4).*cnb ./ggec;
regen   =(1.0d0-u_docc).*u_nitro;
u_nh4   = u_nh4-regen;
u_doc   = u_carbo+u_don.*cnd;
%
%     conversion factors and constants
%
convc = 12.0d0*cnp;
convcd = 12.0d0*cndia;
convc_theta = theta*convc;
convcd_thetad = thetad*convcd;
if COCCOLITHOPHORES == 1
    convcc = 12.d0*cncoc;
    convcc_thetac = thetac*convcc;
end
%
%-----------------------------------------------------------------------
% Calculate the nitrate terms.
%-----------------------------------------------------------------------
%
%        fl1_6 = aj*Q1*phyto            % NO3 uptake by phyto
fl1_6 = aj.*fvalue.*q.*phyto;           % NO3 uptake by phyto Z.W. 5.31.18
fl1_8 = fl1_6.*rho_c.*convc;            % NO3 uptake by chlorophyll
%        fl1_13 = ajd*Qd1*Qd3*dia       % NO3 uptake by diatoms
fl1_13 = ajd.*fvalued.*qd.*dia;         % NO3 uptake by diatoms Z.W. 5.31.18
fl1_17 = fl1_13.*rhod_c.*convcd;        % NO3 uptake by diatom chlorophyll
if COCCOLITHOPHORES == 1
    fl1_coc = ajc*fvaluec*qcoc*cocco;
    fl1_chlc = fl1_coc*rhoc_c*convcc;
end
%
%-----------------------------------------------------------------------
% Calculate the ammonium terms.
%-----------------------------------------------------------------------
%
fl2_1 = anitrf .* nh4;                  % nitrifaction at low light 11/17/2000 yvette                                        
fl2_5 = u_nh4;                          % NH4 uptake and regeneration by bact.                                     
%        fl2_6 = aj*Q2*phyto            % NH4 uptake by phyto
fl2_6 = aj.*(1-fvalue).*q.*phyto;       % NH4 uptake by phyto Z.W. 5.31.18
fl2_8 = fl2_6.*rho_c.*convc;            % NH4 uptake by chlorophyll
%        fl2_13 = ajd*Qd2*Qd3*dia       % NH4 uptake by diatom
fl2_13 = ajd.*(1-fvalued).*qd.*dia;     % NH4 uptake by diatom Z.W. 5.31.18
fl2_17 = fl2_13.*rhod_c.*convcd;        % NH4 uptake by diatom chlorophyll
if COCCOLITHOPHORES == 1
    fl2_coc = ajc*(1-fvaluec)*qcoc*cocco;
    fl2_chlc = fl2_coc*rhoc_c*convcc;
end
%                                         
%-----------------------------------------------------------------------
% Calculate the DON terms.
%-----------------------------------------------------------------------
%
fl3_5 = u_don;                          % uptake of DON by bact.
%
%-----------------------------------------------------------------------
% Calculate the detrital nitrogen terms.
%-----------------------------------------------------------------------
%
% Breakdown to DON
%
fl4_3 =(mud_dom+mub_dom*bact).*detn;
%
% Ingestion by micro/nanozooplankton
%
fl4_7 = g3;
%
%-----------------------------------------------------------------------
% Calculate the bacteria terms.
%-----------------------------------------------------------------------
%
fl5_2 = mub_a * bact;      % bact excretion to NH4
if COCCOLITHOPHORES == 1
    fl5_dic = fl5_2 *cnb;  % bact excretion to DIC
    fl5_alk = fl5_dic;
end
fl5_4 = mub_d * bact;      % bact mortality to DETN
fl5_7 = g2 * betaz_b;      % bact ingestion by zoo
fl5_9 = fl5_4 * cnb;       % bact mortality to DETC
if IRON_BACTERIA == 1
    fl5_18 = fl5_2 * cnb * fec_bac;     % bact excretion to FeD
    fl5_19 = fl5_9 * fec_bac;           % bact mortality to DetFe
end
%
%-----------------------------------------------------------------------
% Calculate the Phytoplankton terms.
%    3.12.18  --> Quadratic mortality
%    5.25.18  --> Minimum respiration threshold
%-----------------------------------------------------------------------
%
if EXUD_MINVAL == 1
% Turn off exudation to DON if biomass is below EXUD_MIN.
%    1.24.19 --> minimum exudation applied to gamma2
%
    phy_exud1 = phyto - exud_min;
    phy_exud2 = max(phy_exud1,0.d0)*gamma1p/phy_exud1;
    phy_exud3 = max(phy_exud1,0.d0)*gamma2p/phy_exud1;
    fl6_3 = phy_exud2*phyto + phy_exud3*sigma*phyto;
else
    % Constant exudation to DON
    fl6_3 = gamma1p * phyto + gamma2p * sigma * phyto;
end

if QUADRATIC_MORT == 1
    % Quadratic mortality to DetN
    fl6_4 = min(phy_mort_max, phy_mort_max * ...
                      (phyto*phyto/(Km_phy+phyto*phyto))) * ...
                       phyto;  
else
    % Linear mortality to DetN
    fl6_4 = mup_d * phyto;
end

% Ingestion by zooplankton
fl6_7 = g1;     % Micro/nano
fl6_11 = g6;    % Meso
%
% Mortality to DetC
%
fl6_9 = fl6_4 * cnp;
%
% Exudation to DOC
%
fl6_10 = fl6_3 * cnp + gamma3p * max(aj./Vpt-q,0.0d0) * aj * phyto * cnp;
%
if IRON == 1
% Exudation to FeD
%
    fl6_18 = fl6_3 * cnp * fec_phy;
%
% Mortality to DetFe
%
    fl6_19 = fl6_9 * fec_phy;
end    
%
%-----------------------------------------------------------------------
% Calculate the microzooplankton terms.
%-----------------------------------------------------------------------
%
detc2detn = max(detc./detn,0.0d0);

if IRON == 1
    detfe2detc = max(detfe/detc,0.0d0);
end
% mortality rate
mort_zoo = muz_a*zoo;
%
% Excretion to NH4
%
if COCCOLITHOPHORES ==1
    fl7_2  = epsilonz_a * mort_zoo +(betaz_p*g1+betaz_b*g2+betaz_d*g3+betaz_c*GS_coc) * muzg_a;
    % Excretion to DIC
    fl7_dic = fl7_2 * cnz;
    fl7_alk = fl7_dic;    
else       
    fl7_2  = epsilonz_a * mort_zoo +(betaz_p * g1 +betaz_b * g2 +betaz_d * g3) * muzg_a;
end

% Mortality to DON
%
fl7_3  = omegaz_d * mort_zoo;
%
% Growth-specific respiration and mortality to detrital pools
%
if COCCOLITHOPHORES == 1
    fl7_4  =(1.0d0-betaz_p) * g1 +(1.0d0-betaz_b) * g2 +(1.0d0-betaz_d) * g3 + (1.d0-betaz_c)*GS_coc;
    fl7_9  =(1.0d0-betaz_p) * g1 * cnp  +(1.0d0-betaz_b) * g2 * cnb  +(1.0d0-betaz_d) * g3 * detc2detn + (1.d0-betal_c)*GS_coc;
else
    fl7_4  =(1.0d0-betaz_p) * g1 +(1.0d0-betaz_b) * g2 +(1.0d0-betaz_d) * g3;
    fl7_9  =(1.0d0-betaz_p) * g1 * cnp  +(1.0d0-betaz_b) * g2 * cnb  +(1.0d0-betaz_d) * g3 * detc2detn;
end
fl7_4  = fl7_4 +(1.d0 - epsilonz_a - omegaz_d) * mort_zoo;
fl7_9  = fl7_9 +(1.d0 - epsilonz_a - omegaz_d) * mort_zoo * cnz;
if COCCOLITHOPHORES == 1
    fl7_pic = fl7_9*PIC2C;
end
%
% Mortality to DOC
%
fl7_10 = fl7_3 .* cnz;
%
% Grazing by mesozooplankton
%
fl7_11 = g5;
%
if IRON == 1 
    if IRON_BACTERIA == 1
% Growth-specific respiration to FeD.
% Fe acquired from bacteria, small phytoplankton and detritus.
%
        fl7_18 = fl7_2 * cnz * fec_zoo;
    else
% Growth-specific respiration to FeD.
% Fe acquired from small phytoplankton and detritus only, not bacteria.
%
        fl7_18 = (fl7_2-betaz_b*g2*muzg_a) * cnz * fec_zoo;
    end
%
% Add in fluxes to FeD from micro/nano-zooplankton-specific mortality.
%
        fl7_18 = fl7_18 + fl7_10*fec_zoo;
%        
    if IRON_BACTERIA == 1
%
% Growth-specific excretion to DetFe (sloppy feeding).
% Fe acquired from bacteria, small phytoplankton, and detritus.
%
        if COCCOLITHOPHORES == 1
            fl7_19 = (1.d0-betaz_p) * g1 * cnp * fec_phy + (1.d0-betaz_b) * ...
                          g2 * cnb * fec_bac + (1.d0-betaz_d) * g3 * detc2detn * detfe2detc + ...
                          (1.d0-betaz_c)*GS_coc*cncoc*fec_coc;
        else
            fl7_19 = (1.d0-betaz_p) * g1 * cnp * fec_phy + (1.d0-betaz_b) * ...            
                           g2 * cnb * fec_bac + (1.d0-betaz_d) * g3 * detc2detn * detfe2detc;
        end
    else
%
% Growth-specific excretion to DetFe (sloppy feeding).
% Fe acquired from small phytoplankton and detritus only, not bacteria.
%
        if COCCOLITHOPHORES == 1
            fl7_19 = (1.d0-betaz_p) * g1 * cnp * fec_phy + (1.d0-betaz_d) * g3 * detc2detn * detfe2detc + ...
                           (1.d0-betaz_c) * GS_coc * cncoc *fec_coc;
        else
            fl7_19 = (1.d0-betaz_p) * g1 * cnp * fec_phy + (1.d0-betaz_d) * g3 * detc2detn * detfe2detc;
        end
    end
%
% Add in fluxes to DetFe of contributions from growth-specific
% excretion and mortality.
        fl7_19 = fl7_19 + (1 - omegaz_d - epsilonz_a) * muz_a * zoo * cnz * fec_zoo;
end
    
%
%-----------------------------------------------------------------------
% Calculate the chlorophyll terms (small phytoplankton).
%-----------------------------------------------------------------------
%
fl8_4  = fl6_4.*convc_theta;    % chl mortality to DETN
fl8_7  = fl6_7.*convc_theta;    % chl ingestion by microzoo
fl8_11 = fl6_11.*convc_theta;   % chl ingestion by mesozoo
if K_LOSS == 1
    if K0 == 1
        fl8_dmg = chlor*chl_dmg;  % constant chl loss [d-1]
    elseif K1 == 1
        fl8_dmg = chlor*(aj/Vpt)*chls_dmg;  % variable loss with production
    elseif K2 == 1
        fl8_dmg = chlor*chl_dmg*theta*parz;  % variable loss with light and theta
    elseif K3 == 1
        fl8_dmg = chlor*chl_dmg*parz;  % variable loss with light intensity
    end
    fl8_dmg_diag(t,k) = fl8_dmg;
    fl8_4 = fl8_4 + fl8_dmg;  % stick on mort. so don't have to include extra term in integration.
end
%
%-----------------------------------------------------------------------
% Calculate the detrital carbon terms.
%-----------------------------------------------------------------------
%
% Ingestion by micro/nanozooplankton
%
fl9_7  = g3.*detc2detn;
%
% Breakdown to DOC
%
fl9_10 =(mud_dom+mub_dom*bact) .* detc;
%
%-----------------------------------------------------------------------
% Calculate the DOC terms.
%-----------------------------------------------------------------------
%
% DOC uptake by bact
fl10_5 = u_doc;
%
%-----------------------------------------------------------------------
% Calculate the Mesozooplankton terms.
%-----------------------------------------------------------------------
%
% mortality rate
if QUAD_MORT_ZOO == 1
    mort_zool = mul*zool*zool;
else
    mort_zool = mul * zool;
end
%
% Mortality and excretion to NH4
%
if COCCOLITHOPHORES ==1
    fl11_2  = epsilonl * mort_zool +(betal_dia*g4+betal_z*g5+betal_p*g6+betal_c*GL_coc) * mulg_a;
    fl11_dic = fl11_2 * cnz;
    fl11_alk = fl11_dic;
else
    fl11_2  = epsilonl * mort_zool +(betal_dia * g4 + betal_z * g5 +betal_p * g6) * mulg_a;
end
%
% Mortality to DON
%
fl11_3  = omegal * mort_zool;
%
% Mortality to DOC
%
fl11_10 = fl11_3 * cnz;
%
if IRON == 1
% Mortality and excretion to FeD
%
    fl11_18 = (fl11_2*cnz + fl11_10) * fec_zoo;
end
%
% Mortality and excretion to LDetN and LDetC (and PIC if applicable)
%
if COCCOLITHOPHORES == 1
    fl11_14=(1.d0 - epsilonl - omegal) * mort_zool + ...
                  (1.d0-betal_dia) * g4 +(1.d0-betal_z) * g5 +(1.d0-betal_p) * g6 + (1.d0-betal_c)*GL_coc;
    fl11_15=(1.d0 - epsilonl - omegal) * mort_zool * cnz + ... 
                  (1.d0-betal_dia) * g4 * cndia + (1.d0-betal_z) * g5 * cnz + (1-betal_p) * g6 * cnp + ...
                   (1.d0-betal_c)*GL_coc *cncoc*fec_coc;
    fl11_pic = fl11_15*PIC2C;
else
    fl11_14=(1.d0 - epsilonl - omegal) * mort_zool + ...
                  (1.d0-betal_dia) * g4 +(1-betal_z) * g5 +(1-betal_p) * g6;
    fl11_15=(1.d0 - epsilonl - omegal) * mort_zool * cnz + ... 
                  (1.d0-betal_dia) * g4 * cndia + (1.d0-betal_z) * g5 * cnz + (1-betal_p) * g6 * cnp;
end
%
% Mortality and excretion to LDetS
%
fl11_16 = (1.d0-epsilonl-omegal)*mort_zool*lpsin + (1-betal_dia)*g4*lpsin;
fl11_16_diag(t,k) = fl11_16;
if IRON == 1
%
% Mortality and excretion to LDetFe
%
if COCCOLITHOPHORES == 1
    fl11_20 = (1.d0 - epsilonl - omegal) * mort_zool*cnz*fec_zoo + ...                         
        (1-betal_dia)*g4*cndia*fec_dia+(1-betal_z)*g5*cnz*fec_zoo+(1-betal_p)*g6*cnp*fec_phy + ...
        (1.d0-betal_c)*GL_coc*cncoc*fec_coc;    
else
    fl11_20 = (1.d0 - epsilonl - omegal) * mort_zool*cnz*fec_zoo + ...                         
        (1-betal_dia)*g4*cndia*fec_dia+(1-betal_z)*g5*cnz*fec_zoo+(1-betal_p)*g6*cnp*fec_phy;
end
end
%
%-----------------------------------------------------------------------
% Calculate the silica terms terms.
% Variable Si:N added -- Z.W. 12.22.17
%-----------------------------------------------------------------------
%
% Fe-independent silica uptake by diatoms
%
if IRON == 1
% Silica uptake by diatom
%
% Fe-dependent Silica uptake by diatoms
%
    fl12_13 = ajd * qd * dia * lpsin;
else
%
% Fe-independent silica uptake by diatoms
%
    fl12_13 = ajd * qd * dia * snpd;
end
fl12_13_diag(t,k) = fl12_13;
%
%-----------------------------------------------------------------------
% Calculate the Diatom  terms.
%    12.22.17 --> Variable Si:N
%    3.12.18  --> Quadratic mortality
%    5.25.18  --> Minimum respiration threshold
%-----------------------------------------------------------------------
%
if EXUD_MINVAL == 1
% Turn off exudation to DON if biomass is below EXUD_MIN.
%    1.24.19 --> minimum exudation applied to gamma2
%
    dia_exud1 = dia - exud_min;
    dia_exud2 = max(dia_exud1,0.d0)*gamma1dia/dia_exud1;
    dia_exud3 = max(dia_exud1,0.d0)*gamma2dia/dia_exud1;
    fl13_3 = dia_exud2*dia + dia_exud3*sigmad*dia;
else
    % Constant exudation to DON
    fl13_3 = gamma1dia.*dia+gamma2dia.*sigmad.*dia;
end

%
% Exudation to DOC
%
fl13_10 = fl13_3.*cndia+gamma3dia * max(ajd./Vptdia-qd,0.0d0).*ajd .*dia.*cndia;
%
% Ingestion by mesozooplankton
%
fl13_11 = g4;

if QUADRATIC_MORT == 1
    % Quadratic mortality to LDetN
    fl13_14 = min(dia_mort_max, dia_mort_max * ...
                           (dia*dia/(Km_dia+dia*dia))) * ...
                           dia;  
else
    % Linear mortality to LDetN
    fl13_14 = amud1 * dia;
end

% Mortality to LDetC
fl13_15 = fl13_14 * cndia;
%
% Constant Si:N ratio (no Fe-dependence) for mortality to LDetS.
%
if IRON == 1
%-----------------------------------------------------------------------
% Variable Si:N ratio (Fe-dependent) for mortality to LDetS.  Mortality
% is modeled as mortality to LDetN (i.e. linear or quadratic).
%-----------------------------------------------------------------------
    fl13_16 = fl13_14 * lpsin;
else
%
% Constant Si:N ratio (no Fe-dependence) for mortality to LDetS.
%
    fl13_16 = fl13_14 * snpd;
end
if IRON == 1
%
% Exudation to FeD
%
        fl13_18 = fl13_3 * cndia * fec_dia;
%
% Mortality to LDetFe
%
        fl13_20 = fl13_15 * fec_dia;
end
%
%-----------------------------------------------------------------------
% Calculate the large detrital nitrogen terms.
%-----------------------------------------------------------------------
%
% LDetN breakdown to DON
fl14_3 = ldetndis.*ldetn;
%
%-----------------------------------------------------------------------
% Calculate the large detrital carbon terms.
%-----------------------------------------------------------------------
%
% LDetC breakdown to DOC
fl15_10 = ldetcdis.*ldetc;
%
%-----------------------------------------------------------------------
% Calculate the large detrital silica terms.
%-----------------------------------------------------------------------
%
% LDetS breakdown to sil
if SIL_REGEN_BAC == 1
    fl16_12 = (ldetsdis +  mub_dom*bact)*ldets;
else
    fl16_12 = ldetsdis.*ldets;
end
fl16_12_diag(t,k) = fl16_12;
%
%-----------------------------------------------------------------------
% Calculate the chlorophyll terms (diatoms).
%-----------------------------------------------------------------------
%
% chl ingestion by mesozoo
fl17_11 = fl13_11.*convcd.*thetad;
% chl mortality to LDETN
fl17_14 = fl13_14.*convcd.*thetad;
if K_LOSS == 1
    if K0 == 1
        fl17_dmg = chl_dmg*chld;  % constant chl loss [d-1]
    elseif K1 == 1
        fl17_dmg = chld_dmg*chld*(ajd/Vptdia);  % variable loss with production
    elseif K2 == 1
        fl17_dmg = chl_dmg*chld*thetad*parz;  % variable loss with light and theta
    elseif K3 == 1
        fl17_dmg = chl_dmg*chld*parz;  % variable chl loss with light
    end
    fl17_dmg_diag(t,k) =fl17_dmg;
    fl17_14 =fl17_14 + fl17_dmg;  % stick on mort. so don't have to include extra term in integration.
end

if IRON == 1
%
%-----------------------------------------------------------------------
% Calculate the dissolved iron terms.
%-----------------------------------------------------------------------
%
    if IRON_BACTERIA == 1
%
% Uptake by bacteria
%
        fl18_5 = (u_don+u_nh4) * cnb * fec_bac;
    end
%
% Uptake by small phytoplankton
%
    fl18_6 = (fl1_6+fl2_6) * cnp * fec_phy;
%
% Uptake by diatoms
%
    fl18_13 = (fl1_13+fl2_13) * cndia * fec_dia;
    
    if COCCOLITHOPHORES == 1
        % uptake by coccolithophores
        fl18_coc = (fl1_13 + fl2_13) *cncoc * fec_coc;
    end
    
%
%-----------------------------------------------------------------------
% Calculate the slow-sinking (small) detrital iron terms.
%-----------------------------------------------------------------------
%
% Ingestion by micro/nanozooplankton
%
        fl19_7 = fl9_7 * detfe2detc;
%
% Dissolution and bacterial remineralization
%
% Variable (with bacteria concentration) remineralization rate.
%
        fl19_18 = (mud_dom+mub_dom*bact) * detfe;
%
%-----------------------------------------------------------------------
% Calculate the fast-sinking (large) detrital iron terms.
%-----------------------------------------------------------------------
%
        fl20_18 = ldetfedis * ldetfe;

	if IRON_SCAVENGING == 1
        %
        %-----------------------------------------------------------------------
        % Calculate the scavenging terms.  Scavenging is assumed only to occur 
        % on Fe'.  Fe' (called Fe_free, here) is solved (as in PISCES-v2) as a 
        % second-order polynomial from the following equations:
        %
        %        LTot = FeL + L'
        %        FeD = FeL + Fe'
        %        KeqFeL = FeL/(Fe'L')
        %-----------------------------------------------------------------------
        %
        if SCAV_DOC == 1 
            %
            % Compute LTot as function of total DOC.  DOC is equal to the sum of the
            % prognostic labile/semi-labile and refractory pools.  A constant value
            % of 40 uM is assumed for the refractory pool.
            %
            LTot = max(0.09*(doc+40)-3,LMin);
        end
        cff = KeqFeL*LTot - KeqFeL*fed + 1;
        Fe_free = (-cff + sqrt(cff*cff + 4*KeqFeL*fed)) / (2*KeqFeL);

        if SCAV_DOC == 1
            % Compute scavenging rate as a function of particulate carbon
            rscav = scav_min + Lam_Fe*(detc+ldetc);
        end
        scav = Fe_free * rscav;
    else
        scav = 0;
    end
end  % iron

if COCCOLITHOPHORES ==1 
    %
    %-----------------------------------------------------------------------
    % Calculate the coccolithophore biomass terms.
    %-----------------------------------------------------------------------
    %
    if EXUD_MINVAL == 1
        % Turn off exudation to DON below biomass threshold EXUD_MIN
        coc_exud1 = cocco - exud_min;
        coc_exud2 = max(coc_exud1,0.d0)*gamma1c/coc_exud1;
        coc_exud3 = max(coc_exud1,0.d0)*gamma2c/coc_exud1;
        flcoc_3 = coc_exud2*cocco + coc_exud3*sigmac*cocco;
    else
        % Constant exudation to DON
        flcoc_3 = gamma1c*cocco + gamma2c*sigmac;
    end

    if QUADRATIC_MORT == 1
        % Quadratic mortality to DetN
        flcoc_4 = min(coc_mort_max, coc_mort_max * ...
                                (cocco*cocco/(Km_coc+cocco*cocco))) * ...
                                cocco;
    else
        % Constant mortality to DetN
        flcoc_4 = mup_coc * cocco;
    end

    % Grazing by zoo
    flcoc_7 = GS_coc;  % Micro/nano
    flcoc_11 = GL_coc; % Meso

    % Mortality to LDetC
    flcoc_9 = flcoc_4 * cncoc;

    % exudation to DOC
    flcoc_10 = flcoc_3 * cncoc + ...
                      gamma3c * max(ajc/Vptc-qcoc,0.d0) * ...
                      ajc * cocco * cncoc;

    if IRON == 1
        % exudation to FeD
        flcoc_18 = flcoc_3 * cncoc * fec_coc;
        % mortality to DetFe
        flcoc_19 = flcoc_9 * fec_coc;
    end

    % mortality to PIC
    flcoc_pic = flcoc_9 * PIC2C;
    
    %
    %-----------------------------------------------------------------------
    % Calculate the chlorophyll terms (coccolithophores).
    %-----------------------------------------------------------------------
    %
    flchlc_4  = flcoc_4.*convcc_thetac;    % chlc mortality to DETN
    flchlc_7  = flcoc_7.*convcc_thetac;    % chlc ingestion by microzoo
    flchlc_11 = flcoc_11.*convcc_thetac;   % chlc ingestion by mesozoo
    if K_LOSS == 1
        if K0 == 1
            flchlc_dmg = chlc*chl_dmg;  % constant chl loss [d-1]
        elseif K1 == 1
            flchlc_dmg = chlc*(ajc/Vptc)*chl_dmg;  % variable loss with production
        elseif K2 == 1
            flchlc_dmg = chlc*chl_dmg*thetac*parz;  % variable loss with light and theta
        elseif K3 == 1
            flchlc_dmg = chlor*chl_dmg*parz;  % variable loss with light intensity
        end
        %fl8_dmg_diag(t,k) = fl8_dmg;
        flchlc_4 = flchlc_4 + flchlc_dmg;  % stick on mort. so don't have to include extra term in integration.
    end    
    
    %
    %-----------------------------------------------------------------------
    % Calculate the particulate inorganic carbon terms
    %-----------------------------------------------------------------------
    %    
    flpic_dic = PIC * dissPIC;  % solubilization to DIC
    flpic_alk = 2 * flpic_dic; % solubilization to ALK
    %
    %-----------------------------------------------------------------------
    % Calculate the dissolved inorganic carbon terms
    %-----------------------------------------------------------------------
    %        
    fldic_5 = fl2_5 * cnb;  % uptake by bacteria
    fldic_6 = (fl1_6+fl2_6) * cnp;  % uptake by phyto
    fldic_13 = (fl1_13+fl2_13) * cndia; % uptake by diatoms
    fldic_coc = (fl1_coc+fl2_coc) * cncoc; % uptake by coccos
    fldic_pic = fldic_coc * PIC2C;  % increase due to calcification
    %
    %-----------------------------------------------------------------------
    % Calculate the carbonate alkalinity terms
    %-----------------------------------------------------------------------
    %        
    flalk_5 = fldic_5;  % Loss/gain via bacteria net NH4 uptake
    flalk_6 = (q1-q2)*sigma*cnp;
    flalk_13 = (qd1-qd2)*sigmad*cndia;
    flalk_coc = (q1c-q2c)*sigmac*cncoc;
    flalk_pic = 2 * fldic_pic;  % Loss to calcification
    flalk_nitr = fl2_1 * 6.625;  % Loss to nitrification

end  %coccolithophores
%-----------------------------------------------------------------------
% Sum up pool fluxes
%-----------------------------------------------------------------------
%
%     [NO3] nitrate
bio(t+1,k,ino3) =(-fl1_6 - fl1_13 +fl2_1) .*dtdays + no3;
if COCCOLITHOPHORES == 1
    bio(t+1,k,ino3) =(-fl1_6 - fl1_13 + - fl1_coc + fl2_1) .*dtdays + no3;
    % [NH4] ammonium
    bio(t+1,k,inh4) =(fl5_2 + fl7_2 + fl11_2 -fl2_5 - fl2_6 - fl2_1 - fl2_13 - fl2_coc) .*dtdays + nh4;
    % [DON] Dissolved Organic Nitrogen
    bio(t+1,k,idon) =(fl4_3 + fl6_3 + fl7_3 + fl11_3 + fl13_3 +fl14_3 + flcoc_3 - fl3_5) .*dtdays + don;
    % [DETN] Detritus Nitrogen
    bio(t+1,k,idetn) =(fl5_4 + fl6_4 + fl7_4 + flcoc_4 -fl4_3 - fl4_7) .*dtdays + detn;
    % [zoo] zooplankton
    bio(t+1,k,izoos) =(fl4_7 + fl5_7 + fl6_7 + flcoc_7 - fl7_2 - fl7_3 - fl7_4 - fl7_11) .*dtdays + zoo;
    % [DETC] Detritus Carbon
    bio(t+1,k,idetc) =(fl5_9 + fl6_9 + fl7_9 +flcoc_9 - fl9_7 - fl9_10) .*dtdays + detc;
    % [DOC] Dissolved Organic Carbon
    bio(t+1,k,idoc) =(fl6_10 + fl7_10 + fl9_10 + fl11_10 + fl13_10 + flcoc_10 - fl10_5) .*dtdays + doc;
    % [ZOOL] Mesozooplankton --> added grazing on phy 11.16.18
    bio(t+1,k,izool) =(fl7_11 + fl13_11 + fl6_11 + flcoc_11 -fl11_2 - fl11_3 - fl11_14) .*dtdays + zool;       
else
    % [NO3] nitrate
    bio(t+1,k,ino3) =(-fl1_6 - fl1_13 +fl2_1) .*dtdays + no3;
    % [NH4] ammonium
    bio(t+1,k,inh4) =(fl5_2 + fl7_2 + fl11_2 -fl2_5 - fl2_6 - fl2_1 - fl2_13) .*dtdays + nh4;
    % [DON] Dissolved Organic Nitrogen
    bio(t+1,k,idon) =(fl4_3 + fl6_3 + fl7_3 + fl11_3 + fl13_3 +fl14_3 - fl3_5) .*dtdays + don;
    % [DETN] Detritus Nitrogen
    bio(t+1,k,idetn) =(fl5_4 + fl6_4 + fl7_4 -fl4_3 - fl4_7) .*dtdays + detn;
    % [zoo] zooplankton
    bio(t+1,k,izoos) =(fl4_7 + fl5_7 + fl6_7 -fl7_2 - fl7_3 - fl7_4 - fl7_11) .*dtdays + zoo;
    % [DETC] Detritus Carbon
    bio(t+1,k,idetc) =(fl5_9 + fl6_9 + fl7_9-fl9_7 - fl9_10) .*dtdays + detc;
    % [DOC] Dissolved Organic Carbon
    bio(t+1,k,idoc) =(fl6_10 + fl7_10 + fl9_10 + fl11_10 + fl13_10 -fl10_5) .*dtdays + doc;
    % [ZOOL] Mesozooplankton --> added grazing on phy 11.16.18
    bio(t+1,k,izool) =(fl7_11 + fl13_11 + fl6_11 -fl11_2 - fl11_3 - fl11_14) .*dtdays + zool;    
end
% [bact] bacteria
bio(t+1,k,ibac) =(fl2_5 + fl3_5 -fl5_2 - fl5_4 - fl5_7) .*dtdays + bact;
% [phyto] phytoplankton --> added mesozoo grazing 11.16.18
bio(t+1,k,iphy) =(fl1_6 + fl2_6 -fl6_3 - fl6_4 - fl6_7 - fl6_11) .*dtdays + phyto;
% [chlor] Chlorophyll a --> added mesozoo grazing 11.16.18
bio(t+1,k,ichl) =(fl1_8 + fl2_8 -fl8_4 - fl8_7 - fl8_11) .*dtdays + chlor;
% [diat] diatoms
bio(t+1,k,idia) =(fl1_13 + fl2_13 -fl13_3 - fl13_11 - fl13_14) .*dtdays + dia;
% [sil] silica
bio(t+1,k,isil) =(fl16_12 -fl12_13) .*dtdays + sil;
% [LDETN] Large detritus Nitrogen
bio(t+1,k,ildetn) =(fl11_14 + fl13_14 -fl14_3) .*dtdays + ldetn;
% [LDETC] Large detritus Carbon
bio(t+1,k,ildetc) =(fl11_15 + fl13_15 -fl15_10) .*dtdays + ldetc;
% [LDETS] Large detritus Silica
bio(t+1,k,ildets) =(fl11_16 + fl13_16 - fl16_12) .*dtdays + ldets;
% [chld] Diatom Chlorophyll a
bio(t+1,k,ichld) =(fl1_17 + fl2_17 -fl17_14 - fl17_11) .*dtdays + chld;

if IRON == 1
%
%-----------------------------------------------------------------------
% Iron model with constant scavenging of Fe'.
% Fe' is computed from the equilibrium equation Fe' + L <--> FeL (Liu
% and Millero, 2002) with a Keq in the tunable range 10e12 - 10e14 and a
% constant total ligand concentration (LTot).  
% Fscav is the fraction of the scavenged iron placed into the fast-sinking 
% detrital pool (LDetFe) to be remineralized at depth, with the remaining 
% 1-Fscav lost from the system (Moore and Braucher 2008). 
%
% Z.W. 11.28.18
%-----------------------------------------------------------------------
%        
    if IRON_BACTERIA == 1
        if COCCOLITHOPHORES == 1
            % [FeD] Dissolved inorganic iron 
            bio(t+1,k,ifed) = (flcoc_18 + fl5_18 + fl6_18 + fl7_18 + fl11_18 + fl13_18 + ...
                                         fl19_18 + fl20_18 - fl18_5 - fl18_6 - fl18_13 - scav) * ...
                                         dtdays + fed;
            % [DetFe] Slow-sinking (small) detrital iron
            bio(t+1,k,idetfe) = (flcoc_19 + fl5_19 + fl6_19 + fl7_19 - fl19_18 - fl19_7) * ...
                                           dtdays + detfe;            
        else
            % [FeD] Dissolved inorganic iron
            bio(t+1,k,ifed) = (fl5_18 + fl6_18 + fl7_18 + fl11_18 + fl13_18 + ...
                                        fl19_18 + fl20_18 - fl18_5 - fl18_6 - fl18_13 - scav) * ...
                                        dtdays + fed;
            % [DetFe] Slow-sinking (small) detrital iron
            bio(t+1,k,idetfe) = (fl5_19 + fl6_19 + fl7_19 - fl19_18 - fl19_7) * ...
                                           dtdays + detfe;
        end
    else
        if COCCOLITHOPHORES == 1
            % [FeD] Dissolved inorganic iron 
            bio(t+1,k,ifed) = (flcoc_18 + fl6_18 + fl7_18 + fl11_18 + fl13_18 + ...
                                         fl19_18 + fl20_18 - fl18_6 - fl18_13 - fl18_coc - scav) * ...
                                         dtdays + fed;
            % [DetFe] Slow-sinking (small) detrital iron
            bio(t+1,k,idetfe) = (flcoc_19 + fl6_19 + fl7_19 - fl19_18 - fl19_7) * ...
                                           dtdays + detfe;            
        else
            % [FeD] Dissolved inorganic iron
            bio(t+1,k,ifed) = (fl6_18 + fl7_18 + fl11_18 + fl13_18 + ...
                                        fl19_18 + fl20_18 - fl18_6 - fl18_13 - scav) * ...
                                        dtdays + fed;
            % [DetFe] Slow-sinking (small) detrital iron
            bio(t+1,k,idetfe) = (fl6_19 + fl7_19 - fl19_18 - fl19_7) * ...
                                           dtdays + detfe;
        end        
    end  % IRON_BACTERIA
    
%     [LDetFe] Fast-sinking (large) detrital iron
    bio(t+1,k,ildetfe) = (fl11_20 + fl13_20 + fscav*scav - fl20_18) * dtdays + ldetfe;
end  % IRON

if COCCOLITHOPHORES == 1
    bio(t+1,k,iCocc) = (fl1_coc + fl2_coc - flcoc_3 - flcoc_4 - flcoc_7 - flcoc_11 - flcoc_pic)*dtdays+cocco;
    bio(t+1,k,iChlC) = (fl1_chlc + fl2_chlc - flchlc_4 - flchlc_7 - flchlc_11) .*dtdays + chlc;
    bio(t+1,k,iPIC) = (flcoc_pic + fldic_pic + fl7_pic + fl11_pic - flpic_dic) * dtdays + PIC;
    bio(t+1,k,iDIC) = (fl5_dic + flpic_dic + fl7_dic + fl11_dic - ...
                                 fldic_5 - fldic_6 - fldic_13 - fldic_coc - fldic_pic) * dtdays + DIC;
    bio(t+1,k,iALK) = (fl5_alk + flpic_alk + fl7_alk + fl11_alk - ...
                                 flalk_5 - flalk_6 - flalk_13 - flalk_coc - flalk_pic - flalk_nitr) * dtdays + ALK;
end
end  % K Loop
%
%-----------------------------------------------------------------------
%  Vertical sinking terms.
%-----------------------------------------------------------------------
%
%  Reconstruct vertical profile of selected biological constituents
%  'Bio(:,:,isink)' in terms of a set of parabolic segments within each
%  grid box. Then, compute semi-Lagrangian flux due to sinking.
%

for isink=1:nsink % SINK_LOOP
    indx=idsink(isink);
%
%  Copy concentration of biological particulates into scratch array
%  'qc' (q-central, restrict it to be positive) which is hereafter
%  interpreted as a set of grid-box averaged values for biogeochemical
%  constituent concentration.
%
% for now let's try vertically advecting the already biologically time-stepped fields.
qc = zeros(1,N);
wr = zeros(1,N);
wl = zeros(1,N);
br = zeros(1,N);
bl = zeros(1,N);
ksource = zeros(1,N);
%
for k=1:N
        qc(k)=max(bio(t+1,k,indx),0.0);
end
%
for k=N-1:-1:1
        fc(k+1)=(qc(k+1)-qc(k)).*hz_inv2(k);
end
for k=2:N-1
        dltr=hz(k).*fc(k+1);
        dltl=hz(k).*fc(k-1+1);
        cff=hz(k-1)+2.0.*hz(k)+hz(k+1);
        cffr=cff.*fc(k+1);
        cffl=cff.*fc(k-1+1);
%
%  Apply PPM monotonicity constraint to prevent oscillations within the grid box.
%
        if((dltr*dltl)<=0.0)
            dltr=0.0;
            dltl=0.0;
        elseif(abs(dltr)>abs(cffl)) 
            dltr=cffl;
        elseif(abs(dltl)>abs(cffr)) 
            dltl=cffr;
        end
%
%  Compute right and left side values (bR,bL) of parabolic segments
%  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
%
%  NOTE: Although each parabolic segment is monotonic within its grid
%        box, monotonicity of the whole profile is not guaranteed,
%        because bL(k+1)-bR(k) may still have different sign than
%        qc(i,k+1)-qc(i,k).  This possibility is excluded,
%        after bL and bR are reconciled using WENO procedure.
%
        cff=(dltr-dltl).*hz_inv3(k);
        dltr=dltr-cff.*hz(k+1);
        dltl=dltl+cff.*hz(k-1);
        br(k)=qc(k)+dltr;
        bl(k)=qc(k)-dltl;
        wr(k)=(2.0.*dltr-dltl).^2;
        wl(k)=(dltr-2.0.*dltl).^2;
end
    cff=1.0e-14;
    for k=2:N-2
            dltl=max(cff,wl(k));
            dltr=max(cff,wr(k+1));
            br(k)=(dltr.*br(k)+dltl.*bl(k+1))./(dltr+dltl);
            bl(k+1)=br(k);
    end
    fc(N+1)=0.0;  % NO-flux boundary condition
    br(N)=qc(N);    % default strictly monotonic
    bl(N)=qc(N);    % conditions
    br(N-1)=qc(N);
    bl(2)=qc(1);    % bottom grid boxes are
    br(1)=qc(1);    % re-assumed to be
    bl(1)=qc(1);    % piecewise constant.
%
%  Apply monotonicity constraint again, since the reconciled interfacial
%  values may cause a non-monotonic behavior of the parabolic segments
%  inside the grid box.
%
for k=1:N
        dltr=br(k)-qc(k);
        dltl=qc(k)-bl(k);
        cffr=2.0.*dltr;
        cffl=2.0.*dltl;
        if((dltr*dltl)<0.0)
            dltr=0.0;
            dltl=0.0;
        elseif(abs(dltr)>abs(cffl)) 
            dltr=cffl;
        elseif(abs(dltl)>abs(cffr)) 
            dltl=cffr;
        end
        br(k)=qc(k)+dltr;
        bl(k)=qc(k)-dltl;
end 
%
%  After this moment reconstruction is considered complete. The next
%  stage is to compute vertical advective fluxes, FC. It is expected
%  that sinking may occurs relatively fast, the algorithm is designed
%  to be free of CFL criterion, which is achieved by allowing
%  integration bounds for semi-Lagrangian advective flux to use as
%  many grid boxes in upstream direction as necessary.
%
%  In the two code segments below, WL is the z-coordinate of the
%  departure point for grid box interface z_w with the same indices;
%  FC is the finite volume flux; ksource(:,k) is index of vertical
%  grid box which contains the departure point (restricted by N(ng)).
%  During the search: also add in content of whole grid boxes
%  participating in FC.
%
cff=dtdays.*abs(wbio(isink));
for k=1:N
        fc(k-1+1)=0.0;
        wl(k)=z(k-1+1)+cff;
        wr(k)=hz(k).*qc(k);
        ksource(k)=k;
end
for k=1:N
    for ks=k:N-1
        if(wl(k)>z(ks+1))
            ksource(k)=ks+1;
            fc(k-1+1)=fc(k-1+1)+wr(ks);
        end
    end 
end 
%
%  Finalize computation of flux: add fractional part.
%
for k=1:N
        ks=ksource(k);
        cu=min(1.0,(wl(k)-z(ks-1+1)).*hz_inv(ks));
        fc(k-1+1)=fc(k-1+1)+hz(ks).*cu.*(bl(ks)+cu.*(0.5.*(br(ks)-bl(ks))-(1.5-cu).*(br(ks)+bl(ks)-2.0.*qc(ks))));
end
for k=1:N
        bio(t+1,k,indx)=qc(k)+(fc(k+1)-fc(k-1+1)).*hz_inv(k);
end 
if FAST_REMIN == 1
    if indx == 16  % LDetS
        cff=fc(1)*hz_inv(1);
        bio(t+1,1,isil) = bio(t+1,1,ildets)+cff;
    end
end
end %  SINK_LOOP 


%for indx=1:nbt
%    for k=1:n
%            bio(t+1,k,indx)=max(bio(t+1,k,indx),0.0);
%    end 
%end 
end %  time loop

if PLOTS == 1
figure;colormap('jet');

dayPerDT = dt/86400;
time=(1/dt:dayPerDT:num_days);

namelist={'NO_3','NH_4','DON','Det_N','Bac','Phyto','Zoo_S',...
                  'ChlS','DetC','DOC','Zoo_L','SiOH_4','Dia','LDet_N',...
                  'LDet_C','LDet_S','ChlD','FeD','DetFe','LDetFe', ...
                  'Cocco','ChlC','PIC','DIC','ALK'};
if COCCOLITHOPHORES == 1
    %vars2plot = [iCocc,iPIC,iDIC,iALK];  % indices of variables to plot
    vars2plot = [iphy, idia, izoos, izool, iCocc, ifed];  % indices of variables to plot
    %vars2plot = [ino3,izool,isil,idia,ildets];  % indices of variables to plot
else
    vars2plot = [ibac,iphy,izoos,izool,idia,ifed];  % indices of variables to plot              
    %vars2plot = [ino3,izool,isil,idia,ildets];  % indices of variables to plot
end

if APRIL_START == 1
    month = {'Apr','May','Jun','Jul','Aug','Sept','Oct','Nov','Dec',...
                     'Jan','Feb','Mar'};
else
    month = {'Jun','Jul','Aug','Sept','Oct','Nov','Dec',...
                     'Jan','Feb','Mar','Apr','May','Jun'};
end
ax_ticks = (0:30:num_days);  % Monthly ticks
ax_ticks_labels = month;  % Monthly ticks
% Plot tracers in time
ax1 = subplot(211);
plot(ax1,time,squeeze(bio(:,end,vars2plot)),'linewidth',2);
ax1.XLim = [0 num_days];
ax1.XTick = ax_ticks;
ax1.XTickLabel = ax_ticks_labels;
xlabel('Month (end of)','fontsize',12);
ax1.Title.String ='Tracer Concentrations';
legend(namelist(vars2plot),'Location','Best');

% Plot Chl:N ratio
ax2 = subplot(212);
plot(ax2,time,squeeze(bio(:,end,ichl)./bio(:,end,iphy)),'linewidth',2);hold on
plot(ax2,time,squeeze(bio(:,end,ichld)./bio(:,end,idia)),'linewidth',2);
if COCCOLITHOPHORES == 1
    plot(ax2,time,squeeze(bio(:,end,iChlC)./bio(:,end,iCocc)),'linewidth',2);
end
ax2.XLim = [0 num_days];
ax2.XTick = ax_ticks;
ax2.XTickLabel = ax_ticks;
if COCCOLITHOPHORES == 1
    leg_chl2n = {'Small Phyto','Diatoms','Coccolithophores'};
else
    leg_chl2n = {'Small Phyto','Diatoms'};
end
legend(leg_chl2n,'Location','Best');
ax2.Title.String = 'Chl:N Ratios';
xlabel('Days','fontsize',12);
ylabel('Chl:N ratio [mg Chl mmol N^{-1}]')

end

