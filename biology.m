
clear all

%-----------------------------------------------------------------------
% Local Parameters
%-----------------------------------------------------------------------
%
plots = 1;  % Flag to turn on plotting
N = 20;     % Number of vertical levels
h = 10; % meters
nbt = 17;   % Number of biological variables
nsink = 5;  % Number of biological variables that sink
dayspersecond = 1./86400; 
bioparfrac = 0.43; %  Fraction of visible light (par0) active in photosynthesis (400nm-700nm wave lengths)
acof = 0.851; % Eppley 1972 % acof = 0.59   ! Norberg 2004
bcof = 0.0633;      % Norberg 2004
topt_phy = 12;  %      Topt_phy = 19 ! Boyd et al. 2013
topt_dia = 8;   %  Topt_dia = 12
w_phy = 13;     %  w_phy = 32    ! Boyd et al. 2013
w_dia = 13;
rho0 = 1025.0d0;
cp = 3985.0d0;
cnp = 6.0d0;        % Carbon:Nitrogen ratio for phytoplankton
cnd = 4.6d0;        % Carbon:Nitrogen ratio for DOM
cnb = 5.38d0;       % Carbon:Nitrogen ratio for bacteria
cndia = 6.625d0;    % Carbon:Nitrogen ratio for diatoms
cnz = 5.02d0;       % Carbon:Nitrogen ratio for zooplankton
vp = 1.2d0;         % Phytoplankton maximum growth rate
vpdia = 1.5d0;      % Diatom maximum growth rate
kdia_no3 = 1.0d0;   % Kdia_NO3 = 3.0d0    ! Jiang
kp_no3 = 0.59d0;    % Half-saturation for phytoplankton NO3 uptake
kp_nh4 = 3.39d-2;
kdia_nh4 = kp_nh4;  % Kdia_NH4 = 0.5d0    ! Jiang
kdia_si  = 2.0d0;   % Kdia_Si   = Kp_NO3
kb_c = 0.51d0;      % Half saturation constant for DOCC
kb_don = 1.34d0;    % Bacteria half saturation constant for DON uptake
kb_a = 1.0d0;       % Bacteria half saturation constant for NH4 uptake
kz_i = 1.2d0;       % Half saturation constant for ingestion
kl_i = 0.69d0;      % Large zooplankton half saturation constant for ingestion
anitrof = 0.05;
snpd = 2.0;         % Silica:Nitrogen ratio
ldetndis = 0.1;
ldetcdis = 0.1;
ldetsdis = 0.1;
psi = 1.48d0;       % Phytoplankton ammonium inhibition parameter
attsw = 0.04d0;     % Light attenuation due to seawater (1/m)
attphyt = 0.02d0;   % Light attenuation by phytoplankton
attphytdia = 0.02d0;    % Light attenuation by diatoms
alpha = 1.227d0;    % Initial slope of the P-I curve
alphadia = 0.75d0;   % Initial slope of the P-I curve for diatoms
mx_chlcr = 4.15d-2; % Maximum chlorophyll-a to Carbon ratio
mx_chlcrd = 5.5d-2; % Maximum chlorophyll-a to Carbon ratio (diatioms)
prefz_b = 1.14d0;   % Zooplankton preference for bacteria
prefz_d = 0.90d0;   % Zooplankton preference for detritus
prefl_z = 0.94d0;   % Large zooplankton preference for small zooplankton
prefl_p = 0.9d0;    % Large zooplankton preference for small phytoplankton
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
muz_a = 0.12d0;     % Zooplankton specific excretion rate
muzg_a = 0.106d0;   % Growth specific excretion rate
mul = 0.1d0;        % Large zooplankton mortality
mulg_a = 0.1d0;     % Large zooplankton growth specific excretion rate
omegal = 0.15d0;    % Large zooplankton fraction to DON
omegaz_d = 0.18d0;  % Detrital fraction of zooplankton mortality
betaz_b = 0.87d0;   % Zooplankton assimilation efficiency of bacteria
betaz_p = 0.61d0;   % Zooplankton assimilation efficiency of phytoplankton
betaz_d = 0.45d0;   % Zooplankton assimilation efficiency of detritus
betal_dia = 0.72d0; % Large zooplankton assimilation efficiency of diatoms
betal_z = 0.44d0;   % Large zooplankton assimilation efficiency of zooplankton
betal_p = 0.72d0;   % Large zooplankton assimilation efficiency of phytoplankton
%gamma1p = 6.6d-2;   % Fraction of total primary production exuded
gamma1p = 2.4d-2;   % Fraction of total primary production exuded
%gamma2p = 0.115d0;  % Fraction of growth-specific primary production exuded
gamma2p = 8.6d-2;  % Fraction of growth-specific primary production exuded
gamma3p = 0.538d0;  % Fraction of nutrient-limited primary production exuded
%gamma1dia = 6.6d-2; % Diatom fraction of total primary production exuded
gamma1dia = 2.4d-2; % Diatom fraction of total primary production exuded
%gamma2dia = 0.112d0;    % Diatom fraction of growth-specific primary production exuded
gamma2dia = 8.6d-2;    % Diatom fraction of growth-specific primary production exuded
gamma3dia = 6.6d-2; % Diatom fraction of nutrient limited primary production exuded
epsilonz_a = 0.6d0; % Ammonium fraction of zooplankton excretion
epsilonl = 0.51d0;  % Large zooplankton fraction to ammonium
%
%  Set time-stepping according to the number of iterations.
%
dt = 900;    % seconds
num_days = 180;
dtdays=dt*dayspersecond;  % [d ts-1]
time_loops = num_days*86400/dt;

% Specify x,y location for 1-D run
x = 250; y = 400;

% Read in temperature data
f1 = '/home/server/homes/pi/zwallace/';
ft = 'ROMS/zw_tools/temp_spinup_MY25.mat';
tempfile = fullfile(f1,ft);
m = matfile(tempfile);
temp = squeeze(m.temp_spinup_MY25(x,y,:,:));
lt = length(temp(1,:));
% Interpolate to 3-day forcing data
ns3dy = 60*60*24*3;  % number of seconds in a 3-day period
temp_dt(time_loops+1,N) = zeros;
for z=1:N
    temp_dt(:,z) = interp1((0:ns3dy:lt*ns3dy-1),temp(z,:),0:dt:num_days*86400);
end

%Read in shortwave data
f2 = 'ROMS/Project_Patagonia/bulk';
fradsw = 'roms_blk_ERA_icecov_RHO_grid_3day_1900.nc';
swfile = fullfile(f1,f2,fradsw);
% swap swrad so that starting date is in winter (~July)
tmp = ncread(swfile,'radsw');
tmp = squeeze(tmp(x,y,:));  % On shelf
ll = length(tmp);
radsw(ll) = zeros;
radsw(1:ll/2) = tmp(ll/2+1:end);
radsw(ll/2+1:end) = tmp(1:ll/2);
% Interpolate to hourly shortwave radiation
%srflx = interp1((0:3600:ll*3600-1),radsw,0:dt:num_days*86400)
% Interpolate to 3-day shortwave radiation
srflx = interp1((0:ns3dy:ll*ns3dy-1),radsw,0:dt:num_days*86400);

% Depth parameters
z = -h:h/N:-h/N;
hz = ones(1,N)*h/N;

%  Indices of biological variables
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
%
% Initial Conditions
%
bio = zeros(time_loops,N,nbt);
%
bio(1,:,ino3) = 30.0d0;
bio(1,:,inh4) = 0.1d0;
bio(1,:,idon) = 0.15d0;
bio(1,:,idetn) = 1.0d0;
bio(1,:,ibac) = 0.1d0;
bio(1,:,iphy) = 100.d0;
bio(1,:,izoos) = 0.01d0;
bio(1,:,ichl) = 1.59d0;
bio(1,:,idetc) = 6.625d0;
bio(1,:,idoc) = 0.9d0;
bio(1,:,izool) = 0.01d0;
bio(1,:,isil) = 30.0d0;
bio(1,:,idia) = 100.d0;
bio(1,:,ildetn) = 0.02d0;
bio(1,:,ildetc) = 0.002d0;
bio(1,:,ildets) = 0.02d0;
bio(1,:,ichld) = 1.59d0;
%
%  Set vertical sinking velocity vector in the same order as the
%  identification vector, idsink.
%
vsink = 1.35d0;
vsinkld = 50.0d0;
wbio(1)=vsink;
wbio(2)=vsink;
wbio(3)=vsinkld;
wbio(4)=vsinkld;
wbio(5)=vsinkld;
%
%  Set vertical sinking indentification vector.
%
idsink(1)=idetn;
idsink(2)=idetc;
idsink(3)=ildetn;
idsink(4)=ildetc;
idsink(5)=ildets;
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
%
for t = 1:time_loops-1
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
%
% Diatoms
%
denomd = dia .* cndia .* 12.0d0;
thetad = max(chld./denomd, 0.0d0);
%  Eppley 1974 Temperature modification of phytoplankton growth
%  which gives Vpt =3.26 at 21C and 5.09 at 28C
%         Vpt = 0.851*(1.066**temp)
%         Vptdia = 0.851*(1.066**temp)
%
% Formulation for growth rate taken from Boyd et al. 2013.
%
         vpt = acof*exp(bcof*temp_dt(t,k))*(1-((temp_dt(t,k)-topt_phy)/(w_phy/2))^2);
         vptdia = acof*exp(bcof*temp_dt(t,k))*(1-((temp_dt(t,k)-topt_dia)/(w_dia/2))^2);
%
% Formulation for growth rate from Yvette 1/11/2019
%       vpt = vp.*exp(0.0633.*temp)./(1.0+exp(temp.*2.0-9.0));
%       vptdia = vpdia.*exp(0.0633.*temp)./(1.0+exp(temp.*2.0-9.0));
%
% Maintain background population by setting growth rate to 0.01 [d-1]
% where it would otherwise be zero.
%
vpt = max(1.0d-2,vpt);
vptdia = max(1.0d-2,vptdia);
%
% Prevent growth rates from exceeding prescribed maxima
%
vpt = min(vpt,vp);
vptdia = min(vptdia,vpdia);
%
%-----------------------------------------------------------------------
% Calculate the nutrient limitation terms on small phytoplankton
% and diatoms.
%-----------------------------------------------------------------------
%
% Small Phytoplankton
%
ak=kp_no3/kp_nh4;
deno=kp_no3+no3+nh4.*ak;
q1=no3.*exp(-psi.*nh4) ./ deno; % nondimensional NO3 limitation factors
q2 =(nh4.*ak) ./ deno;  % nondimensional NH4 limitation factors
q = q1+q2;
if( q == 0.0 )
    fvalue = 0.0;
else
    fvalue=q1./(q1+q2);
end
%
% Diatoms
%
akd=kdia_no3./kdia_nh4;
denod=kdia_no3+no3+nh4.*akd;
qd1=no3.*exp(-psi.*nh4) ./ denod;   % nondimensional NO3 limitation factors
qd2 =(nh4.*akd) ./ denod;   % nondimensional NH4 limitation factors
qd3 = sil ./(sil+kdia_si);  % nondimensional SIL limitation factors
qd = min(qd1+qd2, qd3);
if( qd == 0.0 )
    fvalued = 0.0;
else
    fvalued=qd1./(qd1+qd2);
end
%
%  Calculate Surface Photosynthetically Available Radiation (PAR).  The
%  net shortwave radiation is scaled back to Watts/m2 and multiplied by
%  the fraction that is photosynthetically available, PARfrac.
%
    parsur=bioparfrac*srflx(t);
%
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
else
    kdpar = attsw;
end
%
% Compute biological light attenuation as a function of chlorophyll.
%
parz = parsur*exp(z(k).*(kdpar +(attphyt.*chlor+attphytdia.*chld)));
if (parz < 1.0d-6)
    parz = 0.0;
end
aj = vpt.*theta.*alpha.*parz./sqrt((vpt.*vpt) +(theta.*theta.*alpha.*alpha.*parz.*parz));
ajd = vptdia.*thetad.*alphadia.*parz./sqrt((vptdia.*vptdia)+(thetad.*thetad.*alphadia.*alphadia.*parz.*parz));
%           Check for low light.  If light is low, turn on nitrification.
if(parz >= 0.1 * parsur)
    anitrf = 0.;
else
    anitrf = anitrof;
end
% PAR0 is zero, no light due to earth rotation or earth seaonal variation.
else
    parz  = 0.0;
    aj    = 0.0;
    ajd   = 0.0;
    anitrf = anitrof;   %no light, turn on nitrification
end
sigma=aj.*q;
sigmad=ajd.*qd;
if(parz<=0. || theta<=0 || thetad<=0)
    rho_c =0.;
    rhod_c =0.;
else
    rho_c = mx_chlcr.*sigma./(alpha.*theta.*parz);
    rhod_c = mx_chlcrd.*sigmad./(alphadia.*thetad.*parz);
end
%
%-----------------------------------------------------------------------
% Calculate small and large zooplankton grazing terms.
%-----------------------------------------------------------------------
%
%   coeff. and formulation from BATS paper.
value=((kz_i.*(phyto+(prefz_b.*bact)+(prefz_d .*detn)))+((phyto.*phyto)+(prefz_b .*bact.*bact)+(prefz_d .*detn.*detn)));
if(value > 0.0d0)
    g1=(gmax.*zoo.*phyto.*phyto)./value;
    g2=(gmax.*zoo.*prefz_b .*bact.*bact)./value;
    g3=(gmax.*zoo.*prefz_d .*detn.*detn)./value;
else
    g1=0.0d0;
    g2=0.0d0;
    g3=0.0d0;
end
%  Include mesozooplankton grazing on phy, dia, and microzoo.
%  Z.W. 11.16.18
%
valuel=kl_i.*(dia+prefl_z.*zoo+prefl_p.*phyto) +dia.*dia +prefl_z.*zoo.*zoo +prefl_p.*phyto.*phyto;
if(valuel~=0.)
    g4=grazl .* zool .* dia.*dia ./ valuel;
    g5=grazl .* zool .* prefl_z .* zoo.*zoo ./ valuel;
    g6=grazl .* zool .* prefl_p .* phyto.*phyto ./ valuel;
else
    g4 = 0.0d0;
    g5 = 0.0d0;
    g6 = 0.0d0;
end
%
%-----------------------------------------------------------------------
% Calculate the bacteria terms.
%-----------------------------------------------------------------------
%
docn = don.*cnd;    %DOC_nitro
docc = max(1.0d-12, doc-docn);  %DOC_carbo
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
fl5_4 = mub_d * bact;      % bact mortality to DETN
fl5_7 = g2 * betaz_b;      % bact ingestion by zoo
fl5_9 = fl5_4 * cnb;       % bact mortality to DETC
%
%-----------------------------------------------------------------------
% Calculate the Phytoplankton terms.
%    3.12.18  --> Quadratic mortality
%    5.25.18  --> Minimum respiration threshold
%-----------------------------------------------------------------------
%
% Constant exudation to DON
%
fl6_3 = gamma1p * phyto + gamma2p * sigma * phyto;
%
% Linear mortality to DetN
%
fl6_4 = mup_d * phyto;
%
% Ingestion by zooplankton
%
fl6_7 = g1;     % Micro/nano
fl6_11 = g6;    % Meso
%
% Mortality to DetC
%
fl6_9 = fl6_4 * cnp;
%
% Exudation to DOC
%
fl6_10 = fl6_3 * cnp + gamma3p * max(aj./vpt-q,0.0d0) * aj * phyto * cnp;
%
%-----------------------------------------------------------------------
% Calculate the microzooplankton terms.
%-----------------------------------------------------------------------
%
detc2detn = max(detc./detn,0.0d0);
%
% Excretion to NH4
%
fl7_2  = epsilonz_a * muz_a * zoo +(betaz_p * g1 +betaz_b * g2 +betaz_d * g3) * muzg_a;
%
% Mortality to DON
%
fl7_3  = omegaz_d * muz_a * zoo;
%
% Growth-specific respiration and mortality to detrital pools
%
fl7_4  =(1.0d0-betaz_p) * g1 +(1.0d0-betaz_b) * g2 +(1.0d0-betaz_d) * g3;
fl7_4  = fl7_4 +(1.0d0 - epsilonz_a - omegaz_d) * muz_a * zoo;
fl7_9  =(1.0d0-betaz_p) * g1 * cnp  +(1.0d0-betaz_b) * g2 * cnb  +(1.0d0-betaz_d) * g3 * detc2detn;
fl7_9  = fl7_9 +(1.0d0 - epsilonz_a - omegaz_d) * muz_a * zoo * cnz;
%
% Mortality to DOC
%
fl7_10 = fl7_3 .* cnz;
%
% Grazing by mesozooplankton
%
fl7_11 = g5;
%
%-----------------------------------------------------------------------
% Calculate the chlorophyll terms (small phytoplankton).
%-----------------------------------------------------------------------
%
fl8_4  = fl6_4.*convc_theta;    % chl mortality to DETN
fl8_7  = fl6_7.*convc_theta;    % chl ingestion by microzoo
fl8_11 = fl6_11.*convc_theta;   % chl ingestion by mesozoo
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
mort_zool = mul * zool;
%
% Mortality and excretion to NH4
%
fl11_2  = epsilonl * mort_zool +(betal_dia * g4 + betal_z * g5 +betal_p * g6) * mulg_a;
%
% Mortality to DON
%
fl11_3  = omegal * mort_zool;
%
% Mortality to DOC
%
fl11_10 = fl11_3 * cnz;
%
% Mortality and excretion to LDetN
%
fl11_14  =(1.0d0 - epsilonl - omegal) * mort_zool + (1-betal_dia) * g4 +(1-betal_z) * g5 +(1-betal_p) * g6;
%
% Mortality and excretion to LDetC
%
fl11_15  =(1.0d0 - epsilonl - omegal) * mort_zool * cnz +(1-betal_dia) * g4 * cndia +(1-betal_z) * g5 * cnz + (1-betal_p) * g6 * cnp;
%
%-----------------------------------------------------------------------
% Calculate the silica terms terms.
% Variable Si:N added -- Z.W. 12.22.17
%-----------------------------------------------------------------------
%
% Fe-independent silica uptake by diatoms
%
fl12_13 = ajd * qd * dia * snpd;
%
%-----------------------------------------------------------------------
% Calculate the Diatom  terms.
%    12.22.17 --> Variable Si:N
%    3.12.18  --> Quadratic mortality
%    5.25.18  --> Minimum respiration threshold
%-----------------------------------------------------------------------
%
% Constant exudation to DON
%
fl13_3 = gamma1dia.*dia+gamma2dia.*sigmad.*dia;
%
% Exudation to DOC
%
fl13_10 = fl13_3.*cndia+gamma3dia * max(ajd./vptdia-qd,0.0d0).*ajd .*dia.*cndia;
%
% Ingestion by mesozooplankton
%
fl13_11 = g4;
%
% Linear mortality to LDetN
%
fl13_14 = amud1 * dia;
%
% Mortality to LDetC
%
fl13_15 = fl13_14 * cndia;
%
% Constant Si:N ratio (no Fe-dependence) for mortality to LDetS.
%
fl13_16 = fl13_14 * snpd;
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
fl16_12 = ldetsdis.*ldets;
%
%-----------------------------------------------------------------------
% Calculate the chlorophyll terms (diatoms).
%-----------------------------------------------------------------------
%
% chl ingestion by mesozoo
fl17_11 = fl13_11.*convcd.*thetad;
% chl mortality to LDETN
fl17_14 = fl13_14.*convcd.*thetad;
%
%-----------------------------------------------------------------------
% Sum up pool fluxes
%-----------------------------------------------------------------------
%
%     [NO3] nitrate
bio(t+1,k,ino3) =(-fl1_6 - fl1_13 +fl2_1) .*dtdays + no3;
%     [NH4] ammonium
bio(t+1,k,inh4) =(fl5_2 + fl7_2 + fl11_2 -fl2_5 - fl2_6 - fl2_1 - fl2_13) .*dtdays + nh4;
%     [DON] Dissolved Organic Nitrogen
bio(t+1,k,idon) =(fl4_3 + fl6_3 + fl7_3 + fl11_3 + fl13_3 +fl14_3 - fl3_5) .*dtdays + don;
%     [DETN] Detritus Nitrogen
bio(t+1,k,idetn) =(fl5_4 + fl6_4 + fl7_4 -fl4_3 - fl4_7) .*dtdays + detn;
%     [bact] bacteria
bio(t+1,k,ibac) =(fl2_5 + fl3_5 -fl5_2 - fl5_4 - fl5_7) .*dtdays + bact;
%     [phyto] phytoplankton --> added mesozoo grazing 11.16.18
bio(t+1,k,iphy) =(fl1_6 + fl2_6 -fl6_3 - fl6_4 - fl6_7 - fl6_11) .*dtdays + phyto;
%     [zoo] zooplankton
bio(t+1,k,izoos) =(fl4_7 + fl5_7 + fl6_7 -fl7_2 - fl7_3 - fl7_4 - fl7_11) .*dtdays + zoo;
%     [chlor] Chlorophyll a --> added mesozoo grazing 11.16.18
bio(t+1,k,ichl) =(fl1_8 + fl2_8 -fl8_4 - fl8_7 - fl8_11) .*dtdays + chlor;
%     [DETC] Detritus Carbon
bio(t+1,k,idetc) =(fl5_9 + fl6_9 + fl7_9 -fl9_7 - fl9_10) .*dtdays + detc;
%     [DOC] Dissolved Organic Carbon
bio(t+1,k,idoc) =(fl6_10 + fl7_10 + fl9_10 + fl11_10 + fl13_10-fl10_5) .*dtdays + doc;
%     [ZOOL] Mesozooplankton --> added grazing on phy 11.16.18
bio(t+1,k,izool) =(fl7_11 + fl13_11 + fl6_11 -fl11_2 - fl11_3 - fl11_14) .*dtdays + zool;
%     [diat] diatoms
bio(t+1,k,idia) =(fl1_13 + fl2_13 -fl13_3 - fl13_11 - fl13_14) .*dtdays + dia;
%     [sil] silica
bio(t+1,k,isil) =(fl16_12 -fl12_13) .*dtdays + sil;
%     [LDETN] Large detritus Nitrogen
bio(t+1,k,ildetn) =(fl11_14 + fl13_14 -fl14_3) .*dtdays + ldetn;
%     [LDETC] Large detritus Carbon
bio(t+1,k,ildetc) =(fl11_15 + fl13_15 -fl15_10) .*dtdays + ldetc;
%     [LDETS] Large detritus Silica
bio(t+1,k,ildets) =(fl13_16 -fl16_12) .*dtdays + ldets;
%     [chld] Diatom Chlorophyll a
bio(t+1,k,ichld) =(fl1_17 + fl2_17 -fl17_14 - fl17_11) .*dtdays + chld;
%     
end 
%======================================================================
%Benthic Sub Model
%======================================================================
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
        qc(k)=max(bio(t,k,indx),0.0);
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
end %  SINK_LOOP 


%for indx=1:nbt
%    for k=1:n
%            bio(t+1,k,indx)=max(bio(t+1,k,indx),0.0);
%    end 
%end 
end %  time loop

if plots == 1
figure
%time=dt:dt:num_days*86400;
dayPerDT = dt/86400;
time=(1/dt:dayPerDT:num_days);

vars2plot = [5,6,8,13,17];  % indices of variables to plot
plot(time,squeeze(bio(:,end,vars2plot)),'linewidth',2);hold on
set(gca,'xlim',[0 num_days],'xtick',0:10:num_days,'xticklabel',0:10:num_days);
set(gca,'ylim',[0 30]);
xlabel('Days','fontsize',12);

namelist={'NO_3','NH_4','DON','Det_N','Bac','Phyto','Zoo_S',...
                  'ChlS','DetC','DOC','Zoo_L','SiOH_4','Dia','LDet_N',...
                  'LDet_C','LDet_S','ChlD'};
legend(namelist(vars2plot),'Location','Best');
end