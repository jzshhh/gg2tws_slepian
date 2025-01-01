%% Constant variable
global Aq  Gq gq pw lmax;
Aq=6371000; Gq=6.67259*10^-11;gq=9.8242;pw =1000; lmax=70;
% basis_degree=1:42;
basis_degree=[];

%% GPS data
timespan = [20090101 20230630];
udir_gps='data/gps';form='*.up';site_info_file='data/sites.all';
%% GRACE data
udir_grace='data/grace';mascon_file='CSR_GRACE_GRACE-FO_RL0602_Mascons_all-corrections.nc'; grid_size=2.5;

%% study area
regionName = 'AUS';
integral_edge_file= [pwd '/data/AUS_border_L2.txt']; %Expanded boundary
long_range = 110:0.5:155;
lati_range =-45:0.5:-8;

%% inversion
cutoff = 0.1;
Gauss_radius=250;

%% Plotting and Statistics
boundary_file=[pwd '/data/AUS_border_L1.txt'];
[long_bou, lat_bou]=textread(boundary_file,'%f %f','commentstyle','shell');
boundary=[long_bou lat_bou];

%%
Initialization(long_range,lati_range,integral_edge_file,regionName);

