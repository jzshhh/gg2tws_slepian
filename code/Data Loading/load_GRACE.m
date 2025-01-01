function [GRACE_sites,GRACE_positions,GRACE_data,GNSS_vcd_far]=load_GRACE(udir,mascon_file,boundary_file,grace_alph,gnss_time_mon,GPS_sites,GPS_positions)
%
% Description: Calculate vertical displacements of GRACE virtual stations and remove far-field effects at GNSS stations due to significant coastal signals in Australia.
%
% Input:
%   udir             File directory that stores the GRACE data
%   mascon_file      CSR mascon file
%   boundary_file    Study period, Start and end dates, e.g., timespan=[20060101 20201231]
%   grace_alph       Grid size for generating GRACE virtual stations
%   gnss_time_mon    GNSS datetime
%   GPS_sites        GNSS sites
%   GPS_positions    Location of GNSS sites
% Output:
%   GRACE_sites      Site list for GRACE virtual stations
%   GRACE_positions  Longitude and latitude coordinates of GRACE virtual stations
%   GRACE_data       GRACE-predicted VCD time series at GRACE virtual stations
%   GNSS_vcd_far     Far-field contributions at GNSS stations
%
% Author: Zhongshan Jiang
% Date: 01/01/2025 
% Institution: Sun Yat-Sen University 
% E-mail: jiangzhsh@mail.sysu.edu.cn

if ~exist(fullfile(udir,mascon_file),'file')
    fullURL=['http://download.csr.utexas.edu/outgoing/grace/RL06_mascons/' mascon_file];
    filename=['data/grace/' mascon_file];
    websave(filename,fullURL);
end

time = ncread(fullfile(udir,mascon_file),'time')+datenum([2002 01 01])-1;
grace_time = str2num(datestr(double(time),'yyyymmdd'));
ewh = ncread(fullfile(udir,mascon_file),'lwe_thickness');

for i=1:size(ewh,1)
    for j=1:size(ewh,2)
        ewh0=mean(ewh(i,j,:));
        ewh(i,j,:)=ewh(i,j,:)-ewh0;
    end
end

lon_all=ncread(fullfile(udir,mascon_file),'lon');
lat_all=ncread(fullfile(udir,mascon_file),'lat');

%% Rectangular area
boundary=load(boundary_file);
lon1=floor(min(boundary(:,1)));lon2=ceil(max(boundary(:,1)));
lat1=floor(min(boundary(:,2)));lat2=ceil(max(boundary(:,2)));

alph=0.5;
[lon3,lat3]=meshgrid(lon1:alph:lon2,lat1:alph:lat2);
lon4=reshape(lon3,[],1);lat4=reshape(lat3,[],1);

in=inpolygon(lon4,lat4,boundary(:,1),boundary(:,2));
area_grid=[lon4(in) lat4(in)];

grace_ewh_tmp=nan(size(grace_time,1),size(area_grid,1));
for i=1:size(grace_time)
    grace_ewh_tmp(i,:)=double(interp2(lon_all,flipud(lat_all),rot90(ewh(:,:,i)/100,1),area_grid(:,1),area_grid(:,2)))';
end

time1=datenum(datevec(num2str(grace_time),'yyyymmdd'));
time2=datenum(datevec(num2str(gnss_time_mon),'yyyymmdd'));
% grace_miss_month=[20020615 20020715 20030615 20110115 20110615 20120515 20121015 20130315 20130815 20130915 20140215  ...
%     20140715 20141215 20150615 20151015 20151115 20160415 20160915 20161015 20170215 20170715 20170815  ...
%     20170715 20180915 20171015 20171115 20171215 20180115 20180215 20180315 20180415 20180515 20180915];

grace_long_gap=[20170715 20170815 20170915 20171015 20171115 20171215 20180115 20180215 20180315 20180415 20180515];

grace_ewh=nan(length(time2),size(area_grid,1));
for i=1:size(area_grid,1)
    grace_ewh(:,i)=interp1(time1,grace_ewh_tmp(:,i),time2,'pchip'); %%pchip
    grace_ewh(:,i)=grace_ewh(:,i)-nanmean(grace_ewh(:,i));
end
[~,ok1,ok2]=intersect(gnss_time_mon,grace_long_gap);
grace_ewh(ok1,:)=NaN;


[lon5,lat5]=meshgrid(lon1:grace_alph:lon2,lat1:grace_alph:lat2);
lon6=reshape(lon5,[],1);lat6=reshape(lat5,[],1);

in=inpolygon(lon6,lat6,boundary(:,1),boundary(:,2));

GRACE_positions=[lon6(in) lat6(in)];

for i=1:size(GRACE_positions,1)
    GRACE_sites(i,:)=sprintf('P%03d',i);
end

constants.Aq=6371000;
constants.Gq=6.67259*10^-11;
constants.gq=9.8242;
constants.pw =1000;

Green=is_compute_greens(constants,GRACE_positions,area_grid,GRACE_sites,alph,'grace');

GRACE_data=(Green*grace_ewh'*1000)';

%%% remove contribution from far-field loads.

%load the mask file (global or mainland China), 'mask' has been used in code of Line 101
mask=rot90(ncread('code/Data Loading/CSR_GRACE_GRACE-FO_RL06_Mascons_v02_LandMask.nc','LO_val'));

alph_out1=0.5; % grid size in buffer 0°--5°
alph_out2=1.0; % grid size in buffer 5°--20°
alph_out3=2.0; % grid size in buffer > 20°

%% Generate buffer of boundery
buf_r1=5;
polyout = polybuffer([boundary(:,1) boundary(:,2)],'lines',buf_r1);
out = inpolygon(polyout.Vertices(:,1),polyout.Vertices(:,2),boundary(:,1),boundary(:,2));
edge_points5=[polyout.Vertices(~out,1) polyout.Vertices(~out,2)];

buf_r2=20;
polyout = polybuffer([boundary(:,1) boundary(:,2)],'lines',buf_r2);
out = inpolygon(polyout.Vertices(:,1),polyout.Vertices(:,2),boundary(:,1),boundary(:,2));
edge_points20=[polyout.Vertices(~out,1) polyout.Vertices(~out,2)];

%% Generate global 0.5 grids
long_range_out1=0.25:alph_out1:359.75;
lat_range_out1=89.75:-alph_out1:-89.75;

[along_out1, alat_out1]=meshgrid(long_range_out1,lat_range_out1);
blong_out1=reshape(along_out1',[],1);
blat_out1=reshape(alat_out1',[],1);

%% Generate global 1.0 grids
long_range_out2=0.5:alph_out2:359.5;
lat_range_out2=89.5:-alph_out2:-89.5;

[along_out2, alat_out2]=meshgrid(long_range_out2,lat_range_out2);
blong_out2=reshape(along_out2',[],1);
blat_out2=reshape(alat_out2',[],1);

%% Generate global 2.0 grids
long_range_out3=1:alph_out3:359;
lat_range_out3=89:-alph_out3:-89;

[along_out3, alat_out3]=meshgrid(long_range_out3,lat_range_out3);
blong_out3=reshape(along_out3',[],1);
blat_out3=reshape(alat_out3',[],1);

%% Obtain grids outside study area
% get the grids from boundery to 5°buffer
is_area_grid = inpolygon(blong_out1,blat_out1,boundary(:,1),boundary(:,2));
mid=[blong_out1(~is_area_grid) blat_out1(~is_area_grid)];
is_area_grid=inpolygon(mid(:,1),mid(:,2),edge_points5(:,1),edge_points5(:,2));
area_grid_1=mid(is_area_grid,:);
% get the grids from 5°buffer to 20°buffer
is_area_grid = inpolygon(blong_out2,blat_out2,edge_points5(:,1),edge_points5(:,2));
mid=[blong_out2(~is_area_grid) blat_out2(~is_area_grid)];
is_area_grid=inpolygon(mid(:,1),mid(:,2),edge_points20(:,1),edge_points20(:,2));
area_grid_2=mid(is_area_grid,:);
% get the grids from 20°buffer to global
is_area_grid = inpolygon(blong_out3,blat_out3,edge_points20(:,1),edge_points20(:,2));
area_grid_3=[blong_out3(~is_area_grid) blat_out3(~is_area_grid)];

%% Obtain ewh at grids outside study area
global_ewh=rot90(ewh/100,1).*mask; %% convert cm to m
for i=1:size(global_ewh,3)
    ewh_outside1(:,i)=interp2(lon_all,lat_all,global_ewh(:,:,i),area_grid_1(:,1),area_grid_1(:,2))';
    ewh_outside2(:,i)=interp2(lon_all,lat_all,global_ewh(:,:,i),area_grid_2(:,1),area_grid_2(:,2))';
    ewh_outside3(:,i)=interp2(lon_all,lat_all,global_ewh(:,:,i),area_grid_3(:,1),area_grid_3(:,2))';
end

%% calculate GFs
green1=is_compute_greens(constants,GPS_positions,area_grid_1,GPS_sites,alph_out1,'GRACE1');
green2=is_compute_greens(constants,GPS_positions,area_grid_2,GPS_sites,alph_out2,'GRACE2');
green3=is_compute_greens(constants,GPS_positions,area_grid_3,GPS_sites,alph_out3,'GRACE3');

%% calculate 3-D (n,e,u) displacement at GPS sites due to water loads outside study area
GNSS_VCD_farall=green1*ewh_outside1*10^3+green2*ewh_outside2*10^3+green3*ewh_outside3*10^3; % all far field-induced displacement
grace_decyr = decyear(datevec(num2str(grace_time),'yyyymmdd'));
gnss_decyr = decyear(datevec(num2str(gnss_time_mon),'yyyymmdd'));

GNSS_vcd_far=NaN(length(gnss_time_mon),size(GPS_positions,1));
for i=1:size(GNSS_VCD_farall,1)
    GNSS_vcd_far(:,i)=vcd_due_far([grace_decyr GNSS_VCD_farall(i,:)'],gnss_decyr,ok1);

end
end

function grace_vcd=vcd_due_far(ts_data,ts2,index)
time=ts_data(:,1);data=ts_data(:,2);
A=[ones(size(data,1),1) time-time(1) cos(2*pi*time) sin(2*pi*time) cos(4*pi*time) sin(4*pi*time)];
B=data;
coeff=lscov(A,B);
data2=coeff(1)+(ts2-ts2(1))*coeff(2)+coeff(3)*cos(2*pi*ts2)+coeff(4)*sin(2*pi*ts2)+coeff(5)*cos(4*pi*ts2)+coeff(6)*sin(4*pi*ts2);
grace_vcd=interp1(time,data,ts2,'pchip'); %%pchip
grace_vcd(index)=data2(index);
end




