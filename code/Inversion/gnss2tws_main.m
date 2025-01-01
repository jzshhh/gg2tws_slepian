function gnss2tws_main(GPS_positions,GPS_data,val,vec,cutoff,Gauss_radius,GPS_times,GPS_sites,gnss_time)
%
% Description: Performing PCA decomposition using the 'als' algorithm for missing data
%
% Input:
%   GPS_positions         Raw observation data
%   GPS_data           Number of selected PCs
%   val         Raw observation data
%   vec           Number of selected PCs
%   cutoff         Raw observation data
%   Gauss_radius           Number of selected PCs
%   GPS_times           Number of selected PCs
%   GPS_sites         Raw observation data
%   gnss_time           Number of selected PCs
% 
% Author: Zhongshan Jiang
% Date: 12/09/2022
% Institution: Sun Yat-Sen University
% E-mail: jiangzhsh@mail.sysu.edu.cn

global Aq lmax;
position=GPS_positions;
clat=90-position(:,2);
lon=position(:,1);

J=find(val>=cutoff, 1, 'last' );
ll=1:J;
G=my_galpha(clat,lon,vec,ll,lmax)*Aq;

%%
ddir1='result\Ucoef.mat';%存储位移系数
ddir2='result\TWScoef.mat';%存储等效水高系数
ddir3='result\TWS_site.mat';%存储等效水高
ddir4='result\Udata.mat';%存储等效水高

%%
Data=GPS_data;
coef=datafit(G,Data);
coef1=coef*gnss_time';
save(ddir1,'coef1','GPS_times');
fprintf('Ucoef Compute Finish!\n');

%%
Upre=coef2pot(coef,vec,position,[],[]);
Upre=gnss_time*Upre';
Data=gnss_time*Data';
save(ddir4,'Upre','Data','GPS_times','GPS_positions','GPS_sites');
fprintf('Upre Compute Finish!\n');

%%
TWScoef=VCDcoef2TWScoef(coef,vec);
TWScoef=filter_slep_coff_gauss(TWScoef,vec,Gauss_radius);
save(ddir2,'TWScoef','GPS_times');
fprintf('TWScoef Compute Finish!\n');

%%
TWS=coef2pot(TWScoef,vec,position,[],val);
TWS=gnss_time*TWS';
save(ddir3,'TWS','GPS_times');
fprintf('TWS Compute Finish!\n');
