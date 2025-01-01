clear; clc;
clearvars; close all;
code_dir = [pwd '\code'];
disp('Setting Paths...');
addpath(genpath(code_dir));
load_scenario;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load GPS and GRACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Loading Data...');
[GPS_sites,GPS_positions,GPS_times,GPS_data]=load_GPS(udir_gps,form,site_info_file,timespan);
[GRACE_sites,GRACE_positions,GRACE_data,GNSS_VCD_farall]=load_GRACE(udir_grace,mascon_file,boundary_file,grid_size,GPS_times,GPS_sites,GPS_positions);
disp('Data Loaded.');

save result/gnss_grace_data.mat 
load result/gnss_grace_data.mat 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DECOMPOSITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Decompositing');
all_data=[GPS_data-GNSS_VCD_farall GRACE_data];
all_data=all_data-mean(all_data,'omitnan');
all_positions=[GPS_positions;GRACE_positions];
all_sites=[GPS_sites;GRACE_sites];

[space_func,time_func,each_pc_explained,all_pc_explained,data_recon,rmse_all]=gps_vbpca(all_data,10);
disp('Decomposition Complete');

save result/pca_decomposition.mat 
load result/pca_decomposition.mat 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating eigenvalues and eigenvectors of D matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating eigenvalues and eigenvectors');
cd('result');
[eigenvalues,eigenvectors]=calc_eigen_paraments_main(integral_edge_file,regionName);
cd('..');
disp('Eigen prarameters calculated!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Starting Inversion');
gnss2tws_main(all_positions,space_func,eigenvalues,eigenvectors,cutoff,Gauss_radius,GPS_times,all_sites,time_func);
disp('Inversion Finished');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot/display information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Building Predictions');
[POS_grid,TWS_grid]=plot_gnss2tws_main(eigenvalues,eigenvectors,lmax,basis_degree,regionName,long_range,lati_range,cutoff,GPS_times,all_positions,boundary,time_func,space_func);
TWS_grid=TWS_grid-mean(TWS_grid);

save result/gps_ewh.mat TWS_grid GPS_times POS_grid;
save Inversion_TWS.mat;

% close all;


