function Green=compute_greens(constants,gnss_lonlat,grid_lonlat,sites,alph,data_type)
%
% Description: Calculating Green's functions
%
% Input:
%   gnss_lonlat          Positions (lon, lat) of GNSS sites
%   grid_lonlat          Positions (lon, lat) of grids
%   sites                GNSS sites
%   alph                 Grid interval, e.g., 0.25
% Output:
%   Green                Green's function matrix
%
% Author: Zhongshan Jiang
% Date: 28/10/2021 
% Institution: Southwest Jiaotong University 
% E-mail: jzshhh@my.swjtu.edu.cn

disp('Computing Greens function!');
gpslon=gnss_lonlat(:,1);gpslat=gnss_lonlat(:,2);
loadlon=grid_lonlat(:,1);loadlat=grid_lonlat(:,2);

m=length(gpslon);n=length(loadlon);
Green=zeros(m,n)*NaN;
[love_h,love_l,love_k]=load_LLNs;
h=waitbar(0,'Calculate Green''s function...');
for i=1:m % m: number of sites
    for j=1:n % n: number of grids
        [sinang,cosang,az] = angrad(loadlat(j),loadlon(j),gpslat(i),gpslon(i));
        thet=acosd(cosang);
        alph_disk=rad2deg(sqrt((deg2rad(alph)*constants.Aq)^2*cosd(loadlat(j))/pi)/constants.Aq); %Applying area-weighted scheme to calculate disk radius
        % [ndis,edis,udis]=disk_greens_3d(az,thet,20000,love_h,love_l,alph_disk,constants);
        % Green(3*(i-1)+1,j)=ndis;Green(3*(i-1)+2,j)=edis;Green(3*(i-1)+3,j)=udis;
        Green(i,j)=disk_greens_1d(thet,20000,love_h,alph_disk);
    end
    smg=['Green''s function: '  sites(i,:)];
    waitbar(i/m,h,smg)
end
close(h);
GNSS_bak=gnss_lonlat;
load_bak=grid_lonlat;
save(['result/' data_type '_LGFs.mat'],'GNSS_bak','load_bak','Green'); %Saving Green's functions
end