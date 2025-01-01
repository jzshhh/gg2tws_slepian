clear;clc;

%% Load original boundary data
border=load('AUS_border_L1.txt');

%% Generate a buffer with any radiuses for study area
d=1.0; % 1.5 degree buffer
polyout = polybuffer([border(:,1) border(:,2)],'lines',d);
out = inpolygon(polyout.Vertices(:,1),polyout.Vertices(:,2),border(:,1),border(:,2));
edge_points=[polyout.Vertices(~out,1) polyout.Vertices(~out,2)];

d=0.5; % 1.5 degree buffer
polyout = polybuffer([border(:,1) border(:,2)],'lines',d);
out = inpolygon(polyout.Vertices(:,1),polyout.Vertices(:,2),border(:,1),border(:,2));
edge_points1=[polyout.Vertices(~out,1) polyout.Vertices(~out,2)];

sta_info=importdata('sites.all');


%% Plot original and extended boundary
figure('color',[1 1 1])
plot(border(:,1),border(:,2));
hold on
scatter(sta_info.data(:,1)+360,sta_info.data(:,2),20,'r','filled');
hold on
plot(edge_points(:,1),edge_points(:,2))
hold on 
hold on
plot(edge_points1(:,1),edge_points1(:,2))
set(gca,'xlim',[110 158],'ylim',[-45 -8]);

%% Manually extract extended boundaries
num=0;
while 1
    if strcmpi(get(gcf,'CurrentCharacter'),'q') % Finish when enter 'q'
        break;
    end
    num=num+1;
    aa(num,:) = ginput(1);
    hold on
    plot(aa(num,1),aa(num,2),'o-')
    hold on
end
aa=[aa; aa(1,:)];

save AUS_border_L2.txt aa -ascii % save extended boundary data
