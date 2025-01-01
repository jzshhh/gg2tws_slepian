function [north,east,up]=disk_greens_3d(az,thet,lmax,love_h,love_l,alph_disk,constants)

Pot=Potential(lmax+1,alph_disk,constants);
[Pl,Pl1]=pLegendre(lmax,thet);
up=sum(love_h(1:lmax+1).*Pot(1:lmax+1).*Pl(1:lmax+1))/constants.gq;
away=sum(love_l(1:lmax+1).*Pot(1:lmax+1).*Pl1(1:lmax+1))/constants.gq;
east = -sind(az)*away;
north = -cosd(az)*away;
end

function Pot=Potential(n,alph_disk,constants)
[Pl,~]=pLegendre(n,alph_disk);
%% 重力位扰动
mass=disk_mass(alph_disk,constants);
Pot=NaN(n+1,1);
Pot(1,1)=mass*constants.Gq/constants.Aq;

i=2:n;
Pot(i,1)=mass*constants.Gq.*(Pl(i-1)-Pl(i+1))./(constants.Aq*(1-cosd(alph_disk))*(2*(i-1)+1))';
end

function m=disk_mass(alph,constants)
m=1*constants.pw*2*pi*(1-cosd(alph))*constants.Aq^2;
end

function [P,P1]=pLegendre(n,deg)
myx=cosd(deg);
mydx=-1*sind(deg);
nn=0:1:n+1;
P=zeros(length(nn),1);
P1=zeros(length(nn),1);
for jj=1:length(nn)
    myn=nn(jj);
    if myn==0
        P(jj)=1;
        P1(jj)=0;
    elseif myn==1
        P(jj)=myx;
        P1(jj)=mydx;
    else
        P(jj)=((2*myn-1)/myn)*myx*P(jj-1)-((myn-1)/myn)*P(jj-2);
        P1(jj)=myn*(P(jj)-myx*P(jj-1))./((1-myx).*(1+myx)).*sqrt(1-myx.^2);
    end
end
end
