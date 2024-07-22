tic; clc; clear all; close all

% Pasta para as funções
addpath('./Funcoes'); addpath('./Dados')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constantes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Elipsoide GRS80
a=6378137.0; b=6356752.3141; e2 =0.00669438002290; f=1/298.257222101; m=0.00344978600308; rad=pi/180; G = 6.67259E-11;   % m3 kg-1 s-2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Importação das imagens %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelo_densidade='nao';

H_detailed = double(imread('BRAZ_detalhado.tif'));
H_reference = double(imread('BRAZ_referencia.tif'));

H_detailed(1,:)=[];
H_detailed(:,1)=[];
H_reference(1,:)=[];
H_reference(:,1)=[];

H=H_detailed-H_reference;
[lin, col]=size(H_detailed);

% Correções de densidade
if strcmp('sim',modelo_densidade)
densidade = double(imread('BRAZ_densidades.tif'))*1000;

densidade(1,:)=[];
densidade(:,1)=[];
densidade(densidade<0)=2670;

densidade_media=mean(nonzeros(densidade));
RET=1-1030/densidade_media;
densidade(densidade==2670)=densidade_media;
densidade_RET=densidade(:); H_detailed=H_detailed(:); H_reference=H_reference(:);
densidade_RET(densidade_RET==0)=NaN;
index_1=~isnan(densidade_RET);   % área oceânica zero e área continetal um
index_2=isnan(densidade_RET);    % área oceânica um e área continetal zero
H_detailed=H_detailed.*index_1+H_detailed.*index_2.*RET; % soma de zero na parte oceânica e valor na parte continetal com zero na parte continental e valor*RET na parte oceânica
H_reference=H_reference.*index_1+H_reference.*index_2.*RET; % soma de zero na parte oceânica e valor na parte continetal com zero na parte continental e valor*RET na parte oceânica
H_detailed=reshape(H_detailed,lin,col);
H_reference=reshape(H_reference,lin,col);
densidade(densidade==0)=densidade_media;
else
rho=2670;   % kg/m3
end

% Resolução
res=0.0002777777777777779402; res_rad=res*rad;     % Resolução em ° e em rad, respectivamente
radius_metros=210000;                   % raio de integração em metros

% Coordenadas geodésicas do canto superior esquerdo do primeiro pixel da imagem no canto NW - superior esquerdo - depois passando para o centro geométrico do pixel
minlat=(-13.2468055555555431-res)*rad-(res_rad/2);    % latN
minlon=(-50.5781944444444491+res)*rad+(res_rad/2);   % longW

% Grade de coordenadas geodésicas
gradeMDA_lat=(minlat:-res_rad:(minlat-(lin-1)*res_rad))'.*ones(1,col);
gradeMDA_long=(minlon:res_rad:(minlon+(col-1)*res_rad)).*ones(lin,1);

clear minlat minlon index_1 index_2 densidade_RET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%% Coordenadas geodésicas, cartesianas e locais do ponto de computação e gravidade normal %%%%%%%%%%%%%%%%%

estacao_RTM=load('BRAZ.txt');
lat=deg2rad(estacao_RTM(1));    long=deg2rad(estacao_RTM(2));    h=estacao_RTM(3);
lin_CP=find(abs(lat-gradeMDA_lat(:,1))==min(min(abs(lat-gradeMDA_lat(:,1)))));       % linha do MDS detalhado que contém o ponto de computação
col_CP=find(abs(long-gradeMDA_long(1,:))==min(min(abs(long-gradeMDA_long(1,:)))));   % coluna do MDS detalhado que contém o ponto de computação
Nell_cp=a./sqrt(1 - (e2.*sin(lat).*sin(lat)));
r=Nell_cp+h;  r_2=r.*r; r_3=r.*r.*r;

sin_lat=sin(lat); cos_lat=cos(lat); sin_long=sin(long); cos_long=cos(long);

% Correção harmônica
if H(lin_CP,col_CP)<0
H_reference_cp=interp2(gradeMDA_long,gradeMDA_lat,H_reference,long,lat,'cubic');         % height of computation point on reference surface
% H_detailed_cp=interp2(gradeMDA_long,gradeMDA_lat,H_detailed,long,lat,'cubic');         % height of computation point on detailed surface
if strcmp('sim',modelo_densidade)
densidade_cp=interp2(gradeMDA_long,gradeMDA_lat,densidade,long,lat,'cubic');         % density of computation point
dg_RTM_HC=(4*pi*G*densidade(lin_CP,col_CP)*(h-H_reference_cp))*1D5   % mGal
else
dg_RTM_HC=(4*pi*G*rho*(h-H_reference_cp))*1D5   % mGal
end
else
dg_RTM_HC=0;
end
clear Nell_cp estacao_RTM  H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%% Cálculo dos distúrbios de gravidade RTM no entorno da estação de cálculo %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESSEROIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_lat=res_rad; delta_long=res_rad; cte=G*delta_lat*delta_long;
lat_o=gradeMDA_lat(:); long_o=gradeMDA_long(:); H_detailed_o=H_detailed(:); H_reference_o=H_reference(:);
Nell_o=a./(1 - (e2.*sin(lat_o).*sin(lat_o))).^0.5;
zzp=H_detailed_o+Nell_o;
[x, y, z]=global2local(long,lat,r,long_o,lat_o,zzp); % coordinates of computation point in tesseroid-based Cartesian coordinate system
lo=(x.^2+y.^2+z.^2).^0.5;

% Raio de integração
lo(lo>radius_metros)=NaN; index=~isnan(lo); index=double(index); index(index==0)=NaN;

long_o_tesseroide=long_o.*index; long_o_tesseroide(isnan(long_o_tesseroide))=[];
lat_o_tesseroide=lat_o.*index; lat_o_tesseroide(isnan(lat_o_tesseroide))=[];
H_detailed_o_tesseroide=H_detailed_o.*index; H_detailed_o_tesseroide(isnan(H_detailed_o_tesseroide))=[];
H_reference_o_tesseroide=H_reference_o.*index; H_reference_o_tesseroide(isnan(H_reference_o_tesseroide))=[];
Nell_o_tesseroide=Nell_o.*index; Nell_o_tesseroide(isnan(Nell_o_tesseroide))=[];

sin_lat_o=sin(lat_o_tesseroide); cos_lat_o=cos(lat_o_tesseroide); sin_long_o=sin(long_o_tesseroide); cos_long_o=cos(long_o_tesseroide);
sin_dlong=sin_long_o.*cos_long-cos_long_o.*sin_long; cos_dlong=cos_long_o.*cos_long+sin_long_o.*sin_long;
r1=H_reference_o_tesseroide+Nell_o_tesseroide;
r2=H_detailed_o_tesseroide+Nell_o_tesseroide;
delta_r=r2-r1;
ro=(r1+r2)./2; ro_2=ro.*ro;  ro_3=ro.*ro.*ro;

cos_psi_o=sin_lat_o.*sin_lat + cos_lat_o.*cos_lat.*cos_dlong; sin_psi_o=real((1-cos_psi_o.*cos_psi_o).^0.5);

lo_tesseroide=real((r_2+ro_2-2.*r.*ro.*cos_psi_o).^0.5);
lo_tesseroide_2=lo_tesseroide.*lo_tesseroide;
lo_tesseroide_3=lo_tesseroide.*lo_tesseroide.*lo_tesseroide;
lo_tesseroide_4=lo_tesseroide.*lo_tesseroide.*lo_tesseroide.*lo_tesseroide;
lo_tesseroide_5=lo_tesseroide.*lo_tesseroide.*lo_tesseroide.*lo_tesseroide.*lo_tesseroide;
lo_tesseroide_7=lo_tesseroide.*lo_tesseroide.*lo_tesseroide.*lo_tesseroide.*lo_tesseroide.*lo_tesseroide.*lo_tesseroide;

L000=(ro_2.*(r-ro.*cos_psi_o).*cos_lat_o)./(lo_tesseroide_3);
L002=(ro./lo_tesseroide).^3.*cos_lat.*cos_lat_o.*cos_lat_o.* (cos_dlong-((3.*r)./(lo_tesseroide_2)).*(2.*ro.*cos_lat.*cos_lat_o.*sin_dlong.*sin_dlong+(r-ro.*cos_psi_o).*cos_dlong)+((15.*r_2.*ro)./(lo_tesseroide_4)).*cos_lat.*cos_lat_o.*(r-ro.*cos_psi_o).*sin_dlong.*sin_dlong);
L200=(r.*cos_lat_o./lo_tesseroide_3).*(2-(3.*ro./lo_tesseroide_2).*(5.*ro-(2.*r+3.*ro.*cos_psi_o).*cos_psi_o)+(15.*ro_3./lo_tesseroide_4).*sin_psi_o.*sin_psi_o.*(ro-r.*cos_psi_o));
L020=(ro./lo_tesseroide).^3.*cos_lat.*(1-2.*sin_lat_o.*sin_lat_o).*cos_dlong+...
    ((ro_2)./(lo_tesseroide_5)).*(-r.*(r_2+ro_2).*cos_lat_o+...
    ro.*sin_lat.*(-r.*ro.*(sin_lat.*cos_lat_o-cos_lat.*sin_lat_o.*cos_dlong)+...
    sin_lat_o.*cos_lat_o.*(2.*r_2+4.*ro_2-3.*r.*ro.*sin_lat.*sin_lat_o))+...
      ro_2.*cos_lat.*cos_dlong.*(1-2.*sin_lat_o.*sin_lat_o).*...
      (ro+r.*cos_lat.*cos_lat_o.*cos_dlong)+...
      r.*ro_2.*cos_lat.*sin_lat_o.*cos_lat_o.*cos_dlong.*...
      (3.*sin_lat.*cos_lat_o-4.*cos_lat.*sin_lat_o.*cos_dlong))+...
	  ((5.*r.*ro_3)./(lo_tesseroide_7)).*(-r.*(r_2+ro_2).*sin_lat_o+...
      ro_2.*cos_lat.*sin_lat_o.*cos_lat_o.*cos_dlong.*...
      (ro+r.*cos_lat.*cos_lat_o.*cos_dlong)+...
      ro.*sin_lat.*(2.*r_2-ro_2-r.*ro.*cos_psi_o+sin_lat_o.*sin_lat_o.*...
      (r_2+2.*ro_2-r.*ro.*sin_lat.*sin_lat_o))).*...
      (sin_lat.*cos_lat_o-cos_lat.*sin_lat_o.*cos_dlong);

if strcmp('sim',modelo_densidade)
densidade_o_tesseroide=densidade(:); densidade_o_tesseroide=densidade_o_tesseroide.*index; densidade_o_tesseroide(isnan(densidade_o_tesseroide))=[];
Vz_tesseroide=densidade_o_tesseroide.*cte.*delta_r.*(L000+(1/24).*(L200.*delta_r.^2+L020.*delta_lat.^2+L002.*delta_long.^2));  % m/s2
else
Vz_tesseroide=rho.*cte.*delta_r.*(L000+(1/24).*(L200.*delta_r.^2+L020.*delta_lat.^2+L002.*delta_long.^2));  % m/s2
end
dg_RTM_tesseroide=sum(Vz_tesseroide,"omitnan").*1D5+dg_RTM_HC   % mGal
toc

% dg por fatias de raios para geração de gráfico
dg_RTM_tesseroide_raios=zeros(1,radius_metros/1000);
radius_grafico=zeros(1,radius_metros/1000);
for i=1:radius_metros/1000
lo_tesseroide_grafico=lo_tesseroide; lo_tesseroide_grafico(lo_tesseroide_grafico>(i*1000))=NaN;
index=~isnan(lo_tesseroide_grafico);
Vz_grafico=Vz_tesseroide.*index;    %m/s2
dg_RTM_tesseroide_raios(1,i)=sum(Vz_grafico,"omitnan").*1D5+dg_RTM_HC;   % mGal
radius_grafico(1,i)=i*1000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Geração de imagens e gráficos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
plot(radius_grafico,dg_RTM_tesseroide_raios);
xlabel('Radius of integration (km)','FontName','Times New Roman','FontSize',14)
ylabel('RTM gravity disturbance (mGal)','FontName','Times New Roman','FontSize',14)
set(gca,'FontName','Times New Roman','FontSize',14)
dg_RTM_tesseroide_raios=dg_RTM_tesseroide_raios';