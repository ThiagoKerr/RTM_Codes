tic; clc; clear all; close all

% Pasta para as funções
addpath('./Funcoes'); addpath('./Dados')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constantes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Elipsoide GRS80
a=6378137.0; b=6356752.3141; e2 =0.00669438002290; f=1/298.257222101; m=0.00344978600308; rad=pi/180; G = 6.67259E-11;   % m3 kg-1 s-2
gama_eq = 9.7803267715; gama_po = 9.8321863685;  % gravidade normal no equador e no polo, respectivamente (m/s^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Importação das imagens %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelo_densidade='nao';

H_detailed = double(imread('BRAZ_detalhadoCop.tif'));
H_reference = double(imread('BRAZ_referenciaCop.tif'));

H_detailed(1,:)=[];
H_detailed(:,1)=[];
H_reference(1,:)=[];
H_reference(:,1)=[];

H=H_detailed-H_reference;
[lin, col]=size(H_detailed);

% Correções de densidade
if strcmp('sim',modelo_densidade)
densidade = double(imread('BRAZ_densidades.tif'))*1000;

%CORREÇÃO DOS VALORES NAN DO MODELO DE DENSIDADES
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
res=0.0002777777777777779402; res_rad=res*rad;   % Resolução em ° e em rad, respectivamente
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
col_CP=find(abs(long-gradeMDA_long(1,:))==min(min(abs(long-gradeMDA_long(1,:)))));   % coluna do MDS detalhado que contém o ponto de computaçãogama_elipsoide = (a*gama_eq*cos(lat)*cos(lat) + b*gama_po*sin(lat)*sin(lat))/sqrt(a^2*cos(lat)*cos(lat) + b^2*sin(lat)*sin(lat));
gama_elipsoide = (a*gama_eq*cos(lat)*cos(lat) + b*gama_po*sin(lat)*sin(lat))/sqrt(a^2*cos(lat)*cos(lat) + b^2*sin(lat)*sin(lat));
gama = gama_elipsoide*(1-2/a*(1+f+m-2*f*sin(lat)*sin(lat))*h+3*h^2/a^2);
Nell_cp=a./sqrt(1 - (e2.*sin(lat).*sin(lat)));
r=Nell_cp+h;

% Correção harmônica
if H(lin_CP,col_CP)<0
H_reference_cp=interp2(gradeMDA_long,gradeMDA_lat,H_reference,long,lat,'cubic');         % height of computation point on reference surface
% H_detailed_cp=interp2(gradeMDA_long,gradeMDA_lat,H_detailed,long,lat,'cubic');         % height of computation point on detailed surface
if strcmp('sim',modelo_densidade)
densidade_cp=interp2(gradeMDA_long,gradeMDA_lat,densidade,long,lat,'cubic');         % height of computation point on reference surface
Eta_RTM_HC_Omang=(-4*pi*G*densidade_cp*(h-H_reference_cp)^2)/gama   % m
else
Eta_RTM_HC_Omang=(-4*pi*G*rho*(h-H_reference_cp)^2)/gama   % m
end
else
Eta_RTM_HC_Omang=0;
end

clear Nell_cp estacao_RTM h H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%% Cálculo das anomalias de altura RTM no entorno da estação de cálculo %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRISMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat_o=gradeMDA_lat(:); long_o=gradeMDA_long(:); H_detailed_o=H_detailed(:); H_reference_o=H_reference(:);
Nell_o=a./(1 - (e2.*sin(lat_o).*sin(lat_o))).^0.5;
zzp=H_detailed_o+Nell_o;
[x, y, z]=global2local(long,lat,r,long_o,lat_o,zzp); % coordinates of computation point in prism-based Cartesian coordinate system
lo_prism=(x.^2+y.^2+z.^2).^0.5; lo=(x.^2+y.^2+z.^2).^0.5;

% Raio de integração
lo(lo>radius_metros)=NaN;
index=~isnan(lo); index=double(index);
index(index==0)=NaN;

lo_prism=lo_prism.*index; lo_prism(isnan(lo_prism))=[];
long_o_prism=long_o.*index; long_o_prism(isnan(long_o_prism))=[];
lat_o_prism=lat_o.*index; lat_o_prism(isnan(lat_o_prism))=[];
H_detailed_o_prism=H_detailed_o.*index; H_detailed_o_prism(isnan(H_detailed_o_prism))=[];
H_reference_o_prism=H_reference_o.*index; H_reference_o_prism(isnan(H_reference_o_prism))=[];
Nell_o_prism=Nell_o.*index; Nell_o_prism(isnan(Nell_o_prism))=[];

clear H_detailed H_reference H_reference_o H_reference_cp H_detailed_o lat_o gradeMDA_lat gradeMDA_long lo long_o Nell_o

xx1=long_o_prism-res_rad/2;                       % canto inferior esquerdo do pixel
yy1=lat_o_prism-res_rad/2;          % canto inferior esquerdo do pixel
xx2=long_o_prism+res_rad/2;          % canto superior direito do pixel
yy2=lat_o_prism+res_rad/2;                       % canto superior direito do pixel
zz1=H_detailed_o_prism+Nell_o_prism;
zz2=H_reference_o_prism+Nell_o_prism;

if strcmp('sim',modelo_densidade)
densidade=densidade(:); densidade=densidade.*index; densidade(isnan(densidade))=[];
V=Prism_nosum(xx1,xx2,yy1,yy2,zz1,zz2,long,lat,r,densidade);  % m
else
V=Prism_nosum(xx1,xx2,yy1,yy2,zz1,zz2,long,lat,r,rho);  % m
end
Eta_RTM_prisma=sum(V,"omitnan")/gama+Eta_RTM_HC_Omang   % m
toc

% dg por fatias de raios para geração de gráfico
Eta_RTM_prisma_raios=zeros(1,radius_metros/1000);
radius_grafico=zeros(1,radius_metros/1000);
for i=1:radius_metros/1000
lo_prism_grafico=lo_prism; lo_prism_grafico(lo_prism_grafico>(i*1000))=NaN;
index=~isnan(lo_prism_grafico);
V_grafico=V.*index;   % m
Eta_RTM_prisma_raios(1,i)=sum(V_grafico,"omitnan")/gama+Eta_RTM_HC_Omang;   % m
radius_grafico(1,i)=i*1000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Geração de imagens e gráficos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1);
% plot(radius_grafico,Eta_RTM_prisma_raios);
% xlabel('Radius of integration (km)','FontName','Times New Roman','FontSize',14)
% ylabel('RTM height anomaly (mGal)','FontName','Times New Roman','FontSize',14)
% set(gca,'FontName','Times New Roman','FontSize',14)
Eta_RTM_prisma_raios=Eta_RTM_prisma_raios';