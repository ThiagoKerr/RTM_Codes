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
res=0.0002777777777777779402; res_rad=res*rad; % Resolução em ° e em rad, respectivamente
radius_metros_intervalo=500; radius_metros_total=210000; % m
radius_poliedro_total=radius_metros_total/108000; % °
radius_pixel_poliedro=ceil(radius_poliedro_total/res); % raio de integração em pixels arredondado para cima abordagem poliedro

% Coordenadas geodésicas do canto superior esquerdo do primeiro pixel da imagem no canto NW - superior esquerdo - depois passando para o centro geométrico do pixel
minlat=(-13.2468055555555431-res)*rad-(res_rad/2);    % latN
minlon=(-50.5781944444444491+res)*rad+(res_rad/2); %Lon W

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLIEDRO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Criando as subgrades dos poliedros
sub_grade_lat_NW=(gradeMDA_lat(lin_CP-radius_pixel_poliedro:lin_CP+radius_pixel_poliedro-1 , col_CP-radius_pixel_poliedro:col_CP+radius_pixel_poliedro-1));
sub_grade_lat_NE=(gradeMDA_lat(lin_CP-radius_pixel_poliedro:lin_CP+radius_pixel_poliedro-1 , col_CP-radius_pixel_poliedro+1:col_CP+radius_pixel_poliedro));
sub_grade_lat_SW=(gradeMDA_lat(lin_CP-radius_pixel_poliedro+1:lin_CP+radius_pixel_poliedro , col_CP-radius_pixel_poliedro:col_CP+radius_pixel_poliedro-1));
sub_grade_lat_SE=(gradeMDA_lat(lin_CP-radius_pixel_poliedro+1:lin_CP+radius_pixel_poliedro , col_CP-radius_pixel_poliedro+1:col_CP+radius_pixel_poliedro));

sub_grade_long_NW=(gradeMDA_long(lin_CP-radius_pixel_poliedro:lin_CP+radius_pixel_poliedro-1 , col_CP-radius_pixel_poliedro:col_CP+radius_pixel_poliedro-1));
sub_grade_long_NE=(gradeMDA_long(lin_CP-radius_pixel_poliedro:lin_CP+radius_pixel_poliedro-1 , col_CP-radius_pixel_poliedro+1:col_CP+radius_pixel_poliedro));
sub_grade_long_SW=(gradeMDA_long(lin_CP-radius_pixel_poliedro+1:lin_CP+radius_pixel_poliedro , col_CP-radius_pixel_poliedro:col_CP+radius_pixel_poliedro-1));
sub_grade_long_SE=(gradeMDA_long(lin_CP-radius_pixel_poliedro+1:lin_CP+radius_pixel_poliedro , col_CP-radius_pixel_poliedro+1:col_CP+radius_pixel_poliedro));

sub_grade_H_detailed_NW=(H_detailed(lin_CP-radius_pixel_poliedro:lin_CP+radius_pixel_poliedro-1 , col_CP-radius_pixel_poliedro:col_CP+radius_pixel_poliedro-1));
sub_grade_H_detailed_NE=(H_detailed(lin_CP-radius_pixel_poliedro:lin_CP+radius_pixel_poliedro-1 , col_CP-radius_pixel_poliedro+1:col_CP+radius_pixel_poliedro));
sub_grade_H_detailed_SW=(H_detailed(lin_CP-radius_pixel_poliedro+1:lin_CP+radius_pixel_poliedro , col_CP-radius_pixel_poliedro:col_CP+radius_pixel_poliedro-1));
sub_grade_H_detailed_SE=(H_detailed(lin_CP-radius_pixel_poliedro+1:lin_CP+radius_pixel_poliedro , col_CP-radius_pixel_poliedro+1:col_CP+radius_pixel_poliedro));

sub_grade_H_reference_NW=(H_reference(lin_CP-radius_pixel_poliedro:lin_CP+radius_pixel_poliedro-1 , col_CP-radius_pixel_poliedro:col_CP+radius_pixel_poliedro-1));
sub_grade_H_reference_NE=(H_reference(lin_CP-radius_pixel_poliedro:lin_CP+radius_pixel_poliedro-1 , col_CP-radius_pixel_poliedro+1:col_CP+radius_pixel_poliedro));
sub_grade_H_reference_SW=(H_reference(lin_CP-radius_pixel_poliedro+1:lin_CP+radius_pixel_poliedro , col_CP-radius_pixel_poliedro:col_CP+radius_pixel_poliedro-1));
sub_grade_H_reference_SE=(H_reference(lin_CP-radius_pixel_poliedro+1:lin_CP+radius_pixel_poliedro , col_CP-radius_pixel_poliedro+1:col_CP+radius_pixel_poliedro));

if strcmp('sim',modelo_densidade)
sub_grade_densidade_NW=(densidade(lin_CP-radius_pixel_poliedro:lin_CP+radius_pixel_poliedro-1 , col_CP-radius_pixel_poliedro:col_CP+radius_pixel_poliedro-1));
sub_grade_densidade_NE=(densidade(lin_CP-radius_pixel_poliedro:lin_CP+radius_pixel_poliedro-1 , col_CP-radius_pixel_poliedro+1:col_CP+radius_pixel_poliedro));
sub_grade_densidade_SW=(densidade(lin_CP-radius_pixel_poliedro+1:lin_CP+radius_pixel_poliedro , col_CP-radius_pixel_poliedro:col_CP+radius_pixel_poliedro-1));
sub_grade_densidade_SE=(densidade(lin_CP-radius_pixel_poliedro+1:lin_CP+radius_pixel_poliedro , col_CP-radius_pixel_poliedro+1:col_CP+radius_pixel_poliedro));
densidade_o_poliedro=1/4.*(sub_grade_densidade_NW(:)+sub_grade_densidade_NE(:)+sub_grade_densidade_SW(:)+sub_grade_densidade_SE(:));
end

lat_o_NW=sub_grade_lat_NW(:); lat_o_NE=sub_grade_lat_NE(:); lat_o_SW=sub_grade_lat_SW(:); lat_o_SE=sub_grade_lat_SE(:);
long_o_NW=sub_grade_long_NW(:); long_o_NE=sub_grade_long_NE(:); long_o_SW=sub_grade_long_SW(:); long_o_SE=sub_grade_long_SE(:);
H_detailed_o_NW=sub_grade_H_detailed_NW(:); H_detailed_o_NE=sub_grade_H_detailed_NE(:); H_detailed_o_SW=sub_grade_H_detailed_SW(:); H_detailed_o_SE=sub_grade_H_detailed_SE(:);
H_reference_o_NW=sub_grade_H_reference_NW(:); H_reference_o_NE=sub_grade_H_reference_NE(:); H_reference_o_SW=sub_grade_H_reference_SW(:); H_reference_o_SE=sub_grade_H_reference_SE(:);

sin_lat_o_poliedro_NW=sin(lat_o_NW); sin_lat_o_poliedro_NE=sin(lat_o_NE); sin_lat_o_poliedro_SW=sin(lat_o_SW); sin_lat_o_poliedro_SE=sin(lat_o_SE);

Nell_NW=a./(1 - (e2.*sin_lat_o_poliedro_NW.*sin_lat_o_poliedro_NW)).^0.5; Nell_SE=a./(1 - (e2.*sin_lat_o_poliedro_SE.*sin_lat_o_poliedro_SE)).^0.5; Nell_NE=a./(1 - (e2.*sin_lat_o_poliedro_NE.*sin_lat_o_poliedro_NE)).^0.5; Nell_SW=a./(1 - (e2.*sin_lat_o_poliedro_SW.*sin_lat_o_poliedro_SW)).^0.5;

clear sub_grade_lat_NW sub_grade_lat_NE sub_grade_lat_SW sub_grade_lat_SE sub_grade_long_NW sub_grade_long_NE sub_grade_long_SW sub_grade_long_SE
clear sub_grade_H_detailed_NW sub_grade_H_detailed_NE sub_grade_H_detailed_SW sub_grade_H_detailed_SE
clear sub_grade_H_reference_NW sub_grade_H_reference_NE sub_grade_H_reference_SW sub_grade_H_reference_SE
clear sin_lat_o_poliedro_NW sin_lat_o_poliedro_NE sin_lat_o_poliedro_SW sin_lat_o_poliedro_SE

clear densidade
clear gradeMDA_lat
clear gradeMDA_long
clear H_detailed
clear H_reference
clear lat_o_NE
clear lat_o_SE
clear long_o_SE
clear long_o_SW
clear sub_grade_densidade_NE
clear sub_grade_densidade_NW
clear sub_grade_densidade_SE
clear sub_grade_densidade_SW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Iterações nas camadas

coluna=0;

for k=0:radius_metros_intervalo:(radius_metros_total-radius_metros_intervalo)

radius_metros_1=0+k; radius_metros_2=radius_metros_intervalo+k; % m
coluna=coluna+1;

% define origin of polyhedron-based coordinate system
z1_gridNW=H_detailed_o_NW+Nell_NW; z1_gridSE=H_detailed_o_SE+Nell_SE; z1_gridNE=H_detailed_o_NE+Nell_NE; z1_gridSW=H_detailed_o_SW+Nell_SW;
z1_refNW=H_reference_o_NW+Nell_NW; z1_refSE=H_reference_o_SE+Nell_SE; z1_refNE=H_reference_o_NE+Nell_NE; z1_refSW=H_reference_o_SW+Nell_SW;

lon_o_poliedro=1/2.*(long_o_NW+long_o_NE); % cetro geométrico do poliedro considerando os 4 pixels vizinhos
lat_o_poliedro=1/2.*(lat_o_NW+lat_o_SW); % cetro geométrico do poliedro considerando os 4 pixels vizinhos
zz_o_poliedro=1/4.*(z1_gridNW+z1_gridSE+z1_gridNE+z1_gridSW); % cetro geométrico do poliedro considerando os 4 pixels vizinhos
[x, y, z]=global2local(long,lat,r,lon_o_poliedro,lat_o_poliedro,zz_o_poliedro); % coordinates of computation point in polyhedron-based Cartesian coordinate system

lo_poliedro=(x.^2+y.^2+z.^2).^0.5; lo=(x.^2+y.^2+z.^2).^0.5;

% Raio de integração
lo(lo>radius_metros_2)=NaN; lo(lo<=radius_metros_1)=NaN; index=~isnan(lo); 
index=double(index); index(index==0)=NaN;

lo_poliedro=lo_poliedro.*index; lo_poliedro(isnan(lo_poliedro))=[];
lon_o_poliedro=lon_o_poliedro.*index; lon_o_poliedro(isnan(lon_o_poliedro))=[];
lat_o_poliedro=lat_o_poliedro.*index; lat_o_poliedro(isnan(lat_o_poliedro))=[];

z1_gridNW=z1_gridNW.*index; z1_gridNW(isnan(z1_gridNW))=[];
z1_gridNE=z1_gridNE.*index; z1_gridNE(isnan(z1_gridNE))=[];
z1_gridSW=z1_gridSW.*index; z1_gridSW(isnan(z1_gridSW))=[];
z1_gridSE=z1_gridSE.*index; z1_gridSE(isnan(z1_gridSE))=[];

z1_refNW=z1_refNW.*index; z1_refNW(isnan(z1_refNW))=[];
z1_refNE=z1_refNE.*index; z1_refNE(isnan(z1_refNE))=[];
z1_refSW=z1_refSW.*index; z1_refSW(isnan(z1_refSW))=[];
z1_refSE=z1_refSE.*index; z1_refSE(isnan(z1_refSE))=[];

zz_o_poliedro=zz_o_poliedro.*index; zz_o_poliedro(isnan(zz_o_poliedro))=[];

x=x.*index; x(isnan(x))=[]; y=y.*index; y(isnan(y))=[]; z=z.*index; z(isnan(z))=[];

% coordinates of corners of polyhedron in polyhedron-based Cartesian coordinate system - origem no Computatition point com eixo x para norte e y para leste
dx=zz_o_poliedro*res_rad;
dy=zz_o_poliedro.*cos(lat_o_poliedro).*res_rad;
dx=0.5.*dx; dy=0.5.*dy;

xx1=-dx-x; yy1=-dy-y; zz1=z1_gridNW-zz_o_poliedro-z;
xx2=-dx-x; yy2=dy-y; zz2=z1_gridNE-zz_o_poliedro-z;
xx3=dx-x; yy3=dy-y; zz3=z1_gridSE-zz_o_poliedro-z; 
xx4=dx-x; yy4=-dy-y; zz4=z1_gridSW-zz_o_poliedro-z; 

zz_botom=1/4.*(z1_refNW+z1_refSE+z1_refNE+z1_refSW);
xx5=xx1; yy5=yy1; zz5= zz_botom-zz_o_poliedro-z;
xx6=xx2; yy6=yy2; zz6=zz_botom-zz_o_poliedro-z;
xx7=xx3; yy7=yy3; zz7=zz_botom-zz_o_poliedro-z;
xx8=xx4; yy8=yy4; zz8=zz_botom-zz_o_poliedro-z;

[m, n]=size(xx1);
size1=m*n;
 % irregular polyhedron detecting
f1=(zz1-zz5).*(zz2-zz6)<0; f2=(zz2-zz6).*(zz3-zz7)<0; f3=(zz3-zz7).*(zz4-zz8)<0; f4=(zz1-zz5).*(zz4-zz8)<0;
n1=find( f1& f2 ); n2=find( f3& f2); n3=find(f3& f4 ); n4=find( f1& f4); n5=find( f4& f2); n6=find(f1& f3);
n=[n1;n2; n3; n4; n5; n6]; n=unique(n); [p, l]=size(n); p=p*l; temp=zeros(p,1)-10000;
x_s=[xx1(n) xx2(n)  xx3(n) xx4(n)  xx5(n) xx6(n) xx7(n) xx8(n)];
y_s=[yy1(n) yy2(n)  yy3(n) yy4(n)  yy5(n) yy6(n) yy7(n) yy8(n)];
z_s=[temp temp  temp temp  zz5(n) zz6(n) zz7(n) zz8(n)];
lon_s=lon_o_poliedro(n); lat_s=lat_o_poliedro(n);
  
zz5(n)=-10000; zz6(n)=-10000; zz7(n)=-10000; zz8(n)=-10000;

x_poliedro=[xx1(:) xx2(:)  xx3(:) xx4(:)  xx5(:) xx6(:) xx7(:) xx8(:)];
y_poliedro=[yy1(:) yy2(:)  yy3(:) yy4(:)  yy5(:) yy6(:) yy7(:) yy8(:)];
z_poliedro=[zz1(:) zz2(:)  zz3(:) zz4(:)  zz5(:) zz6(:) zz7(:) zz8(:)];

x_poliedro=[x_poliedro; x_s]; y_poliedro=[y_poliedro; y_s]; z_poliedro=[z_poliedro; z_s];
  
f1=(abs(x_poliedro)<0.05); x_poliedro(f1)=x_poliedro(f1)+0.1;
f1=(abs(y_poliedro)<0.05); y_poliedro(f1)=y_poliedro(f1)+0.1;

lon_o_poliedro=[lon_o_poliedro(:);lon_s]; lat_o_poliedro=[lat_o_poliedro(:);lat_s]; lo_poliedro=[lo_poliedro(:);lo_poliedro(n)];
size1=size1+p;

% call function polyhedron_gx in fortran
g=polyhedron_g_gx(size1,x_poliedro,y_poliedro,z_poliedro);

if strcmp('sim',modelo_densidade)
densidade_o_poliedro_c=densidade_o_poliedro.*index; densidade_o_poliedro_c(isnan(densidade_o_poliedro_c))=[];
densidade_o_poliedro_c=[densidade_o_poliedro_c;densidade_o_poliedro_c(n)];
V_poliedro(1:size1,coluna)=G.*densidade_o_poliedro_c.*g(1:size1,1);
else
V_poliedro(1:size1,coluna)=G.*rho.*g(1:size1,1);
end

% clear lo lon_o_poliedro lat_o_poliedro zz_o_poliedro index x y z xx1 yy1 xx2 yy2 zz1 zz2 z1_gridNW z1_gridNE z1_gridSW z1_gridSE z1_refNW z1_refNE z1_refSW z1_refSE dx dy
% clear xx1 xx2 xx3 xx4 xx5 xx6 xx7 xx8 yy1 yy2 yy3 yy4 yy5 yy6 yy7 yy8 zz1 zz2 zz3 zz4 zz5 zz6 zz7 zz8 zz_botom m n size1 f1 f2 f3 f4 n1 n2 n3 n4 x_s y_s z_s lon_s
% clear x_poliedro y_poliedro z_poliedro g
end

Eta_RTM_poliedro=(-sum(V_poliedro(:)./2,"omitnan"))./gama+Eta_RTM_HC_Omang   % m
toc

% Anomalia de altura por fatias de raios para geração de gráfico
Eta_RTM_poliedro_raios=cumsum((-sum(V_poliedro./2,"omitnan"))./gama)+Eta_RTM_HC_Omang;   % m
radius_grafico=1:radius_metros_total/1000;

%PEGANDO SÓ VALORES DE 1 EM 1KM ATÉ 210KM
j=1;
for i=2:2:420
    Eta_RTM_poliedro_raios1km(j) = Eta_RTM_poliedro_raios(i);
    j=j+1;
end
Eta_RTM_poliedro_raios1km = Eta_RTM_poliedro_raios1km';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Geração de imagens e gráficos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1);
% plot(radius_grafico,Eta_RTM_poliedro_raios);
% xlabel('Radius of integration (km)','FontName','Times New Roman','FontSize',14)
% ylabel('RTM height anomaly (mGal)','FontName','Times New Roman','FontSize',14)
% set(gca,'FontName','Times New Roman','FontSize',14)
% Eta_RTM_poliedro_raios=Eta_RTM_poliedro_raios';