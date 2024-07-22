tic; clc; clear all; close all; format long g

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
res=0.0002777777777777779402;              % Resolução em °
res_rad=res*rad;                       % Resolução em rad
res_metros=res*3600*30;                % Resolução em metros
radius_metros_poliedro=2160; radius_metros_prisma=3240; radius_metros_tesseroide=16200;
radius_poliedro=radius_metros_poliedro/108000;           % raio de integração em ° para poliedro (considerando 1"=30 m e 1°=108000 m) - 2,34 km
radius_metros=210000;                   % raio de integração em metros
radius_pixel_poliedro=ceil(radius_poliedro/res); % raio de integração em pixels arredondado para cima abordagem poliedro

% Coordenadas geodésicas do canto superior esquerdo do primeiro pixel da imagem no canto NW - superior esquerdo - depois passando para o centro geométrico do pixel
minlat=(-13.2468055555555431-res)*rad-(res_rad/2);    % latN
minlon=(-50.5781944444444491+res)*rad+(res_rad/2); %Lon W

% Grade de coordenadas geodésicas
gradeMDA_lat=(minlat:-res_rad:(minlat-(lin-1)*res_rad))'.*ones(1,col);
gradeMDA_long=(minlon:res_rad:(minlon+(col-1)*res_rad)).*ones(lin,1);
% gradeMDA_lat=flipud(gradeMDA_lat); H_detailed=flipud(H_detailed); H_reference=flipud(H_reference);
lat_o=gradeMDA_lat(:); long_o=gradeMDA_long(:); H_detailed_o=H_detailed(:); H_reference_o=H_reference(:);
Nell_o=a./(1 - (e2.*sin(lat_o).*sin(lat_o))).^0.5; zzp=H_detailed_o+Nell_o;

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

lat_o_NW=sub_grade_lat_NW(:); lat_o_NE=sub_grade_lat_NE(:); lat_o_SW=sub_grade_lat_SW(:); lat_o_SE=sub_grade_lat_SE(:);
long_o_NW=sub_grade_long_NW(:); long_o_NE=sub_grade_long_NE(:); long_o_SW=sub_grade_long_SW(:); long_o_SE=sub_grade_long_SE(:);
H_detailed_o_NW=sub_grade_H_detailed_NW(:); H_detailed_o_NE=sub_grade_H_detailed_NE(:); H_detailed_o_SW=sub_grade_H_detailed_SW(:); H_detailed_o_SE=sub_grade_H_detailed_SE(:);
H_reference_o_NW=sub_grade_H_reference_NW(:); H_reference_o_NE=sub_grade_H_reference_NE(:); H_reference_o_SW=sub_grade_H_reference_SW(:); H_reference_o_SE=sub_grade_H_reference_SE(:);

sin_lat_o_poliedro_NW=sin(lat_o_NW); sin_lat_o_poliedro_NE=sin(lat_o_NE); sin_lat_o_poliedro_SW=sin(lat_o_SW); sin_lat_o_poliedro_SE=sin(lat_o_SE);

Nell_NW=a./(1 - (e2.*sin_lat_o_poliedro_NW.*sin_lat_o_poliedro_NW)).^0.5; Nell_SE=a./(1 - (e2.*sin_lat_o_poliedro_SE.*sin_lat_o_poliedro_SE)).^0.5; Nell_NE=a./(1 - (e2.*sin_lat_o_poliedro_NE.*sin_lat_o_poliedro_NE)).^0.5; Nell_SW=a./(1 - (e2.*sin_lat_o_poliedro_SW.*sin_lat_o_poliedro_SW)).^0.5;

% define origin of polyhedron-based coordinate system
z1_gridNW=H_detailed_o_NW+Nell_NW; z1_gridSE=H_detailed_o_SE+Nell_SE; z1_gridNE=H_detailed_o_NE+Nell_NE; z1_gridSW=H_detailed_o_SW+Nell_SW;
z1_refNW=H_reference_o_NW+Nell_NW; z1_refSE=H_reference_o_SE+Nell_SE; z1_refNE=H_reference_o_NE+Nell_NE; z1_refSW=H_reference_o_SW+Nell_SW;

clear sub_grade_lat_NW sub_grade_lat_NE sub_grade_lat_SW sub_grade_lat_SE sub_grade_long_NW sub_grade_long_NE sub_grade_long_SW sub_grade_long_SE
clear sub_grade_H_detailed_NW sub_grade_H_detailed_NE sub_grade_H_detailed_SW sub_grade_H_detailed_SE
clear sub_grade_H_reference_NW sub_grade_H_reference_NE sub_grade_H_reference_SW sub_grade_H_reference_SE
clear sin_lat_o_poliedro_NW sin_lat_o_poliedro_NE sin_lat_o_poliedro_SW sin_lat_o_poliedro_SE
clear Nell_NW Nell_SE Nell_NE Nell_SW

lon_o_poliedro=1/2.*(long_o_NW+long_o_NE); % cetro geométrico do poliedro considerando os 4 pixels vizinhos
lat_o_poliedro=1/2.*(lat_o_NW+lat_o_SW); % cetro geométrico do poliedro considerando os 4 pixels vizinhos
zz_o_poliedro=1/4.*(z1_gridNW+z1_gridSE+z1_gridNE+z1_gridSW); % cetro geométrico do poliedro considerando os 4 pixels vizinhos
[x, y, z]=global2local(long,lat,r,lon_o_poliedro,lat_o_poliedro,zz_o_poliedro); % coordinates of computation point in polyhedron-based Cartesian coordinate system

lo_poliedro=(x.^2+y.^2+z.^2).^0.5; lo=(x.^2+y.^2+z.^2).^0.5;

% Raio de integração
lo(lo>radius_metros_poliedro)=NaN; index=~isnan(lo); index=double(index); index(index==0)=NaN;

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

if strcmp('sim',modelo_densidade)
sub_grade_densidade_NW=(densidade(lin_CP-radius_pixel_poliedro:lin_CP+radius_pixel_poliedro-1 , col_CP-radius_pixel_poliedro:col_CP+radius_pixel_poliedro-1));
sub_grade_densidade_NE=(densidade(lin_CP-radius_pixel_poliedro:lin_CP+radius_pixel_poliedro-1 , col_CP-radius_pixel_poliedro+1:col_CP+radius_pixel_poliedro));
sub_grade_densidade_SW=(densidade(lin_CP-radius_pixel_poliedro+1:lin_CP+radius_pixel_poliedro , col_CP-radius_pixel_poliedro:col_CP+radius_pixel_poliedro-1));
sub_grade_densidade_SE=(densidade(lin_CP-radius_pixel_poliedro+1:lin_CP+radius_pixel_poliedro , col_CP-radius_pixel_poliedro+1:col_CP+radius_pixel_poliedro));
densidade_o_poliedro=1/4.*(sub_grade_densidade_NW(:)+sub_grade_densidade_NE(:)+sub_grade_densidade_SW(:)+sub_grade_densidade_SE(:));
densidade_o_poliedro=densidade_o_poliedro.*index; densidade_o_poliedro(isnan(densidade_o_poliedro))=[];
clear sub_grade_densidade_NW sub_grade_densidade_NE sub_grade_densidade_SW sub_grade_densidade_SE
end

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

lon_o_poliedro=[lon_o_poliedro(:);lon_s];lat_o_poliedro=[lat_o_poliedro(:);lat_s]; lo_poliedro=[lo_poliedro(:);lo_poliedro(n)];
size1=size1+p;

% call function polyhedron_gx in fortran
g=polyhedron_g_gx(size1,x_poliedro,y_poliedro,z_poliedro);
g(1:size1,4)=-g(1:size1,4);
[temp_Vx, temp_Vy, temp_Vz]=g_prism2point(g(1:size1,2),g(1:size1,3),g(1:size1,4),long,lat,lon_o_poliedro(:),lat_o_poliedro(:));  % m/s2

if strcmp('sim',modelo_densidade)
densidade_o_poliedro=[densidade_o_poliedro;densidade_o_poliedro(n)];
Vz_poliedro=G.*densidade_o_poliedro.*temp_Vz;
else
Vz_poliedro=G.*rho.*temp_Vz;
end

dg_RTM_poliedro=sum(Vz_poliedro,"omitnan").*1D5   % mGal

clear lo index x y z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%% Cálculo dos distúrbios de gravidade RTM no entorno da estação de cálculo %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZONA DE CONEXÃO ENTRE POLIEDRO E PRISMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x, y, z]=global2local(long,lat,r,long_o,lat_o,zzp); % coordinates of computation point in prism-based Cartesian coordinate system
lo_pp=(x.^2+y.^2+z.^2).^0.5; lo=(x.^2+y.^2+z.^2).^0.5;

% Raio de integração
lo(lo>radius_metros_poliedro+res_metros/2)=NaN;
lo(lo<=radius_metros_poliedro-res_metros/2)=NaN;
index=~isnan(lo); index=double(index);
index(index==0)=NaN;

lo_pp=lo_pp.*index; lo_pp(isnan(lo_pp))=[];
long_o_pp=long_o.*index; long_o_pp(isnan(long_o_pp))=[];
lat_o_pp=lat_o.*index; lat_o_pp(isnan(lat_o_pp))=[];
H_detailed_o_pp=H_detailed_o.*index; H_detailed_o_pp(isnan(H_detailed_o_pp))=[];
H_reference_o_pp=H_reference_o.*index; H_reference_o_pp(isnan(H_reference_o_pp))=[];
Nell_o_pp=Nell_o.*index; Nell_o_pp(isnan(Nell_o_pp))=[];

xx1=long_o_pp-res_rad/2;                       % canto inferior esquerdo do pixel
yy1=lat_o_pp-res_rad/2;          % canto inferior esquerdo do pixel
xx2=long_o_pp+res_rad/2;          % canto superior direito do pixel
yy2=lat_o_pp+res_rad/2;                       % canto superior direito do pixel
zz1=H_detailed_o_pp+Nell_o_pp;
zz2=H_reference_o_pp+Nell_o_pp;

if strcmp('sim',modelo_densidade)
densidade_o_pp=densidade(:); densidade_o_pp=densidade_o_pp.*index; densidade_o_pp(isnan(densidade_o_pp))=[];
Vz_pp=Prism_g_nosum(xx1,xx2,yy1,yy2,zz1,zz2,long,lat,r,densidade_o_pp);  % m/s2
else
Vz_pp=Prism_g_nosum(xx1,xx2,yy1,yy2,zz1,zz2,long,lat,r,rho);  % m/s2
end
dg_RTM_pp=(-sum(Vz_pp,"omitnan").*1D5)./2   % mGal

clear lo index x y z xx1 yy1 xx2 yy2 zz1 zz2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%% Cálculo dos distúrbios de gravidade RTM no entorno da estação de cálculo %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRISMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x, y, z]=global2local(long,lat,r,long_o,lat_o,zzp); % coordinates of computation point in prism-based Cartesian coordinate system
lo_prism=(x.^2+y.^2+z.^2).^0.5; lo=(x.^2+y.^2+z.^2).^0.5;

% Raio de integração
lo(lo>radius_metros_prisma)=NaN;
lo(lo<=radius_metros_poliedro)=NaN;
index=~isnan(lo); index=double(index);
index(index==0)=NaN;

lo_prism=lo_prism.*index; lo_prism(isnan(lo_prism))=[];
long_o_prism=long_o.*index; long_o_prism(isnan(long_o_prism))=[];
lat_o_prism=lat_o.*index; lat_o_prism(isnan(lat_o_prism))=[];
H_detailed_o_prism=H_detailed_o.*index; H_detailed_o_prism(isnan(H_detailed_o_prism))=[];
H_reference_o_prism=H_reference_o.*index; H_reference_o_prism(isnan(H_reference_o_prism))=[];
Nell_o_prism=Nell_o.*index; Nell_o_prism(isnan(Nell_o_prism))=[];

xx1=long_o_prism-res_rad/2;                       % canto inferior esquerdo do pixel
yy1=lat_o_prism-res_rad/2;          % canto inferior esquerdo do pixel
xx2=long_o_prism+res_rad/2;          % canto superior direito do pixel
yy2=lat_o_prism+res_rad/2;                       % canto superior direito do pixel
zz1=H_detailed_o_prism+Nell_o_prism;
zz2=H_reference_o_prism+Nell_o_prism;

if strcmp('sim',modelo_densidade)
densidade_o_prism=densidade(:); densidade_o_prism=densidade_o_prism.*index; densidade_o_prism(isnan(densidade_o_prism))=[];
Vz_prisma=Prism_g_nosum(xx1,xx2,yy1,yy2,zz1,zz2,long,lat,r,densidade_o_prism);  % m/s2
else
Vz_prisma=Prism_g_nosum(xx1,xx2,yy1,yy2,zz1,zz2,long,lat,r,rho);  % m/s2
end
dg_RTM_prisma=-sum(Vz_prisma,"omitnan").*1D5   % mGal

clear lo index x y z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%% Cálculo dos distúrbios de gravidade RTM no entorno da estação de cálculo %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESSEROIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x, y, z]=global2local(long,lat,r,long_o,lat_o,zzp); % coordinates of computation point in tesseroid-based Cartesian coordinate system
lo=(x.^2+y.^2+z.^2).^0.5;

% Raio de integração
lo(lo>radius_metros_tesseroide)=NaN;
lo(lo<=radius_metros_prisma)=NaN;
index=~isnan(lo); index=double(index);
index(index==0)=NaN;

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

delta_lat=res_rad; delta_long=res_rad; cte=G*delta_lat*delta_long;

cos_psi_o=sin_lat_o.*sin_lat + cos_lat_o.*cos_lat.*cos_dlong;
sin_psi_o=real((1-cos_psi_o.*cos_psi_o).^0.5);

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
dg_RTM_tesseroide=sum(Vz_tesseroide,"omitnan").*1D5   % mGal

clear lo index x y z sin_lat_o cos_lat_o sin_long_o cos_long_o sin_dlong cos_dlong r1 r2 delta_r ro cos_psi_o sin_psi_o L000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%% Cálculo dos distúrbios de gravidade RTM no entorno da estação de cálculo %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PONTO MASSA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x, y, z]=global2local(long,lat,r,long_o,lat_o,zzp); % coordinates of computation point in tesseroid-based Cartesian coordinate system
lo=(x.^2+y.^2+z.^2).^0.5;

% Raio de integração
lo(lo<=radius_metros_tesseroide)=NaN;
index=~isnan(lo); index=double(index);
index(index==0)=NaN;

long_o_ponto_massa=long_o.*index; long_o_ponto_massa(isnan(long_o_ponto_massa))=[];
lat_o_ponto_massa=lat_o.*index; lat_o_ponto_massa(isnan(lat_o_ponto_massa))=[];
H_detailed_o_ponto_massa=H_detailed_o.*index; H_detailed_o_ponto_massa(isnan(H_detailed_o_ponto_massa))=[];
H_reference_o_ponto_massa=H_reference_o.*index; H_reference_o_ponto_massa(isnan(H_reference_o_ponto_massa))=[];
Nell_o_ponto_massa=Nell_o.*index; Nell_o_ponto_massa(isnan(Nell_o_ponto_massa))=[];

sin_lat_o=sin(lat_o_ponto_massa); cos_lat_o=cos(lat_o_ponto_massa); sin_long_o=sin(long_o_ponto_massa); cos_long_o=cos(long_o_ponto_massa);
sin_dlong=sin_long_o.*cos_long-cos_long_o.*sin_long; cos_dlong=cos_long_o.*cos_long+sin_long_o.*sin_long;
r1=H_reference_o_ponto_massa+Nell_o_ponto_massa;
r2=H_detailed_o_ponto_massa+Nell_o_ponto_massa;
delta_r=r2-r1;
ro=(r1+r2)./2; ro_2=ro.*ro;

clear cos_long_o densidade_o_poliedro densidade_o_pp densidade_o_prism densidade_o_tesseroide
clear gradeMDA_lat gradeMDA_long H_detailed H_detailed_o H_detailed_o_NE H_detailed_o_NW H_detailed_o_SE H_detailed_o_SW
clear H_detailed_o_tesseroide H_detailed_o_ponto_massa H_detailed_o_prism H_detailed_o_pp
clear H_reference_o_tesseroide H_reference_o_ponto_massa H_reference_o_prism H_reference_o_pp
clear H_reference H_reference_o H_reference_o_NE H_reference_o_NW H_reference_o_SE H_reference_o_SW dx dy f1 f2 f3 f4 g
clear L002 L020 L200 lat_o lat_o_NE lat_o_NW lat_o_poliedro lat_o_ponto_massa lat_o_pp lat_o_prism lat_o_SE lat_o_SW lat_o_tesseroide
clear lo lo_tesseroide_2 lo_tesseroide_3 lo_tesseroide_4 lo_tesseroide_5 lo_tesseroide_6
clear lo_tesseroide_7 lon_o_poliedro lon_s long_o long_o_NE long_o_SW long_o_ponto_massa long_o_pp long_o_prism long_o_SE long_o_SW
clear long_o_tesseroide n n1 n2 n3 n4 n5 n6 Nell_o Nell_o_ponto_massa Nell_o_pp Nell_o_prism Nell_o_tesseroide r1 r2 temp temp_Vx
clear temp_Vy temp_Vz x x_poliedro x_s xx1 xx2 xx3 xx4 xx5 xx6 xx7 xx8 y y_poliedro y_s yy1 yy2 yy3 yy4 yy5 yy6 yy7 yy8 z z1_gridNE
clear z1_gridNW z1_gridSE z1_gridSW z1_refNE z1_refNW z1_refSE z1_refSW z_poliedro z_2 zz1 zz2 zz3 zz4 zz5 zz6 zz7 zz8 zz_botom
clear zz_o_poliedro zzp lat_s long_o_NW z_s

cos_psi_o=sin_lat_o.*sin_lat + cos_lat_o.*cos_lat.*cos_dlong;
sin_psi_o=real((1-cos_psi_o.*cos_psi_o).^0.5);

lo_ponto_massa=real((r_2+ro_2-2.*r.*ro.*cos_psi_o).^0.5);
lo_ponto_massa_3=lo_ponto_massa.*lo_ponto_massa.*lo_ponto_massa;

L000=(ro_2.*(r-ro.*cos_psi_o).*cos_lat_o)./(lo_ponto_massa_3);

if strcmp('sim',modelo_densidade)
densidade_o_ponto_massa=densidade(:); densidade_o_ponto_massa=densidade_o_ponto_massa.*index; densidade_o_ponto_massa(isnan(densidade_o_ponto_massa))=[];
Vz_ponto_massa=densidade_o_ponto_massa.*cte.*delta_r.*L000;  % m/s2
else
Vz_ponto_massa=rho.*cte.*delta_r.*L000;  % m/s2
end

dg_RTM_ponto_massa=sum(Vz_ponto_massa,"omitnan").*1D5   % mGal

clear index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Total %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dg_RTM_total=  dg_RTM_poliedro  +  dg_RTM_pp  +  dg_RTM_prisma  +  dg_RTM_tesseroide  +  dg_RTM_ponto_massa  +  dg_RTM_HC   % mGal

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dg por fatias de raios para geração de gráfico
lo_todos=[lo_poliedro; lo_pp; lo_prism; lo_tesseroide; lo_ponto_massa];
Vz_todos=[Vz_poliedro; -Vz_pp./2; -Vz_prisma; Vz_tesseroide; Vz_ponto_massa];

dg_RTM_todos_raios=zeros(1,radius_metros/1000);
radius_grafico=zeros(1,radius_metros/1000);
for i=1:radius_metros/1000
lo_todos_grafico=lo_todos;
lo_todos_grafico(lo_todos_grafico>(i*1000))=NaN;
index=~isnan(lo_todos_grafico);
Vz_grafico=Vz_todos.*index;    %m/s2
dg_RTM_todos_raios(1,i)=sum(Vz_grafico,"omitnan").*1D5+dg_RTM_HC;   % mGal
radius_grafico(1,i)=i*1000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Geração de imagens e gráficos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
plot(radius_grafico,dg_RTM_todos_raios);
xlabel('Radius of integration (km)','FontName','Times New Roman','FontSize',14)
ylabel('RTM gravity disturbance (mGal)','FontName','Times New Roman','FontSize',14)
set(gca,'FontName','Times New Roman','FontSize',14)
dg_RTM_todos_raios=dg_RTM_todos_raios';