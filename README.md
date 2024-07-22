# RTM_Codes
Codes for applying the RTM technique using different approaches

These routines will perform the RTM technique through different approaches, with each routine corresponding to a single approach.

In the folder where the routine will be located, place the "Funcoes" folder and a "Dados" folder where the rasters containing density models (if used), the detailed DEM, and the reference DEM will be placed.

In the "Importação de imagens" section, you can choose whether to use the density model or not, simply by writing "sim" or "nao" where requested.

In the "Resolução" section, you must specify the resolution of the detailed raster (res) and the integration radius (radius_metros).

In "minlat" and "minlon," provide the LATITUDE N and LONGITUDE W, respectively.

In "estacao_RTM," load a txt file containing the latitude, longitude, and ellipsoidal altitude of the calculation station.

By taking these precautions, the routine will likely calculate and return the results adequately. If you have any questions, contact "thiago.k.padilha@gmail.com."
