# 1D Spectral analysis based on Welch method in Matlab
Matlab functions to compute time and spatial spectra in 1D based on Welch method 

# Getting started example data:


```MATLAB
clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code to compute 15 first modes of EOF from monthly data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi = ncread('input.nc','data');

lon = ncread('input.nc','longitude');

lat = ncread('input.nc','latitude');

[yy xx] = meshgrid(lat,lon);

time = double(datenum('2000-01-15')+ncread('input.nc','time')*30);

figure;
subplot(5,1,[1 2 3 4]); 
pcolor(xx,yy,squeeze(xi(:,:,1))); shading interp; colormap jet; caxis([-1 1]*3); 
colorbar; title('month 1'); hold on; plot(275,-15,'x');
subplot(5,1,5); 
plot(time,squeeze(xi(37,28,:))); datetick('x','yyyy'); grid on
title('time serie at X');
```
![alt text](https://github.com/fecampos/spectral_analysis_1D/blob/main/example.png)

# Getting started welch_method.m for temporal data:

```MATLAB
[fk Ik] = welch_method(squeeze(xi(20,20,:)),100,0.1,1,0,'t');
```
![alt text](https://github.com/fecampos/spectral_analysis_1D/blob/main/example_time_spectrum.png)

# Getting started welch_method.m for spatial data:
