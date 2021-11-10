function [fk Ik] = welch_method(data,window,overlap,delta,numzeros)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectral analysis based in welch method for spacing of the spetral estimates % 
%                                                                              %
% input:        data,                                                          %
%               window,   (number of partitions, in integer number)            %
%               overlap,  (overlapping data, in percentage)                    %
%               delta,    (spatial variation, in km, m, degrees, etc)          %
%               numzeros, (number of added zeros, in percentage)               %
%                                                                              %
% output:       fk wavenumber, size: nf = round(M/2), line 27                  %
%               Ik modified periodigram, size nf = round(M/2), line 27         %
%                                                                              %
% author:       Fernando Campos (2021)                                         %
%               fcampos@cicese.edu.mx                                          %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = detrend(data);

nl = length(data);

L = round(nl/window);

D = round(overlap*L);

LL = round(numzeros*L);

M = LL+L;

U = sum(hanning(L).^2)/L;

ak = fft(detrend(data(1:L)).*hanning(L),M);

for i = 2:window
  ak = ak + fft(detrend(data(1+L*(i-1)-D:i*L)).*hanning(L+D),M);
end

ak = ak/M;

Ik = abs(ak)*M^2/(L*U);

index = 1:0.5*M;

fk = index/(M*delta); Ik = Ik(index);

return
