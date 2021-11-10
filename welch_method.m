function [fk Ik] = welch_method(data,window,overlap,delta,numzeros,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectral analysis based in welch method for spetral estimates                % 
%                                                                              %
% input:        data,                                                          %
%               window,   (number of points of the windows as an integer)      %
%               overlap,  (overlapping data, [0-1])                            %
%               delta,    (spatial variation, in km, m, degrees, etc)          %
%               numzeros, (number of added zeros, [0-1])                       %
%               type,     (time->'t' or spatial->'s')                          %
%                                                                              %
% output:       fk wavenumber, size: nf = round(M/2)                           %
%               Ik modified periodigram, size nf = round(M/2)                  %
%                                                                              %
% author:       Fernando Campos (2021)                                         %
%               fcampos@cicese.edu.mx                                          %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if type == 's'

  data = detrend(data);
  nl   = length(data);
  L    = window;
  ns   = floor(nl/L);
  D    = round(overlap*L);
  LL   = round(numzeros*L);
  M    = LL+L;
  U    = sum(hanning(L).^2)/L;
  ak   = fft(detrend(data(1:L)).*hanning(L),M);
  for i = 2:ns
    ak = ak + fft(detrend(data(1+L*(i-1)-D:i*L)).*hanning(L+D),M);
  end
  ak   = ak/M;
  Ik   = abs(ak)*M^2/(L*U);
  index = 1:0.5*M;
  fk   = index/(M*delta); Ik = Ik(index);

elseif  type == 't'
  numzeros = 0;
  data = detrend(data);
  nl   = length(data);
  L    = window;
  ns   = floor(nl/L);
  D    = round(overlap*L);
  LL   = round(numzeros*L);
  M    = LL+L;
  U    = sum(hanning(L).^2)/L; 
  ak   = fft(detrend(data(1:L)).*hanning(L),M);
  for i = 2:ns
    ak = ak + fft(detrend(data(1+L*(i-1)-D:i*L)).*hanning(L+D),M);
  end
  ak   = ak/M;
  Ik   = abs(ak)*M^2/(L*U);
  index = 1:0.5*M;
  fk   = index/(M*delta); Ik = Ik(index);

end

return
