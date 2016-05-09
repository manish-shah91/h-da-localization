clc;
close all;
clear;

fc=2.442e9;
load('phase_mag_246_measurements_5x1x5_10cm.mat');
x_amp=amplitude_set;
x_phas=phase_set;
y_amp=reshape([x_amp x_amp(49) x_amp(98) x_amp(147) x_amp(196)],[],5);
y_phas=reshape([x_phas x_phas(49) x_phas(98) x_phas(147) x_phas(196)],[],5);

ma=y_amp.*(0-1);
ph=y_phas;
data=1e-3*10.^(ma/10).*exp(j*ph*pi/180);

plots=true;
N_spatial_smoothing=2;
N_NumSignals=1;
rxsig=data(46,:);
inter_element_distance=65;
lambda=122.5;
N=size(data,2);                         %number of elements in array
samples=size(data,1);
el_space=inter_element_distance/1000;   %spacing between antennas (mm to m)
lambda=lambda/1000;

hula = phased.ULA('NumElements',N,'ElementSpacing',el_space);                              %omnidirectional
%hula.Element.FrequencyRange = [2.2e8 2.6e9];
if plots
    figure(3);
    plotResponse(hula,fc,physconst('LightSpeed'),'RespCut','az','Format','polar');
    figure(4);
    plotResponse(hula,fc,physconst('LightSpeed'),'RespCut','3D','Format','uv');
    figure(6);
    plotResponse(hula,fc,physconst('LightSpeed'),'RespCut','3D','Format','polar');
    figure(5),
    viewArray(hula,'shownormals',true);
end


if plots
    %     figure(92);
    %     plot(1:samples,(rxsig));
    figure(93);
    [theta, rho]=cart2pol(real(rxsig), imag(rxsig));
    polar(theta, rho);
end

%Beamscan
hbeam2 = phased.BeamscanEstimator('SensorArray',hula,...
    'OperatingFrequency',fc,...%'ScanAngles',-90:90,...
    'DOAOutputPort',true,'NumSignals',N_NumSignals,'SpatialSmoothing',N_spatial_smoothing);
[~,ang_Beamscan] = step(hbeam2,rxsig);
if plots
    figure(2);
    plotSpectrum(hbeam2);
end

a=plotSpectrum(hbeam2);
% axesObjs = get(a, 'Children');
% dataObjs = get(axesObjs, 'Children');
% spatial_spectrum_beamscan = get(dataObjs, 'YData');
spatial_spectrum_beamscan=a.YData;

%MVDR
hbeam = phased.MVDREstimator('SensorArray',hula,...
    'OperatingFrequency',fc,'ScanAngles',-90:90,...
    'DOAOutputPort',true,'NumSignals',N_NumSignals,'SpatialSmoothing',N_spatial_smoothing);
[~,ang_MVDR] = step(hbeam,rxsig);
if plots
    figure(1);
    plotSpectrum(hbeam);
end
% a=plotSpectrum(hbeam);
% axesObjs = get(a, 'Children');
% dataObjs = get(axesObjs, 'Children');
% spatial_spectrum_MVDR = get(dataObjs, 'YData')
spatial_spectrum_MVDR=a.YData;

