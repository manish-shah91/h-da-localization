% -----------------------------------------------------------------
% TEST
% -----------------------------------------------------------------
clear;
close all;

%------------------------------------------------------------------
% LOAD DATA
% comment: these data sets have 4 measurements per wavelength.
% however, 2 measurements per wavelength is also sufficient.
%------------------------------------------------------------------

% load('anechoic_array_90_-90_5200MHz_37.mat'); %load 5.2GHz measurements
% fc=5.2e9;

load('anechoic_array_90_-90_2450MHz_37.mat'); %load 2.45GHz measurements
fc=2.45e9;

%1=90º   -   10=45º   -   19=0º   -   28=-45º   -   37=-90º
ma=magn(:,:,20);
ph=phas(:,:,20);


data=1e-3*10.^(ma/10).*exp(1j*ph*pi/180);       %transform magnitude and phase to complex data



%--------------------------
% AOA estimation
%--------------------------
                
plots=true;
N_spatial_smoothing=2;  %number of spatial smoothing operations to perform
N_NumSignals=2;         %number of AOA to determine


%fc=f(frequency_selection)*1e6;          %frequency selection in MHz
rxsig=data;%(:,:,frequency_selection);    %selected measurement data
%rxsig=rxsig-mean(mean(rxsig));

% noise_power=3e-7;
% noise = noise_power*(randn(size(rxsig))+1i*randn(size(rxsig)));
% rxsig=repmat(mean(rxsig,1),size(rxsig,1),1);
% rxsig=rxsig+noise;


N=size(data,2);                         %number of elements in array
samples=size(data,1);                   %number of samples
el_space=inter_element_distance/1000;   %spacing between antennas (mm to m)
lambda=lambda/1000;                     %wavelength (mm to m)
%lambda = physconst('LightSpeed')/fc;    %wavelength

% number_of_periods=10;                                           %obsolete
% noise_power=3e-11;                                              %obsolete
% dist=[2 2.03];                                                  %obsolete
% x_loss=[(lambda/(4*pi*dist(1)))^2, (lambda/(4*pi*dist(2)))^2];  %obsolete
% ang_az=[-50 20];                                                %obsolete
% ang_el=[0 0];                                                   %obsolete


%MVDR en Beamscan

%el_patch=patchfunc();                                                                     % Define patch antenna array element
%hula = phased.ULA('Element',el_patch,'NumElements',N,'ElementSpacing',el_space);          %directional patch antennas
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
    'OperatingFrequency',fc,'ScanAngles',-90:90,...
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

                

%---------------------                
                



