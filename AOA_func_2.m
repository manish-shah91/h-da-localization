
% -------------------------------------------------------------------------
% RELEASE NOTES
% -------------------------------------------------------------------------
% AOA_func_1:  -New function, based on the AOA_2.m file
%              -Introduced el_per_lambda = # elements per wavelength
%              -Tested with AOA_3.m
%              -Try/Catch in Esprit/rootMUSIC
% -------------------------------------------------------------------------


function [ang_MVDR,spatial_spectrum_MVDR, ang_Beamscan, spatial_spectrum_beamscan, ang_Music, ang_Esprit]=AOA_func_2(N,N_spatial_smoothing,N_NumSignals, fc, el_per_lambda, samples, number_of_periods, noise_power, dist, x_loss, ang_az, ang_el,plots)

%fs = 8000; t = (0:1/fs:1).';
%x1 = cos(2*pi*t*300); x2 = cos(2*pi*t*400);
%fc = 1e9;
%lambda = physconst('LightSpeed')/fc;
%x1 = 1*ones(300000,1);
%x2 = 1/3*ones(300000,1);
%ha = phased.ULA('Element',el_patch,'NumElements',10,'ElementSpacing',lambda/2);
%ha.Element.FrequencyRange = [100e6 300e6];

%z = collectPlaneWave(ha,[x1 x2],[10 0;60 0]',fc);   %houdt geen rekening met elementen: enkel voor isotrope elementen
%noise = 0.001*(randn(size(z))+1i*randn(size(z)));
%hdoa = phased.BeamscanEstimator('SensorArray',ha,...
%    'OperatingFrequency',fc,...
%    'DOAOutputPort',true,'NumSignals',2);
%[y,doas] = step(hdoa,z+noise);
%doas = broadside2az(sort(doas),[20 -5]);
%figure(1);
%plotSpectrum(hdoa);

% hsv=phased.SteeringVector('IncludeElementResponse',true,'SensorArray',ha);
% sv=step(hsv,fc,[10;45])
% y1=x1*sv.';
% sv=step(hsv,fc,[60;45]);
% y2=x2*sv.';
% noise = 0.001*(randn(size(y1))+1i*randn(size(y1)));
% hdoa = phased.BeamscanEstimator('SensorArray',ha,...
%     'OperatingFrequency',fc,...
%     'DOAOutputPort',true,'NumSignals',2);
% [y,doas] = step(hdoa,y1+y2+noise)
% doas = broadside2az(sort(doas),[45 45]);
% figure(2);
% plotSpectrum(hdoa);


%MVDR en Beamscan
%N = 10;                                 % # antenna elements        %becomes function input
%N_spatial_smoothing=1;                                              %becomes function input
%N_NumSignals=3;                                                     %becomes function input
%fc = 2.4e9;                             % Carrier frequency

lambda = physconst('LightSpeed')/fc;
el_space=lambda/el_per_lambda;                      % Element spacing

%hula = phased.ULA('NumElements',N,'ElementSpacing',el_space);
el_patch=patchfunc();                   % Define patch antenna array element
hula = phased.ULA('Element',el_patch,'NumElements',N,'ElementSpacing',el_space);
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


%signals

%samples=3e3;                                            %becomes function input
%number_of_periods=10;                                   %becomes function input
%noise_power=3e-11;  % -75 dbm = 3.16e-11 Watt           %becomes function input
% dist_1=3;
% dist_2=3.33;%4.24;
%dist=[3 13.33 8];                                          %becomes function input
%x_loss=[(lambda/(4*pi*dist(1)))^2, (lambda/(4*pi*dist(2)))^2 (lambda/(4*pi*dist(3)))^2];      %becomes function input
%ang_az= [10 5 -40];                                        %becomes function input
%ang_el= [0 0 0];                                          %becomes function input

%material_interface;





distance_difference=min(dist)-dist; %[m]
phase_difference=distance_difference/lambda*2*pi;
x = sqrt(2)*ones(samples,length(phase_difference)).*(cos(repmat((2*pi*(1/samples:1/samples:1)*number_of_periods)',1,length(phase_difference))+repmat(phase_difference,samples,1)));

%x2 = sqrt(2)*ones(samples,1).*transpose(cos(2*pi*(1/samples:1/samples:1)*number_of_periods+phase_difference(2)));




% R_s=0.28%0.5;   32,96
% R_p=R_s;
% R_p=0.2;
x=x.*repmat(x_loss,size(x,1),1);

% if plots
%     figure(91);
%     plot(1:samples,x(:,1),1:samples,x(:,2));
% end

% ang1 = [10;00];
% ang2 = [10;00];

ang=[ang_az;ang_el];
hsv=phased.SteeringVector('IncludeElementResponse',true,'SensorArray',hula); %steering vector voor array definieren
sv=step(hsv,fc,ang);   %release(hsv) before step(hsv,...) with different arguments (more angles)
% y1=x(:,1)*sv(:,1).';
% %sv=step(hsv,fc,ang2);
% y2=x(:,2)*sv(:,2).';
% rng default;
% noise = noise_power*(randn(size(y1))+1i*randn(size(y1)));
% rxsig = y1 + y2 + noise;

% x=reshape(repmat(x,N,1),size(x,1),[]);
% sv=repmat(reshape(sv,1,[]),size(x,1),1);
% y=x.*sv;
% rx_sig=y*repmat(diag(ones(N,1)),size(y,2),1)

x=reshape(repmat(x,N,1),size(x,1),[]);
sv=repmat(reshape(sv,1,[]),size(x,1),1);
y=x.*sv;
rxsig=y*repmat(diag(ones(N,1)),size(y,2)/N,1);
noise = noise_power*(randn(size(rxsig))+1i*randn(size(rxsig)));
rxsig=rxsig+noise;


% figure(91);
% plot(1:samples,x(:,1),1:samples,x(:,2));
if plots
    figure(92);
    plot(1:samples,(rxsig));
end




% x1 = 1*ones(samples,1);
% x2 = 1/3*ones(samples,1);
% ang1 = [10;10];
% ang2 = [60;10];
% hsv=phased.SteeringVector('IncludeElementResponse',true,'SensorArray',hula); %steering vector voor array definieren
% sv=step(hsv,fc,ang1);
% y1=x1*sv.';
% sv=step(hsv,fc,ang2);
% y2=x2*sv.';
% rng default;
% noise = noise_power*(randn(size(y1))+1i*randn(size(y1)));
% rxsig = y1 + y2 + noise;


% hwav = phased.LinearFMWaveform('SampleRate',1e7,'SweepBandwidth',1e5,...
%     'PulseWidth',5e-6,'OutputFormat','Pulses','NumPulses',1);
% sig1 = step(hwav); % returns samples of the linear FM pulse in a column vector Y
% sig2 = sig1;
% ang1 = [30;0];
% ang2 = [-10;0];
% arraysig = collectPlaneWave(hula,[sig1 sig2],[ang1 ang2],fc);   % returns the received signals at the sensor array, H, when the input signals indicated by X arrive at the array from the directions specified in ANG.
% rng default;              %rng('default') puts the settings of the random number generator used by rand, randi, and randn to their default values so that they produce the same random numbers as if you restarted MATLAB.
% npower = 0.01;
% noise = sqrt(npower/2)*...
%     (randn(size(arraysig))+1i*randn(size(arraysig)));
% rxsig = arraysig+noise;



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
%if plots
%figure(1);
%plotSpectrum(hbeam);
%end
a=plotSpectrum(hbeam);
% axesObjs = get(a, 'Children');
% dataObjs = get(axesObjs, 'Children');
% spatial_spectrum_MVDR = get(dataObjs, 'YData')
spatial_spectrum_MVDR=a.YData;



%High resolution DoA

%ha = phased.ULA('NumElements',N,'ElementSpacing',0.5)

% Model the multichannel received signals at the array
%fc = 300e6;                               % Operating frequency
%fs = 8192;                                % Sampling frequency
%lambda = physconst('LightSpeed')/fc;      % Wavelength
%pos = getElementPosition(hula)/lambda;      % Element position in wavelengths

%ang1 = [90;73];
%ang2 = [90;68];           % azimuth, elevation: Direction of the signals
%angs = ang;%[ang1 ang2];
%Nsamp = 1024;                             % Number of snapshots


%rs = rng(2012);                           % Set random number generator
% ----!!!!!---- x = sensorsig(pos,Nsamp,angs,npower);       %Simulate received signal at sensor array



x=rxsig;


%----------------MUSIC------------------
%hdoaMusic.ForwardBackwardAveraging = true;

AOA_error_Music=1;
N_NumSignals_Music=N_NumSignals+1;

while AOA_error_Music
    try
        if N_NumSignals_Music>1
            N_NumSignals_Music=N_NumSignals_Music-1;
            hdoaMusic = phased.RootMUSICEstimator('SensorArray',hula,...
                'OperatingFrequency',fc,...
                'NumSignalsSource','Property','NumSignals',N_NumSignals_Music,'SpatialSmoothing',N_spatial_smoothing);
            ang_Music = step(hdoaMusic,x);
        else
            disp('an error occurred during MUSIC AOA estimation')
            ang_Music=[];
        end
        
        AOA_error_Music=0;
    catch
        AOA_error_Music=1;
    end
end

%figure(3);
%plotSpectrum(hdoaMusic);





%----------------ESPRIT------------------

AOA_error_esprit=1;
N_NumSignals_Esprit=N_NumSignals+1;
while AOA_error_esprit
    try
        if N_NumSignals_Esprit>1
        N_NumSignals_Esprit=N_NumSignals_Esprit-1;
        hdoaEsprit = phased.ESPRITEstimator('SensorArray',hula,...
            'OperatingFrequency',fc,'ForwardBackwardAveraging',true,...
            'NumSignalsSource','Property','NumSignals',N_NumSignals_Esprit,'SpatialSmoothing',N_spatial_smoothing);
        ang_Esprit = step(hdoaEsprit,x);
        else
            ang_Esprit = [];
        end
        AOA_error_esprit=0;
    catch
        disp('an error occurred during ESPRIT AOA estimation')
        AOA_error_esprit=1;
    end
end
% release(hdoaEsprit);
% hdoaEsprit.NumSignalsSource = 'Auto';
% hdoaEsprit.NumSignalsMethod = 'AIC';
% ang_Esprit = step(hdoaEsprit,x);

% hdoaBSEsprit = phased.BeamspaceESPRITEstimator('SensorArray',hula,...
%                'OperatingFrequency',fc,...
%                'NumBeamsSource','Property','NumBeams',3,...
%                'BeamFanCenter',20);
% ang_hdoaBSEsprit = step(hdoaBSEsprit,x);


%-----------------voorlopig_einde
%
%
% hdoaWSF = phased.RootWSFEstimator('SensorArray',hula,...
%           'OperatingFrequency',fc,'MaximumIterationCount',2);
% ang_hdoaWSF = step(hdoaWSF,x)
%
% release(hdoaEsprit);
% hdoaEsprit.RowWeighting = 4
% ang_hdoaEsprit = step(hdoaEsprit,x)
%
% scov = eye(4);
% magratio = [1;0.25;0.5];
% scov(1:3,1:3) = magratio*magratio';
%
% % Incident azimuth
% az_ang = [-23 0 12 40];
% % When the elevation is zero, the azimuth within [-90 90] is the same as
% % the broadside angle.
% el_ang = zeros(1,4);
%
% % The received signals
% %x = sensorsig(pos,Nsamp,[az_ang; el_ang],npower,scov);
% %rng(rs);                                 % Restore random number generator
%
% release(hdoaEsprit);
% hdoaEsprit.NumSignalsSource = 'Auto';
% hdoaEsprit.NumSignalsMethod = 'AIC';
% ang_hdoaEsprit = step(hdoaEsprit,x)
%
% release(hdoaWSF);
% hdoaWSF.NumSignalsSource = 'Property';
% hdoaWSF.NumSignals = 4;
% ang_hdoaWSF = step(hdoaWSF,x)
%
% release(hdoaMusic);
% hdoaMusic.NumSignalsSource = 'Property';
% hdoaMusic.NumSignals = 4;
% hdoaMusic.ForwardBackwardAveraging = true;
% ang_hdoaMusic = step(hdoaMusic,x)
%
% release(hdoaMusic);
% Nr = 2;    % Number of multipath reflections
% hdoaMusic.SpatialSmoothing = Nr
% ang_hdoaMusic_spatialsmoothing = step(hdoaMusic,x)