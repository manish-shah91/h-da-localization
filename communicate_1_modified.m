% -----------------------------------------------------------------
% HEADER title
% -----------------------------------------------------------------
% Revision history
% communicate_anritsu_1:  -setup connection with anritsu VNA
%                         -fetch data &
% -----------------------------------------------------------------
clear all;

% --- list interfaces ---
%out=instrhwinfo
%out.SupportedInterfaces
% -----------------------

fclose('all');                               % Close all connections
visa_object = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::2907::65529::1026026_150_7::0::INSTR', 'Tag', '');
if isempty(visa_object)
    visa_object=visa('ni', 'USB0::2907::65529::1026026_150_7::0::INSTR');
else
    fclose(visa_object);
    visa_object=visa_object(1);
end

set(visa_object, 'InputBufferSize', 10000);   % Configure instrument object, obj1
set(visa_object, 'OutputBufferSize', 512);    % Configure instrument object, obj1
set(visa_object, 'Timeout', 10);               % Configure instrument object, obj1
fopen(visa_object);                           % Open connection

ID = query(visa_object, '*IDN?');
disp(strcat('Established connection with ',ID));

fprintf(visa_object, 'FREQ:START 2443 MHZ');
fprintf(visa_object, 'FREQ:STOP 2446 MHZ');
fprintf(visa_object, 'FORM:READ:DATA ASC');%INT,32  REAL,32 ASC
%query(visa_object, 'FORM:READ:DATA?');
%WERKT SOMS query(visa_object, 'SENSE:SWEEP:POINTS?')

index=1;

prompt='Measurement? [enter/n]';
x=input(prompt,'s');
while(~isempty(x))
data=query(visa_object, 'trac:data?');

    %extract data without header and comma end
    %split on comma
    %make doubles from ASCII
    %sort real & imaginary components
data=reshape(str2double(strsplit(data(8:end-1),',')),2,[]);

    %remove outliers
percentile=2;   %remove upper and lower 1%
p=prctile(data(1,:),[percentile/4 100-percentile/4]);
i=find(data(1,:)<p(1));
data(:,i)=[];
i=find(data(1,:)>p(2));
data(:,i)=[];
p=prctile(data(2,:),[percentile/4 100-percentile/4]);
i=find(data(2,:)<p(1));
data(:,i)=[];
i=find(data(2,:)>p(2));
data(:,i)=[];

%check measurement validity
standard_deviation=std(data');
standard_deviation(isnan(standard_deviation))=0;
data_mean=mean(data,2);
max_rescaled_std=max(abs(standard_deviation'./data_mean));
if (max_rescaled_std>0.1)
    disp(strcat('Possible measurement error! rescaled standard deviation ',num2str(max(standard_deviation))));
else
    disp('Measurement results seems valid');
end

phase=mod(180/pi*(atan(data_mean(2)/data_mean(1))-pi*(data_mean(1)<0))+180,360)-180
amplitude=20*log10(sqrt(data_mean(1)^2+data_mean(2)^2))


phase_set(index)=phase;
amplitude_set(index)=amplitude;


prompt='Next Measurement? [enter/n]';
x=input(prompt,'s');
index=index+1;
end


index=index-1;
prompt='save? [y/n]';
x=input(prompt,'s');
if x=='y'
    prompt='file index?';
    x=input(prompt,'s');
   save(strcat('phase_mag_',num2str(index),'_measurements_',x,'.mat'),'phase_set','amplitude_set');
end

fclose(visa_object)
disp(strcat('Connection with',ID,' closed.'));
