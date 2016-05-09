% -----------------------------------------------------------------
% HEADER title
% -----------------------------------------------------------------
% Revision history
% communicate_anritsu_1:  -setup connection with anritsu VNA
%                         -fetch data
% communicate_mac_2:      -setup connection with prologix & 8753ES
%                         -fetch complex S21 parameters
%                         -detect scaled standard deviation errors
%                         -parameters: f, points, points_multiplier
% -----------------------------------------------------------------


% Measurement reverb 1 & 2:
% large metal box darmstadt: x= 3.5m, y= 3m
%                            rx position: x=1.15, y=0.01
%                            tx position 1: x=1.75, y=1.5
%                            tx position 2: x=0.5, y=2.35

% Measurement 3 & 4:
% large metal box darmstadt: x= 7.9m, y= 4.9m
%                            rx position: x=6.75, y=0.01
%                            tx position 1: x=4.95, y=1.15
%                            tx position 2: x=7.2, y=2.5

clear all;

f=[2400 5230];          %MHz
points=1601;              %Points in sweep: Choose from: 3, 11, 21, 26, 51, 101, 201, 401, 801, 1601
points_multiplier=2;    %Number of sweeps to be performed per measurement


% --- list interfaces ---
%out=instrhwinfo
%out.SupportedInterfaces
% -----------------------

fclose('all');                               % Close all connections
obj1 = instrfind('Type', 'serial', 'Port', '/dev/tty.usbserial-PXWD4NPN', 'Tag', '');
% Create the serial port object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = serial('/dev/tty.usbserial-PXWD4NPN');
else
    fclose(obj1);
    obj1 = obj1(1);
end

% Configure instrument object, obj1.
set(obj1, 'Terminator', {'CR','CR'});
set(obj1, 'baudrate', 9600);
set(obj1, 'Timeout', 10);               % Configure instrument object, obj1
set(obj1, 'InputBufferSize', 100000);   % Configure instrument object, obj1
set(obj1, 'OutputBufferSize', 512);    % Configure instrument object, obj1

% Connect to instrument object, obj1.
fopen(obj1);

%-----------------------------------------------------------------------
% USB-GPIB adapter configuration
%-----------------------------------------------------------------------
fprintf(obj1,'++rst');              %reset GPIB adapter
fprintf(obj1,'++mode 1');           %controller mode
%device_mode=query(obj1,'++mode')
fprintf(obj1,'++addr 16');          %GPIB client address = 16
fprintf(obj1,'++eot_enable 1')      %append character at GPIB eoi
fprintf(obj1,'++eot_char 13')       %appended character= ASCII 13 (CR)
connection=query(obj1,'++ver')      %query adapter version


%-----------------------------------------------------------------------
% configure VNA:
%   -reset
%   -channel 1 for S21 phase
%   -channel 2 for S21 log magnitude
%   -select sweep domain & points
%-----------------------------------------------------------------------
fprintf(obj1,'PRES;')
device=query(obj1,'*IDN?')

fprintf(obj1,'chan1;')      %channel 1
fprintf(obj1,'auxcoff;')    %turn off aux channel
fprintf(obj1,'s21;')        %S21 parameter
fprintf(obj1,'pola;')       %polar measurement
fprintf(obj1,'powe10;')     %output power 10dBm
fprintf(obj1,'auto;')       %autoscale

% fprintf(obj1,'chan2;')
% fprintf(obj1,'auxcoff;')
% fprintf(obj1,'s21;')
% fprintf(obj1,'logm;')       %log magnitude measurement
% fprintf(obj1,'scal15;')     %15dB/div


fprintf(obj1,'duacoff;')     %dual channel display off
fprintf(obj1,strcat('poin',num2str(points),';'))     %Points in sweep: Choose from: 3, 11, 21, 26, 51, 101, 201, 401, 801, 1601


position_index=1
prompt='Measurement? [enter/n]';
x=input(prompt,'s');

while(isempty(x))

    for index_f=[1 : length(f)]
        fprintf(obj1,strcat('star',num2str(f(index_f)),'mhz;'))
        fprintf(obj1,strcat('stop',num2str(f(index_f)),'mhz;'))
        
        repeat_index=points_multiplier;
        multiple_measurements=[];
        
        while(repeat_index >=1)
            
            fprintf(obj1,'sing;')
            fprintf(obj1,'form4;')
            measurement=query(obj1,'outpform;');                       %Fetch data
            measurement=(str2double(strsplit(measurement,{',',char(10)})));   %Split data & convert
            measurement=reshape(measurement(1:end-1),2,[])'*[1; i];           %Make complex
            
            %check measurement validity
            standard_deviation=std(imag(measurement));
            %standard_deviation(standard_deviation==NaN)=0;
            max_rescaled_std=max(abs(standard_deviation'./mean(imag(measurement))));
            if (max_rescaled_std>0.1)
                disp(strcat('Possible measurement error! rescaled standard deviation ',num2str(max(standard_deviation)/mean(imag(measurement)))));
            end
            
            multiple_measurements= [multiple_measurements; measurement];
            repeat_index=repeat_index-1;
        end
        
        data(:,position_index,index_f)=multiple_measurements;
    end
    prompt='Next Measurement? [enter/n]';
    x=input(prompt,'s');
    position_index=position_index+1
end


position_index=position_index-1;
prompt='save? [y/n]';
x=input(prompt,'s');
if x=='y'
    prompt='file index?';
    x=input(prompt,'s');
    save(strcat(num2str(position_index),'_elements_',num2str(points*points_multiplier),'_sweep_points_',x,'.mat'),'data','f');
    disp(strcat(num2str(position_index),'_elements_',num2str(points*points_multiplier),'_sweep_points_',x,'.mat'));
end

fprintf(obj1,'powe0;')                              %output power 0dBm
fprintf(obj1,'++loc')                               %local control again
fclose(obj1);                                       %close connection
disp(strcat('Connection with',device,' closed.'));













% 
%     %extract data without header and comma end
%     %split on comma
%     %make doubles from ASCII
%     %sort real & imaginary components
% data=reshape(str2double(strsplit(data(8:end-1),',')),2,[]); 
% 
%     %remove outliers
% percentile=2;   %remove upper and lower 1%
% p=prctile(data(1,:),[percentile/4 100-percentile/4]);
% i=find(data(1,:)<p(1));
% data(:,i)=[];
% i=find(data(1,:)>p(2));
% data(:,i)=[];
% p=prctile(data(2,:),[percentile/4 100-percentile/4]);
% i=find(data(2,:)<p(1));
% data(:,i)=[];
% i=find(data(2,:)>p(2));
% data(:,i)=[];
% 
% %check measurement validity
% standard_deviation=std(data');
% standard_deviation(standard_deviation==NaN)=0;
% data_mean=mean(data,2)*1e-6;
% max_rescaled_std=max(abs(standard_deviation'./data_mean));
% if (max_rescaled_std>0.1)
%     disp(strcat('Possible measurement error! rescaled standard deviation ',num2str(max(standard_deviation))));
% else
%     disp('Measurement results seems valid');
% end
% 
% phase=mod(180/pi*(atan(data_mean(2)/data_mean(1))-pi*(data_mean(1)<0))+180,360)-180
% amplitude=20*log10(sqrt(data_mean(1)^2+data_mean(2)^2))
% 
% 
% phase_set(index)=phase;
% amplitude_set(index)=amplitude;
% 
% 
% prompt='Next Measurement? [enter/n]';
% x=input(prompt,'s');
% index=index+1;
% end
% 
% 
% index=index-1;
% prompt='save? [y/n]';
% x=input(prompt,'s');
% if x=='y'
%     prompt='file index?';
%     x=input(prompt,'s');
%     save(strcat('phase_mag_',num2str(index),'_measurements_',x,'.mat'),'phase_set','amplitude_set');
% end



