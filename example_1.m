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

fprintf(visa_object, 'FREQ:START 2442 MHZ');
fprintf(visa_object, 'FREQ:STOP 2442 MHZ');
fprintf(visa_object, 'FORM:READ:DATA ASC');%INT,32  REAL,32 ASC
%query(visa_object, 'FORM:READ:DATA?');
%WERKT SOMS query(visa_object, 'SENSE:SWEEP:POINTS?')

data=query(visa_object, 'trac:data?');

    %extract data without header and comma end
    %split on comma
    %make doubles from ASCII
    %sort real & imaginary components
data=reshape(str2double(strsplit(data(8:end-1),',')),2,[]);
