function dpi_saveData(filename,blockFolder)   
% SAVE DPI DATA Saving the collected analog input.
%   PK, 22/10/2019
global DPI
[data,t] = getdata(DPI,DPI.SamplesAvailable);
data = cat(2,t,data);

save([blockFolder,filename],'data')


end

