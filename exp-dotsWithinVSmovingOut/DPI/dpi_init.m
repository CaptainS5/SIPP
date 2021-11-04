function dpi_init(dpi_set, photo)
% DPI INITIALIZATION - configurate analog input module:
%   requires Matlab Data Acquisition Toolbox (only available on Windows)!
%   PK 17/10/2019

if dpi_set.mode
    global DPI
    DPI = analoginput('nidaq','Dev1');                                      % DPI is an analog input, initialize input using data acquisition toolbox
    if photo.mode
        addchannel(DPI, [0 1 3]);                                           % if photdiode us used, we will record from channels 0-2
    else
        addchannel(DPI, 0:1);                                               % if photdiode us not used, we only need channels 0-1
    end
    set(DPI, 'SampleRate', 1000);                                           % set the sampling rate for the DPI and the Photodiode
    set(DPI, 'SamplesPerTrigger', inf);                                     % acquires data to buffer
end 

end

