function instrreset
%INSTRRESET Disconnect and delete all instrument objects.
%
%   INSTRRESET disconnects and deletes all instrument objects. If data 
%   is being written or read asynchronously, the asynchronous operation
%   is stopped.
%
%   An instrument object cannot be reconnected to the instrument after 
%   it has been deleted and should be removed from the workspace with
%   CLEAR.
%
%   See also INSTRHELP, ICINTERFACE/STOPASYNC, ICINTERFACE/FCLOSE,
%   ICINTERFACE/DELETE, ICDEVICE/DISCONNECT.
%

%   Copyright 1999-2017 The MathWorks, Inc.


 
% Delete all instrument objects
try
    instrument.internal.udm.InstrumentManager.getInstance.reset;
catch e
end

% Clear the list of Bluetooth devices on which instrhwinfo has been called.
try
    tempout = com.mathworks.toolbox.instrument.BluetoothDiscovery.clearInstrhwinofCalledOnceList();
catch e
end

% delete IVI-C class complaint objects created via TMTOOL
try
    com.mathworks.toolbox.instrument.browser.ivicWrapper.IviCInstrumentObjectStore.dispose();
catch e
end

try
   % Find all objects.  Return if none were found.   
   obj = instrfind;
   if isempty(obj)
      return;
   end
   
   delete(obj);
catch aException
   rethrow(aException);     
end



 