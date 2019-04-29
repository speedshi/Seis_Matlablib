function wtraces2segy(fname,data,dt,soup,recp)
% This function is used to call the SegyMAT package to output traces in
% segy format. Note the input data is in nt*nrec, i.e. 1st dimension is
% time samples, 2ed dimension is traces. Distance is expressed in meter.
% NOTE here for the coordinate system, Z is vertical down, thus the
% 'ReceiverGroupElevation' is negative.
% INPUT---------------------------------------------------------------------------------
% fname: file name of the output segy file;
% data: seismic traces, nt*nrec;
% dt: time sampling intervals, in second;
% soup: source positions, 1*3, X-Y-Z;
% recp: receiver positions, nrec*3, X-Y-Z, each row indicates a receiver;
% OUTPUT-----------------------------------------------

[nt,nrec]=size(data); % obtain the total number of time samples and traces

itr=1:nrec; % index of the traces

% generate and output the segy file
WriteSegy(fname,data,'ns',nt,'dt',dt,'TraceNumber',itr);

nn=ones(nrec,1);

% edit the trace headers
WriteSegyTraceHeaderValue(fname,nn,'key','TraceIdenitifactionCode'); % type of the data: seismic data
WriteSegyTraceHeaderValue(fname,0*nn,'key','SourceSurfaceElevation'); % surface elevation at source
WriteSegyTraceHeaderValue(fname,soup(1)*nn,'key','SourceX'); % X coordinate of the source
WriteSegyTraceHeaderValue(fname,soup(2)*nn,'key','SourceY'); % Y coordinate of the source
WriteSegyTraceHeaderValue(fname,soup(3)*nn,'key','SourceDepth'); % source depth below surface
WriteSegyTraceHeaderValue(fname,recp(:,1),'key','GroupX'); % X coordinate of the receivers
WriteSegyTraceHeaderValue(fname,recp(:,2),'key','GroupY'); % Y coordinate of the receivers
WriteSegyTraceHeaderValue(fname,-recp(:,3),'key','ReceiverGroupElevation'); % elevations of the receivers
WriteSegyTraceHeaderValue(fname,nn,'key','CoordinateUnits'); % coordinate units: use length

% edit the binary file header
[~,SegyTraceHeaders,SegyHeader]=ReadSegy(fname);
SegyHeader.TraceSorting=1; % set the trace sorting code: as recorded (no sorting)
SegyHeader.MeasurementSystem=1; % set the unit of distance in meter.
SegyHeader.TextualFileHeader=32*ones(3200,1); % clear the textual file header

% ensamble all the headers and data
WriteSegyStructure(fname,SegyHeader,SegyTraceHeaders,data);

end