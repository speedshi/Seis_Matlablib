function [east,north,depth]=geod2cart(latitude,longitude,elevation)
% This function is used to convert geodetic coordinate (latitude,
% longitude, elevation) to Cartesian coordinate (North, East, Vertical down).
% For the axis of Vertical down, the sea-level is 0, and above the
% sea-level is negative, and below the sea-level is positive; can be
% interpret as depth relative to the sea-level.
%
% UTM coordinate system of wgs84Ellipsoid is used by default.
% Unit of latitude and longitude: degree.
% Unit of elevation: meter.
% Unit of output Cartesian coordinate (North, East and Vertical down): meter.
%
% INPUT--------------------------------------------------------------
% latitude: latitude in degree, can be scalar or a vector;
% longitude: longitude in degree, can be scalar or a vector;
% elevation: elevation in meter, can be scalar or a vector.
%
% OUTPUT-------------------------------------------------------------
% east: east component in meter;
% north: north component in meter;
% depth: vertical down component in meter.

latitude=latitude(:);
longitude=longitude(:);
elevation=elevation(:);

% set UTM zone
ut_zone = utmzone(mean(latitude,'omitnan'),mean(longitude,'omitnan')); % first latitude; second longitude

utmstruct = defaultm('utm');
utmstruct.zone = ut_zone;
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);


% convert coordinate. Note transfer elevation to depth by '*-1'
[east,north,depth] = mfwdtran(utmstruct,latitude,longitude,-elevation);
