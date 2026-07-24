function XYZ = Beam2XYZ(Beam1, Beam2, Beam3, Beam4)
% convert beam velocity to XYZ velocity for Signature 1000
% Beam1/2/3/4 - Beam velocity in form [profile * cell#]

% Bofu Zheng
% Sept. 6 2018

theta = 25;
XYZ(:,1,:) = (Beam1 - Beam3)/sin(theta/180*pi)/2;     % from beam velocity to XYZ
XYZ(:,2,:) = (Beam4 - Beam2)/sin(theta/180*pi)/2;     % from beam velocity to XYZ
XYZ(:,3,:) = (Beam1 + Beam2 + Beam3 + Beam4)/cos(theta/180*pi)/4;  % z from beam velocity to XYZ
end