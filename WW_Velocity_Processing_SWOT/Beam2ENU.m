function ENU = Beam2ENU(Beam1,Beam2,Beam3,Beam4, P, R, H)
% convert beam velocity to ENU velocity for Signature 1000
% Beam1/2/3/4 - Beam velocity in form [profile * cell#]
% R - Roll angle (Rotate around X) in degree
% P - Pitch angle (Rotate around Y) in degree
% H - Yaw angle (Rotate around Z) in degree

% Bofu Zheng
% Sept. 6 2018

theta = 25;   % Beam angle

% from beam to ENU
hh = pi*(H - 90)/180;     % Heading-90 degrees for ENU velocity
pp = pi*P/180;
rr = pi*R/180;

% Make xyz matrix
T =[1/sin(theta/180*pi)/2  0  -1/sin(theta/180*pi)/2 0;...
    0 -1/sin(theta/180*pi)/2  0 1/sin(theta/180*pi)/2 ;...
    1/cos(theta/180*pi)/4 1/cos(theta/180*pi)/4 1/cos(theta/180*pi)/4 1/cos(theta/180*pi)/4];

for j = 1:length(P)
    % Make heading matrix
    HM = [cos(hh(j)) sin(hh(j)) 0; -sin(hh(j)) cos(hh(j)) 0; 0 0 1];
    
    % Make tilt matrix
    PM = [cos(pp(j)) -sin(pp(j))*sin(rr(j)) -cos(rr(j))*sin(pp(j));...
        0             cos(rr(j))          -sin(rr(j));  ...
        sin(pp(j)) sin(rr(j))*cos(pp(j))  cos(pp(j))*cos(rr(j))];
    
    % Make resulting transformation matrix
    R = HM*PM*T;
    
    % calculate ENU
    Beam = [Beam1(j,:); Beam2(j,:); Beam3(j,:); Beam4(j,:)];
    velenu = R * Beam;      % Beam velocity to ENU
    ENU(j,1,:) = velenu(1,:);
    ENU(j,2,:) = velenu(2,:);
    ENU(j,3,:) = velenu(3,:);
end
end