function ENU = XYZ2ENU(velx, vely, velz, P, R, H,varargin)
% convert xyz velocity to ENU velocity for Signature 1000
% velx/y/z - xyz velocity in form [profile# * cell#]
% R - Roll angle (Rotate around X) in degree
% P - Pitch angle (Rotate around Y) in degree
% H - Yaw angle (Rotate around Z) in degree
%
% ENU = XYZ2ENU(velx, vely, velz, P, R, H)

% Bofu Zheng
% Sept. 6 2018

% from xyz to ENU
hh = pi*(H - 90)/180;     % Heading-90 degrees for ENU velocity
pp = pi*P/180;
rr = pi*R/180;

for j = 1:length(P)
    % Make heading matrix
    H = [cos(hh(j)) sin(hh(j)) 0; -sin(hh(j)) cos(hh(j)) 0; 0 0 1];
    
    % Make tilt matrix
    P = [cos(pp(j)) -sin(pp(j))*sin(rr(j)) -cos(rr(j))*sin(pp(j));...
        0             cos(rr(j))          -sin(rr(j));  ...
        sin(pp(j)) sin(rr(j))*cos(pp(j))  cos(pp(j))*cos(rr(j))];
    
    % Make resulting transformation matrix
    R = H*P;
    
    % calculate ENU
    if isempty(varargin)
        XYZ = [velx(j,:); vely(j,:); velz(j,:)];
        velenu = R * XYZ;      % XYZ velocity to ENU
    elseif strcmp(varargin{1},'reverse')
         % calculate ENU
        XYZ = [velx(j,:); vely(j,:); velz(j,:)];
        velenu = transpose(transpose(XYZ) * R);      % XYZ velocity to ENU
    end
    
    ENU(j,1,:) = velenu(1,:);  % east-west
    ENU(j,2,:) = velenu(2,:);  % north-south
    ENU(j,3,:) = velenu(3,:);  % up-down
end

end
