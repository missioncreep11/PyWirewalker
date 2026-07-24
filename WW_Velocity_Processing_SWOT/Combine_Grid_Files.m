folder = '/Volumes/NorseTPADS/NOPP2_ADCP/Grid/WW1/D1/';

ADCP_temp = load([folder '/NOPP2_WW1_D1_6.mat'],'ADCP');
fn = fieldnames(ADCP_temp.ADCP);
fn(strcmp(fn,'Nav'))=[];
fn(strcmp(fn,'Notes'))=[];
fn_nav = fieldnames(ADCP_temp.ADCP.Nav);

ADCP = struct();
for f = fn';f=f{:};
    ADCP.(f) = ADCP_temp.ADCP.(f);
end
for f = fn_nav';f=f{:};
    ADCP.Nav.(f) = ADCP_temp.ADCP.Nav.(f);
end

for n = 2:6
    ADCP_temp = load([folder '/NOPP2_WW1_D1_' num2str(n) '.mat']);
    for f = fn';f=f{:};
        ADCP.(f) = [ADCP.(f) ADCP_temp.ADCP.(f)];
    end
    for f = fn_nav';f=f{:};
        ADCP.Nav.(f) = [ADCP.Nav.(f) ADCP_temp.ADCP.Nav.(f)];
    end
end
