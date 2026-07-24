function WWmeta = SetupPath(WWmeta,MainPath,varargin)
%%
    if ~strcmp(MainPath(end),'/')
        MainPath = [MainPath '/'];
    end
    
    % Setup main directories
    WWmeta.matpath = [MainPath 'Combined/'];
     % path to the save processed data - profile (separated into up/down casts)
    WWmeta.propath = [MainPath 'Profile/'];
    % path to the copied profiles - for combining cut-off profiles
    WWmeta.propath_rearrange = [MainPath 'ReOrdered/'];
    % path to the save estimated velocity field
    WWmeta.gridpath = [MainPath 'Grid/'];
    % (optional) path to save the profiling figure, can be commented out
    WWmeta.figpath = [MainPath 'Fig/'];
    
    % Add Sub-Directories
    if nargin>2
        for n = 1:nargin-2
            WWmeta.matpath = [WWmeta.matpath varargin{n} '/'];
            WWmeta.propath = [WWmeta.propath varargin{n} '/'];
            WWmeta.propath_rearrange = [WWmeta.propath_rearrange varargin{n} '/'];
            WWmeta.gridpath = [WWmeta.gridpath varargin{n} '/'];
            WWmeta.figpath = [WWmeta.figpath varargin{n} '/'];
        end
    end
%%    Make directories if they dont already exist
    if ~isfolder(WWmeta.matpath);mkdir(WWmeta.matpath);end
    if ~isfolder(WWmeta.propath);mkdir(WWmeta.propath);end
    if ~isfolder(WWmeta.propath_rearrange);mkdir(WWmeta.propath_rearrange);end
    if ~isfolder(WWmeta.gridpath);mkdir(WWmeta.gridpath);end
    if ~isfolder(WWmeta.figpath);mkdir(WWmeta.figpath);end
end

