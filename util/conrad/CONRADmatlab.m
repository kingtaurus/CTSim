%%
%
% Class to initialize the java software CONRAD for its use in Matlab
%
%
%%
classdef CONRADmatlab < hgsetget
    
    properties
        isInitialized = 0;
        ConradPath = '';
    end
    
    
    methods
        function obj = CONRADmatlab(path)
            obj.ConradPath=path;
        end
        
        function obj = set.isInitialized(obj,value) % Handle class
            obj.isInitialized = value;
        end
        function value = get.isInitialized(obj) % Handle class
            value = obj.isInitialized;
        end
        
        function obj = set.ConradPath(obj,path)
            if (nargin < 2 || ~isdir(path))
                obj.ConradPath = uigetdir;
            else
                obj.ConradPath = path;
            end
        end
        
        function val = get.ConradPath(obj)
            val = obj.ConradPath;
        end
        
        function obj = initialize(obj)
            if (obj.isInitialized)
                clear java;
            end
            p = cd;
            lib = [obj.ConradPath '\CONRAD\lib'];
            cd(lib);
            list = cellstr(ls);
            for i=3:size(list,1)
                javaaddpath([lib '\' list{i}]);
            end
            cd(p);
            javaaddpath([obj.ConradPath '\CONRAD\src']);
            addpath([obj.ConradPath '\CONRAD\src']);
            javaaddpath([obj.ConradPath '\CONRAD']);
            addpath([obj.ConradPath '\CONRAD']);
            obj.isInitialized = 1;
            disp('Initialization done.');
        end
        
        
        function RecoPipeline(obj)
            if (~obj.isInitialized)
                obj.initialize;
            end
            import edu.stanford.rsl.*
            conrad.utils.CONRAD.setup();
            a=apps.gui.ReconstructionPipelineFrame();
            a.setVisible(1);
        end
        
        
    end
    
end