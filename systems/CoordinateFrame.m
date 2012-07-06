classdef CoordinateFrame < handle
% Every input, state, and output in a DynamicalSystem has a coordinate frame
% attached to it.  Many bugs can be avoided by forcing developers to be 
% explicit about these coordinate systems when they make combinations of
% systems.
  
  properties (SetAccess=private,GetAccess=public)
    name;           % string name for this coordinate system
    dim;            % scalar dimension of this coordinate system
    transforms={};  % handles to CoordinateTransform objects

    coordinates={}; % list of coordinate names

    angle_flag=[];  % angle_flag(i)=true iff variable i wraps around 2pi
    poly=[];        % optional msspoly variables for this frame
    prefix;
  end
  
  methods
    function obj=CoordinateFrame(name,dim,prefix)
      typecheck(name,'char');
      obj.name = name;
      
      typecheck(dim,'double');
      sizecheck(dim,[1 1]);
      obj.dim = dim;
      
      if (nargin<3)
        ind = strfind(name,':');
        if isempty(ind)
          prefix = name(1);
        else
          prefix = name(ind(end)+[1:2]);
        end
      else
        typecheck(prefix,'char');
        sizecheck(prefix,[1 1]);
      end
      obj.prefix = prefix;
      
      ind=1;
      function str=coordinateName(~);
        str=[prefix,num2str(ind)];
        ind=ind+1;
      end
      obj.coordinates=cellfun(@coordinateName,cell(dim,1),'UniformOutput',false);
      
      if checkDependency('spot_enabled') && dim>0
        if (prefix=='t') error('oops.  destined for a collision with msspoly representing time'); end
        obj.poly = msspoly(prefix,dim);
      end
      
      obj.angle_flag = repmat(false,dim,1);
    end
    
    function obj=addTransform(obj,transform)
      typecheck(transform,'CoordinateTransform');
      if (getInputFrame(transform) ~= obj || getOutputFrame(transform) == obj)
        error('transform must be from this coordinate frame to another to be added');
      end
      if ~isempty(findTransform(obj,getOutputFrame(transform)))
        error('i already have a transform that gets me to that frame');
      end
      obj.transforms{end+1}=transform;
    end
    
    function obj=updateTransform(obj,newtransform)
      % find a current transform with the same target frame and replace it
      % with the new transform
      typecheck(newtransform,'CoordinateTransform');
      ind=find(cellfun(@(a)getOutputFrame(a)==newtransform.getOutputFrame,obj.transforms));
      if (isempty(ind))
        error('no transform to update');
      elseif length(ind)>1
        error('found multiple transforms');
      else
        obj.transforms{ind}=newtransform;
      end
    end
    
    function tf=findTransform(obj,target,throw_error_if_fail)
      typecheck(target,'CoordinateFrame');
      ind=find(cellfun(@(a)getOutputFrame(a)==target,obj.transforms));
      if (isempty(ind))
        tf=[];
        if (nargin>2 && throw_error_if_fail)
          error(['Could not find any transform between ',obj.name,' and ', target.name]);
        end
      elseif length(ind)>1
        error(['Found multiple transforms between ',obj.name,' and ', target.name]);
      else
        tf=obj.transforms{ind};
      end
    end
    
    function obj=setAngleFlags(obj,flags)
      typecheck(flags,'logical');
      sizecheck(flags,[obj.dim,1]);
      obj.angle_flag = flags;
    end
    
    function obj=setCoordinateNames(obj,cnames)
      if (iscell(cnames) && isvector(cnames) && length(cnames)==obj.dim && all(cellfun(@ischar,cnames)))
        obj.coordinates=cnames;
      else
        error('cnames must be a cell vector of length dim populated with strings'); 
      end
    end
    
    function fr=subFrame(obj,dims)
      if ~isnumeric(dims) || ~isvector(dims) error('dims must be a numeric vector'); end
      if (any(dims>obj.dim | dims<1)) error(['dims must be between 1 and ',obj.dim]); end
      fr = CoordinateFrame([obj.name,mat2str(dims)], length(dims), obj.prefix);
      fr.coordinates = obj.coordinates(dims);
      fr.angle_flag = obj.angle_flag(dims);
      fr.poly = obj.poly(dims);
    end
  end

  % todo: consider putting LCM encode/decode in here
end
