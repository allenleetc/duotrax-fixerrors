classdef SeqTypes 
  properties (Constant)
    type2prettystr = struct(...
      'birth','Track Birth',...
      'death','Track Death',...
      'swap','Match Cost Ambiguity',...
      'touch','Is Touching',...
      'jump','Large Jump',...
      'orientchange','Large Change in Orientation',...
      'orientvelmismatch','Velocity & Orient. Mismatch',...
      'largemajor','Large Major Axis');
%      'User Defined','user',
  end
  methods (Static)
    function prettytype = type2pretty(ty)
      smap = SeqTypes.type2prettystr;
      if isfield(smap,ty)
        prettytype = smap.(ty);
      else
        assert(strncmpi(ty,'user',4));
        usertype = ty(5:end);
        if isempty(usertype)
          usertype = '<unnamed>';
        end
        prettytype = sprintf('User Defined: %s',usertype);
      end
    end
  end
end