classdef SeqStatus
  properties
    tablestr
  end
  enumeration
     UNKNOWN('unreviewed')
     CORRECT('correct') % correct, unknown how
     CORRECTFIXED('fixed') % actively/manually fixed
     CORRECTUNCHANGED('correct') % was correct when checked, unchanged
     CORRECTOTHER('correct other') % became correct thru side-effect of other change
  end
  methods 
    function obj = SeqStatus(str)
      obj.tablestr = str;
    end
    function tf = isCorrect(obj)
      tf = obj~=SeqStatus.UNKNOWN;
    end
  end
end
     
     
     
     
    