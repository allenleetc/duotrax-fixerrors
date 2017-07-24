classdef SeqTable < handle
  
  properties
    hParent
    jtable
    cbkSelectSeq
  end
  
  properties (Dependent)
    nRows
  end
  
  methods
    function v = get.nRows(obj)
      v = size(obj.jtable.Data,1);
    end
  end
  
  methods
    function obj = SeqTable(hParent,pos,posunits,cbk)
      obj.hParent = hParent;
      jt = uiextras.jTable.Table(...
        'parent',hParent,...
        'SelectionMode','discontiguous',...
        'ColumnName',{'type' 'startframe' 'length' 'status'},...
        'ColumnPreferredWidth',[200 200 200 200],...
        'Editable','off');
      jt.Units = posunits;
      jt.Position = pos;
      jt.Units = 'normalized'; % so will confirm to parent
      jt.MouseClickedCallback = @(src,evt)obj.cbkClick(src,evt);
      obj.jtable = jt;
      obj.cbkSelectSeq = cbk;
    end
  end
  
  methods
    
    function setSeqData(obj,seqs)
      t = struct2table(seqs);
      jt = obj.jtable;
      
      startframe = cellfun(@(x)x(1),t.frames);
      startframeF = floor(startframe);
      if ~isequal(startframeF,startframe)
        warningNoTrace('SeqTable:round','Encountered non-integer startframes.');
      end
      startframeF = arrayfun(@num2str,startframeF,'uni',0);
      len = cellfun(@numel,t.frames);
      len = arrayfun(@num2str,len,'uni',0);
      
      dat = [t.type startframeF len {t.status.tablestr}'];
      jt.Data = dat;
    end
    
    function updateStatus(obj,irow,status)
      assert(isa(status,'SeqStatus'));
      jt = obj.jtable;
      dat = jt.Data;
      tf = strcmp(jt.ColumnName,'status');
      assert(nnz(tf)==1);
      dat{irow,tf} = status.tablestr;
      jt.Data = dat;
    end
    
    function irows = getSelectedRows(obj)
      jt = obj.jtable;
      irows = sort(jt.SelectedRows);
    end
    
    function setSelectedRows(obj,irow)
      assert(isscalar(irow));
      
      jt = obj.jtable;
      nrows = size(jt.Data,1);
      if 0<irow && irow<=nrows
        jt.SelectedRows = irow;
      else
        jt.SelectedRows = [];
      end
    end
    
    function cbkClick(obj,src,evt) %#ok<INUSD>
      irow = obj.getSelectedRows();
      if numel(irow)>1
        warning('SeqTable:sel',...
          'Multiple rows of table selected. Using first selected row.');
        irow = irow(1);
      end
      try
        obj.cbkSelectSeq(irow);
      catch ME
        disp(ME.message);
      end
    end
    
  end
  
end