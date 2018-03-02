function verifyTrx(trxfilename,abstol,reltol)
% Perform various checks on a trxfile

if exist('abstol','var')==0
  abstol = 1e-6;
end
if exist('reltol','var')==0
  reltol = 1e-3;
end

trx = load(trxfilename,'-mat');
if ~isfield(trx,'trx')
  error('Trxfile ''%s'' does not contain a ''trx'' variable.');
end
trx = trx.trx;

for fly=1:numel(trx)
  trxEl = trx(fly);
  x = trxEl.x;
  y = trxEl.y;
  a = trxEl.a;
  th = trxEl.theta;
  wal = trxEl.wing_anglel;
  war = trxEl.wing_angler;
  xwl = trxEl.xwingl;
  ywl = trxEl.ywingl;
  xwr = trxEl.xwingr;
  ywr = trxEl.ywingr;

  if any(wal>0)
    nf = nnz(wal>0);    
    warningNoTrace('Fly %d: wing_anglel>0 in %d frames. Maximum val: %.3f deg',...
      fly,nf,max(wal)/pi*180);
  end
  if any(war<0)
    nf = nnz(war<0);
    warningNoTrace('Fly %d: wing_angler<0 in %d frames. Minimum val: %.3f deg',...
      fly,nf,min(war)/pi*180);
  end
  
  waAbsL = th + pi + wal;
  waAbsR = th + pi + war;
  xwl2 = x + 4*a.*cos(waAbsL);
  ywl2 = y + 4*a.*sin(waAbsL);
  xwr2 = x + 4*a.*cos(waAbsR);
  ywr2 = y + 4*a.*sin(waAbsR);

  dxwl = abs(xwl-xwl2);
  dxwlrel = dxwl./abs(xwl);
  dywl = abs(ywl-ywl2);
  dywlrel = dywl./abs(ywl);
  dxwr = abs(xwr-xwr2);
  dxwrrel = dxwr./abs(xwr);
  dywr = abs(ywr-ywr2);
  dywrrel = dywr./abs(ywr);

  tfSusp = dxwl>abstol | dxwlrel>reltol;
  if any(tfSusp)
    nf = nnz(tfSusp);    
    warningNoTrace('Fly %d: wing_anglel/xwingl inconsistency in %d frames.',fly,nf);
  end  
  tfSusp = dywl>abstol | dywlrel>reltol;
  if any(tfSusp)
    nf = nnz(tfSusp);    
    warningNoTrace('Fly %d: wing_anglel/xwingl inconsistency in %d frames.',fly,nf);
  end
  tfSusp = dxwr>abstol | dxwrrel>reltol;
  if any(tfSusp)
    nf = nnz(tfSusp);    
    warningNoTrace('Fly %d: wing_angler/xwingr inconsistency in %d frames.',fly,nf);
  end
  tfSusp = dywr>abstol | dywrrel>reltol;
  if any(tfSusp)
    nf = nnz(tfSusp);    
    warningNoTrace('Fly %d: wing_angler/ywingr inconsistency in %d frames.',fly,nf);
  end
end
  

  