function varargout = fileparts2(name)

% Added by Allen Lee for DTFE in 2017

[p,f,e] = fileparts(name);
varargout{1} = p;
if nargout > 2
  varargout{2} = f;
  varargout{3} = e;
else
  varargout{2} = [f e];
end
