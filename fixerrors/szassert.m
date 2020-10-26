function szassert(x,sz,varargin)

% Added by Allen Lee for DTFE in 2017

assert(isequal(size(x),sz),varargin{:});

