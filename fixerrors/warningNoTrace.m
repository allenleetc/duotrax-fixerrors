function warningNoTrace(varargin)

% Added by Allen Lee for DTFE in 2017

warnst = warning('off','backtrace');
warning(varargin{:});
warning(warnst);