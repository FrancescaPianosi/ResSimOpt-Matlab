function  [Rel_reg, Def_reg, Vul_reg, Rel_down, Qmax_down] = compute_res_perf( Qreg, xreg, Qdown, xdown, varargin )
%
% [Rel_reg, Def_reg, Vul_reg, Rel_down, Qmax_down] = compute_res_perf( Qreg, xreg, Qdown, xdown )
% 
% compute several commonly used indicators of reservoir performance for 
% water supply:
% - Rel_reg: frequency of failing to release the target
% - Def_reg: mean volume of release deficit wrt target
% - Vul_reg: mean *squared* volume of release deficit wrt target
% flood control:
% - Rel_down: frequency of exceeding max downstream flow
% - Qmax_down: maximum downstream flow peak
% 
% Input:
% Qreg: time series of regulated releases                    - vector (T,1)
% xreg: time series of release targets                       - vector (T,1)
%       (we would like Qreg >= xreg)
% Qdown: time series of total flow in downstream river       - vector (T,1) 
%       (does not coincide with Qreg as Qreg may not include spills, or if
%       for example Qreg is directly diverted into a canal/pipe and not
%       released in the downstream river)
% xdown: time series of target flow not to be exceeded       - vector (T,1)
%       (we would like Qdown <= xdown)

T = length(Qreg) ;
idx = Qreg < xreg  ;
Rel_reg = sum(idx)/T ;
Def_reg = sum( xreg(idx) - Qreg(idx) )/T;
Vul_reg = sum( (xreg(idx) - Qreg(idx)).^2 )/T;

idx = Qdown > xdown  ;
Rel_down   = sum(idx)/T ;
Qmax_down = max( Qdown );


if nargin == 5 ; flag = varargin{1}; else flag = 1 ; end
if flag > 0 ;
fprintf('\n')
if ~isempty(Qreg) ; fprintf('Frequency of failing to release the target: %g \n',Rel_reg); end
if ~isempty(Qreg) ; fprintf('Mean volume of release deficit wrt target: %g \n',Def_reg); end
if ~isempty(Qreg) ; fprintf('Mean *squared* volume of release deficit wrt target: %g \n',Vul_reg); end
if ~isempty(Qdown); fprintf('Frequency of exceeding max downstream flow: %g \n',Rel_down); end
if ~isempty(Qdown); fprintf('Maximum downstream flow peak: %g \n',Qmax_down); end
fprintf('\n')
end
