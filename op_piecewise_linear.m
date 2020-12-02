function u = op_piecewise_linear( s, x )
%
% u = op_piecewise_linear( s, x )
%
% Piece-wise linear operating policy
%
% s = storage value (volume)
% u = regulated release (volume/time)
% x = vector of 6 policy parameters as follows:
%       x(1): slope of first linear piece (radiant)
%       x(2): storage at which second linear piece starts (volume)
%       x(3): slope of second linear piece (radiant)
%       x(4): target release (volume)
%       x(5): maximum storage (volume)
%       x(6): length of the simulation time-step (time)
% Note: 
% x(1) and x(2) must vary in (0,pi/2)
% x(2) and x(5) must be in the same units, and x(5)>x(2)
% x(4) must also be in the same volumetric units as x(2) and x(5), so it is
% the overall target release over the simulation time-step
% x(6) is used to convert 'u' into flow units (volume/time)
%
% Example:
% (for a case when volumes are in Ml, flows are in Ml/day,
% and the simulation time step is 1 week) 
% x(1) = pi/8  ; % (rad)
% x(2) = 400   ; % (Ml)
% x(3) = pi/16 ; % (rad)
% x(4) = 14    ; % (Ml/week) - this correspond to a target of 2 Ml/day
% x(5) = 500   ; % (Ml)
% x(6) = 7     ; % (days in a week)  
% s_ = [0:x(5)]; for i=1:length(s_); u_(i) = op_piecewise_linear( s_(i), x ) ; end
% figure; plot(s_,u_,'k'); xlabel('storage (Ml)'); ylabel('flow (Ml/day)')

%  Copyright (c) 2020, Francesca Pianosi
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.

% Recover policy parameters:
ms    = x(5) ;
tr    = x(4) ;
delta = x(6) ;

% Check paramter values are feasible:
if x(1)<=0    ; error('First component of input ''x'' (1st slope angle) must be positive'); end 
if x(1)>=pi/2 ; error('First component of input ''x'' (1st slope angle) must be smaller than pi/2'); end 
if x(2)>ms    ; error('Second component of input ''x'' cannot be larger than fifth component (max storage)'); end 
if x(3)<=0    ; error('Third component of input ''x'' (2nd slope angle) must be positive'); end 
if x(3)>=pi/2 ; error('Third component of input ''x'' (2nd slope angle) must be smaller than pi/2'); end 

% calculate outflow according to policy:
Si = [ 0 tr/tan(x(1)) x(2)                     ms ]' ;
Ui = [ 0           tr  tr  tr+(ms-x(2))*tan(x(3)) ]' ;
u  = interp_lin_scalar(Si,Ui,s) / delta ;

function y = interp_lin_scalar( X , Y , x )
%
% y = interp_lin_scalar( X , Y , x )
%
%            Y(k+1) + Y(k)
% y = Y(k) + ------------- ( x - X(k) ) 
%            X(k+1) - X(k)
%
% with 'k' such that X(k) <= x < X(k+1)
%
% input :
% X = vector of independent variables - ( n , 1 )
% Y = vector of dependent variables   - ( n , 1 )
% x = scalar                          - ( 1 , 1 )
%
% output :
% y = scalar                          - ( 1 , 1 )
%
% Last update : Francesca 09/11/2011

% -------------
% extreme cases
% -------------
%if x <= X( 1 ) ; y = Y( 1 ) ; return ; end
%if x >= X(end) ; y = Y(end) ; return ; end
if x > X(end) ; 
    y = Y(end) + (Y(end)-Y(end-1))/(X(end)-X(end-1))*(x-X(end)) ; 
    return ;
elseif x < X(1)
    y = Y(1) + (Y(2)-Y(1))/(X(2)-X(1))*(x-X(1)) ; 
    return ;
end

% -------------
% otherwise
% -------------

% Find index 'k' of subinterval [ X(k) , X(k+1) ] s.t. X(k) <= x < X(k+1)
[ ignore , i ] = min( abs( X - x ) ) ;

% If X( i ) = x     then   y = Y( i ) :
if X( i ) == x ; y = Y( i ) ; return ; end

% Else :
% if X( i ) < x     then   k = i  
% if X( i ) > x     then   k = i - 1
k = i - ( X( i ) > x )   ;       
% Line joining points ( X(k) , Y(k) ) and ( X(k+1) , Y(k+1) ) 
Dy = Y( k + 1 ) - Y( k ) ;
Dx = X( k + 1 ) - X( k ) ;
m  = Dy / Dx             ; % slope
% Interpolate :
y = Y( k ) +  m * ( x - X( k ) ) ;