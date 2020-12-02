function Op_param = op_piecewise_linear_transform(x,idx,s_max,Qtarget,delta)
%
% Op_param = op_piecewise_linear_transform(x,idx,s_max,Qtarget,delta)
%
% Takes set of decision variables 'x' and transform them into the
% parameters of the piece-wise linear operating policy implemented in the
% function 'op_piecewise_linear'.
% x       - vector of decision variables  - vector (G*3,1)
% idx     - time series of indices to attribute different parameter sets
%           over time                     - vector (T,1)
% s_max   - maximum reservoir storage     - scalar
% Qtarget - time series of target outflow - vector (T,1)
% delta   - scalar
% Op_param: time series of parameters of the operating policy - matrix (T,6)

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

T = length(idx) ; % number of time-steps in the simulation horizon
G = length(unique(idx)) ; % number of temporal groups within which the
% policy parameters are the same
x = reshape( x, 3, G )'; % (G,3)

% % Bring back parameters 'x' from uniform [0,1] (range of variations for all
% % decision variables in the nsga framework) to their actual variation
% % ranges: [0,pi/2) for x(1) and x(3) and [0,s_max] for x(2) 
% x(:,1) = x(:,1)*pi/2  ;
% x(:,2) = x(:,2)*s_max ;
% x(:,3) = x(:,3)*pi/2  ;


% create the matrix of policy parameters:
Op_param = nan(T,6) ;

% fill in the first 3 columns with the (time-varying) parameters
% optimised by the nsga algorithm
for j=1:3 % for each parameter
    for i=1:G % for each time group
        Op_param(idx==i,j) = x(i,j) ;
    end
end

% fill in the other 3 columns with the (constant) parameters
% that are not varied by the nsga algorithm
Op_param(:,4) = Qtarget*delta ;
Op_param(:,5) = s_max    ;
Op_param(:,6) = delta    ;

