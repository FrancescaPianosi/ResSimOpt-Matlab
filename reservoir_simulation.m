function [ s, Qreg, Qspill, E, S ] = reservoir_simulation( I, e, env_min, s0, s_min, s_max, operating_rule, x, delta )
%
% Perform reservoir simulation by implementing the mass balance equation:
%
% s(t+1) = s(t) + delta*( I(t) - E(t) - env_min(t) - Qreg(t) - Qspill(t) ) 
%
% where 's' is the reservoir storage (vol), 'I' the inflow (vol/time) 
% 'E' the evaporation (vol/time), 'env_min' the minimum release for the
% downstream environment (vol/time), 'Qreg' the flow abstracted for human
% use (vol/time), 'Qspill' the flow released through the emergency
% spillways (vol/time).
% The function assumes that all flow variables are expressed in units of
% volume over units of time (for example: m3/sec) and therefore 'delta' is 
% the length of the time-step (for example: delta = 60*60*24 sec if the
% time step is one day)
% 
% Inflow (I) and environmental flow (env_min) are given inputs.
%
% Evaporation (E) is computed as the product of the evaporation rate (e)
% per surface unit, and of the reservoir surface area (S). The evaporation
% rate is a given input, the surface area is derived from the storage via
% the embedded function 'storage_surface_function'. 
%
% Flow throught the emergency spillways (Qspill) is zero if the next-step
% storage estimated based on all other terms in the mass balance does not
% exceed the maximum capacity (s_max), otherwise it is set equal to the
% excess storage. This ensures that throughout the simulation the storage
% (s) never exceeds s_max. 
% 
% The regulated release (Qreg) is either a given input or it is calculated
% from the storage through an operating policy. Either way, its value
% is reduced, when needed, to prevent the next-step storage to go below
% s_min). Also notice that Qreg lumps into one both the flow that is
% abstracted (for example through an outlet tower) to be directly used for
% domestic or industrial use, and the flow that is released into the 
% downstream river (for instance as a way to drawdown the reservoir in 
% anticipation of an incoming flood). The downstream river also receives 
% the flow from the spillways, but that variable (Qspill) is kept separate 
% as it is not flow that is intentionally released, and the operator has no
% control on it.s
%
% USAGE:
%
% [ s, Qreg, Qspill, E, S ] = reservoir_simulation( I, e, env_min, ...
%                             s0, s_min, s_max, operating_rule, x, delta )
%
% Input:
% I = time series of reservoir inflows                       - vector (N,1)
% e = time series of evaporation per unit surface area       - vector (N,1)
% env_min = time series of environmental flows               - vector (N,1)
% s0    = reservoir storage at initial time step             - scalar
% s_min = minimum storage (under which releases are reduced) - scalar
% s_max = maximum storage (at which spills occur)            - scalar
% operating_rule = operating rule function (*)               - string
% x     = parameters of the operating rule (*)               - vector
% delta = length of the simulation time step                 - scalar
%
% Output:
% s      = time series of reservoir storage                - vector (N+1,1)
% Qreg   = time series of regulated release                - vector (N,1)
% Qspill = time series of release through spillways        - vector (N,1)
% E      = time series of evaporation                      - vector (N,1)
% S      = time series of reservoir surface area           - vector (N,1)
%
% (*) There are two options here:
% - set operating_rule to 'rel_seq': it means the simulator does not 
% actually use an operating policy but it uses a given release sequence
% prescribed by the user, and passed over to the function through the
% argument 'x'. Hence 'x' must be a vector (N,1)
% - set operating_rule to the name of the function (for instance: myrule.m)
% that defines the operating policy. Such function must have the form:
%     Qreg = myrule( s, param )
% where param is a vector of size (1,p) holding the parametes of the
% operating policy function. 
% The argument 'x' is then a matrix of size (N,p) holding the parameter
% sets at each time t over the simulation horizon (the matrix will have all
% identical rows if the same parameter set is used at all times).

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

N     = length(I)   ;

s      = nan*ones(N+1,1) ;
Qreg   = nan*ones(N,1) ;
Qspill = nan*ones(N,1) ; 
S      = nan*ones(N,1) ;
E      = nan*ones(N,1) ;

s(1)  = s0 ;

for t = 1:N
    
    % create a temporary variable that represent the available storage
    % at time 't'; we will use this variable to add and subtract inflows,
    % evaporation and environmental flow, before we can compute the amount
    % of water that can be abstracted from the reservoir or that must be
    % spilled
    s_ = s(t)+delta*I(t) ; % (m3) available storage after adding the inflow
    
    
    % ++++++++++++++++++++++++++++++++++++++++++++++
    % Compute and subtract evaporation E
    % ++++++++++++++++++++++++++++++++++++++++++++++
    % compute reservoir surface area:
    S(t) = storage_surface_function(s(t)); % (m2)
    % compute evaporation:
    E(t) = e(t)*S(t) ; % (m3/delta)
    % adjust computed evaporation:
    if s_-delta*E(t)<0 ; E(t)=s_/delta; end % this
    % adjustment is needed to avoid the (rare) case when the 
    % computed evaporation is larger than available storage, 
    % which would generate a negative value for the next storage 
    % even if the release was zero
    
    s_ = s_-delta*E(t) ; % (m3) available storage after subtracting evap
    
    % ++++++++++++++++++++++++++++++++++++++++++++++
    % Subtract minimum environmental flow 
    % ++++++++++++++++++++++++++++++++++++++++++++++
    if s_-delta*env_min(t)<0 ; env_min(t)=s_/delta; end
    s_ = s_-delta*env_min(t) ; % (m3) available storage after (...)
        
    % ++++++++++++++++++++++++++++++++++++++++++++++
    % Compute regulated release Qreg
    % ++++++++++++++++++++++++++++++++++++++++++++++
    switch operating_rule 
        case 'rel_seq' % use the given time series of release decision
            Qreg(t) = x(t) ; % (m3/s)
        otherwise % use the given operating rule
            Qreg(t) = feval(operating_rule,s_,x(t,:)) ; % (m3/s)
    end
    % adjust release to take into account
    % constraints on minimum feasible storage:
    if s_ - delta*Qreg(t)<s_min ; Qreg(t) = max(0, (s_-s_min)/delta ) ; end
    
    s_ = s_ - delta*Qreg(t) ; % (m3) available storage (...)
    
    % ++++++++++++++++++++++++++++++++++++++++++++++
    % Compute flow through spillways Qspill
    % ++++++++++++++++++++++++++++++++++++++++++++++
    
    if s_ > s_max ; Qspill(t) = (s_-s_max)/delta ; else Qspill(t)=0 ; end
    s_ = s_ - delta*Qspill(t) ; % (m3) available storage (...)

    % ++++++++++++++++++++++++++++++++++++++++++++++
    % Storage at next timestep
    % ++++++++++++++++++++++++++++++++++++++++++++++
    s(t+1) = s_ ; % (m3)    
    
end



