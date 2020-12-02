% This script shows an example of how to use the 'reservoir_simulation'
% function. It guides the user through the simulation of a reservoir system
% in three cases: 
% Case 1 - where a very simple operations strategy is assumed, whereby the
% operator always try to simply release the target release. This case
% provides an example of using the 'reservoir_simulation' function using
% a release sequence
% Case 2 - where the release decision is a function of the current storage
% of the reservoir, whereby the target release is returned when the storage
% is in its mid range, and a smaller/larger release is returned when the
% storage is particularly low/high. This case provides an example of using
% the 'reservoir_simulation' function using an operating policy
% Case 3 - an extension of Case 2 to illustrate the case where the
% parameters of the operating policy vary over the year.
% (for a discussion of the difference between using a release sequence and
% an operating policy, see Dobson et al (2019))

%% Load data and define system characteristics
% Load time series of reservoir inflows and
% evaporation rate (per unit surface area):
load Data
T = 365 ; % length of time series (days)
I = Data(1:T,1) ; % (m3/s)
e = Data(1:T,2) ; % (m/s)

% The time series above are daily records for one hydrological year
% from October to September of the next calendar year. So we create now two
% variables that we will use to customise the plot axes and label them with
% the initial of the months
xtick_month = cumsum([1 31 30 31 31 28 31 30 31 30 31 30 31]);
label_month = {'O','N','D','J','F','M','A','M','J','J','A','S'} ;

% Plot inflow and unit evaporation data:
figure; 
subplot(211); plot(I); ylabel('inflow (m3/s)')
set(gca,'XTick',xtick_month,'XTickLabel',label_month)
subplot(212); plot(e); ylabel('unit evap. (m3/s)')
set(gca,'XTick',xtick_month,'XTickLabel',label_month)

% Define minimum environmental flow, initial storage, minimum and maximum
% storage, and the simulation time-step length
env_min = 0.5*ones(T,1)  ; % (m3/s) assumed constant over the time series
s_min = 1*10^7  ; % (m3)
s_max = 1*10^8  ; % (m3)
s0    = 0.9*s_max  ; % (m3)
delta = 60*60*24 ; % (sec/day)

%% CASE 1: 
% First, let's consider a case when we want to always release a constant
% amount over the simulation period. This constant amount is the 'target
% release' that would cover all downstream demand for water, for instance
% for domestic use and/or irrigation

% Let's define the target release 
tr = 7 ; % (m3/s) the target flow we would ideally always release

% Run simulation:
operating_rule = 'rel_seq' ;
op_param = tr*ones(T,1) ; % (m3/s) as we always want to release the same target
% we create a release sequence of the length of the simulation horizon  
[ s, Qreg, Qspill, E ] = reservoir_simulation( I, e, env_min, ...
                             s0, s_min, s_max, operating_rule, op_param, delta ) ;
                         
% Plot results:
flow_label =  'flow (m3/s)' ; stor_label = 'storage (m3)' ; time_label = 'time (days)';
leg = {'Inflow','Evap','Release','Spills','Min Env Flow'} ;
plot_res_timeseries( s, s_min, s_max, I, E, Qreg, Qspill, env_min, leg,flow_label,stor_label , time_label ) ;
set(gca,'XTick',xtick_month,'XTickLabel',label_month)

% Quantify reservoir performance:
[Rel_reg, Def_reg, Vul_reg] = compute_res_perf( Qreg, tr*ones(T,1), [], [] );

%% CASE 2: 
% now let's consider a case when we still want to release the target demand
% but we would also like to (1) apply some hedging (that is, an intentional
% reduction of the release - even if it would still be feasible to release
% the target demand - aimed at saving more water and thus facing smaller
% deficits at later time); and (2) attenuate downstream peak flows for
% flood control purpose. In order to achieve this, instead of releasing a
% constant amount we will relase an amount proportional to the storage
% value: lower than the target release when the storage is very low, higher
% than the target release when the storage is very high, and equal to the
% target only when the storage is in its mid range. We model this operator
% behaviour through an operating policy, that is, a function (in this case
% piece-wise linear) that returns the regulated release depending on the
% storage value. 

% Let's define the parameters of the operating policy and visualise it:
x1 = pi/ 180    ; % (radiant) slope of first linear piece
x2 = 0.9*s_max  ; % (m3) storage at which second linear piece starts
x3 = pi / 32    ; % (radiant) slope of second linear piece
op_param = [ x1, x2, x3, tr*delta, s_max, delta ]  ;
s_test = [0:.25*10^6:s_max];
for i=1:length(s_test); u_test(i) = op_piecewise_linear( s_test(i), op_param ) ; end
figure; plot(s_test,u_test,'k'); xlabel(stor_label); ylabel(flow_label)

Q_max_down = 50 ; % (m3/s) maximum flow that should (ideally) be released 
% downstream (considering both regulated releases and spills)  

% Run simulation:
operating_rule = 'op_piecewise_linear';
Op_param = repmat(op_param,T,1) ; % create a vector of the policy parameters over the 
% simulation horizon (in this case we use the same set of parameters at
% every time-step)
[ s, Qreg, Qspill, E ] = reservoir_simulation( I, e, env_min, ...
                             s0, s_min, s_max, operating_rule, Op_param, delta ) ;

% Plot results:
plot_res_timeseries( s, s_min, s_max, I, E, Qreg, Qspill, env_min, leg,flow_label,stor_label , time_label ) ;
set(gca,'XTick',xtick_month,'XTickLabel',label_month)

% Quantify reservoir performance:
[Rel_reg, Def_reg, Vul_reg, Rel_spill, Qmax_spill] = compute_res_perf( Qreg, tr*ones(T,1), Qspill+Qreg, Q_max_down );

%% CASE 3:
% Now let's make the operating policy vary over time. Specifically, we
% define two seasons: the dry season (from April to August) when the main
% objective is to deliver the target release and the risk of floods is low,
% and the wet season when the risk of flood is higher and so besides
% delivering the target release we want to leave some storage space for
% flood control.

% In the dry season, we can use the same operating policy (with the same
% parameters) as in Case 2:
x1d = pi/ 180    ; % (radiant) slope of first linear piece
x2d = 0.9*s_max  ; % (m3) storage at which second linear piece starts
x3d = pi / 32    ; % (radiant) slope of second linear piece
op_param_d = [ x1d, x2d, x3d, tr*delta, s_max, delta ]  ;
% In the wet season, instead, we start releasing more than the target
% release at much lower storage values:
x1w = pi/ 16    ; % (radiant) 
x2w = 0.7*s_max ; % (m3) 
x3w = pi / 100   ; % (radiant) 
op_param_w = [ x1w, x2w, x3w, tr*delta, s_max, delta ]  ;

% Let's visualise the two piece-wise linear function:
for i=1:length(s_test); 
u_testw(i) = op_piecewise_linear( s_test(i), op_param_w ) ; 
u_testd(i) = op_piecewise_linear( s_test(i), op_param_d ) ; 
end
figure; plot(s_test,u_testd,'color',[102 204 0]/255); hold on; plot(s_test,u_testw,'color',[0 0 255]/255); legend('dry season','wet season')
xlabel(stor_label); ylabel(flow_label)

% Let's put together the two parameterisations over the simulation horizon:
Op_param = repmat(op_param_w,T,1) ;
idxd = 183:335     ;  % Apr-Aug
Op_param(idxd,1) = x1d ;
Op_param(idxd,2) = x2d ;
Op_param(idxd,3) = x3d ;

% Another way to visualise the time-varying operating policy:
for t=1:T
for i=1:length(s_test); U_test(t,i) = op_piecewise_linear( s_test(i), Op_param(t,:) ) ; end
end
figure; 
clrs = gray ; clrs = clrs(end:-1:1,:); colormap(clrs)
imagesc(U_test(:,end:-1:1)'); ylabel(stor_label); xlabel(time_label)
colorbar; 
set(gca, 'YTick',[1,length(s_test)],'YTickLabel',{num2str(s_test(end)),num2str(s_test(1))})
set(gca,'XTick',xtick_month,'XTickLabel',label_month)

% Run simulation:
operating_rule = 'op_piecewise_linear';
[ s, Qreg, Qspill, E, S ] = reservoir_simulation( I, e, env_min, ...
                             s0, s_min, s_max, operating_rule, Op_param, delta ) ;

% Plot results:
plot_res_timeseries( s, s_min, s_max, I, E, Qreg, Qspill, env_min, leg,flow_label,stor_label , time_label ) ;

[Rel_reg, Def_reg, Vul_reg, Rel_spill, Qmax_spill] = compute_res_perf( Qreg, tr*ones(T,1), Qspill+Qreg, Q_max_down );

%% References
%
% Dobson, Wagener, Pianosi (2019), An argument-driven classification and
% comparison of reservoir operation optimization methods, Advances in Water
% Resources, 128, 74-86. doi: 10.1016/j.advwatres.2019.04.012
% https://www.sciencedirect.com/science/article/pii/S0309170818307759
