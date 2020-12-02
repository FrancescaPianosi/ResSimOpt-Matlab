% This script shows an example of how to run an optimisation of the
% reservoir system, using the 'reservoir_simulation' function and a
% multi-objective evolutionary algorithm (NSGA-II). 
% (for a discussion about using multi-objective evol. algorithms in the
% context of reservoir operations optimisation, see Pianosi et al (2011))
% It builds on the previous workflow on reservoir simulation, so refer to
% that workflow for more information about the rationale behind the
% definition of the simulation scenario and reservoir characteristic
% parameters. 
% This workflow uses the NSGA-II algorithm as implemented in a set of
% functions largely based on the code by Aravind Seshadri (2009),
% whose original version is available at:
% www.mathworks.com/matlabcentral/fileexchange/10429-nsga-ii-a-multi-objective-optimization-algorithm/

%% Load data and define system characteristics

% Load time series of reservoir inflows and
% evaporation rate (per unit surface area):
load Data
T = 365 ; % length of time series (days)
I = Data(1:T,1) ; % (m3/s)
e = Data(1:T,2) ; % (m/s)

% Define reservoir characteristics:
s_min = 1*10^7    ; % (m3)
s_max = 1*10^8    ; % (m3)
s0    = 0.9*s_max ; % (m3)
delta = 60*60*24  ; % (sec/day)
env_min        = 0.5*ones(T,1) ; % (m3/s) 
Qtarget        = 7*ones(T,1)   ; % (m3/s) the target release
Q_max_down     = 30 ; % (m3/s)
operating_rule = 'op_piecewise_linear';

%% Choose objective functions
% In this workflow we will use the function 'compute_res_perf' to calculate
% the reservoir performance over the simulation/optimisation period. We can
% then choose two or more of those performance indicators as objectives to
% be minimised by the optimiser.

% Define labels of the performance indicators returned by
% 'compute_res_perf' function (this will be useful later on to costumise
% plots):
res_perf_label = {'Freq. Failure Qtarget', 'Vol. Failure Qtarget', 'Vol.^2 Failure Qtarget', 'Freq. Failure Qdown', 'Max Qdown'};

% Choose which indicators will be used as objective functions:
res_perf_id = [3,4];

%% Define decision variables to be optimised
% In this workflow we will use an optimisation algorithm to find the
% combinations of parameters of the operating policy that optimise the
% objective functions. 
% As operating policy, we will use the piecewise linear function
% 'op_piecewise_linear'. This has 6 parameters, however 3 of them
% correspond to characteristics of the reservoir system and/or simulation
% model (the target release; the maximum storage; and the length of the
% simulation time-step) that cannot be varied by the optimiser. So the
% parameters that can be varied (the decision variables) are only three:
% x(1): the slope of the first linear piece (radiant)
% x(2): the storage at which second linear piece starts (volume)
% x(3): the slope of the second linear piece (radiant)
% If we allow these parameters to vary over time, for example according to
% the season, then the set of decision variables to be optimised comprises
% G*3 parameters, where G is the number of seasons (or sub-periods of the
% year) where we use distinct parameterisations.
% In this example we will use two seasons (G=2), wet and dry.

% Create time series to distinguish wet (1) and dry (2) season
idx  = ones(T,1) ;
idx(183:335) = 2 ; % Apr-Aug

% Define min and max values for the parameters to be optimised: 
x_min  = [    eps;     0;    eps;     eps;     0;    eps ] ;
x_max  = [ pi/2.1; s_max; pi/2.1;  pi/2.1; s_max; pi/2.1 ] ;

% Create a 'reasonable' parameter set for the operating policy 
% that will be passed on to the optimiser as algorithm 'initialisation' 
% (that is, these values will be included in the population evaluated
% at the first iteration; this can significantly speed up
% the convergence of the optimisation, see Pianosi et al 2011).
% A possible set of parameters in the dry season:
x1d = pi/ 180    ; % (radiant) slope of first linear piece
x2d = 0.9*s_max  ; % (m3) storage at which second linear piece starts
x3d = pi / 32    ; % (radiant) slope of second linear piece
% ... and in the wet season:
x1w = pi/ 16     ; % (radiant) 
x2w = 0.7*s_max  ; % (m3) 
x3w = pi / 100   ; % (radiant) 
% Let's put together the two parameterisations:
x_ini  = [    x1w;   x2w;   x3w ;     x1d;   x2d;    x3d ] ;

%% Run optimisation

% Add path to the folder including the NSGA-II functions:
addpath([ pwd '/opt'])
% Save all input data and parameters in a data struct that will be passed
% over to the optimiser:
sys_param.e       = e ;
sys_param.I       = I ;
sys_param.env_min = env_min ;
sys_param.Qtarget = Qtarget ;
sys_param.idx     = idx     ;
sys_param.s0      = s0      ; 
sys_param.s_min   = s_min   ; 
sys_param.s_max   = s_max   ; 
sys_param.delta   = delta   ;
sys_param.Q_max_down = Q_max_down ;
sys_param.res_perf_id    = res_perf_id  ;
sys_param.operating_rule = 'op_piecewise_linear';
sys_param.option = 0 ; % this is a parameter to activate a specific 

% Choose tuning parameters of NSGA-II and run the optimisation:
pop = 20 ; % nmber of decision vectors to be evaluated at each iteration
gen = 5  ; % number of iterations
M   = length(res_perf_id) ; % number of objectives
V   = length(x_ini)       ; % number of decision variables
[ chromosome_0, chromosome, f_intermediate ] = ...
                   nsga_2(pop,gen,M,V,x_min,x_max,x_ini,sys_param);

% plot initial and final population: 
x0 = chromosome_0(:,1:V) ;
xf = chromosome(:,1:V)   ;
figure
for i=1:V ; 
    subplot(1,V,i); plot(x0(:,i),'xk'); hold on; plot(xf(:,i),'or'); plot(x_ini(i),'xb') ; 
    axis([0,pop+1,x_min(i),x_max(i)]); xlabel('population')
end

% plot performances over generations: 
figure
for k=1:gen
    tmp  = reshape(f_intermediate(k,:)',M,pop)'; %(pop,M)
    for j=1:M
    subplot(1,M,j); hold on; plot(k,tmp(:,j),'or'); set(gca,'XLim',[0,gen+1]); box on
    xlabel('generation'); ylabel(res_perf_label{res_perf_id(j)})
    end
end
for j=1:M
    subplot(1,M,j); hold on; plot(0.5,chromosome_0(:,V+j),'xk'); set(gca,'XLim',[0,gen+1]); box on
end

% plot performances of initial and final population: 
f0    = chromosome_0(:,V+1:V+M) ;
ff    = chromosome(:,V+1:V+M)   ;
f_ini = chromosome_0(end,V+1:V+M) ;
if M==2
    figure;
    plot(f0(:,1),f0(:,2),'xk',ff(:,1),ff(:,2),'or',f_ini(end,1),f_ini(end,2),'xb','MarkerSize',10)
    xlabel(res_perf_label{res_perf_id(1)})
    ylabel(res_perf_label{res_perf_id(2)})
elseif M<4
    figure;
    plot_Pareto_Front(f0, {res_perf_label{res_perf_id}},1)
    title('before opt.')
    figure;
    plot_Pareto_Front(ff, {res_perf_label{res_perf_id}},1)
    title('after opt.')
end

%% Choose one optimised operating policy, plot and simulate it

flow_label =  'flow (m3/s)' ; stor_label = 'storage (m3)' ; time_label = 'time (days)';

idx_sol = 5 ;
x_cho = chromosome(idx_sol,1:V)   ;
Op_param_cho = op_piecewise_linear_transform(x_cho,idx,s_max,Qtarget,delta) ;

% visualise the operating policy:
%figure; for i=1:6 ; subplot(1,6,i); plot(Op_param_cho(:,i),'xk'); end
% Let's visualise the two piece-wise linear function:
s_test = [0:.25*10^6:s_max];
for i=1:length(s_test); 
u_testw(i) = op_piecewise_linear( s_test(i), [ x_cho(1:3),delta*Qtarget(1),s_max,delta] ) ; 
u_testd(i) = op_piecewise_linear( s_test(i), [ x_cho(4:6),delta*Qtarget(1),s_max,delta] ) ; 
end
figure; plot(s_test,u_testd,'color',[102 204 0]/255); hold on; plot(s_test,u_testw,'color',[0 0 255]/255); legend('dry season','wet season')
xlabel(stor_label); ylabel(flow_label)

% Simulate the operating policy 
% (this simulation was performed already by the optimiser in order to
% calculate the objective values, but as the simulated time series were not
% saved, we must run it again): 
[ s, Qreg, Qspill, E, S ] = reservoir_simulation( I, e, env_min, s0, s_min, s_max, operating_rule, Op_param_cho, delta ) ;
[Rel_reg, Def_reg, Vul_reg, Rel_spill, Qmax_spill] = compute_res_perf( Qreg, Qtarget, Qspill+Qreg, Q_max_down );
% Plot results:
leg = {'Inflow','Evap','Release','Spills','Min Env Flow'} ;
plot_res_timeseries( s, s_min, s_max, I, E, Qreg, Qspill, env_min, leg,flow_label,stor_label , time_label ) ;

%% References

% Pianosi, Quach Thi, Soncini-Sessa, 2011, Artificial Neural Networks and 
% Multi Objective Genetic Algorithms for water resources management: an
% application to the Hoabinh reservoir, IFAC Proceedings Volumes, 44(1),
% 10579-10584. doi:10.3182/20110828-6-IT-1002.02208.
% http://www.sciencedirect.com/science/article/pii/S1474667016453127


