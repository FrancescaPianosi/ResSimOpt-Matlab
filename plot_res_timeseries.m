function [] = plot_res_timeseries( s, s_min, s_max, I, E, Qreg, Qspill, env_min, leg,flow_label,stor_label , time_label )
%
% plot_res_timeseries( s, s_min, s_max, I, E, Qreg, Qspill, env_min, leg,flow_label,stor_label , time_label )
%
% s       = time series of storage                           - vector (N,1)
% s_min   = minimum storage (under which releases are reduced) - scalar
% s_max   = maximum storage (at which spills occur)            - scalar
% I       = time series of reservoir inflows                 - vector (N,1)
% E       = time series of evaporation                       - vector (N,1)
% Qreg    = time series of regulated release                 - vector (N,1)
% Qspill  = time series of release through spillways         - vector (N,1)
% env_min = time series of environmental flows               - vector (N,1)
% leg     = legend for the flow subplot                  - array of strings
%           for example: {'I','E','Qreg','Qspill','min env'}
% flow_label = label for vertical axis in flow subplot       - string
% stor_label = label for vertical axis in storage subplot    - string
% time_label = label for horizontal axis                     - string

T = length(I) ;

figure

% plot inflow/outflow fluxes:
subplot(211); hold on; box on
plot(I      ,'color',[215,25,28]/255)
plot(E      ,'color',[253,174,97]/255)
plot(Qreg   ,'color',[44,123,182]/255)
plot(Qspill ,'color',[171,217,233]/255)
plot(env_min,'color',[94,60,153]/255)

legend(leg)
ylabel(flow_label)
set(gca,'XLim',[1,T])

% plot storage:
subplot(212); hold on; box on
plot(s,'k')
plot(s_min*ones(T,1),':k'); plot(s_max*ones(T,1),':k')
set(gca,'XLim',[1,T])
ylabel(stor_label)
xlabel(time_label)

