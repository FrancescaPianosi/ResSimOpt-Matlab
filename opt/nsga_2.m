function [ chromosome_0, chromosome, f_intermediate ] = nsga_2(pop,gen,M,V,min_range,max_range,X_ini,sys_param)
%
% [ chromosome_0, chromosome, f_intermediate ] = ...
%                   nsga_2(pop,gen,M,V,min_range,max_range,X_ini,sys_param)
%
% is a multi-objective optimization function implementing the NSGA-II
% algorithm. The function is largely based on the code by Aravind Seshadri
% (2009) available at:
% www.mathworks.com/matlabcentral/fileexchange/10429-nsga-ii-a-multi-objective-optimization-algorithm/
% but with some modifications to make it suitable for reservoir operation
% optimisation. All changes made are marked by comments.
%
% Input arguments
% pop - population size
% gen - total number of generations
% M   - number of objectives 
% V   - number of decision variables
% min_range - minimum feasable value for the decision variables
% max_range - maximum feasable value for the decision variables
% x_ini   - set of decision variables to be included in the population
%             evaluated at first iteration
% sys_param - additional system parameters for simulation/optimisation
%
% This functions is based on evolutionary algorithm for finding the optimal
% solution for multiple objective i.e. pareto front for the objectives. 
% Initially enter only the population size and the stoping criteria or
% the total number of generations after which the algorithm will
% automatically stopped. 
%
% You will be asked to enter the number of objective functions, the number
% of decision variables and the range space for the decision variables.
% Also you will have to define your own objective funciton by editing the
% evaluate_objective() function. A sample objective function is described
% in evaluate_objective.m. Kindly make sure that the objective function
% which you define match the number of objectives that you have entered as
% well as the number of decision variables that you have entered. The
% decision variable space is continuous for this function, but the
% objective space may or may not be continuous.
%
% Original algorithm NSGA-II was developed by researchers in Kanpur Genetic
% Algorithm Labarotary and kindly visit their website for more information
% http://www.iitk.ac.in/kangal/
%
%  Copyright (c) 2009, Aravind Seshadri
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

% Simple error checking
% Number of Arguments
% Check for the number of arguments. The two input arguments are necessary
% to run this function.
if nargin < 2
    error('NSGA-II: Please enter the population size and number of generations as input arguments.');
end
% Both the input arguments need to of integer data type
if isnumeric(pop) == 0 || isnumeric(gen) == 0
    error('Both input arguments pop and gen should be integer datatype');
end
% Minimum population size has to be 20 individuals
if pop < 20
    error('Minimum population for running this function is 20');
end
if gen < 5
    error('Minimum number of generations is 5');
end
% Make sure pop and gen are integers
pop = round(pop);
gen = round(gen);

% infomation to file name

%%%%%%%%%% begin insertion by Francesca Pianosi
% choose the length of the simulation horizon:
N = length(sys_param.I) ;
if ~isfield(sys_param, 'option'); sys_param.option = 0; end
if sys_param.option == 0 % option 1 - always use the entire time series
    h_0 = N ;
else % option 2 - increase the length of the time series at each iteration,
    % starting from the length P specified by user
    P = sys_param.option    ;
    h_0 = min(N,P)          ; %(just in case the user mistakenly set P<N)
end
sys_param.nsga.e       = sys_param.e(1:h_0)       ;
sys_param.nsga.I       = sys_param.I(1:h_0)       ;
sys_param.nsga.env_min = sys_param.env_min(1:h_0) ;
sys_param.nsga.Qtarget = sys_param.Qtarget(1:h_0) ;
sys_param.nsga.idx     = sys_param.idx(1:h_0)     ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end

% Initialize the population
% Population is initialized with random values which are within the
% specified range. Each chromosome consists of the decision variables. Also
% the value of the objective functions, rank and crowding distance
% information is also added to the chromosome vector but only the elements
% of the vector which has the decision variables are operated upon to
% perform the genetic operations like corssover and mutation.


%%% begin modification by Francesca Pianosi (include 'sys_param' to inputs)
%chromosome = initialize_variables(pop, M, V, min_range,max_range);% old
chromosome   = initialize_variables(pop, M, V, min_range, max_range,X_ini, sys_param);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end

chromosome_0 = chromosome ; % save initialization

%
% Sort the initialized population
% Sort the population using non-domination-sort. This returns two columns
% for each individual which are the rank and the crowding distance
% corresponding to their position in the front they belong. At this stage
% the rank and the crowding distance for each chromosome is added to the
% chromosome vector for easy of computation.
chromosome = non_domination_sort_mod(chromosome, M, V);

% Start the evolution process
% The following are performed in each generation
% * Select the parents which are fit for reproduction
% * Perfrom crossover and Mutation operator on the selected parents
% * Perform Selection from the parents and the offsprings
% * Replace the unfit individuals with the fit individuals to maintain a
%   constant population size.
f_intermediate = ones(gen,M*pop) ;
for i = 1 : gen
    
    %%%%%%%%%% begin insertion by Francesca Pianosi
    % length of the simulation horizon to be used at
    % current iteration:
    if sys_param.option == 0 % option 1 - always use the entire time series
        h_k = N ;
    else % option 2 - increase the length of the time series at each iteration
        P = sys_param.option        ;
        h_k = ceil(i/(gen/(N/P)))*P ;
        h_k = min(h_k, S)           ;
    end
    % select disturbance trajectory from the
    % available data set:
    sys_param.nsga.e       = sys_param.e(1:h_k)       ;
    sys_param.nsga.I       = sys_param.I(1:h_k)       ;
    sys_param.nsga.env_min = sys_param.env_min(1:h_k) ;
    sys_param.nsga.Qtarget = sys_param.Qtarget(1:h_k) ;
    sys_param.nsga.idx     = sys_param.idx(1:h_k)     ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end
    
    % Select the parents
    % Parents are selected for reproduction to generate offspring. The
    % original NSGA-II uses a binary tournament selection based on the
    % crowded-comparision operator. The arguments are 
    % pool - size of the mating pool. It is common to have this to be half the
    %        population size.
    % tour - Tournament size. Original NSGA-II uses a binary tournament
    %        selection, but to see the effect of tournament size this is kept
    %        arbitary, to be choosen by the user.
    pool = round(pop/2);
    tour = 2;
    % Selection process
    % A binary tournament selection is employed in NSGA-II. In a binary
    % tournament selection process two individuals are selected at random
    % and their fitness is compared. The individual with better fitness is
    % selected as a parent. Tournament selection is carried out until the
    % pool size is filled. Basically a pool size is the number of parents
    % to be selected. The input arguments to the function
    % tournament_selection are chromosome, pool, tour. The function uses
    % only the information from last two elements in the chromosome vector.
    % The last element has the crowding distance information while the
    % penultimate element has the rank information. Selection is based on
    % rank and if individuals with same rank are encountered, crowding
    % distance is compared. A lower rank and higher crowding distance is
    % the selection criteria.
    parent_chromosome = tournament_selection(chromosome, pool, tour);

    % Perform crossover and Mutation operator
    % The original NSGA-II algorithm uses Simulated Binary Crossover (SBX) and
    % Polynomial  mutation. Crossover probability pc = 0.9 and mutation
    % probability is pm = 1/n, where n is the number of decision variables.
    % Both real-coded GA and binary-coded GA are implemented in the original
    % algorithm, while in this program only the real-coded GA is considered.
    % The distribution indeices for crossover and mutation operators as mu = 20
    % and mum = 20 respectively.
    mu = 20;
    mum = 20;
%%% begin modification by Francesca Pianosi (include 'sys_param' to inputs)
    offspring_chromosome = ...
        genetic_operator(parent_chromosome, ...
        M, V, mu, mum, min_range, max_range, sys_param);
%     offspring_chromosome = ...
%         genetic_operator(parent_chromosome, ...
%         M, V, mu, mum, min_range, max_range); % old %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end

    % Intermediate population
    % Intermediate population is the combined population of parents and
    % offsprings of the current generation. The population size is two
    % times the initial population.
    
    [main_pop,temp] = size(chromosome);
    [offspring_pop,temp] = size(offspring_chromosome);
    % temp is a dummy variable.
    clear temp
    % intermediate_chromosome is a concatenation of current population and
    % the offspring population.
    intermediate_chromosome(1:main_pop,:) = chromosome;
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = ...
        offspring_chromosome;

    % Non-domination-sort of intermediate population
    % The intermediate population is sorted again based on non-domination sort
    % before the replacement operator is performed on the intermediate
    % population.
    intermediate_chromosome = ...
        non_domination_sort_mod(intermediate_chromosome, M, V);
    % Perform Selection
    % Once the intermediate population is sorted only the best solution is
    % selected based on it rank and crowding distance. Each front is filled in
    % ascending order until the addition of population size is reached. The
    % last front is included in the population based on the individuals withchromosome_0
    % least crowding distance
    chromosome = replace_chromosome(intermediate_chromosome, M, V, pop);
    if ~mod(i,100)
        clc
        fprintf('%d generations completed\n',i);
    end
	f = chromosome(:,V+1:V+M)'; 
	f = f(:); 
	f_intermediate(i,:) = f';   
end

%% Result
% % Save the result in ASCII text format.

%eval(['save -ascii ini_', sys_param.filename, ' chromosome_0']) ;
%eval(['save -ascii sol_', sys_param.filename, ' chromosome'  ]) ;

%% Visualize
% % The following is used to visualize the result if objective space
% % dimension is visualizable.
% if M == 2
%     plot(chromosome_0(:,V+1),chromosome_0(:,V+2),'or') 
%     hold on
%     plot(chromosome(:,V+1),chromosome(:,V+2),'*b')
%     legend('before opt.','after opt.')
% elseif M ==3
%     plot3(chromosome_0(:,V+1),chromosome_0(:,V+2),chromosome_0(:,V+3),'or') 
%     hold on
%     plot3(chromosome(:,V+1),chromosome(:,V+2),chromosome(:,V+3),'*')
%     legend('before opt.','after opt.')
% end
