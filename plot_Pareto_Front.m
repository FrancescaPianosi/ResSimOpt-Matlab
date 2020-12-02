function [] = plot_Pareto_Front(M, labs, BW, font_size) 
%
% [] = plot_Pareto_Front(M, labs, BW, font_size) 
%
% Function to plot a 4D pareto front in 2D assigning the remaining
% objectives to the dimension of the circles and to the color (summer
% colormap). 
%
% input:
%   M     = matrix of objective values (to be minimized)
%           It can have from 2 to 4 columns.
%   labs  = cell array containing the objective labels
%           It must have the same number of elements as
%           the columns in M.
%   BW    = black and white option for the plot. 
%           If BW == 1, grey-scale is used, 
%           otherwise the summer colormap is used
%           [default: 0]
% font_size = font size of labels and axes ticks
%           [default: 12]

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

if(nargin<1) 
  error(  'too few arguments'  )
end

if(nargin>4)
  error(  'too many arguments'  )
end

if(nargin<4) % set default 'font_size' value
   font_size = 12 ;
end

if(nargin<3) % set default 'BW' value
    BW = 0 ; 
end

[N,q]=size(M); % check size of 'M'
if q ~= 3; error('''M'' must have at least 3 columns'); end

if(nargin>1) % check size of 'labs'
    q2 = length(labs) ;
    if q ~= q2; error('''labs'' must have as many elements as the number of columns in ''M'''); end
end

% Plot pareto front:

x = M(:,1) ;
y = M(:,2) ;
w = M(:,3) ;
if(BW == 1)
    w = -w ; % the sign of w is changed in order to plot the best alternatives in white and the worst in black
    w_min = min(w) ;
    w_max = max(w) ;
    ww = ( w - w_min )/( w_max - w_min ) ;
    mymap = repmat(ww,1,3);
    M = [x, y, w, mymap]  ;
    M = sortrows(M, 3)    ;
    mymap = M(:,end-2:end);
else
    M = [x, y, w] ;
    M = sortrows(M, 3) ;
    mymap = summer(length(w)) ;
end
labs3 = {} ; for i=1:N; labs3{i} = num2str(M(i,3)) ; end
colormap(mymap); lcolorbar(labs3,'TitleString', labs(3),'FontSize',font_size)
hold on
for i=1:N
    plot(M(i,1),M(i,2),'ko','MarkerSize',10,'MarkerFaceColor',mymap(i,:))
end
xlabel(labs{1},'FontSize',font_size) ; ylabel(labs{2},'FontSize',font_size)
grid on ; box on
axis([min(x)-std(x)/2,max(x)+std(x)/2,min(y)-std(y)/2,max(y)+std(y)/2])
set(gca,'FontSize',font_size)
