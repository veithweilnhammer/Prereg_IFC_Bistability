function [traj, infStates] = tapas_sutton_k1_binary(r, p, varargin)
% Calculates the trajectories of v under the Rescorla-Wagner learning model
%
% This function can be called in two ways:
% 
% (1) tapas_sutton_k1_binary(r, p)
%   
%     where r is the structure generated by tapas_fitModel and p is the parameter vector in native space;
%
% (2) tapas_sutton_k1_binary(r, ptrans, 'trans')
% 
%     where r is the structure generated by tapas_fitModel, ptrans is the parameter vector in
%     transformed space, and 'trans' is a flag indicating this.
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2013 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Transform paramaters back to their native space if needed
if ~isempty(varargin) && strcmp(varargin{1},'trans');
    p = tapas_sutton_k1_binary_transp(r, p);
end

% Unpack parameters
mu     = p(1);
Rhat   = p(2);
vhat_1 = p(3);
h_1    = p(4);

% Number of trials
u = r.u(:,1);
n = length(u);

% Initialize updated quantities
da    = NaN(n,1);
be    = NaN(n+1,1);
al    = NaN(n,1);
h     = NaN(n+1,1);
vhat  = NaN(n+1,1);

% Priors
vhat(1) = vhat_1;
be(1)   = log(Rhat);
h(1)    = h_1;

% Pass through value update loop
for k = 1:n
    if not(ismember(k, r.ign))
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Effect of input u(k)
        %%%%%%%%%%%%%%%%%%%%%%
        
        % Prediction error
        da(k) = u(k)-vhat(k);

        % Beta
        be(k+1) = be(k)+mu*da(k)*h(k);

        % Alpha
        al(k) =  exp(be(k+1))/(Rhat + exp(be(k+1)));
        
        % h
        h(k+1) = (h(k)+al(k)*da(k))*max((1-al(k)),0);

        % Prediction
        vhat(k+1) = vhat(k)+al(k)*da(k);
    else
        da(k)      = 0;
        be(k+1)    = be(k);
        al(k)      = al(k-1);
        h(k+1)     = h(k);
        vhat(k+1)  = vhat(k);
    end
end

% Posterior value
v = vhat;
v(1) = [];

% Remove ends of overlong trajectories
be(end)   = [];
h(end)    = [];
be(end)   = [];
vhat(end) = [];

% Create result data structure
traj = struct;

traj.da    = da;
traj.be    = be;
traj.al    = al;
traj.h     = h;
traj.v     = v;
traj.vhat  = vhat;

% Create matrix (in this case: vector) needed by observation model
infStates = traj.vhat;

return;
