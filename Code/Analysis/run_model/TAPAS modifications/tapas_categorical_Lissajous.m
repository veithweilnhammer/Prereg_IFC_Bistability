function logp = tapas_rs_surprise(r, infStates, ptrans)


% Calculates the log-probability of response speed y (in units of ms^-1) according to the surprise
% model introduced in:
%
% Vossel, S.*, Mathys, C.*, Daunizeau, J., Bauer, M., Driver, J., Friston, K. J., and Stephan, K. E.
% (2013). Spatial Attention, Precision, and Bayesian Inference: A Study of Saccadic Response Speed.
% Cerebral Cortex.
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2013 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Transform zetas to their native space
ze = exp(ptrans(1));





% Initialize returned log-probabilities as NaNs so that NaN is
% returned for all irregualar trials
logp = NaN(length(infStates(:,1,1)),1);


% Weed irregular trials out from inferred states, responses, and inputs
mu1hat = infStates(:,1,1);
mu1hat(r.irr) = [];

x=infStates(:,1,5);
x(r.irr,:)=[];

y = r.y; % response speed
 % preceptual decision
y(r.irr,:) = [];

u = r.u;
u(r.irr,:) = [];





logx = log(x);
log1pxm1 = log1p(x-1);
logx(1-x<1e-4) = log1pxm1(1-x<1e-4);
log1mx = log(1-x);
log1pmx = log1p(-x);
log1mx(x<1e-4) = log1pmx(x<1e-4); 

% % Calculate log-probabilities for non-irregular trials
logp(not(ismember(1:length(logp),r.irr))) = y(:,1).*ze.*(logx -log1mx) +ze.*log1mx -log((1-x).^ze +x.^ze);




return;
