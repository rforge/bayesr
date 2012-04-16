%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* campus *.* lmu *.* de]
%% Time-stamp: <[example.prg] by DSB Die 10/03/2009 10:55 (GMT) on daniel@puc-home>
%%
%% Description:
%% Independent BayesX batch file for running the example model to create the
%% samples data. 
%%
%% History:
%% 10/03/2009   file creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open log file
logopen, replace using log.txt

% create bayesreg object
bayesreg bayesregObject

% specify basename
bayesregObject.outfile = res

% load the response and covariates data
dataset data
data.infile using data.txt

% load graph file for the districts in Tanzania
map tanzania
tanzania.infile using tanzania.bnd
tanzania.reorder

%%% MCMC run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model specification and mcmc options:
bayesregObject.regress y = x1(psplinerw2) + x2(rw1) + x3 + x4 + district(spatial, map=tanzania), family=gaussian iterations=10000 burnin=1000 step=10 predict using data
% Warning:
% Of course, these low numbers for iterations and burnin are just for the sake of illustration!! 

%%% tidy up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the samples
bayesregObject.getsample

% close log
logclose

% and quit
quit


% Local Variables:
% mode: latex
% coding: iso-latin-1-unix
% End:
