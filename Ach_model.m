function [eta1, eta2, scalingfactor, err] = Ach_model(trialparam, fr)
%ACH_MODELFUN   TDRL model fitting on cholinergic spiking data.
%   [ETA, S, ERR] = ACH_MODEL(TRIALPARAM,FR) fits TDRL 'surprise' model on
%   cholinergic spiking data (FR). Trial contingency information is passed
%   in TRIALPARAM. The model is defined by ACH_RL_MODELFUN. The best fit
%   parameters (ETA, S) and the MLE error (ERR) are returned.
%
%   See also BEST_FIT_CHOLINERGIC and ACH_RL_MODELFUN.

%   Panna Hegedus
%   Insititute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   15-Apr-2020

%   Code review: BH 4/15/20

% Initial parameters - eta and S
initparam = [0.5 0.5 30];

% Call 'fminsearchbnd' for optimization
modelfun = @Ach_RL_modelfun;   % goodness-of-fit function
[param, err] = fminsearchbnd(@(p)FitAchModel(p,fr,modelfun,trialparam),initparam,...
    [0, 0, -Inf],[1, 1, Inf],optimset('MaxFunEvals',500));

% Best fit parameters
eta1 = param(1);
eta2 = param(2);
scalingfactor = param(3);