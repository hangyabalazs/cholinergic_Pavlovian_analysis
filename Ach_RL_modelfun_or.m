function y = Ach_RL_modelfun_or(p, x)
%ACH_RL_MODELFUN_OR   TDRL model for salience prediction error.
%   Y = ACH_RL_MODELFUN_OR(P,TRIALPARAM) evaluates a TDRL model for
%   salience prediction error at paramaters P. The first parameter is a
%   discount factor for reward predection (ETA1) and the second parameter
%   is a discount factor for punishment prediction (ETA2). The third
%   parameters maps model outcomes onto spike numbers of firing rates (S).
%   Trial contingency information is passed in TRIALPARAM. The model is a
%   modified TDRL model. Temporal discounting is neglected due to short
%   delays. Salience prediction after cues is the sum of reward and
%   punishment prediction (the latter discounted by ETA). In this version,
%   omission responses were not discarded. Fitted model values for cues,
%   reward and punishment are returned (Y).
%
%   See also ACH_MODEL.

%   Panna Hegedus
%   Insititute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   15-Apr-2020

%   Code review: BH 4/15/20

% Parameters
eta1 = p(1);  % reward prediction discount factor
eta2 = p(2);  % punishment prediction discount factor
S = p(3);  % scaling factor between model outcomes and spike number/firing rate

% Trial contingency information
E_R = x(:, 1);  % reward pred.
E_P = x(:, 2);  % punish pred.
R = x(:, 3);  % reward indicator
P = x(:, 4);  % punishment indicator

% Model outcomes
inx = P==0 & R==0; % indices for cues (note: cues and reinforcements concatenated)
RPE = [R(~inx)-eta1*E_R(~inx)-eta2*E_P(~inx), zeros(length(R(~inx)),1)];  % reward prediction error
PPE = [P(~inx)-eta1*E_R(~inx)-eta2*E_P(~inx), zeros(length(P(~inx)),1)];  % punishment prediction error

y1 = S * (eta1*E_R(inx) + eta2*E_P(inx));  % model outcome for cues
y2 = S * (max(RPE') + max(PPE'));  % model outcome for reinforcements
y = [y1; y2'];   % concatenate cues and outcmes