function [hypo10, hypo25, hypo75, original_model] = create_hypo_model(model, alpha)
%CREATE_HYPO_MODEL Build three manuscript-specific C. reinhardtii nutrient models.
%   model: Mixo COBRA struct with the listed positive nutrient uptake IDs.
%   alpha sets biomass for the lower reference solution. Outputs set uptake
%   upper bounds at 10/25/75% of the lower-to-upper solution-flux range and
%   return the original with upper reference caps. Equal reference fluxes
%   instead use lower flux times the fraction. No general media conversion.

hypo10 = model;
hypo25 = model;
hypo75 = model;
original_model = model;

% find the index for growth rate
bio_index = find(model.c == 1);

% find the optimal biomass value 
opt = optimizeCbModel(model, 'max');
opt_bio = opt.f;

% find index of dissolved nutrient
dis_nu = {'EX_pi_e_REV','EX_nh4_e_REV','EX_so4_e_REV', 'EX_fe2_e_REV', 'EX_mg2_e_REV','EX_na1_e_REV', 'EX_ac_e_REV'};

index = false(size(model.rxns));

for k = 1:numel(dis_nu)
    index = index | contains(model.rxns, dis_nu{k});
end

nu_index = find(index);

% calculate the uptake flux value of optimal growth rate
% max flux
model.lb(bio_index) = opt_bio;
model.ub(bio_index) = opt_bio;
opt1 = optimizeCbModel(model);
max_nu = opt1.v(nu_index);

% min flux
model.lb(bio_index) = opt_bio*alpha;
model.ub(bio_index) = opt_bio*alpha;
opt2 = optimizeCbModel(model);
min_nu = opt2.v(nu_index);

original_model.ub(nu_index) = max_nu;

% hypoosmotic TAP
for k = 1:numel(nu_index)
    if max_nu(k) == min_nu(k)
        hypo10.ub(nu_index(k)) = min_nu(k)*0.1;
        hypo25.ub(nu_index(k)) = min_nu(k)*0.25;
        hypo75.ub(nu_index(k)) = min_nu(k)*0.75;
    else
        hypo10.ub(nu_index(k)) = min_nu(k) + (max_nu(k)-min_nu(k))*0.1;
        hypo25.ub(nu_index(k)) = min_nu(k) + (max_nu(k)-min_nu(k))*0.25;
        hypo75.ub(nu_index(k)) = min_nu(k) + (max_nu(k)-min_nu(k))*0.75;
    end
end
