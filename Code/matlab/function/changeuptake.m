function [auto_model, mixo_model, hetero_model] = changeuptake(model1, model2, model3)
%CHANGEUPTAKE Set C. reinhardtii nutrient caps from growth-constrained solutions.
%   Inputs are Auto/Mixo/Hetero COBRA structs with positive split uptake IDs.
%   Outputs copy these models with selected nutrient upper bounds replaced
%   by uptake fluxes at exactly 99% of optimal biomass. IDs and 0.99 are
%   case-study settings; these are solution fluxes, not uptake FVA extrema.

auto_model = model1;
mixo_model = model2;
hetero_model = model3;

% find the index for growth rate
bio_index1 = find(model1.c == 1);
bio_index2 = find(model2.c == 1);
bio_index3 = find(model3.c == 1);

% find the optimal biomass value 
opt1 = optimizeCbModel(model1, 'max');
opt2 = optimizeCbModel(model2, 'max');
opt3 = optimizeCbModel(model3, 'max');
bio1 = opt1.f;
bio2 = opt2.f;
bio3 = opt3.f;

% find index of dissolved nutrient
dis_nu = {'EX_pi_e_REV','EX_nh4_e_REV','EX_so4_e_REV', 'EX_fe2_e_REV', 'EX_mg2_e_REV','EX_na1_e_REV', 'EX_ac_e_REV'};
index1 = false(size(model1.rxns));
index2 = false(size(model2.rxns));
index3 = false(size(model3.rxns));

for k = 1:numel(dis_nu)
    index1 = index1 | contains(model1.rxns, dis_nu{k});
    index2 = index2 | contains(model2.rxns, dis_nu{k});
    index3 = index3 | contains(model3.rxns, dis_nu{k});
end

nu_index1 = find(index1);
nu_index2 = find(index2);
nu_index3 = find(index3);

% calculate the uptake flux value of optimal growth rate
% max flux
% model1
model1.lb(bio_index1) = bio1*0.99;
model1.ub(bio_index1) = bio1*0.99;
opt1 = optimizeCbModel(model1);
max_nu1 = opt1.v(nu_index1);

auto_model.ub(nu_index1) = max_nu1;

% model2
model2.lb(bio_index2) = bio2*0.99;
model2.ub(bio_index2) = bio2*0.99;
opt2 = optimizeCbModel(model2);
max_nu2 = opt2.v(nu_index2);

mixo_model.ub(nu_index2) = max_nu2;

% model3
model3.lb(bio_index3) = bio3*0.99;
model3.ub(bio_index3) = bio3*0.99;
opt3 = optimizeCbModel(model3);
max_nu3 = opt3.v(nu_index3);

hetero_model.ub(nu_index3) = max_nu3;
