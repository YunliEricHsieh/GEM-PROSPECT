function [CO2model, original_model] = create_CO2_model(model, alpha)
%CREATE_CO2_MODEL Build the manuscript's C. reinhardtii CO2-enriched variant.
%   model: COBRA struct with EX_co2_e_REV positive uptake and one c==1 biomass.
%   alpha: growth fraction for the lower reference solution. Outputs are a
%   CO2 variant and an updated original model. The variant upper bound is
%   max_co2 + 0.03*(max_co2-min_co2); the original cap is max_co2. Reference
%   fluxes come from biomass-constrained solutions, not separate uptake FVA.

CO2model = model;
original_model = model;

% find the index for growth rate
bio_index = find(model.c == 1);

% find the optimal biomass value 
opt = optimizeCbModel(model, 'max');
opt_bio = opt.f;

% find index of uptake CO2 rxn
co2_index = find(ismember(model.rxns,'EX_co2_e_REV'));

% calculate the uptake flux value of optimal growth rate
% max flux 
model.lb(bio_index) = opt_bio;
model.ub(bio_index) = opt_bio;
opt1 = optimizeCbModel(model);
max_co2 = opt1.v(co2_index);

% min flux
model.lb(bio_index) = opt_bio*alpha;
model.ub(bio_index) = opt_bio*alpha;
opt2 = optimizeCbModel(model);
min_co2 = opt2.v(co2_index);

% fix the max CO2 for original model
original_model.ub(co2_index) = max_co2;

% increase 3% of CO2 for CO2 model
CO2model.ub(co2_index) = max_co2 + (max_co2-min_co2)*0.03;
