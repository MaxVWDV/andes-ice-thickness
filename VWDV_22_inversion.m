function [ice_thickness] = VWDV_22_inversion(velocity,slope,inversion_parameters)
%This function inverts for ice thickness based on ice surface velocity
%measurements. Uncertainties in ice motion are accounted for
%probabilistically.
%
%Please see the following publications / submssions for additional details:
%
%Van Wyk de Vries, M., & Wickert, A. D. (2021). Glacier Image Velocimetry:
% an open-source toolbox for easy and rapid calculation of high-resolution
% glacier velocity fields. The Cryosphere, 15(4), 2115-2132.
%
% (Glacier thickness and ice volume of the inner tropical Andes) -
% submitted to Nature Scientific Data.
%
%Max Van Wyk de Vries and David Carchipulla, March 2021


%% Sample distributions of each parameter needed

% Slope
slope_sampled = normrnd(slope,inversion_parameters.slope_uncertainties); %Uncertainty from DEM sources
slope_sampled(slope_sampled<inversion_parameters.slope_minimum_threhsold)=inversion_parameters.slope_minimum_threhsold; %Do not allow slopes lower than minimum threshold

%Velocity
velocity_sampled = normrnd(velocity,inversion_parameters.ice_velocity_uncertainties); %Uncertainty from our own velocity calculations
velocity_sampled(velocity_sampled<0)=0; %Velocities must be positive

%Glen's flow law parameters
glens_flow_law_n_sampled = unifrnd(inversion_parameters.glen_flow_law_n.min,inversion_parameters.glen_flow_law_n.max);

temperature_sampled = unifrnd(inversion_parameters.temperature.min,inversion_parameters.temperature.max); %Note: this is used to calculate Glen's flow A

%Lateral drag correction ('shape factor') - Use info from Gantayat et al 2014
lateral_drag_correction_sampled = unifrnd(inversion_parameters.lateral_drag_correction.min,inversion_parameters.lateral_drag_correction.max);

%Basal sliding correction - E.g. Farinotti et al., 2009
basal_sliding_factor_sampled = unifrnd(inversion_parameters.basal_sliding_factor_b.min,inversion_parameters.basal_sliding_factor_b.max);

%Ice density
ice_density_sampled =...
    unifrnd(inversion_parameters.rho_ice_density.min,inversion_parameters.rho_ice_density.max);%Thouret et al - calculate density for different ice thickness columns; Tamayo and Arias, 2010; Cuffey and Paterson, 2010.


%% Now, iterate a few times with a small number of velocity maps to get a smoothing around 5x max ice thickness
%The idea here is to avoid a-priori assumptions

mean_ice_thick = inversion_parameters.initial_mean_ice_thickness; %Initialize

for loop_outer = 1:inversion_parameters.number_of_thickness_iterations
    %Filter the slope to ~5 times estimated mean ice thickness
    smooth_scale = double(round(5*mean_ice_thick/10));
    
    slope_sampled_sm = nanconv(slope_sampled,fspecial('disk',smooth_scale));    

    ice_thickness=real((((1-basal_sliding_factor_sampled)*velocity_sampled/(60*60*24*365))./...
    (((2*(inversion_parameters.glen_flow_law_A_zero*exp(-(115000/8.31)*((1/(temperature_sampled+273.15))-(1/273.15))))...
    /(glens_flow_law_n_sampled+1))...
    *(lateral_drag_correction_sampled*ice_density_sampled*inversion_parameters.gravity_g).^glens_flow_law_n_sampled)...
    *sin(deg2rad(slope_sampled_sm)).^glens_flow_law_n_sampled)).^(1/(glens_flow_law_n_sampled+1)));

    ice_thickness(ice_thickness>inversion_parameters.threshold_maximum_ice_thickness)=NaN;
    
%     fract(loop_outer)= mean_ice_thick/ nanmean(ice_thickness,'all');
    mean_ice_thick = nanmean(ice_thickness,'all');
    
end

ice_thickness(isnan(slope))=NaN;

