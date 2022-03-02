function [ice_thickness] = GlabTop2_basins_inversion(drainage_stats,inversion_parameters)
%This code file is an implementation of the GlabTop2 ice
%thickness inversion. It also makes use of topographic analysis toolbox
%Topotolbox.
%
%Relevant citations to read are:
%
%RAAJ Ramsankaran, Ankur Pandit & Mohd Farooq Azam (2018)
% Spatially distributed ice-thickness modelling for Chhota Shigri Glacier in western
% Himalayas, India, International Journal of Remote Sensing, 39:10, 3320-3343, DOI:
% 10.1080/01431161.2018.1441563
%
%Haeberli, W., & Hoelzle, M. (1995). Application of inventory data for
% estimating characteristics of and regional climate-change effects on
% mountain glaciers: A pilot study with the European Alps. Annals of
% Glaciology, 21, 206-212. doi:10.3189/S0260305500015834
%
%Please see the publications listed above for full details, and comments
%throughout the code.

%% Step 1: Sample distributions of each parameter needed

% Slope
slope_sampled = normrnd(drainage_stats.slope,inversion_parameters.slope_uncertainties); %Uncertainty from DEM sources
slope_sampled(slope_sampled<inversion_parameters.slope_minimum_threhsold)=inversion_parameters.slope_minimum_threhsold; %Do not allow slopes lower than minimum threshold

%Glen's flow law parameters
glens_flow_law_n_sampled = unifrnd(inversion_parameters.glen_flow_law_n.min,inversion_parameters.glen_flow_law_n.max);

temperature_sampled = unifrnd(inversion_parameters.temperature.min,inversion_parameters.temperature.max); %Note: this is used to calculate Glen's flow A

%Lateral drag correction ('shape factor') - Use info from Gantayat et al 2014
lateral_drag_correction_sampled = unifrnd(inversion_parameters.lateral_drag_correction.min,inversion_parameters.lateral_drag_correction.max);

%Ice density
ice_density_sampled =...
    unifrnd(inversion_parameters.rho_ice_density.min,inversion_parameters.rho_ice_density.max);%Thouret et al - calculate density for different ice thickness columns; Tamayo and Arias, 2010; Cuffey and Paterson, 2010.

basal_shear_stress = NaN(size(drainage_stats.dem));



%% Step 2: Calculate basal shear stress (see Ramanskaran et al., 2018)
%Loop through all basins
for loop = 1:drainage_stats.number_of_basins

    if length(drainage_stats.dem(drainage_stats.basins==loop))>1

        elevation_difference = (max(drainage_stats.dem(drainage_stats.basins==loop),[],'all')-min(drainage_stats.dem(drainage_stats.basins==loop),[],'all'))/1000;

        if elevation_difference<1.6
            basal_shear_stress(drainage_stats.basins==loop) = (0.5+159.8*elevation_difference-43.5*elevation_difference.^2)*1000;
        else
            basal_shear_stress(drainage_stats.basins==loop) =150000;
        end



    end
end



%% Step 3: Calculate the ice thickness

mean_ice_thick = inversion_parameters.initial_mean_ice_thickness; %Initialize with an ice thickness of 100

for loop_outer = 1:inversion_parameters.number_of_thickness_iterations

    %Filter the slope to ~5 times estimated mean ice thickness
    smooth_scale = double(round(5*mean_ice_thick/10));

    slope_sampled_sm = nanconv(slope_sampled,fspecial('disk',smooth_scale));

    ice_thickness = real(basal_shear_stress./...
        (lateral_drag_correction_sampled.*ice_density_sampled.*inversion_parameters.gravity_g.*sin(deg2rad(slope_sampled_sm))));

    ice_thickness(ice_thickness>inversion_parameters.threshold_maximum_ice_thickness)=NaN;

    mean_ice_thick = nanmean(ice_thickness,'all');

end

