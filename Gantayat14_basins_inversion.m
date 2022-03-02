function [ice_thickness] = Gantayat14_basins_inversion(drainage_stats,velocity,slope,inversion_parameters)
%This function inverts for ice thickness based on ice surface velocity
%measurements. Uncertainties in ice motion are accounted for
%probabilistically. We used a slightly modified version  of Gantayat et
%al., 2014's method.
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
%Gantayat, P., Kulkarni, A. V., & Srinivasan, J. (2014). Estimation of ice
% thickness using surface velocities and slope: case study at Gangotri
% Glacier, India. Journal of Glaciology, 60(220), 277-282.
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

%% Step 2: Calculate the velocity and slope in 100m elevation bands

slope_sampled_bands = NaN(size(slope_sampled));
velocity_sampled_bands = NaN(size(slope_sampled));

%Loop through all basins
for loop = 1:drainage_stats.number_of_basins

    if length(drainage_stats.dem(drainage_stats.basins==loop))>1


        elevation_band = drainage_stats.dem(drainage_stats.basins==loop);

        elevation_brackets = floor(min(elevation_band,[],'all')./100)*100 ...
            :100:ceil(max(elevation_band,[],'all')./100)*100;

        %     clear ELA_mismatch

        for loop_ELA = 1:length(elevation_brackets)

            ELA = elevation_brackets(loop_ELA);

            % Calculate mean of this bracket
            slope_band = nanmean(slope_sampled(...
                (double(drainage_stats.basins==loop) +double(drainage_stats.dem>ELA)+double(drainage_stats.dem<=ELA+100))==3 ...
                ),'all');

            velocity_band = nanmean(velocity_sampled(...
                (double(drainage_stats.basins==loop) +double(drainage_stats.dem>ELA)+double(drainage_stats.dem<=ELA+100))==3 ...
                ),'all');

            slope_sampled_bands((double(drainage_stats.basins==loop)+double(drainage_stats.dem>ELA)+double(drainage_stats.dem<=ELA+100))==3)=slope_band;

            velocity_sampled_bands((double(drainage_stats.basins==loop)+double(drainage_stats.dem>ELA)+double(drainage_stats.dem<=ELA+100))==3)=velocity_band;

        end
    end
end


%% Now, iterate a few times with a small number of velocity maps to get a smoothing around 5x max ice thickness
%The idea here is to avoid a-priori assumptions

mean_ice_thick = inversion_parameters.initial_mean_ice_thickness; %Initialize

for loop_outer = 1:inversion_parameters.number_of_thickness_iterations
    %Filter the slope to ~5 times estimated mean ice thickness
    smooth_scale = double(round(5*mean_ice_thick/10));

    slope_sampled_sm = nanconv(slope_sampled_bands,fspecial('disk',smooth_scale));

    ice_thickness= real(((1.5*velocity_sampled_bands/(60*60*24*365))./ (...
        (inversion_parameters.glen_flow_law_A_zero*exp(-(115000/8.31)*((1/(temperature_sampled+273.15))-(1/273.15)))) .* ...
        (lateral_drag_correction_sampled*ice_density_sampled*inversion_parameters.gravity_g * ...
        sin(deg2rad(slope_sampled_sm))).^glens_flow_law_n_sampled))...
        .^(1/(glens_flow_law_n_sampled+1)));

    ice_thickness(ice_thickness>inversion_parameters.threshold_maximum_ice_thickness)=NaN;

    ice_thickness = nanconv(ice_thickness,fspecial('disk',9));  %Smooth with a filter. Gantayat et al. used 3x3, this 9x9 is similar with the oversampled resolution


    mean_ice_thick = nanmean(ice_thickness,'all');

end

ice_thickness(isnan(slope))=NaN;

