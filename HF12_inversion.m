function [ice_thickness] = HF12_inversion(drainage_stats,inversion_parameters)
%This code file is an implementation of the Huss and Farinotti ice
%thickness inversion. It also makes use of topographic analysis toolbox
%Topotolbox.
%
%Relevant citations to read are:
%
%Farinotti, D., Huss, M., Bauder, A., & Funk, M. (2009). An estimate of the
% glacier ice volume in the Swiss Alps. Global and Planetary Change, 68(3),
%225-231.
%
%Huss, M., & Farinotti, D. (2012). Distributed ice thickness and volume of
% all glaciers around the globe. Journal of Geophysical Research: Earth
%Surface, 117(F4).
%
%Farinotti, D., Brinkerhoff, D. J., Clarke, G. K., Fürst, J. J., Frey, H.,
% Gantayat, P., ... & Andreassen, L. M. (2017). How accurate are estimates
% of glacier ice thickness? Results from ITMIX, the Ice Thickness Models
% Intercomparison eXperiment. The Cryosphere, 11(2), 949-970.
%
%Farinotti, D., Huss, M., Fürst, J. J., Landmann, J., Machguth, H.,
% Maussion, F., & Pandit, A. (2019). A consensus estimate for the ice
% thickness distribution of all glaciers on Earth. Nature Geoscience, 12(3)
%, 168-173.
%
%Schwanghart, W., & Scherler, D. (2014). TopoToolbox 2–MATLAB-based
% software for topographic analysis and modeling in Earth surface sciences.
% Earth Surface Dynamics, 2(1), 1-7.
%
%Please see the publications listed above for full details, and comments
%throughout the code.

%% Sample distributions of each parameter needed

% Slope
slope_sampled = normrnd(drainage_stats.slope,inversion_parameters.slope_uncertainties); %Uncertainty from DEM sources
slope_sampled(slope_sampled<inversion_parameters.slope_minimum_threhsold)=inversion_parameters.slope_minimum_threhsold; %Do not allow slopes lower than minimum threshold

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

%Apparent mass balance in the ablation zone - E.g. Farinotti et al., 2009
mass_balance_ablation_sampled = unifrnd(inversion_parameters.mass_balance_ablation.min,inversion_parameters.mass_balance_ablation.max);

%Apparent mass balance in the accumulation zone - E.g. Farinotti et al., 2009
mass_balance_accumulation_sampled = unifrnd(inversion_parameters.mass_balance_accumulation.min,inversion_parameters.mass_balance_accumulation.max);

%% Step 1: Calculate ELA for each sub-basin

%Loop through all basins
for loop = 1:drainage_stats.number_of_basins

    if length(drainage_stats.dem(drainage_stats.basins==loop))>1

        elevation_band = drainage_stats.dem(drainage_stats.basins==loop);

        %         elevation_brackets = ceil(min(elevation_band,[],'all')):inversion_parameters.ELA_resolution:floor(max(elevation_band,[],'all'));
        elevation_brackets = floor(min(elevation_band,[],'all')./inversion_parameters.ELA_resolution)*inversion_parameters.ELA_resolution...
            :inversion_parameters.ELA_resolution:ceil(max(elevation_band,[],'all')./inversion_parameters.ELA_resolution)*inversion_parameters.ELA_resolution;

        clear ELA_mismatch

        for loop_ELA = 1:length(elevation_brackets)

            ELA = elevation_brackets(loop_ELA);

            ELA_mismatch(loop_ELA) = abs(sum((ELA-elevation_band(elevation_band<ELA))*mass_balance_ablation_sampled,'all') - sum((elevation_band(elevation_band>=ELA)-ELA)*mass_balance_accumulation_sampled,'all'));

%             width_temp = max(sum((double(elevation_band>ELA)+double(elevation_band<ELA+inversion_parameters.ELA_resolution))==2),1)*10/tan(pi()/2-deg2rad(...
%                 nanmean(drainage_stats.slope((double(drainage_stats.basins==loop)+double(drainage_stats.dem>ELA)+double(drainage_stats.dem<=ELA+inversion_parameters.ELA_resolution))==3),'all')));

            drainage_stats.width((double(drainage_stats.basins==loop)+double(drainage_stats.dem>ELA)+double(drainage_stats.dem<=ELA+inversion_parameters.ELA_resolution))==3)=max(sum((double(elevation_band>ELA)+double(elevation_band<ELA+inversion_parameters.ELA_resolution))==2),1)*10/tan(pi()/2-deg2rad(...
                nanmean(drainage_stats.slope((double(drainage_stats.basins==loop)+double(drainage_stats.dem>ELA)+double(drainage_stats.dem<=ELA+inversion_parameters.ELA_resolution))==3),'all')));
        end

        [~,idx]=min(ELA_mismatch);

        drainage_stats.ELA(loop)=elevation_brackets(idx);

        %Calculate map-view apparent mass balance
        %Ablation area
        drainage_stats.apparent_mass_balance((double(drainage_stats.basins==loop)+double(drainage_stats.apparent_mass_balance<drainage_stats.ELA(loop)))==2)=...
            ((drainage_stats.apparent_mass_balance((double(drainage_stats.basins==loop)+double(drainage_stats.apparent_mass_balance<drainage_stats.ELA(loop)))==2))-drainage_stats.ELA(loop))*mass_balance_ablation_sampled;


        %Accumulation area
        drainage_stats.apparent_mass_balance((double(drainage_stats.basins==loop)+double(drainage_stats.apparent_mass_balance>=drainage_stats.ELA(loop)))==2)=...
            ((drainage_stats.apparent_mass_balance((double(drainage_stats.basins==loop)+double(drainage_stats.apparent_mass_balance>=drainage_stats.ELA(loop)))==2))-drainage_stats.ELA(loop))*mass_balance_accumulation_sampled;

    end
end

%Clip extreme mass balance endmembers (errors)
drainage_stats.apparent_mass_balance(abs(drainage_stats.apparent_mass_balance)>inversion_parameters.apparent_mass_balance_threshold)=NaN;

%% Step 4: Calculate the upsream flux q

%Loop through all basins
for loop = 1:drainage_stats.number_of_basins

    if length(drainage_stats.dem(drainage_stats.basins==loop))>1


        elevation_band = drainage_stats.dem(drainage_stats.basins==loop);

        %         elevation_brackets = ceil(min(elevation_band,[],'all')):inversion_parameters.ELA_resolution:floor(max(elevation_band,[],'all'));
        elevation_brackets = floor(min(elevation_band,[],'all')./inversion_parameters.ELA_resolution)*inversion_parameters.ELA_resolution...
            :inversion_parameters.ELA_resolution:ceil(max(elevation_band,[],'all')./inversion_parameters.ELA_resolution)*inversion_parameters.ELA_resolution;

        %     clear ELA_mismatch

        for loop_ELA = 1:length(elevation_brackets)

            ELA = elevation_brackets(loop_ELA);

            % Calculate contributing upsream area
%             flux_temp = sum(drainage_stats.apparent_mass_balance(...
%                 double(drainage_stats.basins==loop) + double(drainage_stats.dem>ELA)==2 ...
%                 ),'all')*10*10;

            drainage_stats.flux((double(drainage_stats.basins==loop)+double(drainage_stats.dem>ELA)+double(drainage_stats.dem<=ELA+inversion_parameters.ELA_resolution))==3)=sum(drainage_stats.apparent_mass_balance(...
                double(drainage_stats.basins==loop) + double(drainage_stats.dem>ELA)==2 ...
                ),'all')*10*10;

        end
    end
end



%% Step 5: Calculate the ice thickness

mean_ice_thick = inversion_parameters.initial_mean_ice_thickness; %Initialize with an ice thickness of 100

for loop_outer = 1:inversion_parameters.number_of_thickness_iterations
    %Filter the slope to ~5 times estimated mean ice thickness
    smooth_scale = double(round(5*mean_ice_thick/10));

    slope_sampled_sm = nanconv(slope_sampled,fspecial('disk',smooth_scale));

    ice_thickness = real((((1-basal_sliding_factor_sampled)*(drainage_stats.flux./(60*60*24*365.*drainage_stats.width))./ ...
        (2.*2*(inversion_parameters.glen_flow_law_A_zero*exp(-(115000/8.31)*((1/(temperature_sampled+273.15))-(1/273.15)))))).*...
        ((glens_flow_law_n_sampled+2)./...
        (lateral_drag_correction_sampled*ice_density_sampled.*inversion_parameters.gravity_g.*sin(deg2rad(slope_sampled_sm))).^glens_flow_law_n_sampled))...
        .^(1./(glens_flow_law_n_sampled+2)));

    ice_thickness(ice_thickness>inversion_parameters.threshold_maximum_ice_thickness)=NaN;

    mean_ice_thick = nanmean(ice_thickness,'all');

end

