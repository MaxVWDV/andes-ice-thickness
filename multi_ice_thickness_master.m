%This code file is a the master file for multiple ice thickness
%measurements
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

%% Step 0 - set path
name = 'Cotopaxi';




addpath(genpath('./sub_functions'))
path_to_file = strcat('./velocity_grids/',name,'/');



%% Step 1: Initialize parameters

%Number of Monte-Carlo iterations
inversion_parameters.number_of_iterations = 10;

%Buffer for drainage basin delineation
inversion_parameters.drainage_basin_buffer = 1000/10; % m / resolution = pixels
inversion_parameters.apparent_mass_balance_threshold = 50;
inversion_parameters.ELA_resolution = 10;
inversion_parameters.number_of_thickness_iterations = 5;
inversion_parameters.initial_mean_ice_thickness = 100;
inversion_parameters.threshold_maximum_ice_thickness = 300;

inversion_parameters.slope_uncertainties = 5;
inversion_parameters.slope_minimum_threhsold = 2;

inversion_parameters.gravity_g = 9.79;

%Glen's flow law parameters
inversion_parameters.glen_flow_law_n.min=3;
inversion_parameters.glen_flow_law_n.max=3;

inversion_parameters.glen_flow_law_A_zero = 2.4e-24;

inversion_parameters.temperature.min = -5;
inversion_parameters.temperature.max = -1;

inversion_parameters.rho_ice_density.min = 743; %Mean ice density of 830, Thouret et al - calculate density for different ice thickness columns; Tamayo and Arias, 2010; Cuffey and Paterson, 2010.
inversion_parameters.rho_ice_density.max = 917;

inversion_parameters.lateral_drag_correction.min = 0.8;
inversion_parameters.lateral_drag_correction.max = 1;

inversion_parameters.basal_sliding_factor_b.min = 0;
inversion_parameters.basal_sliding_factor_b.max = 0.4;

%Melt factors for accumulation and ablation zone (Farinotti 09)
inversion_parameters.mass_balance_accumulation.min = 0.4e-2;
inversion_parameters.mass_balance_accumulation.max = 0.6e-2;

inversion_parameters.mass_balance_ablation.min = 0.8e-2;
inversion_parameters.mass_balance_ablation.max = 1e-2;


%% Step 2: Load the ice, topography, and velocity
%Load full dem
dem_full = GRIDobj(strcat(path_to_file,'demdata/dem.tif'));
dem_full.cellsize=10;

%Load glacier mask
glacier_mask = GRIDobj(strcat(path_to_file,'demdata/glacier_mask.tif'));
glacier_mask.Z = flipud (glacier_mask.Z); %mask needs flipping to correct orientation

%Load x and y velocity components
velocity_x_component = GRIDobj(strcat(path_to_file,'Georeferenced Velocity Data/Mean, Standard Deviation and other statistics/Mean Velocity x component.tif'));
velocity_y_component = GRIDobj(strcat(path_to_file,'Georeferenced Velocity Data/Mean, Standard Deviation and other statistics/Mean Velocity y component.tif'));

%Set all to the same scale as the highest resolution endmember (mask)
reference_size = glacier_mask.size;

dem_full.Z = interp2(dem_full.Z,linspace(1,dem_full.size(2),reference_size(2))',linspace(1,dem_full.size(1),reference_size(1)));

velocity_x_component.Z = interp2(velocity_x_component.Z,linspace(1,velocity_x_component.size(2),reference_size(2))',...
    linspace(1,velocity_x_component.size(1),reference_size(1)));

velocity_y_component.Z = interp2(velocity_y_component.Z,linspace(1,velocity_y_component.size(2),reference_size(2))',...
    linspace(1,velocity_y_component.size(1),reference_size(1)));

dem_glacier = dem_full;

%Buffer mask for basins
glacier_mask_buffered=glacier_mask;
glacier_mask_buffered.Z = bwmorph(bwmorph(glacier_mask.Z==1,'thicken',inversion_parameters.drainage_basin_buffer),'close');

%Clip all datasets to the 'active' area (where ice is present)
dem_glacier_buffered = dem_full;
% dem_glacier_buffered.Z(glacier_mask_buffered.Z~=1)=NaN;
dem_glacier.Z(glacier_mask.Z~=1)=NaN;

% glacier_mask.Z = glacier_mask.Z(all(~isnan(nanmedian(glacier_mask.Z,2)),2),:); %Clip NaN in direction 1
% glacier_mask.Z = glacier_mask.Z(:,all(~isnan(nanmedian(glacier_mask.Z,1)),1)); %Clip NaN in direction 2
glacier_mask.size = size(glacier_mask.Z);
% 
% dem_glacier_buffered.Z = dem_glacier_buffered.Z(all(~isnan(nanmedian(glacier_mask_buffered.Z,2)),2),:); %Clip NaN in direction 1
% dem_glacier_buffered.Z = dem_glacier_buffered.Z(:,all(~isnan(nanmedian(glacier_mask_buffered.Z,1)),1)); %Clip NaN in direction 2
dem_glacier_buffered.size = size(dem_glacier_buffered.Z);
% 
% dem_glacier.Z = dem_glacier.Z(all(~isnan(nanmedian(glacier_mask.Z,2)),2),:); %Clip NaN in direction 1
% dem_glacier.Z = dem_glacier.Z(:,all(~isnan(nanmedian(glacier_mask.Z,1)),1)); %Clip NaN in direction 2
dem_glacier.size = size(dem_glacier.Z);

%Create slope
dem_slope = gradient8(dem_glacier,'deg');

%% Step 3: Sub-divide glacier basins into individual drainages
%Note, this solves equation 5 of Farinotti et al, (2009)

%Fill sinks in DEM
dem_filled_sinks = fillsinks(dem_glacier_buffered);

%Calculate a flow direction map
dem_flow_direction = FLOWobj(dem_filled_sinks);

%Divide this flow direction map into drainage basins
dem_drainage_basins = drainagebasins(dem_flow_direction);

%Clip to non-buffered mask
dem_drainage_basins.Z(glacier_mask.Z~=1)=NaN;
dem_filled_sinks.Z(glacier_mask.Z~=1)=NaN;
dem_slope.Z(glacier_mask.Z~=1)=NaN;

dem_drainage_basins.Z = dem_drainage_basins.Z(all(~isnan(nanmedian(glacier_mask.Z,2)),2),:); %Clip NaN in direction 1
dem_drainage_basins.Z = dem_drainage_basins.Z(:,all(~isnan(nanmedian(glacier_mask.Z,1)),1)); %Clip NaN in direction 2
dem_drainage_basins.size = size(dem_drainage_basins.Z);

dem_filled_sinks.Z = dem_filled_sinks.Z(all(~isnan(nanmedian(glacier_mask.Z,2)),2),:); %Clip NaN in direction 1
dem_filled_sinks.Z = dem_filled_sinks.Z(:,all(~isnan(nanmedian(glacier_mask.Z,1)),1)); %Clip NaN in direction 2
dem_filled_sinks.size = size(dem_filled_sinks.Z);

dem_slope.Z = dem_slope.Z(all(~isnan(nanmedian(glacier_mask.Z,2)),2),:); %Clip NaN in direction 1
dem_slope.Z = dem_slope.Z(:,all(~isnan(nanmedian(glacier_mask.Z,1)),1)); %Clip NaN in direction 2
dem_slope.size = size(dem_slope.Z);

dem_drainage_basins = shufflelabel(dem_drainage_basins);

drainage_stats.basins = dem_drainage_basins.Z;

drainage_stats.apparent_mass_balance = dem_filled_sinks.Z;

drainage_stats.upstream_area = dem_filled_sinks.Z;

drainage_stats.flux = dem_filled_sinks.Z*0;

drainage_stats.width = dem_filled_sinks.Z*0;

drainage_stats.dem = dem_filled_sinks.Z;

drainage_stats.slope = dem_slope.Z;

drainage_stats.number_of_basins = nanmax(drainage_stats.basins,[],'all');

drainage_stats.ELA = nan(1,drainage_stats.number_of_basins);

drainage_stats.elevation_range = [nanmin(dem_filled_sinks.Z,[],'all') nanmax(dem_filled_sinks.Z,[],'all')];

%% Step 4: Calculate the ice velocity uncertainties
%Calculate the ice speed
velocity = hypot(velocity_x_component.Z,velocity_y_component.Z);

try
    %Crop to non-ice areas to create uncertainties histogram
    velocity_x_component.Z(ice_mask==1)=NaN;%Crop to stable region

    velocity_y_component.Z(ice_mask==1)=NaN;

    %Create a histogram of x and y components
    [hist_counts_x,bins_x] = histcounts(velocity_x_component.Z,100); %Calculate the histogram counts
    [hist_counts_y,bins_y] = histcounts(velocity_y_component.Z,100); %Calculate the histogram counts

    gauss_data_x = fit(bins_x(1:end-1)',hist_counts_x','gauss1'); %Fit a single gaussian to the data
    gauss_data_y = fit(bins_y(1:end-1)',hist_counts_y','gauss1'); %Fit a single gaussian to the data

    uncertainties.velocity = (gauss_data_x.c1.^2+gauss_data_y.c1.^2).^(0.5); %The standard deviation of the sum of the two gaussians

catch
    uncertainties.velocity = 0.25; %If previous failed (e.g. no velocities, etc) assign 0.25 m/yr uncertainty (reasonable bound for large dataset)
end

if isnan(uncertainties.velocity)
    uncertainties.velocity = 0.25; %If previous failed (e.g. all NaN, etc) assign 0.25 m/yr uncertainty (reasonable bound for large dataset)
end

inversion_parameters.ice_velocity_uncertainties = uncertainties.velocity;

velocity(glacier_mask.Z~=1)=NaN;
velocity = velocity(all(~isnan(nanmedian(glacier_mask.Z,2)),2),:); %Clip NaN in direction 1
velocity = velocity(:,all(~isnan(nanmedian(glacier_mask.Z,1)),1)); %Clip NaN in direction 2


%% Step 5: Calculate an ice thickness using multiple methods

HF = NaN(size(velocity,1),size(velocity,2),inversion_parameters.number_of_iterations);
VWDV = NaN(size(velocity,1),size(velocity,2),inversion_parameters.number_of_iterations);
GT2b = NaN(size(velocity,1),size(velocity,2),inversion_parameters.number_of_iterations);
GT2w = NaN(size(velocity,1),size(velocity,2),inversion_parameters.number_of_iterations);
G14b = NaN(size(velocity,1),size(velocity,2),inversion_parameters.number_of_iterations);
G14w = NaN(size(velocity,1),size(velocity,2),inversion_parameters.number_of_iterations);

for parloop = 1:inversion_parameters.number_of_iterations
    [HF(:,:,parloop)] = HF12_inversion(drainage_stats,inversion_parameters);
    [VWDV(:,:,parloop)] = VWDV_22_inversion(velocity,dem_slope.Z,inversion_parameters);
    [GT2b(:,:,parloop)] = GlabTop2_basins_inversion(drainage_stats,inversion_parameters);
    [GT2w(:,:,parloop)] = GlabTop2_whole_inversion(drainage_stats,inversion_parameters);
    [G14b(:,:,parloop)] = Gantayat14_basins_inversion(drainage_stats,velocity,dem_slope.Z,inversion_parameters);
    [G14w(:,:,parloop)] = Gantayat14_whole_inversion(drainage_stats,velocity,dem_slope.Z,inversion_parameters);
    disp(strcat('Iteration_',num2str(parloop)));
end

save(strcat(name,'_thickness.mat'),"G14w","G14b","GT2w","GT2b","VWDV","HF");