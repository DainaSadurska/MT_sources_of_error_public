%% SCRIPT 3 - Base frequency calculation

% STEP 1a: Load data
%%{
% Set input paths
paths = struct();
paths.root = 'path/to/main/directory';
pahs.output = fullfile(paths.root,"Output_3");

in_paths = struct();
in_dir1 = 'Output_1';
in_dir2 = 'Output_2';
%in_dir = '/Volumes/SSK 1TB/DPhil backup/2022/MATLAB/Data_20220516_1KG_reanalysis/Output_20220601';

in_paths.Basecalls_All = fullfile(paths.root,in_dir1,'Basecalls_ALL.mat');
in_paths.Basecalls_For = fullfile(paths.root,in_dir1,'Basecalls_FOR.mat');
in_paths.Basecalls_Rev = fullfile(paths.root,in_dir1,'Basecalls_REV.mat');
%in_paths.Indels       = fullfile(paths.root,in_dir1,'Indels.mat');
%in_paths.Reads_TOTAL  = fullfile(paths.root,in_dir1,'Reads_TOTAL.mat');
in_paths.SampleIDs     = fullfile(paths.root,in_dir1,"SampleIDs.mat");
in_paths.Individuals   = fullfile(paths.root,in_dir1,"NumberIDs.mat");
in_paths.Positions     = fullfile(paths.root,in_dir1,"Positions.mat");
in_paths.rCRS          = fullfile(paths.root,in_dir2,"rCRS.mat");
in_paths.Reference     = fullfile(paths.root,in_dir2,"G1K_HiCov_Consensus_sequences.mat");

% Create data struct
if ~exist('data', 'var'); data = struct(); end

% Load variables
data.Basecalls.All = load(in_paths.Basecalls_All, 'A', 'C', 'G', 'T');
data.Basecalls.For = load(in_paths.Basecalls_For, 'A', 'C', 'G', 'T');
data.Basecalls.Rev = load(in_paths.Basecalls_Rev, 'A', 'C', 'G', 'T');
%data.Indels        = load(in_paths.Indels, 'Ins', 'Del');
%data.Reads_TOTAL   = load(in_paths.Reads_TOTAL, 'Reads_TOTAL');
data.SampleIDs     = load(in_paths.SampleIDs, 'SampleIDs');
data.Individuals   = load(in_paths.Individuals, 'NumberIDs');
data.Positions     = load(in_paths.Positions, 'Positions');
data.rCRS          = load(in_paths.rCRS, 'rCRS');
data.Reference     = load(in_paths.Reference,'ConcensusSeqs');

% Clear input variables
clear -regexp ^in
%}


% STEP 1b: Set parameters
%%{
if ~exist('parameters', 'var'); parameters = struct(); end

% Whch dimension is individuals vs positions
parameters.dim_positions = 1;
parameters.dim_individuals = 2;

% Number of samples and positions
parameters.n_positions   = size(data.Positions,   parameters.dim_positions);
parameters.n_individuals = size(data.Individuals, parameters.dim_individuals);
parameters.dimensions = [parameters.n_positions parameters.n_individuals];

% rCRS deletion artefact site
parameters.deletion_site = 3107;

% Minimum required directional read depth
parameters.read_depth_min = 500;

% Setting for saving 
settings = struct();
settings.save1 = true;
settings.save2 = true;

%}

%% STEP 2: Calculate total basecall counts
%%{
% Read depths
for temp_dir = ["For", "Rev", "All"]
    % Calculate total number of basecalls made (excluding Ns)
    temp_reads  = data.Basecalls.(temp_dir).A ...
                + data.Basecalls.(temp_dir).C ...
                + data.Basecalls.(temp_dir).G ...
                + data.Basecalls.(temp_dir).T;
    data.Reads.(temp_dir) = temp_reads;

    % Create a mask for positions with sufficient read depth
    temp_mask = true(parameters.dimensions);
    switch temp_dir
        case "For" % Sufficieint forward depth positions
            temp_mask(temp_reads<parameters.read_depth_min) = false;

        case "Rev" % Sufficieint reverse depth positions
            temp_mask(temp_reads<parameters.read_depth_min) = false;

        case "All" % Sufficient in both directions (FALSE if either FOR or REV < MIN) 
            temp_mask(~data.Mask_read_depth.For) = false;
            temp_mask(~data.Mask_read_depth.Rev) = false;
    end
    data.Mask_read_depth.(temp_dir) = temp_mask;
end
% Clear temp
clear -regexp ^temp
%}

%% STEP 3: Calculate base frequencies
%%{
for temp_dir = ["For", "Rev", "All"]
    for temp_base = ["A", "C", "G", "T"]
        % Get base call counts
        temp_reads_base  = data.Basecalls.(temp_dir).(temp_base); % Calls of this base
        temp_reads_total = data.Read_depth.(temp_dir);            % All basecalls
        temp_reads_total(temp_reads_total==0) = NaN;              % replace 0 with NaN

        % Calculate base frequency (NaNs in places where there are no reads)
        data.Frequencies.(temp_dir).(temp_base) = temp_reads_base ./ temp_reads_total;
    end
end
% Clear temp 
clear -regexp ^temp
%}

%% STEP 4: Save basecall, read depth and base frequency structs to MAT file
%%{
if settings.save == true
    out_fields = ["Basecalls","Reads","Mask_read_depth","Frequencies"];
    out_file_ext = ".mat";

    out_path = paths.output;
    if ~exist(out_path, 'dir'); mkdir(out_path); end

    for out_field = out_fields
        out = data.(out_field);
        out_filepath = fullfile(out_path, append(out_field,out_file_ext));
        save(out_filepath, '-struct', 'out')
    end
end; clear -regexp ^out
%}  

%% STEP 5: Find consensus and excess base frequency components and read count fractions  
% Create output variables for base frequency shared and excess components
temp_frequency_shared_component = struct();
temp_frequency_excess_component = struct();
% shared component + excess component = total raw frequency

% Create output variables for read excess and shared fractions
temp_readcount_excess_fraction = struct(); % fraction of reads not supported bi-directionally (between 0 and 1]
temp_readcount_shared_fraction = struct(); % fraction of reads supported bi-directionally

% Calculate consensus/excess frequency components and read count fractions for each base type
for temp_base = ["A", "C", "G", "T"]

    % Get directional base frequencies (values between 0 and 1)
    temp_freq_F = data.Frequencies.For.(temp_base);
    temp_freq_R = data.Frequencies.Rev.(temp_base);

    % Replace NaN (no reads) with zeros
    temp_freq_F(isnan(temp_freq_F)) = 0;
    temp_freq_R(isnan(temp_freq_R)) = 0;

    % Find highest and lowest value in each  F and R base frequency value pair (at each site)
    temp_freq_min = min(temp_freq_F, temp_freq_R,'omitnan');
    temp_freq_max = max(temp_freq_F, temp_freq_R,'omitnan');

    % Find base frequency excess component: difference btween min and max directional frequencies
    temp_freq_excess = temp_freq_max - temp_freq_min; % can be from either of directions

    % Get directional base frequency excess components
    temp_freq_excess_F = temp_freq_F - temp_freq_min; % >0 at sites with excess calls in F direction
    temp_freq_excess_R = temp_freq_R - temp_freq_min; % >0 at sites with excess calls in R direction

    % save into struct
    temp_frequency_shared_component.(temp_base) = temp_freq_min; % same between directions!
    temp_frequency_excess_component.All.(temp_base) = temp_freq_excess;
    temp_frequency_excess_component.For.(temp_base) = temp_freq_excess_F;
    temp_frequency_excess_component.Rev.(temp_base) = temp_freq_excess_R;


    % Calculate F and R read fractions for this base that make up its exces frequency component
    temp_read_excess_fraction_F = temp_freq_excess_F ./ temp_freq_F;
    temp_read_excess_fraction_R = temp_freq_excess_R ./ temp_freq_R;

    temp_read_excess_fraction_F(isnan(temp_read_excess_fraction_F)) = 0; % replace NaNs with 0s
    temp_read_excess_fraction_R(isnan(temp_read_excess_fraction_R)) = 0; % replace NaNs with 0s

    % Get F and R read fractions for this base that make up its shared frequency component
    temp_read_shared_fraction_F = 1-temp_read_excess_fraction_F;
    temp_read_shared_fraction_R = 1-temp_read_excess_fraction_R;

    temp_read_shared_fraction_F(isnan(temp_read_shared_fraction_F)) = 0; % replace NaNs with 0s
    temp_read_shared_fraction_R(isnan(temp_read_shared_fraction_R)) = 0; % replace NaNs with 0s
    

    % Save read fractions into struct
    temp_readcount_excess_fraction.For.(temp_base) = temp_read_excess_fraction_F;
    temp_readcount_excess_fraction.Rev.(temp_base) = temp_read_excess_fraction_R;

    temp_readcount_shared_fraction.For.(temp_base) = temp_read_shared_fraction_F;
    temp_readcount_shared_fraction.Rev.(temp_base) = temp_read_shared_fraction_R;
end

% Calculate the sum excess base frequency components (sum excess components of the 4 bases together)
for temp_dir = ["All", "For", "Rev"]
    temp_frequency_excess_component.(temp_dir).sum = ...
        temp_frequency_excess_component.(temp_dir).A + ...
        temp_frequency_excess_component.(temp_dir).C + ...
        temp_frequency_excess_component.(temp_dir).G + ...
        temp_frequency_excess_component.(temp_dir).T;
end

% Get the sum shared frequency component (sum shared components of all 4 bases together)
temp_frequency_shared_component.sum = ...
    temp_frequency_shared_component.A + ...
    temp_frequency_shared_component.C + ...
    temp_frequency_shared_component.G + ...
    temp_frequency_shared_component.T;

% Sanity check:
% temp_frequency_shared_component.sum + temp_frequency_excess_component.For.sum = 1
% temp_frequency_shared_component.sum + temp_frequency_excess_component.Rev.sum = 1

% Save data into out struct
data.Read_fractions.Shared = temp_readcount_shared_fraction;
data.Read_fractions.Excess = temp_readcount_excess_fraction;

data.Base_frequency_components.Shared = temp_frequency_shared_component;
data.Base_frequency_components.Excess = temp_frequency_excess_component;

% Clear temp
clear -regexp ^temp


%% STEP 6: Calculate the consensus frequencies of each base type
% Consensus base frequencies - % of the shared frequency component sum made up by each base type
temp_consensus_base_frequencies = struct(); % Same for both directions!

for temp_base = ["A", "C", "G", "T"]
    % Get shared base and sum frequency component values
    temp_shared_component_base = data.Base_frequency_components.Shared.(temp_base);
    temp_shared_somponent_sum  = data.Base_frequency_components.Shared.sum;

    % Calculate consensus base frequency
    temp_consensus_freq = temp_shared_component_base ./ temp_shared_somponent_sum;
    temp_consensus_freq(isnan(temp_consensus_freq)) = 0; % replace NaNs (no consensus reads) with 0s

    % Save
    temp_consensus_base_frequencies.(temp_base) = temp_consensus_freq;
end

% Save base consensus frequencies into data struct
data.Consensus_frequencies = temp_consensus_base_frequencies;

% Clear temp
clear -regexp ^temp

%% STEP 7: Calculate the actual consensus and excess read counts
% Input variables: raw base call counts and calculated read count shared and excess fractions
temp_basecalls = data.Basecalls;
temp_fractions = data.Read_fraction;

% Output variables: 
temp_basecalls_shared   = struct(); % Each base calls making up consensus read fraction in each direction
temp_basecalls_excess   = struct(); % Each base calls making up excess read fraction in each direction
temp_total_reads_shared = struct(); % Total reads (A+C+G+T) in the shared fraction
temp_total_reads_excess = struct(); % Total reads (A+C+G+T) in the excess fraction

for temp_dir = ["For", "Rev", "All"]
    for temp_base = ["A", "C", "G", "T"]

        % Calculate shared and excess read counts
        if any(temp_dir == ["For", "Rev"])
           
            % Get raw basecall count and excess and shared read fractions in one direction
            temp_reads = temp_basecalls.(temp_dir).(temp_base);
            temp_excess_fr = temp_fractions.Excess.(temp_dir).(temp_base);
            temp_shared_fr = temp_fractions.Shared.(temp_dir).(temp_base);

            % Calculate shared and excess fractions of the basecall counts in one direction
            temp_reads_shared = temp_reads .* temp_shared_fr;
            temp_reads_excess = temp_reads .* temp_excess_fr;

        elseif temp_dir == "All"
            % Get shared and excess basecall counts in F and R directions
            temp_reads_shared_F = temp_basecalls_shared.For.(temp_base);
            temp_reads_shared_R = temp_basecalls_shared.Rev.(temp_base);
            temp_reads_excess_F = temp_basecalls_excess.For.(temp_base);
            temp_reads_excess_R = temp_basecalls_excess.Rev.(temp_base);

            % Calculate sum shared and sum excess basecall count across both directions
            temp_reads_shared = temp_reads_shared_F + temp_reads_shared_R;
            temp_reads_excess = temp_reads_excess_F + temp_reads_excess_R;
        end

        % Deposit into struct
        temp_basecalls_shared.(temp_dir).(temp_base) = temp_reads_shared;
        temp_basecalls_excess.(temp_dir).(temp_base) = temp_reads_excess;
    end

    % calculate total shared read count by summing consensus calls of each of the 4 bases
    temp_total_reads_shared.(temp_dir) = ...
        temp_basecalls_shared.(temp_dir).A + ...
        temp_basecalls_shared.(temp_dir).C + ...
        temp_basecalls_shared.(temp_dir).G + ...
        temp_basecalls_shared.(temp_dir).T;

    % calculate total excess read count by summing excess calls of each of the 4 bases
    temp_total_reads_excess.(temp_dir) = ...
        temp_basecalls_excess.(temp_dir).A + ...
        temp_basecalls_excess.(temp_dir).C + ...
        temp_basecalls_excess.(temp_dir).G + ...
        temp_basecalls_excess.(temp_dir).T;
end

% Save variables counts into the data struct
data.Consensus_reads = temp_total_reads_shared;
data.Excess_reads = temp_total_reads_excess;
data.Consensus_basecalls = temp_basecalls_shared;
data.Excess_basecalls = temp_basecalls_excess;

% Clear temp
clear -regexp ^temp
%%}

%% STEP 8: Save Consensus and excess read count and frequency variables to MAT file
if settings.save2 == true
    out_path = paths.output;
    if ~exist(out_path, 'dir'); mkdir(out_path); end
    
    out_fields = [...
        "Consensus_frequencies",...
        "Consensus_reads", ...
        "Excess_reads", ...
        "Consensus_basecalls", ...
        "Excess_basecalls",...
        "Read_fractions",...
        "Base_frequency_components"];
    
    % Save eachfield as a MAT variable
    for out_field = out_fields
        out_filepath = fullfile(out_path, append(out_field,".mat"));
        out = data.(temp_fn);
        save(out_filepath, '-struct', 'out');
    end
    clear -regexp ^out
end
    

