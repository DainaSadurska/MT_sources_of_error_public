%% SCRIPT_3_motif_error_rate_calculator
% Script for calculating motif assciated sequencing error rates
%
% INPUTS:
%  > Motif_spaces.mat  - Lists with the motifs to be searched for, motif space IDs: N2-N8
%  > SampleIDs.mat     - Sample IDs
%  > Reads.mat         - Total read count of all base types each mtDNA position
%  > Reads_Ref.mat     - Counts of reference reads at each mtDNA position
%  > Reads_NonRef.mat  - Counts of non-reference reads at each mtDNA position
%
% INPUTS for each motif space N2 to N8: 
%  > [MotifSpaceID]_Motif_locations.mat     < produced by SCRIPT_2_Motif_location_finder
%    >> variables: 'Motif_table', 'L_strand', 'H_strand'
%
%
% OUTPUTS for each motif space N2 to N8:  
%  > [MotifSpaceID]_Motif_table.mat         - Motif space overview with population-level motif error rate values 
%  > [MotifSpaceID]_Motif_instances.mat     - Motif counts analysed in each sample
%  > [MotifSpaceID]_Readcount_totals.mat    - Arrays with sample-level FM/FMM/RM/RMM totals and FER/RER/ERD values
%  > [MotifSpaceID]_Motif_error_rates.mat   - Population-level and Sample-level motif error rate (ERD) tables
%
%{
% Notation used: 
% FM  = Forward Match    - reads from the strand with the motif that match reference base
% FMM = Forward MisMatch - reads from the strand with the motif that do not match reference base
% RM  = Reverse Match    - reads from the opposite strand relative to motif orientation that match reference base
% RMM = Reverse MisMatch - reads from the opposite strand relative to motif orientation that do not match reference base
%
% RER (Reverse Error Rate) - Reverse mismatch fraction relative to motif orientation
% FER (Forward Error Rate) - Forward mismatch fraction relative to motif orientation
% ERD (Error Rate Difference) - Difference between motif forward and reverse mismatch rates
%
% ERD is calculated as:
% RER = RMM / (RM + RMM)
% FER = FMM / (FM + FMM)
% ERD = FER - RER
%
% IMPORTANT: F ('Forward') and R ('Reverse') in FM/FMM/RM/RMM indicate read orientation relative to 
% the MOTIT and do NOT relate to the sequencing direction. For distinguishing between reads arising 
% from forward and reverse sequencing directions L and H strand notation is used instead 
% (L strand > forward reads, H strand > reverse reads).
%}

%% Set paths & Load data
% Paths: outputs
paths = struct();
paths.root = '/path/to/main/directory/';
paths.output = fullfile(paths.root,"Output_2");
if ~exist(paths.output,'dir'); mkdir(paths.output); end

% Paths: inputs
paths.input1 = '/path/to/input/directory/';
paths.input2 = fullfile(paths.root,"Output_1");

in_paths = struct();
in_dir1 = paths.input1;
in_dir2 = paths.input2;

% Filepaths from input1
in_paths.Motif_spaces = fullfile(in_dir1, "Motif_spaces.mat");
in_paths.SampleIDs    = fullfile(in_dir1, "SampleIDs.mat");      
in_paths.Reads        = fullfile(in_dir1, "Reads.mat");
in_paths.Reads_Ref    = fullfile(in_dir1, "Reads_Ref.mat");
in_paths.Reads_NonRef = fullfile(in_dir1, "Reads_NonRef.mat");

in_motif_spaces = ['N2','N3','N4','N5','N6','N7','N8'];
for in_motif_space = in_motif_spaces
    in_paths.(in_motif_space).Locations = fullfile(in_dir2, append(in_motif_space,"_Motif_locations.mat"));
end

% Load data from input 1: Sample IDs and read counts
data = struct();
data.Motif_spaces       = load(in_paths.Motif_spaces); % 'N2','N3','N4','N5','N6','N7','N8'
data.SampleIDs          = load(in_paths.SampleIDs, 'SampleIDs');    
data.Reads              = load(in_paths.Reads, 'For','Rev');       % Total read count
data.Counts_Ref         = load(in_paths.Reads_Ref, 'For','Rev');   % Reference base read count
data.Counts_NonRef      = load(in_paths.Reads_NonRef,'For','Rev'); % Non-reference base read count

% Load files from input 2: Motif locations
for in_motif_space = in_motif_spaces
    data.(in_motif_space).Motif_list      = load(in_paths.(in_motif_space).Locations,'Motif_table');
    data.(in_motif_space).Motif_locations = load(in_paths.(in_motif_space).Locations,'L_strand','H_strand');
end
clear -regexp ^in_

%% Get dataset parameters
parameters = struct();

% All motif space IDs
parameters.motif_spaces = string(fieldnames(data.Motif_spaces))';

% Dataset dimensions
parameters.dim_samples = 1;      % Rows = samples
parameters.dim_positions = 2;    % Cols = positions

parameters.n_positions = size(data.Sequences,parameters.dim_positions);
parameters.n_samples   = size(data.Sequences,parameters.dim_samples);
parameters.dimensions  = [parameters.n_samples parameters.n_positions];

% Sample and position numbers
parameters.samples    = 1:parameters.n_samples;
parameters.positions  = 1:parameters.n_positions;

% Artefact site 'N' at position 3107
parameters.del_artefact_site = 3107;

% Minimum required directional read depth
parameters.min_reads_per_dir = 500;

% Maximum allowed RER - reverse error (mismatch) rate
parameters.max_RER = 0.01;

%% OPTS - set options for this run
opts = struct();
opts.motif_spaces = ["N2", "N3", "N4"]; % Which motif spaces to run (must be subset of parameters.motif_spaces)
opts.save_MAT = false;        % Save outputs to file?  If FALSE, outputs only deposited into DATA struct
opts.save_TXT = false;        % Write output tables to tab-delimited txt file? 
opts.time_check_interval = 4; % How frequenty time remaining is estimated

%% STEP 1 - Get the match and mismatch read count matrices (to be used for motif error rate calculation)
% Assemble reference base match and mismatch count matrices
in_counts_L = struct();
in_counts_H = struct();

% L strand (forward sequencing direction) FM/ FMM/ RM/ RMM count matrices
in_counts_L.FM  = data.Counts_Ref.For;     % FM  = Forward Match
in_counts_L.FMM = data.Counts_NonRef.For;  % FMM = Forward MisMatch
in_counts_L.RM  = data.Counts_Ref.Rev;     % RM  = Reverse Match
in_counts_L.RMM = data.Counts_NonRef.Rev;  % RMM = Reverse MisMatch

% H strand (reverse sequencing direction) FM/ FMM/ RM/ RMM counts matrices
in_counts_H.FM  = data.Counts_Ref.Rev;     % FM  = Forward Match
in_counts_H.FMM = data.Counts_NonRef.Rev;  % FMM = Forward MisMatch
in_counts_H.RM  = data.Counts_Ref.For;     % RM  = Reverse Match
in_counts_H.RMM = data.Counts_NonRef.For;  % RMM = Reverse MisMatch


%% STEP 2 - Filter read count data by (1) minimum read depth and (2) maximum reverse mismatch fraction
% Get filter thresholds
temp_min_reads = parameters.min_reads_per_dir; % Minimum 500 reads in each direction
temp_max_RER   = parameters.max_RER;           % Maximum 1% allowed Reverse Rrror Rate (i.e. reverse mismatch fraction)
fprintf("Filtering data: \n  - Excluding sites with < %d reads per direction\n  - Excluding sites with RER > %.1f%%\n",...
    temp_min_reads, temp_max_RER*100)

% Get total read counts in each direction
temp_reads_L = data.Reads.For;       % L strand total reads
temp_reads_H = data.Reads.Rev;       % H strand total reads

% Find locations where read depth is too low (<500) in either sequencing direction
temp_mask_read_depth_too_low = (temp_reads_L < temp_min_reads | temp_reads_H < temp_min_reads);

% Calculate L and H strand RER (Reverse Error Rate) for every position (RER = reverse mismatch fraction)
temp_RER_L = in_counts_L.RMM ./ (in_counts_L.RM + in_counts_L.RMM);
temp_RER_H = in_counts_H.RMM ./ (in_counts_H.RM + in_counts_H.RMM);

% Find locations in each direction where the appropriate direction specific RER exceeds maximum
temp_mask_L_RER_too_high = temp_RER_L  > temp_max_RER;
temp_mask_H_RER_too_high = temp_RER_H  > temp_max_RER;

% Create combined masks for excluding sites that fail either of the filters
in_mask_L_excluded_locations = temp_mask_L_RER_too_high | temp_mask_read_depth_too_low;
in_mask_H_excluded_locations = temp_mask_H_RER_too_high | temp_mask_read_depth_too_low;

% Apply filter to FORWARD (L strand) reads: insert NaN at sites that fail either of the 2 filter
in_counts_L.FM(in_mask_L_excluded_locations)  = NaN;
in_counts_L.FMM(in_mask_L_excluded_locations) = NaN;
in_counts_L.RM(in_mask_L_excluded_locations)  = NaN;
in_counts_L.RMM(in_mask_L_excluded_locations) = NaN;

% Apply filter to REVERSE (H strand) reads: insert NaN at sites that fail either of the 2 filter
in_counts_H.FM(in_mask_H_excluded_locations)  = NaN;
in_counts_H.FMM(in_mask_H_excluded_locations) = NaN;
in_counts_H.RM(in_mask_H_excluded_locations)  = NaN;
in_counts_H.RMM(in_mask_H_excluded_locations) = NaN;

% Clean-up
clear -regexp ^temp_

%% STEP 3 - Calculate motif counts and FER/RER/ERD values

% Iterate through motif spaces
for in_motif_space = opts.motif_spaces
    %% 3.1 - Get motif space data
    disp(append("Starting motif error rate calculation for motif space ",in_motif_space))    

    % Motif space overview table with count totals and found status
    in_motif_list = data.(in_motif_space).Motif_list;

    % Full motif list
    loop_motif_list = in_motif_list.Motif;
    [loop_n_motifs, loop_motif_size] = size(char(loop_motif_list));
    loop_motif_IDs = in_motif_list.Motif_ID;        % Motif number

    % Start motif space timer
    time_start = tic;
    time_check = opts.time_check_interval; % How frequently time remaining is exsimated    
    
    %% 3.2 - Get motif location variables and lists of motifs found / not found
    % Get motif location variables (onlt contain motifs found - i.e. may not be the full motif list)
    in_L_loc = data.(in_motif_space).Motif_locations.L_strand; % struct, contains only motifs found!
    in_H_loc = data.(in_motif_space).Motif_locations.H_strand; % struct

    %% 4.5 - Make empty arrays for motif-specific read counts during calculations
    loop_n_samples = parameters.n_samples;
    loop_n_positions = parameters.n_positions;

    % Same dimensions as input counts: rows - samples, columns - positions
    temp_dims = [loop_n_samples, loop_n_positions];
    loop_readcount_arrays = struct();
    loop_readcount_arrays.FM  = zeros(temp_dims);
    loop_readcount_arrays.FMM = zeros(temp_dims);
    loop_readcount_arrays.RM  = zeros(temp_dims);
    loop_readcount_arrays.RMM = zeros(temp_dims);
    clear -regexp ^temp_

    %% 3.3 - Create output containers for SAMPLE level motif FM/FMM/RM/RMM & FER/RER/ERD data
    % Output dimensions: rows - motifs, columns - samples
    temp_dims = [loop_n_motifs, parameters.n_samples]; 
    
    % Struct for depositing sample-level motif data outputs
    out_arrays = struct();

    % Arrays for saving motif sample-level FM/FMM/RM/RMM totals
    out_arrays.FM  = zeros(temp_dims);
    out_arrays.FMM = zeros(temp_dims);
    out_arrays.RM  = zeros(temp_dims);
    out_arrays.RMM = zeros(temp_dims);

    % Array for saving motif sample-level ERD (Forward and Reverse mismatch rate difference) values
    out_arrays.FER = zeros(temp_dims);
    out_arrays.RER = zeros(temp_dims);
    out_arrays.ERD = zeros(temp_dims);

    % Arrays for saving umber of motif instances analysed in each sample
    out_arrays.instances_L     = zeros(temp_dims);
    out_arrays.instances_H     = zeros(temp_dims);
    out_arrays.instances_total = zeros(temp_dims);
    clear -regexp ^temp_

    %% 3.4 - Create output table for POPULATION level motif FM/FMM/RM/RMM & FER/RER/ERD data
    % Table columns with the motif list 
    out_motif_table = in_motif_list(:,1:2); % Col 1: Motif ID, Col 2: Motif sequence

    % Motif status column: TRUE = at least 1 motif instance passed filters and was analysed across all samples 
    out_motif_table.Motif_analysed = false(loop_n_motifs,1);  % AFTER read depth and RER filtering

    % Columns for sample count & analysed motif instance totals 
    out_motif_table.Samples_with_motif  = zeros(loop_n_motifs,1);
    out_motif_table.Instances_L_strand  = zeros(loop_n_motifs,1); 
    out_motif_table.Instances_H_strand  = zeros(loop_n_motifs,1);
    out_motif_table.Instances_total     = zeros(loop_n_motifs,1);

    % Columns for population-level FM/FMM/RM/RMM count totals and population-level FER/RER/ERD values
    out_motif_table.Sum_FM   = zeros(loop_n_motifs,1);
    out_motif_table.Sum_FMM  = zeros(loop_n_motifs,1);
    out_motif_table.Sum_RM   = zeros(loop_n_motifs,1);
    out_motif_table.Sum_RMM  = zeros(loop_n_motifs,1);

    % Columns for population-level motif population-level FER/RER/ERD values
    out_motif_table.Sum_FER  = zeros(loop_n_motifs,1);
    out_motif_table.Sum_RER  = zeros(loop_n_motifs,1);
    out_motif_table.Sum_ERD  = zeros(loop_n_motifs,1);

    %% 3.5 - Iterate through motifs and calculate each motif FM/FMM/RM/RMM & FER/RER/ERD values
    % Iterate over all motifs
    for loop_motif_num = 1:loop_n_motifs

        % Start loop timer
        time_loop_start = tic;

        % Get motif sequence
        loop_motif_seq = in_motif_list.Motif(loop_motif_num); % motif sequence
        fprintf('Motif %d of %d: %s',loop_motif_num,loop_n_motifs,loop_motif_seq);

        % Get motif found status PRE-FILTER: allows to save time for larger motif spaces 
        loop_motif_found = in_motif_list.Found(loop_motif_num);
        loop_motif_analysed = loop_motif_found; % will be set to false if no motif instances pass filters

        % Only process motifs that were found in at least one samples
        if loop_motif_found == true % Motif has been found in at least some samples
            %% 3.5.1 - Create motif location masks (rows = samples, cols = positions)
            % Make empty arrays for saving motif locations - positions = cols, samples = rows
            loop_motif_locations_L = false([loop_n_samples, loop_n_positions]); 
            loop_motif_locations_H = false([loop_n_samples, loop_n_positions]); 

            % Combine sample-specific motif location data from all samples into a single motif locations mask
            for temp_sample = 1:loop_n_samples
                % Get motif locations in this sample
                temp_ind_loc_L = in_L_loc.(loop_motif_seq){temp_sample,1};
                temp_ind_loc_H = in_H_loc.(loop_motif_seq){temp_sample,1};

                % Set the corresponding motif location mask sites to TRUE
                loop_motif_locations_L(temp_sample,temp_ind_loc_L) = true;
                loop_motif_locations_H(temp_sample,temp_ind_loc_H) = true;
            end
            clear -regexp ^temp

            % Apply read depth and RER filter mask to exclude motif locations at low depth or high RER sites
            loop_motif_locations_L(in_mask_L_excluded_locations)=false;
            loop_motif_locations_H(in_mask_H_excluded_locations)=false;   

            %% 3.5.2 - Find the number of motif instances in each sample
            % Calculate number of filter-passing motif instances in L and R (per sample) 
            loop_motif_instances_L     = sum(loop_motif_locations_L,2,'omitnan');     % Column vector, dims: [samples 1]
            loop_motif_instances_H     = sum(loop_motif_locations_H,2,'omitnan');             
            loop_motif_instances_total = loop_motif_instances_L+loop_motif_instances_H; 

            % Calculate total number of samples with any filter-passing motif instances
            temp_n_samples_with_motif = numel(find(loop_motif_instances_total>0));     % dims: [1 1]
            
            % T/F indicating whether any instances of this motif passed filters
            loop_motif_analysed = temp_n_samples_with_motif >0;                        % dims: [1 1]            

            % Save number of samples with this motif into output table
            out_motif_table.Motif_analysed(loop_motif_num)     = loop_motif_analysed;
            out_motif_table.Samples_with_motif(loop_motif_num) = temp_n_samples_with_motif;

            % Transpose sample-level motif instance counts into row vectors and place into into output arrays
            out_arrays.Instances_L(loop_motif_num,:)     = loop_motif_instances_L';     % dims: [1 samples]
            out_arrays.Instances_H(loop_motif_num,:)     = loop_motif_instances_H';     % dims: [1 samples]
            out_arrays.Instances_total(loop_motif_num,:) = loop_motif_instances_total'; % dims: [1 samples]

            % Calculate total motif instances counts across all samples and save into output table
            out_motif_table.Instances_L_strand(loop_motif_num) = sum(loop_motif_instances_L,1);     % dims: [1 1]
            out_motif_table.Instances_H_strand(loop_motif_num) = sum(loop_motif_instances_H,1);     % dims: [1 1]
            out_motif_table.Instances_total(loop_motif_num)    = sum(loop_motif_instances_total,1); % dims: [1 1]
            
            
            %% 3.5.3 - If any motif instances pass filters, calculate motif error rates
            if loop_motif_analysed == true
                %% A) Calculate sample-specific motif error rate
                % Grab an empty copy of structs with arrays for motif read counts
                temp_motif_readcounts_L = loop_readcount_arrays; % array of zeros
                temp_motif_readcounts_H = loop_readcount_arrays;
    
                % Extract FM FMM RM & RMM counts from all locations where this motif is found
                for temp_field = ["FM","FMM","RM","RMM"]
                    temp_motif_readcounts_L.(temp_field)(loop_motif_locations_L) = in_counts_L.(temp_field)(loop_motif_locations_L);
                    temp_motif_readcounts_H.(temp_field)(loop_motif_locations_H) = in_counts_H.(temp_field)(loop_motif_locations_H);
                end
    
                % Calculate motif Sample-level FM, FMM, RM & RMM read counts on L strand
                temp_L_FM  = sum(temp_motif_readcounts_L.FM, 2,'omitnan'); % Column vector, dims: [samples, 1]
                temp_L_FMM = sum(temp_motif_readcounts_L.FMM,2,'omitnan'); 
                temp_L_RM  = sum(temp_motif_readcounts_L.RM, 2,'omitnan'); 
                temp_L_RMM = sum(temp_motif_readcounts_L.RMM,2,'omitnan'); 
                
                % Calculate motif Sample-level FM, FMM, RM & RMM read counts on H strand
                temp_R_FM  = sum(temp_motif_readcounts_H.FM, 2,'omitnan'); % Column vector, dims: [samples, 1]
                temp_R_FMM = sum(temp_motif_readcounts_H.FMM,2,'omitnan'); 
                temp_R_RM  = sum(temp_motif_readcounts_H.RM, 2,'omitnan'); 
                temp_R_RMM = sum(temp_motif_readcounts_H.RMM,2,'omitnan'); 
    
                % Calculate motif Sample-level FM, FMM, RM & RMM sums across both strands
                temp_FM  = temp_L_FM  + temp_R_FM;   % Column vector, dims: [samples, 1]
                temp_FMM = temp_L_FMM + temp_R_FMM;
                temp_RM  = temp_L_RM  + temp_R_RM;
                temp_RMM = temp_L_RMM + temp_R_RMM;
    
                % Calculate motif sample-level FER and RER values
                temp_FER = temp_FMM ./ (temp_FM + temp_FMM); % Column vector, dims: [samples, 1]
                temp_RER = temp_RMM ./ (temp_RM + temp_RMM); 
    
                % Calculate motif sample-level ERD value
                temp_ERD = temp_FER - temp_RER;   % Column vector, dims: [samples, 1]

                % SAVE: Transpose sample-level values into row vectors and place into the output arrays,
                out_arrays.FM(loop_motif_num,:)  = temp_FM';
                out_arrays.FMM(loop_motif_num,:) = temp_FMM';
                out_arrays.RM(loop_motif_num,:)  = temp_RM';
                out_arrays.RMM(loop_motif_num,:) = temp_RMM';
    
                out_arrays.FER(loop_motif_num,:) = temp_FER';
                out_arrays.RER(loop_motif_num,:) = temp_RER';
                out_arrays.ERD(loop_motif_num,:) = temp_ERD';
    
                %% B) Calculate the overall population-level motif error rate
                % Calculate motif population-level FM, FMM, RM & RMM read count totals - output is scalar
                temp_SUM_FM  = sum(temp_FM,1);    % dims: [1, 1]
                temp_SUM_FMM = sum(temp_FMM,1);
                temp_SUM_RM  = sum(temp_RM,1);
                temp_SUM_RMM = sum(temp_RMM,1);
    
                % Calculate motif population-level FER and RER values
                temp_SUM_FER = temp_SUM_FMM ./ (temp_SUM_FM + temp_SUM_FMM); % dims: [1, 1]
                temp_SUM_RER = temp_SUM_RMM ./ (temp_SUM_RM + temp_SUM_RMM);
    
                % Calculate motif population-level ERD
                temp_SUM_ERD = temp_SUM_FER - temp_SUM_RER;

                % SAVE: Place population-level match/mismatch and FER/RER/ERD values into output table
                out_motif_table.Sum_FM(loop_motif_num)  = temp_SUM_FM;
                out_motif_table.Sum_FMM(loop_motif_num) = temp_SUM_FMM;
                out_motif_table.Sum_RM(loop_motif_num)  = temp_SUM_RM;
                out_motif_table.Sum_RMM(loop_motif_num) = temp_SUM_RMM;
                out_motif_table.Sum_FER(loop_motif_num) = temp_SUM_FER;
                out_motif_table.Sum_RER(loop_motif_num) = temp_SUM_RER;
                out_motif_table.Sum_ERD(loop_motif_num) = temp_SUM_ERD;
            end
        end
        % Print elapsed time per motif
        if loop_motif_found && loop_motif_analysed
             temp_time_iteration = seconds(toc(time_loop_start));
             temp_time_iteration.Format = 's'; %set duration format as seconds
             fprintf('   \t%s', string(temp_time_iteration));
        else;fprintf(' not found');
        end

        % Estimate time remaining every 25 iterations
        if mod(loop_motif_num,time_check) == 0
            % get total time elapsed
            temp_iteration = loop_motif_num;
            temp_n_iterations = loop_n_motifs;
            temp_time_elapsed = seconds(toc(time_start));
            temp_time_elapsed.Format = 'dd:hh:mm:ss'; %set duration format
    
            % get average iteration length
            temp_time_avg_iteration = temp_time_elapsed./temp_iteration;
            temp_time_avg_iteration.Format = 's'; %set duration format
    
            % get total remaining time (seconds)
            temp_time_remaining = ((temp_n_iterations-temp_iteration)*temp_time_avg_iteration);
            temp_time_remaining.Format = 'dd:hh:mm:ss';   %set duration format
    
            fprintf("\tTime elapsed: %s, Estimated time remaining: %s", temp_time_elapsed, temp_time_remaining) 
        end;fprintf('\n');

        % Clear temp
        clear -regexp ^temp
    end
    disp(append("Motif space ",in_motif_space," completed"))  

    %% 3.6 - Save motif space outputs into DATA struct
    % Save motif summary table with population level counts and totals to data struct
    data.(in_motif_space).Motif_table = out_motif_table(:,[1:4,7:14]); % skip directional instances

    % Save sample- and population-level error rate tables to data struct
    data.(in_motif_space).Population_MER = out_motif_table;  % includes directional instance counts
    temp_table = array2table(out_arrays.ERD);
    temp_table.Properties.VariableNames = data.SampleIDs;
    data.(in_motif_space).Sample_MER = [in_motif_list(:,1:2), temp_table];

    % Save motif instance count arrays to data struct
    data.(in_motif_space).Motif_counts.L_strand = out_arrays.Instances_L;
    data.(in_motif_space).Motif_counts.H_strand = out_arrays.Instances_H;
    data.(in_motif_space).Motif_counts.Total    = out_arrays.Instances_total;

    % Save arrays with motif readcount totals to data struct
    data.(in_motif_space).Readcount_totals.FM  = out_arrays.FM;
    data.(in_motif_space).Readcount_totals.FMM = out_arrays.FMM;
    data.(in_motif_space).Readcount_totals.RM  = out_arrays.RM;
    data.(in_motif_space).Readcount_totals.RMM = out_arrays.RMM;
    data.(in_motif_space).Readcount_totals.FER = out_arrays.FER;
    data.(in_motif_space).Readcount_totals.RER = out_arrays.RER;
    data.(in_motif_space).Readcount_totals.ERD = out_arrays.ERD;

    % Clean-up
    clear -regexp ^out;
    clear -regexp ^loop;
    clear -regexp ^temp;
    clear -regexp ^time;    
    %}


    %% 3.7 - Save outputs to file 
    % Save output variables to MAT files
    if opts.save_MAT == true
        disp(append("Saving motif space ",in_motif_space," output variables to .mat file"))

        % output path
        out_path = paths.output;
        if ~exist(out_path,'dir'),mkdir(out_path),end
        
        % Set MAT variable file names
        temp_filepath_instances   = fullfile(out_path,append(in_motif_space,'_Motif_instances.mat'));
        temp_filepath_readcounts  = fullfile(out_path,append(in_motif_space,'_Readcount_totals.mat'));
        temp_filepath_error_rates = fullfile(out_path,append(in_motif_space,'_Motif_error_rates.mat'));
        temp_filepath_motif_table = fullfile(out_path,append(in_motif_space,'_Motif_table.mat'));

        % Save motif instance count variables
        out = data.(in_motif_space).Motif_counts;
        out.Motif_list = data.(in_motif_space).Motif_table(:,1:2);
        out.SampleIDs = data.SampleIDs;
        save(temp_filepath_instances, '-struct', 'out');

        % Save basecount totals arrays (includes RER/FER/ERD arrays)
        out = data.(in_motif_space).Readcount_totals;
        out.Motif_list = data.(in_motif_space).Motif_table(:,1:2);
        out.SampleIDs = data.SampleIDs;
        save(temp_filepath_readcounts, '-struct', 'out');

        % Save motif population- and sample-lever error rate tables
        out = struct();
        out.Population_MER = data.(in_motif_space).Population_MER;
        out.Sample_MER     = data.(in_motif_space).Sample_MER;
        save(temp_filepath_error_rates, '-struct', 'out');

        % Save motif table
        out = struct();
        out.Motif_table = data.(in_motif_space).Motif_table;
        save(temp_filepath_motif_table, '-struct', 'out');

        % Clear temp out
        clear -regexp ^temp
        clear -regexp ^out
    end

    % Write output tables as TXT files
    if opts.save_TXT == true
        disp(append("Writing motif space ",in_motif_space," outputs to .txt file"))

        % output path
        out_path = paths.output;
        if ~exist(out_path,'dir'),mkdir(out_path),end

        % Save motif space overview table with population level counts and error rates  
        temp_filepath_motif_table = fullfile(out_path,append(in_motif_space,'_Motif_table.txt'));
        writetable(out_motif_table, temp_filepath_motif_table, 'Delimiter','\t');

        % Make sample level data tables (motifs x samples) with sample IDs as variable names
        out_tables = struct();
        temp_arrays = data.(in_motif_space).Readcount_totals;
        temp_arrays.Instances_L     = data.(in_motif_space).Motif_counts.L_strand;
        temp_arrays.Instances_H     = data.(in_motif_space).Motif_counts.H_strand;
        temp_arrays.Instances_total = data.(in_motif_space).Motif_counts.Total;
        for temp_field = string(fieldnames(temp_arrays))'
            temp_table = array2table(temp_arrays.(temp_field));
            temp_table.Properties.VariableNames = data.SampleIDs;
            out_tables.(temp_field) = [data.(in_motif_space).Motif_table(:,1:2), temp_table];
        end

        % Wrtie Tables with counts of motif instances in each sample
        temp_filepath_counts_L     = fullfile(out_path,append(in_motif_space,'_Motif_counts_L_strand.txt'));
        temp_filepath_counts_H     = fullfile(out_path,append(in_motif_space,'_Motif_counts_H_strand.txt'));
        temp_filepath_counts_total = fullfile(out_path,append(in_motif_space,'_Motif_counts_total.txt'));

        writetable(out_tables.Instances_L,    temp_filepath_counts_L,    'Delimiter','\t')
        writetable(out_tables.Instances_H,    temp_filepath_counts_H,    'Delimiter','\t')
        writetable(out_tables.Instances_total,temp_filepath_counts_total,'Delimiter','\t')

        % Write Tables with sample-level FM FMM RM RMM FER RER and ERD data
        temp_filepath_FM  = fullfile(out_path,append(in_motif_space,'_FM.txt' ));
        temp_filepath_FMM = fullfile(out_path,append(in_motif_space,'_FMM.txt'));
        temp_filepath_RM  = fullfile(out_path,append(in_motif_space,'_RM.txt' ));
        temp_filepath_RMM = fullfile(out_path,append(in_motif_space,'_RMM.txt'));
        temp_filepath_FER = fullfile(out_path,append(in_motif_space,'_FER.txt'));
        temp_filepath_RER = fullfile(out_path,append(in_motif_space,'_RER.txt'));
        temp_filepath_ERD = fullfile(out_path,append(in_motif_space,'_ERD.txt'));

        writetable(out_tables.FM, temp_filepath_FM, 'Delimiter','\t');
        writetable(out_tables.FMM,temp_filepath_FMM,'Delimiter','\t');
        writetable(out_tables.RM, temp_filepath_RM, 'Delimiter','\t');
        writetable(out_tables.RMM,temp_filepath_RMM,'Delimiter','\t');
        writetable(out_tables.FER,temp_filepath_FER,'Delimiter','\t');
        writetable(out_tables.RER,temp_filepath_RER,'Delimiter','\t');
        writetable(out_tables.ERD,temp_filepath_ERD,'Delimiter','\t');

        % Clear temp out
        clear -regexp ^temp
        clear -regexp ^out
    end
    
end
disp("All motif spaces completed")
clear -regexp ^in_









