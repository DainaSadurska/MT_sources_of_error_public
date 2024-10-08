%% SCRIPT 7: Heteroplamy_filtering_CONSENSUS_METHOD.m
% Script for filtering heteroplasmies based on the consensus read counts and frequencies.
% This script filters heteroplasmies according to CONSENSUS METHOD principles, and the equivalent 
% script for standard method is SCRIPT 5.
% Inputs for this script were generated by SCRIPT 4
%
% FILTER SETTINGS USED:
% Filter 1 - Allele supported by at least 1 consensus read
% Filter 2 - Allele supported by at least 1 consensus read in each sequencing direction
% Filter 3 - Appele supported by minimum 5 consensus reads in each sequencing direction
% Filter 4 - Minimum 5 consensus reads in each diection + overall allele consensus frequency must be at least 0.1%
% Filter 5 - Minimum 5 consensus reads in each diection + allele consensus frequency must be at least 0.1% in each sequencing direction
% Filter 6 - Minimum 5 consensus reads in each diection + overall allele consensus frequency must be at least 1%
% Filter 7 - Minimum 5 consensus reads in each diection + allele consensus frequency must be at least 1% in each sequencing direction
%
% INPUTS from Output_4/Consensus_method:
%   > Samples           - Sample numbers and IDs
%   > Positions         - Position numbers and rCRS bases
%   > Sequences         - Sample consensus sequences
%   > Excess_reads      - Sum excess reads across all alleles
%   > Consensus_reads   - Sum consensus reads across all alleles
%   > Allele_order.mat                  - Allele order by consensus reads
%   > Allele_count.mat                  - Number of alleles detected at each position
%   > Allele_bases.mat                  - Sorted allele bases 
%   > Allele_consensus_reads.mat        - Sorted allele consensus reads 
%   > Allele_excess_reads.mat           - Sorted allele excess reads
%   > Allele_consensus_frequencies.mat  - Sorted allele frequencies
%   > Parameters.mat   
%
% OUTPUTS 
% in Output_4/Consensus_method
%   > Heteroplasmy_filters
%   > Site_masks
%   > Excess_reads
%   > Consensus_reads
%
% in each of the subdrectories Filter_1 to Filter_7:
%   > Allele_count        - Number of alleles passing the specific filter
%   > Allele_masks        - indices for extracting the filter-passing allele data
%   > Allele_bases        - filter passing allelle bases
%   > Allele_consensus_reads        - filter passing allele consensus reads
%   > Allele_excess_reads           - filter passing allele excess reads
%   > Allele_consensus_frequencies  - filter passing allele frequencies
%   > Heteroplasmic_sites - index of sites with 2 or more filter-passing allelels
%   > Homoplasmic_sites   - index of sites with 1 filter-passing allelele
%
% in sub-directory NoFilter:
%   > Allele_count        - Any alleles observed, with no other filters
%   > Allele_masks        - indices for all non-filtered alleles
%   > Allele_bases        - bases of all non-filtered alleles
%   > Allele_consensus_reads        - non-filtered allele consensus reads
%   > Allele_excess_reads           - non-filtered allele excess reads
%   > Allele_consensus_frequencies  - non-filtered allele frequencies
%   > Heteroplasmic_sites - index of sites with 2 or more alleles
%   > Homoplasmic_sites   - index of sites with 1 allelee
%


%% STEP 1: Paths and Load
% 1) Load data
paths.root = '/path/to/main/direcory/';
paths.output = fullfile(paths.root,'Output_5/Consensus_method');

% Set input paths
in_paths = struct();
in_dir = 'Output_4/Standard_method';

% from Output_4/Consensus_method
in_paths.Parameters         = fullfile(paths.root,in_dir,'Parameters.mat');
in_paths.Samples            = fullfile(paths.root,in_dir,'Samples.mat');
in_paths.Positions          = fullfile(paths.root,in_dir,'Positions.mat');
in_paths.Sequences          = fullfile(paths.root,in_dir,'Sequences.mat');
in_paths.Consensus_reads    = fullfile(paths.root,in_dir,'Consensus_reads.mat');
in_paths.Excess_reads       = fullfile(paths.root,in_dir,'Excess_reads.mat');
in_paths.Allele_count       = fullfile(paths.root,in_dir,'Allele_count.mat');
in_paths.Allele_bases       = fullfile(paths.root,in_dir,'Allele_bases.mat');
in_paths.Allele_consensus_reads       = fullfile(paths.root,in_dir,'Allele_consensus_reads.mat');
in_paths.Allele_excess_reads          = fullfile(paths.root,in_dir,'Allele_excess_reads.mat');
in_paths.Allele_consensus_frequencies = fullfile(paths.root,in_dir,'Allele_consensus_frequencies.mat');
in_paths.Allele_excess_frequencies    = fullfile(paths.root,in_dir,'Allele_excess_frequencies.mat');


% Create data struct
if ~exist('data', 'var'); data = struct(); end

% Load variables
data.Mask_del        = load(in_paths.Mask_del,'Mask_del');
data.Mask_read_depth = load(in_paths.Mask_read_depth,'Mask_read_depth');
data.SampleIDs       = load(in_paths.Samples, 'SampleIDs');
data.Individuals     = load(in_paths.Samples, 'Individuals');
data.Positions       = load(in_paths.Positions, 'Positions');
data.Reference       = load(in_paths.Sequences, 'ConcensusSeqs');
data.rCRS            = load(in_paths.Sequences, 'rCRS');
data.Consensus_reads    = load(in_paths.Consensus_reads,'For', 'Rev', 'All');
data.Excess_reads       = load(in_paths.Excess_reads,'For', 'Rev', 'All');
data.Allele_count       = load(in_paths.Allele_count, 'Allele_count');
data.Allele_bases       = load(in_paths.Allele_bases, 'Major','Minor_1','Minor_2','Minor_3');
data.Allele_consensus_reads       = load(in_paths.Allele_consensus_reads,'For', 'Rev', 'All');
data.Allele_excess_reads          = load(in_paths.Allele_excess_reads,'For', 'Rev', 'All');
data.Allele_consensus_frequencies = load(in_paths.Allele_consensus_frequencies, 'Major','Minor_1','Minor_2','Minor_3');
data.Allele_excess_frequencies    = load(in_paths.Allele_excess_frequencies, 'Major','Minor_1','Minor_2','Minor_3');
parameters              = load(in_paths.Parameters);

% Clear input variables
clear -regexp ^in

% 2) Set parameters
parameters = parameters.Parameters;

% rCRS deletion artefact site
parameters.deletion_site = 3107;

% Minimum required directional read depth
parameters.min_read_depth = 500;

%% STEP 2: Make masks for excluding problematic sites: position 3107 / insufficient depth / deletions 
% sites with low depth
temp_idx_low_depth = (...
    data.Consensus_reads.For < parameters.min_read_depth | ...
    data.Consensus_reads.Rev < parameters.min_read_depth );

% Sites with reference deletions or artefacts
temp_idx_deletions = (data.Reference == '-' | data.Reference == 'N');
temp_idx_deletions(parameters.deletion_site,:) = true;

% Assemble mask for sites that pass all criteria
temp_mask_pass = ~temp_idx_low_depth & ~temp_idx_deletions;

% Save masks to Data struct
data.Site_masks.FAIL_low_depth = temp_idx_low_depth;
data.Site_masks.FAIL_deletions = temp_idx_deletions;
data.Site_masks.FAIL = ~temp_mask_pass;
data.Site_masks.PASS = temp_mask_pass;
clear -regexp ^temp_

%% STEP 3: Set heteroplasmic site filter settings
% Create a table with different filter settings
temp_filter_settings = struct();
temp_filter_settings.Filter = compose("Filter_%d",1:7)';

% MINIMUM THRESHOLDS SET FOR FILTERS: F1   F2    F3     F4     F5     F6     F7
temp_filter_settings.Consensus_reads_overall  =[ 1;   2;   10;    10;    10;   10;    10];
temp_filter_settings.Consensus_reads_each_dir =[ 0;   1;    5;     5;     5;    5;     5];
temp_filter_settings.Freq_overall   =[ 0;   0;    0; 0.001; 0.001; 0.01;  0.01];
temp_filter_settings.Freq_each_dir  =[ 0;   0;    0;     0; 0.001;    0;  0.01];

temp_filter_settings = struct2table(temp_filter_settings);
data.Heteroplasmy_filters = temp_filter_settings;

% Clear temp 
clear -regexp ^temp



%% STEP 4: Identify alleles passing each set of filter settings
% Initialise struct for heteroplasmy masks and allele counts
temp_allele_masks = struct();  % logical, 1 = allel passes filter
temp_allele_counts = struct(); % number of alleles per position that pass each filter

% Set allele names to analyse
temp_alleles = ["Major", "Minor_1", "Minor_2", "Minor_3"];

% Get allele filter settings
temp_filters = data.Heteroplasmy_filters;

% Iterate through heteroplasmy filter settings
for temp_filter_num = 1:numel(temp_filters.Filter)

    % Get filter settings
    temp_filter            = temp_filters.Filter(temp_filter_num);
    temp_min_reads_total   = temp_filters.Consensus_reads_overall(temp_filter_num);
    temp_min_reads_per_dir = temp_filters.Consensus_reads_each_dir(temp_filter_num);
    temp_min_freq_total    = temp_filters.Freq_overall(temp_filter_num);
    temp_min_freq_per_dir  = temp_filters.Freq_each_dir(temp_filter_num);

    % Iterate though alleles
    for temp_allele = temp_alleles
        % Get allele frequency and read matrices
        temp_cons_reads  = data.Allele_consensus_reads.(temp_allele);
        temp_allele_freqs = data.Allele_consensus_frequencies.(temp_allele);
    
        % Find indices of alleles passing each filter
        temp_mask = ... % 1 = allele passes filter, 0 = Allele fails filter
            temp_cons_reads.All >= temp_min_reads_total   & ... % Apply min reads overall
            temp_cons_reads.For >= temp_min_reads_per_dir & ... % Apply min reads in F dir
            temp_cons_reads.Rev >= temp_min_reads_per_dir & ... % Apply min reads in R dir
            temp_allele_freqs.All >= temp_min_freq_total    & ... % Apply min frequency overall
            temp_allele_freqs.For >= temp_min_freq_per_dir  & ... % Apply min frequency in F dir
            temp_allele_freqs.Rev >= temp_min_freq_per_dir;       % Apply min frequency in R dir

        % Put filter passing allele mask  into struct
        temp_allele_masks.(temp_filter).(temp_allele) = temp_mask;
    end

    % Get number of alleles at each position passing the filter
    temp_allele_counts.(temp_filter) = ...
        temp_allele_masks.(temp_filter).Major + ...
        temp_allele_masks.(temp_filter).Minor_1 + ...
        temp_allele_masks.(temp_filter).Minor_2 + ...
        temp_allele_masks.(temp_filter).Minor_3;
end

% Save variables into data struct
data.Filtered_allele_masks  = temp_allele_masks;
data.Filtered_allele_counts = temp_allele_counts;

% Clear temp 
clear -regexp ^temp

%% STEP 5: Exclude low depth and deletion sites
temp_alleles = ["Major", "Minor_1", "Minor_2", "Minor_3"];
temp_filters = data.Heteroplasmy_filters.Filter; % Get filter names

% Retrieve heteroplasmy masks
temp_allele_masks  = data.Filtered_allele_masks;
temp_allele_counts = data.Filtered_allele_counts;

% Get the mask for sites that PASS depth and deletion filters
temp_mask_pass = data.Site_masks.PASS; 

for temp_filter = temp_filters'
    % Apply filter to each allele-specific mask
    for temp_allele = temp_alleles
        temp_allele_masks.(temp_filter).(temp_allele) = ...
            temp_allele_masks.(temp_filter).(temp_allele) & temp_mask_pass;
    end

    % Apply site filter to allele count matrices
    temp_allele_counts.(temp_filter)(~temp_mask_pass) = NaN;
end

% Save variables into data struct
data.Filtered_allele_masks  = temp_allele_masks;
data.Filtered_allele_counts = temp_allele_counts;

% Clear temp 
clear -regexp ^temp

%% STEP 6: Apply read depth & deletion filtration to Allele reads, Allele freq, Allele bases, Allele count
% Get mask for all sites pasing min depth and deletion filtering
temp_mask_pass = data.Site_masks.PASS;

% Filter allele count
data.Allele_count(~temp_mask_pass) = NaN;

for temp_allele = ["Major", "Minor_1", "Minor_2", "Minor_3"]
    % Filter allele bases
    data.Allele_bases.(temp_allele)(~temp_mask_pass) = 'N';

    % Filter allele reads and allele frequencies
    for temp_dir = ["All", "For", "Rev"]
        data.Allele_consensus_reads.(temp_allele).(temp_dir)(~temp_mask_pass) = NaN;
        data.Allele_excess_reads.(temp_allele).(temp_dir)(~temp_mask_pass) = NaN;
        data.Allele_consensus_frequencies.(temp_allele).(temp_dir)(~temp_mask_pass) = NaN;
    end
end

% Filter total reads
for temp_dir = ["All", "For", "Rev"]
    data.Consensus_reads.(temp_allele).(temp_dir)(~temp_mask_pass) = NaN;
    data.Excess_reads.(temp_allele).(temp_dir)(~temp_mask_pass) = NaN;
end

% Clear temp 
clear -regexp ^temp

%% STEP 7: Extract allele reads and frequencies for each filter and put into appropriate struct
for temp_filter = data.Heteroplasmy_filters.Filter'
    % Create output struct
    data.(temp_filter) = struct();

    % Save correct filter masks and allele counts
    data.(temp_filter).Allele_count = data.Filtered_allele_counts.(temp_filter);
    data.(temp_filter).Allele_masks  = data.Filtered_allele_masks.(temp_filter);
    disp(temp_filter)
    for temp_allele = ["Major", "Minor_1", "Minor_2", "Minor_3"]
        disp(temp_allele)
        % Get mask for where this allele passes the filter
        temp_allele_mask = data.Filtered_allele_masks.(temp_filter).(temp_allele);

        % Get filter passing allele bases & save into appropriate FILTER struct
        temp_base = data.Allele_bases.(temp_allele);
        temp_base(~temp_allele_mask) = ' ';
        data.(temp_filter).Allele_bases.(temp_allele) = temp_base;

        % Filter allele reads frequencies
        for temp_dir = ["All", "For", "Rev"]
            % Get filter-passing allele reads
            temp_cons_reads = data.Allele_consensus_reads.(temp_allele).(temp_dir);
            temp_cons_reads(~temp_allele_mask) = NaN;

            temp_excs_reads = data.Allele_excess_reads.(temp_allele).(temp_dir);
            temp_excs_reads(~temp_allele_mask) = NaN;
        
            % Get filter-passing allele frequencies
            temp_frequency = data.Allele_consensus_frequencies.(temp_allele).(temp_dir);
            temp_frequency(~temp_allele_mask) = NaN;

            % Save filtered reads & frequencies into appropriate FILTER struct
            data.(temp_filter).Allele_consensus_reads.(temp_allele).(temp_dir)       = temp_cons_reads;
            data.(temp_filter).Allele_excess_reads.(temp_allele).(temp_dir)          = temp_excs_reads;
            data.(temp_filter).Allele_consensus_frequencies.(temp_allele).(temp_dir) = temp_frequency;
        end
    end

    % Save correct filter heteroplasmic & homoplasmic site masks
    data.(temp_filter).Heteroplasmic_sites = data.(temp_filter).Allele_count >=2;
    data.(temp_filter).Homoplasmic_sites   = data.(temp_filter).Allele_count ==1;
end

% Clear temp 
clear -regexp ^temp

%% STEP 8: Calculate filtered minor alelle sum read counts and frequencies
for temp_filter = data.Heteroplasmy_filters.Filter'
    for temp_dir = ["All", "For", "Rev"]
        for temp_var = ["Allele_consensus_reads","Allele_excess_reads","Allele_consensus_frequencies"]
        temp_minor_1 = data.(temp_filter).(temp_var).Minor_1.(temp_dir);
        temp_minor_2 = data.(temp_filter).(temp_var).Minor_2.(temp_dir);
        temp_minor_3 = data.(temp_filter).(temp_var).Minor_3.(temp_dir);

        temp_NaN_mask1 = isnan(temp_minor_1);
        temp_NaN_mask2 = isnan(temp_minor_2);
        temp_NaN_mask3 = isnan(temp_minor_3);

        temp_minor_1(temp_NaN_mask1) = 0;
        temp_minor_2(temp_NaN_mask2) = 0;
        temp_minor_3(temp_NaN_mask3) = 0;
        temp_minor_sum = temp_minor_1+temp_minor_2+temp_minor_3;
        temp_minor_sum(temp_NaN_mask1) = NaN;
        
        data.(temp_filter).(temp_var).Minor_sum.(temp_dir) = temp_minor_sum;
        end
    end
end
clear -regexp ^temp

%% STEP 9: Calculate unfiltered minor alelle sum read counts and frequencies
for temp_dir = ["All", "For", "Rev"]
    for temp_var = ["Allele_consensus_reads","Allele_excess_reads","Allele_consensus_frequencies"]
        temp_minor_1 = data.(temp_var).Minor_1.(temp_dir);
        temp_minor_2 = data.(temp_var).Minor_2.(temp_dir);
        temp_minor_3 = data.(temp_var).Minor_3.(temp_dir);

        temp_NaN_mask1 = isnan(temp_minor_1);
        temp_NaN_mask2 = isnan(temp_minor_2);
        temp_NaN_mask3 = isnan(temp_minor_3);

        temp_minor_1(temp_NaN_mask1) = 0;
        temp_minor_2(temp_NaN_mask2) = 0;
        temp_minor_3(temp_NaN_mask3) = 0;
        temp_minor_sum = temp_minor_1+temp_minor_2+temp_minor_3;
        temp_minor_sum(temp_NaN_mask1) = NaN;
        
        data.(temp_var).Minor_sum.(temp_dir) = temp_minor_sum;
    end
end
clear -regexp ^temp

%% STEP 10: Assemble a 'NoFilter struct'
% Assemble a struct for pre-filter data
temp_allele_masks = struct();
temp_allele_bases = data.Allele_bases;

for temp_dir = ["All", "For", "Rev"] 
    for temp_allele = ["Major", "Minor_1", "Minor_2", "Minor_3"]
        temp_mask_N = temp_allele_bases.(temp_allele) =='N';
        temp_allele_bases.(temp_allele)(temp_mask_N) = ' ';
        temp_allele_masks.(temp_allele).(temp_dir) = ~temp_mask_N;
    end
end

data.NoFilter = struct();
data.NoFilter.Allele_count = data.Allele_count;
data.NoFilter.Allele_masks = temp_allele_masks;
data.NoFilter.Allele_bases = temp_allele_bases;
data.NoFilter.Allele_consensus_reads = data.Allele_consensus_reads;
data.NoFilter.Allele_excess_reads = data.Allele_excess_reads;
data.NoFilter.Allele_consensus_frequencies = data.Allele_consensus_frequencies;
data.NoFilter.Heteroplasmic_sites = data.Allele_count >=2;
data.NoFilter.Homoplamic_sites = data.Allele_count == 1;
clear -regexp ^temp

%% STEP 11: SAVE variables
% 1) Save filter settigns and site masks
out_path = paths.output;
if ~exist(out_path,'dir'); mkdir(out_path); end

% Heteroplasmy filters
out_field = "Heteroplasmy_filters";
out_filepath = fullfile(out_path, append(out_field,".mat"));
save(out_filepath, '-struct', 'data', out_field);

% Site masks
out_field = "Site_masks";
out_filepath = fullfile(out_path, append(out_field,".mat"));
out = struct();
out.Sites_included = data.Site_masks.PASS;
out.Sites_excluded = data.Site_masks.FAIL;
out.Mask_FAIL_low_depth = data.Site_masks.FAIL_low_depth;
out.Mask_FAIL_deletions = data.Site_masks.FAIL_deletions;
save(out_filepath, '-struct', 'out');

% Consensus_reads, Excess_reads
for out_field = ["Consensus_reads","Excess_reads"]
    out_filepath = fullfile(out_path, append(out_field,".mat"));
    out = data.out_field;
    save(out_filepath, '-struct', 'out');
end


% 2) Save data from each filter in separate directory
out_sets = data.Heteroplasmy_filters.Filter';
for out_set = out_sets
    out_path = fullfile(paths.output,out_set);
    if ~exist(out_path,'dir'); mkdir(out_path); end
    
    out = data.(out_set);
    for out_field = string(fieldnames(out))'
        out_filepath = fullfile(out_path, append(out_field,".mat"));
        if isstruct(out.(out_field))
             out2 = out.(out_field);
             save(out_filepath, '-struct', 'out2');
        else;save(out_filepath, '-struct', 'out','out_field');
        end
    end
end

% 3) Save unfiltered allele data
out_set = "NoFilter";
out_path = fullfile(paths.output,out_set);
if ~exist(out_path,'dir'); mkdir(out_path); end

out = data.(out_set);
for out_field = string(fieldnames(out))'
    out_filepath = fullfile(out_path, append(out_field,".mat"));
    if isstruct(out.(out_field))
         out2 = out.(out_field);
         save(out_filepath, '-struct', 'out2');
    else;save(out_filepath, '-struct', 'out','out_field');
    end
end

clear -regexp ^out




