%% SCRIPT_6_Allele_sorting_by_consensus_frequency.m
% Script for sorting alleles by consensus read frequencies > CONSENSUS METHOD
% 
% INPUTS
% from Output_1
%   > SampleIDs.mat
%   > Individuals.mat
%   > Positions.mat
%
% from Output_2
%   > Reference.mat
%   > rCRS.mat
%
% from Output_3
%   > Consensus_reads.mat
%   > Excess_reads.mat
%   > Consensus_basecalls.mat
%   > Excess_basecalls.mat
%   > Consensus_frequency.mat
%
% OUTPUTS saved in Output_4/Consensus_method
% 1) Sorted allele data:
%   > Allele_count.mat                 - Number of alleles detected at each position
%   > Allele_order.mat                 - Index for sorting 4 bases in order by allele consensus frequency
%   > Allele_bases.mat                 - Sorted allele bases
%   > Allele_consensus_reads.mat       - Sorted allele consensus reads 
%   > Allele_excess_reads.mat          - Sorted allele excess reads 
%   > Allele_consensus_frequencies.mat - Sorted allele consensus frequencies 
%   > Consensus_reads.mat              - Total consensus reads across all alleles
%   > Excess_reads.mat                 - Total excess reads across all alleles
%   
% 2) Metadata:
%   > Parameters     - Parameters for re-use in other scripts
%   > Positions      - Position info: Position number and rCRS base
%   > Samples        - Sample info: Sample number and G1K string ID
%   > Sequences      - Sample consensus sequences + sample & position IDs

%% STEP 1: Paths and inputs
%STEP 1a: Load variables
paths = struct();
paths.root = 'path/to/main/directory';
paths.output = fullfile(paths.root,'Output_4');
if ~exist(paths.output, 'dir'); mkdir(paths.output); end

% Set input paths
in_paths = struct();
in_dir1 = 'Output_1';
in_dir2 = 'Output_2';
in_dir3 = 'Output_3';

% from Output_1
in_paths.SampleIDs     = fullfile(paths.root,in_dir1,"SampleIDs.mat");
in_paths.Individuals   = fullfile(paths.root,in_dir1,"NumberIDs.mat");
in_paths.Positions     = fullfile(paths.root,in_dir1,"Positions.mat");

% from Output_2
in_paths.Reference     = fullfile(paths.root,in_dir2,"G1K_HiCov_Consensus_sequences.mat");
in_paths.rCRS          = fullfile(paths.root,in_dir2,"rCRS.mat");

% from Output_3
in_paths.Consensus_reads       = fullfile(paths.root,in_dir3,'Consensus_reads.mat');
in_paths.Excess_reads          = fullfile(paths.root,in_dir3,'Excess_reads.mat');
in_paths.Consensus_basecalls   = fullfile(paths.root,in_dir3,'Consensus_basecalls.mat');
in_paths.Excess_basecalls      = fullfile(paths.root,in_dir3,'Excess_basecalls.mat');
in_paths.Consensus_frequencies = fullfile(paths.root,in_dir3,'Consensus_frequencies.mat');

% Create data struct
if ~exist('data', 'var'); data = struct(); end

% Load variables
data.SampleIDs    = load(in_paths.SampleIDs,  'SampleIDs');
data.Individuals  = load(in_paths.Individuals,'NumberIDs');
data.Positions    = load(in_paths.Positions,  'Positions');
data.Reference    = load(in_paths.Reference,  'ConcensusSeqs');    
data.rCRS         = load(in_paths.rCRS,       'rCRS');    
data.Consensus_reads       = load(in_paths.Consensus_reads,    'For', 'Rev', 'All');
data.Excess_reads          = load(in_paths.Excess_reads,       'For', 'Rev', 'All');
data.Consensus_basecalls   = load(in_paths.Consensus_basecalls,'For', 'Rev', 'All');
data.Excess_basecalls      = load(in_paths.Excess_basecalls,   'For', 'Rev', 'All');
data.Consensus_frequencies = load(in_paths.Consensus_frequencies,'A','C','G','T');

% Clear input variables
clear -regexp ^in
            
% STEP 1b: Set parameters
if ~exist('parameters', 'var'); parameters = struct(); end

% Which dimension is individuals vs positions
parameters.dim_positions = 1;   % Rows
parameters.dim_individuals = 2; % Columns

% Number of samples and positions
parameters.n_positions   = size(data.Positions,   parameters.dim_positions);
parameters.n_individuals = size(data.Individuals, parameters.dim_individuals);
parameters.dimensions    = [parameters.n_positions parameters.n_individuals];

% Sample and position metadata
parameters.Individuals = data.Individuals;
parameters.SampleIDs   = data.SampleIDs;
parameters.Positions   = data.Positions;
parameters.rCRS_bases  = data.rCRS;

% Base, allele and direction strings
parameters.directions = ["All", "For", "Rev"];
parameters.bases      = ['A', 'C', 'G', 'T'];
parameters.alleles    = ["Major", "Minor_1", "Minor_2", "Minor_3"];



%% STEP 2: Sort alleles by total (F+R) number of consensus reads
temp_bases = parameters.bases;     % ['A', 'C', 'G', 'T'];
temp_alleles = parameters.alleles; % ["Major", "Minor_1", "Minor_2", "Minor_3"];

% Variables for outputs
temp_allele_count = struct();
temp_allele_bases = struct();
temp_allele_order = struct();
temp_allele_con_reads = struct();

% Use consensus read count across BOTH directions for sorting alleles
for temp_dir = parameters.directions
    switch temp_dir
        case "All" %% Allele order determined by consensus read count across BOTH directions
            % Combine base-specific consensus read count 2D matrices into single 3D matrix
            temp_reads_3D        = data.Consensus_basecalls.(temp_dir).(temp_bases(1)); % A
            temp_reads_3D(:,:,2) = data.Consensus_basecalls.(temp_dir).(temp_bases(2)); % C
            temp_reads_3D(:,:,3) = data.Consensus_basecalls.(temp_dir).(temp_bases(3)); % G
            temp_reads_3D(:,:,4) = data.Consensus_basecalls.(temp_dir).(temp_bases(4)); % T
            
            % Get number of alternatove alleles at each position
            temp_allele_count = sum(temp_reads_3D>0,3);
                            
            % Sort alleles by consensus read count & get sorted index
            [temp_reads_sorted, temp_alt_sort_idx] = maxk(temp_reads_3D,4,3);

            % Make 3D matrix of bases arranged in decreasing read count order along dim3
            temp_bases_sorted = temp_bases(temp_alt_sort_idx);
            
            % Base sort order: find locations where each base type is 1st, 2ns, 3rd or 4th allele  
            for temp_base = temp_bases
                temp_allele_order.(temp_base) = zeros(parameters.dimensions); 
                for temp_allele_num = 1:numel(temp_alleles)
                    temp_base_idx = temp_bases_sorted(:,:,temp_allele_num)==temp_base; 
                    temp_allele_order.(temp_base)(temp_base_idx) = temp_allele_num;
                end
            end
            
            % Find indices of 0 read count values and replace all bases with 0 reads with N
            temp_idx_zeros = temp_reads_sorted == 0;
            temp_bases_sorted(temp_idx_zeros) = 'N';
            
            % Slice sorted allele and read count 3D matrices into 4 2D slices
            for temp_allele_num = 1:numel(temp_alleles)
                temp_allele = temp_alleles(temp_allele_num);
                temp_allele_bases.(temp_allele) = temp_bases_sorted(:,:,temp_allele_num);
                temp_allele_cons_reads.(temp_allele).(temp_dir) = temp_reads_sorted(:,:,temp_allele_num);
            end
            
        otherwise %% FOR and REV consensus read count data rearranged in the same order 
            for temp_allele = temp_alleles                                                            
                temp_reads = zeros(temp_dims);
                
                for temp_base = temp_bases
                    % Identify positions where this base is allele number X
                    temp_base_idx = temp_allele_bases.(temp_allele)==temp_base;

                    % Extract base specific directional reads and add to allele read matrix
                    temp_base_reads = data.Consensus_basecalls.(temp_dir).(temp_base);
                    temp_reads(temp_base_idx) = temp_base_reads(temp_base_idx);   
                end
                % Save directional read counts to struct
                temp_allele_cons_reads.(temp_allele).(temp_dir) = temp_reads;
            end
    end
end

% Save variables to DATA struct
data.Sorted_alleles.Allele_count = temp_allele_count;
data.Sorted_alleles.Allele_order = temp_allele_order;
data.Sorted_alleles.Allele_bases = temp_allele_bases;
data.Sorted_alleles.Allele_consensus_reads = temp_allele_cons_reads;

% Clear temp
clear -regexp ^temp


%% STEP 3: Rearrange excess reads and concensus frequencies into the now sorted allele order
% Initiate result structs
temp_allele_cons_freqs = struct();
temp_allele_excs_reads = struct();

% Iterate through sorted alleles
for temp_allele = parameters.alleles

    % Make empty comtainers for sorted allele consensus frequency and excess read data
    temp_allele_cons_freqs.(temp_allele) = zeros(parameters.dimensions);
    for temp_dir = parameters.directions
        temp_allele_excs_reads.(temp_allele).(temp_dir) = zeros(parameters.dimensions);
    end

    % Get sorted allele number (1 = major alelle, 1-3 = minor alleles in decr. frequency order)
    switch temp_allele
        case "Major";   temp_allele_num = 1;
        case "Minor_1"; temp_allele_num = 2;
        case "Minor_2"; temp_allele_num = 3;
        case "Minor_3"; temp_allele_num = 4;
    end

    % Iterate through base types & rearrange base-specific variables into the sorted allele order
    for temp_base = parameters.bases

        % Get locations where the specific allele (1st, 2nd, 3rd or 4th) is the base of interest
        temp_idx = data.Sorted_alleles.Allele_order.(temp_base) == temp_allele_num;

        % Use location idx to extract & write to output the relevant base-specific consensus frequencies
        temp_allele_cons_freqs.(temp_allele)(temp_idx) = data.Consensus_frequencies.(temp_base)(temp_idx);

        % Use location idx to extract & write to output the relevant base-specific excess reads (directional)
        for temp_dir = parameters.directions
            temp_allele_excs_reads.(temp_allele).(temp_dir)(temp_idx) = data.Excess_basecalls.(temp_dir).(temp_base)(temp_idx);
        end
    end
end

% Save outputs to DATA struct
data.Sorted_alleles.Allele_excess_reads = temp_allele_excs_reads;
data.Sorted_alleles.Allele_consensus_frequencies = temp_sorted_consensus_freqs;

% Clear temp
clear -regexp ^temp


%% STEP 4: Save outputs to MAT file
out_dir = "Consensus_method";
out_path = fullfile(paths.output,out_dir);
if ~exist(out_path, 'dir'); mkdir(out_path); end

% 4.1  Save sorted allele variables
out_fields = [...
    "Allele_count",...
    "Allele_order",...
    "Allele_bases",...
    "Allele_consensus_reads",...
    "Allele_excess_reads",...
    "Allele_consensus_frequencies"];

out_str = data.Sorted_alleles;
for out_field = out_fields
    out_filepath = fullfile(out_path, append(out_field,".mat"));
    if isstruct(out_str.(out_field))
         out = out_str.(out_field);
         save(out_filepath, '-struct', 'out')
    else;save(out_filepath, '-struct', 'out_str', out_field);
    end
end

% 4.2 Save sum read count data 
out_fields = ["Consensus_reads", "Excess_reads"];
for out_field = out_fields
    out_filepath = fullfile(out_path, append(out_field,".mat"));
    out = data.(out_field);
    save(out_filepath, '-struct', 'out')
end

% 4.3  Save Individuals, SampleIDs, Positions, rCRS, Reference, Parameters 
out_files = ["Samples", "Positions", "Sequences", "Parameters"];
for out_file = out_files
    out_filepath = fullfile(out_path, append(out_file,".mat"));
    out = struct();
    switch out_file
        case "Samples"   
            % File for sample metadata
            out.Individuals = data.Individuals;
            out.SampleIDs   = data.SampleIDs;
        case "Positions" 
            % File for position metadata
            out.Positions   = data.Positions;
            out.rCRS        = data.rCRS;
        case "Sequences" 
            % File for consensus sequences
            out.Individuals = data.Individuals;
            out.SampleIDs   = data.SampleIDs;
            out.Positions   = data.Positions;
            out.rCRS        = data.rCRS;
            out.ConsensusSeqs = data.Reference;
        case "Parameters"
            out = parameters;
    end
    save(out_filepath, '-struct', 'out')
end
clear -regexp ^out

