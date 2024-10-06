%% SCRIPT_4_Allele_sorting_by_raw_frequency.m
% Script for sorting alleles by total raw frequencies > STANDARD METHOD
%
% Inputs from Output_1
%   > SampleIDs.mat
%   > Individuals.mat
%   > Positions.mat
% 
% Inputs from Output_2
%   > Reference.mat
%   > rCRS.mat
% 
% Inputs from Output_3
%   > Reads.mat
%   > Basecalls.mat
%   > Frequencies.mat
%
% OUTPUTS saved in Output_4/Standard_method
% 1) Sorted allele data:
%   > Allele_count.mat      - Number of alleles detected at each position
%   > Allele_order.mat      - Index of each base type in the sorted alele frequency order
%   > Allele_bases.mat      - Sorted allele bases (alleles sorted by total reads across both directions)
%   > Allele_reads.mat      - Sorted allele reads (alleles sorted by total reads across both directions)
%   > Allele_frequencies.mat  - Sorted allele frequencies (STANDARD METHOD)
%   > Reads.mat             - Total reads across all alleles
%
% 2) Metadata:
%   > Parameters     - Parameters for re-use in other scripts
%   > Positions      - Position info: Position number and rCRS base
%   > Samples        - Sample info: Sample number and G1K string ID
%   > Sequences      - Sample consensus sequences + sample & position IDs

%% STEP 1: Load data
%STEP 1a: Load variables
% Set paths
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
in_paths.Basecalls     = fullfile(paths.root,in_dir3,'Basecalls.mat');
in_paths.Frequencies   = fullfile(paths.root,in_dir3,'Frequencies.mat');

% Create data struct
if ~exist('data', 'var'); data = struct(); end

% Load variables
data.SampleIDs    = load(in_paths.SampleIDs,   'SampleIDs');
data.Individuals  = load(in_paths.Individuals, 'NumberIDs');
data.Positions    = load(in_paths.Positions,   'Positions');
data.Reference    = load(in_paths.Reference,   'ConcensusSeqs');
data.rCRS         = load(in_paths.rCRS,        'rCRS');
data.Reads        = load(in_paths.Reads,       'For', 'Rev', 'All');
data.Basecalls    = load(in_paths.Basecalls,   'For', 'Rev', 'All');
data.Frequencies  = load(in_paths.Frequencies, 'For', 'Rev', 'All');

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
parameters.dimensions = [parameters.n_positions parameters.n_individuals];

% Sample and position metadata
parameters.Individuals = data.Individuals;
parameters.SampleIDs   = data.SampleIDs;
parameters.Positions   = data.Positions;
parameters.rCRS_bases  = data.rCRS;

% Base, allele and direction strings
parameters.directions = ["All", "For", "Rev"];
parameters.bases      = ['A', 'C', 'G', 'T'];
parameters.alleles    = ["Major", "Minor_1", "Minor_2", "Minor_3"];

%% STEP 2: Sort all alleles by total (F+R) number of reads
temp_dims = parameters.dimensions;
temp_dirs = ["All","For","Rev"];
temp_bases = ['A','C','G','T'];
temp_alleles = ["Major","Minor_1","Minor_2","Minor_3"];

% Initiate result structs
temp_alt_num = struct();
temp_alt_bases = struct();
temp_alt_counts = struct();

temp_alt_base_order = struct();
temp_alt_base_order.A = zeros(temp_dims);
temp_alt_base_order.C = zeros(temp_dims);
temp_alt_base_order.G = zeros(temp_dims);
temp_alt_base_order.T = zeros(temp_dims);

for temp_dir = temp_dirs
    switch temp_dir
        case "All" % Sort 4 alleles in order by TOTAL read count
            % Combine base total count matrices as sheets into 3D matrix
            temp_all        = data.Basecalls.(temp_dir).A;
            temp_all(:,:,2) = data.Basecalls.(temp_dir).C;
            temp_all(:,:,3) = data.Basecalls.(temp_dir).G;
            temp_all(:,:,4) = data.Basecalls.(temp_dir).T;
            
            % Get number of alternatove alleles at each position
            temp_alt_num = sum(temp_all>0,3);
                            
            % Find 1st, 2nd, 3rd and 4th max read count (as sheets in 3D matrix) and its index
            [temp_alt_reads_sorted, temp_alt_sort_idx] = maxk(temp_all,4,3);

            % Make 3D matrix of bases arranged in decreasing read count order along dim3
            temp_alt_bases_sorted = temp_bases(temp_alt_sort_idx);
            
            % Get sorted allele order position for each base type               
            for temp_allele = 1:4
                for temp_base = temp_bases
                    temp_base_idx = temp_alt_bases_sorted(:,:,temp_allele)==temp_base;            
                    temp_alt_base_order.(temp_base)(temp_base_idx) = temp_allele;
                end
            end

            % Find indices of 0 read count values
            temp_idx_zeros = temp_alt_reads_sorted==0;
            
            % Replace all bases with 0 reads with N
            temp_alt_bases_sorted(temp_idx_zeros) = 'N';

            % Split the sorted alternative allele read count 3D matrix  into 4 per-allele 2D matrices
            temp_alt_counts.Major.(temp_dir) = temp_alt_reads_sorted(:,:,1);
            temp_alt_counts.Minor_1.(temp_dir) = temp_alt_reads_sorted(:,:,2);
            temp_alt_counts.Minor_2.(temp_dir) = temp_alt_reads_sorted(:,:,3);
            temp_alt_counts.Minor_3.(temp_dir) = temp_alt_reads_sorted(:,:,4);

            % Split the sorted alternative allele base 3D matrix into 4 per-allele 2D matrices
            temp_alt_bases.Major = temp_alt_bases_sorted(:,:,1);
            temp_alt_bases.Minor_1 = temp_alt_bases_sorted(:,:,2);
            temp_alt_bases.Minor_2 = temp_alt_bases_sorted(:,:,3);
            temp_alt_bases.Minor_3 = temp_alt_bases_sorted(:,:,4);
                        
        otherwise % Rearrange sorted alelle direction specific read counts 
            for temp_allele = temp_alleles                    
                
                % Make a matrix for allele X direction specific read counts
                temp_allele_reads = zeros(temp_dims);

                for temp_base = temp_bases
                    % positions where this base is allele number X
                    temp_base_idx = temp_alt_bases.(temp_allele)==temp_base;

                    % Extract base specific directional reads and add to allele read matrix
                    temp_base_reads = data.Basecalls.(temp_dir).(temp_base);
                    temp_allele_reads(temp_base_idx) = temp_base_reads(temp_base_idx);  

                end
                % Save directional read counts to struct
                temp_alt_counts.(temp_allele).(temp_dir) = temp_allele_reads;
            end
    end
end

% Save to data struct
data.Sorted_alleles.Allele_count = temp_alt_num;
data.Sorted_alleles.Allele_reads = temp_alt_counts;
data.Sorted_alleles.Allele_bases = temp_alt_bases;
data.Sorted_alleles.Allele_order = temp_alt_base_order;

% Clear temp vars
clear -regexp ^temp

%% STEP 3: Rearrange base frequencies in the sorted allele order
temp_dirs = ["All","For","Rev"];
temp_bases = ['A','C','G','T'];
temp_alleles = ["Major","Minor_1","Minor_2","Minor_3"];

% output dta acontainer
temp_freq_sorted = struct();

for temp_dir = temp_dirs
    for temp_allele = temp_alleles
        % Get correct order number for each sorted allele
        switch temp_allele
            case "Major";   temp_allele_idx = 1;
            case "Minor_1"; temp_allele_idx = 2;
            case "Minor_2"; temp_allele_idx = 3;
            case "Minor_3"; temp_allele_idx = 4;
        end
        % empty comtainer for frequency values
        temp_allele_freq = zeros(parameters.dimensions);

        % Take requred frequency values from each base type and save into specific sorted allele  
        for temp_base = temp_bases
            temp_idx = data.Sorted_alleles.Allele_order.(temp_base) == temp_allele_idx;
            temp_allele_freq(temp_idx) = data.Frequencies.(temp_dir).(temp_base)(temp_idx);
        end
        temp_freq_sorted.(temp_allele).(temp_dir) = temp_allele_freq;
    end
end
% Place into DATA struct
data.Sorted_alleles.Allele_frequencies = temp_freq_sorted;

% Clear temp
clear -regexp ^temp

%% STEP 4: Save all required outputs for the next steps to file in Output_4/Standard_method
out_dir = "Standard_method";
out_path = fullfile(paths.output,out_dir);
if ~exist(out_path, 'dir'); mkdir(out_path); end

% 4.1  Save sorted allele variables
out_fields = [...
    "Allele_count",...
    "Allele_order",...
    "Allele_bases",...
    "Allele_reads",...
    "Allele_frequencies"];
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
out_field = ["Reads", "Excess_reads"];
out_filepath = fullfile(out_path, append(out_field,".mat"));
out = data.(out_field);
save(out_filepath, '-struct', 'out')


% 4.3  Save Individuals, SampleIDs, Positions, rCRS, Reference, Parameters
out_files = ["Samples", "Positions", "Sequences", "Parameters"];
for out_file = out_files
    out_filepath = fullfile(out_path, append(out_file, ".mat"));
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

