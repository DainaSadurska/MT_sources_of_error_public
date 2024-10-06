%% SCRIPT 2: Consensus seq import
% Script for importing 3202 sample consensus sequences from FASTA file 
% Input file: HiCov_3202_MT_consensus_sequences.fa
% Output files in directory Output_2: 
% > G1K_HiCov_Consensus_sequences.mat 
% > Mask_del.mat
% > Mask_ref.mat

%% STEP 1: Load rCRS and sample consensus sequences and save as MAT variables

% Set paths
paths = struct();
paths.root = '/path/to/main/directory';
paths.input = "/path/to/input/directory";
paths.output = fullfile(paths.root, 'Output_2');

% Create data struct
data = struct();

% Read rCRS sequence from FASTA file
in_filename = "rCRS.fa";
in_filepath = fullfile(paths.input,in_filename);
temp_rCRS = fastaread(in_filepath);
temp_rCRS = string({temp_rCRS.Sequence}.');
data.rCRS = (char(temp_rCRS))';
clear -regexp ^in

% Read consensus sequences from FASTA file
in_filename = "HiCov_3202_MT_consensus_sequences_20220520.fa";
in_filepath = fullfile(paths.input,in_filename);
temp_FASTAData = fastaread(in_filepath);
clear -regexp ^in

% Parse FASTAData struct
temp_n_inds = numel(temp_FASTAData);
temp_ind_nums = (1:temp_n_inds)';
temp_IDs  = string({temp_FASTAData.Header}.');
temp_seqs = string({temp_FASTAData.Sequence}.');
temp_seqs_char = char(temp_seqs);

% Modify consensus sequences
% 1) Replace deletions marked as "*" as "-"
temp_seqs_char(temp_seqs_char == '*') = '-';

% 2) Set 3107N as '-' in all samples
temp_seqs_char(:,3107) = '-';

% 3) Set all to upper case
temp_seqs_char(temp_seqs_char == 'a') = 'A';
temp_seqs_char(temp_seqs_char == 'c') = 'C';
temp_seqs_char(temp_seqs_char == 'g') = 'G';
temp_seqs_char(temp_seqs_char == 't') = 'T';
temp_seqs_char(temp_seqs_char == 'n') = 'N';

% 4) Replace any non characters other than [ACGTN-] with N
temp_seqs_char(...
    temp_seqs_char ~= 'A'& ...
    temp_seqs_char ~= 'C'& ...
    temp_seqs_char ~= 'G'& ...
    temp_seqs_char ~= 'T'& ...
    temp_seqs_char ~= '-'& ...
    temp_seqs_char ~= 'N') = 'N';
temp_seqs = string(temp_seqs_char);

% Save sample IDs and consensus sequences in data struct - TRANSPOSED
data.Positions = (1:16569)'; % positions = rows
data.SampleIDs = temp_IDs';  % samples = columns
data.ConsensusSeqs = (temp_seqs_char)'; % TRANSPOSED so that samples = columns

% Create 3202 HiCov sample overview table with IDs and consensus sequences
temp_table = struct();
temp_table.Sample = (1:length(temp_IDs))';
temp_table.G1K_ID = temp_IDs;
temp_table.ConsensusSeq = temp_seqs;
data.Sample_overview = struct2table(temp_table);
clear -regexp ^temp

% Make outout directory
if ~exist(paths.output,'dir'); mkdir(paths.output); end

% Save rCRS to MAT file
out_filepath = fullfile(paths.output, "rCRS.mat");
save(out_filepath, '-struct', 'data', 'rCRS');

% Save consensus sequences into MAT file
out_filepath = fullfile(paths.output, "G1K_HiCov_Consensus_sequences.mat");
out = rmfield(data, 'rCRS');
save(out_filepath, '-struct', 'out');
clear -regexp ^out

%% STEP 2: Make sample consensus sequence mask variables
%) Base type mask: locations where consensus sequences have each type of base
data.Mask_ref = struct();
for temp_base = ['A','C','G','T']
    data.Mask_ref.(temp_base) = data.Reference==temp_base;
end

% 2) Deletion and 'N' site mask 
temp_mask_del = (data.Reference=='N' | data.Reference=='-');
data.Mask_del = temp_mask_del;
clear -regexp ^temp

% Save reference mask to MAT file
out_filepath = fullfile(paths.output, "Mask_ref.mat");
out = data.Mask_ref;
save(out_filepath, '-struct', 'out')

% Save deletion / N site mask to MAT file
out_filepath = fullfile(paths.output,"Mask_del.mat");
save(out_filepath, '-struct', 'data', 'Mask_del');
clear -regexp ^out


