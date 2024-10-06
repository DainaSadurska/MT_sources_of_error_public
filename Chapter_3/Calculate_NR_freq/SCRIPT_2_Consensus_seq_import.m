%% SCRIPT 2: Consensus seq import
% Script for importing 2535 sample consensus sequences from FASTA file, reordering to match  
% the order of SampleIDs list, and removing any for sequences of samples not in the SampleIDs list.
% As SampleIDs list does not contains samples HG00102 and NA19031, sequences for these samples are 
% not included in the output MAT files
%
% Input FASTA files: 
% > G1K_phase3_MT_consensus_sequences.fa
% > rCRS.fa
%
% Input MAF file:
% > SampleIDs.mat
%
% Output files saved in directory Output_1: 
% > MT_consensus_sequences.mat 
% > rCRS.mat 


%% STEP 1: Load rCRS and sample consensus sequences

% Set paths
paths = struct();
paths.root  = '/path/to/main/directory';
paths.input = '/path/to/input/directory';
paths.output = fullfile(paths.root, 'Output_1');

% Create data struct
if ~exist('data','var');data = struct();end

% Read rCRS sequence from FASTA file
in_filename = "rCRS.fa";
in_filepath = fullfile(paths.input,in_filename);
temp_rCRS = fastaread(in_filepath);
temp_rCRS = string({temp_rCRS.Sequence}.');
data.rCRS = (char(temp_rCRS));
clear -regexp ^in

% Read consensus sequences from FASTA file
in_filename = "G1K_phase3_MT_consensus_sequences.fa";
in_filepath = fullfile(paths.input,in_filename);
temp_FASTAData = fastaread(in_filepath);
clear -regexp ^in

% Parse FASTAData struct
temp_n_inds = numel(temp_FASTAData);
temp_ind_nums = (1:temp_n_inds)';
temp_IDs  = string({temp_FASTAData.Header}.');
temp_seqs = string({temp_FASTAData.Sequence}.');
temp_seqs_char = char(temp_seqs);

%% Apply consensus sequence modifications
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

% 4) Replace any characters other than [ACGTN-] with N
temp_seqs_char(...
    temp_seqs_char ~= 'A'& ...
    temp_seqs_char ~= 'C'& ...
    temp_seqs_char ~= 'G'& ...
    temp_seqs_char ~= 'T'& ...
    temp_seqs_char ~= '-'& ...
    temp_seqs_char ~= 'N') = 'N';
temp_seqs = string(temp_seqs_char);

%% Sort consensus sequences to be in the same order as Sample IDs from BAM processing
% Load sampleIDs list produced in SCRIPT_1_Basecount_import
in_dir = paths.output; % from Output_1
in_filename = "SampleIDs.mat";
in = load(fullfile(in_dir,in_filename),'SampleIDs');
data.SampleIDs = in.SampleIDs;
clear -regexp ^in_

% Make an empty array fopr sorted consensus seqeunces
temp_seqs_sorted = strings(size(data.SampleIDs));

% Pass through sample IDs and grab the corresponding consensus sequence
for temp_sample = 1:numel(data.SampleIDs)
    temp_ID_str = data.SampleIDs(temp_sample);

    % Find which sample in consensus sequences file is the one I am looking for
    temp_idx = find(temp_IDs == temp_ID_str);

    %throw error if no match found. If no match, check IDs from FASTA for leading/trailing spaces
    if isempty(temp_idx); error(append("No match for Sample ",temp_ID_str)); end 

    % Grab that sample consensus sequence and put into the correct row of the sorted seqs array
    temp_seqs_sorted(temp_sample) = temp_seqs(temp_idx);
end

% Save sample IDs and sorted consensus sequences into data struct
data.Positions = 1:16569;   % positions = columns
data.SampleIDs = temp_IDs;  % samples = rows
data.ConsensusSeqs = char(temp_seqs_sorted); 

% Make a sample overview table with IDs and consensus sequences
temp_table = struct();
temp_table.Sample = (1:numel(data.SampleIDs))';
temp_table.G1K_ID = data.SampleIDs;
temp_table.ConsensusSeq = temp_seqs_sorted;
data.Sample_overview = struct2table(temp_table);
clear -regexp ^temp

%% Save sample consensus sequences as MAT variable
% Make outout directory
if ~exist(paths.output,'dir'); mkdir(paths.output); end

% Save rCRS to its own individual MAT file called "rCRS.mat")
out_filepath = fullfile(paths.output, "rCRS.mat");
save(out_filepath, '-struct', 'data', 'rCRS');

% Save all DATA struct fields as variables into MAT file called "MT_consensus_sequences"
out_filepath = fullfile(paths.output, "MT_consensus_sequences.mat");
save(out_filepath, '-struct', 'data');

%{
% Save consensus sequences to their own individual MAT file
out_filepath = fullfile(paths.output, "MT_consensus_sequences.mat");
out = rmfield(data, 'rCRS');
save(out_filepath, '-struct', 'out');
%}
clear -regexp ^out


