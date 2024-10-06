%% SCRIPT_2_Motif_location_finder
% Script for identifying all motif-of-interest locations in sample consensus sequences. Outputs 
% contain lists of mtDNA positions denoting the 3'-most base of each motif of interest instance.
%
% INPUTS:
%  > Motif_spaces            - Lists with the motifs to be searched for, motif space IDs: N2-N8
%  > MT_consensus_sequences  - Sample consensus sequences
%  > SampleIDs               - List of sample IDs
%
% OUTPUTS for each motif space N2 to N8:
%  > [MotifSpaceID]_Motif_locations.mat        
%  > [MotifSpaceID]_Motif_counts.mat
%
%  > [MotifSpaceID]_Motif_counts_L_strand.txt
%  > [MotifSpaceID]_Motif_counts_H_strand.txt
%  > [MotifSpaceID]_Motif_counts_total.txt
%  > [MotifSpaceID]_Motif_table.txt
%   
%
%% SET PATHS
% Paths: outputs
paths = struct();
paths.root = '/path/to/main/directory/';
paths.output = fullfile(paths.root,"Output");
if ~exist(paths.output,'dir'); mkdir(paths.output); end

% Paths: inputs
paths.inputs = '/path/to/input/directory/';

%% LOAD DATA
% Load data
in_dir = paths.inputs;
in_paths = struct();
in_paths.Motif_spaces = fullfile(in_dir, "Motif_spaces.mat");
in_paths.Sequences    = fullfile(in_dir, "MT_consensus_sequences.mat");    
in_paths.SampleIDs    = fullfile(in_dir, "SampleIDs.mat");      

data = struct();
data.Motif_spaces       = load(in_paths.Motif_spaces,'N2','N3','N4','N5','N6','N7','N8');
data.SampleIDs          = load(in_paths.SampleIDs,'SampleIDs');
data.Sequences          = load(in_paths.Sequences,'ConsensusSeqs');

clear -regexp ^in_


%% Get dataset parameters
parameters = struct();

% All motif space IDs
parameters.motif_spaces = string(fieldnames(data.Motif_spaces))';% ["N2", "N3", "N4", "N5", "N6", "N7", "N8"];

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

%% OPTS - set options for this run
opts = struct();
opts.motif_spaces = ["N2", "N3", "N4"]; % Which motif spaces to run (must be subset of parameters.motif_spaces)
opts.save = false;  % Save outputs to file?  If FALSE, outputs only deposited into DATA struct


%% STEP 1 - Get sample IDs and reference sequences with position 3107 removed
% Get position numbers for positions to be included only
temp_col_removed = parameters.del_artefact_site;
temp_cols = parameters.positions;
temp_cols(temp_col_removed)=[];

% Get sample reference sequences and generate their complements 
temp_seqs = data.Sequences(:,temp_cols);
temp_seqs_compl = repmat('N',size(temp_seqs));
for temp_row = 1:size(temp_seqs,1)
    temp_seqs_compl(temp_row,:) = seqcomplement(temp_seqs(temp_row,:));
end

% Save reference sequences for use in scripts
in_refseqs_L = temp_seqs;
in_refseqs_H = temp_seqs_compl; % has not been flipped yet

% Included positions only for use in the script
in_positions = parameters.positions(:,temp_cols);

% Get sample ID list
in_sampleIDs = data.SampleIDs;

% clear temp
clear -regexp ^temp

%% STEP 2 - Loop through motif spaces and find locations of each motif in every sample
for in_motif_space = opts.motif_spaces
    %% 2.1 - Get the motif list
    disp(append("Starting motif space '",in_motif_space))
    
    loop_motif_list = data.Motif_spaces.(in_motif_space);
    [loop_n_motifs, loop_motif_size] = size(loop_motif_list);

    % Start global timer
    time_start = tic;
    time_check = 25; % How frequently time remaining is extimated
    
    
    %% 2.2 - Extend reference sequences to handle motifs crossing 16569/1 split
    % Get included position numbers list
    loop_positions = in_positions;    % excludes 3107!!!
    
    % Set sequence and position list extension
    temp_ext_posnums = 1:(loop_motif_size-1);
    temp_ext_refseqs_L = in_refseqs_L(:,1:(loop_motif_size-1));
    temp_ext_refseqs_H = in_refseqs_H(:,1:(loop_motif_size-1));
    
    % Make extended positions list
    loop_posnums_extended_L = [loop_positions, temp_ext_posnums];
    loop_posnums_extended_H = [loop_positions, temp_ext_posnums];
    
    % Make extended ref sequences to cover 1/16569 split
    loop_refseqs_extended_L = [in_refseqs_L, temp_ext_refseqs_L];
    loop_refseqs_extended_H = [in_refseqs_H, temp_ext_refseqs_H];
    
    % Flip H strand sequence and position list
    loop_posnums_extended_H = flip(loop_posnums_extended_H,2);
    loop_refseqs_extended_H = flip(loop_refseqs_extended_H,2);
    
    % clear temp
    clear -regexp ^temp
    
    %% 2.3 - Make output variables for this motif space
    loop_n_samples = parameters.n_samples; %2533
    
    % 1) Table with the motif list
    out_motif_table = struct();
    out_motif_table.Motif_ID = (1:loop_n_motifs)';     % Motif numbers
    out_motif_table.Motif = string(loop_motif_list);   % Motif sequences
    out_motif_table.Found              = false(loop_n_motifs,1); % Placeholder values
    out_motif_table.Instances_total    = zeros(loop_n_motifs,1); % Placeholder values
    out_motif_table.Instances_L_strand = zeros(loop_n_motifs,1); % Placeholder values
    out_motif_table.Instances_H_strand = zeros(loop_n_motifs,1); % Placeholder values
    out_motif_table.Found_in_samples   = zeros(loop_n_motifs,1); % Placeholder values
    out_motif_table = struct2table(out_motif_table);
    
    % 2) Tables for L, H and total motif counts in each sample
    temp_table = array2table(zeros(loop_n_motifs,loop_n_samples));
    temp_table.Properties.VariableNames = in_sampleIDs;

    out_motif_counts_L     = [out_motif_table(:,1:2), temp_table];
    out_motif_counts_H     = [out_motif_table(:,1:2), temp_table];
    out_motif_counts_total = [out_motif_table(:,1:2), temp_table];

    % Sctruct for the motif hit locations in each samples
    out_locations = struct();
    
    %% 2.4 - Iterate through motifs and find each motif locations in every samples
    for loop_motif = 1:loop_n_motifs
        % Start loop timer
        time_loop_start = tic;
    
        % get motif sequence
        loop_motif_seq = loop_motif_list(loop_motif,:);
    
        % make motif hit output containers: 
        loop_motif_locations_L = {[loop_n_samples,1]}; % Position numbers for all hits
        loop_motif_locations_H = {[loop_n_samples,1]};
    
        loop_motif_hit_count_per_sample_L = zeros(loop_n_samples,1); % number of hits found in each sample
        loop_motif_hit_count_per_sample_H = zeros(loop_n_samples,1);
        loop_motif_hit_count_per_sample = zeros(loop_n_samples,1);
    
        %disp(append("Motif ",string(loop_motif),": ", loop_motif_seq))
        fprintf('Motifs %d: %s',loop_motif,loop_motif_seq);
        
        % Initialise counters
        loop_counter_hits_total = 0;
        loop_counter_hits_L_total = 0;
        loop_counter_hits_H_total = 0;
        loop_counter_motifs_found_in_samples = 0;
    
        % Iterate through samples and find all motif locations
        for temp_ind = 1:loop_n_samples
    
            % Get a single reference sequence
            temp_refseq_L = loop_refseqs_extended_L(temp_ind,:);
            temp_refseq_H = loop_refseqs_extended_H(temp_ind,:);
        
            temp_posnums_L = loop_posnums_extended_L;
            temp_posnums_H = loop_posnums_extended_H;
    
            % Search L sequence for motif instances - STRFIND returns the index of leftmost letter of each match
            temp_hit_idx_L_START = strfind(temp_refseq_L,loop_motif_seq);
            temp_hit_idx_L_END = temp_hit_idx_L_START + (loop_motif_size - 1);
    
            % Search H sequence for motif instances 
            temp_hit_idx_H_START = strfind(temp_refseq_H,loop_motif_seq);
            temp_hit_idx_H_END = temp_hit_idx_H_START + (loop_motif_size - 1);
    
            % Get each motif instance LAST base locations - i.e. position with error rate
            temp_hit_positions_L = temp_posnums_L(temp_hit_idx_L_END);
            temp_hit_positions_H = temp_posnums_H(temp_hit_idx_H_END);
    
            % Save this motif location data into motif space output containers
            loop_motif_locations_L{temp_ind,1} = temp_hit_positions_L; % Position numbers for all hits
            loop_motif_locations_H{temp_ind,1} = temp_hit_positions_H;
    
            temp_hit_count_L = numel(temp_hit_positions_L);
            temp_hit_count_H = numel(temp_hit_positions_H);
    
            loop_motif_hit_count_per_sample_L(temp_ind,1) = temp_hit_count_L; % number of hits
            loop_motif_hit_count_per_sample_H(temp_ind,1) = temp_hit_count_H;
            loop_motif_hit_count_per_sample(temp_ind,1) = temp_hit_count_L+temp_hit_count_H;
    
            % Update motif count totals
            loop_counter_hits_L_total = loop_counter_hits_L_total+temp_hit_count_L;
            loop_counter_hits_H_total = loop_counter_hits_H_total+temp_hit_count_H;
            loop_counter_hits_total   = loop_counter_hits_total+temp_hit_count_L+temp_hit_count_H;
    
            % Increment temp_found_in_samples counter
            if temp_hit_count_L+temp_hit_count_H >0
                loop_counter_motifs_found_in_samples = loop_counter_motifs_found_in_samples+1;
            end        
        end
    
        % Save outputs for the current motif ONLY if motif was found at all
        if loop_counter_hits_total > 0
            % Update Found status
            out_motif_table.Found(loop_motif,1) = true;
            out_motif_table.Found_in_samples(loop_motif,1) = loop_counter_motifs_found_in_samples;
    
            % Save total instances counters
            out_motif_table.Instances_total(loop_motif,1)    = loop_counter_hits_total;
            out_motif_table.Instances_L_strand(loop_motif,1) = loop_counter_hits_L_total;
            out_motif_table.Instances_H_strand(loop_motif,1) = loop_counter_hits_H_total;
            
            % Save motif locations
            out_locations.Motif_locations_L.(loop_motif_seq) = loop_motif_locations_L;
            out_locations.Motif_locations_H.(loop_motif_seq) = loop_motif_locations_H;
    
            % Save Motif counts in each sample
            out_motif_counts_L{loop_motif,3:end} = loop_motif_hit_count_per_sample_L';
            out_motif_counts_H{loop_motif,3:end} = loop_motif_hit_count_per_sample_H';
            out_motif_counts_total{loop_motif,3:end} = loop_motif_hit_count_per_sample';
        end
    
        % Display loop time elapsed per motif
        temp_time_iteration = seconds(toc(time_loop_start));
        temp_time_iteration.Format = 's'; %set duration format as seconds
        fprintf('    %s\n', string(temp_time_iteration))
        
        % Estimate time remaining every 25 loops
        if mod(loop_motif,time_check) == 0
            % get total time elapsed
            temp_iteration = loop_motif;
            temp_n_iterations = loop_n_motifs;
            temp_time_elapsed = seconds(toc(time_start));
            temp_time_elapsed.Format = 'dd:hh:mm:ss'; %set duration format
    
            % get average iteration length
            temp_time_avg_iteration = temp_time_elapsed./temp_iteration;
            temp_time_avg_iteration.Format = 's'; %set duration format
    
            % get total remaining time (seconds)
            temp_time_remaining = ((temp_n_iterations-temp_iteration)*temp_time_avg_iteration);
            temp_time_remaining.Format = 'dd:hh:mm:ss';   %set duration format
    
            % get % of loops remaining
            temp_iterations_remaining = ((temp_n_iterations-temp_iteration)./temp_n_iterations*100);
    
            % Display % of motif space remaining
            fprintf(' > %0.2f%% of motifs remaining\n > Avg time per motif %s\n',...
                temp_iterations_remaining, string(temp_time_avg_iteration));
    
            % Display time elapsed and remaining
            fprintf(" > Time elapsed: %s, Estimated time remaining: %s \n", temp_time_elapsed, temp_time_remaining)  
        end
    
        % clear temp
        clear -regexp ^temp
    end

    %% 2.5 - SAVE motif location variables to file
    if opts.save == true
        disp(append("Saving locations for motif space ",in_motif_space))
    
        % output path
        out_path = paths.output;
        if ~exist(out_path,'dir'),mkdir(out_path),end
        
        % Set MAT variable file names
        temp_filepath_locations = fullfile(out_path,append(in_motif_space,'_Motif_locations.mat'));
        temp_filepath_counts = fullfile(out_path,append(in_motif_space,'_Motif_counts.mat'));

        % Save locations variables
        out_locations.Motif_table = out_motif_table;
        save(temp_filepath_locations, '-struct', 'out_locations');
        
        % Save counts variables
        out_counts = struct();
        out_counts.Motif_table = out_motif_table;
        out_counts.Motif_counts.L_strand = out_motif_counts_L;
        out_counts.Motif_counts.H_strand = out_motif_counts_H;
        out_counts.Motif_counts.Total = out_motif_counts_total;
        save(temp_filepath_counts, '-struct', 'out_counts');
        
        % Write count tables to txt file
        temp_filepath_counts_L = fullfile(out_path,append(in_motif_space,'_Motif_counts_L_strand.txt'));
        temp_filepath_counts_H = fullfile(out_path,append(in_motif_space,'_Motif_counts_H_strand.txt'));
        temp_filepath_counts_total = fullfile(out_path,append(in_motif_space,'_Motif_counts_total.txt'));
        temp_filepath_motif_table = fullfile(out_path,append(in_motif_space,'_Motif_table.txt'));
        
        writetable(out_motif_table,temp_filepath_motif_table,'Delimiter','\t');
        writetable(out_motif_counts_L,temp_filepath_counts_L,'Delimiter','\t')
        writetable(out_motif_counts_H,temp_filepath_counts_H,'Delimiter','\t')
        writetable(out_motif_counts_total,temp_filepath_counts_total,'Delimiter','\t')
        
        % clear temp
        clear -regexp ^temp
    end

    %% 2.6 - Put loop outputs into data struct and clear variables
    % save outputs into DATA struct
    data.(in_motif_space).Motif_table = out_motif_table;
    data.(in_motif_space).Motif_counts.L_strand = out_motif_counts_L;
    data.(in_motif_space).Motif_counts.H_strand = out_motif_counts_H;
    data.(in_motif_space).Motif_counts.Total    = out_motif_counts_total;
    data.(in_motif_space).Motif_locations.L_strand = out_locations.Motif_locations_L;
    data.(in_motif_space).Motif_locations.H_strand = out_locations.Motif_locations_H;


    % Clear higher level LOOP variables
    clear -regexp ^out
    clear -regexp ^loop
    clear -regexp ^temp
    clear -regexp ^time

    disp(append("Motif space '",in_motif_space,"' completed"))
end
% Clear input variables
clear -regexp ^in


%% FOR TESTING & DESCRIPTION
%{
% TEST SEQUENCE: 
% position   1  2  3  4  5  6  7  8  9  0 11 12 13 14 15 16 17 18 19 20 21    
% sequence   T  T  A  A  T  T  C  C  C  C  A  A  -  T  T  C  A  A  T  A  A
%
% TEST MOTIF:
% motif      A  A  T 

% Dummy motif and sequence
test_motif_seq = 'AAT';
test_motif_size = length(test_motif_seq);

test_refseq_orig = 'TTAATTCCCCAA-TTCAATAA';
test_posnums_orig = 1:length(test_refseq_orig);

% Remove position '-' - replicates situation with 3107
test_pos_to_remove = strfind(test_refseq_orig,'-');
test_refseq = test_refseq_orig;
test_refseq(test_pos_to_remove)=[];
test_posnums = test_posnums_orig;
test_posnums(test_pos_to_remove)=[];

% Extend refseq to capture motifs crossing over the split
test_ext_seq = test_refseq(1:(test_motif_size-1));
test_ext_pos = 1:(test_motif_size-1); 

% Make extended L strand sequence and position list
test_refseq_L = [test_refseq,  test_ext_seq]; % add extensions after last base
test_posnums_L= [test_posnums, test_ext_pos]; 

% Get extended H strand sequence by flipping complementary sequence to extended L seq
test_refseq_H = flip(seqcomplement(test_refseq_L),2);
test_posnums_H = flip(test_posnums_L,2); % also flip position list

%
% MODIFIED SEQUENCE: extended & pos 13 removed
% pos        1  2  3  4  5  6  7  8  9  0 11 12 14 15 16 17 18 19 20 21     1  2   
% H seq: <<< A  A  T  T  A  A  G  G  G  G  T  T  A  A  G  T  T  A  T  T  |  A  A
% L seq: >>> T  T  A  A  T  T  C  C  C  C  A  A  T  T  C  A  A  T  A  A  |  T  T
%
% MOTIF OF INTEREST
% on L strand:  >>>  A  A  T 
% on H strand:  <<<  T  A  A
%
%
% MOTIFS PRESENT:
% H hits              <------                 <------                 <---------
% H seq: <<< A  A  T  T  A  A  G  G  G  G  T  T  A  A  G  T  T  A  T  T  |  A  A
%            1  2  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 20 21     1  2    
% L seq: >>> T  T  A  A  T  T  C  C  C  C  A  A  T  T  C  A  A  T  A  A  |  T  T
% L hits           ------>                 ------>        ------>  --------->            
%
%
% FINDING MOTIFS IN L STRAND
% HITS:            ------>                 ------>        ------>  ------|-->
% L seq: >>> T  T  A  A  T  T  C  C  C  C  A  A  T  T  C  A  A  T  A  A  |  T  T
% pos L      1  2  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 20 21  |  1  2   
% idx        1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20  | 21 22
% hit idx          3-----5                11----13       16----18 19-----|-21    
% hit pos          3-----5                11----14       17----19 20-----|--1    
%
% L idx START: 3,11,16,19  < motif 5' end index given by strfind()
% L idx END:   5,13,18,21  < 5' end index + motif legth - 1 
%
% L pos START: 3,11,17,20  < Actual position number of motif start base at 3'
% L pos END:   5,14,19,1   < Actual position number of motif end base at 5' end
%
%
% FINDING MOTIFS IN L STRAND
% HITS:      ------|-->                 ------>                 ------>
% H flip >>> A  A  |  T  T  A  T  T  G  A  A  T  T  G  G  G  G  A  A  T  T  A  A
% pos H      2  1  | 21 20 19 18 17 16 15 14 12 11 10  9  8  7  6  5  4  3  2  1
% idx        1  2  |  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
% hit idx    1-----|--3                 9----11                17----19   
% hit pos    2-----|-21                15----12                 6-----4   
%
% H idx START: 1,9,17    < output from strfind()
% H idx END:   3,11,19   < idx END = idx STRAT + (motif length - 1)
%
% H pos START: 2,15,6    < use IDX to extract position from position numbers list 
% H pos END:   21,12,4
%
%
% FINAL OUTPUT - POSITION NUMBERS OF MOTIF 3' BASE IN THE ORIGINAL SEQUENCE
% L hit locations: 5,14,19,1
% H hit locations: 21,12,4
%
       
% Search L sequence for motif instances - STRFIND returns the index of leftmost letter of each match
test_hit_idx_L_START = strfind(test_refseq_L,test_motif_seq);
test_hit_idx_L_END = test_hit_idx_L_START + (test_motif_size - 1);

% Search H sequence for motif instances 
test_hit_idx_H_START = strfind(test_refseq_H,test_motif_seq);
test_hit_idx_H_END = test_hit_idx_H_START + (test_motif_size - 1);

% Get each motif instance LAST base locations - i.e. position with error rate
test_hit_locations_L = test_posnums_L(test_hit_idx_L_END);
test_hit_locations_H = test_posnums_H(test_hit_idx_H_END);

disp("Expected H locations: 21,12,4;  Calculated: ")
disp(test_hit_locations_H)

disp("Expected L locations: 5,14,19,1;  Calculated: ")
disp(test_hit_locations_L)

%disp(test_refseq_L(test_hit_idx_L_START(1):test_hit_idx_L_END(1)))

%clear -regexp ^test_
%}

