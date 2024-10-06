%% SCRIPT_3_Ref_and_NonRef_counts
% Script loads base count variables, calculates base frequencies and identified reference and
% non-reference reads
% 
% INPUTS from Output_1:
% > Basecalls_All.mat 
% > Basecalls_For.mat
% > Basecalls_Rev.mat
% > SampleIDs.mat
% > Positions.mat
% > MT_consensus_sequences.mat
%
% OUTPUTS saved into Output_2:
% > Sample IDs.mat 
% > Positions.mat 
% > Basecalls.mat 
% > Reads.mat 
% > Frequency.mat 
% > Reads_Ref.mat 
% > Reads_NonRef.mat 
% > Frequency_Ref.mat 
% > Frequency_NonRef.mat 
% > Excluded_sites.mat 

%% Load data
% Set input paths
paths = struct();
paths.root = 'path/to/main/directory';
pahs.output = fullfile(paths.root,"Output_2");

in_paths = struct();
in_dir = 'Output_1';

in_paths.Basecalls_All = fullfile(paths.root,in_dir,'Basecalls_All.mat');
in_paths.Basecalls_For = fullfile(paths.root,in_dir,'Basecalls_For.mat');
in_paths.Basecalls_Rev = fullfile(paths.root,in_dir,'Basecalls_Rev.mat');
in_paths.SampleIDs     = fullfile(paths.root,in_dir,"SampleIDs.mat");
in_paths.Positions     = fullfile(paths.root,in_dir,"Positions.mat");
in_paths.Consensus     = fullfile(paths.root,in_dir,"MT_consensus_sequences.mat");

% Create data struct
if ~exist('data', 'var'); data = struct(); end

% Load variables
data.Basecalls.All = load(in_paths.Basecalls_All, 'A', 'C', 'G', 'T');
data.Basecalls.For = load(in_paths.Basecalls_For, 'A', 'C', 'G', 'T');
data.Basecalls.Rev = load(in_paths.Basecalls_Rev, 'A', 'C', 'G', 'T');
data.SampleIDs     = load(in_paths.SampleIDs, 'SampleIDs');
data.Positions     = load(in_paths.Positions, 'Positions');
data.Consensus     = load(in_paths.Consensus,'ConcensusSeqs');

% Clear input variables
clear -regexp ^in


%% Set parameters
if ~exist('parameters', 'var'); parameters = struct(); end

% Whch dimension is individuals vs positions
parameters.dim_samples = 1;   % rows = samples
parameters.dim_positions = 2; % columns = positions

% Number of samples and positions
parameters.n_positions = size(data.Positions, parameters.dim_positions);
parameters.n_samples   = size(data.SampleIDs, parameters.dim_samples);
parameters.dimensions = [parameters.n_samples parameters.n_positions];

% rCRS deletion artefact site
parameters.deletion_site = 3107;

% Minimum required directional read depth
parameters.read_depth_min = 500;


%% OPTS: settings for the run
opts = struct();
opts.filter_depth = true;       % Remove positions with read depth below set minimum (500)?
opts.filter_nonACGT_ref = true; % Remove positions where consensus sequences have 'N', deletions, etc?
opts.save = false;              % Save outputs to MAF file? If FALSE, outputs only placed into workspace DATA struct
opts.write_txt = false;         % Write output tables to TXT file?

%% Calculate sum A+C+G+T read totals
%%{
% Read depths
for temp_dir = ["For", "Rev", "All"]
    % Calculate total number of basecalls made (excluding Ns)
    temp_reads  = data.Basecalls.(temp_dir).A ...
                + data.Basecalls.(temp_dir).C ...
                + data.Basecalls.(temp_dir).G ...
                + data.Basecalls.(temp_dir).T;
    data.Reads.(temp_dir) = temp_reads;
end
% Clear temp
clear -regexp ^temp
%}

%% Calculate base frequencies (fraction of A+C+G+T totals, ignoring Ns and indels
for temp_dir = ["For", "Rev", "All"]
    for temp_base = ["A", "C", "G", "T"]
        % Get base call counts
        temp_reads_base  = data.Basecalls.(temp_dir).(temp_base); % Calls of this base
        temp_reads_total = data.Reads.(temp_dir);            % All basecalls
        temp_reads_total(temp_reads_total==0) = NaN;         % Put NaN where there are 0 total reads

        % Calculate base frequency (NaNs in places where there are no reads at all)
        data.Frequency.(temp_dir).(temp_base) = temp_reads_base ./ temp_reads_total;
    end
end
% Clear temp 
clear -regexp ^temp

%% Assemble reference and non-reference read matrices
% Initialise empty arrays for reference and non-reference read counts
temp_ref_reads = struct();
temp_nonref_reads = struct();
for temp_dir = ["All","For","Rev"]
    temp_ref_reads.(temp_dir)    = zeros(parameters.dimensions);
    temp_nonref_reads.(temp_dir) = zeros(parameters.dimensions);
end

% Make 4 masks, one for each base type where it the reference base in each sample
temp_ref_seqs = data.Consensus;
temp_ref_mask = struct();
for temp_base = ['A','C','G','T']
    % Get index of where each base is the reference
    temp_ref_mask.(temp_base) = temp_ref_seqs == temp_base;
end

% Add A reads to ref array in places where A is reference
for temp_base = ['A','C','G','T']
    % Get locations where this base is reference vs non-reference
    temp_base_ref_idx = temp_ref_mask.(temp_base);
    temp_base_nonref_idx = ~temp_ref_mask.(temp_base);

    for temp_dir = ["All","For","Rev"]
        % get read array of this base type
        temp_base_reads = data.Basecalls.(temp_dir).(temp_base);

        % Make empty arrays of ref and nonref reads only
        temp_base_ref_reads    = zeros(parameters.dimensions);
        temp_base_nonref_reads = zeros(parameters.dimensions);

        % Place reads in ref and nonref arrays according to ref and nonref locations
        temp_base_ref_reads(temp_base_ref_idx)       = temp_base_reads(temp_base_ref_idx);
        temp_base_nonref_reads(temp_base_nonref_idx) = temp_base_reads(temp_base_nonref_idx);

        % Add the base-specific ref and nponref reads to the overall arrays
        temp_ref_reads.(temp_dir)    = temp_ref_reads.(temp_dir)   +temp_base_ref_reads;
        temp_nonref_reads.(temp_dir) = temp_nonref_reads.(temp_dir)+temp_base_nonref_reads;
    end
end

% Save ref and nonref counts into the Data struct
data.Reads_Ref    = temp_ref_reads;
data.Reads_NonRef = temp_nonref_reads;

% Clear temp
clear -regexp ^temp_

%% Calculate reference and non-reference read frequencies (as fraction of all ACGT reads)
temp_freq_ref = struct();
temp_freq_nonref = struct();

%get idx of nonACGT ref positions
temp_mask_nonACGT = (data.Consensus=='A'|data.Consensus=='C'|data.Consensus=='G'|data.Consensus=='T');

for temp_dir = ["All","For","Rev"]
    % Calculate reference and nonreference frequncies
    temp_freq_ref.(temp_dir)    = data.Reads_Ref.(temp_dir)    ./ data.Reads.(temp_dir);
    temp_freq_nonref.(temp_dir) = data.Reads_NonRef.(temp_dir) ./ data.Reads.(temp_dir);

    % place NaN in locations with nonACGT reference
    temp_freq_ref.(temp_dir)(temp_mask_nonACGT)    = NaN;
    temp_freq_nonref.(temp_dir)(temp_mask_nonACGT) = NaN;
end

% Save ref and nonref frequencies into Data struct
data.Frequency_Ref = temp_freq_ref;
data.Frequency_NonRef = temp_freq_nonref;

% clear temp
clear -regexp ^temp_

%% Apply filter to exclude data from position 3107, sites with low depth & or non-ACGT consensus base
% Make a mask for sites to exclude
temp_mask_to_exclude = false(parameters.dimensions);

% Filter position 3107
temp_mask_to_exclude(:,3107) = true; 

% Filter low depth sites
if opts.filter_depth
    temp_idx_low_depth = data.Reads.For < parameters.read_depth_min ...
                       | data.Reads.Rev < parameters.read_depth_min;

    temp_mask_to_exclude(temp_idx_low_depth) = true;
end

% Filter sites with nonACGT in consensus sequences
if opts.filter_nonACGT
    temp_idx_nonACGT_reference = ~(data.Consensus=='A' ...
                                 | data.Consensus=='C' ...
                                 | data.Consensus=='G' ...
                                 | data.Consensus=='T');
    temp_mask_to_exclude(temp_idx_nonACGT_reference) = true;
end

% Use mask to replace excluded site values with NaN
for temp_dir = ["All","For","Rev"]
    for temp_base = ['A','C','G','T']
        data.Basecalls.(temp_dir).(temp_base)(temp_mask_to_exclude) = NaN;
        data.Frequency.(temp_dir).(temp_base)(temp_mask_to_exclude) = NaN;
    end
    data.Reads.(temp_dir)(temp_mask_to_exclude)            = NaN;
    data.Reads_Ref.(temp_dir)(temp_mask_to_exclude)        = NaN;
    data.Reads_NonRef.(temp_dir)(temp_mask_to_exclude)     = NaN;
    data.Frequency_Ref.(temp_dir)(temp_mask_to_exclude)    = NaN;
    data.Frequency_NonRef.(temp_dir)(temp_mask_to_exclude) = NaN;
end

% Add the excluded position mask to data struct
data.Excluded_sites = temp_mask_to_exclude;

% clear temp
clear -regexp ^temp_

%% Save outputs to MAT file
if opts.save
    % Make outout directory
    out_dir = paths.output;
    if ~exist(out_dir,'dir'); mkdir(out_dir); end
    
    % Save read count and frequency variables
    out_fields = ["Basecalls","Reads","Frequency","Reads_Ref","Reads_NonRef","Frequency_Ref","Frequency_NonRef"];
    for out_field = out_fields
        out = data.(out_field);
        out_filepath = fullfile(out_dir,append(out_field,".mat"));
        save(out_filepath, '-struct', 'out');
    end
    
    % Save excluded site mask, sample IDs and positions
    out_fields = ["Excluded_sites","SampleIDs","Positions"];
    for out_field = out_fields
        out_filepath = fullfile(out_dir,append(out_field,".mat"));
        save(out_filepath, '-struct', 'data', out_field);
    end
    
    % No need to re-save consensus sequences as these have not been changed in any way
    clear -regexp ^out_
end

%% Write outputs to TXT file
if opts.save_txt
    out_dir = fullfile(paths.output,"txt");
    if ~exist(out_dir,'dir'); mkdir(out_dir); end

    out_varnames = string(1:16569);
    out_samples = array2table(data.SampleIDs,'VariableNames',"Sample");
    
    for out_dir = ["All","For","Rev"]
        for out_field = ["Basecalls","Frequency"]
            for out_base = ['A','C','G','T']
                out_filename = append(out_field,"_",out_dir,"_",out_base,".txt");
                out_filepath = fullfile(out_dir,out_filename);
                out_array = data.(out_field).(out_dir).(out_base);

                out_table = [out_samples, array2table(out_array, 'VariableNames',out_varnames)];
                writetable(out_table, out_filepath, 'Delimiter','\t');
            end
        end

        for out_field = ["Reads","Reads_Ref","Reads_NonRef","Frequency_Ref","Frequency_NonRef"]
            out_filename = append(out_field,"_",out_dir,".txt");
            out_filepath = fullfile(out_dir,out_filename);
            out_array = data.(out_field).(out_dir);

            out_table = [out_samples, array2table(out_array, 'VariableNames',out_varnames)];
            writetable(out_table, out_filepath, 'Delimiter','\t');
        end
    end

    out_field = "Excluded_sites";
    out_filename = append(out_field,".txt");
    out_filepath = fullfile(out_dir,out_filename);
    out_array = data.(out_field);
    out_table = [out_samples, array2table(out_array, 'VariableNames',out_varnames)];
    writetable(out_table, out_filepath, 'Delimiter','\t');

    clear -regexp ^out_
end


