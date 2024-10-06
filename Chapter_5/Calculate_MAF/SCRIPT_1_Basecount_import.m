%% SCRIPT 1: basecount_import
% Script for reading .txt files with raw basecount tables generated after parsing pileups.
% Each table is transposed so that samples = columns and positions = rows, then saved to .MAT file

% Set paths
paths = struct();
paths.root   = '/path/to/main/directory';
paths.output = fullfile(paths.root,append('Output_1'));
paths.input  = "/path/to/input/directory";

% Create data struct
data = struct(); 

% Set input directory and file names
temp_files = ["Sample_IDs", "Reads",...
    "All_As", "All_Cs", "All_Gs", "All_Ts", "All_Ns",...
    "For_As", "For_Cs", "For_Gs", "For_Ts", "For_Ns",...
    "Rev_As", "Rev_Cs", "Rev_Gs", "Rev_Ts", "Rev_Ns",...
    "Insertions", "Deletions"];

% Load txt files
for temp_file = temp_files
    temp_filepath = fullfile(paths.input, append(temp_file,".txt"));

    % Load each table, trim, format & transpose so that samples = columns, positions = rows
    temp_table = readtable(temp_filepath, 'delimiter', '\t', 'ReadVariableNames', true);
    
    if temp_file == "Sample_ID"
        temp_sampleIDs = string(temp_table{:,:})); % Get sample ID strings
        data.SampleIDs = temp_sampleIDs'; % transpose to columns
        data.NumberIDs = 1:length(temp_sampleIDs); % columns
        data.Positions = (1:16569)'; % rows
    else
        temp_datamatrix = temp_table{1:end, 3:end}; % basecount matrix only
        temp_datamatrix = temp_datamatrix'; % transpose so that samples = columns
        switch temp_file
            case "All_As"; data.Basecalls_ALL.A = temp_datamatrix;
            case "All_Cs"; data.Basecalls_ALL.C = temp_datamatrix;
            case "All_Gs"; data.Basecalls_ALL.G = temp_datamatrix;
            case "All_Ts"; data.Basecalls_ALL.T = temp_datamatrix;
            case "All_Ns"; data.Basecalls_ALL.N = temp_datamatrix;
            case "For_As"; data.Basecalls_FOR.A = temp_datamatrix;
            case "For_Cs"; data.Basecalls_FOR.C = temp_datamatrix;
            case "For_Gs"; data.Basecalls_FOR.G = temp_datamatrix;
            case "For_Ts"; data.Basecalls_FOR.T = temp_datamatrix;
            case "For_Ns"; data.Basecalls_FOR.N = temp_datamatrix;
            case "Rev_As"; data.Basecalls_REV.A = temp_datamatrix;
            case "Rev_Cs"; data.Basecalls_REV.C = temp_datamatrix;
            case "Rev_Gs"; data.Basecalls_REV.G = temp_datamatrix;
            case "Rev_Ts"; data.Basecalls_REV.T = temp_datamatrix;
            case "Rev_Ns"; data.Basecalls_REV.N = temp_datamatrix;
            case "Insertions"; data.Indels.Ins = temp_datamatrix;
            case "Deletions";  data.Indels.Del = temp_datamatrix;
            case "Reads"; data.Reads_TOTAL = temp_datamatrix;
        end
    end
end
% Clear temp variables
clear -regexp ^temp 

% Save basecount matrices as MAT files
out_dir = paths.output;
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

for out_field = ["Basecalls_ALL", "Basecalls_FOR", "Basecalls_REV", "Indels"]
    out = data.(out_field);
    out_filepath = fullfile(out_dir, append(out_field,".mat"));
    save(out_filepath, '-struct', 'out')
end

for out_field = ["Reads_TOTAL","SampleIDs","NumberIDs","Positions"]
    out = struct();
    out.(out_field) = data.(out_field);
    out_filepath = fullfile(out_dir, append(out_field,".mat"));
    save(out_filepath, '-struct', 'out')
end
clear -regexp ^out
