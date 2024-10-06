%% SCRIPT 1: basecount_import
% Script for reading .txt files with raw basecount tables generated after parsing pileups
% with outout dimensions samples = rows and positions = columns. Script also removes rows 
% with data associated with samples "HG00102" and "NA19031" (due to partially missing 
% alignment data). All data is tehn saved as Matlab variables to .MAT file
%
% OUTPUTS:
% > Basecalls_All.mat
% > Basecalls_For.mat
% > Basecalls_Rev.mat
% > Indels.mat
% > Reads_TOTAL.mat
% > SampleIDs.mat
% > NumberIDs.mat
% > Positions.mat
%

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

    % Load each table, trim, format. Dims: samples = rows, positions = columns
    temp_table = readtable(temp_filepath, 'delimiter', '\t', 'ReadVariableNames', true);
    
    if temp_file == "Sample_ID"
        % Get sample ID strings
        data.SampleIDs = string(temp_table{:,:})); 
        
        % Find row numbers for data from the truncated samples "HG00102" and "NA19031" 
        temp_rows_to_remove = sort([...
            find(data.SampleIDs == "NA19031") ...
            find(data.SampleIDs == "HG00102")]);

        % Remove samples "HG00102" and "NA19031" from SampleIDs list
        data.SampleIDs(temp_rows_to_remove)=[];

        % Make sample number and position number variables
        data.NumberIDs = (1:length(data.SampleIDs))'; % row numbers
        data.Positions = 1:16569; % column numbers

    else
        % Get basecount matrix only, samples = rows, positions = columns
        temp_datamatrix = temp_table{1:end, 3:end};  

        % Remove rows associated with samples "HG00102" and "NA19031" 
        temp_datamatrix(temp_rows_to_remove,:) = [];
        switch temp_file
            case "All_As"; data.Basecalls.All.A = temp_datamatrix;
            case "All_Cs"; data.Basecalls.All.C = temp_datamatrix;
            case "All_Gs"; data.Basecalls.All.G = temp_datamatrix;
            case "All_Ts"; data.Basecalls.All.T = temp_datamatrix;
            case "All_Ns"; data.Basecalls.All.N = temp_datamatrix;
            case "For_As"; data.Basecalls.For.A = temp_datamatrix;
            case "For_Cs"; data.Basecalls.For.C = temp_datamatrix;
            case "For_Gs"; data.Basecalls.For.G = temp_datamatrix;
            case "For_Ts"; data.Basecalls.For.T = temp_datamatrix;
            case "For_Ns"; data.Basecalls.For.N = temp_datamatrix;
            case "Rev_As"; data.Basecalls.For.A = temp_datamatrix;
            case "Rev_Cs"; data.Basecalls.For.C = temp_datamatrix;
            case "Rev_Gs"; data.Basecalls.For.G = temp_datamatrix;
            case "Rev_Ts"; data.Basecalls.For.T = temp_datamatrix;
            case "Rev_Ns"; data.Basecalls.For.N = temp_datamatrix;
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

out_field = "Basecalls";
for out_field2 = ["All", "For", "Rev"]
    out = data.(out_field).(out_field2);
    out_filepath = fullfile(out_dir, append(out_field,"_"out_field2,".mat"));
    save(out_filepath, '-struct', 'out')
end

% Save indels
out_field = "Indels";
out = data.(out_field);
out_filepath = fullfile(out_dir, append(out_field,".mat"));
save(out_filepath, '-struct', 'out')

% Save total reads, sample IDS, numbers and position numbers
for out_field = ["Reads_TOTAL","SampleIDs","NumberIDs","Positions"]
    out = struct();
    out.(out_field) = data.(out_field);
    out_filepath = fullfile(out_dir, append(out_field,".mat"));
    save(out_filepath, '-struct', 'out')
end
clear -regexp ^out
