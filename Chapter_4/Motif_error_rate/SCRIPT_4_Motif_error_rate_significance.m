%% SCRIPT_4_Motif_error_rate_significance
% Script for statistical significance testing of observed motif-assciated sequencing error rates 
% using Fisher's exact and Chi squared tests. 
%
% INPUTS for each motif space N2 to N8:     < produced by SCRIPT_3_Motif_error_rate_calculator
%  > [MotifSpaceID]_Motif_error_rates.mat   - Population-level and Sample-level motif error rate (ERD) tables
%  > [MotifSpaceID]_Readcount_totals.mat    - Arrays with sample-level FM/FMM/RM/RMM count totals
%
% OUTPUTS for each motif space N2 to N8:  
%  > [MotifSpaceID]_MER_stats_summary.mat         
%  > [MotifSpaceID]_MER_Sample_level_stats.mat - Test results using sample-level contingency tables
%  > [MotifSpaceID]_MER_Global_stats.mat       - Test results using population-level contingency tables
%
% Tables saved as tab-deloimited text files:
%  > [MotifSpaceID]_MER_stats_summary.txt
%  > [MotifSpaceID]Sample_MER_significance_Result.txt
%  > [MotifSpaceID]Sample_MER_significance_P_value.txt
%  > [MotifSpaceID]Sample_MER_significance_Test_type.txt
%  > [MotifSpaceID]Population_MER_significance_test_results.txt
%
% Notation: MER - Motif Error Rate, equivalent to ERD (Error Rate Difference)
%{
% Notes:
%
% 1) Statistical test choice: 
% Fishers exact test (left-tailed) is used preferentially, except where contingency table values 
% exceed 10^7. Due to limitation of FISHERTEST function, Chi squared test is used for contingency
% tables above 10^7 instead.
%
% 2) Skipped motifs:
% To reduce computation time script skips over motifs where by definition the result can not be 
% significant. Testing is skipped in two types of scenarios: 
%  - where no instances of the motif were analysed and therefore FM/FMM/RM/RMM counts are 0
%  - where ERD is negative (i.e. reverse mismatch rate exceeds motif forward mismatch rate) and
%    therefore by definition motif error rate is indistingushable from 0. 
% In these cases test is skipped and the following default values written into outputs: 
% Motif not found: test type = "n/a",  h = NaN,  p = NaN
% ERD is negative: test type = "skip", h = 0,    p = NaN 
%
% 3) Correction for multiple comparions:
% Bonferroni correction is applied, adjusting sample level and population level significance 
% threhsolds by a factor of number_of_samples x number_of_motifs for sample level tests and a factor 
% of number_of_motifs for population level tests. 
%
%
% Function used for performong Chi square test on contingency tables is based on CHISQUARECONT by
% Panagiotis Moulos, 2008.
%}


%% Set paths & Load data
% Paths: outputs
paths = struct();
paths.root = '/path/to/main/directory/';
paths.output = fullfile(paths.root,"Output_3");
if ~exist(paths.output,'dir'); mkdir(paths.output); end

% Paths: inputs
paths.input1 = '/path/to/input/directory/';     % For the list of samples
paths.input2 = fullfile(paths.root,"Output_2"); % For error rate and readcount totals matrices

in_paths = struct();

% Assemble filepaths
in_dir = paths.input1;
in_paths.Samples      = fullfile(in_dir, "Samples.mat"); 
in_paths.Motif_spaces = fullfile(in_dir, "Motif_spaces.mat");

in_dir = paths.input2;
in_motif_spaces = ['N2','N3','N4','N5','N6','N7','N8'];
for in_motif_space = in_motif_spaces
    in_paths.(in_motif_space).Motif_error_rates = fullfile(in_dir, append(in_motif_space,"_Motif_error_rates.mat"));
    in_paths.(in_motif_space).Readcount_totals  = fullfile(in_dir, append(in_motif_space,"_Readcount_totals.mat"));
end

% Create DATA struct
if ~exist('data','var'); data = struct(); end

% Load files
data.SampleIDs    = load(in_paths.Samples, 'SampleIDs');
data.Motif_spaces = load(in_paths.Motif_spaces);
for in_motif_space = in_motif_spaces
    data.(in_motif_space).Population_MER   = load(in_paths.(in_motif_space).Motif_error_rates,'Population_MER');
    data.(in_motif_space).Sample_MER       = load(in_paths.(in_motif_space).Motif_error_rates,'Sample_MER');
    data.(in_motif_space).Readcount_totals = load(in_paths.(in_motif_space).Readcount_totals,'FM','FMM','RM','RMM');
end
clear -regexp ^in_

%% Get dataset parameters
parameters = struct();

% Motif space IDs
parameters.motif_spaces = string(fieldnames(data.Motif_spaces))';

% Dataset dimensions
parameters.dim_motifs = 1;  % Rows = motifs
parameters.dim_samples = 2; % Cols = samples 

% Sample count and number IDs
parameters.n_samples  = numel(data.SampleIDs);
parameters.samples    = 1:parameters.n_samples;

% Significance trheshold
parameters.alpha = 0.05;

% maximum contingency table value for fisher test (Matlab cannot handle larger)
parameters.FT_conttab_max = 10000000;

%% OPTS - set options for this run
opts = struct();
opts.motif_spaces = ["N2", "N3", "N4"]; % Which motif spaces to run (must be subset of parameters.motif_spaces)
opts.save_MAT = false;        % Save outputs to file?  If FALSE, outputs only deposited into DATA struct
opts.save_TXT = false;        % Write output tables to tab-delimited txt file? 
opts.time_check_interval = 4; % How frequenty time remaining is estimated

%% MAIN Significance testing of sample-level motif error rates
% Loop through motif spaces, load FM/FMM/RM/RMM arrays & perform significance testing on contingency tables
for in_motif_space = parameters.motif_spaces
    %% Get motif space data
    disp(append("Starting motif sapce ",in_motif_space))
    
    loop_motif_list       = data.(in_motif_space).Population_MER(:,1:3); % Motif number, sequence and analysed status only!
    loop_global_ERD_table = data.(in_motif_space).Population_MER;        % Full motif table
    loop_sample_ERD_table = data.(in_motif_space).Sample_MER; % Rows - motifs, Cols - samples
    
    loop_sampleIDs = string(loop_sample_ERD_table.Properties.VariableNames(3:end));
    loop_n_samples = numel(loop_sampleIDs);
    
    loop_motif_seqs = loop_motif_list.Motif;
    loop_n_motifs   = numel(loop_motif_seqs);
    loop_motif_IDs  = loop_motif_list.Motif_ID;
    loop_idx_analysed_motifs = loop_motif_list.Motif_analysed;
    
    % Input varables: Column vectors of population-level count totals
    loop_global_FM  = loop_global_ERD_table.Sum_FM;
    loop_global_FMM = loop_global_ERD_table.Sum_FMM;
    loop_global_RM  = loop_global_ERD_table.Sum_RM;
    loop_global_RMM = loop_global_ERD_table.Sum_RMM;
    
    % Calculate global ERD and readcount total
    loop_global_ERD = (loop_global_FMM./(loop_global_FMM+loop_global_FM))-(loop_global_RMM./(loop_global_RMM+loop_global_RM));
    loop_global_readcount_total = loop_global_FM + loop_global_FMM + loop_global_RM + loop_global_RMM;
    
    % Input variables: Arrays of sample-level count totals
    loop_sample_FM  = data.(in_motif_space).Readcount_totals.FM;
    loop_sample_FMM = data.(in_motif_space).Readcount_totals.FMM;
    loop_sample_RM  = data.(in_motif_space).Readcount_totals.RM;
    loop_sample_RMM = data.(in_motif_space).Readcount_totals.RMM;
    
    % Calculate sample ERD and readcount total
    loop_sample_ERD = (loop_sample_FMM./(loop_sample_FMM+loop_sample_FM))-(loop_sample_RMM./(loop_sample_RMM+loop_sample_RM));
    loop_sample_readcount_total = loop_sample_FM + loop_sample_FMM + loop_sample_RM + loop_sample_RMM;


    %% Set significance level adjusted for multiple comparisons (Bonferroni correction)
    loop_alpha_sample = parameters.alpha ./(loop_n_motifs .* loop_n_samples);
    loop_alpha_global = parameters.alpha ./ loop_n_motifs;
    
    
    %% Make empty containers for outputs
    % Sample level outputs
    temp_dims = [loop_n_motifs loop_n_samples];
    loop_out_sample_p = NaN(temp_dims);    % P value
    loop_out_sample_h = NaN(temp_dims);    % Test decision
    
    % Population level outputs
    loop_out_global_p = NaN(loop_n_motifs,1);
    loop_out_global_h = NaN(loop_n_motifs,1);
    clear -regexp ^temp_
    
    %% Sample-level data: Set which sample-level cells to skip and which test to run
    %{
    % For performace reasons we can skip cells where result cannot be significant by definition:
    % 1) cells where no motif instances were found in the particular sample (i.e. readcount total = 0)
    % 2) cells where sample ERD is negative (i.e. FER does not exceed RER)
    % When cell is skipped, the following default values are written into the output arrays:
    % Scenario 1 (Motif not found): test type = "n/a",  h = NaN,  p = NaN
    % Scenario 2 (ERD negative):    test type = "skip", h = 0,    p = NaN %Inf
    %}
    % Get the index of cells to skip that cannot be significant by definition
    temp_idx_skip1 = loop_sample_readcount_total == 0;
    temp_idx_skip2 = loop_sample_ERD < 0;
    loop_idx_skip = temp_idx_skip1 | temp_idx_skip2;
    
    % Insert default values into P and H sample-level output arrays for skipped cells
    loop_out_sample_p(loop_idx_skip) = NaN;
    loop_out_sample_h(temp_idx_skip1) = NaN;
    loop_out_sample_h(temp_idx_skip2) = 0;
    
    % Determine whether fisher test can be run at each cell(contingency table values must be <10000000)
    temp_idx_FT = loop_sample_FM  < parameters.FT_conttab_max ...
                & loop_sample_FMM < parameters.FT_conttab_max ...
                & loop_sample_RM  < parameters.FT_conttab_max ...
                & loop_sample_RMM < parameters.FT_conttab_max;
    temp_idx_FT(loop_idx_skip) = false;
    
    % For any contingency tables with values >10^7 chi2 test must be used
    temp_idx_X2 = ~temp_idx_FT;
    temp_idx_X2(loop_idx_skip) = false;
    
    % Sample-level: Make a array with information of which test is run in which cell
    loop_sample_test_type = strings(size(loop_out_sample_p));
    loop_sample_test_type(temp_idx_FT) = "FT";  
    loop_sample_test_type(temp_idx_X2) = "X2";
    loop_sample_test_type(temp_idx_skip1) = "n/a";   % skipped cells due to motif not found
    loop_sample_test_type(temp_idx_skip2) = "skip"; % skipped cells due to ERD <0
    
    loop_sample_analysed = ~temp_idx_skip1; % all cells with some readcounts
    
    %% Population-level data: Set which motifs to skip and which tests to run
    % Get the index of cells to skip that cannot be significant by definition
    temp_global_idx_skip1 = loop_global_readcount_total == 0;
    temp_global_idx_skip2 = loop_global_ERD < 0;
    loop_global_idx_skip = temp_global_idx_skip1 | temp_global_idx_skip2;
    
    % Global: Insert default values into P and H output arrays for skipped cells
    loop_out_global_p(loop_global_idx_skip) = NaN;
    loop_out_global_h(temp_global_idx_skip1) = NaN;
    loop_out_global_h(temp_global_idx_skip2) = 0;
    
    % Sample: Determine whether fisher test can be run at each cell(contingency table values must be <10000000)
    temp_global_idx_FT = loop_global_FM  < parameters.FT_conttab_max ...
                       & loop_global_FMM < parameters.FT_conttab_max ...
                       & loop_global_RM  < parameters.FT_conttab_max ...
                       & loop_global_RMM < parameters.FT_conttab_max;
    temp_global_idx_FT(loop_global_idx_skip) = false;
    
    % Sample: For any contingency tables with values >10^7 chi2 test must be used
    temp_global_idx_X2 = ~temp_global_idx_FT;
    temp_global_idx_X2(loop_global_idx_skip) = false;
    
    % Global: Make an array with information of which test is run for which motif
    loop_global_test_type = strings(size(loop_out_global_p));
    loop_global_test_type(temp_global_idx_FT) = "FT";  
    loop_global_test_type(temp_global_idx_X2) = "X2";
    loop_global_test_type(temp_global_idx_skip1) = "n/a"; % skipped cells due to motif not found
    loop_global_test_type(temp_global_idx_skip2) = "skip";% skipped cells due to ERD <0
    
    loop_global_analysed = ~temp_global_idx_skip1; % all cells with some readcounts
    
    clear -regexp ^temp_
    
    %% Handle motifs that will be fully skipped due to no samples containing them
    % Get lists of motifs analaysed versus those skipped entirely with no instances found in any samples
    loop_motif_list_analysed = loop_motif_list(loop_idx_analysed_motifs,:);
    loop_motif_list_ignored  = loop_motif_list(~loop_idx_analysed_motifs,:);
    
    % List of skipped motif IDs - use to unsert NaNs into corresponding output array rows
    temp_idx_rows_ignored = loop_motif_list_ignored.Motif_ID;
    loop_out_sample_p(temp_idx_rows_ignored,:) = NaN;
    loop_out_sample_h(temp_idx_rows_ignored,:) = NaN;
    
    loop_out_global_p(temp_idx_rows_ignored,:) = NaN;
    loop_out_global_h(temp_idx_rows_ignored,:) = NaN;
    clear -regexp ^temp_
    
    % List of motif IDs that were analysed: loop through THIS set only
    loop_motif_IDs_analysed = loop_motif_list_analysed.Motif_ID;
    
    %% Run significnce testing
    % Start timer
    time_start = tic;
    time_check = opts.time_check_interval;
    
    % Initialise iteration counters
    time_iterations_total = numel(loop_motif_IDs_analysed);
    time_iteration = 0;
       
    % iterate over the analysed motifs (rows)
    for temp_motif = loop_motif_IDs_analysed'%1:loop_n_motifs
        % Start loop timer
        time_loop_start = tic;
        time_iteration = time_iteration+1; %Increment iteration number
    
        % Get motif sequence and number
        temp_motif_seq = loop_motif_list.Motif(temp_motif);
        fprintf('Motif %d: %s  [%d of %d]', temp_motif, temp_motif_seq, time_iteration, time_iterations_total)
    
        %% POPULATION-LEVEL: Run test on the motif population level contingency table
        temp_skip = loop_global_idx_skip(temp_motif);
        if ~temp_skip
            temp_cont_table = [...
                loop_global_FM(temp_motif)   ...
                loop_global_FMM(temp_motif); ...
                loop_global_RM(temp_motif)   ...
                loop_global_RMM(temp_motif)];
        
            switch loop_global_test_type(temp_motif)
                case "FT" % Run fishers test (left tailed)
                    [temp_h,temp_p] = fishertest(temp_cont_table,'Tail','left','Alpha',loop_alpha_global);
        
                case "X2" % Run chi square test
                    [temp_h,temp_p] = chi2cont(temp_cont_table, loop_alpha_global); 
                otherwise
                    error('Test type missing')
            end 
    
            % Save p value and test decision into population-leveloutput arrays
            loop_out_global_p(temp_motif) = temp_p;
            loop_out_global_h(temp_motif) = temp_h;
        end
        
        %% SAMPLE-LEVEL: Iterate over samples (columns) and run test on each sample contingency table
        for temp_sample = 1:loop_n_samples
    
            % get SKIP or NOT SKIP decision
            temp_skip = loop_idx_skip(temp_motif,temp_sample); 
    
            % If cell NOT skipped, run tests
            if ~temp_skip 
                % Assemble contingency table with sample FM/FMM/RM/RMM totals
                temp_cont_table = [...
                    loop_sample_FM(temp_motif,temp_sample)   ...
                    loop_sample_FMM(temp_motif,temp_sample); ...
                    loop_sample_RM(temp_motif,temp_sample)   ...
                    loop_sample_RMM(temp_motif,temp_sample)];
    
                % Get which test to run
                temp_test = loop_sample_test_type(temp_motif,temp_sample);
                switch temp_test
                    case "FT" % Run fishers test (left tailed)
                        [temp_h,temp_p] = fishertest(temp_cont_table,'Tail','left','Alpha',loop_alpha_sample);
    
                    case "X2" % Run chi square test
                        [temp_h,temp_p] = chi2cont(temp_cont_table, loop_alpha_sample); 
                    otherwise
                        error('Test type missing')
                end
    
                % Save p value and test decision into output arrays
                loop_out_sample_p(temp_motif,temp_sample) = temp_p;
                loop_out_sample_h(temp_motif,temp_sample) = temp_h;
            end
            
        end
    
        %% Loop timing
        % Print elapsed time per motif
        temp_time_this_iteration = seconds(toc(time_loop_start));
        temp_time_this_iteration.Format = 's'; %set duration format as seconds
        fprintf('   \t%s', string(temp_time_this_iteration));
    
        % Estimate time remaining every 4 motifs
        if mod(time_iteration,time_check) == 0
            % get total time elapsed
            temp_time_elapsed = seconds(toc(time_start));
            temp_time_elapsed.Format = 'dd:hh:mm:ss'; %set duration format
    
            % get average iteration length
            temp_time_avg_iteration = temp_time_elapsed./time_iteration;
            temp_time_avg_iteration.Format = 's'; %set duration format
    
            % get total remaining time (seconds)
            temp_time_remaining = ((time_iterations_total-time_iteration)*temp_time_avg_iteration);
            temp_time_remaining.Format = 'dd:hh:mm:ss';   %set duration format
    
            % Display time elapsed and remaining
            fprintf("\tTime elapsed: %s, Estimated time remaining: %s", temp_time_elapsed, temp_time_remaining) 
        end;fprintf('\n');
    
        % Clear temp
        clear -regexp ^temp_
            
    end
    %% Save loop outputs into DATA struct
    % Sample level stats
    temp_out_sample = struct();
    temp_out_sample.Result    = loop_out_sample_h;
    temp_out_sample.P_values  = loop_out_sample_p;
    temp_out_sample.Test_type = loop_sample_test_type;
    temp_out_sample.Analysed  = loop_sample_analysed;
    temp_out_sample.Alpha     = loop_alpha_sample;
    data.(in_motif_space).Sample_level_stats = temp_out_sample;
    
    % Global stats
    temp_out_global = struct();
    temp_out_global.Result    = loop_out_global_h;
    temp_out_global.P_values  = loop_out_global_p;
    temp_out_global.Test_type = loop_global_test_type;
    temp_out_global.Analysed  = loop_global_analysed;
    temp_out_global.Alpha     = loop_alpha_global;
    data.(in_motif_space).Global_stats = temp_out_global;
    
    % Create summary table for the motif space
    % Get counts
    temp_n_samples_significant     = sum((loop_out_sample_h == 1),2);
    temp_n_samples_not_significant = sum((loop_out_sample_h == 0),2);
    temp_n_samples_with_motif      = temp_n_samples_significant + temp_n_samples_not_significant;
    
    temp_table = struct();
    temp_table.Motif_ID        = loop_global_ERD_table.Motif_ID;
    temp_table.Motif           = loop_global_ERD_table.Motif;
    temp_table.Motif_analysed  = loop_global_ERD_table.Motif_analysed;
    temp_table.Samples_with_motif     = temp_n_samples_with_motif;
    temp_table.Total_motif_instances  = loop_global_ERD_table.Instances_total;
    temp_table.Global_MER      = loop_global_ERD_table.Sum_ERD;
    
    % Population level significance
    temp_table.Global_MER_significant = loop_out_global_h;
    temp_table.Global_p_value     = loop_out_global_p;
    temp_table.Global_test_type   = loop_global_test_type;
    
    % Sample-level ERD summary statistics
    temp_table.Sample_MER_mean   = mean(loop_sample_ERD,2,'omitnan');
    temp_table.Sample_MER_stdev  = std(loop_sample_ERD,1,2,'omitnan');
    temp_table.Sample_MER_median = median(loop_sample_ERD,2,'omitnan');
    temp_table.Sample_MER_Q1     = quantile(loop_sample_ERD,0.25,2);
    temp_table.Sample_MER_Q3     = quantile(loop_sample_ERD,0.75,2);
    
    % Sample-level significance overview
    temp_table.Samples_significant    = temp_n_samples_significant;
    temp_table.Samples_not_significant= temp_n_samples_not_significant;
    temp_table.Signif_sample_fraction = temp_n_samples_significant ./ temp_n_samples_with_motif;
    temp_table.Samples_FT_tests       = sum(loop_sample_test_type=="FT",2);
    temp_table.Samples_X2_tests       = sum(loop_sample_test_type=="X2",2);
    temp_table.Samples_ERD_negative   = sum(loop_sample_ERD<0,2);
    
    % Significance trhesholds used
    temp_table.Alpha                  = rempmat(parameters.alpha,[loop_n_motifs 1]);
    temp_table.Global_alpha_adjusted  = rempmat(loop_alpha_global,[loop_n_motifs 1]);
    temp_table.Sample_alpha_adjusted  = rempmat(loop_alpha_sample,[loop_n_motifs 1]);
    
    data.(in_motif_space).MER_stats_summary = struct2table(temp_table);
    
    %% Save loop outputs to MAT file
    if opts.save_MAT
        disp(append("Saving motif space ",in_motif_space," output variables to file"))

        % Output directory
        out_path = paths.output;
        if ~exist(out_path,'dir'); mkdir(out_path); end
    
        % Save sample level test results
        out_field = "Sample_level_stats";
        out_filepath = fullfile(out_path,append(in_motif_space,"_MER_Sample_level_stats.mat"));
        out = data.(in_motif_space).(out_field);
        save(out_filepath, '-struct', 'out');
    
        % Save global test results
        out_field = "Global_stats";
        out_filepath = fullfile(out_path,append(in_motif_space,"_MER_Global_stats.mat"));
        out = data.(in_motif_space).(out_field);
        save(out_filepath, '-struct', 'out');
    
        % Save stats summary table
        out_field = "MER_stats_summary";
        out_filepath = fullfile(out_path,append(in_motif_space,"_MER_stats_summary.mat"));
        out = data.(in_motif_space).(out_field);
        save(out_filepath, '-struct', 'out');
    
        % Clear temp
        clear -regexp ^out_
    end
    
    
    %% Write loop outputs to TXT file
    if opts.save_TXT
        disp(append("Writing motif space ",in_motif_space," outputs to txt file"))

        % Output directory
        out_path = paths.output;
        if ~exist(out_path,'dir'); mkdir(out_path); end
    
        % Get motif list and sample IDs for tables
        out_motif_list = loop_sample_ERD_table(:,1:2);
        out_sample_IDs = loop_sample_ERD_table.Properties.VariableNames(2:end);
    
        % 1) Save sample level stats tables as txt files
        out_fields = ["Result", "P_value", "Test_type"];
        for out_field = out_fields
    
            % File path
            out_filepath = fullfile(out_path,append(in_motif_space,"_Sample_MER_significance_",out_field,".txt"));
    
            % Assemble output table
            out_table = [out_motif_list, ...
                array2table(data.(in_motif_space).Sample_level_stats.(out_field), 'VariableNames', out_sample_IDs)];
    
            % Write to file
            writetable(out_table, out_filepath, 'Delimiter','\t');
        end
    
        % 2) Save population level stats table
        % File path
        out_filepath = fullfile(out_path,append(in_motif_space,"_Population_MER_significance_test_results.txt"));
    
        % Assemble output table
        out_table = struct();
        out_table.Result    = data.(in_motif_space).Sample_level_stats.Result;
        out_table.P_value   = data.(in_motif_space).Sample_level_stats.P_value;
        out_table.Alpha     = data.(in_motif_space).MER_stats_summary.Global_alpha_adjusted;
        out_table.Test_type = data.(in_motif_space).Sample_level_stats.Test_type;
        out_table = [out_motif_list, struct2table(out_table)];
    
        % Write to file
        writetable(out_table, out_filepath, 'Delimiter','\t');
    
        % 3) Write MER stats summary table to txt file
        out_table = data.(in_motif_space).MER_stats_summary.Global_alpha_adjusted;
        out_filepath = fullfile(out_path,append(in_motif_space,"_MER_stats_summary.txt"));
        writetable(out_table, out_filepath, 'Delimiter','\t');
    
        % Clear temp
        clear -regexp ^out_
    end
    
    %% Cleanup
    clear -regexp ^temp_
    clear -regexp ^loop_
    clear -regexp ^time_

    disp(append("Mptif space ",in_motif_space," completed"))
end
disp("All mptif spaces completed")
clear -regexp ^in






%% LOCAL FUNCTIONS
function [h,p,x2,alpha] = chi2cont(conttab, alpha)
% CHI2CONT takes as input a 2x2 matrix that represents a 2x2 contingency table and
% calculates the probability of obtaining the observed and each of the more extreme tables
% based on the pearson chi square test of independence. 
% Chi square test may become unreliable with low contingency table values (total < 20 or 
% any one value < 5), and Fishers exact test should be used instead.
% Requires the Statistics Toolbox. 
%
% Usage
% p = chi2cont(conttab)
% [h,p,x2,alpha] = chi2cont(conttab)
% [__] = chi2cont(conttab, alpha)
%
% Inputs
% CONTTAB - a 2x2 contigency table arraycreated from the frequency data
% ALPHA   - significance threshold (optional), default value 0.05
%
% Outputs 
% P       - p-value of the test (probability of observing as or more extreme table)
% X2      - The chi square statistic value
% H       - Test decision at specific alpha
% ALPHA   - Significance threshold used
% 
%
% Based on function CHISQUARECONT by Panagiotis Moulos:
% Panagiotis Moulos (2008). Chi-square test - contingency tables 
% (https://www.mathworks.com/matlabcentral/fileexchange/18705-chi-square-test-contingency-tables), 
% MATLAB Central File Exchange. Retrieved October 15, 2019.
%
    %% Argument validation
    % Check input arguments
    switch nargin
        case 0; error('CHI2CONT:BadInput','Input argument missing.')
        case 1; alpha = 0.05;
        case 2; if alpha<0||alpha>1; error('CHI2CONT:BadInput','Alpha must between 0 and 1.'); end
        otherwise; error('CHI2CONT:BadInput','Too many input arguments.')
    end

    % Check input contingency table shape
    if ~all(size(conttab)==[2 2])
        error('CHI2CONT:BadInput','Contingency table dimensions must be [2 2].')
    end
    
    % Issue warning if input vaues too low (total frequencies are < 20 or any cell value is < 5).
    if sum(sum(conttab))<20 || ~all(all(conttab>5))
        warning('CHI2CONT:Unreliable','CHI2CONT unreliable with low values.');
    end

    %% Chi square test calculation
    OBS = reshape(conttab',1,4);

    % Calculate expected values
    EXP = zeros(1,4);

    EXP(1) = ((OBS(1)+OBS(2)).*(OBS(1)+OBS(3)))./sum(OBS);   
    EXP(2) = ((OBS(1)+OBS(2)).*(OBS(2)+OBS(4)))./sum(OBS); 
    EXP(3) = ((OBS(3)+OBS(4)).*(OBS(1)+OBS(3)))./sum(OBS);   
    EXP(4) = ((OBS(3)+OBS(4)).*(OBS(2)+OBS(4)))./sum(OBS);

    % Calculate chi square statistic
    x2=sum((OBS-EXP).^2./EXP);

    % Find p-value (with 1 degree of freedom)
    p=1-chi2cdf(x2,1);
    h = p<=alpha;
end
