%% SCRIPT_1_Motif_space_generator
% Script for generating different length DNA sequence motif spaces
% Motif space = all possible sequence permutations of a given lenghts.
% 
% OUTPUT: 
% Motif_spaces   - Struct containging 8 motif spaces from 2bp (N2) to 8bp (N8) lengths
% Saved as individuals variables in file "Motif_spaces.mat"
%
% Script used PERMN function by Jos van der Geest
% Citeation: Jos van der Geest (10584) (2019). PERMN.
% (https://www.mathworks.com/matlabcentral/fileexchange/7147-permn), 
% MATLAB Central File Exchange. Retrieved March 6, 2019

%% Set parameters
opts = struct();
opts.space_min = 2;                      % smallest motif space to generate 
opts.space_max = 8;                      % largest motif space to generate 
opts.save = false;                       % Save outputs to MAT file? 
opts.path = '/path/to/output/directory'; % Path to directory for outputs
opts.filename = "Motif_spaces";      % File name for out pouts

%% Generate motif spaces
% Set output path
out_filepath = fullfile(opts.path, append(opts.filename,".mat"));

% Create Motif spaces struct
Motif_spaces = struct(); 

% Generate motif spaces
for motif_length = 2:4%8
    % Motif space ID
    motif_space_ID = sprintf('N%d',motif_length);
    disp(append("Generating motif space ",motif_space_ID))

    % Generate motif space
    Motif_spaces.(motif_space_ID) = flip(permn('ACGT', motif_length),2);

    % OPTIONAL: Save motif space as variable into MAT file
    if opts.save == true
        if ~exist(out_filepath,'file')
             save(out_filepath, '-struct', 'Motif_spaces', motif_space_ID)
        else;save(out_filepath, '-struct', 'Motif_spaces', motif_space_ID,'-append')
        end; disp(append("Motif space ",motif_space_ID," saved"))
    end
end

% Clear 
clear -regexp ^out
clear motif_length motif_space_ID
clear opts


%% Local functions
function [M, I] = permn(V, N, K)
% PERMN - permutations with repetition
%
%   M = PERMN(V,N) returns all permutations of N elements taken from the 
%   vector V, with repetitions. M has the size numel(V).^N-by-N.
%
% version 6.2 (jan 2019)
% (c) Jos van der Geest
% Matlab File Exchange Author ID: 10584
% email: samelinoa@gmail.com

narginchk(2, 3) ;

if fix(N) ~= N || N < 0 || numel(N) ~= 1
    error('permn:negativeN','Second argument should be a positive integer') ;
end
nV = numel(V) ;

if nargin==2 
    %% PERMN(V,N) - return all permutations
    if nV == 0 || N == 0
        M = zeros(nV, N) ;
        I = zeros(nV, N) ;
    elseif N == 1
        % return column vectors
        M = V(:) ;
        I = (1:nV).' ;
    else
        % this is faster than the math trick used with 3 inputs below
        [Y{N:-1:1}] = ndgrid(1:nV) ;
        I = reshape(cat(N+1, Y{:}), [], N) ;
        M = V(I) ;
    end
else
    %% PERMN(V,N,K) - return a subset of all permutations
    nK = numel(K) ;
    if nV == 0 || N == 0 || nK == 0
        M = zeros(numel(K), N) ;
        I = zeros(numel(K), N) ;
    elseif nK < 1 || any(K<1) || any(K ~= fix(K))
        error('permn:InvalidIndex','Third argument should contain positive integers.') ;
    else
        V = reshape(V, 1, []) ; % v1.1 make input a row vector
        nV = numel(V) ;
        Npos = nV^N ;
        if any(K > Npos)
            warning('permn:IndexOverflow', ...
                'Values of K exceeding the total number of combinations are saturated.')
            K = min(K, Npos) ;
        end
             
        % The engine is based on version 3.2 with the correction
        % suggested by Roger Stafford. This approach uses a single matrix
        % multiplication.
        B = nV.^(1-N:0) ;
        I = ((K(:)-.5) * B) ; % matrix multiplication
        I = rem(floor(I), nV) + 1 ;
        M = V(I) ;
    end
end
end
