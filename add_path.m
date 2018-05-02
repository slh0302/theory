%% 首先运行该脚本，再运行solution中的脚本
function paths = add_path()
% toolkit_path Appends all toolkit directories to the Matlab path
%
% If no output argument is given, the function appends all toolkit directories
% to the Matlab/Octave path, otherwise it returns a list as a cell array.
%
% Output:
% - paths (cell): A list of all paths.

script_directory = fileparts(mfilename('fullpath'));
include_dirs = cellfun(@(x) fullfile(script_directory, x), {'', 'basic', 'linner', 'solution','level', 'solution-2', 'parallel_solution_2'}, 'UniformOutput', false);

if exist(fullfile(script_directory, 'native'), 'dir')
   include_dirs{end+1} = fullfile(script_directory, 'native');
end

if nargout > 0
    paths = include_dirs;
else
    addpath(include_dirs{:});
end
