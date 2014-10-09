function defaults = parseOptions(opts, defaults, nowarning)
% parses options
% opts = parseOptions(opts, defaults, nowarning)
% INPUTS
% 	opts 	 [struct]
% 	defaults [struct]
% 	nowarning  [bool]
% OUTPUTS
% 	opts 	[struct]

if isempty(opts)
	opts = struct();
end

assert(isstruct(opts), 'opts must be a struct')
assert(isstruct(defaults), 'defaults must be a struct')

if nargin < 3
	nowarning = 0;
end

requiredOptions = fieldnames(defaults);
existingOptions = fieldnames(opts);

for ii = 1:numel(existingOptions)
	if any(cellfun(@(x) strcmp(x, existingOptions{ii}), requiredOptions, 'UniformOutput', true))
		defaults.(existingOptions{ii}) = opts.(existingOptions{ii});
	elseif ~nowarning
		warning(sprintf('%s is not an option I know\n', existingOptions{ii}))
	end
end
