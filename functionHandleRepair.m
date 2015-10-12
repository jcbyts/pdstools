function obj=functionHandleRepair(obj, verbose)
% fix broken function handles in a loaded struct or handle class
% obj=functionHandleRepair(obj)
if nargin==1
    verbose=true;
end

if isstruct(obj)
    fields=fieldnames(obj);
else
    fields=properties(obj);
end

for kField=1:numel(fields)
    if isnumeric(obj.(fields{kField})) || ischar(obj.(fields{kField}))
        continue
    end
    n=numel(obj.(fields{kField}));
    if n>1
        for k=1:n
            if isa(obj.(fields{kField})(k), 'struct') && isfield(obj.(fields{kField})(k), 'workspace')
                obj.(fields{kField})(k)=reconstructFcn(obj.(fields{kField})(k));
            elseif isa(obj.(fields{kField})(k), 'handle') || (isa(obj.(fields{kField})(k), 'struct') && ~isfield(obj.(fields{kField})(k), 'workspace'))
                obj.(fields{kField})(k)=functionHandleRepair(obj.(fields{kField})(k), verbose);
            end
        end
    else
        if isa(obj.(fields{kField}), 'struct') && isfield(obj.(fields{kField}), 'workspace')
            obj.(fields{kField})=reconstructFcn(obj.(fields{kField}));
            if verbose
                fprintf('Function [%s] reconstructed\n', (fields{kField}))
            end
            %             if strcmp(S.type, 'anonymous')
            %                 keyboard
            %
            %                 fstr=func2str(obj.(fields{kField}));
            %                 if strncmp(fstr,'@',1)
            %                     obj.(fields{kField})=eval(fstr);
            %                     fprintf('function handle [%s] fixed\n', fields{kField})
            %                 end
            %             end
        elseif isa(obj.(fields{kField}), 'handle') || (isa(obj.(fields{kField}), 'struct') && ~isfield(obj.(fields{kField}), 'workspace'))
            obj.(fields{kField})=functionHandleRepair(obj.(fields{kField}), verbose);
        elseif isa(obj.(fields{kField}), 'function_handle')
            f=functions(obj.(fields{kField}));
            if strcmp(f.type, 'anonymous')
                obj.(fields{kField})=f;
                if verbose
                    fprintf('Function handle [%s] meta data included\n', fields{kField})
                end
            end
        end
    end
end
if verbose
    fprintf('Done\n')
end

function out = reconstructFcn(s)

for iWks = 1:numel(s.workspace)
    wkspace = s.workspace{iWks};
    varnames = fieldnames(wkspace);
    if ~isempty(varnames)
        for i = 1:numel(varnames)
            try
                tmp = wkspace.(varnames{i});
                eval([varnames{i} ' = tmp;']);
            end
        end
    end
end

fcn_str = s.function;
fcn = eval(fcn_str);
out = fcn;


