function output = getOptPolicySaveStrAppend(inp)
% Flexible function that gives the appended string indicating the set of
% parameters used
% if input is the param structure, output the appended string
% if inpus is the appended string, output the params


%Get string identifier for the set of parameters
if isstruct(inp)
    params = inp;
    if ~isfield(params,'aGamma_a')
        output = sprintf('_zsRes(%.04f)_c(%.04f)_sig2(%.04f)_sig2_z(%.04f)_dt(%.04f)_T(%.04f)_aG(%.04f)_c_s(%.04f)_Cs(%.04f)_t_s(%.04f)',...
            params.zsRes,params.c,params.sig2,params.sig2_z,params.dt,params.T,params.aGamma,params.c_s,params.Cs,params.t_s);
    else
        % Make distinction when using flexible aGamma, since the attended
        % gamma will be 1-aGamma
        output = sprintf('_zsRes(%.04f)_c(%.04f)_sig2(%.04f)_sig2_z(%.04f)_dt(%.04f)_T(%.04f)_aG(%.04f)_aGamma_a(%.04f)_aGamma_n(%.04f)_c_s(%.04f)_Cs(%.04f)_t_s(%.04f)',...
            params.zsRes,params.c,params.sig2,params.sig2_z,params.dt,params.T,params.aGamma,params.aGamma_a,params.aGamma_n,params.c_s,params.Cs,params.t_s);
    end
elseif isstr(inp)
    i_parentheses = strfind(inp,')');
    allparamsearchstr = {'c(','sig2(','sig2_z(','aG(','c_s(','zsRes(','dt(','T(','t_s('};
    paramstr = {'c','sig2','sig2_z','aGamma','c_s','zsRes','dt','T','t_s'};
    params = struct;
    for p = 1:length(paramstr)
        i_paramstr = strfind(inp,allparamsearchstr{p});
        parentheses_after = i_parentheses(i_parentheses>i_paramstr);
        params.(paramstr{p}) = str2num(inp(i_paramstr+length(allparamsearchstr{p}) : parentheses_after(1)-1));
    end
    output = params;
end



end

