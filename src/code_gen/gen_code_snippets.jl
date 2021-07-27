# Code snippet handling
struct CodeSnippet
    code_lines::Vector{String}
    lang::Any
end
function init_code(lang)
    return CodeSnippet(Vector{String}(undef, 0), lang)
end
function push_code!(code, str; ind_lvl = 1, ind_str = "    ")
    # Indent only non-empty lines.
    indentation = repeat(ind_str, ind_lvl)
    indented_string = isempty(str) ? "" : indentation * str
    return push!(code.code_lines, indented_string)
end
function push_code_verbatim_string!(code, str)
    for line in split(str, "\n")
        push!(code.code_lines, line)
    end
end
function push_comment!(code, str; ind_lvl = 1, ind_str = "    ")
    # Convert empty comments to empty lines.
    return push_code!(
        code,
        isempty(str) ? "" : comment(code.lang, str),
        ind_lvl = ind_lvl,
        ind_str = ind_str,
    )
end
function to_string(code)
    return join(code.code_lines, "\n")
end
