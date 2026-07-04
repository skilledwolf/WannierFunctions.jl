using Documenter
using WannierFunctions

makedocs(;
    sitename = "WannierFunctions.jl",
    modules = [WannierFunctions],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true",
                             edit_link = "main"),
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
        "File formats" => "file-formats.md",
        "Migrating from Wannier90" => "migrating-from-wannier90.md",
        "API" => "api.md",
    ],
    warnonly = [:missing_docs],   # not every internal helper carries a docstring yet
)

# deploydocs is a no-op until the repository has a public URL; fill in `repo` when publishing:
# deploydocs(; repo = "github.com/<org>/WannierFunctions.jl", devbranch = "main")
