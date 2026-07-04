using Documenter
using WannierFunctions

# The landing page is the README with repo-relative links rewritten for the docs tree.
readme = read(joinpath(@__DIR__, "..", "README.md"), String)
readme = replace(readme,
    "docs/theory.md" => "theory.md",
    "docs/file-formats.md" => "file-formats.md",
    "docs/migrating-from-wannier90.md" => "migrating-from-wannier90.md")
write(joinpath(@__DIR__, "src", "index.md"), readme)

makedocs(;
    sitename = "WannierFunctions.jl",
    modules = [WannierFunctions],
    # remotes/edit links need a public repository URL; set both when publishing.
    remotes = nothing,
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true",
                             edit_link = nothing),
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
        "File formats" => "file-formats.md",
        "Migrating from Wannier90" => "migrating-from-wannier90.md",
        "API" => "api.md",
    ],
    warnonly = [:missing_docs, :cross_references],   # README keeps a few repo-relative links
)

# deploydocs is a no-op until the repository has a public URL; fill in `repo` when publishing:
# deploydocs(; repo = "github.com/<org>/WannierFunctions.jl", devbranch = "main")
