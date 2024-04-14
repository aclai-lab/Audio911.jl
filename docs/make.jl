using Audio911
using Documenter

DocMeta.setdocmeta!(Audio911, :DocTestSetup, :(using Audio911); recursive = true)

makedocs(;
    modules = [Audio911],
    authors = "Mauro Milella, Giovanni Pagliarini, Alberto Paparella, Eduard I. Stan",
    repo=Documenter.Remotes.GitHub("aclai-lab", "Audio911.jl"),
    sitename = "Audio911.jl",
    format = Documenter.HTML(;
        size_threshold = 4000000,
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://aclai-lab.github.io/Audio911.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting-started.md",
    ],
    # NOTE: warning
    warnonly = :true,
)

@info "`makedocs` has finished running. "

deploydocs(;
    repo = "github.com/aclai-lab/Audio911.jl",
    target = "build",
    branch = "gh-pages",
    versions = ["main" => "main", "stable" => "v^", "v#.#", "dev" => "dev"],
)
