# Readers for the Wannier90 file formats: .win (input), .amn (projections),
# .mmn (overlaps), .eig (band energies). Writers live in output.jl.

using StaticArrays

"Parse a Fortran-style float, tolerating `d`/`D` exponent markers (e.g. \"17.0d0\")."
parse_f64(s::AbstractString) = parse(Float64, replace(strip(s), r"[dD]" => "e"))

# ---------------------------------------------------------------------------
# Input validation. Three tiers:
#   SUPPORTED    — consumed by this package (silent)
#   IGNORED_OK   — recognised, genuinely irrelevant to what we compute (silent)
#   known        — in the reference parser's catalogue but not supported here (warn once)
#   anything else — a typo: error, with a did-you-mean suggestion
# ---------------------------------------------------------------------------

const SUPPORTED_KEYWORDS = Set{String}([
    "num_wann", "num_bands", "mp_grid", "num_iter", "kmesh_tol", "search_shells",
    "use_ws_distance", "dis_win_min", "dis_win_max", "dis_froz_min", "dis_froz_max",
    "dis_num_iter", "dis_mix_ratio", "conv_tol", "conv_window", "trial_step", "num_cg_steps",
    "bands_plot", "bands_num_points", "write_hr", "hr_plot", "write_tb", "guiding_centres",
    "postproc_setup", "exclude_bands", "gamma_only", "spinors",
    "wannier_plot", "wannier_plot_list", "wannier_plot_supercell", "wannier_plot_format",
])
const IGNORED_KEYWORDS = Set{String}([
    "wvfn_formatted", "num_print_cycles", "iprint", "timing_level",
])
const SUPPORTED_BLOCKS = Set{String}([
    "unit_cell_cart", "kpoints", "kpoint_path", "atoms_frac", "atoms_cart", "projections",
])

"Levenshtein edit distance (for did-you-mean suggestions)."
function _levenshtein(a::AbstractString, b::AbstractString)
    m, n = length(a), length(b)
    d = collect(0:n)
    for (i, ca) in enumerate(a)
        prev = d[1]
        d[1] = i
        for (j, cb) in enumerate(b)
            cur = d[j+1]
            d[j+1] = min(d[j+1] + 1, d[j] + 1, prev + (ca == cb ? 0 : 1))
            prev = cur
        end
    end
    return d[n+1]
end

function _closest(name::String, pool)
    best, bestd = "", typemax(Int)
    for cand in pool
        dist = _levenshtein(name, cand)
        dist < bestd && ((best, bestd) = (cand, dist))
    end
    return bestd <= max(2, length(name) ÷ 3) ? best : ""
end

"Validate parsed keys/blocks against the reference catalogue; error on typos."
function _validate_win(raw::Dict{String,String}, blocks::Dict{String,Vector{String}},
                       path::AbstractString; strict::Bool=true)
    for key in keys(raw)
        key in SUPPORTED_KEYWORDS && continue
        key in IGNORED_KEYWORDS && continue
        if key in W90_KNOWN_KEYWORDS
            @warn "$path: keyword `$key` is recognised (a wannier90 keyword) but not supported " *
                  "by WannierFunctions.jl — ignoring it" _id = Symbol(:unsup_, key) maxlog = 1
        else
            sugg = _closest(key, union(W90_KNOWN_KEYWORDS, SUPPORTED_KEYWORDS))
            msg = "$path: unknown keyword `$key`" * (isempty(sugg) ? "" : " — did you mean `$sugg`?")
            strict ? error(msg) : @warn msg
        end
    end
    for name in keys(blocks)
        name in SUPPORTED_BLOCKS && continue
        if name in W90_KNOWN_BLOCKS || name in W90_KNOWN_KEYWORDS
            @warn "$path: block `$name` is recognised but not supported by WannierFunctions.jl — " *
                  "ignoring it" _id = Symbol(:unsupblk_, name) maxlog = 1
        else
            sugg = _closest(name, union(W90_KNOWN_BLOCKS, SUPPORTED_BLOCKS))
            msg = "$path: unknown block `$name`" * (isempty(sugg) ? "" : " — did you mean `$sugg`?")
            strict ? error(msg) : @warn msg
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# .win input file
# ---------------------------------------------------------------------------

"""
    WinInput

Parsed `.win` input. Scalar parameters not promoted to typed fields are kept in `raw`
(lower-cased keys), and block bodies (raw lines, comments stripped) in `blocks`.
"""
struct WinInput
    num_wann::Int
    num_bands::Int
    mp_grid::NTuple{3,Int}
    unit_cell::SMatrix{3,3,Float64,9}   # columns = lattice vectors a₁,a₂,a₃ (Å)
    kpoints::Vector{SVector{3,Float64}}
    num_iter::Int
    kmesh_tol::Float64
    search_shells::Int
    use_ws_distance::Bool
    dis_win_min::Float64
    dis_win_max::Float64
    dis_froz_min::Float64
    dis_froz_max::Float64
    dis_num_iter::Int
    dis_mix_ratio::Float64
    gamma_only::Bool
    raw::Dict{String,String}
    blocks::Dict{String,Vector{String}}
end

"Strip a trailing `!`/`#` comment and surrounding whitespace."
function strip_comment(line::AbstractString)
    for (i, c) in pairs(line)
        (c == '!' || c == '#') && return strip(line[1:prevind(line, i)])
    end
    return strip(line)
end

"Split the raw text of a `.win` file into scalar key/value pairs and named blocks."
function _parse_win_structure(text::AbstractString)
    raw = Dict{String,String}()
    blocks = Dict{String,Vector{String}}()
    lines = split(text, '\n')
    i = 1
    while i <= length(lines)
        line = strip_comment(lines[i])
        if isempty(line)
            i += 1
            continue
        end
        low = lowercase(line)
        if startswith(low, "begin ")
            name = strip(low[7:end])
            body = String[]
            i += 1
            while i <= length(lines)
                inner = strip_comment(lines[i])
                if startswith(lowercase(inner), "end ")
                    break
                end
                isempty(inner) || push!(body, String(inner))
                i += 1
            end
            blocks[name] = body
        else
            # scalar assignment: key = value  OR  key : value
            m = match(r"^([A-Za-z0-9_]+)\s*[:=]\s*(.*)$", line)
            if m !== nothing
                raw[lowercase(m.captures[1])] = String(strip(m.captures[2]))
            end
        end
        i += 1
    end
    return raw, blocks
end

_getint(raw, key, default) = haskey(raw, key) ? parse(Int, split(raw[key])[1]) : default
_getfloat(raw, key, default) = haskey(raw, key) ? parse_f64(split(raw[key])[1]) : default
function _getbool(raw, key, default)
    haskey(raw, key) || return default
    v = lowercase(raw[key])
    return occursin("t", v) && !occursin("f", v)
end

"Parse a `unit_cell_cart` block into a 3×3 matrix whose columns are the lattice vectors (Å)."
function _parse_cell(body::Vector{String})
    scale = 1.0
    rows = SVector{3,Float64}[]
    for ln in body
        toks = split(ln)
        if length(toks) == 1
            u = lowercase(toks[1])
            (u == "bohr" || u == "b") && (scale = BOHR)
            continue                     # "ang"/"angstrom" ⇒ scale stays 1
        end
        length(toks) >= 3 || continue
        push!(rows, SVector{3,Float64}(parse_f64.(toks[1:3])))
    end
    length(rows) == 3 || error("unit_cell_cart must have 3 lattice-vector rows")
    # w90 lists each lattice vector as a ROW; store as columns.
    A = hcat((scale .* rows)...)         # columns = a₁,a₂,a₃
    return SMatrix{3,3,Float64}(A)
end

"Parse a `kpoints` block (fractional coordinates, one k-point per row)."
function _parse_kpoints(body::Vector{String})
    kpts = SVector{3,Float64}[]
    for ln in body
        toks = split(ln)
        length(toks) >= 3 || continue
        push!(kpts, SVector{3,Float64}(parse_f64.(toks[1:3])))
    end
    return kpts
end

"""
    read_win(path; strict=true) -> WinInput

Parse a Wannier90 `.win` input file. Only the parameters needed for wannierisation and
interpolation are promoted to typed fields; the remainder stay in `raw`/`blocks`.

Unknown keywords (not in the reference wannier90 parser's catalogue) are an **error** with a
did-you-mean suggestion — a silently ignored typo like `num_itre` is worse than a hard stop.
Recognised-but-unsupported keywords warn once and are ignored. Pass `strict=false` to downgrade
unknown-keyword errors to warnings.
"""
function read_win(path::AbstractString; strict::Bool=true)
    text = read(path, String)
    raw, blocks = _parse_win_structure(text)
    _validate_win(raw, blocks, path; strict=strict)

    haskey(raw, "num_wann") || error("$path: num_wann is required")
    num_wann = parse(Int, split(raw["num_wann"])[1])
    num_bands = _getint(raw, "num_bands", num_wann)

    haskey(raw, "mp_grid") || error("$path: mp_grid is required")
    mp = parse.(Int, split(raw["mp_grid"]))
    length(mp) == 3 || error("$path: mp_grid needs 3 integers")

    haskey(blocks, "unit_cell_cart") || error("$path: unit_cell_cart block is required")
    cell = _parse_cell(blocks["unit_cell_cart"])

    haskey(blocks, "kpoints") || error("$path: kpoints block is required")
    kpoints = _parse_kpoints(blocks["kpoints"])

    return WinInput(
        num_wann, num_bands, (mp[1], mp[2], mp[3]), cell, kpoints,
        _getint(raw, "num_iter", 100),
        _getfloat(raw, "kmesh_tol", KMESH_TOL_DEFAULT),
        _getint(raw, "search_shells", 36),
        _getbool(raw, "use_ws_distance", true),
        _getfloat(raw, "dis_win_min", -Inf),
        _getfloat(raw, "dis_win_max", Inf),
        _getfloat(raw, "dis_froz_min", -Inf),
        _getfloat(raw, "dis_froz_max", -Inf),
        _getint(raw, "dis_num_iter", 200),
        _getfloat(raw, "dis_mix_ratio", 0.5),
        _getbool(raw, "gamma_only", false),
        raw, blocks,
    )
end

# ---------------------------------------------------------------------------
# .amn — projection overlaps A[m,n,k] = ⟨ψ_{m,k}|g_n⟩
# ---------------------------------------------------------------------------

"""
    read_amn(path) -> (A, num_bands, num_kpts, num_wann)

`A[m,n,k]` is (num_bands × num_wann × num_kpts) complex. File layout: comment line;
`num_bands num_kpts num_wann`; then num_bands·num_wann·num_kpts records `m n k Re Im`.
"""
function read_amn(path::AbstractString)
    open(path, "r") do io
        readline(io)                                   # comment
        nb, nk, nw = parse.(Int, split(readline(io)))
        A = Array{ComplexF64,3}(undef, nb, nw, nk)
        for _ in 1:(nb * nw * nk)
            t = split(readline(io))
            m, n, k = parse(Int, t[1]), parse(Int, t[2]), parse(Int, t[3])
            A[m, n, k] = complex(parse_f64(t[4]), parse_f64(t[5]))
        end
        return A, nb, nk, nw
    end
end

# ---------------------------------------------------------------------------
# .mmn — overlaps M[m,n,b,k] = ⟨u_{m,k}|u_{n,k+b}⟩ + neighbour connectivity
# ---------------------------------------------------------------------------

"""
    read_mmn(path) -> (M, kpb, gpb, num_bands, num_kpts, nntot)

`M[m,n,b,k]` is (num_bands × num_bands × nntot × num_kpts). `kpb[b,k]` is the neighbour
k-index and `gpb[:,b,k]` its reciprocal-lattice fold. File layout: comment line;
`num_bands num_kpts nntot`; then num_kpts·nntot blocks, each a line `k kb g1 g2 g3`
followed by num_bands² records `Re Im` with the row index m varying fastest.
"""
function read_mmn(path::AbstractString)
    open(path, "r") do io
        readline(io)                                   # comment
        nb, nk, nntot = parse.(Int, split(readline(io)))
        M = Array{ComplexF64,4}(undef, nb, nb, nntot, nk)
        kpb = Matrix{Int}(undef, nntot, nk)
        gpb = Array{Int,3}(undef, 3, nntot, nk)
        slot = zeros(Int, nk)                          # per-k running neighbour counter
        for _ in 1:(nk * nntot)
            h = split(readline(io))
            k = parse(Int, h[1])
            b = (slot[k] += 1)
            kpb[b, k] = parse(Int, h[2])
            gpb[1, b, k] = parse(Int, h[3])
            gpb[2, b, k] = parse(Int, h[4])
            gpb[3, b, k] = parse(Int, h[5])
            for n in 1:nb, m in 1:nb                   # m fastest
                t = split(readline(io))
                M[m, n, b, k] = complex(parse_f64(t[1]), parse_f64(t[2]))
            end
        end
        return M, kpb, gpb, nb, nk, nntot
    end
end

# ---------------------------------------------------------------------------
# .eig — band energies eig[m,k] (eV)
# ---------------------------------------------------------------------------

"""
    read_eig(path) -> eig

`eig[m,k]` (num_bands × num_kpts), eV. Each line: `band_index kpt_index energy`.
"""
function read_eig(path::AbstractString)
    entries = Tuple{Int,Int,Float64}[]
    nb = nk = 0
    for ln in eachline(path)
        t = split(ln)
        isempty(t) && continue
        m, k = parse(Int, t[1]), parse(Int, t[2])
        nb = max(nb, m); nk = max(nk, k)
        push!(entries, (m, k, parse_f64(t[3])))
    end
    eig = Matrix{Float64}(undef, nb, nk)
    for (m, k, e) in entries
        eig[m, k] = e
    end
    return eig
end
