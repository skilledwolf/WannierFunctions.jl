# The Wannier90 checkpoint file (`seedname.chk`) — Fortran sequential-unformatted records.
#
# This is the full-precision interchange point with `wannier90.x`/`postw90.x`: the final gauge
# U (and U_opt + window data when disentangled), the final rotated overlaps M, centres, and
# spreads. Reading lets us consume an existing Wannier90 run; writing lets `wannier90.x`
# `restart = plot` (and postw90) consume ours. Record sequence per
# w90_readwrite_read_chkpt{_header,_matrices} (readwrite.F90:2296-2532); every `write(unit)` is
# one record framed by int32 byte-length markers (gfortran convention).

"Everything a Wannier90 checkpoint stores."
struct Checkpoint
    header::String                              # 33-char date stamp
    exclude_bands::Vector{Int}
    real_lattice::Matrix{Float64}               # 3×3, columns = a_i (Å)
    recip_lattice::Matrix{Float64}              # 3×3, columns = b_i (Å⁻¹)
    kpt_latt::Matrix{Float64}                   # 3 × num_kpts, fractional
    mp_grid::NTuple{3,Int}
    nntot::Int
    checkpoint::String                          # position tag: "postdis" | "postwann"
    have_disentangled::Bool
    omega_invariant::Float64                    # Ω_I (disentangled only, else NaN)
    lwindow::Union{Nothing,Matrix{Bool}}        # num_bands × num_kpts
    ndimwin::Union{Nothing,Vector{Int}}
    u_matrix_opt::Union{Nothing,Array{ComplexF64,3}}  # num_bands × num_wann × num_kpts
    u_matrix::Array{ComplexF64,3}               # num_wann × num_wann × num_kpts
    m_matrix::Array{ComplexF64,4}               # num_wann × num_wann × nntot × num_kpts
    centres::Matrix{Float64}                    # 3 × num_wann (Å)
    spreads::Vector{Float64}                    # num_wann (Ų)
end

num_wann(c::Checkpoint) = size(c.u_matrix, 1)
num_bands(c::Checkpoint) = c.u_matrix_opt === nothing ? size(c.u_matrix, 1) : size(c.u_matrix_opt, 1)

function Base.show(io::IO, ::MIME"text/plain", c::Checkpoint)
    print(io, "Checkpoint \"", c.checkpoint, "\": ", num_bands(c), " bands → ", num_wann(c),
          " WF, ", size(c.kpt_latt, 2), " k-points",
          c.have_disentangled ? @sprintf(", disentangled (Ω_I = %.9f)", c.omega_invariant) : "")
end

# --- Fortran sequential-unformatted record framing --------------------------------------------

"Read one record's payload, verifying the int32 length markers."
function _frec(io::IO)
    n = read(io, Int32)
    buf = read(io, n)
    n2 = read(io, Int32)
    n == n2 || error("corrupt Fortran record: length markers $n ≠ $n2")
    return buf
end

_frec(io::IO, ::Type{T}) where {T} = reinterpret(T, _frec(io))

"Write one record with framing."
function _wrec(io::IO, payload::AbstractVector{UInt8})
    write(io, Int32(length(payload)))
    write(io, payload)
    write(io, Int32(length(payload)))
    return nothing
end
_wrec(io::IO, x::AbstractArray) = _wrec(io, Vector{UInt8}(reinterpret(UInt8, vec(x))))
_wrec(io::IO, s::AbstractString, len::Int) = _wrec(io, Vector{UInt8}(rpad(s, len)[1:len]))

# --- reader ------------------------------------------------------------------------------------

"""
    read_chk(path) -> Checkpoint

Read a Wannier90 binary checkpoint (`seedname.chk`).
"""
function read_chk(path::AbstractString)
    open(path, "r") do io
        header = String(_frec(io))
        nb = Int(only(_frec(io, Int32)))
        nexcl = Int(only(_frec(io, Int32)))
        excl = Int.(_frec(io, Int32))                 # empty record when nexcl == 0
        length(excl) == nexcl || error("exclude_bands record length mismatch")
        A = reshape(Vector(_frec(io, Float64)), 3, 3)
        B = reshape(Vector(_frec(io, Float64)), 3, 3)
        nk = Int(only(_frec(io, Int32)))
        mp = Int.(_frec(io, Int32))
        kl = reshape(Vector(_frec(io, Float64)), 3, nk)
        nntot = Int(only(_frec(io, Int32)))
        nw = Int(only(_frec(io, Int32)))
        tag = strip(String(_frec(io)))
        havedis = only(_frec(io, Int32)) != 0         # Fortran default logical = 4 bytes

        ωI = NaN
        lwindow = nothing; ndimwin = nothing; uopt = nothing
        if havedis
            ωI = only(_frec(io, Float64))
            lwindow = reshape(_frec(io, Int32) .!= 0, nb, nk)
            ndimwin = Int.(_frec(io, Int32))
            uopt = reshape(Vector(_frec(io, ComplexF64)), nb, nw, nk)
        end
        u = reshape(Vector(_frec(io, ComplexF64)), nw, nw, nk)
        m = reshape(Vector(_frec(io, ComplexF64)), nw, nw, nntot, nk)
        centres = reshape(Vector(_frec(io, Float64)), 3, nw)
        spreads = Vector(_frec(io, Float64))

        # NB: the reference reads real/recip lattice as ((i,j),i=1,3),j=1,3) with
        # real_lattice(row=i, col=j) where ROW i is lattice vector a_i; our Lattice stores
        # vectors as COLUMNS, so transpose on the way in.
        return Checkpoint(header, excl, Matrix(transpose(A)), Matrix(transpose(B)), kl,
                          (mp[1], mp[2], mp[3]), nntot, String(tag), havedis, ωI,
                          lwindow, ndimwin, uopt, u, m, centres, spreads)
    end
end

# --- writer ------------------------------------------------------------------------------------

"""
    write_chk(path, chk::Checkpoint)

Write a Wannier90 binary checkpoint readable by `wannier90.x` (`restart = plot/wannierise`)
and `postw90.x`.
"""
function write_chk(path::AbstractString, c::Checkpoint)
    open(path, "w") do io
        _wrec(io, c.header, 33)
        _wrec(io, Int32[num_bands(c)])
        _wrec(io, Int32[length(c.exclude_bands)])
        _wrec(io, Int32.(c.exclude_bands))
        _wrec(io, Matrix(transpose(c.real_lattice)))   # rows = a_i on disk
        _wrec(io, Matrix(transpose(c.recip_lattice)))
        _wrec(io, Int32[size(c.kpt_latt, 2)])
        _wrec(io, Int32[c.mp_grid...])
        _wrec(io, c.kpt_latt)
        _wrec(io, Int32[c.nntot])
        _wrec(io, Int32[num_wann(c)])
        _wrec(io, c.checkpoint, 20)
        _wrec(io, Int32[c.have_disentangled ? 1 : 0])
        if c.have_disentangled
            _wrec(io, [c.omega_invariant])
            _wrec(io, Int32.(c.lwindow))
            _wrec(io, Int32.(c.ndimwin))
            _wrec(io, c.u_matrix_opt)
        end
        _wrec(io, c.u_matrix)
        _wrec(io, c.m_matrix)
        _wrec(io, c.centres)
        _wrec(io, c.spreads)
    end
    return path
end

"""
    Checkpoint(model, win, result) -> Checkpoint

Assemble a checkpoint from a completed `run_wannier` result, in the state `wannier90.x` labels
`"postwann"`. For disentangled runs the window bookkeeping and the rectangular `U_opt` are
reconstructed from the result's `DisentangleResult`.
"""
function Checkpoint(model::Model, win::WinInput, result::WannierResult)
    nw = model.num_wann
    nb = model.num_bands
    nk = nkpt(model.kgrid)
    kl = Matrix{Float64}(undef, 3, nk)
    for k in 1:nk
        kl[:, k] = model.kgrid.frac[k]
    end
    lwindow = nothing; ndimwin = nothing; uopt = nothing
    ωI = NaN
    if result.disentangled
        dis = result.dis
        eig = model.eig
        wd = dis_windows(eig, nw;
                         win_min=win.dis_win_min, win_max=win.dis_win_max,
                         froz_min=win.dis_froz_min,
                         froz_max=(win.dis_froz_max == -Inf ? nothing : win.dis_froz_max))
        lwindow = falses(nb, nk) |> Matrix{Bool}
        for k in 1:nk, i in 1:wd.ndimwin[k]
            lwindow[wd.winbands[k][i], k] = true
        end
        ndimwin = copy(wd.ndimwin)
        uopt = zeros(ComplexF64, nb, nw, nk)
        for k in 1:nk
            uopt[1:wd.ndimwin[k], :, k] = dis.Uopt[k]   # window-local rows, zero-padded above
        end
        ωI = dis.omega_I
    end
    # Γ-only: the model was expanded to the closed b-set on load; the checkpoint (like the
    # .mmn/.nnkp) stores only the file half, which is the first nntot/2 slots by construction.
    nntot = model.bvectors.nntot
    Mout = result.Mrot
    if win.gamma_only
        nntot = nntot ÷ 2
        Mout = Mout[:, :, 1:nntot, :]
    end
    return Checkpoint("written by WannierFunctions.jl", Int[],
                      Matrix(model.lattice.A), Matrix(model.lattice.B), kl,
                      model.kgrid.mp_grid, nntot, "postwann",
                      result.disentangled, ωI, lwindow, ndimwin, uopt,
                      result.U, Mout, result.spread.centres, result.spread.spreads)
end

# --- formatted variant (.chk.fmt) — the cross-platform transport format of w90chk2chk.x -------
# Line-oriented, list-directed: same record sequence as the binary, floats as "re im" pairs,
# logicals as 0/1 integers (conv_read_chkpt_fmt / conv_write_chkpt_fmt in w90chk2chk.F90).

"""
    read_chk_fmt(path) -> Checkpoint

Read a formatted checkpoint (`seedname.chk.fmt`, as produced by `w90chk2chk.x -u2f`).
"""
function read_chk_fmt(path::AbstractString)
    it = eachline(path)
    st = iterate(it)
    nextline() = (l = st[1]; st = iterate(it, st[2]); l)
    ints(l) = parse.(Int, split(l))
    flts(l) = parse.(Float64, replace.(split(l), r"[dD]" => "e"))

    header = nextline()
    nb = ints(nextline())[1]
    nexcl = ints(nextline())[1]
    excl = [ints(nextline())[1] for _ in 1:nexcl]
    A = reshape(flts(nextline()), 3, 3)
    B = reshape(flts(nextline()), 3, 3)
    nk = ints(nextline())[1]
    mp = ints(nextline())
    kl = Matrix{Float64}(undef, 3, nk)
    for k in 1:nk
        kl[:, k] = flts(nextline())
    end
    nntot = ints(nextline())[1]
    nw = ints(nextline())[1]
    tag = strip(nextline())
    havedis = ints(nextline())[1] != 0
    ωI = NaN
    lwindow = nothing; ndimwin = nothing; uopt = nothing
    cplx() = (v = flts(nextline()); complex(v[1], v[2]))
    if havedis
        ωI = flts(nextline())[1]
        lwindow = Matrix{Bool}(undef, nb, nk)
        for k in 1:nk, i in 1:nb
            lwindow[i, k] = ints(nextline())[1] != 0
        end
        ndimwin = [ints(nextline())[1] for _ in 1:nk]
        uopt = Array{ComplexF64,3}(undef, nb, nw, nk)
        for k in 1:nk, j in 1:nw, i in 1:nb
            uopt[i, j, k] = cplx()
        end
    end
    u = Array{ComplexF64,3}(undef, nw, nw, nk)
    for k in 1:nk, j in 1:nw, i in 1:nw
        u[i, j, k] = cplx()
    end
    m = Array{ComplexF64,4}(undef, nw, nw, nntot, nk)
    for k in 1:nk, n in 1:nntot, j in 1:nw, i in 1:nw
        m[i, j, n, k] = cplx()
    end
    centres = Matrix{Float64}(undef, 3, nw)
    for j in 1:nw
        centres[:, j] = flts(nextline())
    end
    spreads = [flts(nextline())[1] for _ in 1:nw]
    return Checkpoint(String(strip(header)), excl, Matrix(transpose(A)), Matrix(transpose(B)),
                      kl, (mp[1], mp[2], mp[3]), nntot, String(tag), havedis, ωI,
                      lwindow, ndimwin, uopt, u, m, centres, spreads)
end

"""
    write_chk_fmt(path, chk::Checkpoint)

Write a formatted checkpoint readable by `w90chk2chk.x -f2u` (and any Wannier90 build).
"""
function write_chk_fmt(path::AbstractString, c::Checkpoint)
    g(x) = @sprintf("%25.17E", x)
    open(path, "w") do io
        println(io, rpad(c.header, 33)[1:33])
        println(io, num_bands(c))
        println(io, length(c.exclude_bands))
        for b in c.exclude_bands
            println(io, b)
        end
        At = transpose(c.real_lattice)                 # rows = a_i on disk
        println(io, join((g(At[i, j]) for j in 1:3 for i in 1:3), ""))
        Bt = transpose(c.recip_lattice)
        println(io, join((g(Bt[i, j]) for j in 1:3 for i in 1:3), ""))
        println(io, size(c.kpt_latt, 2))
        println(io, join(c.mp_grid, " "))
        for k in 1:size(c.kpt_latt, 2)
            println(io, g(c.kpt_latt[1, k]), g(c.kpt_latt[2, k]), g(c.kpt_latt[3, k]))
        end
        println(io, c.nntot)
        println(io, num_wann(c))
        println(io, c.checkpoint)
        println(io, c.have_disentangled ? 1 : 0)
        if c.have_disentangled
            println(io, g(c.omega_invariant))
            for k in 1:size(c.kpt_latt, 2), i in 1:num_bands(c)
                println(io, c.lwindow[i, k] ? 1 : 0)
            end
            for k in 1:size(c.kpt_latt, 2)
                println(io, c.ndimwin[k])
            end
            for k in 1:size(c.kpt_latt, 2), j in 1:num_wann(c), i in 1:num_bands(c)
                println(io, g(real(c.u_matrix_opt[i, j, k])), g(imag(c.u_matrix_opt[i, j, k])))
            end
        end
        nw = num_wann(c)
        for k in 1:size(c.kpt_latt, 2), j in 1:nw, i in 1:nw
            println(io, g(real(c.u_matrix[i, j, k])), g(imag(c.u_matrix[i, j, k])))
        end
        for k in 1:size(c.kpt_latt, 2), n in 1:c.nntot, j in 1:nw, i in 1:nw
            println(io, g(real(c.m_matrix[i, j, n, k])), g(imag(c.m_matrix[i, j, n, k])))
        end
        for j in 1:nw
            println(io, g(c.centres[1, j]), g(c.centres[2, j]), g(c.centres[3, j]))
        end
        for j in 1:nw
            println(io, g(c.spreads[j]))
        end
    end
    return path
end

"""
    gauge_v_windows(chk, nb) -> (vs, winidx)

Per-q windowed gauge matrices `v(q) = U_opt·U` (ndimwin(q) × nw rows of the window) and the
window band indices, shared by every postw90-style operator assembly (A/B/C, spin, sHu/sIu).
For non-disentangled runs `v = U` on all `nb` bands.
"""
function gauge_v_windows(chk::Checkpoint, nb::Int)
    nk = size(chk.kpt_latt, 2)
    vs = Vector{Matrix{ComplexF64}}(undef, nk)
    winidx = Vector{Vector{Int}}(undef, nk)
    for q in 1:nk
        if chk.have_disentangled
            nd = chk.ndimwin[q]
            vs[q] = chk.u_matrix_opt[1:nd, :, q] * chk.u_matrix[:, :, q]
            winidx[q] = findall(@view chk.lwindow[:, q])
        else
            vs[q] = chk.u_matrix[:, :, q]
            winidx[q] = collect(1:nb)
        end
    end
    return vs, winidx
end
