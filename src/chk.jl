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
            lwindow[wd.nfirstwin[k]+i-1, k] = true
        end
        ndimwin = copy(wd.ndimwin)
        uopt = zeros(ComplexF64, nb, nw, nk)
        for k in 1:nk
            uopt[1:wd.ndimwin[k], :, k] = dis.Uopt[k]   # window-local rows, zero-padded above
        end
        ωI = dis.omega_I
    end
    return Checkpoint("written by WannierFunctions.jl", Int[],
                      Matrix(model.lattice.A), Matrix(model.lattice.B), kl,
                      model.kgrid.mp_grid, model.bvectors.nntot, "postwann",
                      result.disentangled, ωI, lwindow, ndimwin, uopt,
                      result.U, result.Mrot, result.spread.centres, result.spread.spreads)
end
