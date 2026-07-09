# Ballistic (Landauer–Büttiker) quantum conductance of a periodic chain — wannier90's
# `transport_mode = bulk` (transport.F90): the 1D-periodic Wannier Hamiltonian is cut into
# principal layers H00/H01, lead surface Green functions are built by the López-Sancho
# decimation (tran_transfer), and T(E) = Tr[Γ_L G^r Γ_R G^a] (Fisher–Lee) plus the
# Green-function DOS are scanned over an energy window. The `_htB.dat` read/write formats are
# drop-in compatible, so a wannier90.x `tran_read_ht` run on the same file is byte-comparable.

using LinearAlgebra

const TRAN_ETA = 0.0005im       # small imaginary part of E for the retarded GF (transport.F90:90)
const TRAN_NTERX = 50           # max decimation iterations

"""
    tran_transfer(H00, H01, ecmp) -> (T, T̄)

López-Sancho/Rubio decimation for the lead transfer matrices at complex energy `ecmp`
(port of `tran_transfer`; converges quadratically, errors after $(TRAN_NTERX) iterations).
"""
function tran_transfer(H00::Matrix{Float64}, H01::Matrix{Float64}, ecmp::ComplexF64)
    n = size(H00, 1)
    Id = Matrix{ComplexF64}(I, n, n)
    t11 = (ecmp * I - H00) \ Id                            # (e − H00)⁻¹
    tau = t11 * H01'
    taut = t11 * H01
    tot = copy(tau)
    tsum = copy(taut)
    tott = copy(taut)
    tsumt = copy(tau)
    for _ in 1:TRAN_NTERX
        s1 = I - tau * taut - taut * tau
        s2 = s1 \ Id
        tau2 = s2 * (tau * tau)
        taut2 = s2 * (taut * taut)
        tot += tsum * tau2
        tsum = tsum * taut2
        tott += tsumt * taut2
        tsumt = tsumt * tau2
        tau, taut = tau2, taut2
        conver = sum(abs, tau2)
        conver2 = sum(abs, taut2)
        (conver < 1e-7 && conver2 < 1e-7) && return tot, tott
    end
    error("tran_transfer: transfer matrix not converged")
end

"""
    tran_green(tot, tott, H00, H01, e; igreen) -> g

Green function at real energy `e` from the transfer matrices (port of `tran_green`):
`igreen = 0` bulk g_nn, `1` surface g00 (right lead), `-1` dual surface ḡ00 (left lead).
"""
function tran_green(tot::Matrix{ComplexF64}, tott::Matrix{ComplexF64},
                    H00::Matrix{Float64}, H01::Matrix{Float64}, e::Float64; igreen::Int=0)
    n = size(H00, 1)
    ehinv = if igreen == 1
        e * I - H00 - H01 * tot
    elseif igreen == -1
        e * I - H00 - H01' * tott
    elseif igreen == 0
        e * I - H00 - H01 * tot - H01' * tott
    else
        error("tran_green: igreen must be -1, 0 or 1")
    end
    return Matrix{ComplexF64}(ehinv) \ Matrix{ComplexF64}(I, n, n)
end

"""
    transport_bulk(H00, H01; win_min=-3.0, win_max=3.0, energy_step=0.01)
        -> (energies, qc, dos)

Landauer conductance T(E) (in units of 2e²/h) and Green-function DOS of the periodic chain
with principal-layer blocks `H00`/`H01` (Fermi level already subtracted from `H00`'s
diagonal), scanned over `win_min:energy_step:win_max` — the reference's `tran_bulk`.
"""
function transport_bulk(H00::Matrix{Float64}, H01::Matrix{Float64};
                        win_min::Float64=-3.0, win_max::Float64=3.0,
                        energy_step::Float64=0.01)
    n_e = floor(Int, (win_max - win_min) / energy_step) + 1
    energies = [win_min + (n - 1) * energy_step for n in 1:n_e]
    qc = Vector{Float64}(undef, n_e)
    dos = Vector{Float64}(undef, n_e)
    for (i, e) in enumerate(energies)
        tot, tott = tran_transfer(H00, H01, e + TRAN_ETA)
        gB = tran_green(tot, tott, H00, H01, e; igreen=0)
        sLr = H01' * tott                       # Σ_L^r
        sRr = H01 * tot                         # Σ_R^r
        ΓL = im * (sLr - sLr')
        ΓR = im * (sRr - sRr')
        c1 = ΓL * gB * ΓR * gB'
        qc[i] = real(tr(c1))
        dos[i] = -imag(tr(gB)) / π
    end
    return energies, qc, dos
end

"""
    transport_from_tb(lattice, irvec, Hr, centres, mp_grid;
                      one_dim_axis=:z, fermi_energy, dist_cutoff=1500.0,
                      dist_cutoff_mode="one_dim", hr_cutoff=0.0) -> (H00, H01, num_pl)

Assemble the principal-layer blocks from a real-space Wannier Hamiltonian (the reference's
`tran_reduce_hr` + `tran_cut_hr_one_dim` + `tran_get_ht`): keep H(R) with R purely along the
1D axis, zero elements beyond `dist_cutoff` (WF-centre distances, `one_dim` or full 3D mode),
choose the principal-layer size from `hr_cutoff` decay, and tile H00/H01. The Fermi energy is
subtracted from the H00 diagonal.
"""
function transport_from_tb(lattice::Lattice, irvec::Vector{NTuple{3,Int}},
                           Hr::Array{ComplexF64,3}, centres::Matrix{Float64},
                           mp_grid::NTuple{3,Int};
                           one_dim_axis::Symbol=:z, fermi_energy::Float64,
                           dist_cutoff::Float64=1000.0,
                           dist_cutoff_mode::AbstractString="three_dim",
                           hr_cutoff::Float64=0.0)
    nw = size(Hr, 1)
    dir = one_dim_axis === :x ? 1 : one_dim_axis === :y ? 2 : 3
    # lattice vector parallel to the transport axis
    onedv = 0
    for i in 1:3
        abs(abs(lattice.A[dir, i]) - norm(lattice.A[:, i])) < 1e-8 && (onedv = i)
    end
    onedv != 0 || error("transport: no lattice vector is parallel to one_dim_axis $one_dim_axis")

    irvec_max = maximum(r -> max(abs(r[1]), abs(r[2]), abs(r[3])), irvec) + 1
    hr1 = zeros(nw, nw, 2irvec_max + 1)                    # index n1 + irvec_max + 1
    two = filter(!=(onedv), 1:3)
    for n1 in -irvec_max:irvec_max
        for (ir, r) in enumerate(irvec)
            if mod(n1 - r[onedv], mp_grid[onedv]) == 0 && r[two[1]] == 0 && r[two[2]] == 0
                hr1[:, :, n1+irvec_max+1] = real.(Hr[:, :, ir])
                break
            end
        end
    end

    # distance cutoff between WF centres shifted by n1 lattice vectors
    maxdist = mp_grid[onedv] * abs(lattice.A[dir, onedv]) / 2
    dcut = min(dist_cutoff, maxdist)
    for i in 1:nw, j in 1:nw, n1 in -irvec_max:irvec_max
        d = if occursin("one_dim", dist_cutoff_mode)
            abs(centres[dir, i] - centres[dir, j] + n1 * lattice.A[dir, onedv])
        else
            norm(centres[:, i] - centres[:, j] + n1 .* lattice.A[:, onedv])
        end
        d > dcut && (hr1[j, i, n1+irvec_max+1] = 0.0)
    end

    # principal layer size from the H(R) decay
    num_pl = 0
    for n1 in -irvec_max:irvec_max
        hmax = maximum(abs, @view hr1[:, :, n1+irvec_max+1])
        if hmax > hr_cutoff
            num_pl = max(num_pl, abs(n1))
        else
            hr1[:, :, n1+irvec_max+1] .= 0.0
        end
    end
    for n1 in -num_pl:num_pl, i in 1:nw, j in 1:nw
        abs(hr1[j, i, n1+irvec_max+1]) < hr_cutoff && (hr1[j, i, n1+irvec_max+1] = 0.0)
    end

    nbb = num_pl * nw
    H00 = zeros(nbb, nbb)
    H01 = zeros(nbb, nbb)
    for j in 0:num_pl-1, i in 0:num_pl-1
        H00[j*nw+1:(j+1)*nw, i*nw+1:(i+1)*nw] = hr1[:, :, i-j+irvec_max+1]
    end
    for j in 1:num_pl, i in 0:j-1
        H01[(j-1)*nw+1:j*nw, i*nw+1:(i+1)*nw] = hr1[:, :, i-j+1+num_pl+irvec_max+1]
    end
    for i in 1:nbb
        H00[i, i] -= fermi_energy
    end
    return H00, H01, num_pl
end

"""
    translate_centres_home(centres, lattice, atoms_frac_mean) -> centres

Translate Cartesian WF centres into the home cell centred on the mean atomic position
(the reference's `internal_translate_centres` with `automatic_translation`): fractional
coordinates are shifted into `[c̄ − 0.5, c̄ + 0.5)`.
"""
function translate_centres_home(centres::Matrix{Float64}, lattice::Lattice,
                                cfrac::AbstractVector{Float64})
    A = Matrix(lattice.A)
    Ainv = inv(A)
    out = similar(centres)
    for n in 1:size(centres, 2)
        rf = Ainv * centres[:, n]
        rf .+= -floor.(rf .- (cfrac .- 0.5))
        out[:, n] = A * rf
    end
    return out
end

"""
    run_transport(model, win, res; seedname) -> (energies, qc, dos)

The `.win`-driven `transport_mode = bulk` flow: build H(R) from the wannierised result,
translate the WF centres home, assemble the principal layers and scan the energy window,
writing `_qc.dat` / `_dos.dat` (and `_htB.dat` with `tran_write_ht`).
"""
function run_transport(model::Model, win::WinInput, res::WannierResult;
                       seedname::AbstractString=model.seedname)
    occursin("bulk", lowercase(get(win.raw, "transport_mode", "bulk"))) ||
        error("transport: only transport_mode = bulk is implemented (tran_lcr is not)")
    haskey(win.raw, "fermi_energy") ||
        error("transport requires fermi_energy to be set")
    op = hamiltonian_operator(model, res)
    # mean atomic position (fractional) for the automatic home-cell translation
    atoms = get(win.blocks, "atoms_frac", nothing)
    cfrac = zeros(3)
    if atoms !== nothing
        nat = 0
        for ln in atoms
            t = split(ln)
            length(t) >= 4 || continue
            nat += 1
            cfrac .+= parse_f64.(t[2:4])
        end
        nat > 0 && (cfrac ./= nat)
    end
    centres = translate_centres_home(res.spread.centres, model.lattice, cfrac)
    axis = Symbol(lowercase(strip(get(win.raw, "one_dim_axis", "z"))))
    H00, H01, _ = transport_from_tb(model.lattice, op.irvec, op.data[:, :, :, 1], centres,
                                    model.kgrid.mp_grid;
                                    one_dim_axis=axis,
                                    fermi_energy=_getfloat(win.raw, "fermi_energy", 0.0),
                                    dist_cutoff=_getfloat(win.raw, "dist_cutoff", 1000.0),
                                    dist_cutoff_mode=get(win.raw, "dist_cutoff_mode", "three_dim"),
                                    hr_cutoff=_getfloat(win.raw, "hr_cutoff", 0.0))
    _getbool(win.raw, "tran_write_ht", false) && write_ht(seedname * "_htB.dat", H00, H01)
    energies, qc, dos = transport_bulk(H00, H01;
                                       win_min=_getfloat(win.raw, "tran_win_min", -3.0),
                                       win_max=_getfloat(win.raw, "tran_win_max", 3.0),
                                       energy_step=_getfloat(win.raw, "tran_energy_step", 0.01))
    write_transport(seedname, energies, qc, dos)
    return energies, qc, dos
end

"""
    write_ht(path, H00, H01; header) — write a `seedname_htB.dat` (reference format).
"""
function write_ht(path::AbstractString, H00::Matrix{Float64}, H01::Matrix{Float64};
                  header::AbstractString="written by WannierFunctions.jl")
    n = size(H00, 1)
    open(path, "w") do io
        println(io, " ", header)
        println(io, lpad(n, 6))
        _write_f126(io, H00)
        println(io, lpad(n, 6))
        _write_f126(io, H01)
    end
    return path
end

"Column-major 6F12.6 block writer (the reference's `(6F12.6)` format)."
function _write_f126(io::IO, M::Matrix{Float64})
    vals = vec(M)                                # column-major, (j,i) order as the reference
    for start in 1:6:length(vals)
        stop = min(start + 5, length(vals))
        println(io, join((@sprintf("%12.6f", v) for v in vals[start:stop])))
    end
end

"""
    read_ht(path) -> (H00, H01)  — read a `seedname_htB.dat`.
"""
function read_ht(path::AbstractString)
    lines = readlines(path)                      # line 1 is a free-text header
    p = 2
    n = parse(Int, strip(lines[p]))
    vals = Float64[]
    p += 1
    while length(vals) < n * n
        append!(vals, parse.(Float64, split(lines[p])))
        p += 1
    end
    H00 = reshape(vals[1:n*n], n, n)
    n2 = parse(Int, strip(lines[p]))
    n2 == n || error("read_ht: inconsistent block sizes $n vs $n2")
    vals2 = Float64[]
    p += 1
    while length(vals2) < n * n && p <= length(lines)
        append!(vals2, parse.(Float64, split(lines[p])))
        p += 1
    end
    H01 = reshape(vals2[1:n*n], n, n)
    return Matrix(H00), Matrix(H01)
end

"""
    write_transport(seedname, energies, qc, dos) — write `_qc.dat` / `_dos.dat` (reference format).
"""
function write_transport(seedname::AbstractString, energies, qc, dos)
    open(seedname * "_qc.dat", "w") do io
        println(io, " ## written by WannierFunctions.jl")
        for (e, q) in zip(energies, qc)
            @printf(io, "%15.9f%18.9f\n", e, q)
        end
    end
    open(seedname * "_dos.dat", "w") do io
        println(io, " ## written by WannierFunctions.jl")
        for (e, d) in zip(energies, dos)
            @printf(io, "%15.9f%18.9f\n", e, d)
        end
    end
    return seedname * "_qc.dat"
end
