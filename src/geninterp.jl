# Generalised interpolation (postw90's geninterp module): band energies — and optionally their
# Cartesian k-derivatives (band velocities, Hellmann–Feynman: dE_n/dk_α = Re[U†∂_αH U]_nn) —
# at an arbitrary user-supplied k-point list, written in the seedname_geninterp.dat format.

using Printf
using StaticArrays

"Light per-k interpolation: eigenvalues and (optionally) their Cartesian derivatives (eV·Å)."
function eig_deleig(bm::BerryModel, kfrac::AbstractVector; deriv::Bool=true)
    nw = num_wann(bm)
    H = zeros(ComplexF64, nw, nw)
    dH = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    kf = SVector{3,Float64}(kfrac...)
    for ir in 1:length(bm.irvec)
        fac = cis(TWOPI * dot(kf, SVector{3,Float64}(bm.irvec[ir]...))) / bm.ndegen[ir]
        @views H .+= fac .* bm.Hr[:, :, ir]
        if deriv
            R = bm.Rcart[ir]
            @views for c in 1:3
                dH[c] .+= (fac * im * R[c]) .* bm.Hr[:, :, ir]
            end
        end
    end
    F = eigen(Hermitian((H + H') / 2))
    E, U = F.values, F.vectors
    dE = zeros(3, nw)
    if deriv
        for c in 1:3
            dHh = U' * dH[c] * U
            for n in 1:nw
                dE[c, n] = real(dHh[n, n])
            end
        end
    end
    return E, dE
end

"""
    read_geninterp_kpt(path) -> (comment, idx, kpts)

Read a `seedname_geninterp.kpt` file: comment line, coordinate flag (`crystal`/`frac`), count,
then `index kx ky kz` per line.
"""
function read_geninterp_kpt(path::AbstractString)
    lines = readlines(path)
    comment = strip(lines[1])
    flag = lowercase(strip(lines[2]))
    startswith(flag, "crystal") || startswith(flag, "frac") ||
        error("geninterp kpt: unsupported coordinate flag `$flag` (crystal/frac only)")
    n = parse(Int, strip(lines[3]))
    idx = Vector{Int}(undef, n)
    kpts = Vector{SVector{3,Float64}}(undef, n)
    for i in 1:n
        t = split(lines[3+i])
        idx[i] = parse(Int, t[1])
        kpts[i] = SVector{3,Float64}(parse_f64.(t[2:4])...)
    end
    return String(comment), idx, kpts
end

"""
    geninterp(bm, seedname; alsofirstder=true) -> path

Run generalised interpolation: read `seedname_geninterp.kpt`, interpolate E_n(k) (and dE/dk
when `alsofirstder`), write `seedname_geninterp.dat` in the postw90 format.
"""
function geninterp(bm::BerryModel, seedname::AbstractString; alsofirstder::Bool=true)
    comment, idx, kpts = read_geninterp_kpt(seedname * "_geninterp.kpt")
    out = seedname * "_geninterp.dat"
    open(out, "w") do io
        println(io, "# Written by WannierFunctions.jl")
        println(io, "# Input file comment: ", comment)
        if alsofirstder
            println(io, "#  Kpt_idx  K_x (1/ang)       K_y (1/ang)        K_z (1/ang)       " *
                        "Energy (eV)      EnergyDer_x       EnergyDer_y       EnergyDer_z")
        else
            println(io, "#  Kpt_idx  K_x (1/ang)       K_y (1/ang)        K_z (1/ang)       Energy (eV)")
        end
        for (i, kf) in enumerate(kpts)
            E, dE = eig_deleig(bm, kf; deriv=alsofirstder)
            kcart = bm.lattice.B * kf
            for n in 1:length(E)
                @printf(io, "%10d  ", idx[i])
                @printf(io, "%.10g      %.10g      %.10g      ", kcart[1], kcart[2], kcart[3])
                if alsofirstder
                    @printf(io, "%.10g      %.10g      %.10g      %.10g    \n",
                            E[n], dE[1, n], dE[2, n], dE[3, n])
                else
                    @printf(io, "%.10g    \n", E[n])
                end
            end
        end
    end
    return out
end
