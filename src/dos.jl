# Density of states on the Wannier-interpolated bands (postw90's dos module), with the same
# adaptive Gaussian smearing as the Kubo path: per-band width η_n = min(fac·|∇ε_n|·Δk, η_max).
#
#   DOS(ε) = (1/N_k) Σ_{n,k} g_s w_n(k) δ_{η_nk}(ε − ε_nk)      [states/eV per cell]
#
# w_n(k) = Σ_{j∈project} |U_jn(k)|² for a WF-projected DOS (1 otherwise); with a `SpinModel`
# the DOS is split into up/down channels with weights (1 ± ⟨σ·n̂⟩_nk)/2 — the reference does
# NOT apply the g_s electron factor to the spin channels (and forces g_s = 1 for spinors).

using LinearAlgebra
using StaticArrays

"""
    density_of_states(bm; energies, kmesh, adaptive=true, adpt_fac=√2, adpt_max=1.0,
                      smr_width=0.0, spin=nothing, project=nothing, polar=0.0, azimuth=0.0,
                      elec_per_state=2) -> (energies, dos[, dos_up, dos_dn])

DOS on the energy grid (eV). With `spin` (a `SpinModel`), also returns the spin-decomposed
DOS (pass `elec_per_state = 1` for spinor calculations, as postw90 enforces). With `project`
(WF indices), the DOS is projected onto those Wannier functions.
"""
function density_of_states(bm::BerryModel; energies::AbstractVector{<:Real},
                           kmesh::NTuple{3,Int}=(25, 25, 25),
                           adaptive::Bool=true, adpt_fac::Float64=sqrt(2.0),
                           adpt_max::Float64=1.0, smr_width::Float64=0.0,
                           spin::Union{Nothing,SpinModel}=nothing,
                           project::Union{Nothing,Vector{Int}}=nothing,
                           polar::Float64=0.0, azimuth::Float64=0.0,
                           elec_per_state::Int=2)
    es = collect(Float64, energies)
    ne = length(es)
    nktot = prod(kmesh)
    Δk = kmesh_spacing(bm.lattice, kmesh)
    kl = [SVector(i / kmesh[1], j / kmesh[2], k / kmesh[3])
          for i in 0:kmesh[1]-1 for j in 0:kmesh[2]-1 for k in 0:kmesh[3]-1]
    nch = spin === nothing ? 1 : 3            # total / up / down
    acc = [zeros(ne, nch) for _ in 1:nktot]
    Threads.@threads for idx in 1:nktot
        _dos_kpoint!(acc[idx], bm, kl[idx], es, Δk, adaptive, adpt_fac, adpt_max, smr_width,
                     spin, project, polar, azimuth, elec_per_state)
    end
    tot = sum(acc) ./ nktot
    return spin === nothing ? (es, tot[:, 1]) : (es, tot[:, 1], tot[:, 2], tot[:, 3])
end

function _dos_kpoint!(out::Matrix{Float64}, bm::BerryModel, kf::SVector{3,Float64},
                      es::Vector{Float64}, Δk::Float64, adaptive::Bool, adpt_fac::Float64,
                      adpt_max::Float64, smr_width::Float64, spin, project,
                      polar::Float64, azimuth::Float64, elec_per_state::Int)
    E, dE, U = eig_deleig_vec(bm, kf; deriv=adaptive)
    nw = length(E)
    binwidth = length(es) > 1 ? es[2] - es[1] : 1.0
    spn = spin === nothing ? nothing : _spin_diag(spin, kf, U; polar=polar, azimuth=azimuth)
    for n in 1:nw
        w = project === nothing ? 1.0 : sum(abs2(U[j, n]) for j in project)
        αsq = spn === nothing ? 0.0 : (1.0 + spn[n]) / 2.0
        η = adaptive ?
            min(adpt_fac * norm(SVector(dE[1, n], dE[2, n], dE[3, n])) * Δk, adpt_max) : smr_width
        # When the width is small compared to the grid bin, the reference switches to a
        # histogram: the whole state goes into the nearest bin with weight 1/binwidth
        # (dos.F90:696-716, min_smearing_binwidth_ratio = 2).
        if η / binwidth < 2.0
            ie = round(Int, (E[n] - es[1]) / binwidth) + 1
            if 1 <= ie <= length(es)                  # states outside the window are dropped
                rdum = w / binwidth
                out[ie, 1] += rdum * elec_per_state
                if spn !== nothing
                    out[ie, 2] += rdum * αsq
                    out[ie, 3] += rdum * (1.0 - αsq)
                end
            end
        else
            for (ie, ε) in enumerate(es)
                abs(ε - E[n]) > 10.0 * η && continue        # smearing_cutoff
                rdum = w * _w0gauss((ε - E[n]) / η) / η
                out[ie, 1] += rdum * elec_per_state
                if spn !== nothing
                    out[ie, 2] += rdum * αsq
                    out[ie, 3] += rdum * (1.0 - αsq)
                end
            end
        end
    end
    return out
end

"Write `seedname-dos.dat` in the postw90 format (energy, dos[, up, down])."
function write_dos(path::AbstractString, es::Vector{Float64}, cols::Vector{Float64}...)
    open(path, "w") do io
        for i in 1:length(es)
            print(io, fortran_e(es[i], 16, 8))
            for c in cols
                print(io, fortran_e(c[i], 16, 8))
            end
            println(io)
        end
    end
    return path
end
