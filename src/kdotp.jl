# k·p expansion coefficients by Löwdin quasi-degenerate perturbation theory around a
# reference k-point (postw90's berry_task = kdotp; exact conventions in
# docs/reference-notes/kdotp.md). Only H(R) is needed.
#
#   H^{kp}(κ) ≈ T0 + Σ_a T1_a κ_a + Σ_ab T2_ab κ_a κ_b        [eV, eV·Å, eV·Å²]
#
#   T0(n,m)    = δ_nm ε_bn
#   T1_a(n,m)  = [U†∂_aH U](bn,bm)
#   T2_ab(n,m) = ½[U†∂²_abH U](bn,bm)
#                + ½ Σ_{r∉A} [U†∂_aH U](bn,r)[U†∂_bH U](r,bm)·(1/(ε_bn−ε_r) + 1/(ε_bm−ε_r))
#
# with A = `bands` (the quasi-degenerate set). The stored T2 already contains the ½ — contract
# with κ_a κ_b directly.

using Printf
using StaticArrays

"""
    kdotp(bm; kpoint=[0,0,0], bands) -> (; T0, T1, T2)

k·p coefficient matrices around fractional `kpoint` for the interpolated-band indices `bands`.
`T0` is nA×nA (eV), `T1[a]` nA×nA (eV·Å), `T2[a,b]` nA×nA (eV·Å²), all in the H(k) eigenbasis.
"""
function kdotp(bm::BerryModel; kpoint::AbstractVector=[0.0, 0.0, 0.0], bands::Vector{Int})
    nw = num_wann(bm)
    all(b -> 1 <= b <= nw, bands) || error("kdotp: bands out of range 1:$nw")
    kf = SVector{3,Float64}(kpoint...)
    H = zeros(ComplexF64, nw, nw)
    dH = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    d2H = [zeros(ComplexF64, nw, nw) for _ in 1:3, _ in 1:3]
    for ir in 1:length(bm.irvec)
        fac = cis(TWOPI * dot(kf, SVector{3,Float64}(bm.irvec[ir]...))) / bm.ndegen[ir]
        Rc = bm.Rcart[ir]
        @views H .+= fac .* bm.Hr[:, :, ir]
        @views for a in 1:3
            dH[a] .+= (fac * im * Rc[a]) .* bm.Hr[:, :, ir]
            for b in 1:3
                d2H[a, b] .+= (-fac * Rc[a] * Rc[b]) .* bm.Hr[:, :, ir]
            end
        end
    end
    F = eigen(Hermitian((H + H') / 2))
    E, U = F.values, F.vectors
    dHh = [U' * dH[a] * U for a in 1:3]
    d2Hh = [U' * d2H[a, b] * U for a in 1:3, b in 1:3]

    nA = length(bands)
    inA = falses(nw)
    inA[bands] .= true
    T0 = ComplexF64.(Diagonal(E[bands]))
    T1 = [ComplexF64[dHh[a][bn, bm] for bn in bands, bm in bands] for a in 1:3]
    T2 = Matrix{Matrix{ComplexF64}}(undef, 3, 3)
    for a in 1:3, b in 1:3
        T = zeros(ComplexF64, nA, nA)
        for (j, bmj) in enumerate(bands), (i, bni) in enumerate(bands)
            s = 0.5 * d2Hh[a, b][bni, bmj]
            for r in 1:nw
                inA[r] && continue
                s += 0.5 * dHh[a][bni, r] * dHh[b][r, bmj] *
                     (1.0 / (E[bni] - E[r]) + 1.0 / (E[bmj] - E[r]))
            end
            T[i, j] = s
        end
        T2[a, b] = T
    end
    return (; T0, T1, T2)
end

"""
    write_kdotp(seedname, res) -> seedname

Write `seedname-kdotp_{0,1,2}.dat` in the postw90 layout: one complex value per line in
`(2E18.8E3)`, column-major within each block; order-1 has three blocks (a = x,y,z), order-2
nine blocks (a,b) with b fastest.
"""
function write_kdotp(seedname::AbstractString, res)
    e18(z) = fortran_e(real(z), 18, 8; edigits=3) * fortran_e(imag(z), 18, 8; edigits=3)
    open(seedname * "-kdotp_0.dat", "w") do io
        for z in vec(res.T0)
            println(io, e18(z))
        end
    end
    open(seedname * "-kdotp_1.dat", "w") do io
        for a in 1:3, z in vec(res.T1[a])
            println(io, e18(z))
        end
    end
    open(seedname * "-kdotp_2.dat", "w") do io
        for a in 1:3, b in 1:3, z in vec(res.T2[a, b])
            println(io, e18(z))
        end
    end
    return seedname
end
