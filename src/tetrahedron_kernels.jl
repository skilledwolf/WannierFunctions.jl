# tetrahedron_kernels.jl
#
# Faithful, verbatim Julia port of the pure-numerical tetrahedron-integration
# functions from Wannier90's postw90/tetrahedron.F90.
#
# Reference: Minsu Ghim and Cheol-Hwan Park, PRB 106, 075126.
#            Kawamura correction, PRB 89, 094515.
#
# These are self-contained math kernels (Float64 arrays in, Float64 out) with NO
# external dependencies. Fortran is 1-based and column-major; Julia is 1-based and
# column-major too, so index bases are kept identical to the source (a(i) -> a[i],
# t(i,k) -> t[i,k]).
#
# No `module` wrapper — bare functions intended to be `include`d.

using LinearAlgebra

# ---------------------------------------------------------------------------
# tetrahedron_P_matrix_init  ->  tet_p_matrix() :: Matrix{Float64} (4x20)
# The Kawamura P-matrix: literal integer entries divided by 1260.0.
# ---------------------------------------------------------------------------
function tet_p_matrix()::Matrix{Float64}
    P = Matrix{Float64}(undef, 4, 20)
    # correction, Kawamura PRB 89 094515
    P[1, 1:4] = Float64[1440, 0, 30, 0]
    P[2, 1:4] = Float64[0, 1440, 0, 30]
    P[3, 1:4] = Float64[30, 0, 1440, 0]
    P[4, 1:4] = Float64[0, 30, 0, 1440]
    #
    P[1, 5:8] = Float64[-38, 7, 17, -28]
    P[2, 5:8] = Float64[-28, -38, 7, 17]
    P[3, 5:8] = Float64[17, -28, -38, 7]
    P[4, 5:8] = Float64[7, 17, -28, -38]
    #
    P[1, 9:12] = Float64[-56, 9, -46, 9]
    P[2, 9:12] = Float64[9, -56, 9, -46]
    P[3, 9:12] = Float64[-46, 9, -56, 9]
    P[4, 9:12] = Float64[9, -46, 9, -56]
    #
    P[1, 13:16] = Float64[-38, -28, 17, 7]
    P[2, 13:16] = Float64[7, -38, -28, 17]
    P[3, 13:16] = Float64[17, 7, -38, -28]
    P[4, 13:16] = Float64[-28, 17, 7, -38]
    #
    P[1, 17:20] = Float64[-18, -18, 12, -18]
    P[2, 17:20] = Float64[-18, -18, -18, 12]
    P[3, 17:20] = Float64[12, -18, -18, -18]
    P[4, 17:20] = Float64[-18, 12, -18, -18]
    #
    P[1:4, 1:20] = P[1:4, 1:20] ./ 1260.0
    return P
end

# ---------------------------------------------------------------------------
# tetrahedron_array_init  ->  tet_array() :: Matrix{Int} (6x20)
# The 6x20 integer stencil indices (with Kawamura correction).
# ---------------------------------------------------------------------------
function tet_array()::Matrix{Int}
    ta = Matrix{Int}(undef, 6, 20)
    # with correction, Kawamura PRB 89 094515
    ta[1, 1:4] = [22, 38, 39, 43]
    ta[2, 1:4] = [22, 42, 38, 43]
    ta[3, 1:4] = [22, 26, 42, 43]
    ta[4, 1:4] = [22, 27, 26, 43]
    ta[5, 1:4] = [22, 23, 27, 43]
    ta[6, 1:4] = [22, 39, 23, 43]
    #
    ta[1, 5:20] = [6, 37, 35, 64, 5, 33, 56, 48, 1, 54, 40, 47, 59, 23, 42, 18]
    ta[2, 5:20] = [2, 46, 33, 64, 6, 41, 54, 44, 1, 62, 34, 48, 63, 18, 47, 17]
    ta[3, 5:20] = [18, 10, 41, 64, 2, 9, 62, 60, 1, 30, 58, 44, 47, 38, 27, 21]
    ta[4, 5:20] = [17, 28, 9, 64, 18, 11, 30, 59, 1, 32, 25, 60, 48, 21, 44, 5]
    ta[5, 5:20] = [21, 19, 11, 64, 17, 3, 32, 63, 1, 24, 31, 59, 44, 26, 39, 6]
    ta[6, 5:20] = [5, 55, 3, 64, 21, 35, 24, 47, 1, 56, 7, 63, 60, 6, 59, 2]
    return ta
end

# ---------------------------------------------------------------------------
# tetrahedron_log1p  ->  tet_log1p(x) :: Float64
# Numerically-stable log(1+x).
# ---------------------------------------------------------------------------
function tet_log1p(x::Float64)::Float64
    if abs(x) > 0.5
        return log(abs(1.0 + x))
    else
        y = 1.0 + x
        z = y - 1.0
        if z == 0
            return x
        else
            return x * log(y) / z
        end
    end
end

# ---------------------------------------------------------------------------
# tetrahedron_sort  ->  tet_sort!(a, b1, b2, t)
# Size-4 bubble sort of `a` ascending, permuting payloads b1, b2 identically,
# and rebuilding the 3x3 edge matrix t so that after sorting
#   t[:, i-1] = t_temp[:, ref[i]] - t_temp[:, ref[1]]   for i = 2..4
# where t_temp has column 1 = 0 and columns 2..4 = original t columns 1..3.
#
# Returns the permuted (a, b1, b2, t) as new arrays.
# ---------------------------------------------------------------------------
function tet_sort!(a::Vector{Float64}, b1::Vector{Float64}, b2::Vector{Float64},
                   t::Matrix{Float64})
    # Work on copies so caller arrays are not aliased/mutated unexpectedly.
    a = copy(a)
    b1 = copy(b1)
    b2 = copy(b2)
    t = copy(t)

    reference = Vector{Int}(undef, 4)
    for i in 1:4
        reference[i] = i
    end
    b1_temp = copy(b1)
    b2_temp = copy(b2)

    t_temp = Matrix{Float64}(undef, 3, 4)
    for j in 1:3
        t_temp[j, 1] = 0.0
    end
    for i in 2:4
        for j in 1:3
            t_temp[j, i] = t[j, i-1]
        end
    end

    # bubble sort (ascending), tracking permutation in `reference`
    for i in 4:-1:1
        for j in 1:(i-1)
            if a[j] > a[j+1]
                temp = a[j]
                a[j] = a[j+1]
                a[j+1] = temp
                temp2 = reference[j]
                reference[j] = reference[j+1]
                reference[j+1] = temp2
            end
        end
    end

    # rearrange t(j,i)
    for i in 2:4
        for j in 1:3
            t[j, i-1] = t_temp[j, reference[i]] - t_temp[j, reference[1]]
        end
    end
    # rearrange b(j,i)
    for i in 1:4
        b1[i] = b1_temp[reference[i]]
        b2[i] = b2_temp[reference[i]]
    end

    return a, b1, b2, t
end

# ---------------------------------------------------------------------------
# utility_inv3 (from utility.F90) — inlined helper.
# Returns b = adjoint of 3x3 matrix a, and det.  inverse(a) = b/det.
# ---------------------------------------------------------------------------
function _tet_inv3(a::Matrix{Float64})
    b = Matrix{Float64}(undef, 3, 3)
    b[1, 1] = a[2, 2] * a[3, 3] - a[3, 2] * a[2, 3]
    b[1, 2] = a[2, 3] * a[3, 1] - a[3, 3] * a[2, 1]
    b[1, 3] = a[2, 1] * a[3, 2] - a[3, 1] * a[2, 2]
    b[2, 1] = a[3, 2] * a[1, 3] - a[1, 2] * a[3, 3]
    b[2, 2] = a[3, 3] * a[1, 1] - a[1, 3] * a[3, 1]
    b[2, 3] = a[3, 1] * a[1, 2] - a[1, 1] * a[3, 2]
    b[3, 1] = a[1, 2] * a[2, 3] - a[2, 2] * a[1, 3]
    b[3, 2] = a[1, 3] * a[2, 1] - a[2, 3] * a[1, 1]
    b[3, 3] = a[1, 1] * a[2, 2] - a[2, 1] * a[1, 2]
    det = a[1, 1] * b[1, 1] + a[1, 2] * b[1, 2] + a[1, 3] * b[1, 3]
    return b, det
end

# ---------------------------------------------------------------------------
# tetrahedron_jacobian  ->  tet_jacobian(t, x, type_) :: Float64
# Jacobian part of surface integrations (lines 543-594).
# ---------------------------------------------------------------------------
function tet_jacobian(t::Matrix{Float64}, x::Vector{Float64}, type_::Int)::Float64
    J = Matrix{Float64}(undef, 3, 2)
    if type_ == 1
        J[1, 1] = -x[1]; J[1, 2] = -x[1]
        J[2, 1] = x[2];  J[2, 2] = 0.0
        J[3, 1] = 0.0;   J[3, 2] = x[3]
    elseif type_ == 2
        J[1, 1] = x[1];  J[1, 2] = 0.0
        J[2, 1] = 0.0;   J[2, 2] = x[2]
        J[3, 1] = 1.0 - x[1] - x[3]; J[3, 2] = -x[3]
    elseif type_ == 3
        y = x[1] * (x[2] - 1.0) * x[3] / (-x[2] + x[1] * x[2] + x[2] * x[3] - x[1] * x[3])
        J[1, 1] = x[1] - y;   J[1, 2] = -y
        J[2, 1] = y - 1.0;    J[2, 2] = y - 1.0 + x[2]
        J[3, 1] = 1.0 - x[1]; J[3, 2] = 0.0
    else
        J[1, 1] = x[1];  J[1, 2] = 0.0
        J[2, 1] = 0.0;   J[2, 2] = x[2]
        J[3, 1] = 1.0 - x[1] - x[3]; J[3, 2] = 1.0 - x[2] - x[3]
    end
    Ans = 0.0
    for j_ in 1:3
        for k in 1:3
            for a in 1:3
                for b in 1:3
                    for c in 1:3
                        for d in 1:3
                            Ans = Ans + t[j_, a] * t[j_, b] * t[k, c] * t[k, d] *
                                        J[a, 1] * J[c, 2] * (J[b, 1] * J[d, 2] - J[b, 2] * J[d, 1])
                        end
                    end
                end
            end
        end
    end
    return sqrt(abs(Ans))
end

# ---------------------------------------------------------------------------
# tetrahedron_integral  ->  tet_integral(F, D, t, hw, type_, cutoff, avoid_deg)
# The core analytic single-tetrahedron integral (lines 270-486).
# F, D are length-4 Float64; t is 3x3.
# ---------------------------------------------------------------------------
function tet_integral(F_in::Vector{Float64}, D_in::Vector{Float64}, t_in::Matrix{Float64},
                      hw::Float64, type_::Int, tet_cutoff::Float64, avoid_deg::Float64)::Float64
    # copy inputs (Fortran intent(in): D = D_in; F = F_in; t = t_in)
    D = copy(D_in)
    F = copy(F_in)
    t = copy(t_in)

    dummy = zeros(Float64, 4)
    D, F, dummy, t = tet_sort!(D, F, dummy, t)
    Ans = 0.0

    dd = zeros(Float64, 3)
    ll = zeros(Float64, 3)
    ff = 1.0
    Det_t = 0.0

    # case 1 and 3: nondissipative part, case 2: dissipative part
    if type_ == 1 || type_ == 3
        if type_ == 3
            # treatment for accidental small band splitting (but degenerate actually)
            for j in 1:4
                if abs(D[j]) < avoid_deg
                    D[j] = avoid_deg * (abs(D[j]) / D[j])
                    fill!(F, 0.0)
                end
            end
        end

        # cutoff treatment, hw == 0.0 for case 3
        DAV = (D[2] + D[3]) / 2.0
        if abs((D[2] - D[3]) / (DAV + hw)) < tet_cutoff
            D_small_prev = D[2]
            D_large_prev = D[3]
            D[3] = DAV + 0.5 * abs(DAV + hw) * tet_cutoff
            D[2] = DAV - 0.5 * abs(DAV + hw) * tet_cutoff
            if D[1] > D[2]
                D[1] = D[1] + (D[2] - D_small_prev)
            end
            if D[3] > D[4]
                D[4] = D[4] + (D[3] - D_large_prev)
            end
        end
        DAV = (D[1] + D[2]) / 2.0
        if abs((D[1] - D[2]) / (DAV + hw)) < tet_cutoff
            if D[2] > 0
                D[1] = D[2] * (2.0 - tet_cutoff) / (2.0 + tet_cutoff)
            else
                D[1] = D[2] * (2.0 + tet_cutoff) / (2.0 - tet_cutoff)
            end
        end
        DAV = (D[3] + D[4]) / 2.0
        if abs((D[3] - D[4]) / (DAV + hw)) < tet_cutoff
            if D[3] > 0
                D[4] = D[3] * (2.0 + tet_cutoff) / (2.0 - tet_cutoff)
            else
                D[4] = D[3] * (2.0 - tet_cutoff) / (2.0 + tet_cutoff)
            end
        end

        # intermediate variables for case 1 and 3
        for i in 1:3
            dd[i] = (D[4] - D[i]) / (D[i] + hw)
            ll[i] = tet_log1p(dd[i])
        end
        ff = 1.0

        # determinant factor from parametrisation (tetrahedron volume)
        Det_t = abs(t[1, 1] * t[2, 2] * t[3, 3] + t[1, 2] * t[2, 3] * t[3, 1] + t[1, 3] * t[2, 1] * t[3, 2]
                    - t[1, 1] * t[2, 3] * t[3, 2] - t[1, 3] * t[2, 2] * t[3, 1] - t[1, 2] * t[2, 1] * t[3, 3])
    end

    bb = zeros(Float64, 4)
    cc = zeros(Float64, 4, 3)

    if type_ == 1
        for i in 1:3
            a = i
            b = mod(i, 3) + 1
            c = mod(i + 1, 3) + 1
            cc[a, a] = -(1.0 + dd[a]) * (3.0 * dd[a]^2 - 2.0 * (dd[b] + dd[c]) * dd[a] + dd[b] * dd[c]) *
                       ((dd[b] - dd[c]) * dd[b] * dd[c])^2
            cc[b, a] = -dd[a] * (1.0 + dd[b]) * (dd[c] - dd[a]) * ((dd[b] - dd[c]) * dd[b] * dd[c])^2
            cc[c, a] = dd[a] * (1.0 + dd[c]) * (dd[a] - dd[b]) * ((dd[b] - dd[c]) * dd[b] * dd[c])^2
            cc[4, a] = -(dd[a] - dd[b]) * (dd[c] - dd[a]) * ((dd[b] - dd[c]) * dd[b] * dd[c])^2
            bb[a] = cc[4, a] * dd[a]
            ff = ff * (1.0 + dd[a]) / (dd[a] * (dd[a] - dd[b]))^2
        end
        bb[4] = -dd[1] * dd[2] * dd[3] * ((dd[1] - dd[2]) * (dd[2] - dd[3]) * (dd[3] - dd[1]))^2
        ff = -ff / 6.0

        for i in 1:4
            Ans = Ans + F[i] * (cc[i, 1] * ll[1] + cc[i, 2] * ll[2] + cc[i, 3] * ll[3] + bb[i])
        end
        Ans = Ans * ff / (D[4] + hw)

        return Ans * Det_t

    elseif type_ == 2
        # calculate integ d2k / |grad D| * F(k)

        # obtaining |grad D|
        t_inverse, Det_t = _tet_inv3(t)
        t_inverse = t_inverse ./ Det_t

        GradD = 0.0
        for i in 1:3
            for j in 1:3
                for k in 1:3
                    GradD = GradD + t_inverse[i, k] * t_inverse[j, k] *
                                    (D[i+1] - D[1]) * (D[j+1] - D[1])
                end
            end
        end
        GradD = sqrt(abs(GradD))

        x = zeros(Float64, 3)
        # F_uv is DIMENSION(0:2) in Fortran (0-based); map F_uv(0),F_uv(1),F_uv(2)
        F_uv0 = 0.0
        F_uv1 = 0.0
        F_uv2 = 0.0

        if hw < D[1]
            Ans = 0.0
        elseif hw < D[2]
            # parametrization
            x[1] = (hw - D[1]) / (D[2] - D[1])
            x[2] = (hw - D[1]) / (D[3] - D[1])
            x[3] = (hw - D[1]) / (D[4] - D[1])
            # Jacobian factor
            Jac = tet_jacobian(t, x, 1)
            # integration formula
            F_uv0 = F[1] + (F[2] - F[1]) * x[1]
            F_uv1 = (F[3] - F[1]) * x[2] - (F[2] - F[1]) * x[1]
            F_uv2 = (F[4] - F[1]) * x[3] - (F[2] - F[1]) * x[1]
            Ans = Jac * (F_uv0 / 2.0 + (F_uv1 + F_uv2) / 6.0) / GradD
        elseif hw < D[3]
            # parametrization
            x[1] = (hw - D[4]) / (D[2] - D[4])
            x[2] = (hw - D[1]) / (D[3] - D[1])
            x[3] = (hw - D[1]) / (D[4] - D[1])

            ## triangle 1
            Jac = tet_jacobian(t, x, 2)
            F_uv0 = F[1] + (F[4] - F[1]) * x[3]
            F_uv1 = (F[2] - F[1]) * x[1] + (F[4] - F[1]) * (1.0 - x[1] - x[3])
            F_uv2 = (F[3] - F[1]) * x[2] - (F[4] - F[1]) * x[3]
            Ans = Jac * (F_uv0 / 2.0 + (F_uv1 + F_uv2) / 6.0)

            ## triangle 2
            y = (hw - D[3]) / (D[2] - D[3])
            Jac = tet_jacobian(t, x, 3)
            F_uv0 = F[1] + (F[2] - F[1]) * y + (F[3] - F[1]) * (1.0 - y)
            F_uv1 = (F[2] - F[1]) * (x[1] - y) +
                    (F[3] - F[1]) * (y - 1.0) + (F[4] - F[1]) * (1.0 - x[1])
            F_uv2 = -(F[2] - F[1]) * y + (F[3] - F[1]) * (y - 1.0 + x[2])
            Ans = Ans + Jac * (F_uv0 / 2.0 + (F_uv1 + F_uv2) / 6.0)

            Ans = Ans / GradD
        elseif hw < D[4]
            # parametrization
            x[1] = (hw - D[4]) / (D[2] - D[4])
            x[2] = (hw - D[4]) / (D[3] - D[4])
            x[3] = (hw - D[1]) / (D[4] - D[1])
            Jac = tet_jacobian(t, x, 4)
            F_uv0 = F[1] + (F[4] - F[1]) * x[3]
            F_uv1 = (F[2] - F[1]) * x[1] + (F[4] - F[1]) * (1.0 - x[1] - x[3])
            F_uv2 = (F[3] - F[1]) * x[2] + (F[4] - F[1]) * (1.0 - x[2] - x[3])
            Ans = Jac * (F_uv0 / 2.0 + (F_uv1 + F_uv2) / 6.0) / GradD
        else
            Ans = 0.0
        end
        return Ans

    elseif type_ == 3
        for i in 1:3
            a = i
            b = mod(i, 3) + 1
            c = mod(i + 1, 3) + 1
            cc[a, a] = -(1.0 + dd[a]) * (2.0 * dd[a]^3 + (3.0 - dd[b] - dd[c]) * dd[a]^2
                                         - 2.0 * (dd[b] + dd[c]) * dd[a] + dd[b] * dd[c]) *
                       ((dd[b] - dd[c]) * dd[b] * dd[c])^2
            cc[b, a] = -(1.0 + dd[a]) * dd[a] * (1.0 + dd[b]) * (dd[c] - dd[a]) *
                       ((dd[b] - dd[c]) * dd[b] * dd[c])^2
            cc[c, a] = (1.0 + dd[a]) * dd[a] * (1.0 + dd[c]) * (dd[a] - dd[b]) *
                       ((dd[b] - dd[c]) * dd[b] * dd[c])^2
            cc[4, a] = -(1.0 + dd[a]) * (dd[a] - dd[b]) * (dd[c] - dd[a]) *
                       ((dd[b] - dd[c]) * dd[b] * dd[c])^2
            bb[a] = cc[4, a] * dd[a]
            ff = ff * (1.0 + dd[a]) / (dd[a] * (dd[a] - dd[b]))^2
        end
        bb[4] = -dd[1] * dd[2] * dd[3] * ((dd[1] - dd[2]) * (dd[2] - dd[3]) * (dd[3] - dd[1]))^2
        ff = ff / 2.0
        for i in 1:4
            Ans = Ans + F[i] * (cc[i, 1] * ll[1] + cc[i, 2] * ll[2] + cc[i, 3] * ll[3] + bb[i])
        end
        Ans = Ans * ff / (D[4] + hw)^2

        return Ans * Det_t
    else
        return 0.0
    end
end

# ---------------------------------------------------------------------------
# tetrahedron_fermidirac  ->  tet_fermidirac(F, E_ref, E2, t, hw, Ef, type_, cutoff, avoid_deg)
# The 5-case Fermi-level split (lines 169-267).
# ---------------------------------------------------------------------------
function tet_fermidirac(F::Vector{Float64}, E_ref::Vector{Float64}, E2::Vector{Float64},
                        t::Matrix{Float64}, hw::Float64, Ef::Float64, type_::Int,
                        tet_cutoff::Float64, avoid_deg::Float64)::Float64
    Ans = 0.0
    # sorting vertices according to E_ref
    F_s = copy(F)
    E1_s = copy(E_ref)
    E2_s = copy(E2)
    t_s = copy(t)
    E1_s, E2_s, F_s, t_s = tet_sort!(E1_s, E2_s, F_s, t_s)
    D = E1_s .- E2_s

    x = zeros(Float64, 3)
    F_small = zeros(Float64, 4)
    D_small = zeros(Float64, 4)
    t_small = Matrix{Float64}(undef, 3, 3)

    # case 1,2,3,4,5
    if Ef < E1_s[1]                 # case 1: zero
        Ans = Ans + 0.0
    elseif Ef < E1_s[2]             # case 2: a small tet.

        x[1] = (Ef - E1_s[1]) / (E1_s[2] - E1_s[1])
        x[2] = (Ef - E1_s[1]) / (E1_s[3] - E1_s[1])
        x[3] = (Ef - E1_s[1]) / (E1_s[4] - E1_s[1])

        F_small[1] = F_s[1]
        D_small[1] = D[1]
        for i in 1:3
            F_small[i+1] = F_s[1] + (F_s[i+1] - F_s[1]) * x[i]
            D_small[i+1] = D[1] + (D[i+1] - D[1]) * x[i]
            t_small[:, i] = t_s[:, i] .* x[i]
        end
        Ans = Ans + tet_integral(F_small, D_small, t_small, hw, type_, tet_cutoff, avoid_deg)

    elseif Ef < E1_s[3]             # case 3: two tet.'s with cases 2 and 4

        x[1] = (Ef - E1_s[4]) / (E1_s[2] - E1_s[4])
        x[2] = (Ef - E1_s[1]) / (E1_s[3] - E1_s[1])
        x[3] = (Ef - E1_s[1]) / (E1_s[4] - E1_s[1])
        y = (Ef - E1_s[3]) / (E1_s[2] - E1_s[3])

        F_small = copy(F_s)
        D_small = copy(D)
        t_small = copy(t_s)
        F_small[4] = F_s[1] + (F_s[4] - F_s[1]) * x[3]
        D_small[4] = D[1] + (D[4] - D[1]) * x[3]
        t_small[:, 3] = t_s[:, 3] .* x[3]
        Ans = Ans + tet_integral(F_small, D_small, t_small, hw, type_, tet_cutoff, avoid_deg)

        F_small[1] = F_s[1] + (F_s[3] - F_s[1]) * x[2]
        F_small[2] = F_s[3] + (F_s[2] - F_s[3]) * y
        D_small[1] = D[1] + (D[3] - D[1]) * x[2]
        D_small[2] = D[3] + (D[2] - D[3]) * y
        t_small[:, 1] = t_s[:, 1] .* y .+ t_s[:, 2] .* (1 - y - x[2])
        t_small[:, 2] = t_s[:, 2] .* (1 - x[2])
        t_small[:, 3] = t_s[:, 3] .* x[3] .- t_s[:, 2] .* x[2]
        Ans = Ans - tet_integral(F_small, D_small, t_small, hw, type_, tet_cutoff, avoid_deg)

        F_small[1] = F_s[4] + (F_s[2] - F_s[4]) * x[1]
        F_small[3] = F_small[2]
        F_small[2] = F_s[2]
        D_small[1] = D[4] + (D[2] - D[4]) * x[1]
        D_small[3] = D_small[2]
        D_small[2] = D[2]
        t_small[:, 1] = (t_s[:, 1] .- t_s[:, 3]) .* (1 - x[1])
        t_small[:, 2] = t_s[:, 1] .* (y - x[1]) .+ t_s[:, 2] .* (1 - y) .+ t_s[:, 3] .* (x[1] - 1)
        t_small[:, 3] = -t_s[:, 1] .* x[1] .+ t_s[:, 3] .* (x[1] + x[3] - 1)
        Ans = Ans + tet_integral(F_small, D_small, t_small, hw, type_, tet_cutoff, avoid_deg)

    elseif Ef < E1_s[4]             # case 4: a large tet. - a small tet.

        x[1] = (Ef - E1_s[4]) / (E1_s[2] - E1_s[4])
        x[2] = (Ef - E1_s[4]) / (E1_s[3] - E1_s[4])
        x[3] = (Ef - E1_s[1]) / (E1_s[4] - E1_s[1])

        F_small[1] = F_s[4]
        F_small[2] = F_s[1] + (F_s[4] - F_s[1]) * x[3]
        F_small[3] = F_s[4] + (F_s[2] - F_s[4]) * x[1]
        F_small[4] = F_s[4] + (F_s[3] - F_s[4]) * x[2]
        D_small[1] = D[4]
        D_small[2] = D[1] + (D[4] - D[1]) * x[3]
        D_small[3] = D[4] + (D[2] - D[4]) * x[1]
        D_small[4] = D[4] + (D[3] - D[4]) * x[2]
        t_small[:, 1] = -t_s[:, 3] .* (1 - x[3])
        t_small[:, 2] = (t_s[:, 1] .- t_s[:, 3]) .* x[1]
        t_small[:, 3] = (t_s[:, 2] .- t_s[:, 3]) .* x[2]
        Ans = Ans + tet_integral(F_s, D, t_s, hw, type_, tet_cutoff, avoid_deg) -
              tet_integral(F_small, D_small, t_small, hw, type_, tet_cutoff, avoid_deg)

    else                            # case 5: a large tet.
        Ans = Ans + tet_integral(F_s, D, t_s, hw, type_, tet_cutoff, avoid_deg)
    end

    return Ans
end

# ---------------------------------------------------------------------------
# tetrahedron_spinhall  ->  tet_spinhall(F, E1, E2, t, hw, Ef, type_, cutoff, avoid_deg)
# flag1/flag2 quick-zero, then
#   tet_fermidirac(F, E1, E2, ...) - tet_fermidirac(F, E2, E1, ...)   (lines 125-167).
# ---------------------------------------------------------------------------
function tet_spinhall(F::Vector{Float64}, E1::Vector{Float64}, E2::Vector{Float64},
                      t::Matrix{Float64}, hw::Float64, Ef::Float64, type_::Int,
                      tet_cutoff::Float64, avoid_deg::Float64)::Float64
    # fnk-fmk = 0 then quickly returns zero
    occ1 = zeros(Float64, 4)
    occ2 = zeros(Float64, 4)
    flag1 = true
    flag2 = true
    for i in 1:4
        if E1[i] < Ef
            occ1[i] = 1.0
        end
        if E2[i] < Ef
            occ2[i] = 1.0
        end
        if occ1[i] != 1.0 || occ2[i] != 1.0
            flag1 = flag1 && false
        end
        if occ1[i] != 0.0 || occ2[i] != 0.0
            flag2 = flag2 && false
        end
    end
    if flag1 || flag2
        return 0.0
    else
        return tet_fermidirac(F, E1, E2, t, hw, Ef, type_, tet_cutoff, avoid_deg) -
               tet_fermidirac(F, E2, E1, t, hw, Ef, type_, tet_cutoff, avoid_deg)
    end
end
