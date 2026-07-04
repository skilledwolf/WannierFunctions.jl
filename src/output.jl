# Wannier90 output-file writers, byte-format compatible with the reference
# (hamiltonian.F90 / plot.F90; see docs/reference-notes/file-formats.md §6-12).
#
# Precision by file (gotcha 2): _hr.dat = F12.6; _tb.dat = E15.8; _band.dat = E16.8.
# H(R) is stored UNDIVIDED by ndegen (gotcha 4). Matrix rows print `j i` with j fast
# (row index) and value ham_r[j,i,irpt].

using Printf

# --- Fortran E-format helper -------------------------------------------------
# Fortran `Ew.d` normalises the mantissa to `0.dddddE±nn` (leading `0.`), unlike
# C/Julia `%E` which gives `d.ddddE±nn`. We reproduce the Fortran convention so
# _tb.dat / _band.dat byte-match the reference.
#
# `fortran_e(x, w, d)` renders `x` in `(Ew.d)`, i.e. field width `w`, `d` mantissa
# digits after the point, exponent always signed with (at least) two digits.
function fortran_e(x::Real, w::Int, d::Int)
    if !isfinite(x)
        return lpad(string(x), w)
    end
    if x == 0.0
        body = "0." * lpad("", d, "0") * "E+00"   # 0.000...E+00
        return lpad(body, w)
    end
    neg = x < 0
    ax = abs(float(x))
    e = floor(Int, log10(ax)) + 1          # so that ax = m * 10^e with 0.1<=m<1
    m = ax / 10.0^e
    # round mantissa to d digits; may carry to 1.0 → renormalise
    mr = round(m; digits=d)
    if mr >= 1.0
        mr /= 10.0
        e += 1
    end
    if mr != 0.0 && mr < 0.1
        mr *= 10.0
        e -= 1
    end
    mant_str = @sprintf("%.*f", d, mr)      # "0.dddddddd"
    frac = mant_str[3:end]                  # digits after "0."
    esign = e < 0 ? "-" : "+"
    eabs = abs(e)
    estr = eabs < 100 ? @sprintf("%02d", eabs) : string(eabs)
    body = (neg ? "-" : "") * "0." * frac * "E" * esign * estr
    return lpad(body, w)
end

# =============================================================================
# _hr.dat  (hamiltonian.F90:631-692)
# =============================================================================
"""
    write_hr(path, num_wann, irvec, ndegen, Hr; header)

Write `seedname_hr.dat`. Layout:
- line 1: header string (list-directed, leading space)
- line 2: num_wann (list-directed)
- line 3: nrpts    (list-directed)
- ndegen block: `(15I5)` — 15 ints per line, width 5, wrapping
- matrix rows, nested `irpt{ i{ j }}`: `(5I5,2F12.6)` = Rx Ry Rz j i Re Im,
  value `Hr[j,i,irpt]` (j fast/row).

`Hr` is (num_wann × num_wann × nrpts), stored undivided by ndegen.
"""
function write_hr(path::AbstractString, num_wann::Integer,
                  irvec::Vector{NTuple{3,Int}}, ndegen::Vector{Int},
                  Hr::Array{ComplexF64,3};
                  header::AbstractString="written by Wannier90.jl")
    nrpts = length(irvec)
    open(path, "w") do io
        # list-directed writes: leading space, matches Fortran write(unit,*)
        println(io, " ", header)
        println(io, "          ", num_wann)
        println(io, "          ", nrpts)
        # ndegen, (15I5)
        for (i, d) in enumerate(ndegen)
            print(io, lpad(d, 5))
            (i % 15 == 0 || i == nrpts) && print(io, "\n")
        end
        # matrix rows
        for irpt in 1:nrpts
            R1, R2, R3 = irvec[irpt]
            for i in 1:num_wann, j in 1:num_wann
                h = Hr[j, i, irpt]
                @printf(io, "%5d%5d%5d%5d%5d%12.6f%12.6f\n",
                        R1, R2, R3, j, i, real(h), imag(h))
            end
        end
    end
    return path
end

"""
    read_hr(path) -> (num_wann, irvec, ndegen, Hr)

Inverse of [`write_hr`](@ref). Parses by tokenizing (never by column slicing):
header line skipped, then num_wann, nrpts, the `(15I5)` ndegen block, then the
matrix rows. Returns `Hr[j,i,irpt] = Re + im*Im` (undivided by ndegen).
"""
function read_hr(path::AbstractString)
    lines = readlines(path)
    idx = 1
    idx += 1                                   # skip header
    num_wann = parse(Int, strip(lines[idx])); idx += 1
    nrpts    = parse(Int, strip(lines[idx])); idx += 1
    # ndegen: read nrpts integers across (15I5) wrapped lines
    ndegen = Int[]
    while length(ndegen) < nrpts
        toks = split(strip(lines[idx]))
        append!(ndegen, parse.(Int, toks))
        idx += 1
    end
    length(ndegen) == nrpts ||
        error("read_hr: parsed $(length(ndegen)) ndegen, expected $nrpts")
    irvec = Vector{NTuple{3,Int}}(undef, nrpts)
    Hr = zeros(ComplexF64, num_wann, num_wann, nrpts)
    for irpt in 1:nrpts
        for i in 1:num_wann, j in 1:num_wann
            toks = split(strip(lines[idx])); idx += 1
            R1 = parse(Int, toks[1]); R2 = parse(Int, toks[2]); R3 = parse(Int, toks[3])
            jj = parse(Int, toks[4]); ii = parse(Int, toks[5])
            re = parse(Float64, toks[6]); im = parse(Float64, toks[7])
            (j == 1 && i == 1) && (irvec[irpt] = (R1, R2, R3))
            Hr[jj, ii, irpt] = complex(re, im)
        end
    end
    return num_wann, irvec, ndegen, Hr
end

# =============================================================================
# _tb.dat  (hamiltonian.F90:862-994)
# =============================================================================
"""
    write_tb(path, lattice, num_wann, irvec, ndegen, Hr; header)

Write `seedname_tb.dat`: lattice + `<0n|H|Rm>` + `<0n|r|Rm>`.
- line 1: header (list-directed)
- lines 2-4: a1,a2,a3 in Å (list-directed) — rows are lattice vectors
- num_wann, nrpts (list-directed)
- ndegen `(15I5)`
- H part: per R, a blank line then `(3I5)` irvec, then rows `(2I5,3x,2(E15.8,1x))`
  = j i Re(H) Im(H), value `Hr[j,i,irpt]`
- r part: same block structure; the position operator `<0n|r|Rm>` is written as
  zeros (TODO: r-matrices not yet computed) but the block layout is valid.

`lattice.A` has lattice vectors as COLUMNS (Å); the reference stores
`real_lattice(k,:) = a_k` (rows), so we emit the columns of `lattice.A`.
"""
function write_tb(path::AbstractString, lattice, num_wann::Integer,
                  irvec::Vector{NTuple{3,Int}}, ndegen::Vector{Int},
                  Hr::Array{ComplexF64,3};
                  header::AbstractString="written by Wannier90.jl")
    nrpts = length(irvec)
    A = lattice.A          # columns are a_1,a_2,a_3 (Å)
    open(path, "w") do io
        println(io, " ", header)
        # a_1, a_2, a_3 (each is a column of A). List-directed E-format-ish output;
        # reference uses write(*,*) which prints full precision. Use E15.8-style.
        for k in 1:3
            @printf(io, "  %s  %s  %s\n",
                    fortran_e(A[1, k], 15, 8),
                    fortran_e(A[2, k], 15, 8),
                    fortran_e(A[3, k], 15, 8))
        end
        println(io, "          ", num_wann)
        println(io, "          ", nrpts)
        for (i, d) in enumerate(ndegen)
            print(io, lpad(d, 5))
            (i % 15 == 0 || i == nrpts) && print(io, "\n")
        end
        # H part
        for irpt in 1:nrpts
            R1, R2, R3 = irvec[irpt]
            print(io, "\n")                    # blank line before each R block
            @printf(io, "%5d%5d%5d\n", R1, R2, R3)
            for i in 1:num_wann, j in 1:num_wann
                h = Hr[j, i, irpt]
                @printf(io, "%5d%5d   %s %s \n", j, i,
                        fortran_e(real(h), 15, 8), fortran_e(imag(h), 15, 8))
            end
        end
        # r part: <0n|r|Rm>. TODO: r-matrices (position operator) not yet computed;
        # written as zeros so downstream parsers see a structurally valid block.
        for irpt in 1:nrpts
            R1, R2, R3 = irvec[irpt]
            print(io, "\n")
            @printf(io, "%5d%5d%5d\n", R1, R2, R3)
            for i in 1:num_wann, j in 1:num_wann
                z = fortran_e(0.0, 15, 8)
                @printf(io, "%5d%5d   %s %s %s %s %s %s \n", j, i,
                        z, z, z, z, z, z)
            end
        end
    end
    return path
end

# =============================================================================
# _band.dat  (plot.F90:1157-1172)
# =============================================================================
"""
    write_band_dat(path, xvals, energies)

Write `seedname_band.dat`. `energies` is (nb × nk). For each band b, for each
k-point, a line `(2E16.8)` = xval energy, with a blank line between bands.
"""
function write_band_dat(path::AbstractString, xvals::AbstractVector,
                        energies::AbstractMatrix)
    nb, nk = size(energies)
    length(xvals) == nk ||
        error("write_band_dat: length(xvals)=$(length(xvals)) != nk=$nk")
    open(path, "w") do io
        for b in 1:nb
            for k in 1:nk
                print(io, fortran_e(xvals[k], 16, 8),
                          fortran_e(energies[b, k], 16, 8), "\n")
            end
            b < nb && print(io, " \n")         # blank line between bands
        end
    end
    return path
end

# =============================================================================
# _band.kpt  (plot.F90:712-717)
# =============================================================================
"""
    write_band_kpt(path, kpts; weight=1.0)

Write `seedname_band.kpt`: first line the count, then `(3f12.6,3x,a)` =
k1 k2 k3 weight (weight rendered as "1.0" by default, matching the reference).
"""
function write_band_kpt(path::AbstractString, kpts::Vector{<:AbstractVector};
                        weight="1.0")
    open(path, "w") do io
        println(io, lpad(length(kpts), 12))
        for k in kpts
            @printf(io, "%12.6f%12.6f%12.6f   %s\n",
                    k[1], k[2], k[3], string(weight))
        end
    end
    return path
end

# =============================================================================
# _band.labelinfo.dat  (plot.F90:721-741)
# =============================================================================
"""
    write_labelinfo(path, labels, idxs, xvals, kpts)

Write `seedname_band.labelinfo.dat`. One line per special point, format
`(a,3x,I10,3x,4f18.10)` = label, point-index, xval, k1 k2 k3. The reference
label variable is `character(len=20)`, so the label field is left-justified in
a 20-char field (matches the reference byte spacing).
"""
function write_labelinfo(path::AbstractString, labels::AbstractVector,
                         idxs::AbstractVector{<:Integer},
                         xvals::AbstractVector,
                         kpts::Vector{<:AbstractVector};
                         label_width::Int=20)
    n = length(labels)
    (length(idxs) == n && length(xvals) == n && length(kpts) == n) ||
        error("write_labelinfo: mismatched lengths")
    open(path, "w") do io
        for i in 1:n
            k = kpts[i]
            # '(a,3x,I10,3x,4f18.10)': label in a len=20 field, then index, xval, coords.
            @printf(io, "%s   %10d   %18.10f%18.10f%18.10f%18.10f\n",
                    rpad(string(labels[i]), label_width),
                    idxs[i], xvals[i], k[1], k[2], k[3])
        end
    end
    return path
end
