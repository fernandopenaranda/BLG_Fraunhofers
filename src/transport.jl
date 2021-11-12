######## Hybrid kernel polynomial method ########
# Functions:

# fraunhofer
#     methods fullKPM, hybridKPM, exact

# Here we compute: <Ic> = <Ic>KPM - <Ic>KPM_proj_ABS + Ic_exact_ABS
#                  <Ic> = Ic_KPM_continuum + Ic_exact_ABS

const σ0τ0 = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
const σ0τz = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1]
const ħoe = ustrip(u"eV/nA",ħ/e)
const a0 = 0.246 

fraunhofer_kpm(Bperprange::StepRangeLen, p; Δϕpoints = 3, kw...) = 
    fraunhofer_kpm(collect(Bperprange), p, collect(-1.:2/Δϕpoints:1); kw...)

"""
"the adaptive maximum value search for CPRs is enabled with the kw `adaptive`.
Disabled by default
"""
fraunhofer_kpm(Bperp::Array{T,1}, p; kw...) where {T} =
    fraunhofer_kpm(Bperp, p, missing; kw...) 
function fraunhofer_kpm(Bperp::Array{T,1}, p, Δϕ::Union{Array{T,1}, Missing}; 
    adaptive = false, abssubspace = missing, kw...) where {T}
    println("computing current matrix...")
    Bperp, ic = icϕ_kpm(Bperp, p, Δϕ, abssubspace; adaptive = adaptive, kw...)
    if adaptive == false
        icmax = similar(Bperp)
        for i in 1:length(Bperp)
            println("it: ", i/length(Bperp))
            icmax[i] = maximum(abs.(ic[i,:])) 
        end
        return Bperp, icmax, ic
    else 
        return Bperp, ic, abs.(ic)
    end
end

"computes the fraunhofer pattern only considering the states inside the parent gap using 
exact diagonalization"
fraunhofer_abs_exact(Bperp::Array{T,1}, p; kw...) where {T} = 
    fraunhofer_abs_exact(Bperp, p, missing; kw...) 
function fraunhofer_abs_exact(Bperp::Array{T,1}, p,
        Δϕ::Union{Array{T,1}, Missing}; kw...) where {T}
    println("computing current matrix...")
    icmax = SharedArray(similar(Bperp))
    ic = icϕ_exactdiag(Bperp, p, Δϕ; kw...)
    for i in 1:length(Bperp)
        println("it: ", i/length(Bperp))
        icmax[i] = maximum(abs.(ic[:, i]))
    end
    return Bperp, icmax, ic
end

"""
    `icϕ_exactdiag(Bperplist::Array{T,1}, p, Δϕ::Union{Array{T,1}, Missing}; kw...)`
supercurrent sweep with both Bperp and Δϕ. See: `supercurrent_exactdiag()`
"""
function icϕ_exactdiag(Bperplist::Array{T,1}, p, 
    Δϕ::Union{Array{T,1}, Missing}; kw...) where {T}
    Δϕlist = ifelse(isa(Δϕ, Missing) == true, -0.5:0.25:π+0.5, Δϕ)
    I = zeros(Float64, length(Δϕlist), length(Bperplist))  
    [I[:, i] = supercurrent_exactdiag(collect(Δϕlist), 
        reconstruct(p, B_perp = Bperplist[i]); kw...) 
            for i in 1:length(Bperplist)]
    return I 
end

"""
    `supercurrent_exactdiag(Δϕlist, p = Params() ; nev = 10)` 
Computes the supercurrent for a fixed Bperp value specified in p. If `nev::Missing` it 
iteratively finds all ABS (i.e. with energies inside a lengthscale which we set us the gap).
It returns the supercurrent contribution of these subgap subspace. Note that this subspace
is chosen to always capture all the ABS contribution although it may contain also continuum
states, because the parent gap will be larger that the induced gap. However, once the
subspace is specified, we can always substract this contribution for the full KPM 
contribution. 
    Two different methods implemented 
        `method = :free_energy`
        `method = :current_cut` (`sum(j_ik)` along a cut)
"""
function supercurrent_exactdiag(Δϕlist, p = Params(); nev = 10, 
    method = :free_energy, kw...)
    println("we are using the following method: ", method)
    @unpack Δ = p
    hpar, = paramhams(p; flat = false)
    if method == :free_energy
        f = SharedArray(zeros(Float64, length(Δϕlist)))
        @sync @distributed for i in 1:length(Δϕlist)
            println(i/length(Δϕlist))
            f[i] = f_e(hpar, Δϕlist[i], nev, Δ; kw...) 
        end
        ip = interpolate((Δϕlist,), f, Gridded(Linear()));
        deriv =[Interpolations.gradient(ip, i)[1] for i in Δϕlist]
        return @. 2/ħoe .* deriv
    else
        I = SharedArray(zeros(Float64, length(Δϕlist)))
        @sync @distributed for i in 1:length(Δϕlist) 
            println(i/length(Δϕlist))
             st = cut_st(hpar, Δϕlist[i], nev, Δ) 
            I[i] = supercurrent_cut(p, st, hpar())
        end
        return I     
    end
end


"computes the supercurrent matrix (Bperp, Δϕ) if adaptive = false or the critical
supercurrent vs magnetic field if adaptive = true, using a KPM routine
see: adaptive_max_finder documentation. If adaptive is true and a ketmodel is given
for the ABSsubspace it will use the hybridKPM approach, otherwise it will
run a  KPM calculation over the whole spectrum"
function icϕ_kpm(Bperp::Array{T,1}, p, Δϕ::Union{Array{T,1}, Missing},
        abssubspace::Union{Quantica.KetModel, Missing}; 
        adaptive = false, rk = 10, kw...) where {T}    
    # method = (isa(abssubspace, Missing) ? :fullKPM, :hybridKPM)
    if adaptive == false
        if isa(Δϕ, Missing)
            @warn "provide a valid Δϕ range" 
        else
            # ic = SharedArray(Array{Float64, 2}(undef, length(Bperp), length(Δϕ)))
            # @sync @distributed for j in 1:length(Δϕ) 
            ic = Array{Float64, 2}(undef, length(Bperp), length(Δϕ))
            for j in 1:length(Δϕ) 
                println("phase: ", Δϕ[j])
                [ic[i, j] = supercurrent_kpm(rk, abssubspace, 
                    reconstruct(p, B_perp = Bperp[i]), Δϕ = Δϕ[j]; kw...)[2]
                 for i in 1:length(Bperp)]
            end
        end
        return Bperp, ic
    else
        icmax = SharedArray(similar(Bperp))
        # local definition of generating functions for the adaptive calculation path
        gen_model(x::Number; it = 1, kw...) =
            supercurrent_kpm(rk, abssubspace, 
                reconstruct(p, B_perp = Bperp[it]), Δϕ = x; kw...)[2]
            # /rk comes from normalization of the stochastic trace
        gen_model(x::Array; kw...) = 
            [gen_model(x[i]; kw ...) for i in 1:length(x)]
        @sync @distributed for i in 1:length(Bperp) 
            @unpack bandrange = kw
            brange = _bandrange(reconstruct(p, B_perp = Bperp[i]), bandrange)
            kw1 = [key=> ifelse(key === :bandrange, brange, value) for (key, value) in kw] 
            # update kw with bandrange kwarg so there is no need to compute it for each ϕ, 
            # it only varies with Bperp
            icmax[i] = adaptive_max_finder(gen_model, fourier_model, it = i; kw1...)[1]
        end
        # I could have used pmap for parallelization but since each task requires
        # the same effort
        return Bperp, icmax
    end
end

"""
it should be called icϕ_continuum_kpm
"""
function icϕ_hybrid(Bperp::Array{T,1}, p, Δϕ::Union{Array{T,1}, Missing};
    adaptive = false, nev = 10, rk = 10, bandrange = missing, kw...) where {T}  
    method =(isa(nev, Missing) ? :adaptive_nev : :nev)
    @unpack Δ = p
    if adaptive == false
        if isa(Δϕ, Missing)
            @warn "provide a valid Δϕ range" 
        else 
            ic_continuum = Array{Float64, 2}(undef, length(Bperp),
            length(Δϕ))
            @sync @distributed for i in 1:length(Bperp)
                p = reconstruct(p, B_perp = Bperp[i])
                hpar, = paramhams(p; flat = false)
                brange = _bandrange(p, bandrange)
                for j in 1:length(Δϕ)
                    println("phase: ", Δϕ[j])
                    st = eigen_st(hpar, Δϕ[j], nev, Δ, 
                        method = method)[2]
                    abssubspace = subspace_kets(st, hpar())
                    ic_continuum[i, j] = 
                        supercurrent_kpm(rk, abssubspace, p, Δϕ = Δϕ[j], 
                        bandrange = brange; kw...)[2]
                end
            end
        end
    else 
        
        icmax = SharedArray(similar(Bperp))
        # local definition of generating functions for the adaptive
        # calculation path
        gen_model(x::Number; it = 1, kw...) =
            supercurrent_kpm(rk, abssubspace, 
                reconstruct(p, B_perp = Bperp[it]), Δϕ = x; kw...)[2]
            # /rk comes from normalization of the stochastic trace
        gen_model(x::Array; kw...) = 
            [gen_model(x[i]; kw ...) for i in 1:length(x)]
        @sync @distributed for i in 1:length(Bperp) 
            @unpack bandrange = kw
            brange = _bandrange(reconstruct(p, B_perp = Bperp[i]), 
                bandrange)
            kw1 = [key=> ifelse(key === :bandrange, brange, value) 
                for (key, value) in kw] 
            # update kw with bandrange kwarg so there is no need 
            # to compute it for each ϕ, 
            # it only varies with Bperp
            icmax[i] = adaptive_max_finder(gen_model, fourier_model,
            it = i; kw1...)[1]
        end
        # I could have used pmap for parallelization but since each
        # task requires the same effort
        # it wouldn't
        return Bperp, icmax
    end

    return Bperp, ic_continuum
end


"""
Computes the supercurrent for a hamiltonian built with parameters `p::Params` using `rk` 
randomvectors with non-zero values only in the right lead (this results in a better
stochastic trace estimation with less resources). If `rk` given by the user do not ensure
the convergence it will be overwritten see: stochastictracecriteria. It uses distributed
with a loop over the phases instead of intel mklsparse
""" 
function supercurrent_kpm(rk, abssubspace, p = Params(); Δϕ = 1,
    bandrange = missing, kw...)
    @unpack Δ, Ls, scale = p
    # bandranges = _bandrange(p, bandrange)
    hfunc, dhfunc = paramhams(p; flat = false)
    # rk = criteria(p, rk, dhfunc) 
    if isa(abssubspace, Missing)
        return [Δϕ, av(π*Δϕ, Δ, bandrange, hfunc, dhfunc; kw...)]
    else 
        return [Δϕ, 
        av(π*Δϕ, Δ, bandrange, hfunc, dhfunc, rk = rk; kw...) - 
        av(π*Δϕ, Δ, bandrange, hfunc, dhfunc, abssubspace; kw...)]
    end
end

"""
computes the thermodynamic average of the current operator using KPM. If
abssubspace is not an argument it computes the whole KPM, else it 
substracts the contribution from ABS states
"""
av(ϕ, Δ, bandranges, hfunc, dhfunc; rk = 10, kw...) = 2/ħoe *  
    real(averageKPM(hfunc(phi = ϕ),dhfunc(phi = ϕ); flat = Val(false), 
        ket = ketsmodel(rk), bandrange = bandranges, kw...))/rk
av(ϕ, Δ, bandranges, hfunc, dhfunc, abssubspace::Vector; kw...) = 
    sum(x -> av(ϕ, Δ, bandranges, hfunc, dhfunc, x; kw...), abssubspace)
av(ϕ, Δ, bandranges, hfunc, dhfunc, abssubspace::Quantica.Ket; kw...) = 
    2/ħoe * real(averageKPM(hfunc(phi = ϕ),dhfunc(phi = ϕ); 
    flat = Val(false), ket = abssubspace, bandrange = bandranges, kw...))

"""
random ket model for the stochastic evaluation of the two KPM traces
"""
ketsmodel(numkets) = randomkets(numkets,  r -> cis(2π * rand()) * σ0τ0)
"""
    subspace_kets(sp::Quantica.Spectrum, h::Quantica.Hamiltonian)
builds a ket array of the states contained in spectrum that should
be filtered out in a hybrid KPM calculation.
"""
subspace_kets(sp::Quantica.Spectrum, h::Quantica.Hamiltonian) = 
    subspace_kets(sp.states, h)
function subspace_kets(states::Matrix{ComplexF64},
    h::Quantica.Hamiltonian)
    st = [reshape(states[:,1], 4, :) for i in 1:size(states,2)]
    ket_0 = ket(h)
    basis_ket = [ket_0 for i in 1:size(st)[1]]
    [basis_ket[i].amplitudes[j]= st[i][:,j] for i in 1:length(basis_ket)
    for j in 1:length(st[1][1,:])]
    return basis_ket
end

"""
returns the free energy computed using the exact eigenvalues for the exact evaluation of
the supercurrent over a finite set of eigenvectors
""" 
function f_e(hpar::Quantica.ParametricHamiltonian, Δϕ, 
    nev::Union{Int64, Missing}, Δ, dacp = false; kw...)
    method =(isa(nev, Missing) ? :adaptive_nev : :nev)
    if dacp == false
        λ = negative_eigen(hpar, Δϕ, nev, Δ, method = method; kw...)
    else
        λ = negative_eigen_dacp(hpar, Δϕ, gap; maxdeg = 2, kw...)
    end
    f = -sum(λ) #this minus comes from the tanh(1/β) when T->0
    return f
end

"""
returns the free energy computed using the exact eigenvalues and the andreev eigenstates 
associated
"""
function f_e_st(hpar::Quantica.ParametricHamiltonian, Δϕ,
    nev::Union{Int64, Missing}, Δ)
    method =(isa(nev, Missing) ? :adaptive_nev : :nev)
    λ, st = negative_eigen(hpar, Δϕ, nev, Δ, method = method)
    f = -sum(λ) #this minus comes from the tanh(1/β) when T->0
    return f, st 
end

function cut_st(hpar::Quantica.ParametricHamiltonian, Δϕ,
    nev::Union{Int64, Missing}, Δ)
    method =(isa(nev, Missing) ? :adaptive_nev : :nev)
    λ, st = eigen_st(hpar, Δϕ, nev, Δ/2, method = method, 
    negative = true)
    # return convert(Matrix{ComplexF64}, qr(st).Q)  
    # QR factorization to avoid degeneracies 
    return st
end
"""
computes the eigenvalues and eigenvectors withing the gap using an adaptive procedure for 
setting the number of nevs calculated with the shift invert method
"""
function eigen_st(hpar, Δϕ, nev, gap; method = :nev, nev0 = 4, 
        negative = false)
    sp = spectrum(hpar(phi = Δϕ), method =
         ArpackPackage(nev = ifelse(method == :nev, nev, nev0), 
        sigma = -0.001im))
    λ = real(sp.energies)
    st = sp.states
    if method == :adaptive_nev
        if maximum(abs.(λ)) < gap
            nev0 *= 2
            println(nev0)
            if nev0 < 401
                λ, st = eigen_st(hpar, Δϕ, missing, nev0 = nev0,
                    gap, method = :adaptive_nev)
            else @warn "we are missing subgap states!" end
        else nothing end
    else nothing end
    return ifelse(method == :nev, λ, λ[λ .< gap]), 
        ifelse(negative == false, st, st[:,findall(x -> x < 0, λ)])
end

"""
returns the NEGATIVE eigenvalues withing the gap using for the exact supercurrent
calculation an adaptive procedure for setting the number of nevs calculated with the shift
invert method
"""
function negative_eigen(hpar, Δϕ, nev, gap; method = :nev, exact = false, nev0 = 5)
    if exact == true
        sp = spectrum(hpar(phi = Δϕ))
    else
        sp = spectrum(hpar(phi = Δϕ), method = 
            ArpackPackage(nev = ifelse(method == :nev, nev, nev0), 
            sigma = -0.001im))
    end
    λ = real(sp.energies)
    λneg = λ[λ.<=0];
    if method == :adaptive_nev
        if maximum(abs.(λneg)) < gap
            nev0 *= 2
            println(nev0)
            if nev0 < 401
                λneg = negative_eigen(hpar, Δϕ, missing, 
                    nev0 = nev0, gap, method =:adaptive_nev)
            else @warn "we are missing subgap states!" end
        else nothing end
    else nothing end
    return ifelse(method == :nev, λneg, λneg[λneg .< gap])
end

"""
    `negative_eigen_dacp(hpar, Δϕ, gap; maxdeg = 2, kw...)`
returns the negative eigenvalues inside the interval (-a, a) using the DACP method,
which is suited to deal with large matrices
"""
function negative_eigen_dacp(hpar, Δϕ, gap; maxdeg = 2, kw...)
    a = gap * 2 # we select both gap and continuum spectrum 
    eigs_dacp = DACPdiagonaliser(hpar(phi = Δϕ), a, numkets = maxdeg; kw...)
    λ = real(eigs_dacp)
    return λ[λ.<=0];
end

"""
current_cut. I do not set very precisely the boundaries but it doesn't matter since non 
interface sites will have zero matrix elements H_ab
"""
edge_positions(h, myregion) = 
    [pos for pos in sitepositions(h, region = myregion)] 
edge_indices_l(h, scale, a0) = 
    [ind for ind in siteindices(h, region = (r -> - 3 * scale * a0/2 < r[1] < 0))]
edge_indices_r(h, scale, a0) = 
    [ind for ind in siteindices(h, region = (r -> 0 <= r[1] < 9*scale*a0/√3))]

"""
computes the supercurrent between two sites i and j
"""
function supercurrent_cut(p, states, h)
    eigenstates = subspace_kets(states, h)
    dn = h.harmonics[1].dn
    lat_positions = Quantica.allsitepositions(h.lattice)
    I = 0
    for i in 1:length(eigenstates)
        I += supercurrent_e(h, eigenstates[i], dn, lat_positions, p.scale, p.a0)
    end
    return I
end

"""
computes the supercurrent for a given energy state, psi, along a vertical (x = 0) cut 
(this could have been made more efficiently selecting only those sites conected by hoppings
with a hopping selector and a region - future optimization)
"""
function supercurrent_e(h, psi, dn, lat_positions, scale, a0)
    I_e = 0
    for i in edge_indices_l(h, scale, a0) 
        for j in edge_indices_r(h, scale, a0)
           I_e += supercurrent_ij(h, psi, i, j, dn, lat_positions)
        end
    end
    return I_e
end

"computes the supercurrent due to a given state psi across the i - j link"
supercurrent_ij(h, psi, i, j, dn, lat_positions) = 
    im*(psi[j]' * h[dn][i, j]' * σ0τz  * psi[i] -
        psi[i]' * σ0τz * h[dn][i, j] * psi[j])

