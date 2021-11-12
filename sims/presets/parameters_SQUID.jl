#######
# Here we store the parameters corresponding to the SQUID-like geometry (`mask = true`)
######

####### Flux to B conversion
Φ0 = 2067.8338 # T*nm^2
function Bprange(p, Φrange) 
    Bperp = collect(Φrange./(p.Ln * p.W)) #conversion from flux to Bfield (size dependent)
    Bcriteria(Bperp, p.scale)
    return Bperp
end

"We need to make sure that we are not in the Hofstadter regime, for the scaled and
non-scaled system"
function Bcriteria(B::Array, scale; threshold = 50.)
    ħc = ustrip(u"eV*nm",ħ*c_0)
    a0 = 0.246    #nm
    Bmax = findmax(abs.(B))[1]
    l_B = ħc/Bmax
    l_scale = scale * a0    
    println(l_B/ l_scale)
    println(l_B/a0)
    l_B/ l_scale > threshold ? nothing : @warn("lB ≈ a0*scale, decrease scaling")
    l_B/ a0 > threshold ? nothing : @warn("lB ≈ a0, increase BLG area (we are in t
        he Hofstadter regime)")
end

#################
# Params
#################
p = Params(scale = 10, Δ=.0, Ln=200, Ls=0, W = 320, U = 0.00005, μSL = .8, μSR=.8,
    μN = 0*-0.0034, V_perp =  0., V_par = 1e-6, B_par = 0, B_perp = 0., Vu = 0.,
    edges = :armchair_zigzag, smooth = false, dir = :x, U_dimer = 0.0, α = 0.,
    smoothμ = false, localselfenergy = false, analyticalselfenergy = false, normal=true);

p_inv = reconstruct(p, λ = 0.02,
    scale = 10, Ln = 70, W = 74,
    smoothμ = false, mask = true, Δx_mask = 7, Δy_mask = 5,
    analyticalselfenergy = false, normal = true)

p_inv_sc = reconstruct(p_inv,
    localselfenergy = true, normal = false, Δ = 0.1)

p_inv_sc_smooth = reconstruct(p_inv_sc, 
    smoothμ = true, Δxμ = 7, λμ = 1, λμbarrier = 1, Δμbarrier = -3)


    ############################################################
p_norm_sc_smooth = reconstruct(p_inv_sc_smooth, 
    Δx_mask = 14,  μN = -0.34)

##############################
p_inv_2 = reconstruct(p, λ = 0.02,
    scale = 10, Ln = 50, W = 50,
    smoothμ = false, mask = true, Δx_mask = 7, Δy_mask = 5,
    analyticalselfenergy = false, normal = true)

p_inv_sc_2 = reconstruct(p_inv_2,
    localselfenergy = true, normal = false, Δ = 0.1)



####################

p_inv_3 = reconstruct(p, λ = 0.02,
    scale = 20, Ln = 50, W = 80,
    smoothμ = false, mask = true, Δx_mask = 5, Δy_mask = 8,
    analyticalselfenergy = false, normal = true)

p_inv_sc_3 = reconstruct(p_inv_3,
    localselfenergy = true, normal = false, Δ = 0.1)


#####################

p_nodeg = reconstruct(p_inv_sc_smooth, V_par = 1e-4, α = 1e-4, Ex = 1e-4, Ey = 1e-4, 
    Ez = 1e-4)