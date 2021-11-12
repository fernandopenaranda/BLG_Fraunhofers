#######
# Here we store the parameters corresponding to the non-SQUID-like geometry (`mask = false`)
######

p = Params(scale = 10, Δ=.0, Ln=200, Ls=0, W = 320, U = 0.0005, μSL = .8, μSR=.8,
    μN = 0*-0.0034, V_perp =  0., V_par = 1e-6, B_par = 0, B_perp = 0., Vu = 0., 
    edges = :armchair_zigzag, smooth = false, dir = :x, U_dimer = 0.0, α = 0., smoothμ = false,
    localselfenergy = false, analyticalselfenergy = false, normal=true);

p_λarm  = reconstruct(p, λ = 0.03, edges = :armchair_zigzag) # edge regime
p_λzig  = reconstruct(p, λ = 0.03, edges = :zigzag_armchair) # edge regime

p_λarm_gap  = reconstruct(p, λ = 0.03, α = 0.02, V_par = 0.005, 
    edges = :armchair_zigzag) # gapped regime
p_λzig_gap  = reconstruct(p, λ = 0.03, α = 0.02, V_par = 0.005, 
    edges = :zigzag_armchair) # gapped regime

p_λarm_bulk = reconstruct(p, λ = 0.03, α = 0.02, V_par = 0.005, 
    edges = :armchair_zigzag, μN = .1) # bulk regime
p_λzig_bulk = reconstruct(p, λ = 0.03, α = 0.02, V_par = 0.005, 
    edges = :zigzag_armchair, μN = .1) # bulk regime

p_λarm_gap_ls  = reconstruct(p_λarm_gap, μSL = 0.00001, μSR = 0.00001,
    Ls = 100, Δ = 1, normal = false) # gapped regime
p_λarm_bulk_ls  = reconstruct(p_λarm_bulk, μSL = 8, μSR = 8, 
    Ls = 100, Ln = 0, Δ = 1, normal = false) # bulk regime
p_λarm_ls  = reconstruct(p_λarm, μSL = 0.00001, μSR = 0.00001, 
    Ls = 100, Δ = 1, normal = false) # edge regime
pa_inv_smoothμ = reconstruct(p, smoothμ = true, λ = 0.02, 
    Ls = 0, Ln = 200, W = 200, scale = 20, Δxμ = 5, λμ = 5, λμbarrier = 2, Δμbarrier = -2)

p_test2 = reconstruct(p_λarm_gap_ls, scale = 1, Ls = 1, Ln = 2, W = 5)
p_test4 = reconstruct(p_test2, scale = 1, Ls = 2, Ln = 4, W = 20, μN = .5);
p_test5 = reconstruct(p_test2, scale = 1, Ls = 4, Ln = 20, W = 20, μN = .5);

pa_inv_smoothμ_ls = reconstruct(pa_inv_smoothμ, Ls = 50, μSL = .8,
     μSR = .8, Δ = 0.004, normal = false)
pa_inv_smoothμ_ls_trivgap = reconstruct(pa_inv_smoothμ, smoothμ = false, U = 0.04)
pa_inv_smoothμ_ls_bulkconduction = 
    reconstruct(pa_inv_smoothμ_ls, smoothμ = false, μN = 0.015)

####### Flux to B conversion
Φ0 = 2067.8338 # T*nm^2
function Bprange(p, Φrange) 
    Bperp = collect(Φrange./(p.Ln * p.W)) #conversion from flux to Bfield (size dependent)
    Bcriteria(Bperp, p.scale)
    return Bperp
end

"We need to make sure that we are not in the Hofstadter regime, 
    for the scaled and non scaled system"
function Bcriteria(B::Array, scale; threshold = 50.)
    ħc = ustrip(u"eV*nm",ħ*c_0)
    a0 = 0.246    #nm
    Bmax = findmax(abs.(B))[1]
    l_B = ħc/Bmax
    l_scale = scale * a0    
    println(l_B/ l_scale)
    println(l_B/a0)
    l_B/ l_scale > threshold ? nothing : @warn("lB ≈ a0*scale, decrease scaling")
    l_B/ a0 > threshold ? nothing : @warn("lB ≈ a0, increase BLG area 
        (we are in the Hofstadter regime)")
end