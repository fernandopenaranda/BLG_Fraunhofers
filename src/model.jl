#####################################################################################
# Model
#####################################################################################

const ħoec = ustrip(u"T*nm^2",ħ/e)
const μB = ustrip(u"T^-1*meV",μ_B)   
const σ0´ = @SMatrix[1 0; 0 1]
const σx´ = @SMatrix[0 1; 1 0]
const σy´ = @SMatrix[0 -im; im 0]
const σz´ = @SMatrix[1 0; 0 -1]
const σ0τ0´ = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
const σ0τz´ = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1]
const σxτ0´ = @SMatrix[0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]
const σyτ0´ = @SMatrix[0 -im 0 0; im 0 0 0; 0 0 0 -im; 0 0 im 0]
const σzτ0´ = @SMatrix[1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 -1]
const σxτz´ = @SMatrix[0 1 0 0; 1 0 0 0; 0 0 0 -1; 0 0 -1 0]
const σyτz´ = @SMatrix[0 -im 0 0; im 0 0 0; 0 0 0 im; 0 0 -im 0]
const σzτz´ = @SMatrix[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
const σyτy´ = @SMatrix[0 0 0 -1; 0 0 1 0; 0 1 0 0; -1 0 0 0]

"""
Note if localselfenergy == True || analyticalselfenergy == True, LS = 0 as a preset
Note as it is written you should pass missing args corresponding to the fields: a0, dinter, dintra t0   

To do:
  1) fix kfs
"""

# Scaling model: t1 belongs to the low energy sector (high energy sector commented _model1) info:
# https://iopscience.iop.org/article/10.1088/0034-4885/76/5/056503/pdf

@with_kw struct Params @deftype Float64
                                        # Params AB bernal bilayer obtained from DFT calculations:
                                        # params ref -> https://arxiv.org/pdf/1511.06706.pdf (flag <-:dft)
    normal::Bool = false                # determines whether the system is normal or Nambu 
    smooth::Bool = false                # U smooth profile if true otherwise constant
    smoothμ::Bool = false               # Chemical potential profile if true otherwise steplike
    mask::Bool = false                  # Determines whether a rectangular hole is made on the N region
    dir::Symbol = :y                    # Determines the direction in which U is smooth. Possible values: :x, :y or :xy    
    edges::Symbol = :armchair_zigzag    # edges kwarg sets the nature of the SN and Nvacuum interfaces: 
                                        # edges =:typeofxinterface_typeofyinterfaces
                                        # Default :armchair_zigzag. 
                                        # Other options :zigzag_armchair, :zigzag_zigzag
    car::Bool = false                   # Interlayer hoppings in the SC leads are by default suppressed
    localselfenergy::Bool = false       # Determines whether we implement a localselfenergy ∝Δσyτy or not 
    analyticalselfenergy::Bool = false  # Determines whether the self-energy due to a 3d unbounded SC
                                        # (local + non local terms) is implented or not. See notes
    scale = 1.                          # scaling parameter
    a0 = 0.246 
    dinter = 1.36 * a0                  # Note that we leave dinter unaltered by scale
    dintra = a0 /sqrt(3) 
    Ls = 3.                             # Length of the SC leads
    Ln = 3.                             # Length of the N region
    W = 3.                              # Width of the JJ junction
    t0 =  2.6                           # <-:dft 
    t1 = .36                            # <-:dft
    τL = 1.                             # Transparency left
    τR = 1.                             # Transparency right
    t3 = .28                            # <-:dft #WARN number 1) I NEED TO RETHINK THE SCALING
    μN = 0. * t0
    μSL = t0 / 3.25                     # t0/3 > μS > 2 * t1
    μSR = t0 / 3.25                     # t0/3 > μS > 2 * t1
    U_dimer = 0.0* 0.15                 # Energy difference between dimer and non dimer sites: https://arxiv.org/pdf/1901.01332.pdf
    g1_on = 1.3e-3                      # Comparable to the gap
    g1_hop = 1.3e-3                     # Comparable to the gap
    kfs = 1.                            # Fermi wavelength of the superconductor
    U = 0.0*1e-3                        # Interlayer bias. If smooth = true is the maximal value that is reached at the edges,
                                        # else is the homogeneous value
    Δ = 1.3e-3                          # That of SC MoRe: ref ->  https://arxiv.org/abs/2006.05522
    B_par = 0*0.5                       # Applied on the outof plane
    B_perp = 0.                         # In-plane field and perpendicular to the system. 
    λ = 0*0.5e-3                        # Between [-2.5,2.5] meV: ref -> https://arxiv.org/abs/1905.08688
    α = 0.0                             # Rashba coef
    g = 2.                              # g-factor
    V_par = 0.0                         # Zeeman field in y direction
    V_perp = 0.0                        # Zeeman fiel in z direction. You might want to uncomment them and compute it
                                        # with B playing attention to the scaling.
    Vu = 0.                             # modulus of the effective onsite energy Vuvector due to the e-e interactions (spin canting)
    θ = 0.                              # AF canting with respect to the y direction.
    Ex = 0                              # Electric field
    Ey = 0
    Ez = 0
    Δxμ = 10                            # Center of the sigmoidal μSL/R-μN
    λμ = 5                              # Smooth parameter
    λμbarrier = 1                        
    Δμbarrier = -10                     # Height of the leads barrier (tunnel regime) -> increasing normal reflections
    Δx_mask = 2                         # width of the mask on x
    Δy_mask = 2                         # width of the mask on y
end                                     # Units eV,  T, nm

"hamiltonian models given `p::Params`"
function models_device(p = Params())
    @unpack a0, dinter, dintra, Ls, Ln, W, μN, μSL, μSR, t0, t1, τL, τR, normal,
            t3, Δ, B_par, B_perp, λ, g, g1_on, g1_hop, kfs, U, U_dimer, α, V_par, V_perp,
            Vu, θ, smooth, smoothμ, scale, edges, dir, car, localselfenergy, analyticalselfenergy, Ex, Ey, Ez,
            Δxμ, λμ, λμbarrier, Δμbarrier = p

    (σ0τ0, σ0τz, σxτz, σyτz, σzτz, σyτy), σvect = sigmas(Val(normal))#normalorNambu(p)
    
    a0, Ln, Ls, W, regs = modelregs(p)    
    edgereg, Redgereg, normalreg, scRreg, scLreg, screg, screg2 = [regs[i] for i in 1:length(regs)]
    t0 = t0 / scale
    t3 = t3 / scale
    dintra = dintra * scale
    dinter_warp = sqrt(dintra^2 + dinter^2)

    # t0, t1, and t3 hoppings
    modelintra = hopping(-t0 * σ0τz, range = dintra,
            sublats = (:At=>:Bt,:Ab=>:Bb,:Bb=>:Ab,:Bt=>:At))
    modelinter = hopping(0*t1 * σ0τz, range = dinter, sublats = (:At=>:Bb,:Bb=>:At)) +
                 hopping(0*t3 * σ0τz, range = dinter_warp, sublats = (:Bt=>:Ab, :Ab=>:Bt))        

    ## SOC
    # ISING SOC (different sign of λ in different layers) - "Valley Zeeman"
    # modelKM = λ * (hopKM(:At)-hopKM(:Bt)-hopKM(:Ab)+hopKM(:Bb)) 
    # we don't have a KM but an Ising with different lambda on each layer
    ϕ0 = π/5
    hopKM(s) = hopping((r,dr) -> im *
        sign(cos(3*atan(dr[2],dr[1]) - ϕ0)) * ifelse(normal === true, σzτz, σzτ0´), sublats = (s => s),
         range = a0)
    
    # alternatively cos(3*atan(dr[2],dr[1])). Note that in order this to be independent on the
    # lattice orientation phi0 should not coincide with the direction of any NNN in the rotated lattice.
    # The following renormalization of λ comes from the relation between the SOC
    # gap in the k-Hamiltonian and the TB hopping parameter, t2: λ_KM = 3 √3 t2.
    # See: http://www.phys.ufl.edu/~pjh/teaching/phz7427/7427notes/ch6.pdf
    # see also https://journals.aps.org/prb/pdf/10.1103/PhysRevB.82.245412
    # https://journals.aps.org/prb/pdf/10.1103/PhysRevB.100.085412 This is Andor's et al paper 
    
    λ /= (2*3*sqrt(3))
    modelIsing =  hopping((r, dr) ->λ * im * sign(r[3]) *
        sign(cos(3*atan(dr[2],dr[1]) - ϕ0)) * ifelse(normal === true, σzτz, σzτ0´),
        sublats = (:At => :At, :Ab => :Ab, :Bt => :Bt, :Bb => :Bb), range = a0)

    ising! = @hopping!((t, r, dr; λ = λ) -> t + λ * im * sign(r[3]) *
    sign(cos(3*atan(dr[2],dr[1]) - ϕ0)) * ifelse(normal === true, σzτz, σzτ0´);
    sublats = (:At => :At, :Ab => :Ab, :Bt => :Bt, :Bb => :Bb), range = a0)

    # RASHBA SOC (differnt sign of α in different layers)
    hopRashba() = hopping((r,dr) -> sign(r[3])*(dr[2]/dintra * ifelse(normal === true, σxτz, σxτ0´) - dr[1]/dintra * σyτz), 
         range = dintra, sublats = (:At=>:Bt, :Ab=>:Bb, :Bt=>:At, :Bb=>:Ab))
    modelRashba = 2*im/3 * α * hopRashba()
    
    model00 = onsite(0I) + modelintra + modelinter + modelIsing 
    model0 = onsite(0I) + modelintra + modelinter + modelIsing + modelRashba 

    ## ZEEMAN and canting due to e-e interactions
        # Vu is spatially local (Hartree level, no Fock) and enters as an onsite
        # then it does not contribute to the flux
        # V_par = -1/2*g*μB*B_par
        # V_perp = -1/2*g*μB*B_perp # Uncomment this in params and think about the correct scaling

    Vzeem = @SVector[0, V_par, V_perp]
    Vuvect = Rz([0,Vu,0], θ)
    Vvect = Vzeem + Vuvect
    
    ##  UNBOUNDED TERMS

    modelNtop = onsite((- μN + U/2 ) * σ0τz + sum(σvect .* Vzeem), sublats= (:At, :Bt)) +
                onsite(sum(σvect .* Vuvect), sublats= (:At)) -
                onsite(sum(σvect .* Vuvect), sublats= (:Bt))  

    modelNbot = onsite((- μN - U/2) * σ0τz + sum(σvect .* Vzeem), sublats= (:Ab, :Bb)) +
                onsite(sum(σvect .* Vuvect), sublats= (:Bb)) -
                onsite(sum(σvect .* Vuvect), sublats= (:At))  

    modelS = onsite(- μSL * σ0τz - Δ * σyτy) 

    ## Peierls phases

    piecewise(x) = ifelse(edges == :armchair_zigzag, clamp(x, 0*-Ln/2, Ln), 
        clamp(x, -Ln/2 - a0/(2*√3) , Ln/2 + a0/(2*√3)))
    
    A(r, B_par, B_perp) = SA[0, B_perp, B_par]/ħoec * piecewise(r[1] + Ln/2)
    eφ(r, dr, B_par, B_perp) = diagphi(dot(A(r, B_par, B_perp), (dr + Δz(r))))
    diagphi(φ) = normal ? Diagonal(SA[cis(φ), cis(φ)]) : Diagonal(SA[cis(φ), cis(φ), cis(-φ), cis(-φ)])
    eφ_hoppings(t, r, dr, B_par, B_perp) = t * eφ(r, dr, B_par, B_perp) 

    ## Pairing phases
    #We put φSC/2 to the left and -φSC/2 to the right. This is why we have a 4 factor
    eφ_onsite(o, r, φSC, Bpar, Bperp) =
        eφ(r, r + [0, -W/2, 0], Bpar, Bperp) * diagphi(sign(r[1])*φSC/4) * o * diagphi(sign(r[1])*φSC/4)' *
        eφ(r, r + [0,  -W/2, 0], Bpar, Bperp)'
    deφ_onsite(o, r, φSC, Bpar, Bperp) =
        sign(r[1]) * im/4 * eφ(r,  r + [0, -W/2, 0], Bpar, Bperp) * diagphi(sign(r[1])*φSC/4) *
        (σ0τz * o - o * σ0τz) * diagphi(sign(r[1])*φSC/4)' * eφ(r, r + [0,  -W/2, 0], Bpar, Bperp)'
    
    # ALTERNATIVE: Now we put φSC on the left and 0 on the right 
    # eφ_onsite(o, r, φSC, Bpar, Bperp) =
    #     eφ(r, r, Bpar, Bperp) * diagphi(ifelse(sign(r[1])<0, φSC/2, 0)) * o * 
    #     diagphi(ifelse(sign(r[1])<0, φSC/2, 0))' * eφ(r, r, Bpar, Bperp)' 
    # deφ_onsite(o, r, φSC, Bpar, Bperp) =
    #     sign(r[1]) * im/4 * eφ(r, r, Bpar, Bperp) * diagphi(ifelse(sign(r[1])<0, φSC/2, 0)) *
    #     (σ0τz * o - o * σ0τz) * diagphi(ifelse(sign(r[1])<0, φSC/2, 0))' * eφ(r, r, Bpar, Bperp)'

    ## Ripples: length scale of 10*a0
    Random.seed!(5); 
    kx, ky = rand(0:1e-5:2*2π/10,10,2) 
    Δz(r) = [0,0,0*0.5*a0/length(kx) * sum(sin.((kx.*r[1]+ky.*r[2])/a0))] #disabled
    
    # smooth function describing the band bending due to the SC. λN is the width of the tanh,
    # xμ is the center of the sigmoidal. We want xμ to be inside the normal region
     smoothmu(r) = μN/2*(tanh((r[1]+(Ln/2 - Δxμ))/λμ)-tanh((r[1]-(Ln/2 - Δxμ))/λμ)) + 
         ifelse(r[1] < 0, μSL, μSR)/2*(2-tanh((r[1]+(Ln/2 - Δxμ))/λμ)+tanh((r[1]-(Ln/2 - Δxμ))/λμ)) + 
         Δμbarrier * ((sigmoidal(r[1] - Ln/2, λμbarrier, λμbarrier/2) * sigmoidal(-r[1] + Ln/2, λμbarrier, λμbarrier/2)) +
                     (sigmoidal(r[1] + Ln/2, λμbarrier, λμbarrier/2) * sigmoidal(-r[1] - Ln/2, λμbarrier, λμbarrier/2)))

    ## CAR processes on or not. Interlayer hoppings suppressed in the leads if CAR == false
    hopSCcondition(r,dr) = ifelse(screg(r - dr/2) == true || screg(r + dr/2) == true , 
        ifelse(car == false, 0, 1), 1) 
    ## Transparencies between hopping sites at the SC-N left and right interfaces
    function transparency(r, dr, τL, τR)
        if xor(scRreg(r - dr/2), scRreg(r +  dr/2))
            τR
        elseif xor(scLreg(r - dr/2), scLreg(r +  dr/2))
            τL
        else 1.0 end
    end

    ## Different Δdimer energy (top and bottom):
    Δdimer_top! = @onsite!((o,r) -> o + U_dimer * σ0τz; region = normalreg, sublats = (:At))
    Δdimer_bot! = @onsite!((o,r) -> o + U_dimer * σ0τz; region = normalreg, sublats = (:Bb)) 
 
    ###  BOUNDED MODIFIERS 
    onsiteNtop! = @onsite!((o, r; μN = μN, U = U) -> (-ifelse(smoothμ == false, μN, smoothmu(r)) + U/2 * ifelse(smooth == true, 
        smoothU(r, W, Ln, dir),1)) * σ0τz + sum(σvect .* Vzeem);
        region = normalreg, sublats = (:At, :Bt)) 
    onsiteNbot! = @onsite!((o, r; μN = μN, U = U) -> (- ifelse(smoothμ == false, μN, smoothmu(r)) - U/2 * 
        ifelse(smooth == true, smoothU(r, W, Ln, dir),1)) * σ0τz + 
        sum(σvect .* Vzeem); region = normalreg, sublats = (:Ab, :Bb))

    cantingNAtBb! = @onsite!((o) -> o + sum(σvect .* Vuvect); 
        region = normalreg, sublats= (:At,:Bb))
    cantingNAbBt! = @onsite!((o) -> o - sum(σvect .* Vuvect); 
        region = normalreg, sublats= (:Ab,:Bt))
    field! = @onsite!((o, r) -> o + [Ex, Ey, Ez]'r * σ0τz)
    onsiteS! = @onsite!((o,r; Δ=Δ) -> - ifelse(smoothμ == false, ifelse(scLreg(r), μSL, μSR), smoothmu(r)) * σ0τz - Δ  * σyτy; region = screg)
    peierls! =@hopping!((t, r, dr; B_par=B_par, B_perp=B_perp) -> 
        eφ_hoppings(t,r,dr,B_par,B_perp))

    scphases! = @onsite!((o, r; phi=0, B_par=B_par, B_perp=B_perp) ->
         eφ_onsite(o,r,phi,B_par,B_perp);  region = ifelse(localselfenergy == false && 
         analyticalselfenergy == false, screg, edgereg))

    dpeierls! = @onsite!((o, r; phi=0, B_par=B_par, B_perp=B_perp) -> 
        deφ_onsite(o,r,phi,B_par,B_perp); region = screg)

    transparency! = @hopping!((t,r,dr) -> t * transparency(r, dr, τL, τR); range = 4a0)

    pairing! = @onsite!((o, r) -> -Δ * σyτy; region = screg)

    localselfenergy! = @onsite!((o,r) -> ifelse(localselfenergy == false, o, o - Δ * σyτy); region = edgereg)
    # Self energy in the static and Andreev limit:
    Σon() = g1_on * σyτy # The Fermi level is renormalized by E -> E - δμ_N term 
    Σhop(dr) = 0 * σ0τz + g1_hop * sin(dr * kfs)/dr * σyτy 
    # I've neglected the self energy corrections to the normal hopping above
    analyticalselfenergyon! = @onsite!((o,r) -> 
        ifelse(analyticalselfenergy == false, o, o + Σon()); region = edgereg)

    analyticalselfenergyhop! = @hopping!((t, r, dr) -> 
        ifelse(analyticalselfenergy == false, t, t + Σhop(norm(dr))); 
         sublats = (:At =>:At, :Ab=>:Ab, :Bt=>:Bt, :Bb=>:Bb, :At=>:Bb,:Bb=>:At), range = dinter)

    ## Bounded modifiers required for the unbounded phams (note that for infinite y Bperp has to be 0)
    # in order to have a periodic sys.
    onsiteNtop_ribbon! = @onsite!((o, r; μN=μN,U=U, V_par = V_par, V_perp = V_perp) -> 
        (-ifelse(smoothμ == false, μN, smoothmu(r))  + U/2 * ifelse(smooth == true, smoothU(r, W, Ln, dir),1) ) * σ0τz 
        + V_par * ifelse(normal === true, σyτz, σyτ0´) + V_perp * σzτz; sublats = (:At, :Bt))
        #+ sum(σvect .* Vvect); sublats = (:At, :Bt))
    onsiteNbot_ribbon! = @onsite!((o, r; μN=μN,U=U , V_par = V_par, V_perp = V_perp)  -> 
        (-ifelse(smoothμ == false, μN, smoothmu(r)) - U/2 * ifelse(smooth == true, smoothU(r, W, Ln, dir),1) ) * σ0τz 
        + V_par * ifelse(normal === true, σyτz, σyτ0´) + V_perp * σzτz; sublats = (:Ab, :Bb))
        #+ sum(σvect .* Vvect); sublats = (:Ab,:Bb))       
        
    onsiteNtop_inf! = @onsite!((o; μN=μN,U=U) -> o + (- μN + U/2 ) * σ0τz +
        sum(σvect .* Vvect); sublats = (:At, :Bt))
    onsiteNbot_inf! = @onsite!((o; μN=μN,U=U)  -> o + (- μN - U/2 ) * σ0τz +
        sum(σvect .* Vvect); sublats = (:Ab,:Bb))       

    rashba! = @hopping!((t, r, dr; α = α) -> t + 2*im/3 * α * ifelse(screg(r),0.,1.)* 
        sign(r[3])*( dr[2]/dintra * ifelse(normal === true, σxτz, σxτ0´) - dr[1]/dintra * σyτz); 
        range = dintra, sublats = (:At=>:Bt, :Ab=>:Bb, :Bt=>:At, :Bb=>:Ab))
  
    scphasesnoB! = @onsite!((o, r; phi=0) -> eφ_onsite(o,r,phi,0,0); 
        region = ifelse(localselfenergy == false && analyticalselfenergy == false, screg, edgereg))
    scphasesnoBperp! = @onsite!((o, r; phi=0, B_par=B_par) -> 
        eφ_onsite(o,r,phi,B_par,0); region = ifelse(localselfenergy == false && 
        analyticalselfenergy == false, screg, edgereg))
        
    interhopSCt1! = @hopping!((t, r, dr) -> t + hopSCcondition(r,dr) * t3 * σ0τz; 
        range = dinter_warp, sublats = (:Bt=>:Ab, :Ab=>:Bt))
    interhopSCt3! = @hopping!((t, r, dr) -> t + hopSCcondition(r,dr) * t1 * σ0τz; 
        range = dinter, sublats = (:At=>:Bb,:Bb=>:At))         
    peierlsnoBperp! = @hopping!((t, r, dr; B_par=B_par) -> eφ_hoppings(t,r,dr,B_par,0))

    ### Build models N-Tuple
    models = (; model0, model00, onsiteNtop!, onsiteNbot!, onsiteNtop_inf!, onsiteNbot_inf!,
              onsiteS!, scphases!, scphasesnoB!, peierls!, dpeierls!,
              scLreg, scRreg, modelNtop, modelNbot, modelS, transparency!, Δdimer_top!, pairing!,
              Δdimer_bot!,scphasesnoBperp!, peierlsnoBperp!, cantingNAtBb!, cantingNAbBt!,
              onsiteNtop_ribbon!, localselfenergy!, analyticalselfenergyon!,
              analyticalselfenergyhop!, rashba!, ising!, onsiteNbot_ribbon!, 
              interhopSCt1!,  interhopSCt3!, field!)
    return models
end

"shape of the smooth function f(x,y) s.t U(x, y) = U * f(r) if smooth == true"
function smoothU(r, W, Ln, dir)
    ξ = 10
    if dir == :x
        1 - sigmoidal(r[1], Ln/4, ξ) * sigmoidal(-r[1], Ln/4, ξ)
    elseif dir == :y
        1 - sigmoidal(-r[2], -W/4, ξ) * sigmoidal(r[2], 3W/4, ξ) 
    elseif dir == :xy
        1 - sigmoidal(r[1], Ln/4, ξ) * sigmoidal(-r[1], Ln/4, ξ) *
        sigmoidal(-r[2], -W/4, ξ) * sigmoidal(r[2], 3W/4, ξ)
    end
end

sigmoidal(r, length, ξ) = 1/(1 + exp((r-length)/ξ))

"Rotation matrix of an angle θ in around x axis"
Rx(Rvect, θ) = @SMatrix[1 0 0; 0 cos(θ) -sin(θ); 0 sin(θ) cos(θ)] * Rvect
"Rotation matrix of an angle θ in around z axis"
Rz(Rvect, θ) = @SMatrix[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] * Rvect

"Unbounded lattice takes `p:Params`"
function lat_infinite(p = Params())
    @unpack a0, dinter, scale, edges = p
    a0 = a0 * scale
    if edges == :armchair_zigzag
        sAbot = sublat((0.0,-1.0a0/sqrt(3.0),  - dinter/2); name = :Ab)
        sBbot = sublat((0.0, 0.0a0/sqrt(3.0),  - dinter/2); name = :Bb)
        sAtop = sublat((0.0, 0.0a0/sqrt(3.0),  + dinter/2); name = :At)
        sBtop = sublat((0.0, 1.0a0/sqrt(3.0),  + dinter/2); name = :Bt)
        br = a0 * SA[cos(pi/3) sin(pi/3) 0; -cos(pi/3) sin(pi/3) 0]'
    elseif edges == :zigzag_armchair
        sAbot = sublat((-1.0a0/sqrt(3.0),0.0, - dinter/2); name = :Ab)
        sBbot = sublat((0.0a0/sqrt(3.0),0.0,  - dinter/2); name = :Bb)
        sAtop = sublat((0.0a0/sqrt(3.0),0.0,  + dinter/2); name = :At)
        sBtop = sublat((1.0a0/sqrt(3.0),0.0,  + dinter/2); name = :Bt)
        br = a0 * Rz(SA[cos(pi/3) sin(pi/3) 0; -cos(pi/3) sin(pi/3) 0]',π/2)
    end
    return lattice(sAtop, sBtop, sAbot, sBbot; bravais = br)
end

# Be careful below I have changed the names to force BdG symmetry but it is not
# intuitive. Simplify this 
sigmas(::Val{true}) = [σ0´, σ0´, σx´, σy´, σz´, 0*σ0´], @SVector[σx´, σy´, σz´] 
sigmas(::Val{false}) = [σ0τ0´, σ0τz´, σxτz´, σyτz´, σzτz´, σyτy´], 
    @SVector[σxτz´, σyτ0´, σzτz´] 