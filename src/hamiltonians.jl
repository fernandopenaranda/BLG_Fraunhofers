#####################################################################################
# Hamiltonians 
#####################################################################################

"Build the unbounded superconducting hamiltonian (leads)"
function hS(p = Params())
    modelS = models_device(p).modelS
    lat = lat_infinite(p)
    return hamiltonian(lat, modelS, orbitals= (:eup,:edown,:hup,:hdown))
end

function nanoribbon_irregular_edge(p = Params())
    @unpack edges, normal = p
    m = models_device(p)
    hterms =  m.modelNtop + m.modelNbot + m.model0
    lat = unitcell(lat_infinite(p), (1, 3); region = r -> 0 < r[2] < 20)
    return lat |> hamiltonian(hterms, orbitals=  normal ? (:eup, :edown) : (:eup,:edown,:hup,:hdown)) |> 
        parametric(m.onsiteNtop_ribbon!, m.onsiteNbot_ribbon!) 
end

" Constructor returns a SNS junction y infinite (infinite SC-N armchair interface)
and a normal unbounded parametric hamiltonian as 2D ham in 3D space. 
IMPORTANT SMOOTH MUST BE SET TO FALSE or dir = :y "
function nanoribbonperiodic(p = Params())
    @unpack normal = p
    lat = lat_infinite(p)
    m = models_device(p)
    devicereg, = xyregs(p)
    br1, br2, br3 = nanoribbonbravais(p)
    hterms =  m.modelNtop + m.modelNbot + m.model00
    hamN = lat |> 
        hamiltonian(hterms, orbitals = normal ? (:eup, :edown) : (:eup,:edown,:hup,:hdown))
    phamN = hamN |> parametric(m.rashba!, m.onsiteNtop_inf!, m.onsiteNbot_inf!, 
        m.Δdimer_top!, m.Δdimer_bot!,m.cantingNAtBb!, m.cantingNAbBt!)
    hamSNS = hamN |> unitcell(br1, br2) |> unitcell(br3, region = devicereg)
    phamSNS = hamSNS |> 
        parametric(m.onsiteNtop_ribbon!, m.onsiteNbot_ribbon!, m.field!,
            m.Δdimer_top!,m.Δdimer_bot!, m.cantingNAtBb!, m.cantingNAbBt!, m.onsiteS!,
            m.rashba!, m.interhopSCt1!, m.interhopSCt3!, m.peierlsnoBperp!,
            m.transparency!, m.localselfenergy!, m.analyticalselfenergyon!, 
            m.scphasesnoBperp!, m.analyticalselfenergyhop!)
    return phamSNS, phamN
end

"Build a nanoribon of width W with params P::Params and infinite edge (armchair or zigzag)
 to the vacuum"
function nanoribbon(p = Params())
    @unpack edges, normal = p
    m = models_device(p)
    hterms =  m.modelNtop + m.modelNbot + m.model0
    nanoribreg = xyregs(p)[2]
    lat = unitcell(lat_infinite(p), ifelse(edges == :armchair_zigzag,(-1,1),(1,1)), region = nanoribreg) 
    return lat |> hamiltonian(hterms, orbitals=  normal ? (:eup, :edown) : (:eup,:edown,:hup,:hdown)) |> 
        parametric(m.onsiteNtop_ribbon!, m.onsiteNbot_ribbon!) 
end

"Build the bounded parametric hamiltonian given `p::Params`, free of dangling bonds.
if `mask = true` we create a rectangular hole in the center of the device whose dimensions
are specified by `Δx_mask, Δy_mask`"
function paramhams(p = Params(); flat = true)
    @unpack t0, t1, t3, λ, normal, mask , Δx_mask, Δy_mask, Ln, W, a0, scale = p
    m = models_device(p)
    lat = lat_infinite(p)
    resetonsites! = @onsite!(o -> 0*o)
    resethoppings! = @hopping!(t -> 0*t)
    reg = Quantica.combine(xyregs(p)[1], xyregs(p)[2])
    mask == false ?  devicereg = reg : devicereg = xor(reg, RegionPresets.rectangle((Ln-2*Δx_mask, W-2*Δy_mask), (0., W/2 + a0 *scale)))
    ham = lat |> hamiltonian(m.model0, orbitals =  normal ? (:eup, :edown) : (:eup,:edown,:hup,:hdown)) |>
            unitcell(mincoordination = coordinationfunc(t0,t1,t3,λ), region = devicereg)
            #note when λ = 0 coordinationfunc stops working correctly, put an infinitesimal but finite value
    pham = ham |> 
        parametric(m.onsiteNtop!, m.onsiteNbot!, m.Δdimer_top!, m.Δdimer_bot!, m.cantingNAtBb!,
        m.cantingNAbBt!, m.onsiteS!, m.interhopSCt1!, m.interhopSCt3!, m.peierls!, 
        m.transparency!, m.localselfenergy!, m.analyticalselfenergyon!, m.scphases!,
        m.analyticalselfenergyhop!)
    pdham = ham |> parametric(resetonsites!, resethoppings!, m.pairing!, m.dpeierls!) # old: m.onsiteS!, before dpeierls
    return pham, pdham
end

"Implementation of a monolayer graphene model with Kane-Mele SOC + SC -> Kitaev instead of 
Bilayer + Rashba + Zeeman + SC -> Oreg"
function Fukane end

"Determines the minimum connectivity in the lattice so that no dangling bonds form. It reads
how many nonzero hopping elements there are (up to a given precision) to return the minimum
coordination number in such a way we avoid dangling sites in the corners"
coordinationfunc(t0::Number, t1::Number, t3::Number, λ::Number) = coordinationfunc([t0, t1, t3, λ])  
coordinationfunc(hopvect::Array) = length(hopvect[abs.(hopvect) .> 1e-10]) + 1