#####################################################################################
# Regions required to build the Hamiltonians in hamiltonian.jl
#####################################################################################

"New method to Quantica.combine that accepts Regions"
Quantica.combine(regtup::Tuple{Quantica.RegionPresets.Region,Quantica.RegionPresets.Region}) =
    Quantica.combine(regtup[1], regtup[2])
Quantica.combine(reg1::Quantica.RegionPresets.Region, 
    reg2::Quantica.RegionPresets.Region) = 
 Quantica.RegionPresets.Region{2}(r-> reg1.f(r) && reg2.f(r))

"Bravais vectors for the periodic hamiltonian"
function nanoribbonbravais(p)
    @unpack edges = p
    if edges == :armchair_zigzag
        br1 = (1, 1)
        br2 = (1, -1)
        br3 = (1, 0)
    elseif edges == :zigzag_armchair
        br1 = (1, 1)
        br2 = (1, -1)
        br3 = (0, 1)
    else 
        nothing
    end
    return br1, br2, br3
end

"x and y boundaries for all types of edges :armchair_zigzag, :zigzag_armchair used in hamiltonians.jl"
function xyregs(p)
    @unpack a0, Ln, Ls, W, scale, edges = p
    ϵ = a0/1000
    a0 = a0 * scale
    Ln, Ls, W = edgecondition(Ln, Ls, W, a0, edges)
    if edges == :armchair_zigzag
        xreg = Quantica.RegionPresets.Region{2}(r -> -Ln/2 -Ls < r[1] < Ln/2 + Ls)
        yreg =  Quantica.RegionPresets.Region{2}(r -> ifelse(sign(r[3]) < 0,
            - ϵ <= r[2] <= W + a0/(2*√3)  + ϵ,
            -a0/(2*√3) - ϵ <= r[2] <= W + a0/(2*√3)  + ϵ))
    elseif edges == :zigzag_armchair 
        xreg = Quantica.RegionPresets.Region{2}(r -> ifelse(sign(r[3]) < 0,
            -Ln/2 -Ls -a0/(2*√3) - ϵ <= r[1] <= Ln/2 + a0/(2*√3) + Ls + ϵ,
            -Ln/2 - Ls -a0/(2*√3) - ϵ <= r[1] <= Ln/2 + a0/(2*√3) + Ls + ϵ))
        yreg = Quantica.RegionPresets.Region{2}(r ->  0 < r[2] < W -a0/2 + ϵ) 
    end
    return xreg, yreg
end

"regions used in model.jl"
function modelregs(p)
    @unpack a0, dinter, dintra, Ls, Ln, W, scale, edges = p

    a0 = a0 * scale
    Ln, Ls, W = edgecondition(Ln, Ls, W, a0, edges)
    W = ifelse(edges == :armchair_zigzag, W + a0/sqrt(3), W-a0/2)
    ϵ = a0*1e-3
    edgereg(r) = ifelse(edges == :zigzag_armchair, 
    -Ln/2 -a0/(2*√3) - ϵ <= r[1] <=  -Ln/2 - ϵ  || Ln/2 <= r[1] <=  Ln/2 + a0/(2*√3) + ϵ,
    -Ln/2 < r[1] < -Ln/2 + 2a0/(2*√3) || Ln/2 - 2a0/(2*√3) < r[1] < Ln/2)

    Redgereg(r) = ifelse(edges == :zigzag_armchair, 
    Ln/2 <= r[1] <=  Ln/2 + a0/(2*√3) + ϵ,
    Ln/2 - 2a0/(2*√3) < r[1] < Ln/2)

    normalreg(r) = ifelse(edges == :zigzag_armchair, 
    -Ln/2 -a0/(2*√3) - ϵ <= r[1] <= Ln/2 + a0/(2*√3) + ϵ,
    -Ln/2 < r[1] < Ln/2)
    #- Ln/2 -a0/(2*sqrt(3)) <= r[1] <= Ln/2 +a0/(2*sqrt(3))
    scRreg(r) = r[1] > ifelse(edges == :zigzag_armchair, Ln/2 + a0/(2*√3) + ϵ, Ln/2)
    scLreg(r) = r[1] < ifelse(edges == :zigzag_armchair,-Ln/2 -a0/(2*√3) - ϵ, -Ln/2)
    scRreg2(r) = r[1] > ifelse(edges == :zigzag_armchair, Ln/2 + a0/(2*√3) + ϵ, Ln/2 - a0/(2*√3))
    scLreg2(r) = r[1] < ifelse(edges == :zigzag_armchair,-Ln/2 -a0/(2*√3) - ϵ, -Ln/2 + a0/(2*√3))
    screg(r) = scRreg(r) || scLreg(r)
    screg2(r) = scRreg2(r) || scLreg2(r)
    return a0, Ln, Ls, W, (edgereg, Redgereg, normalreg, scRreg, scLreg, screg, screg2)
end

"condition to prevent dangling bonds in both the zigzag and armchair direction"
edgecondition(Ln, Ls, W, a0, edges) = ifelse(edges == :zigzag_armchair,
    (zigzagcond(Ln, a0), zigzagcond(Ls, a0), armchaircond(W, a0)), 
    (armchaircond(Ln, a0), armchaircond(Ls, a0), zigzagcond(W, a0)))
"Condition to prevent dangling bonds in the zigzag direction"
zigzagcond(length, a0) = 6a0/sqrt(3) * round((length)/(6a0/sqrt(3)))
"Condition to prevent dangling bonds in the armchair direction"
armchaircond(length, a0) = 4a0 * round(length/(4a0))