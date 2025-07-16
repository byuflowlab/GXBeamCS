
# ------------  CLT for Laminates -----------------

"""
    Qbar(lamina)

Computes the lamina stiffness matrix at some arbitrary orientation theta
Transforms a specially orthotropic lamina from principal axis to
an arbitrary axis defined by the ply orientation.
"""
function Qbar(lamina)

    mat = lamina.material
    E1 = mat.E1
    E2 = mat.E2
    nu12 = mat.nu12
    nu21 = nu12*E2/E1
    G12 = mat.G12
    delta = 1.0/(1 - nu12*nu21)

    s, c = sincos(lamina.theta)
    # s = sin(lamina.theta)
    c2 = c*c
    s2 = s*s
    cs = c*s

    Q = Symmetric([
        E1*delta  nu12*E2*delta  0.0;
        nu12*E2*delta  E2*delta  0.0;
        0.0  0.0  G12])

    # Q = @SMatrix [E1*delta  nu12*E2*delta  0.0;
    #               nu12*E2*delta  E2*delta  0.0;
    #               0.0  0.0  G12] #todo: These three allocations are a significant portion of allocations. 
    
    # Q = Symmetric(Q) #Note: Speedup, but creates a typing problem later on for inv(A)

    # @show typeof(Q)

    # Ts = [c2 s2 -2*cs;
    #       s2 c2 2*cs;
    #       cs -cs c2-s2]

    # Te = [c2 s2 -cs;
    #       s2 c2 cs;
    #       2*cs -2*cs c2-s2]

    # Q = @SMatrix [E1*delta  nu12*E2*delta  0.0;
    #               nu12*E2*delta  E2*delta  0.0;
    #               0.0  0.0  G12] #Note: For some reason I wasn't able to use a SMatrix here... 

    Ts = @SMatrix [c2 s2 -2*cs;
                   s2 c2 2*cs;
                   cs -cs c2-s2]

    Te = @SMatrix [c2 s2 -cs;
                   s2 c2 cs;
                   2*cs -2*cs c2-s2]

    Qbar = Symmetric(Ts * Q * Ts')

    # @show typeof(Qbar)

    return Qbar, Ts, Te
end

"""
compute the z locations and overall height for each lamina in the laminate.
"""
function zspacing(laminate)

    n = length(laminate)
    TF = eltype(laminate[1])

    # compute z vector
    z = zeros(TF, n+1)
    for i = 1:n
        z[i+1] = z[i] + laminate[i].t
    end
    h = z[end] - z[1]
    mid = h/2.0  # compute midpoint
    z .-= mid  # recenter at midpoint

    return z, h
end

"""
Same z locations, but doubles them up on the interior so that strain can be computed at either side of the ply.
"""
function zspacingdouble(laminate)
    z, _ = zspacing(laminate)

    # setup new z vector at top and bottom of each ply
    nz = 2*(length(z)-1)
    TF = eltype(laminate[1])
    zvec = zeros(TF, nz)
    zvec[1] = z[1]
    zvec[end] = z[end]
    j = 2
    for i = 2:length(z)-1
        zvec[j] = z[i]
        zvec[j+1] = z[i]
        j += 2
    end

    return zvec
end

"""
Compute the A, B, D stiffness matrices for a thin laminate
"""
function laminatestiffnessmatrix(laminate, z)
    A = Symmetric(zeros(3, 3))
    B = Symmetric(zeros(3, 3))
    D = Symmetric(zeros(3, 3))

    n = length(laminate)
    for k = 1:n
        Q, _, _ = Qbar(laminate[k])
        A += Q*(z[k+1] - z[k])
        B += Q*(z[k+1]^2 - z[k]^2)/2
        D += Q*(z[k+1]^3 - z[k]^3)/3
    end

    return A, B, D
end

"""
    inv3x3(A::AbstractMatrix{T}) where T

Analytically computes the inverse of a 3x3 matrix A.

Returns the inverse matrix if A is invertible.
"""
function inv3x3(A::AbstractMatrix{T}) where T #Todo: WAY WAY slower than inv(A)
    @assert size(A) == (3, 3) "Matrix must be 3x3"

    a, b, c = A[1,1], A[1,2], A[1,3]
    d, e, f = A[2,1], A[2,2], A[2,3]
    g, h, i = A[3,1], A[3,2], A[3,3]

    detA = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
    @assert detA != 0 "Matrix is singular"

    # invA = similar(A) #Todo: double check that this math is correct. -> Error is 2e-13... which is higher than expected. -> The math matches a source online. 
    # invA[1,1] =  (e*i - f*h) / detA
    # invA[1,2] = -(b*i - c*h) / detA
    # invA[1,3] =  (b*f - c*e) / detA
    # invA[2,1] = -(d*i - f*g) / detA
    # invA[2,2] =  (a*i - c*g) / detA
    # invA[2,3] = -(a*f - c*d) / detA
    # invA[3,1] =  (d*h - e*g) / detA
    # invA[3,2] = -(a*h - b*g) / detA
    # invA[3,3] =  (a*e - b*d) / detA

    invA11 =  (e*i - f*h) / detA
    invA12 = -(b*i - c*h) / detA
    invA13 =  (b*f - c*e) / detA
    invA21 = -(d*i - f*g) / detA
    invA22 =  (a*i - c*g) / detA
    invA23 = -(a*f - c*d) / detA
    invA31 =  (d*h - e*g) / detA
    invA32 = -(a*h - b*g) / detA
    invA33 =  (a*e - b*d) / detA

    invA = @SMatrix [invA11 invA12 invA13;
                   invA21 invA22 invA23;
                   invA31 invA32 invA33]
    
    return invA
end

"""
Compute the alpha, beta, delta compliance matrices for a thin laminate given the stiffness matrices
"""
function laminatecompliancematrix(A, B, D)
    # @show size(A), size(B), size(D)
    Ainv = inv(A)
    Hinv = Symmetric(inv(D - B*Ainv*B))
    alpha = Symmetric(Ainv + Ainv*B*Hinv*B*Ainv)
    beta = -Ainv*B*Hinv
    delta = Hinv

    return alpha, beta, delta
end

"""
convenience method to go from laminate definition to its compliance matrices
"""
function laminatecompliance(laminate)
    z, h = zspacing(laminate)
    A, B, D = laminatestiffnessmatrix(laminate, z)
    alpha, beta, delta = laminatecompliancematrix(A, B, D)
    return alpha, beta, delta
end

"""
strains for a thin laminate with forces: N1, N2, N12, M1, M2, M12
"""
function laminatestrains(alpha, beta, delta, laminate, forces)
    C = [alpha beta; beta' delta]

    alleps = C*forces
    epsilonbar = alleps[1:3]
    kappa = alleps[4:6]

    # setup new z vector at top and bottom of each ply
    zvec = zspacingdouble(laminate)
    nz = length(zvec)

    # @show typeof(alpha), typeof(beta), typeof(delta)

    TF = promote_type(typeof(forces[1]), typeof(zvec[1]), typeof(alpha[1]), typeof(beta[1]), typeof(delta[1]))

    epsilonp = zeros(TF, 3, nz)
    for i = 1:3
        epsilonp[i, :] = epsilonbar[i] .+ kappa[i]*zvec
    end

    return epsilonbar, kappa, zvec, epsilonp
end

"""
stress for a thin laminate given stresses.
"""
function laminatestresses(laminate, epsilonp)

    TF = promote_type(typeof(epsilonp[1]), typeof(laminate[1].t)) #todo: Somehow this requires a bunch of allocations. 

    n = length(laminate)
    sigmap = zeros(TF, 3, 2*n)
    sigma = zeros(TF, 3, 2*n)
    epsilon = zeros(TF, 3, 2*n)

    i = 1
    for k = 1:n
        Q, Ts, Te = Qbar(laminate[k])
        # compute for top and bottom
        sigmap[:, i] = Q*epsilonp[:, i]
        sigmap[:, i+1] = Q*epsilonp[:, i+1]
        sigma[:, i] = Te' * sigmap[:, i]
        sigma[:, i+1] = Te' * sigmap[:, i+1]
        epsilon[:, i] = Ts' * epsilonp[:, i]
        epsilon[:, i+1] = Ts' * epsilonp[:, i+1]
        i += 2
    end

    return sigmap, sigma, epsilon
end



function clt(laminate, forces)

    z, h = zspacing(laminate)
    A, B, D = laminatestiffnessmatrix(laminate, z)
    alpha, beta, delta = laminatecompliancematrix(A, B, D)
    epsilonbar, kappa, zvec, epsilonp = laminatestrains(alpha, beta, delta, laminate, forces)
    sigmap, sigma, epsilon = laminatestresses(laminate, epsilonp)
    failure = tsai_wu_plane(sigma, laminate)

    return sigma, epsilon, failure
end



function tsai_hill(sigma, laminate)
    
    n = length(laminate)
    # failure = zeros(2*n)  # failure if > 1
    T = eltype(sigma)
    failure = Vector{T}(undef, 2*n)
    i = 1
    for k = 1:n
        (; S1t, S1c, S2t, S2c, S12) = laminate[k].material
        for j = i:i+1
            S1 = sigma[1, j] >= 0.0 ? S1t : S1c
            S2 = sigma[2, j] >= 0.0 ? S2t : S2c
            
            failure[j] = sigma[1, j]^2/S1^2 + sigma[2, j]^2/S2^2 + sigma[3, j]^2/S12^2 - sigma[1, j]*sigma[2, j]/S1^2
        end
        i += 2
    end

    return failure
end

function tsai_wu_plane(sigma, laminate)
    
    n = length(laminate)
    # failure = zeros(2*n)  # failure if > 1
    T = eltype(sigma)
    failure = Vector{T}(undef, 2*n)
    i = 1
    for k = 1:n
        (; S1t, S1c, S2t, S2c, S12) = laminate[k].material
        for j = i:i+1
            # failure[j] = sigma[1, j]^2/(S1t*S1c) + sigma[2, j]^2/(S2t*S2c) - sigma[1, j]*sigma[2, j]/sqrt(S1t*S1c*S2t*S2c) +
            #     sigma[1, j]*(1/S1t - 1/S1c) + sigma[2, j]*(1/S2t - 1/S2c) + sigma[3, j]^2/S12^2 #Dr. Ning's math. 
            failure[j] = sigma[1, j]^2/(S1t*S1c) + sigma[2, j]^2/(S2t*S2c) - 2*sigma[1, j]*sigma[2, j]/sqrt(S1t*S1c*S2t*S2c) +
                sigma[1, j]*(1/S1t - 1/S1c) + sigma[2, j]*(1/S2t - 1/S2c) + sigma[3, j]^2/S12^2 #My math (Sympy)
        end
        i += 2
    end

    return failure
end

function full_tsai_wu(sigma, laminate)
    # S1t = strength.S1t
    # S1c = strength.S1c
    # S2t = strength.S2t
    # S2c = strength.S2c
    # S12 = strength.S12

    # _, n = size(sigma)
    n = length(laminate)
    # failure = zeros(2*n)  # failure if > 1
    T = eltype(sigma)
    failure = Vector{T}(undef, 2*n)
    i = 1
    for k = 1:n
        (; S1t, S1c, S2t, S2c, S12) = laminate[k].material
        S3t = S2t
        S3c = S2c
        S13 = S12
        S23 = S12
        for j = i:i+1
            failure[j] = sigma[1, j]^2/(S1t*S1c) +
                    sigma[2, j]^2/(S2t*S2c) +
                    sigma[3, j]^2/(S3t*S3c) +
                    sigma[4, j]^2/S12^2 +
                    sigma[5, j]^2/S13^2 +
                    sigma[6, j]^2/S23^2 +
                    sigma[1, j]*(1/S1t - 1/S1c) +
                    sigma[2, j]*(1/S2t - 1/S2c) +
                    sigma[3, j]*(1/S3t - 1/S3c) -
                    sigma[1, j]*sigma[2, j]/sqrt(S1t*S1c*S2t*S2c) -
                    sigma[1, j]*sigma[3, j]/sqrt(S1t*S1c*S3t*S3c) -
                    sigma[2, j]*sigma[3, j]/sqrt(S2t*S2c*S3t*S3c)
        end
        i += 2
    end

    return failure
end

function max_stress(sigma, laminate)
    
    n = length(laminate)
    T = eltype(sigma)
    failure = Array{T, 2}(undef, 3, 2*n)  # failure if > 1
    i = 1
    for k = 1:n
        (; S1t, S1c, S2t, S2c, S12) = laminate[k].material
        for j = i:i+1
            S1 = sigma[1, j] >= 0.0 ? S1t : -S1c
            S2 = sigma[2, j] >= 0.0 ? S2t : -S2c
            
            failure[1, j] = sigma[1, j]/S1
            failure[2, j] = sigma[2, j]/S2
            failure[3, j] = sqrt((sigma[3, j]/S12)^2)
        end
        i += 2
    end

    return failure
end

function mybuckling_strain(epsilon, laminate, b, E_axial)
    
    n = length(laminate)
    T = eltype(epsilon)
    failure = Vector{T}(undef, 2*n)  # failure if > 1
    i = 1
    for k = 1:n # iterate over the layers
        z, _ = zspacing(laminate)
        _, _, D = laminatestiffnessmatrix(laminate, z)
        #Todo: I'm not sure that this should be the laminate stiffness matrix or something else. 
        Ncrit = 3.6*(pi/b)^2*D[1, 1]
        t = laminate[k].t 
        for j = i:i+1 # iterate over the top and bottom of the layer
            
            epsilon_crit = -Ncrit/(E_axial*t)
            
            failure[j] = epsilon[1, j]/epsilon_crit
        end
        i += 2
    end


    return failure
end

function buckling_strain(epsilon, laminate, b)
    

    ### Based on RotorSE
    #This is what they use for Eaxial: https://github.com/WISDEM/RotorSE/blob/d497e96a612b56aa37b57d35c64aa5049bb9136f/src/rotorse/precomp.py

    # A, B, D, totalHeight = self.compositeMatrices(sector)

    # S = np.vstack((np.hstack((A, B)), np.hstack((B, D))))

    # # E_eff_x = N_x/h/eps_xx and eps_xx = S^{-1}(0,0)*N_x (approximately)
    # detS = np.linalg.det(S)
    # Eaxial = detS/np.linalg.det(S[1:, 1:])/totalHeight


    #This is their main approach
    #https://github.com/WISDEM/RotorSE/blob/d497e96a612b56aa37b57d35c64aa5049bb9136f/src/rotorse/rotor.py
    z, _ = zspacing(laminate)
    h = z[end] - z[1]
    # @show h, sum(l.t for l in laminate) #True
    A, B, D = laminatestiffnessmatrix(laminate, z)

    S = [A B; B D]
    E_axial = det(S)/det(S[2:end, 2:end])/h

    dterm = sqrt(D[1, 1]*D[2, 2]) + D[1, 2] + 2*D[3,3]
    Ncrit = 3.6*((pi/b)^2)*dterm

    eps_crit = -Ncrit/(E_axial*h)

    failure = epsilon[1,1]/eps_crit


    return failure
end

# ------------  CLT for Beams of Laminates -----------------

struct BeamSection{VL, VF} #Note: I suggest that we rename this to Region, or SectionRegion, or something like that. Beam section makes it sound like it is a section of a beam (a cross section), not a section of a cross-section. -> Segment might also be a good name. -> segments are passed in... so maybe not. 
    laminate::VL  #Vector{Lamina}
    y::VF  # Vector{Float}
    z::VF  # Vector{Float}
end

function get_section_floattype(section)
    return promote_type(typeof(section.y[1]), typeof(section.z[1]))
end


"""
    get_beam_sections(x, y, chord, twist, paxis, xbreak, weblocs, segments, web_segments)

Creates a vector of BeamSections from the information for the airfoil mesh. 

**Arguments**
- x::Vector{Float64}: normalized x-coordinates of the airfoil
- y::Vector{Float64}: normalized y-coordinates of the airfoil
- chord::Float64: chord length of the airfoil (meters)
- twist::Float64: twist angle of the airfoil (radians)
- paxis::Float64: position of the reference axis as a fraction of the chord. Twisting about this axis. 
- xbreak::Vector{Float64}: normalized x-coordinates of the different regions of the airfoil
- weblocs::Vector{Float64}: normalized x-coordinates of the webs
- segments::Vector{Lamina}: laminas for each region
- web_segments::Vector{Lamina}: laminas for each web

**Returns**

"""
function get_beam_sections(x, y, chord, twist, paxis, xbreak, weblocs, segments, web_segments; fit=Akima, le_idx=argmin(x))
    ### Check that inputs are good. 
    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end
    #Todo: Check that xbreak starts and ends with 0 and 1, respectively.

    if !isapprox(minimum(x), 0) || !isapprox(maximum(x), 1)
        throw(ArgumentError("x must start at 0 and end at 1"))
    end
    

    # le_idx = argmin(x) #Todo. This might not capture the LE of the airfoil. -> Sticking it as an optional argument so the user can specify it. 
    xtop = reverse(x[1:le_idx])
    ytop = reverse(y[1:le_idx]) 
    xbot = x[le_idx:end] #todo: I probably need to know if the LE is repeated. 
    ybot = y[le_idx:end]
    topfit = fit(xtop, ytop)
    botfit = fit(xbot, ybot)

    nw = length(web_segments) #todo: Maybe add multiple points for web? -> Don't know if it changes anything. 
    nb = length(xbreak)
    ns = 2*(nb - 1) + nw #Number of regions

    layertype = typeof(segments[1])
    floattype = promote_type(typeof(x[1]), typeof(chord), typeof(twist), typeof(paxis))


    vectortype = Vector{floattype}
    sections = Vector{BeamSection{layertype, vectortype}}(undef, ns) 

    
    # @show vectortype
    # println("")

    # sections = Vector{BeamSection{layertype, vectortype}}(undef, ns) 

    s, c = sincos(twist) #I think applied twist correctly. 
    xc = paxis * chord
    idx_outer_top = 1
    idx_outer_bot = 1
    # @show xc, s, c
    
    for i = 1:(nb-1) #Iterate over the regions
        ### Top
        idx_inner = findfirst(x -> x >= xbreak[i+1], xtop)
        idx_inner = xtop[idx_inner] != xbreak[i+1] ? idx_inner - 1 : idx_inner

        #Scale the coordinates
        xi = xtop[idx_outer_top:idx_inner]*chord 
        yi = ytop[idx_outer_top:idx_inner]*chord

        #Rotate about the pitch axis.
        x_i = @. (xi - xc)*c + yi*s + xc
        y_i = @. -(xi - xc)*s + yi*c
        
        sections[i] = BeamSection(segments[i], reverse(x_i), reverse(y_i))

        #Update the outer index
        idx_outer_top = idx_inner



        ### Bottom
        idx_inner = findfirst(x -> x >= xbreak[i+1], xbot)
        idx_inner = xbot[idx_inner] != xbreak[i+1] ? idx_inner - 1 : idx_inner

        xi = xbot[idx_outer_bot:idx_inner]*chord
        yi = ybot[idx_outer_bot:idx_inner]*chord

        x_i = @. (xi - xc)*c + yi*s + xc
        y_i = @. -(xi - xc)*s + yi*c
        
        sections[i+nb-1] = BeamSection(segments[i], x_i, y_i)

        idx_outer_bot = idx_inner
    end


    ### Add in the webs
    idx = 2*(nb - 1)
    # @show xc, s, c
    for i in eachindex(weblocs)
        #Get the upper and lower coordinates of the web.
        ytop = topfit(weblocs[i])
        ybot = botfit(weblocs[i])
        xi = [weblocs[i], weblocs[i]].*chord
        yi = [ybot, ytop].*chord

        #Rotate about the pitch axis. 
        x_i = @. (xi - xc)*c + yi*s + xc 
        y_i = @. -(xi - xc)*s + yi*c

        # @show weblocs[i], ybot, ytop, x_i, y_i
        sections[idx + i] = BeamSection(web_segments[i], x_i, y_i)
        # println("")
    end

    return sections
end

function make_top_section(xtop, ytop, chord, s, c, xc, xbreak, segments, i)
    
    ### Top
    idx_inner = findfirst(x -> x >= xbreak[i+1], xtop)
    idx_inner = xtop[idx_inner] != xbreak[i+1] ? idx_inner - 1 : idx_inner

    if i == 1
        idx_outer_top = 1
    else
        idx_outer_top = findfirst(x -> x >= xbreak[i], xtop)
        idx_outer_top = xtop[idx_outer_top] != xbreak[i] ? idx_outer_top - 1 : idx_outer_top
    end

    #Scale the coordinates
    xi = xtop[idx_outer_top:idx_inner]*chord 
    yi = ytop[idx_outer_top:idx_inner]*chord

    #Rotate about the pitch axis.
    x_i = @. (xi - xc)*c + yi*s + xc
    y_i = @. -(xi - xc)*s + yi*c
    
    section = BeamSection(segments[i], reverse(x_i), reverse(y_i))

    return section
end

function make_bot_section(xbot, ybot, chord, s, c, xc, xbreak, segments, i)
    ### Bottom
    idx_inner = findfirst(x -> x >= xbreak[i+1], xbot)
    idx_inner = xbot[idx_inner] != xbreak[i+1] ? idx_inner - 1 : idx_inner

    if i == 1
        idx_outer_bot = 1
    else
        idx_outer_bot = findfirst(x -> x >= xbreak[i], xbot)
        idx_outer_bot = xbot[idx_outer_bot] != xbreak[i] ? idx_outer_bot - 1 : idx_outer_bot
    end

    xi = xbot[idx_outer_bot:idx_inner]*chord
    yi = ybot[idx_outer_bot:idx_inner]*chord

    x_i = @. (xi - xc)*c + yi*s + xc
    y_i = @. -(xi - xc)*s + yi*c
    
    section = BeamSection(segments[i], x_i, y_i)

    return section
end

function make_web_section(topfit, botfit, webloc, chord, s, c, xc, web_segment)
    
        #Get the upper and lower coordinates of the web.
    ytop = topfit(webloc)
    ybot = botfit(webloc)
    xi = [webloc, webloc].*chord
    yi = [ybot, ytop].*chord

    #Rotate about the pitch axis. 
    x_i = @. (xi - xc)*c + yi*s + xc 
    y_i = @. -(xi - xc)*s + yi*c

    section = BeamSection(web_segment, x_i, y_i)   
    return section
end

"""
An out of place verion of `get_beam_sections()` for ReverseDiff compatibility.
"""
function get_beam_sections_oop(x, y, chord, twist, paxis, xbreak, weblocs, segments, web_segments; fit=Akima, le_idx=argmin(x))
    ### Check that inputs are good. 
    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end
    #Todo: Check that xbreak starts and ends with 0 and 1, respectively.

    if !isapprox(minimum(x), 0) || !isapprox(maximum(x), 1)
        throw(ArgumentError("x must start at 0 and end at 1"))
    end
    

    # le_idx = argmin(x) #Todo. This might not capture the LE of the airfoil. -> Sticking it as an optional argument so the user can specify it. 
    xtop = reverse(x[1:le_idx])
    ytop = reverse(y[1:le_idx]) 
    xbot = x[le_idx:end] #todo: I probably need to know if the LE is repeated. 
    ybot = y[le_idx:end]
    topfit = fit(xtop, ytop)
    botfit = fit(xbot, ybot)

    nw = length(web_segments) #todo: Maybe add multiple points for web? -> Don't know if it changes anything. 
    nb = length(xbreak)
    ns = 2*(nb - 1) + nw #Number of regions

    layertype = typeof(segments[1])
    floattype = promote_type(typeof(x[1]), typeof(chord), typeof(twist), typeof(paxis))


    vectortype = Vector{floattype}
    sections = Vector{BeamSection{layertype, vectortype}}(undef, ns) 

    
    # @show vectortype
    # println("")

    # sections = Vector{BeamSection{layertype, vectortype}}(undef, ns) 

    s, c = sincos(twist) #I think applied twist correctly. 
    xc = paxis * chord
    

    top_sections = [make_top_section(xtop, ytop, chord, s, c, xc, xbreak, segments, i) for i in 1:(nb-1)]
    bot_sections = [make_bot_section(xbot, ybot, chord, s, c, xc, xbreak, segments, i) for i in 1:(nb-1)]


    # i in eachindex(weblocs)
    web_sections = [make_web_section(topfit, botfit, weblocs[i], chord, s, c, xc, web_segments[i]) for i in eachindex(weblocs)]

    sections = vcat(top_sections, bot_sections, web_sections)

    return sections
end




struct CLTCache{TM1, TM2, TM3, TM4}
    F::TM1
    L::TM2
    Wbar::TM3
    S::TM4
end

struct CLT{TL, TF, TM1, TM2, TM3, TM4} <: CompositeSectionAnalysis 
    sections::Vector{BeamSection{TL, TF}}  # a vector of beam sections
    closed_section::Bool
    cache::CLTCache{TM1, TM2, TM3, TM4}
end

function CLT(sections, closed_section, TF::DataType)  # initialize empty cache

    F = zeros(TF, 2, 2)
    L = zeros(TF, 2, 4)
    Wbar = Symmetric(zeros(TF, 4, 4))
    S = Symmetric(zeros(TF, 4, 4))
    return CLT(sections, closed_section, CLTCache(F, L, Wbar, S))
end

function CLT(sections, closed_section)
    TF = promote_type(typeof(sections[1].laminate[1].material.E1), typeof(sections[1].y[1]))

    return CLT(sections, closed_section, TF)
end

CLT(sections) = CLT(sections, true)  # default to closed section


"""
    count_clt_elements(clt::CLT)
Count the number of elements in a CLT object.

**Arguments**
- clt::CLT: the CLT object

**Returns**
- Int: the number of elements in the CLT object
"""
function count_clt_elements(clt::CLT)
    n = 0
    for i in eachindex(clt.sections)
        sec = clt.sections[i]
        n_seg = length(sec.y)-1
        n += n_seg*length(sec.laminate)
    end #End section loop
    return n
end

"""
    find_section_element(clt::CLT, i::Int)

Find the section of a CLT object that corresponds to the ith element of the layup.

**Arguments**
- clt::CLT: the CLT object
- i::Int: the index of the element

**Returns**
- Int: the index of the section containing the ith element
"""
function find_section_element(clt::CLT, i::Int)
    idx = 1
    for j in eachindex(clt.sections)
        sec = clt.sections[j]
        n_seg = length(sec.y) - 1
        n_elem = n_seg * length(sec.laminate)

        if i <= idx + n_elem - 1
            return j
        end

        idx += n_elem
    end

    throw(ArgumentError("Element index $i is out of range for the given CLT object."))
end

"""
    find_section_layer(clt::CLT, i::Int)

Find the layer of a section corresponding to the ith element of a CLT layup.

**Arguments**
- clt::CLT: the CLT object
- i::Int: the index of the element

**Returns**
- Tuple{Int, Int}: the section index and the layer index within the section
"""
function find_section_layer(clt::CLT, i::Int)
    idx = 1
    for j in eachindex(clt.sections)
        sec = clt.sections[j]
        n_seg = length(sec.y) - 1
        n_elem = n_seg * length(sec.laminate)

        if i <= idx + n_elem - 1
            local_idx = i - idx
            layer_idx = div(local_idx, n_seg) + 1
            return j, layer_idx
        end

        idx += n_elem
    end

    throw(ArgumentError("Element index $i is out of range for the given CLT object."))
end


"""
    get_clt_element_position(clt, i)
Grab the position of the ith layer element in the CLT.

**Arguments**
- clt::CLT: the CLT object
- i::Int: the index of the element

**Returns**
- Tuple{Vector{Float64}, Vector{Float64}}: the y and z coordinates of the element
"""
function get_clt_element_position(clt::CLT, i)

    idx = 1
    for j in eachindex(clt.sections)
        sec = clt.sections[j]
        n_seg = length(sec.y)-1

        for k in 1:n_seg
            lam = sec.laminate
            T = 0.0
            for l in eachindex(lam)
                if idx == i
                    x = sec.y[k+1] - sec.y[k]
                    y = sec.z[k+1] - sec.z[k]
                    L = sqrt(x^2 + y^2) #Element length
                    
                    #Normal vector to the element
                    nx = -y/L
                    ny = x/L

                    xm = (sec.y[k] + sec.y[k+1])/2
                    ym = (sec.z[k] + sec.z[k+1])/2

                    xp = xm + nx*(T + lam[l].t/2)
                    yp = ym + ny*(T + lam[l].t/2)
                    
                    return xp, yp
                end
                idx += 1
                T += lam[l].t
            end #End laminate (element) loop
        end #End segment loop
    end #End section loop
    
    @warn("layer index out of range ($i > $idx)")
    return 0
end

"""
    get_section_indices(clt::CLT, i)
Get the element indices for the ith section from a CLT layup.

**Arguments**
- clt::CLT: the CLT object
- i::Int: the index of the section

**Returns**
- Vector{Int}: the indices of the elements in the ith section
"""
function get_section_indices(clt::CLT, i)
    idx = 1
    n = length(clt.sections)
    for j in 1:(i-1)
        sec = clt.sections[j]
        n_seg = length(sec.y)-1
        idx += n_seg*length(sec.laminate)

    end #End section loop

    seci = clt.sections[i]
    n_segi = length(seci.y)-1
    m_segi = length(seci.laminate)
    idxs = idx+1:idx+n_segi*m_segi

    return idxs
end

"""
    get_section_layer_indices(section::BeamSection, indices, i_layer)
Get the layer indices from a section. 

**Arguments**
- section::BeamSection: the section object
- indices::Vector{Int}: the indices of the elements in the section
- i_layer::Int: the index of the layer

**Returns**
- Vector{Int}: the indices of the elements in the i_layerth layer
"""
function get_section_layer_indices(section::BeamSection, indices, i_layer)

    n_idx = length(indices)
    num_layers = length(section.laminate)
    n_seg = length(section.y)-1

    if n_idx != n_seg*num_layers
        @warn("Number of indices does not match number of elements in section")
    end

    return indices[i_layer:num_layers:end]
end


"""
    get_section_layer_indices(clt::CLT, i, i_layer)
Get the element indices for the ith section and the i_layerth layer from a CLT layup. 

**Arguments**
- clt::CLT: the CLT object
- i::Int: the index of the section
- i_layer::Int: the index of the layer

**Returns**
- Vector{Int}: the indices of the elements in the ith section and i_layerth layer
"""
function get_section_layer_indices(clt::CLT, i, i_layer)
    idxs = get_section_indices(clt, i)
    sec = clt.sections[i]
    return get_section_layer_indices(sec, idxs, i_layer)
end

"""
    get_section_thickness(clt::CLT, i)
Get the thickness of the ith section in a CLT layup.

**Arguments**
- clt::CLT: the CLT object
- i::Int: the index of the section

**Returns**
- Float64: the thickness of the ith section
"""
function get_section_thickness(clt::CLT, i)
    
    sec = clt.sections[i]
    return sum(l.t for l in sec.laminate)
end



function fullstiffnessmatrix(P, S)
    K = zeros(6, 6)
    K[2:3, 2:3] .= S
    idx = [1, 5, 6, 4]
    j = 1
    for i in idx
        K[i, idx] .= P[j, :]
        j += 1
    end

    return K
end

"""
"""
function rotatestiffnessmatrix(K, r, theta)
    R = [1.0 0 0;
        0 cos(theta) -sin(theta);
        0 sin(theta) cos(theta)]
    p = [0.0 -r[3] r[2];
        r[3] 0 -r[1];
        -r[2] r[1] 0]
    HinvT = [R p*R; zeros(3, 3) R]
    Hinv = transpose(HinvT)

    Kp = Hinv*K*HinvT
    return Kp
end

function rotatecompliancematrix(K, r, theta)
    R = [1.0 0 0;
        0 cos(theta) -sin(theta);
        0 sin(theta) cos(theta)]
    p = [0.0 -r[3] r[2];
        r[3] 0 -r[1];
        -r[2] r[1] 0]
    HinvT = [R p*R; zeros(3, 3) R]
    # Hinv = transpose(HinvT)

    HT = inv(HinvT)
    H = transpose(HT)

    Kp = HT*K*H
    return Kp
end


"""
    rotate_compliance(S, theta)

Rotate the compliance matrix S by theta radians.

**Arguments**
- S::Symmetric{Float64, Matrix{Float64}}: the compliance matrix
- theta::Float64: the angle to rotate by
"""
function rotate_compliance(S, theta)
    #Todo: I don't know if the elements are in the correct order. 
    s, c = sincos(theta)
    s2, c2 = s^2, c^2
    sc = s*c

    Tsigma = [c2   s2   0  -2sc  0  0
              s2   c2   0   2sc  0  0
               0    0   1    0  0  0
              sc  -sc   0 c2-s2 0  0
               0    0   0    0  c -s
               0    0   0    0  s  c] #TODO: I could mathematically invert this and code it in. 
    
    Tepsilon = [c2   s2   0  -sc  0  0
                s2   c2   0   sc  0  0
                 0    0   1    0  0  0
                2sc  -2sc   0 c2-s2 0  0
                 0    0   0    0  c -s
                 0    0   0    0  s  c]

    Srot = Tepsilon*S*inv(Tsigma)

    return Srot
end



function compliance_matrix(clt::CLT, shear_center=true)

    m = length(clt.sections) #number of sections

    Pbar = zeros(4, 4)
    Im = zeros(2, 4)
    F = zeros(2, 2)
    A = 0.0

    for i = 1:m #Iterate over the number of sections
        sec = clt.sections[i]
        yp = sec.y #The coordinates of the section
        zp = sec.z
        n = length(yp) #The number of points in the section

        alpha, beta, delta = laminatecompliance(sec.laminate)
        a = Symmetric([alpha[1, 1] beta[1, 1] beta[1, 3];
                 beta[1, 1] delta[1, 1] delta[1, 3];
                 beta[1, 3] delta[1, 3] delta[3, 3]])
        atinv = (a[2, 2]*a[3, 3] - a[2, 3]^2) / det(a)
        # Atilde = inv(atilde)  # TODO: we only need (1, 1) component so don't need to invert everything

        for k = 2:n #Iterate over the "elements" in the section
            ybar = (yp[k-1] + yp[k])/2 #The element midpoint
            zbar = (zp[k-1] + zp[k])/2
            b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2) #The length of the element
            ca = (yp[k] - yp[k-1])/b #The cosine of the angle of the element #Todo: Will this correctly orient elements? Check the math.
            sa = (zp[k] - zp[k-1])/b #The sine of the angle of the element
            #Add the element area to the cross section area #TODO: Wait... A is reset at the beginning of the function, so this is like the cumulative area of all the sections? 
            A += 0.5*(zp[k-1] + zp[k]) * (yp[k-1] - yp[k])  # if closed section (trapezoid formula for polygon area: https://en.wikipedia.org/wiki/Shoelace_formula)

            Rk = [1.0 zbar ybar 0.0;
                0.0 ca -sa 0.0;
                0.0 sa ca 0.0;
                0 0 0 1] #Todo: What is this matrix? -> It looks like a rotation matrix, but it is a 4x4 and also has the midpoints on it. 
            omega = 1.0/b*Symmetric([
                    alpha[1, 1] beta[1, 1] 0.0 -beta[1, 3]/2.0;
                    beta[1, 1] delta[1, 1] 0.0 -delta[1, 3]/2.0;
                    0 0 12.0/(atinv*b^2) 0; #Note: b changes every k iteration, so this must be created every iteration. 
                    -beta[1, 3]/2.0 -delta[1, 3]/2.0 0 delta[3, 3]/4.0])
            # omegainv = inv(omega)
            Pbar += Rk'*(omega\Rk)

            I1 = [alpha[1, 3] beta[3, 1] 0.0 -beta[3, 3]/2.0;
                  beta[1, 2] delta[1, 2] 0.0 -delta[2, 3]/2.0]
            Im += I1*(omega\Rk)  # repeated, could cache

            F1 = [alpha[3, 3] beta[3, 2];
                 beta[3, 2] delta[2, 2]]
            F += b*F1 - I1*(omega\I1')
        end

        # if i <= 10   #TODO: temporary hack
        #     A += abs(Asub)/2.0
        # end
        # A += Asub/2.0
    end
    L = -Im
    L[1, 4] += 2*A
    if clt.closed_section
        Pbar += L'*(F\L)
    end
    Wbar = inv(Symmetric(Pbar))
    cent = -Symmetric([Wbar[2, 2] Wbar[2, 3]; Wbar[2, 3] Wbar[3, 3]]) \ [Wbar[1, 2]; Wbar[1, 3]]
    zc = cent[1]; yc = cent[2]

    Rb = [1.0 0 0 0;
          zc 1 0 0;
          yc 0 1 0;
          0 0 0 1]
    S = Symmetric(Rb'*Wbar*Rb)  # compliance
    # K = inv(S)
    sc = [0.0, 0.0]  #Shear Center #Todo: 
    tc = [yc, zc] #Tension Center

    TF = eltype(S)
    Sfull = zeros(TF, 6, 6)
    idx = [1, 5, 6, 4]
    for i = 1:4
        for j = i:4
            Sfull[idx[i], idx[j]] = S[i, j]
        end
    end
 
    clt.cache.S .= S #Todo. Does this S need the effect from the shear flow? -> It looks like shear flow doesn't replace any of the indices in S, so it should be fine.


    if clt.closed_section
        s, S, ysc, zsc = shearflow_general_attempt(clt.sections, Wbar, F, L, yc, zc) 
        # s, S, ysc, zsc = shearflow(clt.sections, S, yc, zc; closedsection=clt.closed_section) #Super mega slow. Todo: I'm not sure on the S input. 
        Sfull[2, 2] = s[1, 1]
        Sfull[3, 3] = s[2, 2]

        # sc = [ysc, zsc] #Note: I don't think ysc and zsc are actual the shear center coordinates. 
        # Sfull[2, 2] = s[2, 2]
        # Sfull[3, 3] = s[1, 1]
    end


    Sfull = Symmetric(Sfull)

    # move to sc
    if shear_center
        ysc = sc[1]; zsc = sc[2]
        P = [0 zsc -ysc; -zsc 0 0; ysc 0 0]
        Hinv = [I P; zeros(3, 3) I]
        HinvT = [I zeros(3, 3); transpose(P) I]
        Sfull = Hinv * Sfull * HinvT
    end

    # save entries in cache for strain evaluation
    clt.cache.F .= F
    clt.cache.L .= L
    clt.cache.Wbar .= Wbar
    # clt.cache.S .= S

    # s, S, ysc, zsc = shearflow(sections, Wbar, F, L, yc, zc)
    return Sfull, sc, tc
end


function mass_matrix_clt(clt::CLT; reference=nothing)
    return mass_matrix(clt.sections; reference=reference)
end


"""
    mass_matrix(sections::Vector{BeamSection}; reference=nothing)

Finding the mass matrix of the composite section.

**Arguments**
- sections::Vector{BeamSection}: a vector of beam sections, each with a `laminate` field containing the laminates in the section.
- reference::Vector{Float64}: the reference frame for the mass matrix: [x, y, theta]. x, y are the origin of the reference frame (often some percentage of the chord), and theta is the twist angle.

**Returns**
- M: the mass matrix
- [xm, ym]: the center of mass
"""
function mass_matrix(sections; reference=nothing)

    mass = 0.0 #TODO: Typing? -> Dr. Ning uses this. 
    #Center of mass
    xm = 0.0 
    ym = 0.0

    for i in eachindex(sections)
        sec = sections[i]

        ns = length(sec.y)

        for j in 1:ns-1 #Iterate over the number of elements in the section
            x = sec.y[j+1] - sec.y[j] #The x distance between the element end points
            y = sec.z[j+1] - sec.z[j] #The y distance between the element end points
            L = sqrt(x^2 + y^2)

            #Midpoint of the element
            xmid = (sec.y[j+1] + sec.y[j])/2
            ymid = (sec.z[j+1] + sec.z[j])/2

            #Normal vector to the element
            nx = -y/L
            ny = x/L

            # @show nx, ny

            T = 0.0 #The total thickness travelled so far. 

            for k in eachindex(sec.laminate) #Iterate over the laminates in the section
                t = sec.laminate[k].t #The thickness of the laminate
                rho = sec.laminate[k].material.rho #The density of the laminate
                A = L*t #Area of the section

                #Centroid of the element lamina
                xc = xmid + nx*(T + t/2)
                yc = ymid + ny*(T + t/2)

                # @show xc, yc

                dmass = A*rho
                mass += dmass
                xm  += dmass*xc
                ym  += dmass*yc

                T += t #Increment the thickness
            end #End looping over laminates
        end #end looping over segment elements
    end #End looping over sections

    xm /= mass
    ym /= mass

    Ixx = 0.0
    Iyy = 0.0
    Ixy = 0.0

    for i in eachindex(sections)
        sec = sections[i]

        ns = length(sec.y)

        for j in 1:ns-1 #Iterate over the number of elements in the section
            x = sec.y[j+1] - sec.y[j] #The x distance between the element end points
            y = sec.z[j+1] - sec.z[j] #The y distance between the element end points
            L = sqrt(x^2 + y^2)

            #Midpoint of the element
            xmid = (sec.y[j+1] + sec.y[j])/2
            ymid = (sec.z[j+1] + sec.z[j])/2

            #Normal vector to the element
            nx = -y/L
            ny = x/L

            #Angles
            ctheta = x/L
            stheta = y/L

            T = 0.0 #The total thickness travelled so far.

            for k in eachindex(sec.laminate) #Iterate over the laminates in the section
                t = sec.laminate[k].t #The thickness of the laminate
                rho = sec.laminate[k].material.rho #The density of the laminate
                A = L*t #Area of the section

                #Centroid of the element lamina
                xc = xmid + nx*(T + t/2)
                yc = ymid + ny*(T + t/2)

                # Ixx += rho*A*(yc - ym)^2
                # Iyy += rho*A*(xc - xm)^2
                # Ixy += rho*A*(xc - xm)*(yc - ym)

                ## Mass moment of inertia about the centroid
                Ixx_c = rho*L*(t^3)/12
                Iyy_c = rho*t*(L^3)/12
                # Ixy_c = 0.0

                ## Rotate the inertia
                Ixx_cp = Ixx_c*(ctheta^2) + Iyy_c*(stheta^2) #-2*Ixy_c*stheta*ctheta
                Iyy_cp = Ixx_c*(stheta^2) + Iyy_c*(ctheta^2) #+2*Ixy_c*stheta*ctheta
                Ixy_cp = (Ixx_c - Iyy_c)*stheta*ctheta #+ Ixy_c*(ctheta^2 - stheta^2)

                ## Shift the inertia (parallel axis theorem)
                # Ixx += Ixx_cp + rho*A*((yc - ym)^2)
                # Iyy += Iyy_cp + rho*A*((xc - xm)^2)
                # Ixy += Ixy_cp + rho*A*((xc - xm)*(yc - ym))
                Ixx_ = Ixx_cp + rho*A*((yc - ym)^2)
                Iyy_ = Iyy_cp + rho*A*((xc - xm)^2)
                Ixy_ = Ixy_cp + rho*A*((xc - xm)*(yc - ym))

                Ixx += Ixx_
                Iyy += Iyy_
                Ixy += Ixy_

                # @show t, L, A, rho, xc, yc, xm, ym, Ixx_, Iyy_, Ixy_
                T += t #Increment the thickness
            end #End looping over laminates
        end #end looping over segment elements
    end #End looping over sections

    if !isnothing(reference)
        #Shift to the reference point
        x = [xm-reference[1], ym-reference[2]] 
        #Rotate about the reference point into reference frame
        R = [cos(reference[3]) -sin(reference[3]); sin(reference[3]) cos(reference[3])]
        x = R*x
        #Rename 
        xm = x[1]
        ym = x[2]
    end

    M = Symmetric([
        mass 0.0 0 0 mass*ym -mass*xm
        0 mass 0 -mass*ym 0 0
        0 0 mass mass*xm 0 0
        0 -mass*ym mass*xm Ixx+Iyy 0 0
        mass*ym 0 0 0 Ixx -Ixy
        -mass*xm 0 0 0 -Ixy Iyy
    ])

    return M, [xm, ym]
end

"""
    shift_mass_matrix(M, cm, sc)

Shift the mass matrix from the center of mass to the shear center (or any other point). 

**Arguments**
- M: the mass matrix
- cm: the center of mass
- sc: the shear center

**Returns**
- M: the shifted mass matrix
- rvec: the vector from the center of mass to the shear center
"""
function shift_mass_matrix(M, cm, sc; check=true)
    mass = M[1, 1]
    xm_m = -M[1, 6]/mass
    ym_m = M[1, 5]/mass

    if check
        if !isapprox(cm[1], xm_m) || !isapprox(cm[2], ym_m)
            @warn("The given center of mass does not match the mass matrix center of mass.")
        end
    end

    rvec = sc .- cm

    xm = sc[1]
    ym = sc[2]

    Ixx = M[5, 5] + mass*(rvec[2]^2)
    Iyy = M[6, 6] + mass*(rvec[1]^2)
    Ixy = -M[5, 6] + mass*rvec[1]*rvec[2]

    return Symmetric([
        mass 0.0 0 0 mass*ym -mass*xm
        0 mass 0 -mass*ym 0 0
        0 0 mass mass*xm 0 0
        0 -mass*ym mass*xm Ixx+Iyy 0 0
        mass*ym 0 0 0 Ixx -Ixy
        -mass*xm 0 0 0 -Ixy Iyy
    ]), rvec
end


function laminate_loads(clt::CLT, forces, moments)
    m = length(clt.sections) #number of sections

    Pbar = zeros(4, 4)
    I = zeros(2, 4)
    F = zeros(2, 2)
    A = 0.0

    for i = 1:m #Iterate over the number of sections
        sec = clt.sections[i]
        yp = sec.y #The coordinates of the section
        zp = sec.z
        n = length(yp) #The number of points in the section

        alpha, beta, delta = laminatecompliance(sec.laminate)
        a = Symmetric([alpha[1, 1] beta[1, 1] beta[1, 3];
                 beta[1, 1] delta[1, 1] delta[1, 3];
                 beta[1, 3] delta[1, 3] delta[3, 3]]) #mu_k
        atinv = (a[2, 2]*a[3, 3] - a[2, 3]^2) / det(a)
        # Atilde = inv(atilde)  # TODO: we only need (1, 1) component so don't need to invert everything

        for k = 2:n #Iterate over the "elements" in the section
            ybar = (yp[k-1] + yp[k])/2 #The element midpoint
            zbar = (zp[k-1] + zp[k])/2
            b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2) #The length of the element
            ca = (yp[k] - yp[k-1])/b #The cosine of the angle of the element #Todo: Will this correctly orient elements? Check the math.
            sa = (zp[k] - zp[k-1])/b #The sine of the angle of the element
            #Add the element area to the cross section area #TODO: Wait... A is reset at the beginning of the function, so this is like the cumulative area of all the sections? 
            A += 0.5*(zp[k-1] + zp[k]) * (yp[k-1] - yp[k])  # if closed section (trapezoid formula for polygon area: https://en.wikipedia.org/wiki/Shoelace_formula)

            Rk = [1.0 zbar ybar 0.0;
                0.0 ca -sa 0.0;
                0.0 sa ca 0.0;
                0 0 0 1] #Todo: What is this matrix? -> It looks like a rotation matrix, but it is a 4x4 and also has the midpoints on it. 
            omega = 1.0/b*Symmetric([
                    alpha[1, 1] beta[1, 1] 0.0 -beta[1, 3]/2.0;
                    beta[1, 1] delta[1, 1] 0.0 -delta[1, 3]/2.0;
                    0 0 12.0/(atinv*b^2) 0; #Note: b changes every k iteration, so this must be created every iteration. 
                    -beta[1, 3]/2.0 -delta[1, 3]/2.0 0 delta[3, 3]/4.0])
            # omegainv = inv(omega)
            Pbar += Rk'*(omega\Rk)

            I1 = [alpha[1, 3] beta[3, 1] 0.0 -beta[3, 3]/2.0;
                  beta[1, 2] delta[1, 2] 0.0 -delta[2, 3]/2.0]
            I += I1*(omega\Rk)  # repeated, could cache

            F1 = [alpha[3, 3] beta[3, 2];
                 beta[3, 2] delta[2, 2]]
            F += b*F1 - I1*(omega\I1')

        end

    end
    L = -I
    L[1, 4] += 2*A
    if clt.closed_section
        Pbar += L'*(F\L)
    end
    Wbar = inv(Symmetric(Pbar))
    cent = -Symmetric([Wbar[2, 2] Wbar[2, 3]; Wbar[2, 3] Wbar[3, 3]]) \ [Wbar[1, 2]; Wbar[1, 3]]
    zc = cent[1]; yc = cent[2]

    Finv = inv(F)

    applied_loads = [forces[1], moments[2], moments[3], moments[1]]

    loads = ()
    

    for i in 1:m
        sec = clt.sections[i]
        yp = sec.y #The coordinates of the section
        zp = sec.z
        n = length(yp) #The number of points in the section

        alpha, beta, delta = laminatecompliance(sec.laminate)
        a = Symmetric([alpha[1, 1] beta[1, 1] beta[1, 3];
                 beta[1, 1] delta[1, 1] delta[1, 3];
                 beta[1, 3] delta[1, 3] delta[3, 3]]) #mu_k

        nu_k = [ alpha[1, 3] beta[1, 2];
                 beta[3, 1] delta[1, 2];
                 beta[3, 3] delta[2, 3] ]


        mu_inv = inv(a)

        section_loads = zeros(6, n-1)
        for k in 2:n
            Nen, Mn = Finv*L*Wbar*applied_loads

            ybar = (yp[k-1] + yp[k])/2 #The element midpoint
            zbar = (zp[k-1] + zp[k])/2
            b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2) #The length of the element
            ca = (yp[k] - yp[k-1])/b #The cosine of the angle of the element #Todo: Will this correctly orient elements? Check the math.
            sa = (zp[k] - zp[k-1])/b #The sine of the angle of the element 

            Rk = [1.0 zbar ybar 0.0;
                0.0 ca -sa 0.0;
                0.0 sa ca 0.0;
                0 0 0 1]

            # eta = zbar
            eta = 0.0
            Reta = [1.0 0.0 eta 0.0;
                    0.0 1.0 0.0 0.0;
                    0.0 0.0 0.0 -2.0]


            Ne, Me, Men = mu_inv*(Reta*Rk - nu_k*Finv*L)*Wbar*applied_loads

            section_loads[:, k-1] = [0.0, Ne, Nen, Mn, Me, Men]
        end

        loads = (loads..., section_loads)
    end
    
    return loads
end


# ---- alternative (simpler) methods for orthotropic, no longer used ------------
function centroid(beamsections)

    m = length(beamsections)

    ynum = 0.0
    znum = 0.0
    den = 0.0
    for i = 1:m
        sec = beamsections[i]
        yp = sec.y
        zp = sec.z
        n = length(yp)

        alpha, beta, delta = laminatecompliance(sec.laminate)
        D = alpha[1, 1]*delta[1, 1] - beta[1, 1]^2

        s = 0.0
        ybar = 0.0
        zbar = 0.0
        ca = 0.0
        sa = 0.0
        for j = 2:n
            ds = sqrt((yp[j] - yp[j-1])^2 + (zp[j] - zp[j-1])^2)
            ym = (yp[j] + yp[j-1])/2
            zm = (zp[j] + zp[j-1])/2
            s += ds
            ca += (yp[j] - yp[j-1])  # / ds * ds
            sa += (zp[j] - zp[j-1])  # / ds * ds
            ybar += ym * ds
            zbar += zm * ds
        end

        ynum += delta[1, 1]/D*ybar + beta[1, 1]/D*sa
        znum += delta[1, 1]/D*zbar - beta[1, 1]/D*ca
        den += delta[1, 1]/D*s
    end

    yc = ynum/den
    zc = znum/den


    return yc, zc
end

function beamstiffnessold(beamsections)

    yc, zc = centroid(beamsections)

    m = length(beamsections)

    EA = 0.0
    EIyy = 0.0
    EIzz = 0.0
    EIyz = 0.0
    A = 0.0
    GJden = 0.0
    for i = 1:m
        sec = beamsections[i]
        yp = sec.y
        zp = sec.z
        n = length(yp)

        alpha, beta, delta = laminatecompliance(sec.laminate)  # different sections should generally have different laminates for a closed section (otherwise combine into one section), so no loss in efficiency here
        a11hat = alpha[1, 1] - beta[1, 2]^2/delta[2, 2]
        b11hat = beta[1, 1] - beta[1, 2]*delta[1, 2]/delta[2, 2]
        d11hat = delta[1, 1] - delta[1, 2]^2/delta[2, 2]
        Dhat = a11hat*d11hat - b11hat^2
        D = alpha[1, 1]*delta[1, 1] - beta[1, 1]^2
        alpha_nu = alpha[3, 3] - beta[3, 3]^2/delta[3, 3]

        s = 0.0
        z2 = 0.0
        y2 = 0.0
        yz = 0.0
        zca = 0.0
        zsa = 0.0
        yca = 0.0
        ysa = 0.0
        c2 = 0.0
        s2 = 0.0
        cs = 0.0
        Asub = 0.0
        for j = 2:n
            ds = sqrt((yp[j] - yp[j-1])^2 + (zp[j] - zp[j-1])^2)
            ca = (yp[j] - yp[j-1])/ds
            sa = (zp[j] - zp[j-1])/ds
            ym = (yp[j] + yp[j-1])/2
            zm = (zp[j] + zp[j-1])/2
            y = ym - yc
            z = zm - zc
            s += ds
            z2 += z^2 * ds + ds^3/12*sa^2
            y2 += y^2 * ds + ds^3/12*ca^2
            yz += y*z * ds + ds^3/12*sa*ca
            zca += z*ca * ds
            zsa += z*sa * ds
            yca += y*ca * ds
            ysa += y*sa * ds
            c2 += ca*ca * ds
            s2 += sa*sa * ds
            cs += ca*sa * ds
            Asub += (zp[j-1] + zp[j]) * (yp[j-1] - yp[j])  # if closed section
        end

        EA += d11hat/Dhat * s  # for a closed section
        EIyy += delta[1, 1]/D*z2 - 2*beta[1, 1]/D*zca + alpha[1, 1]/D*c2
        EIzz += delta[1, 1]/D*y2 + 2*beta[1, 1]/D*ysa + alpha[1, 1]/D*s2
        EIyz += delta[1, 1]/D*yz + beta[1, 1]/D*(zsa - yca) - alpha[1, 1]/D*cs
        GJden += s*alpha_nu
        A += abs(Asub)/2.0
    end

    GJ = 4*A^2/GJden

    return EA, EIyy, EIzz, EIyz, GJ
end

# --------------------------------------




# --------- strains ----------------




"""
from double prime (e1, e2, e3) to primed c.s.
"""
function clt_Ttheta(theta)

    s, c = sincos(theta)

    c2 = c^2
    s2 = s^2
    sc = s*c

    Tsigma = @SMatrix [
        c2    s2    0  -2*sc  0  0;
        s2    c2    0   2*sc  0  0;
        0     0     1   0     0  0;
        sc   -sc    0  c2-s2  0  0;
        0     0     0   0     c -s;
        0     0     0   0     s  c;
    ]

    Teps = @SMatrix [
        c2    s2    0  -sc   0  0;
        s2    c2    0   sc   0  0;
        0     0     1   0    0  0;
        2*sc -2*sc  0  c2-s2 0  0;
        0     0     0   0    c -s;
        0     0     0   0    s  c;
    ]

    return Tsigma, Teps
end


"""
from primed to unprimed
"""
function clt_Talpha(c, s)

    # s, c = sincos(alpha)

    c2 = c^2
    s2 = s^2
    sc = s*c

    Tsigma = @SMatrix [
        1  0    0  0  0    0;
        0  c2  s2  0  0  -2*sc;
        0  s2  c2  0  0   2*sc;
        0  0    0  c  -s   0;
        0  0    0  s  c    0;
        0  sc  -sc 0  0   c2-s2
    ]

    Teps = @SMatrix [
        1  0      0   0  0    0;
        0  c2    s2   0  0   -sc;
        0  s2    c2   0  0    sc;
        0  0      0   c  -s   0;
        0  0      0   s  c    0;
        0  2*sc -2*sc 0  0   c2-s2
    ]

    return Tsigma, Teps
end


# """
# from ply to beam
# """
# function rotate_stress_and_strains(sigmap, epsilonp, theta, cosalpha, sinalpha)
#     Tsigma1, Teps1 = clt_Ttheta(theta)
#     Tsigma2, Teps2 = clt_Talpha(cosalpha, sinalpha)
#     println(sigmap)
#     println(epsilonp)
#     sigmab = Tsigma2 * Tsigma1 * sigmap
#     epsilonb = Teps2 * Teps1 * epsilonp

#     return sigmab, epsilonb
# end

"""
from laminate to beam
"""
function rotate_stress_and_strains(sigmaprime, epsilonprime, cosalpha, sinalpha)
    Tsigma2, Teps2 = clt_Talpha(cosalpha, sinalpha)
    sigmab = Tsigma2 * sigmaprime
    epsilonb = Teps2 * epsilonprime

    return sigmab, epsilonb
end
# TODO: create general functions for both methods (CLT and FEA)

function num_strain_locs(clt::CLT)
    ntotal = 0
    for i in eachindex(clt.sections)
        ntotal += (length(clt.sections[i].y) - 1) * 2*length(clt.sections[i].laminate)
    end
    return ntotal
end


"""
    strains_and_stresses(F, M, clt::CLT)

Compute strains and stresses in the beam sections.

**Arguments**
- F::Vector{TF}: forces in the beam
- M::Vector{TF}: moments in the beam
- clt::CLT: composite section

**Returns**
- `strain_b::Vector(6, ne)`: strains in beam coordinate system for each element. order: xx, yy, zz, xy, xz, yz
- `stress_b::Vector(6, ne)`: stresses in beam coordinate system for each element. order: xx, yy, zz, xy, xz, yz
- `strain_p::Vector(6, ne)`: strains in ply coordinate system for each element. order: 11, 22, 33, 12, 13, 23
- `stress_p::Vector(6, ne)`: stresses in ply coordinate system for each element. order: 11, 22, 33, 12, 13, 23
"""
function strains_and_stresses(F, M, clt::CLT)
    # map GXBeam forces to internal order
    Nxbar = F[1]  # deformations due to shear neglected in this method
    Txbar, Mybar, Mzbar = M
    Mzbar *= -1  # opposite sign convention used internally
    FMvec = [Nxbar; Mybar; Mzbar; Txbar]

    # rename for convenience
    cc = clt.cache

    # number of sections
    m = length(clt.sections)

    # count how many locations I have to compute strain at
    ntotal = 0
    for i = 1:m
        ntotal += (length(clt.sections[i].y) - 1) * 2*length(clt.sections[i].laminate)
    end

    ### Preallocate
    #Get the typing for the arrays
    TF = promote_type(typeof(F[1]), typeof(M[1]), typeof(clt.sections[1].y[1]))

    strain_p = zeros(TF, 6, ntotal) 
    stress_p = zeros(TF, 6, ntotal)
    strain_b = zeros(TF, 6, ntotal)
    stress_b = zeros(TF, 6, ntotal)

    idx = 1
    for sec in clt.sections #Iterate across the sectinos
        
        yp = sec.y
        zp = sec.z
        n = length(yp)

        alpha, beta, delta = laminatecompliance(sec.laminate)
        muk = Symmetric([alpha[1, 1] beta[1, 1] beta[1, 3];
            beta[1, 1] delta[1, 1] delta[1, 3]
            beta[1, 3] delta[1, 3] delta[3, 3]])
        nuk = [alpha[1, 3] beta[1, 2]
            beta[3, 1] delta[1, 2]
            beta[3, 3] delta[2, 3]]

        for k = 2:n #Iterate across the elements in the section
            ybar = (yp[k-1] + yp[k])/2
            zbar = (zp[k-1] + zp[k])/2
            b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2)
            ca = (yp[k] - yp[k-1])/b
            sa = (zp[k] - zp[k-1])/b

            # Rk = [1.0 zbar ybar 0.0;
            #     0.0 ca -sa 0.0;
            #     0.0 sa ca 0.0;
            #     0 0 0 1] #todo: a significant number of allocations

            eta = 0.0  # just compute at midpoint
            # Reta = [1.0 0 eta 0;
            #         0 1 0 0;
            #         0 0 0 -2] #todo: a significant number of allocations

            RetaRk = @SMatrix [1.0 zbar+(sa*eta) ybar+(ca*eta) 0.0;
                       0.0 ca -sa 0.0;
                       0.0 0.0 0.0 -2.0] #todo: a significant number of allocations

            # pg 270
            if clt.closed_section
                NM1 = cc.F\cc.L*cc.Wbar*FMvec
                # NM2 = (muk\(Reta*Rk - nuk*(cc.F\cc.L)))*cc.Wbar*FMvec
                NM2 = (muk\(RetaRk - nuk*(cc.F\cc.L)))*cc.Wbar*FMvec #todo: A significant number of allocations
                forces = [NM2[1]; 0.0; NM1[1]; NM2[2]; NM1[2]; NM2[3]]  # N1, N2, N12, M1, M2, M12 #Note: coud probably just directly multiply out the elements and form the forces vector directly. #todo: A significant number of allocations. 
            else
                # NM = muk\Reta*Rk*cc.Wbar*FMvec
                NM = muk\RetaRk*cc.Wbar*FMvec
                forces = [NM[1]; 0.0; 0.0; NM[2]; 0.0; NM[3]]
            end

            _, _, _, epsilonprime = laminatestrains(alpha, beta, delta, sec.laminate, forces)
            sigmaprime, sigma_p, epsilon_p = laminatestresses(sec.laminate, epsilonprime)


            # remap from internal representation to common representation
            # epsilonp is 11, 22, 12 (other 3 components are zero)  TODO: 33 is actually not zero
            # strain_p::Vector(6, nloc)`: strains in ply coordinate system for each element. order: 11, 22, 33, 12, 13, 23
            _, nz = size(epsilonprime)
            strain_p[1, idx:idx+nz-1] = epsilon_p[1, :]
            strain_p[2, idx:idx+nz-1] = epsilon_p[2, :]
            strain_p[4, idx:idx+nz-1] = epsilon_p[3, :]

            stress_p[1, idx:idx+nz-1] = sigma_p[1, :]
            stress_p[2, idx:idx+nz-1] = sigma_p[2, :]
            stress_p[4, idx:idx+nz-1] = sigma_p[3, :]

            # iz = [1, 1]
            # for iii = 2:length(sec.laminate)*2
            #     iz = [iz; iii; iii]
            # end

            # for ir = idx:idx+nz-1
            #     thetak = sec.laminate[iz[ir - idx + 1]].theta
            #     stress_b[:, ir], strain_b[:, ir] = rotate_stress_and_strains(stress_p[:, ir], strain_p[:, ir], thetak, ca, sa)
            # end

            sigma_temp = zeros(TF, 6)
            epsilon_temp = zeros(TF, 6)
            for ir = idx:idx+nz-1
                sigma_temp[[1, 2, 4]] .= sigmaprime[:, ir-idx+1]
                epsilon_temp[[1, 2, 4]] .= epsilonprime[:, ir-idx+1]
                stress_b[:, ir], strain_b[:, ir] = rotate_stress_and_strains(sigma_temp, epsilon_temp, ca, sa)
            end

            idx += nz

        end
    end


    return strain_b, stress_b, strain_p, stress_p
end


"""
    tsai_wu(stress, clt::CLT)

Compute Tsai-Wu failure criterion for a composite section.

**Arguments**
- stress::Matrix{TF}: 6 stresses in the ply coordinate system. Stresses should be calculated at the top and bottom of each ply (as in `strains_and_stresses`).
- clt::CLT: the composite section

**Returns**
- failure::Vector{TF}: Tsai-Wu failure criterion for the top and bottom of each ply (where the stresses are calculated)
"""
function tsai_wu(stress, clt::CLT)

    m = length(clt.sections)

    # count how many locations I have to compute failure at
    ntotal = 0
    for i = 1:m
        ntotal += (length(clt.sections[i].y) - 1) * 2*length(clt.sections[i].laminate)
    end

    T = eltype(stress)
    failure = Vector{T}(undef, ntotal)

    stress_idx = [1, 2, 4] #Stresses map (1, 2, 4) -> (1, 2, 3)
    idx = 1
    for i = 1:m
        sec = clt.sections[i]
        for j = 1:length(sec.y)-1
            nl = length(sec.laminate)
            idxs = idx:idx+2*nl-1

            # failure[idxs] = tsai_wu_plane(stress[stress_idx, idxs], sec.laminate)
            failure[idxs] = full_tsai_wu(stress[:, idxs], sec.laminate)

            idx += 2*nl
        end
    end

    return failure
end

"""
    tsai_hill(stress, clt::CLT)

Compute Tsai-Hill failure criterion for a composite section.

**Arguments**
- stress::Matrix{TF}: 6 stresses in the ply coordinate system. Stresses should be calculated at the top and bottom of each ply (as in `strains_and_stresses`).
- clt::CLT: the composite section

**Returns**
- failure::Vector{TF}: Tsai-Wu failure criterion for the top and bottom of each ply (where the stresses are calculated)
"""
function tsai_hill(stress, clt::CLT)

    m = length(clt.sections)

    # count how many locations I have to compute failure at
    ntotal = 0
    for i = 1:m
        ntotal += (length(clt.sections[i].y) - 1) * 2*length(clt.sections[i].laminate)
    end

    T = eltype(stress)
    failure = Vector{T}(undef, ntotal)

    stress_idx = [1, 2, 4] #Stresses map (1, 2, 4) -> (1, 2, 3)
    idx = 1
    for i = 1:m
        sec = clt.sections[i]
        for j = 1:length(sec.y)-1
            nl = length(sec.laminate)
            idxs = idx:idx+2*nl-1

            failure[idxs] = tsai_hill(stress[stress_idx, idxs], sec.laminate)

            idx += 2*nl
        end
    end

    return failure
end

"""
    max_stress(stress, clt::CLT)

Compute the maximum stress failure criterion for a composite section.

**Arguments**
- stress::Matrix{TF}: 6 stresses in the ply coordinate system. Stresses should be calculated at the top and bottom of each ply (as in `strains_and_stresses`).
- clt::CLT: the composite section

**Returns**
- failure::Vector{TF}: Tsai-Wu failure criterion for the top and bottom of each ply (where the stresses are calculated)
"""
function max_stress(stress, clt::CLT)

    m = length(clt.sections)

    # count how many locations I have to compute failure at
    ntotal = 0
    for i = 1:m
        ntotal += (length(clt.sections[i].y) - 1) * 2*length(clt.sections[i].laminate)
    end

    T = eltype(stress)
    failure = Array{T, 2}(undef, 3, ntotal)

    stress_idx = [1, 2, 4] #Stresses map (1, 2, 4) -> (1, 2, 3)
    idx = 1
    for i = 1:m
        sec = clt.sections[i]
        for j = 1:length(sec.y)-1
            nl = length(sec.laminate)
            idxs = idx:idx+2*nl-1

            failure[:, idxs] = max_stress(stress[stress_idx, idxs], sec.laminate)

            idx += 2*nl
        end
    end

    return failure
end

"""
    buckling(clt::CLT, strain, E_axial)

Compute the buckling failure criterion for a composite section.

**Arguments**
- clt::CLT: the composite section
- strain::Matrix{TF}: 6 strains in the ply coordinate system. Strains should be calculated at the top and bottom of each ply (as in `strains_and_stresses`)
- E_axial::TF: axial modulus of elasticity

**Returns**
- failure::Vector{TF}: buckling failure criterion for the top and bottom of each ply (where the stresses are calculated)
"""
function buckling(clt::CLT, strain, E_axial)

    m = length(clt.sections)

    # count how many locations I have to compute failure at
    ntotal = 0
    for i = 1:m
        # ntotal += (length(clt.sections[i].y) - 1) * 2*length(clt.sections[i].laminate)
        ntotal += (length(clt.sections[i].y) - 1) #* 2*length(clt.sections[i].laminate)
    end

    T = eltype(strain)
    failure = Vector{T}(undef, ntotal)

    strain_idx = [1, 2, 4] #Strain map (1, 2, 4) -> (1, 2, 3)
    idx = 1
    for i = 1:m #Iterate over the sections
        sec = clt.sections[i]
        for j = 1:length(sec.y)-1 #Iterate over the elements in the section
            nl = length(sec.laminate)
            idxs = idx:idx+2*nl-1

            x = sec.y[j+1] - sec.y[j] #The x distance between the element end points
            y = sec.z[j+1] - sec.z[j] #The y distance between the element end points
            L = sqrt(x^2 + y^2) #The length of the element

            # failure[idxs] = buckling_strain(strain[strain_idx, idxs], sec.laminate, L, E_axial)
            failure[idx] = buckling_strain(strain[strain_idx, idxs], sec.laminate, L) #RotorSE approach
            #https://github.com/WISDEM/RotorSE/blob/d497e96a612b56aa37b57d35c64aa5049bb9136f/src/rotorse/rotor.py

            # idx += 2*nl
            idx += 1
        end
    end

    return failure
end

function max_strain(clt::CLT, strain, strain_ult)
end



"""
    interpolate_load(load, clt::CLT, x, y)

Find the strain, stress, or failure at a specific point in the cross section.

**Arguments**
- load::Vector{TF}: the load (strain, stress, or failure) at the top and bottom of each ply (as in `strains_and_stresses` or `tsai_wu`)
- clt::CLT: the composite section
- x::TF: x-coordinate of the point
- y::TF: y-coordinate of the point
- fit::Function: interpolation function (default: linear)

**Returns**
- load_new: the interpolated value at the point
"""
function interpolate_load(load, clt::CLT, x, y; fit=linear)


    #= 
        - I need to find the section that the point is in
        - How do I guarantee that the point is in the cross section? (Not floating in empty space or outside the shape)
            - Floating outside the empy shape should be easy... right? 
            - Well I guess that the defining points are hidden away in the sections of CLT. I could probably iterate through, but I don't know if that would tell me anything about the actual shape. How do I know what order they should be in. 
            - I could probably just make a 2D unstructured grid then interpolate based on that? 

            - Probably just easier to visually inspect at this point. 
    =#
    
end


# ---------- shear flow --------------------

function qopen(eta, y1, y2, z1, z2, m1, m2)
    ys = y1 + eta*(y2 - y1)
    zs = z1 + eta*(z2 - z1)
    ybarsub = (y1 + ys)/2
    zbarsub = (z1 + zs)/2
    bsub = sqrt((ys - y1)^2 + (zs - z1)^2)
    qo = m1*zbarsub*bsub + m2*ybarsub*bsub
    return qo
end

function shearflow(sections, P, yc, zc; closedsection=true, npts=20)

    EIyy = P[2, 2]
    EIzz = P[3, 3]
    EIyz = P[2, 3]

    eta = range(0.0, 1.0, length=npts)

    function factors(Vy, Vz)
        factor1 = (-EIzz*Vz + EIyz*Vy)/(EIyy*EIzz - EIyz^2)
        factor2 = (EIyz*Vz - EIyy*Vy)/(EIyy*EIzz - EIyz^2)
        return factor1, factor2
    end
    factor1y, factor2y = factors(1.0, 0.0)
    factor1z, factor2z = factors(0.0, 1.0)


    m = length(sections)


    qcy = 0.0
    qcz = 0.0
    if closedsection
        # ------ compute qc -------------
        qoprevy = 0.0
        qoprevz = 0.0
        qcy_num = 0.0
        qcy_den = 0.0
        qcz_num = 0.0
        qcz_den = 0.0
        for i = 1:m
            sec = sections[i]
            yp = sec.y .- yc  # relative to centroid
            zp = sec.z .- zc
            n = length(yp)

            alpha, beta, delta = laminatecompliance(sec.laminate)
            a11 = alpha[1, 1]
            a_nu = alpha[3, 3] - beta[3, 3]^2/delta[3, 3]

            m1y = factor1y/a11
            m2y = factor2y/a11
            m1z = factor1z/a11
            m2z = factor2z/a11

            for k = 2:n
                b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2)

                qo(s, qoprev, m1, m2) = qoprev + qopen(s, yp[k-1], yp[k], zp[k-1], zp[k], m1, m2)

                # integraly, _ = quadgk(s -> qo(s, qoprevy, m1y, m2y), 0.0, b)
                # integralz, _ = quadgk(s -> qo(s, qoprevz, m1z, m2z), 0.0, b)
                integraly = b*trapz(eta, qo.(eta, qoprevy, m1y, m2y))
                integralz = b*trapz(eta, qo.(eta, qoprevz, m1z, m2z))
                qcy_num += a_nu*integraly
                qcy_den += a_nu*b
                qcz_num += a_nu*integralz
                qcz_den += a_nu*b

                qoprevy = qo(1.0, qoprevy, m1y, m2y)
                qoprevz = qo(1.0, qoprevz, m1z, m2z)
            end
        end

        qcy = -qcy_num/qcy_den
        qcz = -qcz_num/qcz_den
    end


    # ------ compute ysc, zsc, s, S -------------
    qoprevy = 0.0
    qoprevz = 0.0
    ysc = 0.0
    zsc = 0.0
    syy = 0.0
    szz = 0.0
    syz = 0.0

    for i = 1:m
        sec = sections[i]
        yp = sec.y .- yc  # relative to centroid
        zp = sec.z .- zc
        n = length(yp)

        alpha, beta, delta = laminatecompliance(sec.laminate)
        a11 = alpha[1, 1]
        a_nu = alpha[3, 3] - beta[3, 3]^2/delta[3, 3]

        m1y = factor1y/a11
        m2y = factor2y/a11
        m1z = factor1z/a11
        m2z = factor2z/a11

        for k = 2:n
            b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2)

            qo(s, qoprev, m1, m2) = qoprev + qopen(s, yp[k-1], yp[k], zp[k-1], zp[k], m1, m2)
            qy(s) = qcy + qo(s, qoprevy, m1y, m2y)
            qz(s) = qcz + qo(s, qoprevz, m1z, m2z)
            function p(s)
                ys = yp[k-1] + s*(yp[k] - yp[k-1])
                zs = zp[k-1] + s*(zp[k] - zp[k-1])
                return sqrt(ys^2 + zs^2)
            end

            # integraly, _ = quadgk(s -> qy(s)*p(s), 0.0, b)
            # integralz, _ = quadgk(s -> qz(s)*p(s), 0.0, b)
            integraly = b*trapz(eta, @. qy(eta)*p(eta))
            integralz = b*trapz(eta, @. qz(eta)*p(eta))

            zsc -= integraly
            ysc += integralz


            # integral, _ = quadgk(s -> qy(s)^2, 0.0, b)
            integral = b*trapz(eta, @. qy(eta)^2)
            syy += a_nu*integral
            # integral, _ = quadgk(s -> qz(s)^2, 0.0, b)
            integral = b*trapz(eta, @. qz(eta)^2)
            szz += a_nu*integral
            # integral, _ = quadgk(s -> qy(s)*qz(s), 0.0, b)
            integral = b*trapz(eta, @. qy(eta)*qz(eta))
            syz += a_nu*integral

            qoprevy = qo(1.0, qoprevy, m1y, m2y)
            qoprevz = qo(1.0, qoprevz, m1z, m2z)
        end
    end

    s = Symmetric([syy syz;
                   syz szz])
    S = inv(s)

    return s, S, ysc, zsc
end



function shearflow_general_attempt(sections, Wbar, F, L, yc, zc)

    m = length(sections)

    TF = get_section_floattype(sections[1])

    qoy = 0.0
    qoz = 0.0
    qcy = 0.0
    qcz = 0.0
    eta = 0.0
    AMy = zeros(TF, 2, 2)
    bvy = zeros(TF, 2)
    AMz = zeros(TF, 2, 2)
    bvz = zeros(TF, 2)
    for i = 1:m
        sec = sections[i]
        yp = sec.y
        zp = sec.z
        n = length(yp)

        alpha, beta, delta = laminatecompliance(sec.laminate)
        muk = Symmetric([alpha[1, 1] beta[1, 1] beta[1, 3];
                            beta[1, 1] delta[1, 1] delta[1, 3]
                            beta[1, 3] delta[1, 3] delta[3, 3]])
        nuk = [alpha[1, 3] beta[1, 2]
            beta[3, 1] delta[1, 2]
            beta[3, 3] delta[2, 3]]

        for k = 2:n
            ybar = (yp[k-1] + yp[k])/2
            zbar = (zp[k-1] + zp[k])/2
            b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2)
            ca = (yp[k] - yp[k-1])/b
            sa = (zp[k] - zp[k-1])/b

            # open section shear flow
            eta = b/2.0
            Rk = [1.0 zbar ybar 0.0;
                0.0 ca -sa 0.0;
                0.0 sa ca 0.0;
                0 0 0 1]
            Reta = [1.0 0 eta 0;  # TODO: this should be local eta, so just b/2 at the center I think
                    0 1 0 0;
                    0 0 0 -2]

            M = (muk\(Reta*Rk - nuk*(F\L)))*Wbar

            # --------
            # @show typeof(M) #Matrix{ForwardDiff.Dual{ForwardDiff.Tag{var"#objwrap#92"{Int64}, Float64}, Float64, 12}}
            # @show typeof(Vy) #Defined below
            # @show typeof(Vz) #Defined below
            # @show b # Dual
            # @show alpha #Matrix{Floats}
            # @show beta #Matrix{Floats}
            # @show delta #Matrix{Floats}
            # @show qoy #TF
            # @show AMy #Matrix{Floats}
            # @show bvy #Vector{Floats}


            Vy = 1.0; Vz = 0.0
            qoy = shearflowsub!(M, Vy, Vz, b, alpha, beta, delta, qoy, AMy, bvy)

            Vy = 0.0; Vz = 1.0
            qoz = shearflowsub!(M, Vy, Vz, b, alpha, beta, delta, qoz, AMz, bvz)
        end
    end

    qcy = 0.0  #(AMy\bvy)[1]
    qcz = 0.0  #(AMz\bvz)[1]

    eta = 0.0
    qoy = 0.0
    qoz = 0.0
    ysc = 0.0
    zsc = 0.0
    syy = 0.0
    szz = 0.0
    syz = 0.0
    for i = 1:m
        sec = sections[i]
        yp = sec.y
        zp = sec.z
        n = length(yp)

        alpha, beta, delta = laminatecompliance(sec.laminate)
        muk = Symmetric([alpha[1, 1] beta[1, 1] beta[1, 3];
                            beta[1, 1] delta[1, 1] delta[1, 3]
                            beta[1, 3] delta[1, 3] delta[3, 3]])
        nuk = [alpha[1, 3] beta[1, 2]
            beta[3, 1] delta[1, 2]
            beta[3, 3] delta[2, 3]]

        alpha_nu = alpha[3, 3] - beta[3, 3]^2/delta[3, 3]

        for k = 2:n
            ybar = (yp[k-1] + yp[k])/2
            zbar = (zp[k-1] + zp[k])/2
            b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2)
            ca = (yp[k] - yp[k-1])/b
            sa = (zp[k] - zp[k-1])/b

            # open section shear flow
            eta = b/2.0
            Rk = [1.0 zbar ybar 0.0;
                0.0 ca -sa 0.0;
                0.0 sa ca 0.0;
                0 0 0 1]
            Reta = [1.0 0 eta 0;  # TODO: this should be local eta, so just b/2 at the center I think
                    0 1 0 0;
                    0 0 0 -2]

            M = (muk\(Reta*Rk - nuk*(F\L)))*Wbar

            Vy = 1.0; Vz = 0.0
            a = M*[0.0; Vz; Vy; 0]
            dN = a[1]
            qoy -= dN*b
            qy = qoy + qcy

            Vy = 0.0; Vz = 1.0
            a = M*[0.0; Vz; Vy; 0]
            dN = a[1]
            qoz -= dN*b
            qz = qoz + qcz

            # shear center
            p = sqrt((ybar - yc)^2 + (zbar - zc)^2)
            ysc += qz*p*b
            zsc -= qy*p*b

            # stiffnesses
            syy += alpha_nu*qy^2*b
            szz += alpha_nu*qz^2*b
            syz += alpha_nu*qy*qz*b
        end
    end

    s = Symmetric([syy syz;
                   syz szz])
    S = inv(s)

    return s, S, ysc, zsc
end

function shearflowsub!(M, Vy, Vz, b, alpha, beta, delta, qo, AM, bv)

    a = M*[0.0; Vz; Vy; 0]
    dN = a[1]
    qo -= dN*b

    # closed section due to Vy/Vz
    AM .-= b*Symmetric([alpha[3, 3] beta[2, 3]
                      beta[2, 3] delta[2, 2]])
    bv .+= qo*b*[alpha[3, 3]; beta[2, 3]]

    return qo
end


# -------------------------------------------------------------

function get_laminate_points(laminate::Array{TL, 1}, x, y) where {TL<:Layer}
    nhat = [y[1] - y[2], x[1] - x[2]] # Normal to the layup. 
    nhat = nhat./norm(nhat)

    nl = length(laminate)
    np = 2nl + 2
    xp = zeros(np)
    yp = zeros(np)

    xp[1:2] = x[1:2]
    yp[1:2] = y[1:2]

    idx = 1
    for i in 3:2:np
        xp[i] = xp[i-2] + nhat[1]*laminate[idx].t
        xp[i+1] = xp[i-1] + nhat[1]*laminate[idx].t

        yp[i] = yp[i-2] + nhat[2]*laminate[idx].t
        yp[i+1] = yp[i-1] + nhat[2]*laminate[idx].t
        idx += 1
    end

    return xp, yp
end

@recipe function plot_laminate(laminate::Array{TL, 1}, x, y) where {TL<:Layer} #Todo: I think I'm going to move the majority of this into a function so it doesn't get repeated so much. 

    xp, yp = get_laminate_points(laminate, x, y)

    # Loop through each layer and plot. 
    for i in 1:nl
        @series begin
            j = 2*(i-1)+1
            idxs = [j, j+1, j+3, j+2, j]
        
            xr = xp[idxs]
            yr = yp[idxs]

            label --> false
            seriescolor --> :black

            xr, yr
        end
    end
end # End recipe


@recipe function plot_clt(clt::CLT; highlight_idxs=[], highlight=:red)


    idx = 1 #Element counter
    for i in eachindex(clt.sections)
        sec = clt.sections[i]
        ns = length(sec.y)

        for j in 1:ns-1 #Iterate over the number of elements in the section
            x = sec.y[j+1] - sec.y[j] #The x distance between the element end points
            y = sec.z[j+1] - sec.z[j] #The y distance between the element end points
            L = sqrt(x^2 + y^2)


            #Normal vector to the element
            nx = -y/L
            ny = x/L

            # @show nx, ny

            T = 0.0 #The total thickness travelled so far. 

            for k in eachindex(sec.laminate) #Iterate over the laminates in the section
                t = sec.laminate[k].t #The thickness of the laminate
                
                xp = [sec.y[j] + nx*T, sec.y[j+1] + nx*T, sec.y[j+1] + nx*(T+t), sec.y[j] + nx*(T+t), sec.y[j] + nx*T]
                yp = [sec.z[j] + ny*T, sec.z[j+1] + ny*T, sec.z[j+1] + ny*(T+t), sec.z[j] + ny*(T+t), sec.z[j] + ny*T]

                @series begin
                    label --> false
                    if in(idx, highlight_idxs)
                        seriescolor --> highlight
                        markershape --> :circle
            
                    else
                        seriescolor --> :black
                    end

                    xp, yp
                end
                
                idx += 1 #Increment the element counter
                T += t #Increment the thickness
            end #End looping over laminates
        end #end looping over segment elements
    end #End looping over sections

end #End recipe


"""

Notes: Check https://docs.juliaplots.org/latest/generated/colorschemes/#matplotlib for color schemes. 
"""
@recipe function plot_clt_solution(clt::CLT, solution::Array{TF, 1}, cgcolor; double_sol=true, highlight_idxs=[], highlight=:red, loval=minimum(solution), hival=maximum(solution), diffun=nothing) where {TF}

    colorbar --> true
    fill --> true
    clim --> (loval, hival)

    if isnothing(diffun)
        diffun = x -> (x[1] + x[2])/2
    end

    idx = 1 #Element counter
    sol_idx = 1
    for i in eachindex(clt.sections)
        sec = clt.sections[i]
        ns = length(sec.y)

        for j in 1:ns-1 #Iterate over the number of elements in the section
            x = sec.y[j+1] - sec.y[j] #The x distance between the element end points
            y = sec.z[j+1] - sec.z[j] #The y distance between the element end points
            L = sqrt(x^2 + y^2)


            #Normal vector to the element
            nx = -y/L
            ny = x/L

            # @show nx, ny

            T = 0.0 #The total thickness travelled so far. 

            for k in eachindex(sec.laminate) #Iterate over the laminates in the section
                t = sec.laminate[k].t #The thickness of the laminate
                
                xp = [sec.y[j] + nx*T, sec.y[j+1] + nx*T, sec.y[j+1] + nx*(T+t), sec.y[j] + nx*(T+t), sec.y[j] + nx*T]
                yp = [sec.z[j] + ny*T, sec.z[j+1] + ny*T, sec.z[j+1] + ny*(T+t), sec.z[j] + ny*(T+t), sec.z[j] + ny*T]
                if double_sol
                    # val = sum(solution[sol_idx:sol_idx+1])/2
                    val = diffun(solution[sol_idx:sol_idx+1])
                else
                    val = solution[idx]
                end

                @series begin
                    label --> false
                    if in(idx, highlight_idxs)
                        seriescolor --> highlight
                        markershape --> :circle
            
                    else
                        seriescolor --> :black
                    end

                    fc --> cgcolor
                    fill_z --> val

                    xp, yp
                end
                
                idx += 1 #Increment the element counter
                sol_idx += 2
                T += t #Increment the thickness
            end #End looping over laminates
        end #end looping over segment elements
    end #End looping over sections

end #End recipe