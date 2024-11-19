
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

    c = cos(lamina.theta)
    s = sin(lamina.theta)
    c2 = c*c
    s2 = s*s
    cs = c*s

    Q = Symmetric([E1*delta  nu12*E2*delta  0.0;
                  nu12*E2*delta  E2*delta  0.0;
                  0.0  0.0  G12])
    Ts = [c2 s2 -2*cs;
          s2 c2 2*cs;
          cs -cs c2-s2]

    Te = [c2 s2 -cs;
          s2 c2 cs;
          2*cs -2*cs c2-s2]

    Qbar = Symmetric(Ts * Q * Ts')

    return Qbar, Ts, Te
end

"""
compute the z locations and overall height for each lamina in the laminate.
"""
function zspacing(laminate)

    n = length(laminate)

    # compute z vector
    z = zeros(n+1)
    for i = 1:n
        z[i+1] = z[i] + laminate[i].t
    end
    h = z[end] - z[1]
    mid = h/2.0  # compute midpoint
    z .-= mid  # recenter at midpoint

    return z, h
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
Compute the alpha, beta, delta compliance matrices for a thin laminate given the stiffness matrices
"""
function laminatecompliancematrix(A, B, D)
    Ainv = inv(A)
    Hinv = Symmetric(inv(D - B*Ainv*B))
    alpha = Symmetric(Ainv + Ainv*B*Hinv*B*Ainv)
    beta = -Ainv*B*Hinv
    delta = Hinv

    return alpha, beta, delta
end

"""
convenicne method to go from laminate definition to its compliance matrices
"""
function laminatecompliance(laminate)
    z, h = zspacing(laminate)
    A, B, D = laminatestiffnessmatrix(laminate, z)
    alpha, beta, delta = laminatecompliancematrix(A, B, D)
    return alpha, beta, delta
end

"""
strains for a thin laminate with forces: N1, N1, N12, M1, M2, M12
"""
function strains(alpha, beta, delta, z, forces)
    C = [alpha beta; beta' delta]

    alleps = C*forces
    epsilonbar = alleps[1:3]
    kappa = alleps[4:6]

    # setup new z vector at top and bottom of each ply
    nz = 2*(length(z)-1)
    zvec = zeros(nz)
    zvec[1] = z[1]
    zvec[end] = z[end]
    j = 2
    for i = 2:length(z)-1
        zvec[j] = z[i]
        zvec[j+1] = z[i]
        j += 2
    end

    epsilonp = zeros(3, nz)
    for i = 1:3
        epsilonp[i, :] = epsilonbar[i] .+ kappa[i]*zvec
    end

    return epsilonbar, kappa, zvec, epsilonp
end

"""
stress for a thin laminate given stresses.
"""
function stresses(laminate, epsilonp)

    n = length(laminate)
    sigmap = zeros(3, 2*n)
    sigma = zeros(3, 2*n)
    epsilon = zeros(3, 2*n)

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
    epsilonbar, kappa, zvec, epsilonp = strains(alpha, beta, delta, z, forces)
    sigmap, sigma, epsilon = stresses(laminate, epsilonp)
    failure = tsai_wu_plane(sigma, laminate)

    return sigma, epsilon, failure
end



function tsai_hill(sigma, strength)
    S1t = strength.S1t
    S1c = strength.S1c
    S2t = strength.S2t
    S2c = strength.S2c
    S12 = strength.S12

    _, n = size(sigma)
    failure = zeros(n)  # fails if > 1
    for i = 1:n
        if sigma[1, i] >= 0.0
            S1 = S1t
        else
            S1 = S1c
        end
        if sigma[2, i] >= 0.0
            S2 = S2t
        else
            S2 = S2c
        end
        failure[i] = sigma[1, i]^2/S1^2 + sigma[2, i]^2/S2^2 + sigma[3, i]^2/S12^2 - sigma[1, i]*sigma[2, i]/S1^2
    end

    return failure
end

function tsai_wu_plane(sigma, laminate)
    # S1t = strength.S1t
    # S1c = strength.S1c
    # S2t = strength.S2t
    # S2c = strength.S2c
    # S12 = strength.S12

    # _, n = size(sigma)
    n = length(laminate)
    failure = zeros(2*n)  # failure if > 1
    i = 1
    for k = 1:n
        (; S1t, S1c, S2t, S2c, S12) = laminate[k].material
        for j = i:i+1
            failure[j] = sigma[1, j]^2/(S1t*S1c) + sigma[2, j]^2/(S2t*S2c) - sigma[1, j]*sigma[2, j]/sqrt(S1t*S1c*S2t*S2c) +
                sigma[1, j]*(1/S1t - 1/S1c) + sigma[2, j]*(1/S2t - 1/S2c) + sigma[3, j]^2/S12^2
        end
        i += 2
    end

    return failure
end


# ------------  CLT for Beams of Laminates -----------------

struct BeamSection{VL, VF} #Note: I suggest that we rename this to Region, or SectionRegion, or something like that. Beam section makes it sound like it is a section of a beam (a cross section), not a section of a cross-section.
    laminate::VL  #Vector{Lamina}
    y::VF  # Vector{Float}
    z::VF  # Vector{Float}
end


"""
    get_beam_sections(x, y, chord, twist, paxis, xbreak, weblocs, segments, web_segments)

Creates a vector of BeamSections from the information for the airfoil mesh. 

**Arguments**

"""
function get_beam_sections(x, y, chord, twist, paxis, xbreak, weblocs, segments, web_segments; fit=Akima)
    ### Check that inputs are good. 
    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end
    #Todo: Check that xbreak starts and ends with 0 and 1, respectively.

    if !isapprox(minimum(x), 0) || !isapprox(maximum(x), 1)
        throw(ArgumentError("x must start at 0 and end at 1"))
    end
    
    

    ns = length(x) - 1
    # n = length(xbreak) - 1 #Number of regions

    sections = Vector{BeamSection}(undef, ns) #Todo: Typing

    s, c = sincos(twist) #I think applied twist correctly. 
    xc = paxis * chord
    for i = 1:ns #Iterate over the af coordinates and create a sections
        xbar = (x[i] + x[i+1])/2 #Midpoint of the section
        idx = findfirst(x ->  x >= xbar, xbreak) - 1 #Find which region of the cross-section the section is. 

        y_i = @. (x[i:i+1]*chord - xc)*c + y[i:i+1]*s + xc
        z_i = @. -(x[i:i+1] - xc)*s + y[i:i+1].*chord*c
        # x = nodes[i].x
        # y = nodes[i].y
        # nodes[i] = Node((x - xc)*c + y*s + xc, -(x - xc)*s + y*c)
        sections[i] = BeamSection(segments[idx], y_i, z_i)
    end

    ### Add in the webs
    idx = argmin(x) #todo: This might not split the airfoils well. 
    xtop = reverse(x[1:idx])
    ytop = reverse(y[1:idx])
    xbot = x[idx+1:end]
    ybot = y[idx+1:end]
    topfit = fit(xtop, ytop)
    botfit = fit(xbot, ybot)

    websections = Vector{BeamSection}(undef, length(weblocs))
    for i in eachindex(weblocs)
        y_i = [weblocs[i], weblocs[i]].*chord
        ztop = topfit(weblocs[i])
        zbot = botfit(weblocs[i])
        z_i = [zbot, ztop].*chord
        y_i = @. (y_i - xc)*c + z_i*s + xc
        z_i = @. -(y_i - xc)*s + z_i*c
        websections[i] = BeamSection(web_segments[i], y_i, z_i)
    end

    return vcat(sections, websections)
end


struct CLT <: CompositeSectionAnalysis
    sections::Vector{BeamSection}  # a vector of beam sections
    closed_section::Bool
end

CLT(sections) = CLT(sections, true)  # default to closed section


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



function compliance_matrix(clt::CLT, shear_center=true)

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
            I += I1*(omega\Rk)  # repeated, could cache

            F1 = [alpha[3, 3] beta[3, 2];
                 beta[3, 2] delta[2, 2]]
            F += b*F1 - I1*(omega\I1')

        end

        # if i <= 10   #TODO: temporary hack
        #     A += abs(Asub)/2.0
        # end
        # A += Asub/2.0
    end
    L = -I
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
    sc = [0.0, 0.0]  # TODO
    tc = [yc, zc]

    # TODO: move to sc or tc.

    # return S, sc, tc #What does S look like before he nerfs it? -> 4x4 symmetric... could have looked above to figure that one out. 
    

    Sfull = zeros(6, 6)
    idx = [1, 5, 6, 4]
    for i = 1:4
        for j = i:4
            Sfull[idx[i], idx[j]] = S[i, j]
        end
    end

    s, S, ysc, zsc = shearflow_general_attempt(clt.sections, Wbar, F, L, yc, zc) 
    Sfull[2, 2] = s[1, 1]
    Sfull[3, 3] = s[2, 2]
    # Sfull[2, 2] = s[2, 2]
    # Sfull[3, 3] = s[1, 1]


    Sfull = Symmetric(Sfull)

    
    
    return Sfull, sc, tc
end

function mass_matrix(clt::CLT, shear_center=true) #Todo: Finish me!
    mu = 0.0 #Todo: Typing
    xm2 = 0.0
    xm3 = 0.0
    i22 = 0.0
    i33 = 0.0
    i23 = 0.0
    
    for i in eachindex(clt.sections)
        sec = clt.sections[i]
        x = abs(diff(sec.y)[1])
        y = abs(diff(sec.z)[1])
        L = sqrt(x^2 + y^2)
        
        for j in eachindex(sec.laminate)

            # A = 0.5*(zp[k-1] + zp[k]) * (yp[k-1] - yp[k])  # if closed section (trapezoid formula for polygon area: https://en.wikipedia.org/wiki/Shoelace_formula) -> Stolen from below. 
            A = L*sec.laminate[j].t

            mu += A*sec.laminate[j].material.rho
        end
    end


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

    qoy = 0.0
    qoz = 0.0
    qcy = 0.0
    qcz = 0.0
    eta = 0.0
    AMy = zeros(2, 2)
    bvy = zeros(2)
    AMz = zeros(2, 2)
    bvz = zeros(2)
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
            eta += b
            Rk = [1.0 zbar ybar 0.0;
                0.0 ca -sa 0.0;
                0.0 sa ca 0.0;
                0 0 0 1]
            Reta = [1.0 0 eta 0;
                    0 1 0 0;
                    0 0 0 -2]

            M = (muk\(Reta*Rk - nuk*(F\L)))*Wbar

            # --------
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
            eta += b
            Rk = [1.0 zbar ybar 0.0;
                0.0 ca -sa 0.0;
                0.0 sa ca 0.0;
                0 0 0 1]
            Reta = [1.0 0 eta 0;
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

