using GXBeamCS
using Test
using FLOWMath

function stiffness_in_internal_order(S)
    S2 = S[[1, 4, 5, 6], [1, 4, 5, 6]]  # remove zeros for shear flow
    K = inv(S2)
    K = K[[1, 3, 4, 2], [1, 3, 4, 2]]  # change from gxbeam order Faxial {Fy, Fz omitted}, torsion, Myy, Mzz to Faxial, Myy, Mzz, torsion
    return K, S2[[1, 3, 4, 2], [1, 3, 4, 2]]
end

# function stiffness_in_gxbeam_order(S)
#     S2 = S[[1, 4, 5, 6], [1, 4, 5, 6]]  # remove zeros for shear flow
#     K2 = inv(S2)
#     K = zeros(6, 6)
#     K[[1, 4, 5, 6], [1, 4, 5, 6]] .= K2

#     return K
# end

#------- Foundations of Classical Laminate Theory, Andreas Öchsner -----

function lam_test_material()

    E1 = 129000.0e6  # note must be typo on this book
    E2 = 11000.0e6
    G12 = 6600.0e6
    nu12 = 0.28
    rho = 1.0
    S1t = 1950.0e6
    S1c = 1480e6
    S2t = 48e6
    S2c = 200e6
    S12 = 79e6
    mat = MaterialPlane(E1, E2, G12, nu12, rho, S1t, S1c, S2t, S2c, S12)

    t = 8e-3

    return mat, t
end

@testset "Lamiante CLT Test 1" begin
# -------- 4.2. Problem 1 ------------

mat, t = lam_test_material()

theta = [45 -45 0 90 90 0 -45 45]*pi/180
laminate = Layer.(Ref(mat), t/8, theta)

Q1, _, _ = GXBeamCS.Qbar(laminate[1])
@test isapprox(Q1[1, 1]/1e6, 43385.924, atol=1e-3)
@test isapprox(Q1[1, 2]/1e6, 30185.924, atol=1e-3)
@test isapprox(Q1[1, 3]/1e6, 29698.543, atol=1e-3)
@test isapprox(Q1[2, 2]/1e6, 43385.924, atol=1e-3)
@test isapprox(Q1[2, 3]/1e6, 29698.543, atol=1e-3)
@test isapprox(Q1[3, 3]/1e6, 33685.195, atol=1e-3)
@test Q1[1, 2] == Q1[2, 1]
@test Q1[1, 3] == Q1[3, 1]
@test Q1[2, 3] == Q1[3, 2]

Q2, _, _ = GXBeamCS.Qbar(laminate[2])
@test isapprox(Q2[1, 1]/1e6, 43385.924, atol=1e-3)
@test isapprox(Q2[1, 2]/1e6, 30185.924, atol=1e-3)
@test isapprox(Q2[1, 3]/1e6, -29698.543, atol=1e-3)
@test isapprox(Q2[2, 2]/1e6, 43385.924, atol=1e-3)
@test isapprox(Q2[2, 3]/1e6, -29698.543, atol=1e-3)
@test isapprox(Q2[3, 3]/1e6, 33685.195, atol=1e-3)
@test Q2[1, 2] == Q2[2, 1]
@test Q2[1, 3] == Q2[3, 1]
@test Q2[2, 3] == Q2[3, 2]

Q3, _, _ = GXBeamCS.Qbar(laminate[3])
@test isapprox(Q3[1, 1]/1e6, 129868.204, atol=1e-3)
@test isapprox(Q3[1, 2]/1e6, 3100.729, atol=1e-3)
@test isapprox(Q3[1, 3]/1e6, 0.0, atol=1e-3)
@test isapprox(Q3[2, 2]/1e6, 11074.033, atol=1e-3)
@test isapprox(Q3[2, 3]/1e6, 0.0, atol=1e-3)
@test isapprox(Q3[3, 3]/1e6, 6600.0, atol=1e-3)

Q4, _, _ = GXBeamCS.Qbar(laminate[4])
@test isapprox(Q4[1, 1]/1e6, 11074.033, atol=1e-3)
@test isapprox(Q4[1, 2]/1e6, 3100.729, atol=1e-3)
@test isapprox(Q4[1, 3]/1e6, 0.0, atol=1e-3)
@test isapprox(Q4[2, 2]/1e6, 129868.204, atol=1e-3)
@test isapprox(Q4[2, 3]/1e6, 0.0, atol=1e-3)
@test isapprox(Q4[3, 3]/1e6, 6600.0, atol=1e-3)

z, h = GXBeamCS.zspacing(laminate)
A, B, D = GXBeamCS.laminatestiffnessmatrix(laminate, z)

@test isapprox(A[1, 1]/1e3, 455428.170, atol=1e-3)
@test isapprox(A[1, 2]/1e3, 133146.612, atol=1e-3)
@test isapprox(A[1, 3]/1e3, 0.0, atol=1e-3)
@test isapprox(A[2, 2]/1e3, 455428.170, atol=1e-3)
@test isapprox(A[2, 3]/1e3, 0.0, atol=1e-3)
@test isapprox(A[3, 3]/1e3, 161140.779, atol=1e-3)

@test isapprox(B[1, 1], 0.0, atol=1e-10)
@test isapprox(B[1, 2], 0.0, atol=1e-10)
@test isapprox(B[1, 3], 0.0, atol=1e-10)
@test isapprox(B[2, 2], 0.0, atol=1e-10)
@test isapprox(B[2, 3], 0.0, atol=1e-10)
@test isapprox(B[3, 3], 0.0, atol=1e-10)

@test isapprox(D[1, 1]*1e3, 2233175.466, atol=1e-3)
@test isapprox(D[1, 2]*1e3, 1143478.381, atol=1e-3)
@test isapprox(D[1, 3]*1e3, 356382.514, atol=1e-3)
@test isapprox(D[2, 2]*1e3, 1757998.781, atol=1e-3)
@test isapprox(D[2, 3]*1e3, 356382.514, atol=1e-3)
@test isapprox(D[3, 3]*1e3, 1292780.601, atol=1e-3)

alpha, beta, delta = GXBeamCS.laminatecompliancematrix(A, B, D)

@test isapprox(alpha[1, 1]/1e-10, 24.009482, atol=1e-6)
@test isapprox(alpha[1, 2]/1e-10, -7.019287, atol=1e-6)
@test isapprox(alpha[1, 3]/1e-10, 0.0, atol=1e-6)
@test isapprox(alpha[2, 2]/1e-10, 24.009482, atol=1e-6)
@test isapprox(alpha[2, 3]/1e-10, 0.0, atol=1e-6)
@test isapprox(alpha[3, 3]/1e-10, 62.057538, atol=1e-6)

@test isapprox(beta[1, 1], 0.0, atol=1e-10)
@test isapprox(beta[1, 2], 0.0, atol=1e-10)
@test isapprox(beta[1, 3], 0.0, atol=1e-10)
@test isapprox(beta[2, 2], 0.0, atol=1e-10)
@test isapprox(beta[2, 3], 0.0, atol=1e-10)
@test isapprox(beta[3, 3], 0.0, atol=1e-10)

@test isapprox(delta[1, 1]/1e-4, 6.771890, atol=1e-6)
@test isapprox(delta[1, 2]/1e-4, -4.264613, atol=1e-6)
@test isapprox(delta[1, 3]/1e-4, -0.691184, atol=1e-6)
@test isapprox(delta[2, 2]/1e-4, 8.710637, atol=1e-6)
@test isapprox(delta[2, 3]/1e-4, -1.225641, atol=1e-6)
@test isapprox(delta[3, 3]/1e-4, 8.263678, atol=1e-6)

Nx = 1000e3
forces = [Nx; 0.0; 0.0; 0.0; 0.0; 0.0]

epsilonbar, kappa, zvec, epsilonp = GXBeamCS.strains(alpha, beta, delta, z, forces)

@test isapprox(epsilonbar[1]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilonbar[2]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilonbar[3]*1e3, 0.0, atol=1e-10)
@test isapprox(kappa[1], 0.0, atol=1e-10)
@test isapprox(kappa[2], 0.0, atol=1e-10)
@test isapprox(kappa[3], 0.0, atol=1e-10)


sigmap, sigma, epsilon = GXBeamCS.stresses(laminate, epsilonp)

# using PyPlot
# close("all"); pygui(true)

# figure()
# plot(sigmap[1, :]/1e6, zvec)
# figure()
# plot(epsilonp[1, :]/1e6, zvec)

@test isapprox(sigma[1, 1]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 2]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 3]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 4]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 5]/1e6, 309.630, atol=1e-3)
@test isapprox(sigma[1, 6]/1e6, 309.630, atol=1e-3)
@test isapprox(sigma[1, 7]/1e6, -83.714, atol=1e-3)
@test isapprox(sigma[1, 8]/1e6, -83.714, atol=1e-3)
@test isapprox(sigma[1, 9]/1e6, -83.714, atol=1e-3)
@test isapprox(sigma[1, 10]/1e6, -83.714, atol=1e-3)
@test isapprox(sigma[1, 11]/1e6, 309.630, atol=1e-3)
@test isapprox(sigma[1, 12]/1e6, 309.630, atol=1e-3)
@test isapprox(sigma[1, 13]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 14]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 15]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 16]/1e6, 112.958, atol=1e-3)

@test isapprox(sigma[2, 1]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 2]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 3]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 4]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 5]/1e6, -0.328, atol=1e-3)
@test isapprox(sigma[2, 6]/1e6, -0.328, atol=1e-3)
@test isapprox(sigma[2, 7]/1e6, 24.412, atol=1e-3)
@test isapprox(sigma[2, 8]/1e6, 24.412, atol=1e-3)
@test isapprox(sigma[2, 9]/1e6, 24.412, atol=1e-3)
@test isapprox(sigma[2, 10]/1e6, 24.412, atol=1e-3)
@test isapprox(sigma[2, 11]/1e6, -0.328, atol=1e-3)
@test isapprox(sigma[2, 12]/1e6, -0.328, atol=1e-3)
@test isapprox(sigma[2, 13]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 14]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 15]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 16]/1e6, 12.042, atol=1e-3)

@test isapprox(sigma[3, 1]/1e6, -20.479, atol=1e-3)
@test isapprox(sigma[3, 2]/1e6, -20.479, atol=1e-3)
@test isapprox(sigma[3, 3]/1e6, 20.479, atol=1e-3)
@test isapprox(sigma[3, 4]/1e6, 20.479, atol=1e-3)
@test isapprox(sigma[3, 5]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 6]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 7]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 8]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 9]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 10]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 11]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 12]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 13]/1e6, 20.479, atol=1e-3)
@test isapprox(sigma[3, 14]/1e6, 20.479, atol=1e-3)
@test isapprox(sigma[3, 15]/1e6, -20.479, atol=1e-3)
@test isapprox(sigma[3, 16]/1e6, -20.479, atol=1e-3)

@test isapprox(epsilon[1, 1]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 2]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 3]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 4]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 5]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[1, 6]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[1, 7]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[1, 8]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[1, 9]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[1, 10]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[1, 11]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[1, 12]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[1, 13]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 14]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 15]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 16]*1e3, 0.849510, atol=1e-6)

@test isapprox(epsilon[2, 1]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 2]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 3]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 4]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 5]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[2, 6]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[2, 7]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[2, 8]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[2, 9]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[2, 10]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[2, 11]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[2, 12]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[2, 13]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 14]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 15]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 16]*1e3, 0.849510, atol=1e-6)

@test isapprox(epsilon[3, 1]*1e3, -3.102877, atol=1e-6)
@test isapprox(epsilon[3, 2]*1e3, -3.102877, atol=1e-6)
@test isapprox(epsilon[3, 3]*1e3, 3.102877, atol=1e-6)
@test isapprox(epsilon[3, 4]*1e3, 3.102877, atol=1e-6)
@test isapprox(epsilon[3, 5]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 6]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 7]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 8]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 9]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 10]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 11]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 12]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 13]*1e3, 3.102877, atol=1e-6)
@test isapprox(epsilon[3, 14]*1e3, 3.102877, atol=1e-6)
@test isapprox(epsilon[3, 15]*1e3, -3.102877, atol=1e-6)
@test isapprox(epsilon[3, 16]*1e3, -3.102877, atol=1e-6)

forces = [0.0; 0.0; 0.0; 1000; 0.0; 0.0]

epsilonbar, kappa, zvec, epsilonp = GXBeamCS.strains(alpha, beta, delta, z, forces)

@test isapprox(epsilonbar[1], 0.0, atol=1e-6)
@test isapprox(epsilonbar[2], 0.0, atol=1e-6)
@test isapprox(epsilonbar[3], 0.0, atol=1e-10)
@test isapprox(kappa[1], 0.677189, atol=1e-6)
@test isapprox(kappa[2], -0.426461, atol=1e-6)
@test isapprox(kappa[3], -0.069118, atol=1e-6)

sigmap, sigma, epsilon = GXBeamCS.stresses(laminate, epsilonp)


# figure()
# plot(sigmap[1, :]/1e6, zvec*1e3)
# figure()
# plot(epsilonp[1, :]*1e3, zvec*1e3)

@test isapprox(sigma[1, 1]/1e6, -49.154, atol=1e-3)
@test isapprox(sigma[1, 2]/1e6, -36.866, atol=1e-3)
@test isapprox(sigma[1, 3]/1e6, -63.151, atol=1e-3)
@test isapprox(sigma[1, 4]/1e6, -42.101, atol=1e-3)
@test isapprox(sigma[1, 5]/1e6, -173.246, atol=1e-3)
@test isapprox(sigma[1, 6]/1e6, -86.623, atol=1e-3)
@test isapprox(sigma[1, 7]/1e6, 53.284, atol=1e-3)
@test isapprox(sigma[1, 8]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[1, 9]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[1, 10]/1e6, -53.284, atol=1e-3)
@test isapprox(sigma[1, 11]/1e6, 86.623, atol=1e-3)
@test isapprox(sigma[1, 12]/1e6, 173.246, atol=1e-3)
@test isapprox(sigma[1, 13]/1e6, 42.101, atol=1e-3)
@test isapprox(sigma[1, 14]/1e6, 63.151, atol=1e-3)
@test isapprox(sigma[1, 15]/1e6, 36.866, atol=1e-3)
@test isapprox(sigma[1, 16]/1e6, 49.154, atol=1e-3)

@test isapprox(sigma[2, 1]/1e6, -8.210, atol=1e-3)
@test isapprox(sigma[2, 2]/1e6, -6.158, atol=1e-3)
@test isapprox(sigma[2, 3]/1e6, -4.504, atol=1e-3)
@test isapprox(sigma[2, 4]/1e6, -3.003, atol=1e-3)
@test isapprox(sigma[2, 5]/1e6, 5.246, atol=1e-3)
@test isapprox(sigma[2, 6]/1e6, 2.623, atol=1e-3)
@test isapprox(sigma[2, 7]/1e6, -6.177, atol=1e-3)
@test isapprox(sigma[2, 8]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[2, 9]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[2, 10]/1e6, 6.177, atol=1e-3)
@test isapprox(sigma[2, 11]/1e6, -2.623, atol=1e-3)
@test isapprox(sigma[2, 12]/1e6, -5.246, atol=1e-3)
@test isapprox(sigma[2, 13]/1e6, 3.003, atol=1e-3)
@test isapprox(sigma[2, 14]/1e6, 4.504, atol=1e-3)
@test isapprox(sigma[2, 15]/1e6, 6.158, atol=1e-3)
@test isapprox(sigma[2, 16]/1e6, 8.210, atol=1e-3)

@test isapprox(sigma[3, 1]/1e6, 29.136, atol=1e-3)
@test isapprox(sigma[3, 2]/1e6, 21.852, atol=1e-3)
@test isapprox(sigma[3, 3]/1e6, -21.852, atol=1e-3)
@test isapprox(sigma[3, 4]/1e6, -14.568, atol=1e-3)
@test isapprox(sigma[3, 5]/1e6, 0.912, atol=1e-3)
@test isapprox(sigma[3, 6]/1e6, 0.456, atol=1e-3)
@test isapprox(sigma[3, 7]/1e6, -0.456, atol=1e-3)
@test isapprox(sigma[3, 8]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 9]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 10]/1e6, 0.456, atol=1e-3)
@test isapprox(sigma[3, 11]/1e6, -0.456, atol=1e-3)
@test isapprox(sigma[3, 12]/1e6, -0.912, atol=1e-3)
@test isapprox(sigma[3, 13]/1e6, 14.568, atol=1e-3)
@test isapprox(sigma[3, 14]/1e6, 21.852, atol=1e-3)
@test isapprox(sigma[3, 15]/1e6, -21.852, atol=1e-3)
@test isapprox(sigma[3, 16]/1e6, -29.136, atol=1e-3)
end

# ----- 4.3. Problem 2 --------
@testset "Lamiante CLT Test 2" begin

mat, t = lam_test_material()
theta = [45 -45 0 90 0 90 45 -45]*pi/180
laminate = Layer.(Ref(mat), t/8, theta)

forces = [1000*1e3; 0.0; 0.0; 0.0; 0.0; 0.0]

z, h = GXBeamCS.zspacing(laminate)
A, B, D = GXBeamCS.laminatestiffnessmatrix(laminate, z)
alpha, beta, delta = GXBeamCS.laminatecompliancematrix(A, B, D)
epsilonbar, kappa, zvec, epsilonp = GXBeamCS.strains(alpha, beta, delta, z, forces)
sigmap, sigma, epsilon = GXBeamCS.stresses(laminate, epsilonp)

@test isapprox(A[1, 1]/1e3, 455428.170, atol=1e-3)
@test isapprox(A[1, 2]/1e3, 133146.612, atol=1e-3)
@test isapprox(A[1, 3]/1e3, 0.0, atol=1e-3)
@test isapprox(A[2, 2]/1e3, 455428.170, atol=1e-3)
@test isapprox(A[2, 3]/1e3, 0.0, atol=1e-3)
@test isapprox(A[3, 3]/1e3, 161140.779, atol=1e-3)

@test isapprox(B[1, 1], -118794.171, atol=1e-3)
@test isapprox(B[1, 2], 1.455e-11, atol=1e-3)
@test isapprox(B[1, 3], -59397.086, atol=1e-3)
@test isapprox(B[2, 2], 118794.171, atol=1e-3)
@test isapprox(B[2, 3], -59397.086, atol=1e-3)
@test isapprox(B[3, 3], 0.0, atol=1e-10)

@test isapprox(D[1, 1]*1e3, 1995587.124, atol=1e-3)
@test isapprox(D[1, 2]*1e3, 1143478.381, atol=1e-3)
@test isapprox(D[1, 3]*1e3, 0.0, atol=1e-3)
@test isapprox(D[2, 2]*1e3, 1995587.124, atol=1e-3)
@test isapprox(D[2, 3]*1e3, 0.0, atol=1e-3)
@test isapprox(D[3, 3]*1e3, 1292780.601, atol=1e-3)

@test isapprox(alpha[1, 1]/1e-10, 24.562273, atol=1e-6)
@test isapprox(alpha[1, 2]/1e-10, -6.911749, atol=1e-6)
@test isapprox(alpha[1, 3]/1e-10, 0.445254, atol=1e-6)
@test isapprox(alpha[2, 2]/1e-10, 24.562273, atol=1e-6)
@test isapprox(alpha[2, 3]/1e-10, -0.445254, atol=1e-6)
@test isapprox(alpha[3, 3]/1e-10, 62.948045, atol=1e-6)

@test isapprox(beta[1, 1]/1e-7, 1.834321, atol=1e-6)
@test isapprox(beta[1, 2]/1e-7, -0.626374, atol=1e-6)
@test isapprox(beta[1, 3]/1e-7, 0.810957, atol=1e-6)
@test isapprox(beta[2, 2]/1e-7, -1.834321, atol=1e-6)
@test isapprox(beta[2, 3]/1e-7, 0.810957, atol=1e-6)
@test isapprox(beta[3, 3]/1e-7, 0.0, atol=1e-6)

@test isapprox(delta[1, 1]/1e-4, 7.677865, atol=1e-6)
@test isapprox(delta[1, 2]/1e-4, -4.400777, atol=1e-6)
@test isapprox(delta[1, 3]/1e-4, 0.113057, atol=1e-6)
@test isapprox(delta[2, 2]/1e-4, 7.677865, atol=1e-6)
@test isapprox(delta[2, 3]/1e-4, -0.113057, atol=1e-6)
@test isapprox(delta[3, 3]/1e-4, 7.809784, atol=1e-6)

@test isapprox(epsilonbar[1]*1e3, 2.456227, atol=1e-6)
@test isapprox(epsilonbar[2]*1e3, -0.691175, atol=1e-6)
@test isapprox(epsilonbar[3]*1e3, 0.044525, atol=1e-6)
@test isapprox(kappa[1], 0.183432, atol=1e-6)
@test isapprox(kappa[2], -0.062637, atol=1e-6)
@test isapprox(kappa[3], 0.081096, atol=1e-6)

@test isapprox(sigma[1, 1]/1e6, 67.486, atol=1e-3)
@test isapprox(sigma[1, 2]/1e6, 80.657, atol=1e-3)
@test isapprox(sigma[1, 3]/1e6, 105.854, atol=1e-3)
@test isapprox(sigma[1, 4]/1e6, 108.745, atol=1e-3)
@test isapprox(sigma[1, 5]/1e6, 269.587, atol=1e-3)
@test isapprox(sigma[1, 6]/1e6, 293.215, atol=1e-3)
@test isapprox(sigma[1, 7]/1e6, -74.580, atol=1e-3)
@test isapprox(sigma[1, 8]/1e6, -82.146, atol=1e-3)
@test isapprox(sigma[1, 9]/1e6, 316.843, atol=1e-3)
@test isapprox(sigma[1, 10]/1e6, 340.470, atol=1e-3)
@test isapprox(sigma[1, 11]/1e6, -89.711, atol=1e-3)
@test isapprox(sigma[1, 12]/1e6, -97.277, atol=1e-3)
@test isapprox(sigma[1, 13]/1e6, 146.513, atol=1e-3)
@test isapprox(sigma[1, 14]/1e6, 159.684, atol=1e-3)
@test isapprox(sigma[1, 15]/1e6, 123.199, atol=1e-3)
@test isapprox(sigma[1, 16]/1e6, 126.090, atol=1e-3)

@test isapprox(sigma[2, 1]/1e6, 10.201, atol=1e-3)
@test isapprox(sigma[2, 2]/1e6, 10.734, atol=1e-3)
@test isapprox(sigma[2, 3]/1e6, 9.149, atol=1e-3)
@test isapprox(sigma[2, 4]/1e6, 10.328, atol=1e-3)
@test isapprox(sigma[2, 5]/1e6, 0.212, atol=1e-3)
@test isapprox(sigma[2, 6]/1e6, 0.087, atol=1e-3)
@test isapprox(sigma[2, 7]/1e6, 23.220, atol=1e-3)
@test isapprox(sigma[2, 8]/1e6, 25.057, atol=1e-3)
@test isapprox(sigma[2, 9]/1e6, -0.038, atol=1e-3)
@test isapprox(sigma[2, 10]/1e6, -0.163, atol=1e-3)
@test isapprox(sigma[2, 11]/1e6, 26.894, atol=1e-3)
@test isapprox(sigma[2, 12]/1e6, 28.731, atol=1e-3)
@test isapprox(sigma[2, 13]/1e6, 13.398, atol=1e-3)
@test isapprox(sigma[2, 14]/1e6, 13.931, atol=1e-3)
@test isapprox(sigma[2, 15]/1e6, 16.225, atol=1e-3)
@test isapprox(sigma[2, 16]/1e6, 17.405, atol=1e-3)
end


# ------- 4.4. Problem 3 --------

@testset "Lamiante CLT Test 3" begin

mat, t = lam_test_material()
theta = [45 0 30 -45]*pi/180
laminate = Layer.(Ref(mat), t/4, theta)

forces = [1000*1e3; 500e3; 0.0; 0.0; 0.0; 0.0]

# S1t = 1950.0e6
# S1c = 1480.0e6
# S2t = 48.0e6
# S2c = 200.0e6
# S12 = 79.0e6
# strength = CompositeStrength(S1t, S1c, S2t, S2c, S12)

sigma, epsilon, wu = GXBeamCS.clt(laminate, forces)

@test isapprox(sigma[1, 1]/1e6, 18.776, atol=1e-3)
# @test isapprox(sigma[1, 2]/1e6, 18.776, atol=1e-3)  # typo in text, repeated value
@test isapprox(sigma[1, 3]/1e6, 238.813, atol=1e-3)
@test isapprox(sigma[1, 4]/1e6, 229.609, atol=1e-3)
@test isapprox(sigma[1, 5]/1e6, 149.134, atol=1e-3)
@test isapprox(sigma[1, 6]/1e6, 229.851, atol=1e-3)
@test isapprox(sigma[1, 7]/1e6, 214.113, atol=1e-3)
@test isapprox(sigma[1, 8]/1e6, 13.857, atol=1e-3)

@test isapprox(epsilon[3, 1]*1e3, 2.663312, atol=1e-6)
@test isapprox(epsilon[3, 2]*1e3, 1.790924, atol=1e-6)
@test isapprox(epsilon[3, 3]*1e3, -4.138210, atol=1e-6)
@test isapprox(epsilon[3, 4]*1e3, -1.996388, atol=1e-6)
@test isapprox(epsilon[3, 5]*1e3, -0.202719, atol=1e-6)
@test isapprox(epsilon[3, 6]*1e3, 0.112682, atol=1e-6)
@test isapprox(epsilon[3, 7]*1e3, -0.046148, atol=1e-6)
@test isapprox(epsilon[3, 8]*1e3, 0.826241, atol=1e-6)


# failure = tsai_hill(sigma, strength)
# R = sqrt.(1.0 ./ failure)

# @test isapprox(R[1], 0.683, atol=1e-3)
# @test isapprox(R[2], 0.880, atol=1e-3)
# @test isapprox(R[3], 1.001, atol=1e-3)
# @test isapprox(R[4], 1.346, atol=1e-3)
# @test isapprox(R[5], 1.214, atol=1e-3)
# @test isapprox(R[6], 1.999, atol=1e-3)
# @test isapprox(R[7], 1.929, atol=1e-3)
# @test isapprox(R[8], 1.828, atol=1e-3)

R = [0.681; 0.901; 1.041; 1.467; 1.296; 2.325; 2.216; 1.836]
sigmaR = copy(sigma)
sigmaR[1, :] .*= R
sigmaR[2, :] .*= R
sigmaR[3, :] .*= R
failure = GXBeamCS.tsai_wu_plane(sigmaR, laminate)

@test isapprox(failure[1], 1.0, atol=1e-2)
@test isapprox(failure[2], 1.0, atol=1e-2)
@test isapprox(failure[3], 1.0, atol=1e-2)
@test isapprox(failure[4], 1.0, atol=1e-2)
@test isapprox(failure[5], 1.0, atol=1e-2)
@test isapprox(failure[6], 1.0, atol=1e-2)
@test isapprox(failure[7], 1.0, atol=1e-2)
@test isapprox(failure[8], 1.0, atol=1e-2)

forces = [0.0; 0.0; 0.0; 1000; 500; 0.0]

sigma, epsilon, wu = GXBeamCS.clt(laminate, forces)

# failure = tsai_hill(sigma, strength)
# R = sqrt.(1.0 ./ failure)

# @test isapprox(R[1], 5.324, atol=1e-2)
# @test isapprox(R[2], 9.990, atol=3e-3)
# # @test isapprox(R[3], 5.973, atol=1e-3)  # typo?  all stresses and strains match
# @test isapprox(R[4], 10.557, atol=1e-3)
# @test isapprox(R[5], 20.539, atol=1e-3)
# @test isapprox(R[6], 9.977, atol=1e-3)
# @test isapprox(R[7], 4.698, atol=1e-3)
# @test isapprox(R[8], 2.523, atol=1e-3)

R = [6.727; 10.909; 9.786; 14.021; 19.030; 10.114; 4.463; 2.602]
sigmaR = copy(sigma)
sigmaR[1, :] .*= R
sigmaR[2, :] .*= R
sigmaR[3, :] .*= R
failure = GXBeamCS.tsai_wu_plane(sigmaR, laminate)

@test isapprox(failure[1], 1.0, atol=1e-2)
@test isapprox(failure[2], 1.0, atol=1e-2)
# @test isapprox(failure[3], 1.0, atol=1e-2)  # again there appears an error in the text for this station.  stresses checkout.
@test isapprox(failure[4], 1.0, atol=1e-2)
@test isapprox(failure[5], 1.0, atol=1e-2)
@test isapprox(failure[6], 1.0, atol=1e-2)
@test isapprox(failure[7], 1.0, atol=1e-2)
@test isapprox(failure[8], 1.0, atol=1e-2)
end

# ------- Kollar example 6.2 ------

function beam_test_material()
    E1 = 148e9
    E2 = 9.65e9
    G12 = 4.55e9
    nu12 = 0.3
    rho = 1.0
    m0 = MaterialPlane(E1, E2, G12, nu12, rho)
    t0 = 0.1e-3

    E1 = 16.39e9
    E2 = 16.39e9
    G12 = 38.19e9
    nu12 = 0.801
    rho = 1.0
    m45 = MaterialPlane(E1, E2, G12, nu12, rho)
    t45 = 0.2e-3

    return m0, t0, m45, t45
end

@testset "Beam CLT Test 6.2" begin

m0, t0, m45, t45 = beam_test_material()


t = [t45*ones(2); t0*ones(12); t45*ones(2)]
theta = [0.0; 0.0; zeros(12); 0.0; 0.0]*pi/180
mat = [fill(m45, 2); fill(m0, 12); fill(m45, 2)]
laminate = Layer.(mat, t, theta)

# alpha, beta, delta = laminatecompliance(laminate)

# b1 = BeamSection(laminate, [0.0; 50e-3], [61e-3; 61e-3])
# b2 = BeamSection(laminate, [25e-3; 25e-3], [0; 60e-3])
b1 = BeamSection(laminate, [0.0; 50e-3], [0.0; 0])
b2 = BeamSection(laminate, [25e-3; 25e-3], [-1e-3; -61e-3])
sections = [b1; b2]

yc, zc = GXBeamCS.centroid(sections)

@test isapprox(yc, 25e-3, atol=1e-4)
@test isapprox(61e-3+zc, 0.0441, atol=1e-4)

EA, EIyy, EIzz, EIyz, GJ = GXBeamCS.beamstiffnessold(sections)
closed_section = false
clt = CLT(sections, closed_section)
shear_center = false
S, _, _ = compliance_matrix(clt, shear_center)
K, _ = stiffness_in_internal_order(S)

@test isapprox(EA/1e6, 21.22, atol=1e-2)  # same for open/closed if B = 0
@test isapprox(EIyy/1e3, 8.530, atol=1e-3)

@test isapprox(K[1, 1]/1e6, 21.22, atol=1e-2)  # same for open/closed if B = 0
@test isapprox(K[2, 2]/1e3, 8.530, atol=1e-2)
end


# ----------- example 6.3 ------------
@testset "Beam CLT Test 6.3" begin

m0, t0, m45, t45 = beam_test_material()


t = [t45*ones(2); t0*ones(12); t45*ones(2)]
theta = [0.0; 0.0; zeros(12); 0.0; 0.0]*pi/180
mat = [fill(m45, 2); fill(m0, 12); fill(m45, 2)]
laminate = Layer.(mat, t, theta)

b1 = BeamSection(laminate, [0.0; 52e-3], [0.0; 0.0])
b2 = BeamSection(laminate, [51e-3; 51e-3], [1e-3; 69e-3])
b3 = BeamSection(laminate, [52e-3; 0.0], [70e-3; 70e-3])
b4 = BeamSection(laminate, [1e-3; 1e-3], [69e-3; 1e-3])
sections = [b1; b2; b3; b4]
EA, EIyy, EIzz, EIyz, GJ = GXBeamCS.beamstiffnessold(sections)

@test isapprox(EIyy/1e3, 34.692, atol=1e-3)
@test isapprox(EIzz/1e3, 20.924, atol=1e-3)

# the above gives exact with book for EI because we defined geometry consistent with the formulas, but is approximate with GJ.  to match GJ exactly need closed path below.
b = BeamSection(laminate, [0.0; 50e-3; 50e-3; 0.0; 0.0], [0.0; 0.0; 70e-3; 70e-3; 0.0])
EA, EIyy, EIzz, EIyz, GJ = GXBeamCS.beamstiffnessold([b])
@test isapprox(EIyy/1e3, 34.692, atol=1e-1)  # looser tolerances for these - exact above.
@test isapprox(EIzz/1e3, 20.924, atol=1e-1)  # looser tolerances for these
@test isapprox(GJ/1e3, 7.352, atol=1e-3)

# S, K, _ = beamstiffness(sections)
clt = CLT(sections)
shear_center = false
S, _, _ = compliance_matrix(clt, shear_center)
K, _ = stiffness_in_internal_order(S)

@test isapprox(K[1, 1]/1e6, 46.303, atol=1e-3)
@test isapprox(K[2, 2]/1e3, 34.692, atol=1e-3)
@test isapprox(K[3, 3]/1e3, 20.924, atol=1e-3)
@test isapprox(K[1, 2]/1e3, 0.0, atol=1e-6)
@test isapprox(K[1, 3]/1e3, 0.0, atol=1e-6)
@test isapprox(K[1, 4]/1e3, 0.0, atol=1e-6)
@test isapprox(K[2, 1]/1e3, 0.0, atol=1e-6)
@test isapprox(K[2, 3]/1e3, 0.0, atol=1e-6)
@test isapprox(K[2, 4]/1e3, 0.0, atol=1e-6)
@test isapprox(K[3, 4]/1e3, 0.0, atol=1e-6)

# S, K, _ = beamstiffness([b])
clt = CLT([b])
S, _, _ = compliance_matrix(clt, shear_center)
K, _ = stiffness_in_internal_order(S)

@test isapprox(K[4, 4]/1e3, 7.3, atol=1e-1)  # looser tolerance - Area is a bit ambiguous
end

# ---------- example 6.5 ------------------

@testset "Beam CLT Test 6.5" begin

m0, t0, m45, t45 = beam_test_material()

t = [t0*ones(10); t0*ones(10)]
theta = [zeros(10); 45*ones(10)]*pi/180
mat = [fill(m0, 10); fill(m0, 10)]
laminate = Layer.(mat, t, theta)

b = BeamSection(laminate, [0.0; 50e-3; 50e-3; 0.0; 0.0], [0.0; 0.0; 70e-3; 70e-3; 0.0])
# b = BeamSection(laminate, [-25; 25; 25; -25; -25]*1e-3, [-35; -35; 35; 35; -35]*1e-3)

# S, K = beamstiffness([b])
clt = CLT([b])
shear_center = false
S, _, _ = compliance_matrix(clt, shear_center)
_, S = stiffness_in_internal_order(S)
@test isapprox(S[1, 1]/1e-6, 0.02576, atol=1e-5)
@test isapprox(S[1, 4]/1e-6, -0.4237, atol=1e-4)
@test isapprox(S[2, 2]/1e-6, 33.91, atol=1e-2)
@test isapprox(S[3, 3]/1e-6, 55.63, atol=1e-2)
@test isapprox(S[4, 4]/1e-6, 250.56, atol=1e-2)

# TODO: add checks for yc/zc - although these are explicitly checked by EA


# # (just an internal check on shear flow)
# # -------- example 6.2 -------
# t = [t45*ones(2); t0*ones(12); t45*ones(2)]
# theta = [0; 0; zeros(12); 0; 0]*pi/180
# mat = [fill(m45, 2); fill(m0, 12); fill(m45, 2)]
# laminate = Lamina.(t, theta, mat)

# b1 = BeamSection(laminate, [25e-3; 25e-3], [-60e-3; 0.0])
# b2 = BeamSection(laminate, [0.0; 50e-3], [0.0; 0.0])
# sections = [b1; b2]
# W, P, yc, zc = beamstiffness(sections)
# s, S, ysc, zsc = shearflow(sections, P, yc, zc)
end


# ----- example 6.6 / Appendix A1 and A.8 ------

@testset "Beam shear flow test 6.6 and A7" begin

m0, t0, m45, t45 = beam_test_material()

t = [t45*ones(2); t0*ones(12); t45*ones(2)]
theta = [0; 0; zeros(12); 0; 0]*pi/180
mat = [fill(m45, 2); fill(m0, 12); fill(m45, 2)]
laminate = Layer.(mat, t, theta)

df = 49e-3
d = 62e-3
b = BeamSection(laminate, [0.0; df; df; 0.0], [0.0; 0.0; d; d])
sections = [b]
# W, P, yc, zc = beamstiffness(sections, closedsection=false)
closed_section = false
clt = CLT(sections, closed_section)
shear_center = false
W, sc, tc = compliance_matrix(clt, shear_center)
K, _ = stiffness_in_internal_order(W)

@test isapprox(K[1, 1]/1e6, 30.87, atol=0.01)
@test isapprox(K[2, 2]/1e3, 22.015, rtol=0.001)
@test isapprox(K[3, 3]/1e3, 8.188, rtol=0.001)
@test isapprox(K[4, 4], 13.19, atol=0.01)

yc = tc[1]; zc = tc[2]
s, S, ysc, zsc = GXBeamCS.shearflow(sections, K, yc, zc, closedsection=false)
# S, K, yc, zc, Wbar, F, L = beamstiffness(sections, closedsection=false)
# s2, S2, ysc2, zsc2 = GXBeamCS.shearflow_general_attempt(sections, clt.cache.Wbar, clt.cache.F, clt.cache.L, yc, zc)

a11 = 5.18e-9
bf = 50e-3
e = 3*bf^2/a11 / ( (6*bf + d)/a11 )
ysca = e + df - yc  # TODO: shear center does not seem to agree.  book may be using the approximate formula though...

alpha, beta, delta = GXBeamCS.laminatecompliance(laminate)
alpha_nu = alpha[3, 3] - beta[3, 3]^2/delta[3, 3]
deltay = (df - yc)/yc
rhoy = 3/5*(8 - 9*deltay + 3*deltay^2)/(2 - deltay)^2
qyu = 1/(2*df)*(3 - 3*deltay)/(2 - deltay)
gamma = 1 + 1/6*d/df
syy = rhoy/2*alpha_nu/df + d/3*qyu^2*alpha_nu
szz = alpha_nu/d + 2/3*alpha_nu*df/(d*gamma)^2
syz = 0.0

@test isapprox(syy/1e-7, s[1, 1]/1e-7, atol=0.01)
@test isapprox(szz/1e-7, s[2, 2]/1e-7, atol=0.1)
@test isapprox(syz/1e-7, s[1, 2]/1e-7, atol=1e-6)



# ---- appendix A.7


t = [t45*ones(2); t0*ones(12); t45*ones(2)]
theta = [0; 0; zeros(12); 0; 0]*pi/180
mat = [fill(m45, 2); fill(m0, 12); fill(m45, 2)]
laminate = Layer.(mat, t, theta)

df = 2.0
d = 4.0
b = BeamSection(laminate, [0.0; df; df; 0.0; 0.0], [0.0; 0.0; d; d; 0.0])
# W, P, yc, zc = beamstiffness([b])
closed_section = false
clt = CLT([b], closed_section)
shear_center = false
W, _, tc = compliance_matrix(clt, shear_center)
P, _ = stiffness_in_internal_order(W)
yc = tc[1]; zc = tc[2]
s, S, ysc, zsc = GXBeamCS.shearflow([b], P, yc, zc)

gammaz = 1.0 + 1/3*d/df
gammay = 1.0 + 1/3*df/d
szz = alpha_nu/(2*d) + alpha_nu*df/(6*d^2*gammaz^2)
syy = alpha_nu/(2*df) + alpha_nu*d/(6*df^2*gammay^2)
syz = 0.0

@test isapprox(syy/1e-8, s[1, 1]/1e-8, atol=0.01)
@test isapprox(szz/1e-8, s[2, 2]/1e-8, atol=0.05)
@test isapprox(syz/1e-8, s[1, 2]/1e-8, atol=1e-6)


end


# Loss of Accuracy Using Smeared Properties in Composite Beam Modeling
# Ning Liu, Purdue University

@testset "CLT stress. rectangular c.s" begin

# ----- rectangular cross-section ------
E1 = 41.5e9
E2 = 7.83e9
nu12 = 0.3
G12 = 3.15e9
rho = 1.0
mat1 = MaterialPlane(E1, E2, G12, nu12, rho)

t = fill(0.25, 8)  # they said mm in thesis, but units only work out if using meters
theta = [25.0, 25.0, 50.0, 0, 50, 0, 25, 25]*pi/180
laminate1 = Layer.(Ref(mat1), t, theta)

b = BeamSection(laminate1, [-2.0, 2], [0.0, 0.0])

closed_section = false
clt = CLT([b], closed_section)

shear_center = true
S, _, tc = compliance_matrix(clt, shear_center)
K = clt_stiffness_matrix(S)

# F = [0.0; 0; 0]
# M = [-1000.0; 0; 0]
# # reorder to GXBeam order
# F = F[[3, 1, 2]]
# M = M[[3, 1, 2]]
F = zeros(3)
M = [0.0, -1e2, 0.0]
epsilon_b, sigma_b, epsilon_p, sigma_p = strains_and_stresses(F, M, clt)

y = [0.0]
for i = 1:7
    y = [y; 0.25*i - 1e-6; 0.25*i + 1e-6]
end
y= [y; 0.25*8]

data1 = [
0.028975009054690637 38.14473513346067
0.1267656646142703 34.628242119317385
0.22383194494748315 31.153110802240203
0.27888446215139484 29.167330035493514
0.3766751177109746 25.692183743145858
0.47301702281782 22.217067401339154
0.5287939152480986 3.6098894974791307
0.6258601955813113 3.111722149687239
0.7229264759145239 2.654901523690981
0.7787033683448025 24.36077737056926
0.8757696486780152 15.593265663651561
0.9728359290112282 6.90844740032513
1.0286128214415062 0.9947130877232127
1.1256791017747196 0.4965457399313209
1.2227453821079322 -0.0016216078605637563
1.27852227453821 -20.75882904509333
1.3748641796450563 -29.52632577674057
1.4719304599782688 -38.29383748365823
1.5277073524085476 -16.712001802166853
1.6247736327417601 -20.228479841039665
1.7218399130749735 -23.786304601708117
1.7776168055052515 -25.81344706552091
1.8746830858384644 -29.288578382598097
1.9724737413980438 -32.846418118536995
]

# using PyPlot
# close("all"); pygui(true)
# figure()
# plot(y, sigma_b[1, :])
# plot(data1[:, 1], data1[:, 2], "o")

s11mine = linear(y, sigma_b[1, :], data1[:, 1])
n = length(data1[:, 1])
for i = 1:n
    @test isapprox(s11mine[i], data1[i, 2], atol=0.4)
end


data2 = [
    0.026831036983321233 1.8824542507358295
    0.1254532269760696 1.765169773493156
    0.2211747643219727 1.6625805144392833
    0.27701232777374907 1.618668152540212
    0.3741841914430747 1.530790107921348
    0.4713560551124004 1.4429120633024866
    0.5279187817258886 -8.064232926673204
    0.6250906453952142 -6.166816853645008
    0.7222625090645389 -4.357636074734456
    0.7788252356780276 -0.11478106471014016
    0.8745467730239305 0.047335558588925686
    0.971718636693256 0.19475163161712317
    1.0282813633067436 1.128783112229672
    1.125453226976069 2.8717874205519816
    1.2218999274836835 4.673612592245025
    1.2777374909354606 0.7399943479930133
    1.3741841914430741 0.9241724608625264
    1.4720812182741114 1.1230617881670533
    1.5271936185641768 -0.16350029859658566
    1.6250906453952139 -0.28814038305676615
    1.722262509064539 -0.40543019238151157
    1.7781000725163156 -0.4714013778099968
    1.87454677302393 -0.6033997355287175
    1.9724437998549675 -0.7280398199888971
]

# figure()
# plot(y, sigma_b[2, :])
# plot(data2[:, 1], data2[:, 2], "o")

s22mine = linear(y, sigma_b[2, :], data2[:, 1])
n = length(data2[:, 1])
for i = 1:n
    @test isapprox(s22mine[i], data2[i, 2], atol=0.4)
end

data3 = [
    0.02900994112760258 2.4119771295010715
0.1268271911530225 2.265637550457141
0.22247097118809128 2.126634931192692
0.2782636572284641 2.053515759766518
0.3760817952906523 1.9218820696052106
0.47245065734701713 1.7902377230026847
0.5283174944575406 -9.05493972672444
0.6248142338085367 -7.068936074229217
0.7227566970182837 -5.1417453208232615
0.7787567396440485 -1.7810394381568955
0.8751775518513565 -1.05238928512596
0.97232211402481 -0.338439692757035
1.028911369046324 0.7796234635298873
1.1254018921399418 2.6626858938467493
1.2218950793438648 4.589865990811488
1.2782578849894701 1.9579274820294676
1.374679141215163 2.6939305795017146
1.4718250354437676 3.4299390051945755
1.528122570386904 -0.2828823364602542
1.6252142943726444 -0.44393313260741785
1.7223060183583843 -0.6049839287545806
1.7780969283252204 -0.7075148779460028
1.8751886523109609 -0.8685656740931655
1.971556182312173 -1.022268854019627
]

# figure()
# plot(y, sigma_b[4, :])
# plot(data3[:, 1], data3[:, 2], "o")

s12mine = linear(y, sigma_b[4, :], data3[:, 1])
n = length(data3[:, 1])
for i = 1:n
    @test isapprox(s12mine[i], data3[i, 2], atol=0.4)
end

end



# -------- multi-layer composite pipe ----------

@testset "CLT stress.  multi-layer pipe" begin

E1 = 20.59e6
E2 = 1.42e6
nu12 = 0.42
G12 = 0.87e6
rho = 1.0
mat1 = MaterialPlane(E1, E2, G12, nu12, rho)

# t = [2.54e-3, 2.54e-3]
t = [0.1, 0.1]
theta = [0.0, 90.0]*pi/180
laminate1 = Layer.(Ref(mat1), t, theta)

theta = [-45, 45]*pi/180
laminate2 = Layer.(Ref(mat1), t, theta)

# s = 50.8e-3
# r = 10.16e-3
s = 2.0
r = 0.4
b1 = BeamSection(laminate1, [-s/2, s/2], [-r, -r])
# figure(); plot([-s/2, s/2], [-r, -r])
theta = range(-pi/2, pi/2, length=20)
y = r*cos.(theta)
z = r*sin.(theta)
y .+= s/2
# plot(y, z)
b2 = BeamSection(laminate2, y, z)
b3 = BeamSection(laminate1, [s/2, -s/2], [r, r])
# plot([s/2, -s/2], [r, r])
theta = range(pi/2, 3*pi/2, length=20)
y = r*cos.(theta)
z = r*sin.(theta)
y .-= s/2
b4 = BeamSection(laminate2, y, z)
# plot(y, z)

beams = [b1, b2, b3, b4]
closed_section = true
clt = CLT(beams, closed_section)

shear_center = false
S, _, tc = compliance_matrix(clt, shear_center)
K = clt_stiffness_matrix(S)

@test isapprox(K[1, 1], 1.03890e7, rtol=0.005)
@test isapprox(K[4, 4], 6.87060e5, rtol=0.005)
@test isapprox(K[5, 5], 1.88227e6, rtol=0.005)
@test isapprox(K[6, 6], 5.38148e6, rtol=0.005)


# println("K11 = ", round((K[1, 1]/1.03890e7 - 1)*100, digits=2), "%")
# println("K44 = ", round((K[4, 4]/6.87060e5 - 1)*100, digits=2), "%")
# println("K55 = ", round((K[5, 5]/1.88227e6 - 1)*100, digits=2), "%")
# println("K66 = ", round((K[6, 6]/5.38148e6 - 1)*100, digits=2), "%")
# println("K14 = ", round((K[1, 4]/9.83566e4 - 1)*100, digits=2), "%")


F = [0.0; 0; 0]
M = [0.0; -1000.0; 0]
epsilon_b, sigma_b, epsilon_p, sigma_p = strains_and_stresses(F, M, clt)

# extract portion on top
start = 4 + 19*4 + 1
s11 = sigma_b[1, start+4-1:-1:start]
s22 = sigma_b[2, start+4-1:-1:start]
y = [0.0, 0.1-1e-6, 0.1+1e-6, 0.2]



data1 = [
0.00007028112449799367 -0.22500000000000087
0.025115749192054414 -0.2388888888888897
0.05001511094884587 -0.2583333333333342
0.07506054521114758 -0.2750000000000008
0.0999827931253634 -0.2916666666666675
0.10003296012332155 -4.413888888888886
0.12495206414885132 -4.688888888888886
0.1499213351723392 -4.963888888888886
0.17489060619582708 -5.238888888888888
0.20000591591957056 -5.513888888888888
]


# figure()
# plot(y, s11)
# plot(data1[:, 1], data1[:, 2]*1e3, "o")

s11mine = linear(y, s11, data1[:, 1])
for i = 1:length(s11mine)
    @test isapprox(s11mine[i]/1e3, data1[i, 2], atol=0.02)
end

data2 = [
    0.00005049521816367225 0.005689229309128241
0.024988981632448892 0.06940620969508204
0.050000540778775575 0.133123266128051
0.07515811019449115 0.1964241932932616
0.09980430568061074 0.2601408694911542
0.10027153854265528 -0.10119074416567281
0.12496910801244873 -0.10428715048791759
0.14995903600775481 -0.10717511194120638
0.17494889640571387 -0.11027121407538987
0.20008490226775583 -0.11336716411554265

]

# figure()
# plot(y, s22)
# plot(data2[:, 1], data2[:, 2]*1e3, "o")

# s22mine = linear(y, s22, data1[:, 1])
# for i = 1:length(s22mine)
#     @test isapprox(s22mine[i]/1e3, data2[i, 2], atol=0.02)
# end

end

# -------- airfoil ----------

@testset "CLT stress. airfoil" begin

chord = 1.9
web = [0.15, 0.5]
nodes = [0.0, 0.0016, 0.0041, 0.1147, 0.5366, 1.0]
mat = Vector{Material{Float64}}(undef, 5)
mat[1] = MaterialPlane(3.7e10, 9.0e9, 4.0e9, 0.28, 1.0)
mat[2] = MaterialPlane(1.03e10, 1.03e10, 8.0e9, 0.30, 1.0)
mat[3] = MaterialPlane(10.0, 10.0, 1.0, 0.30, 1.0)
mat[4] = MaterialPlane(1.03e10, 1.03e10, 8.0e9, 0.30, 1.0)
mat[5] = MaterialPlane(1.0e7, 1.0e7, 2.0e5, 0.30, 1.0)


t = [0.000381, 0.00051, 18*0.00053]
theta = [0, 0, 20]*pi/180
segment1 = Layer.([mat[3], mat[4], mat[2]], t, theta)
t = [0.000381, 0.00051, 33*0.00053]
theta = [0, 0, 20]*pi/180
segment3 = Layer.([mat[3], mat[4], mat[2]], t, theta)
t = [0.000381, 0.00051, 17*0.00053, 38*0.00053, 1*0.003125, 37*0.00053, 16*0.00053]
theta = [0, 0, 20, 30, 0, 30, 20]*pi/180
idx = [3, 4, 2, 1, 5, 1, 2]
segment4 = Layer.(mat[idx], t, theta)
t = [0.000381, 0.00051, 17*0.00053, 0.003125, 16*0.00053]
theta = [0, 0, 20, 0, 0]*pi/180
idx = [3, 4, 2, 5, 2]
segment5 = Layer.(mat[idx], t, theta)
t = [38*0.00053, 0.003125, 38*0.00053]
theta = [0, 0, 0]*pi/180
idx = [1, 5, 1]
webs = Layer.(mat[idx], t, theta)

yn = [1.00000000, 0.99619582, 0.98515158, 0.96764209, 0.94421447, 0.91510964, 0.88074158, 0.84177999, 0.79894110, 0.75297076, 0.70461763, 0.65461515, 0.60366461, 0.55242353, 0.50149950, 0.45144530, 0.40276150, 0.35589801, 0.31131449, 0.26917194, 0.22927064, 0.19167283, 0.15672257, 0.12469599, 0.09585870, 0.07046974, 0.04874337, 0.03081405, 0.01681379, 0.00687971, 0.00143518, 0.00053606, 0.00006572, 0.00001249, 0.00023032, 0.00079945, 0.00170287, 0.00354717, 0.00592084, 0.01810144, 0.03471169, 0.05589286, 0.08132751, 0.11073805, 0.14391397, 0.18067874, 0.22089879, 0.26433734, 0.31062190, 0.35933893, 0.40999990, 0.46204424, 0.51483073, 0.56767889, 0.61998250, 0.67114514, 0.72054815, 0.76758733, 0.81168064, 0.85227225, 0.88883823, 0.92088961, 0.94797259, 0.96977487, 0.98607009, 0.99640466, 1.00000000]
zn = [0.00000000, 0.00017047, 0.00100213, 0.00285474, 0.00556001, 0.00906779, 0.01357364, 0.01916802, 0.02580144, 0.03334313, 0.04158593, 0.05026338, 0.05906756, 0.06766426, 0.07571157, 0.08287416, 0.08882939, 0.09329359, 0.09592864, 0.09626763, 0.09424396, 0.09023579, 0.08451656, 0.07727756, 0.06875796, 0.05918984, 0.04880096, 0.03786904, 0.02676332, 0.01592385, 0.00647946, 0.00370956, 0.00112514, -0.00046881, -0.00191488, -0.00329201, -0.00470585, -0.00688469, -0.00912202, -0.01720842, -0.02488211, -0.03226730, -0.03908459, -0.04503763, -0.04986836, -0.05338180, -0.05551392, -0.05636585, -0.05605816, -0.05472399, -0.05254383, -0.04969990, -0.04637175, -0.04264894, -0.03859653, -0.03433153, -0.02996944, -0.02560890, -0.02134397, -0.01726049, -0.01343567, -0.00993849, -0.00679919, -0.00402321, -0.00180118, -0.00044469, 0.00000000]


i0 = argmin(yn)
yu = yn[1:i0]
zu = zn[1:i0]

ynew = yu
znew = zu
for i = 2:5
    ix = findlast(ynew .> nodes[i])

    ynew = [ynew[1:ix]; nodes[i]; ynew[ix+1:end]]
    znew = [znew[1:ix]; linear(reverse(yu), reverse(zu), nodes[i]); znew[ix+1:end]]
end

idx = zeros(Int64, 6)
idx[1] = length(ynew)
for i = 2:5
    idx[i] = findfirst(ynew .== nodes[i])
end
idx[6] = 1

ynew *= chord
znew *= chord

b1 = BeamSection(segment5, ynew[idx[6]:idx[5]], znew[idx[6]:idx[5]])
b2 = BeamSection(segment4, ynew[idx[5]:idx[4]], znew[idx[5]:idx[4]])
b3 = BeamSection(segment3, ynew[idx[4]:idx[3]], znew[idx[4]:idx[3]])
b4 = BeamSection(segment1, ynew[idx[3]:idx[2]], znew[idx[3]:idx[2]])
b5 = BeamSection(segment1, ynew[idx[2]:idx[1]], znew[idx[2]:idx[1]])

# figure()
# plot(ynew[idx[6]:idx[5]], znew[idx[6]:idx[5]])
# plot(ynew[idx[5]:idx[4]], znew[idx[5]:idx[4]])
# plot(ynew[idx[4]:idx[3]], znew[idx[4]:idx[3]])
# plot(ynew[idx[3]:idx[2]], znew[idx[3]:idx[2]])
# plot(ynew[idx[2]:idx[1]], znew[idx[2]:idx[1]])

yl = yn[i0:end]
zl = zn[i0:end]

ynew = yl
znew = zl
for i = 2:5
    ix = findlast(ynew .< nodes[i])

    ynew = [ynew[1:ix]; nodes[i]; ynew[ix+1:end]]
    znew = [znew[1:ix]; linear(yl, zl, nodes[i]); znew[ix+1:end]]
end

idx[1] = 1
for i = 2:5
    idx[i] = findfirst(ynew .== nodes[i])
end
idx[6] = length(ynew)

ynew *= chord
znew *= chord

# figure()
# plot(ynew[idx[1]:idx[2]], znew[idx[1]:idx[2]])
# plot(ynew[idx[2]:idx[3]], znew[idx[2]:idx[3]])
# plot(ynew[idx[3]:idx[4]], znew[idx[3]:idx[4]])
# plot(ynew[idx[4]:idx[5]], znew[idx[4]:idx[5]])
# plot(ynew[idx[5]:idx[6]], znew[idx[5]:idx[6]])




b6 = BeamSection(segment1, ynew[idx[1]:idx[2]], znew[idx[1]:idx[2]])
b7 = BeamSection(segment1, ynew[idx[2]:idx[3]], znew[idx[2]:idx[3]])
b8 = BeamSection(segment3, ynew[idx[3]:idx[4]], znew[idx[3]:idx[4]])
b9 = BeamSection(segment4, ynew[idx[4]:idx[5]], znew[idx[4]:idx[5]])
b10 = BeamSection(segment5, ynew[idx[5]:idx[6]], znew[idx[5]:idx[6]])

webseg = Vector{BeamSection}(undef, 2)
for i = 1:2
    wzu = linear(reverse(yu), reverse(zu), web[i])
    wzl = linear(yl, zl, web[i])
    webseg[i] = BeamSection(webs, chord*[web[i], web[i]], chord*[wzl, wzu])
    webseg[i] = BeamSection(webs, chord*[web[i], web[i]], chord*[wzl, wzu])
end

profile = [b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, webseg[1], webseg[2]]

closed_section = true
clt = CLT(profile, closed_section)

shear_center = false
S, _, tc = compliance_matrix(clt, shear_center)
Kp = clt_stiffness_matrix(S)

# W, P, yc, zc = beamstiffness(profile)

# K = fullstiffnessmatrix(P, zeros(2, 2))

# theta = 0.0
# r = [0.0, 0.031, 0.040]
# Kp = rotatestiffnessmatrix(K, r, theta)

abs(Kp[1, 1]/2.389e9 - 1)*100
# abs(Kp[2, 2]/8.252e6 - 1)*100
# abs(Kp[3, 3]/2.444e6 - 1)*100
abs(Kp[4, 4]/2.167e7 - 1)*100
abs(Kp[5, 5]/1.97e7 - 1)*100
abs(Kp[6, 6]/4.406e8 - 1)*100
abs(Kp[1, 4]/-3.382e7 - 1)*100
abs(Kp[1, 5]/-2.627e7 - 1)*100
abs(Kp[1, 6]/-4.736e8 - 1)*100
abs(Kp[4, 5]/-6.279e4 - 1)*100
abs(Kp[4, 6]/1.430e6 - 1)*100
abs(Kp[5, 6]/1.209e7 - 1)*100

F = [0.0; 0; 0]
M = [0.11e6; 0.0; 0]  # 1 in-lb in N-m
epsilon_b, sigma_b, epsilon_p, sigma_p = strains_and_stresses(F, M, clt)

x2 = 0.55
# segment 4 upper surface
# idxs = 140 + 1
# idxf = 140 + 154

# ny = 12
# nl = 7
# (ny-1)*nl*2 = 154
idxs = 140 + 5*7*2 + 1
idxf = 140 + 6*7*2

sigma_b[:, idxs:idxf]

t = [0.000381, 0.00051, 17*0.00053, 38*0.00053, 1*0.003125, 37*0.00053, 16*0.00053]
x3 = [sum(t)]
loc = sum(t)
for i = 1:length(t)-1
    loc -= t[i]
    x3 = [x3; loc+1e-6; loc-1e-6]
end
x3 = [x3; 0.0]

figure()
plot(x3 * 39.3701, sigma_b[1, idxs:idxf] * 0.000145038)  # convert m to in and N/m^2 to psi

figure()
plot(x3 * 39.3701, sigma_b[2, idxs:idxf] * 0.000145038)  # convert m to in and N/m^2 to psi

end