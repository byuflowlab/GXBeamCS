using GXBeamCS
using Test

@testset "construct box" begin
    
    ### Basic test
    ncells = 3
    xyb, nnb = GXBeamCS.create_edge(1, 0.1, ncells)
    xy_gold = [0 0; 0 0.1; 1/3 0; 1/3 0.1; 2/3 0; 2/3 0.1; 1 0;1 0.1] 
    nn_gold = [1 3 4 2; 3 5 6 4; 5 7 8 6]
    @test xyb ≈ xy_gold
    @test nnb == nn_gold


    ### Test more cells
    ncells = 4
    xyb, nnb = GXBeamCS.create_edge(1, 0.1, ncells)
    xy_gold = [0 0; 0 0.1; 1/4 0; 1/4 0.1; 2/4 0; 2/4 0.1; 3/4 0; 3/4 0.1; 1 0; 1 0.1] 
    nn_gold = [1 3 4 2; 3 5 6 4; 5 7 8 6; 7 9 10 8]
    @test xyb ≈ xy_gold
    @test nnb == nn_gold


    ### Test an odd length
    xy, nn = GXBeamCS.create_edge(1.3, 0.1, 5)
    xy_gold = [0 0; 0 0.1; 0.26 0; 0.26 0.1; 0.52 0; 0.52 0.1; 0.78 0; 0.78 0.1; 1.04 0; 1.04 0.1; 1.3 0; 1.3 0.1]
    nn_gold = [1 3 4 2; 3 5 6 4; 5 7 8 6; 7 9 10 8; 9 11 12 10]
    @test xy ≈ xy_gold
    @test nn == nn_gold


    ### Test a swapped edge
    ncells = 4
    xyb, nnb = GXBeamCS.create_edge(1, 0.1, ncells; swap=true)
    xy_gold = [0 0; 0 0.1; 1/4 0; 1/4 0.1; 2/4 0; 2/4 0.1; 3/4 0; 3/4 0.1; 1 0; 1 0.1] 
    xy_gold[:, 1], xy_gold[:, 2] = xy_gold[:, 2], xy_gold[:, 1]
    nn_gold = [2 4 3 1; 4 6 5 3; 6 8 7 5; 8 10 9 7]
    @test xyb ≈ xy_gold
    @test nnb == nn_gold


    # plt = scatter(xyb[:,1], xyb[:,2], legend = false, xlims=(-0.1, 1.1), ylims=(-0.1, 1.1))
    # for i in 1:size(xyb, 1) 
    #     annotate!(plt, xyb[i, 1], xyb[i, 2], text(string(i), :left))
    # end
    # display(plt)


    ### Test a shifted edge
    xyb, nnb = GXBeamCS.create_edge(1, 0.1, 4)
    xyb[:,1] .+= 0.5
    @test all(xyb[:,1] .>= 0.5)
    @test xyb[1,1] == 0.5
    @test xyb[end,1] == 1.5

    ### Insert extra points.
    xy, nn = GXBeamCS.create_edge(1, 0.1, 4; xextra=[0.1, 0.9])
    xy_gold = [0 0; 0 0.1; 0.1 0; 0.1 0.1; 0.25 0; 0.25 0.1; 0.5 0; 0.5 0.1; 0.75 0; 0.75 0.1; 0.9 0; 0.9 0.1; 1.0 0.0; 1.0 0.1]
    nn_gold = [1 3 4 2; 3 5 6 4; 5 7 8 6; 7 9 10 8; 9 11 12 10; 11 13 14 12]

    @test xy ≈ xy_gold
    @test nn == nn_gold

    ### Insert repeated points
    xy, nn = GXBeamCS.create_edge(1, 0.1, 4; xextra=[0.25, 0.9])
    xy_gold = [0 0; 0 0.1; 0.25 0; 0.25 0.1; 0.5 0; 0.5 0.1; 0.75 0; 0.75 0.1; 0.9 0; 0.9 0.1; 1.0 0.0; 1.0 0.1]
    nn_gold = [1 3 4 2; 3 5 6 4; 5 7 8 6; 7 9 10 8; 9 11 12 10]

    @test xy ≈ xy_gold
    @test nn == nn_gold


    ### Insert close points
    xy, nn = GXBeamCS.create_edge(1, 0.1, 4; xextra=[(0.25-1e-9), 0.9])
    xy_gold = [0 0; 0 0.1; 0.25 0; 0.25 0.1; 0.5 0; 0.5 0.1; 0.75 0; 0.75 0.1; 0.9 0; 0.9 0.1; 1.0 0.0; 1.0 0.1]
    nn_gold = [1 3 4 2; 3 5 6 4; 5 7 8 6; 7 9 10 8; 9 11 12 10]

    @test xy ≈ xy_gold
    @test nn == nn_gold

    ### Todo: Rotate an edge (about the center of the edge)
    # xyb, nnb = GXBeamCS.create_edge(1, 0.1, 4)

    # xyb_trans = xyb .- [0.5 0.05]
    # # plt = scatter(xyb_trans[:,1], xyb_trans[:,2], legend = false, aspect_ratio = :equal)
    # # display(plt)
    # R = [cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)]
    # n, _ = size(xyb_trans)
    # xyb_rot = reduce(hcat, [R*xyb_trans[i,:] for i in 1:n])'
    # plt = scatter(xyb_rot[:,1], xyb_rot[:,2], legend = false, aspect_ratio = :equal)
    # display(plt)
    #Todo: How to translate it back so it rotated about the center of the edge in the original frame? -> Does matter right now. 


    # n, _ = size(xyb)
    # xyb_rot = reduce(hcat, [GXBeamCS.rotate_points(xyb[i,:], pi/4, [0.5, -0.05]) for i in 1:n])'
    # display(xyb_rot)
    # plt = scatter(xyb_rot[:,1], xyb_rot[:,2], legend = false, aspect_ratio = :equal)
    # display(plt)


    ### Construct a box from edges
    w = 1
    h = 1
    t = 0.1
    nw = 4
    nh = 4

    xy, nn = GXBeamCS.mesh_box(w, h, t, t, t, t, nw, nh)
    nn_gold = [1 3 4 2; 3 5 6 4; 5 7 8 6; 7 9 10 8; 9 11 12 10; 11 13 14 12;
     4 16 15 2; 16 18 17 15; 18 20 19 17; 20 29 27 19;
      14 22 21 12; 22 24 23 21; 24 26 25 23; 26 39 37 25;
       27 29 30 28; 29 31 32 30; 31 33 34 32; 33 35 36 34; 35 37 38 36; 37 39 40 38]

    # for i in 1:size(nn, 1)
    #     println("test $i: ", nn[i, :] == nn_gold[i, :])
    # end
    @test nn == nn_gold #Check the connectivity graph
    @test size(unique(xy, dims=1))==(40,2) #Check that the points are unique
    #check a couple of random points
    @test xy[1,:] == [0.0, 0.0]
    @test xy[7,:] == [0.5, 0.0]
    @test xy[8,:] == [0.5, t]
    @test xy[13,:] == [1.0, 0.0]
    @test xy[14,:] == [1.0, t]
    @test xy[17, :] == [0.0, 0.5]
    @test xy[18, :] == [t, 0.5]
    @test xy[23, :] == [1-t, 0.5]
    @test xy[24, :] == [1.0, 0.5]
    @test xy[28,:] == [0.0, 1.0]
    @test xy[40,:] == [1.0, 1.0]



    # plt = scatter(xy[:,1], xy[:,2], legend = false, xlims=(-0.1, 1.1), ylims=(-0.1, 1.1))
    # for i in 1:size(xy, 1) 
    #     annotate!(plt, xy[i, 1], xy[i, 2], text(string(i), :left))
    # end
    # display(plt)


    ### Test a box with different a different number of cells
    xy, nn = GXBeamCS.mesh_box(w, h, t, t, t, t, 5, 6)
    nn_gold = [1 3 4 2; 3 5 6 4; 5 7 8 6; 7 9 10 8; 9 11 12 10; 11 13 14 12; 13 15 16 14;
     4 18 17 2; 18 20 19 17; 20 22 21 19; 22 24 23 21; 24 26 25 23; 26 39 37 25;
      16 28 27 14; 28 30 29 27; 30 32 31 29; 32 34 33 31; 34 36 35 33; 36 51 49 35;
       37 39 40 38; 39 41 42 40; 41 43 44 42; 43 45 46 44; 45 47 48 46; 47 49 50 48; 49 51 52 50]

    # for i = 1:size(nn, 1)
    #     println("test $i: ", nn[i, :] == nn_gold[i, :])
    # end
    @test nn == nn_gold #Check the connectivity graph
    @test size(unique(xy, dims=1))==(52,2) #Check that the points are unique
    #Todo: Test some random points

    n = size(xy, 1)

    # plt = scatter(xy[:,1], xy[:,2], leg=false, seriescolor=:black)
    # for i in 1:n
    #     annotate!(plt, xy[i, 1], xy[i, 2], text(string(i), :left))
    # end
    # display(plt)
    

    ### Todo:Test the wrapper
    xy, nn = GXBeamCS.mesh_box(w, h, t)

    # plt = scatter(xy[:,1], xy[:,2], leg=false, seriescolor=:black)
    # display(plt)


    ### Todo:Test a box with different thicknesses
    tb = 0.2
    tt = 0.1
    tl = tr = 0.05
    xy, nn = GXBeamCS.mesh_box(w, h, tt, tb, tl, tr, 5, 5)

    # plt = scatter(xy[:,1], xy[:,2], leg=false, seriescolor=:black)
    # display(plt)



    ### Todo:Test a box with thicknesses that coincide with discretization
    

end