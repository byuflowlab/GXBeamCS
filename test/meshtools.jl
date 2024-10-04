
using GXBeamCS, Test

@testset "Mesh Deformation" begin
    @testset "Mesh Stretching" begin
        xpoints = collect(0:2)
        ypoints = collect(0:2)

        nodes = Vector{Node{Float64}}(undef, length(xpoints)*length(ypoints))
        k = 1
        for j in eachindex(ypoints)
            # global k
            for i in eachindex(xpoints)
                nodes[k] = Node(xpoints[i], ypoints[j])
                k += 1
            end
        end

        ### Test unidirectional mesh shifting
        O = (2.0, 1.0)
        k = (0.5, 1.0)
        nodes_new = GXBeamCS.stretch_mesh.(nodes, Ref(O), Ref(k)) #Note: by testing `stretch_mesh` we are also testing `stretch_mesh!`.
        @test nodes_new[7:9] == [Node(1.0, 2.0), Node(1.5, 2.0), Node(2.0, 2.0)]

        #Test that the orthogonal origin position has no effect
        O = (2.0, 10.0)
        k = (0.5, 1.0)
        nodes_new = GXBeamCS.stretch_mesh.(nodes, Ref(O), Ref(k))
        @test nodes_new[7:9] == [Node(1.0, 2.0), Node(1.5, 2.0), Node(2.0, 2.0)]

        ### Test contraction mesh shifting
        O = (1.0, 1.0)
        k = (0.5, 0.5)
        nodes_new = GXBeamCS.stretch_mesh.(nodes, Ref(O), Ref(k))
        @test nodes_new[7:9] == [Node(0.5, 1.5), Node(1.0, 1.5), Node(1.5, 1.5)]

        ### Test expansion mesh shifting (by undoing the contraction mesh shifting)
        O = (1.0, 1.0)
        k = (2.0, 2.0)
        nodes_new = GXBeamCS.stretch_mesh.(nodes_new, Ref(O), Ref(k))
        @test nodes_new[7:9] == [Node(0.0, 2.0), Node(1.0, 2.0), Node(2.0, 2.0)]

        ### Test bi-directional mesh shifting
        O = (2, 2)
        k = (0.4, 0.6)
        nodes_new = GXBeamCS.stretch_mesh.(nodes, Ref(O), Ref(k))
        @test nodes_new[1:3] == [Node(1.2, 0.8), Node(1.6, 0.8), Node(2.0, 0.8)]
    end

    @testset "Mesh Translation" begin
        ### Test find direction
        @test isapprox(GXBeamCS.find_direction([Node(0.0, 0.0), Node(1.0, 0.0)]), [1.0, 0.0])
        @test isapprox(GXBeamCS.find_direction([Node(0.0, 0.0), Node(0.0, 1.0)]), [0.0, 1.0])
        @test isapprox(GXBeamCS.find_direction([Node(0.0, 0.0), Node(1.0, 1.0)]), [1/sqrt(2), 1/sqrt(2)])
        @test isapprox(GXBeamCS.find_direction([Node(0.0, 0.0), Node(1.0, -1.0)]), [1/sqrt(2), -1/sqrt(2)])

        ### Find k
        ## Node and Origin are in different locations, same amount in each direction
        n1 = Node(0.0, 0.0) 
        O = [1.0, 1.0]

        #Negative stretch
        d = [-.1, -.1]
        k = GXBeamCS.find_k(n1, O, d)
        nn = GXBeamCS.stretch_mesh(n1, O, k) #Tested previously
        @test isapprox(nn.x, -0.1)&&isapprox(nn.y, -0.1)

        #Positive stretch
        d = [0.1, 0.1]
        k = GXBeamCS.find_k(n1, O, d)
        nn = GXBeamCS.stretch_mesh(n1, O, k)
        @test isapprox(nn.x, 0.1)&&isapprox(nn.y, 0.1)

        ## Node and Origin are in different locations, different amount in each direction
        n1 = Node(0.1, 0.2)
        O = [1.1, 2.1]
        d = [0.1, -0.1]
        k = GXBeamCS.find_k(n1, O, d)
        nn = GXBeamCS.stretch_mesh(n1, O, k)
        @test isapprox(nn.x, 0.2)&&isapprox(nn.y, 0.1)

        ## Node and Origin have a common distance. 
        n = Node(1.0, 1.0)
        O = [3, 1] 
        d = [0.1, 0.1]
        k = GXBeamCS.find_k(n, O, d)
        nn = GXBeamCS.stretch_mesh(n, O, k)
        @test isapprox(nn.x, 1.1)
        @test isapprox(nn.y, 1.0) #Note: the y value should not change because it lines up with the origin. #Todo: Failing.  


        ### Mesh Translation
        nodes = [Node(0.0, 0.0), Node(1.0, 0.0), Node(1.0, 1.0), Node(0.0, 1.0)]
        nhat = [1.0, 1.0]./sqrt(2)
        t = 0.1
        d = t*nhat
        nn = GXBeamCS.translate_mesh.(nodes, Ref(d))
        @test isapprox(nn[1].x, t/sqrt(2))&&isapprox(nn[1].y, t/sqrt(2))
        @test isapprox(nn[2].x, 1.0 + t/sqrt(2))&&isapprox(nn[2].y, t/sqrt(2))
        @test isapprox(nn[3].x, 1.0 + t/sqrt(2))&&isapprox(nn[3].y, 1.0 + t/sqrt(2))
        @test isapprox(nn[4].x, t/sqrt(2))&&isapprox(nn[4].y, 1.0 + t/sqrt(2))
    end

    @testset "Layer Translation" begin
        xpoints = collect(0:2)
        ypoints = collect(0:2)
        nodes = Vector{Node{Float64}}(undef, length(xpoints)*length(ypoints))

        k = 1
        for j in eachindex(ypoints)
            for i in eachindex(xpoints)
                nodes[k] = Node(xpoints[i], ypoints[j])
                k += 1
            end
        end

        l1 = view(nodes, 1:3)
        l2 = view(nodes, 4:6)
        l3 = view(nodes, 7:9)

        layers = [l1, l2, l3]
        t = [-0.1, -0.15] #Note: the translation vector is in the opposite direction of the mesh direction.
        GXBeamCS.thicken_layers!(layers, t)

        @test nodes[1].x == 0.0 && nodes[1].y == 0.0
        @test nodes[2].x == 1.0 && nodes[2].y == 0.0
        @test nodes[3].x == 2.0 && nodes[3].y == 0.0
        @test nodes[4].x == 0.0 && nodes[4].y == 0.9
        @test nodes[5].x == 1.0 && nodes[5].y == 0.9
        @test nodes[6].x == 2.0 && nodes[6].y == 0.9
        @test nodes[7].x == 0.0 && nodes[7].y == 1.85
        @test nodes[8].x == 1.0 && nodes[8].y == 1.85
        @test nodes[9].x == 2.0 && nodes[9].y == 1.85


        ### Todo: Test if the layers are upside down
        xpoints = collect(0:2)
        ypoints = collect(0:2)
        nodes = Vector{Node{Float64}}(undef, length(xpoints)*length(ypoints))

        k = 1
        for j in eachindex(ypoints)
            for i in eachindex(xpoints)
                nodes[k] = Node(xpoints[i], ypoints[j])
                k += 1
            end
        end

        l1 = view(nodes, 1:3)
        l2 = view(nodes, 4:6)
        l3 = view(nodes, 7:9)

        layers = reverse([l1, l2, l3])
        t = [-0.1, -0.15]

        GXBeamCS.thicken_layers!(layers, t)

        @test nodes[1].x == 0.0 && nodes[1].y == 0.15
        @test nodes[2].x == 1.0 && nodes[2].y == 0.15
        @test nodes[3].x == 2.0 && nodes[3].y == 0.15
        @test nodes[4].x == 0.0 && nodes[4].y == 1.1
        @test nodes[5].x == 1.0 && nodes[5].y == 1.1
        @test nodes[6].x == 2.0 && nodes[6].y == 1.1
        @test nodes[7].x == 0.0 && nodes[7].y == 2.0
        @test nodes[8].x == 1.0 && nodes[8].y == 2.0
        @test nodes[9].x == 2.0 && nodes[9].y == 2.0


        ### Todo: Test if the layers are rotated

        ### Todo: Test if the user puts the layers in backwards

        ### Todo: Test if the user puts the wrong number of layers in. 
    end

    @testset "Mesh Selection" begin
        ### Test the Winding Algorithm
        polygon = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
        
        #Inside the polygon
        @test GXBeamCS.is_point_in_polygon(polygon, (2.0, 2.0)) 
        @test GXBeamCS.is_point_in_polygon(polygon, (3.0, 1.0))
        @test GXBeamCS.is_point_in_polygon(polygon, (1.2, 3.99))

        #Outside the polygon
        @test !GXBeamCS.is_point_in_polygon(polygon, (2.5, 5.0)) #Top 
        @test !GXBeamCS.is_point_in_polygon(polygon, (2.5, -1.0)) #Bottom
        @test !GXBeamCS.is_point_in_polygon(polygon, (0.0, 5.0)) #Right
        @test !GXBeamCS.is_point_in_polygon(polygon, (-0.25, 1.0)) #Left
        @test !GXBeamCS.is_point_in_polygon(polygon, (5.0, 5.0)) #Top right

        #Polygon edge
        @test GXBeamCS.is_point_in_polygon(polygon, (2.0, 4.0)) #Top
        @test GXBeamCS.is_point_in_polygon(polygon, (4.0, 2.0)) #Right
        @test GXBeamCS.is_point_in_polygon(polygon, (2.0, 0.0)) #Bottom
        @test GXBeamCS.is_point_in_polygon(polygon, (0.0, 2.0)) #Left 

        #Polygon vertex
        @test GXBeamCS.is_point_in_polygon(polygon, (0.0, 0.0)) #Bottom left
        @test GXBeamCS.is_point_in_polygon(polygon, (4.0, 0.0)) #Bottom right
        @test GXBeamCS.is_point_in_polygon(polygon, (4.0, 4.0)) #Top right
        @test GXBeamCS.is_point_in_polygon(polygon, (0.0, 4.0)) #Top left

        #Displaced Polygon
        polygon = [(1.0, 1.0), (5.0, 1.0), (5.0, 5.0), (1.0, 5.0)]
        @test GXBeamCS.is_point_in_polygon(polygon, (2.0, 2.0))

        #Clockwise polygon
        polygon = reverse(polygon)
        @test GXBeamCS.is_point_in_polygon(polygon, (2.0, 2.0))

        #Polygons with more sides #Note: more sides means more tests, means longer function call. 
        polygon = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (2.0, 6.0), (0.0, 4.0)]
        @test GXBeamCS.is_point_in_polygon(polygon, (2.0, 2.0)) #todo: I could add more shapes here. 

        #Polygon with 3 sides
        polygon = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0)]
        @test GXBeamCS.is_point_in_polygon(polygon, (3.0, 2.0))

        #Polygon with 2 sides
        polygon = [(0.0, 0.0), (4.0, 0.0)]
        try 
            GXBeamCS.is_point_in_polygon(polygon, (3.0, 2.0))
            @test false
        catch e
            @test e == ErrorException("winding_number: Polygon must have at least 3 vertices")
        end

        #Todo: Test if having the first and last point the same will break it. (i.e. a quadrilateral that is really a triangle)


        ### Test the Mesh selection function
        xpoints = collect(0:2)
        ypoints = collect(0:2)

        nodes = Vector{Node{Float64}}(undef, length(xpoints)*length(ypoints))
        k = 1
        for j in eachindex(ypoints)
            # global k
            for i in eachindex(xpoints)
                nodes[k] = Node(xpoints[i], ypoints[j])
                k += 1
            end
        end

        points = [(-0.1, 1.1), (2.1, 1.1), (2.1, 3.1), (-0.1, 3.1)]
        selected_nodes, idxs = GXBeamCS.select_nodes(nodes, points)
        @test sort(idxs) == [7, 8, 9]
    end

end

@testset "Mesh Quality" begin

    @testset "Cell Area" begin
        N1 = Node(0.0, 0.0)
        N2 = Node(1.0, 0.0)
        N3 = Node(1.0, 1.0)
        N4 = Node(0.0, 1.0)

        @test isapprox(GXBeamCS.signed_area([N1, N2, N3, N4]), 1.0) #Counter-clockwise
        @test isapprox(GXBeamCS.signed_area([N1, N4, N3, N2]), -1.0) #clockwise
    end

    @testset "Cell self-intersection" begin #Checking self-intersection 
        p1, p2, p3, p4 = (0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)  # Square

        @test !GXBeamCS.does_quadrilateral_hourglass([p1, p2, p3, p4]) #False
        @test GXBeamCS.does_quadrilateral_hourglass([p1, p2, p4, p3]) #True
        @test GXBeamCS.does_quadrilateral_hourglass([p1, p3, p2, p4]) #True
        @test GXBeamCS.does_quadrilateral_hourglass([p1, p3, p4, p2]) #True
        @test GXBeamCS.does_quadrilateral_hourglass([p1, p4, p2, p3]) #True
        @test !GXBeamCS.does_quadrilateral_hourglass([p1, p4, p3, p2]) #False
    end

    @testset "Cell-Cell Intersection" begin
        polygon1 = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]  # Square
        polygon2 = [(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)]  # Overlapping square
        polygon3 = [(5.0, 5.0), (9.0, 5.0), (9.0, 9.0), (5.0, 9.0)]  # Non-overlapping square
        polygon4 = [(0.0, 4.0), (4.0, 4.0), (4.0, 8.0), (0.0, 8.0)]  # Touching square
        polygon5 = [(2.0, 2.0), (6.0, 0.0), (7.5, 2.0), (4.0, 4.0)]  # Intersecting quad with shared node. 

        @test GXBeamCS.polygons_intersect(polygon1, polygon2)
        @test !GXBeamCS.polygons_intersect(polygon1, polygon3)
        @test !GXBeamCS.polygons_intersect(polygon1, polygon4)
        @test GXBeamCS.polygons_intersect(polygon1, polygon5)
    end

    @testset "Mesh test" begin
        @test true #Todo: 
    end

end