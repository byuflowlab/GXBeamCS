#=
Tools to work with the mesh:
- mesh deformation
- mesh quality checking
=#


### ----- Mesh deformation ----- ###
#=
General mesh deformation approach is to select some nodes to deform and select some nodes to respond. 
Deformation is done by deciding a direction (either by user or by the find direction function)
and then using either the translate or expand functions. Resonse is done by the contract function. 
=#

"""
    is_left(p1, p2, point)

Helper function to determine if a point is to the left, right, or on an edge
"""
function is_left(p1, p2, point)
    return (p2[1] - p1[1]) * (point[2] - p1[2]) - (point[1] - p1[1]) * (p2[2] - p1[2])
end

"""
    winding_number(polygon, point)

Function to calculate the winding number of a point with respect to a polygon
"""
function winding_number(polygon, point)
    n = length(polygon)
    if n < 3
        error("winding_number: Polygon must have at least 3 vertices")
    end
    wn = 0  # Winding number counter

    for i in 1:n
        j = i % n + 1  # Wrap around to the first vertex
        # print(i, ": ")

        if polygon[i][2] <= point[2] #Above the ith point
            # And below the jth (next) point && 
            # print("a")
            il = is_left(polygon[i], polygon[j], point)
            
            if polygon[j][2] >= point[2] && is_left(polygon[i], polygon[j], point) >= 0
                wn += 1  # Upward crossing
                # print("b")
            end
            # print(" il: ", il)
        else #Below the ith point
            # and above the jth point &&
            # print("c")
            il = is_left(polygon[i], polygon[j], point)
            
            if polygon[j][2] <= point[2] && is_left(polygon[i], polygon[j], point) < 0
                wn -= 1  # Downward crossing
                # print("d")
            end
            # print(" il: ", il)
        end
        # print("\n")
    end

    return wn
end

# Main function to check if a point is inside a polygon using the winding number algorithm
function is_point_in_polygon(polygon, point) #Todo: Does the polygon have to close? I don't think it does. 
    return winding_number(polygon, point) != 0
end

"""
    select_nodes(nodes, points)
Given four points that define a quadrilateral, return the nodes that are inside the quadrilateral. 

**Arguments**
- `nodes::Array{Node}`: The array of all the nodes in the mesh.
- `points::Array{Node}`: The four points that define the quadrilateral starting from top left and going clockwise.
"""
function select_nodes(nodes, points)
    #Check if the node is inside the quadrilateral
    flags = [is_point_in_polygon(points, [node.x, node.y]) for node in nodes]
    
    n = sum(Int, flags)
    idxs = Vector{Int}(undef, n)
    k = 1
    for i in eachindex(nodes)
        if flags[i]
            idxs[k] = i
            k += 1
        end
    end
    return nodes[idxs], idxs
end



"""
    find_direction(nodes)
Given two nodes, find the direction of the layer.
"""
function find_direction(nodes)
    nx = nodes[2].x - nodes[1].x
    ny = nodes[2].y - nodes[1].y
    n = sqrt(nx^2 + ny^2)

    return [nx/n, ny/n]
end

"""
    find_k(node, O, d)
Given a node in the global reference frame, find the vector `k` that will stretch the node by a distance `d` with respect to the origin `O`.

**Arguments**
- `node::Node`: The node to be stretched.
- `O::Array{Float64}`: The origin of the local reference frame, length 2. 
- `d::Array{Float64}`: The distance to stretch the node, length 2.
"""
function find_k(node, O, d; verbose=false)
    P = [node.x, node.y]
    Pp = P - O #Map the point to the local reference frame.

    kx = d[1]/Pp[1] + 1 #Find the shift. 
    ky = d[2]/Pp[2] + 1

    if isapprox(Pp[1], 0)
        kx /= kx #Divide by kx to get kx=1 with type stability. 
        verbose ? println("X distance to origin is zero, returning kx=1") : nothing
    end

    if isapprox(Pp[2], 0)
        ky /= ky
        verbose ? println("Y distance to origin is zero, returning ky=1") : nothing
    end

    return [kx, ky]
end



"""
    stretch_mesh(nodes, O, k)

Deform the nodes about an origin `O` by a factor `k`. `k`` is a 2D vector.
The origin `O` is the origin of a local reference frame (parallel with the global frame),
i.e. a point in the mesh that is not moved.
"""
function stretch_mesh(node, O, k)

    
    x = O[1] + k[1] * (node.x - O[1])
    y = O[2] + k[2] * (node.y - O[2])
    #Inside the parenthesis maps the node to the local reference frame.
    #Then we stretch the node in the local reference frame by multiplying by kx and ky.
    #Then we map the node back to the global reference frame.

    return Node(x, y)
end



"""
    translate_mesh(node, d)

Translate the node by a distance `d`.

**Arguments**
- `node::Node`: The node to be translated.
"""
function translate_mesh(node, d)
    #Note: I can start with a single direction for all the nodes in the layer, then move to individual directions for each node. 

    x = node.x + d[1]
    y = node.y + d[2]
    return Node(x, y)
end


"""
thicken_layers!(layers, t)

Thicken each layer in the mesh by a distance in the vector `t`.

**Arguments**
- `t::Array{Float64}`: The translation vector.
- `layer_nodes::Array{Array{Node}}`: The layers of nodes in the mesh. Note that the length of layers must be one more than the length of `t` and the length of each individual layer must be the same. The layers should be in order of outside to inside. 

**Notes**
- The first layer of nodes does not move. 
"""
function thicken_layers!(layer_nodes, t)
    nl = length(layer_nodes) #Number of layer_nodes
    layer_length = length(layer_nodes[1]) #Number of nodes in each layer #Todo: Why not do a matrix of nodes.... This will give us the size constraints automatically. 

    # Function input checks. 
    for i in eachindex(layer_nodes)
        if length(layer_nodes[i]) != layer_length
            error("thicken_layers!(): All layer_nodes must have the same number of nodes.")
        end
    end
    if length(t) != nl - 1
        error("thicken_layers!(): The number of layer_nodes must be one more than the number of thicknesses.")
    end

    #Translate the layer_nodes individually. 
    for i in nl:-1:2 #Iterate backwards through the layer_nodes so you don't move the layer_nodes that you're using to find the unit normal. 
        for j in 1:layer_length
            # Find the direction to move the node.
            d = t[i-1] .* find_direction((layer_nodes[i-1][j], layer_nodes[i][j]))  #todo: This might have extra allocations. 

            # Move the node in the direction.
            layer_nodes[i][j] = translate_mesh(layer_nodes[i][j], d)
        end
    end
end





### ----- Mesh quality checking ----- ###
"""
    element_center(element)

Return the midpoint of an element. 
"""
function element_center(nodes, element)
    nodes_interest = nodes[element.nodenum]
    xbar = 0. #Todo: Typing
    ybar = 0.

    for node in nodes_interest
        xbar += node.x
        ybar += node.y
    end
    return xbar, ybar
end


"""
    signed_area(nodes)
Use the shoelace formula to calculate the signed area of a mesh element. 

**Arguments**
- `nodes::Vector{Node}`: A vector of the `Node`s of a cell to calculate the area of.
"""
function signed_area(nodes) #Todo: Test me!
    n = length(nodes)
    area = 0.0

    for i in 1:n
        j = (i % n) + 1  # Wrap around to the first vertex
        area += nodes[i].x * nodes[j].y - nodes[j].x * nodes[i].y
    end

    return area / 2.0
end

"""
    cell_sizes(nodes, element)
Calculate the size of all the cells. 
"""
function cell_sizes(nodes, elements)
    ncells = length(elements)
    A = zeros(ncells)

    for i in eachindex(elements)
        # element = elements[i]
        # nodes_interest = nodes[element.nodenum]
        A[i] = signed_area(nodes[elements[i].nodenum])
    end
    # @show A
    return A
end

"""
    cell_sizes(fem)
Calculate the size of all the cells in a mesh.
"""
function cell_sizes(fem)
    return cell_sizes(fem.nodes, fem.elements)
end


function orientation(p, q, r)
    val = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])
    if val == 0
        return 0  # collinear
    elseif val > 0
        return 1  # clockwise
    else
        return 2  # counterclockwise
    end
end

# Helper function to check if a point is on a line segment
function on_segment(p, q, r)
    return min(p[1], r[1]) <= q[1] <= max(p[1], r[1]) && min(p[2], r[2]) <= q[2] <= max(p[2], r[2])
end

# Helper function to check if two line segments intersect
function segments_intersect(p1, q1, p2, q2)
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    if o1 != o2 && o3 != o4
        return true
    end

    if o1 == 0 && on_segment(p1, p2, q1)
        return true
    end

    if o2 == 0 && on_segment(p1, q2, q1)
        return true
    end

    if o3 == 0 && on_segment(p2, p1, q2)
        return true
    end

    if o4 == 0 && on_segment(p2, q1, q2)
        return true
    end

    return false
end

# Function to check if a quadrilateral hourglasses
function does_quadrilateral_hourglass(quadrilateral)
    if length(quadrilateral) != 4
        error("Input must be a quadrilateral with 4 vertices.")
    end

    (p1, p2, p3, p4) = quadrilateral

    # Check intersections between non-adjacent sides
    if segments_intersect(p1, p2, p3, p4) || segments_intersect(p2, p3, p4, p1)
        return true
    end

    return false
end

"""
    check_mesh_hourglass(nodes, elements)

Check if the mesh has any hourglassing elements.
"""
function check_mesh_hourglass(nodes, elements)
    vec = zeros(Bool, length(elements))

    for j in eachindex(elements)
        element = elements[j]
        points = [[nodes[i].x, nodes[i].y] for i in element.nodenum]

        vec[j] = does_quadrilateral_hourglass(points)
    end

    return vec
end

"""
    check_mesh_hourglass(fem)

Check if the mesh has any hourglassing elements.
"""
function check_mesh_hourglass(fem)
    return check_mesh_hourglass(fem.nodes, fem.elements)
end


####--------------- ChatGPT 4.0 code ----------- ###
# using LinearAlgebra

# Helper function to check if two line segments intersect
function segments_intersect2(p1, p2, p3, p4)
    # Vector cross product
    function cross(v1, v2)
        return v1[1] * v2[2] - v1[2] * v2[1]
    end

    # Check if point q lies on line segment pr
    function on_segment(p, q, r)
        return min(p[1], r[1]) ≤ q[1] ≤ max(p[1], r[1]) && min(p[2], r[2]) ≤ q[2] ≤ max(p[2], r[2])
    end

    # Directions for the segments
    d1 = cross(p4 - p3, p1 - p3)
    d2 = cross(p4 - p3, p2 - p3)
    d3 = cross(p2 - p1, p3 - p1)
    d4 = cross(p2 - p1, p4 - p1)

    if d1 * d2 < 0 && d3 * d4 < 0
        return true
    elseif d1 == 0 && on_segment(p3, p1, p4)
        return true
    elseif d2 == 0 && on_segment(p3, p2, p4)
        return true
    elseif d3 == 0 && on_segment(p1, p3, p2)
        return true
    elseif d4 == 0 && on_segment(p1, p4, p2)
        return true
    else
        return false
    end
end

# Main function to check if a polygon is self-intersecting
function is_self_intersecting(polygon)
    n = length(polygon)
    for i in 1:n
        for j in i+2:n
            if j != i % n + 1 && segments_intersect2(polygon[i], polygon[i % n + 1], polygon[j], polygon[j % n + 1])
                return true
            end
        end
    end
    return false
end

# Example usage
# polygon = [[0, 0], [4, 0], [4, 4], [0, 4], [2, 2]]
# println(is_self_intersecting(polygon))  # Should print true for a self-intersecting polygon



# using DataStructures

# # Define a structure for an event
# struct Event
#     point::Tuple{Float64, Float64}
#     edge_index::Int
#     event_type::Symbol  # :left or :right
# end

# # Compare events by their x-coordinates, breaking ties by y-coordinates
# function Base.isless(e1::Event, e2::Event)
#     if e1.point[1] == e2.point[1]
#         return e1.point[2] < e2.point[2]
#     else
#         return e1.point[1] < e2.point[1]
#     end
# end

# # Helper function to check if two line segments intersect
# function segments_intersect(p1, p2, p3, p4)
#     function orientation(p, q, r)
#         val = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])
#         return val == 0 ? 0 : (val > 0 ? 1 : 2)
#     end

#     o1 = orientation(p1, p2, p3)
#     o2 = orientation(p1, p2, p4)
#     o3 = orientation(p3, p4, p1)
#     o4 = orientation(p3, p4, p2)

#     if o1 != o2 && o3 != o4
#         return true
#     end

#     return false
# end

# # Main function to check if a polygon is self-intersecting using Bentley-Ottmann algorithm
# function is_self_intersecting(polygon)
#     n = length(polygon)
#     events = PriorityQueue{Event, Float64}()

#     # Insert all segment endpoints as events
#     for i in 1:n
#         p1 = polygon[i]
#         p2 = polygon[i % n + 1]
#         if p1[1] <= p2[1]
#             push!(events, Event(p1, i, :left) => p1[1])
#             push!(events, Event(p2, i, :right) => p2[1])
#         else
#             push!(events, Event(p2, i, :left) => p2[1])
#             push!(events, Event(p1, i, :right) => p1[1])
#         end
#     end

#     status = SortedDict{Int, Int}()
    
#     while !isempty(events)
#         event, _ = dequeue(events)
#         p, i, etype = event.point, event.edge_index, event.event_type
#         if etype == :left
#             status[i] = p[2]
#             for (j, _) in status
#                 if j != i && segments_intersect(polygon[i], polygon[i % n + 1], polygon[j], polygon[j % n + 1])
#                     return true
#                 end
#             end
#         else
#             delete!(status, i)
#         end
#     end

#     return false
# end

# # Example usage
# polygon = [[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0], [2.0, 2.0]]
# println(is_self_intersecting(polygon))  # Should print true for a self-intersecting polygon
