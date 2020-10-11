module MeshSignedDistance

using GeometryBasics
using LinearAlgebra


mag2(x) = sum(v^2 for v in x)

# find distance x0 is from segment x1-x2
function point_segment_distance(x0, x1, x2)

   dx = x2-x1
   m2 = mag2(dx)
   # find parameter value of closest point on segment
   s12 = dot(x2-x0, dx)/m2
   if s12<0
      s12 = zero(s12)
   elseif s12>1
      s12 = one(s12)
   end
   # and find the distance
   return norm(x0 - s12*x1+(1-s12)*x2)

end

# find distance x0 is from triangle x1-x2-x3
function point_triangle_distance(p, face)

    x0 = p
    x1 = face.points[1]
    x2 = face.points[2]
    x3 = face.points[3]

    # first find barycentric coordinates of closest point on infinite plane
    x13 = x1-x3
    x23 = x2-x3
    x03 = x0-x3
    m13 = mag2(x13)
    m23 = mag2(x23)
    d = dot(x13,x23)
    v = m13*m23-d*d
    invdet = 1/max(v, eps(v))
    a = dot(x13,x03)
    b = dot(x23,x03)
    # the barycentric coordinates themselves
    w23=invdet*(m23*a-d*b)
    w31=invdet*(m13*b-d*a)
    w12=1-w23-w31
    if w23>=0 && w31>=0 && w12>=0 # if we're inside the triangle
        return norm(x0 - w23*x1+w31*x2+w12*x3)
    else # we have to clamp to one of the edges
        if w23>0 # this rules out edge 2-3 for us
            return min(point_segment_distance(x0,x1,x2), point_segment_distance(x0,x1,x3))
        elseif w31>0 # this rules out edge 1-3
            return min(point_segment_distance(x0,x1,x2), point_segment_distance(x0,x2,x3))
        else # w12 must be >0, ruling out edge 1-2
            return min(point_segment_distance(x0,x1,x3), point_segment_distance(x0,x2,x3))
        end
    end
end





# http://fileadmin.cs.lth.se/cs/personal/tomas_akenine-moller/code/raytri_tam.pdf
function triangle_ray_intersection(orig, dir, face)

    # find vectors for two edges sharing vert0 
    edge1 = face.points[2] - face.points[1]
    edge2 = face.points[3] - face.points[1]

    # begin calculating determinant - also used to calculate U parameter
    pvec = cross(dir, edge2)

    # if determinant is near zero, ray lies in plane of triangle
    det = dot(edge1, pvec)

    if det ≈ 0
        return false, zero(det)
    end

    inv_det = 1/det

    # calculate distance from vert0 to ray origin
    tvec = orig - vert0

    # calculate U parameter and test bounds 
    u = dot(tvec, pvec) * inv_det
    if u < 0 || u > 1
        return false, zero(det)
    end

    # prepare to test V parameter
    qvec = cross(tvec, edge1)

    # calculate V parameter and test bounds
    v = dot(dir, qvec) * inv_det
    if v < 0 || u + v > 1
        return false, zero(det)
    end

    # calculate t, ray intersects triangle
    t = dot(edge2, qvec) * inv_det

    return true, t
end
    


function signed_distance(point, trimesh::TriangleMesh)

    @assert length(point)==3

    num_intersect = 0
    distance = 1/eps(eltype(point))

    for face in trimesh

        distance = min(distance, point_triangle_distance(point, face))

        intersect, dist = triangle_ray_intersection(point, -point, face)

        if intersect && dist >= 0
            num_intersect += 1
        end
        
    end

    if num_intersect%2 == 1 
        return -distance
    else
        return distance
    end

end



end
