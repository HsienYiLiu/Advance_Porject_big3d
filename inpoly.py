import numpy as np
from astropy.io.fits import file

def sign(x):
    """
    Return 1 if x is positive, -1 if it's negative, and 0 if it's zero.
    """
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0


def vertex_sign(P, O):
    """
    Sign of the vertex P with respect to O, as defined above.
    """
    #print("0   ", sign(P[0] - O[0]))
    #print("1   ", sign(P[1] - O[1]))
    #print("2   ", sign(P[2] - O[2]))
    result = sign(P[0] - O[0]) or sign(P[1] - O[1]) or sign(P[2] - O[2])
    print(result)
    
    if not result:
       raise ValueError("vertex coincides with origin")
    return result


def edge_sign(P, Q, O):
    """
    Sign of the edge PQ with respect to O, as defined above.
    """
    result = (
        sign((P[1] - O[1]) * (Q[0] - O[0]) - (P[0] - O[0]) * (Q[1] - O[1])) or
        sign((P[2] - O[2]) * (Q[0] - O[0]) - (P[0] - O[0]) * (Q[2] - O[2])) or
        sign((P[2] - O[2]) * (Q[1] - O[1]) - (P[1] - O[1]) * (Q[2] - O[2]))
    )
    #if not result:
        # with origin")
    return result


def triangle_sign(P, Q, R, O):
    """
    Sign of the triangle PQR with respect to O, as defined above.
    """
    m1_0 = P[0] - O[0]
    m1_1 = P[1] - O[1]
    m2_0 = Q[0] - O[0]
    m2_1 = Q[1] - O[1]
    m3_0 = R[0] - O[0]
    m3_1 = R[1] - O[1]
    result = sign(
        (m1_0 * m2_1 - m1_1 * m2_0) * (R[2] - O[2]) +
        (m2_0 * m3_1 - m2_1 * m3_0) * (P[2] - O[2]) +
        (m3_0 * m1_1 - m3_1 * m1_0) * (Q[2] - O[2]))
    if not result:
        raise ValueError("vertices coplanar with origin")
    return result


def triangle_chain(v1, v2, v3, origin):
    """
    Return the contribution of this triangle to the winding number.
    Raise ValueError if the face contains the origin.
    """
    v1sign = vertex_sign(v1, origin)
    v2sign = vertex_sign(v2, origin)
    v3sign = vertex_sign(v3, origin)

    face_boundary = 0
    if v1sign != v2sign:
        face_boundary += edge_sign(v1, v2, origin)
    if v2sign != v3sign:
        face_boundary += edge_sign(v2, v3, origin)
    if v3sign != v1sign:
        face_boundary += edge_sign(v3, v1, origin)
    if not face_boundary:
        return 0

    return triangle_sign(v1, v2, v3, origin)

def winding_number(point):
            x, y, z = point
            if -1 < x < 1 and -1 < y < 1 and -1 < z < 1:
                return "inside"
            if -1 <= x <= 1 and -1 <= y <= 1 and -1 <= z <= 1:
                return "boundary"
            return "outside"

class Polyhedron(object):
    def __init__(self, triangles, vertex_positions):
        """
        Initialize from list of triangles and vertex positions.
        """
        # Validate: check the combinatorial data.
        edges = set()
        vertices = set()
        for triangle in triangles:
            vertices.update(triangle)
            P, Q, R = triangle
            for edge in ((P, Q), (Q, R), (R, P)):
                if edge[0] == edge[1]:
                    raise ValueError("Self edge: {!r}".format(edge))
                if edge in edges:
                    raise ValueError("Duplicate edge: {!r}".format(edge))
                edges.add(edge)

        # For each edge that appears, the reverse edge should also appear.
        for P, Q in edges:
            if not (Q, P) in edges:
                raise ValueError("Unmatched edge: {!r}".format((P, Q)))

        # Vertex set should match indices in vertex_positions.
        if vertices != set(range(len(vertex_positions))):
            raise ValueError("Vertex set doesn't match position indices.")

        # Vertex positions in R^3.
        self.vertex_positions = vertex_positions
        # Indices making up each triangle, counterclockwise
        # around the outside of the face.
        self.triangles = triangles

    def triangle_positions(self):
        """
        Triples of vertex positions.
        """
        for triangle in self.triangles:
            yield tuple(self.vertex_positions[vx] for vx in triangle)

    def volume(self):
        """
        Return the volume of this polyhedron.
        """
        acc = 0
        for p1, p2, p3 in self.triangle_positions():
            # Twice the area of the projection onto the x-y plane.
            det = ((p2[1] - p3[1]) * (p1[0] - p3[0]) -
                   (p2[0] - p3[0]) * (p1[1] - p3[1]))
            # Three times the average height.
            height = p1[2] + p2[2] + p3[2]
            acc += det * height
        return acc / 6.0

    def winding_number(self, point):
        """Determine the winding number of *self* around the given point.
        """
        return sum(
            triangle_chain(v1, v2, v3, point)
            for v1, v2, v3 in self.triangle_positions()) // 2
"""
def read_off(file):
    #if 'OFF' != file.readline().strip():
     #   raise('Not a valid OFF header')
    n_verts, n_faces, n_dontknow = tuple([int(s) for s in file.readline().strip().split(' ')])
    verts = [[float(s) for s in file.readline().strip().split(' ')] for i_vert in range(n_verts)]
    faces = [[int(s) for s in file.readline().strip().split(' ')][1:] for i_face in range(n_faces)]
    return verts, faces
"""

import os
import sys
 
def read_off(file):
    """
    Reads vertices and faces from an off file.
 
    :param file: path to file to read
    :type file: str
    :return: vertices and faces as lists of tuples
    :rtype: [(float)], [(int)]
    """
 
    assert os.path.exists(file)
 
    with open(file, 'r') as fp:
        
        lines = fp.readlines()
        lines = [line.strip() for line in lines]
        assert lines[0] == 'OFF'
 
        parts = lines[1].split(' ')
        assert len(parts) == 3
 
        num_vertices = int(parts[0])
        assert num_vertices > 0
 
        num_faces = int(parts[1])
        assert num_faces > 0
 
        vertices = []
        for i in range(num_vertices):
            
            vertex = lines[2 + i+1].split(' ')
            
            vertex = [float(point) for point in vertex]
            assert len(vertex) == 3
 
            vertices.append(vertex)
 
        faces = []
        for i in range(num_faces):
            
            face = lines[2 + num_vertices + i + 1].split(' ')
            face = [int(index) for index in face]
            
            #print(face[0], "  " , len(face) - 1)
            #assert face[0] == len(face) - 1
            for index in face:
                assert index >= 0 and index < num_vertices
 
            assert len(face) > 1
 
            faces.append(face)
 
        return vertices, faces
if __name__ == "__main__":
    file = '0.off'
    vertices, faces = read_off(file)
    point = [0,0,0.25]
    print(Polyhedron(faces,vertices).winding_number(point))
    #return verts, faces
    		


            