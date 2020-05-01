import numpy as np
import time
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


# Defines the Vertex class
class Vertex:
    def __init__(self, x, y, nextV=None, prevV=None, intersect=False, inside=False, inbound=False):
        self.x = x
        self.y = y
        self.nextV = nextV
        self.prevV = prevV
        self.intersect = intersect
        self.inside = inside
        self.inbound = inbound
        
    def link(self, v):
        self.neighbor = v
        v.neighbor = self

# Defines the Polygon_ class 
class Polygon_:
    def __init__(self, listPoints, closed):
        self.closed = closed
        self.listPoints = listPoints
        self.vertices = []
        
        # Set vertices if necessary
        if len(self.listPoints) > 0:
            for x, y in listPoints:
                self.vertices.append(Vertex(x, y))
            
        self.set_next_prev()
       
    def add(self, v):
        if (v.x, v.y) in self.listPoints:
            self.closed = True
        else:
            self.listPoints.append((v.x, v.y))
            self.vertices.append(v)
            self.set_next_prev()
    
    def sort(self, v, V0, V1):
        V1_index = self.vertices.index(V1)
        
        self.listPoints.insert(V1_index, (v.x, v.y))
        self.vertices.insert(V1_index, v)
            
        self.set_next_prev()
    
    def set_next_prev(self):
        if len(self.listPoints) > 0:
            self.vertices[0].prevV = self.vertices[len(self.vertices) - 1]
            self.vertices[len(self.vertices) - 1].nextV = self.vertices[0]
            
            # Set nextV and prevV for all vertices
            for i in range(len(self.vertices)):
                v = self.vertices[i]
                if i > 0: v.prevV = self.vertices[i - 1]
                if i < len(self.vertices) - 1: v.nextV = self.vertices[i + 1]
    
    def set_interior(self):
        for v in self.vertices:
            # Horizontal edges are considered both left and right
            if v.nextV.y == v.y: v.interior = 'both'
            
            v.interior = 'left'
            v.interior = 'right'
    
    def inside_or_on(self, v):
        P = Polygon(self.listPoints)
        p = Point(v.x, v.y)
        if P.contains(p) or P.touches(p): return True
        else: return False
    
    def intersect_vertex(self, processed):
        for v, p in zip(self.vertices, processed):
            if v.intersect and not p: return v
        return None
                
# Determines if line segment S0-S1 intersects with C0-C1
def intersect(S0, S1, C0, C1):    
    SdeltaX = float(S1.x - S0.x)
    SdeltaY = float(S1.y - S0.y)
    
    CdeltaX = float(C1.x - C0.x)
    CdeltaY = float(C1.y - C0.y)
    
    den = CdeltaY * SdeltaX - CdeltaX * SdeltaY
    # Slopes are equal
    if den == 0: return None
    
    alphaS = (CdeltaX * (S0.y - C0.y) - CdeltaY * (S0.x - C0.x)) / den
    alphaC = (SdeltaX * (S0.y - C0.y) - SdeltaY * (S0.x - C0.x)) / den
    
    if (0 <= alphaS <= 1) and (0 <= alphaC <= 1):
        x = S0.x + alphaS * SdeltaX
        y = S0.y + alphaS * SdeltaY
        return [(x, y), alphaS, alphaC]
    else:
        return None
  
    
    
    
# weiler_atherton
#   DESCRIPTION:    Implements the Weiler_Atherton polygon clipping algorithm
#   INPUTS:         
#                   C: Clip polygon represented as a Polygon_ object
#                   S: Subject polygon represented as a Polygon_ object
def weiler_atherton(C, S):
    # Loop through S vertices, label as inside or outside clipping polygon
    for v in S.vertices:
        if C.inside_or_on(v):
            v.inside = True
    
    # Find all intersections and insert them into both lists with a link
    #   (same as Greiner-Hormann phase 1)
    S_list = np.copy(S.vertices)
    C_list = np.copy(C.vertices)
    
    num_intersections = 0
    for S0 in S_list:
        S1 = S0.nextV
        
        for C0 in C_list:
            C1 = C0.nextV
            
            intersect_out = intersect(S0, S1, C0, C1)
            
            if intersect_out != None:
                (x, y), alphaS, alphaC = intersect_out
                I1 = Vertex(x, y, intersect=True)
                I2 = Vertex(x, y, intersect=True)
                I1.link(I2)
                S.sort(I1, S0, S1)
                C.sort(I2, C0, C1)
                num_intersections += 1
    
    # There are only three possibilities if there are no intersections
    if num_intersections == 0:
        # S is inside C
        if C.inside_or_on(S.vertices[0]):
            return [S]
        # C is inside S
        elif S.inside_or_on(C.vertices[0]):
            return [C]
        else: return []
    
    # Generate a list of inbound intersections
    for v in S.vertices:
        if v.intersect and v.nextV.inside:
                v.inbound = True
    
    # Like phase three of Greiner-Hormann
    clipped_polygons = []
    processed = [False] * len(S.vertices)
    current = S.intersect_vertex(processed)
    
    # While unprocessed intersecting points in S
    while current != None:
        processed[S.vertices.index(current)] = True
        new_polygon = Polygon_([], False)
        new_polygon.add(Vertex(current.x, current.y))
        
        # Until polygon closed
        while True:
            # Check if inbound
            if current.inbound:
                # Until we reach another intersection
                while True:
                    current = current.nextV
                    p = (current.x, current.y)
                    new_polygon.add(Vertex(p[0], p[1]))
                    if p in S.listPoints:
                        processed[S.listPoints.index(p)] = True
                    if current.intersect: break
            else:
                # Until current.intersect
                while True:
                    current = current.prevV
                    p = (current.x, current.y)
                    new_polygon.add(Vertex(p[0], p[1]))
                    if p in S.listPoints:
                        processed[S.listPoints.index(p)] = True
                    if current.intersect: break
                   
            # Hop to other polygon
            current = current.neighbor
            if new_polygon.closed: break
        
        clipped_polygons.append(new_polygon)
        current = S.intersect_vertex(processed)
        if current == None: break
            
    return clipped_polygons


if __name__=="__main__":
    
    
    times = []
    for i in range(5000):
        # Initialize the clip polygon
        C = Polygon_([(-2, 1), (-5, 6), (0, 1)], True)
        
        # Initialize the subject polygon
        S = Polygon_([(10, 10), (-5, 10), (-10, 5), (-10, -10)], True)
    
        start = time.time()
        # Run Weiler-Atherton polygon clipping algorithm
        R_g = weiler_atherton(C, S)
        times.append(time.time() - start)
        
    print("Average time=", sum(times) / len(times))
    
    # Display results
    C_p = Polygon(C.listPoints)
    plt.plot(*C_p.exterior.xy, color='red', label='Clip Polygon')
    
    S_p = Polygon(S.listPoints)
    plt.plot(*S_p.exterior.xy, color='blue', label='Subject Polygon')
    
    for R_i in R_g:
        R_p = Polygon(R_i.listPoints)
        print(R_i.listPoints)
        plt.plot(*R_p.exterior.xy, color='black', label='Resulting Polygon')
        
    plt.title('Results of Weiler-Atherton')
    plt.legend()
    plt.show()
    
    