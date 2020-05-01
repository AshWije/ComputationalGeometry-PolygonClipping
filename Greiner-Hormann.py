import numpy as np
import time
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon



# Defines the Vertex class
class Vertex:
    def __init__(self, x, y, nextV=None, prevV=None, intersect=False, entry_exit=False, neighbor=None, alpha=0.0):
        self.x = x
        self.y = y
        self.nextV = nextV
        self.prevV = prevV
        self.intersect = intersect
        self.entry_exit = entry_exit
        self.neighbor = neighbor
        self.alpha = alpha
        
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
    
    def inside(self, v):
        if Polygon(self.listPoints).contains(Point(v.x, v.y)):
            return True
        else: return False
    
    def intersect_vertex(self, processed):
        for v, p in zip(self.vertices, processed):
            if v.intersect and not p: return v
        return None
                
        

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
  
    
    
    
# greiner_hormann
#   DESCRIPTION:    Implements the Greiner-Hormann polygon clipping algorithm
#   INPUTS:         
#                   C: Clip polygon represented as a Polygon_ object
#                   S: Subject polygon represented as a Polygon_ object
def greiner_hormann(C, S):
    # Phase one
    S_list = np.copy(S.vertices)
    C_list = np.copy(C.vertices)
    for S0 in S_list:
        S1 = S0.nextV
        
        for C0 in C_list:
            C1 = C0.nextV
            
            intersect_out = intersect(S0, S1, C0, C1)
            
            if intersect_out != None:
                (x, y), alphaS, alphaC = intersect_out
                I1 = Vertex(x, y, intersect=True, alpha=alphaS)
                I2 = Vertex(x, y, intersect=True, alpha=alphaC)
                I1.link(I2)
                S.sort(I1, S0, S1)
                C.sort(I2, C0, C1)
    
    # Phase two
    # Set exit_entry for C intersections
    if S.inside(C.vertices[0]): status = 'exit'
    else: status = 'entry'
    
    for Ci in C.vertices:
        if Ci.intersect:
            Ci.entry_exit = status
            
            # Toggle status
            if status == 'exit': status = 'entry'
            else: status = 'exit'
    
    
    # Set exit_entry for S intersections
    if C.inside(S.vertices[0]): status = 'exit'
    else: status = 'entry'
    
    for Si in S.vertices:
        if Si.intersect:
            Si.entry_exit = status
            
            # Toggle status
            if status == 'exit': status = 'entry'
            else: status = 'exit'
            
    # Phase three
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
            if current.entry_exit == 'entry':
                
                # Until current.intersect
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
    for i in range(1):
        # Initialize the clip polygon
        C = Polygon_([(-2, -10), (-5, 6), (0, 1)], True)
        
        # Initialize the subject polygon
        S = Polygon_([(10, 10), (-5, 10), (-10, 5), (-10, -10)], True)
    
        start = time.time()
        # Run Greiner-Hormann polygon clipping algorithm
        R_g = greiner_hormann(C, S)
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
        
    plt.title('Results of Greiner-Hormann')
    plt.legend()
    plt.show()
    
    