import time
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon


# Return True or False depending on whether the point p is on the interior side of the edge (C0, C1)
def inside(p, C0, C1):
    return (C1[0]-C0[0])*(p[1]-C0[1]) > (C1[1]-C0[1])*(p[0]-C0[0])

# Returns intersection between S0, S1, C0, and C1
def intersect(C0, C1, S0, S1):
    SdeltaX = S0[0] - S1[0]
    SdeltaY = S0[1] - S1[1]
    
    CdeltaX = C0[0] - C1[0]
    CdeltaY = C0[1] - C1[1]
    
    n1 = C0[0] * C1[1] - C0[1] * C1[0]
    n2 = S0[0] * S1[1] - S0[1] * S1[0] 
    n3 = 1.0 / (CdeltaX * SdeltaY - CdeltaY * SdeltaX)
    
    # Return intersection point
    return ((n1*SdeltaX - n2*CdeltaX) * n3, (n1*SdeltaY - n2*CdeltaY) * n3)
    
# sutherland_hodgman
#   DESCRIPTION:    Implements the Sutherland_Hodgman polygon clipping algorithm
#   INPUTS:         
#                   C: Clip polygon represented as a Polygon_ object
#                   S: Subject polygon represented as a Polygon_ object
def sutherland_hodgman(S, C):
    clippedPolygon = S
    C0 = C[-1]
 
    # For each edge in clip polygon
    for c in C:
        C1 = c
        inputList = clippedPolygon
        clippedPolygon = []
        
        # S0 is previous point
        if len(inputList) > 0:
            S0 = inputList[-1]
        
        # For each point in clipped subject polygon
        for s in inputList:
            S1 = s
            # S1 is inside clip polygon, add S1
            if inside(S1, C0, C1):
                
                # S0 is outside clip polygon, add intersection
                if not inside(S0, C0, C1):
                    clippedPolygon.append(intersect(C0, C1, S0, S1))
                    
                # Add S1
                clippedPolygon.append(S1)
                
            # S0 is inside clip polygon, add intersection
            elif inside(S0, C0, C1):
                clippedPolygon.append(intersect(C0, C1, S0, S1))
            S0 = S1
        C0 = C1
        
    # Return clipped polygon
    return clippedPolygon


if __name__=="__main__":    
    times = []
    for i in range(5000):
        # Initialize the clip polygon
        C = [(-2, 1), (-5, 6), (0, 1)]
        
        # Initialize the subject polygon
        S = [(10, 10), (-5, 10), (-10, 5), (-10, -10)]
        
        start = time.time()
        # Run Sutherland-Hodgman clipping algorithm
        R_s = sutherland_hodgman(C, S)
        times.append(time.time() - start)
        
    print("Average time=", sum(times) / len(times))
    
    C_p = Polygon(C)
    plt.plot(*C_p.exterior.xy, color='red', label='Clip Polygon')
    
    S_p = Polygon(S)
    plt.plot(*S_p.exterior.xy, color='blue', label='Subject Polygon')
    
    print(R_s)
    if len(R_s) > 0:
        R_p = Polygon(R_s)
        plt.plot(*R_p.exterior.xy, color='black', label='Resulting Polygon')
        
    plt.title('Results of Sutherland-Hodgman')
    plt.legend()
    plt.show()
    
    