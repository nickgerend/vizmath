# vizmath
# Author: Nick Gerend

from math import sin, cos, pi, exp, atan2, sqrt
    
def angle_by_two_points(x1, y1, x2, y2):
    a = atan2(x2-x1, y2-y1)*180/pi
    if a < 0:
        a += 360. 
    return a

def distance_between_two_points(x1, y1, x2, y2):
    return sqrt((x2-x1)**2+(y2-y1)**2)

def circle(xo, yo, r=1., points=50, spread=360., end_cap=False):
    if points % 2 != 0:
        points += 1
    list_xy = []
    angle = 0. #-90.
    path_i = 1
    points_to_draw = points
    if end_cap:
        points_to_draw += 1
    for i in range(points_to_draw):
        list_xy.append((xo+r*sin(angle*pi/180.), yo+r*cos(angle*pi/180.), path_i))
        angle += 1./points*spread
        path_i += 1
    return list_xy

def rescale(x, xmin, xmax, newmin, newmax):
    rescaled = (newmax-newmin)*((x-xmin)/(xmax-xmin))+newmin
    return rescaled

def sigmoid(x1, y1, x2, y2, points, orientation = 'h', limit = 6):
    x_1 = x1
    y_1 = y1
    x_2 = x2
    y_2 = y2
    if orientation == 'v':
        x1 = y_1
        y1 = x_1
        x2 = y_2
        y2 = x_2
    x = []
    y = []
    amin = 1./(1.+exp(limit))
    amax = 1./(1.+exp(-limit))
    da = amax-amin
    for i in range(points):
        i += 1
        xi = (i-1.)*((2.*limit)/(points-1.))-limit
        yi = ((1.0/(1.0+exp(-xi)))-amin)/da
        x.append((xi-(-limit))/(2.*limit)*(x2-x1)+x1)
        y.append((yi-(0.))/(1.)*(y2-y1)+y1)
    return { 'h': list(zip(x,y)), 'v': list(zip(y,x))}.get(orientation, None)

def polarize(x, max_x, y=0.):
    angle = (2.*pi)*(((x)%(max_x))/(max_x))
    angle_deg = angle * 180./pi
    angle_rotated = (abs(angle_deg-360.)+90.) % 360. 
    angle_new = angle_rotated * pi/180.
    xp = (y)*cos(angle_new)
    yp = (y)*sin(angle_new)
    return xp, yp

def point_along_path(x1, y1, x2, y2, d):
    u = distance_between_two_points(x1, y1, x2, y2)
    s = d/u
    x = x1+s*(x2-x1)
    y = y1+s*(y2-y1)
    return x, y

def two_lines_intercept(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y, p3_x, p3_y):
    s1_x = p1_x - p0_x     
    s1_y = p1_y - p0_y
    s2_x = p3_x - p2_x
    s2_y = p3_y - p2_y
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y)
    i_x = p0_x + (t * s1_x)
    i_y = p0_y + (t * s1_y)
    return i_x, i_y