# Copyright (c) 2023-2024, Nick Gerend
# This file is part of the vizmath library, distributed under a Dual License: Non-Commercial Use and Commercial Use. See LICENSE-NC and LICENSE-COM for details.

from math import atan2, sqrt, sin, cos, tan, pi, exp, atan, radians
import scipy.optimize as optimize
import numpy as np
import pandas as pd

#region basic functions

def distance_between_two_points(x1, y1, x2, y2):
    return sqrt((x2-x1)**2+(y2-y1)**2)

def circle(xo, yo, r=1., points=50, spread=360., end_cap=False):
    if points % 2 != 0:
        points += 1
    list_xy = []
    angle = 0. #-90.
    path_i = 1
    points_to_draw = points
    for i in range(points_to_draw):
        list_xy.append((xo+r*sin(angle*pi/180.), yo+r*cos(angle*pi/180.), path_i))
        angle += 1./points*spread
        path_i += 1
    if end_cap:
        list_xy.append(list_xy[0])
    return list_xy

def eq_triangle(xo, yo, side_length=1., area=None, end_cap=False):
    if area is not None:
        side_length = sqrt(4 * area / sqrt(3))
    height = sqrt(3) / 2 * side_length
    top_vertex = (xo, yo + height * 2 / 3)
    right_vertex = (xo + side_length / 2, yo - height / 3)
    left_vertex = (xo - side_length / 2, yo - height / 3)
    list_xy = [top_vertex, right_vertex, left_vertex]
    if end_cap:
        list_xy.append(list_xy[0])
    return list_xy

def regular_polygon(xo, yo, n_sides, area=1, end_cap=False):
    angle = pi / n_sides
    radius = sqrt((2 * area) / (n_sides * sin(2 * angle)))
    list_xy = []
    for k in range(n_sides):
        theta = -2 * pi * k / n_sides + pi / 2
        xk = xo + radius * cos(theta)
        yk = yo + radius * sin(theta)
        list_xy.append((xk, yk))
    if end_cap:
        list_xy.append(list_xy[0])
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

def polarize(x, max_x, y=0., y_offset = 0.):
    angle = (2.*pi)*(((x)%(max_x))/(max_x))
    angle_deg = angle * 180./pi
    angle_rotated = (abs(angle_deg-360.)+90.) % 360. 
    angle_new = angle_rotated * pi/180.
    xp = (y_offset+y)*cos(angle_new)
    yp = (y_offset+y)*sin(angle_new)
    return xp, yp

def point_along_path(x1, y1, x2, y2, d):
    u = distance_between_two_points(x1, y1, x2, y2)
    s = d/u
    x = x1+s*(x2-x1)
    y = y1+s*(y2-y1)
    return x, y

def two_lines_intercept(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y, p3_x, p3_y):
    """
    Compute the intersection point of two lines, if it exists.
    Line 1 is represented by points (p0_x, p0_y) and (p1_x, p1_y).
    Line 2 is represented by points (p2_x, p2_y) and (p3_x, p3_y).
    """
    xdiff = (p0_x - p1_x, p2_x - p3_x)
    ydiff = (p0_y - p1_y, p2_y - p3_y)

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
        return None # lines do not intersect

    d = (det((p0_x, p0_y), (p1_x, p1_y)), det((p2_x, p2_y), (p3_x, p3_y)))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y

def line_to_point_distance(x0, y0, x1, y1, x2, y2):
    # point: x0, y0
    # line: x1, y1, x2, y2
    n = abs((y1-y2)*x0+(x2-x1)*y0+x1*y2-x2*y1)
    d = sqrt((x2-x1)**2+(y2-y1)**2)
    return n/d

def distance_between_two_points_3d(x1, y1, z1, x2, y2, z2):
    return sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

def rotate(x, y, angledeg, x_offset=0, y_offset=0):
    if angledeg == 0:
        xa = x
        ya = y
        xa -= x_offset
        ya -= y_offset
    xa = x*cos(angledeg*pi/180) + y*sin(angledeg*pi/180)
    ya = -x*sin(angledeg*pi/180) + y*cos(angledeg*pi/180)
    xa -= x_offset
    ya -= y_offset
    return xa, ya

def angle_of_two_points(x1, y1, x2, y2):
    return atan2(x2-x1, y2-y1)*180/pi

def point_along_path_3d(x1, y1, z1, x2, y2, z2, d):
    u = distance_between_two_points_3d(x1, y1, z1, x2, y2, z2)
    s = d/u
    x = x1+s*(x2-x1)
    y = y1+s*(y2-y1)
    z = z1+s*(z2-z1)
    return x,y,z

def ellipse(x1, x2, height, points, sign=1.):
    width = abs(x2-x1)
    a = width/2
    b = height/2
    x = np.linspace(-a, a, num=points)
    angle = pi-(x+a)/(width)*pi
    u = np.tan(angle/2.)
    xt = a*(1-u**2)/(u**2+1)
    y = sign*2*b*u/(u**2+1)
    xt += x1 + (x2-x1)/2.
    return list(zip(xt,y))

def radial_projection(x1, y1, length, xo=0, yo=0):
    angle = atan2(y1 - yo, x1 - xo)
    x2 = x1 + length * cos(angle)
    y2 = y1 + length * sin(angle)
    return x2, y2

def is_point_on_line(x0, y0, x1, y1, x2, y2):
    """
    Check if point (x0, y0) is on the line segment between points (x1, y1) and (x2, y2).
    """
    crossproduct = (y0 - y1) * (x2 - x1) - (x0 - x1) * (y2 - y1)
    if abs(crossproduct) > 1e-7:
        return False
    dotproduct = (x0 - x1) * (x2 - x1) + (y0 - y1) * (y2 - y1)
    if dotproduct < 0:
        return False
    squaredlengthba = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)
    if dotproduct > squaredlengthba:
        return False
    return True

#endregion

#region binning functions
def range_group(df, item_col='item', min_col='min', max_col='max', group_name='group'):
    # Nick's range group function
    dfi = df[[item_col, min_col, max_col]]
    records = dfi.to_records(index=False)
    list_c = list(records)
    group = 0
    groups = {}
    groups[group] = [list_c[0]]
    pop_list = []
    minx = list_c[0][1]
    maxx = list_c[0][2]
    list_c.pop(0)
    while len(list_c) > 0:
        for i in range(len(list_c)):
            minx2 = list_c[i][1]
            maxx2 = list_c[i][2]
            c1 = minx <= minx2 and maxx >= minx2
            c2 = minx <= maxx2 and maxx >= maxx2
            c3 = minx2 <= minx and maxx2 >= minx
            c4 = minx2 <= maxx and maxx2 >= maxx
            if c1 or c2 or c3 or c4:
                if minx2 <= minx:
                    minx = minx2
                if maxx2 >= maxx:
                    maxx = maxx2
                groups[group].append(list_c[i])
                pop_list.append(i)
        if len(pop_list) == 0:
            group += 1
            minx = list_c[0][1]
            maxx = list_c[0][2]
            groups[group] = [list_c[0]]
            list_c.pop(0)
        else:
            pop_list = sorted(pop_list, reverse=True)
            for idx in pop_list:
                if idx < len(list_c):
                    list_c.pop(idx)
            pop_list = []

    rows = []
    for key in groups:
        for i in range(len(groups[key])):
            rows.append((key,groups[key][i][0]))

    df_join = pd.DataFrame(rows, columns=[group_name, item_col])
    df_out = df.merge(df_join, on=[item_col], how='left')
    return df_out
#endregion

#region chord functions

def chord(x0, y0, x1, y1, x2, y2, points, h_override=0.):
    # Nick's chord function
    h = line_to_point_distance(x0, y0, x1, y1, x2, y2)
    w = distance_between_two_points(x1, y1, x2, y2)
    if h_override == 0.:
        new_h = (1.-(h/w)/10.)*h
        if new_h < h*0.01:
            h = h*0.01
        else:
            h = new_h
    else:
        h = h*h_override
    a = angle_of_two_points(x1, y1, x2, y2) - 90
    xr = []
    yr = []
    for i in range(points+1):
        arc_percent = i/(points/2.)
        if i > points/2.:
            arc_percent = (points-i)/(points/2.)
        if i == 0 or i == points:
            arc = 0.
        else:
            arc = sqrt((h/2.)**2-((h/2.)-(h/2.)/((points)/2.)*i)**2.)
        percent = arc/(h/2.)
        y_1 = -percent*arc+(1-percent)*arc_percent
        y_2 = percent*arc+(1-percent)*arc_percent
        xr_1, yr_1 = rotate(i/points*w, y_1, a, -x1, -y1)
        xr_2, yr_2 = rotate(i/points*w, y_2, a, -x1, -y1)
        d1 =  distance_between_two_points(x0, y0, xr_1, yr_1)
        d2 =  distance_between_two_points(x0, y0, xr_2, yr_2)
        if d1 < d2:
            xr.append(xr_1)
            yr.append(yr_1)
        else:
            xr.append(xr_2)
            yr.append(yr_2)
    return list(zip(xr, yr))

def __chord_center(x0, y0, z0, x1, y1, z1, x2, y2, z2):
    p_to_p_dist = distance_between_two_points_3d(x1, y1, z1, x2, y2, z2)
    p_to_center_dist = distance_between_two_points_3d(x1, y1, z1, x0, y0, z0)
    a = p_to_p_dist/2
    b=sqrt(p_to_center_dist**2-a**2)
    x1 = -a
    y1 = b
    x2 = a
    y2 = b
    return x1, y1, x2, y2

def chord_3d(x0, y0, z0, x1, y1, z1, x2, y2, z2, points, h_override=0.):
    # Nick's 3d chord function
    x0 *= 100
    y0 *= 100
    z0 *= 100
    x1 *= 100
    y1 *= 100
    z1 *= 100
    x2 *= 100
    y2 *= 100
    z2 *= 100
    xc1, yc1, xc2, yc2 = __chord_center(x0, y0, z0, x1, y1, z1, x2, y2, z2)
    chord_list = chord(0., 0., xc1, yc1, xc2, yc2, points, h_override)
    chord_list_3d = []
    for i in range(len(chord_list)):
        xc = chord_list[i][0]
        yc = chord_list[i][1]
        xp, yp = two_lines_intercept(xc1, yc1, xc2, yc2, 0., 0., xc, yc)
        d1 = distance_between_two_points(xc1, yc1, xp, yp)
        d2 = distance_between_two_points(xp, yp, xc, yc)
        xp1, yp1, zp1 = point_along_path_3d(x1, y1, z1, x2, y2, z2, d1)
        xp2, yp2, zp2 = point_along_path_3d(xp1, yp1, zp1, x0, y0, z0, d2)
        chord_list_3d.append((xp2/100., yp2/100., zp2/100., i))
    return chord_list_3d

#endregion

#region heart functions

def __zheart(x, y, sign):
    s = 1.
    if sign == 'neg':
        s = -1.
    return s * np.sqrt(.5 * (-x**2 - y**2 + 1) + (2**(1/3) * x**3)/(-86400 * x**5 + 1641600 * x**3 * y**2 + 86400 * x**3 + np.sqrt((-86400 * x**5 + 1641600 * x**3 * y**2 + 86400 * x**3)**2 - 55296000 * x**9))**(1/3) + (-86400 * x**5 + 1641600 * x**3 * y**2 + 86400 * x**3 + np.sqrt((-86400 * x**5 + 1641600 * x**3 * y**2 + 86400 * x**3)**2 - 55296000 * x**9))**(1/3)/(240 * 2**(1/3)))

def __zheart_x0(x, y, sign):
    s = 1.
    if sign == 'neg':
        s = -1.
    return s * np.sqrt(1 - y**2)/np.sqrt(2)

def __zheart_y0(x, y, sign):
    s = 1.
    if sign == 'neg':
        s = -1.
    return s * np.sqrt(1/2 * (1 - x**2) + x**3/(4 * 15**(1/3) * (-45 * x**5 + 45 * x**3 + np.sqrt(15) * np.sqrt(135 * x**10 - x**9 - 270 * x**8 + 135 * x**6))**(1/3)) + (-45 * x**5 + 45 * x**3 + np.sqrt(15) * np.sqrt(135 * x**10 - x**9 - 270 * x**8 + 135 * x**6))**(1/3)/(4 * 15**(2/3)))

def heart_3d(N=1.5, n=500):
    x = np.linspace(-N,N,n)
    y = np.linspace(-N,N,n)
    Xgrid1, Ygrid1 = np.meshgrid(x, y)
    Xgrid2, Ygrid2 = np.meshgrid(x, y)
    Zgrid1 = __zheart(Xgrid1, Ygrid1, 'pos')
    Zgrid2 = __zheart(Xgrid2, Ygrid2, 'neg')
    Xout1 = np.reshape(Xgrid1, -1)
    Yout1 = np.reshape(Ygrid1, -1)
    Zout1 = np.reshape(Zgrid1, -1)
    Xout2 = np.reshape(Xgrid2, -1)
    Yout2 = np.reshape(Ygrid2, -1)
    Zout2 = np.reshape(Zgrid2, -1)
    Xout = np.concatenate([Xout1, Xout2])
    Yout = np.concatenate([Yout1, Yout2])
    Zout = np.concatenate([Zout1, Zout2])

    Zgrid1_x0 = __zheart_x0(Xgrid1, Ygrid1, 'pos')
    Zgrid1_y0 = __zheart_y0(Xgrid1, Ygrid1, 'pos')
    Zout1_x0 = np.reshape(Zgrid1_x0, -1)
    Zout1_y0 = np.reshape(Zgrid1_y0, -1)
    Zgrid2_x0 = __zheart_x0(Xgrid2, Ygrid2, 'neg')
    Zgrid2_y0 = __zheart_y0(Xgrid2, Ygrid2, 'neg')
    Zout2_x0 = np.reshape(Zgrid2_x0, -1)
    Zout2_y0 = np.reshape(Zgrid2_y0, -1)

    Zout_x0 = np.concatenate([Zout1_x0, Zout2_x0])
    Zout_y0 = np.concatenate([Zout1_y0, Zout2_y0])

    return zip(Yout,Xout,Zout,Zout_x0,Zout_y0)

#endregion

#region image functions

def laplacian_of_gaussian(x, y, sigma, zscale=100000):
    ratio = (x**2+y**2)/(2*sigma**2)
    return (-1/(np.pi*sigma**4)*(1-ratio)*np.exp(-ratio))*zscale

def laplacian_of_gaussian_grid(N, n, x, y, sigma, zscale=100000):
    N = 49 # odd number for even distribution around 0
    n = 500
    x = np.linspace(0.0,N,n)
    y = np.linspace(0.0,N,n)
    Xgrid, Ygrid = np.meshgrid(x, y)
    Zgrid = -laplacian_of_gaussian(Xgrid-N//2, Ygrid-N//2, sigma=sigma)
    Xout = np.reshape(Xgrid, -1)
    Yout = np.reshape(Ygrid, -1)
    Zout = np.reshape(Zgrid, -1)
    return Xout,Yout,Zout

#endregion

#region artistic functions

def feather(points, sf = 2.0, ef = 2.0, sx = 0.0, ex = 2.0, scale = 25.):
    # Nick's feather function
    xi = np.linspace(1.5, 3.5, num=points)
    x = np.linspace(sx, ex, num=points)
    a = np.linspace(0., 40., num=points)
    f = np.linspace(sf, ef, num=points)
    y = (1./xi**f)*np.sin(a*np.pi/180.)*scale
    return x, y

def __teardrop_angle_function(a, y):
    at = -np.cos(a*np.pi/180,)
    at -= y
    return at

def __teardrop_angle(a1, a2, y, maxiter=1000):
    a = optimize.bisect(__teardrop_angle_function, a1, a2, args=(y), maxiter=maxiter)
    return a

def __teardrop_yz(y, miny, maxy, minz, maxz, m, f=1.):
    zs = rescale(y, miny, maxy, -1., 1.)
    a1 = 0.
    a2 = 180.
    a = __teardrop_angle(a1, a2, zs)
    theta = a*pi/180.
    z = -cos(theta)
    zs =  rescale(z, -1., 1., minz, maxz)
    yt = np.sin(theta)*np.sin(0.5*theta)**m*f
    yts = rescale(yt, 0., 1.0, miny, maxy)
    return yts, zs, a

def teardrop(list_xy, miny, maxy, minz, maxz):
    # Nick's teardrop function
    for i in range(len(list_xy)):
        y = list_xy[i].y
        if i==20:
            stop = 0
        yt, z, alast = __teardrop_yz(y, miny, maxy, minz, maxz, 1.5, f=1.2) # balloon shape: y, 0., 10., 0., 8., 1.5, f=1.2
        list_xy[i].y = yt
        list_xy[i].z = z
    return list_xy

#endregion

#region geospatial functions

def lat_long_to_3d(latitude, longitude, radius, out='xyz'):
    ls = atan((1-0)**2 * tan(radians(latitude)))
    r = sqrt((radius**2)/((1 + (1/((1-0)**2)-1) * (sin(ls)**2))))
    x = r * cos(ls) * cos(radians(longitude)) + 1 * cos(radians(latitude)) * cos(radians(longitude))
    y = r * cos(ls) * sin(radians(longitude)) + 1 * cos(radians(latitude)) * sin(radians(longitude))
    z = r * sin(ls) + 1 * sin(radians(latitude))
    output = x,y,z
    if out=='x':
        output = x
    if out=='y':
        output = y
    if out=='z':
        output = z
    return output

def earth():
    Earth_radius = 6371. # km
    coefs = (1, 1, 1)  
    rx, ry, rz = [Earth_radius/np.sqrt(coef) for coef in coefs]

    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)

    x = rx * np.outer(np.cos(u), np.sin(v))
    y = ry * np.outer(np.sin(u), np.sin(v))
    z = rz * np.outer(np.ones_like(u), np.cos(v))

    x = np.reshape(x, -1)
    y = np.reshape(y, -1)
    z = np.reshape(z, -1)

    return list(zip(x,y,z))

#endregion

#region distribution functions

def sphere(points=1000):
    num_pts = points
    indices = np.arange(0, num_pts, dtype=float) + 0.5
    phi = np.arccos(1 - 2*indices/num_pts)
    theta = np.pi * (1 + 5**0.5) * indices
    x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
    return x,y,z

def __sunflower_radius(i,n,b):
    if i > n-b:
        r = 1
    else:
        r = sqrt(i-1/2)/sqrt(n-(b+1)/2)
    return r

def sunflower(n, alpha, extent):
    list_xy = []
    b = round(alpha*sqrt(n))
    phi = (sqrt(5)+1)/2
    for i in range(n):
        i += 1
        r = __sunflower_radius(i,n,b)
        theta = 2*pi*i/phi**2
        x = r*cos(theta)*extent
        y = r*sin(theta)*extent
        list_xy.append([x,y])
    return list_xy

def sunflower_donut(population, count, start, extent):
    xy = sunflower(population, 1., extent)
    xy_donut = xy[start:start+count]
    return xy_donut

def concentric_spread(radius, arc_length, points, xo=0, yo=0, rotate=0, style=''):
    # Nick's concentric spread function
    rr = radians(rotate)
    angles =  [radians(90) - i * arc_length / radius for i in range(points)]
    avg_a = rr
    if style == 'gravity':
        avg_a += sum(angles) / points + pi/2
    elif style == 'antigravity':
        avg_a += sum(angles) / points - pi/2 
    xy_cs = [[xo + radius * cos(a-avg_a), yo + radius * sin(a-avg_a)] for a in angles]
    return xy_cs

#endregion

#region physics functions

def gravity(e_rad, e_mass):
    e_m = 5.96e24
    e_r = 6.37e6
    G = 6.673e-11
    return G*(e_mass*e_m)/((e_rad*e_r)**2)

def orbital(apogee, eccentricity=0, inclination=0, ascension=0, perigee=0):
    inc = inclination * pi / 180.
    M1 = np.matrix([ [0, cos(inc), -sin(inc)], [0, sin(inc), cos(inc)], [1, 0, 0] ])

    rotation = (ascension + perigee) * pi/180
    M2 = np.matrix([ [0, 0, 1], [cos(rotation), -sin(rotation), 0], [sin(rotation), cos(rotation), 0] ])    
    angle = np.linspace(0,2*pi, 182)
    r = (apogee * (1-eccentricity**2)) / (1 + eccentricity*cos(angle))

    xr = r*cos(angle)
    yr = r*sin(angle)
    zr = 0.

    pts = np.matrix(list(zip(xr,yr,zr)))
    pts =  (M1 * M2 * pts.T).T
    xr,yr,zr = pts[:,0].A.flatten(), pts[:,1].A.flatten(), pts[:,2].A.flatten()

    return list(zip(xr,yr,zr))

def effective_solar_flux(teff):
    if teff < 2600 or teff > 7200:
        return 0, 0, 0, 0
    # credit: arXiv:1404.5292
    # Coeffcients to be used in the analytical expression to calculate habitable zone flux boundaries
    seff = [0,0,0,0,0,0]
    seffsun  = [1.776,1.107, 0.356, 0.320, 1.188, 0.99] 
    a = [2.136e-4, 1.332e-4, 6.171e-5, 5.547e-5, 1.433e-4, 1.209e-4]
    b = [2.533e-8, 1.580e-8, 1.698e-9, 1.526e-9, 1.707e-8, 1.404e-8]
    c = [-1.332e-11, -8.308e-12, -3.198e-12, -2.874e-12, -8.968e-12, -7.418e-12]
    d = [-3.097e-15, -1.931e-15, -5.575e-16, -5.011e-16, -2.084e-15, -1.713e-15]
    tstar = teff - 5780.0
    for i in range(len(a)):
        seff[i] = seffsun[i] + a[i]*tstar + b[i]*tstar**2 + c[i]*tstar**3 + d[i]*tstar**4
    recentVenus = seff[0] # optimistic too close
    runawayGreenhouse = seff[1] # conservative too close
    maxGreenhouse = seff[2] # conservative too far
    earlyMars = seff[3] # optimistic too far
    #fivemeRunaway = seff[4] # earth close limit
    #tenthmeRunaway = seff[5] # closer limit 
    #c_close, c_far, o_close, o_far
    return runawayGreenhouse, maxGreenhouse, recentVenus, earlyMars

#endregion

#region collision functions

def circle_collide(x1, y1, x2, r1, r2, place='top'):
    y2 = None
    # if (x1-r1 <= x2-r2 and x1+r1 >= x2-r2) or (x1-r1 <= x2+r2 and x1+r1 >= x2+r2):
    x1_left, x1_right = x1 - r1, x1 + r1
    x2_left, x2_right = x2 - r2, x2 + r2
    if not (x2_right < x1_left or x2_left > x1_right):
        yc = -x2**2 + 2*x2*x1 + r2**2 + 2*r2*r1 + r1**2 - x1**2
        if yc < 0:
            yc = 0
        if place == 'bottom':
            y2 = y1 - sqrt(yc)
        if place == 'top':
            y2 = sqrt(yc) + y1
    return y2

def circle_collided(x1, y1, x2, y2, r1, r2, tol=10, r_tol=0.0):
    r1 = r1-r1*r_tol
    A = (x2-x1)**2 + (y2-y1)**2
    B = (r1+r2)**2
    result = False
    if round(A, tol) < round(B, tol):
        result = True
    return result

#endregion

#region polygon functions

def sort_vertices(polygon):
    centroid = np.mean(polygon, axis=0)
    def sort_key(point):
        return np.arctan2(point[1] - centroid[1], point[0] - centroid[0])
    return sorted(polygon, key=sort_key)

def is_point_in_convex_polygon(x, y, polygon):
        def cross_product(o, a, b):
            return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
        def is_left_turn(p, q, r):
            return cross_product(p, q, r) > 0
        n = len(polygon)
        if n < 3:
            return False  # not a polygon
        sorted_polygon = sort_vertices(polygon)
        # check for each edge if the point is on the same side
        for i in range(n):
            if not is_left_turn(sorted_polygon[i], sorted_polygon[(i + 1) % n], (x, y)):
                return False
        return True

def line_polygon_intercepts(x1, y1, x2, y2, polygon):
    """
    Find intersection points between a line (represented by points (x1, y1) and (x2, y2))
    and a polygon (represented as a list of (x, y) points).
    """
    sorted_polygon = sort_vertices(polygon)
    intersections = []
    for i in range(len(sorted_polygon)):
        p1 = sorted_polygon[i]
        p2 = sorted_polygon[(i + 1) % len(sorted_polygon)]
        try:
            intersect_x, intersect_y = two_lines_intercept(x1, y1, x2, y2, p1[0], p1[1], p2[0], p2[1])
            # check if the intercept is on both the line segment of the polygon edge and the given line
            if is_point_on_line(intersect_x, intersect_y, p1[0], p1[1], p2[0], p2[1]):
                intersections.append((intersect_x, intersect_y))
        except:
            pass  # lines are parallel, no intersection
    return intersections

def rectangle(width, height, rotate=0, area=50, x_offset=0, y_offset=0):
    # adjust width and height according to the specified area
    scale_factor = (area / (width * height)) ** 0.5
    width *= scale_factor
    height *= scale_factor
    # rectangle vertices before rotation
    rectangle = np.array([
        [-width / 2, -height / 2],
        [-width / 2, height / 2],
        [width / 2, height / 2],
        [width / 2, -height / 2]
    ])
    rectangle += np.array([x_offset, y_offset])
    if rotate == 0:
        return rectangle
    else:
        # rotation matrix
        theta = np.radians(rotate)
        rotation_matrix = np.array([
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)]
        ])
        # rotated rectangle
        rotated_rectangle = np.dot(rectangle, rotation_matrix)
        return rotated_rectangle

def rotate_polygon(polygon, rotate):
    angle_radians = radians(rotate)
    rotated_polygon = []
    for x, y in polygon:
        # rotate point
        rotated_x = x * cos(angle_radians) - y * sin(angle_radians)
        rotated_y = x * sin(angle_radians) + y * cos(angle_radians)
        rotated_polygon.append((rotated_x, rotated_y))
    return rotated_polygon

#endregion