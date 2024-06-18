# Copyright (c) 2023-2024, Nick Gerend
# This file is part of the vizmath library, distributed under a Dual License: Non-Commercial Use and Commercial Use. See LICENSE-NC and LICENSE-COM for details.

# Nick's Path Swarm Chart Algorithm
#%%
import pandas as pd
from math import pi, radians, cos, sin
import numpy as np
import random
import matplotlib.pyplot as plt
import string
from scipy.interpolate import CubicSpline
from matplotlib.patches import Circle
# from IPython.display import display # for continued plotting

from . import functions as vf
from .draw import points as dp

# import functions as vf
# from draw import points as dp

#%%

class pathswarm:
    
    def __init__(self, df, id_field, position_field, 
        min=None, max=None, path=[],
        size_field=None, order_field=None, direction_field=None,
        buffer=0., size_override=None, direction_override=None,
        rotation=0, tol_overlap=10, tol_r=1e-10, interp='cubic_spline', df_only=False, 
        draw_points=True, kwargs={}):

        self.df = df
        self.id_field = id_field
        self.position_field = position_field
        self.min = min
        self.max = max
        self.path = path
        self.size_field = size_field
        self.order_field = order_field
        self.direction_field = direction_field
        self.buffer = buffer
        self.size_override = size_override
        self.direction_override = direction_override
        self.rotation = rotation
        self.tol_overlap = tol_overlap
        self.tol_r = tol_r
        self.interp = interp
        self.df_only = df_only
        self.draw_points = draw_points
        self.kwargs = kwargs
        
        self.i_path = None
        self.i_path_length = None
        self.i_path_progress = None

        self.nodes = None
        self.o_pathswarm = dp()

        if not df_only:
            self.__path_generate()
            self.swarm_setup()
            self.path_swarm()

    @classmethod
    def __generate_path(cls):
        num_points = random.randint(5, 10)
        x_values = sorted(np.random.rand(num_points) * 100)
        y_values = np.random.rand(num_points) * 20
        return list(zip(x_values, y_values))

    @classmethod
    def random_pathswarm(cls, size, min=None, max=None,
        buffer=0., size_override=None, direction_override=None,
        rotation=0, tol_overlap=10, interp='cubic_spline', kwargs={}):
        data = [[''.join(random.choices(string.ascii_letters, k=5)), random.randint(1, 100)/50, random.randint(1, 100)/50] for _ in range(size)]
        df = pd.DataFrame(data, columns=['id', 'position', 'size'])
        path = cls.__generate_path()
        return cls(df, 'id', 'position', min=min, max=max, path=path,
        size_field='size', order_field=None, direction_field=None,
        buffer=buffer, size_override=size_override, direction_override=direction_override,
        rotation=rotation, tol_overlap=tol_overlap, interp=interp, df_only=False, kwargs=kwargs)
    
    @classmethod
    def random_path(cls):
        return cls.__generate_path()

    def __path_x_intersections(self, line_points, path_points, known_point=None, tol=1e-5):

        def calc_line_eq(point1, point2):
            slope = (point2[1] - point1[1]) / (point2[0] - point1[0]) if point2[0] - point1[0] != 0 else float('inf')
            intercept = point1[1] - slope * point1[0] if slope != float('inf') else point1[0]
            return slope, intercept
        
        def calc_intersect(line_slope, line_intercept, segment_point1, segment_point2):
            segment_slope, segment_intercept = calc_line_eq(segment_point1, segment_point2)
            
            if line_slope == segment_slope:  # Parallel or identical
                return None
            
            if line_slope == float('inf'):  # Line is vertical
                x = line_intercept
                y = segment_slope * x + segment_intercept
            elif segment_slope == float('inf'):  # Segment is vertical
                x = segment_intercept
                y = line_slope * x + line_intercept
            else:
                x = (segment_intercept - line_intercept) / (line_slope - segment_slope)
                y = line_slope * x + line_intercept
            
            # Check if the intersection is within the segment bounds
            if min(segment_point1[0], segment_point2[0]) <= x <= max(segment_point1[0], segment_point2[0]) and \
            min(segment_point1[1], segment_point2[1]) <= y <= max(segment_point1[1], segment_point2[1]):
                return x, y
            return None

        def remove_known_intersect(known_point, intersection_list, tolerance=1e-5):
            return [p for p in intersection_list if ((p[0] - known_point[0])**2 + (p[1] - known_point[1])**2)**0.5 > tolerance] or None

        intersections = []
        line_slope, line_intercept = calc_line_eq(line_points[0], line_points[1])
        for i in range(len(path_points) - 1):
            intersection = calc_intersect(line_slope, line_intercept, path_points[i], path_points[i + 1])
            if intersection:
                intersections.append(intersection)

        if known_point is not None:
            intersections = remove_known_intersect(known_point, intersections, tol)
            
        return intersections
    
    def __path_x_bounds(self, x_point, p1, p2, x_intersects):
        
        def distance(point1, point2):
            return ((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)**0.5
        
        A = p2[1] - p1[1]
        B = p1[0] - p2[0]
        C = p2[0]*p1[1] - p1[0]*p2[1]
        g1 = []
        g2 = []
        for point in x_intersects:
            x, y = point
            v = A*x + B*y + C
            if v >= 0: # one side or on line
                g1.append(point)
            else: # other side
                g2.append(point)
        b1 = min(g1, key=lambda point: distance(point, x_point), default=None)
        b2 = min(g2, key=lambda point: distance(point, x_point), default=None)
        return b1, b2
    
    def __path_cubic_spline(self, points, closed=False, res=int(1e4)):
        bc_type = 'natural'
        if points[0][0] == points[-1][0] and points[0][1] == points[-1][1] and closed:
            bc_type = 'periodic'
        points_np = np.array(points)
        x, y = points_np[:, 0], points_np[:, 1]
        units = np.linspace(0, 1, len(points_np))
        cs_x = CubicSpline(units, x, bc_type=bc_type)
        cs_y = CubicSpline(units, y, bc_type=bc_type)
        cs_res = np.linspace(0, 1, res)
        cs_xy = np.vstack((cs_x(cs_res), cs_y(cs_res))).T
        cumdist = np.insert(np.cumsum(np.sqrt(np.sum(np.diff(cs_xy, axis=0)**2, axis=1))), 0, 0)
        self.i_path_length = cumdist[-1]
        self.i_path = cs_xy.tolist()
        self.i_path_progress = cumdist

    def __path_linear(self, points):
        points_np = np.array(points)
        cumdist = np.insert(np.cumsum(np.sqrt(np.sum(np.diff(points_np, axis=0)**2, axis=1))), 0, 0)
        self.i_path_length = cumdist[-1]
        self.i_path = points
        self.i_path_progress = cumdist

    def __path_generate(self):
        
        if not self.path:
            self.path = [(self.df[self.position_field].min(),0.),(self.df[self.position_field].max(),0.)]
            self.interp = 'linear'

        if self.interp =='cubic_spline':
            defaults = {'closed': False, 'res':int(1e5)}
            params = {**defaults, **self.kwargs}
            closed = params['closed']
            res = params['res']
            self.__path_cubic_spline(self.path, closed, res)
        elif self.interp =='linear':
            self.__path_linear(self.path)
    
    class __node:

        def __init__(self, id, path_x, path_y, path_xo, path_yo, path_fraction, path_value,
            node_vector,node_vector_angle, node_intersects, node_bounds, 
            node_radius, node_area, node_order, node_size, node_size_by, node_radius_buffer):
            
            self.id = id
            self.path_x = path_x
            self.path_y = path_y
            self.path_xo = path_xo
            self.path_yo = path_yo
            self.path_value = path_value
            self.path_fraction = path_fraction

            self.node_vector = node_vector
            self.node_vector_angle = node_vector_angle
            self.node_intersects = node_intersects
            self.node_bounds = node_bounds

            self.node_radius = node_radius
            self.node_area = node_area
            self.node_order = node_order
            self.node_size = node_size
            self.node_size_by = node_size_by
            self.node_radius_buffer = node_radius_buffer

            self.loop_x = None
            self.loop_y = None

            self.node_x = None
            self.node_y = None
            self.node_prospects = []
            self.node_placed = False
            self.node_neighbor = None
            
    def __node_direction(self, x_point, angle, length):
        north_angle = (360 - angle + 90) % 360
        rads = radians(north_angle)
        dx = (length / 2) * cos(rads)
        dy = (length / 2) * sin(rads)
        point1 = (x_point[0] - dx, x_point[1] - dy)
        point2 = (x_point[0] + dx, x_point[1] + dy)
        return point1, point2

    def swarm_setup(self):

        #region initialize

        #region defaults
        defaults = {'closed': False, 'perp_draw_length': 1., 'calc_intersects': False, 
            'calc_bounds': False, 'sort':'desc', 'size_by':'area', 'order_by':'size'}
        params = {**defaults, **self.kwargs} # default if not initialized
        closed = params['closed']
        perp_draw_length = params['perp_draw_length']
        calc_intersects = params['calc_intersects']
        calc_bounds = params['calc_bounds']
        sort = params['sort']
        size_by = params['size_by']
        order_by = params['order_by']
        #endregion

        #region get initial positions
        i_positions = self.df[self.position_field].values
        #endregion

        #region min max clip
        if self.min == None:
            min_pos = min(i_positions)
        else:
            min_pos = self.min
        if self.max == None:
            max_pos = max(i_positions)
        else:
            max_pos = self.max
        df = self.df[(self.df[self.position_field] >= min_pos) & ( self.df[self.position_field] <= max_pos)].copy(deep=True)
        if self.df_only:
            return df
        #endregion

        #region get positions
        ids = df[self.id_field].values
        positions = df[self.position_field].values
        #endregion

        #region get sizing
        size = None
        if self.size_field is not None and self.size_override is None:
            sizes = df[self.size_field].values
        elif self.size_override is not None:
            size = self.size_override 
            sizes = [size for _ in ids]
        else:
            size = (max_pos-min_pos)/len(positions)*10
            sizes = [size for _ in ids]
        #endregion

        #region get directions
        if self.direction_field:
            angles = df[self.direction_field].values
        #endregion

        #region get ordering
        if self.order_field is not None:
            orders = df[self.order_field].values
        else:
            sort_by = None
            if order_by == 'size':
                sort_by = sizes
            if order_by == 'position':
                sort_by = positions
            sb_set = set(sort_by)
            if sort == 'asc':
                unique_sorted_numbers = sorted(sb_set, reverse=False)
            elif sort == 'desc':
                unique_sorted_numbers = sorted(sb_set, reverse=True)
            elif sort == 'random':
                unique_sorted_numbers = random.sample(list(sb_set), len(sb_set))
            rank_dict = {num: idx + 1 for idx, num in enumerate(unique_sorted_numbers)}
            orders = [rank_dict[num] for num in sort_by]
        #endregion

        #region setup path details
        fractions = [(x - min_pos) / (max_pos - min_pos) for x in positions]
        points_np = np.array(self.i_path)
        total_length = self.i_path_length
        cumdist = self.i_path_progress
        #endregion

        #endregion

        nodes = []

        for i in range(len(ids)):

            #region find item's point along path along with before and after points
            fraction = fractions[i]
            target_length = fraction * self.i_path_length
            if closed and (target_length == 0 or target_length == total_length): # closed loop starting or ending point
                x_point = points_np[0]
                point_before = points_np[-2]
                point_after = points_np[1]
            elif target_length == 0: # starting point
                x_point = points_np[0]
                point_before = points_np[0]
                point_after = points_np[1]
            elif target_length == total_length: #ending point
                x_point = points_np[-1]
                point_before = points_np[-2]
                point_after = points_np[-1]
            else: # in between starting and ending points
                for j in range(len(cumdist) - 1):
                    if cumdist[j] <= target_length <= cumdist[j + 1]:
                        if cumdist[j] == target_length: # is at a point
                            x_point = points_np[j]
                            point_before = points_np[j-1]
                            point_after = points_np[j+1]
                            break
                        elif cumdist[j+1] == target_length: # is at the next point
                            x_point = points_np[j+1]
                            point_before = points_np[j]
                            point_after = points_np[j+2]
                            break
                        else: # is in between points
                            seg_adj = (target_length - cumdist[j]) / (cumdist[j + 1] - cumdist[j])
                            x_point = (1 - seg_adj) * points_np[j] + seg_adj * points_np[j + 1]
                            point_before = points_np[j]
                            point_after = points_np[j + 1]
                            break
            #endregion -> x_point (point on path), point_before, point_after

            #region get tangent
            if self.direction_override is not None:
                x_start, x_end = self.__node_direction(x_point, self.direction_override, perp_draw_length)
            elif self.direction_field:
                x_start, x_end = self.__node_direction(x_point, angles[i], perp_draw_length)
            else:            
                tangent_vector = point_after - point_before
                tangent_vector_normalized = tangent_vector / np.linalg.norm(tangent_vector)
                perpendicular_vector = np.array([-tangent_vector_normalized[1], tangent_vector_normalized[0]])
                x_start = x_point - perpendicular_vector * (perp_draw_length / 2)
                x_end = x_point + perpendicular_vector * (perp_draw_length / 2)
            #endregion -> x_start, x_end

            #region get tangent angle
            x1, y1 = x_start
            x2, y2 = x_end
            node_vector_angle = vf.angle_of_two_points(x1,y1,x2,y2)
            #endregion

            #region find intersections
            x_intersects = []
            if calc_intersects:
                x_intersects = self.__path_x_intersections([x_start,x_end], self.i_path, x_point)
            #endregion -> x_intersects
            
            #region filter intersections to bounds
            b1, b2 = None, None
            if x_intersects and calc_bounds:
                b1, b2 = self.__path_x_bounds(x_point, point_before, point_after, x_intersects)
            #endregion -> b1, b2
                
            #region size
            if size is None:
                i_size = sizes[i]
            else:
                i_size = size
            if size_by == 'radius':
                i_radius = i_size + abs(self.buffer)
                i_area = pi*i_radius**2
            if size_by == 'area':
                i_area = i_size
                i_radius = (i_area/pi)**0.5 + abs(self.buffer)
                i_area = pi*i_radius**2
            #endregion
            
            #region order
            i_order = orders[i]
            #endregion
            
            #region create node
            nodes.append(
                self.__node(id = ids[i],
                path_x = x_point[0], path_y = x_point[1], path_xo = x_point[0], path_yo = x_point[1], 
                path_fraction = fraction, path_value = positions[i],
                node_vector = (x_start, x_end), node_vector_angle = node_vector_angle,
                node_intersects = x_intersects, node_bounds = (b1, b2),
                node_radius = i_radius, node_area = i_area, node_order = i_order,
                node_size = size, node_size_by = size_by, node_radius_buffer = self.buffer
                ))
            #endregion

        self.nodes = nodes

    def path_plot(self, plot=True, index=None):
        path_points = self.i_path
        node_points = [(o.path_x, o.path_y) for o in self.nodes]
        node_vectors = [o.node_vector for o in self.nodes]
        node_vector_angles = [o.node_vector_angle for o in self.nodes]
        node_bounds = [i for o in self.nodes if o.node_bounds is not None for i in o.node_bounds if i is not None]
        if index is not None:
            xya = [(x,y,a) for (x,y), a in zip(node_points, node_vector_angles)]
            rotate = xya[index][2]
            path_points = [vf.rotate(x, y, -rotate) for (x,y) in self.i_path]
            node_points = [vf.rotate(x, y, -rotate) for (x,y,a) in xya]
            node_vectors = [(vf.rotate(x1, y1, -rotate), vf.rotate(x2, y2, -rotate)) for((x1,y1),(x2,y2)) in node_vectors]
            node_bounds = [vf.rotate(x, y, -rotate) for (x,y) in node_bounds]
        fig, ax = plt.subplots()
        path_xs, path_ys = zip(*path_points) 
        ax.plot(path_xs, path_ys, 'k-', label='Path')
        data_xs, data_ys = zip(*node_points)
        ax.plot(data_xs, data_ys, 'ro', label='Data Points')
        for pair in node_vectors:
            xs, ys = zip(*pair)
            ax.plot(xs, ys, 'b-', label='Line Pairs')
        if node_bounds:
            extra_xs, extra_ys = zip(*node_bounds)
            ax.plot(extra_xs, extra_ys, 'bo', label='Extra Points')
        ax.axis('equal')
        if plot:
            plt.show()
        return fig, ax
    
    def node_filter(self, node, type='v'):
        a = -node.node_vector_angle
        x, y = vf.rotate(node.path_x, node.path_y, a)
        r = node.node_radius
        nodes_in_shadow = []
        left_bound = x - r
        right_bound = x + r
        upper_bound = y + r
        lower_bound = y - r
        for n in [o for o in self.nodes if o.node_placed]: # check placed nodes
            x_n, y_n = vf.rotate(n.node_x, n.node_y, a)
            x_left, x_right = x_n - n.node_radius, x_n + n.node_radius
            y_up, y_down = y_n + n.node_radius, y_n - n.node_radius
            if not (x_right < left_bound or x_left > right_bound) and n.id != node.id and type == 'v':
                n.loop_x, n.loop_y = x_n, y_n
                nodes_in_shadow.append(n)
            if not (y_down > upper_bound or y_up < lower_bound) and n.id != node.id and type == 'h':
                n.loop_x, n.loop_y = x_n, y_n
                nodes_in_shadow.append(n)
            if not (x_right < left_bound or x_left > right_bound) and \
                not (y_down > upper_bound or y_up < lower_bound) and \
                n.id != node.id and type == 'b':
                n.loop_x, n.loop_y = x_n, y_n
                nodes_in_shadow.append(n)
        return x, y, nodes_in_shadow
    
    def place_node(self, node, x_rotated, y_rotated, radius, horizon, mode, placed_nodes_in_shadow):
        for nis in placed_nodes_in_shadow:
            x1, y1, r1 = nis.loop_x, nis.loop_y, nis.node_radius
            y_t_rotated = vf.circle_collide(x1, y1, x_rotated, r1, radius, place='top')
            y_b_rotated = vf.circle_collide(x1, y1, x_rotated, r1, radius, place='bottom')
            cb_t = False
            cb_b = False
            # filter for the horizon
            if horizon is not None:
                if horizon == 'top':
                    if y_t_rotated <= y_rotated:
                        cb_t = True
                    if y_b_rotated <= y_rotated:
                        cb_b = True
                elif horizon == 'bottom':
                    if y_t_rotated >= y_rotated:
                        cb_t = True
                    if y_b_rotated >= y_rotated:
                        cb_b = True
            node.node_prospects.append((x_rotated, y_t_rotated, cb_t, nis.id))
            node.node_prospects.append((x_rotated, y_b_rotated, cb_b, nis.id))
        # filter out prospects with collisions
        node.node_prospects = [(x, y,
            not c and not any(vf.circle_collided(x, y, nis.loop_x, nis.loop_y, radius, nis.node_radius, self.tol_overlap, self.tol_r) 
            for nis in placed_nodes_in_shadow),
            i) for x, y, c, i in node.node_prospects]
        if mode == 'closest':
            x_r, y_r, _, nbr = min((p for p in node.node_prospects if p[2]), key=lambda p: abs(p[1] - y_rotated))
        if mode == 'random':
            x_r, y_r, _, nbr = random.choice([p for p in node.node_prospects if p[2]])
        node.node_x, node.node_y = vf.rotate(x_r, y_r, node.node_vector_angle)
        node.node_placed = True
        node.node_neighbor = nbr

    def __set_boundary(self, node, horizon='top', offset=0.):
        a = -node.node_vector_angle
        x, y = vf.rotate(node.path_x, node.path_y, a)
        r = node.node_radius
        yb = y
        if horizon == None:
            yb = y + offset
        elif horizon == 'top':
            yb = y + r + offset
        elif horizon == 'bottom':
            yb = y - r + offset
        xsb, ysb = vf.rotate(x, yb, -a)
        node.path_x, node.path_y = xsb, ysb

    def __swarm_xy(self, circle_resolution=50):
        if self.draw_points:
            p = 0
            for x, y in self.i_path:
                # add path
                self.o_pathswarm.append(p, x, y, p, type='path')
                p += 1
            for n in self.nodes:
                # add circles
                n_circle = vf.circle(n.node_x, n.node_y, r=n.node_radius-n.node_radius_buffer, end_cap=True, points=circle_resolution)
                for c in n_circle:
                    self.o_pathswarm.append(id=n.id, x=c[0], y=c[1], path=c[2], type='node')
                # add path to circle lines
                self.o_pathswarm.append(id=n.id, x=n.path_xo, y=n.path_yo, path=0)
                self.o_pathswarm.append(id=n.id, x=n.node_x, y=n.node_y, path=0)
            
    def swarm_rotate(self, degrees):
        # path
        for i in range(len(self.i_path)):
            self.i_path[i] = vf.rotate(self.i_path[i][0], self.i_path[i][1], degrees)
        # rotate nodes
        for n in self.nodes:
            n.path_x, n.path_y = vf.rotate(n.path_x, n.path_y, degrees)
            n.path_xo, n.path_yo = vf.rotate(n.path_xo, n.path_yo, degrees)
            n.node_x, n.node_y = vf.rotate(n.node_x, n.node_y, degrees)
        # rotate list_xy:
        if self.draw_points:
            for o in self.o_pathswarm.viz:
                o.x, o.y = vf.rotate(o.x, o.y, degrees)

    def path_swarm(self):
        defaults = {'horizon':None, 'offset':0., 'mode':'closest'}
        params = {**defaults, **self.kwargs} # default if not initialized
        horizon = params['horizon']
        offset = params['offset']
        mode = params['mode']
        # adjust the initial points as needed
        if horizon is not None or offset != 0.:
            [self.__set_boundary(node, horizon, offset) for node in self.nodes]
        nodes = sorted(self.nodes, key=lambda o: o.node_order)
        swarm = []
        for node in nodes:
            radius = node.node_radius
            # check the current path location as a first check
            _, _, placed_nodes_in_place = self.node_filter(node, type='b')
            if placed_nodes_in_place:
                x1, y1, r1 = node.path_x, node.path_y, radius
                place = True
                for nip in placed_nodes_in_place:
                    x2, y2, r2 = nip.node_x, nip.node_y, nip.node_radius
                    if vf.circle_collided(x1, y1, x2, y2, r1, r2, self.tol_overlap, self.tol_r):
                        place = False
                if place:
                    #place the node
                    node.node_x, node.node_y, node.node_placed = node.path_x, node.path_y, True
                    swarm.append(node)
                    continue
            else:
                node.node_x, node.node_y, node.node_placed = node.path_x, node.path_y, True
                swarm.append(node)
                continue
            # check all placed node surfaces along tangent for placement (in rotated plane)
            x_rotated, y_rotated, placed_nodes_in_shadow = self.node_filter(node)
            if placed_nodes_in_shadow:
                self.place_node(node, x_rotated, y_rotated, radius, horizon, mode, placed_nodes_in_shadow)
            else:
                node.node_x, node.node_y, node.node_placed = node.path_x, node.path_y, True
            swarm.append(node)
        self.nodes = swarm
        self.__swarm_xy()
        if self.rotation != 0:
            self.swarm_rotate(self.rotation)

    def plot_path_swarm(self, plot_lines=True,  plot=True):
        # plot the path
        fig, axs = plt.subplots()
        # plot path
        path_x, path_y = zip(*self.i_path)
        axs.plot(path_x, path_y, 'k-', label='Path', zorder=1)
        # plot circles
        for n in self.nodes:
            
            # plot tangent lines
            if plot_lines:
                xs, ys = zip(*((n.path_xo, n.path_yo), (n.node_x, n.node_y))) 
                axs.plot(xs, ys, 'b-', linewidth=0.5, zorder=2)

            circle = Circle((n.node_x, n.node_y), n.node_radius - n.node_radius_buffer, fc='white', edgecolor='black', linewidth=1, alpha=0.9, zorder=3)
            axs.add_patch(circle)
        axs.axis('equal')
        if plot:
            plt.show()
        return fig, axs

class superswarm:

    def __init__(self, df, id_field, position_field, 
        min=None, max=None, size_field=None, order_field=None, 
        buffer=0., size_override=None, rotation=0, tol_overlap=10, tol_r=1e-10,
        shape='c', draw_points=True, kwargs={}):
        
        self.df = df
        self.id_field = id_field
        self.position_field = position_field
        self.min = min
        self.max = max
        self.size_field = size_field
        self.order_field = order_field
        self.buffer = buffer
        self.size_override = size_override
        self.rotation = rotation
        self.tol_overlap = tol_overlap
        self.tol_r = tol_r
        self.shape = shape
        self.draw_points = draw_points
        self.kwargs = kwargs
        
        self.pathswarm = None
        self.path_area = None
        self.super_swarm(df=df, id_field=id_field, position_field=position_field,
            min=min, max=max, size_field=size_field, order_field=order_field,
            buffer=buffer, size_override=size_override, rotation=rotation,
            tol_overlap=tol_overlap, tol_r=tol_r, draw_points=draw_points, kwargs=kwargs)

    @classmethod
    def random_superswarm(cls, size, min=None, max=None, buffer=0., size_override=None, 
        rotation=0, tol_overlap=10, shape='c', draw_points=True, kwargs={}):
        data = [[''.join(random.choices(string.ascii_letters, k=5)), random.randint(1, 100)/10, random.randint(1, 100)/10] for _ in range(size)]
        df = pd.DataFrame(data, columns=['id', 'position', 'size'])
        return cls(df, 'id', 'position', min=min, max=max, size_field='size', order_field=None, buffer=buffer, 
            size_override=size_override, rotation=rotation, tol_overlap=tol_overlap, shape=shape,
            draw_points=draw_points, kwargs=kwargs)

    def shape_path(self, area, curve):
        path = []
        interp = ''
        if self.shape == 'c': # circle
            r = (area/pi)**(1/2)
            path = [(x,y) for x,y,_ in vf.circle(0., 0., r=r, points=20, end_cap=True)]
            interp = 'cubic_spline'
        elif self.shape[0] == 'p': # regular (equal sided) polygon
            n_sides = int(self.shape[1:])
            path = vf.regular_polygon(0., 0., n_sides, area=area, end_cap=True)
            interp = 'linear'
        if curve:
            interp = 'cubic_spline'
        return path, interp

    def super_swarm(self, df, id_field, position_field, min, max, size_field, order_field, 
        buffer, size_override, rotation, tol_overlap, tol_r, draw_points, kwargs):
        # defaults
        defaults = {'sort':'desc', 'size_by':'area', 'order_by':'size',
            'offset':None, 'mode':'closest', 'curve':False}
        params = {**defaults, **kwargs} # default if not initialized
        sort = params['sort']
        size_by = params['size_by']
        order_by = params['order_by']
        offset = params['offset']
        mode = params['mode']
        curve = params['curve']
        # get area of for the closed path
        o_ps_df = pathswarm(df=df, id_field=id_field, position_field=position_field,
            min=min, max=max, df_only=True)
        df_clip = o_ps_df.swarm_setup()
        area = 1.
        size = None
        ids = df_clip[id_field].values
        if self.size_field is not None and self.size_override is None:
            sizes = df_clip[self.size_field].values
        elif self.size_override is not None:
            size = self.size_override 
            sizes = [size for _ in ids]
        else:
            size = (df_clip[position_field].max()-df_clip[position_field].min())/len(df_clip)*10
            sizes = [size for _ in ids]
        df_clip['__size'] = sizes
        if size_by == 'radius':
            area = sum(pi * df_r**2 for df_r in df_clip['__size'].values) # by radius
        elif size_by == 'area':
            area = df_clip['__size'].sum() # by area
        if offset is None:
            offset = area/40.
        # provide path from shape option:
        path, interp = self.shape_path(area, curve)
        # get the pathswarm object
        o_ps = pathswarm(df=df, id_field=id_field, position_field=position_field,
            min=min, max=max, path=path, size_field=size_field, order_field=order_field,
            buffer=buffer, size_override=size_override, rotation=rotation,
            tol_overlap=tol_overlap, tol_r=tol_r, interp=interp, draw_points=draw_points,
            kwargs={'closed':True, 'horizon':'top', 'sort':sort, 'order_by':order_by,
                'offset':offset, 'mode':mode, 'size_by': size_by})
        self.pathswarm = o_ps
        self.path_area = area

    def plot_super_swarm(self, plot_lines=True,  plot=True):
        return self.pathswarm.plot_path_swarm(plot_lines=plot_lines,  plot=plot)

class beeswarm:

    def __init__(self, df, id_field, position_field, min=None, max=None,
        size_field=None, order_field=None, buffer=0., size_override=None, 
        rotation=0, tol_overlap=10, tol_r=1e-10, center_clusters=True, scale_to_data=True, 
        draw_points=True, kwargs={}):
        
        self.df = df
        self.id_field = id_field
        self.position_field = position_field
        self.min = min
        self.max = max
        self.size_field = size_field
        self.order_field = order_field
        self.buffer = buffer
        self.size_override = size_override
        self.rotation = rotation
        self.tol_overlap = tol_overlap
        self.tol_r=tol_r
        self.center_clusters = center_clusters
        self.scale_to_data = scale_to_data
        self.draw_points = draw_points
        self.kwargs = kwargs
        
        self.pathswarm = None
        self.bee_swarm(df=df, id_field=id_field, position_field=position_field,
            min=min, max=max, size_field=size_field, order_field=order_field,
            buffer=buffer, size_override=size_override, rotation=rotation,
            tol_overlap=tol_overlap, tol_r=tol_r, draw_points=draw_points, kwargs=kwargs)
        
    @classmethod
    def random_beeswarm(cls, size, min=None, max=None, buffer=0., size_override=None, 
        rotation=0, tol_overlap=10, center_clusters=True, scale_to_data=True, 
        draw_points=True, kwargs={}):
        data = [[''.join(random.choices(string.ascii_letters, k=5)), random.randint(1, 1000)/10, random.randint(1, 100)/50] for _ in range(size)]
        df = pd.DataFrame(data, columns=['id', 'position', 'size'])
        return cls(df, 'id', 'position', min=min, max=max, size_field='size', order_field=None, 
            buffer=buffer, size_override=size_override, rotation=rotation, tol_overlap=tol_overlap, 
            center_clusters=center_clusters, scale_to_data=scale_to_data, draw_points=draw_points, kwargs=kwargs)
        
    def bee_swarm(self, df, id_field, position_field, min, max, size_field, 
        order_field, buffer, size_override, rotation, tol_overlap, tol_r, 
        draw_points, kwargs):
         # defaults
        defaults = {'sort':'desc', 'size_by':'area', 'order_by':'size',
            'mode':'closest'}
        params = {**defaults, **kwargs} # default if not initialized
        sort = params['sort']
        size_by = params['size_by']
        order_by = params['order_by']
        mode = params['mode']
        # get df_clip
        o_ps_df = pathswarm(df=df, id_field=id_field, position_field=position_field,
            min=min, max=max, df_only=True)
        df_clip = o_ps_df.swarm_setup()
        # define path:
        interp = 'linear'
        path = [(0,0),(10,0)]
        if min is not None or max is not None:
            path = []
            if min is not None:
                path.append((min, 0))
            else:
                path.append((0, 0))
            if max is not None:
                path.append((max, 0))
            else:
                path.append((10, 0))
        if self.scale_to_data:
            path = [(df_clip[position_field].min(), 0.),(df_clip[position_field].max(), 0.)]
        # get the pathswarm object
        o_ps = pathswarm(df=df, id_field=id_field, position_field=position_field,
            min=min, max=max, path=path, size_field=size_field, order_field=order_field,
            buffer=buffer, size_override=size_override, rotation=0,
            tol_overlap=tol_overlap, tol_r=tol_r, interp=interp, draw_points=draw_points,
            kwargs={'closed':False, 'horizon':None, 'sort':sort, 'order_by':order_by,
                'offset':0., 'mode':mode, 'size_by': size_by})
        if self.center_clusters:
            # input setup for clustering
            size = None
            ids = df_clip[id_field].values
            if self.size_field is not None and self.size_override is None:
                sizes = df_clip[self.size_field].values
            elif self.size_override is not None:
                size = self.size_override 
                sizes = [size for _ in ids]
            else:
                size = (df_clip[position_field].max()-df_clip[position_field].min())/len(df_clip)*10
                sizes = [size for _ in ids]
            df_clip['__size'] = sizes
            df_clip_center = df_clip[[id_field, position_field, '__size']].copy(deep=True)
            df_clip_center['__r'] = df_clip_center.apply(lambda row: row['__size'] 
                if size_by == 'radius' else (row['__size'] / np.pi) ** 0.5, axis=1)
            df_clip_center['__min'] = -df_clip_center['__r'] + df_clip_center[position_field]
            df_clip_center['__max'] = df_clip_center['__r'] + df_clip_center[position_field]
            # group by range clustering
            df_clip_center = vf.range_group(df_clip_center, item_col=id_field, min_col='__min', max_col='__max', group_name='__group')
            df_grouped = df_clip_center.groupby(['__group'])
            df_nodes = pd.DataFrame([{'id': obj.id, 'node_y': obj.node_y, 'node_radius': obj.node_radius} for obj in o_ps.nodes])
            group_ids = {id: group_key for group_key, group_df in df_grouped for id in group_df['id'].unique()}
            df_nodes['group'] = df_nodes['id'].map(group_ids)
            df_nodes['max_value'] = df_nodes.groupby('group')['node_y'].transform(lambda x: (x + df_nodes.loc[x.index, 'node_radius']).max())
            df_nodes['min_value'] = df_nodes.groupby('group')['node_y'].transform(lambda x: (x - df_nodes.loc[x.index, 'node_radius']).min())
            df_nodes['midpoint'] =  df_nodes['min_value'] + (df_nodes['max_value'] - df_nodes['min_value']) / 2
            df_nodes['new_path_y'] = df_nodes['node_y'] - df_nodes['midpoint']
            # get y=0 nodes for isolation adjustment and adjust
            path_nodes = [n for n in o_ps.nodes if n.node_y == n.path_yo]
            for node in o_ps.nodes:
                node.node_y = df_nodes.loc[df_nodes['id'] == node.id, 'new_path_y'].values[0]
            isolated_nodes = [n for n in path_nodes if n.id not in {o.node_neighbor for o in o_ps.nodes}]
            for node in isolated_nodes:
                radius = node.node_radius
                x_rotated, y_rotated, placed_nodes_in_shadow = o_ps.node_filter(node)
                if placed_nodes_in_shadow:
                    o_ps.place_node(node, x_rotated, y_rotated, radius, False, mode, placed_nodes_in_shadow)
            # rotation
            if rotation != 0:
                o_ps.swarm_rotate(rotation)
        self.pathswarm = o_ps

    def plot_bee_swarm(self, plot_lines=True,  plot=True):
        return self.pathswarm.plot_path_swarm(plot_lines=plot_lines,  plot=plot)
