#%%
import pandas as pd
from math import sqrt, inf
import copy
import matplotlib.pyplot as plt
import random
import string
import numpy as np
from scipy.spatial import ConvexHull

from . import functions as vf
from .draw import points as dp

# import functions as vf
# from draw import points as dp

class quadtile:

    def __init__(self, df, id_field, value_field, xo=0., yo=0., packing='auto', overflow=6, buffer=0.5, rotate=45., constraints=None):
        self.df = df
        self.id_field = id_field
        self.value_field = value_field
        self.xo = xo
        self.yo = yo
        self.packing = packing
        self.overflow = overflow
        self.buffer = buffer
        self.rotate = rotate
        if constraints is not None:
            self.constraints = vf.sort_vertices(constraints)

        self.o_quadtile_chart = None
        self.o_squares = dp()
        self.quadtile_chart()
    
    @classmethod
    def random_quadtile(cls, size, buffer=0.5, rotate=45.):
        data = [[''.join(random.choices(string.ascii_letters, k=5)), random.randint(1, 1000)] for _ in range(size)]
        df = pd.DataFrame(data, columns=['id', 'value'])
        return cls(df, 'id', 'value', buffer=buffer, rotate=rotate)

    class __line:
        def __init__(self, index, x1, y1, x2, y2): 
            self.index = index
            self.x1 = x1
            self.y1 = y1
            self.x2 = x2
            self.y2 = y2
            self.len = sqrt((x2-x1)**2+(y2-y1)**2)
        def update_len(self):
            self.len = sqrt((self.x2-self.x1)**2+(self.y2-self.y1)**2)

    def __add_square(self, list_xy, index, i, x, y, w, a, side):
        list_xy.append(id=index, x=x, y=y, path=0, index=index, item=i, w=w, a=a, side=side)
        index += 1
        list_xy.append(id=index, x=x, y=y+w, path=1, index=index, item=i, w=w, a=a, side=side)
        index += 1
        list_xy.append(id=index, x=x+w, y=y+w, path=2, index=index, item=i, w=w, a=a, side=side)
        index += 1
        list_xy.append(id=index, x=x+w, y=y, path=3, index=index, item=i, w=w, a=a, side=side)
        index += 1
        list_xy.append(id=index, x=x, y=y, path=4, index=index, item=i, w=w, a=a, side=side)
        index += 1
        self.o_squares.append(id=i, x=x, y=y, path=0, w=w, a=a, side=side)

    def __max_packing(self, side, i, overflow, list_xy, w1):
        # find max and pick another
        if side == 'top':
            if i > overflow:
                side_max = max([o.y for o in list_xy.viz if o.side == side]) - w1/2.
                next_side_max = max([o.x for o in list_xy.viz if o.side == 'right']) - w1/2.
                if side_max > next_side_max:
                    return 'right'
        if side == 'right':
            if i > overflow:
                side_max = max([o.x for o in list_xy.viz if o.side == side]) - w1/2.
                next_side_max = abs(min([o.y for o in list_xy.viz if o.side == 'bottom']) - w1/2.)
                if side_max > next_side_max:
                    return 'bottom'
        if side == 'bottom':
            if i > overflow:
                side_max = abs(min([o.y for o in list_xy.viz if o.side == side]) - w1/2.)
                next_side_max = abs(min([o.x for o in list_xy.viz if o.side == 'left']) - w1/2.)
                if side_max > next_side_max:
                    return 'left'
        if side == 'left':
            if i > overflow:
                side_max = abs(min([o.x for o in list_xy.viz if o.side == side]) - w1/2.)
                next_side_max = max([o.y for o in list_xy.viz if o.side == 'top']) - w1/2.
                if side_max > next_side_max:
                    return 'top'
        return side

    def __min_packing(self, side, i, overflow, list_xy, w1):
        # find min and stay
        if side == 'top':
            if i > overflow:
                side_min = min([o.y for o in list_xy.viz if o.side == side]) - w1/2.
                next_side_min = min([o.x for o in list_xy.viz if o.side == 'right']) - w1/2.
                if side_min > next_side_min:
                    return 'right'
        if side == 'right':
            if i > overflow:
                side_min = min([o.x for o in list_xy.viz if o.side == side]) - w1/2.
                next_side_min = abs(max([o.y for o in list_xy.viz if o.side == 'bottom']) - w1/2.)
                if side_min > next_side_min:
                    return 'bottom'
        if side == 'bottom':
            if i > overflow:
                side_min = abs(max([o.y for o in list_xy.viz if o.side == side]) - w1/2.)
                next_side_min = abs(max([o.x for o in list_xy.viz if o.side == 'left']) - w1/2.)
                if side_min > next_side_min:
                    return 'left'
        if side == 'left':
            if i > overflow:
                side_min = abs(max([o.x for o in list_xy.viz if o.side == side]) - w1/2.)
                next_side_min = min([o.y for o in list_xy.viz if o.side == 'top']) - w1/2.
                if side_min > next_side_min:
                    return 'top'
        return side

    def __num_packing(self, side, i, overflow, list_xy):
        # find max and pick another
        if side == 'top':
            if i > overflow:
                side_num = len([o for o in list_xy.viz if o.side == side])
                next_side_num = len([o for o in list_xy.viz if o.side == 'right'])
                if side_num > next_side_num:
                    return 'right'
        if side == 'right':
            if i > overflow:
                side_num = len([o for o in list_xy.viz if o.side == side])
                next_side_num = len([o for o in list_xy.viz if o.side == 'bottom'])
                if side_num > next_side_num:
                    return 'bottom'
        if side == 'bottom':
            if i > overflow:
                side_num = len([o for o in list_xy.viz if o.side == side])
                next_side_num = len([o for o in list_xy.viz if o.side == 'left'])
                if side_num > next_side_num:
                    return 'left'
        if side == 'left':
            if i > overflow:
                side_num = len([o for o in list_xy.viz if o.side == side])
                next_side_num = len([o for o in list_xy.viz if o.side == 'top'])
                if side_num > next_side_num:
                    return 'top'
        return side

    def __inc_packing(self, side, side_list, w, polygon):
        x = 0.
        y = 0.
        next_side = ''
        if side == 'top':
            for o in sorted(side_list, key=lambda o: (o.y1, -o.len)):
                lengnth = o.len
                if w <= lengnth:
                    x = o.x1
                    y = o.y1
                    break
            next_side = 'right'
        if side == 'right':
            for o in sorted(side_list, key=lambda o: (o.x1, -o.len)):
                lengnth = o.len
                if w <= lengnth:
                    x = o.x1
                    y = o.y1 - w
                    break
            next_side = 'bottom'
        if side == 'bottom':
            for o in sorted(side_list, key=lambda o: (-o.y1, -o.len)):
                lengnth = o.len
                if w <= lengnth:
                    x = o.x1 - w
                    y = o.y1 - w
                    break
            next_side = 'left'
        if side == 'left':
            for o in sorted(side_list, key=lambda o: (-o.x1, -o.len)):
                lengnth = o.len
                if w <= lengnth:
                    x = o.x1 - w
                    y = o.y1
                    break
            next_side = 'top'
        point_in = vf.is_point_in_convex_polygon(x, y, polygon)
        if not point_in:
            return next_side
        return side

    def quadtile_chart(self):

        input = self.df[[self.id_field, self.value_field]].values.tolist()
        squares_sorted = sorted(input, key=lambda x: x[1], reverse=True)
        widths_sorted = copy.deepcopy(squares_sorted)
        xo = self.xo
        yo = self.yo
        packing = self.packing
        overflow = self.overflow
        buffer = self.buffer

        for i in range(len(widths_sorted)):
            val = widths_sorted[i][1]
            val = sqrt(val) + buffer*2
            widths_sorted[i][1] = val

        #region initialize
        list_xy = dp()
        index = 0
        w1 = widths_sorted[0][1]
        side = 'top'
        side_dict = {'top' : [self.__line(1, xo, yo+w1, xo+w1, yo+w1)],
            'right' : [self.__line(1, xo+w1, yo+w1, xo+w1, yo)],
            'bottom' : [self.__line(1, xo+w1, yo, xo, yo)],
            'left' : [self.__line(1, xo, yo, xo, yo+w1)]}
        item = widths_sorted[0][0]
        a = squares_sorted[0][1]
        self.__add_square(list_xy, index, item, xo, yo, w1, a, 'center')
        x = 0
        y = 0
        #endregion

        #region loop algorithm
        for i in range(1,len(widths_sorted)):

            if self.packing == 'auto':
                if i % 3 == 0:
                    packing = 'max'
                elif i % 5 == 0:
                    packing = 'min'
                else:
                    packing = 'num'

            item = widths_sorted[i][0]
            w = widths_sorted[i][1]
            a = squares_sorted[i][1]
            side_list = side_dict[side]

            p = 0
            while p < 4: # try all 4 sides
                side_p = side
                try:
                    if packing == 'inc':
                        side = self.__inc_packing(side, side_list, w, self.constraints)
                    if packing == 'num':
                        side = self.__num_packing(side, i, overflow, list_xy)
                    if packing == 'max':
                        side = self.__max_packing(side, i, overflow, list_xy, w1)
                    if packing == 'min':
                        side = self.__min_packing(side, i, overflow, list_xy, w1)
                    side_list = side_dict[side]
                except:
                    side = side_p
                if side_p == side:
                    break
                p += 1

            #region algorithm (build up for each side until adjacent side helps build out)
            if side == 'top':
                # inspect lowest y in list with longest length :)
                for o in sorted(side_list, key=lambda o: (o.y1, -o.len)):
                    lengnth = o.len
                    if w <= lengnth:
                        x = o.x1
                        y = o.y1
                        self.__add_square(list_xy, index, item, x, y, w, a, side)
                        # replace segment with new segments :)
                        o.x1 += w
                        o.update_len()
                        idx = max([o.index for o in side_list]) + 1
                        side_dict[side].append(self.__line(idx, x, y+w, x+w, y+w))
                        break
                if x == xo:
                    # blend segment on left side :)
                    side_list_adj = side_dict['left']
                    for o in sorted(side_list_adj, key=lambda o: (-o.y2)):
                        o.y2 = y+w
                        o.update_len()
                        break
                side = 'right'
                continue
            if side == 'right':
                # inspect lowest x in list with longest length :)
                for o in sorted(side_list, key=lambda o: (o.x1, -o.len)):
                    lengnth = o.len
                    if w <= lengnth:
                        x = o.x1
                        y = o.y1 - w
                        self.__add_square(list_xy, index, item, x, y, w, a, side)
                        # replace segment with new segments :)
                        o.y1 -= w
                        o.update_len()
                        idx = max([o.index for o in side_list]) + 1
                        side_dict[side].append(self.__line(idx, x+w, y+w, x+w, y))
                        break
                if y+w == yo+w1:
                    # blend segment on left side :)
                    side_list_adj = side_dict['top']
                    for o in sorted(side_list_adj, key=lambda o: (-o.x2)):
                        o.x2 = x+w
                        o.update_len()
                        break
                side = 'bottom'
                continue
            if side == 'bottom':
                # inspect highest y in list with longest length :)
                for o in sorted(side_list, key=lambda o: (-o.y1, -o.len)):
                    lengnth = o.len
                    if w <= lengnth:
                        x = o.x1 - w
                        y = o.y1 - w
                        self.__add_square(list_xy, index, item, x, y, w, a, side)
                        # replace segment with new segments :)
                        o.x1 -= w
                        o.update_len()
                        idx = max([o.index for o in side_list]) + 1
                        side_dict[side].append(self.__line(idx, x+w, y, x, y))
                        break
                if x+w == xo+w1:
                    # blend segment on left side :)
                    side_list_adj = side_dict['right']
                    for o in sorted(side_list_adj, key=lambda o: (o.y2)):
                        o.y2 = y
                        o.update_len()
                        break
                side = 'left'
                continue
            if side == 'left':
                # inspect highest x in list with longest length :)
                for o in sorted(side_list, key=lambda o: (-o.x1, -o.len)):
                    lengnth = o.len
                    if w <= lengnth:
                        x = o.x1 - w
                        y = o.y1
                        self.__add_square(list_xy, index, item, x, y, w, a, side)
                        # replace segment with new segments :)
                        o.y1 += w
                        o.update_len()
                        idx = max([o.index for o in side_list]) + 1
                        side_dict[side].append(self.__line(idx, x, y, x, y+w))
                        break
                if y == yo:
                    # blend segment on left side
                    side_list_adj = side_dict['bottom']
                    for o in sorted(side_list_adj, key=lambda o: (o.x2)):
                        o.x2 = x
                        o.update_len()
                        break
                side = 'top'
            #endregion
        #endregion

        #region rotate
        path_counter = 1
        for o in sorted(list_xy.viz, key=lambda o: (o.item, o.path)):
            if path_counter == 6:
                path_counter = 1
            if path_counter == 1 or path_counter == 5:
                o.x += buffer
                o.y += buffer
            if path_counter == 2:
                o.x += buffer
                o.y -= buffer
            if path_counter == 3:
                o.x -= buffer
                o.y -= buffer
            if path_counter == 4:
                o.x -= buffer
                o.y += buffer
            xa, ya = vf.rotate(o.x, o.y, self.rotate, 0., 0.)
            o.x = xa
            o.y = ya
            path_counter += 1
        #endregion

        self.o_quadtile_chart = list_xy
        self.o_quadtile_chart.to_dataframe()

    def quadtile_plot(self, opacity=0.5, show_constraints=False, polygon=None):
        colors = {}
        sides = {0:'center', 1:'top', 2:'right', 3:'bottom', 4:'left'}
        for i in range(5):
            r = random.random()
            b = random.random()
            g = random.random()
            color = (r, g, b)
            colors[sides[i]] = color
        df_lvl_group = self.o_quadtile_chart.df.groupby(['item'])
        fig, axs = plt.subplots()
        axs.set_aspect('equal', adjustable='box')
        for group, rows in df_lvl_group:
            x = rows['x'].values
            y = rows['y'].values
            c = colors[rows['side'].values[0]]
            set_linewidth=0.5
            axs.fill(x, y, alpha=opacity, fc=c)
            plt.plot(x, y, 'k-', linewidth=set_linewidth)
        if show_constraints:
            sorted_polygon = self.constraints
            x_sorted_polygon, y_sorted_polygon = zip(*sorted_polygon)
            x_sorted_polygon += (x_sorted_polygon[0],) # close the polygon
            y_sorted_polygon += (y_sorted_polygon[0],)
            plt.plot(x_sorted_polygon, y_sorted_polygon, 'b-', label='Sorted Polygon')
            plt.fill(x_sorted_polygon, y_sorted_polygon, 'skyblue', alpha=0.3)
        if polygon is not None:
            sorted_polygon = vf.sort_vertices(polygon)
            x_sorted_polygon, y_sorted_polygon = zip(*sorted_polygon)
            x_sorted_polygon += (x_sorted_polygon[0],) # close the polygon
            y_sorted_polygon += (y_sorted_polygon[0],)
            plt.plot(x_sorted_polygon, y_sorted_polygon, 'b-', label='Sorted Polygon')
            plt.fill(x_sorted_polygon, y_sorted_polygon, 'skyblue', alpha=0.3)
        plt.show(block=True)

    def to_df(self):
        return self.o_quadtile_chart.df

    def to_csv(self, file_name):
        self.o_quadtile_chart.dataframe_to_csv(file_name)

class polyquadtile:

    def __init__(self, df, id_field, value_field, xo=0., yo=0., 
            buffer=0.05, rotate=45., constraints=None, sides=['top','right','bottom','left'],
            collapse=False, auto=True, auto_max_iter=20, auto_min_val=0.0001, auto_max_val=20.):
        
        self.df = df
        self.id_field = id_field
        self.value_field = value_field
        self.xo = xo
        self.yo = yo
        self.buffer = buffer
        self.rotate = rotate
        self.sides = sides
        self.collapse = collapse
        if constraints is not None:
            self.constraints = vf.sort_vertices(constraints)

        self.o_polyquadtile_chart = None
        self.o_polysquares = dp()
        self.multiplier = 1.
        if len(df) < 3:
            qt = quadtile(df, id_field, value_field, xo, yo, buffer=buffer,
                rotate=rotate, constraints=constraints)
            self.o_polyquadtile_chart = qt.o_quadtile_chart
            self.o_polysquares = qt.o_squares
        else:
            if auto:
                self.auto_polyquadtile_chart(auto_max_iter, auto_min_val, auto_max_val)
            else:
                self.polyquadtile_chart(df, value_field, rotate)

    @classmethod
    def __angle_between(cls, p1, p2, p3):
        v1, v2 = p1 - p2, p3 - p2
        angle = np.degrees(np.math.atan2(np.linalg.det([v1,v2]),np.dot(v1,v2)))
        return abs(angle)

    @classmethod
    def __generate_random_convex_polygon(cls):
        # Randomly choose the number of sides between 4 and 8
        num_sides = np.random.randint(4, 9)

        # Generate points on a circle with random radii
        angles = np.sort(np.random.uniform(0, 2*np.pi, num_sides))
        radii = np.random.uniform(0.1, np.sqrt(200/np.pi), num_sides)
        points = np.array([radii * np.cos(angles), radii * np.sin(angles)]).T

        # Create the polygon using ConvexHull
        hull = ConvexHull(points)
        polygon = points[hull.vertices]

        # Check for minimum angle and area condition
        area = ConvexHull(polygon).volume
        num_polygon_sides = len(polygon)
        min_angle = min([cls.__angle_between(polygon[i-1], polygon[i], polygon[(i+1) % num_polygon_sides]) for i in range(num_polygon_sides)])
        
        if area >= 200 or min_angle <= 20:
            return cls.__generate_random_convex_polygon()

        # Adjust the centroid to be at (0,0)
        centroid = np.mean(polygon, axis=0)
        adjusted_polygon = polygon - centroid

        return adjusted_polygon

    @classmethod
    def random_polyquadtile(cls, size, xo=0., yo=0., buffer=0.05, rotate=45., sides=['top','right','bottom','left'],
            collapse=False, auto=True, auto_max_iter=20, auto_min_val=0.0001, auto_max_val=20.):
        data = [[''.join(random.choices(string.ascii_letters, k=5)), random.randint(1, 1000)] for _ in range(size)]
        df = pd.DataFrame(data, columns=['id', 'value'])
        poly = cls.__generate_random_convex_polygon()
        return cls(df, 'id', 'value', xo=xo, yo=yo, buffer=buffer, rotate=rotate, constraints=poly, sides=sides,
            collapse=collapse, auto=auto, auto_max_iter=auto_max_iter, auto_min_val=auto_min_val, auto_max_val=auto_max_val)

    def __seg_len_to_ploy(self, x, y, poly):
        x0 = -1
        y0 = y
        intercepts = vf.line_polygon_intercepts(x0, y0, x, y, poly)
        intercept = max(intercepts, key=lambda p: p[0]) # max along the x-axis
        segment_length = vf.distance_between_two_points(x, y, intercept[0], intercept[1])
        return segment_length

    class __segment:
        def __init__(self, x, y, length, height=inf, active=True): 
            self.x = x
            self.y = y
            self.length = length
            self.height = height
            self.active = active
        def __eq__(self, other):
            return self.x == other.x and self.y == other.y and self.length == other.length and self.height == other.height
        def update_length(self, poly):
            self.length = self.__seg_len_to_ploy(self.x, self.y, poly)
    
    def __quad_vertices(self, offset):
        # assume a center of 0, 0 and shift the x forward and y down by the half-width
        polygon = self.constraints
        quad = {
            'top' : [vf.rotate(x, y, 0, -offset, offset) for x, y in polygon],
            'right' : [vf.rotate(x, y, -90, -offset, offset) for x, y in polygon],
            'bottom' : [vf.rotate(x, y, -180, -offset, offset) for x, y in polygon],
            'left' : [vf.rotate(x, y, -270, -offset, offset) for x, y in polygon]
        }
        return quad
    
    def __place(self, segments, segment, seg_idx, w, poly, lvl, lvl_h):
        placed = False
        # check the top side of the square to see if the points are within the polygon
        top1 = vf.is_point_in_convex_polygon(segment.x, segment.y+w, poly)
        top2 = vf.is_point_in_convex_polygon(segment.x+w, segment.y+w, poly)
        # split the segment into 2 new ones after setting the square (or 1 new one if there's complete overlap)
        if top1 and top2 and w <= segment.length and segment.y+w <= segment.height: # checking segment height for under-bridge segments
            placed = True
            if segment.y == lvl:
                lvl_h = segment.y+w
            # child segments shares the same height threshold as the parent in case the parent is an under-bridge segment
            if w == segment.length: # 1 segment (full overlap)
                segments[seg_idx] = self.__segment(segment.x, segment.y+w, w, segment.height, False)
            else: # 2 segments
                segments[seg_idx:seg_idx+1] = [
                    self.__segment(segment.x, segment.y+w, w, segment.height, False),
                    self.__segment(segment.x+w, segment.y, segment.length-w, segment.height, False)]
        else:
            if seg_idx is not None: 
                segments[seg_idx].active = False
        return placed, lvl_h
    
    def __add_square(self, list_xy, index, i, x, y, w, a, side):
        list_xy.append(id=index, x=x, y=y, path=0, index=index, item=i, w=w, a=a, side=side)
        index += 1
        list_xy.append(id=index, x=x, y=y+w, path=1, index=index, item=i, w=w, a=a, side=side)
        index += 1
        list_xy.append(id=index, x=x+w, y=y+w, path=2, index=index, item=i, w=w, a=a, side=side)
        index += 1
        list_xy.append(id=index, x=x+w, y=y, path=3, index=index, item=i, w=w, a=a, side=side)
        index += 1
        list_xy.append(id=index, x=x, y=y, path=4, index=index, item=i, w=w, a=a, side=side)
        index += 1
        self.o_polysquares.append(id=i, x=x, y=y, path=0, w=w, a=a, side=side)

    def __next_side(self, side='top'):
        if side == 'top':
            if 'right' in self.sides:
                return 'right'
            if 'bottom' in self.sides:
                return 'bottom'
            if 'left' in self.sides:
                return 'left'
            return 'top'
        if side == 'right':
            if 'bottom' in self.sides:
                return 'bottom'
            if 'left' in self.sides:
                return 'left'
            if 'top' in self.sides:
                return 'top'
            return 'right'
        if side == 'bottom':
            if 'left' in self.sides:
                return 'left'
            if 'top' in self.sides:
                return 'top'
            if 'right' in self.sides:
                return 'right'
            return 'bottom'
        if side == 'left':
            if 'top' in self.sides:
                return 'top'
            if 'right' in self.sides:
                return 'right'
            if 'bottom' in self.sides:
                return 'bottom'
            return 'left'

    def __squares_to_centroids(self, offset, rotate):
        for o in self.o_polysquares.viz:
            o.x += o.w/2.
            o.y += o.w/2.
            if o.side == 'center':
                revert_x = offset
                revert_y = offset
                revert_r = 0.
            if o.side == 'top':
                revert_x = offset
                revert_y = -offset
                revert_r = 0.
            if o.side == 'right':
                revert_x = -offset
                revert_y = -offset
                revert_r = 90.
            if o.side == 'bottom':
                revert_x = -offset
                revert_y = offset
                revert_r = 180.
            if o.side == 'left':
                revert_x = offset
                revert_y = offset
                revert_r = 270.
            xa, ya = vf.rotate(o.x, o.y, revert_r, revert_x, revert_y)
            o.x = xa
            o.y = ya
        for o in self.o_polysquares.viz:
            xa, ya = vf.rotate(o.x, o.y, rotate)
            o.x = xa+self.xo
            o.y = ya+self.yo
        self.o_polysquares.to_dataframe()

    def __no_room(self, side_dict):
        top = False
        if 'top' in self.sides:
            if side_dict['top']['no_room']:
                top = True
        else:
            top = True
        right = False
        if 'right' in self.sides:
            if side_dict['right']['no_room']:
                right = True
        else:
            right = True
        bottom = False
        if 'bottom' in self.sides:
            if side_dict['bottom']['no_room']:
                bottom = True
        else:
            bottom = True
        left = False
        if 'left' in self.sides:
            if side_dict['left']['no_room']:
                left = True
        else:
            left = True
        no_room = False
        if top and right and bottom and left:
            no_room = True
        return no_room

    def polyquadtile_chart(self, df, value_field, rotate):

        #region initialize
        list_xy = dp()
        buffer = self.buffer
        # sort the values and calculate the width (including the buffer)
        df_sorted = df.sort_values(by=value_field, ascending=False)
        df_sorted['width'] = (df_sorted[value_field].apply(sqrt)) + buffer * 2
        df_center = df_sorted.iloc[0:1] # first record
        df_quad = df_sorted.iloc[1:] # quad records
        # calculate the length of the first segment between the origin and the polygon boundary
        center_w = df_center['width'].values[0]
        offset = center_w/2. # offset by the half-width
        poly_quad = self.__quad_vertices(offset)
        # center details
        center_id = df_center[self.id_field].values[0]
        center_a = df_center[value_field].values[0]
        # retrieve values for looping
        widths = df_quad['width'].values
        ids = df_quad[self.id_field].values
        areas = df_quad[value_field].values
        self.__add_square(list_xy, 0, center_id, 0., 0., center_w, center_a, 'center')
        

        side_dict = {
            'top' : {
                'segments' : [self.__segment(0., 0., length=self.__seg_len_to_ploy(0., 0., poly_quad['top']))],
                'lvl' : 0.,
                'lvl_h' : 0.,
                'no_room' : False
            },
            'right' : {
                'segments' : [self.__segment(0., 0., length=self.__seg_len_to_ploy(0., 0., poly_quad['right']))],
                'lvl' : 0.,
                'lvl_h' : 0.,
                'no_room' : False
            },
            'bottom' : {
                'segments' : [self.__segment(0., 0., length=self.__seg_len_to_ploy(0., 0., poly_quad['bottom']))],
                'lvl' : 0.,
                'lvl_h' : 0.,
                'no_room' : False
            },
            'left' : {
                'segments' : [self.__segment(0., 0., length=self.__seg_len_to_ploy(0., 0., poly_quad['left']))],
                'lvl' : 0.,
                'lvl_h' : 0.,
                'no_room' : False
            }}

        #endregion

        #region algorithm
        side = self.sides[0]
        i = 0

        #region debugging
        # fig, axs = plt.subplots()
        # axs.set_aspect('equal', adjustable='box')
        # x_sorted_polygon, y_sorted_polygon = zip(*poly_quad[side])
        # x_sorted_polygon += (x_sorted_polygon[0],) # close the polygon
        # y_sorted_polygon += (y_sorted_polygon[0],)
        # axs.plot(x_sorted_polygon, y_sorted_polygon, 'b-')
        # axs.fill(x_sorted_polygon, y_sorted_polygon, 'skyblue', alpha=0.3)
        # display(fig)
        # plt.show()
        #endregion

        while i < len(widths):

            #region debugging
            # if i == 3:
            #     print('debug')

            # if next((s for s in side_dict[side]['segments'] if s.x == 4.717317018950091 and s.y ==1.635081752286866), None) is not None:
            #     print('debug')
            #endregion

            w = widths[i]
            id = ids[i]
            a = areas[i]
            # get the last segment/bridge at the end of the line
            last_lvl_segment = max((s for s in side_dict[side]['segments'] if s.y == side_dict[side]['lvl']), key=lambda s: s.x)

            # segments sorted by smallest height (for bridge segments) then smallest x that are active
            sorted_segments = copy.deepcopy(sorted([s for s in side_dict[side]['segments'] if s.active], key=lambda s: (s.height, s.x)))
            for j in range(len(sorted_segments)):
                
                segment = sorted_segments[j]
                seg_idx = next((i for i, s in enumerate(side_dict[side]['segments']) if s == segment), None)
                placed, side_dict[side]['lvl_h'] = self.__place(side_dict[side]['segments'], segment, seg_idx, w, poly_quad[side], side_dict[side]['lvl'], side_dict[side]['lvl_h']) # place square if possible

                if segment == last_lvl_segment and not placed and side_dict[side]['lvl'] != side_dict[side]['lvl_h']: # create a bridge segment
                    last_lvl_segment.height = side_dict[side]['lvl_h'] # set a threshold for the under-bridge segment
                    side_dict[side]['lvl'] = side_dict[side]['lvl_h'] # reset the level
                    # check if a segment exists at the side_dict[side]['lvl_h'], if so, extend it, otherwise create a new one
                    lvl_h_seg = max((s for s in side_dict[side]['segments'] if s.y == side_dict[side]['lvl_h']), key=lambda s: s.x, default=None)
                    if lvl_h_seg is not None: # extend it
                        lvl_h_seg.length += self.__seg_len_to_ploy(lvl_h_seg.x, side_dict[side]['lvl_h'], poly_quad[side])
                    else: # create new
                        side_dict[side]['segments'].append(self.__segment(segment.x, side_dict[side]['lvl_h'], length=self.__seg_len_to_ploy(segment.x, side_dict[side]['lvl_h'], poly_quad[side])))
                    [setattr(s, 'active', True) for s in side_dict[side]['segments']] # reset all to active
                    break
                elif placed:
                    # add square
                    self.__add_square(list_xy, i+1, id, segment.x, segment.y, w, a, side)
                    side_dict[side]['no_room'] = False

                    #region degbugging:
                    # square_plot = [s for s in list_xy.viz if s.item == id]
                    # xp = [s.x for s in square_plot]
                    # yp = [s.y for s in square_plot]
                    # r = random.random()
                    # b = random.random()
                    # g = random.random()
                    # cp = (r, g, b)
                    # set_linewidth=0.5
                    # axs.plot(xp, yp, 'k-', linewidth=set_linewidth)
                    # axs.fill(xp, yp, 'skyblue', alpha=0.3)
                    # display(fig)
                    # plt.show()
                    #endregion

                    if j+1 == len(sorted_segments): # last segment
                        [setattr(s, 'active', True) for s in side_dict[side]['segments']] # reset all to active
                    i += 1 # move on to the next square
                    # adjust the side
                    side = self.__next_side(side)
                    break
                else:
                    if j+1 == len(sorted_segments): # last segment reached so add more bridges if possible:
                        # if (x,y) of the last_lvl_segment match the lower right corner of an existing square
                        #   then add another bridge at the top right corner of that square
                        #   (or extend a segment if one is already there)
                        upper_square = next((s for s in self.o_polysquares.viz if s.x +s.w == last_lvl_segment.x and s.y == last_lvl_segment.y), None)
                        if upper_square is not None:
                            # set height restriction on old bridge
                            last_lvl_segment.height = upper_square.w
                            side_dict[side]['lvl'] += upper_square.w
                            side_dict[side]['lvl_h'] = side_dict[side]['lvl']
                            # check if a segment exists at the side_dict[side]['lvl_h'], if so, extend it, otherwise create a new one
                            lvl_h_seg = max((s for s in side_dict[side]['segments'] if s.y == side_dict[side]['lvl']), key=lambda s: s.x, default=None)
                            if lvl_h_seg is not None: # extend it
                                lvl_h_seg.length += self.__seg_len_to_ploy(lvl_h_seg.x, side_dict[side]['lvl'], poly_quad[side])
                            else: # create new
                                usx = upper_square.x + upper_square.w
                                side_dict[side]['segments'].append(self.__segment(usx, side_dict[side]['lvl'], length=self.__seg_len_to_ploy(usx, side_dict[side]['lvl'], poly_quad[side])))
                            [setattr(s, 'active', True) for s in side_dict[side]['segments']] # reset all to active
                            break
                        # no fit found, reset and check other sides
                        side_dict[side]['no_room'] = True
                        [setattr(s, 'active', True) for s in side_dict[side]['segments']] # reset all to active
                        side = self.__next_side(side)
                        # check if all sides haven't found a fit, if so then finished
                        # if side_dict['top']['no_room'] and side_dict['right']['no_room'] and side_dict['bottom']['no_room'] and side_dict['left']['no_room']:
                        if self.__no_room(side_dict):
                            i = len(widths) # no segments work
        #endregion

        #region rotate
        path_counter = 1
        for o in sorted(list_xy.viz, key=lambda o: (o.item, o.path)):
            if path_counter == 6:
                path_counter = 1
            if path_counter == 1 or path_counter == 5:
                o.x += buffer
                o.y += buffer
            if path_counter == 2:
                o.x += buffer
                o.y -= buffer
            if path_counter == 3:
                o.x -= buffer
                o.y -= buffer
            if path_counter == 4:
                o.x -= buffer
                o.y += buffer
            if o.side == 'center':
                revert_x = offset
                revert_y = offset
                revert_r = 0.
            if o.side == 'top':
                revert_x = offset
                revert_y = -offset
                revert_r = 0.
            if o.side == 'right':
                revert_x = -offset
                revert_y = -offset
                revert_r = 90.
            if o.side == 'bottom':
                revert_x = -offset
                revert_y = offset
                revert_r = 180.
            if o.side == 'left':
                revert_x = offset
                revert_y = offset
                revert_r = 270.
            xa, ya = vf.rotate(o.x, o.y, revert_r, revert_x, revert_y)
            o.x = xa
            o.y = ya
            path_counter += 1
        for o in list_xy.viz:
            xa, ya = vf.rotate(o.x, o.y, rotate)
            o.x = xa+self.xo
            o.y = ya+self.yo
        rotated_polygon = []
        for o in self.constraints:
            xc, yc = vf.rotate(o[0], o[1], rotate)
            rotated_polygon.append((xc+self.xo, yc+self.yo))
        self.constraints = rotated_polygon
        #endregion

        self.o_polyquadtile_chart = list_xy
        self.o_polyquadtile_chart.to_dataframe()
        self.__squares_to_centroids(offset, rotate)

    def auto_polyquadtile_chart(self, max_iter=20, min_value=0.0001, max_value=20.):
        df = self.df
        val_f = self.value_field
        df['apqt_norm_i'] = df[val_f].apply(lambda x: 0.01 + 0.99 * (x - df[val_f].min()) / (df[val_f].max() - df[val_f].min()))
        squares = len(df)
        # bisect to find the right multiplier
        i = 1
        auto_pqt = True
        mid_point = 1.
        while auto_pqt:
            mid_point = (min_value + max_value) / 2
            df['apqt_norm'] = df['apqt_norm_i']
            df['apqt_norm'] = df['apqt_norm'] * mid_point
            squares_fit = 0
            try:
                self.o_polyquadtile_chart = None
                self.o_polysquares = dp()
                self.polyquadtile_chart(df, 'apqt_norm', 0)
                squares_fit = len(self.o_polysquares.viz)
            except:
                pass
            if squares == squares_fit:
                if i >= max_iter:
                    break
                min_value = mid_point
            else:
                max_value = mid_point
            if i > max_iter + 5:
                auto_pqt = False
            i += 1
        self.multiplier = mid_point
        if self.rotate != 0:
            self.o_polyquadtile_chart = None
            self.o_polysquares = dp()
            self.polyquadtile_chart(df, 'apqt_norm', self.rotate)
            squares_fit = len(self.o_polysquares.viz)

    def polyquadtile_plot(self, opacity=0.5, show_constraints=False, polygon=None):
        colors = {}
        sides = {0:'center', 1:'top', 2:'right', 3:'bottom', 4:'left'}
        for i in range(5):
            r = random.random()
            b = random.random()
            g = random.random()
            color = (r, g, b)
            colors[sides[i]] = color
        df_lvl_group = self.o_polyquadtile_chart.df.groupby(['item'])
        fig, axs = plt.subplots()
        axs.set_aspect('equal', adjustable='box')
        for group, rows in df_lvl_group:
            x = rows['x'].values
            y = rows['y'].values
            c = colors[rows['side'].values[0]]
            set_linewidth=0.5
            axs.fill(x, y, alpha=opacity, fc=c)
            plt.plot(x, y, 'k-', linewidth=set_linewidth)
        if show_constraints:
            sorted_polygon = self.constraints
            x_sorted_polygon, y_sorted_polygon = zip(*sorted_polygon)
            x_sorted_polygon += (x_sorted_polygon[0],) # close the polygon
            y_sorted_polygon += (y_sorted_polygon[0],)
            plt.plot(x_sorted_polygon, y_sorted_polygon, 'b-', label='Sorted Polygon')
            plt.fill(x_sorted_polygon, y_sorted_polygon, 'skyblue', alpha=0.3)
        if polygon is not None:
            sorted_polygon = vf.sort_vertices(polygon)
            x_sorted_polygon, y_sorted_polygon = zip(*sorted_polygon)
            x_sorted_polygon += (x_sorted_polygon[0],) # close the polygon
            y_sorted_polygon += (y_sorted_polygon[0],)
            plt.plot(x_sorted_polygon, y_sorted_polygon, 'b-', label='Sorted Polygon')
            plt.fill(x_sorted_polygon, y_sorted_polygon, 'skyblue', alpha=0.3)
        plt.show(block=True)

    def to_df(self):
        return self.o_polyquadtile_chart.df

    def to_csv(self, file_name):
        self.o_polyquadtile_chart.dataframe_to_csv(file_name)
