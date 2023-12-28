#%%
import pandas as pd
from math import sqrt
import copy
import matplotlib.pyplot as plt
import random
import string

from . import functions as vf
from .draw import points as dp

# import functions as vf
# from draw import points as dp

class quadtile:

    def __init__(self, df, id_field, value_field, xo=0., yo=0., packing='min', overflow=6, buffer=0.5):
        self.df = df
        self.id_field = id_field
        self.value_field = value_field
        self.xo = xo
        self.yo = yo
        self.packing = packing
        self.overflow = overflow
        self.buffer = buffer

        self.o_quadtile_chart = None
        self.quadtile_chart()
    
    @classmethod
    def random_quadtile(cls, size):
        data = [[''.join(random.choices(string.ascii_letters, k=5)), random.randint(1, 1000)] for _ in range(size)]
        df = pd.DataFrame(data, columns=['id', 'value'])
        return cls(df, 'id', 'value')

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

            if packing == 'num':
                #region find max and pick another
                if side == 'top':
                    if i > overflow:
                        side_max = len([o for o in list_xy.viz if o.side == side])
                        next_side_max = len([o for o in list_xy.viz if o.side == 'right'])
                        if side_max > next_side_max:
                            side = 'right'
                            side_list = side_dict[side]
                if side == 'right':
                    if i > overflow:
                        side_max = len([o for o in list_xy.viz if o.side == side])
                        next_side_max = len([o for o in list_xy.viz if o.side == 'bottom'])
                        if side_max > next_side_max:
                            side = 'bottom'
                            side_list = side_dict[side]
                if side == 'bottom':
                    if i > overflow:
                        side_max = len([o for o in list_xy.viz if o.side == side])
                        next_side_max = len([o for o in list_xy.viz if o.side == 'left'])
                        if side_max > next_side_max:
                            side = 'left'
                            side_list = side_dict[side]
                if side == 'left':
                    if i > overflow:
                        side_max = len([o for o in list_xy.viz if o.side == side])
                        next_side_max = len([o for o in list_xy.viz if o.side == 'top'])
                        if side_max > next_side_max:
                            side = 'top'
                            side_list = side_dict[side]
                #endregion
            
            if packing == 'max':
                #region find max and pick another
                if side == 'top':
                    if i > overflow:
                        side_max = max([o.y for o in list_xy.viz if o.side == side]) - w1/2.
                        next_side_max = max([o.x for o in list_xy.viz if o.side == 'right']) - w1/2.
                        if side_max > next_side_max:
                            side = 'right'
                            side_list = side_dict[side]
                if side == 'right':
                    if i > overflow:
                        side_max = max([o.x for o in list_xy.viz if o.side == side]) - w1/2.
                        next_side_max = abs(min([o.y for o in list_xy.viz if o.side == 'bottom']) - w1/2.)
                        if side_max > next_side_max:
                            side = 'bottom'
                            side_list = side_dict[side]
                if side == 'bottom':
                    if i > overflow:
                        side_max = abs(min([o.y for o in list_xy.viz if o.side == side]) - w1/2.)
                        next_side_max = abs(min([o.x for o in list_xy.viz if o.side == 'left']) - w1/2.)
                        if side_max > next_side_max:
                            side = 'left'
                            side_list = side_dict[side]
                if side == 'left':
                    if i > overflow:
                        side_max = abs(min([o.x for o in list_xy.viz if o.side == side]) - w1/2.)
                        next_side_max = max([o.y for o in list_xy.viz if o.side == 'top']) - w1/2.
                        if side_max > next_side_max:
                            side = 'top'
                            side_list = side_dict[side]
                #endregion

            if packing == 'min':
                #region find min and stay
                if side == 'top':
                    if i > overflow:
                        side_min = min([o.y for o in list_xy.viz if o.side == side]) - w1/2.
                        next_side_min = min([o.x for o in list_xy.viz if o.side == 'right']) - w1/2.
                        if side_min > next_side_min:
                            side = 'right'
                            side_list = side_dict[side]
                if side == 'right':
                    if i > overflow:
                        side_min = min([o.x for o in list_xy.viz if o.side == side]) - w1/2.
                        next_side_min = abs(max([o.y for o in list_xy.viz if o.side == 'bottom']) - w1/2.)
                        if side_min > next_side_min:
                            side = 'bottom'
                            side_list = side_dict[side]
                if side == 'bottom':
                    if i > overflow:
                        side_min = abs(max([o.y for o in list_xy.viz if o.side == side]) - w1/2.)
                        next_side_min = abs(max([o.x for o in list_xy.viz if o.side == 'left']) - w1/2.)
                        if side_min > next_side_min:
                            side = 'left'
                            side_list = side_dict[side]
                if side == 'left':
                    if i > overflow:
                        side_min = abs(max([o.x for o in list_xy.viz if o.side == side]) - w1/2.)
                        next_side_min = min([o.y for o in list_xy.viz if o.side == 'top']) - w1/2.
                        if side_min > next_side_min:
                            side = 'top'
                            side_list = side_dict[side]
                #endregion

            #region algorithm
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
            xa, ya = vf.rotate(o.x, o.y, 45., 0., 0.)
            o.x = xa
            o.y = ya
            path_counter += 1
        #endregion

        self.o_quadtile_chart = list_xy
        self.o_quadtile_chart.to_dataframe()

    def quadtile_plot(self, opacity=0.5):
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
        plt.show(block=True)

    def to_df(self):
        return self.o_quadtile_chart.df

    def to_csv(self, file_name):
        self.o_quadtile_chart.dataframe_to_csv(file_name)
