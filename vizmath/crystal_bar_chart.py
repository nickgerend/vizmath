# Nick's Crystal Bar Chart Algorithm

from math import inf
import random
import matplotlib.pyplot as plt
# from IPython.display import display # for continued plotting

from . import functions as vf
from .draw import points as dp

# import functions as vf
# from draw import points as dp

class crystals:

    def __init__(self, df, id_field, height_field, height_range, width_field = None, 
        bottom_up = False, width_override = None, offset=0., reset_origin=False, rotation=0.):

        self.df = df
        self.id_field = id_field
        self.height_field = height_field
        self.height_range = height_range
        self.width_field = width_field
        self.bottom_up = bottom_up
        self.width_override = width_override
        self.offset = offset
        self.reset_origin = reset_origin
        self.rotation = rotation

        self.o_crystal_bar_chart = None
        self.crystal_bar_chart()

    def __place_width(self, side, width, mw, rw, lw, level=0, count=1, middle_taken=False):
        if side == 'm':
            xl = -width/2.
            xr = width/2.
        elif side == 'r':
            if (count == 1 or middle_taken == False) and level == 0 :
                xl = width/2. - width/2.
                xr = width/2. + width/2.
            elif middle_taken:
                xl = rw + mw/2. + width/2. - width/2.
                xr = rw + mw/2. + width/2. + width/2.
            else:
                xl = rw + width/2. - width/2.
                xr = rw + width/2. + width/2.
        elif side == 'l':
            if (count == 1 or middle_taken == False) and level == 0:
                xl = -width/2. - width/2.
                xr = -width/2. + width/2.
            elif middle_taken:
                xl = -lw - mw/2. - width/2. - width/2.
                xr = -lw - mw/2. - width/2. + width/2.
            else:
                xl = -lw - width/2. - width/2.
                xr = -lw - width/2. + width/2.
        return xl, xr

    def __nextside(self, side, group_i=0, count=1):
        next_side = ''
        if side == 'm':
            next_side = 'r'
        elif side == 'r':
            next_side = 'l'
        elif side == 'l':
            if group_i == count-1:
                next_side = 'm'
            else:
                next_side = 'r'
        return next_side

    def __draw_top(self, list_xy, id, xmin, xmax, ymin, ymax, value, group, height, width):
        y = (ymax-ymin)/2.+ymin
        x = (xmax-xmin)/2.+xmin
        list_xy.append(id=id,x=xmin,y=y,path=0,side=0,value=value,group=group,height=height,width=width)
        list_xy.append(id=id,x=x,y=ymax,path=1,side=0,value=value,group=group,height=height,width=width)
        list_xy.append(id=id,x=xmax,y=y,path=2,side=0,value=value,group=group,height=height,width=width)
        list_xy.append(id=id,x=x,y=ymin,path=3,side=0,value=value,group=group,height=height,width=width)
        list_xy.append(id=id,x=xmin,y=y,path=4,side=0,value=value,group=group,height=height,width=width)

    def __draw_left(self, list_xy, id, xmin, xmax, ymin, ymax, value, group, height, width):
        y = (ymax-ymin)/2.+ymin
        x = (xmax-xmin)/2.+xmin

        if x >= 0:
            if y >= 0:
                list_xy.append(id=id,x=xmin,y=y,path=0,side=1,value=value,group=group,height=height,width=width)
                list_xy.append(id=id,x=x,y=ymin,path=1,side=1,value=value,group=group,height=height,width=width)
                list_xy.append(id=id,x=0.,y=0.,path=2,side=1,value=value,group=group,height=height,width=width)
                list_xy.append(id=id,x=xmin,y=y,path=3,side=1,value=value,group=group,height=height,width=width)
            else:
                list_xy.append(id=id,x=xmin,y=y,path=0,side=1,value=value,group=group,height=height,width=width)
                list_xy.append(id=id,x=x,y=ymax,path=1,side=1,value=value,group=group,height=height,width=width)
                list_xy.append(id=id,x=0.,y=0.,path=2,side=1,value=value,group=group,height=height,width=width)
                list_xy.append(id=id,x=xmin,y=y,path=3,side=1,value=value,group=group,height=height,width=width)
        else:
            if y >= 0:
                slope1 = inf
                if x != 0:
                    slope1 = abs(ymin)/abs(x)
                slope2 = abs(y)/abs(xmin)
                if slope1 > slope2 and ymin > 0:
                    list_xy.append(id=id,x=xmin,y=y,path=0,side=1,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=x,y=ymin,path=1,side=1,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=0.,y=0.,path=2,side=1,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=xmin,y=y,path=3,side=1,value=value,group=group,height=height,width=width)
                else:
                    list_xy.append(id=id,x=x,y=ymax,path=0,side=1,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=xmax,y=y,path=1,side=1,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=0.,y=0.,path=2,side=1,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=x,y=ymax,path=3,side=1,value=value,group=group,height=height,width=width)
            else:
                slope1 = inf
                if x != 0:
                    slope1 = abs(ymax)/abs(x)
                slope2 = abs(y)/abs(xmin)
                if slope1 > slope2 and ymax < 0:
                    list_xy.append(id=id,x=xmin,y=y,path=0,side=1,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=x,y=ymax,path=1,side=1,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=0.,y=0.,path=2,side=1,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=xmin,y=y,path=3,side=1,value=value,group=group,height=height,width=width)
                else:
                    list_xy.append(id=id,x=x,y=ymin,path=0,side=1,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=xmax,y=y,path=1,side=1,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=0.,y=0.,path=2,side=1,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=x,y=ymin,path=3,side=1,value=value,group=group,height=height,width=width)

    def __draw_right(self, list_xy, id, xmin, xmax, ymin, ymax, value, group, height, width):
        y = (ymax-ymin)/2.+ymin
        x = (xmax-xmin)/2.+xmin
        if x >= 0:
            if y >= 0:
                slope1 = inf
                if x != 0:
                    slope1 = abs(ymin)/abs(x)
                slope2 = abs(y)/abs(xmax)
                if slope1 > slope2 and ymin > 0:
                    list_xy.append(id=id,x=x,y=ymin,path=0,side=2,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=xmax,y=y,path=1,side=2,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=0.,y=0.,path=2,side=2,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=x,y=ymin,path=3,side=2,value=value,group=group,height=height,width=width)
                else:
                    list_xy.append(id=id,x=x,y=ymax,path=0,side=2,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=xmin,y=y,path=1,side=2,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=0.,y=0.,path=2,side=2,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=x,y=ymax,path=3,side=2,value=value,group=group,height=height,width=width)
            else:
                slope1 = inf
                if x != 0:
                    slope1 = abs(ymax)/abs(x)
                slope2 = abs(y)/abs(xmax)
                if slope1 > slope2 and ymax < 0:
                    list_xy.append(id=id,x=x,y=ymax,path=0,side=2,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=xmax,y=y,path=1,side=2,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=0.,y=0.,path=2,side=2,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=x,y=ymax,path=3,side=2,value=value,group=group,height=height,width=width)
                else:
                    list_xy.append(id=id,x=x,y=ymin,path=0,side=2,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=xmin,y=y,path=1,side=2,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=0.,y=0.,path=2,side=2,value=value,group=group,height=height,width=width)
                    list_xy.append(id=id,x=x,y=ymin,path=3,side=2,value=value,group=group,height=height,width=width)
        else:
            if y >= 0:
                list_xy.append(id=id,x=x,y=ymin,path=0,side=2,value=value,group=group,height=height,width=width)
                list_xy.append(id=id,x=xmax,y=y,path=1,side=2,value=value,group=group,height=height,width=width)
                list_xy.append(id=id,x=0.,y=0.,path=2,side=2,value=value,group=group,height=height,width=width)
                list_xy.append(id=id,x=x,y=ymin,path=3,side=2,value=value,group=group,height=height,width=width)
            else:
                list_xy.append(id=id,x=x,y=ymax,path=0,side=2,value=value,group=group,height=height,width=width)
                list_xy.append(id=id,x=xmax,y=y,path=1,side=2,value=value,group=group,height=height,width=width)
                list_xy.append(id=id,x=0.,y=0.,path=2,side=2,value=value,group=group,height=height,width=width)
                list_xy.append(id=id,x=x,y=ymax,path=3,side=2,value=value,group=group,height=height,width=width)

    def __crystal_bar(self, df_in, id_field, height_field, height_range, width_field = None, 
        bottom_up = False, width_override = None):
        if width_override is not None:
            df = df_in.rename({id_field: 'item', height_field: 'x'}, axis=1)
            df['width'] = width_override
        else:
            df = df_in.rename({id_field: 'item', height_field: 'x', width_field: 'width'}, axis=1)
        df = df[['item', 'x', 'width']]
        df = df.sort_values('x', ascending=bottom_up)
        df['min'] = df['x'] - height_range/2.
        df['max'] = df['x'] + height_range/2.
        df = vf.range_group(df)
        df_groups = df.groupby('group')
        list_xy = dp()
        side = 'm'
        for group, rows in df_groups:
            ll = 0
            rl = 0
            mw = 0.
            lw = 0.
            rw = 0.
            xmin = 0.
            xmax = 0.
            count = len(rows)
            middle_first = False
            group_i = 0
            for i, row in rows.sort_values('x', ascending=bottom_up).iterrows():
                width = row['width']
                id = row['item']
                ymin = row['min']
                ymax = row['max']
                value = row['x']
                if side == 'm':
                    mw += width
                    xmin,xmax = self.__place_width(side, width, mw, rw, lw, count=count)
                    side = self.__nextside(side)
                    if group_i == 0:
                        middle_first = True
                elif side == 'r':
                    xmin,xmax = self.__place_width(side, width, mw, rw, lw, level=rl, count=count, middle_taken=middle_first)
                    side = self.__nextside(side, group_i, count)
                    rl += 1
                    rw += width
                elif side == 'l':
                    xmin,xmax = self.__place_width(side, width, mw, rw, lw, level=ll, count=count, middle_taken=middle_first)
                    side = self.__nextside(side, group_i, count)
                    ll += 1
                    lw += width
                self.__draw_top(list_xy, id, xmin, xmax, ymin, ymax, value, group, height_range, width)
                self.__draw_left(list_xy, id, xmin, xmax, ymin, ymax, value, group, height_range, width)
                self.__draw_right(list_xy, id, xmin, xmax, ymin, ymax, value, group, height_range, width)
                group_i += 1
        self.o_crystal_bar_chart = list_xy

    def crystal_bar_chart(self):
        
        self.df[self.height_field] -= self.offset
        self.__crystal_bar(self.df, self.id_field, self.height_field, self.height_range,
            self.width_field, self.bottom_up, self.width_override)
        if not self.reset_origin:
            [o.__setattr__('y', o.y + self.offset) for o in self.o_crystal_bar_chart.viz]
        if self.rotation != 0:
            [o.__setattr__('x', new_x) or o.__setattr__('y', new_y) for o, (new_x, new_y) in zip(self.o_crystal_bar_chart.viz, [vf.rotate(o.x, o.y, self.rotation) for o in self.o_crystal_bar_chart.viz])]
        self.df[self.height_field] += self.offset
        self.o_crystal_bar_chart.to_dataframe()

    def cbc_plot(self, opacity=0.8, legend=True, alternate_color=False, color=True, show=True):
        fig, axs = plt.subplots()
        axs.set_aspect('equal', adjustable='box')
        df_lvl_group = self.o_crystal_bar_chart.df.groupby(['id','side'])
        colors = []
        for i in range(self.o_crystal_bar_chart.df['group'].max()+1):
            r = random.random()
            b = random.random()
            g = random.random()
            c = (r, g, b)
            colors.append(c)
            if alternate_color:
                if i % 2 == 0:
                    if color:
                        colors[i] = colors[0]
                    else:
                        colors[i] = (1,1,1)
                else:
                    if color:
                        colors[i] = colors[1]
                    else:
                        colors[i] = (.75,.75,.75)
            if legend:
                axs.plot([], [], color=colors[i], label=i+1)
        nested_values = {name: group['group'].iloc[0] for name, group in df_lvl_group}
        sorted_group_names = sorted(nested_values, key=nested_values.get, reverse=self.bottom_up)
        if legend:
            axs.legend(bbox_to_anchor=(1, 1))
        for name in sorted_group_names:
            rows = df_lvl_group.get_group(name)
            x = rows['x'].values
            y = rows['y'].values
            i = rows['group'].values[0]
            axs.fill(x, y, alpha=opacity, fc=colors[i], linewidth=0.5, edgecolor='black')
        if show:
            plt.show()
        return fig, axs

    def to_df(self):
        return self.o_crystal_bar_chart.df

    def to_csv(self, file_name):
        self.o_crystal_bar_chart.dataframe_to_csv(file_name)
