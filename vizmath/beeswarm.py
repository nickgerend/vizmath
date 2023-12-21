
#%%
# Nick's (modified) Beeswarm Chart

from math import inf, pi, sqrt
import numpy as np
import random
import matplotlib.pyplot as plt

from . import functions as vf
from .draw import points as dp

# import functions as vf
# from draw import points as dp

class swarm():

    def __init__(self, df, id_field, position_field, size_field=None, 
        buffer=0., size_override=None, rotation=0):

        self.df = df
        self.id_field = id_field
        self.position_field = position_field
        self.size_field = size_field
        self.buffer = buffer
        self.size_override = size_override
        self.rotation = rotation

        self.o_beeswarm = None
        self.beeswarm()

    def beeswarm(self):
        self.df['item'] = self.df[self.id_field]
        self.df['x'] = self.df[self.position_field]
        if self.size_override is not None:
            self.df['A'] = self.size_override
        else:
            if self.size_field is None:
                max_x = self.df['x'].max()
                self.df['A'] = max_x/20.
            else:
                self.df['A'] = self.df[self.size_field]
        self.df['Ao'] = self.df['A']
        self.df['r'] = np.sqrt(self.df['A']/pi) + self.buffer
        self.df['A'] = pi*self.df['r']**2
        self.df['min'] = -self.df['r'] + self.df['x']
        self.df['max'] = self.df['r'] + self.df['x']
        self.df = vf.range_group(self.df)
        df_grouped = self.df.groupby(['group'])
        list_xy = dp()
        dict_groups = {}
        group_i = 0
        for group, rows in df_grouped:
            rows_sort = rows.sort_values(by=['A'], ascending=False)
            for i, row in rows_sort.iterrows():
                item = row['item']
                x = row['x']
                r = row['r']
                A = row['A']
                y = 0.
                if group_i == 0:
                    dict_groups[group] = [(x,0.,r)]
                    # list_xy.append(point(item, x, 0., A, group))
                    list_xy.append(id=item, x=x, y=0., path=0, area=A, group=group)
                else:
                    yi = 0.
                    for j in range(len(dict_groups[group])):
                        #try top and bottom of every prior-placed item
                        #place if no collision and minimum y-extent
                        xo = dict_groups[group][j][0]
                        yo = dict_groups[group][j][1]
                        ro = dict_groups[group][j][2]
                        yt = vf.circle_collide(xo, yo, x, ro, r, place='top')
                        if yt == None:
                            if j == 0:
                                ylast = 0.
                            continue
                        yb = vf.circle_collide(xo, yo, x, ro, r, place='bottom')
                        # check collion against every other placed item
                        ytest = [yt,yb]
                        collided = False
                        for k in range(len(dict_groups[group])):
                            xk = dict_groups[group][k][0]
                            yk = dict_groups[group][k][1]
                            rk = dict_groups[group][k][2]
                            if vf.circle_collided(xk, yk, x, yt, rk, r):
                                ytest[0] = None
                            if vf.circle_collided(xk, yk, x, yb, rk, r):
                                ytest[1] = None
                            if ytest[0] == None and ytest[1] == None:
                                collided = True
                                break
                        if collided:
                            continue
                        if ytest[0] == None:
                            yi = yb
                        elif ytest[1] == None:
                            yi = yt
                        else:
                            if abs(yt) > abs(yb):
                                yi = yb
                            else:
                                yi = yt
                        if j == 0:
                            ylast = yi
                            y = yi
                        elif abs(yi) < abs(ylast) or ylast == 0.:
                            y = yi
                        ylast = yi
                    dict_groups[group].append((x,y,r))
                    list_xy.append(id=item, x=x, y=y, path=0, area=A, group=group)
                ylast = 0.
                group_i += 1
            group_i = 0

        self.o_beeswarm = dp()
        for i in range(len(list_xy.viz)):
            item = list_xy.viz[i].id
            x = list_xy.viz[i].x
            y = list_xy.viz[i].y
            A = list_xy.viz[i].area
            group = list_xy.viz[i].group
            r = sqrt(A/pi)
            r -= self.buffer
            A = pi*r**2
            circle_i = vf.circle(x, y, r=r, end_cap=True)
            for j in circle_i:
                self.o_beeswarm.append(id=item, x=j[0], y=j[1], path=j[2], area=A, group=group, xo=x, yo=y)
        if self.rotation != 0:
            [o.__setattr__('x', new_x) or o.__setattr__('y', new_y) for o, (new_x, new_y) in zip(self.o_beeswarm.viz, [vf.rotate(o.x, o.y, self.rotation) for o in self.o_beeswarm.viz])]
        self.o_beeswarm.to_dataframe()

    def beeswarm_plot(self, opacity=0.5, color=True):
        fig, axs = plt.subplots()
        axs.set_aspect('equal', adjustable='box')
        df_group = self.o_beeswarm.df.groupby(['id'])
        c = None
        for group, rows in df_group:
            x = rows['x'].values
            y = rows['y'].values
            r = random.random()
            b = random.random()
            g = random.random()
            if color:
                c = (r, g, b)
            else:
                c = (1, 1, 1)
            axs.fill(x, y, alpha=opacity, fc=c)
            plt.plot(x, y, 'k-', linewidth=1.2)
        plt.show()
