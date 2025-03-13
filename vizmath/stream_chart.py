# Nick's Stream Chart Algorithm
#%%
import random
import matplotlib.pyplot as plt

from . import functions as vf
from .draw import points as dp

# import functions as vf
# from draw import points as dp

#%%
class stream:

    def __init__(self, df, x_field, item_field, value_field, buffer=0., step_points=10, step_factor=4.):
        
        self.df = df
        self.x_field = x_field
        self.item_field = item_field
        self.value_field = value_field
        self.buffer = buffer
        self.step_points = step_points
        self.step_factor = step_factor

        self.o_stream = dp()
        self.stream_chart()
        self.o_stream.to_dataframe()

    def stream_chart(self):

        item_dict = {}
        value_add = 0.
        df_x = self.df.groupby([self.x_field], sort=True)
        item_field = self.item_field
        value_field = self.value_field
        buffer = self.buffer
        points = self.step_points
        factor = self.step_factor
        iter = 1

        for name, group in df_x:
            order = 0
            value_add = 0.
            num_items = len(group)
            y_offset = (group[value_field].sum() + buffer * (num_items-1))/2
            x = name[0]
            for i, row in group.sort_values(by=[value_field], ascending=True).iterrows():
                item = row[item_field]
                value = row[value_field]
                path = 1
                value += value_add - y_offset
                rank=num_items-order
                if iter > 1:
                    if item in item_dict:
                        path = item_dict[item][1]
                        path += 1
                    else:
                        item_dict[item] = [value, 1]
                    if value == item_dict[item][0]:
                        self.o_stream.append(item, x, value, path, value=row[value_field], rank=rank)
                    else:
                        sigmoid = vf.sigmoid(last_x, item_dict[item][0], x, value, points, limit=factor)
                        for i in range(len(sigmoid)):
                            self.o_stream.append(item, sigmoid[i][0], sigmoid[i][1], path, value=row[value_field], rank=rank)
                            path += 1
                else:
                    self.o_stream.append(item, x, value, 1, value=row[value_field], rank=rank)
                item_dict[item] = [value, path]
                value += buffer
                order += 1
                value_add = value + y_offset
            last_x = x
            iter += 1

        item_dict_2 = {}
        iter = 1
        for name, group in df_x:
            order = 0
            value_add = 0.
            num_items = len(group)
            y_offset = (group[value_field].sum() + buffer * (num_items-1))/2
            x = name[0]
            for i, row in group.sort_values(by=[value_field], ascending=True).iterrows():
                item = row[item_field]
                value = row[value_field]
                path = 0
                value += value_add - y_offset
                rank=num_items-order
                if iter > 1:
                    if item not in item_dict_2:
                        item_dict_2[item] = [value - row[value_field], 0]
                    if value - row[value_field] == item_dict_2[item][0]:
                        #add same y, new x
                        self.o_stream.append(item, x, value - row[value_field], path, value=row[value_field], rank=rank)
                    else:
                        #sigmoid
                        sigmoid = vf.sigmoid(last_x, item_dict_2[item][0], x, value - row[value_field], points, limit=factor)
                        for i in range(len(sigmoid)):
                            self.o_stream.append(item, sigmoid[i][0], sigmoid[i][1], path, value=row[value_field], rank=rank)
                else:
                    self.o_stream.append(item, x, value - row[value_field], 1, value=row[value_field], rank=rank)
                
                item_dict_2[item] = [value - row[value_field], path]
                value += buffer
                value_add = value + y_offset
            last_x = x
            iter += 1

        for i in range(len(self.o_stream.viz)-1, -1, -1):
            item = self.o_stream.viz[i].id
            if self.o_stream.viz[i].path != 0:
                break
            path = item_dict[item][1]+1
            item_dict[item][1] = path
            self.o_stream.viz[i].path = path
    
    def stream_plot(self, opacity=0.5, show=True):
        fig, axs = plt.subplots()
        # axs.set_aspect('equal', adjustable='box')

        df_stream = self.o_stream.df.groupby(['id'])
        for group, rows in df_stream:
            rows = rows.sort_values(by='path')
            x = rows['x'].values
            y = rows['y'].values
            r = random.random()
            b = random.random()
            g = random.random()
            color = (r, g, b)
            axs.fill(x, y, alpha=opacity, fc=color)
            axs.plot(x, y, 'k-', linewidth=0.5)

        if show:
            plt.show()
        return fig, axs
