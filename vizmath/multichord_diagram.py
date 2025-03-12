# Nick's Multi-Chord Algorithm

import pandas as pd
import numpy as np
from math import isnan
import random
import string
import matplotlib.pyplot as plt

from . import functions as vf
from .draw import points as dp

# import functions as vf
# from draw import points as dp

class multichord:
    
    def __init__(self, df, multiset_field, value_field, order=None,    
        percent = 100., set_radial_buffer = 25., 
        multiset_offset = 0.04, multiset_height = 0.04, 
        set_offset = 0.08, set_height = 0.04, rotate_deg=0.0,
        upset_output=True):
        
        self.df = df
        self.multiset_field = multiset_field
        self.value_field = value_field
        self.order = order
        self.percent = percent
        self.set_radial_buffer = set_radial_buffer
        self.multiset_offset = multiset_offset
        self.multiset_height = multiset_height
        self.set_offset = set_offset
        self.set_height = set_height
        self.rotate_deg = rotate_deg
        self.multiplier = None

        self.o_multichord = dp()
        self.multichord_diagram()
        self.o_multichord.to_dataframe()
        if upset_output:
            self.upset_df = self.__multichord_upset()

    @classmethod
    def __random_multiset(cls, num_sets, num_multisets, id_length):
        # Ensure the parameters allow generating unique multi-sets covering all sets.
        if num_multisets < num_sets:
            raise ValueError("num_multisets must be at least num_sets")
        if num_multisets > (2**num_sets - 1):
            raise ValueError("num_multisets exceeds possible unique combinations for the given num_sets")
        
        # Generate unique set identifiers
        sets_list = set()
        while len(sets_list) < num_sets:
            sets_list.add(''.join(random.choices(string.ascii_letters, k=id_length)))
        sets_list = list(sets_list)  # Convert to list for indexing

        multisets = set()  # Track unique multi-set combos
        mutiset_list = []
        
        # Ensure every set is represented in at least one multi-set.
        for s in sets_list:
            subset = {s} | {x for x in sets_list if x != s and random.random() < 0.5}
            combo = ','.join(sorted(subset))
            if combo not in multisets:
                multisets.add(combo)
                mutiset_list.append([combo, round(random.uniform(0, 100), 2)])
        
        # Generate additional unique multi-sets until reaching the desired count.
        while len(mutiset_list) < num_multisets:
            subset = set(random.sample(sets_list, random.randint(1, num_sets)))
            combo = ','.join(sorted(subset))
            if combo not in multisets:
                multisets.add(combo)
                mutiset_list.append([combo, round(random.uniform(0, 100), 2)])
        
        return mutiset_list

    @classmethod
    def random_multichord(cls, num_sets=4, num_multisets=5, id_length=2,    
        percent = 100., set_radial_buffer = 25., 
        multiset_offset = 0.04, multiset_height = 0.04, 
        set_offset = 0.08, set_height = 0.04, rotate_deg=0.0,
        upset_output=True):
        multisets = cls.__random_multiset(num_sets, num_multisets, id_length)
        df = pd.DataFrame(multisets, columns = ['multiset', 'value'])
        return cls(df, 'multiset', 'value', order=None,    
        percent=percent, set_radial_buffer=set_radial_buffer, 
        multiset_offset=multiset_offset, multiset_height=multiset_height, 
        set_offset=set_offset, set_height=set_height, rotate_deg=rotate_deg,
        upset_output=upset_output)

    class __group:
        def __init__(self, multiset, set, value, position, offset, x1, x2): 
            self.multiset = multiset
            self.set = set
            self.value = value
            self.position = position
            self.offset = offset
            self.x1 = x1
            self.x2 = x2
        def to_dict(self):
            return {
                'multiset' : self.multiset,
                'set' : self.set,
                'value' : self.value,
                'position' : self.position,
                'offset' : self.offset,
                'x1' : self.x1, 
                'x2' : self.x2 }

    def __polarize(self, x, max_x, y_offset = 0.):
        y = max_x/4.+y_offset
        x_out, y_out = vf.polarize(x, max_x, y, y_offset)
        return x_out, y_out

    def __multichord_calc(self):
        
        #region initialize element list
        element_list = []
        for i, row in self.df.iterrows():
            multiset_i = row[self.multiset_field]
            value_i = row[self.value_field] * self.multiplier
            sets = multiset_i.split(',')
            for j in range(len(sets)):
                element_list.append(self.__group(multiset_i, sets[j], value_i, 0, 0., 0., 0.))
        df_elements = pd.DataFrame.from_records([s.to_dict() for s in element_list])
        #endregion

        #region set up group propperties
        df_group = df_elements.groupby(['set'])['value'].sum().reset_index()
        set_radial_buffer = self.set_radial_buffer
        df_group['offset'] = 0.

        order_list = self.order.split(',')
        df_group['order'] = [order_list.index(x)+1 for x in df_group['set']]
        
        offset = 0.
        for i, row in df_group.sort_values(by=['order'], ascending=True).iterrows():
            df_group.at[i,'offset'] = offset
            offset += row['value']
        dict_group = df_group.set_index('set').T.to_dict('list')
        #endregion

        #region endpoints
        value_offset = 0.
        last_group = order_list[0]
        df_elements['order'] = [order_list.index(x)+1 for x in df_elements['set']]
        for i, row in df_elements.sort_values(by=['order', 'value'], ascending=True).iterrows():
            if row['set'] != last_group:
                value_offset = 0.
                last_group = row['set']
            
            order_i = row['order']
            offset = dict_group[row['set']][1]

            df_elements.at[i,'offset'] = value_offset
            df_elements.at[i,'x1'] = offset + value_offset + set_radial_buffer*(order_i-1)
            df_elements.at[i,'x2'] = offset + value_offset + row['value'] + set_radial_buffer*(order_i-1)

            value_offset += row['value']
        #endregion

        #region setup virtual points
        span = df_elements['x2'].max()+set_radial_buffer
        yest = span/4.

        if self.percent < 100.:
            span = df_elements['x2'].max()*(1./(self.percent/100.))

        df_elements_prior = df_elements.copy()
        df_elements_next = df_elements.copy()
        df_elements_prior['x1'] = df_elements_prior['x1'] - span
        df_elements_prior['x2'] = df_elements_prior['x2'] - span
        df_elements_prior['position'] = -1
        df_elements_next['x1'] = df_elements_next['x1'] + span
        df_elements_next['x2'] = df_elements_next['x2'] + span
        df_elements_next['position'] = 1
        dfs = [df_elements_prior, df_elements, df_elements_next]
        df_combined = pd.concat(dfs, axis=0)
        #endregion

        #region point dictionary
        chord_dict = {}
        for name, multiset_i in df_combined.groupby(['multiset', 'position']):
            x_list = []
            key = name
            for i, row in multiset_i.sort_values(by=['order'], ascending=True).iterrows():
                e = row['set']
                x1 = row['x1']
                x2 = row['x2']
                x_list.append([e,x1,x2])
            chord_dict[key] = x_list
        #endregion

        #region select point path
        chord_arc_dict = {}
        for i, row in self.df.iterrows():
            list_i = chord_dict[row[self.multiset_field],0]
            list_neg = chord_dict[row[self.multiset_field],-1]
            list_pos = chord_dict[row[self.multiset_field],1]

            if len(list_i) > 1:
                list_new = []
                p1 = list_i[0][2]
                p2 = 0.
                for j in range(len(list_i)):

                    p1_0 = p1
                    if j == len(list_i)-1:
                        j = -1

                    p2_i = list_i[j+1][1]
                    p2_neg = list_neg[j+1][1]
                    p2_pos = list_pos[j+1][1]

                    i_delta = abs(p1-p2_i)
                    neg_delta = abs(p1-p2_neg)
                    pos_delta = abs(p1-p2_pos)

                    if i_delta < neg_delta:
                        if i_delta < pos_delta:
                            p2 = p2_i
                            p1 = list_i[j+1][2]
                        else:
                            p2 = p2_pos
                            p1 = list_pos[j+1][2]
                    elif pos_delta < neg_delta:
                        p2 = p2_pos
                        p1 = list_pos[j+1][2]
                    else:
                        p2 = p2_neg
                        p1 = list_neg[j+1][2]

                    list_new.append([row[self.multiset_field],p1_0, p2])
                chord_arc_dict[row[self.multiset_field]] = list_new
            else:
                chord_arc_dict[row[self.multiset_field]] = list_i
        #endregion

        #region draw point path
        for i, row in self.df.iterrows():
            multiset_i = row[self.multiset_field]
            value = row[self.value_field]
            count = len(multiset_i.split(','))
            list_path = chord_arc_dict[row[self.multiset_field]]

            if count == 1:
                #region 1 element
                p1 = list_path[0][1]
                p2 = list_path[0][2]
                points = int(abs(p2-p1)/(span)/2.*720)+10
                s1 = np.linspace(p2, p1, num=points)
                path = 1
                for k in range(len(s1)):
                    x, y = self.__polarize(s1[k], span)
                    self.o_multichord.append(multiset_i, x, y, path, value=value, count=count, type='chord')
                    path += 1
                x1, y1 = self.__polarize(p1, span)
                x2, y2 = self.__polarize(p2, span)
                chord_list = vf.chord(0., 0., x1, y1, x2, y2, 20, 0.25)
                for k in range(len(chord_list)):
                    self.o_multichord.append(multiset_i, chord_list[k][0], chord_list[k][1], path, value=value, count=count, type='chord')
                    path += 1
                #endregion
            
            else:
                #region multiple elements

                #region first straight segment
                path = 1
                p1 = chord_dict[row[self.multiset_field],0][0][1]
                p2 = list_path[0][1]
                points = int(abs(p2-p1)/(span)/2.*720)
                #straight between p1 to p2
                s1 = np.linspace(p1, p2, num=points)
                for k in range(len(s1)):
                        x, y = self.__polarize(s1[k], span)
                        self.o_multichord.append(multiset_i, x, y, path, value=value, count=count, type='chord')
                        path += 1        
                #endregion
                
                for j in range(len(list_path)):

                    #region curves
                    p2 = list_path[j][1]
                    x1, y1 = self.__polarize(p2, span)
                    p3 = list_path[j][2]
                    x2, y2 = self.__polarize(p3, span)
                    chord_list = vf.chord(0., 0., x1, y1, x2, y2, 50)
                    for k in range(len(chord_list)):
                            self.o_multichord.append(multiset_i, chord_list[k][0], chord_list[k][1], path, value=value, count=count, type='chord')
                            path += 1
                    #endregion

                    #region lagging straight segments
                    if j < len(list_path)-1:
                        p4 = list_path[j+1][1]
                        #straight between p3 to p4
                        points = int(abs(p4-p3)/(span)/2.*720)
                        s2 = np.linspace(p3, p4, num=points)
                        for k in range(len(s2)):
                            x, y = self.__polarize(s2[k], span)
                            self.o_multichord.append(multiset_i, x, y, path, value=value, count=count, type='chord')
                            path += 1
                    #endregion

                #endregion
        #endregion

        #region rescale to unit radius
        for i in range(len(self.o_multichord.viz)):
            self.o_multichord.viz[i].x = vf.rescale(self.o_multichord.viz[i].x, -span/4., span/4., -1., 1.)
            self.o_multichord.viz[i].y = vf.rescale(self.o_multichord.viz[i].y, -span/4., span/4., -1., 1.)
        #endregion

        #region add legends
        set_offset = self.set_offset
        for i, row in df_elements.iterrows():
            multiset_i =  row.multiset
            value = row.value
            count = len(multiset_i.split(','))
            set = row.set

            p1 = row.x1
            p2 = row.x2
            points = int(abs(p2-p1)/(span)/2.*720)+3
            
            path = 1
            top = np.linspace(p1, p2, num=points)
            x1 = 0
            y1 = 0
            initial = True
            for k in range(len(top)):
                top_k = vf.rescale(top[k], -span/4., span/4., -1., 1.)
                x, y = self.__polarize(top_k, 4., self.multiset_offset)
                if initial:
                    x1 = x
                    y1 = y
                    initial = False
                self.o_multichord.append(set + ' - ' + multiset_i, x, y, path, value=value, count=count, type='multiset')
                path += 1
            bottom = np.linspace(p2, p1, num=points)
            for k in range(len(bottom)):
                bottom_k = vf.rescale(bottom[k], -span/4., span/4., -1., 1.)
                x, y = self.__polarize(bottom_k, 4., self.multiset_offset + self.multiset_height)
                self.o_multichord.append(set + ' - ' + multiset_i, x, y, path, value=value, count=count, type='multiset')
                path += 1
            self.o_multichord.append(set + ' - ' + multiset_i, x1, y1, path, value=value, count=count, type='multiset')

        for i, row in df_group.iterrows():
            set_i =  row.set
            value = row.value
            count = None

            p1 = row.offset + (row.order-1)*set_radial_buffer
            p2 = row.offset + row.value + (row.order-1)*set_radial_buffer
            points = int(abs(p2-p1)/(span)/2.*720)+3
            
            path = 1
            top = np.linspace(p1, p2, num=points)
            x1 = 0
            y1 = 0
            initial = True
            for k in range(len(top)):
                top_k = vf.rescale(top[k], -span/4., span/4., -1., 1.)
                x, y = self.__polarize(top_k, 4., set_offset)
                if initial:
                    x1 = x
                    y1 = y
                    initial = False
                self.o_multichord.append(set_i, x, y, path, value=value, count=count, type='set')
                path += 1
            bottom = np.linspace(p2, p1, num=points)
            for k in range(len(bottom)):
                bottom_k = vf.rescale(bottom[k], -span/4., span/4., -1., 1.)
                x, y = self.__polarize(bottom_k, 4., self.set_offset + self.set_height)
                self.o_multichord.append(set_i, x, y, path, value=value, count=count, type='set')
                path += 1
            self.o_multichord.append(set_i, x1, y1, path, value=value, count=count, type='set')
        #endregion

        #region rotate
        if self.rotate_deg != 0:
            for i in range(len(self.o_multichord.viz)):
                x, y = self.o_multichord.viz[i].x, self.o_multichord.viz[i].y
                xr, yr = vf.rotate(x, y, self.rotate_deg)
                self.o_multichord.viz[i].x, self.o_multichord.viz[i].y = xr, yr
        #endregion

    def multichord_diagram(self):
        
        val_max = self.df[self.value_field].max()
        self.multiplier = 100/val_max

        if self.order == None:
            order_dict = {}
            for i, row in self.df.iterrows():
                items = row[self.multiset_field].split(',')
                value = row[self.value_field]
                for j in range(len(items)):
                    if items[j] in order_dict:
                        order_dict[items[j]] += value
                    else:
                        order_dict[items[j]] = value
            df_order = pd.DataFrame(order_dict.items(), columns=[self.multiset_field, self.value_field])
            order = df_order.sort_values(by=[self.value_field], ascending=False)[self.multiset_field].values.tolist()
            self.order = ",".join(order)
        
        self.__multichord_calc()

    def multichord_plot(self, level = None, transparency = 0.5):
        max_count = int(max(self.o_multichord.df['count'].unique()))
        colors = []
        for i in range(max_count):
            r = random.random()
            b = random.random()
            g = random.random()
            color = (r, g, b)
            colors.append(color)
        df_lvl_group = self.o_multichord.df.groupby(['type','id'])
        fig, axs = plt.subplots()
        axs.set_aspect('equal', adjustable='box')
        for group, rows in df_lvl_group:
            x = rows['x'].values
            y = rows['y'].values
            count_check = rows['count'].values[0]
            color = ()
            set_linewidth=0.5
            if isinstance(count_check, (int, float)) and not isnan(count_check):
                color = colors[int(count_check)-1]
                if int(count_check) == level:
                    set_linewidth=3
            else:
                r = random.random()
                b = random.random()
                g = random.random()
                color = (r, g, b)
            axs.fill(x, y, alpha=transparency, fc=color)
            plt.plot(x, y, 'k-', linewidth=set_linewidth)
        plt.show(block=True)

    def __multichord_upset(self):
        # df_mc = self.o_multichord.df[self.o_multichord.df['type']=='chord'][['id','value']]
        # df_mc_cc = pd.concat([df_mc, df_mc['id'].str.get_dummies(sep=',')], axis = 1)
        # df_mc_upset = pd.melt(df_mc_cc, id_vars=['id', 'value'], var_name='id_2', value_name='value_2')
        # return df_mc_upset[df_mc_upset['id_2'] != 'path'].drop_duplicates()

        df_multiset = pd.concat([self.df, self.df[self.multiset_field].str.get_dummies(sep=',')], axis = 1)
        df_upset = pd.melt(df_multiset, id_vars=[self.multiset_field, self.value_field], var_name = 'set', value_name = 'link')
        return df_upset
