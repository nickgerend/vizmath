# Nick's Radial Treeamp Algorithm

import pandas as pd
import numpy as np
from math import pi, cos, sin, sqrt, log
import random
import matplotlib.pyplot as plt
import copy
# from IPython.display import display # for continued plotting

from . import functions as vf
from .draw import points as dp

# import functions as vf
# from draw import points as dp

class rad_treemap:

    def __init__(self, df, groupers, value=None, r1=1, r2=2, a1=0, a2=360, points=200,
        default_sort=False, default_sort_override=True, default_sort_override_reversed=True,
        mode='smart', no_groups=False, rotate_deg=0, full=False, rectangular=False):
        
        self.df = df
        self.groupers = groupers
        self.value = value
        self.r1 = r1
        self.r2 = r2
        self.a1 = a1
        self.a2 = a2
        self.points = points
        self.default_sort = default_sort
        self.default_sort_override = default_sort_override
        self.default_sort_override_reversed = default_sort_override_reversed
        self.mode = mode
        self.no_groups = no_groups
        self.rotate_deg = rotate_deg
        self.full = full
        self.rectangular = rectangular

        self.o_rad_treemap = None
        self.rad_treemap()

    @classmethod
    def __generate_hierarchical_data(cls, num_levels=3, num_top_level_items=10, items_range=(2,5), value_range=(1,100), 
        outlier_fraction=0.2, use_log=True, sig=100):
    
        data = []
        min_val, max_val = value_range
        min_items, max_items = items_range
        # adjust the value range for each top-level item to add more variability
        top_level_ranges = [(random.uniform(min_val, max_val/2), random.uniform(max_val/2, max_val)) for _ in range(num_top_level_items)]

        def generate_value(current_range):
            # decide if the value should be an outlier
            if random.random() < outlier_fraction:
                value_range = current_range
            else:
                value_range = (current_range[0], current_range[0] + (current_range[1] - current_range[0]) * (1 - outlier_fraction))
            # generate value based on specified distribution
            if use_log:
                mean = (log(value_range[0]) + log(value_range[1])) / 2
                sigma = (log(value_range[1]) - mean) / sig
                return random.lognormvariate(mean, sigma)
            else:
                return random.uniform(*value_range)

        def create_rows(prefix, level_idx):
            if level_idx == num_levels:
                # determine the value range based on the top level
                current_range = top_level_ranges[int(prefix[0][1:]) - 1] if level_idx > 1 else value_range
                value = generate_value(current_range)
                data.append(prefix + [value])
            else:
                for item in range(random.randint(min_items, max_items) if level_idx > 0 else num_top_level_items):
                    create_rows(prefix + [f'{chr(97 + level_idx)}{item + 1}'], level_idx + 1)

        create_rows([], 0)
        headers = [chr(97 + i) for i in range(num_levels)] + ['value']
        df = pd.DataFrame(data, columns = headers)
        return df

    @classmethod
    def random_rad_treemap(cls, num_levels=3, num_top_level_items=10, items_range=(2,5), value_range=(1,100), 
        outlier_fraction=0.2, use_log=True, sig=100,
        r1=1, r2=2, a1=0, a2=360, points=200,
        default_sort=False, default_sort_override=True, default_sort_override_reversed=True,
        mode='smart', no_groups=False, rotate_deg=0, full=False, rectangular=False,
        data_only=False):

        groupers = [chr(97 + i) for i in range(num_levels)]
        df = cls.__generate_hierarchical_data(num_levels=num_levels, num_top_level_items=num_top_level_items, 
            items_range=items_range, value_range=value_range, outlier_fraction=outlier_fraction, use_log=use_log, sig=sig)
        
        if data_only:
            return df

        return cls(df=df, groupers=groupers, value='value', r1=r1, r2=r2, a1=a1, a2=a2, points=points,
            default_sort=default_sort, default_sort_override=default_sort_override, default_sort_override_reversed=default_sort_override_reversed,
            mode=mode, no_groups=no_groups, rotate_deg=rotate_deg, full=full, rectangular=rectangular)

    def __corner_dictionary(self, group, start, count, value=None):
        df_dict = None
        if value:
            df_dict = group.agg({value: 'sum'})
        else:
            df_dict = group.size()
        df_dict = df_dict.to_dict()
        if value:
            df_dict = df_dict[value]
        end = 0
        count_ = 0
        dfc_dict = group.size()
        dfc_dict = dfc_dict.to_dict()
        for key in df_dict:
            count_ = df_dict[key]
            int_count = dfc_dict[key]
            end = df_dict[key] + start
            df_dict[key] = []
            df_dict[key].append(end/count - start/count) # [0]
            df_dict[key].append(0.) # r1 -> v1 [1]
            df_dict[key].append(0.) # r2 -> v2 [2]
            df_dict[key].append(0.) # a1 -> h1 [3]
            df_dict[key].append(0.) # a2 -> h2 [4]
            df_dict[key].append(count_) # count or value [5]
            df_dict[key].append(int_count) # integer count [6]
            df_dict[key].append(0.) # level rank [7]
            df_dict[key].append(0.) # overall rank [8]
            start = end
        return df_dict

    def __next_corner_calc(self, r1, r2, a1, a2, area, orientation='v', fill_angle=360., rectangular=False):
        out = 0., 0., 0., 0.
        if rectangular:
            if orientation == 'h':
                r3 = r1+area/(a2-a1)
                out = r1, r3, a1, a2 
            else:
                a3 = a1+area/(r2-r1)
                out = r1, r2, a1, a3
        else:
            if orientation == 'h':
                f = (a2-a1)/fill_angle
                r3 = sqrt(((area/(f))+pi*r1**2)/pi) # inner radius
                out = r1, r3, a1, a2 
            else:
                a3 = (area*fill_angle)/(pi*r2**2-pi*r1**2) # inner angle
                out = r1, r2, a1, a1+a3
        return out

    def __polygon_compare(self, r1, r2, a1, a2, area, rectangular=False):
        r1v, r2v, a1v, a2v = self.__next_corner_calc(r1, r2, a1, a2, area, orientation='v', rectangular=rectangular)
        r1h, r2h, a1h, a2h = self.__next_corner_calc(r1, r2, a1, a2, area, orientation='h', rectangular=rectangular)
        if rectangular:
            hlenv = a2v-a1v
            vlenv = r2v-r1v
            vmax = max(hlenv,vlenv)
            hlenh = a2h-a1h
            vlenh = r2h-r1h
            hmax = max(hlenh,vlenh)
        else:
            arclenv = (pi*r2v*2)*((a2v-a1v)/360)
            radlenv = r2v-r1v
            vmax = max(arclenv,radlenv)
            arclenh = (pi*r2h*2)*((a2h-a1h)/360)
            radlenh = r2h-r1h
            hmax = max(arclenh,radlenh)
        return vmax, hmax, r1v, r2v, a1v, a2v, r1h, r2h, a1h, a2h

    def __rad_treemap(self, df, groupers, value, r1, r2, a1, a2, points, default_sort, 
        default_sort_override, default_sort_override_reversed, mode, no_groups, rotate_deg,
        rectangular=False):

        area = (pi*r2**2-pi*r1**2)*((a2-a1)/360.)
        if rectangular:
            area = (r2-r1)*(a2-a1)
        rad_tm = dp()

        #region algorithm

        #region initialize
        list_group = []
        for i in range(len(groupers)):
            group_i = groupers[0:i+1]
            df_group = df.groupby(group_i, sort=default_sort)
            total = len(df)
            if value:
                total = df[value].sum()
            df_dict = self.__corner_dictionary(df_group, 0., total, value=value)
            if default_sort_override:
                df_dict = dict(sorted(df_dict.items(), key=lambda item: item[1][0], reverse=default_sort_override_reversed))
            list_group.append(df_dict)
        loop = True
        overall_rank = 1
        for g in range(len(groupers)):
            r1o, r2o, a1o, a2o = r1, r2, a1, a2
            level_rank = 1
            for i in list_group[g]:
                list_group[g][i][7] = level_rank
                list_group[g][i][8] = overall_rank
                if loop:
                    ai = list_group[g][i][0]*area
                    vmax, hmax, r1v, r2v, a1v, a2v, r1h, r2h, a1h, a2h = self.__polygon_compare(r1o, r2o, a1o, a2o, ai, rectangular=rectangular)
                    if (vmax > hmax or mode == 'out') and mode != 'around' and mode != 'alternate' and mode != 'legend': # draw out by radius (h)
                        list_group[g][i][1], list_group[g][i][2], list_group[g][i][3], list_group[g][i][4] = r1h, r2h, a1h, a2h
                        r1o = r2h
                    else: # draw around by angle (v)
                        list_group[g][i][1], list_group[g][i][2], list_group[g][i][3], list_group[g][i][4] = r1v, r2v, a1v, a2v
                        a1o = a2v
                level_rank += 1
                overall_rank += 1
            if not no_groups:
                loop = False
        #endregion

        #region point grid new
        for i in range(1, len(groupers), 1):
            lg_copy = copy.deepcopy(list_group[i-1]) # store a copy to return after overwriting
            for j in list_group[i]:
                if not no_groups:
                    # parent box
                    key = j[0:i]
                    if i == 1:
                        key = key[0]
                    rp1 = list_group[i-1][key][1]
                    rp2 = list_group[i-1][key][2]
                    ap1 = list_group[i-1][key][3]
                    ap2 = list_group[i-1][key][4]
                    # child box
                    ac = list_group[i][j][0]*area
                    vmax, hmax, r1v, r2v, a1v, a2v, r1h, r2h, a1h, a2h = self.__polygon_compare(rp1, rp2, ap1, ap2, ac, rectangular=rectangular)
                    if ((vmax > hmax or mode == 'out') and mode != 'around' and mode != 'alternate') or (mode == 'alternate' and i % 2 != 0): # draw out by radius (h)
                        list_group[i][j][1], list_group[i][j][2], list_group[i][j][3], list_group[i][j][4] = r1h, r2h, a1h, a2h
                        list_group[i-1][key][1] = r2h
                    else: # draw around by angle (v)
                        list_group[i][j][1], list_group[i][j][2], list_group[i][j][3], list_group[i][j][4] = r1v, r2v, a1v, a2v
                        list_group[i-1][key][3] = a2v
            list_group[i-1] = lg_copy # restore
        #endregion

        #region draw
        ix = 1
        level = 1
        for i in range(len(groupers)):
            for j in list_group[i]:
                # starting point calc
                x1, y1 = 0., 0.
                ad1, ad2 = list_group[i][j][3], list_group[i][j][4]
                a_1, a_2 = (ad1-90.)*pi/180., (ad2-90.)*pi/180.
                r_1, r_2 = list_group[i][j][1], list_group[i][j][2]
                if rectangular:
                    x1, y1 = ad1, r_1
                    angles1 = [ad1, ad2]
                    angles2 = [ad2, ad1]
                else:
                    x1, y1 = r_1*cos(a_1), r_1*sin(a_1)
                    # angles calc
                    point_frac = (ad2-ad1)/365.
                    points1 = int(point_frac*points)+5
                    points2 = int(point_frac*points*r_1/r_2)+5
                    angles1 = np.linspace(a_1, a_2, num=points1)
                    angles2 = np.linspace(a_2, a_1, num=points2)
                # values
                count, int_count = list_group[i][j][5], list_group[i][j][6]
                level_rank, overall_rank = list_group[i][j][7], list_group[i][j][8]
                # calc points
                x, y = 0., 0.
                if rectangular:
                    x, y = x1, y1
                else:
                    x, y = x1, -y1
                if vf.rotate != 0:
                    x, y = vf.rotate(x, y, rotate_deg)
                rad_tm.append(id=ix, x=x, y=y, path=0, level=level, group=j, count=int_count, value=count, level_rank=level_rank, overall_rank=overall_rank)
                for k in range(len(angles1)):
                    if rectangular:
                        x, y = angles1[k], r_2
                    else:
                        x, y = r_2*cos(angles1[k]), -r_2*sin(angles1[k])
                    if vf.rotate != 0:
                        x, y = vf.rotate(x, y, rotate_deg)
                    rad_tm.append(id=ix, x=x, y=y, path=k+1, level=level, group=j, count=int_count, value=count, level_rank=level_rank, overall_rank=overall_rank)
                for k in range(len(angles2)):
                    if rectangular:
                        x, y = angles2[k], r_1
                    else:
                        x, y = r_1*cos(angles2[k]), -r_1*sin(angles2[k])
                    if vf.rotate != 0:
                        x, y = vf.rotate(x, y, rotate_deg)
                    rad_tm.append(id=ix, x=x, y=y, path=len(angles1)+k+1, level=level, group=j, count=int_count, value=count, level_rank=level_rank, overall_rank=overall_rank)
                ix += 1
            level += 1
        #endregion

        #endregion

        rad_tm.to_dataframe()
        return rad_tm

    def __rad_treemap_legends(self, df, groupers, r1, a1, a2, points, chart_width, legend_width, buffer):
        dfs = []
        dfi = self.__rad_treemap(df, groupers, None, r1, r1+chart_width, 
            a1, a2, points, mode='legend', no_groups=False, rotate_deg=0, rectangular=self.rectangular)
        for i in range(len(groupers)):
            if i > 0:
                df_i = dfi.copy(deep=True)
                df_i[['x', 'y']] = df_i.apply(lambda row: vf.radial_projection(row['x'], row['y'],
                    (chart_width + legend_width + buffer * 2) * (i)), axis=1, result_type='expand')
                df_i['type'] = 'chart_' + str(i+1)
                dfs.append(df_i.loc[df_i['level'].isin([i,i+1])])
            dfl = self.__rad_treemap(df, groupers, None, r1+chart_width+buffer, r1+chart_width+buffer+legend_width,
                a1, a2, points, mode='legend', no_groups=True, rotate_deg=0, rectangular=self.rectangular)
            dfl['type'] = 'legend_' + str(i+1)
            dfs.append(dfl.loc[dfl['level'] == i+1])
            r1 += chart_width + legend_width + buffer * 2
        dfi['type'] = 'chart_0'
        dfs.append(dfi.loc[dfi['level'] == len(groupers)])
        df_out = pd.concat(dfs, axis=0)
        return df_out

    def __rad_tm_plot(self, pie_tree_df, level, opacity, line_level, mult, show=True):
        fig, axs = plt.subplots()
        axs.set_aspect('equal', adjustable='box')
        df_lvl = pie_tree_df.loc[pie_tree_df['level'] == level]
        df_lvl_group = df_lvl.groupby(['group'])
        for group, rows in df_lvl_group:
            x = rows['x'].values
            y = rows['y'].values
            r = random.random()
            b = random.random()
            g = random.random()
            color = (r, g, b)
            axs.fill(x, y, alpha=opacity, fc=color)
            axs.plot(x, y, 'k-', linewidth=0.5)
        if line_level > 0:
            df_lvl = pie_tree_df.loc[pie_tree_df['level'] == line_level]
            df_lvl_group = df_lvl.groupby(['group'])
            for group, rows in df_lvl_group:
                x = rows['x'].values
                y = rows['y'].values
                axs.plot(x, y, 'k-', linewidth=3)
        elif line_level < 0:
            mult += level
            colors = [str(i / (mult-1)) for i in range(mult)]
            for i in range(level-1, 0, -1):
                df_lvl = pie_tree_df.loc[pie_tree_df['level'] == i]
                df_lvl_group = df_lvl.groupby(['group'])
                for group, rows in df_lvl_group:
                    x = rows['x'].values
                    y = rows['y'].values
                    axs.plot(x, y, color = colors[i-1], linewidth=3)
        if show:
            plt.show()
        return fig, axs

    def __rad_tm_plot2(self, pie_tree_df, level, opacity, title, axis, limit, fill, show=True):
        fig, axs = plt.subplots()
        axs.set_aspect('equal', adjustable='box')
        df_lvl = pie_tree_df.loc[pie_tree_df['level'] == level]
        df_lvl_group = df_lvl.groupby(['group'])
        val_max = pie_tree_df['value'].max()
        for group, rows in df_lvl_group:
            x = rows['x'].values
            y = rows['y'].values
            v = rows['value'].values
            color = 'w'
            if fill == 'gray':
                color = str(v[0]/val_max)
            if fill == 'random':
                r = random.random()
                b = random.random()
                g = random.random()
                color = (r, g, b)
            axs.fill(x, y, alpha=opacity, fc=color)
        for i in range(level):
            df_lvl = pie_tree_df.loc[pie_tree_df['level'] == i+1]
            df_lvl_group = df_lvl.groupby(['group'])
            lw = (level - i)*2.5
            colors = [str(i / (level)) for i in range(level+1)]
            for group, rows in df_lvl_group:
                x = rows['x'].values
                y = rows['y'].values
                if limit == i+1 or limit == 0:
                    if limit == i+1:
                        axs.plot(x, y, color='black', linewidth=7)
                        axs.plot(x, y, color='grey', linewidth=4)
                    else:
                        axs.plot(x, y, color=colors[i], linewidth=lw)
        axs.set_title(title)
        axs.axis(axis)
        if show:
            plt.show()
        return fig, axs

    def rad_treemap(self, chart_width=1, legend_width=0.1, buffer=.01):
        if not self.full:
            self.o_rad_treemap = self.__rad_treemap(self.df, self.groupers, self.value, self.r1, self.r2, self.a1, self.a2,
                self.points, self.default_sort, self.default_sort_override, self.default_sort_override_reversed,
                self.mode, self.no_groups, self.rotate_deg, self.rectangular)
            self.df_rad_treemap = self.o_rad_treemap.df
        else:
            self.o_rad_treemap = self.__rad_treemap_legends(self.df, self.groupers, self.r1, self.a1, self.a2,
                self.points, chart_width=chart_width, legend_width=legend_width, buffer=buffer)
            self.df_rad_treemap = self.o_rad_treemap.df
            
    def plot_level(self, level=1, opacity=0.5, line_level=0, mult=0, show=True):
        fig, axs = self.__rad_tm_plot(self.df_rad_treemap, level=level, opacity=opacity, line_level=line_level, mult=mult, show=show)
        return fig, axs

    def plot_levels(self, level=1, opacity=0.5, title='', axis='on', limit=0, fill='gray', show=True):
        fig, axs = self.__rad_tm_plot2(self.df_rad_treemap, level=level, opacity=opacity, title=title, axis=axis, limit=limit, fill=fill, show=show)
        return fig, axs

    def to_df(self):
        return self.o_rad_treemap.df

    def to_csv(self, file_name):
        self.o_rad_treemap.dataframe_to_csv(file_name)
