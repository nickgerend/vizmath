# vizmath
# Author: Nick Gerend

import pandas as pd
import os
import matplotlib.pyplot as plt
from . import functions as vf

class points:

    def __init__(self):
            self.viz=[]
            self.df = pd.DataFrame()

    class point:

        def __init__(self, id=None, x=0.0, y=0.0, path=0, kwargs=None):
            self.id = id
            self.x = x
            self.y = y
            self.path = path
            for key, value in kwargs.items():
                self.__setattr__(key, value)
        
    def append(self, id, x, y, path, **kwargs):
        self.viz.append(self.point(id=id,x=x,y=y,path=path,kwargs=kwargs))

    def to_dataframe(self):
        self.df = pd.DataFrame([{attr: getattr(p,attr) for attr in dir(p) if not attr.startswith('_')} for p in self.viz])
    
    def dataframe_rescale(self, xmin=None, xmax=None, ymin=None, ymax=None, nxmin=-1, nxmax=1, nymin=-1, nymax=1):
        if not self.df.empty:
            if xmin is None:
                xmin = min(self.df['x'])
            if xmax is None:
                xmax = max(self.df['x'])
            if ymin is None:
                ymin = min(self.df['y'])
            if ymax is None:
                ymax = max(self.df['y'])
            self.df['x'] = [vf.rescale(x,xmin,xmax,nxmin,nxmax) for x in self.df['x']]
            self.df['y'] = [vf.rescale(y,ymin,ymax,nymin,nymax) for y in self.df['y']]

    def dataframe_to_csv(self, file_name, directory=None):
        if not self.df.empty:
            target_directory = directory if directory else os.getcwd()
            file_path = os.path.join(target_directory, file_name + '.csv')
            self.df.to_csv(file_path, encoding='utf-8', index=False)

    def plot_xy(self):
        x = [o.x for o in self.viz]
        y = [o.y for o in self.viz]
        plt.scatter(x, y)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()
