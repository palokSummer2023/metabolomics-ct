import pandas as pd
import numpy as np
from matplotlib import pyplot


struct_df_file = {
    'cec_m1' : 'a_cecal_mouse_1.hdf5',
    'cec_m2' : 'a_cecal_mouse_2.hdf5',
    'ser_m1' : 'a_serum_mouse_1.hdf5',
    'ser_m2' : 'a_serum_mouse_2.hdf5'
}


class Experiment:
    def __init__(self, d7, d15, df, type=None, dfname=None, mouse=1):
        self.d0 = "con"
        self.d7 = d7
        self.d15 = d15

        self.df = df

        self.df0 = self.df[f'{type}_0_{self.d0}_mouse_{mouse}'].astype(float)
        self.df7 = self.df[f'{type}_7_{self.d7}_mouse_{mouse}'].astype(float)
        if self.d7 == "con" and self.d15 == "con":
            self.df15 = self.df[f'{type}_15_con_mouse_{mouse}'].astype(float)
        else:
            self.df15 = self.df[f'{type}_15_{self.d7}_{self.d15}_mouse_{mouse}'].astype(float)
        self.shortname = f"E_{self.d7}_{self.d15}"
        self.name = f"E_{self.shortname}___{dfname}__mouse{mouse}"
        self.type = type
        self.cat = dfname

    def getmstability(self):
        if np.mean(self.df0) <= np.mean(self.df7) and np.mean(self.df7) <= np.mean(self.df15):
            return False
        if np.mean(self.df0) >= np.mean(self.df7) and np.mean(self.df7) >= np.mean(self.df15):
            return False
        return True

    def getoutlierstability(self):
        outliers = []
        for i in range(len(self.df0)):
            l = self.df0.iloc[i] <= self.df7.iloc[i] and self.df7.iloc[i] <= self.df15.iloc[i]
            h = self.df0.iloc[i] >= self.df7.iloc[i] and self.df7.iloc[i] >= self.df15.iloc[i]
            if (l or h) != self.getmstability():
                outliers.append(self.df['Metabolite Name'].iloc[i])
        return outliers, len(outliers)/len(self.df0)

    def getconvexityavg(self):
        jensen = (np.mean(self.df15) - np.mean(self.df7*2) + np.mean(self.df0))/(np.mean(self.df15) - np.mean(self.df0))
        return jensen
    
    def getconvexity(self):
        jensen = (self.df15 - self.df7*2 + self.df0)/(self.df15 - self.df0)
        return jensen

    def getconvexitynum(self):
        c = self.getconvexity()
        c = np.array(c.dropna())
        cmask = np.isinf(c)
        c = c[~cmask]
        return np.mean(c)

    def plot(self):
        pyplot.plot(["0", "7", "15"], [np.mean(self.df0), np.mean(self.df7), np.mean(self.df15)])
        pyplot.title(self.name)
        pyplot.show()

    def getmetabolites(self):
        absent_mb = (self.df0 == 0) & (self.df7 == 0) & (self.df15 == 0)
        return self.df[~absent_mb]['Metabolite Name']

class ExperimentException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class Experiment2(Experiment):
    def __init__(self, d7, d15, df, type=None, dfname=None, mouse=None):
        self.d0 = "con"
        self.d7 = d7
        self.d15 = d15

        self.df = df

        try:
            self.df0 = self.df[f'0_{self.d0}'].astype(float)
            p0 = f'0_{self.d0}'
            self.df7 = self.df[f'7_{self.d7}'].astype(float)
            p7 = f'7_{self.d7}'
            if self.d7 == "con" and self.d15 == "con":
                self.df15 = self.df[f'15_con'].astype(float)
                p15 = f'15_con'
            else:
                self.df15 = self.df[f'15_{self.d7}_{self.d15}'].astype(float)
                p15 = f'15_{self.d7}_{self.d15}'
            self.shortname = f"E_{self.d7}_{self.d15}"
            self.name = f"E_{self.shortname}___{dfname}__mouse{mouse}"
            self.type = type
            self.cat = dfname
            self.points = [p0, p7, p15]
        except:
            raise(ExperimentException(f'experiment does not exist with con, {self.d7}, {self.d15}'))

def exper_t(dfname, d7, d15):
        return Experiment(d7, d15, globals()[dfname], type=struct_types[dfname][0], mouse=struct_types[dfname][2], dfname=struct_types[dfname][1])





