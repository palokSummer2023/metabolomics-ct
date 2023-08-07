from utils import *


anom_f_df = {}
for dfn in ["cec_m1", "ser_m1", "cec_m2", "ser_m2"]:
    df = pd.read_hdf(f'generated/{struct_df_file[dfn]}', 'cvfiltered')
    e = Experiment2('con', 'con', df)
    S_0_7 = set(df[e.df0 != 0]['Metabolite Name']) ^ set(df[e.df7 != 0]['Metabolite Name'])
    S_7_15 = set(df[e.df7 != 0]['Metabolite Name']) ^ set(df[e.df15 != 0]['Metabolite Name'])
    S_0_15 = set(df[e.df0 != 0]['Metabolite Name']) ^ set(df[e.df15 != 0]['Metabolite Name'])
    S = S_0_7 | S_7_15 | S_0_15
    # print(f'{dfn} has {len(S_0_7)}, {len(S_0_15)}, {len(S_7_15)} differences. total: {len(S)}')
    anom_f_df[dfn] = list(S)

cldf = {}
for dfn in ["cec_m1", "ser_m1", "cec_m2", "ser_m2"]:
    df = pd.read_hdf(f'generated/{struct_df_file[dfn]}', 'cvfiltered')
    cldf[dfn] = df[~df['Metabolite Name'].isin(anom_f_df[dfn])]

differences = {}
for dfn in ["cec_m1", "ser_m1", "cec_m2", "ser_m2"]:
    e0 = Experiment2('con', 'con', cldf[dfn])
    mb0 = set(e0.df[(e0.df0 != 0)]['Metabolite Name'])
    sets_mb = {}
    for j, i in enumerate(['dss', 'ab', 'dab']):
        e = Experiment2(i, 'h2o', cldf[dfn])
        sets_mb[i] = set(e.df[(e.df7 != 0) & (e.df15 != 0)]['Metabolite Name'])
    diff = {k : v ^ mb0 for k, v in sets_mb.items()}
    differences[dfn] = diff