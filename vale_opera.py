# -*- coding: utf-8 -*-
"""
Created on Thu Oct 05 14:50:18 2017

@author: Ambidados
"""

import glob
import pandas as pd
import datetime
#import pessoal
import numpy as np
import matplotlib.pyplot as plt
#%%
#path =r'C:\Users\Ambidados\Documents\anem\dados' # use your path
allFiles = glob.glob("*.csv")
list_ = []; frame = pd.DataFrame()
b10 = pd.read_csv(allFiles[0], sep=',', header=0)
#b10 = b10.drop(b10.index[len(b10)-1])
b10.Data = pd.to_datetime(b10.Data[:],format ='%Y-%m-%d %H:%M:%S')
b10.index = b10.Data
b10.Hs = b10.Hs.astype(float)  
#b10.loc[lambda x : x.loc[:,'Hs']>3.5,'Hs'] =  np.NaN
#tmp =  b10.loc[lambda x : x.loc[:,'Hs']>3.5,'Hs']
b10 = b10.Hs[b10.Hs<3.5]
b10.loc[lambda x : x.loc[:]==0]=  np.NaN
b10 = b10[(b10.index.year>2013)]
b10_2014 = b10[(b10.index.year==2014)].astype(float)
b10_2015 = b10[(b10.index.year==2015)].astype(float)
b10_2016 = b10[(b10.index.year==2016)].astype(float)
b10_2017 = b10[(b10.index.year==2017)].astype(float)
b10_2014_1h = pd.DataFrame();b10_2015_1h = pd.DataFrame();b10_2016_1h = pd.DataFrame();b10_2017_1h = pd.DataFrame()    
#%%
files  = [b10_2014,b10_2015,b10_2016,b10_2017]
files_n  = ['b10_2014','b10_2015','b10_2016','b10_2017']
files2 = [b10_2014_1h,b10_2015_1h,b10_2016_1h,b10_2017_1h]
for ix in range(0,4):
        files2[ix] = files[ix].groupby(pd.TimeGrouper('H')).mean()
        print files_n[ix]
        print files2[ix] .count()/float(len(files2[ix]))
        pass

#%%
#path =r'C:\Users\Ambidados\Documents\anem\dados' # use your path
cd ..
cd Boia4_historico/
#%%
allFiles = glob.glob("*.csv")
list_ = []; frame = pd.DataFrame()
b4 = pd.read_csv(allFiles[0], sep=',', header=0)
#b4 = b4.drop(b4.index[len(b4)-1])
b4.Data = pd.to_datetime(b4.Data[:],format ='%Y-%m-%d %H:%M:%S')
b4.index = b4.Data
b4  =  b4.Hs.astype(float)
    
b4.loc[lambda x : x.loc[:]>3.5]=  np.NaN
b4.loc[lambda x : x.loc[:]==0]=  np.NaN
b4_2014 = b4[(b4.index.year==2014)].astype(float)
b4_2015 = b4[(b4.index.year==2015)].astype(float)
b4_2016 = b4[(b4.index.year==2016)].astype(float)
b4_2017 = b4[(b4.index.year==2017)].astype(float)
b4_2014_1h = pd.DataFrame();b4_2015_1h = pd.DataFrame();b4_2016_1h = pd.DataFrame();b4_2017_1h = pd.DataFrame()    
#%%
files  = [b4_2014,b4_2015,b4_2016,b4_2017]
files_n  = ['b4_2014','b4_2015','b4_2016','b4_2017']
files2 = [b4_2014_1h,b4_2015_1h,b4_2016_1h,b4_2017_1h]
for ix in range(0,4):
        files2[ix] = files[ix].groupby(pd.TimeGrouper('H')).mean()
        print files_n[ix]
        print files2[ix] .count()/float(len(files2[ix]))
        pass
    
#%%
b4=pd.DataFrame(b4.dropna())
b10=pd.DataFrame(b10.dropna())
bt  =  b4.combine_first(b10)
bt.loc[lambda x : x.loc[:,'Hs']>3.5,'Hs']=  np.NaN
bt.loc[lambda x : x.loc[:,'Hs']==0,'Hs']=  np.NaN
bt_2014 = bt[(bt.index.year==2014)]['Hs'].astype(float)
bt_2015 = bt[(bt.index.year==2015)]['Hs'].astype(float)
bt_2016 = bt[(bt.index.year==2016)]['Hs'].astype(float)
bt_2017 = bt[(bt.index.year==2017)]['Hs'].astype(float)
bt_2014_1h = pd.DataFrame();bt_2015_1h = pd.DataFrame();bt_2016_1h = pd.DataFrame();bt_2017_1h = pd.DataFrame()    
#%%
files  = [bt_2014,bt_2015,bt_2016,bt_2017]
files_n  = ['bt_2014','bt_2015','bt_2016','bt_2017']
files2 = [bt_2014_1h,bt_2015_1h,bt_2016_1h,bt_2017_1h]
for ix in range(0,4):
        files2[ix] = files[ix].groupby(pd.TimeGrouper('H')).mean()
        print files_n[ix]
        print files2[ix] .count()/float(len(files2[ix]))
        pass

#%% figuras altura significativa
plt.style.use('seaborn-whitegrid')
st = datetime.datetime(2014, 01, 01)
en = datetime.datetime(2018, 01, 01)
fig = plt.figure(figsize=(12,9))
ax1=fig.add_subplot(311)
ax1.plot(b4, '.k',linewidth=0.4); ax1.set_title(u'Dados Bóia 4')
ax1.set_xlim(st,en);
ax1.set_ylim(0,3.5); ax1.set_ylabel('Altura Significativa (m)')
ax2=fig.add_subplot(312)
ax2.plot(b10, '.k',linewidth=0.4); ax2.set_title(u'Dados Bóia 10')
ax2.set_ylim(0,3.5); ax2.set_ylabel('Altura Significativa (m)')
ax2.set_xlim(st,en); 
ax3=fig.add_subplot(313)
ax3.plot(bt, '.k',linewidth=0.4); ax3.set_title(u'Dados Bóias 4 e 10 Unificados')
ax3.set_ylim(0,3.5); ax3.set_xlabel('Tempo'); ax3.set_ylabel('Altura Significativa (m)')
ax3.set_xlim(st,en)
fig.tight_layout()
fig.savefig('altura_sig.png')
