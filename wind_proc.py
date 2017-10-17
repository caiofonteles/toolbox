# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:47:08 2017

@author: Ambidados
"""

import glob
import pandas as pd
import datetime
import pessoal
import numpy as np
#path =r'C:\Users\Ambidados\Documents\anem\dados' # use your path
allFiles = glob.glob("dados/rajada_2016_*.txt")
list_ = []; frame = pd.DataFrame()
for file_ in allFiles:
    df = pd.read_csv(file_, sep='\t+', header=0)
    df.columns = ['data','Spd','Dir']
    df.data = pd.to_datetime(df.data[:],format ='%d/%m/%Y %H:%M:%S')
    frame = frame.append(df)
    
    
#frame.index = range(0,len(frame.data))
#%%
frame.index = frame.data
U = np.cos(frame.Dir*np.pi/180)*frame.Spd
V = np.sin(frame.Dir*np.pi/180)*frame.Spd
U_h =  U.groupby(pd.TimeGrouper('H')).mean()
V_h =  V.groupby(pd.TimeGrouper('H')).mean()

frame.index= frame.data
temp = frame.groupby(pd.TimeGrouper('H')).mean()
x = CompToDirVel(U_h,V_h)
temp.Spd = x['vel']
temp.Dir = x['dir']
raj_h = temp
temp.index.names =  ['Data']
temp = temp.rename(columns={'Spd':u'Velocidade (m/s)','Dir':u'Direção (°)'})
temp.to_excel('rajada_2016_h.xlsx')

U_h =  U.groupby(pd.TimeGrouper('900S')).mean()
V_h =  V.groupby(pd.TimeGrouper('900S')).mean()

frame.index= frame.data
temp = frame.groupby(pd.TimeGrouper('900S')).mean()
x = CompToDirVel(U_h,V_h)
temp.Spd = x['vel']
temp.Dir = x['dir']
raj_h = temp
temp.index.names =  ['Data']
temp = temp.rename(columns={'Spd':u'Velocidade (m/s)','Dir':u'Direção (°)'})
temp.to_excel('rajada_2016_15m.xlsx')
#%%  vento TPM 
allFiles = glob.glob("dados/vento_2016_*.txt") # listar nomes dos arquivos
list_ = []; frame = pd.DataFrame()
for file_ in allFiles:
    df = pd.read_csv(file_, sep='\t+', header=0)
    df.columns = ['data','Spd','Dir']
    df.data = pd.to_datetime(df.data[:],format ='%d/%m/%Y %H:%M:%S')
    frame = frame.append(df)

frame.index= frame.data
#tmp1 = [(str(aux)).replace(',','.').replace('9.9.2','9.') for aux in frame.loc[:,'Spd'] ]
tmp = frame.Spd.convert_objects(convert_numeric=True)
#frame.Spd =pd.to_numeric(frame.Spd)  # mostra o que não consegue ser convertido

frame.Spd = tmp
#%%
U = np.cos(frame.Dir*np.pi/180)*frame.Spd
V = np.sin(frame.Dir*np.pi/180)*frame.Spd
U_h =  U.groupby(pd.TimeGrouper('H')).mean()
V_h =  V.groupby(pd.TimeGrouper('H')).mean()

temp = frame.groupby(pd.TimeGrouper('H')).mean()
x = CompToDirVel(U_h,V_h)
temp.Spd = x['vel']
temp.Dir = x['dir']
vento_h = temp
temp.index.names =  ['Data']
temp = temp.rename(columns={'Spd':u'Velocidade (m/s)','Dir':u'Direção (°)'})
temp.to_excel('vento_2016_h.xlsx')

U_h =  U.groupby(pd.TimeGrouper('900S')).mean()
V_h =  V.groupby(pd.TimeGrouper('900S')).mean()

temp = frame.groupby(pd.TimeGrouper('900s')).mean()
x = CompToDirVel(U_h,V_h)
temp.Spd = x['vel']
temp.Dir = x['dir']
vento_h = temp
temp.index.names =  ['Data']
temp = temp.rename(columns={'Spd':u'Velocidade (m/s)','Dir':u'Direção (°)'})
temp.to_excel('vento_2016_15m.xlsx')
#%% teste escrevendo hora
HH = temp.index.strftime('%H:%M:%S')
DT = temp.index.strftime('%Y/%m/%d')
temp = temp.assign(Data=DT,Hora=HH)
temp.to_excel('vento_2016_15m_teste.xlsx',columns=['Hora',u'Velocidade (m/s)',u'Direção (°)'])
