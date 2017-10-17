# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 13:10:25 2017

@author: Ambidados
"""

## figuras vale
import matplotlib as mpl
import matplotlib.pyplot as plt
import datetime
import fnmatch
tmp =  dir(ax1)
[aux for aux in tmp if fnmatch.fnmatch(aux,'*line*')]
#########################################
# figuras vento
#########################################
plt.style.use('seaborn-whitegrid')
st = datetime.datetime(2015, 01, 01)
en = datetime.datetime(2016, 01, 01)
fig = plt.figure()
ax1=fig.add_subplot(211)
ax1.plot(vento_h['Velocidade (m/s)'], '-k',linewidth=0.4); ax1.set_title('Vento TPM 2015')
ax1.set_xlim(st,en)
ax1.set_ylim(0,50); ax1.set_xlabel('Tempo'); ax1.set_ylabel('Velocidade (m/s)')
ax2=fig.add_subplot(212)
ax2.plot(vento_h['Direcao (graus)'], '-k',linewidth=0.4)
ax2.set_ylim(0,360)
ax2.set_xlabel('Tempo'); ax2.set_ylabel('Direcao (graus)')
ax2.set_xlim(st,en)
fig.savefig('figuras/vento_TPM_2015_1h.png')

#########################################
## figuras rajadas
#########################################
plt.style.use('seaborn-whitegrid')
st = datetime.datetime(2015, 01, 01)
en = datetime.datetime(2016, 01, 01)
fig = plt.figure()
ax1=fig.add_subplot(211)
ax1.plot(raj_h['Velocidade (m/s)'], '-k',linewidth=0.4); ax1.set_title('Vento TPM 2015')
ax1.set_xlim(st,en)
ax1.set_ylim(0,50); ax1.set_xlabel('Tempo'); ax1.set_ylabel('Velocidade (m/s)')
ax2=fig.add_subplot(212)
ax2.plot(raj_h['Direcao (graus)'], '-k',linewidth=0.4)
ax2.set_ylim(0,360)
ax2.set_xlabel('Tempo'); ax2.set_ylabel('Direcao (graus)')
ax2.set_xlim(st,en)
fig.savefig('figuras/rajada_TPM_2015_1h.png')

#%%  figuras mensais vento

plt.style.use('seaborn-whitegrid')
st = datetime.datetime(2015, 01, 01)
en = datetime.datetime(2016, 01, 01)
d_name=['Jan','Fev','Mar','Abr','Mai','Jun','Jul','Ago','Set','Out','Nov','Dez']
for aux in range(3,13):
            tmp_v=vento_h[dates[aux-1]:dates[aux]]
            fig = plt.figure()
            ax1=fig.add_subplot(211)
            ax1.plot(tmp_v['Velocidade (m/s)'], '-k',linewidth=0.4); ax1.set_title('Vento TPM - '+d_name[aux-1]+'/2015')
            ax1.set_xlim(dates[aux-1],dates[aux])
            ax1.set_ylim(0,50); ax1.set_xlabel('Tempo'); ax1.set_ylabel('Velocidade (m/s)')
            ax2=fig.add_subplot(212)
            ax2.plot(tmp_v['Direcao (graus)'], '-k',linewidth=0.4)
            ax2.set_ylim(0,360)
            ax2.set_xlabel('Tempo'); ax2.set_ylabel('Direcao (graus)')
            ax2.set_xlim(dates[aux-1],dates[aux])
            fig.savefig('figuras/vento_TPM_2015_1h_'+d_name[aux-1]+'.png')
            plt.close()
            pass
        
#%%  figuras mensais Rajada

plt.style.use('seaborn-whitegrid')
st = datetime.datetime(2015, 01, 01)
en = datetime.datetime(2016, 01, 01)
d_name=['Jan','Fev','Mar','Abr','Mai','Jun','Jul','Ago','Set','Out','Nov','Dez']
for aux in range(2,13):
            tmp_v=raj_h[dates[aux-1]:dates[aux]]
            fig = plt.figure()
            ax1=fig.add_subplot(211)
            ax1.plot(tmp_v['Velocidade (m/s)'], '-k',linewidth=0.4); ax1.set_title('Rajada TPM - '+d_name[aux-1]+'/2015')
            ax1.set_xlim(dates[aux-1],dates[aux])
            ax1.set_ylim(0,50); ax1.set_xlabel('Tempo'); ax1.set_ylabel('Velocidade (m/s)')
            ax2=fig.add_subplot(212)
            ax2.plot(tmp_v['Direcao (graus)'], '-k',linewidth=0.4)
            ax2.set_ylim(0,360)
            ax2.set_xlabel('Tempo'); ax2.set_ylabel('Direcao (graus)')
            ax2.set_xlim(dates[aux-1],dates[aux])
            fig.savefig('figuras/rajada_TPM_2015_1h_'+d_name[aux-1]+'.png')
            plt.close()
            pass