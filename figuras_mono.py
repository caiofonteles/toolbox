# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 08:42:20 2017

@author: Ambidados
"""
import matplotlib.pyplot as plt
#%%
fig = plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(re_sst,'-k')
ax1.plot(pd_sst,'-b')
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel('Temperatura (Celsius)');
ax1.set_title('Séries de TSM')