# -*- coding: utf-8 -*-
"""
Created on Mon Oct 02 15:28:04 2017

@author: Ambidados
"""

y=a.ppt[a.sst<22.5]
y[y==0].count()
y[y>0].count()
y=a.ppt[a.sst>25]
y[y==0].count()
y[y>0].count()
y=ppt_mc
y[y==0].count()
y[y>0].count()

results, edges = np.histogram(np.array(ppt_n.ppt[ppt_n.ppt[:]>=0]),bins=range(0,150,10), normed=1)
binWidth = edges[1] - edges[0]
plt.bar(edges[:-1], results*binWidth, binWidth)
