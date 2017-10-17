# -*- coding: cp1252 -*-
"""
Created on Tue Sep 12 08:39:42 2017

@author: Ambidados
"""
import datetime
import numpy as np
import pandas as pd
#import pessoal as p
import matplotlib.pyplot as plt
import netCDF4 as nc
import glob
#import xarray as xr
#%% chuva macae
temp = pd.read_csv('chuva_macae.csv',sep=';',header=1)
temp.columns=['Data','Hora','ppt']

temp.replace(',','.')
temp['ppt']=temp['ppt'].astype('float64')
tmp=pd.to_datetime(temp['Data'])+pd.to_timedelta(temp['Hora'],unit='h')
temp.index=tmp
ppt_mc=temp['ppt'].groupby(pd.TimeGrouper('D')).sum()
ppt_mc = ppt_mc.to_frame()
#plt.plot(ppt_mc,'.')

#%%
temp = pd.read_csv('chuva_rj.csv',sep=';',header=1)
temp = temp.drop(temp.columns[[range(3,41)]],1)
temp.columns=['Data','Hora','ppt']
#tmp1 =temp.loc[:,'ppt'].values.astype(str).replace(',','.')
tmp1 = [(str(aux)).replace(',','.').replace('////','NaN') for aux in temp.loc[:,'ppt'] ]
temp['ppt'] = tmp1
temp['ppt']=temp['ppt'].astype('float64')
tmp=pd.to_datetime(temp['Data'])+pd.to_timedelta(temp['Hora'],unit='h')
temp.index=tmp
ppt_rj=temp['ppt'].groupby(pd.TimeGrouper('D')).sum() 
ppt_rj=ppt_rj.to_frame()


#%% importando e calculando médias de tsm
nc_list = glob.glob("*.nc")
i=0;
pd_sst=[];sst_date=[];re_sst =[]
for aux in nc_list[:]:
                Dset = nc.Dataset(aux)
                sst = Dset.variables["analysed_sst"][0,:,:] -273
                sst_r = Dset.variables["analysed_sst"][0,100:151,199:250] -273 #recortando a area de ressurgencia
                sst_units = Dset.variables["analysed_sst"].units
                sst_d = datetime.datetime.utcfromtimestamp(Dset.variables['time'][:]-32400+346809600)
                sst_date.append(sst_d)
                pd_sst.append(np.mean(sst))
                re_sst.append(np.mean(sst_r))
                i=i+1
                pass

#%%
df_sst=pd.DataFrame(-np.array(pd_sst)+np.array(re_sst))
df_sst.columns=['sst']     
df_sst.index=[sst_date]
pd_sst=pd.DataFrame(np.array(pd_sst))
pd_sst.columns=['sst']     
pd_sst.index=[sst_date]
re_sst=pd.DataFrame(np.array(re_sst))
re_sst.columns=['sst']     
re_sst.index=[sst_date]

#%% figuras ppt macae
fig = plt.figure()
ax1=fig.add_subplot(111)
l1, =plt.plot(ppt_mc,'-k',linewidth=0.6)
l2, =plt.plot(ppt_rj,'-b',linewidth=0.6)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Precipitação (mm)');
ax1.set_title(u'Precipitação Acumulada Diária'); ax1.legend([l1,l2],[u'Macaé','Rio de Janeiro'],loc='upper right')
ax1.set_ylim(0,220)
fig.savefig('figuras/series_ppt.png')
#%%
fig = plt.figure()
ax1=fig.add_subplot(211)
ax1.plot(ppt_mc,'-k',linewidth=0.6); ax1.set_title(u'Precipitação Acumulada Diária - Macaé');
ax1.set_ylim(0,220);  ax1.set_ylabel(u'Precipitação (mm)')
ax2=fig.add_subplot(212)
ax2.plot(ppt_rj,'-k',linewidth=0.6); ax2.set_title(u'Rio de Janeiro');
ax2.set_ylim(0,220)
ax2.set_xlabel('Tempo'); ax2.set_ylabel(u'Precipitação (mm)')
fig.savefig('figuras/series_ppt.png')
#%%
nc_list2 = pd.DataFrame(nc_list)
nc_list2.columns =  ['files']
nc_list2.index = re_sst.index
tst = nc_list2[re_sst.sst.between(22.71,27, inclusive=False)]
sst_3d=list(range(len(tst))); i=0
for aux in tst['files'][:]:
    #print aux
    Dset = nc.Dataset(aux)
    sst = Dset.variables["analysed_sst"][:]
    sst_3d[i]=np.array(sst-273)
    i=i+1
    pass

tt= np.array(sst_3d)
tt = tt.reshape(len(tst),300,401)
comp_ps = np.mean(tt,0)
comp_ps[np.where(comp_ps==-32768)]=np.nan
#%%
tst = nc_list2[re_sst.sst<21.62]
sst_3d=list(range(len(tst))); i=0
for aux in tst['files'][:]:
    #print aux
    Dset = nc.Dataset(aux)
    sst = Dset.variables["analysed_sst"][:]
    sst_3d[i]=np.array(sst-273)
    i=i+1
    pass

tt= np.array(sst_3d)
tt = tt.reshape(len(tst),300,401)
comp_ng = np.mean(tt,0)
comp_ng[np.where(comp_ng==-32768)]=np.nan
#run composites.py
#%% figura tsm costeira e offshore
fig = plt.figure()
ax1=fig.add_subplot(111)
l1, =plt.plot(re_sst,'-k',linewidth=0.6)
l2, =plt.plot(pd_sst,'-b',linewidth=0.6)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Temperatura (Celsius)');
ax1.set_title(u'Séries de TSM'); ax1.legend([l1,l2],['TSM costeira','TSM offshore'],loc='lower left')
ax1.set_ylim(18,30)
fig.savefig('figuras/series_tsm.png')

#%% figura de diferença entre a tsm costeira e offshore
fig = plt.figure()
ax1=fig.add_subplot(111)
l1, =plt.plot(df_sst,'-k',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Temperatura (Celsius)');
ax1.set_title(u'Diferença entre TSM (Costeira - Offshore)');
ax1.set_ylim(-4,4)
fig.savefig('figuras/diff_tsm.png')

#%% histogramas rio de janeiro

ppt_n = ppt_rj[red.sst<-2];ac_ng = ppt_n.sum();ac_ng = ac_ng/len(ppt_n) # media-1dp
ppt_p = ppt_rj[re_sst.sst>24.71];ac_ps = ppt_p.sum();ac_ps = ac_ps/len(ppt_p) # media+1dp
plt.figure();plt.plot(ppt_n,'.');plt.ylim(0,220);plt.title(u'Água fria')
plt.figure();plt.plot(ppt_p,'.');plt.ylim(0,220);plt.title(u'Água Quente')
#dm = ppt_mc.merge(df_sst, left_index=True, right_index=True, how='inner')
plt.figure();plt.hist(np.array(ppt_rj[ppt_rj.ppt[:]>0]),50,normed=1, facecolor='green'); plt.title(u'PPT - Macaé (2008-2016)');plt.ylim(0,0.7);plt.legend(ppt_mc.sum()/len(ppt_mc))
plt.savefig('figuras/hist_rj.png');plt.close()
plt.figure();plt.hist(np.array(ppt_n[ppt_n.ppt[:]>=0]), 50,normed=1, facecolor='green');plt.title(u'Água Fria');plt.ylim(0,0.7);plt.legend(ac_ng)
plt.savefig('figuras/hist_rj_anom_neg.png');plt.close()
plt.figure();plt.hist(np.array(ppt_p[ppt_p.ppt[:]>=0]), 50,normed=1, facecolor='green');plt.title(u'Água Quente');plt.ylim(0,0.7);plt.legend(ac_ps)
plt.savefig('figuras/hist_rj_anom_pos.png');plt.close()
#%% histogramas macae

ppt_n = ppt_mc[re_sst.sst<21.62];ac_ng = ppt_n.sum();ac_ng = ac_ng/len(ppt_n) # media-1dp
ppt_p = ppt_mc[re_sst.sst>24.71];ac_ps = ppt_p.sum();ac_ps = ac_ps/len(ppt_p) # media+1dp
#plt.figure();plt.plot(ppt_n,'.');plt.ylim(0,220);plt.title(u'Água fria')
#plt.figure();plt.plot(ppt_p,'.');plt.ylim(0,220);plt.title(u'Água Quente')
#dm = ppt_mc.merge(df_sst, left_index=True, right_index=True, how='inner')
plt.figure();plt.hist(np.array(ppt_mc[ppt_mc.ppt[:]>0]),50,normed=1, facecolor='green'); plt.title(u'PPT - Macaé (2008-2016)');plt.ylim(0,0.4);plt.legend(ppt_mc.sum()/len(ppt_mc))
plt.savefig('figuras/hist_macae.png');plt.close()
plt.figure();plt.hist(np.array(ppt_n[ppt_n.ppt[:]>=0]), 50,normed=1, facecolor='green');plt.title(u'Água Fria');plt.ylim(0,0.4);plt.legend(ac_ng)
plt.savefig('figuras/hist_mc_anom_neg.png');plt.close()
plt.figure();plt.hist(np.array(ppt_p[ppt_p.ppt[:]>=0]), 50,normed=1, facecolor='green');plt.title(u'Água Quente');plt.ylim(0,0.4);plt.legend(ac_ps)
plt.savefig('figuras/hist_mc_anom_pos.png');plt.close()

#%% analise de gradiente para encontrar os dias de ressurgencia

an_d = np.diff(re_sst.sst) # anomalia dia a dia
plt.figure();plt.plot(re_sst.index[1:],np.cumsum(an_d),'.-')
gd_sst = np.gradient(re_sst.sst)
plt.figure();plt.plot(re_sst.index[:],gd_sst,'.-')

#%%
import scipy.signal  as sg

b, a = sg.butter(4, 0.1, 'high')
y = sg.filtfilt(b, a, re_sst.sst, padlen=150)
plt.plot(re_sst.index[:],np.cumsum(y),'.-')

red =  pd.DataFrame(y) # teste para teste de spearman
red =  pd.DataFrame(np.cumsum(y))

red.index = re_sst.index
red.columns=['sst']

#%% filtro de passa baixa para mostrar variabilidade de longo periodo
import scipy.signal  as sg

b, a = sg.butter(4, 1/365., 'low')
y = sg.filtfilt(b, a, re_sst.sst, padlen=150)
plt.plot(re_sst.index[:],y,'.-')
a = ppt_mc.groupby(pd.TimeGrouper('W')).mean()

tsm_low =  pd.DataFrame(y) # teste para teste de spearman

tsm_low.index = re_sst.index
tsm_low.columns=['sst']
#%% plot da série de tsm diária e baixa frequencia
fig = plt.figure()
ax1=fig.add_subplot(111)
l1, =plt.plot(re_sst,'-k',linewidth=0.8)
l2, =plt.plot(tsm_low,'-b',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Temperatura (Celsius)');
ax1.set_title(u'Série de TSM');ax1.legend([l1,l2],[u'TSM diária',u'TSM baixa frequência'])
ax1.set_ylim(18,30);ax1.set_xlim('20071201','20170201')
fig.savefig('figuras/serie_tsm_res_low.png')
#%%
nc_list2 = pd.DataFrame(nc_list)
nc_list2.columns =  ['files']
nc_list2.index = re_sst.index
tst = a[a.sst>25]
sst_3d=list(range(len(tst))); i=0
for aux in tst['files'][:]:
    #print aux
    Dset = nc.Dataset(aux)
    sst = Dset.variables["analysed_sst"][:]
    sst_3d[i]=np.array(sst-273)
    i=i+1
    pass

tt= np.array(sst_3d)
tt = tt.reshape(len(tst),300,401)
comp_ps = np.mean(tt,0)
comp_ps[np.where(comp_ps==-32768)]=np.nan
#%%
tst = a[a.sst<22.5]
sst_3d=list(range(len(tst))); i=0
for aux in tst['files'][:]:
    #print aux
    Dset = nc.Dataset(aux)
    sst = Dset.variables["analysed_sst"][:]
    sst_3d[i]=np.array(sst-273)
    i=i+1
    pass

tt= np.array(sst_3d)
tt = tt.reshape(len(tst),300,401)
comp_ng = np.mean(tt,0)
comp_ng[np.where(comp_ng==-32768)]=np.nan
#run composites.py

#%% teste de spearman
import scipy.stats as est
import scipy.signal as sg

a = re_sst[(re_sst.index.month>10) | (re_sst.index.month<5)]


dm = ppt_mc.merge(a, left_index=True, right_index=True, how='inner')
dm  =  dm.dropna()
rank = est.spearmanr(dm.sst[:-15],dm.ppt[15:]); print rank
rank = sg.correlate(dm.sst[:],dm.ppt[:],mode='full'); print rank
plt.plot(dm.sst,dm.ppt,'.')

t_r =  red[red.sst>2]
plt.plot(dm.ppt,'.')

plt.figure();plt.hist(np.array(dm.ppt), 50,normed=1, facecolor='green'); plt.title(u'Ressurgência');plt.ylim(0,0.7);#plt.legend(ac_ng)
plt.savefig('figuras/hist_rj_anom_neg.png');plt.close()

#%% merging dados macae
dm = nc_list2.merge(re_sst, left_index=True, right_index=True, how='inner')
dm =  dm.merge(ppt_mc, left_index=True, right_index=True, how='inner')
a = dm[(dm.index.month>10) | (dm.index.month<5)]
a  =  a.dropna()
rank = est.spearmanr(a.sst[:-10],a.ppt[10:]); print rank

#%% histogramas macae

ppt_n = a[a.sst<22.5];ac_ng = ppt_n.ppt.sum();ac_ng = str(ac_ng/len(ppt_n)) # media-1dp
ppt_p = a[a.sst>25];ac_ps = ppt_p.ppt.sum();ac_ps = str(ac_ps/len(ppt_p)) # media+1dp
#plt.figure();plt.plot(ppt_n,'.');plt.ylim(0,220);plt.title(u'Água fria')
#plt.figure();plt.plot(ppt_p,'.');plt.ylim(0,220);plt.title(u'Água Quente')
#dm = ppt_mc.merge(df_sst, left_index=True, right_index=True, how='inner')
plt.figure();plt.hist(np.array(ppt_mc[ppt_mc.ppt[:]>=0]),bins=range(0,150),normed=1, facecolor='green'); plt.title(u'PPT - Macaé (2008-2016)');plt.ylim(0,1);plt.legend(ppt_mc.sum()/len(ppt_mc))
plt.savefig('figuras/hist_macae.png');plt.close()
plt.figure();plt.hist(np.array(ppt_n.ppt[ppt_n.ppt[:]>=0]), bins=range(0,150),normed=1, facecolor='green');plt.title(u'Ressurgência');plt.ylim(0,1);plt.legend([ac_ng])
plt.savefig('figuras/hist_mc_anom_neg.png');plt.close()
plt.figure();plt.hist(np.array(ppt_p.ppt[ppt_p.ppt[:]>=0]), bins=range(0,150),normed=1, facecolor='green');plt.title(u'Água Quente');plt.ylim(0,1);plt.legend([ac_ps])
plt.savefig('figuras/hist_mc_anom_pos.png');plt.close()
#%% figura de chuva e tsm durante ressurgencia e agua quente macae
fig = plt.figure()
ax1=fig.add_subplot(211)
l1, =plt.plot(re_sst,'-k',linewidth=0.8)
l1, =plt.plot(a.sst[a.sst<22.5],'.b',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Temperatura (Celsius)');
ax1.set_title(u'Série TSM');
ax1.set_ylim(18,30);ax1.set_xlim('20071201','20170201')
ax1=fig.add_subplot(212); 
#l1, =plt.plot(re_sst,'-k',linewidth=0.8)
l1, =plt.plot(a.ppt[a.sst<23],'.b',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Precipitação (mm)');
ax1.set_title(u'Precipitação durante ressurgência');
ax1.set_ylim(0,120); ax1.set_xlim('20071201','20170201')

fig.savefig('figuras/serie_mc_tsm_res.png')

fig = plt.figure()
ax1=fig.add_subplot(211)
l1, =plt.plot(re_sst,'-k',linewidth=0.8)
l1, =plt.plot(a.sst[a.sst>25],'.b',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Temperatura (Celsius)');
ax1.set_title(u'Série de TSM');
ax1.set_ylim(18,30);ax1.set_xlim('20071201','20170201')
ax1=fig.add_subplot(212); 
#l1, =plt.plot(re_sst,'-k',linewidth=0.8)
l1, =plt.plot(a.ppt[a.sst>25],'.b',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Precipitação (mm)');
ax1.set_title(u'Precipitação durante água quente');
ax1.set_ylim(0,120); ax1.set_xlim('20071201','20170201')

fig.savefig('figuras/serie_mc_tsm_pos.png')
#%% merging dados rj
dm = nc_list2.merge(re_sst, left_index=True, right_index=True, how='inner')
dm =  dm.merge(ppt_rj, left_index=True, right_index=True, how='inner')
a = dm[(dm.index.month>10) | (dm.index.month<5)]
a  =  a.dropna()
rank = est.spearmanr(a.sst[:-10],a.ppt[10:]); print rank

plt.plot(a.ppt[a.ppt>30],a.sst[a.ppt>30],'.')
plt.plot(a.ppt[(a.ppt>0.1) & (a.ppt<30)],a.sst[(a.ppt>0.1) & (a.ppt<30)],'.')
#%% histogramas rio de janeiro

ppt_n = a[a.sst<22.5];ac_ng = ppt_n.ppt.sum();ac_ng = ac_ng/len(ppt_n) # media-1dp
ppt_p = a[a.sst>25];ac_ps = ppt_p.ppt.sum();ac_ps = ac_ps/len(ppt_p) # media+1dp
#plt.figure();plt.plot(ppt_n,'.');plt.ylim(0,220);plt.title(u'Água fria')
#plt.figure();plt.plot(ppt_p,'.');plt.ylim(0,220);plt.title(u'Água Quente')
#dm = ppt_mc.merge(df_sst, left_index=True, right_index=True, how='inner')
plt.figure();plt.hist(np.array(ppt_rj[ppt_rj.ppt[:]>=0]),bins=range(0,150),normed=1, facecolor='green'); plt.title(u'PPT - Rio de Janeiro (2008-2016)');plt.ylim(0,1);plt.legend(ppt_rj.sum()/len(ppt_mc))
plt.savefig('figuras/hist_rj.png');plt.close()
plt.figure();plt.hist(np.array(ppt_n.ppt[ppt_n.ppt[:]>=0]), bins=range(0,150), normed=1,  facecolor='green');plt.title(u'Água Fria');plt.ylim(0,1);plt.legend([ac_ng])
plt.savefig('figuras/hist_rj_anom_neg.png');plt.close()
plt.figure();plt.hist(np.array(ppt_p.ppt[ppt_p.ppt[:]>=0]), bins=range(0,150), normed=1,  facecolor='green');plt.title(u'Água Quente');plt.ylim(0,1);plt.legend([ac_ps])
plt.savefig('figuras/hist_rj_anom_pos.png');plt.close()

#%% figura de chuva e tsm durante ressurgencia e agua quente RJ
fig = plt.figure()
ax1=fig.add_subplot(211)
l1, =plt.plot(re_sst,'-k',linewidth=0.8)
l1, =plt.plot(a.sst[a.sst<22.5],'.b',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Temperatura (Celsius)');
ax1.set_title(u'Série TSM');
ax1.set_ylim(18,30);ax1.set_xlim('20071201','20170201')
ax1=fig.add_subplot(212); 
#l1, =plt.plot(re_sst,'-k',linewidth=0.8)
l1, =plt.plot(a.ppt[a.sst<23],'.b',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Precipitação (mm)');
ax1.set_title(u'Precipitação durante ressurgência');
ax1.set_ylim(0,120); ax1.set_xlim('20071201','20170201')

fig.savefig('figuras/serie_rj_tsm_res.png')

fig = plt.figure()
ax1=fig.add_subplot(211)
l1, =plt.plot(re_sst,'-k',linewidth=0.8)
l1, =plt.plot(a.sst[a.sst>25],'.b',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Temperatura (Celsius)');
ax1.set_title(u'Série de TSM');
ax1.set_ylim(18,30);ax1.set_xlim('20071201','20170201')
ax1=fig.add_subplot(212); 
#l1, =plt.plot(re_sst,'-k',linewidth=0.8)
l1, =plt.plot(a.ppt[a.sst>25],'.b',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Precipitação (mm)');
ax1.set_title(u'Precipitação durante água quente');
ax1.set_ylim(0,120); ax1.set_xlim('20071201','20170201')

fig.savefig('figuras/serie_rj_tsm_pos.png')
#%% figura tsm diaria e low frew + precipitação macae
dm = re_sst.merge(ppt_mc, left_index=True, right_index=True, how='inner')
tsm_mo = re_sst.groupby(pd.TimeGrouper('W')).mean()
tsm_mo = tsm_mo.dropna()
a = ppt_mc.groupby(pd.TimeGrouper('W')).mean()
a = a.dropna()
fig = plt.figure()
ax1=fig.add_subplot(211)
l1, =plt.plot(a.ppt/a.ppt.max()*6+26,'-k',linewidth=0.8)
l2, =plt.plot(tsm_mo,'-b',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Temperatura (Celsius)');
ax1.set_title(u'Série de TSM');ax1.legend([l1,l2],[u'TSM diária',u'TSM baixa frequência'])
ax1.set_ylim(18,36);ax1.set_xlim('20071201','20170201')

ax1=fig.add_subplot(212); 
l1, =plt.plot(a.ppt,'.-r',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Precipitação (mm)');
ax1.set_title(u'Precipitação');
ax1.set_ylim(0,1000); ax1.set_xlim('20071201','20170201')

#fig.savefig('figuras/serie_tsm_prec.png')
#%% filtro de passa baixa para mostrar variabilidade de longo periodo
import scipy.signal  as sg

b, a = sg.butter(4, 1/180., 'low')
y = sg.filtfilt(b, a, re_sst.sst, padlen=150)
b, a = sg.butter(4, 1/365., 'low')
y1 = sg.filtfilt(b, a, re_sst.sst, padlen=150)
#plt.plot(re_sst.index[:],y,'.-')
a = ppt_mc.groupby(pd.TimeGrouper('6M')).sum()
c = re_sst.groupby(pd.TimeGrouper('6M')).mean()

fig = plt.figure()
ax1=fig.add_subplot(211)
l1, =plt.plot(re_sst.index[:],y,'-k',linewidth=0.8)
l2, =plt.plot(re_sst.index[:],y1,'-b',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Temperatura (Celsius)');
ax1.set_title(u'Série de TSM');ax1.legend([l1,l2],[u'TSM diária',u'TSM baixa frequência'])
ax1.set_ylim(20,26);ax1.set_xlim('20071201','20170201')

ax1=fig.add_subplot(212); 
l1, =plt.plot(a.ppt,'.-r',linewidth=0.8)
ax1.set_xlabel('Tempo (Ano)');ax1.set_ylabel(u'Precipitação (mm)');
ax1.set_title(u'Precipitação');
ax1.set_ylim(0,1500); ax1.set_xlim('20071201','20170201')
fig.savefig('figuras/serie_tsm_prec_long.png')
