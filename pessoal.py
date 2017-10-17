# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:57:51 2017

@author: VPNUser
"""

import numpy as np
import os
from datetime import datetime
import time
import fnmatch

#import threading
#import subprocess

#==============================================================================
# 
#==============================================================================
def DirVelToComp(tdir,tvel):
    if type(tdir) is np.ndarray or list:
        dirRad = [np.deg2rad(tdir[ix]) for ix,aux in enumerate(tdir)]
        U = [np.cos(dirRad[ix])*tvel[ix] for ix,aux in enumerate(dirRad)]
        V = [np.sin(dirRad[ix])*tvel[ix] for ix,aux in enumerate(dirRad)]
        return dict(U=U,V=V)
            
    else:
        dirRad = np.deg2rad(tdir)
        U = np.cos(dirRad)*tvel
        V = np.sin(dirRad)*tvel
        return dict(U=U,V=V)
    pass
    

#import pessoal as p

#dirr = [random.random()*360 for aux in xrange(100)]
#vell = [random.random()*10 for aux in xrange(100)]
#
#comp = [p.DirVelToComp(dirr[ix],vell[ix]) for ix,aux in enumerate(vell)]
#comp[0]['U']
#comp[0]['V']

#ou 

#p.DirVelToComp(dirr,vell)

#%%
#==============================================================================
# 
#==============================================================================
def CompToDirVel(u,v):
    
    if type(u) is np.ndarray or list:
        print('np.ndarray or list')
        tdir = [np.rad2deg(np.arctan2(v[ix],u[ix])) for ix,aux in enumerate(u)]
        for ix,aux in enumerate(tdir):
            if tdir[ix] < 0:
                tdir[ix] += 360
            pass
        tvel = [pow(pow(u[ix],2) + pow(v[ix],2),0.5) for ix,aux in enumerate(u)]            
        return dict(dir=tdir,vel=tvel)
    
    else:
        print('float')
        tdir = np.rad2deg(np.arctan2(v,u)) #faz a soma vetoria e acha a direção resultante
        if tdir < 0:
            tdir += 360
            pass
        
        tvel = pow(pow(u,2) + pow(v,2),0.5) # Acha a velocidade vetorial resultante
        return dict(dir=tdir,vel=tvel)
        pass


#%%
#if type(np.asarray(U)) is np.ndarray:
#    print('ok')
#%%
#type(U[0])

#==============================================================================
#     
#==============================================================================
def estatSerie(velocidade,direcao): #entrar com lista d2, posVel = -2, posDir = -1 pra primeira camada.
    velMax = np.nanmax(velocidade)
    velMin = np.nanmin(velocidade)
    velMean = np.nanmean(velocidade)
    velStd = np.nanstd(velocidade)
         
    v = np.nanmean(np.sin(np.deg2rad(direcao)))
    u = np.nanmean(np.cos(np.deg2rad(direcao)))
    dirMean = np.arctan2(v,u)*(180/np.pi) 
    if dirMean < 0:
        dirMean += 360
        pass
    
    d = dict(velMax=velMax,
             velMin=velMin,
             velMean=velMean,
             velStd=velStd,
             dirMean=dirMean,           
            )
    
    return d 
    pass

    
def finder(path,what):
    
    dirr= []
    filee = []
    found = []
    try:
        for aux in os.listdir(path):
            if os.path.isdir(path+aux+'/'):
                dirr.append(path+aux+'/')
                fileTemp,dirTemp,foundTemp = finder(path+aux+'/',what)
                filee = filee + fileTemp
                dirr = dirr + dirTemp 
                found = found + foundTemp
                if (path+aux).find(what) >= 0:
                    print ('pasta: ' + aux)
                    found.append(path+aux)
                    pass
                pass
            elif os.path.isfile(path+aux):
                filee.append(path+aux)
                if (path+aux).find(what) >= 0:
                    datec = os.path.getctime(path+aux)
                    datec = datetime.fromtimestamp(datec)
                    datem = os.path.getmtime(path+aux)
                    datem = datetime.fromtimestamp(datem)
                    found.append([path+aux,datec,datem])
                    print ('arquivo: ' + path + aux)

                    pass
                pass
            else:
                pass
            pass
        found.sort(key = lambda temp: temp[2], reverse = True)
        
    except:
        print('acesso negado: ' + aux)
        pass
    
        
    return(filee,dirr,found)   


def finderEx(what,where):
    matches = []
    for root, dirnames, filenames in os.walk(where):
        for filename in fnmatch.filter(filenames, what):
            matches.append(os.path.join(root, filename))
    return matches



#%%

def finderPlus():
    import os
    from fnmatch import fnmatch
    
    d = []
    for root,dirr,files in os.walk(r"c:\Google Drive"):
        for aux in files:
            for aux2 in xrange(4,11,1):
                if fnmatch(aux,'*(%d)*'%aux2):
        #            print (os.path.join(root,aux))
                    d.append(os.path.join(root,aux))
                    pass
                pass
            pass
        pass
    
    import shutil as sh
    for aux in d:
        sh.move(aux,'c:/Google Drive/ZZZ/' + os.path.basename(aux))
        pass



#%%









if __name__ == '__main__':
    
    tem0 = finderEx('*ascapaf*','c:/google drive/')
    tem = [os.path.realpath(aux) for aux in tem0]
    



    
