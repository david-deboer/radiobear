#! /usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import copy
import sqlite3 as sqlite

class Look():
    def __init__(self,fn='resp.db',setName='Set1HI',**par):
        """Look is a class that reads in the sql database and generates the plots.
           Right now, there is only set1 (the first run of the 2000+ models provided)
                      but in two tables (low freq (Set1) and hi freq (Set1HI))
           Parameters are:
              use_T:  which T value to use (only one).  For set1, these are a polar limb profile at
                             0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9 (use_T=0...6] and disc ave [7]
              use_limb:  which T values to use for limb plot.
              use_model = list of numbered models in database, alphabetical by model name
              use_freq = list of frequency numbers.  
                              For set1:  1.46, 2.44, 3.04, 3.44, 4.42, 5.45,
                              5.9, 6.45, 7.45, 8.58, 9.61, 10.38, 11.41, 13.18, 14.21, 15.18, 16.21,
                              17.38, 22.45, 23.45, 24.45, 25.45 [0...21]
                              For Set1HI:  29.98, 41.81, 58.31, 81.31, 113.40, 158.14, 220.54, 307.56, 428.92, 598.16
              switch_to_one_plot_at = don't split them out if below this number of plots"""
        
        if fn is not None:  ###If given a name, this will read in the database on creating class instance
            self.readData(fn,setName)
            self.fn = fn
        self.allowedPar = {'use_T':int,  ###These are the allowed parameters in the plot outputs and their types
                           'use_model':list,
                           'use_freq':list,
                           'use_limb':list,
                           'switch_to_one_plot_at':int}
        self.setpar(use_T=7, use_limb=range(0,7), switch_to_one_plot_at=8)  ###Set some defaults
        try:
            numModels = len(self.model)   ###self.model contains all of the model names
            srn = range(numModels)
            self.setpar(use_model=srn)   ###This selects all models
        except AttributeError:
            numModels = 10
        try:
            frn = range(len(self.freq))   ###self.freq contains all of the frequencies
            self.setpar(use_freq=frn)     ###This selects all frequencies
        except AttributeError:
            pass
        self.setpar(**par)    ###This sets any other par you've fed it
        self.axesRange = [0,numModels,110,450]
        
    def readData(self,fn,setName):
        """Read in data and set arrays.  For the specified set it populates:
                 self.model:  names of models e.g. am07tt1nc2cd3gg3
                 self.freq:   list of frequencies
                 self.Tdata:  T data for that model/freq -- limb and disc ave"""
        db = sqlite.connect(fn)   ###This just reads in the sql and makes the arrays described above.
        self.set = setName
        c = db.cursor()
        c_exec = "PRAGMA TABLE_INFO('%s')" % (self.set)
        c.execute(c_exec)
        self.col_headers = []
        for ch in c.fetchall():
            self.col_headers.append(str(ch[1]))
        c_exec = 'SELECT * from '+self.set
        c.execute(c_exec)
        data = c.fetchall()
        model = []
        freq = []
        for mf in data:
            mn = str(mf[0])
            if mn not in model:
                model.append(mn)
            if mf[1] not in freq:
                freq.append(mf[1])
        Tdata = []
        ctr = 0
        for ii, mn in enumerate(model):
            Tdata.append([])
            for jj, f in enumerate(freq):
                Tdata[ii].append(data[ctr][2:])
                ctr+=1
        self.Tdata = np.array(Tdata)
        self.model = model
        self.freq = np.array(freq)
        del(data)
    def setpar(self, show=False, **par):
        """Sets the parameters used for plotting"""
        if show is True:
            print 'Parameters:'
            for pk in self.allowedPar.keys():
                try:
                    v = getattr(self,pk)
                    if type(v)==list and len(v)>10:
                        v = '[Many...]'
                except AttributeError:
                    v = 'No value'
                print '\t'+str(pk)+'  ',self.allowedPar[pk],v
            return
        for pk in par.keys():
            if pk in self.allowedPar.keys() and type(par[pk])==self.allowedPar[pk]:
                setattr(self,pk,par[pk])
                #print 'Set ',pk
            elif pk in self.allowedPar.keys() and self.allowedPar[pk] == list:
                setattr(self,pk,[par[pk]])
            else:
                print 'Parameter error: ',pk,par[pk],
        try:
            if type(self.use_model[0]) == str:
                self.changeModName2Num()
        except AttributeError:
            self.use_model=None
    def changeModName2Num(self):
        nums = []
        for umn in self.use_model:
            nums.append(self.model.index(umn)) #since each model should only be there once...
        self.use_model=nums
    def info(self,setName='Set1HI'):
        """Reads and displays the sql database information"""
        db = sqlite.connect(self.fn)
        c = db.cursor()
        c.execute("PRAGMA TABLE_INFO('SetInfo')")
        prag = c.fetchall()
        c.execute("SELECT * FROM SetInfo")
        sii = c.fetchall()
        for ii,p in enumerate(prag):
            print '============'+p[1]+'============'
            print sii[0][ii]            

    def plotSpectra(self,**par):
        """Plots the spectra based on parameter values"""
        self.setpar(**par)
        indiv_plots = True
        if len(self.use_model)>self.switch_to_one_plot_at:
            plt.figure('Spectra')
            indiv_plots = False
        for umn in self.use_model:
            spectrum = self.Tdata[umn,:,self.use_T]
            if indiv_plots:
                plt.figure(self.model[umn])
            plt.plot(self.freq,spectrum)
        plt.xlabel('Freq [GHz]')
        plt.ylabel(r'$T_B$ [K]')

    def plotLimb(self,**par):
        """Plots limb(s) for given models and frequencies"""
        self.setpar(**par)
        plt.figure('Limb')
        for umn in self.use_model:
            for ufn in self.use_freq:
                limb = self.Tdata[umn,ufn,self.use_limb]
                plt.plot(self.use_limb,limb)
        plt.xlabel('b (index)')
        plt.ylabel(r'T$_B$')

    def plotModel(self,**par):
        """Set function below for 'default' plotModel version"""
        self.plotModel_NoSort(**par)

    def plotModel_NoSort(self,**par):
        """Plots the models based on parameter values, does not sort by value (so essentially alphabetical)"""
        self.setpar(**par)
        plt.figure('Raw')
        for f in self.use_freq:
            models = self.Tdata[:,f,self.use_T]
            g = '%.2f GHz' % (self.freq[f])
            plt.plot(models,label=g)
            x = 0.8*len(self.Tdata)/len(self.freq)*f
            plt.text(x,self.Tdata[int(x)][f][self.use_T],g)
        plt.title('db order')
        plt.xlabel('Model #')
        plt.ylabel(r'$T_B$ [K]')
        plt.axis(self.axesRange)

    def plotModel_AllSort(self,**par):
        """Plots..., sorts all spectra to be monotonically increasing in T"""
        self.setpar(**par)
        #Get individually sorted model results, keyed on freq
        plt.figure('T_sort_individually')
        for f in self.use_freq:
            T_sort_individually = copy.deepcopy(self.Tdata[:,f,self.use_T])
            T_sort_individually.sort()
            g = '%.2f GHz' % (self.freq[f])
            plt.plot(T_sort_individually,label=g)
            x = 0.8*len(self.Tdata)/len(self.freq)*f
            plt.text(x,T_sort_individually[int(x)],g)
        plt.title('Sort all')
        plt.xlabel('Model #')
        plt.ylabel(r'$T_B$ [K]')
        plt.axis(self.axesRange)

    def plotModel_IndivSort(self,**par):
        """Lots of plots.  For each used freq, it sorts on that and keeps that sort order for the rest"""
        self.setpar(**par)
        #Sort one-by-one, keyed on freq
        subplotprefix='22'
        sbn = 0
        for f in self.use_freq:
            TsortVals = []
            for i,ts in enumerate(self.Tdata[:,f,self.use_T]):
                TsortVals.append((ts,i))
            Tsortd = sorted(TsortVals,key=lambda x: x[0])
            if sbn == 0:
                plt.figure(str(self.freq[f])+'+')
            plt.subplot(subplotprefix+str(sbn+1))
            sbn=(sbn+1) % 4
            Tsort_on_one = []
            for iif in self.use_freq:
                Tsort_on_one.append([])
                for sv in Tsortd:
                    Tsort_on_one[-1].append(self.Tdata[sv[1]][iif][self.use_T])
                g = '%.2f GHz' % (self.freq[iif])
                plt.plot(Tsort_on_one[iif],label=g)
                x = 0.8*len(self.Tdata)/len(self.freq)*iif
                plt.text(x,Tsort_on_one[iif][int(x)],g)
            plt.title('Sort on '+str(self.freq[f])+' GHz')
            plt.xlabel('Model #')
            plt.ylabel(r'$T_B$ [K]')
            plt.axis(self.axesRange)

###Script that runs to make the plots
if __name__ == "__main__":
    l = Look('resp.db','Set1HI')
    l.info()
    l.plotSpectra()
    l.setpar(use_freq=range(len(l.freq)))
    l.plotModel_NoSort()
    l.plotModel_AllSort()
    l.plotModel_IndivSort()
    plt.show()

