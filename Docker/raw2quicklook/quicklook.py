from scipy.fftpack import fft,fftshift,ifft  #VERSIONE APPLY CHE MODIFICA GLI ATTREIBUTI DELL'OGGETTO
from joblib import Parallel, delayed

def container(obj, dataBlock):
            return obj.multiRangeCompression(dataBlock)

# observations: the fastest way seem to be np.apply on axis
# convert data to np with padding if needed (after compression)
# use np apply inside the block crucher

class quicklook: 
    import logging 
    from contextlib import closing 
    import multiprocessing as mp 
    import pandas as pd 
    import numpy as np 
    from logging import raiseExceptions
    from scipy.constants import speed_of_light as c 
    import scipy.ndimage.morphology as scmorph 
    import math
    import sys
    import plotly as pl
    import plotly.express as px
    from plotly.subplots import make_subplots
    import os
    import rimeOST
    import pdfgen
    # from tqdm import tqdm
    # from  parallel_pandas import ParallelPandas
    # from multiprocesspandas import applyparallel 
    OPENING_MIN = 4 #SNR

    sel_pulse_width = {
        0: 50e-6, 
        1: 100e-6, 
        2: 150e-6, 
        3: 200e-6, 
        4: 250e-6, 
        5: 75e-6, 
        6: 125e-6, 
        7: 175e-6 
    }
    sel_sampling_freq = {
        0: 3e6,
        1: 1.2e6
    }
    sel_B = {
        0: 2.8e6,
        1: 1e6
    }
    modeDict = {
        'active': 1,
        'passive': 0
    }
    RxGain = -7 +83 -0.9 -5

    def setOutputFolder(self,outputFolder): 
        self.outputFolder = outputFolder

    def __init__(self,rime=None, mode=None, processors=1,exec='0|',ost=None): 
        # self.ParallelPandas.initialize(n_cpu=processors,split_factor=4,disable_pr_bar=True)
        self.exID = exec
        self.logger = self.logging.getLogger('quicklook')
        self.processors = processors

        if mode is None:
            self.mode = not rime.DataFrame['Rx'][0]
        else:
            self.mode = self.modeDict[mode]

        if rime is not None:
            # self.rime = rime
            self.rawData   = rime.DataFrame['data'].values
            self.bitNumber = rime.DataFrame['Bit_Num'][0]
            self.B           = self.sel_B[rime.DataFrame['B'][0]]
            self.samplingFr  = self.sel_sampling_freq[rime.DataFrame['B'][0]]
            self.N_OST = rime.DataFrame['N_OST'].to_numpy()
            self.OSTId = rime.DataFrame.iloc[0].OST_Id
            # self.freqAxis = self.np.flip(9e6+self.np.fft.fftshift(self.np.fft.fftfreq(len(self.rawData[0]),1/self.samplingFr)))

            """self.bitNumber = rime.DataFrame['Bit_Num']
            self.B           = self.np.asarray([self.sel_B[letter] for letter in rime.DataFrame['B']])
            self.samplingFr = self.np.asarray([self.sel_sampling_freq[letter] for letter in rime.DataFrame['B']])"""
            if self.mode: 
                self.SWST        = rime.DataFrame['SWST']
                self.pulse_width = self.sel_pulse_width[rime.DataFrame['Pulse_sel'][0]]
                self.wnd_samples = int(self.pulse_width*self.samplingFr) # they are not constat TODO
                #self.pulse_width = [self.sel_pulse_width[letter] for letter in rime.DataFrame['Pulse_sel']]
                #self.wnd_samples = (self.np.asarray(self.pulse_width)*self.np.asarray(self.samplingFr)).astype('int')
                self.rank = rime.DataFrame['Rank']
        else: 
            raise SystemError('Obsolete, RIMERaw object needed as input')

        #save ost 
        if ost is not None: 
            self.OST = ost
            self.attenuation = ost.get_attenuation()

        self.logger.debug(self.exID+'Quicklook generation: beginning')
        if self.mode: 
            #costruttore caso attivo 
            self.subcostructorPassive()
            self.subcostructorActive()

        else: 
            self.subcostructorPassive()


    def cruncher(self,input, output, function, step=None):
            if  step is None: 
                step = len(input)//self.processors
            if isinstance(input, list):
                in2 = input
            else:
                in2 = input.tolist()
            def funIterator(values):
                return [function(xx) for xx in values]
            temp = Parallel(n_jobs=self.processors, prefer='threads' )(delayed(funIterator)(in2[jj : jj+step])  for jj in range(0, len(input), step))
            output[:] = [i for elem in temp for i in elem]

    def cruncherDf(self,input, output, function, step=None):
            if  step is None: 
                step = len(input)//self.processors
            return Parallel(n_jobs=self.processors, prefer='threads' )(delayed(function)(input.loc[jj : jj+step])  for jj in range(0, len(input), step))
    
    def subcostructorPassive(self): 
        self.maxDim = self.np.max([xx.shape[0] for xx in self.rawData])
        self.power_ADC = self.ADCRMSPower() 
        self.power_Antenna = self.AntennaPower()
        self.RMS = self.RMSAmplitude()
        self.logger.debug(self.exID+'Power completed')

        self.saturation_values_real = []
        self.saturation_values_imag = []
        self.cruncher(self.rawData, self.saturation_values_real,self.saturationReal )
        self.cruncher(self.rawData, self.saturation_values_imag,self.saturationImag )
        self.saturation_values_real =  self.np.stack(self.saturation_values_real)
        self.saturation_values_imag =  self.np.stack(self.saturation_values_imag)
        self.logger.debug(self.exID+'Saturation completed')

        self.spectrogramValues = []
        self.cruncher(self.rawData, self.spectrogramValues,self.spectogram )
        self.calibratedSpectogramValues = []
        self.cruncher(self.rawData, self.calibratedSpectogramValues,self.calibratedSpectogram )
        self.logger.debug(self.exID+'Spectrogram completed')

        self.meanSpectrum = 20*self.np.log10(self.np.mean(self.np.abs(self.spectrogramValues)))
    
    def subcostructorActive(self): 
        self.compressed=[]
        self.cruncher(self.rawData, self.compressed,self.rangeCompressionRow )
        self.logger.debug(self.exID+'Range Compression completed')

        self.snr = []
        self.cruncher(self.compressed, self.snr,self.SNR )
        self.logger.debug(self.exID+'SNR completed')

        self.shifted = []
        self.cruncher(self.compressed, self.shifted,self.shiftFunction )
        self.logger.debug(self.exID+'shift completed')

        self.altitude = self.height(self.compressed)
        self.logger.debug(self.exID+'height completed')

        self.compressed = self.np.stack(self.compressed)

        self.powerColumns = self.powerClm() # TODO change to calibrated power
        self.logger.debug(self.exID+'power completed')

    def get_mask(self, rowId):
        return self.N_OST == rowId
    
    def rangeCompression(self,row): 
        # N_fft = row.count()
        N_fft = len(row.values)
        return self.pd.Series(fftshift(ifft(fft(row.fillna(0).values,N_fft)*self.np.conj(fft(row.fillna(0).values,N_fft)))))
    
    def rangeSelfCompressionRow(self,row):
        result = self.np.zeros(self.maxDim, dtype ='complex128')
        N_fft = len(row)
        temp = self.np.nan_to_num(row)
        result[:N_fft] = fftshift(ifft(fft(temp,N_fft)*self.np.conj(fft(temp,N_fft))))
        return result

    def rangeCompressionRow(self,row):
        result = self.np.zeros(self.maxDim, dtype ='complex128')
        N_fft = len(row)
        temp = self.np.nan_to_num(row)
        t = self.np.arange(start=-self.pulse_width/2, stop = (self.pulse_width/2)+(1/self.samplingFr),step=1/self.samplingFr)
        rate = self.B/self.pulse_width
        s_ref = self.np.exp(+1j*self.math.pi*rate*self.np.power(t,2))
        result[:N_fft] = ifft(fft(temp,N_fft)*fft(s_ref,N_fft))
        # result[:N_fft] = fftshift(ifft(fft(temp,N_fft)*fft(s_ref,N_fft)))
        return result

    def multiRangeCompression(self, values):
        # temp = self.np.nan_to_num(self.np.vstack(values))
        # N_fft = temp.shape[1]
        # return fftshift(ifft(fft(temp,N_fft)*self.np.conj(fft(temp,N_fft))))
        result = []
        for row in values:
            N_fft = len(row)
            temp = self.np.nan_to_num(row)
            result.append(fftshift(ifft(fft(temp,N_fft)*self.np.conj(fft(temp,N_fft)))))
        return result

    def multiRangeCompressionDf(self, df):
        return df.loc[:,'data'].apply(self.rangeCompressionRow)
    
    def power_adjusting_factor(self): 
        # Conversion factor from integer power to dBm
        full_scale = 4 #dBm TODO not correct, 4dBm is the average power for saturation, considering a gaussian signal the relative Vrms is std dev
        levels = 128
        adjusting_factor = full_scale-20*self.np.log10(levels) 
        return adjusting_factor
    
    def power(self,adjusting_factor): 
        # powerValue = (self.rawData.abs()).pow(2)
        # return (powerValue.mean(axis=1,skipna=True))
        # powerValue = 10*self.np.log10(self.np.abs(self.rawData/128)**2)

        powerValue = self.np.abs(self.rawData)**2
        result =  self.np.vectorize(self.np.nanmean)(powerValue)
        
        #dB + shift values
        result[result <= 0] = self.sys.float_info.epsilon
        result = 10*self.np.log10(result) + adjusting_factor
        return result
    
    def RMSAmplitude(self): 
        powerValue = self.np.abs(self.rawData)**2
        result =  self.np.vectorize(self.np.nanmean)(powerValue)
        result = self.np.sqrt(result)
        return result

    def ADCRMSVoltage(self):
        Vfs = 0.3544 # Peak to peak
        powerValue = self.np.abs(self.rawData * (Vfs) / 255)**2
        result =  self.np.vectorize(self.np.nanmean)(powerValue)
        result = self.np.sqrt(result)
        return result
    
    def ADCRMSPower(self):
        Z = 50 
        result = (self.ADCRMSVoltage()**2)/Z
        result[result <= 0] = self.sys.float_info.epsilon
        return result
    
    def AntennaPower(self):
        return self.power_ADC * 10**( (self.attenuation - self.RxGain)/10)

    def AntennaPowerOld(self, ostRowId = None):
        '''
        Calibrated power at antenna port in dBm
        '''
        mask = self.get_mask(ostRowId)
        antPower = 30 + 10*self.np.log10(self.power_ADC[mask])+ self.attenuation[mask] - self.RxGain
        return antPower

    def powerClm(self,adjusting_factor = 0): #TODO aggiungere adjusting_factor
        # powerValue = (self.rawData.abs()).pow(2)
        # return (powerValue.mean(axis=0,skipna=True))
        powerValue = self.np.abs(self.shifted)**2
        result =  self.np.nanmean(self.np.vstack(powerValue), axis=0)

        #dB + shift values
        result[result <= 0] = self.sys.float_info.epsilon
        result = 10*self.np.log10(result) 
        return result

    def powerData(self,data): 
        n_samples = self.np.sum(~self.np.isnan(data))
        total = self.np.nansum((abs(data))**2)
        return self.np.float32(total/n_samples)

    def SNR(self,data): #data = np.array
        LIMIT = 2*self.wnd_samples #TODO check that the window lenght could be different for each PRI
        max = self.np.nanargmax(data)

        if LIMIT <= max: 
            noise = (data[0:self.wnd_samples])
        else: 
            noise = (self.np.concatenate((data[-(LIMIT-max):],data[:(max-self.wnd_samples)])))
        
        if not self.np.any(noise):
            return 0 
            # self.logger.fatal(self.exID+'Error identifying the noise')
            # self.raiseExceptions('Unable to distinguish noise from signal')

        threshold = self.np.nanmean(noise)+3*self.np.nanstd(noise)
        noise_mask = self.np.abs(data) > threshold
        noise_mask = self.scmorph.binary_opening(noise_mask,structure=self.np.ones(self.OPENING_MIN,dtype=int))
        index = self.np.argwhere(noise_mask == True)
        noise = (self.np.abs(data))

        if len(index) != 0: 
            noise[int(self.np.nanmin(index)):int(self.np.nanmax(index))+1] = self.np.NaN

        SNR = ((abs(data[max]))**2) / self.powerData(noise)
        return self.np.float32(20*self.np.log10(SNR))

    def multiSNR(self, input):
        return [self.SNR(xx) for xx in input]

    def saturationReal(self,row): 
        hist,bins = self.np.histogram(self.np.real(row[~self.np.isnan(row)]), bins=self.np.arange(-128.5, 128.5, 1.0))
        return self.pd.Series(hist)

    def saturationImag(self,row): 
        hist,bins = self.np.histogram(self.np.imag(row[~self.np.isnan(row)]), bins=self.np.arange(-128.5, 128.5, 1.0))
        return self.pd.Series(hist)

    def spectogram(self,row):
        result = self.np.zeros(self.maxDim) 
        result[:len(row)] = self.np.abs(fftshift(fft(row,len(row))))
        return result   
   
    def instVoltage(self, row):
        Vfs = 0.3544
        return Vfs * row  / 255

    def calibratedSpectogram(self,row):
        # PowerSpectrum reported in W/Hz at ADC port
        # TODO check dfreq for each row and put back the real spectral power density
        # ref: https://www.sjsu.edu/people/burford.furman/docs/me120/FFT_tutorial_NI.pdf
        Z = 50
        # dfreq = self.freqAxis[1] - self.freqAxis[0]
        dfreq = 1 
        # dfreq = self.OST.get_dfreq(row) 
        result = self.np.zeros(self.maxDim)
        fftres = fft(self.instVoltage(row),len(row))
        result[:len(row)] = self.np.abs(fftshift( (fftres * self.np.conjugate(fftres))/len(row) )) / (dfreq *Z)
        return result   

    def height(self,row): 
        # time = ((max.divide(self.samplingFr)).add(self.SWST)).multiply(10**-6)
        max = self.np.argmax(row)
        time = (max / self.samplingFr) + self.SWST * 1e-6 + self.OST.get_PRI()*self.rank
        return self.np.float32(time / (2*self.c))
        # return self.np.float32((time.divide(self.c)).divide(2))

    def shiftFunction(self,row): 
        # dim = row.shape
        dim = (self.maxDim,)
        peak_position = dim[0] - 1 
        shift_amount = (peak_position - self.np.argmax(row))
        new_dim = dim
        new_dim = (2*dim[0],)
        result = self.np.empty_like(row, shape=new_dim)
        result[:shift_amount] = self.np.nan
        result[shift_amount:shift_amount+row.shape[0]] = row
        result[shift_amount+row.shape[0]:] = self.np.nan
        return result

    # def spectogram(self,row): 
    #     # row = row.to_numpy()
    #     return self.pd.Series(self.np.abs(fftshift(fft(row,len(row)))))


    #####PLOTS 

    def rowMeanPower_plot(self,silent=True,pdf=None): 
        fig = self.px.line(10*self.np.log10(self.power_Antenna)+30,
                           labels={'index':'Packet','value':'Power [dBm]'},
                           title='Calibrated Mean Power at antenna port [dBm]')
        fig.add_hline(y=-77.8,line_color = 'red')
        fig.add_hline(y=-83,line_color = 'red')
        fig.update_layout(showlegend = True)

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename= self.outputFolder+"/meanPowerAtAntenna.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)
    
    def rowRMSAmplitude_plot(self,silent=True,pdf=None): 
        fig = self.px.line(self.RMS,labels={'index':'Packet','value':'RMS Amplitude (DN)'},title='RMS Amplitude for each packet (DN)')

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename= self.outputFolder+"/RMSAmplitude.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def rowMeanPowerADC_plot(self,silent=True,pdf=None): 
        fig = self.px.line(10*self.np.log10(self.power_ADC) +30,
                        labels={'index':'Packet','value':'Power (dBm)'},
                        title='Mean Power for each packet at ADC port (dBm)')
        fig.add_hline(y=4,line_color = 'red')
        fig.update_layout(showlegend = True)

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename= self.outputFolder+"/meanPowerADC.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def bitSaturation(self,silent=True,pdf=None): 
        figI = self.px.imshow(self.saturation_values_imag.T,aspect='auto',labels=dict(x='Packet',y='Bit'),title='Bit saturation: imaginary part',y = self.np.arange(-128,128,1.0))
        figR = self.px.imshow(self.saturation_values_real.T,aspect='auto',labels=dict(x='Packet',y='Bit'),title='Bit saturation: real part',y = self.np.arange(-128,128,1.0))

        fig = self.pl.subplots.make_subplots(rows=1, cols=2,subplot_titles=('Bit saturation: imaginary part','Bit saturation: real part'))
        fig.add_trace(figR.data[0],row=1,col=1)
        fig.add_trace(figI.data[0],row=1,col=2)

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/saturation.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def spectogram_plot(self,silent=True,pdf=None): #TODO
        # result = self.np.abs(self.spectrogramValues[::1+len(self.spectrogramValues)//1000])
        result = self.np.abs(self.spectrogramValues[2500:5000])
        result[result == 0] = self.sys.float_info.epsilon
        fig = self.px.imshow(20*self.np.log10(result),
                             x= self.OST.get_freqAxis(2), #TODO
                             aspect='auto',
                             labels=dict(x='Frequency',y='Packets'),
                             title='Spectogram',
                            #  binary_string =True, 
                             binary_format="jpeg",
                             zmax = self.np.max(20*self.np.log10(result)),
                             zmin =  self.np.max(20*self.np.log10(result)) - 42,
                             binary_compression_level=2)

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/spectogram.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def columnMeanPower_plot(self,silent=True,pdf=None): 
        fig = self.px.line(self.powerColumns,labels={'index':'Samples','value':'Power [dBm]'},title='Mean Power for each column')
        fig.add_hline(y=(self.powerColumns).mean())
        fig.update_layout(showlegend = True)

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/meanPowerColumns.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def snr_plot(self,silent=True,pdf=None): 
        fig = self.px.line(self.snr,labels={'index':'packet','value':'SNR [dB]'},title='Signal to Noise Ratio')
        fig.add_hline(y=self.np.mean(self.snr))
        fig.update_traces(showlegend=True)

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/SNR.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def altitude_plot(self,silent=True,pdf=None): 
        fig = self.px.line(self.altitude,labels={'index':'packet','value':'altitude [m]'},title='Altitude')
        
        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/Altitude.html",auto_open=False)
        
        if pdf is not None:
            pdf.addPlot(fig)


    def saturatedSamples(self, ostRow):
        mask =self.get_mask(ostRow)
        counts,bins = self.np.histogram(self.np.concatenate(self.rawData[mask]).real,bins=self.np.arange(-128.5, 128.5, 1.0))
        counts2,bins = self.np.histogram(self.np.concatenate(self.rawData[mask]).imag,bins=bins)
        counts += counts2

        return 100 * (counts[0] + counts[-1]) / self.np.sum(counts)

    def bitHistogram_plot(self,silent=True,pdf=None): 
        counts,bins = self.np.histogram(self.np.concatenate(self.rawData).real,bins=self.np.arange(-128.5, 128.5, 1.0))
        counts2,bins = self.np.histogram(self.np.concatenate(self.rawData).imag,bins=bins)
        counts += counts2

        saturated = 100 * (counts[0] + counts[-1]) / self.np.sum(counts)
        bins = 0.5 * (bins[:-1] + bins[1:])
        fig = self.px.bar(x=bins,y=counts,labels={'x':'DN','y':'Counts'},title=f'RAW data histogram ({saturated:.2f}% saturated samples)')

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/HistogramAbsValues.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def mean_spectrum_plot(self,silent=True,pdf=None):  #OBSOLETE
        # freq = self.n p.flip(9e6+self.np.fft.fftshift(self.B+self.np.fft.fftfreq(len(self.spectrogramValues[0]),1/self.samplingFr)))
        fig = self.px.line( x= self.freqAxis, #TODO
                            y=20*self.np.log10(self.np.mean(self.np.abs(self.spectrogramValues), axis=0)),
                            title='Amplitude average spectrum',
                            labels={'x':'Frequency','y':'.'})

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/MeanSpectrum.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def sum_by_group(self,values, groups):
        order = self.np.argsort(groups)
        groups = groups[order]
        values = values[order]
        values.cumsum(out=values)
        index = self.np.ones(len(groups), 'bool')
        index[:-1] = groups[1:] != groups[:-1]
        values = values[index]
        groups = groups[index]
        values[1:] = values[1:] - values[:-1]
        return values, groups

    def mean_spectrum_plot_by_row(self,silent=True,pdf=None,csv=None):
        eps = self.np.finfo(self.calibratedSpectogramValues[0].dtype).eps
        sps = self.pd.DataFrame()
        for i,ostrow in self.OST.OST.iterrows():
            if ostrow.RIME_OST_Wait_Flag == 'Science_Data':
                # print(ostrow.RIME_OST_Entry_Id)000
                dfreq = self.OST.get_dfreq(ostrow.RIME_OST_Entry_Id)
                rxPathAtt = self.np.mean( self.attenuation[self.get_mask(ostrow.RIME_OST_Entry_Id)]) - self.RxGain

                # linSpectrum = self.np.mean(self.np.abs(self.np.array(self.calibratedSpectogramValues)[self.get_mask(ostrow.RIME_OST_Entry_Id)] /dfreq), axis=0)
                # sp = rxPathAtt + 10*self.np.log10(linSpectrum,out=self.np.zeros_like(linSpectrum), where=(linSpectrum!=0))
                spectrMatrix = self.np.abs( self.np.array(self.calibratedSpectogramValues)[self.get_mask(ostrow.RIME_OST_Entry_Id)]) / dfreq
                linSpectrum = self.np.mean(spectrMatrix, axis=0)
                sp = rxPathAtt +30+ 10*self.np.log10(linSpectrum,out=self.np.zeros_like(linSpectrum), where=(linSpectrum!=0))

                minSpectrum = self.np.min(spectrMatrix, axis=0)
                spmin = rxPathAtt +30+ 10*self.np.log10(minSpectrum,out=self.np.zeros_like(linSpectrum), where=(linSpectrum!=0))

                maxSpectrum = self.np.max(spectrMatrix, axis=0)
                spmax = rxPathAtt +30+ 10*self.np.log10(maxSpectrum,out=self.np.zeros_like(linSpectrum), where=(linSpectrum!=0))

                stdSpectrum = self.np.std(spectrMatrix, axis=0)
                spstd = 10*self.np.log10(stdSpectrum,out=self.np.zeros_like(linSpectrum), where=(linSpectrum!=0))
                
                spdf = self.pd.DataFrame()
                spdf['Spectrum'] = sp[:len(self.OST.get_freqAxis(ostrow.RIME_OST_Entry_Id))]
                spdf['SpectrumMin'] = spmin[:len(self.OST.get_freqAxis(ostrow.RIME_OST_Entry_Id))]
                spdf['SpectrumMax'] = spmax[:len(self.OST.get_freqAxis(ostrow.RIME_OST_Entry_Id))]
                spdf['SpectrumStd'] = spstd[:len(self.OST.get_freqAxis(ostrow.RIME_OST_Entry_Id))]
                spdf['Freqency'] = self.OST.get_freqAxis(ostrow.RIME_OST_Entry_Id)
                spdf['EntryId'] = ostrow.RIME_OST_Entry_Id
                spdf['BitNumber'] = ostrow.RIME_OST_Bit_Number
                spdf['RxAtt'] = ostrow.RIME_OST_Rx_Att
                spdf['Band'] = ostrow.RIME_OST_Band
                spdf['dFreq'] = dfreq
                sps = self.pd.concat([sps, spdf])
                
                from PIL import Image, ImageColor
                from plotly.express.colors import sample_colorscale
                spectrMatrix[spectrMatrix <=0] = eps
                specMatDB = 10*self.np.log10(self.np.fliplr(spectrMatrix[:,:len(self.OST.get_freqAxis(ostrow.RIME_OST_Entry_Id))]))
                vmax = self.np.percentile(specMatDB, 99) # -50
                vmin = self.np.percentile(specMatDB, 1) # -100 #
                specMatDB[specMatDB > vmax] = vmax
                specMatDB[specMatDB < vmin] = vmin
                exportSectrMat = ( 255* (specMatDB-vmin)/(vmax-vmin) ).astype(self.np.uint8)
                im = Image.fromarray(exportSectrMat)
                x = self.np.linspace(0, 1, 256)
                im.putpalette((self.np.array([ ImageColor.getcolor(c, "RGB") for c in sample_colorscale(self.px.colors.sequential.Magma, list(x))])).astype(self.np.uint8))
                im.save(self.outputFolder+f"/Spectrogram_{vmin}-{vmax}dB_row{i}.png")

                #TODO compute histogram of each column and save to numpy file and colormappe png, for each rowù
                # same method could be helpful for saving spectrograms
                # 1v stack values and verify that everithing still works
                # 2v apply colormap to spectrogam and save to png (maybe indexed png is better)
                # 3- compute histrogram of each column and save to numpy and colormapped png
                #   -35 -140dBW/Hz
                # move everything in generation function, not in the plot ones, as fields
                #  add export functions to save spectrograms and spetral histogram, as well to plot them
                bins = self.np.arange(-140,-35, 0.1)
                specHist  = self.np.zeros((len(bins)-1, len(self.OST.get_freqAxis(ostrow.RIME_OST_Entry_Id)) ))
                for s in range(specHist.shape[1]):
                    specHist[:,s] = self.np.histogram(10*self.np.log10(1e-100+spectrMatrix[:,s]), bins)[0]
                fig = self.px.imshow(specHist, y=bins[:-1], x=self.OST.get_freqAxis(ostrow.RIME_OST_Entry_Id), aspect='Auto',origin='lower')
                self.pl.offline.plot(fig,filename=self.outputFolder+f"/SpectrumHistogram_row{i}.html",auto_open=False)

        fig = self.px.line( sps,
                            x= 'Freqency',
                            y= 'Spectrum',
                            labels={'Freqency':'Frequency [Hz]','Spectrum':'Spectrum [dBW/Hz]'},
                            color = 'EntryId',
                            color_discrete_sequence= self.px.colors.sequential.Cividis,
                            title='Average Power Spectrum at antenna port (for each OST row)')
        # fig.update_yaxes(range = [-170,-120])
        if csv is not None:
            sps.to_csv(self.outputFolder+"/MeanSpectrum_byRow.csv")

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/MeanSpectrum_byRow.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def mean_spectrum_plot_by_row_attenuation(self, att =0, bit=None, yRange = None, customTitle = None, silent=True,pdf=None, png=None, csv=None):
        '''
        Power spectrum in dbmW, not normalized by DeltaFreq to evaluate complaince to EMI requirements
        '''
        sps = self.pd.DataFrame()
        
        if bit is None:
            OSTFilt = self.OST.OST.query("RIME_OST_Rx_Att == @att and RIME_OST_Wait_Flag == 'Science_Data'")
        else:
            OSTFilt = self.OST.OST.query("RIME_OST_Rx_Att == @att and RIME_OST_Wait_Flag == 'Science_Data' and RIME_OST_Bit_Number == @bit")

        if OSTFilt.empty:
            self.logger.warning(self.exID+'mean_spectrum_plot_by_row_attenuation: no rows with the selected paramters')
            return
        
        for _,ostrow in OSTFilt.iterrows():
            rxPathAtt = self.np.mean( self.attenuation[self.get_mask(ostrow.RIME_OST_Entry_Id)]) - self.RxGain
            # sp = rxPathAtt +30+ 10*self.np.log10(self.np.mean(self.np.abs( self.np.array(self.calibratedSpectogramValues)[self.get_mask(ostrow.RIME_OST_Entry_Id)]), axis=0))
            spectrMatrix = self.np.abs( self.np.array(self.calibratedSpectogramValues)[self.get_mask(ostrow.RIME_OST_Entry_Id)])
            linSpectrum = self.np.mean(spectrMatrix, axis=0)
            sp = rxPathAtt +30+ 10*self.np.log10(linSpectrum,out=self.np.zeros_like(linSpectrum), where=(linSpectrum!=0))

            minSpectrum = self.np.min(spectrMatrix, axis=0)
            spmin = rxPathAtt +30+ 10*self.np.log10(minSpectrum,out=self.np.zeros_like(linSpectrum), where=(linSpectrum!=0))

            maxSpectrum = self.np.max(spectrMatrix, axis=0)
            spmax = rxPathAtt +30+ 10*self.np.log10(maxSpectrum,out=self.np.zeros_like(linSpectrum), where=(linSpectrum!=0))

            spdf = self.pd.DataFrame()
            spdf['Spectrum'] = sp[:len(self.OST.get_freqAxis(ostrow.RIME_OST_Entry_Id))]
            spdf['SpectrumMin'] = spmin[:len(self.OST.get_freqAxis(ostrow.RIME_OST_Entry_Id))]
            spdf['SpectrumMax'] = spmax[:len(self.OST.get_freqAxis(ostrow.RIME_OST_Entry_Id))]
            spdf['Freqency'] = self.OST.get_freqAxis(ostrow.RIME_OST_Entry_Id)
            spdf['EntryId'] = ostrow.RIME_OST_Entry_Id
            sps = self.pd.concat([sps, spdf])

        fig = self.px.line( sps,
                            x= 'Freqency',
                            y= 'SpectrumMin',
                            labels={'Freqency':'Frequency [Hz]','Spectrum':'Spectrum [dBm]'},
                            color = 'EntryId',
                            color_discrete_sequence= self.px.colors.qualitative.Safe,
                            title=f'Average Power Spectrum at antenna port (OST rows with {att}dB attenuation) ' + customTitle)

        # fig.add_hline(y=-89,
        #             line_color="gray", 
        #             annotation_text="Limit 1MHz", 
        #             annotation_position="top right",
        #             annotation_font_color="black" )
        # fig.add_hline(y=-84,
        #             line_color="red", 
        #             annotation_text="Limit 3MHz", 
        #             annotation_position="top right",
        #             annotation_font_color="red" )

        if csv is not None:
            sps.to_csv(self.outputFolder+f"/MeanSpectrum_byRow_{att}_{bit}b.csv")

        if yRange is not None:
            fig.update_yaxes(range = yRange)

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+f"/MeanSpectrum_byRow_att{att}dB.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

        if png is not None:
            fig.write_image(self.outputFolder+f"/MeanSpectrum_byRow_att{att}dB.png",format="png", width=1920, height=1080, scale=1)

    def single_spectrum_plot(self,idx,silent=True,pdf=None): 
        fig = self.px.line(20*self.np.log10(self.np.abs(self.spectrogramValues[idx])),title='Amplitude spectrum')

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/SingleSpectrumPacket"+str(idx)+".html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def single_raw(self,idx,silent=True,pdf=None): 
        fig= self.px.line(self.np.real(self.rawData[idx]),labels={'index':'Samples','value':'Absolute Value'},title='Raw data: index '+str(idx))
        fig.update_traces(line_color='orange')
    
        fig2 = self.px.line(self.np.imag(self.rawData[idx]))
    
        fig3 = self.px.line(self.np.abs(self.rawData[idx]))
        fig3.update_traces(line_color='purple')
    
        fig.add_traces(list(fig2.select_traces()))
        fig.add_traces(list(fig3.select_traces()))

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/RawData["+str(idx)+"].html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def single_comp(self,idx,silent=True,pdf=None): 
        fig = self.px.line(20*self.np.log10(self.np.abs(self.compressed[idx])),layout_yaxis_range=[0,90])

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/CompressedData["+str(idx)+"].html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def rangeCompression_plot(self,silent=True,pdf=None): 
        #RANGE COMPRESSION
        image = self.np.abs(self.compressed[::1+len(self.compressed)//5000])
        fig = self.px.imshow(image,aspect='auto',labels=dict(x='Packet',y='Sample'),title='Range compression')

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/RangeCompression.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def rangeCompression2_plot(self,silent=True,pdf=None): 
        #RANGE COMPRESSION
        result = self.np.abs(self.compressed[::1+len(self.compressed)//5000])
        result[result == 0] = self.sys.float_info.epsilon
        image = 20*self.np.log10(result)
        fig = self.px.imshow(image,aspect='auto',labels=dict(x='Packet',y='Sample'),title='Range compression')

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/RangeCompressiondB.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def rangeCompressionShifted_plot(self,silent=True,pdf=None): 
        #RANGE COMPRESSION
        image = self.np.abs(self.shifted[::1+len(self.compressed)//5000])
        fig = self.px.imshow(image,aspect='auto',labels=dict(x='Packet',y='Sample'),title='Range compression')

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/RangeCompressionShiftedData.html",auto_open=False)

        if pdf is not None:
            pdf.addPlot(fig)

    def rangeRaw_plot(self,silent=True,pdf=None): 
        #RANGE COMPRESSION
        image = self.np.real(self.np.stack(self.rawData[::1+len(self.compressed)//5000]))
        fig = self.px.imshow(image,zmin=-128,zmax=128,aspect='auto',labels=dict(x='Packet',y='Sample'),title='Range compression')

        if not silent: 
            fig.show()
        else: 
            self.pl.offline.plot(fig,filename=self.outputFolder+"/RangeRawData.html",auto_open=False)
        if pdf is not None:
            pdf.addPlot(fig)

    def table_report(self,pdf,csv=False): 
        pdf.addBlankPage()

        pdf.addText('Statistics',100,150,fontSize=100,alignment='center')
        #Parameters
        

        OSTFilt = self.OST.OST.query("RIME_OST_Wait_Flag == 'Science_Data'")     
        stats = [['Row ID', "RxAtt [dB]", "Band", "SWL [us]", "BitNum", "Pre summing","Mean DN RMS", "Ant. Mean Power [dBm]", "Saturation [%]"]]
        for _,ostrow in OSTFilt.iterrows():
            mask =self.get_mask(ostrow.RIME_OST_Entry_Id)
            rxAtt = self.np.mean( self.attenuation[mask])

            # mean_powerADC = 30+ 10*self.np.log10( self.np.mean(self.power_ADC[mask]) )
            dnRMS = self.np.mean( self.RMS[mask])
            mean_powerAntenna = 30+ 10*self.np.log10(self.np.mean(  self.power_Antenna[mask]))
            sat = self.saturatedSamples(ostrow.RIME_OST_Entry_Id)
            stats.append([ostrow.RIME_OST_Entry_Id,
                        str(rxAtt),
                        ostrow.RIME_OST_Band,
                        str(ostrow.RIME_OST_SWL),
                        str(ostrow.RIME_OST_Bit_Number if ostrow.RIME_OST_Bit_Number>0 else 8),
                        str(2**(ostrow.RIME_OST_PresNum)),
                        f"{dnRMS:.2f}",
                        f"{mean_powerAntenna:.2f}",
                        f"{sat:.2f}"])
        
        pdf.addTable(stats,100,1850)

        if csv:
            df = self.pd.DataFrame(stats[1:], columns =stats[0])
            df.to_csv(self.outputFolder+"/stats_byRow.csv")

        if self.mode: #active
            max_alt = self.altitude.max()
            min_alt = self.altitude.min()
            max_pow_col = self.np.nanmax(self.powerColumns)
            max_snr = self.np.max(self.snr)
            min_snr = self.np.min(self.snr)

            stats = [['', "Value"],
            # ["Mean Power at ADC", f"{mean_powerADC:.2f} dBm"],
            ["Mean Power at Antenna ", f"{mean_powerAntenna:.2f} dBm"],
            ["MAX alt", str(max_alt)],
            ["MIN alt", str(min_alt)],
            ["MAX col pow", str(max_pow_col) + ' dBm'],
            ["MAX SNR", str(max_snr)+' dB'],
            ["MIN SNR", str(min_snr)+' dB']]
            
            pdf.addTable(stats,100,850)

        
        

