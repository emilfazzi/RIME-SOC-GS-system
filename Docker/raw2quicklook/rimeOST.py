class rimeOST:
    """
    Parsing of RIME Telemetries
    JUI-TASR-RIM-ICD-002 issue 09
    copy of Alberto's version 18/03/2023 
    from commit:
        0cdf1d5
        RIME flyby programmer, optimisation Tested on trajectory "crema_5_1_150lb_23_1" version "v422"
        Branch,develop
        2023-02-20
    Then edited by Max
    """
    
    import re
    import os
    import pandas as pd
    import numpy as np
    from jinja2 import Environment, FileSystemLoader

    re_OST_id = re.compile(r'LEP10026\t(\d*)')
    re_OST_TotalEntry = re.compile(r'LEP10090\t(\d*)')
    re_OST_TCEntry = re.compile(r'LEP10091\t(\d*)')

    sel_sampling_freq = {
        'High_Band': 3e6,
        'Low_Band': 1.2e6
    }

    sel_gain_selector = {
        'Nom_Attenuation' : 0,
        'Max_Attenuation' : 19
    }

    paramList = ['LEP10092','LEP10093','LEP10094','LEP10096','LEP10095','LEP10027','LEP10028','LEP10029','LEP10030','LEP10036','LEP10031','LEP10032','LEP10034','LEP10037','LEP10098','LEP10038','LEP10039','LEP10040','LEP10033','LEP10097']
    
    first_entry = {"RIME_OST_Entry_Id": 1,
                "RIME_OST_PRI_D": 1,
                "RIME_OST_PRI_I": 1,
                "RIME_OST_Time_Margin": 20,
                "RIME_OST_PresNum": 0,
                "RIME_OST_Blk_Len": 7264,
                "RIME_OST_Wait_Flag": "No_Science_Data",
                "RIME_OST_Gain_Selector": "Nom_Attenuation",
                "RIME_OST_PRI": 4000,
                "RIME_OST_Band": "High_Band",
                "RIME_OST_Traker_Mode": 0,
                "RIME_OST_SWST": 500,
                "RIME_OST_Frequency_Shift": 0,
                "RIME_OST_Rank": 0,
                "RIME_OST_Rx_Att": 0,
                "RIME_OST_SWST_Rate": 0,
                "RIME_OST_Bit_Number": 0,
                "RIME_OST_SWL": 300,
                "RIME_OST_Pulse_Selection": "100_Microseconds",
                "RIME_OST_SWL_Tck": 300}

    last_entry = {"RIME_OST_Entry_Id": 0,
                "RIME_OST_PRI_D": 1,
                "RIME_OST_PRI_I": 1,
                "RIME_OST_Time_Margin": 20,
                "RIME_OST_PresNum": 0,
                "RIME_OST_Blk_Len": 2000,
                "RIME_OST_Wait_Flag": "No_Science_Data",
                "RIME_OST_Gain_Selector": "Nom_Attenuation",
                "RIME_OST_PRI": 4000,
                "RIME_OST_Band": "High_Band",
                "RIME_OST_Traker_Mode": 0,
                "RIME_OST_SWST": 500,
                "RIME_OST_Frequency_Shift": 0,
                "RIME_OST_Rank": 0,
                "RIME_OST_Rx_Att": 0,
                "RIME_OST_SWST_Rate": 0,
                "RIME_OST_Bit_Number": 0,
                "RIME_OST_SWL": 300,
                "RIME_OST_Pulse_Selection": "100_Microseconds",
                "RIME_OST_SWL_Tck": 300}

    def __init__(self, OST, ostId=None, TC=False):
        if isinstance(OST, str):
            if TC:
                # Read telecomand from csv
                self.OST = self.parseOSTTC(self.pd.read_csv(self.os.path.join(OST), sep=";"))
                self.cast_type()
            else:
                # Read OST from csv
                assert ostId!=None, "Check OST constructor, 'ostId' must be specified!"
                self.ostId = ostId
                self.OST = self.parseOST(OST)
                # self.OST = self.OST.drop([0,self.OST.shape[0]-1], axis=0)
                self.cast_type()
                # self.OST["RIME_OST_Entry_Id"] = self.OST["RIME_OST_Entry_Id"]-1
        else:
            if TC:
                # Read telecomand from dataframe
                self.OST = self.parseOSTTC(OST)
                self.cast_type()
            else:
                # Read OST from dataframe
                #assert ostId!=None,"Check OST constructor, 'ostId' must be specified!"
                self.ostId = ostId if ostId is not None else 0
                self.OST = OST
                self.cast_type()
        assert self.OST.shape[0]<=24, "Too many OST rows!"
        assert self.OST.shape[1]==20, "Incorrect column names or number!"
        
    def cast_type(self):
        self.OST = self.OST.astype({"RIME_OST_Entry_Id": 'int64',
                                    "RIME_OST_PRI_D": 'int64',
                                    "RIME_OST_PRI_I": 'int64',
                                    "RIME_OST_Time_Margin": 'int64',
                                    "RIME_OST_PresNum": 'int64',
                                    "RIME_OST_Blk_Len": 'int64',
                                    "RIME_OST_Wait_Flag": 'string',
                                    "RIME_OST_Gain_Selector": 'string',
                                    "RIME_OST_PRI": 'int64',
                                    "RIME_OST_Band": 'string',
                                    "RIME_OST_Traker_Mode": 'int64',
                                    "RIME_OST_SWST": 'int64',
                                    "RIME_OST_Frequency_Shift": 'int64',
                                    "RIME_OST_Rank": 'int64',
                                    "RIME_OST_Rx_Att": 'int64',
                                    "RIME_OST_SWST_Rate": 'int64',
                                    "RIME_OST_Bit_Number": 'int64',
                                    "RIME_OST_SWL": 'int64',
                                    "RIME_OST_Pulse_Selection": 'string',
                                    "RIME_OST_SWL_Tck": 'int64'})
    
    def getParam(self, code, line):
        mo = self.re.search(code+r'\t(.*)\t(.*)\n', line)
        return (mo.group(2).strip(), mo.group(1).strip())

    def checkOST(self):
        OSTContraints = {"RIME_OST_PRI_I":          [0, 7],
                         "RIME_OST_Time_Margin":    [0, 63],
                         "RIME_OST_PresNum":        [0, 7],
                         "RIME_OST_Blk_Len":        [1, 1048575],
                        #  "RIME_OST_Gain_Selector":  [0, 1],
                         "RIME_OST_PRI":            [1000, 10000],
                        #  "RIME_OST_Band":           [0, 1],
                         "RIME_OST_Traker_Mode":    [0, 1],
                         "RIME_OST_SWST":           [0, 10000],
                         "RIME_OST_Frequency_Shift":[-4096, 4095],
                         "RIME_OST_Rank":           [0, 7],
                         "RIME_OST_Rx_Att":         [0, 31],
                         "RIME_OST_SWST_Rate":      [-64, 63],
                         "RIME_OST_Bit_Number":     [0, 7],
                         "RIME_OST_SWL":            [0, 2560],
                        #  "RIME_OST_Pulse_Selection":[0, 7],
                         "RIME_OST_SWL_Tck":        [0, 2047]}
        for k in OSTContraints:
            assert all(self.OST[k] >= OSTContraints[k][0]), "Invalid Parameter in OST: lower than expected " + k
            assert all(self.OST[k] <= OSTContraints[k][1]), "Invalid Parameter in OST: higher than expected " + k
        assert all(self.OST['RIME_OST_SWL']+self.OST["RIME_OST_SWST"] <= self.OST["RIME_OST_PRI"]), "Invalid Parameter in OST: SWST must be lower than PRI"
        assert all( (self.OST['RIME_OST_Blk_Len'] % (2**self.OST['RIME_OST_PRI_I'] * 2**self.OST['RIME_OST_PRI_D']))==0 ), "Invalid Parameter in OST: Block_Len must be integer multiple of  2^PRI_I_TCK * 2^PRI_D_TCK"
        assert all( (self.OST['RIME_OST_Blk_Len'] % (2**self.OST['RIME_OST_PresNum'] ))==0 ), "Invalid Parameter in OST: Block_Len must be integer multiple of  2^PRI_I_TCK * 2^PRI_D_TCK"

        assert all( self.OST['RIME_OST_SWST']+ \
                    self.OST['RIME_OST_Blk_Len']*self.OST['RIME_OST_PRI']*1e-6*self.OST['RIME_OST_SWST_Rate']> 315), \
                        'SWST<0 at block end' 
        
        assert all( self.OST['RIME_OST_SWL']+ \
                    self.OST['RIME_OST_SWST']+ \
                    self.OST['RIME_OST_Blk_Len']*self.OST['RIME_OST_PRI']*1e-6*self.OST['RIME_OST_SWST_Rate'] < self.OST['RIME_OST_PRI']), \
                        'SWST+SWL>PRI at block end' 
        return True

    def parseOST(self, path):
        """
        Parse the OST in a given path selected by ostNumber
        OST files are like OST1_1.csv
        Result is returned as pandas dataframe
        """

        totalEntries = 0
        df  = self.pd.DataFrame()
        for OSTsubN in range(1,4):
            with open(self.os.path.join(path+'OST'+str(self.ostId)+'_'+str(OSTsubN)+'.csv')) as f:
                contents = f.read()
            fileSplit = contents.split('#')
            line = fileSplit[3]

            RIME_OST_Entry_Id = int(self.re_OST_id.search(line).group(1))
            RIME_OST_TotalEntry = int(self.re_OST_TotalEntry.search(line).group(1) )
            RIME_OST_TCEntry = int(self.re_OST_TCEntry.search(line).group(1) )

            rows = []
            for i in range(4, len(fileSplit)-1):
                line = fileSplit[i]
                dictRow = {}
                for lep in self.paramList:
                    x = self.getParam(lep, line)
                    dictRow[x[0]] = x[1]
                rows.append(dictRow)
            df = self.pd.concat([df, self.pd.DataFrame(rows)])

            totalEntries += RIME_OST_TCEntry
            if totalEntries >= RIME_OST_TotalEntry:
                break
        df = df.reset_index(drop=True)
        return df

    def wordSlice(self, msg, line, start=0, lenght=16):
        """lenght bits starting at start from the N-th Word to int"""
        if ((start+lenght)>16):
            raise Exception("Unvalid paramters: start+lenght >16")
        word = msg[(line-1)*4:(line-1+1)*4]
        # shift to the right then mask the unnecessary left part 
        print(word)
        return ( int(str(word),16) >> (16-(lenght+start)) ) & ((2**lenght)-1)

    def twosComplement(self, val, bits):
        # bits = 16 # Number of bits in a hexadecimal number format
        # val = int(hexval, bits)
        if val & (1 << (bits-1)):
            val -= 1 << bits
        return val

    def parseOSTTC(self, inputrow):
        '''
        TC (207,132) - Load Operation Sequence Table (OST) (LEC10024)
        '''
        pkt = inputrow.PRIMARY_HEADER+inputrow.DATAFIELD_HEADER+inputrow.DATAFIELD
        print(pkt)
        OST = self.pd.DataFrame()
        self.ostId = self.wordSlice( pkt, 6)
        OST['Id'] = [ self.ostId  ] 
        OST['TotalEntries'] = [ self.wordSlice( pkt, 7, 0, 8) ]
        OST['TC_Defined_Entries'] = [ self.wordSlice( pkt, 7, 8, 8) ] 
        OST['Table'] = [ self.pd.DataFrame() ]

        for N in range(OST['TC_Defined_Entries'].values[0]):
            row = {}
            row['RIME_OST_Entry_Id'] = self.wordSlice( pkt, 8+9*N, 2, 6)
            row['RIME_OST_PRI_D'] = self.wordSlice( pkt, 8+9*N, 9, 3)
            row['RIME_OST_PRI_I'] = self.wordSlice( pkt, 8+9*N, 13, 3)
            row['RIME_OST_Time_Margin'] = self.wordSlice( pkt, 9+9*N, 2, 6)
            row['RIME_OST_PresNum'] = self.wordSlice( pkt, 9+9*N, 9, 3)
            RIME_OST_Blk_LenMS = self.wordSlice( pkt, 9+9*N, 12, 4)
            RIME_OST_Blk_LenLS = self.wordSlice( pkt, 10+9*N)
            row['RIME_OST_Blk_Len'] = RIME_OST_Blk_LenMS*(2**16) + RIME_OST_Blk_LenLS
            row['RIME_OST_Wait_Flag'] = self.wordSlice( pkt, 11+9*N, 0, 1)
            row['RIME_OST_Gain_Selector'] = self.wordSlice( pkt, 11+9*N, 1, 1)
            row['RIME_OST_PRI'] = self.wordSlice( pkt, 11+9*N, 2, 14)
            row['RIME_OST_Band'] = self.wordSlice( pkt, 12+9*N, 0, 1)
            row['RIME_OST_Traker_Mode'] = self.wordSlice( pkt, 12+9*N, 1, 1)
            row['RIME_OST_SWST'] = self.wordSlice( pkt, 12+9*N, 2, 14)
            row['RIME_OST_Frequency_Shift'] = self.twosComplement( self.wordSlice( pkt, 13+9*N, 3, 13), 13)
            row['RIME_OST_Rank'] = self.wordSlice( pkt, 14+9*N, 0, 3)
            row['RIME_OST_Rx_Att'] = self.wordSlice( pkt, 14+9*N, 3, 5)
            row['RIME_OST_SWST_Rate'] = self.twosComplement( self.wordSlice( pkt, 14+9*N, 9, 7), 7)
            row['RIME_OST_Bit_Number'] = self.wordSlice( pkt, 15+9*N, 1, 3)
            row['RIME_OST_SWL'] = self.wordSlice( pkt, 15+9*N, 4, 12)
            row['RIME_OST_Pulse_Selection'] = self.wordSlice( pkt, 16+9*N, 1, 3)
            row['RIME_OST_SWL_Tck'] = self.wordSlice( pkt, 16+9*N, 4, 12)
            rowDf = self.pd.DataFrame(row, index=[N])
            OST['Table'] = self.pd.concat([OST['Table'][0], rowDf])
        
        return OST['Table']

    def mergeOST(self, df):
        dfout = self.pd.DataFrame()
        for i in range(len(df)):
            pass
        
    def exportOST(self, out_path=None, addFirstLast=False):
        """
        Export the OST in CSV format to be used on RIME EGSE
        WARNING: the line TC_TAG		{{TC_TAG}} has 2 tabs in the middle
        """

        if addFirstLast:
            N_tot_entries = self.OST.shape[0]+2
            self.last_entry["RIME_OST_Entry_Id"] = N_tot_entries
            
            # self.OST["RIME_OST_Entry_Id"] = self.OST["RIME_OST_Entry_Id"]+1
            self.OST = self.pd.concat( \
                [self.pd.DataFrame([self.first_entry]),\
                self.OST,\
                self.pd.DataFrame([self.last_entry])])
        else:            
            N_tot_entries = self.OST.shape[0]
                                    
        N_files = (self.OST.shape[0])//11
        #TODO check, probably not needed, add empty file at the end
        if (self.OST.shape[0])%11 != 0:
            N_files +=1

        if out_path==None:
            # Export string variables
            data = []
            result_filenames = []
            for file in range(N_files):
                entries = self.OST.iloc[11*file:11*(file+1)]
                result_filenames.append("OST"+str(self.ostId)+"_"+str(file+1)+".csv")
                environment = self.Environment(loader=self.FileSystemLoader(self.os.path.join(self.os.path.dirname(__file__), "templates/")))
                template = environment.get_template("OST_template.txt")
                context = {
                    "TC_TAG" : "LEC10024",
                    "RIME_OST_Id" : str(self.ostId),
                    "RIME_OST_Total_Entry" : N_tot_entries,
                    "RIME_OST_TC_Entry" : entries.shape[0],
                    "entries" : entries}
                data.append(template.render(context))
            return [data, result_filenames]
        else:
            # Export and save to csv
            for file in range(N_files):
                entries = self.OST.iloc[11*file:11*(file+1)]
                results_filename = self.os.path.join(out_path,"OST"+str(self.ostId)+"_"+str(file+1)+".csv")
                environment = self.Environment(loader=self.FileSystemLoader(self.os.path.join(self.os.path.dirname(__file__), "templates/")))
                template = environment.get_template("OST_template.txt")
                context = {
                    "TC_TAG" : "LEC10024",
                    "RIME_OST_Id" : str(self.ostId),
                    "RIME_OST_Total_Entry" : N_tot_entries,
                    "RIME_OST_TC_Entry" : entries.shape[0],
                    "entries" : entries}
                with open(results_filename, mode="w", encoding="utf-8") as results:
                    results.write(template.render(context))

    def exportOST2SPOT(self, out_path=None):
        """
        Export the OST in CSV format to be ingested by ESA SPOT tool
        WARNING: 1st and last line should be already in
        """
        rowsPerBlock = 3
                    
        N_tot_entries = self.OST.shape[0]
        self.last_entry["RIME_OST_Entry_Id"] = N_tot_entries
              
        N_files = (self.OST.shape[0])//rowsPerBlock

        environment = self.Environment(loader=self.FileSystemLoader(self.os.path.join(self.os.path.dirname(__file__), "templates/")))
        template = environment.get_template("OST_SPOT_template.txt")

        data = []
        for file in range(N_files):
            entries = self.OST.iloc[rowsPerBlock*file:rowsPerBlock*(file+1)]
            context = {
                "TC_TAG" : "LEC10024",
                "RIME_OST_Id" : str(self.ostId),
                "RIME_OST_Total_Entry" : N_tot_entries,
                "RIME_OST_TC_Entry" : entries.shape[0],
                "entries" : entries}
            data.append(template.render(context))

        results_filename = self.os.path.join(out_path,"RIM_PDOR_OST"+str(self.ostId)+".csv")
        if out_path==None:
            # Export string variables
            return ["\n".join(data), results_filename]
        else:
            # Export and save to csv
            with open(results_filename, mode="w", encoding="utf-8") as results:
                results.write("\n".join(data))
                # for row in data:
                #     results.write(row+"\n")
        # print(data)
        # return data

    def estimateResources(self, out_path=None):

        mask = self.np.double([1 if x else 0 for x in self.OST.RIME_OST_Wait_Flag == 'Science_Data'])
        undersampling = self.np.double([4 if x else 10 for x in self.OST.RIME_OST_Band == 'High_Band'])
        presum = 2** self.np.double(self.OST.RIME_OST_PresNum)
        SWL = 1e-6 * self.np.double(self.OST.RIME_OST_SWL)
        bitNum = self.np.double([8 if x==0 else x for x in self.OST.RIME_OST_Bit_Number] )
        blkLen = self.np.double(self.OST.RIME_OST_Blk_Len)
        PRI = 1e-6 * self.np.double(self.OST.RIME_OST_PRI)

        SW_bitSize = mask *( (16*16)+ (12e6 * SWL * 2 * bitNum ) / (undersampling * presum) )
        blkDatavolume = SW_bitSize * blkLen
        totalDatavolume = self.np.sum(blkDatavolume)
        blkDatarate = SW_bitSize / PRI
        totalOSTtime = self.np.sum(blkLen * PRI)

        msg =  '- - - - - - - - - - - - - - - - - - \n'
        msg += f'Total OST time    [sec] = {totalOSTtime:.2f} \n'
        msg += f'Total Datavolume [Gbit] = {totalDatavolume*1e-9:.3f} \n'
        msg += f'Max Datarate     [Mbps] = {max(blkDatarate*1e-6):.2f} \n'
        msg += f'Mean Datarate    [Mbps] = {1e-6*totalDatavolume/totalOSTtime:.2f} \n'
        msg += f'OST Row, datarate and Datavolume [id s Mbps Mbit DC%] = \n'    
        for k,v1,v2,v3, v4 in zip(range(len(blkDatarate)), blkLen * PRI, blkDatarate, self.np.array(blkDatavolume), 100*100e-6/(PRI)):
            msg += f'{k: 4d}{v1: 10.2f}{v2/1e6: 10.2f}{v3/1e6: 10.1f}{v4: 10.1f} \n'

        print(msg)
        if out_path is not None:
            with open(self.os.path.join(out_path,"resources.txt"), "w") as text_file:
                text_file.write(msg)


    def get_attenuation(self, start=0, stop=None):
        extraAttenuation = self.OST['RIME_OST_Gain_Selector'].map({'Nom_Attenuation':0,'Max_Attenuation':19})
        attenuationValue = (self.OST['RIME_OST_Rx_Att'] + extraAttenuation)[self.OST['RIME_OST_Wait_Flag'] == 'Science_Data']
        numberOccurences = (self.OST['RIME_OST_Blk_Len']//2**self.OST['RIME_OST_PresNum'])[self.OST['RIME_OST_Wait_Flag'] == 'Science_Data']
        out = self.np.repeat(attenuationValue.tolist(),numberOccurences.tolist())
        return out[start:] if stop is None else out[start:stop+1]
            
    
    def get_rowId(self):
        attenuationValue = (self.OST['RIME_OST_Entry_Id'])[self.OST['RIME_OST_Wait_Flag'] == 'Science_Data']
        numberOccurences = (self.OST['RIME_OST_Blk_Len']//2**self.OST['RIME_OST_PresNum'])[self.OST['RIME_OST_Wait_Flag'] == 'Science_Data']
        return self.np.repeat(attenuationValue.tolist(),numberOccurences.tolist())
    
    def get_mask(self, id):
        return self.get_rowId() == id
        
    def get_PRI(self): 
        numberOccurences = (self.OST['RIME_OST_Blk_Len']//2**self.OST['RIME_OST_PresNum'])[self.OST['RIME_OST_Wait_Flag'] == 'Science_Data']
        PRI = self.OST['RIME_OST_PRI'][self.OST['RIME_OST_Wait_Flag'] == 'Science_Data']
        return self.np.repeat(PRI.tolist(),numberOccurences.tolist())
    
    def get_freqAxis(self, row):
        swl = 1e-6*self.OST.iloc[row-1]["RIME_OST_SWL"]
        timeStep = 1/self.sel_sampling_freq[self.OST.iloc[row-1]["RIME_OST_Band"]]
        numSamples =  round(swl * self.sel_sampling_freq[self.OST.iloc[row-1]["RIME_OST_Band"]])
        return self.np.flip(9e6+(1/swl)+self.np.fft.fftshift(self.np.fft.fftfreq(numSamples,timeStep)))
    
    def get_dfreq(self, row):
        freqAx = self.get_freqAxis(row)
        return abs(freqAx[0] - freqAx[1])

    def plot_SWST(self):
        import plotly.graph_objects as go
        import plotly as pl
        from plotly.subplots import make_subplots
        import bisect

        OST = self.OST

        blkLen = self.np.double(self.OST.RIME_OST_Blk_Len)
        PRI = 1e-6 * self.np.double(self.OST.RIME_OST_PRI)
        blkTime = PRI*blkLen
        block_edges = self.np.concatenate((0, self.np.cumsum(blkTime)), axis=None)

        slowTime = self.np.arange(start=0, stop=block_edges[-1], step=1.0 )
        # SWST_profile = data.Profile.values[0]["SWST"]
        SWST_linear_profile = self.np.zeros_like(slowTime)
        SWL = self.np.zeros_like(slowTime)
        PulseLenght = self.np.zeros_like(slowTime)
        PRI = self.np.zeros_like(slowTime)
        rxAtt = self.np.zeros_like(slowTime)
        bits = self.np.zeros_like(slowTime)

        for i in range(len(slowTime)):
            if i<=0:
                OST_row_index = 0
            elif i>=block_edges[-1]:
                OST_row_index = OST.shape[0]-1
            else:
                OST_row_index = bisect.bisect_right(block_edges, i)-1
            SWST_linear_profile[i] = round(OST.iloc[OST_row_index]["RIME_OST_SWST"] + OST.iloc[OST_row_index]["RIME_OST_SWST_Rate"]*(i-block_edges[OST_row_index]))
            if self.OST.iloc[OST_row_index]['RIME_OST_Wait_Flag'] == 'Science_Data':
                SWL[i] = OST.iloc[OST_row_index]["RIME_OST_SWL"]
                zz = OST.iloc[OST_row_index]["RIME_OST_Bit_Number"]
                bits[i] = 8 if zz == 0 else zz
            else:
                SWL[i] = 0
                bits[i] = 0
            PulseLenght[i] = int(OST.iloc[OST_row_index]["RIME_OST_Pulse_Selection"].split("_")[0])
            PRI[i] = OST.iloc[OST_row_index]["RIME_OST_PRI"]
            rxAtt[i] = OST.iloc[OST_row_index]["RIME_OST_Rx_Att"] + self.sel_gain_selector[OST.iloc[OST_row_index]["RIME_OST_Gain_Selector"]]
            
        
        # max_SWST = data.Profile.values[0]["PRI_profile"]*1e6 - 2 - SWL
        
        DC = 100*PulseLenght / PRI
        
        fig1 = go.Figure()
        fig1 = pl.subplots.make_subplots(rows=2, cols=2)

        fig1.add_trace(go.Scatter(x=slowTime, y=SWST_linear_profile, mode='lines', line= {"shape": 'hv'}, name='SWST', line_width=2, line_color="red"), row =1, col=1)
        fig1.add_trace(go.Scatter(x=slowTime, y=SWST_linear_profile+SWL, mode='lines', line= {"shape": 'hv'}, name='Rx Window', line_width=0, line_color="red", fill='tonexty'), row =1, col=1)
        fig1.add_trace(go.Scatter(x=slowTime, y=PulseLenght, mode='lines', line= {"shape": 'hv'}, name='Tx Pulse', line_width=2, line_color="Blue", fill='tozeroy'), row =1, col=1)
        fig1.add_trace(go.Scatter(x=slowTime, y=PRI, mode='lines', line= {"shape": 'hv'}, name='PRI', line_width=2, line_color="Orange"), row =1, col=1)

        fig1.update_yaxes(range=[0, 10000], row =1, col=1)
        
        fig1.add_trace(go.Scatter(x=slowTime, y=DC, line= {"shape": 'hv'}, mode='lines', name='Ducty Cycle [%]', line_width=2, line_color="red"), row =2, col=1)
      
        fig1.update_yaxes(range=[0, 12], row =2, col=1)

        fig1.add_trace(go.Scatter(x=slowTime, y=rxAtt, line= {"shape": 'hv'}, mode='lines', name='Rx Att [dB]', line_width=2, line_color="purple"), row =1, col=2)
      
        fig1.update_yaxes(range=[0, 50], row =1, col=2)
       
        fig1.add_trace(go.Scatter(x=slowTime, y=bits, line= {"shape": 'hv'}, mode='lines', name='Bit Number', line_width=2, line_color="green"), row =2, col=2)
      
        fig1.update_yaxes(range=[0, 10], row =2, col=2)
        
        fig1['layout']['xaxis']['title']='Slow time [s]'
        fig1['layout']['xaxis2']['title']='Slow time [s]'
        fig1['layout']['xaxis3']['title']='Slow time [s]'
        fig1['layout']['xaxis4']['title']='Slow time [s]'
        fig1['layout']['yaxis']['title']='Fast time [us]'
        fig1['layout']['yaxis2']['title']='Attenuation [dB]'
        fig1['layout']['yaxis3']['title']='DC [%]'
        fig1['layout']['yaxis4']['title']='Bit Number'
        # fig1.show()

        # fig1 = go.Figure()
        # fig1 = pl.subplots.make_subplots(rows=1, cols=1)

        # fig1.add_trace(go.Scatter(x=slowTime, y=SWST_linear_profile, mode='lines', line= {"shape": 'hv'}, name='SWST', line_width=2, line_color="red"), row =1, col=1)
        # fig1.add_trace(go.Scatter(x=slowTime, y=SWST_linear_profile+SWL, mode='lines', line= {"shape": 'hv'}, name='Rx Window', line_width=0, line_color="red", fill='tonexty'), row =1, col=1)
        # fig1.add_trace(go.Scatter(x=slowTime, y=PulseLenght, mode='lines', line= {"shape": 'hv'}, name='Tx Pulse', line_width=2, line_color="Blue", fill='tozeroy'), row =1, col=1)
        # fig1.add_trace(go.Scatter(x=slowTime, y=PRI, mode='lines', line= {"shape": 'hv'}, name='PRI', line_width=2, line_color="Orange"), row =1, col=1)

        # fig1.update_yaxes(range=[0, 10000], row =1, col=1)
        # fig1['layout']['xaxis']['title']='Slow time [s]'
        # fig1['layout']['yaxis']['title']='Fast time [us]'


        # SWST_error = (SWST_linear_profile-SWST_profile*1e6)
        # SWST_error[SWST_error<0] = 0
        # SWST_margin = 299792458*SWST_error*1e-6
        
        # fig2 = make_subplots(rows=1, cols=1, specs=[[{"secondary_y": True}]])
        # fig2.add_trace(go.Scatter(x=time_to_ca, y=SWST_error, mode='lines', name='Time error'), row=1, col=1, secondary_y=False)
        # fig2.add_trace(go.Scatter(x=time_to_ca, y=SWST_margin/1e3, mode='lines', name='Margin', line_color="gray"), row=1, col=1, secondary_y=True)

        # fig2.update_yaxes(title_text="Margin [km]", secondary_y=True)
        # fig2.update_layout(
        #     height=300,
        #     yaxis=dict(color="#636EFA"),
        #     xaxis_title="Time to Closest Approach [s]",
        #     yaxis_title="SWST error [us]",
        #     margin=dict(t=25),)
        
        return fig1