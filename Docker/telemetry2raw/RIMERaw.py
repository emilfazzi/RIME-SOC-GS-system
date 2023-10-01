from joblib import Parallel, delayed
import xml.etree.ElementTree as ET

class RIMERaw:
    import os
    from io import BytesIO
    from bitstring import BitArray
    # import multiprocessing as mp 

    import numpy as np   
    import pandas as pd
    import logging
    import datetime as dt
    """
    Parsing of RIME Telemetries
    JUI-TASR-RIM-ICD-002 issue 09
    """

    parsers = {'3' : {'25':'HKReport'},
             '207' : {'167':'Science'}}

    
    OBTSYNCMASK = (1<<31)-1 # '0x7fffffff'
    # FINETIMECONST = 1.01105
    FINETIMECONST = 1.01

    PEC_TEST = False

    def __init__(self,file='',processors=1,exec="0|"): ## costructor 
        self.exID = exec

        self.Table = self.initTable()
        self.logger = self.logging.getLogger("RIME")

        file = str(file)

        if len(file)>0:
            if not self.PEC_TEST:
                self.logger.debug(self.exID+'PEC field check disabled')
            if (file.endswith('.txt')):
                self.logger.debug(self.exID+'Reading text file')
                self.DataFrame = self.subcostructorTxt(file,processors)
            else: 
                self.logger.debug(self.exID+'Reading binary file')
                self.DataFrame = self.subcostructorBin(file,processors)

    def initTable(self): 
        table = self.np.zeros(256,dtype=self.np.uintc)
        for i in range(0,256,1):
            tmp = 0
            if((i & 1) != 0): 
                tmp = tmp ^ 0x1021
            if ((i & 2) != 0): 
                tmp = tmp ^ 0x2042
            if ((i & 4) != 0): 
                tmp = tmp ^ 0x4084
            if ((i & 8) != 0): 
                tmp = tmp ^ 0x8108
            if ((i & 16) != 0): 
                tmp = tmp ^ 0x1231
            if ((i & 32) != 0): 
                tmp = tmp ^ 0x2462
            if ((i & 64) != 0): 
                tmp = tmp ^ 0x48C4
            if ((i & 128) != 0): 
                tmp = tmp ^ 0x9188
            table[i] = tmp
        return table 

    def subcostructorBin(self,file,processors): 

        #binary file
        tlmList = []
        try:
            with open(file,'rb') as f: 
                buf = self.BytesIO(f.read())
                if self.os.stat(file).st_size == 0:

                    self.logger.error(self.exID+'File ' + file +  ' is empty')
                    raise Exception('The current file ' +  file +' is empty')

                self.logger.debug(self.exID+'Binary file reading start')

                while pkt_header := buf.read(6):
                    header = self.HeaderBin(pkt_header)
                    pkt_length = header['PacketLength'] #pkt length is the 3rd word 

                    if pkt_length == 0: 
                        self.logger.error(self.exID+'Packet dimension is incrorect')
                        raise ValueError('Packet length is 0, packet dimension is incorrect')

                    pkt_data = buf.read(pkt_length+1) #add PEC
                    msg = pkt_header + pkt_data

                    if self.PEC_TEST: 
                        if not self.CrcOptBin(msg[:-2],msg[-2:]):
                            self.logger.error(self.exID+'PEC field is not valid')
                            raise ValueError('Error: PEC field not valid')

                    tlmList.append((header,pkt_data)) #TEST read txt
                    # tlmList.append((header,msg))
                
                self.logger.debug(self.exID+'Binary file reading end')
                
        except: 
            self.logger.fatal(self.exID+'Error opening the file')
            raise RuntimeError('File opening has failed')

        df = self.pd.DataFrame(tlmList)
        self.logger.info(self.exID+str(processors) + ' processors available for multiprocessing')
        seq = df.values.tolist()
        # result = [self.ScienceBin(seq[jj]) for jj in range(len(seq))]

        try: 
            self.logger.debug(self.exID+'Dataframe processing start')
            result = Parallel(n_jobs=processors)(delayed(self.ScienceBin)(seq[jj]) for jj in range(len(seq)))
            # with self.closing(self.mp.Pool(processes=processors)) as pool: 
            #     result = pool.map(self.ScienceBin,seq)
        except: 
            self.logger.fatal(self.exID+'System has not been able to implement multiproccessing')
            raise SystemError('Multiprocessing implementation has failed')

        self.logger.debug(self.exID+'Dataframe processing end')
        return self.pd.DataFrame(result)
    
    def subcostructorTxt(self,file,processors): 
        #text file
        try: 
            if self.os.stat(file).st_size == 0:
                    self.logger.error(self.exID+'File ' + file +  ' is empty')
                    raise Exception('The current file ' +  file +' is empty')
            
            self.logger.debug(self.exID+'Text file reading start')
            df = self.pd.read_table(file,header=None,dtype=str)
            self.logger.debug(self.exID+'Text file reading end')
        except: 
            self.logger.error(self.exID+'Error reading the file')
            raise RuntimeError('File reading has failed')

        self.logger.info(self.exID+str(processors) + ' processors available for multiprocessing')
        seq = df.values.tolist()
        try: 
            self.logger.debug(self.exID+'Dataframe processing start')
            result = Parallel(n_jobs=processors)(delayed(self.ScienceTxt)(seq[jj]) for jj in range(len(seq)))

            # with self.closing(self.mp.Pool(processes=processors)) as pool: 
            #     result = pool.map(self.ScienceTxt,seq)
        except: 
            self.logger.fatal(self.exID+'System has not been able to implement multiproccessing')
            raise SystemError('Multiprocessing implementation has failed')
        
        self.logger.debug(self.exID+'Dataframe processing end')
        return self.pd.DataFrame(result)
        

    def dataInt(self,data,nb): 
        #MSB
        return (data<<(8-nb)) | ((int(not(data>>(nb-1)))*(pow(2,8-nb)-1)))

    ##FOR BINARY FILES 
    def CrcOptBin(self,DataIN,SyndromeIN):
        Syndrome = 0xFFFF
        for j in range(0,len(DataIN)): 
            Data = DataIN[j]
            Syndrome =  (((Syndrome <<8) & 0xFF00) ^ (self.Table[(((Syndrome >> 8) ^ Data) & 0x00FF)]))
            Data = Data << 1
        return Syndrome == int(self.hexGeneration(SyndromeIN[0],SyndromeIN[1]))

    def CrcControlBin(self, DataIN, SyndromeIN):
        '''PEC control function, return TRUE if syndrome is correct'''
        # translated from C example here: ECSS-E-70-41A 30 January 2003 page 210
        # could be improved with the optimized version reported in the same document
        # unsigned char Data; /* Byte to be encoded */
        # unsigned Syndrome; /* Original CRC syndrome */
        '''PEC control function, return TRUE if syndrome is correct'''
        Syndrome = 0xFFFF
        for j in range(0, len(DataIN)):
            Data = DataIN[j]
            for i in range(0,8):
                if ((Data & 0x80) ^ ((Syndrome & 0x8000) >> 8)):
                    Syndrome = ((Syndrome << 1) ^ 0x1021) & 0xFFFF
                else:
                    Syndrome = (Syndrome << 1) & 0xFFFF
                Data = Data << 1 
        return Syndrome == int(self.hexGeneration(SyndromeIN[0],SyndromeIN[1]))
    
    def wordSliceBin(self,n1,n2,start=0,dataLength=16): 
        word = self.hexGeneration(n1,n2)
        return (word >> (16-dataLength-start)) & ((2**dataLength)-1)
    
    def hexToBin(self,line,pos1,pos2):
        int_value = int.from_bytes(line[pos1:pos2+1],"big")
        return bin(int_value)[2:].zfill(16)

    def hexGeneration(self,n1,n2): 
        return n1*(2**8)+n2
   
    def CUCTimeStamp(self, coarse, fine):
        return float(coarse & self.OBTSYNCMASK) + self.FINETIMECONST*float(fine)*(2**-16)

    def HeaderBin(self,msgHeader): 
        out = {}
        out['VN'] = self.wordSliceBin(n1=msgHeader[0],n2=msgHeader[1],start=0,dataLength=3)
        out['T'] = self.wordSliceBin(n1=msgHeader[0],n2=msgHeader[1],start=3,dataLength=1)
        out['D'] = self.wordSliceBin(n1=msgHeader[0],n2=msgHeader[1],start=4,dataLength=1)
        out['ProcessId'] = self.wordSliceBin(n1=msgHeader[0],n2=msgHeader[1],start=5,dataLength=7)
        out['Category'] = self.wordSliceBin(n1=msgHeader[0],n2=msgHeader[1],start=12,dataLength=4)

        out['GF'] = self.wordSliceBin(n1=msgHeader[2],n2=msgHeader[3],start=0,dataLength=2)
        out['SequenceCounter'] = self.wordSliceBin(n1=msgHeader[2],n2=msgHeader[3],start=2,dataLength=14)
        out['PacketLength'] = self.wordSliceBin(n1=msgHeader[4],n2=msgHeader[5])
        return out

    def ScienceBin(self,row): 
        msg = row[1]
        out = row[0]
        
        #TEST
        # return self.ScienceTxt([row[1].hex()])

        #data field header
        out['P'] = self.wordSliceBin(n1=msg[0],n2=msg[1],start=3,dataLength=1)
        out['ServiceType'] = self.wordSliceBin(n1=msg[0],n2=msg[1],start=8,dataLength=8)
        out['ServiceSubType'] = self.wordSliceBin(n1=msg[2],n2=msg[3],start=0,dataLength=8)
        out['Destination'] = self.wordSliceBin(n1=msg[2],n2=msg[3],start=8,dataLength=8)
        
        #read time stamp 
        word1 = self.hexGeneration(n1=msg[4],n2=msg[5])
        word2 = self.hexGeneration(n1=msg[6],n2=msg[7])
        out['TimeStampCoarse'] = word1*(2**16)+word2
        out['TimeStampFine'] = self.hexGeneration(n1=msg[8],n2=msg[9])
        out['TimeStamp'] = self.CUCTimeStamp(out['TimeStampCoarse'],out['TimeStampFine']) #with mask
        
        out['priCounter'] = (self.wordSliceBin(n1=msg[10],n2=msg[11]))*(2**16)+(self.wordSliceBin(n1=msg[12],n2=msg[13]))
        out['SWST'] = self.wordSliceBin(n1=msg[14],n2=msg[15],start=1,dataLength=15) 
        out['Df'] = self.wordSliceBin(n1=msg[14],n2=msg[15],start=0,dataLength=1) 
        out['N_OST'] = self.wordSliceBin(n1=msg[16],n2=msg[17],start=0,dataLength=8) 
        out['Tx'] = self.wordSliceBin(n1=msg[16],n2=msg[17],start=8,dataLength=2)
        out['Rx'] = self.wordSliceBin(n1=msg[16],n2=msg[17],start=10,dataLength=2)
        out['SubMode_Id'] = self.wordSliceBin(n1=msg[16],n2=msg[17],start=12,dataLength=4)
        out['Pulse_sel'] = self.wordSliceBin(n1=msg[18],n2=msg[19],start=1,dataLength=3)
        out['Bit_Num'] = self.wordSliceBin(n1=msg[18],n2=msg[19],start=5,dataLength=3)
        if out['Bit_Num'] == 0: 
            out['Bit_Num'] = 8
        out['B'] = self.wordSliceBin(n1=msg[18],n2=msg[19],start=9,dataLength=1)
        out['M'] = self.wordSliceBin(n1=msg[18],n2=msg[19],start=11,dataLength=1)
        out['Rank'] = self.wordSliceBin(n1=msg[18],n2=msg[19],start=13,dataLength=3)
        out['OST_Id'] = self.wordSliceBin(n1=msg[20],n2=msg[21])
        out['N_Acquisition'] = (self.wordSliceBin(n1=msg[22],n2=msg[23]))*(2**16)+(self.wordSliceBin(n1=msg[24],n2=msg[25]))
        
        #check data length and data read 
        # dataLength = int((out['PacketLength']-1)/2-13)
        dataLength = (out['PacketLength']+1-(14*2))*8//(2*out['Bit_Num'])

        nb = out['Bit_Num'] 
        s = self.BitArray(msg[26:(len(msg)-2)])
        temp = self.np.int8(s.unpack(str(dataLength*2)+'*int'+str(out['Bit_Num']))) <<(8-nb)
        # same speed as bistring
        # temp = self.np.int8(self.bitstruct.unpack(dataLength*2*('s'+str(out['Bit_Num'])), msg[26:(len(msg)-2)])) <<(8-nb)

        out['data'] = self.np.complex64(temp[::2] + 1j* temp[1::2] )

        # if dataLength != int(len(msg[26:(len(msg)-2)])/2): 
        #     self.logger.warning(self.exID+'Packet length value from packet is not equal to the real packet length')
        #     # ValueError("Error in the pkt length ")

        # word_result = 0
        # for i in range(26,round(dataLength)+26,2):
        #     word = self.hexGeneration(msg[i],msg[i+1])
        #     if word_result == 0: 
        #         word_result = word
        #     else: 
        #         word_result = word_result*(2**16) + word

        # nb = out['Bit_Num']
        # j = range(1,dataLength*2,2)
        # data = self.np.zeros(round(dataLength),self.np.complex128)
        # try: 
        #     for i in range(0,dataLength): 
        #         dataI = self.np.int8(((word_result >> (dataLength*(2*nb)-j[i]*nb)) & (pow(2,nb)-1))<<(8-nb))
        #         dataQ = self.np.int8(((word_result >> (dataLength*(2*nb)-j[i]*nb-nb)) & (pow(2,nb)-1))<<(8-nb))
        #         data[i] = dataI+1j*dataQ
        # except: 
        #     self.logger.fatal(self.exID+'Error in reading I and Q data')
        #     raise RuntimeError('Error reading data')
        # out['data'] = data
        return self.pd.Series(out)


    ##FOR TEXT FILES
    def CrcOpt(self,DataIN,SyndromeIN):
        Syndrome = 0xFFFF
        for j in range(0,round(len(DataIN)/2)): 
            Data = int(DataIN[2*j:2*j+2],16)
            Syndrome =  (((Syndrome <<8) & 0xFF00) ^ (self.Table[(((Syndrome >> 8) ^ Data) & 0x00FF)]))
        return Syndrome == int(SyndromeIN,16)

    def CrcControlTxt(self,DataIN, SyndromeIN):
        '''PEC control function, return TRUE if syndrome is correct'''
        # translated from C example here: ECSS-E-70-41A 30 January 2003 page 210
        # could be improved with the optimized version reported in the same document
        # unsigned char Data; /* Byte to be encoded */
        # unsigned Syndrome; /* Original CRC syndrome */
        Syndrome = 0xFFFF
        for j in range(0, round(len(DataIN)/2)):
            Data = int(DataIN[2*j: 2*j+2],16)
            for i in range(0,8):
                if ((Data & 0x80) ^ ((Syndrome & 0x8000) >> 8)):
                    Syndrome = ((Syndrome << 1) ^ 0x1021) & 0xFFFF
                else:
                    Syndrome = (Syndrome << 1) & 0xFFFF
                Data = Data << 1  
        return Syndrome == int(SyndromeIN, 16)

    def hexSlice16( self,msg, start, length):
        """N-th to N+length Word to int"""
        return int(msg[(start-1)*4:(start-1+length)*4],16)

    def hexSlice8( self,msg, start, offset):
        """half of N-th Word to int, offset 0 for left 1 for right"""
        return int(msg[(start-1)*4+(2*offset):(start-1)*4+2+(2*offset)],16)

    def wordSliceTxt(self, msg, line, start=0, length=16):
        """lenght bits starting at start from the N-th Word to int"""
        if ((start+length)>16): 
            raise Exception("Invalid parameters: start+lenght >16")

        word = msg[(line-1)*4:(line)*4]
        # shift to the right then mask the unnecessary left part 
        return ( int(word,16) >> (16-(length+start)) ) & ((2**length)-1)
    

    def HeaderTxt(self,msg): 
        msg = msg.strip('\n')
        if self.PEC_TEST: 
            if not self.CrcOpt(msg[:-4],msg[-4:]): 
                self.logger.error(self.exID+'Invalid PEC field')
                raise ValueError('Error: PEC field not valid')
        out = {}
        out['VN'] =  self.wordSliceTxt(msg, line=1, start=0, length=3)
        out['T'] =  self.wordSliceTxt(msg, line=1, start=3, length=1)
        out['D'] =  self.wordSliceTxt(msg, line=1, start=4, length=1)
        out['ProcessId'] =  self.wordSliceTxt(msg, line=1, start=5, length=7)
        out['Category'] =  self.wordSliceTxt(msg, line=1, start=12, length=4)
        out['GF'] =  self.wordSliceTxt(msg, line=2, start=0, length=2)
        out['SequenceCounter'] =  self.wordSliceTxt(msg, line=2, start=2, length=14) 
        out['PacketLength'] =  self.wordSliceTxt(msg, line=3) 
       
        return out 

    def DataFieldHeaderTxt(self, msg, out):
        #data field header
        out['P'] =  self.wordSliceTxt(msg,line=4,start=3,length=1)
        out['ServiceType'] =  self.wordSliceTxt(msg, line=4, start=8, length=8)
        out['ServiceSubType'] =  self.wordSliceTxt(msg, line=5, start=0, length=8)
        out['Destination'] =  self.wordSliceTxt(msg, line=5, start=8, length=8)

        #read time stamp 
        out['TimeStampCoarse'] =  self.hexSlice16(msg,start=6,length=2)
        out['TimeStampFine'] =  self.wordSliceTxt(msg,line=8)
        out['TimeStamp'] =  self.CUCTimeStamp(out['TimeStampCoarse'],out['TimeStampFine'])
        out['TimeSync'] = (out['TimeStampCoarse'] & ~self.OBTSYNCMASK) == 0
        return out

    def ScienceTxt(self,row): 
        msg = row[0]
        out = self.HeaderTxt(msg)
        out = self.DataFieldHeaderTxt(msg, out)

        # #data field header
        # out['P'] =  self.wordSliceTxt(msg,line=4,start=3,length=1)
        # out['ServiceType'] =  self.wordSliceTxt(msg, line=4, start=8, length=8)
        # out['ServiceSubType'] =  self.wordSliceTxt(msg, line=5, start=0, length=8)
        # out['Destination'] =  self.wordSliceTxt(msg, line=5, start=8, length=8)
        
        # #read time stamp 
        # out['TimeStampCoarse'] =  self.hexSlice16(msg,start=6,length=2)
        # out['TimeStampFine'] =  self.wordSliceTxt(msg,line=8)
        # out['TimeStamp'] =  self.CUCTimeStampTxt(out['TimeStampCoarse'],out['TimeStampFine'])
        
        out['priCounter'] =  self.hexSlice16(msg,start=9,length=2)
        out['SWST'] =  self.hexSlice16(msg,start=11,length=1) & (1<<15)-1 # mask 1st bit with '0x7fff'
        out['Df'] =  self.hexSlice8(msg, start=11,offset=0)>>7 #== 1
        out['N_OST'] =  self.hexSlice8(msg, start=12,offset=0)
        out['Tx'] =  self.hexSlice8(msg, start=12,offset=1) >>6 
        out['Rx'] = ( self.hexSlice8(msg, start=12,offset=1) & int('0x3f',16)) >>4 
        out['SubMode_Id'] = ( self.hexSlice8(msg, start=12,offset=1) & int('0x0f',16))
        out['Pulse_sel'] = ( self.hexSlice8(msg, start=13,offset=0) & int('0x7f',16)) >>4 
        out['Bit_Num'] = ( self.hexSlice8(msg, start=13,offset=0) & int('0x07',16)) 
        if out['Bit_Num'] == 0: 
            out['Bit_Num'] = 8
        out['B'] = ( self.hexSlice8(msg,start=13,offset=1) & int('0x7f',16)) >>6
        out['M'] = ( self.hexSlice8(msg, start=13,offset=1) & int('0x1f',16)) >> 4 
        out['Rank'] = ( self.hexSlice8(msg, start=13,offset=1) & int('0x07',16)) 
        out['OST_Id'] =  self.hexSlice16(msg,start=14,length=1)
        out['N_Acquisition'] =  self.hexSlice16(msg,start=15,length=2)

        #DATA READ
        # dataLength = int((out['PacketLength']-1)/2-13)
        # dataLength = int( (len(msg)//4 ) -1 -16 )
        # fix TAS 15/2/23
        dataLength = (out['PacketLength']+1-(14*2))*8//(2*out['Bit_Num'])

        #check data length 
       
        nb = out['Bit_Num'] 
        data= self.np.zeros(dataLength,self.np.complex128)
        s = self.hexSlice16(msg,start=17,length=dataLength)
        j = range(1,dataLength*2,2)

        try: 
            for i in range(0,dataLength):
                # dataI = self.np.int8(self.dataInt(((s >> (dataLength*16-j[i]*nb)) & (pow(2,nb)-1)),nb))
                # dataQ = self.np.int8(self.dataInt(((s >> (dataLength*16-j[i]*nb-nb)) & (pow(2,nb)-1)),nb))
                dataI = self.np.int8(((s >> (dataLength*(2*nb)-j[i]*nb)) & (pow(2,nb)-1))<<(8-nb))
                dataQ = self.np.int8(((s >> (dataLength*(2*nb)-j[i]*nb-nb)) & (pow(2,nb)-1))<<(8-nb))
                data[i] = dataI+1j*dataQ
        except: 
            self.logger.fatal(self.exID+'Error in reading I and Q data')
            raise RuntimeError('Error reading data')

        out['data'] = data 
        return self.pd.Series(out)

    def CheckErrorAtTheBeginning(self,TSListCoarse):
        diff = self.np.diff(TSListCoarse)
        start = min(self.np.argmax(diff > 1), self.np.argmax(diff<0))
        sel = diff[diff != 0]
        if (sel.shape)[0] == 0: 
            return TSListCoarse
        check = self.np.argmax(sel!=1)
        end = self.np.argmax(diff[start:] == 1) + start
        if start == end: 
            self.logger.debug(self.exID+'Wrong time stamp up to the last packet')
            first_group = diff[start:]
        else: 
            first_group = diff[start:end]
        first_group = first_group[first_group != 0]
        if len(first_group) != self.np.sum(first_group): 
            if check != 0:
                return TSListCoarse
            if start == end: 
                return TSListCoarse
            self.logger.debug(self.exID+'Wrong time stamp at the beginning')
            ref = TSListCoarse[end]
            count = len(first_group)
            TSListCoarse[0:start+1] = ref-count 
        return TSListCoarse

    def TlcReport(self, msg): 
        out= self.HeaderTxt(msg) 
        out['TelecommandPktId'] = self.wordSliceTxt(msg, 9) 
        out['PacketSequenceControl'] = self.wordSliceTxt(msg, 10)         
        return out 

    def EventReport(self, msg):
        out= self.HeaderTxt(msg)
        out['RIDId1'] = self.wordSliceTxt(msg,9)
        out['RIDId2'] = self.wordSliceTxt(msg,10)

    def applyTimeStampCorrection(self):
        TSListCoarse = self.DataFrame['TimeStampCoarse']
        TSListFine = self.DataFrame['TimeStampFine']
        self.DataFrame['TimeStamp'] = self.TimeStampCorrection(TSListCoarse,TSListFine)

    def TimeStampCorrection(self,TSListCoarse,TSListFine): 
        TSListFine = TSListFine.to_numpy()
        TSListCoarse = TSListCoarse.to_numpy()
        TSListCoarse = self.CheckErrorAtTheBeginning(TSListCoarse)
        TSListCorrected = [] 
        flag = False #correct or not
        refCT = 0
        correctCT = 0
        refError = 0
        for coarseTime,fineTime in zip(TSListCoarse,TSListFine): 
            if refCT == 0: 
                refCT = coarseTime
                correctCT = coarseTime
            else: 
                deltaCT = coarseTime-refCT
                if not flag: 
                    if (0 <= deltaCT <= 1): 
                        correctCT = coarseTime
                        refCT = correctCT
                        refError = 0 
                        flag = False
                    elif deltaCT > 1 or deltaCT < 0:
                        self.logger.debug(self.exID+'First error in time stamp')
                        refError = coarseTime
                        correctCT = refCT+1
                        refCT = correctCT
                        flag = True
                else: 
                    if 0 <= deltaCT <= 1: 
                        refCT = coarseTime
                        correctCT = coarseTime
                        refError = 0
                        flag = False
                    elif deltaCT >= 2 or deltaCT < 0: 
                        if refError != coarseTime: 
                            self.logger.debug(self.exID+'New error in the time stamp')
                            refError = coarseTime
                            correctCT = refCT+1
                            refCT = correctCT
                            flag = True
                        else: 
                            self.logger.debug(self.exID+'Same error in the time stamp')
                            refError = coarseTime
                            correctCT = refCT
                            flag = True
            ts = self.CUCTimeStamp(correctCT,fineTime)
            TSListCorrected.append(ts)
        return self.pd.Series(TSListCorrected)


    def pdslabel_generation(self,xmlFileName, info): 
        tree = ET.parse('pds4label_template.xml')
        root = tree.getroot()

        ns = root.tag.split('}')[0]+'}'

        #Date update
        title = root.find(f'.//*{ns}title')
        date = self.dt.datetime.now()
        title.text = 'RIME raw data'

        modification = root.find(f'.//*/*{ns}modification_date')
        modification.text = date.strftime('%Y-%m-%d')
       
        #Time update
        start_date_time = root.find(f'.//*/*{ns}start_date_time')
        start_date_time.text = (self.dt.datetime(2000, 1, 1)+ self.dt.timedelta(seconds=self.DataFrame['TimeStamp'].iloc[0]) ).isoformat()[:-3]+'Z'
        stop_date_time = root.find(f'.//*/*{ns}stop_date_time')
        stop_date_time.text = (self.dt.datetime(2000, 1, 1)+ self.dt.timedelta(seconds=self.DataFrame['TimeStamp'].iloc[-1]) ).isoformat()[:-3]+'Z'

        spacecraft_clock_start = root.find('.//*/*{http://psa.esa.int/psa/v1}spacecraft_clock_start_count')
        spacecraft_clock_start.text = (self.dt.datetime(2000, 1, 1)+ self.dt.timedelta(seconds=self.DataFrame['TimeStamp'].iloc[0]) ).isoformat()[:-3]+'Z'
        spacecraft_clock_stop = root.find('.//*/*{http://psa.esa.int/psa/v1}spacecraft_clock_stop_count')
        spacecraft_clock_stop.text = (self.dt.datetime(2000, 1, 1)+ self.dt.timedelta(seconds=self.DataFrame['TimeStamp'].iloc[-1]) ).isoformat()[:-3]+'Z'

        #Mission phase update
        mission_phase_name = root.find('.//*/*{http://psa.esa.int/psa/v1}mission_phase_name')
        mission_phase_name.text = info['mission_phase']

        #Input file name update
        file_name_input = root.find('.//*/*{http://psa.esa.int/psa/v1}file_name')
        file_name_input.text = 'NOME FILE INPUT'

        #Ouptut file information update
        file_name_output = root.find(f'.//*/*{ns}file_name')
        file_name_output.text = 'NOME FILE OUTPUT'
        file_size = root.find(f'.//*/*{ns}file_size')
        file_size.text = 'DIMENSIONE IN BYTE'
        checksum = root.find(f'.//*/*{ns}md5_checksum')
        checksum.text = 'CHECKSUM'


        #### SAVE XML FILE 
        tree.write(str(xmlFileName)+'.xml')
