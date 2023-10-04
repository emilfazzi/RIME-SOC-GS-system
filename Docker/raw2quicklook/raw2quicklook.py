import os
import argparse
import logging
import pathlib
import time
import logging.config
import pickle

import multiprocessing as mp
import pandas as pd

from RIMERaw import RIMERaw
import quicklook
import rimeOST
from pdfgen import pdfGen

ostFiles ={
    0:"ost/FFT_00/OST3.xlsx", 
    5:"ost/RIM_01.2_RxOnly#5/OST5.xlsx",
    6:"ost/RIM_01.3_RxOnly#6/OST6.xlsx",
    7:"ost/RIM_02.5_TxOps#7/OST7.xlsx",
    8:"ost/RIM_02.6_RxPresumming#8/OST8.xlsx",
    14:"ost/Hmix#14/OST14.xlsx",
    15:"ost/Lmix#15/OST15.xlsx",
    16:"ost/Hatt#16/OST16.xlsx",
    17:"ost/Latt#17/OST17.xlsx",
    18:"ost/RIM_IIFC_RX#18/OST18.xlsx",
    19:"ost/RIM_IIFC_TX#19/OST19.xlsx"
}

def main():
   
    my_parser = argparse.ArgumentParser(description='Telemetry to raw', epilog='Example: python script.py -i inputFolder -o outputFolder -l outputFolderLogFiles -vid productVersionID')

    my_parser.add_argument('-i',
                       '--input',
                       action='store',
                       required=True,
                       help='Input folder containing the data-pack')

    # my_parser.add_argument('-s',
    #                     '--ost', 
    #                     action='store', 
    #                     required=True, 
    #                     help='OST file')

    my_parser.add_argument('-o',
                        '--output',
                        action='store',
                        required=True,
                        help='Output folder')

    my_parser.add_argument('-l',
                        '--log',
                        action='store',
                        required=True,
                        help='Output folder for the log files')

    my_parser.add_argument('-vid',
                        '--version-id',
                        action='store',
                        default='PDS4',
                        help='Product version identifier (not mandatory)')
 
    args = (my_parser.parse_args())

    inPath = pathlib.Path(args.input)
    outPath = pathlib.Path(args.output)
    logPath = pathlib.Path(args.log)

    outPath.mkdir(parents=True, exist_ok=True)
    logPath.mkdir(parents=True, exist_ok=True)

    #check input and log paths
    if not inPath.is_dir(): 
        loo.error(executableID+"Input folder does not exist")
        exit()

    #log file name definition
    executableName = pathlib.Path(__file__).stem #"tm2raw"
    instrumentName = 'RIME'
    executableID = str(int(time.time())) + '|'
    logFileName = logPath  / (instrumentName+'-'+executableName+'.log')

    #logger generation
    logging.config.fileConfig('logger.conf',defaults={'logfilename' : logFileName})
    loo = logging.getLogger('general')
    
    num_processors = mp.cpu_count()

   #extract all files in the input folder 
    try:
        files = list(inPath.glob('*.pkl'))
    except: 
        loo.error(executableID+'Zero files have been found in the input folder')
        exit()

    #INFO
    loo.info(executableID+'Start execution')
    loo.info(executableID+'Input folder: '+ inPath.as_posix())
    loo.info(executableID+'Output folder: '+ outPath.as_posix())
    loo.info(executableID+'Output folder for log files: '+ logPath.as_posix())
    for fileName in files: 
        loo.info(executableID+'Input file: '+ fileName.as_posix())

    #Read files
    for fileName in files: 

        loo.debug(executableID+'Reading file: ' + fileName.as_posix())
        rime = pickle.load( open( fileName, "rb" ) )

        ostId = rime.DataFrame.OST_Id.iloc[0]
        ost = rimeOST.rimeOST(pd.read_excel( ostFiles[ostId], usecols='B:U'), ostId=ostId)

        loo.info(executableID+'End of file reading')

        #QUICKLOOK IMPLEMENTATION 
        logging.config.fileConfig('logger_quicklook.conf',defaults={'logfilename' : logFileName})
        logger2 = logging.getLogger('quicklook')

        logger2.info(executableID+'Quicklook start')
        q = quicklook.quicklook(rime=rime, mode='passive',processors=num_processors,exec=executableID,ost=ost)

        # Output folder definition
        logger2.debug(executableID+'Quicklook export')
        q.file_input = fileName

        q.setOutputFolder((outPath / fileName.name).as_posix())
        os.makedirs(q.outputFolder,exist_ok=True)
        # pickle.dump( q, open( file_name+"_ql.pkl", "wb" ) )

        silent = True

        p = pdfGen(fileName,q.outputFolder,imagesPerPage=3,columns=1)
        #Plot generation
        q.single_raw(1,silent,p)
        q.bitHistogram_plot(silent,p)
        q.rowRMSAmplitude_plot(silent,p)
        q.rowMeanPowerADC_plot(silent,p)
        q.rowMeanPower_plot(silent,p)
        q.mean_spectrum_plot_by_row(silent,p, csv=True)
        q.bitSaturation(silent,p)
        # q.spectogram_plot(silent,p)
        # q.mean_spectrum_plot_by_row_attenuation(att = 15, customTitle = 'Antenna fully deployed', silent=silent, pdf=p)

        if q.mode:  # only for active mode
            q.columnMeanPower_plot(silent,p)
            q.altitude_plot(silent,p)
            q.snr_plot(silent,p)
            q.rangeCompression2_plot(silent,p)    
        q.table_report(p, csv=True)
        p.pdfFileGenerator()

    
    logger2.info(executableID+'Quicklook end')
    

    logging.shutdown()


if __name__ == '__main__': 
    main()