import os 
import pathlib
import argparse 
import time
import struct
import logging, logging.config
import pickle

import multiprocessing as mp 
import numpy as np

from RIMERaw import RIMERaw

def main():
    
    my_parser = argparse.ArgumentParser(description='Telemetry to raw', epilog='Example: python script.py -i inputFolder -o outputFolder -l outputFolderLogFiles -vid productVersionID')

    my_parser.add_argument('-i',
                       '--input',
                       action='store',
                       required=True,
                       help='Input folder containing the data-pack')

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

    processors = mp.cpu_count()

    #extract all files in the input folder 
    try:
        files = list(inPath.glob('RIM1_38*')) + list(inPath.glob('*SCIENCE_OST*.txt'))
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
    #TODO how labels should be generated for multiple files?
    for fileName in files: 

        loo.debug(executableID+'Reading file: ' + fileName.as_posix())
        rawData = RIMERaw(fileName.as_posix(),processors=processors,exec=executableID)

        pickle.dump( rawData, open( fileName+"_raw.pkl", "wb" ) )
        continue
        #Time stamp correction 
        loo.debug(executableID+'Start time stamp coreection')
        # rawData.applyTimeStampCorrection()
        loo.debug(executableID+'End time stamp correction')

        #Write output file
        loo.debug(executableID+'Writing output file')

        #pad all rangeline to max lenght and stack in 2d array
        rawVal = rawData.DataFrame['data'].values
        maxLen = max([len(x) for x in rawData.DataFrame['data'].values])
        for i, x in enumerate(rawVal):
            rawVal[i] = np.pad(x, (0, maxLen-x.shape[0]), mode='constant')
        rawVal = np.stack(rawVal)

        info = {}
        info['file_I'] = outPath /  (fileName.stem + '_I.dat')
        info['file_Q'] = outPath / (fileName.stem +'_Q.dat')
        info['file_Extra'] = outPath / (fileName.stem +'_extra.csv')

        with open(info['file_I'], 'wb') as f:
            intVal = np.real(rawVal).astype('int8')
            np.save(f, intVal)

        with open(info['file_Q'], 'wb') as f:
            intVal =np.imag(rawVal).astype('int8')
            np.save(f, intVal)

        rawData.DataFrame.drop(columns=['data']).to_csv(info['file_Extra'])

        with open(info['file_I'], 'rb') as f:
            # compute <offset unit="byte"> 128 </offset> always multiple of 64
            # valid for v1.0, header len as  unsigned short little endian
            f.seek(8)
            info['headerSize'] = 9+ struct.unpack('<H', f.read(2))[0] +1

        info['mission_phase'] = 'NECP'

        # write PDS Label
        loo.debug(executableID+'Writing PDS Label')
        rawData.pdslabel_generation(outPath / fileName.stem, info)

        loo.debug(executableID+'Output file available in the output folder')

    loo.info(executableID+'End of execution')
    logging.shutdown()


if __name__ == '__main__': 
    main()
