import os
import platform
import sys
import time
import threading

from FastqUtils import FastqHandler
from controlutils import HelpTheMinion
import configargparse
from watchdog.observers.polling import PollingObserver as Observer

from concurrent.futures import ThreadPoolExecutor
from requests_futures.sessions import FuturesSession

session = FuturesSession(executor=ThreadPoolExecutor(max_workers=10))



# noinspection PyGlobalUndefined
if __name__ == '__main__':

    global OPER

    OPER = platform.system()
    if OPER is 'Windows':  # MS
        OPER = 'windows'
    else:
        OPER = 'linux'  # MS
    #print(OPER)  # MS

    if OPER is 'linux':
        config_file = os.path.join(os.path.sep, \
                                   os.path.dirname(os.path.realpath('__file__' \
                                                                    )), 'minfq_posix.config')
    if OPER is 'windows':
        config_file = os.path.join(os.path.sep, sys.prefix,
                                   'minfq_windows.config')

    parser = \
        configargparse.ArgParser(description='minFQ: A program to analyse minION fastq files in real-time or post-run and monitor the activity of MinKNOW.' \
                                 , default_config_files=[config_file])
    parser.add(
        '-w',
        '--watch-dir',
        type=str,
        required=True,
        default=None,
        help='The path to the folder containing the downloads directory with fast5 reads to analyse - e.g. C:\\data\\minion\\downloads (for windows).',
        dest='watchdir'
    )

    parser.add(
        '-ip',
        '--ip-address',
        type=str,
        dest='ip',
        required=True,
        default=None,
        help='The IP address of the minKNOW machine.',
    )

    parser.add(
        '-v',
        '--verbose',
        action='store_true',
        help='Display debugging information.',
        default=False,
        dest='verbose',
    )
    parser.add(
        '-nf',
        '--no_fastq',
        action='store_true',
        help='Run minFQ without monitoring fastq files.',
        default=False,
        dest='noFastQ',
    )

    parser.add(
        '-nm',
        '--no_minKNOW',
        action='store_true',
        help='Run minFQ without monitoring minKNOW for live activity.',
        default=False,
        dest='noMinKNOW',
    )

    parser.add(
        '-n',
        '--name',
        type=str,
        required=True,
        default=None,
        help='The run name you wish to provide.',
        dest='run_name',
    )

    parser.add(
        '-f',
        '--is_flowcell',
        action='store_true',
        help='If you add this flag, all runs added here will be considered as a single flow cell with the name set by the name flag.',
        dest='is_flowcell',
    )

    parser.add(
        '-k',
        '--key',
        type=str,
        required=True,
        default=None,
        help='The api key for uploading data.',
        dest='api_key',
    )

    parser.add(
        '-p',
        '--port',
        type=int,
        #required=True,
        default=8100,
        help='The port number for the local server.',
        dest='port_number',
    )

    parser.add(
        '-hn',
        '--hostname',
        type=str,
        #required=True,
        default='127.0.0.1',
        help='The run name you wish to provide.',
        dest='host_name',
    )

    parser.add(
        '-c',
        '--custom_barcode',
        type=str,
        required=False,
        default=None,
        help='Optionally split reads based on odd/even channel description. Not a standard option.',
        dest='cust_barc'
    )

    parser.add(
        '-s',
        '--skip_sequence',
        action='store_true',
        required=False,
        help='If selected only read metrics, not sequence, will be uploaded to the databse.',
        dest='skip_sequence'
    )

    args = parser.parse_args()
    # Makes no sense to run if both no minKNOW and no FastQ is set:
    if args.noFastQ and args.noMinKNOW:
        print ("You must monitor either FastQ or MinKNOW.")
        print ("This program will now exit.")
        os._exit(0)


    args.full_host = "http://" + args.host_name + ":" + str(args.port_number) + "/"

    if args.is_flowcell:
        print ("Flowcell option is set.")
    #GLobal creation of header (needs fixing)
    global header
    header = {'Authorization': 'Token ' + args.api_key, 'Content-Type': 'application/json'}


    if not args.noFastQ:
        event_handler = FastqHandler(args,header)
        # This block handles the fastq
        observer = Observer()
        observer.schedule(event_handler, path=args.watchdir, recursive=True)
        observer.daemon = True

    if not args.noMinKNOW:
        # this block is going to handle the running of minControl.
        minwsip = "ws://" + args.ip + ":9500/"
        helper = HelpTheMinion(minwsip, args)
        helper.connect()

    try:
        if not args.noFastQ:
            observer.start()
        if not args.noMinKNOW:
            t = threading.Thread(target=helper.process_minion())
            t.daemon = True
            t.start()
        while 1:
            time.sleep(1)
    except KeyboardInterrupt:
        print(": ctrl-c event")
        if not args.noMinKNOW:
            helper.mcrunning=False
            helper.hang_up()
        if not args.noFastQ:
            observer.stop()
            observer.join()
        os._exit(0)

