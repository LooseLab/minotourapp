#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import re
import time
import errno
from socket import error as socket_error
import urllib2
import threading
import json
import multiprocessing
import copy
import platform
import hashlib
import MySQLdb
import configargparse
from ws4py.client.threadedclient import WebSocketClient
from thrift import Thrift
from thrift.transport import TTransport
from thrift.protocol import TCompactProtocol
from dictdiffer import diff, patch, swap, revert
import numpy as np

from twisted.internet.protocol import ReconnectingClientFactory

from autobahn.twisted.websocket import WebSocketClientProtocol, \
    WebSocketClientFactory



lock = threading.Lock()

# import memcache



# Unbuffered IO
# sys.stdin = os.fdopen(sys.stdin.fileno(), 'w', 0) # MS
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)  # MS
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)  # MS

global OPER
OPER = platform.system()
if OPER == 'Windows':  # MS
    OPER = 'windows'
elif OPER == 'Darwin':
    OPER = 'osx'
else:
    OPER = 'linux'  # MS
print OPER  # MS


config_file = script_dir = os.path.dirname(os.path.realpath('__file__' \
        )) + '/' + 'minup_posix.config'
parser = \
    configargparse.ArgParser(description= \
    'interaction: A program to provide real time interaction for minION runs.' \
    , default_config_files=[config_file])
parser.add(
    '-ws',
    '--websocket-host',
    type=str,
    dest='wshost',
    required=True,
    default='localhost',
    help="The location of the minotour server being connected to.",
)
parser.add(
    '-dbh',
    '--mysql-host',
    type=str,
    dest='dbhost',
    required=False,
    default='localhost',
    help="The location of the MySQL database. default is 'localhost'.",
    )
parser.add(
    '-dbu',
    '--mysql-username',
    type=str,
    dest='dbusername',
    required=True,
    default=None,
    help='The MySQL username with create & write privileges on MinoTour.'
    ,
    )
parser.add(
    '-dbp',
    '--mysql-port',
    type=int,
    dest='dbport',
    required=False,
    default=3306,
    help="The MySQL port number, else the default port '3306' is used."
    ,
    )
parser.add(
    '-pw',
    '--mysql-password',
    type=str,
    dest='dbpass',
    required=True,
    default=None,
    help='The password for the MySQL username with permission to upload to MinoTour.'
    ,
    )
parser.add(
    '-pin',
    '--security-pin',
    type=str,
    dest='pin',
    required=True,
    default=None,
    help='This is a security feature to prevent unauthorised remote control of a minION device. \
    You need to provide a four digit pin number which must be entered on the website to remotely \
    control the minION.'
    ,
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
global args
args = parser.parse_args()

version = '0.3'  # 9th January 2017

### test which version of python we're using

###Machine to connect to address
global ipadd
ipadd = args.ip


def bytesto(bytes, to, bsize=1024):
    """convert bytes to megabytes, etc.
       sample code:
           print('mb= ' + str(bytesto(314575262000000, 'm')))
       sample output:
           mb= 300002347.946
    """

    a = {'k' : 1, 'm': 2, 'g' : 3, 't' : 4, 'p' : 5, 'e' : 6}
    r = float(bytes)
    for i in range(a[to]):
        r = r / bsize

    return r

def checkpin(dbname):
    #if check_db_connection() == 1:
    try:
        checkpin = "select message from %s.messages where message = 'pin'" % (dbname)
        print checkpin
        cursor=db.query(checkpin)
        #db.commit()
        rows = cursor.fetchall()
        return len(rows)
    except Exception,err:
        print "line 173 Error",err
        check_db_connection()
        print "Database problem"
        return 0

def checkwarning(dbname):
        #if check_db_connection() == 1:
        try:
            checkwarn = "select name from %s.alerts where name = 'minKNOWwarning'" % (dbname)
            cursor=db.query(checkwarn)
            #db.commit()
            rows = cursor.fetchall()
            return len(rows)
        except Exception,err:
            print "line 187 Error",err
            check_db_connection()
            print "Database problem"
            return 0



def check_email_twitter():
    #if check_db_connection() == 1:

    try:
        getuserinfo = "select twitnote,emailnote from Gru.users where user_name = '%s'" %(args.dbusername)
        #print getuserinfo
        cursor=db.query(getuserinfo)
        #print cursor
        #db.commit()
        row=cursor.fetchone()
    except Exception,err:
        print "line 202 Error",err


def save_data():
    print "saving data"
    clean_up()
    datatosave=sendingminIONdict
    for thing in datatosave["DETAILS"]:
        for minion in datatosave["DETAILS"][thing]:
            if "detailsdata" in datatosave["DETAILS"][thing][minion]:
                if "status" in datatosave["DETAILS"][thing][minion]["detailsdata"]["engine_states"]:
                    insert_data = 'insert into Gru.json_store (user_name,date_time,json_string,minion,runname) VALUES ("%s","%s","%s","%s","%s")' % (args.dbusername,time.time(),(MySQLdb.escape_string(str(json.dumps(datatosave)))),minion,str(datatosave["DETAILS"][thing][minion]["livedata"]["dataset"]["result"]))
                    try:
                        cursor=db.query(insert_data)
                        #db.commit()
                    except Exception,err:
                            check_db_connection()
                            print "line 219 Error",err


def check_databases():
    print "CHECK DATABASES"
    for minION in minIONdict: # if yes then we check each minION to see what is happening
        #print minION, minIONdict[minION]["state"]
        if minIONdict[minION]["state"] == "active":
            #look and see if there are any active runs with this user for this minION.
            #print "MinION is active - fetching databases"
            fetchruns = "SELECT runname FROM Gru.minIONruns where activeflag = 1 and flowcellid = '%s' and user_name = '%s'" % (minION,args.dbusername)
            print fetchruns
            #if check_db_connection() == 1:
            try:
                cursor=db.query(fetchruns)
                #db.commit()
                rows = cursor.fetchall()
                print rows
                for row in rows:
                    print "ROW", row[0]
                    print "CHECKPIN", checkpin(row[0])
                    if checkpin(row[0]) < 1:
                        print "TRYING PIN INSERT"
                        #pininsert = \
                        #"insert into %s.messages (message,target,param1,complete) VALUES ('pin','all','%s','1')" \
                        #% (row[0],hashlib.md5(args.pin).hexdigest())
                        pininsert = \
                        "insert into %s.messages (message,target,param1,complete) VALUES ('pin','all','%s','1')" \
                        % (row[0],args.pin)

                        print "pininsert",pininsert
                        cursor=db.query(pininsert)
                        #db.commit()
                    if checkwarning(row[0]) < 1:
                        #print "Checking Drive Warnings from MinKNOW",type(minIONdict[minION]["livedata"]['disk_space']['result'][0]['recommend_alert'])
                        if minIONdict[minION]["livedata"]['disk_space']['result'][0]['recommend_alert'] == True:
                            warninginsert = "insert into %s.alerts (name,type,complete) VALUES ('minKNOWwarning','%s',0)" % (row[0],MySQLdb.escape_string(minIONdict[minION]["livedata"]["machine_id"]['result'])) # add in minion.livedata.machine_id.result as "type"
                            print warninginsert
                            cursor=db.query(warninginsert)
                            #db.commit()
                            print "DANGER DISK WARNING"
                            globalwarningdict[minION]=1
                        if minION in globalwarningdict.keys() and minIONdict[minION]["livedata"]['disk_space']['result'][0]['recommend_alert'] == False:
                            deletewarnings = "delete * from %s.alerts where name='minKNOWwarning'" %(row[0])
                            cursor=db.query(deletewarnings)
                            print "Disk Warning Removed"
                            globalwarningdict.pop(minION)
                    fetchalerts = "select * from %s.alerts" % (row[0])
                    #print fetchalerts
                    cursor=db.query(fetchalerts)
                    #db.commit()
                    rows2 = cursor.fetchall()
                    for row2 in rows2:
                        #print row2
                        if row2[1] == "disknotify":
                            #print bytesto(minIONdict[minION]["livedata"]['disk_space']['result'][0]['bytes_available'],"g")
                            if bytesto(minIONdict[minION]["livedata"]['disk_space']['result'][0]['bytes_available'],"g") < row2[7]:
                                try:
                                    update = "update %s.alerts set end = 1, type = '%s' where alert_index = %s" % (row[0],MySQLdb.escape_string(minIONdict[minION]["livedata"]["machine_id"]['result']),row2[0])
                                    #print update
                                    cursor=db.query(update)
                                    #db.commit()
                                except Exception,err:
                                    print "line 267 Error",err
                    #check for things to stop the run...
                    fetchinteractions = "select * from %s.interaction where complete != 1" %(row[0])
                    cursor=db.query(fetchinteractions)
                    #db.commit()
                    rows3 = cursor.fetchall()
                    for row3 in rows3:
                        if row3[1] == "stop":
                            #print "THING TO STOP"
                            stoprun(minIONdict[minION]["port"])
                            try:
                                update = "update %s.interaction set complete = 1 where job_index = %s" % (row[0],row3[0])
                                #print update
                                cursor=db.query(update)
                                #db.commit()
                            except Exception,err:
                                print "Error",err
            except Exception,err:
                check_db_connection()
                print "Error",err




    #sys.exit()
    pass
    ##First get the list of running minIONs


def _urlopen(url, *args):
    """Open a URL, without using a proxy for localhost.
    While the no_proxy environment variable or the Windows "Bypass proxy
    server for local addresses" option should be set in a normal proxy
    configuration, the latter does not affect requests by IP address. This
    is apparently "by design" (http://support.microsoft.com/kb/262981).
    This method wraps urllib2.urlopen and disables any set proxy for
    localhost addresses.
    """
    try:
        host = url.get_host().split(':')[0]
    except AttributeError:
        host = urlparse.urlparse(url).netloc.split(':')[0]
    import socket
    # NB: gethostbyname only supports IPv4
    # this works even if host is already an IP address
    addr = socket.gethostbyname(host)
    #if addr.startswith('127.'):
    #    return _no_proxy_opener.open(url, *args)
    #else:
    return urllib2.urlopen(url, *args)


def execute_command_as_string(data, host=None, port=None):
    host_name = host
    port_number = port
    #print host,port,data
    url = 'http://%s:%s%s' % (host_name, port_number, '/jsonrpc')
    #req = urllib2.Request(url, data=data,headers={'Authorization': 'Bearer abc','Content-Length': str(len(data)),'Content-Type': 'application/json'})
    req = urllib2.Request(url, data=data,headers={'Content-Length': str(len(data)),'Content-Type': 'application/json'})
    #print req
    try:
        f = _urlopen(req)
        json_respond = json.loads(f.read())
        f.close()
        return json_respond
    except Exception, err:
        print err
        err_string = \
            'Fail to initialise mincontrol. Likely reasons include minKNOW not running, the wrong IP address for the minKNOW server or firewall issues.'

    return



def send_message_port(message,port):
    message_to_send = \
        '{"id":"1", "method":"user_message","params":{"content":"%s"}}' \
        % message
    results=""
    try:
        results = execute_command_as_string(message_to_send, ipadd, port)
    except Exception, err:
        print "message send fail", err
    return results

def send_message(message):
    message_to_send = \
        '{"id":"1", "method":"user_message","params":{"content":"%s"}}' \
        % message
    results = execute_command_as_string(message_to_send, ipadd, 8000)
    return results

def startstop(command,minION):
    if OPER == "osx":
        p = os.popen('/Applications/MinKNOW.app/Contents/Resources/bin/mk_manager_client -i ' + minION + ' --' + command,"r")
        while 1:
            line = p.readline()
            if not line: break
            #print line
    elif OPER == "windows":
        p = os.popen('C:\\grouper\\binaries\\bin\\mk_manager_client.exe -i ' + minION + ' --' + command, "r")
        #print 'C:\\grouper\\binaries\\bin\\mk_manager_client.exe -i ' + minION + ' --' + command, "r"
        while 1:
            line = p.readline()
            if not line: break
            #print line
    elif OPER == "linux":
        print "!!!!!!!!!!!!!! Sorry cannot handle linux yet."
    else:
        print "!!!!!!!!!!!!!! Sorry - cannot recognise your operating system."

def commands(command):
    return {
        'getstaticdata' : '{"id":1,"method":"get_static_data","params":null}',
        'initialization_status' : '{"params": "null", "id": 5, "method": "initialization_status"}',
        'get_analysis_configuration' : '{"id":1,"method":"get_analysis_configuration","params":null}',
        'initialiseminion' : '{"params": {"command": "init_main_board", "parameters": []}, "id": 0, "method": "board_command_ex"}',
        'shutdownminion' : '{"params": {"state_id": "status", "value": "stop"}, "id": 0, "method": "set_engine_state"}',
        'startmessagenew' : '{"id":"1", "method":"user_message","params":{"content":"minoTour is now interacting with your run. This is done at your own risk. To stop minoTour interaction with minKnow disable upload of read data to minoTour."}}',
        'status':'{"id":"1", "method":"get_engine_state","params":{"state_id":"status"}}',
        'dataset':'{"id":"1", "method":"get_engine_state","params":{"state_id":"data_set"}}',
        'startrun' :'{"id":"1", "method":"start_script","params":{"name":"MAP_Lambda_Burn_In_Run_SQK_MAP005.py"}}',
        'stoprun' : '{"id":"1", "method":"stop_experiment","params":"null"}',
        'stopprotocol' : '{"id":"1", "method":"stop_script","params":{"name":"MAP_48Hr_Sequencing_Run.py"}}',
        'biasvoltageget' : '{"id":"1","method":"board_command_ex","params":{"command":"get_bias_voltage"}}',
        'bias_voltage_gain' : '{"id":"1","method":"get_engine_state","params":{"state_id":"bias_voltage_gain"}}',
        'bias_voltage_set' :'{"id":"1","method":"board_command_ex","params":{"command":"set_bias_voltage","parameters":"-120"}}',
        'machine_id' :'{"id":"1","method":"get_engine_state","params":{"state_id":"machine_id"}}',
        'machine_name' :'{"id":"1","method":"get_engine_state","params":{"state_id":"machine_name"}}',
        'sample_id' :'{"id":"1","method":"get_engine_state","params":{"state_id":"sample_id"}}',
        'flow_cell_id' : '{"id":"1","method":"get_engine_state","params":{"state_id":"flow_cell_id"}}',
        'user_error' :'{"id":"1","method":"get_engine_state","params":{"state_id":"user_error"}}',
        'sequenced_res' :'{"id":"1","method":"get_engine_state","params":{"state_id":"sequenced"}}',
        'yield_res' : '{"id":"1","method":"get_engine_state","params":{"state_id":"yield"}}',
        'current_script' : '{"id":"1","method":"get_engine_state","params":{"state_id":"current_script"}}',
        'get_scripts' : '{"id":"1", "method":"get_script_info_list","params":{"state_id":"status"}}',
        'disk_space' : '{"id":1,"method":"get_disk_space_info","params":null}',
        'sinc_delay' : '{"id":1,"method":"sinc_delay","params":null}',
        'get_seq_metrics' : '{"id":1,"method": "get_seq_metrics","params":null}',
        'get_tracking_id': '{"id":1,"method": "get_tracking_id","params":null}',
        'get_statistics': '{"id":1,"method": "get_statistics","params":null}',

    }[command]





class HelpTheMinion(WebSocketClient):

#    def __init__(self,minIONdict):
#        self.minIONdict = minIONdict

    def opened(self):
        print "Connected to Master MinION Controller!"

    def initialiseminion():
        result = self.send(json.dumps({"params": "null", "id": 5, "method": "initialization_status"}))
        #print result


    def received_message(self, m):
        ##print "message received!"
        for thing in ''.join(map(chr,map(ord,str(m)))).split('\n'):
            #if len(thing) > 5 and "2L" not in thing and "2n" not in thing:
            if len(thing) > 5 and thing[1] == "M":
                if thing[1:8] not in minIONdict:
                    #print "INITIALISING MINIONDICT"
                    minIONdict[thing[1:8]]=dict()
                print "minION ID:", thing[1:8]
                ##print map(ord,thing)
                minIONports =  map(lambda x:x-192+8000,filter(lambda x:x>190,map(ord,thing)))
                ##print minIONports
                if len(minIONports) > 0:
                    minIONdict[thing[1:8]]["state"]="active"
                    port = minIONports[0]
                    ws_longpoll_port = minIONports[1]
                    ws_event_sampler_port = minIONports[2]
                    ws_raw_data_sampler_port = minIONports[3]
                    minIONdict[thing[1:8]]["port"]=port
                    minIONdict[thing[1:8]]["ws_longpoll_port"]=ws_longpoll_port
                    minIONdict[thing[1:8]]["ws_event_sampler_port"]=ws_event_sampler_port
                    minIONdict[thing[1:8]]["ws_raw_data_sampler_port"]=ws_raw_data_sampler_port
                else:
                    minIONdict[thing[1:8]]["state"]="inactive"
                    minIONdict[thing[1:8]]["port"]=""
                    minIONdict[thing[1:8]]["ws_longpoll_port"]=""
                    minIONdict[thing[1:8]]["ws_event_sampler_port"]=""
                    minIONdict[thing[1:8]]["ws_raw_data_sampler_port"]=""


class CallingHome(WebSocketClient):

    def opened(self):
        print "Connection Success - Calling Home."
        self.connected=True
        self.messagesent=False

    def setup(self, timeout=1):
        try:
            self.__init__(self.url)
            self.connect()
            self.connected=True
            self.timeout=2
            self.run_forever()
        except KeyboardInterrupt:
            self.connected=False
            self.messagesent=False
            self.close()
        except:
            self.connected=False
            self.messagesent=False
            newTimeout = timeout + 1
            print("Timing out for %i seconds. . ." % newTimeout)
            time.sleep(newTimeout)
            print("Attempting reconnect. . .")
            self.setup(newTimeout)


    def received_message(self, m):
        #print "message received!",m
        if str(m) != "Connected - waiting for messages.":
            #print "we're off"
            try:
                message = json.loads(str(m))
                if type(message) is dict:
                    ##print "message",m,
                    if message["JOB"] == "testmessage":
                        send_message_port("minoTour is checking communication status with "+str(message["minion"])+".",minIONdict[message["minion"]]["port"])
                    if message["JOB"] =="biasvoltageinc":
                        biasvoltsmod("inc",minIONdict[message["minion"]]["port"])
                    if message["JOB"] =="biasvoltagedec":
                        biasvoltsmod("dec",minIONdict[message["minion"]]["port"])
                    if message["JOB"] == "rename":
                        renamerun(message["NAME"],minIONdict[message["minion"]]["port"])
                    if message["JOB"] == "nameflowcell":
                        renameflowcell(message["NAME"],minIONdict[message["minion"]]["port"])
                    if message["JOB"] == "stopminion":
                        stoprun(minIONdict[message["minion"]]["port"])
                    if message["JOB"] == "startminion":
                        startrun(message["SCRIPT"],minIONdict[message["minion"]]["port"])
                    if message["JOB"] == "initialiseminion":
                        #print "Trying to initialise minION"
                        #result = helper.initialiseminion()
                        startstop('start',message["minion"])
                        print "Try to run the initialiseminion command now."
                        #execute_command_as_string(commands('initialiseminion'), ipadd,minIONdict[message["minion"]]["port"])
                    if message["JOB"] == "shutdownminion":
                        #First we need to stop the minion if it is running...
                        #stoprun(minIONdict[message["minion"]]["port"])
                        #execute_command_as_string(commands('shutdownminion'), ipadd,minIONdict[message["minion"]]["port"])
                        startstop('stop',message["minion"])
                        print "Try to shut down the minION now."
            except Exception, err:
                print "line 508",err,str(m)
                pass
            #if type(message) is dict:
        #    #print "message"
            #try:
            #    #print message["JOB"]
            #except:
            #    pass
        #self.send(json.dumps("message received!", ensure_ascii=False))

    def closed(self, code, reason=None):
        print self.sock.fileno()
        print "Closed down Home", code, reason
        self.connected=False
        self.messagesent=False
        #self.runanalysis = True
        #minIONdict=[]
        time.sleep(3)
        print ("Reconnecting. . .")
        # self.sock.shutdown(socket.SHUT_RDWR)
        #self.sock.close()
        self.setup()


def startrun(script,port):
    #print 'Starting the minION device.'
    try:
        startruncustom = \
            '{"id":1, "method":"start_script","params":{"name":"' \
            + script + '"}}'
        print startruncustom,ipadd,port
        startresult = \
            execute_command_as_string(startruncustom,
                ipadd, port)
        startrunmessage = 'minoTour sent a remote run start command.'
        startresultmessage = send_message_port(startrunmessage,port)
    except Exception, err:
        print >> sys.stderr, err

def stoprun_message(port,message):
    try:
        stopresult = execute_command_as_string(commands('stoprun'),
                ipadd, port)
        stopprotocolresult = \
            execute_command_as_string(commands('stopprotocol'), ipadd,
                port)
        stoprunmessage = 'minoTour sent a remote run stop command for %s.' % (message)
        stopresultmessage = send_message_port(stoprunmessage,port)
    except Exception, err:
        print >> sys.stderr, err

def send_custom_message(port,message):
    try:
        custommessage = "minoTour: %s" % (message)
        custommessageresult=send_message_port(custommessage,port)
    except Exception,err:
        print >> sys.stderr, err

def stoprun(port):
    try:
        stopresult = execute_command_as_string(commands('stoprun'),
                ipadd, port)
        stopprotocolresult = \
            execute_command_as_string(commands('stopprotocol'), ipadd,
                port)
        stoprunmessage = 'minoTour sent a remote run stop command.'
        stopresultmessage = send_message_port(stoprunmessage,port)
    except Exception, err:
        print >> sys.stderr, err


def renameflowcell(name,port):
    try:
        set_flowcell_id = \
            '{"id":1,"method":"set_engine_state","params":{"state_id":"flow_cell_id","value":"'+name+'"}}'
        set_sample_id_result = \
            execute_command_as_string(set_flowcell_id,ipadd,port)
        startresultmessage = \
            send_message_port('minoTour renamed the flowcell to ' + name,port)
    except Exception, err:
        print >> sys.stderr, err

def renamerun(name,port):
    try:
        set_sample_id = \
            '{"id":"1","method":"set_engine_state","params":{"state_id":"sample_id","value":"' \
            + name + '"}}'
        set_sample_id_result = \
            execute_command_as_string(set_sample_id, ipadd,
                port)
        startresultmessage = \
            send_message_port('minoTour renamed the run to '
                + name,port)
    except Exception, err:
        print >> sys.stderr, err


def biasvoltsmod(direction,port):
    if direction == "inc":
        #print 'Incrementing Bias Voltage'
        try:
            biasresultmessage = \
                execute_command_as_string(commands('biasvoltageget'),
                    ipadd, port)
            biasvoltageoffset = \
                execute_command_as_string(commands('bias_voltage_gain'),
                    ipadd, port)
            curr_voltage = int(biasresultmessage['result'
                    ]['bias_voltage']) \
                * int(biasvoltageoffset['result']) + 10
            bias_voltage_inc = \
                '{"id":"1","method":"board_command_ex","params":{"command":"set_bias_voltage","parameters":"%s"}}' \
                % curr_voltage
            biasvoltagereset = \
                execute_command_as_string(bias_voltage_inc,
                    ipadd, port)
            biasresultmessage = \
                execute_command_as_string(commands('biasvoltageget'),
                    ipadd, port)
            biasvoltageoffset = \
                execute_command_as_string(commands('bias_voltage_gain'),
                    ipadd, port)
            curr_voltage = int(biasresultmessage['result'
                    ]['bias_voltage']) \
                * int(biasvoltageoffset['result'])
            incmessage = 'minoTour shifted the bias voltage by +10 mV.'
            incresultmessage = send_message_port(incmessage,port)
        except Exception, err:
            print >> sys.stderr, err
    if direction == "dec":
        #print 'Decreasing Bias Voltage'
        try:
            biasresultmessage = \
                execute_command_as_string(commands('biasvoltageget'),
                    ipadd, port)
            biasvoltageoffset = \
                execute_command_as_string(commands('bias_voltage_gain'),
                    ipadd, port)
            curr_voltage = int(biasresultmessage['result'
                    ]['bias_voltage']) \
                * int(biasvoltageoffset['result']) - 10
            bias_voltage_dec = \
                '{"id":"1","method":"board_command_ex","params":{"command":"set_bias_voltage","parameters":"%s"}}' \
                % curr_voltage
            biasvoltagereset = \
                execute_command_as_string(bias_voltage_dec,
                    ipadd, port)
            biasresultmessage = \
                execute_command_as_string(commands('biasvoltageget'),
                    ipadd, port)
            biasvoltageoffset = \
                execute_command_as_string(commands('bias_voltage_gain'),
                    ipadd, port)
            curr_voltage = int(biasresultmessage['result'
                    ]['bias_voltage']) \
                * int(biasvoltageoffset['result'])
            decmessage = 'minoTour shifted the bias voltage by -10 mV.'
            decresultmessage = send_message_port(decmessage,port)
        except Exception, err:
            print >> sys.stderr, err

class MessagesClient(WebSocketClient):
    def __init__(self, *args,**kwargs):
        super(MessagesClient, self).__init__(*args,**kwargs)
        print "MessagesClient established!"
        self.detailsdict=dict()
        self.daemon=True

    def init_minion(self,minion):
        self.minion=minion

    def opened(self):
        print "MessagesClient Opened!"

    def closed(self, code, reason="MessagesClient disconnected for some reason"):
        print "socket",self.sock.fileno()
        print "Closed down", code, reason

    def received_message(self,m):
        if not m.is_binary:
            print m
            try:
                if self.minion not in minIONdict.keys():
                    minIONdict[self.minion]=dict()
                if "messages" not in minIONdict[self.minion].keys():
                    minIONdict[self.minion]["messages"]=[]
                minIONdict[self.minion]["messages"].append(json.loads(str(m)))
                send_data()
            except:
                pass


class DummyClient(WebSocketClient):

    def __init__(self, *args,**kwargs):
        super(DummyClient, self).__init__(*args,**kwargs)
        print "Client established!"
        self.detailsdict=dict()
        self.daemon=True

    def init_minion(self,minion):
        """ ejecteject has three states:
            0=initisalised or sequencing finished
            1=sequencing happening
            2=run finish detected."""
        self.minion=minion
        self.ejecteject=0
        self.runname=""


    def opened(self):
        self.send(json.dumps({'engine_states':'1','channel_states':'1','channel_info':'1'}))

    def closed(self, code, reason="client disconnected for some reason"):
        print "socket",self.sock.fileno()
        print "Closed down", code, reason

    def received_message(self, m):
        if not m.is_binary:

            #ejecteject=0
            json_object = json.loads(str(m))
            #print json_object
            for element in json_object:
                if element == "statistics" and json_object[element] != "null":
                    for element2 in json_object[element]:
                        if  element2 == "channel_info": #element2 == "read_statistics" or
                            print "found ",element2
                            print json_object[element][element2]
                        if element2 != "read_statistics" and element2 != "channel_info" and json_object[element][element2] != "null":# and element2 != "read_statistics" and element2 != "completed_samples":
                            if element not in self.detailsdict:
                                self.detailsdict[element]=dict()
                            self.detailsdict[element][element2]=json_object[element][element2]
                elif element == "engine_states" and json_object[element] != "null":
                    for element2 in json_object[element]:
                        if json_object[element][element2] != "null":
                            if element not in self.detailsdict:
                                self.detailsdict[element]=dict()
                            if json_object[element][element2] is not dict:
                                self.detailsdict[element][element2]=json_object[element][element2]
                    if "status" in json_object[element]:
                        print "engine states",json_object[element]["status"]
                        if json_object[element]["status"]== "finishing":
                            print "Looks as though a run is finishing!!!!"
                            self.ejecteject=2
                        if json_object[element]["status"]== "starting":
                            print "Looks as though a run is starting!!!!"
                            if self.ejecteject==0:
                                print "resetting history"
                                if self.minion in minIONdict and "yield_history" in minIONdict[self.minion]:
                                    minIONdict[self.minion]["yield_history"]=[]
                                    minIONdict[self.minion]["temp_history"]["asictemp"]=[]
                                    minIONdict[self.minion]["temp_history"]["heatsinktemp"]=[]
                                    minIONdict[self.minion]["temp_history"]["voltage"]=[]
                                    minIONdict[self.minion]["pore_history"]["strand"]=[]
                                    minIONdict[self.minion]["pore_history"]["percent"]=[]
                                    minIONdict[self.minion]["pore_history"]["single"]=[]
                                    if "details" in minIONdict[self.minion]["pore_history"]:
                                        minIONdict[self.minion]["pore_history"].pop("details")
                        if json_object[element]["status"]=="processing" and self.ejecteject==0:
                            print "we're off - resetting history"
                            if self.minion in minIONdict and "yield_history" in minIONdict[self.minion]:
                                minIONdict[self.minion]["yield_history"]=[]
                                minIONdict[self.minion]["temp_history"]["asictemp"]=[]
                                minIONdict[self.minion]["temp_history"]["heatsinktemp"]=[]
                                minIONdict[self.minion]["temp_history"]["voltage"]=[]
                                minIONdict[self.minion]["pore_history"]["strand"]=[]
                                minIONdict[self.minion]["pore_history"]["percent"]=[]
                                minIONdict[self.minion]["pore_history"]["single"]=[]
                                if "details" in minIONdict[self.minion]["pore_history"]:
                                    minIONdict[self.minion]["pore_history"].pop("details")
                            self.ejecteject=1
                elif element == "channel_info" and json_object[element] != "null":
                    ### Getting and updating channel information.
                    if "channels" in json_object[element].keys():
                        for channel in json_object[element]["channels"]:
                            state = "unknown"
                            if "state" in channel.keys():
                                state = channel["state"]
                            logitem(int(channel["name"]),self.minion,state)
                    for element2 in json_object[element]:
                        #pass
                        #turns out we use this data elsewhere - so need to clean it out at send point.
                        if json_object[element][element2] != "null":
                            if element not in self.detailsdict:
                                self.detailsdict[element]=dict()
                            if json_object[element][element2] is not dict:
                                self.detailsdict[element][element2]=json_object[element][element2]
                else:
                    if json_object[element] != "null":
                        #print element
                        self.detailsdict[element]=json_object[element]
            if self.ejecteject==2 and len(self.runname)>0:
                save_data()
                print "resetting history 773"
                try:
                    if self.minion in minIONdict and "yield_history" in minIONdict[self.minion]:
                        minIONdict[self.minion]["yield_history"]=[]
                        minIONdict[self.minion]["temp_history"]["asictemp"]=[]
                        minIONdict[self.minion]["temp_history"]["heatsinktemp"]=[]
                        minIONdict[self.minion]["temp_history"]["voltage"]=[]
                        minIONdict[self.minion]["pore_history"]["strand"]=[]
                        minIONdict[self.minion]["pore_history"]["percent"]=[]
                        minIONdict[self.minion]["pore_history"]["single"]=[]
                        if "details" in minIONdict[self.minion]["pore_history"]:
                            minIONdict[self.minion]["pore_history"].pop("details")
                    self.ejecteject=0
                except Exception, err:
                    print "Problem",err
                    print "Probably already popped"

def get_run_scripts(port):
    #print "trying to get runscripts on port %s" %(port)
    get_scripts = \
        '{"id":"1", "method":"get_script_info_list","params":{"state_id":"status"}}'
    results = execute_command_as_string(get_scripts, ipadd, port)
    #print type(results['result'])
    #print len(results['result'])
    #for thing in results['result']:
    #    print thing
    return results['result']
    for key in results.keys():
        print "mincontrol:", key, results[key]
    scriptlist = list()
    for element in results['result']:
        for item in results['result'][element]:
            scriptlist.append("('runscript','all','" + item['name']
                              + "',1)")

class RepeatedTimer(object):
    def __init__(self,interval,function,*args,**kwargs):
        self._timer = None
        self.interval = interval
        self.function = function
        self.args = args
        self.kwargs = kwargs
        self.is_running = False
        self.start()

    def _run(self):
        self.is_running = False
        self.start()
        self.function(*self.args,**self.kwargs)

    def start(self):
        if not self.is_running:
            self._timer = threading.Timer(self.interval,self._run)
            self._timer.start()
            self.is_running = True

    def stop(self):
        self._timer.cancel()
        self.is_running = False



def process_tracked_yield():
    """
    This function maintains a history of all yield values over the time of a run - it can be configured to store yields
    at given intervals. It is used by the main web page to provide a history of cumulative yield.
    """
    #print "Process_tracked_yield Running"
    #print minIONdict
    for minion in minIONdict:
        if "detailsdata" in minIONdict[minion]:
            if "livedata" in minIONdict[minion]:
                if "yield_res" in minIONdict[minion]["livedata"]:
                    if "result" in minIONdict[minion]["livedata"]["yield_res"]:
                        try:
                            asictemp = float(minIONdict[minion]["detailsdata"]["engine_states"]["minion_asic_temperature"])
                        except:
                            asictemp = 0
                        try:
                            heatsinktemp = float(minIONdict[minion]["detailsdata"]["engine_states"]["minion_heatsink_temperature"])
                        except:
                            heatsinktemp = 0
                        if minIONdict[minion]["livedata"]["bias_voltage_gain"]["result"] != "null":
                            #print "!!!!!!!!STUFF!!!!!", minIONdict[minion]["livedata"]["bias_voltage_gain"]["result"]
                            biasvoltage = int(minIONdict[minion]["livedata"]["bias_voltage_gain"]["result"])
                        else:
                            biasvoltage = 0
                        if minIONdict[minion]["livedata"]["biasvoltageget"]["result"] != "null":
                            voltage_val = int(minIONdict[minion]["livedata"]["biasvoltageget"]["result"]["bias_voltage"])
                        else:
                            voltage_val = 0
                        voltage_value = voltage_val * biasvoltage
                        #print "volts:",voltage_value
                        try:
                            if minIONdict[minion]["livedata"]["yield_res"]["result"] != "null":
                                yieldval = int(minIONdict[minion]["livedata"]["yield_res"]["result"])
                            else:
                                yieldval =  minIONdict[minion]["detailsdata"]["statistics"]["read_event_count"]
                        except:
                            yieldval=0
                        try:
                            meanratio = minIONdict[minion]["detailsdata"]["meanratio"]
                            openpore = minIONdict[minion]["detailsdata"]["openpore"]
                            instrand = minIONdict[minion]["detailsdata"]["instrand"]
                        except:
                            meanratio = 0
                            openpore = 0
                            instrand = 0

                        if "yield_history" not in minIONdict[minion]:
                            minIONdict[minion]["yield_history"]=[]
                        if "temp_history" not in minIONdict[minion]:
                            minIONdict[minion]["temp_history"]=dict()
                        if "pore_history" not in minIONdict[minion]:
                            minIONdict[minion]["pore_history"]=dict()
                        if "meanratio_history" not in minIONdict[minion]["pore_history"]:
                            minIONdict[minion]["pore_history"]["meanratio_history"]=[]
                            minIONdict[minion]["pore_history"]["openpore_history"] = []
                            minIONdict[minion]["pore_history"]["instrand_history"] = []
                        if "strand" not in minIONdict[minion]["pore_history"]:
                            minIONdict[minion]["pore_history"]["strand"]=[]
                            minIONdict[minion]["pore_history"]["percent"]=[]
                            minIONdict[minion]["pore_history"]["single"]=[]
                        if "asictemp" not in minIONdict[minion]["temp_history"]:
                            minIONdict[minion]["temp_history"]["asictemp"]=[]
                            minIONdict[minion]["temp_history"]["heatsinktemp"]=[]
                            minIONdict[minion]["temp_history"]["voltage"]=[]
                        strand = 0
                        single_pore = 0
                        if "simplesummary" in minIONdict[minion]:
                            #for state in minIONdict[minion]["simplesummary"].keys():
                            #    print state,minIONdict[minion]["simplesummary"][state]
                            if "strand" in minIONdict[minion]["simplesummary"].keys():
                                strand = float(minIONdict[minion]["simplesummary"]["strand"])
                            if "single_pore" in minIONdict[minion]["simplesummary"].keys():
                                single_pore = float(minIONdict[minion]["simplesummary"]["single_pore"])
                            elif "good_single" in minIONdict[minion]["simplesummary"].keys():
                                single_pore = float(minIONdict[minion]["simplesummary"]["good_single"])
                        if (strand+single_pore)> 0:
                            percent = (strand/(strand+single_pore)*100)
                        else:
                            percent = 0
                        #print "testing length"
                        if len(minIONdict[minion]["yield_history"]) >1:
                            if yieldval > minIONdict[minion]["yield_history"][-1][1]:
                                minIONdict[minion]["yield_history"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,yieldval))
                            minIONdict[minion]["temp_history"]["asictemp"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,asictemp))
                            minIONdict[minion]["temp_history"]["voltage"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,voltage_value))
                            minIONdict[minion]["temp_history"]["heatsinktemp"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,heatsinktemp))
                            minIONdict[minion]["pore_history"]["strand"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,strand))
                            minIONdict[minion]["pore_history"]["percent"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,percent))
                            minIONdict[minion]["pore_history"]["single"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,single_pore))
                            minIONdict[minion]["pore_history"]["meanratio_history"].append(
                                (minIONdict[minion]["detailsdata"]["timestamp"] * 1000, meanratio))
                            minIONdict[minion]["pore_history"]["openpore_history"].append(
                                (minIONdict[minion]["detailsdata"]["timestamp"] * 1000, openpore))
                            minIONdict[minion]["pore_history"]["instrand_history"].append(
                                (minIONdict[minion]["detailsdata"]["timestamp"] * 1000, instrand))

                            if "simplesummary" in minIONdict[minion]:
                                if "details" not in minIONdict[minion]["pore_history"].keys():
                                    minIONdict[minion]["pore_history"]["details"]=dict()
                                for state in minIONdict[minion]["simplesummary"].keys():
                                    if state not in minIONdict[minion]["pore_history"]["details"].keys():
                                        minIONdict[minion]["pore_history"]["details"][state]=[]
                                    #print "appending at 899"
                                    minIONdict[minion]["pore_history"]["details"][state].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,minIONdict[minion]["simplesummary"][state]))
                        else:
                            minIONdict[minion]["yield_history"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,yieldval))
                            minIONdict[minion]["temp_history"]["asictemp"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,asictemp))
                            minIONdict[minion]["temp_history"]["heatsinktemp"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,heatsinktemp))
                            minIONdict[minion]["temp_history"]["voltage"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,voltage_value))
                            minIONdict[minion]["pore_history"]["strand"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,strand))
                            minIONdict[minion]["pore_history"]["percent"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,percent))
                            minIONdict[minion]["pore_history"]["single"].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,single_pore))
                            minIONdict[minion]["pore_history"]["meanratio_history"].append(
                                (minIONdict[minion]["detailsdata"]["timestamp"] * 1000, meanratio))
                            minIONdict[minion]["pore_history"]["openpore_history"].append(
                                (minIONdict[minion]["detailsdata"]["timestamp"] * 1000, openpore))
                            minIONdict[minion]["pore_history"]["instrand_history"].append(
                                (minIONdict[minion]["detailsdata"]["timestamp"] * 1000, instrand))
                            if "simplesummary" in minIONdict[minion]:
                                if "details" not in minIONdict[minion]["pore_history"].keys():
                                    minIONdict[minion]["pore_history"]["details"]=dict()
                                for state in minIONdict[minion]["simplesummary"].keys():
                                    if state not in minIONdict[minion]["pore_history"]["details"].keys():
                                        minIONdict[minion]["pore_history"]["details"][state]=[]
                                    #print "appending at 915"
                                    minIONdict[minion]["pore_history"]["details"][state].append((minIONdict[minion]["detailsdata"]["timestamp"]*1000,minIONdict[minion]["simplesummary"][state]))
                        #if "details" in minIONdict[minion]["pore_history"]:
                            #print minIONdict[minion]["pore_history"]["details"]
#print "YO DOG"
                else:
                    minIONdict[minion]["yield_history"]=[]
                    minIONdict[minion]["temp_history"]["asictemp"]=[]
                    minIONdict[minion]["temp_history"]["heatsinktemp"]=[]
                    minIONdict[minion]["temp_history"]["voltage"]=[]
                    minIONdict[minion]["pore_history"]["strand"]=[]
                    minIONdict[minion]["pore_history"]["percent"]=[]
                    minIONdict[minion]["pore_history"]["single"]=[]
                    minIONdict[minion]["pore_history"]["openpore_history"] = []
                    minIONdict[minion]["pore_history"]["instrand_history"] = []
                    minIONdict[minion]["pore_history"]["meanratio_history"] = []
                    #minIONdict[minion]["pore_history"]["details"]=dict()
                    if "details" in minIONdict[minion]["pore_history"]:
                        minIONdict[minion]["pore_history"].pop("details")


def logitem(channel,minion,state):
    if minion not in channel_data:
        channel_data[minion]=dict()
        for i in chanlookup:
            channel_data[minion][i]=dict()
    channel_data[minion][channel]["state"]=state

def get_state_summary(minion):
        state_dict = dict()
        for key, value in channel_data[minion].items():
            try:
                if value["state"] in state_dict:
                    state_dict[value["state"]] += 1
                else:
                    state_dict[value["state"]] = 1
            except Exception, err:
                print "line 956",err
                pass
        return state_dict


def send_data():
    ### If the dictionary contains DETAILS data we are going to send a partial update
    ### Otherwise we will send a full update.
    global olddict
    #print "BEFORE",et.messagesent
    if et.messagesent is False:
        print "resetting old dict"
        olddict=dict()
        et.messagesent=True
    lock.acquire()
    #print "AFTER",et.messagesent
    try:
        for job in sendingminIONdict:
            if job == "DETAILS":
                for user in sendingminIONdict[job].keys():
                    deepcopydict=copy.deepcopy(sendingminIONdict[job][user])
                    # here we try and identify if we are sending unnecessary data and clean it out
                    #for bit in deepcopydict:
                    #    if "detailsdata" in deepcopydict[bit].keys():
                    #        deepcopydict[bit]["detailsdata"]["channel_info"].pop("channels")


                    difference = diff(olddict,deepcopydict)
                    difference=list(difference)
                    stufftosend=dict()
                    stufftosend[job]=dict()
                    stufftosend[job][user]=list()
                    stufftosend[job][user]=difference
                    olddict=copy.deepcopy(deepcopydict)
                    #print "XXXXXXXX_new",sys.getsizeof(json.dumps(stufftosend))
                    et.connection.sendMessage(json.dumps(stufftosend))
            else:
                    print "FIRST SEND"
                    #print "XXXXXXXX",sys.getsizeof(json.dumps(sendingminIONdict))
                    et.connection.sendMessage(json.dumps(sendingminIONdict))
    except Exception, err:
        print "Message Send Problem", err
    finally:
        lock.release()
    return(olddict)

def process_channel_information():
    """
    This function maintains the state of the channels in a simple to read format
    """
    colourlookup={}
    for minion in minIONdict:
        if "detailsdata" in minIONdict[minion]:
            if "channelstuff" in minIONdict[minion]:
                for item in minIONdict[minion]["channelstuff"]:
                    if 'style' in minIONdict[minion]["channelstuff"][item]:
                        colourlookup[minIONdict[minion]["channelstuff"][item]['style']['label']]=minIONdict[minion]["channelstuff"][item]['style']['colour']
                if "cust_chan_dict" not in minIONdict[minion]:
                    minIONdict[minion]["cust_chan_dict"]=dict()
                    for i in chanlookup:
                        minIONdict[minion]["cust_chan_dict"][i]=dict()
                for channel in minIONdict[minion]["detailsdata"]["channel_info"]["channels"]:
                    if minion not in statedict:
                        statedict[minion]=dict()
                    if "state_group" in channel:
                        if channel["state"] in colourlookup:
                            statedict[minion].update({channel["name"]:{'state':channel["state"],"state_group":channel["state_group"],"colour":colourlookup[channel["state"]]}})
                            minIONdict[minion]["cust_chan_dict"][int(channel["name"])]=channel["state"]
                        elif channel["state_group"] in colourlookup:
                            statedict[minion].update({channel["name"]:{'state':channel["state"],"state_group":channel["state_group"],"colour":colourlookup[channel["state_group"]]}})
                            minIONdict[minion]["cust_chan_dict"][int(channel["name"])]=channel["state"]
                        else:
                            pass
                statesummarydict[minion]=dict()
                statesummarydict[minion]=get_state_summary(minion)
                if "simplechanstats" not in minIONdict[minion]:
                    minIONdict[minion]["simplechanstats"]={}
                try:
                    minIONdict[minion]["simplesummary"]=statesummarydict[minion]
                except Exception, err:
                    print "line 1081",err
                    pass

def clean_up():
    if minIONdict == minIONdict_test: ## The dictionary is unchanged since the last cycle
        active = 0
        inactive = 0
        for minION in minIONdict:
            if minIONdict[minION]["state"]=="active":
                active += 1
                if minION not in minIONclassdict:
                    minIONclassdict[minION]=dict()
                if "connected" not in minIONclassdict[minION]:
                    connectip = "ws://"+ args.ip + ":"+str(minIONdict[minION]["ws_longpoll_port"])+"/"
                    minIONclassdict[minION]["class"]=DummyClient(connectip)
                    connectip2 = "ws://"+ args.ip + ":"+str(minIONdict[minION]["ws_longpoll_port"])+"/user_messages"
                    minIONclassdict[minION]["class2"]=MessagesClient(connectip2,job="ws_longpoll_port/user_messages")
                try:
                    if "connected" not in minIONclassdict[minION]:
                        try:
                            minIONclassdict[minION]["class"].connect()
                            minIONclassdict[minION]["class2"].connect()
                            print "GETTING THE GOOD STUFF"
                            minIONclassdict[minION]["class"].init_minion(minION)
                            minIONclassdict[minION]["class2"].init_minion(minION)
                            results = execute_command_as_string(commands('get_analysis_configuration'), ipadd,minIONdict[minION]["port"])
                            minIONdict[minION]["channelstuff"]=results["result"]["channel_states"]
                        except Exception, err:
                            print "Connection failed", err
                        minIONclassdict[minION]["connected"]="True"
                    if minION in minIONdict_test:
                        if minIONdict[minION]["state"]==minIONdict_test[minION]["state"]:
                            try: #To catch mysterious bugs
                                livedata=dict() # collect a dictionary of useful data - might change this
                                for query in ('status','dataset','biasvoltageget','bias_voltage_gain','machine_id','machine_name','sample_id','user_error','sequenced_res','yield_res','current_script','disk_space','flow_cell_id'):#,'getstaticdata','get_analysis_configuration'):
                                    results = execute_command_as_string(commands(query), ipadd,minIONdict[minION]["port"])
                                    livedata[query]=results
                                    if query == "disk_space":
                                        #print results
                                        #print results["result"][0]["recommend_stop"]
                                        #print results["result"][0]["recommend_alert"]
                                        check_email_twitter()
                                    #if query == "biasvoltageget":
                                    #    print results
                                minIONdict[minION]["livedata"]=livedata #append results to data stream
                            except Exception, err:
                                print "line 1141",err
                                pass
                        else:
                            if "scripts" in minIONdict[minION]:
                                if minIONdict[minION]["scripts"]!=get_run_scripts(minIONdict[minION]["port"]):
                                    minIONdict[minION]["scripts"]=get_run_scripts(minIONdict[minION]["port"])
                            else:
                                minIONdict[minION]["scripts"]=get_run_scripts(minIONdict[minION]["port"])
                                results = execute_command_as_string(commands('startmessagenew'), ipadd,minIONdict[minION]["port"])
                                try: #To catch mysterious bugs
                                    livedata=dict() # collect a dictionary of useful data - might change this
                                    for query in ('status','dataset','biasvoltageget','bias_voltage_gain','machine_id','machine_name','sample_id','user_error','sequenced_res','yield_res','current_script','disk_space','flow_cell_id'):#,'getstaticdata','get_analysis_configuration'):
                                        results = execute_command_as_string(commands(query), ipadd,minIONdict[minION]["port"])
                                        livedata[query]=results
                                        if query == "disk_space":
                                            #print results
                                            #print results["result"][0]["recommend_stop"]
                                            #print results["result"][0]["recommend_alert"]
                                            check_email_twitter()
                                        #if query == "biasvoltageget":
                                        #    print results
                                    minIONdict[minION]["livedata"]=livedata #append results to data stream
                                except Exception, err:
                                    print "line 1172",err
                                    pass
                except Exception, err:
                    print "line 1177",err
                    print "Connection Error"
                else:
                    if "scripts" not in minIONdict[minION]:
                        minIONdict[minION]["scripts"]=get_run_scripts(minIONdict[minION]["port"])
            else:
                inactive += 1
    else:
        #print "Dictionary has changed"
        for minION in minIONdict: # if yes then we check each minION to see what is happening
            if minIONdict[minION]["state"]=="active": #we have a sequencing minION - what is it doing?
                #print minION, "ACTIVE"
                if minION not in minIONclassdict:
                    connectip = "ws://"+ args.ip + ":"+str(minIONdict[minION]["ws_longpoll_port"])+"/" #Connect to the minIONdict
                    minIONclassdict[minION]=dict() #Add minION to the dictionary
                    connectip2 = "ws://"+ args.ip + ":"+str(minIONdict[minION]["ws_longpoll_port"])+"/user_messages"

                if "connected" not in minIONclassdict[minION]:
                    print "Connecting to this minION",minION
                    minIONclassdict[minION]["class"]=DummyClient(connectip)
                    minIONclassdict[minION]["class2"]=MessagesClient(connectip2,job="ws_longpoll_port/user_messages")
                if "connected" not in minIONclassdict[minION]:
                    try:
                        minIONclassdict[minION]["class"].connect()
                        minIONclassdict[minION]["class2"].connect()
                    except Exception, err:
                        print "Connection failed", err
                    try:
                        print "GETTING THE GOOD STUFF"
                        minIONclassdict[minION]["class"].init_minion(minION)
                        minIONclassdict[minION]["class2"].init_minion(minION)
                        results = execute_command_as_string(commands('get_analysis_configuration'), ipadd,minIONdict[minION]["port"])
                        minIONdict[minION]["channelstuff"]=results["result"]["channel_states"]
                    except Exception, err:
                        print "Connection failed", err
                    minIONclassdict[minION]["connected"]="True"
                livedata=dict() # collect a dictionary of useful data - might change this
                minIONdict[minION]["scripts"]=get_run_scripts(minIONdict[minION]["port"])
                for query in ('status','dataset','biasvoltageget','bias_voltage_gain','machine_id','machine_name','sample_id','user_error','sequenced_res','yield_res','current_script','disk_space','flow_cell_id'):#,'getstaticdata','get_analysis_configuration'):
                    results = execute_command_as_string(commands(query), ipadd,minIONdict[minION]["port"])
                    livedata[query]=results
                    if query == "disk_space":
                        #print results
                        #print results["result"][0]["recommend_stop"]
                        #print results["result"][0]["recommend_alert"]
                        check_email_twitter()
                    #if query == "biasvoltageget":
                        #print results
                results2 = execute_command_as_string(commands('get_analysis_configuration'), ipadd,minIONdict[minION]["port"])
                minIONdict[minION]["channelstuff"]=results2["result"]["channel_states"]
                minIONdict[minION]["livedata"]=livedata #append results to data stream
                if "class" in minIONclassdict[minION]:
                    for element in minIONclassdict[minION]["class"].detailsdict:
                        if "detailsdata" not in minIONdict[minION]:
                            minIONdict[minION]["detailsdata"]=dict()
                        if element == "statistics" and minIONclassdict[minION]["class"].detailsdict[element] != "null":
                            if element not in minIONdict[minION]["detailsdata"]:
                                minIONdict[minION]["detailsdata"][element]=dict()
                            context = dict(list(minIONdict[minION]["detailsdata"][element].items()) + list(minIONclassdict[minION]["class"].detailsdict[element].items()))
                            minIONdict[minION]["detailsdata"][element]=context
                        if element == "engine_states" and minIONclassdict[minION]["class"].detailsdict[element] != "null":
                            for element2 in minIONclassdict[minION]["class"].detailsdict[element]:
                                if minIONclassdict[minION]["class"].detailsdict[element][element2] != "null":
                                    if element not in minIONdict[minION]["detailsdata"]:
                                        minIONdict[minION]["detailsdata"][element]=dict()
                                    minIONdict[minION]["detailsdata"][element][element2]=minIONclassdict[minION]["class"].detailsdict[element][element2]
                        else:
                            if minIONclassdict[minION]["class"].detailsdict[element] != "null":
                                minIONdict[minION]["detailsdata"][element]=minIONclassdict[minION]["class"].detailsdict[element]
            else: #minION state is inactive
                #print minION, "INACTIVE"
                if minION in minIONclassdict:
                    if "connected" in minIONclassdict[minION]:
                        keys = minIONclassdict[minION].keys()
                        for key in keys:
                            print key
                            minIONclassdict[minION].pop(key, None)

    sendingminIONdict["DETAILS"][args.dbusername]=minIONdict


def check_db_connection():
    #global db
    #health = 1
    #print "DATABASE STATE",db.open
    #if db.open==1:
    #    pass
    #else:
    try:
        db = MySQLdb.connect(host=args.dbhost, user=args.dbusername,
                             passwd=args.dbpass, port=args.dbport)
        cursor = db.cursor()
    except Exception, err:
        print >> sys.stderr, "Can't connect to MySQL: %s" % err
        health = 0

class MyClientProtocol(WebSocketClientProtocol):

    """def __init__(self, *args,**kwargs):
        print "I'm initialised"
        self.messagesent=False
    """
    def onConnect(self, response):
        print("Server connected: {0}".format(response.peer))
        self.factory.resetDelay()
        self.counter=0

    def onOpen(self):
        print("WebSocket connection open.")
        self.factory.messagesent=False
        self.factory.connection=self

        def hello():
            self.sendMessage(u"Hello, world!".encode('utf8'))
            self.sendMessage(b"\x00\x01\x03\x04", isBinary=True)
            self.factory.reactor.callLater(1, hello)
            #print "check my bad self"

        def process_shizzle():
            print "sending data to ws"
            if minIONdict == minIONdict_test: ## The dictionary is unchanged since the last cycle
                active = 0
                inactive = 0
                for minION in minIONdict:
                    if minIONdict[minION]["state"]=="active":
                        active += 1
                        if minION not in minIONclassdict:
                            minIONclassdict[minION]=dict()
                        if "connected" not in minIONclassdict[minION]:
                            connectip = "ws://"+ args.ip + ":"+str(minIONdict[minION]["ws_longpoll_port"])+"/"
                            minIONclassdict[minION]["class"]=DummyClient(connectip)
                            connectip2 = "ws://"+ args.ip + ":"+str(minIONdict[minION]["ws_longpoll_port"])+"/user_messages"
                            minIONclassdict[minION]["class2"]=MessagesClient(connectip2)
                        try:
                            if "connected" not in minIONclassdict[minION]:
                                try:
                                    minIONclassdict[minION]["class"].connect()
                                    minIONclassdict[minION]["class2"].connect()
                                    print "GETTING THE GOOD STUFF"
                                    minIONclassdict[minION]["class"].init_minion(minION)
                                    minIONclassdict[minION]["class2"].init_minion(minION)
                                    results = execute_command_as_string(commands('get_analysis_configuration'), ipadd,minIONdict[minION]["port"])
                                    minIONdict[minION]["channelstuff"]=results["result"]["channel_states"]
                                except Exception, err:
                                    print "Connection failed", err
                                minIONclassdict[minION]["connected"]="True"
                            ##Check if connction is new or old:
                            if minION in minIONdict_test:
                                if minIONdict[minION]["state"]==minIONdict_test[minION]["state"]:
                                    try: #To catch mysterious bugs
                                        livedata=dict() # collect a dictionary of useful data - might change this
                                        for query in ('status','dataset','biasvoltageget','bias_voltage_gain','machine_id','machine_name','sample_id','user_error','sequenced_res','yield_res','current_script','disk_space','flow_cell_id','get_seq_metrics'):#,'getstaticdata','get_analysis_configuration'):
                                            results = execute_command_as_string(commands(query), ipadd,minIONdict[minION]["port"])
                                            livedata[query]=results
                                            if query == "disk_space":
                                                #print results
                                                #print results["result"][0]["recommend_stop"]
                                                #print results["result"][0]["recommend_alert"]
                                                check_email_twitter()
                                            #if query == "biasvoltageget":
                                            #    print results
                                            #if query == "get_seq_metrics":
                                            #    print query
                                            print results
                                        minIONdict[minION]["livedata"]=livedata #append results to data stream
                                    except Exception, err:
                                        print "line 1407",err
                                        pass
                                else:
                                    if "scripts" in minIONdict[minION]:
                                        if minIONdict[minION]["scripts"]!=get_run_scripts(minIONdict[minION]["port"]):
                                            minIONdict[minION]["scripts"]=get_run_scripts(minIONdict[minION]["port"])
                                    else:
                                        minIONdict[minION]["scripts"]=get_run_scripts(minIONdict[minION]["port"])
                                        results = execute_command_as_string(commands('startmessagenew'), ipadd,minIONdict[minION]["port"])
                                        try: #To catch mysterious bugs
                                            livedata=dict() # collect a dictionary of useful data - might change this
                                            for query in ('status','dataset','biasvoltageget','bias_voltage_gain','machine_id','machine_name','sample_id','user_error','sequenced_res','sequenced_res','yield_res','current_script','disk_space','flow_cell_id', 'get_seq_metrics'):#,'getstaticdata','get_analysis_configuration'):
                                                #print query
                                                results = execute_command_as_string(commands(query), ipadd,minIONdict[minION]["port"])
                                                livedata[query]=results
                                                if query == "disk_space":
                                                    #print results
                                                    #print results["result"][0]["recommend_stop"]
                                                    #print results["result"][0]["recommend_alert"]
                                                    check_email_twitter()
                                                #if query == "get_seq_metrics":
                                                print query
                                                print results
                                                #if query == "biasvoltageget":
                                                #    print results
                                            minIONdict[minION]["livedata"]=livedata #append results to data stream
                                        except Exception, err:
                                            print "line 1438",err
                                            pass
                        except Exception, err:
                            print "line 1443",err
                            print "Connection Error"
                        else:
                            if "scripts" not in minIONdict[minION]:
                                minIONdict[minION]["scripts"]=get_run_scripts(minIONdict[minION]["port"])
                    else:
                        inactive += 1
            else:
                #print "Dictionary has changed"
                for minION in minIONdict: # if yes then we check each minION to see what is happening
                    if minIONdict[minION]["state"]=="active": #we have a sequencing minION - what is it doing?
                        #print minION, "ACTIVE"
                        try:
                            if minION not in minIONclassdict:
                                connectip = "ws://"+ args.ip + ":"+str(minIONdict[minION]["ws_longpoll_port"])+"/" #Connect to the minIONdict
                                connectip2 = "ws://"+ args.ip + ":"+str(minIONdict[minION]["ws_longpoll_port"])+"/user_messages" #Connect to the minIONdict
                                minIONclassdict[minION]=dict() #Add minION to the dictionary
                            if "connected" not in minIONclassdict[minION]:
                                print "Connecting to this minION",minION
                                minIONclassdict[minION]["class"]=DummyClient(connectip)
                                minIONclassdict[minION]["class2"]=MessagesClient(connectip2)
                            if "connected" not in minIONclassdict[minION]:
                                try:
                                    minIONclassdict[minION]["class"].connect()
                                    minIONclassdict[minION]["class2"].connect()
                                except Exception, err:
                                    print "Connection failed", err
                                try:
                                    print "GETTING THE GOOD STUFF"
                                    minIONclassdict[minION]["class"].init_minion(minION)
                                    minIONclassdict[minION]["class2"].init_minion(minION)
                                    results = execute_command_as_string(commands('get_analysis_configuration'), ipadd,minIONdict[minION]["port"])
                                    minIONdict[minION]["channelstuff"]=results["result"]["channel_states"]
                                except Exception, err:
                                    print "Connection failed", err
                                minIONclassdict[minION]["connected"]="True"
                            livedata=dict() # collect a dictionary of useful data - might change this
                            minIONdict[minION]["scripts"]=get_run_scripts(minIONdict[minION]["port"])
                            for query in (
                            'status', 'dataset', 'biasvoltageget', 'bias_voltage_gain', 'machine_id', 'machine_name',
                            'sample_id', 'user_error', 'sequenced_res', 'sequenced_res', 'yield_res', 'current_script',
                            'disk_space', 'flow_cell_id','get_tracking_id',
                            'get_seq_metrics'):  # ,'getstaticdata','get_analysis_configuration'):
                                #print query
                                results = execute_command_as_string(commands(query), ipadd, minIONdict[minION]["port"])
                                livedata[query] = results
                                if query == "disk_space":
                                    # print results
                                    # print results["result"][0]["recommend_stop"]
                                    # print results["result"][0]["recommend_alert"]
                                    check_email_twitter()
                                if query == "get_seq_metrics":
                                    print query
                                    print results
                            results=execute_command_as_string(commands('get_statistics'), ipadd, minIONdict[minION]["port"])
                            print "Got Stats"
                            meanratio = list()
                            openpore = list()
                            instrand = list()
                            datatofetch = ('seq_pore_level', 'seq_strand_delta','tba_pore_level','tba_event_delta')
                            for item in results['result']['stats_for_channel_name']:
                                #for value in datatofetch:
                                #    print item, value, results['result']['stats_for_channel_name'][item][value]
                                if float(results['result']['stats_for_channel_name'][item]['seq_strand_delta']) != 0:
                                    meanratio.append((float(results['result']['stats_for_channel_name'][item]['seq_pore_level'])/float(results['result']['stats_for_channel_name'][item]['seq_strand_delta'])))
                                    openpore.append(float(results['result']['stats_for_channel_name'][item]['seq_pore_level']))
                                    instrand.append(float(results['result']['stats_for_channel_name'][item]['seq_strand_delta']))
                                    #print item, (float(results['result']['stats_for_channel_name'][item]['seq_pore_level'])/float(results['result']['stats_for_channel_name'][item]['seq_strand_delta']))

                            results2 = execute_command_as_string(commands('get_analysis_configuration'), ipadd,minIONdict[minION]["port"])
                            minIONdict[minION]["channelstuff"]=results2["result"]["channel_states"]
                            minIONdict[minION]["livedata"]=livedata #append results to data stream
                            if "class" in minIONclassdict[minION]:
                                for element in minIONclassdict[minION]["class"].detailsdict:
                                    if "detailsdata" not in minIONdict[minION]:
                                        minIONdict[minION]["detailsdata"]=dict()
                                    if np.mean > 0:
                                        print np.mean(meanratio)
                                        print np.mean(openpore)
                                        print np.mean(instrand)
                                        minIONdict[minION]["detailsdata"]["meanratio"]= np.mean(meanratio)
                                        minIONdict[minION]["detailsdata"]["openpore"] = np.mean(openpore)
                                        minIONdict[minION]["detailsdata"]["instrand"] = np.mean(instrand)
                                    if element == "statistics" and minIONclassdict[minION]["class"].detailsdict[element] != "null":
                                        if element not in minIONdict[minION]["detailsdata"]:
                                            minIONdict[minION]["detailsdata"][element]=dict()
                                        context = dict(list(minIONdict[minION]["detailsdata"][element].items()) + list(minIONclassdict[minION]["class"].detailsdict[element].items()))
                                        minIONdict[minION]["detailsdata"][element]=context
                                    if element == "engine_states" and minIONclassdict[minION]["class"].detailsdict[element] != "null":
                                        for element2 in minIONclassdict[minION]["class"].detailsdict[element]:
                                            if minIONclassdict[minION]["class"].detailsdict[element][element2] != "null":
                                                if element not in minIONdict[minION]["detailsdata"]:
                                                    minIONdict[minION]["detailsdata"][element]=dict()
                                                minIONdict[minION]["detailsdata"][element][element2]=minIONclassdict[minION]["class"].detailsdict[element][element2]
                                    else:
                                        if minIONclassdict[minION]["class"].detailsdict[element] != "null":
                                            minIONdict[minION]["detailsdata"][element]=minIONclassdict[minION]["class"].detailsdict[element]
                        except Exception, err:
                            minIONclassdict.pop(minION)
                            print "line 1411",err
                            print "Connection Error"
                    else: #minION state is inactive
                        #print minION, "INACTIVE"
                        if minION in minIONclassdict:
                            if "connected" in minIONclassdict[minION]:
                                keys = minIONclassdict[minION].keys()
                                for key in keys:
                                    print key
                                    minIONclassdict[minION].pop(key, None)

            sendingminIONdict["DETAILS"][args.dbusername]=minIONdict
            process_tracked_yield()

            #print "tick"
            if self.counter==1:
                #print "Yeah!"
                self.counter=0
            process_channel_information()
            check_databases()
            #if et.connected:
            send_data()
            #else:
            #    print "not connected"
            #    et.messagesent=False
            self.factory.reactor.callLater(5, process_shizzle)

        # start sending messages every 5 seconds..
        process_shizzle()
        #hello()

    def onMessage(self, m, isBinary):
        if isBinary:
            print("Binary message received: {0} bytes".format(len(m)))
        else:
            #print("Text message received: {0}".format(m.decode('utf8')))
            #print "message received!",m
            if str(m) != "Connected - waiting for messages.":
                #print "we're off"
                try:
                    message = json.loads(str(m))
                    if type(message) is dict:
                        ##print "message",m,
                        if message["JOB"] == "testmessage":
                            send_message_port("minoTour is checking communication status with "+str(message["minion"])+".",minIONdict[message["minion"]]["port"])
                        if message["JOB"] =="biasvoltageinc":
                            biasvoltsmod("inc",minIONdict[message["minion"]]["port"])
                        if message["JOB"] =="biasvoltagedec":
                            biasvoltsmod("dec",minIONdict[message["minion"]]["port"])
                        if message["JOB"] == "rename":
                            renamerun(message["NAME"],minIONdict[message["minion"]]["port"])
                        if message["JOB"] == "nameflowcell":
                            renameflowcell(message["NAME"],minIONdict[message["minion"]]["port"])
                        if message["JOB"] == "stopminion":
                            stoprun(minIONdict[message["minion"]]["port"])
                        if message["JOB"] == "startminion":
                            startrun(message["SCRIPT"],minIONdict[message["minion"]]["port"])
                        if message["JOB"] == "custommessage":
                            send_custom_message(minIONdict[message["minion"]]["port"],message["SCRIPT"])
                        if message["JOB"] == "initialiseminion":
                            #print "Trying to initialise minION"
                            #result = helper.initialiseminion()
                            startstop('start',message["minion"])
                            print "Try to run the initialiseminion command now."
                            #execute_command_as_string(commands('initialiseminion'), ipadd,minIONdict[message["minion"]]["port"])
                        if message["JOB"] == "shutdownminion":
                            #First we need to stop the minion if it is running...
                            #stoprun(minIONdict[message["minion"]]["port"])
                            #execute_command_as_string(commands('shutdownminion'), ipadd,minIONdict[message["minion"]]["port"])
                            startstop('stop',message["minion"])
                            print "Try to shut down the minION now."
                except Exception, err:
                    print "line 508",err,str(m)
                    pass

    def onClose(self, wasClean, code, reason):
        print("WebSocket connection closed: {0}".format(reason))
        self.messagesent=False


class MyClientFactory(WebSocketClientFactory, ReconnectingClientFactory):

    protocol = MyClientProtocol

    def clientConnectionFailed(self, connector, reason):
        print("Client connection failed .. retrying ..")
        self.retry(connector)

    def clientConnectionLost(self, connector, reason):
        print("Client connection lost .. retrying ..")
        self.retry(connector)


class DB:
  conn = None

  def connect(self):
    #print "trying to connect to ",args.dbhost
    self.conn = MySQLdb.connect(host=args.dbhost, user=args.dbusername,
                         passwd=args.dbpass, port=args.dbport)
    #print "yay"

  def query(self, sql):
    try:
      #print "trying out ",sql
      cursor = self.conn.cursor()
      #print "Got here!"
      cursor.execute(sql)
      self.conn.commit()
      #print "here?"
    except (AttributeError, MySQLdb.OperationalError):
      self.connect()
      cursor = self.conn.cursor()
      #print "cursor type",type(cursor)
      cursor.execute(sql)
      self.conn.commit()
    #print "Return",type(cursor)
    return cursor

if __name__ == '__main__':
    global minIONdict
    global minIONdict_test
    global minIONclassdict
    global statedict
    global statesummarydict
    global channel_data
    channel_data = dict()
    global chanlookup
    chanlookup ={1:(31,0),2:(31,1),3:(31,2),4:(31,3),5:(31,4),6:(31,5),7:(31,6),8:(31,7),9:(30,0),10:(30,1),11:(30,2),12:(30,3),13:(30,4),14:(30,5),15:(30,6),16:(30,7),17:(29,0),18:(29,1),19:(29,2),20:(29,3),21:(29,4),22:(29,5),23:(29,6),24:(29,7),25:(28,0),26:(28,1),27:(28,2),28:(28,3),29:(28,4),30:(28,5),31:(28,6),32:(28,7),33:(31,15),34:(31,14),35:(31,13),36:(31,12),37:(31,11),38:(31,10),39:(31,9),40:(31,8),41:(30,15),42:(30,14),43:(30,13),44:(30,12),45:(30,11),46:(30,10),47:(30,9),48:(30,8),49:(29,15),50:(29,14),51:(29,13),52:(29,12),53:(29,11),54:(29,10),55:(29,9),56:(29,8),57:(28,15),58:(28,14),59:(28,13),60:(28,12),61:(28,11),62:(28,10),63:(28,9),64:(28,8),65:(3,0),66:(3,1),67:(3,2),68:(3,3),69:(3,4),70:(3,5),71:(3,6),72:(3,7),73:(2,0),74:(2,1),75:(2,2),76:(2,3),77:(2,4),78:(2,5),79:(2,6),80:(2,7),81:(1,0),82:(1,1),83:(1,2),84:(1,3),85:(1,4),86:(1,5),87:(1,6),88:(1,7),89:(0,0),90:(0,1),91:(0,2),92:(0,3),93:(0,4),94:(0,5),95:(0,6),96:(0,7),97:(3,15),98:(3,14),99:(3,13),100:(3,12),101:(3,11),102:(3,10),103:(3,9),104:(3,8),105:(2,15),106:(2,14),107:(2,13),108:(2,12),109:(2,11),110:(2,10),111:(2,9),112:(2,8),113:(1,15),114:(1,14),115:(1,13),116:(1,12),117:(1,11),118:(1,10),119:(1,9),120:(1,8),121:(0,15),122:(0,14),123:(0,13),124:(0,12),125:(0,11),126:(0,10),127:(0,9),128:(0,8),129:(7,0),130:(7,1),131:(7,2),132:(7,3),133:(7,4),134:(7,5),135:(7,6),136:(7,7),137:(6,0),138:(6,1),139:(6,2),140:(6,3),141:(6,4),142:(6,5),143:(6,6),144:(6,7),145:(5,0),146:(5,1),147:(5,2),148:(5,3),149:(5,4),150:(5,5),151:(5,6),152:(5,7),153:(4,0),154:(4,1),155:(4,2),156:(4,3),157:(4,4),158:(4,5),159:(4,6),160:(4,7),161:(7,15),162:(7,14),163:(7,13),164:(7,12),165:(7,11),166:(7,10),167:(7,9),168:(7,8),169:(6,15),170:(6,14),171:(6,13),172:(6,12),173:(6,11),174:(6,10),175:(6,9),176:(6,8),177:(5,15),178:(5,14),179:(5,13),180:(5,12),181:(5,11),182:(5,10),183:(5,9),184:(5,8),185:(4,15),186:(4,14),187:(4,13),188:(4,12),189:(4,11),190:(4,10),191:(4,9),192:(4,8),193:(11,0),194:(11,1),195:(11,2),196:(11,3),197:(11,4),198:(11,5),199:(11,6),200:(11,7),201:(10,0),202:(10,1),203:(10,2),204:(10,3),205:(10,4),206:(10,5),207:(10,6),208:(10,7),209:(9,0),210:(9,1),211:(9,2),212:(9,3),213:(9,4),214:(9,5),215:(9,6),216:(9,7),217:(8,0),218:(8,1),219:(8,2),220:(8,3),221:(8,4),222:(8,5),223:(8,6),224:(8,7),225:(11,15),226:(11,14),227:(11,13),228:(11,12),229:(11,11),230:(11,10),231:(11,9),232:(11,8),233:(10,15),234:(10,14),235:(10,13),236:(10,12),237:(10,11),238:(10,10),239:(10,9),240:(10,8),241:(9,15),242:(9,14),243:(9,13),244:(9,12),245:(9,11),246:(9,10),247:(9,9),248:(9,8),249:(8,15),250:(8,14),251:(8,13),252:(8,12),253:(8,11),254:(8,10),255:(8,9),256:(8,8),257:(15,0),258:(15,1),259:(15,2),260:(15,3),261:(15,4),262:(15,5),263:(15,6),264:(15,7),265:(14,0),266:(14,1),267:(14,2),268:(14,3),269:(14,4),270:(14,5),271:(14,6),272:(14,7),273:(13,0),274:(13,1),275:(13,2),276:(13,3),277:(13,4),278:(13,5),279:(13,6),280:(13,7),281:(12,0),282:(12,1),283:(12,2),284:(12,3),285:(12,4),286:(12,5),287:(12,6),288:(12,7),289:(15,15),290:(15,14),291:(15,13),292:(15,12),293:(15,11),294:(15,10),295:(15,9),296:(15,8),297:(14,15),298:(14,14),299:(14,13),300:(14,12),301:(14,11),302:(14,10),303:(14,9),304:(14,8),305:(13,15),306:(13,14),307:(13,13),308:(13,12),309:(13,11),310:(13,10),311:(13,9),312:(13,8),313:(12,15),314:(12,14),315:(12,13),316:(12,12),317:(12,11),318:(12,10),319:(12,9),320:(12,8),321:(19,0),322:(19,1),323:(19,2),324:(19,3),325:(19,4),326:(19,5),327:(19,6),328:(19,7),329:(18,0),330:(18,1),331:(18,2),332:(18,3),333:(18,4),334:(18,5),335:(18,6),336:(18,7),337:(17,0),338:(17,1),339:(17,2),340:(17,3),341:(17,4),342:(17,5),343:(17,6),344:(17,7),345:(16,0),346:(16,1),347:(16,2),348:(16,3),349:(16,4),350:(16,5),351:(16,6),352:(16,7),353:(19,15),354:(19,14),355:(19,13),356:(19,12),357:(19,11),358:(19,10),359:(19,9),360:(19,8),361:(18,15),362:(18,14),363:(18,13),364:(18,12),365:(18,11),366:(18,10),367:(18,9),368:(18,8),369:(17,15),370:(17,14),371:(17,13),372:(17,12),373:(17,11),374:(17,10),375:(17,9),376:(17,8),377:(16,15),378:(16,14),379:(16,13),380:(16,12),381:(16,11),382:(16,10),383:(16,9),384:(16,8),385:(23,0),386:(23,1),387:(23,2),388:(23,3),389:(23,4),390:(23,5),391:(23,6),392:(23,7),393:(22,0),394:(22,1),395:(22,2),396:(22,3),397:(22,4),398:(22,5),399:(22,6),400:(22,7),401:(21,0),402:(21,1),403:(21,2),404:(21,3),405:(21,4),406:(21,5),407:(21,6),408:(21,7),409:(20,0),410:(20,1),411:(20,2),412:(20,3),413:(20,4),414:(20,5),415:(20,6),416:(20,7),417:(23,15),418:(23,14),419:(23,13),420:(23,12),421:(23,11),422:(23,10),423:(23,9),424:(23,8),425:(22,15),426:(22,14),427:(22,13),428:(22,12),429:(22,11),430:(22,10),431:(22,9),432:(22,8),433:(21,15),434:(21,14),435:(21,13),436:(21,12),437:(21,11),438:(21,10),439:(21,9),440:(21,8),441:(20,15),442:(20,14),443:(20,13),444:(20,12),445:(20,11),446:(20,10),447:(20,9),448:(20,8),449:(27,0),450:(27,1),451:(27,2),452:(27,3),453:(27,4),454:(27,5),455:(27,6),456:(27,7),457:(26,0),458:(26,1),459:(26,2),460:(26,3),461:(26,4),462:(26,5),463:(26,6),464:(26,7),465:(25,0),466:(25,1),467:(25,2),468:(25,3),469:(25,4),470:(25,5),471:(25,6),472:(25,7),473:(24,0),474:(24,1),475:(24,2),476:(24,3),477:(24,4),478:(24,5),479:(24,6),480:(24,7),481:(27,15),482:(27,14),483:(27,13),484:(27,12),485:(27,11),486:(27,10),487:(27,9),488:(27,8),489:(26,15),490:(26,14),491:(26,13),492:(26,12),493:(26,11),494:(26,10),495:(26,9),496:(26,8),497:(25,15),498:(25,14),499:(25,13),500:(25,12),501:(25,11),502:(25,10),503:(25,9),504:(25,8),505:(24,15),506:(24,14),507:(24,13),508:(24,12),509:(24,11),510:(24,10),511:(24,9),512:(24,8)}
    global olddict
    global et
    minIONdict=dict()
    minIONdict_test=dict()
    minIONclassdict=dict()
    statedict=dict()
    statesummarydict=dict()
    minwsip = "ws://"+ args.ip + ":9500/"
    olddict=minIONdict.copy()
    global db
    global globalwarningdict
    globalwarningdict=dict()
    #global cursor

    #trying to implement a reconnecting mySQL connector
    db = DB()


    global helper
    helper = HelpTheMinion(minwsip)
    helper.connect()

    global sendingminIONdict
    sendingminIONdict=dict()
    sendingminIONdict["DETAILS"]=dict()


    import sys

    from twisted.python import log
    from twisted.internet import reactor

    log.startLogging(sys.stdout)

    et = MyClientFactory(u"ws://"+ args.wshost + ":8080")

    reactor.connectTCP(args.wshost, 8080, et)
    try:
        reactor.run()
    except (KeyboardInterrupt,Exception) as err:
        print "ctrl-c detected at top level",err
        print "bye bye"
