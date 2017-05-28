import platform
import configargparse
import os,sys
import json
import requests
import datetime
from urlparse import urlparse

import time
from ws4py.client.threadedclient import WebSocketClient


###API shizzle

class MinControlAPI():
    """
    Notes for requests
    header = {'Authorization':'Token e45c142b457121278f9b67d713285a7e10382b36', 'Content-Type':'application/json'}
    r=requests.post(args.full_host+'api/v1/runs/', headers=header, json={"run_name": "20170612_1630_matt", "run_id": "hj78yy9o-217e-4335-9451-66a7288a9dd5", "barcode": "barcode10"})
    r.text #returns result
    """

    def __init__(self):
        print ("MinControlAPI")
        #self.header=self.header()

    def header (self):
        return ({'Authorization': 'Token ' + args.api_key, 'Content-Type': 'application/json'})


    def create_minion (self,minION):
        print("Seen a new run")
        # Test to see if the run exists
        #header = {'Authorization': 'Token ' + args.api_key, 'Content-Type': 'application/json'}
        r = requests.get(args.full_host + 'api/v1/minions', headers=self.header())
        if minION not in r.text:
            print ('Need to create minION:', minION)
            creatminION = requests.post(args.full_host + 'api/v1/minions/', headers=self.header(),
                                      json={"minION_name": minION})
            print(creatminION.text)
            if creatminION.status_code != 201:
                print(creatminION.status_code)
                print(creatminION.text)
                print("Houston - we have a problem!")
            else:
                print(json.loads(creatminION.text))

    def identify_minion(self,minION):
        r = requests.get(args.full_host + 'api/v1/minions', headers=self.header())
        minionidlink = ""
        for minion in json.loads(r.text):
            if minion["minION_name"] == minION:
                print(minion)
                minionidlink = minion["url"]
        return(minionidlink)

    def update_minion_status (self,minION,computer,status):
        print ("Setting the status of {} to {} on {}".format(minION,status,computer))
        #First get id of minION we are updating
        #r = requests.get(args.full_host + 'api/v1/minions', headers=self.header())
        minionidlink=self.identify_minion(minION)
        #for minion in json.loads(r.text):
        #    if minion["minION_name"] == minION:
        #        print(minion)
        #        minionidlink = minion["url"]
        print (minionidlink)
        #Second get id of status we are looking for:
        statusidlink=""
        r = requests.get(args.full_host + 'api/v1/events', headers=self.header())
        print (r.text)
        for info in json.loads(r.text):
            if info["name"] == status:
                print(info)
                statusidlink= info["url"]
        print ("status",status, statusidlink)
        print ({"computer_name": computer, "datetime": str(datetime.datetime.now()), "event": str(urlparse(statusidlink).path),"minION": str(urlparse(minionidlink).path)})
        updatestatus = requests.post(minionidlink + "events/", headers=self.header(),
                                     json={"computer_name": computer, "datetime": str(datetime.datetime.now()), "event": str(urlparse(statusidlink).path),"minION": str(urlparse(minionidlink).path)})
        print (updatestatus.text)

    def check_scripts(self,minION):
        minionidlink=self.identify_minion(minION)
        print (minionidlink)
        #print(minionidlink+'scripts')
        r=requests.get(minionidlink + 'scripts/', headers = self.header())
        print (r.text)
        self.scripts=r.text

    def update_script(self,minION,script):
        print ("Updating script registry for {} with {}".format(minION,script))
        minionidlink=self.identify_minion(minION)
        print (script['name'])
        print script["tags"]
        #Check if we have seen the script before
        if script["name"] not in self.scripts:
            #so script is not present.
            payload=dict()
            payload["minION"]=minionidlink
            payload["identifier"]= script["identifier"]
            payload["name"]= script["name"]
            if "experiment type" in script["tags"].keys():
                payload["experiment_type"] = script["tags"]["experiment type"]
            if "flow cell" in script["tags"].keys():
                payload["flow_cell"]=script["tags"]["flow cell"]
            if "base calling" in script["tags"].keys():
                payload["base_calling"]=script["tags"]["base calling"]
            if "kit" in script["tags"].keys():
                payload["kit"]=script["tags"]["kit"]
            print ("PAYLOAD", payload)
            updatestate = requests.post(minionidlink + 'scripts/', json=payload, headers=self.header())
        else:
            #so script exists, but is the identifier the same or different?
            if script["identifier"] not in self.scripts:
                print ("path has changed")
                for element in json.loads(self.scripts):
                    if script["name"] in str(element):
                        print (element)
                        payload = dict()
                        payload["minION"] = minionidlink
                        payload["identifier"] = script["identifier"]
                        payload["name"] = script["name"]
                        if "experiment type" in script["tags"].keys():
                            payload["experiment_type"] = script["tags"]["experiment type"]
                        if "flow cell" in script["tags"].keys():
                            payload["flow_cell"] = script["tags"]["flow cell"]
                        if "base calling" in script["tags"].keys():
                            payload["base_calling"] = script["tags"]["base calling"]
                        if "kit" in script["tags"].keys():
                            payload["kit"] = script["tags"]["kit"]
                        print("PAYLOAD", payload)
                        update = requests.put(minionidlink + "scripts/" + str(element["id"]) + "/", json=payload, headers=self.header() )
                        print (update.text)
            else:
                print ("path is the same")






class HelpTheMinion(WebSocketClient):

#    def __init__(self,minIONdict):
#        self.minIONdict = minIONdict

    def opened(self):
        print ("Connected to Master MinION Controller!")

    def initialiseminion():
        result = self.send(json.dumps({"params": "null", "id": 5, "method": "initialization_status"}))
        print (result)


    def received_message(self, m):
        print ("message received!")
        print ("message:", type(m))
        print (m.data)
        #tuple_of_data = struct.unpack("icc", m)
        #print(tuple_of_data)

        for thing in ''.join(map(chr,map(ord,(m.data).decode('windows-1252')))).split('\n'):
            print (thing)
            if len(thing) > 5 and "2L" not in thing and "2n" not in thing:
                if len(thing) > 5 and thing[1] == "M":
                    if thing[1:8] not in minIONdict:
                        print ("INITIALISING MINIONDICT", thing[1:8])
                        minIONdict[thing[1:8]]=dict()
                        APIHelp=MinControlAPI()
                        APIHelp.create_minion(thing[1:8])
                        APIHelp.update_minion_status(thing[1:8],'UNKNOWN','connected')
                    #print ("minION ID:", thing[1:8])
                    minIONports =  list(map(lambda x:x-192+8000,filter(lambda x:x>190,map(ord,thing))))
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
                        APIHelp = MinControlAPI()
                        results = execute_command_as_string(commands("machine_name"), ipadd,
                                                            minIONdict[thing[1:8]]["port"])
                        print (results["result"])
                        APIHelp.update_minion_status(thing[1:8], str(results["result"]), 'active')
                    else:
                        minIONdict[thing[1:8]]["state"]="inactive"
                        minIONdict[thing[1:8]]["port"]=""
                        minIONdict[thing[1:8]]["ws_longpoll_port"]=""
                        minIONdict[thing[1:8]]["ws_event_sampler_port"]=""
                        minIONdict[thing[1:8]]["ws_raw_data_sampler_port"]=""
                        APIHelp = MinControlAPI()
                        APIHelp.update_minion_status(thing[1:8], 'UNKNOWN', 'inactive')

class DummyClient(WebSocketClient):

    def __init__(self, *args,**kwargs):
        super(DummyClient, self).__init__(*args,**kwargs)
        print ("Client established!")
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

"""
    def opened(self):
        self.send(json.dumps({'engine_states':'1','channel_states':'1','channel_info':'1'}))

    def closed(self, code, reason="client disconnected for some reason"):
        print ("socket",self.sock.fileno())
        print ("Closed down", code, reason)

    def received_message(self, m):
        if not m.is_binary:

            #ejecteject=0
            json_object = json.loads(str(m))
            #print json_object
            for element in json_object:
                if element == "statistics" and json_object[element] != "null":
                    for element2 in json_object[element]:
                        if  element2 == "channel_info": #element2 == "read_statistics" or
                            print ("found ",element2)
                            print (json_object[element][element2])
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
                        print ("engine states",json_object[element]["status"])
                        if json_object[element]["status"]== "finishing":
                            print ("Looks as though a run is finishing!!!!")
                            self.ejecteject=2
                        if json_object[element]["status"]== "starting":
                            print ("Looks as though a run is starting!!!!")
                            if self.ejecteject==0:
                                print ("resetting history")
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
                            print ("we're off - resetting history")
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
                print ("resetting history 773")
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
                except Exception as err:
                    print ("Problem",err)
                    print ("Probably already popped")
                    
    """


def execute_command_as_string(data, host=None, port=None):
    """
    Notes for requests
    header = {'Authorization':'Token e45c142b457121278f9b67d713285a7e10382b36', 'Content-Type':'application/json'}
    r=requests.post(args.full_host+'api/v1/runs/', headers=header, json={"run_name": "20170612_1630_matt", "run_id": "hj78yy9o-217e-4335-9451-66a7288a9dd5", "barcode": "barcode10"})
    r.text #returns result
    """
    host_name = host
    port_number = port
    url = 'http://%s:%s%s' % (str(host_name), str(port_number), '/jsonrpc')
    #print (url)
    r=requests.post(url, data=data,headers={'Content-Length': str(len(data)),'Content-Type': 'application/json'})
    try:
        json_respond = json.loads(r.text)
        return json_respond
    except Exception as err:
        print (err)
        err_string = \
            'Fail to initialise mincontrol. Likely reasons include minKNOW not running, the wrong IP address for the minKNOW server or firewall issues.'

    return

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
        'get_engine_states' : '{"id":1,"method": "get_engine_states","params":null}',

    }[command]

class MessagesClient(WebSocketClient):
    def __init__(self, *args,**kwargs):
        super(MessagesClient, self).__init__(*args,**kwargs)
        print ("MessagesClient established!")
        self.detailsdict=dict()
        self.daemon=True

    def init_minion(self,minion):
        self.minion=minion

    def opened(self):
        print ("MessagesClient Opened!")

    def closed(self, code, reason="MessagesClient disconnected for some reason"):
        print ("socket",self.sock.fileno())
        print ("Closed down", code, reason)

    def received_message(self,m):
        if not m.is_binary:
            print (m)
            try:
                if self.minion not in minIONdict.keys():
                    minIONdict[self.minion]=dict()
                if "messages" not in minIONdict[self.minion].keys():
                    minIONdict[self.minion]["messages"]=[]
                minIONdict[self.minion]["messages"].append(json.loads(str(m)))
                #send_data()
            except:
                pass

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
        print ("mincontrol:", key, results[key])
    scriptlist = list()
    for element in results['result']:
        for item in results['result'][element]:
            scriptlist.append("('runscript','all','" + item['name']
                              + "',1)")


if __name__ == '__main__':
    global OPER
    OPER = platform.system()
    if OPER == 'Windows':  # MS
        OPER = 'windows'
    elif OPER == 'Darwin':
        OPER = 'osx'
    else:
        OPER = 'linux'  # MS
    print (OPER)  # MS

    global minIONdict
    minIONdict=dict()

    global minIONdict_test
    minIONdict_test = dict()

    global minIONclassdict
    minIONclassdict=dict()

    config_file = script_dir = os.path.dirname(os.path.realpath('__file__' \
                                                                )) + '/' + 'minup_posix.config'
    parser = \
        configargparse.ArgParser(description= \
                                     'interaction: A program to provide real time interaction for minION runs.' \
                                 , default_config_files=[config_file])

    """
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
    """
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
        '-k',
        '--key',
        type=str,
        required=True,
        default='caf6c3c35db714bacc92e206a75e402b3eec2e05',
        help='The api key for uploading data.',
        dest='api_key',
    )

    parser.add(
        '-p',
        '--port',
        type=int,
        # required=True,
        default=8100,
        help='The port number for the local server.',
        dest='port_number',
    )

    parser.add(
        '-hn',
        '--hostname',
        type=str,
        # required=True,
        default='127.0.0.1',
        help='The run name you wish to provide.',
        dest='host_name',
    )

    global args
    args = parser.parse_args()
    args.full_host = "http://" + args.host_name + ":" + str(args.port_number) + "/"
    version = '0.4'  # 26th May 2017

    ### test which version of python we're using

    ###Machine to connect to address
    global ipadd
    ipadd = args.ip

    minwsip = "ws://" + args.ip + ":9500/"
    global helper
    helper = HelpTheMinion(minwsip)
    helper.connect()

    try:
        while True:
            print (".")
            print (minIONdict)
            if minIONdict == minIONdict_test: ## The dictionary is unchanged since the last cycle
                active = 0
                inactive = 0
                for minION in minIONdict:
                    if minIONdict[minION]["state"] == "active":
                        active += 1
                        print ("minION is active")
                        if minION not in minIONclassdict:
                            print ("adding minION to dict")
                            minIONclassdict[minION] = dict()
                        if "connected" not in minIONclassdict[minION]:
                            connectip = "ws://" + args.ip + ":" + str(minIONdict[minION]["ws_longpoll_port"]) + "/"
                            minIONclassdict[minION]["class"] = DummyClient(connectip)
                            connectip2 = "ws://" + args.ip + ":" + str(
                                minIONdict[minION]["ws_longpoll_port"]) + "/user_messages"
                            minIONclassdict[minION]["class2"] = MessagesClient(connectip2)
                            print ("initial handshake")
                        try:
                            if "connected" not in minIONclassdict[minION]:
                                try:
                                    minIONclassdict[minION]["class"].connect()
                                    minIONclassdict[minION]["class2"].connect()
                                    minIONclassdict[minION]["class"].init_minion(minION)
                                    minIONclassdict[minION]["class2"].init_minion(minION)
                                    results = execute_command_as_string(commands('get_analysis_configuration'), ipadd,
                                                                        minIONdict[minION]["port"])
                                    print ("getting the good stuff")
                                    print(results)
                                    minIONdict[minION]["channelstuff"] = results["result"]["channel_states"]
                                    print ("connection good")
                                except Exception as err:
                                    print ("380 Connection failed", err)
                                minIONclassdict[minION]["connected"] = "True"
                            ##Check if connction is new or old:
                            if minION in minIONdict_test:
                                if minIONdict[minION]["state"] == minIONdict_test[minION]["state"]:
                                    try:  # To catch mysterious bugs
                                        livedata = dict()  # collect a dictionary of useful data - might change this
                                        for query in (
                                        'status', 'dataset', 'biasvoltageget', 'bias_voltage_gain', 'machine_id',
                                        'machine_name', 'sample_id', 'user_error', 'sequenced_res', 'yield_res',
                                        'current_script', 'disk_space', 'flow_cell_id',
                                        'get_seq_metrics', 'get_engine_states'):  # ,'getstaticdata','get_analysis_configuration'):
                                            results = execute_command_as_string(commands(query), ipadd,
                                                                                minIONdict[minION]["port"])
                                            livedata[query] = results
                                            print (query,results)
                                        minIONdict[minION]["livedata"] = livedata  # append results to data stream
                                        print ("data appended")
                                        sys.exit()
                                    except Exception as err:
                                        print ("line 1407", err)
                                        pass
                                    print ("OK")
                                else:
                                    if "scripts" in minIONdict[minION]:
                                        if minIONdict[minION]["scripts"] != get_run_scripts(minIONdict[minION]["port"]):
                                            minIONdict[minION]["scripts"] = get_run_scripts(minIONdict[minION]["port"])
                                    else:
                                        minIONdict[minION]["scripts"] = get_run_scripts(minIONdict[minION]["port"])
                                        results = execute_command_as_string(commands('startmessagenew'), ipadd,
                                                                            minIONdict[minION]["port"])
                                        try:  # To catch mysterious bugs
                                            livedata = dict()  # collect a dictionary of useful data - might change this
                                            for query in (
                                            'status', 'dataset', 'biasvoltageget', 'bias_voltage_gain', 'machine_id',
                                            'machine_name', 'sample_id', 'user_error', 'sequenced_res', 'sequenced_res',
                                            'yield_res', 'current_script', 'disk_space', 'flow_cell_id',
                                            'get_seq_metrics'):  # ,'getstaticdata','get_analysis_configuration'):
                                                # print query
                                                results = execute_command_as_string(commands(query), ipadd,
                                                                                    minIONdict[minION]["port"])
                                                livedata[query] = results
                                                #print (query, results)
                                                if query == "disk_space":
                                                    # print results
                                                    # print results["result"][0]["recommend_stop"]
                                                    # print results["result"][0]["recommend_alert"]
                                                    check_email_twitter()
                                                #if query == "get_seq_metrics":
                                                # if query == "biasvoltageget":
                                                #    print results

                                            minIONdict[minION]["livedata"] = livedata  # append results to data stream
                                        except Exception as err:
                                            print ("line 1438", err)
                                            pass
                        except Exception as err:
                            print ("line 1443", err)
                            print ("Connection Error")
                        else:
                            print ("here")
                            if "scripts" not in minIONdict[minION]:
                                minIONdict[minION]["scripts"] = get_run_scripts(minIONdict[minION]["port"])
                                #print (get_run_scripts(minIONdict[minION]["port"]))
                                APIHelp = MinControlAPI()
                                APIHelp.check_scripts(minION)
                                for script in minIONdict[minION]["scripts"]:
                                    #print (script)
                                    APIHelp.update_script(minION,script)
                                #time.sleep(5)
                                #sys.exit()
                    else:
                        inactive += 1

                print ("at the end")
            else:
                print ("Dictionary has changed")
                for minION in minIONdict:  # if yes then we check each minION to see what is happening
                    if minIONdict[minION]["state"] == "active":  # we have a sequencing minION - what is it doing?
                        # print minION, "ACTIVE"
                        try:
                            if minION not in minIONclassdict:
                                connectip = "ws://" + args.ip + ":" + str(
                                    minIONdict[minION]["ws_longpoll_port"]) + "/"  # Connect to the minIONdict
                                connectip2 = "ws://" + args.ip + ":" + str(minIONdict[minION][
                                                                               "ws_longpoll_port"]) + "/user_messages"  # Connect to the minIONdict
                                minIONclassdict[minION] = dict()  # Add minION to the dictionary
                            if "connected" not in minIONclassdict[minION]:
                                print ("Connecting to this minION", minION)
                                minIONclassdict[minION]["class"] = DummyClient(connectip)
                                minIONclassdict[minION]["class2"] = MessagesClient(connectip2)
                            if "connected" not in minIONclassdict[minION]:
                                try:
                                    minIONclassdict[minION]["class"].connect()
                                    minIONclassdict[minION]["class2"].connect()
                                except Exception as err:
                                    print ("480 Connection failed", err)
                                try:
                                    print ("GETTING THE GOOD STUFF")
                                    minIONclassdict[minION]["class"].init_minion(minION)
                                    minIONclassdict[minION]["class2"].init_minion(minION)
                                    results = execute_command_as_string(commands('get_analysis_configuration'), ipadd,
                                                                        minIONdict[minION]["port"])
                                    minIONdict[minION]["channelstuff"] = results["result"]["channel_states"]
                                except Exception as err:
                                    print ("Connection failed", err)
                                minIONclassdict[minION]["connected"] = "True"
                            livedata = dict()  # collect a dictionary of useful data - might change this
                            minIONdict[minION]["scripts"] = get_run_scripts(minIONdict[minION]["port"])
                            for query in (
                                    'status', 'dataset', 'biasvoltageget', 'bias_voltage_gain', 'machine_id',
                                    'machine_name',
                                    'sample_id', 'user_error', 'sequenced_res', 'sequenced_res', 'yield_res',
                                    'current_script',
                                    'disk_space', 'flow_cell_id', 'get_tracking_id',
                                    'get_seq_metrics'):  # ,'getstaticdata','get_analysis_configuration'):
                                # print query
                                results = execute_command_as_string(commands(query), ipadd, minIONdict[minION]["port"])
                                livedata[query] = results
                                if query == "disk_space":
                                    # print results
                                    # print results["result"][0]["recommend_stop"]
                                    # print results["result"][0]["recommend_alert"]
                                    check_email_twitter()
                                if query == "get_seq_metrics":
                                    print (query)
                                    print (results)
                            results = execute_command_as_string(commands('get_statistics'), ipadd,
                                                                minIONdict[minION]["port"])
                            print ("Got Stats")
                            meanratio = list()
                            openpore = list()
                            instrand = list()
                            datatofetch = ('seq_pore_level', 'seq_strand_delta', 'tba_pore_level', 'tba_event_delta')
                            for item in results['result']['stats_for_channel_name']:
                                # for value in datatofetch:
                                #    print item, value, results['result']['stats_for_channel_name'][item][value]
                                if float(results['result']['stats_for_channel_name'][item]['seq_strand_delta']) != 0:
                                    meanratio.append((float(
                                        results['result']['stats_for_channel_name'][item]['seq_pore_level']) / float(
                                        results['result']['stats_for_channel_name'][item]['seq_strand_delta'])))
                                    openpore.append(
                                        float(results['result']['stats_for_channel_name'][item]['seq_pore_level']))
                                    instrand.append(
                                        float(results['result']['stats_for_channel_name'][item]['seq_strand_delta']))
                                    # print item, (float(results['result']['stats_for_channel_name'][item]['seq_pore_level'])/float(results['result']['stats_for_channel_name'][item]['seq_strand_delta']))

                            results2 = execute_command_as_string(commands('get_analysis_configuration'), ipadd,
                                                                 minIONdict[minION]["port"])
                            minIONdict[minION]["channelstuff"] = results2["result"]["channel_states"]
                            minIONdict[minION]["livedata"] = livedata  # append results to data stream
                            if "class" in minIONclassdict[minION]:
                                for element in minIONclassdict[minION]["class"].detailsdict:
                                    if "detailsdata" not in minIONdict[minION]:
                                        minIONdict[minION]["detailsdata"] = dict()
                                    if np.mean > 0:
                                        print (np.mean(meanratio))
                                        print (np.mean(openpore))
                                        print (np.mean(instrand))
                                        minIONdict[minION]["detailsdata"]["meanratio"] = np.mean(meanratio)
                                        minIONdict[minION]["detailsdata"]["openpore"] = np.mean(openpore)
                                        minIONdict[minION]["detailsdata"]["instrand"] = np.mean(instrand)
                                    if element == "statistics" and minIONclassdict[minION]["class"].detailsdict[
                                        element] != "null":
                                        if element not in minIONdict[minION]["detailsdata"]:
                                            minIONdict[minION]["detailsdata"][element] = dict()
                                        context = dict(list(minIONdict[minION]["detailsdata"][element].items()) + list(
                                            minIONclassdict[minION]["class"].detailsdict[element].items()))
                                        minIONdict[minION]["detailsdata"][element] = context
                                    if element == "engine_states" and minIONclassdict[minION]["class"].detailsdict[
                                        element] != "null":
                                        for element2 in minIONclassdict[minION]["class"].detailsdict[element]:
                                            if minIONclassdict[minION]["class"].detailsdict[element][
                                                element2] != "null":
                                                if element not in minIONdict[minION]["detailsdata"]:
                                                    minIONdict[minION]["detailsdata"][element] = dict()
                                                minIONdict[minION]["detailsdata"][element][element2] = \
                                                minIONclassdict[minION]["class"].detailsdict[element][element2]
                                    else:
                                        if minIONclassdict[minION]["class"].detailsdict[element] != "null":
                                            minIONdict[minION]["detailsdata"][element] = \
                                            minIONclassdict[minION]["class"].detailsdict[element]
                        except Exception as err:
                            minIONclassdict.pop(minION)
                            print ("line 1411", err)
                            print ("Connection Error")
                    else:  # minION state is inactive
                        # print minION, "INACTIVE"
                        if minION in minIONclassdict:
                            if "connected" in minIONclassdict[minION]:
                                keys = minIONclassdict[minION].keys()
                                for key in keys:
                                    print (key)
                                    minIONclassdict[minION].pop(key, None)
            minIONdict_test = minIONdict
            time.sleep(1)


    except (KeyboardInterrupt, Exception) as err:
        print ("ctrl-c detected at top level", err)
        print ("bye bye")