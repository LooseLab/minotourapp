"""
A collection of classes used to interrogate MinKNOW
"""

import requests
import json
import datetime
import time
import sys,os
from urllib.parse import urlparse


from ws4py.client.threadedclient import WebSocketClient

from channelmaps import chanlookup

global channel_data
channel_data=dict()



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
        print ("!!!!!!!!!!!!!! Sorry cannot handle linux yet.")
    else:
        print ("!!!!!!!!!!!!!! Sorry - cannot recognise your operating system.")

def send_message_port(message,ipadd,port):
    message_to_send = \
        '{"id":"1", "method":"user_message","params":{"content":"%s"}}' \
        % message
    results=""
    try:
        results = execute_command_as_string(message_to_send, ipadd, port)
    except Exception as err:
        print ("message send fail", err)
    return results

def startrun(script,ipadd,port):
    #print 'Starting the minION device.'
    try:
        startruncustom = \
            '{"id":1, "method":"start_script","params":{"name":"' \
            + script + '"}}'
        #print (startruncustom,ipadd,port)
        startresult = \
            execute_command_as_string(startruncustom,
                ipadd, port)
        startrunmessage = 'minoTour sent a remote run start command.'
        startresultmessage = send_message_port(startrunmessage, ipadd, port)
    except Exception as err:
        print >> sys.stderr, err

def stoprun_message(port,ipadd,message):
    try:
        stopresult = execute_command_as_string(commands('stoprun'),
                ipadd, port)
        stopprotocolresult = \
            execute_command_as_string(commands('stopprotocol'), ipadd,
                port)
        stoprunmessage = 'minoTour sent a remote run stop command for %s.' % (message)
        stopresultmessage = send_message_port(stoprunmessage,port)
    except Exception as err:
        print >> sys.stderr, err

def send_custom_message(port,message):
    try:
        custommessage = "minoTour: %s" % (message)
        custommessageresult=send_message_port(custommessage,port)
    except Exception as err:
        print >> sys.stderr, err

def stoprun(ipadd,port):
    try:
        stopresult = execute_command_as_string(commands('stoprun'),
                ipadd, port)
        stopprotocolresult = \
            execute_command_as_string(commands('stopprotocol'), ipadd,
                port)
        stoprunmessage = 'minoTour sent a remote run stop command.'
        stopresultmessage = send_message_port(stoprunmessage,ipadd,port)
    except Exception as err:
        print >> sys.stderr, err


def renameflowcell(name,ipadd,port):
    try:
        set_flowcell_id = \
            '{"id":1,"method":"set_engine_state","params":{"state_id":"flow_cell_id","value":"'+name+'"}}'
        set_sample_id_result = \
            execute_command_as_string(set_flowcell_id,ipadd,port)
        startresultmessage = \
            send_message_port('minoTour renamed the flowcell to ' + name,ipadd,port)
    except Exception as err:
        print >> sys.stderr, err

def renamerun(name,ipadd,port):
    try:
        set_sample_id = \
            '{"id":"1","method":"set_engine_state","params":{"state_id":"sample_id","value":"' \
            + name + '"}}'
        set_sample_id_result = \
            execute_command_as_string(set_sample_id, ipadd,
                port)
        startresultmessage = \
            send_message_port('minoTour renamed the run to ' + name,ipadd,port)
    except Exception as err:
        print >> sys.stderr, err




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
    def __init__(self, connectip,args):
        super(MessagesClient, self).__init__(connectip)
        #print ("MessagesClient established!")
        self.detailsdict=dict()
        self.daemon=True
        self.args=args

    def init_minion(self,minion):
        self.minion=minion

    def opened(self):
        print ("MessagesClient Opened!")

    def closed(self, code, reason="MessagesClient disconnected for some reason"):
        print ("socket",self.sock.fileno())
        print ("Closed down", code, reason)

    def received_message(self,m):
        if not m.is_binary:
            #print (m)
            try:
                if self.minion not in minIONdict.keys():
                    minIONdict[self.minion]=dict()
                if "messages" not in minIONdict[self.minion].keys():
                    minIONdict[self.minion]["messages"]=[]
                #minIONdict[self.minion]["messages"].append(json.loads(str(m)))
                minIONdict[self.minion]["APIHelp"].update_message(json.loads(str(m)))
                #send_data()
            except:
                pass

def get_run_scripts(ipadd,port):
    #print ("trying to get runscripts on port %s" %(port))
    get_scripts = \
        '{"id":"1", "method":"get_script_info_list","params":null}'
        #'{"id":"1", "method":"get_script_info_list","params":{"state_id":"status"}}'
    results = execute_command_as_string(get_scripts, ipadd, port)
    #print type(results['result'])
    #print len(results['result'])
    #for thing in results['result']:
    #    print thing
    #print (results)
    return results['result']
    #for key in results.keys():
    #    print ("mincontrol:", key, results[key])
    #scriptlist = list()
    #for element in results['result']:
    #    for item in results['result'][element]:
    #        scriptlist.append("('runscript','all','" + item['name']
    #                          + "',1)")
    #print (scriptlist)



def process_tracked_yield(minIONdict):
    """
    This function maintains a history of all yield values over the time of a run - it can be configured to store yields
    at given intervals. It is used by the main web page to provide a history of cumulative yield.
    """
    #print ("Process_tracked_yield Running")
    #print minIONdict
    for minion in minIONdict:
        #print (minion)
        if "detailsdata" in minIONdict[minion]:
            if "livedata" in minIONdict[minion]:
                #print ("in livedata")
                if "yield_res" in minIONdict[minion]["livedata"]:
                    if "result" in minIONdict[minion]["livedata"]["yield_res"]:
                        try:
                            asictemp = float(minIONdict[minion]["detailsdata"]["engine_states"]["minion_asic_temperature"])
                        except:
                            asictemp = 0
                        #print ("asictemp",asictemp)
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
        for i in chanlookup():
            channel_data[minion][i]=dict()
    channel_data[minion][channel]["state"]=state


def process_channel_information(minIONdict,statedict, statesummarydict):
    """
    This function maintains the state of the channels in a simple to read format
    """
    colourlookup={}
    for minion in minIONdict:
        if "detailsdata" in minIONdict[minion]:
            if "channelstuff" in minIONdict[minion]:
                #print ("channelstuff",minIONdict[minion]["channelstuff"])
                for item in minIONdict[minion]["channelstuff"]:
                    if 'style' in minIONdict[minion]["channelstuff"][item]:
                        colourlookup[minIONdict[minion]["channelstuff"][item]['style']['label']]=minIONdict[minion]["channelstuff"][item]['style']['colour']
                if "cust_chan_dict" not in minIONdict[minion]:
                    minIONdict[minion]["cust_chan_dict"]=dict()
                    for i in chanlookup():
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
                statesummarydict[minion]=get_state_summary(minion,channel_data)
                if "simplechanstats" not in minIONdict[minion]:
                    minIONdict[minion]["simplechanstats"]={}
                try:
                    minIONdict[minion]["simplesummary"]=statesummarydict[minion]
                except Exception as err:
                    print ("line 1081",err)
                    pass

def get_state_summary(minion,channel_data):
    state_dict = dict()
    for key, value in channel_data[minion].items():
        try:
            if value["state"] in state_dict:
                state_dict[value["state"]] += 1
            else:
                state_dict[value["state"]] = 1
        except Exception as err:
            print
            "line 956", err
            pass
    #print (state_dict)
    return state_dict


class DummyClient(WebSocketClient):

    def __init__(self, *args,**kwargs):
        super(DummyClient, self).__init__(*args,**kwargs)
        #print ("Client established!")
        self.detailsdict=dict()
        self.daemon=True

    def init_minion(self,minion,minIONdict):
        """ ejecteject has three states:
            0=initisalised or sequencing finished
            1=sequencing happening
            2=run finish detected."""
        self.minion=minion
        self.ejecteject=0
        self.runname=""
        self.minIONdict=minIONdict


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
                        #if  element2 == "channel_info": #element2 == "read_statistics" or
                        #    print ("found ",element2)
                        #    print (json_object[element][element2])
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
                        #print ("engine states",json_object[element]["status"])
                        if json_object[element]["status"]== "finishing":
                            #print ("Looks as though a run is finishing!!!!")
                            self.ejecteject=2
                        if json_object[element]["status"]== "starting":
                            print ("Looks as though a run is starting!!!!")
                            if self.ejecteject==0:
                                #print ("resetting history")
                                if self.minion in self.minIONdict and "yield_history" in self.minIONdict[self.minion]:
                                    self.minIONdict[self.minion]["yield_history"]=[]
                                    self.minIONdict[self.minion]["temp_history"]["asictemp"]=[]
                                    self.minIONdict[self.minion]["temp_history"]["heatsinktemp"]=[]
                                    self.minIONdict[self.minion]["temp_history"]["voltage"]=[]
                                    self.minIONdict[self.minion]["pore_history"]["strand"]=[]
                                    self.minIONdict[self.minion]["pore_history"]["percent"]=[]
                                    self.minIONdict[self.minion]["pore_history"]["single"]=[]
                                    if "details" in self.minIONdict[self.minion]["pore_history"]:
                                        self.minIONdict[self.minion]["pore_history"].pop("details")
                        if json_object[element]["status"]=="processing" and self.ejecteject==0:
                            print ("we're off - resetting history")
                            if self.minion in self.minIONdict and "yield_history" in self.minIONdict[self.minion]:
                                self.minIONdict[self.minion]["yield_history"]=[]
                                self.minIONdict[self.minion]["temp_history"]["asictemp"]=[]
                                self.minIONdict[self.minion]["temp_history"]["heatsinktemp"]=[]
                                self.minIONdict[self.minion]["temp_history"]["voltage"]=[]
                                self.minIONdict[self.minion]["pore_history"]["strand"]=[]
                                self.minIONdict[self.minion]["pore_history"]["percent"]=[]
                                self.minIONdict[self.minion]["pore_history"]["single"]=[]
                                if "details" in self.minIONdict[self.minion]["pore_history"]:
                                    self.minIONdict[self.minion]["pore_history"].pop("details")
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
                #print ("resetting history 773")
                try:
                    if self.minion in self.minIONdict and "yield_history" in self.minIONdict[self.minion]:
                        self.minIONdict[self.minion]["yield_history"]=[]
                        self.minIONdict[self.minion]["temp_history"]["asictemp"]=[]
                        self.minIONdict[self.minion]["temp_history"]["heatsinktemp"]=[]
                        self.minIONdict[self.minion]["temp_history"]["voltage"]=[]
                        self.minIONdict[self.minion]["pore_history"]["strand"]=[]
                        self.minIONdict[self.minion]["pore_history"]["percent"]=[]
                        self.minIONdict[self.minion]["pore_history"]["single"]=[]
                        if "details" in self.minIONdict[self.minion]["pore_history"]:
                            self.minIONdict[self.minion]["pore_history"].pop("details")
                    self.ejecteject=0
                except Exception as err:
                    print ("Problem",err)
                    print ("Probably already popped")


###API shizzle

class MinControlAPI():
    """
    Notes for requests
    header = {'Authorization':'Token e45c142b457121278f9b67d713285a7e10382b36', 'Content-Type':'application/json'}
    r=requests.post(args.full_host+'api/v1/runs/', headers=header, json={"run_name": "20170612_1630_matt", "run_id": "hj78yy9o-217e-4335-9451-66a7288a9dd5", "barcode": "barcode10"})
    r.text #returns result
    """

    def __init__(self, minion, args,statedict,summarystatedict,minIONdict):
        #print ("MinControlAPI")
        self.args = args
        self.minIONdict=minIONdict
        self.minion = minion
        self.create_minion(self.minion)
        self.minionidlink=self.identify_minion(self.minion)
        self.runidlink = ""
        self.status_summary=dict()
        self.minstatexist=False
        self.current_run_id=""
        self.computer=""
        self.runid = 0
        self.flowcelllink = ""
        self.statedict = statedict
        self.summarystatedict = summarystatedict

    def header (self):
        return ({'Authorization': 'Token ' + self.args.api_key, 'Content-Type': 'application/json'})

    def check_jobs (self,minION):
        #print ("checking jobs")
        r = requests.get(self.minionidlink + 'control/', headers=self.header())
        #print (r.text)
        jobs = json.loads(r.text)
        for job in jobs:
            #print (job['job'])
            if job["job"] == "testmessage":
                send_message_port("minoTour is checking communication status with " + str(self.minion) + ".",
                                  self.args.ip,self.minIONdict[self.minion]["port"])
                r = requests.post(self.minionidlink + 'control/' + str(job["id"]) + '/', headers=self.header())
                #print (r.text)
            if job["job"] == "custommessage":
                send_message_port(str(job["custom"]),
                                  self.args.ip,self.minIONdict[self.minion]["port"])
                r = requests.post(self.minionidlink + 'control/' + str(job["id"]) + '/', headers=self.header())
                #print (r.text)
            if job["job"] == "stopminion":
                stoprun(self.args.ip,self.minIONdict[self.minion]["port"])
                r = requests.post(self.minionidlink + 'control/' + str(job["id"]) + '/', headers=self.header())
                #print(r.text)

            if job["job"] == "rename":
                renamerun(job["custom"], self.args.ip,self.minIONdict[self.minion]["port"])
                r = requests.post(self.minionidlink + 'control/' + str(job["id"]) + '/', headers=self.header())
                #print(r.text)
            if job["job"] == "nameflowcell":
                renameflowcell(job["custom"], self.args.ip,self.minIONdict[self.minion]["port"])
                r = requests.post(self.minionidlink + 'control/' + str(job["id"]) + '/', headers=self.header())
                #print(r.text)
            if job["job"] == "stopminion":
                stoprun(self.args.ip,self.minIONdict[self.minion]["port"])
                r = requests.post(self.minionidlink + 'control/' + str(job["id"]) + '/', headers=self.header())
                #print(r.text)
            if job["job"] == "startminion":
                startrun(job["custom"],self.args.ip, self.minIONdict[self.minion]["port"])
                r = requests.post(self.minionidlink + 'control/' + str(job["id"]) + '/', headers=self.header())
                #print(r.text)
            if job["job"] == "initialiseminion":
                # print "Trying to initialise minION"
                # result = helper.initialiseminion()
                startstop('start', self.minion)
                r = requests.post(self.minionidlink + 'control/' + str(job["id"]) + '/', headers=self.header())
                #print(r.text)
                # execute_command_as_string(commands('initialiseminion'), self.args.ip,self.minIONdict[message["minion"]]["port"])
            if job["job"] == "shutdownminion":
                # First we need to stop the minion if it is running...
                # stoprun(self.minIONdict[message["minion"]]["port"])
                # execute_command_as_string(commands('shutdownminion'), self.args.ip,self.minIONdict[message["minion"]]["port"])
                startstop('stop', self.minion)
                r = requests.post(self.minionidlink + 'control/' + str(job["id"]) + '/', headers=self.header())
                #print(r.text)



    def create_minion (self,minION):
        #print("Seen a new run")
        # Test to see if the run exists
        #header = {'Authorization': 'Token ' + self.args.api_key, 'Content-Type': 'application/json'}
        r = requests.get(self.args.full_host + 'api/v1/minions', headers=self.header())
        if minION not in r.text:
            #print ('Need to create minION:', minION)
            createminION = requests.post(self.args.full_host + 'api/v1/minions/', headers=self.header(),
                                      json={"minION_name": minION})
            #print(createminION.text)
            if createminION.status_code != 201:
                print(createminION.status_code)
                print(createminION.text)
                print("Houston - we have a problem!")
            else:
                print(json.loads(createminION.text))

    def create_flowcell(self, name):
        #print ("Creating a new flowcell")
        #Test to see if the flowcell exists.
        r = requests.get(self.args.full_host+'api/v1/flowcells', headers=self.header())
        flowcellname=name
        if flowcellname not in r.text:
            #print ("we need to create this flowcell")
            createflowcell = requests.post(self.args.full_host+'api/v1/flowcells/', headers=self.header(), json={"name": flowcellname})
            self.flowcelllink = json.loads(createflowcell.text)["url"]
        else:
            #print (json.loads(r.text))
            for flowcell in json.loads(r.text):
                #print (flowcell["name"])
                #print (flowcell["url"])
                if flowcell["name"] == flowcellname:
                    self.flowcelllink = flowcell["url"]
                    break

    def create_flowcell_run(self):
        #print ("Adding run to flowcell")
        createflowcellrun = requests.post(self.flowcelllink,headers=self.header(),json={"flowcell": self.flowcelllink, "run": self.runidlink})


    def create_run(self, runid):
        #print ("create run called")
        # This needs to be optimised further - it should search the specific run, not all runs...
        r = requests.get(self.args.full_host+'api/v1/runs', headers=self.header())
        if runid not in r.text:
            is_barcoded = False
            barcoded = "unclassified"
            createrun = requests.post(self.args.full_host+'api/v1/runs/', headers=self.header(), json={"run_name": self.status_summary['run_name'], "run_id": runid, "barcode": barcoded, "is_barcoded":is_barcoded, "minION":self.minionidlink})
            if createrun.status_code != 201:
                print (createrun.status_code)
                print (createrun.text)
                print ("Houston - we have a problem!")
            else:
                #print ("Now we will try and create a runid")
                self.runidlink = json.loads(createrun.text)["url"]
                self.runid = json.loads(createrun.text)["id"]
                #print ("now we will try and create a flowcell:{}".format(self.status_summary['flow_cell_id']))
                self.create_flowcell(self.status_summary['flow_cell_id'])
                self.create_flowcell_run()
        else:
            for run in json.loads(r.text):
                if run["run_id"] == runid:
                    self.runidlink = run["url"]
        self.update_minion_run_stats()

    def update_minion_run_stats (self):
        payload = {"minION": self.minionidlink,
                   #"minKNOW_status": self.status_summary['status'],
                   "minKNOW_current_script": self.status_summary['current_script'],
                   "minKNOW_sample_name": self.status_summary['sample_name'],
                   "minKNOW_exp_script_purpose": self.status_summary['exp_script_purpose'],
                   "minKNOW_flow_cell_id": self.status_summary['flow_cell_id'],
                   "minKNOW_run_name": self.status_summary['run_name'],
                   "run_id": self.runidlink,
                   "minKNOW_version": self.status_summary['minKNOW_version'],
                   "minKNOW_hash_run_id": self.status_summary['hash_run_id'],
                   "minKNOW_script_run_id": self.status_summary['script_run_id'],
                   "minKNOW_real_sample_rate": self.status_summary['real_sample_rate'],
                   "minKNOW_asic_id": self.status_summary['asic_id'],
                   "minKNOW_start_time": self.status_summary['start_time'],
                   "minKNOW_colours_string": self.status_summary['colour_data'],
                   "minKNOW_computer": self.computer,
                   #"minKNOW_histogram_values": self.status_summary['histogram'],
                   #"minKNOW_histogram_bin_width": self.status_summary['bin_width']
                   #"minKNOW_total_drive_space": self.status_summary['bytes_capacity'],
                   #"minKNOW_disk_space_till_shutdown": self.status_summary['bytes_available'],
                   #"minKNOW_warnings": self.status_summary['recommend_alert'],
                   }
        #print (payload)
        #sys.exit()
        createminIONRunStatus = requests.post(self.runidlink + 'rundetails/', headers=self.header(),
                                           json=payload)
        #print("tried to create")
        #print(createminIONRunStatus.status_code)
        #print(createminIONRunStatus.text)

    def update_message(self,messagestring):
        #r = requests.post(self.minionidlink + 'messages/', headers=self.header())
        #print (messagestring)
        payload = {"minION":self.minionidlink,
                   "run_id":self.runidlink,
                   "minKNOW_message":messagestring["message"],
                   "minKNOW_identifier": messagestring["identifier"],
                   "minKNOW_severity": messagestring["severity"],
                   "minKNOW_message_timestamp": messagestring["timestamp"],
        }
        createminIONStatus = requests.post(self.minionidlink + 'messages/', headers=self.header(),
                                           json=payload)
        #print("tried to create")
        #print(createminIONStatus.status_code)
        #print(createminIONStatus.text)


    def identify_minion(self,minION):
        r = requests.get(self.args.full_host + 'api/v1/minions', headers=self.header())
        minionidlink = ""
        for minion in json.loads(r.text):
            if minion["minION_name"] == minION:
                #print(minion)
                minionidlink = minion["url"]
        return(minionidlink)

    def update_minion_current_stats (self,livedata,detailsdata,simplesummary,channelstuff):
        #print ("update minion current stats called")
        #print (livedata)
        #sys.exit()
        #print(livedata["get_seq_metrics"]["result"]["read_event_count_weighted_hist_bin_width"])
        #for key in livedata:
        #    print (key)
        #    print (livedata[key])
        #print (simplesummary)
        try:
            self.status_summary['colour_data']= json.dumps(channelstuff)
            self.status_summary['start_time'] = datetime.datetime.fromtimestamp(
        int(livedata['get_engine_states']['result']['daq_start_time'])
    ).strftime('%Y-%m-%d %H:%M:%S')
            self.status_summary['status'] = livedata['status']['result']
            self.status_summary['current_script'] = livedata['current_script']['result']
            self.status_summary['sample_name']=livedata['get_engine_states']['result']['sample_id']
            self.status_summary['minKNOW_version']=livedata['get_engine_states']['result']['version_string']
            #print ("MinKNOW Version", self.status_summary['minKNOW_version'])
            self.status_summary['flow_cell_id'] = livedata['get_engine_states']['result']['flow_cell_id']
            self.status_summary['run_name'] = livedata['get_engine_states']['result']['data_set']
            self.status_summary['hash_run_id'] = livedata['get_engine_states']['result']['hash_run_id']
            self.status_summary['exp_script_purpose'] = livedata['get_engine_states']['result']['exp_script_purpose']
            self.status_summary['script_run_id'] = livedata['get_engine_states']['result']['script_run_id']
            self.status_summary['real_sample_rate'] = livedata['get_engine_states']['result']['real_sample_rate']
            self.status_summary['asic_id'] = livedata['get_engine_states']['result']['asic_id']
            #print (livedata['disk_space']['result'][0]['bytes_capacity'])
            self.status_summary['bytes_capacity'] = livedata['disk_space']['result'][0]['bytes_capacity']
            self.status_summary['bytes_available'] = livedata['disk_space']['result'][0]['bytes_available']
            self.status_summary['bytes_when_alert_issued'] = livedata['disk_space']['result'][0]['bytes_when_alert_issued']
            self.status_summary['recommend_alert'] = livedata['disk_space']['result'][0]['recommend_alert']
            self.status_summary['histogram']=str(livedata["get_seq_metrics"]["result"]["read_event_count_weighted_hist"])
            self.status_summary['bin_width']=livedata["get_seq_metrics"]["result"]["read_event_count_weighted_hist_bin_width"]
            self.status_summary['read_count'] = livedata["get_seq_metrics"]["result"]["selected_completed_count"]

            #print (self.status_summary)
            #Check if minION exists in database...
            payload = {"minION": self.minionidlink,
                    "minKNOW_status": self.status_summary['status'],
                    "minKNOW_current_script": self.status_summary['current_script'],
                    "minKNOW_sample_name": self.status_summary['sample_name'],
                    "minKNOW_exp_script_purpose": self.status_summary['exp_script_purpose'],
                    "minKNOW_flow_cell_id": self.status_summary['flow_cell_id'],
                    "minKNOW_run_name": self.status_summary['run_name'],
                    "minKNOW_hash_run_id": self.status_summary['hash_run_id'],
                    "minKNOW_script_run_id": self.status_summary['script_run_id'],
                    "minKNOW_real_sample_rate": self.status_summary['real_sample_rate'],
                    "minKNOW_asic_id": self.status_summary['asic_id'],
                    "minKNOW_total_drive_space": self.status_summary['bytes_capacity'],
                    "minKNOW_disk_space_till_shutdown": self.status_summary['bytes_when_alert_issued'],
                    "minKNOW_disk_available": self.status_summary['bytes_available'],
                    "minKNOW_warnings": self.status_summary['recommend_alert'],
                    #"minKNOW_histogram_values": self.status_summary['histogram'],
                    #"minKNOW_histogram_bin_width": self.status_summary['bin_width']
                }
            #print (type(payload))
            if self.minstatexist is False:
                r = requests.get(self.minionidlink + 'status/', headers=self.header())
                #print (r.text)
                #print (r.status_code)
                if r.status_code == 200:
                    #print ('minION status exists')
                    self.minstatexist = True
                else:
                    #print ('need to create a status')
                    createminIONStatus = requests.post(self.minionidlink + 'status/', headers=self.header(),
                                                      json=payload)
                    #print ("tried to create")
                    #print(createminIONStatus.status_code)
            else:
                createminIONStatus = requests.put(self.minionidlink + 'status/', headers=self.header(),
                                             json=payload)
                #print (createminIONStatus.text)
            if self.status_summary['status'] == 'processing':
                if self.current_run_id != self.status_summary['hash_run_id']:
                    #print ("need to create run")
                    self.create_run(self.status_summary['hash_run_id'])
                    self.current_run_id = self.status_summary['hash_run_id']


                #else:
                #    print ('run exists')
                self.update_minion_stats(livedata,detailsdata,simplesummary)
        except Exception as err:
            print ("Problem",err)

    def update_minion_stats (self,livedata,detailsdata,simplesummary):
        #print ("Trying to update minion stats")
        #print ("LIVEDATA")
        #print (livedata["detailsdata"])
        if "yield_res" in livedata:
            if "result" in livedata["yield_res"]:
                #print (detailsdata)
                try:
                    asictemp = float(detailsdata["engine_states"]["minion_asic_temperature"])
                except:
                    asictemp = 0
                #print("asictemp", asictemp)
                try:
                    heatsinktemp = float(
                        detailsdata["engine_states"]["minion_heatsink_temperature"])
                except:
                    heatsinktemp = 0
                if livedata["bias_voltage_gain"]["result"] != "null":
                    # print "!!!!!!!!STUFF!!!!!", self.minIONdict[minion]["livedata"]["bias_voltage_gain"]["result"]
                    biasvoltage = int(livedata["bias_voltage_gain"]["result"])
                else:
                    biasvoltage = 0
                if livedata["biasvoltageget"]["result"] != "null":
                    voltage_val = int(livedata["biasvoltageget"]["result"]["bias_voltage"])
                else:
                    voltage_val = 0
                voltage_value = voltage_val * biasvoltage
                # print "volts:",voltage_value
                try:
                    if livedata["yield_res"]["result"] != "null":
                        yieldval = int(livedata["yield_res"]["result"])
                    else:
                        yieldval = detailsdata["statistics"]["read_event_count"]
                except:
                    yieldval = 0
                try:
                    meanratio = detailsdata["meanratio"]
                    openpore = detailsdata["openpore"]
                    instrand = detailsdata["instrand"]
                except:
                    meanratio = 0
                    openpore = 0
                    instrand = 0

                payload = {"minION": self.minionidlink,
                   "run_id": self.runidlink,
                   "sample_time": str(datetime.datetime.now()),
                   "event_yield": int(livedata['get_seq_metrics']['result']['read_event_count']),
                   "asic_temp": asictemp,
                   "heat_sink_temp": heatsinktemp,
                   "voltage_value": voltage_value,
                   "mean_ratio": meanratio,
                   "open_pore":openpore,
                   "in_strand":instrand,
                   "minKNOW_histogram_values": self.status_summary['histogram'],
                   "minKNOW_histogram_bin_width": self.status_summary['bin_width'],
                    "minKNOW_read_count": self.status_summary['read_count']
                   }

                for category in simplesummary:
                    #print (category)
                    payload[str(category)]=simplesummary[category]
                #print (payload)
                createminionstat = requests.post(self.runidlink + 'runstats/', headers=self.header(), json=payload)
                #print (createminionstat.text)
                #print(createminionstat.status_code)


    def update_minion_status (self,minION,computer,status):
        #print ("Setting the status of {} to {} on {}".format(minION,status,computer))
        #First get id of minION we are updating
        #r = requests.get(self.args.full_host + 'api/v1/minions', headers=self.header())
        #minionidlink=self.identify_minion(minION)
        #for minion in json.loads(r.text):
        #    if minion["minION_name"] == minION:
        #        print(minion)
        #        minionidlink = minion["url"]
        self.computer = computer
        #print (self.minionidlink)
        #Second get id of status we are looking for:
        statusidlink=""
        r = requests.get(self.args.full_host + 'api/v1/events', headers=self.header())
        #print (r.text)
        for info in json.loads(r.text):
            if info["name"] == status:
                #print(info)
                statusidlink= info["url"]
        #print ("status",status, statusidlink)
        #print ({"computer_name": computer, "datetime": str(datetime.datetime.now()), "event": str(urlparse(statusidlink).path),"minION": str(urlparse(self.minionidlink).path)})
        updatestatus = requests.post(self.minionidlink + "events/", headers=self.header(),
                                     json={"computer_name": computer, "datetime": str(datetime.datetime.now()), "event": str(urlparse(statusidlink).path),"minION": str(urlparse(self.minionidlink).path)})
        #print (updatestatus.text)

    def check_scripts(self,minION):
        #minionidlink=self.identify_minion(minION)
        #print ("checking scripts")
        #print (self.minionidlink)
        #print(minionidlink+'scripts')
        r=requests.get(self.minionidlink + 'scripts/', headers = self.header())
        #print (r.text)
        self.scripts=r.text

    def update_script(self,minION,script):
        #print ("Updating script registry for {} with {}".format(minION,script))
        #minionidlink=self.identify_minion(minION)
        #print (script['name'])
        #print (script["tags"])
        #Check if we have seen the script before
        if script["name"] not in self.scripts:
            #so script is not present.
            payload=dict()
            payload["minION"]=self.minionidlink
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
            #print ("PAYLOAD", payload)
            updatestate = requests.post(self.minionidlink + 'scripts/', json=payload, headers=self.header())
        else:
            #so script exists, but is the identifier the same or different?
            #print ("so script exists, but is the identifier the same or different?")
            #print (script["identifier"])
            #print (self.scripts)
            if str(script["identifier"]) not in self.scripts:
                #print ("path has changed")
                for element in json.loads(self.scripts):
                    if script["name"] in str(element):
                        #print (element)
                        payload = dict()
                        payload["minION"] = self.minionidlink
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
                        update = requests.put(self.minionidlink + "scripts/" + str(element["id"]) + "/", json=payload, headers=self.header() )
                        #print (update.text)
            else:
                pass


class HelpTheMinion(WebSocketClient):
    def __init__(self, minswip, args):
        WebSocketClient.__init__(self, minswip)
        self.args = args
        self.minIONdict = dict()
        self.minIONdict_test = dict()
        self.minIONclassdict = dict()
        self.statedict = dict()
        self.summarystatedict = dict()
        self.mcrunning = True

    def opened(self):
        print ("Connected to Master MinION Controller!")

    def initialiseminion():
        result = self.send(json.dumps({"params": "null", "id": 5, "method": "initialization_status"}))
        #print (result)

    def received_message(self, m):
        #print ("message received!")
        #print (m.data)
        for thing in ''.join(map(chr,map(ord,(m.data).decode('windows-1252')))).split('\n'):
            #print (thing)
            if len(thing) > 5 and "2L" not in thing and "2n" not in thing:
                if len(thing) > 5 and (thing[1] == "M" or thing[1] == "G"):
                    if thing[1:8] not in self.minIONdict:
                        #print ("INITIALISING self.minIONdict", thing[1:8])
                        self.minIONdict[thing[1:8]]=dict()
                        self.minIONdict[thing[1:8]]["APIHelp"]=MinControlAPI(thing[1:8],self.args,self.statedict, self.summarystatedict, self.minIONdict)
                        self.minIONdict[thing[1:8]]["APIHelp"].update_minion_status(thing[1:8],'UNKNOWN','connected')
                    minIONports = list(map(lambda x:x-192+8000,filter(lambda x:x>190,map(ord,thing))))
                    if len(minIONports) > 0:
                        self.minIONdict[thing[1:8]]["state"]="active"
                        port = minIONports[0]
                        ws_longpoll_port = minIONports[1]
                        ws_event_sampler_port = minIONports[2]
                        ws_raw_data_sampler_port = minIONports[3]
                        self.minIONdict[thing[1:8]]["port"]=port
                        self.minIONdict[thing[1:8]]["ws_longpoll_port"]=ws_longpoll_port
                        self.minIONdict[thing[1:8]]["ws_event_sampler_port"]=ws_event_sampler_port
                        self.minIONdict[thing[1:8]]["ws_raw_data_sampler_port"]=ws_raw_data_sampler_port
                        results = execute_command_as_string(commands("machine_name"), self.args.ip,
                                                            self.minIONdict[thing[1:8]]["port"])
                        #print (results["result"])
                        self.minIONdict[thing[1:8]]["APIHelp"].update_minion_status(thing[1:8], str(results["result"]), 'active')
                    else:
                        self.minIONdict[thing[1:8]]["state"]="inactive"
                        self.minIONdict[thing[1:8]]["port"]=""
                        self.minIONdict[thing[1:8]]["ws_longpoll_port"]=""
                        self.minIONdict[thing[1:8]]["ws_event_sampler_port"]=""
                        self.minIONdict[thing[1:8]]["ws_raw_data_sampler_port"]=""
                        self.minIONdict[thing[1:8]]["APIHelp"].update_minion_status(thing[1:8], 'UNKNOWN', 'inactive')

    def process_minion(self):
        print ("We're good to go!")
        while self.mcrunning:
            #print (".")
            #print (self.mcrunning)
            #print (minIONdict)
            if self.minIONdict == self.minIONdict_test: ## The dictionary is unchanged since the last cycle
                active = 0
                inactive = 0
                for minION in self.minIONdict:
                    if self.minIONdict[minION]["state"] == "active":
                        active += 1
                        #print ("minION is active")
                        if minION not in self.minIONclassdict:
                            #print ("adding minION to dict")
                            self.minIONclassdict[minION] = dict()
                        if "connected" not in self.minIONclassdict[minION]:
                            connectip = "ws://" + self.args.ip + ":" + str(self.minIONdict[minION]["ws_longpoll_port"]) + "/"
                            self.minIONclassdict[minION]["class"] = DummyClient(connectip)
                            connectip2 = "ws://" + self.args.ip + ":" + str(
                                self.minIONdict[minION]["ws_longpoll_port"]) + "/user_messages"
                            self.minIONclassdict[minION]["class2"] = MessagesClient(connectip2,self.args)
                            #print ("initial handshake")
                        try:
                            if "connected" not in self.minIONclassdict[minION]:
                                try:
                                    self.minIONclassdict[minION]["class"].connect()
                                    self.minIONclassdict[minION]["class2"].connect()
                                    self.minIONclassdict[minION]["class"].init_minion(minION,self.minIONdict)
                                    self.minIONclassdict[minION]["class2"].init_minion(minION)
                                    results = execute_command_as_string(commands('get_analysis_configuration'), self.args.ip,
                                                                        self.minIONdict[minION]["port"])
                                    #print ("getting the good stuff")

                                    self.minIONdict[minION]["channelstuff"] = results["result"]["channel_states"]
                                    #print ("connection good")
                                except Exception as err:
                                    print ("380 Connection failed", err)
                                self.minIONclassdict[minION]["connected"] = "True"
                            ##Check if connction is new or old:
                            if minION in self.minIONdict_test:
                                if self.minIONdict[minION]["state"] == self.minIONdict_test[minION]["state"]:
                                    try:  # To catch mysterious bugs
                                        livedata = dict()  # collect a dictionary of useful data - might change this
                                        for query in (
                                        'status', 'dataset', 'biasvoltageget', 'bias_voltage_gain', 'machine_id',
                                        'machine_name', 'sample_id', 'user_error', 'sequenced_res', 'yield_res',
                                        'current_script', 'disk_space', 'flow_cell_id',
                                        'get_seq_metrics', 'get_engine_states'):  # ,'getstaticdata','get_analysis_configuration'):
                                            results = execute_command_as_string(commands(query), self.args.ip,
                                                                                self.minIONdict[minION]["port"])
                                            livedata[query] = results
                                            #print (query,results)
                                        self.minIONdict[minION]["livedata"] = livedata  # append results to data stream
                                        #print ("data appended")
                                        #sys.exit()
                                    except Exception as err:
                                        print ("line 1407", err)
                                        pass
                                    #print ("OK")
                                else:
                                    if "scripts" in self.minIONdict[minION]:
                                        #print ("running get run scripts line 1110")
                                        if self.minIONdict[minION]["scripts"] != get_run_scripts(self.args.ip,self.minIONdict[minION]["port"]):
                                            #print("running get run scripts 1111")
                                            self.minIONdict[minION]["scripts"] = get_run_scripts(self.args.ip,self.minIONdict[minION]["port"])
                                            #print ("calling check scripts 1113")
                                            self.minIONdict[minION]["APIHelp"].check_scripts(minION)
                                            for script in self.minIONdict[minION]["scripts"]:
                                                # print (script)
                                                self.minIONdict[minION]["APIHelp"].update_script(minION, script)
                                    else:
                                        #print("running get run scripts 1114")
                                        self.minIONdict[minION]["scripts"] = get_run_scripts(self.args.ip,self.minIONdict[minION]["port"])
                                        #print("calling check scripts 1121")
                                        self.minIONdict[minION]["APIHelp"].check_scripts(minION)
                                        for script in self.minIONdict[minION]["scripts"]:
                                            # print (script)
                                            self.minIONdict[minION]["APIHelp"].update_script(minION, script)
                                        results = execute_command_as_string(commands('startmessagenew'), self.args.ip,
                                                                            self.minIONdict[minION]["port"])
                                        try:  # To catch mysterious bugs
                                            livedata = dict()  # collect a dictionary of useful data - might change this
                                            for query in (
                                            'status', 'dataset', 'biasvoltageget', 'bias_voltage_gain', 'machine_id',
                                            'machine_name', 'sample_id', 'user_error', 'sequenced_res', 'sequenced_res',
                                            'yield_res', 'current_script', 'disk_space', 'flow_cell_id',
                                            'get_seq_metrics'):  # ,'getstaticdata','get_analysis_configuration'):
                                                # print query
                                                results = execute_command_as_string(commands(query), self.args.ip,
                                                                                    self.minIONdict[minION]["port"])
                                                livedata[query] = results
                                            self.minIONdict[minION]["livedata"] = livedata  # append results to data stream
                                        except Exception as err:
                                            print ("line 1438", err)
                                            pass
                        except Exception as err:
                            print ("line 1443", err)
                            print ("Connection Error")
                        else:
                            #print ("here")
                            if "scripts" not in self.minIONdict[minION]:
                                #print ("running get run scripts 1147")
                                self.minIONdict[minION]["scripts"] = get_run_scripts(self.args.ip,self.minIONdict[minION]["port"])
                                #print("calling check scripts 1161")
                                self.minIONdict[minION]["APIHelp"].check_scripts(minION)
                                for script in self.minIONdict[minION]["scripts"]:
                                    # print (script)
                                    self.minIONdict[minION]["APIHelp"].update_script(minION, script)
                                #print (get_run_scripts(self.args.ip,self.minIONdict[minION]["port"]))
                                #APIHelp = MinControlAPI()
                                #print("calling check scripts 1168")
                                #self.minIONdict[minION]["APIHelp"].check_scripts(minION)
                                for script in self.minIONdict[minION]["scripts"]:
                                    #print (script)
                                    self.minIONdict[minION]["APIHelp"].update_script(minION,script)
                                #time.sleep(5)
                                #sys.exit()
                    else:
                        inactive += 1
                    process_tracked_yield(self.minIONdict)
                    process_channel_information(self.minIONdict, self.statedict, self.summarystatedict)
                    self.minIONdict[minION]["APIHelp"].update_minion_current_stats(self.minIONdict[minION]['livedata'],self.minIONdict[minION]['detailsdata'],self.minIONdict[minION]["simplesummary"],self.minIONdict[minION]["channelstuff"])
                    self.minIONdict[minION]["APIHelp"].check_jobs(minION)
            else:
                #print ("Dictionary has changed")
                for minION in self.minIONdict:  # if yes then we check each minION to see what is happening
                    if "state" in self.minIONdict[minION].keys():
                        if self.minIONdict[minION]["state"] == "active":  # we have a sequencing minION - what is it doing?
                            #APIHelp =
                            #print (minION, "ACTIVE")
                            try:
                                if minION not in self.minIONclassdict:
                                    connectip = "ws://" + self.args.ip + ":" + str(
                                        self.minIONdict[minION]["ws_longpoll_port"]) + "/"  # Connect to the self.minIONdict
                                    connectip2 = "ws://" + self.args.ip + ":" + str(self.minIONdict[minION][
                                                                                   "ws_longpoll_port"]) + "/user_messages"  # Connect to the self.minIONdict
                                    self.minIONclassdict[minION] = dict()  # Add minION to the dictionary
                                if "connected" not in self.minIONclassdict[minION]:
                                    #print ("Connecting to this minION", minION)
                                    self.minIONclassdict[minION]["class"] = DummyClient(connectip)
                                    self.minIONclassdict[minION]["class2"] = MessagesClient(connectip2,self.args)
                                if "connected" not in self.minIONclassdict[minION]:
                                    try:
                                        self.minIONclassdict[minION]["class"].connect()
                                        self.minIONclassdict[minION]["class2"].connect()
                                    except Exception as err:
                                        print ("480 Connection failed", err)
                                    try:
                                        #print ("GETTING THE GOOD STUFF")
                                        self.minIONclassdict[minION]["class"].init_minion(minION,self.minIONdict)
                                        self.minIONclassdict[minION]["class2"].init_minion(minION)
                                        results = execute_command_as_string(commands('get_analysis_configuration'), self.args.ip,
                                                                            self.minIONdict[minION]["port"])
                                        self.minIONdict[minION]["channelstuff"] = results["result"]["channel_states"]
                                    except Exception as err:
                                        print ("Connection failed", err)
                                    self.minIONclassdict[minION]["connected"] = "True"
                                livedata = dict()  # collect a dictionary of useful data - might change this
                                #print ("running get run scripts 1199")
                                self.minIONdict[minION]["scripts"] = get_run_scripts(self.args.ip,self.minIONdict[minION]["port"])
                                #print("calling check scripts 1216")
                                self.minIONdict[minION]["APIHelp"].check_scripts(minION)
                                for script in self.minIONdict[minION]["scripts"]:
                                    # print (script)
                                    self.minIONdict[minION]["APIHelp"].update_script(minION, script)
                                for query in (
                                        'status', 'dataset', 'biasvoltageget', 'bias_voltage_gain', 'machine_id',
                                        'machine_name',
                                        'sample_id', 'user_error', 'sequenced_res', 'sequenced_res', 'yield_res',
                                        'current_script',
                                        'disk_space', 'flow_cell_id', 'get_tracking_id',
                                        'get_seq_metrics', 'get_engine_states'):  # ,'getstaticdata','get_analysis_configuration'):
                                    results = execute_command_as_string(commands(query), self.args.ip, self.minIONdict[minION]["port"])
                                    livedata[query] = results
                                results = execute_command_as_string(commands('get_statistics'), self.args.ip,
                                                                    self.minIONdict[minION]["port"])
                                #print ("Got Stats")
                                meanratio = list()
                                openpore = list()
                                instrand = list()
                                datatofetch = ('seq_pore_level', 'seq_strand_delta', 'tba_pore_level', 'tba_event_delta')
                                for item in results['result']['stats_for_channel_name']:
                                    if float(results['result']['stats_for_channel_name'][item]['seq_strand_delta']) != 0:
                                        meanratio.append((float(
                                            results['result']['stats_for_channel_name'][item]['seq_pore_level']) / float(
                                            results['result']['stats_for_channel_name'][item]['seq_strand_delta'])))
                                        openpore.append(
                                            float(results['result']['stats_for_channel_name'][item]['seq_pore_level']))
                                        instrand.append(
                                            float(results['result']['stats_for_channel_name'][item]['seq_strand_delta']))

                                results2 = execute_command_as_string(commands('get_analysis_configuration'), self.args.ip,
                                                                     self.minIONdict[minION]["port"])
                                self.minIONdict[minION]["channelstuff"] = results2["result"]["channel_states"]
                                self.minIONdict[minION]["livedata"] = livedata  # append results to data stream
                                if "class" in self.minIONclassdict[minION]:
                                    for element in self.minIONclassdict[minION]["class"].detailsdict:
                                        if "detailsdata" not in self.minIONdict[minION]:
                                            print ("creating detailsdata")
                                            self.minIONdict[minION]["detailsdata"] = dict()
                                        try:
                                            self.minIONdict[minION]["detailsdata"]["meanratio"] = np.mean(meanratio)
                                            self.minIONdict[minION]["detailsdata"]["openpore"] = np.mean(openpore)
                                            self.minIONdict[minION]["detailsdata"]["instrand"] = np.mean(instrand)
                                        except:
                                            self.minIONdict[minION]["detailsdata"]["meanratio"] = 0
                                            self.minIONdict[minION]["detailsdata"]["openpore"] = 0
                                            self.minIONdict[minION]["detailsdata"]["instrand"] = 0
                                        if element == "statistics" and self.minIONclassdict[minION]["class"].detailsdict[
                                            element] != "null":
                                            if element not in self.minIONdict[minION]["detailsdata"]:
                                                self.minIONdict[minION]["detailsdata"][element] = dict()
                                            context = dict(list(self.minIONdict[minION]["detailsdata"][element].items()) + list(
                                                self.minIONclassdict[minION]["class"].detailsdict[element].items()))
                                            self.minIONdict[minION]["detailsdata"][element] = context
                                        if element == "engine_states" and self.minIONclassdict[minION]["class"].detailsdict[
                                            element] != "null":
                                            for element2 in self.minIONclassdict[minION]["class"].detailsdict[element]:
                                                if self.minIONclassdict[minION]["class"].detailsdict[element][
                                                    element2] != "null":
                                                    if element not in self.minIONdict[minION]["detailsdata"]:
                                                        self.minIONdict[minION]["detailsdata"][element] = dict()
                                                    self.minIONdict[minION]["detailsdata"][element][element2] = \
                                                    self.minIONclassdict[minION]["class"].detailsdict[element][element2]
                                        else:
                                            if self.minIONclassdict[minION]["class"].detailsdict[element] != "null":
                                                self.minIONdict[minION]["detailsdata"][element] = \
                                                self.minIONclassdict[minION]["class"].detailsdict[element]
                            except Exception as err:
                                self.minIONclassdict.pop(minION)
                                print ("line 1411", err)
                                print ("Connection Error")
                            process_tracked_yield(self.minIONdict)
                            process_channel_information(self.minIONdict, self.statedict, self.summarystatedict)
                            if "livedata" in self.minIONdict[minION].keys():
                                self.minIONdict[minION]["APIHelp"].update_minion_current_stats(self.minIONdict[minION]['livedata'],
                                                                                      self.minIONdict[minION]['detailsdata'],
                                                                                      self.minIONdict[minION]["simplesummary"],
                                                                                      self.minIONdict[minION]["channelstuff"])
                                self.minIONdict[minION]["APIHelp"].check_jobs(minION)

                        else:  # minION state is inactive
                            # print minION, "INACTIVE"
                            if minION in self.minIONclassdict:
                                if "connected" in self.minIONclassdict[minION]:
                                    keys = self.minIONclassdict[minION].keys()
                                    for key in keys:
                                        #print (key)
                                        self.minIONclassdict[minION].pop(key, None)


            time.sleep(5)

    def hang_up(self):
        print ("ctrl-c detected at top level")
        print ("Disconnecting MinIONs from minoTour control.")
        for minION in self.minIONclassdict:
            self.minIONdict[minION]["APIHelp"].update_minion_status(minION, 'UNKNOWN', 'unplugged')
        print ("bye bye")