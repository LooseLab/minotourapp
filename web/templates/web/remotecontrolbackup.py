{% extends "web/template_private.html" %}
{% block content %}
    <div class="content-wrapper">
        <!-- Content Header (Page header) -->
        <section class="content-header">
            <h1>
                <span>remoteControl <small>- remote interaction with minKNOW.</small></span>
            </h1>
            <ol class="breadcrumb">
                <li><a href="#"><i class="fa fa-cloud-upload"></i> remoteControl</a></li>
                <li class="active">Here</li>
            </ol>
        </section>
        <!-- Main content -->
        <section class="content">
            <div class="box">
                <div class="box-header">
                    <h3 class="box-title">minIONs Available</h3>
                </div><!-- /.box-header -->
                <div class="box-body">
                    <div class="row">
                        <div class="col-lg-12">
                            <p>This page allows you to remotely control minKNOW. You use it at your own risk.</p>
                            <p>Once connected you will see the minION id below. If you have multiple minIONs connected
                                to one account and all active, they can all be controlled from this page. Further
                                options to implement 'Run Until' are available within the Current Sequencing Run
                                folder.</p>
                            <p>Note that data can take 10-20 seconds to be populated at busy times.</p>
                            <p>The pore history chart is only updated every 5 minutes to reduce page loading times.</p>
                            <div id="app">
                                {% verbatim %}
                                <!-- Nav tabs -->
                                <ul class="nav nav-tabs" role="tablist">
                                    <!--  <li v-for="(key,minion) in minIONdict "><a v-bind:href="'#' + minion" role = "tab" data-toggle="tab">{{ minion }}</a></li>-->
                                    <li v-for="(key,minion) in minions " role="presentation"><a
                                            v-bind:href="'#' + minion.name" role="tab" data-toggle="tab">
                                        <div v-if='minion.status=="processing"'>
                                            Sequencing:{{minion.livedata.machine_id.result}}/{{minion.name}}
                                        </div>
                                        <div v-else>
                                            <div v-if='minion.livedata.machine_id.result!="unknown"'>
                                                On:{{minion.livedata.machine_id.result}}/{{minion.name}}
                                            </div>
                                            <div v-else>Off:{{minion.name}}</div>
                                        </div>
                                    </a></li>
                                </ul>
                                <div class="tab-content">
                                    <div v-for="(key,minion) in minions" role="tabpanel" class="tab-pane" :id="minion.name">
                                        <div v-if='minion.state=="active"'>
                                            <div class="panel panel-info">
                                                <div class="panel-heading">
                                                    <h3 class="panel-title">minKNOW Details -
                                                        {{minion.livedata.dataset.result}}</h3>
                                                </div>
                                                <div class="panel-body">
                                                    <div v-if="minion.livedata.current_script.result.length>0">
                                                        <div class="row">
                                                            <div class="col-md-5">
                                                                <div class="row">
                                                                    <div class="col-md-8"><p><b></i>Experiment
                                                                        Started</i>: {{minion.start_time}}</b></p></div>
                                                                    <div class="col-md-4"><p><i>(Last Update</i>: {{
                                                                        minion.last_update }}</p></div>
                                                                </div>
                                                                <div class="row">
                                                                    <div class="col-md-3"><p><i>MinKNOW version</i>:
                                                                        {{minion.minKNOW_version}}</p></div>
                                                                    <div class="col-md-3"><p><i>Flow Cell ID</i>:
                                                                        {{minion.flow_cell_id}}</p></div>
                                                                    <div class="col-md-3"><p><i>minION ID</i>:
                                                                        {{minion.name}}</p></div>
                                                                    <div class="col-md-3"><p><i>ASIC ID</i>:
                                                                        {{minion.asic_id}}</p></div>
                                                                </div>
                                                                <div class="row">
                                                                    <div class="col-md-8"><p><i>Run Name</i>:
                                                                        {{minion.run_name}}</p></div>
                                                                    <div class="col-md-4"><p><i>Status</i>:
                                                                        {{minion.status}}</p></div>
                                                                </div>
                                                                <div class="row">
                                                                    <div class="col-md-3"><p><i>Yield</i>:
                                                                        {{minion.livedata.event_yield}}</p></div>
                                                                    <!--<div class="col-md-3"><p><i>Channels with Reads</i>: {{minion.statistics.channels_with_read_event_count}}</p></div>-->
                                                                    <!--<div class="col-md-3"><p><i>Read Event Count</i>: {{minion.statistics.read_event_count}}</p></div>-->
                                                                    <div class="col-md-3"><p><i>Completed Read Count</i>:
                                                                        {{minion.read_count}}</p></div>
                                                                </div>
                                                                <div class="row">
                                                                    <div class="col-md-3" :id="minion.name">
                                                                        <div is="container-avg" :title="minion.name"
                                                                             :key="key"
                                                                             :datain="minion.livedata.event_yield"
                                                                             :datain2="minion.read_count"></div>
                                                                    </div>
                                                                    <div class="col-md-3" :id="minion.name">
                                                                        <div is="container-strand" :title="minion.name"
                                                                             :key="key"
                                                                             :datain="minion.currstrand"></div>
                                                                    </div>
                                                                    <div class="col-md-3" :id="minion.name">
                                                                        <div is="container-perc" :title="minion.name"
                                                                             :key="key"
                                                                             :datain="minion.currpercentage"></div>
                                                                    </div>
                                                                </div>
                                                            </div>
                                                            <div class="col-md-4">
                                                                <div class="col-md-12" :id="minion.name">
                                                                    <div is="chartporehist" :title="minion.name"
                                                                         :key="key" :datain="minion.colour_stats"
                                                                         :datain2="minion.pore_history"></div>
                                                                </div>
                                                            </div>
                                                            <div class="col-md-3">
                                                                <h5><b>Messages from MinKNOW:</b></h5>
                                                                <div class="pre-scrollable">
                                                                    <div v-for="message in minion.messages | reverse">
                                                                        <!--<div class="alert alert-{{message.severity}} alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>{{message.message}}<br>{{message.timestamp}}</div>-->
                                                                        <span v-bind:class="'label label-' + message.severity">{{message.severity}}</span>
                                                                        {{message.message}}<br><i>{{message.timestamp}}</i>
                                                                    </div>
                                                                </div>
                                                            </div>
                                                        </div>
                                                        <div class="row">
                                                            <div class="col-md-6">
                                                                <p>The values below are estimated from the run based on
                                                                    a scaling factor from events determined by
                                                                    sequencing speed. The projected yield for 24 hours
                                                                    and 48 hours will update dynamically throuhgout the
                                                                    run.
                                                            </div>
                                                        </div>
                                                        <div class="row">
                                                            <div class="col-md-2">
                                                                <p> Set sequencing speed: </p>
                                                                <select v-model="seqspeed">
                                                                    <option>MegaCrazy Runs</option>
                                                                    <option selected="selected">450 b/s</option>
                                                                    <option>250 b/s</option>
                                                                    <option>70 b/s</option>
                                                                </select>
                                                                <p><span>Selected: {{ seqspeed }}</span></p>
                                                            </div>
                                                            <div is="predictedvals" :seqspeed="seqspeed"
                                                                 :currentyield="minion.livedata.event_yield"
                                                                 :compreads="minion.read_count"
                                                                 :calcurrenttime="minion.start_time"
                                                                 :calcstarttime="minion.start_time"></div>
                                                        </div>
                                                        <div class="row">
                                                        </div>
                                                        <hr>
                                                        <div class="row">
                                                            <div class="col-lg-6" :id="minion.name">
                                                                <div is="chartyield" :title="minion.name" :key="key"
                                                                     :datain2="minion.yield_history"></div>
                                                            </div>
                                                            <div class="col-lg-6" :id="minion.name">
                                                                <div is="porehistory" :title="minion.name" :key="key"
                                                                     :datain2="minion.strand"
                                                                     :datain3="minion.good_single"></div>
                                                            </div>
                                                        </div>
                                                        <div class="row">
                                                            <div class="col-lg-12" :id="minion.name">
                                                                <div is="projchartyield" :start_time="minion.start_time"
                                                                     :seqspeed="seqspeed" :title="minion.name"
                                                                     :key="key" :datain="minion.livedata.event_yield"
                                                                     :datain2="minion.yield_history"></div>
                                                            </div>
                                                        </div>
                                                        <div class="row">
                                                            <div class="col-lg-6" :id="minion.name">
                                                                <div is="perchistory" :title="minion.name" :key="key"
                                                                     :datain2="minion.percentage"></div>
                                                            </div>
                                                            <div class="col-lg-6" :id="minion.name">
                                                                <div is="chartporehistdetails" :title="minion.name"
                                                                     :key="key" :datain="minion.colour_stats"
                                                                     :datain2="minion.pore_history"></div>
                                                            </div>
                                                        </div>
                                                        <div class="row">
                                                            <div class="col-lg-12" :id="minion.name">
                                                                <div is="chartreadhist" :title="minion.name" :key="key"
                                                                     :datain="minion.histogram_values"
                                                                     :datain2="minion.histogram_bin_width"
                                                                     :totalyield="minion.livedata.event_yield"
                                                                     :seqspeed="seqspeed"
                                                                     :readcount="minion.read_count"></div>
                                                            </div>
                                                        </div>
                                                        <div class="row">
                                                            <div class="col-lg-6" :id="minion.name">
                                                                <div is="temphistory" :title="minion.name" :key="key"
                                                                     :datain2="minion.asic_temp"
                                                                     :datain3="minion.heatsinktemp"></div>
                                                            </div>
                                                            <div class="col-lg-6" :id="minion.name">
                                                                <div is="volthistory" :title="minion.name" :key="key"
                                                                     :datain2="minion.voltage"></div>
                                                            </div>
                                                        </div>
                                                        <div class="row">
                                                            <div class="col-lg-6" :id="minion.name">
                                                                <div is="porecurrents" :title="minion.name" :key="key"
                                                                     :datain="minion.instrand_history"
                                                                     :datain2="minion.openpore_history"></div>
                                                            </div>
                                                            <div class="col-lg-6" :id="minion.name">
                                                                <div is="meanratio" :title="minion.name" :key="key"
                                                                     :datain="minion.meanratio_history"></div>
                                                            </div>
                                                        </div>
                                                    </div>
                                                    <div v-else>
                                                        <div class="row">
                                                            <div class="col-md-7">
                                                                <p>This minION is not currently running.</p>
                                                                <button id='inactivateminion'
                                                                        class='btn btn-danger btn-sm'
                                                                        data-toggle='modal'
                                                                        v-bind:data-target="'#' + minion.name + 'offminionmodal'">
                                                                    <i class='fa fa-stop'></i> Switch Off minION
                                                                </button>
                                                            </div>
                                                            <div class="col-md-5">
                                                                <h5><b>Messages from MinKNOW:</b></h5>
                                                                <div class="pre-scrollable">
                                                                    <!--<div v-for="message in minion.messages | reverse" >-->
                                                                    <div v-for="message in minion.messages">
                                                                        <!--<div class="alert alert-{{message.severity}} alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>{{message.message}}<br>{{message.timestamp}}</div>-->
                                                                        <!--<span v-bind:class="'label label-' + message.severity">{{message.severity}}</span>  {{message.message}}<br><i>{{message.timestamp | date "%c"}}</i>-->
                                                                        <span v-bind:class="'label label-' + message.severity">{{message.severity}}</span>
                                                                        {{message.message}}<br><i>{{message.timestamp}}</i>
                                                                    </div>
                                                                </div>
                                                            </div>
                                                        </div>
                                                        <!-- Modal -->
                                                        <div class='modal fade'
                                                             v-bind:id="minion.name + 'offminionmodal'" tabindex='-1'
                                                             role='dialog' aria-labelledby='myModalLabel'
                                                             aria-hidden='true'>
                                                            <div class='modal-dialog'>
                                                                <div class='modal-content'>
                                                                    <div class='modal-header'>
                                                                        <button type='button' class='close'
                                                                                data-dismiss='modal' aria-hidden='true'>
                                                                            &times;
                                                                        </button>
                                                                        <h4 class='modal-title' id='myModalLabel'>Stop
                                                                            your minION</h4>
                                                                    </div>
                                                                    <div class='modal-body'>
                                                                        <div v-bind:id="minion.name + 'offminioninfo'"
                                                                        <p>This action will switch the minION to an
                                                                            inactive state. It should be possible to
                                                                            reactivate the minION remotely but software
                                                                            crashes on the minION controlling device may
                                                                            cause problems. You should only inactivate
                                                                            your minION device remotely if you are
                                                                            certain you wish to do so and <strong> at
                                                                                your own risk</strong>.</p>

                                                                        <p>If you are sure you wish to do this, click
                                                                            'Inactivate minION' below. Otherwise close
                                                                            this window.</p>
                                                                    </div>
                                                                    <div class='modal-footer'>
                                                                        <button type='button' class='btn btn-default'
                                                                                data-dismiss='modal'>Close
                                                                        </button>
                                                                        <button v-on:click="inactivateminion"
                                                                                :id='minion.name' type='button'
                                                                                class='btn btn-danger'
                                                                                data-dismiss='modal'>Switch Off minION
                                                                        </button>
                                                                    </div>
                                                                </div><!-- /.modal-content -->
                                                            </div><!-- /.modal-dialog -->
                                                        </div><!-- /.modal -->
                                                    </div>
                                                </div>

                                            </div>
                                        </div>

                                        <div class="col-md-4">
                                            <div class="panel panel-warning">
                                                <div class="panel-heading">
                                                    <h3 class="panel-title">minKNOW Control options</h3>
                                                </div>
                                                <div class="panel-body">
                                                    <h5>To test if you have a connection to minKNOW:</h5>
                                                    <button v-on:click="testmessage" :id='minion.url' type='button'
                                                            class='btn btn-info btn-sm'><i class='fa fa-magic'></i> Test
                                                        Communication
                                                    </button>
                                                    <br>


                                                    <h5>Log Custom Message:</h5>
                                                    <!--<button id='renamerun' type='button' class='btn btn-info btn-sm'><i class='fa fa-magic'></i> Rename Run</button>-->
                                                    <!-- Indicates a dangerous or potentially negative action -->
                                                    <!-- Button trigger modal -->
                                                    <button id='custommessage' class='btn btn-info btn-sm'
                                                            data-toggle='modal'
                                                            v-bind:data-target="'#' + minion.name + 'custommessage'">
                                                        <i class='fa fa-magic'></i> Send Message
                                                    </button>
                                                    <!-- Modal -->
                                                    <div class='modal fade' v-bind:id="minion.name+'custommessage'"
                                                         tabindex='-1' role='dialog' aria-labelledby='myModalLabel'
                                                         aria-hidden='true'>
                                                        <div class='modal-dialog'>
                                                            <div class='modal-content'>
                                                                <div class='modal-header'>
                                                                    <button type='button' class='close'
                                                                            data-dismiss='modal' aria-hidden='true'>
                                                                        &times;
                                                                    </button>
                                                                    <h4 class='modal-title' id='myModalLabel'>Send a
                                                                        custom message to MinKNOW</h4>
                                                                </div>
                                                                <div class='modal-body'>
                                                                    <div v-bind:id="minion.name + 'custommessage'">
                                                                        <p>You can log a custom message into minKNOW to
                                                                            document things such as addition of library
                                                                            or other event. This will then be written
                                                                            into the log files for this run.</p>
                                                                        <input type="text"
                                                                               v-bind:id="minion.name+'custommessagefield'"
                                                                               class="form-control" placeholder="">
                                                                    </div>
                                                                    <div class='modal-footer'>
                                                                        <button type='button' class='btn btn-default'
                                                                                data-dismiss='modal'>Close
                                                                        </button>
                                                                        <button v-on:click="custommessage"
                                                                                :id='minion.url' :name='minion.name'
                                                                                type='button' class='btn btn-danger'
                                                                                data-dismiss='modal'>Send Message
                                                                        </button>
                                                                    </div>
                                                                </div><!-- /.modal-content -->
                                                            </div><!-- /.modal-dialog -->
                                                        </div><!-- /.modal -->
                                                    </div>


                                                    <h5>Rename Your Run:</h5>
                                                    <!--<button id='renamerun' type='button' class='btn btn-info btn-sm'><i class='fa fa-magic'></i> Rename Run</button>-->
                                                    <!-- Indicates a dangerous or potentially negative action -->
                                                    <!-- Button trigger modal -->
                                                    <button id='renamerun' class='btn btn-info btn-sm'
                                                            data-toggle='modal'
                                                            v-bind:data-target="'#' + minion.name + 'renamemodal'">
                                                        <i class='fa fa-magic'></i> Rename Run
                                                    </button>
                                                    <!-- Modal -->
                                                    <div class='modal fade' v-bind:id="minion.name+'renamemodal'"
                                                         tabindex='-1' role='dialog' aria-labelledby='myModalLabel'
                                                         aria-hidden='true'>
                                                        <div class='modal-dialog'>
                                                            <div class='modal-content'>
                                                                <div class='modal-header'>
                                                                    <button type='button' class='close'
                                                                            data-dismiss='modal' aria-hidden='true'>
                                                                        &times;
                                                                    </button>
                                                                    <h4 class='modal-title' id='myModalLabel'>Rename
                                                                        Your Run</h4>
                                                                </div>
                                                                <div class='modal-body'>
                                                                    <div v-bind:id="minion.url+'renameinfo'">
                                                                        <p>You can rename a run if you wish to do so.
                                                                            Note there is no need to do this unless you
                                                                            wish to change the sample ID for some
                                                                            reason.</p>
                                                                        <input type="text"
                                                                               v-bind:id="minion.name+'newname'"
                                                                               class="form-control"
                                                                               placeholder="New Run Name">
                                                                        <p>If you are sure you wish to do this enter
                                                                            your new name above and click 'Rename Run'
                                                                            below. Otherwise close this window.</p>
                                                                        <p> We dont recommend doing this when a run is
                                                                            in progress!</p>
                                                                    </div>
                                                                    <div class='modal-footer'>
                                                                        <button type='button' class='btn btn-default'
                                                                                data-dismiss='modal'>Close
                                                                        </button>
                                                                        <button v-on:click="renamenow" :id='minion.url'
                                                                                :name='minion.name' type='button'
                                                                                class='btn btn-danger'
                                                                                data-dismiss='modal'>Rename Run
                                                                        </button>
                                                                    </div>
                                                                </div><!-- /.modal-content -->
                                                            </div><!-- /.modal-dialog -->
                                                        </div><!-- /.modal -->
                                                    </div>
                                                    <h5>Rename Flow Cell:</h5>
                                                    <!--<button id='renamerun' type='button' class='btn btn-info btn-sm'><i class='fa fa-magic'></i> Rename Run</button>-->
                                                    <!-- Indicates a dangerous or potentially negative action -->
                                                    <!-- Button trigger modal -->
                                                    <button id='renameflowcell' class='btn btn-info btn-sm'
                                                            data-toggle='modal'
                                                            v-bind:data-target="'#' + minion.name + 'renameflowcell'">
                                                        <i class='fa fa-magic'></i> Rename Flowcell
                                                    </button>
                                                    <!-- Modal -->
                                                    <div class='modal fade' v-bind:id="minion.name+'renameflowcell'"
                                                         tabindex='-1' role='dialog' aria-labelledby='myModalLabel'
                                                         aria-hidden='true'>
                                                        <div class='modal-dialog'>
                                                            <div class='modal-content'>
                                                                <div class='modal-header'>
                                                                    <button type='button' class='close'
                                                                            data-dismiss='modal' aria-hidden='true'>
                                                                        &times;
                                                                    </button>
                                                                    <h4 class='modal-title' id='myModalLabel'>Rename
                                                                        Your Flowcell</h4>
                                                                </div>
                                                                <div class='modal-body'>
                                                                    <div v-bind:id="minion.name+'renameflowcellinfo'">
                                                                        <p>You can rename a flowcell if you wish to do
                                                                            so. Note there is no need to do this unless
                                                                            you wish to change the Flowcell ID for some
                                                                            reason.</p>
                                                                        <input type="text"
                                                                               v-bind:id="minion.name+'newflowcellname'"
                                                                               class="form-control"
                                                                               placeholder="New Flowcell ID">
                                                                        <p>If you are sure you wish to do this enter
                                                                            your new Flowcell ID above and click 'Rename
                                                                            Run' below. Otherwise close this window.</p>
                                                                        <p> We dont recommend doing this when a run is
                                                                            in progress!</p>
                                                                    </div>
                                                                    <div class='modal-footer'>
                                                                        <button type='button' class='btn btn-default'
                                                                                data-dismiss='modal'>Close
                                                                        </button>
                                                                        <button v-on:click="renameflowcellnow"
                                                                                :id='minion.url' :name='minion.name'
                                                                                type='button' class='btn btn-danger'
                                                                                data-dismiss='modal'>Rename Flowcell
                                                                        </button>
                                                                    </div>
                                                                </div><!-- /.modal-content -->
                                                            </div><!-- /.modal-dialog -->
                                                        </div><!-- /.modal -->
                                                    </div>


                                                    <br>
                                                    <h5>Remote Start/Stop Sequencing:</h5>

                                                    <!-- Indicates a dangerous or potentially negative action -->
                                                    <!-- Button trigger modal -->
                                                    <button id='stopminion' class='btn btn-danger btn-sm'
                                                            data-toggle='modal'
                                                            v-bind:data-target="'#' + minion.name + 'stopminionmodal'">
                                                        <i class='fa fa-stop'></i> Stop minION
                                                    </button>

                                                    <!-- Modal -->
                                                    <div class='modal fade' v-bind:id="minion.name+'stopminionmodal'"
                                                         tabindex='-1' role='dialog' aria-labelledby='myModalLabel'
                                                         aria-hidden='true'>
                                                        <div class='modal-dialog'>
                                                            <div class='modal-content'>
                                                                <div class='modal-header'>
                                                                    <button type='button' class='close'
                                                                            data-dismiss='modal' aria-hidden='true'>
                                                                        &times;
                                                                    </button>
                                                                    <h4 class='modal-title' id='myModalLabel'>Stop your
                                                                        minION</h4>
                                                                </div>
                                                                <div class='modal-body'>
                                                                    <div v-bind:id="minion.name+'stopminioninfo'">
                                                                        <p>This will attempt to stop your minION
                                                                            sequencer remotely. It should be possible to
                                                                            restart sequencing remotely but software
                                                                            crashes on the minION controlling device may
                                                                            cause problems. You should only stop your
                                                                            minION device remotely if you are certain
                                                                            you wish to do so and <strong> at your own
                                                                                risk</strong>.</p>

                                                                        <p>If you are sure you wish to do this, click
                                                                            'Stop minION' below. Otherwise close this
                                                                            window.</p>
                                                                    </div>
                                                                    <div class='modal-footer'>
                                                                        <button type='button' class='btn btn-default'
                                                                                data-dismiss='modal'>Close
                                                                        </button>
                                                                        <button v-on:click="stopminion" :id='minion.url'
                                                                                type='button' class='btn btn-danger'
                                                                                data-dismiss='modal'>Stop minION
                                                                        </button>
                                                                    </div>
                                                                </div><!-- /.modal-content -->
                                                            </div><!-- /.modal-dialog -->
                                                        </div><!-- /.modal -->
                                                    </div>


                                                    <!-- Indicates a dangerous or potentially negative action -->
                                                    <!-- Button trigger modal -->
                                                    <button id='startminion' class='btn btn-success btn-sm'
                                                            data-toggle='modal'
                                                            v-bind:data-target="'#' + minion.name + 'startminionmodal'">
                                                        <i class='fa fa-play'></i> Start minION
                                                    </button>

                                                    <!-- Modal -->
                                                    <div class='modal fade' v-bind:id="minion.name+'startminionmodal'"
                                                         tabindex='-1' role='dialog' aria-labelledby='myModalLabel'
                                                         aria-hidden='true'>
                                                        <div class='modal-dialog'>
                                                            <div class='modal-content'>
                                                                <div class='modal-header'>
                                                                    <button type='button' class='close'
                                                                            data-dismiss='modal' aria-hidden='true'>
                                                                        &times;
                                                                    </button>
                                                                    <h4 class='modal-title' id='myModalLabel'>Start your
                                                                        minION</h4>
                                                                </div>
                                                                <div class='modal-body'>
                                                                    <div v-bind:id="minion.name+'startminioninfo'">
                                                                        <p>This will attempt to restart your minION
                                                                            sequencer remotely.</p>

                                                                        <p>If you are sure you wish to do this select an
                                                                            available run script and click 'Start
                                                                            minION' below. Otherwise close this
                                                                            window.</p>
                                                                        <div v-for="script in minion.scripts"
                                                                             class='radio'>
                                                                            <label>
                                                                                <input type='radio' name='scriptRadios'
                                                                                       :id='script.name'
                                                                                       :value='script.identifier'>{{script.name}}.py</label>
                                                                        </div>
                                                                    </div>
                                                                    <div class='modal-footer'>
                                                                        <button type='button'
                                                                                class='btn btn-default btn-sm'
                                                                                data-dismiss='modal'>Close
                                                                        </button>
                                                                        <button v-on:click="startminion"
                                                                                :id='minion.url' type='button'
                                                                                class='btn btn-success'
                                                                                data-dismiss='modal'>Start minION
                                                                        </button>
                                                                    </div>
                                                                </div><!-- /.modal-content -->
                                                            </div><!-- /.modal-dialog -->
                                                        </div><!-- /.modal -->
                                                    </div>
                                                    <br><br>
                                                    <button class="btn btn-primary" type="button" data-toggle="collapse"
                                                            v-bind:data-target="'#collapseExample' + minion.name"
                                                            aria-expanded="false" aria-controls="collapseExample">
                                                        Available Scripts:
                                                    </button>
                                                    <div class="collapse" v-bind:id="'collapseExample' + minion.name">
                                                        <div class="well">
                                                            <div v-for="script in minion.scripts">{{script.name}}</div>
                                                        </div>
                                                    </div>


                                                </div>
                                            </div>
                                        </div>

                                        <div class="col-md-4">
                                            <div class='panel panel-info'>
                                                <div class='panel-heading'>
                                                    <h3 class='panel-title'>minKNOW real time data</h3>
                                                </div>
                                                <div class='panel-body'>
                                                    <div class='table-responsive'>
                                                        <table class='table table-condensed'>
                                                            <tr>
                                                                <th>Category
                                                                </td>
                                                                <th>Info
                                                                </td>
                                                            </tr>
                                                            <tr>
                                                                <td>minKNOW computer name</td>
                                                                <td>{{minion.livedata.machine_id.result}}</td>
                                                            </tr>
                                                            <tr>
                                                                <td>minKNOW Status</td>
                                                                <td>{{minion.status}}</td>
                                                            </tr>
                                                            <tr>
                                                                <td>Current Script</td>
                                                                <td>{{minion.livedata.current_script.result}}</td>
                                                            </tr>
                                                            <tr>
                                                                <td>Sample Name</td>
                                                                <td>{{minion.sample_id}}</td>
                                                            </tr>
                                                            <tr>
                                                                <td>Flow Cell ID</td>
                                                                <td>{{minion.flow_cell_id}}</td>
                                                            </tr>
                                                            <tr>
                                                                <td>Run Name</td>
                                                                <td>{{minion.run_name}}</td>
                                                            </tr>
                                                            <tr>
                                                                <td>Voltage Offset</td>
                                                                <td>{{minion.livedata.voltage_value}} mV</td>
                                                            </tr>
                                                            <tr>
                                                                <td>Yield</td>
                                                                <td>{{minion.livedata.event_yield}}</td>
                                                            </tr>
                                                        </table>
                                                    </div>
                                                </div>

                                            </div>


                                        </div>

                                        <div class="col-lg-4" :id="minion.name">
                                            <div is="diskusage" :title="minion.livedata.machine_id.result" :key="key"
                                                 :datain="minion.livedata.disk_space.result"></div>
                                        </div>


                                        <div class="col-md-12">
                                            <br>

                                        </div>
                                    </div>

                                   <!-- <div v-else>-->
                                    <div>
                                        <h5>You are interacting with minION: {{minion.name}}</h5>
                                        It is currently <i>inactive</i>.

                                        <!-- Indicates a dangerous or potentially negative action -->
                                        <!-- Button trigger modal -->
                                        <button id='initminion' class='btn btn-warning btn-sm' data-toggle='modal'
                                                v-bind:data-target="'#' + minion.name + 'initminionmodal'">
                                            <i class='fa fa-stop'></i> Initialise minION
                                        </button>

                                        <!-- Modal -->
                                        <div class='modal fade' v-bind:id="minion.name+'initminionmodal'" tabindex='-1'
                                             role='dialog' aria-labelledby='myModalLabel' aria-hidden='true'>
                                            <div class='modal-dialog'>
                                                <div class='modal-content'>
                                                    <div class='modal-header'>
                                                        <button type='button' class='close' data-dismiss='modal'
                                                                aria-hidden='true'>&times;
                                                        </button>
                                                        <h4 class='modal-title' id='myModalLabel'>Start your minION</h4>
                                                    </div>
                                                    <div class='modal-body'>
                                                        <div v-bind:id="minion.name+'initminioninfo'">
                                                            <p>This action will switch the minION to the active
                                                                state.</p>

                                                            <p>If you are sure you wish to do this, click 'Initialise
                                                                minION' below. Otherwise close this window.</p>
                                                        </div>
                                                        <div class='modal-footer'>
                                                            <button type='button' class='btn btn-default'
                                                                    data-dismiss='modal'>Close
                                                            </button>
                                                            <button v-on:click="initminion" :id='minion.name'
                                                                    type='button' class='btn btn-warning'
                                                                    data-dismiss='modal'>Initialise minION
                                                            </button>
                                                        </div>
                                                    </div><!-- /.modal-content -->
                                                </div><!-- /.modal-dialog -->
                                            </div><!-- /.modal -->
                                        </div>


                                    </div>
                                </div>

                            </div>


                            <!-- Tab panes-->
                            {% endverbatim %}

                        </div>
                    </div>
                </div>
            </div>


        </section><!-- /.content -->
    </div><!-- /.content-wrapper -->

    <script>

        $(document).ready(function () {
            function add(a, b) {
                return a + b;
            }

            function numberWithCommas(x) {
                return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
            }

            function geteighthours(data, runstart) {
                return (28800 * 1000 - (data[0][0] - (runstart * 1000)));
            }

            function projectresults(syntheticdata, scalingfactor, steps, difference, runstart) {
                testdata = syntheticdata.slice(-2);
                var lastval = testdata[1][1];
                var lasttime = testdata[1][0];
                var timingdiff = testdata[1][0] - testdata[0][0];
                //var valdiff = testdata[1][1]-testdata[0][1];
                var valdiff = difference;
                var newresults = [];
                var muxphase = 0
                //console.log(lastval,timingdiff);
                for (var i = 0; i < steps; i++) {
                    templastval = lastval;
                    lastval = lastval + (valdiff * scalingfactor);
                    valdiff = lastval - templastval;
                    lasttime = lasttime + timingdiff;
                    if (muxphase != Math.floor((lasttime - (runstart * 1000)) / 1000 / 28800)) {
                        difference = difference * 0.9;
                        valdiff = difference;
                        muxphase = Math.floor((lasttime - (runstart * 1000)) / 1000 / 28800);
                    }
                    //remainder = lasttime-(runstart*1000) - (lasttime-(runstart*1000) % (28800*1000))/(28800*1000);
                    //console.log(Math.floor((lasttime-(runstart*1000))/1000/28800));
                    newresults.push([lasttime, Math.ceil(lastval)]);
                }
                //console.log(newresults);
                return newresults;
            }

            function round(value, decimals) {
                return Number(Math.round(value + 'e' + decimals) + 'e-' + decimals);
            }

            function converttobases(data, seqspeed) {
                switch (seqspeed) {
                    case "MegaCrazy Runs":
                        scaling = 3.5;
                        break;
                    case "450 b/s":
                        scaling = 1.8;
                        break;
                    case "250 b/s":
                        scaling = 1.1;
                        break;
                    case "70 b/s":
                        scaling = 1.0;
                        break;
                }
                var scaleddata = [];
                for (var i = 0; i < data.length; i++) {
                    scaleddata.push([data[i][0], data[i][1] * scaling]);
                }
                return scaleddata;
            }

            function projectdata(data) {
                var results = [];
                var holder = [];
                var diffholder = 0;
                var meanholder = 0;
                if (data.length > 3000) {
                    data = data.slice(-3000);
                }


                for (var i = 1; i < data.length; i++) {
                    var diff = data[i][1] - data[i - 1][1];
                    holder.push(diff);
                    meanholder = meanholder + diff;
                }
                for (var i = 2; i < holder.length; i++) {
                    var ratio = holder[i] / holder[i - 1];
                    //if (ratio > 2){
                    //    ratio = 2;
                    //}
                    diffholder = diffholder + ratio;
                    //if (meanholder > 1) {
                    //    meanholder = 1;
                    //}
                }
                //console.log(diffholder/(holder.length - 1));
                if (diffholder / (holder.length - 1) > 1) {
                    return [1, meanholder / (holder.length - 1 )];
                } else if (diffholder / (holder.length - 1) < 0.999) {
                    return ([0.999, meanholder / (holder.length - 1 )]);
                } else {
                    return ([diffholder / (holder.length - 1), meanholder / (holder.length - 1 )]);
                }

            }

            function scaleyield(firstelement, data) {
                var results = [];
                for (var i = 0; i < data.length; i++) {
                    //console.log(data[i]);
                    //console.log((data[i][0]-firstelement[0])/1000);
                    //console.log((data[i][1]/1000000));
                    results.push((((data[i][0] - firstelement[0]) / 1000), Math.ceil((data[i][1] / 1000000))));
                }
                return results;
            }

            function parseporehist(descriptions, counts) {
                var results = [];
                var colors = [];
                var categories = [];
                var datam = [];
                var colorlookup = [];
                //console.log(descriptions);
                //console.log(counts);
                //times=gettimelist(counts);
                //console.log(times);
                for (var thing in descriptions) {
                    //console.log(thing);
                    if (descriptions.hasOwnProperty(thing)) {
                        //console.log('here?');
                        if (descriptions[thing].hasOwnProperty("style")) {
                            //    console.log('and?');
                            //console.log(descriptions[thing]);
                            //console.log(descriptions[thing]["name"]);
                            //console.log(descriptions[thing]["style"]["colour"]);
                            colorlookup[descriptions[thing]["name"]] = descriptions[thing]["style"]["colour"];
                            //        console.log(descriptions[thing]["style"]["label"]);
                            //        console.log(descriptions[thing]["style"]["colour"]);
                        }
                    }
                }
                //console.log(colorlookup);
                for (var pore in counts) {
                    //console.log(pore);
                    //console.log(counts[pore]);
                    //console.log(colorlookup[pore]);
                    //results.push({"name":descriptions[thing]["style"]["label"], "data":[{"y":porenumber}],"color":"#"+descriptions[thing]["style"]["colour"]});
                    results.push({"name": pore, "color": "#" + colorlookup[pore], "data": counts[pore]})//,"color":"#121212"]});
                    //results.push({"name":pore,"data":[{100000,1},{200000,2}]});
                    //break;
                    //console.log(results);
                }
                /*for (var pore in counts){
                 //    console.log(pore);
                 //for (var timepoint in counts[pore]){
                 //console.log(counts[pore]);
                 //}
                 }*/
                return results
            }



            function formatBytes(bytes, decimals) {
                if (bytes == 0) return '0 Byte';
                var k = 1000;
                var dm = decimals + 1 || 3;
                var sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB'];
                var i = Math.floor(Math.log(bytes) / Math.log(k));
                return parseFloat((bytes / Math.pow(k, i)).toFixed(dm)) + ' ' + sizes[i];
            }

            function tohistogram(readeventcountweightedhist, readeventcountweightedhistbinwidth, totalyield) {
                var results = [];
                var categories = [];
                //var counter = 0;
                readeventcountweightedhist = readeventcountweightedhist.replace(/u(?=[^:]+')/g, "").replace(/'/g, "");
                readeventcountweightedhist = JSON.parse(readeventcountweightedhist);
                //console.log(readeventcountweightedhist);
                var n50count = 0;
                var n50index = 0;
                var check = 0;
                //console.log(readeventcountweightedhist);
                for (i in readeventcountweightedhist) {
                    //if (readeventcountweightedhist[i] > 0){
                    //counter+=1;
                    //console.log(readeventcountweightedhistbinwidth);
                    //console.log(i);

                    //console.log(i*readeventcountweightedhistbinwidth, readeventcountweightedhist[i]);
                    n50count += parseInt(readeventcountweightedhist[i]);
                    if (n50count >= (parseInt(totalyield) / 2)) {
                        //console.log('n50',(i+1)*readeventcountweightedhistbinwidth, n50count);
                        check += 1;
                    }
                    //console.log(i);
                    //console.log(parseInt(i)+1);
                    var category = String((parseInt(i)) * readeventcountweightedhistbinwidth) + " - " + String((parseInt(i) + 1) * readeventcountweightedhistbinwidth) + " ev";
                    categories.push(category);
                    if (check == 1) {
                        n50index = i;
                        results.push({"name": category, "y": parseInt(readeventcountweightedhist[i]), "color": 'red'});
                        check += 1;
                    } else {
                        results.push({"name": category, "y": parseInt(readeventcountweightedhist[i]), "color": 'blue'});
                    }

                    //}
                }
                categories.push(">> max ev");
                var missed = 0;
                //var missed = totalyield - readeventcountweightedhist.reduce(add,0);
                results.push({"name": ">> max ev", "y": missed});
                //console.log(n50index);
                return [results, categories, n50index];
            }

            var master_minion_list = [];

            console.log(master_minion_list);
            function start() {
                $.get("{% url 'activeminION_list' %}", function (data) {
                    //$( ".result" ).html( data
                    var data2 = [];
                    var ien = data.length;
                    //var counter = 0;
                    //console.log(data);
                    for (var i = 0; i < ien; i++) {
                        var minIONname = data[i].minION_name;
                        if ($.inArray(data[i].minION_name, master_minion_list) == -1) {
                            minionsthings.minions.push({
                                scripts: '',
                                status: data[i].run_status,
                                state: data[i].status,
                                sample_id: data[i].sample_name,
                                run_name: data[i].run_name,
                                asic_id: '',
                                histogram_bin_width: '',
                                histogram_values: [],
                                read_count: '',
                                start_time: '',
                                last_update: '',
                                flow_cell_id: data[i].flow_cell_id,
                                name: data[i].minION_name,
                                livedata: {
                                    voltage_value: data[i].voltage_value,
                                    event_yield: data[i].event_yield,
                                    machine_id: {
                                        result: data[i].computer
                                    },
                                    current_script: {
                                        result: data[i].currentscript
                                    },
                                    disk_space: {
                                        result: {
                                            bytes_capacity: data[i].total_drive_space,
                                            bytes_to_alert: data[i].space_till_shutdown,
                                            bytes_available: data[i].space_available,
                                            recommend_alert: data[i].warnings,
                                        }
                                    }
                                },
                                minKNOW_version: data[i].minKNOW_version,
                                voltage: [],
                                asictemp: [],
                                heatsinktemp: [],
                                strand: [],
                                good_single: [],
                                percentage: [],
                                currpercentage: [],
                                currstrand: [],
                                yield_history: [],
                                meanratio_history: [],
                                messages: [],
                                instrand_history: [],
                                openpore_history: [],
                                colour_stats: [],
                                pore_history: {
                                    above: [],
                                    adapter: [],
                                    below: [],
                                    good_single: [],
                                    strand: [],
                                    inrange: [],
                                    multiple: [],
                                    pending_mux_change: [],
                                    saturated: [],
                                    unavailable: [],
                                    unblocking: [],
                                    unclassified: [],
                                    unknown: [],
                                },
                                url: data[i].url,
                                last_counter: 0,
                            });
                            master_minion_list.push(data[i].minION_name);
                            $.get(data[i].url + "scripts", function (scriptstuff) {
                                //console.log(scriptstuff);
                                for (obj in minionsthings.minions) {
                                    for (j in data) {
                                        if (minionsthings.minions[obj].name === data[j].minION_name) {
                                            minionsthings.minions[obj].scripts = scriptstuff;
                                        }
                                    }

                                }

                            })
                        } else {
                            //console.log(data[i].minION_name);
                            for (obj in minionsthings.minions) {
                                if (minionsthings.minions[obj].name == data[i].minION_name) {
                                    minionsthings.minions[obj].status = data[i].run_status;
                                    minionsthings.minions[obj].state = data[i].status;
                                    minionsthings.minions[obj].sample_id = data[i].sample_name;
                                    minionsthings.minions[obj].run_name = data[i].run_name;
                                    minionsthings.minions[obj].flow_cell_id = data[i].flow_cell_id;
                                    minionsthings.minions[obj].name = data[i].minION_name;
                                    minionsthings.minions[obj].livedata.voltage_value = data[i].voltage_value;
                                    minionsthings.minions[obj].livedata.event_yield = data[i].event_yield;
                                    minionsthings.minions[obj].livedata.machine_id.result = data[i].computer;
                                    minionsthings.minions[obj].livedata.current_script.result = data[i].currentscript;
                                    minionsthings.minions[obj].livedata.disk_space.result.bytes_capacity = data[i].total_drive_space;
                                    minionsthings.minions[obj].livedata.disk_space.result.bytes_to_alert = data[i].space_till_shutdown;
                                    minionsthings.minions[obj].livedata.disk_space.result.bytes_available = data[i].space_available;
                                    minionsthings.minions[obj].livedata.disk_space.result.recommend_alert = data[i].warnings;
                                }
                            }

                        }
                        //}
                        //console.log(master_minion_list);
                        //data2.push([data[i].computer, data[i].minION_name, data[i].run_status]);
                        //get scripts for this minION

                        //console.log(i + " ival1");
                        if (data[i].last_run != "undefined") {
                            //console.log(data[i].last_run);
                            //console.log(i + " ival2");
                            for (obj in minionsthings.minions) {
                                for (j in data) {
                                    if (minionsthings.minions[obj].name === data[j].minION_name) {
                                        //console.log("found it");
                                        //console.log(minionsthings.minions[obj].last_counter);
                                        counter = minionsthings.minions[obj].last_counter;
                                    }
                                }
                            }
                            //console.log("http://" + window.location.hostname + ":8100/api/v1/runs/" + data[i].last_run + "/runstats/" + counter);
                            //$.get("http://" + window.location.hostname + ":8100/api/v1/runs/" + data[i].last_run + "/runstats/" + counter, function ( data3 ) {
                            $.get("/api/v1/runs/" + data[i].last_run + "/runstats/" + counter, function (data3) {
                                //console.log("http://" + window.location.hostname + ":8100/api/v1/runs/" + data[i].last_run + "/runstats/");
                                //console.log("sausage");
                                for (obj in minionsthings.minions) {
                                    for (j in data) {
                                        if (minionsthings.minions[obj].name === data[j].minION_name) {
                                            //console.log("updating historgram");
                                            minionsthings.minions[obj].histogram_bin_width = data3.slice(-1)[0].minKNOW_histogram_bin_width;
                                            //console.log(minionsthings.minions[obj].histogram_values);
                                            minionsthings.minions[obj].histogram_values = data3.slice(-1)[0].minKNOW_histogram_values;
                                            //console.log(minionsthings.minions[obj].histogram_values);
                                            minionsthings.minions[obj].read_count = data3.slice(-1)[0].minKNOW_read_count;
                                            minionsthings.minions[obj].last_update = Date.parse(data3.slice(-1)[0].sample_time);
                                            minionsthings.minions[obj].last_counter = data3.slice(-1)[0].id;
                                        }
                                    }
                                    for (element in data3) {
                                        minionsthings.minions[obj].voltage.push([Date.parse(data3[element].sample_time), data3[element].voltage_value]);
                                        minionsthings.minions[obj].asictemp.push([Date.parse(data3[element].sample_time), data3[element].asic_temp]);
                                        minionsthings.minions[obj].heatsinktemp.push([Date.parse(data3[element].sample_time), data3[element].heat_sink_temp]);
                                        minionsthings.minions[obj].strand.push([Date.parse(data3[element].sample_time), data3[element].strand]);
                                        minionsthings.minions[obj].good_single.push([Date.parse(data3[element].sample_time), data3[element].good_single]);
                                        minionsthings.minions[obj].currpercentage = data3[element].occupancy;
                                        minionsthings.minions[obj].currstrand = data3[element].strand;
                                        minionsthings.minions[obj].percentage.push([Date.parse(data3[element].sample_time), data3[element].occupancy]);
                                        minionsthings.minions[obj].yield_history.push([Date.parse(data3[element].sample_time), data3[element].event_yield]);
                                        minionsthings.minions[obj].meanratio_history.push([Date.parse(data3[element].sample_time), data3[element].mean_ratio]);
                                        minionsthings.minions[obj].instrand_history.push([Date.parse(data3[element].sample_time), data3[element].in_strand]);
                                        minionsthings.minions[obj].openpore_history.push([Date.parse(data3[element].sample_time), data3[element].open_pore]);
                                        var myStringArray = ["above", "adapter", "below", "good_single", "strand", "inrange", "multiple", "pending_mux_change", "saturated", "unavailable", "unblocking", "unclassified", "unknown"];
                                        var arrayLength = myStringArray.length;
                                        for (var i = 0; i < arrayLength; i++) {
                                            minionsthings.minions[obj].pore_history[myStringArray[i]].push([Date.parse(data3[element].sample_time), data3[element][myStringArray[i]]]);
                                        }
                                    }
                                }

                            });
                            $.get("/api/v1/runs/" + data[i].last_run + "/rundetails/", function (data4) {
                                for (obj in minionsthings.minions) {
                                    for (j in data) {
                                        if (minionsthings.minions[obj].name === data[j].minION_name) {
                                            //console.log(data[j].minION_name);
                                            minionsthings.minions[obj].asic_id = data4.slice(-1)[0].minKNOW_asic_id;
                                            minionsthings.minions[obj].start_time = Date.parse(data4.slice(-1)[0].minKNOW_start_time);
                                            minionsthings.minions[obj].colour_stats = JSON.parse("[" + data4.slice(-1)[0].minKNOW_colours_string.replace("u'", "'") + "]");
                                        }
                                    }

                                }
                            });
                            $.get(data[i].url + 'recentmessages', function (data5) {
                                for (obj in minionsthings.minions) {
                                    for (j in data5) {
                                        if (minionsthings.minions[obj].name === minIONname) {
                                            minionsthings.minions[obj].messages.push({
                                                severity: data5[j].minKNOW_severity,
                                                timestamp: data5[j].minKNOW_message_timestamp,
                                                message: data5[j].minKNOW_message
                                            });
                                        }
                                    }
                                }

                            })
                        }
                    }
                    console.log(minionsthings);
                });


            }

            start();
            setInterval(start, 10000);


            {% verbatim %}


            Vue.filter('currencyDisplay', {
                // model -> view
                // formats the value when updating the input element.
                read: function (val) {
                    return '$' + val.toFixed(2)
                },
                // view -> model
                // formats the value when writing to the data.
                write: function (val, oldVal) {
                    var number = +val.replace(/[^\d.]/g, '')
                    return isNaN(number) ? 0 : parseFloat(number.toFixed(2))
                }
            })

            Vue.component('projchartyield', {
                template: '<div id="projcontaineryield{{title}}" style="margin: 0 auto"</div>',
                props: ['title', 'key', 'datain', 'datain2', 'seqspeed', 'start_time'],
                data: function () {
                    //var d = new Date();
                    //var t = d.getTime();
                    return {
                        opts: {
                            chart: {
                                renderTo: 'projcontaineryield' + this.title,
                                type: 'spline',
                                zoomType: 'x',
                                height: 350,
                            },

                            title: {
                                text: 'Projected Yield over time '
                            },
                            xAxis: {
                                type: 'datetime',
                                tickPixelInterval: 150
                            },
                            yAxis: {
                                title: {
                                    text: 'Cumulative Predicted Bases'
                                },
                                plotLines: [{
                                    value: 0,
                                    width: 1,
                                    color: '#808080'
                                }],
                                min: 0,
                            },
                            credits: {
                                enabled: false
                            },
                            series: [{
                                //regression: true,
                                //regressionSettings: {
                                //type: 'polynomial',
                                //order: 4,
                                //color: 'rgba(223, 83, 83, .9)',
                                //extrapolate: 5
                                //},
                                name: 'Original Data',
                                data: []

                            },
                                {
                                    name: 'Current Projected Data',
                                    dashStyle: 'longdash',
                                    data: []
                                },
                                {
                                    name: 'First Hour Projected Data',
                                    dashStyle: 'shortdash',
                                    data: []
                                },
                                {
                                    name: 'Ideal Results',
                                    dashStyle: 'Dot',
                                    data: []
                                }
                            ]

                        }
                    }
                }
                ,


                ready: function () {
                    this.$nextTick(function () {
                        this.chart = new Highcharts.Chart(this.opts);
                        //this.chart = new Highcharts.stockChart(this.opts);
                        //this.chart.series[0].setData(this.datain2.slice(-15));
                        this.chart.series[0].setData(converttobases(this.datain2, this.seqspeed));
                        this.chart.redraw();

                        setInterval(function () {
                            var timeleft = geteighthours(this.datain2.slice(-1), this.start_time);
                            var firsthour = 120;
                            [scalingfactor, difference] = projectdata(this.datain2);
                            [scalingfactor2, difference2] = projectdata(this.datain2.slice(0, firsthour));
                            var syntheticdata = this.datain2;
                            newarray = projectresults(syntheticdata, scalingfactor, 4000, difference, this.start_time);
                            newarray1 = projectresults(syntheticdata.slice(0, firsthour), scalingfactor2, 4000, difference2, this.start_time);
                            newarray2 = projectresults(syntheticdata.slice(0, firsthour), 1, 4000, difference2, this.start_time);
                            this.chart.series[0].setData(converttobases(this.datain2, this.seqspeed));
                            this.chart.series[1].setData(converttobases(newarray, this.seqspeed));
                            this.chart.series[2].setData(converttobases(newarray1, this.seqspeed));
                            this.chart.series[3].setData(converttobases(newarray2, this.seqspeed));
                        }.bind(this), 15000);
                    });
                }
            })

            Vue.component('chartporehistdetails', {
                template: '<div id="container-porehist{{title}}" style="margin: 0 auto"</div>',
                props: ['title', 'key', 'datain', 'datain2'],
                data: function () {
                    return {
                        opts: {
                            chart: {
                                renderTo: 'container-porehist' + this.title,
                                type: 'area',
                                //zoomType: 'xy',
                                //animation: false
                            },
                            rangeSelector: {
                                enabled: false
                            },
                            title: {
                                text: 'Pore States'
                            },
                            xAxis: {
                                range: 1 * 3600 * 1000, //set range to last hour of data
                            },
                            //type: 'datetime',
                            //tickPixelInterval: 150
                            //},
                            //colors:[],
                            yAxis: {
                                //max: 512,
                                endOnTick: false,
                                title: {
                                    text: 'Channel Classifications'
                                }
                            },
                            legend: {
                                enabled: true
                            },
                            plotOptions: {
                                area: {
                                    stacking: 'percent',
                                },
                                series: {
                                    showInNavigator: true,
                                    dataLabels: {
                                        enabled: false,
                                        formatter: function () {
                                            return this.y;
                                        }
                                    }
                                }
                            },


                            credits: {
                                enabled: false
                            },
                            series: [],

                        }
                    }
                }
                ,


                created: function () {
                },
                ready: function () {
                    //var returndata=parsechanstats(this.datain,this.datain2);
                    this.$nextTick(function () {
                        this.chart = new Highcharts.stockChart(this.opts);
                        var returndata = parseporehist(this.datain[0], this.datain2);
                        while (this.chart.series.length > 0)
                            this.chart.series[0].remove(true);
                        for (var i = 0; i < returndata.length; i++) {
                            this.chart.addSeries(returndata[i]);
                        }
                        setInterval(function () {
                            var returndata = parseporehist(this.datain[0], this.datain2);
                            while (this.chart.series.length > 0)
                                this.chart.series[0].remove(true);
                            for (var i = 0; i < returndata.length; i++) {
                                this.chart.addSeries(returndata[i]);
                            }
                        }.bind(this), 10000);
                    });
                }
            })


            Vue.component('porecurrents', {
                template: '<div id="porecurrents{{title}}" style="margin: 0 auto"</div>',
                props: ['title', 'key', 'datain', 'datain2'],
                data: function () {
                    //var d = new Date();
                    //var t = d.getTime();
                    return {
                        opts: {
                            chart: {
                                renderTo: 'porecurrents' + this.title,
                                type: 'spline',
                                zoomType: 'x',
                                height: 350,
                            },
                            rangeSelector: {
                                enabled: true,
                                buttons: [{
                                    type: 'minute',
                                    count: 1,
                                    text: '1min'
                                }, {
                                    type: 'minute',
                                    count: 5,
                                    text: '5min'
                                }, {
                                    type: 'minute',
                                    count: 30,
                                    text: '1/2hr'
                                }, {
                                    type: 'minute',
                                    count: 60,
                                    text: '1hr'
                                }, {
                                    type: 'day',
                                    count: 0.5,
                                    text: '12hrs'
                                }, {
                                    type: 'day',
                                    count: 1,
                                    text: '1day'
                                }, {
                                    type: 'all',
                                    text: 'All'
                                }]
                            },
                            title: {
                                text: 'Pore State Currents'
                            },
                            xAxis: {
                                type: 'datetime',
                                tickPixelInterval: 150
                            },
                            yAxis: {
                                title: {
                                    text: 'Current'
                                },
                                plotLines: [{
                                    value: 0,
                                    width: 1,
                                    color: '#808080'
                                }],
                                startOnTick: false,
                                //min: 0,
                            },
                            credits: {
                                enabled: false
                            },
                            series: [
                                {
                                    name: 'In Strand',
                                    dataGrouping: {
                                        enabled: true
                                    },
                                    data: []
                                }, {
                                    name: 'Open Pore',
                                    dataGrouping: {
                                        enabled: true
                                    },
                                    data: []
                                }
                            ]
                        }
                    }
                }
                ,
                ready: function () {
                    this.$nextTick(function () {
                        this.chart = new Highcharts.stockChart(this.opts);
                        this.chart.series[0].setData(this.datain);
                        this.chart.series[1].setData(this.datain2);
                        setInterval(function () {
                            //console.log(this.datain);
                            this.chart.series[0].setData(this.datain);
                            this.chart.series[1].setData(this.datain2);
                            this.chart.redraw();
                        }.bind(this), 15000);
                    });
                }
            })

            Vue.component('meanratio', {
                template: '<div id="meanratio{{title}}" style="margin: 0 auto"</div>',
                props: ['title', 'key', 'datain'],
                data: function () {
                    //var d = new Date();
                    //var t = d.getTime();
                    return {
                        opts: {
                            chart: {
                                renderTo: 'meanratio' + this.title,
                                type: 'spline',
                                zoomType: 'x',
                                height: 350,
                            },
                            rangeSelector: {
                                enabled: true,
                                buttons: [{
                                    type: 'minute',
                                    count: 1,
                                    text: '1min'
                                }, {
                                    type: 'minute',
                                    count: 5,
                                    text: '5min'
                                }, {
                                    type: 'minute',
                                    count: 30,
                                    text: '1/2hr'
                                }, {
                                    type: 'minute',
                                    count: 60,
                                    text: '1hr'
                                }, {
                                    type: 'day',
                                    count: 0.5,
                                    text: '12hrs'
                                }, {
                                    type: 'day',
                                    count: 1,
                                    text: '1day'
                                }, {
                                    type: 'all',
                                    text: 'All'
                                }]
                            },
                            title: {
                                text: 'Current Ratio - In Strand/Open Pore'
                            },
                            xAxis: {
                                type: 'datetime',
                                tickPixelInterval: 150
                            },
                            yAxis: {
                                title: {
                                    text: 'Current Ratio'
                                },
                                plotLines: [{
                                    value: 0,
                                    width: 1,
                                    color: '#808080'
                                }],
                                startOnTick: false,
                                //min: 0,
                            },
                            credits: {
                                enabled: false
                            },
                            series: [
                                {
                                    name: 'Current Ratio',
                                    dataGrouping: {
                                        enabled: true
                                    },
                                    data: []
                                }
                            ]
                        }
                    }
                }
                ,
                ready: function () {
                    this.$nextTick(function () {
                        this.chart = new Highcharts.stockChart(this.opts);
                        this.chart.series[0].setData(this.datain);
                        setInterval(function () {
                            //console.log(this.datain);
                            this.chart.series[0].setData(this.datain);
                            this.chart.redraw();
                        }.bind(this), 15000);
                    });
                }
            })


            Vue.component('chartyield', {
                template: '<div id="containeryield{{title}}" style="margin: 0 auto"</div>',
                props: ['title', 'key', 'datain2'],
                data: function () {
                    //var d = new Date();
                    //var t = d.getTime();
                    return {
                        opts: {
                            chart: {
                                renderTo: 'containeryield' + this.title,
                                type: 'spline',
                                zoomType: 'xy'
                            },
                            rangeSelector: {
                                enabled: true,
                                buttons: [{
                                    type: 'minute',
                                    count: 1,
                                    text: '1min'
                                }, {
                                    type: 'minute',
                                    count: 5,
                                    text: '5min'
                                }, {
                                    type: 'minute',
                                    count: 30,
                                    text: '1/2hr'
                                }, {
                                    type: 'minute',
                                    count: 60,
                                    text: '1hr'
                                }, {
                                    type: 'day',
                                    count: 0.5,
                                    text: '12hrs'
                                }, {
                                    type: 'day',
                                    count: 1,
                                    text: '1day'
                                }, {
                                    type: 'all',
                                    text: 'All'
                                }]
                            },
                            title: {
                                text: 'Yield over time '
                            },
                            xAxis: {
                                type: 'datetime',
                                tickPixelInterval: 150
                            },
                            yAxis: {
                                title: {
                                    text: 'Cumulative Events'
                                },
                                plotLines: [{
                                    value: 0,
                                    width: 1,
                                    color: '#808080'
                                }],
                                min: 0,
                            },
                            credits: {
                                enabled: false
                            },
                            series: [{
                                name: 'Event Counts',
                                data: []
                            }]
                        }
                    }
                }
                ,


                ready: function () {
                    this.$nextTick(function () {
                        this.chart = new Highcharts.stockChart(this.opts);
                        //this.chart = new Highcharts.Chart(this.opts);
                        this.chart.series[0].setData(this.datain2);
                        setInterval(function () {
                            //console.log(this.datain2);
                            this.chart.series[0].setData(this.datain2);
                            this.chart.redraw();
                        }.bind(this), 15000);
                    });
                }
            })

            Vue.component('perchistory', {
                template: '<div id="perchistory{{title}}" style="margin: 0 auto"</div>',
                props: ['title', 'key', 'datain2'],
                data: function () {
                    //var d = new Date();
                    //var t = d.getTime();
                    return {
                        opts: {
                            chart: {
                                renderTo: 'perchistory' + this.title,
                                type: 'spline',
                                zoomType: 'xy'
                            },
                            rangeSelector: {
                                enabled: true,
                                buttons: [{
                                    type: 'minute',
                                    count: 1,
                                    text: '1min'
                                }, {
                                    type: 'minute',
                                    count: 5,
                                    text: '5min'
                                }, {
                                    type: 'minute',
                                    count: 30,
                                    text: '1/2hr'
                                }, {
                                    type: 'minute',
                                    count: 60,
                                    text: '1hr'
                                }, {
                                    type: 'day',
                                    count: 0.5,
                                    text: '12hrs'
                                }, {
                                    type: 'day',
                                    count: 1,
                                    text: '1day'
                                }, {
                                    type: 'all',
                                    text: 'All'
                                }]
                            },
                            title: {
                                text: '% Occupancy Over Time'
                            },
                            xAxis: {
                                type: 'datetime',
                                tickPixelInterval: 150
                            },
                            yAxis: {
                                title: {
                                    text: '% Occupancy'
                                },
                                plotLines: [{
                                    value: 0,
                                    width: 1,
                                    color: '#808080'
                                }],
                                min: 0,
                            },
                            credits: {
                                enabled: false
                            },
                            series: [{
                                name: '% Occupancy',
                                data: []
                            }
                            ]
                        }
                    }
                }
                ,
                ready: function () {
                    this.$nextTick(function () {
                        this.chart = new Highcharts.stockChart(this.opts);
                        //this.chart = new Highcharts.Chart(this.opts);
                        this.chart.series[0].setData(this.datain2);
                        setInterval(function () {
                            //console.log(this.datain2);
                            this.chart.series[0].setData(this.datain2);
                            this.chart.redraw();
                        }.bind(this), 15000);
                    });
                }
            })

            Vue.component('porehistory', {
                template: '<div id="porehistory{{title}}" style="margin: 0 auto"</div>',
                props: ['title', 'key', 'datain2', 'datain3'],
                data: function () {
                    //var d = new Date();
                    //var t = d.getTime();
                    return {
                        opts: {
                            chart: {
                                renderTo: 'porehistory' + this.title,
                                type: 'spline',
                                zoomType: 'xy'
                            },
                            rangeSelector: {
                                enabled: true,
                                buttons: [{
                                    type: 'minute',
                                    count: 1,
                                    text: '1min'
                                }, {
                                    type: 'minute',
                                    count: 5,
                                    text: '5min'
                                }, {
                                    type: 'minute',
                                    count: 30,
                                    text: '1/2hr'
                                }, {
                                    type: 'minute',
                                    count: 60,
                                    text: '1hr'
                                }, {
                                    type: 'day',
                                    count: 0.5,
                                    text: '12hrs'
                                }, {
                                    type: 'day',
                                    count: 1,
                                    text: '1day'
                                }, {
                                    type: 'all',
                                    text: 'All'
                                }]
                            },
                            title: {
                                text: 'In Strand Counts'
                            },
                            xAxis: {
                                type: 'datetime',
                                tickPixelInterval: 150
                            },
                            yAxis: {
                                title: {
                                    text: 'Number of Pores In Strand/Single'
                                },
                                plotLines: [{
                                    value: 0,
                                    width: 1,
                                    color: '#808080'
                                }],
                                min: 0,
                            },
                            credits: {
                                enabled: false
                            },
                            series: [
                                {
                                    name: 'In Strand',
                                    data: []
                                }, {
                                    name: 'Single Pore',
                                    data: []
                                }
                            ]
                        }
                    }
                }
                ,
                ready: function () {
                    this.$nextTick(function () {
                        //this.chart = new Highcharts.Chart(this.opts);
                        this.chart = new Highcharts.stockChart(this.opts);
                        this.chart.series[0].setData(this.datain2);
                        this.chart.series[1].setData(this.datain3);
                        setInterval(function () {
                            //console.log(this.datain2);
                            this.chart.series[0].setData(this.datain2);
                            this.chart.series[1].setData(this.datain3);
                            this.chart.redraw();
                        }.bind(this), 15000);
                    });
                }
            })

            Vue.component('temphistory', {
                template: '<div id="temphistory{{title}}" style="margin: 0 auto"</div>',
                props: ['title', 'key', 'datain2', 'datain3'],
                data: function () {
                    //var d = new Date();
                    //var t = d.getTime();
                    return {
                        opts: {
                            chart: {
                                renderTo: 'temphistory' + this.title,
                                type: 'spline',
                                zoomType: 'xy'
                            },
                            rangeSelector: {
                                enabled: true,
                                buttons: [{
                                    type: 'minute',
                                    count: 1,
                                    text: '1min'
                                }, {
                                    type: 'minute',
                                    count: 5,
                                    text: '5min'
                                }, {
                                    type: 'minute',
                                    count: 30,
                                    text: '1/2hr'
                                }, {
                                    type: 'minute',
                                    count: 60,
                                    text: '1hr'
                                }, {
                                    type: 'day',
                                    count: 0.5,
                                    text: '12hrs'
                                }, {
                                    type: 'day',
                                    count: 1,
                                    text: '1day'
                                }, {
                                    type: 'all',
                                    text: 'All'
                                }]
                            },
                            title: {
                                text: 'Temperature over time '
                            },
                            xAxis: {
                                type: 'datetime',
                                tickPixelInterval: 150
                            },
                            yAxis: {
                                title: {
                                    text: 'Cumulative Events'
                                },
                                plotLines: [{
                                    value: 0,
                                    width: 1,
                                    color: '#808080'
                                }],
                                min: 0,
                            },
                            credits: {
                                enabled: false
                            },
                            series: [{
                                name: 'Asic Temperature',
                                data: []
                            }, {
                                name: 'Heat Sink Temperature',
                                data: []
                            }
                            ]
                        }
                    }
                }
                ,
                ready: function () {
                    this.$nextTick(function () {
                        //this.chart = new Highcharts.Chart(this.opts);
                        this.chart = new Highcharts.stockChart(this.opts);
                        this.chart.series[0].setData(this.datain2);
                        this.chart.series[1].setData(this.datain3);
                        setInterval(function () {
                            //console.log(this.datain2["asictemp"]);
                            this.chart.series[0].setData(this.datain2);
                            this.chart.series[1].setData(this.datain3);
                            this.chart.redraw();
                        }.bind(this), 15000);
                    });
                }
            })

            Vue.component('volthistory', {
                template: '<div id="volthistory{{title}}" style="margin: 0 auto"</div>',
                props: ['title', 'key', 'datain2'],
                data: function () {
                    //var d = new Date();
                    //var t = d.getTime();
                    return {
                        opts: {
                            chart: {
                                renderTo: 'volthistory' + this.title,
                                type: 'spline',
                                zoomType: 'x',
                                //height: 350,
                            },
                            rangeSelector: {
                                enabled: true,
                                buttons: [{
                                    type: 'minute',
                                    count: 1,
                                    text: '1min'
                                }, {
                                    type: 'minute',
                                    count: 5,
                                    text: '5min'
                                }, {
                                    type: 'minute',
                                    count: 30,
                                    text: '1/2hr'
                                }, {
                                    type: 'minute',
                                    count: 60,
                                    text: '1hr'
                                }, {
                                    type: 'day',
                                    count: 0.5,
                                    text: '12hrs'
                                }, {
                                    type: 'day',
                                    count: 1,
                                    text: '1day'
                                }, {
                                    type: 'all',
                                    text: 'All'
                                }]
                            },
                            title: {
                                text: 'Global Voltage over time '
                            },
                            xAxis: {
                                type: 'datetime',
                                tickPixelInterval: 150
                            },
                            yAxis: {
                                title: {
                                    text: 'Voltage'
                                },
                                plotLines: [{
                                    value: 0,
                                    width: 1,
                                    color: '#808080'
                                }],
                                //min: 0,
                            },
                            credits: {
                                enabled: false
                            },
                            series: [{
                                name: 'Global Voltage',
                                dataGrouping: {
                                    enabled: true
                                },
                                data: []
                            }
                            ]
                        }
                    }
                }
                ,
                ready: function () {
                    this.$nextTick(function () {
                        //console.log(this.title);
                        //console.log(this.datain2);
                        //console.log(this.key);
                        this.chart = new Highcharts.stockChart(this.opts);
                        this.chart.series[0].setData(this.datain2["voltage"]);
                        //this.chart.series[1].setData(this.datain2["heatsinktemp"]);
                        setInterval(function () {
                            //console.log(this.datain2["asictemp"]);
                            //console.log(this.datain2);
                            this.chart.series[0].setData(this.datain2);
                            //this.chart.series[1].setData(this.datain2["heatsinktemp"]);
                            this.chart.redraw();
                        }.bind(this), 15000);
                    });
                }
            })

            Vue.component('chartporehist', {
                template: '<div id="container-pore{{title}}" style="margin: 0 auto"</div>',
                props: ['title', 'key', 'datain', 'datain2'],
                data: function () {
                    return {
                        opts: {
                            chart: {
                                renderTo: 'container-pore' + this.title,
                                type: 'column',
                                zoomType: 'xy',
                                height: 350,
                                animation: false
                            },
                            title: {
                                text: 'Pore States'
                            },
                            xAxis: {
                                categories: []
                            },
                            //colors:[],
                            yAxis: {
                                title: {
                                    //text: 'Pore Type'
                                }
                            },
                            legend: {
                                enabled: false
                            },
                            plotOptions: {
                                series: {
                                    dataLabels: {
                                        enabled: true,
                                        formatter: function () {
                                            return this.y;
                                        }
                                    }
                                }
                            },
                            credits: {
                                enabled: false
                            },
                            series: [{
                                name: 'porestates',
                                data: [],
                            }
                            ],

                        }
                    }
                }
                ,


                created: function () {
                },
                ready: function () {
                    //var returndata=parsechanstats(this.datain,this.datain2);
                    //console.log(returndata);
                    //var returndata = tohistogram(this.datain,parseInt(this.datain2));
                    //this.chart.series[0].setData(returndata[3]);
                    //console.log(returndata[2]);
                    //this.chart.xAxis[0].setCategories(returndata[2]);
                    //var returndata=parsechanstats(this.datain,this.datain2);
                    this.$nextTick(function () {
                        this.chart = new Highcharts.Chart(this.opts);
                        //minion=this.key;
                        setInterval(function () {
                            //console.log(this.datain);
                            //console.log(this.datain2);
                            var returndata = parsechanstats(this.datain, this.datain2);
                            //console.log(returndata);
                            //var returndata = tohistogram(this.datain,parseInt(this.datain2));
                            this.chart.series[0].setData(returndata[3]);
                            //console.log(returndata[2]);
                            this.chart.xAxis[0].setCategories(returndata[2]);
                            //while(this.chart.series.length > 0)
                            //    this.chart.series[0].remove(true);
                            //for (var i = 0; i < returndata[0].length; i++) {
                            //    this.chart.addSeries(returndata[0][i]);
                            //}
                            //this.chart.colors=returndata[1];
                            //this.chart.redraw();
                            //console.log(returndata[1]);
                        }.bind(this), 15000);
                    });
                }
            })

            function parsechanstats(descriptions, counts) {
                var results = [];
                var colors = [];
                var categories = [];
                var datam = [];
                //console.log(descriptions);
                //console.log(counts);

                //for (var i = 0; i < descriptions.length; i++) {
                for (var thing in descriptions) {
                    //console.log("in decsriptions");
                    if (descriptions.hasOwnProperty(thing)) {
                        //console.log(thing);
                        //console.log(descriptions[thing]);
                        for (var element in descriptions[thing]) {
                            if (descriptions[thing][element].hasOwnProperty("style")) {
                                //console.log(descriptions[thing][element]["style"]["label"]);
                                //console.log(descriptions[thing][element]["style"]["colour"]);

                                if (counts.hasOwnProperty(descriptions[thing][element]["name"])) {
                                    var porenumber = counts[descriptions[thing][element]["name"]].slice(-1)[0][1];
                                    //console.log(counts[descriptions[thing][element]["name"]].slice(-1)[0][1]);
                                } else {
                                    var porenumber = 0;
                                    //console.log("0");
                                }

                                results.push({
                                    "name": descriptions[thing][element]["style"]["label"],
                                    "data": [{"y": porenumber}],
                                    "color": "#" + descriptions[thing][element]["style"]["colour"]
                                });
                                colors.push("#" + descriptions[thing][element]["style"]["colour"]);
                                categories.push(descriptions[thing][element]["style"]["label"]);
                                datam.push({
                                    "y": porenumber,
                                    "color": "#" + descriptions[thing][element]["style"]["colour"]
                                });

                            }
                        }
                    }
                }
                return [results, colors, categories, datam];
            }

            Vue.component('container-avg', {
                template: '<div id="container-avg{{title}}" style="height: 140px; margin: 0 auto"</div>',
                props: ['title', 'key', 'datain', 'datain2'],
                data: function () {
                    return {

                        opts: {
                            chart: {
                                renderTo: 'container-avg' + this.title,
                                type: 'solidgauge'
                            },

                            title: "Average Read Length (Events)",

                            pane: {
                                center: ['50%', '85%'],
                                size: '100%',
                                startAngle: -90,
                                endAngle: 90,
                                background: {
                                    backgroundColor: (Highcharts.theme && Highcharts.theme.background2) || '#EEE',
                                    innerRadius: '60%',
                                    outerRadius: '100%',
                                    shape: 'arc'
                                }
                            },

                            tooltip: {
                                enabled: false
                            },

                            // the value axis
                            yAxis: {
                                type: 'logarithmic',
                                stops: [
                                    [0.3, '#0000FF'], // blue
                                    [0.37, '#DDDF0D'], // green
                                    [0.43, '#DF5353'], // red
                                ],
                                lineWidth: 0,
                                minorTickInterval: null,
                                tickPixelInterval: 400,
                                tickWidth: 0,
                                title: {
                                    y: -70
                                },
                                labels: {
                                    y: 16
                                },
                                min: 0.1,
                                max: 100000,
                                title: {
                                    text: null
                                }
                            },

                            plotOptions: {
                                solidgauge: {
                                    dataLabels: {
                                        y: -30,
                                        borderWidth: 0,
                                        useHTML: true
                                    }
                                }
                            },


                            credits: {
                                enabled: false
                            },

                            series: [{
                                name: 'Events',
                                data: [0],
                                dataLabels: {
                                    format: '<div style="text-align:center"><span style="font-size:15px;color:' +
                                    ((Highcharts.theme && Highcharts.theme.contrastTextColor) || 'black') + '">{y}</span><br/>' +
                                    '<span style="font-size:12px;color:silver"> Avg event len</span></div>'
                                },
                                tooltip: {
                                    valueSuffix: ' events'
                                }
                            }],

                            plotOptions: {
                                solidgauge: {
                                    dataLabels: {
                                        y: 30,
                                        borderWidth: 0,
                                        useHTML: true
                                    }
                                }
                            }
                        }
                    }
                }


                ,
                ready: function () {
                    //alert(this.datain);

                    this.$nextTick(function () {
                        this.chart = new Highcharts.Chart(this.opts);
                        //this.chart.series[0].setData(this.datain);
                        if (this.chart) {
                            //point = this.chart.series[0].points[0];

                            this.chart.series[0].points[0].update(round(parseFloat(this.datain / this.datain2), 0));
                            //point.update(this.datain);
                            //alert("camel");
                        }

                        setInterval(function () {
                            //    this.chart.series[0].setData(this.datain);
                            //    this.chart.redraw();
                            if (this.chart) {
                                point = this.chart.series[0].points[0];
                                point.update(round(parseFloat(this.datain / this.datain2), 0));
                            }
                        }.bind(this), 5000);
                    });
                }
            }),

                Vue.component('container-chan', {
                    template: '<div id="container-chan{{title}}" style="height: 140px; margin: 0 auto"</div>',
                    props: ['title', 'key', 'datain'],
                    data: function () {
                        return {
                            opts: {
                                chart: {
                                    renderTo: 'container-chan' + this.title,
                                    type: 'solidgauge'
                                },

                                title: null,

                                pane: {
                                    center: ['50%', '85%'],
                                    size: '100%',
                                    startAngle: -90,
                                    endAngle: 90,
                                    background: {
                                        backgroundColor: (Highcharts.theme && Highcharts.theme.background2) || '#EEE',
                                        innerRadius: '60%',
                                        outerRadius: '100%',
                                        shape: 'arc'
                                    }
                                },

                                tooltip: {
                                    enabled: false
                                },

                                // the value axis
                                yAxis: {
                                    stops: [


                                        [0.5, '#DF5353'], // red
                                        [0.75, '#DDDF0D'], // yellow
                                        [0.9, '#55BF3B'], // green
                                    ],
                                    lineWidth: 0,
                                    minorTickInterval: null,
                                    tickPixelInterval: 400,
                                    tickWidth: 0,
                                    title: {
                                        y: -70
                                    },
                                    labels: {
                                        y: 16
                                    },
                                    min: 0,
                                    max: 512,
                                    title: {
                                        text: null
                                    }
                                },

                                plotOptions: {
                                    solidgauge: {
                                        dataLabels: {
                                            y: -30,
                                            borderWidth: 0,
                                            useHTML: true
                                        }
                                    }
                                },


                                credits: {
                                    enabled: false
                                },

                                series: [{
                                    name: 'Used Channels',
                                    data: [0],
                                    dataLabels: {
                                        format: '<div style="text-align:center"><span style="font-size:15px;color:' +
                                        ((Highcharts.theme && Highcharts.theme.contrastTextColor) || 'black') + '">{y}</span><br/>' +
                                        '<span style="font-size:12px;color:silver"> Used Channels</span></div>'
                                    },
                                    tooltip: {
                                        valueSuffix: ' Channel Count'
                                    }
                                }],

                                plotOptions: {
                                    solidgauge: {
                                        dataLabels: {
                                            y: 30,
                                            borderWidth: 0,
                                            useHTML: true
                                        }
                                    }
                                }
                            }
                        }
                    }
                    ,
                    ready: function () {
                        //alert(this.datain);

                        this.$nextTick(function () {
                            this.chart = new Highcharts.Chart(this.opts);
                            //this.chart.series[0].setData(this.datain);
                            if (this.chart) {
                                //point = this.chart.series[0].points[0];
                                this.chart.series[0].points[0].update(this.datain);
                                //point.update(this.datain);
                                //alert("camel");
                            }

                            setInterval(function () {
                                //    this.chart.series[0].setData(this.datain);
                                //    this.chart.redraw();
                                if (this.chart) {
                                    point = this.chart.series[0].points[0];

                                    point.update(parseFloat(this.datain));
                                    //this.chart.redraw();


                                }
                            }.bind(this), 5000);
                        });
                    }
                }),


                Vue.component('container-strand', {
                    template: '<div id="container-strand{{title}}" style="height: 140px; margin: 0 auto"</div>',
                    props: ['title', 'key', 'datain'],
                    data: function () {
                        return {
                            opts: {
                                chart: {
                                    renderTo: 'container-strand' + this.title,
                                    type: 'solidgauge'
                                },

                                title: null,

                                pane: {
                                    center: ['50%', '85%'],
                                    size: '100%',
                                    startAngle: -90,
                                    endAngle: 90,
                                    background: {
                                        backgroundColor: (Highcharts.theme && Highcharts.theme.background2) || '#EEE',
                                        innerRadius: '60%',
                                        outerRadius: '100%',
                                        shape: 'arc'
                                    }
                                },

                                tooltip: {
                                    enabled: false
                                },

                                // the value axis
                                yAxis: {
                                    stops: [


                                        [0.5, '#DF5353'], // red
                                        [0.75, '#DDDF0D'], // yellow
                                        [0.9, '#55BF3B'], // green
                                    ],
                                    lineWidth: 0,
                                    minorTickInterval: null,
                                    tickPixelInterval: 400,
                                    tickWidth: 0,
                                    title: {
                                        y: -70
                                    },
                                    labels: {
                                        y: 16
                                    },
                                    min: 0,
                                    max: 512,
                                    title: {
                                        text: null
                                    }
                                },

                                plotOptions: {
                                    solidgauge: {
                                        dataLabels: {
                                            y: -30,
                                            borderWidth: 0,
                                            useHTML: true
                                        }
                                    }
                                },


                                credits: {
                                    enabled: false
                                },

                                series: [{
                                    name: 'In Strand',
                                    data: [0],
                                    dataLabels: {
                                        format: '<div style="text-align:center"><span style="font-size:15px;color:' +
                                        ((Highcharts.theme && Highcharts.theme.contrastTextColor) || 'black') + '">{y}</span><br/>' +
                                        '<span style="font-size:12px;color:silver"> In Strand</span></div>'
                                    },
                                    tooltip: {
                                        valueSuffix: ' In Strand'
                                    }
                                }],

                                plotOptions: {
                                    solidgauge: {
                                        dataLabels: {
                                            y: 30,
                                            borderWidth: 0,
                                            useHTML: true
                                        }
                                    }
                                }
                            }
                        }
                    }
                    ,
                    ready: function () {
                        //alert(this.datain);

                        this.$nextTick(function () {
                            this.chart = new Highcharts.Chart(this.opts);
                            if (this.chart) {
                                this.chart.series[0].points[0].update(this.datain);
                            }

                            setInterval(function () {
                                if (this.chart) {
                                    point = this.chart.series[0].points[0];
                                    point.update(parseFloat(this.datain));
                                }
                            }.bind(this), 5000);
                        });
                    }
                }),

                Vue.component('container-perc', {
                    template: '<div id="container-perc{{title}}" style="height: 140px; margin: 0 auto"</div>',
                    props: ['title', 'key', 'datain'],
                    data: function () {
                        return {
                            opts: {
                                chart: {
                                    renderTo: 'container-perc' + this.title,
                                    type: 'solidgauge'
                                },

                                title: null,

                                pane: {
                                    center: ['50%', '85%'],
                                    size: '100%',
                                    startAngle: -90,
                                    endAngle: 90,
                                    background: {
                                        backgroundColor: (Highcharts.theme && Highcharts.theme.background2) || '#EEE',
                                        innerRadius: '60%',
                                        outerRadius: '100%',
                                        shape: 'arc'
                                    }
                                },

                                tooltip: {
                                    enabled: false
                                },

                                // the value axis
                                yAxis: {
                                    stops: [


                                        [0.5, '#DF5353'], // red
                                        [0.75, '#DDDF0D'], // yellow
                                        [0.9, '#55BF3B'], // green
                                    ],
                                    lineWidth: 0,
                                    minorTickInterval: null,
                                    tickPixelInterval: 400,
                                    tickWidth: 0,
                                    title: {
                                        y: -70
                                    },
                                    labels: {
                                        y: 16
                                    },
                                    min: 0,
                                    max: 100,
                                    title: {
                                        text: null
                                    }
                                },

                                plotOptions: {
                                    solidgauge: {
                                        dataLabels: {
                                            y: -30,
                                            borderWidth: 0,
                                            useHTML: true
                                        }
                                    }
                                },


                                credits: {
                                    enabled: false
                                },

                                series: [{
                                    name: '% Occupancy',
                                    data: [0],
                                    dataLabels: {
                                        format: '<div style="text-align:center"><span style="font-size:15px;color:' +
                                        ((Highcharts.theme && Highcharts.theme.contrastTextColor) || 'black') + '">{y}</span><br/>' +
                                        '<span style="font-size:12px;color:silver"> % Occupancy</span></div>'
                                    },
                                    tooltip: {
                                        valueSuffix: ' % Occupancy'
                                    }
                                }],

                                plotOptions: {
                                    solidgauge: {
                                        dataLabels: {
                                            y: 30,
                                            borderWidth: 0,
                                            useHTML: true
                                        }
                                    }
                                }
                            }
                        }
                    }
                    ,
                    ready: function () {
                        //alert(this.datain);

                        this.$nextTick(function () {
                            this.chart = new Highcharts.Chart(this.opts);
                            //this.chart.series[0].setData(this.datain);
                            if (this.chart) {
                                //point = this.chart.series[0].points[0];
                                this.chart.series[0].points[0].update(this.datain);
                                //point.update(this.datain);
                                //alert("camel");
                            }

                            setInterval(function () {
                                //    this.chart.series[0].setData(this.datain);
                                //    this.chart.redraw();
                                if (this.chart) {
                                    point = this.chart.series[0].points[0];
                                    //console.log(this.datain);
                                    //var single = 0;
                                    //if (parseFloat(this.datain["good_single"]) > 0) {
                                    //    single = parseFloat(this.datain["good_single"]);
                                    //}
                                    //round()
                                    point.update(round(parseFloat(this.datain), 2));
                                    //point.update(round(((parseFloat(this.datain["strand"])/(parseFloat(this.datain["strand"])+single))*100),2));
                                    //this.chart.redraw();


                                }
                            }.bind(this), 5000);
                        });
                    }
                }),

                Vue.component('predictedvals', {
                    template: '<div class="col-md-2"><p>Estimated Current Yield (bases) <strong>  {{ speedresult }}  </strong></p></div><div class="col-md-2"><p>Estimated Average Read Length (bases) <strong>{{averageresult}}</strong></p></div><div class="col-md-2"><p>Theoretical Predicted Yield at 24 hours (bases) <strong>{{startbit}}</strong></p></div><div class="col-md-2"><p> Theoretical Predicted Yield at 48 hours (bases) <strong>{{endbit}}</strong></p></div>',
                    props: ['seqspeed', 'currentyield', 'compreads', 'calcurrenttime', 'calcstarttime'],
                    computed: {
                        speedresult: function () {
                            switch (this.seqspeed) {
                                case "MegaCrazy Runs":
                                    scaling = 3.5;
                                    break;
                                case "450 b/s":
                                    scaling = 1.8;
                                    break;
                                case "250 b/s":
                                    scaling = 1.1;
                                    break;
                                case "70 b/s":
                                    scaling = 1.0;
                                    break;
                            }

                            return numberWithCommas(round(parseFloat(scaling * this.currentyield), 0));
                        },
                        startbit: function () {
                            switch (this.seqspeed) {
                                case "MegaCrazy Runs":
                                    scaling = 3.5;
                                    break;
                                case "450 b/s":
                                    scaling = 1.8;
                                    break;
                                case "250 b/s":
                                    scaling = 1.1;
                                    break;
                                case "70 b/s":
                                    scaling = 1.0;
                                    break;
                            }
                            var startthing = parseInt(this.calcstarttime);
                            var endthing = new Date();
                            //console.log(endthing-startthing);
                            //return round(parseFloat(scaling * this.currentyield),0)/(endthing-startthing);
                            return numberWithCommas(round(round(parseFloat(scaling * this.currentyield), 0) / ((endthing - startthing) / 1000) * (24 * 60 * 60), 0));
                        },
                        endbit: function () {
                            switch (this.seqspeed) {
                                case "MegaCrazy Runs":
                                    scaling = 3.5;
                                    break;
                                case "450 b/s":
                                    scaling = 1.8;
                                    break;
                                case "250 b/s":
                                    scaling = 1.1;
                                    break;
                                case "70 b/s":
                                    scaling = 1.0;
                                    break;
                            }
                            var startthing = parseInt(this.calcstarttime);
                            var endthing = new Date();

                            return numberWithCommas(round(round(parseFloat(scaling * this.currentyield), 0) / ((endthing - startthing) / 1000) * (48 * 60 * 60), 0));
                        },

                        averageresult: function () {
                            switch (this.seqspeed) {
                                case "MegaCrazy Runs":
                                    scaling = 3.5;
                                    break;
                                case "450 b/s":
                                    scaling = 1.8;
                                    break;
                                case "250 b/s":
                                    scaling = 1.1;
                                    break;
                                case "70 b/s":
                                    scaling = 1.0;
                                    break;
                            }
                            return numberWithCommas(round(parseFloat(scaling * (this.currentyield / this.compreads)), 0));
                        }
                    }
                })

            {% endverbatim %}


            function getCookie(name) {
                var cookieValue = null;
                if (document.cookie && document.cookie !== '') {
                    var cookies = document.cookie.split(';');
                    for (var i = 0; i < cookies.length; i++) {
                        var cookie = jQuery.trim(cookies[i]);
                        // Does this cookie string begin with the name we want?
                        if (cookie.substring(0, name.length + 1) === (name + '=')) {
                            cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                            break;
                        }
                    }
                }
                return cookieValue;
            }



            var minionsthings = new Vue({
                el: '#app',
                data: {
                    minions: [],
                    message: '',
                },
                computed: {
                    speedresult: function () {
                        switch (this.seqspeed) {
                            case "MegaCrazy Runs":
                                scaling = 3.5;
                                break;
                            case "450 b/s":
                                scaling = 1.8;
                                break;
                            case "250 b/s":
                                scaling = 1.1;
                                break;
                            case "70 b/s":
                                scaling = 1.0;
                                break;
                        }
                        return this.currentyield;
                    }
                },
                methods: {
                    testmessage: function (event) {
                        var instructionmessage = {
                            "INSTRUCTION": {
                                "USER": "<?php echo $_SESSION['user_name'];?>",
                                "minion": event.target.id,
                                "JOB": "testmessage"
                            }
                        };
                        var minionid = event.target.id.split("/");
                        //var csrftoken = getCookie('csrftoken');
                        function csrfSafeMethod(method) {
                            return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
                        }

                        $.ajaxSetup({
                            beforeSend: function (xhr, settings) {
                                if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                                    xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
                                }
                            }
                        });
                        var url_mask = "{% url 'minIONcontrol-list' pk=12345 %}".replace(/12345/, minionid[minionid.length - 2].toString());
                        $.ajax({
                            "type": "POST",
                            "dataType": "json",
                            "url": url_mask,
                            "data": {
                                "custom": "",
                                "job": "testmessage",
                                "minION": event.target.id,
                                "complete": 'False'
                            },
                            "beforeSend": function (xhr, settings) {
                                console.log("before send");
                                $.ajaxSettings.beforeSend(xhr, settings);
                            },
                            "success": function (result) {
                                console.log(result);
                            },
                            error: function (XMLHttpRequest, textStatus, errorThrown) {
                                console.log("Status: " + textStatus);
                                alert("Error: " + errorThrown);
                            }

                        })
                        //ws.send(JSON.stringify(instructionmessage));
                    },

                    renamenow: function (event) {
                        var minionid = event.target.id.split("/");

                        function csrfSafeMethod(method) {
                            return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
                        }

                        $.ajaxSetup({
                            beforeSend: function (xhr, settings) {
                                if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                                    xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
                                }
                            }
                        });
                        var url_mask = "{% url 'minIONcontrol-list' pk=12345 %}".replace(/12345/, minionid[minionid.length - 2].toString());
                        var stringthing = event.target.name;
                        var message = $("#" + stringthing + "newname").val();
                        $.ajax({
                            "type": "POST",
                            "dataType": "json",
                            "url": url_mask,
                            "data": {
                                "custom": message,
                                "job": "rename",
                                "minION": event.target.id,
                                "complete": 'False'
                            },
                            "beforeSend": function (xhr, settings) {
                                console.log("before send");
                                $.ajaxSettings.beforeSend(xhr, settings);
                            },
                            "success": function (result) {
                                console.log(result);
                            },
                            error: function (XMLHttpRequest, textStatus, errorThrown) {
                                console.log("Status: " + textStatus);
                                alert("Error: " + errorThrown);
                            }

                        })
                        //ws.send(JSON.stringify(instructionmessage));
                        //$('#renamemodal').modal('hide');
                    },
                    renameflowcellnow: function (event) {
                        //var instructionmessage={"INSTRUCTION":{"USER":"<?php echo $_SESSION['user_name'];?>","minion":event.target.id,"JOB":"nameflowcell","NAME":$("#"+event.target.id+"newflowcellname").val()}};
                        var minionid = event.target.id.split("/");

                        function csrfSafeMethod(method) {
                            return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
                        }

                        $.ajaxSetup({
                            beforeSend: function (xhr, settings) {
                                if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                                    xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
                                }
                            }
                        });
                        var url_mask = "{% url 'minIONcontrol-list' pk=12345 %}".replace(/12345/, minionid[minionid.length - 2].toString());
                        var stringthing = event.target.name;
                        var message = $("#" + stringthing + "newflowcellname").val();
                        $.ajax({
                            "type": "POST",
                            "dataType": "json",
                            "url": url_mask,
                            "data": {
                                "custom": message,
                                "job": "nameflowcell",
                                "minION": event.target.id,
                                "complete": 'False'
                            },
                            "beforeSend": function (xhr, settings) {
                                console.log("before send");
                                $.ajaxSettings.beforeSend(xhr, settings);
                            },
                            "success": function (result) {
                                console.log(result);
                            },
                            error: function (XMLHttpRequest, textStatus, errorThrown) {
                                console.log("Status: " + textStatus);
                                alert("Error: " + errorThrown);
                            }

                        })


                    },
                    startminion: function (event) {
                        var script = $("input[type='radio'][name='scriptRadios']:checked").val();
                        //alert(script);
                        var instructionmessage = {
                            "INSTRUCTION": {
                                "USER": "<?php echo $_SESSION['user_name'];?>",
                                "minion": event.target.id,
                                "JOB": "startminion",
                                "SCRIPT": script
                            }
                        };
                        var minionid = event.target.id.split("/");

                        function csrfSafeMethod(method) {
                            return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
                        }

                        $.ajaxSetup({
                            beforeSend: function (xhr, settings) {
                                if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                                    xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
                                }
                            }
                        });
                        var url_mask = "{% url 'minIONcontrol-list' pk=12345 %}".replace(/12345/, minionid[minionid.length - 2].toString());
                        $.ajax({
                            "type": "POST",
                            "dataType": "json",
                            "url": url_mask,
                            "data": {
                                "custom": script,
                                "job": "startminion",
                                "minION": event.target.id,
                                "complete": 'False'
                            },
                            "beforeSend": function (xhr, settings) {
                                console.log("before send");
                                $.ajaxSettings.beforeSend(xhr, settings);
                            },
                            "success": function (result) {
                                console.log(result);
                            },
                            error: function (XMLHttpRequest, textStatus, errorThrown) {
                                console.log("Status: " + textStatus);
                                alert("Error: " + errorThrown);
                            }

                        })
                    },
                    custommessage: function (event) {
                        //var script = $("input[type='radio'][name='scriptRadios']:checked").val();
                        //alert($("#"+event.target.id+"custommessagefield").val());
                        //var instructionmessage={"INSTRUCTION":{"USER":"<?php echo $_SESSION['user_name'];?>","minion":event.target.id,"JOB":"custommessage","SCRIPT":$("#"+event.target.id+"custommessagefield").val()}};
                        var minionid = event.target.id.split("/");

                        function csrfSafeMethod(method) {
                            return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
                        }

                        $.ajaxSetup({
                            beforeSend: function (xhr, settings) {
                                if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                                    xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
                                }
                            }
                        });
                        var url_mask = "{% url 'minIONcontrol-list' pk=12345 %}".replace(/12345/, minionid[minionid.length - 2].toString());
                        var stringthing = event.target.name;
                        var message = $("#" + stringthing + "custommessagefield").val();
                        $.ajax({
                            "type": "POST",
                            "dataType": "json",
                            "url": url_mask,
                            "data": {
                                "custom": message,
                                "job": "custommessage",
                                "minION": event.target.id,
                                "complete": 'False'
                            },
                            "beforeSend": function (xhr, settings) {
                                console.log("before send");
                                $.ajaxSettings.beforeSend(xhr, settings);
                            },
                            "success": function (result) {
                                console.log(result);
                            },
                            error: function (XMLHttpRequest, textStatus, errorThrown) {
                                console.log("Status: " + textStatus);
                                alert("Error: " + errorThrown);
                            }

                        })

                        //$('#startminionmodal').modal('hide');
                    },
                    stopminion: function (event) {
                        //$('#stopminionmodal').modal('hide');
                        var minionid = event.target.id.split("/");

                        function csrfSafeMethod(method) {
                            return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
                        }

                        $.ajaxSetup({
                            beforeSend: function (xhr, settings) {
                                if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                                    xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
                                }
                            }
                        });
                        var url_mask = "{% url 'minIONcontrol-list' pk=12345 %}".replace(/12345/, minionid[minionid.length - 2].toString());
                        $.ajax({
                            "type": "POST",
                            "dataType": "json",
                            "url": url_mask,
                            "data": {"custom": "", "job": "stopminion", "minION": event.target.id, "complete": 'False'},
                            "beforeSend": function (xhr, settings) {
                                console.log("before send");
                                $.ajaxSettings.beforeSend(xhr, settings);
                            },
                            "success": function (result) {
                                console.log(result);
                            },
                            error: function (XMLHttpRequest, textStatus, errorThrown) {
                                console.log("Status: " + textStatus);
                                alert("Error: " + errorThrown);
                            }

                        })
                    },
                    inactivateminion: function (event) {
                        //alert("hello");
                        //$('#stopminionmodal').modal('hide');
                        var minionid = event.target.id.split("/");

                        function csrfSafeMethod(method) {
                            return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
                        }

                        $.ajaxSetup({
                            beforeSend: function (xhr, settings) {
                                if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                                    xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
                                }
                            }
                        });
                        var url_mask = "{% url 'minIONcontrol-list' pk=12345 %}".replace(/12345/, minionid[minionid.length - 2].toString());
                        $.ajax({
                            "type": "POST",
                            "dataType": "json",
                            "url": url_mask,
                            "data": {
                                "custom": "",
                                "job": "shutdownminion",
                                "minION": event.target.id,
                                "complete": 'False'
                            },
                            "beforeSend": function (xhr, settings) {
                                console.log("before send");
                                $.ajaxSettings.beforeSend(xhr, settings);
                            },
                            "success": function (result) {
                                console.log(result);
                            },
                            error: function (XMLHttpRequest, textStatus, errorThrown) {
                                console.log("Status: " + textStatus);
                                alert("Error: " + errorThrown);
                            }

                        })

                    },
                    initminion: function (event) {
                        var minionid = event.target.id.split("/");

                        function csrfSafeMethod(method) {
                            return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
                        }

                        $.ajaxSetup({
                            beforeSend: function (xhr, settings) {
                                if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                                    xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
                                }
                            }
                        });
                        var url_mask = "{% url 'minIONcontrol-list' pk=12345 %}".replace(/12345/, minionid[minionid.length - 2].toString());
                        $.ajax({
                            "type": "POST",
                            "dataType": "json",
                            "url": url_mask,
                            "data": {
                                "custom": "",
                                "job": "initialiseminion",
                                "minION": event.target.id,
                                "complete": 'False'
                            },
                            "beforeSend": function (xhr, settings) {
                                console.log("before send");
                                $.ajaxSettings.beforeSend(xhr, settings);
                            },
                            "success": function (result) {
                                console.log(result);
                            },
                            error: function (XMLHttpRequest, textStatus, errorThrown) {
                                console.log("Status: " + textStatus);
                                alert("Error: " + errorThrown);
                            }

                        })

                    },

                }


            });

            {% verbatim %}

            Vue.component('diskusage', {
                template: "<div class='panel panel-info'><div class='panel-heading'><h3 class='panel-title'>Disk Space</h3></div><div class='panel-body'><div class='table-responsive'><table class='table table-condensed' ><tr><th>Category</td><th>Info</td></tr><tr><td>minKNOW computer name</td><td>{{title}}</td></tr><tr><td>Total Drive Capacity</td><td>{{capacity}}</td></tr><tr><td>Free Drive Space</td><td>{{space}} / {{percent}}%</td></tr><tr><td>Disk Space till Shutdown</td><td>{{bytealert}}</td></tr><tr><td>Warnings?</td><td>{{recalert}}</td></tr></table></div></div></div>",
                props: ['title', 'key', 'datain'],

                data: function () {
                    //var bytes_available = formatBytes(this.datain[0].bytes_available);
                    //console.log(this.datain[0].bytes_available);
                    var bytes_available = formatBytes(this.datain.bytes_available);

                    var drive_capacity = formatBytes(this.datain.bytes_capacity);
                    var percentage = (this.datain.bytes_available / this.datain.bytes_capacity * 100).toFixed(2);
                    //var percentage = "n/a";
                    var bytes_to_alert = formatBytes(this.datain.bytes_available - this.datain.bytes_to_alert);
                    var recommend_alert = this.datain.recommend_alert;
                    return {
                        space: bytes_available,
                        capacity: drive_capacity,
                        percent: percentage,
                        bytealert: bytes_to_alert,
                        recalert: recommend_alert
                    }

                },
                ready: function () {
                    this.$nextTick(function () {
                        //console.log(this.datain.bytes_available);
                        setInterval(function () {
                            //console.log(this.datain);
                            //this.space = formatBytes(this.datain.bytes_available);
                            this.capacity = formatBytes(this.datain.bytes_capacity);
                            //this.percent = (this.datain.bytes_available/this.datain.bytes_capacity * 100).toFixed(2);
                            this.bytealert = formatBytes(this.datain.bytes_available - this.datain.bytes_to_alert);
                            this.recalert = this.datain.recommend_alert;
                        }.bind(this), 15000);
                    });
                }
            })

            Vue.component('chartreadhist', {
                template: '<div id="container{{title}}" style="margin: 0 auto"</div>',
                props: ['title', 'key', 'datain', 'datain2', 'totalyield', 'seqspeed', 'readcount'],
                data: function () {
                    return {
                        opts: {
                            chart: {
                                renderTo: 'container' + this.title,
                                type: 'column',
                                zoomType: 'x'
                            },
                            title: {
                                text: 'Read length Histograms'
                            },
                            xAxis: {
                                categories: [],
                            },
                            yAxis: {
                                title: {
                                    text: 'Total Event Length'
                                }
                            },
                            credits: {
                                enabled: false
                            },
                            series: [{
                                name: 'Read Histogram',
                                //data: this.datain
                            }]
                        }
                    }
                }
                ,


                created: function () {
                },
                ready: function () {

                    this.$nextTick(function () {
                        this.chart = new Highcharts.Chart(this.opts);
                        //minion=this.key;

                        setInterval(function () {
                            var returndata = tohistogram(this.datain, parseInt(this.datain2), this.totalyield);
                            this.chart.series[0].setData(returndata[0]);
                            this.chart.xAxis[0].setCategories(returndata[1]);
                            console.log('updating highcharts histogram');
                            console.log(this.datain);
                            //console.log(returndata[2]-0.5);
                            //console.log(typeof(returndata[2]));
                            var N50 = parseInt(returndata[2]);
                            this.chart.xAxis[0].removePlotBand('plot-band-1');
                            this.chart.xAxis[0].addPlotBand({
                                from: N50 - 0.5,
                                to: N50 + 0.5,
                                color: '#FCFFC5',
                                id: 'plot-band-1',
                            });
                            this.chart.xAxis[0].removePlotBand('plot-band-2');
                            this.chart.xAxis[0].addPlotBand({
                                color: 'black',
                                width: 2,
                                dashStyle: 'longdashdot',
                                value: returndata[2],
                                label: {
                                    text: 'Estimated Read N50',
                                    align: 'left',
                                    rotation: 0,
                                    x: +10 // Amount of pixels the label will be repositioned according to the alignment.
                                },
                                id: 'plot-band-2',
                            });
                            this.chart.xAxis[0].removePlotBand('plot-band-3');
                            this.chart.xAxis[0].addPlotBand({
                                color: 'black',
                                width: 2,
                                dashStyle: 'longdashdot',
                                value: (Math.floor(this.totalyield / this.readcount / this.datain2)),
                                label: {
                                    text: 'Estimated Read Average - ' + Math.round(this.totalyield / this.readcount / 1000 * 100) / 100 + ' K events',
                                    align: 'left',
                                    rotation: 0,
                                    x: +10,
                                    y: +30, // Amount of pixels the label will be repositioned according to the alignment.
                                },
                                id: 'plot-band-3',
                            });
                        }.bind(this), 15000);
                    });
                }
            })

            Vue.component('container-curr-yield-prediction', {
                template: '<div id="container-curr-yield-prediction{{title}}" style="height: 140px; margin: 0 auto"</div>',
                props: ['title', 'key', 'script', 'eventcount'],
                data: function () {
                    return {

                        opts: {
                            chart: {
                                renderTo: 'container-curr-yield-prediction' + this.title,
                                type: 'solidgauge'
                            },

                            title: "Predicted Current Yield",

                            pane: {
                                center: ['50%', '85%'],
                                size: '100%',
                                startAngle: -90,
                                endAngle: 90,
                                background: {
                                    backgroundColor: (Highcharts.theme && Highcharts.theme.background2) || '#EEE',
                                    innerRadius: '60%',
                                    outerRadius: '100%',
                                    shape: 'arc'
                                }
                            },

                            tooltip: {
                                enabled: false
                            },

                            // the value axis
                            yAxis: {
                                type: 'logarithmic',
                                stops: [
                                    [0.3, '#0000FF'], // blue
                                    [0.37, '#DDDF0D'], // green
                                    [0.43, '#DF5353'], // red
                                ],
                                lineWidth: 0,
                                minorTickInterval: null,
                                tickPixelInterval: 400,
                                tickWidth: 0,
                                title: {
                                    y: -70
                                },
                                labels: {
                                    y: 16
                                },
                                min: 0.1,
                                max: 100000,
                                title: {
                                    text: null
                                }
                            },

                            plotOptions: {
                                solidgauge: {
                                    dataLabels: {
                                        y: -30,
                                        borderWidth: 0,
                                        useHTML: true
                                    }
                                }
                            },


                            credits: {
                                enabled: false
                            },

                            series: [{
                                name: 'Events',
                                data: [0],
                                dataLabels: {
                                    format: '<div style="text-align:center"><span style="font-size:15px;color:' +
                                    ((Highcharts.theme && Highcharts.theme.contrastTextColor) || 'black') + '">{y}</span><br/>' +
                                    '<span style="font-size:12px;color:silver"> Predicted Yield</span></div>'
                                },
                                tooltip: {
                                    valueSuffix: ' events'
                                }
                            }],

                            plotOptions: {
                                solidgauge: {
                                    dataLabels: {
                                        y: 30,
                                        borderWidth: 0,
                                        useHTML: true
                                    }
                                }
                            }
                        }
                    }
                }


                ,
                ready: function () {


                    this.$nextTick(function () {

                        switch (this.script) {
                            case "MegaCrazy Runs":
                                scaling = 3.5;
                                break;
                            case "450 b/s":
                                scaling = 1.8;
                                break;
                            case "250 b/s":
                                scaling = 1.1;
                                break;
                            case "70 b/s":
                                scaling = 1.0;
                                break;
                        }
                        //alert(scaling);
                        this.chart = new Highcharts.Chart(this.opts);
                        //this.chart.series[0].setData(this.datain);
                        if (this.chart) {
                            //point = this.chart.series[0].points[0];
                            this.chart.series[0].points[0].update(round(parseFloat((this.eventcount) * scaling), 0));
                            //point.update(this.datain);
                            //alert("camel");
                        }

                        setInterval(function () {
                            //    this.chart.series[0].setData(this.datain);
                            //    this.chart.redraw();
                            if (this.chart) {
                                switch (this.script) {
                                    case "MegaCrazy Runs":
                                        scaling = 3.5;
                                        break;
                                    case "450 b/s":
                                        scaling = 1.8;
                                        break;
                                    case "250 b/s":
                                        scaling = 1.1;
                                        break;
                                    case "70 b/s":
                                        scaling = 1.0;
                                        break;
                                }
                                point = this.chart.series[0].points[0];
                                point.update(round(parseFloat((this.eventcount) * scaling), 0));
                                //this.chart.redraw();


                            }
                        }.bind(this), 15000);
                    });
                }
            })
        })


        {% endverbatim %}

    </script>
{% endblock %}