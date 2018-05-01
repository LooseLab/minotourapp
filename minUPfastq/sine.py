import math
import sys
import threading
import time
from random import randint

import numpy

import minknow.data_reader as data_reader
import pyaudio
from minknow.engine_client import EngineClient

_e = EngineClient()

def sine(frequency, length, rate):
    length = int(length * rate)
    factor = float(frequency) * (math.pi * 2) / rate
    return numpy.sin(numpy.arange(length) * factor)


def play_tone(stream, frequency=440, length=1, rate=128000):
    chunks = []
    chunks.append(sine(frequency*4, length, rate))
    chunk = numpy.concatenate(chunks) * 0.25
    #print chunk
    asign = numpy.sign(chunk)
    signchange = ((numpy.roll(asign, 1) - asign) != 0).astype(int)
    #print signchange
    occurences = numpy.where(signchange == signchange.max())
    #print occurences[-1][-1]
    #for index, item in enumerate(chunks):
    #    if (index >= occurences[-1][-1]):
    #        items[index] = 0
    chunk[occurences[-1][-1]:] = 0
    #print chunk
    stream.write(chunk.astype(numpy.float32).tostring())


class ThreadingExample(object):
    """ Threading example class
    The run() method will be started and it will run in the background
    until the application exits.
    """

    def __init__(self, interval=0):
        """ Constructor
        :type interval: int
        :param interval: Check interval, in seconds
        """
        self.interval = interval
	self.datadict = list()
	self.p = pyaudio.PyAudio()
    	self.stream = self.p.open(format=pyaudio.paFloat32,
                    channels=1, rate=44100, output=1)
        thread = threading.Thread(target=self.run, args=())
        thread.daemon = True                            # Daemonize thread
        thread.start()                                  # Start the execution


    def add_freqs(self,chunks):
	self.datadict.append(chunks)

    def run(self):
        """ Method that runs forever """
        while True:
            # Do something
            #print('Doing something imporant in the background')
	    if len(self.datadict) > 0:
	        for bit in self.datadict.pop(0):
		    play_tone(self.stream,frequency=bit,length=0.01)

            time.sleep(self.interval)


if __name__ == '__main__':
    #p = pyaudio.PyAudio()
    #stream = p.open(format=pyaudio.paFloat32,
    #               channels=1, rate=44100, output=1)
    example = ThreadingExample()
    while True:
        _data = data_reader.latest_raw_data(_e, 800)
	print sys.argv[0]
	chunks = numpy.mean(_data.data[int(sys.argv[1])].reshape(-1, 5), axis=1)
	print len(chunks)
	example.add_freqs(chunks)
        #for bit in chunks:
	#    play_tone(stream,frequency=bit,length=0.0025)
    #for x in range (1000):	
    #	play_tone(stream,frequency=randint(110,660),length=0.0025)
    	#play_tone(stream,frequency=440,length=0.01)
    #for x in range (1000):
    #	play_tone(stream,frequency=randint(660,1000),length=0.0025)

    #stream.close()
    #p.terminate()
