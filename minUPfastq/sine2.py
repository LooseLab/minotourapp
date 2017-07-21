import numpy as np
import sounddevice as sd
import math

def playArray(pts, time):

    samplerate = 44100.0
    sd.default.samplerate = samplerate
    volume = 10000
    total = len(pts)
    count = 0
    phase = 0.0
    samples = np.zeros(math.ceil(samplerate * time  * total), dtype = np.float)

    for i, val in enumerate(samples):
        freq = pts[count]
        freqRads = 2 * np.pi * freq / samplerate
        phase = phase + freqRads
        sampleValue = volume * np.sin(phase)
        samples[i] = sampleValue
        if ( i >  0 and i % (samplerate * time) == 0):
            count = count + 1

    wav_wave = np.array(samples, dtype=np.int16)
    sd.play(wav_wave, blocking=True)
