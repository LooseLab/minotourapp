import math

import numpy

import pyaudio


def sine(frequency, length, rate):
    length = int(length * rate)
    factor = (float(frequency) * (math.pi * 2) / rate)
    return numpy.sin(numpy.arange(length) * factor)


def play_tone(stream, frequency, length, rate=44100):
    chunks = [sine(frequency, length, rate)]

    chunk = numpy.concatenate(chunks) * 0.25

    fade = 200.

    fade_in = numpy.arange(0., 1., 1/fade)
    fade_out = numpy.arange(1., 0., -1/fade)

    chunk[:fade] = numpy.multiply(chunk[:fade], fade_in)
    chunk[-fade:] = numpy.multiply(chunk[-fade:], fade_out)

    stream.write(chunk.astype(numpy.float32).tostring())


def test():
    test_freqs = [50, 100, 200, 400, 800, 1200, 2000, 3200]

    for i in range(2):
        for freq in test_freqs:
            play_tone(stream, freq, 1)


if __name__ == '__main__':
    p = pyaudio.PyAudio()
    stream = p.open(format=pyaudio.paFloat32,
                    channels=1, rate=44100, output=1)


test()
