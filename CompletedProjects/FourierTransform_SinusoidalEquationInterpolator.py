import numpy
from matplotlib import pyplot
from numpy import fft

def series(x, a_hat):
    final = 0
    for k in range(0, len(a_hat)):
        phase = numpy.arctan2(numpy.imag(a_hat[k]), numpy.real(a_hat[k]))
        amplitude = numpy.absolute(a_hat[k])
        final += amplitude*numpy.cos(k*x+phase)

    return final

def seriesprinted(a_hat):
    final = ""
    for k in range(0, len(a_hat)):
        phase = numpy.arctan2(numpy.imag(a_hat[k]), numpy.real(a_hat[k]))
        amplitude = numpy.absolute(a_hat[k])
        final += "" + str(amplitude) + "*cos(" + str(k) + "*" + str(phase) + "*t) + "

    return final

xv = [646,646,646,646,653,666,686,708,734,760,788,817,860,886,902,915,922,923,923,917,907,897,877,856,837,824,811,802,797,791,784,779,772,766,761,755,749,740,730,719,708,696,680,671,660,655,649,646,643,643,643]
yv = [495,496,502,513,528,547,569,587,601,611,615,615,604,587,574,559,546,528,516,495,477,465,447,442,444,451,459,464,468,470,471,471,464,448,432,418,402,390,383,380,380,381,391,400,417,427,439,448,456,462,466]

xv_hat = fft.fft(xv)[:int(len(xv)/2)]
yv_hat = fft.fft(yv)[:int(len(yv)/2)]

x_vec = numpy.linspace(0, 6.28, 400)
x_out = numpy.zeros(len(x_vec))
for x in range(0, len(x_out)):
    x_out[x] = series(x_vec[x], xv_hat)

y_vec = numpy.linspace(0, 6.28, 400)
y_out = numpy.zeros(len(y_vec))
for y in range(0, len(y_out)):
    y_out[y] = series(y_vec[y], yv_hat)

print(seriesprinted(xv_hat))

def distance(x1, x2, y1, y2):
    return ((x2-x1)**2 + (y2-y1)**2)**(1/2)

dist = 0
for i in range(1, len(x_out)):
    pyplot.plot(x_out[i], y_out[i], 'bo')
    dist += distance(x_out[i], x_out[i-1], y_out[i], y_out[i-1])
    pyplot.pause(.001)
print(dist)
pyplot.show()




