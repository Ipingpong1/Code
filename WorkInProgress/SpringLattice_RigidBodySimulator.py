import numpy
from matplotlib import pyplot

class node:
    def __init__(self, x, y, k):
        self.xpos = x
        self.ypos = y
        self.k = k
    connections = numpy.array([])
    connection_distances_0 = numpy.array([])
    yvel = 0
    xvel = 0
    yvel_0 = 0
    xvel_0 = 0

def distance(node, onode):
    return ((onode.xpos - node.xpos)**2 + (onode.ypos - node.ypos)**2) ** (1/2)

dt = .009
g = 9.8

x0 = 50
y0 = 20
size = 25
k = 50

Lbo = node(0 + x0, 0 + y0, k)
Rbo = node(size + x0, 0 + y0 + 10, k)
Lup = node(0 + x0 - 10, size + y0, k)
Rup = node(size + x0, size + y0 + 10, k)

Lbo.connections = numpy.array([Lup, Rbo])
Rbo.connections = numpy.array([Lbo, Rup])
Lup.connections = numpy.array([Lbo, Rup])
Rup.connections = numpy.array([Lup, Rbo])

nodes = numpy.array([Lbo, Rbo, Lup, Rup])

for node in nodes:
    for i in range(0, len(node.connections)):
        node.connection_distances_0 = numpy.append(node.connection_distances_0, distance(node, node.connections[i]))

edges = []
for i in range(0, len(nodes)):
    for edge in nodes[i].connections:
        edges.append([nodes[i], edge])

for t in range(0, 200):
    for node in nodes:
        # update velocities
        if(node.ypos <= 0 or node.ypos + node.yvel <= 0):
            node.yvel *= -1
            node.xvel *= 1
            node.ypos = 0
        else:
            node.yvel -= g*dt
        node.yvel *= .96
        node.xvel *= .96
        # apply displacements
        node.xpos += node.xvel
        node.ypos += node.yvel

    nodesn = numpy.copy(nodes)
    for node in nodes:
        for i, other_node in enumerate(nodes):
            if(other_node in node.connections):
                dist = distance(node, other_node)
                dist_0 = node.connection_distances_0[numpy.where(node.connections == other_node)]
                if(dist > dist_0):
                    angle = numpy.degrees(numpy.arcsin((other_node.ypos - node.ypos)/dist))
                    x = dist - dist_0
                    other_node.xvel -= k*x*numpy.cos(numpy.radians(angle))*dt
                    other_node.yvel -= k*x*numpy.sin(numpy.radians(angle))*dt
                elif (dist < dist_0):
                    angle = numpy.degrees(numpy.arcsin((other_node.ypos - node.ypos) / dist))
                    x = dist_0 - dist
                    other_node.xvel += k*x*numpy.cos(numpy.radians(angle))*dt
                    other_node.yvel += k*x*numpy.sin(numpy.radians(angle))*dt
    nodes = nodesn

    pyplot.xlim(0, 150)
    pyplot.ylim(-5, 150)
    for j in range(0, len(edges)):
        pyplot.plot([edges[j][0].xpos, edges[j][1].xpos], [edges[j][0].ypos, edges[j][1].ypos], marker='o', linestyle='-')
    pyplot.pause(.01)
    pyplot.close()
    pyplot.show()












