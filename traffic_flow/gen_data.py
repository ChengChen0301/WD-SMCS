from trafficSimulator import *
import random


def initialize(v_max, a_max, b_max, T, s0, delta):
    sim = Simulation()
    sim.create_roads([
        ((-10, 106), (290, 106)),  # [0]
        ((-10, 102), (290, 102)),  # [1]
        ((290, 98), (-10, 98)),  # [2]
        ((290, 94), (-10, 94)),  # [3]
    ])
    sim.create_gen({
        'vehicle_rate': 180,
        'vehicles': [
            [3, {"path": [0]}, v_max, a_max, b_max, T, s0, delta],
            [6, {"path": [1]}, v_max, a_max, b_max, T, s0, delta],
            [3, {"path": [2]}, v_max, a_max, b_max, T, s0, delta],
            [6, {"path": [3]}, v_max, a_max, b_max, T, s0, delta],
        ]
    })
    sim.dt = 1.0 / 10
    return sim


if __name__ == '__main__':
    v_max = 8.33
    a_max = 1.44
    b_max = 4.61
    T = 1.6
    s0 = 2
    delta = 4

    sim = initialize(v_max, a_max, b_max, T, s0, delta)
    for i in range(100):
        sim.run(10)
        x1 = []
        x2 = []
        for j in range(len(sim.roads)):
            for k in range(len(sim.roads[j].vehicles)):
                x1.append((sim.roads[j].vehicles[k].x, sim.roads[j].start[1], sim.roads[j].vehicles[k].v, sim.roads[j].vehicles[k].alpha))
                noiseFlag = 1
                u = random.gauss(0, 0.1)
                x2.append((sim.roads[j].vehicles[k].x + u, sim.roads[j].start[1], sim.roads[j].vehicles[k].v, sim.roads[j].vehicles[k].alpha))
        np.savetxt('observation/vehicle_small_noise/truth' + str(i+1) + '.txt', x1)
        np.savetxt('observation/vehicle_large_noise/obsv' + str(i+1) + '.txt', x2)







