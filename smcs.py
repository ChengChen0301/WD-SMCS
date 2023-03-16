from trafficSimulator import *
import random
import ot
import time
import os


def initialize(v_max, a_max, b_max, T, s0, delta):
    sim = Simulation()
    sim.create_roads([
        ((-10, 106), (290, 106)),  # [0]
        ((-10, 102), (290, 102)),  # [1]
        ((290, 98), (-10, 98)),  # [2]
        ((290, 94), (-10, 94)),  # [3]
    ])
    sim.create_gen({
        'vehicle_rate': 60,
        'vehicles': [
            [3, {"path": [0]}, v_max, a_max, b_max, T, s0, delta],
            [6, {"path": [1]}, v_max, a_max, b_max, T, s0, delta],
            [3, {"path": [2]}, v_max, a_max, b_max, T, s0, delta],
            [6, {"path": [3]}, v_max, a_max, b_max, T, s0, delta],
        ]
    })
    sim.dt = 1.0 / 10
    # data = np.loadtxt('observation/round2/obsv24.txt')
    # for i in range(len(data)):
    #     for j in range(len(sim.roads)):
    #         if data[i][1] == sim.roads[j].start[1]:
    #             config = {"path": [j]}
    #             theVehicle = Vehicle(config, v_max, a_max, b_max, T, s0)
    #             theVehicle.x = data[i, 0]
    #             theVehicle.v = data[i, 2]
    #             theVehicle.alpha = data[i, 3]
    #             sim.roads[j].vehicles.append(theVehicle)
    return sim


def disturbPara(para, a, b, level):
    flag = True
    while flag:
        u = np.random.normal(0, 1)
        temp = para + level * u
        if a <= temp <= b:
            flag = False
            para = temp
    return para


def compute_weights(d):
    return 1.0 / np.sqrt(2*np.pi) * np.exp(-(d*d)/2.0)


def resample(weights):
    pick = np.zeros(parNum)
    for j in range(parNum):
        u = random.uniform(0, 1)
        sumWeights = 0.0
        for i in range(parNum):
            sumWeights += weights[i]
            if sumWeights > u:
                pick[j] = i
                break
    return pick


def weighted_std(weights, values):
    means = np.sum(weights*values)
    return np.sqrt(np.sum(np.square(values - means)*weights))


def makedirs(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


if __name__ == '__main__':
    parNum = 500
    staDim = 6
    dataScale = 60
    Nth = 0.5 * parNum
    v_max = 16.6
    a_max = 1.44
    b_max = 4.61
    T = 1.6
    s0 = 2
    delta = 4

    a1 = 8.33
    b1 = 33.33
    a2 = 0.5
    b2 = 5
    a3 = 0.5
    b3 = 5
    a4 = 0.5
    b4 = 4
    a5 = 0.5
    b5 = 5
    # v_max = 16.6
    # a_max = 0.73
    # b_max = 1.67
    a6 = 1
    b6 = 5
    saveDir = 'result/round2_11'
    makedirs(saveDir)

    predSamp = np.zeros((parNum, staDim))
    nextSamp = np.zeros((parNum, staDim))
    save_each_wass = np.zeros((dataScale, parNum, 2*dataScale))
    save_ensemble = np.zeros((dataScale, parNum, 15))
    save_eff = np.zeros(dataScale)
    flag0 = np.zeros(parNum)
    accp = np.zeros(parNum)
    levels = np.zeros(staDim) + 1
    p1s = np.zeros(parNum)
    minwass = np.zeros((dataScale, 2))
    save_wass = np.zeros((dataScale, parNum, dataScale))

    predSamp[:, 0] = np.random.uniform(a1, b1, [parNum, 1]).squeeze()  # v
    predSamp[:, 1] = np.random.uniform(a2, b2, [parNum, 1]).squeeze()  # a
    # predSamp[:, 2] = np.random.uniform(a3, b3, [parNum, 1]).squeeze()  # b
    predSamp[:, 2] = np.zeros(parNum) + b_max
    predSamp[:, 3] = np.random.uniform(a4, b4, [parNum, 1]).squeeze()  # T
    # predSamp[:, 4] = np.random.uniform(a5, b5, [parNum, 1]).squeeze()  # s0
    predSamp[:, 4] = np.zeros(parNum) + s0
    # predSamp[:, 5] = np.random.uniform(a6, b6, [parNum, 1]).squeeze()  # delta
    predSamp[:, 5] = np.zeros(parNum) + delta

    # predSamp[0, 0] = v_max
    # predSamp[0, 1] = a_max
    # predSamp[0, 3] = T
    # predSamp[0, 5] = delta

    weights = np.zeros(parNum) + 1.0 / parNum
    save_ensemble[0, :, 0] = predSamp[:, 0]
    save_ensemble[0, :, 1] = predSamp[:, 1]
    save_ensemble[0, :, 2] = predSamp[:, 2]
    save_ensemble[0, :, 3] = predSamp[:, 3]
    save_ensemble[0, :, 4] = predSamp[:, 4]
    save_ensemble[0, :, 5] = predSamp[:, 5]
    save_ensemble[0, :, 6] = weights
    save_eff[0] = parNum
    levels[0] = np.std(predSamp[:, 0])
    levels[1] = np.std(predSamp[:, 1])
    levels[2] = np.std(predSamp[:, 2])
    levels[3] = np.std(predSamp[:, 3])
    levels[4] = np.std(predSamp[:, 4])
    levels[5] = np.std(predSamp[:, 5])
    print(levels)

    t = 1
    t0 = time.time()
    while t < dataScale:
        t1 = time.time()
        minDist = np.zeros(dataScale) + 100.0
        maxDist = np.zeros(dataScale)
        for i in range(parNum):
            nextSamp[i, 0] = disturbPara(predSamp[i, 0], a1, b1, levels[0])
            nextSamp[i, 1] = disturbPara(predSamp[i, 1], a2, b2, levels[1])
            # nextSamp[i, 2] = disturbPara(predSamp[i, 2], a3, b3, levels[2])
            nextSamp[i, 2] = predSamp[i, 2]
            nextSamp[i, 3] = disturbPara(predSamp[i, 3], a4, b4, levels[3])
            # nextSamp[i, 4] = disturbPara(predSamp[i, 4], a5, b5, levels[4])
            nextSamp[i, 4] = predSamp[i, 4]
            # nextSamp[i, 5] = disturbPara(predSamp[i, 5], a6, b6, levels[5])
            nextSamp[i, 5] = predSamp[i, 5]

            for tt in range(1, t+1, 1):
                x0 = np.loadtxt('observation/round2/obsv' + str(tt+0) + '.txt').reshape(-1, 4)[:, 0:2]
                w0 = np.ones((len(x0),)) / len(x0)

                if tt < t:
                    if flag0[i] == 0:
                        emd1 = save_each_wass[t-1][i][tt]
                    else:
                        emd1 = save_each_wass[t-1][i][tt+dataScale]
                else:
                    sim = initialize(predSamp[i, 0], predSamp[i, 1], predSamp[i, 2], predSamp[i, 3], predSamp[i, 4], predSamp[i, 5])
                    steps = int(1.0 / sim.dt)
                    sim.run(steps * tt)
                    x1 = []
                    for j in range(len(sim.roads)):
                        for k in range(len(sim.roads[j].vehicles)):
                            x1.append((sim.roads[j].vehicles[k].x, sim.roads[j].start[1]))
                    x1 = np.array(x1)
                    np.savetxt(saveDir + '/' + str(t) + '_' + str(i) + '_' + str(tt) + '.txt', x1)
                    # print(x1)
                    w1 = np.ones((len(x1),)) / len(x1)
                    M = np.sqrt(ot.dist(x0, x1))
                    emd1 = ot.emd2(w0, w1, M, numItermax=1e7)
                    if emd1 < np.exp(-20):
                        emd1 = np.exp(-20)
                save_each_wass[t, i, tt] = emd1
                if emd1 < minDist[tt]:
                    if emd1 > 0:
                        minDist[tt] = emd1
                        if tt == t:
                            minwass[t, 0] = i
                            minwass[t, 1] = emd1

                sim = initialize(nextSamp[i, 0], nextSamp[i, 1], nextSamp[i, 2], nextSamp[i, 3], nextSamp[i, 4], nextSamp[i, 5])
                steps = int(1.0 / sim.dt)
                sim.run(steps * tt)
                x2 = []
                for j in range(len(sim.roads)):
                    for k in range(len(sim.roads[j].vehicles)):
                        x2.append((sim.roads[j].vehicles[k].x, sim.roads[j].start[1]))
                x2 = np.array(x2)
                # np.savetxt('test/' + str(t) + '_' + str(i) + '_' + str(tt+dataScale) + '.txt', x2)
                # print(x2)
                w2 = np.ones((len(x2),)) / len(x2)
                M = np.sqrt(ot.dist(x0, x2))
                emd2 = ot.emd2(w0, w2, M, numItermax=1e7)
                if emd2 < np.exp(-20):
                    emd2 = np.exp(-20)
                save_each_wass[t, i, tt+dataScale] = emd2
                if emd2 < minDist[tt]:
                    if emd2 > 0:
                        minDist[tt] = emd2

        maxDist = np.max(save_each_wass[t], axis=0)

        for i in range(parNum):
            p1 = 1.0
            p2 = 1.0
            for tt in range(1, t+1, 1):
                # p1 = p1 * compute_weights(save_each_wass[t, i, tt] / minDist[tt])
                # p2 = p2 * compute_weights(save_each_wass[t, i, tt+dataScale] / minDist[tt])
                p1 = p1 * compute_weights(save_each_wass[t, i, tt] / maxDist[tt])
                p2 = p2 * compute_weights(save_each_wass[t, i, tt + dataScale] / maxDist[tt])

            if p1 >= p2:
                if p1 == 0:
                    accptance = 1
                    print(str(i) + 'yes')
                else:
                    accptance = p2 / p1
            # if p1 > p2:
            #     accptance = p2 / p1
            else:
                accptance = 1.0
            p1s[i] = p1

            u = random.uniform(0, 1)
            if u < accptance:
                predSamp[i, 0] = nextSamp[i, 0]
                predSamp[i, 1] = nextSamp[i, 1]
                predSamp[i, 2] = nextSamp[i, 2]
                predSamp[i, 3] = nextSamp[i, 3]
                predSamp[i, 4] = nextSamp[i, 4]
                predSamp[i, 5] = nextSamp[i, 5]
                flag0[i] = 1
            else:
                flag0[i] = 0

            # weights[i] = weights[i] * compute_weights(save_each_wass[t, i, t] / minDist[t])
            weights[i] = weights[i] * compute_weights(save_each_wass[t, i, t] / maxDist[t])
            accp[i] = accptance

        if np.sum(weights) == 0:
            print("yes")
            weights = np.zeros(parNum) + 1.0 / parNum
        else:
            weights = weights / np.sum(weights)

        save_ensemble[t, :, 0] = predSamp[:, 0]
        save_ensemble[t, :, 1] = predSamp[:, 1]
        save_ensemble[t, :, 2] = predSamp[:, 2]
        save_ensemble[t, :, 3] = predSamp[:, 3]
        save_ensemble[t, :, 4] = predSamp[:, 4]
        save_ensemble[t, :, 5] = predSamp[:, 5]
        save_ensemble[t, :, 6] = weights
        save_ensemble[t, :, 7] = accp
        save_ensemble[t, :, 8] = p1s
        save_ensemble[t, :, 9] = nextSamp[:, 0]
        save_ensemble[t, :, 10] = nextSamp[:, 1]
        save_ensemble[t, :, 11] = nextSamp[:, 2]
        save_ensemble[t, :, 12] = nextSamp[:, 3]
        save_ensemble[t, :, 13] = nextSamp[:, 4]
        save_ensemble[t, :, 14] = nextSamp[:, 5]
        save_eff[t] = 1.0 / np.sum(np.square(weights))
        save_wass[t] = save_each_wass[t, :, 0:dataScale]

        # levels[0] = weighted_std(weights, predSamp[:, 0])
        # levels[1] = weighted_std(weights, predSamp[:, 1])
        # levels[2] = weighted_std(weights, predSamp[:, 2])
        # levels[3] = weighted_std(weights, predSamp[:, 3])
        # levels[4] = weighted_std(weights, predSamp[:, 4])
        # levels[5] = weighted_std(weights, predSamp[:, 5])

        if save_eff[t] < Nth:
            pick = resample(weights)
            newSamp = predSamp[pick.astype(int), :]
            predSamp = newSamp
            weights = np.zeros(parNum) + 1.0 / parNum
            newWass = save_each_wass[t, pick.astype(int), :]
            save_each_wass[t] = newWass
            newflag = flag0[pick.astype(int)]
            flag0 = newflag
        t2 = time.time()

        print("t = " + str(t) + ", Neff = " + str(save_eff[t]) + ", time = " + str(t2 - t1))
        t += 1
        np.save(saveDir + '/ensembles.npy', save_ensemble)
        np.save(saveDir+'/wass.npy', save_wass)
        np.save(saveDir + '/eff.npy', save_eff)
        np.save(saveDir + '/minwass.npy', minwass)

    print("total_time = " + str(t2 - t0))













