import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import ot
import ot.plot
import scipy.stats as st
import matplotlib.patches as mpatches


Nt = 30
v_max = 8.33
a_max = 1.44
b_max = 4.61
T = 1.6
s0 = 2
delta = 4

a = 2000.0
b = 0.08
v = 1.0

#------------------------------------draw filter result--------------------------------------/


def weighted_std(weights, values):
    means = np.sum(weights*values)
    return means, np.sqrt(np.sum(np.square(values - means)*weights))


def draw_idm(i, j):
    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    ax1 = plt.subplot(gs[0], frameon=True)
    ax2 = plt.subplot(gs[1], frameon=True)
    ax3 = plt.subplot(gs[2], frameon=True)
    ax4 = plt.subplot(gs[3], frameon=True)

    dir3 = 'result/round' + str(i) + '_' + str(j)

    pred = np.zeros((Nt, 4))
    t = np.linspace(0, Nt-1, Nt)
    ensembles = np.load(dir3 + '/ensembles.npy')

    for i in range(Nt):
        ensemble = ensembles[i]
        pred[i, 0], std = weighted_std(ensemble[:, 6], ensemble[:, 0])
        ax1.errorbar(i, pred[i, 0], yerr=std, mfc='r', fmt='-', ecolor='black')
        pred[i, 1], std = weighted_std(ensemble[:, 6], ensemble[:, 1])
        ax2.errorbar(i, pred[i, 1], yerr=std, mfc='r', fmt='-', ecolor='black')
        pred[i, 2], std = weighted_std(ensemble[:, 6], ensemble[:, 3])
        ax3.errorbar(i, pred[i, 2], yerr=std, mfc='r', fmt='-', ecolor='black')
        pred[i, 3], std = weighted_std(ensemble[:, 6], ensemble[:, 5])
        ax4.errorbar(i, pred[i, 3], yerr=std, mfc='r', fmt='-', ecolor='black')

    ax1.plot(t, pred[:, 0])
    ax2.plot(t, pred[:, 1])
    ax3.plot(t, pred[:, 2])
    ax4.plot(t, pred[:, 3])
    ax1.plot(t, np.zeros(Nt)+v_max)
    ax2.plot(t, np.zeros(Nt)+a_max)
    ax3.plot(t, np.zeros(Nt)+T)
    ax4.plot(t, np.zeros(Nt)+delta)
    ax1.set_ylabel('$v_0$')
    ax2.set_ylabel('a')
    ax3.set_ylabel('T')
    ax4.set_ylabel('$\delta$')
    fig.text(0.5, 0, 't(s)', ha='center')
    plt.savefig(dir3 + '/idm.png', bbox_inches='tight', pad_inches=0.1)


def draw_idm2(i, j):
    fig = plt.figure(figsize=(21, 4))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
    ax1 = plt.subplot(gs[0], frameon=True)
    ax2 = plt.subplot(gs[1], frameon=True)
    ax3 = plt.subplot(gs[2], frameon=True)
    fs = 15

    dir3 = 'result/round' + str(i) + '_' + str(j)
    pred = np.zeros((Nt, 4))
    t = np.linspace(0, Nt - 1, Nt)
    ensembles = np.load(dir3 + '/ensembles.npy')

    for i in range(Nt):
        ensemble = ensembles[i]
        pred[i, 0], std = weighted_std(ensemble[:, 6], ensemble[:, 0])
        ax1.errorbar(i, pred[i, 0], yerr=std, mfc='r', fmt='-', ecolor='black')
        pred[i, 1], std = weighted_std(ensemble[:, 6], ensemble[:, 1])
        ax2.errorbar(i, pred[i, 1], yerr=std, mfc='r', fmt='-', ecolor='black')
        pred[i, 2], std = weighted_std(ensemble[:, 6], ensemble[:, 3])
        ax3.errorbar(i, pred[i, 2], yerr=std, mfc='r', fmt='-', ecolor='black')

    ax1.plot(t, pred[:, 0])
    ax2.plot(t, pred[:, 1], label='estimation')
    ax3.plot(t, pred[:, 2])
    ax1.plot(t, np.zeros(Nt) + v_max)
    ax2.plot(t, np.zeros(Nt) + a_max, label='truth')
    ax3.plot(t, np.zeros(Nt) + T)

    ax1.set_ylabel('$v_0$', fontsize=fs)
    ax2.set_ylabel('a', fontsize=fs)
    ax3.set_ylabel('$T_s$', fontsize=fs)
    ax1.tick_params(labelsize=fs)
    ax2.tick_params(labelsize=fs)
    ax3.tick_params(labelsize=fs)
    ax2.legend(fontsize=fs)
    fig.text(0.5, -0.05, 't', ha='center', fontsize=fs)
    # plt.savefig(dir3 + '/idm.eps', format='eps', bbox_inches='tight', pad_inches=0.1)
    plt.savefig(dir3 + '/idm.png', bbox_inches='tight', pad_inches=0.1)


def draw_sfm(i, j):
    fig = plt.figure(figsize=(21, 4))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
    ax1 = plt.subplot(gs[0], frameon=True)
    ax2 = plt.subplot(gs[1], frameon=True)
    ax3 = plt.subplot(gs[2], frameon=True)
    fs = 15

    dir3 = '/Users/cc/CLionProjects/PF/result/round' + str(i) + '_' + str(j)
    pred = np.zeros((Nt, 3))
    t = np.linspace(0, Nt-1, Nt)
    ensembles = np.loadtxt(dir3 + '/ensembles.txt')
    ensembles = ensembles.reshape(Nt, 500, -1)
    for i in range(Nt):
        ensemble = ensembles[i]
        pred[i, 0], std = weighted_std(ensemble[:, 3], ensemble[:, 0])
        # ax1.errorbar(i, pred[i, 0], yerr=std, mfc='r', fmt='-', ecolor='black')
        pred[i, 1], std = weighted_std(ensemble[:, 3], ensemble[:, 1])
        # ax2.errorbar(i, pred[i, 1], yerr=std, mfc='r', fmt='-', ecolor='black')
        pred[i, 2], std = weighted_std(ensemble[:, 3], ensemble[:, 2])
        # ax3.errorbar(i, pred[i, 2], yerr=std, mfc='r', fmt='-', ecolor='black')

    ax1.plot(t, pred[:, 0])
    ax2.plot(t, pred[:, 1], label='estimation')
    ax3.plot(t, pred[:, 2])
    ax1.plot(t, np.zeros(Nt) + a)
    ax2.plot(t, np.zeros(Nt) + b, label='truth')
    ax3.plot(t, np.zeros(Nt) + v)
    ax1.set_ylabel('A', fontsize=fs)
    ax2.set_ylabel('B', fontsize=fs)
    ax3.set_ylabel('$v^p$', fontsize=fs)
    ax1.tick_params(labelsize=fs)
    ax2.tick_params(labelsize=fs)
    ax3.tick_params(labelsize=fs)
    ax2.legend(fontsize=fs)
    fig.text(0.5, -0.05, 't', ha='center', fontsize=fs)
    plt.savefig(dir3 + '/sfm.png', bbox_inches='tight', pad_inches=0.1)


def draw_sfm2(i, j):
    fig = plt.figure(figsize=(15, 3))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
    ax1 = plt.subplot(gs[0], frameon=True)
    ax2 = plt.subplot(gs[1], frameon=True)
    ax3 = plt.subplot(gs[2], frameon=True)

    dir3 = '/Users/cc/CLionProjects/PF/result/round' + str(i) + '_' + str(j)
    pred = np.zeros((Nt, 3))
    t = np.linspace(0, Nt-1, Nt)
    ensembles = np.loadtxt(dir3 + '/resamples.txt')
    ensembles = ensembles.reshape(Nt, 500, -1)
    for i in range(Nt):
        ensemble = ensembles[i]
        pred[i, 0], std = weighted_std(np.zeros(500)+1.0/500, ensemble[:, 0])
        ax1.errorbar(i, pred[i, 0], yerr=std, mfc='r', fmt='-', ecolor='black')
        pred[i, 1], std = weighted_std(np.zeros(500)+1.0/500, ensemble[:, 1])
        ax2.errorbar(i, pred[i, 1], yerr=std, mfc='r', fmt='-', ecolor='black')
        pred[i, 2], std = weighted_std(np.zeros(500)+1.0/500, ensemble[:, 2])
        ax3.errorbar(i, pred[i, 2], yerr=std, mfc='r', fmt='-', ecolor='black')

    ax1.plot(t, pred[:, 0])
    ax2.plot(t, pred[:, 1])
    ax3.plot(t, pred[:, 2])
    ax1.plot(t, np.zeros(Nt) + a)
    ax2.plot(t, np.zeros(Nt) + b)
    ax3.plot(t, np.zeros(Nt) + v)
    ax1.set_ylabel('A')
    ax2.set_ylabel('B')
    ax3.set_ylabel('$v^p$')
    fig.text(0.5, 0, 't', ha='center')
    plt.savefig(saveDir + '/sfm2.png', bbox_inches='tight', pad_inches=0.1)


def draw_wass(i, j):
    fig = plt.figure()
    dir3 = '/Users/cc/CLionProjects/PF/result/round' + str(i) + '_' + str(j)
    wass2 = np.loadtxt(dir3 + '/savewass.txt')
    emd = np.std(wass, axis=1)
    eff = np.loadtxt(dir3 + '/eff.txt')
    t = np.linspace(0, Nt - 1, Nt)
    plt.plot(t, emd)
    plt.plot(t, eff/10000)
    plt.xlabel('t')
    plt.ylabel('emd')
    plt.savefig(dir3 + '/eff.png', bbox_inches='tight', pad_inches=0.1)


test_point = np.zeros((100, 2))
for i in range(100):
    index_i = int(i / 10)
    index_j = int(i % 10)
    test_point[i][0] = -5.0 + 0.5 * (2 * index_i + 1) * 10.0 / 10
    test_point[i][1] = -5.0 + 0.5 * (2 * index_j + 1) * 10.0 / 10


def compute_emd(array1, array2):
    # 画格子计算
    count1 = np.zeros(100)
    count2 = np.zeros(100)
    for j in range(100):
        index_i = int((array1[j][0] + 5.0) / 1)
        index_j = int((array1[j][1] + 5.0) / 1)
        if index_i < 10:
            count1[index_i * 10 + index_j] += 1
        index_i = int((array2[j][0] + 5.0) / 1)
        index_j = int((array2[j][1] + 5.0) / 1)
        if index_i < 10:
            count2[index_i * 10 + index_j] += 1

    index1 = np.where(count1 > 0)[0]
    index2 = np.where(count2 > 0)[0]
    xs = test_point[index1]
    xt = test_point[index2]
    a = count1[index1] / np.sum(count1)
    b = count2[index2] / np.sum(count2)

    M = np.sqrt(ot.dist(xs, xt))
    emd = ot.emd2(a, b, M, numItermax=1e7)
    # print(emd)
    return emd


def compute_emd2(array1, array2):
    # 直接根据坐标计算
    index1 = np.where(array1[:, 0] < 5.0)
    index2 = np.where(array2[:, 0] < 5.0)

    xs = array1[index1]
    xt = array2[index2]
    a = np.zeros(len(xs)) + 1.0 / len(xs)
    b = np.zeros(len(xt)) + 1.0 / len(xt)

    M = np.sqrt(ot.dist(xs, xt))
    emd = ot.emd2(a, b, M, numItermax=1e7)
    # print(emd)
    return emd


def draw_kde2D(i, j, k):
    xmin, xmax = -5.0, 5.0
    ymin, ymax = -5.0, 5.0
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])

    dir3 = '/Users/cc/CLionProjects/PF/result/round' + str(i) + '_' + str(j)
    obsv = np.loadtxt('/Users/cc/CLionProjects/PF/observation/round' + str(i) + '/obsv' + str(k * 100) + '.txt')
    x = obsv[:, 0]
    y = obsv[:, 1]
    index = np.where(x < 5.0)[0]
    x = x[index]
    y = y[index]
    values = np.vstack([x, y])
    kernel_1 = st.gaussian_kde(values)
    f1 = np.reshape(kernel_1(positions).T, xx.shape)

    pos = np.loadtxt(dir3 + '/positions.txt')
    pos = pos.reshape(Nt, 100, -1)
    pred = pos[k]
    x = pred[:, 0]
    y = pred[:, 1]
    index = np.where(x < 5.0)[0]
    x = x[index]
    y = y[index]
    values = np.vstack([x, y])
    kernel_2 = st.gaussian_kde(values)
    f2 = np.reshape(kernel_2(positions).T, xx.shape)

    wass = np.loadtxt(dir3 + '/savewass.txt')
    emd = np.min(wass, axis=1)

    fig = plt.figure(figsize=(5, 8))
    ax1 = fig.add_subplot(211, frameon=False)
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_title('obsv')
    cfset = ax1.contourf(xx, yy, f1)
    plt.axis('equal')

    ax2 = fig.add_subplot(212, frameon=False)
    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title('pred')
    plt.axis('equal')
    cfset = ax2.contourf(xx, yy, f2)

    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8)
    cb_ax = fig.add_axes([0.9, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(cfset, cb_ax)
    plt.title('t = ' + str(k) + ', emd = ' + str(emd[k]))
    plt.savefig(dir3 + '/density/' + str(k) + '.png')


def compute_kde2D(i, j, k):
    # 对应compute_emd
    xmin, xmax = -5.0, 5.0
    ymin, ymax = -5.0, 5.0
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])

    dir3 = '/Users/cc/CLionProjects/PF/result/round' + str(i) + '_' + str(j)
    obsv = np.loadtxt('/Users/cc/CLionProjects/PF/observation/round' + str(i) + '/obsv' + str(k * 100) + '.txt')
    x = obsv[:, 0]
    y = obsv[:, 1]
    index = np.where(x < 5.0)[0]
    x = x[index]
    y = y[index]
    values = np.vstack([x, y])
    kernel_1 = st.gaussian_kde(values)
    f1 = np.reshape(kernel_1(positions).T, xx.shape)

    pos = np.loadtxt(dir3 + '/positions.txt')
    pos = pos.reshape(Nt, 100, -1)
    pred = pos[k]
    x = pred[:, 0]
    y = pred[:, 1]
    index = np.where(x < 5.0)[0]
    x = x[index]
    y = y[index]
    values = np.vstack([x, y])
    kernel_2 = st.gaussian_kde(values)
    f2 = np.reshape(kernel_2(positions).T, xx.shape)

    # truth = np.loadtxt('/Users/cc/CLionProjects/PF/observation/round' + str(i) + '/truth' + str(k * 100) + '.txt')
    # x = truth[:, 0]
    # y = truth[:, 1]
    # index = np.where(x < 5.0)[0]
    # x = x[index]
    # y = y[index]
    # values = np.vstack([x, y])
    # kernel_3 = st.gaussian_kde(values)
    # f3 = np.reshape(kernel_3(positions).T, xx.shape)

    prior = np.loadtxt('/Users/cc/CLionProjects/PF/prior/round' + str(i) + '/truth' + str(k * 100) + '.txt')
    x = prior[:, 0]
    y = prior[:, 1]
    index = np.where(x < 5.0)[0]
    x = x[index]
    y = y[index]
    values = np.vstack([x, y])
    kernel_3 = st.gaussian_kde(values)
    f3 = np.reshape(kernel_3(positions).T, xx.shape)

    emd0 = compute_emd(prior, obsv)
    emd1 = compute_emd(obsv, pred)

    return f1, f2, f3, emd0, emd1


def compute_kde2D2(i, j, k):
    # 对应compute_emd2
    xmin, xmax = -5.0, 5.0
    ymin, ymax = -5.0, 5.0
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])

    dir3 = '/Users/cc/CLionProjects/PF/result/round' + str(i) + '_' + str(j)
    obsv = np.loadtxt('/Users/cc/CLionProjects/PF/observation/round' + str(i) + '/obsv' + str(k * 100) + '.txt')
    x = obsv[:, 0]
    y = obsv[:, 1]
    index = np.where(x < 5.0)[0]
    x = x[index]
    y = y[index]
    values = np.vstack([x, y])
    kernel_1 = st.gaussian_kde(values)
    f1 = np.reshape(kernel_1(positions).T, xx.shape)

    # pos = np.loadtxt(dir3 + '/positions.txt')
    # pos = pos.reshape(Nt, 100, -1)
    # pred = pos[k]
    pred = np.loadtxt('/Users/cc/CLionProjects/PF/posterior/round' + str(i) + '_' + str(k) + '/truth' + str(k * 100) + '.txt')
    x = pred[:, 0]
    y = pred[:, 1]
    index = np.where(x < 5.0)[0]
    x = x[index]
    y = y[index]
    values = np.vstack([x, y])
    kernel_2 = st.gaussian_kde(values)
    f2 = np.reshape(kernel_2(positions).T, xx.shape)

    # truth = np.loadtxt('/Users/cc/CLionProjects/PF/observation/round' + str(i) + '/truth' + str(k * 100) + '.txt')
    # x = truth[:, 0]
    # y = truth[:, 1]
    # index = np.where(x < 5.0)[0]
    # x = x[index]
    # y = y[index]
    # values = np.vstack([x, y])
    # kernel_3 = st.gaussian_kde(values)
    # f3 = np.reshape(kernel_3(positions).T, xx.shape)

    prior = np.loadtxt('/Users/cc/CLionProjects/PF/prior/round' + str(i) + '/truth' + str(k * 100) + '.txt')
    x = prior[:, 0]
    y = prior[:, 1]
    index = np.where(x < 5.0)[0]
    x = x[index]
    y = y[index]
    values = np.vstack([x, y])
    kernel_3 = st.gaussian_kde(values)
    f3 = np.reshape(kernel_3(positions).T, xx.shape)

    emd0 = compute_emd2(prior, obsv)
    emd1 = compute_emd2(obsv, pred)

    return f1, f2, f3, emd0, emd1


def compare_kde2D(i, j):
    xmin, xmax = -5.0, 5.0
    ymin, ymax = -5.0, 5.0
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]

    dir3 = '/Users/cc/CLionProjects/PF/result/round' + str(i) + '_' + str(j)
    wass = np.loadtxt(dir3 + '/savewass.txt')
    emd = np.min(wass, axis=1)
    fs = 15

    fig = plt.figure(figsize=(15, 12))
    gs = gridspec.GridSpec(3, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1, 1])
    axes = []
    for k in range(3*3):
        axes.append(plt.subplot(gs[k], frameon=False))
    l = 0
    for k in [1, 20, 36]:
        f1, f2, f3, emd0, emd1 = compute_kde2D2(i, j, k)
        axes[l].set_xlim(xmin, xmax)
        axes[l].set_ylim(ymin, ymax)
        axes[l].set_xticks([])
        axes[l].set_yticks([])
        if l == 0:
            axes[l].set_ylabel('observed', labelpad=8.5, fontsize=fs)
        cfset = axes[l].contourf(xx, yy, f1)
        axes[l + 3].set_title('WD1 = ' + str(format(emd0, '.3f')), loc='center', fontsize=fs)
        # axes[l+3].set_title('emd0 = ' + str(format(emd0, '.3f')) + ', emd1 = ' + str(format(emd1, '.3f')), loc='center', fontsize=fs)

        axes[l+6].set_xlim(xmin, xmax)
        axes[l+6].set_ylim(ymin, ymax)
        axes[l+6].set_xticks([])
        axes[l+6].set_yticks([])
        if l == 0:
            axes[l+6].set_ylabel('posterior', labelpad=8.5, fontsize=fs)
        cfset = axes[l+6].contourf(xx, yy, f2)
        axes[l+6].set_title('WD2 = ' + str(format(emd1, '.3f')), loc='center', fontsize=fs)

        axes[l+3].set_xlim(xmin, xmax)
        axes[l+3].set_ylim(ymin, ymax)
        axes[l+3].set_xticks([])
        axes[l+3].set_yticks([])
        if l == 0:
            axes[l+3].set_ylabel('prior', labelpad=8.5, fontsize=fs)
        cfset = axes[l+3].contourf(xx, yy, f3)
        axes[l].set_title('t = ' + str(k), loc='center', fontsize=fs)

        l = l + 1

    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8)
    # plt.axis('equal')
    cb_ax = fig.add_axes([0.85, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(cfset, cb_ax)
    cbar.ax.tick_params(labelsize=fs)
    plt.savefig(dir3 + '/density.png',  bbox_inches='tight', pad_inches=0.1)


def draw_room(k):

    left = -5.0
    right = 5.0
    up = 5.0
    down = -5.0
    width = 1.0
    r = 0.3
    door = np.array([right, 0])
    wall = np.array([[right, up, left, up],  # 上
                     [left, up, left, down],  # 左
                     [left, down, right, down],  # 下
                     [right, down, door[0], -width],
                     [door[0], width, right, up]])

    fig = plt.figure(figsize=(11, 10), facecolor='white')
    # fig = plt.figure()
    ax = fig.add_subplot(111, frameon=False)

    for j in range(len(wall)):
        ax.plot([wall[j][0], wall[j][2]], [wall[j][1], wall[j][3]], 'black')

    ax.plot([wall[4][0], wall[4][0] + width], [wall[4][1], wall[4][1]], 'black')
    ax.plot([wall[3][2], wall[3][2] + width], [wall[3][3], wall[3][3]], 'black')

    data = np.loadtxt('/Users/cc/CLionProjects/PF/observation/round9/obsv' + str(100*k) + '.txt')

    for i in range(100):
        if data[i, 0] < right:
            circle = mpatches.Circle(xy=(data[i, 0], data[i, 1]), color='grey', radius=r, alpha=1)
            ax.add_patch(circle)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(5, -0.2, 'exit', fontsize=30, color='red')
    plt.savefig('room.png', bbox_inches='tight', pad_inches=0)


def draw_hist():

    x1 = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 0.05, 0.05, 0.05, 0.05, 0.05]
    x2 = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 0.15, 0.15, 0.15, 0.15, 0.15]
    x3 = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95]

    fig = plt.figure(figsize=(21, 5))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
    ax1 = plt.subplot(gs[0], frameon=True)
    ax2 = plt.subplot(gs[1], frameon=True)
    ax3 = plt.subplot(gs[2], frameon=True)

    ax1.hist(x1, 10, density=1, facecolor='blue', alpha=0.5, edgecolor='black')
    ax2.hist(x2, 10, density=1, facecolor='blue', alpha=0.5, edgecolor='black')
    ax3.hist(x3, 10, density=1, facecolor='blue', alpha=0.5, edgecolor='black')
    ax1.set_ylabel('frequency', fontsize=15)
    ax1.set_title('U', fontsize=15)
    ax2.set_title('V', fontsize=15)
    ax3.set_title('W', fontsize=15)
    fig.text(0.5, 0.01, 'x', ha='center', fontsize=15)
    ax1.tick_params(labelsize=15)
    ax2.tick_params(labelsize=15)
    ax3.tick_params(labelsize=15)
    plt.savefig('TV_illustration.png', bbox_inches='tight', pad_inches=0.1)


i = 9
j = 2
# draw_sfm(i, j)
# draw_idm(i, j)
# draw_idm2(i, j)
compare_kde2D(i, j)

# for k in range(1, 30):
#     compare_kde2D(i, j, k)
    # draw_density(i, j, k)
    # draw_density2(i, j, k)
# draw_idm()



