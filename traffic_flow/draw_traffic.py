from trafficSimulator import *
import pygame
from pygame import gfxdraw
import numpy as np
import ot


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
    # data = np.loadtxt('samples8/24.txt')
    # for i in range(len(data)):
    #     for j in range(len(sim.roads)):
    #         if data[i][1] == sim.roads[j].start[1]:
    #             config = {"path": [j]}
    #             theVehicle = Vehicle(config, v_max, a_max, b_max, T, s0)
    #             theVehicle.x = data[i, 0]
    #             theVehicle.v = data[i, 2]
    #             sim.roads[j].vehicles.append(theVehicle)
    return sim


class Window:
    def __init__(self, sim, config={}):
        # Simulation to draw
        self.sim = sim
        self.set_default_config()

        # Update configurations
        for attr, val in config.items():
            setattr(self, attr, val)

    def set_default_config(self):
        """Set default configuration"""
        self.width = 794 #794
        self.height = 200 #200
        self.bg_color = (255, 255, 255)
        self.fps = 60
        self.zoom = (5, 5)
        self.offset = (0, 0)

    def convert(self, x, y=None):
        """Converts simulation coordinates to screen coordinates"""
        if isinstance(x, list):
            return [self.convert(e[0], e[1]) for e in x]
        if isinstance(x, tuple):
            return self.convert(*x)
        return (
            int(self.width/2 + (x + self.offset[0])*self.zoom[0]),
            int(self.height/2 + (y + self.offset[1])*self.zoom[1])
        )

    def background(self, r, g, b):
        """Fills screen with one color."""
        self.screen.fill((r, g, b))

    def polygon(self, vertices, color, filled=True):
        gfxdraw.aapolygon(self.screen, vertices, color)
        if filled:
            gfxdraw.filled_polygon(self.screen, vertices, color)

    def rotated_box(self, pos, size, angle=None, cos=None, sin=None, centered=True, color=(0, 0, 255), filled=True):
        """Draws a rectangle center at *pos* with size *size* rotated anti-clockwise by *angle*."""
        x, y = pos
        l, h = size

        if angle:
            cos, sin = np.cos(angle), np.sin(angle)

        vertex = lambda e1, e2: (
            x + (e1 * l * cos + e2 * h * sin) / 2,
            y + (e1 * l * sin - e2 * h * cos) / 2
        )

        if centered:
            vertices = self.convert(
                [vertex(*e) for e in [(-1, -1), (-1, 1), (1, 1), (1, -1)]]
            )
        else:
            vertices = self.convert(
                [vertex(*e) for e in [(0, -1), (0, 1), (2, 1), (2, -1)]]
            )

        self.polygon(vertices, color, filled=filled)

    def arrow(self, pos, size, angle=None, cos=None, sin=None, color=(255, 255, 255)):
        if angle:
            cos, sin = np.cos(angle), np.sin(angle)

        self.rotated_box(
            pos,
            size,
            cos=(cos - sin) / np.sqrt(2),
            sin=(cos + sin) / np.sqrt(2),
            color=color,
            centered=False
        )

        self.rotated_box(
            pos,
            size,
            cos=(cos + sin) / np.sqrt(2),
            sin=(sin - cos) / np.sqrt(2),
            color=color,
            centered=False
        )

    def draw_roads(self):
        jj = 0
        for road in self.sim.roads:
            jj = jj + 1
            # Draw road background
            if jj < 5:
                self.rotated_box(
                    road.start,
                    (road.length, 3.7),
                    cos=road.angle_cos,
                    sin=road.angle_sin,
                    color=(128, 128, 128),
                    centered=False
                )

            # Draw road arrow
                if road.length > 5:
                    for i in np.arange(-0.5*road.length, 0.5*road.length, 10):
                        pos = (
                            road.start[0] + (road.length/2 + i + 3) * road.angle_cos,
                            road.start[1] + (road.length/2 + i + 3) * road.angle_sin
                        )

                        self.arrow(
                            pos,
                            (-1.25, 0.2),
                            cos=road.angle_cos,
                            sin=road.angle_sin
                        )

    def draw_vehicle(self, vehicle, road, color, filled):
        l, h = vehicle.l,  2
        sin, cos = road.angle_sin, road.angle_cos

        x = road.start[0] + cos * vehicle.x
        y = road.start[1] + sin * vehicle.x

        self.rotated_box((x, y), (l, h), cos=cos, sin=sin, centered=True, color=color, filled=filled)

    def draw_vehicles(self, color, filled):
        for road in self.sim.roads:
            # Draw vehicles
            for vehicle in road.vehicles:
                self.draw_vehicle(vehicle, road, color, filled)

    def draw_status(self, emd1, emd2):
        text_fps = self.text_font.render(f'WD1={emd1:.3}', False, (0, 0, 0))
        text_frc = self.text_font.render(f'WD2={emd2:.3}', False, (0, 0, 0))
        text_leg1 = self.text_font.render('prior', False, (0, 0, 0))
        text_leg2 = self.text_font.render(f'observed', False, (0, 0, 0))
        text_leg3 = self.text_font.render(f'posterior', False, (0, 0, 0))

        self.screen.blit(text_fps, (180, 12))
        self.screen.blit(text_frc, (520, 12))
        self.screen.blit(text_leg1, (80, 12))
        self.screen.blit(text_leg2, (400, 12))
        self.screen.blit(text_leg3, (720, 12))

    def draw(self):

        self.screen = pygame.display.set_mode((self.width, self.height))
        pygame.display.flip()
        clock = pygame.time.Clock()

        self.background(*self.bg_color)
        self.draw_roads()
        # self.rotated_box((, 98), (5, 2), angle=np.pi, color=(255, 0, 0))
        # self.rotated_box((190, 102), (5, 2), angle=np.pi, color=(255, 0, 0))

        self.draw_vehicles((0, 0, 255), True)
        self.rotated_box((100, 90), (5, 2), angle=np.pi, color=(255, 0, 0), filled=False)
        self.rotated_box((140, 90), (5, 2), angle=np.pi, color=(0, 0, 255), filled=True)
        self.rotated_box((180, 90), (5, 2), angle=np.pi, color=(0, 255, 0), filled=False)

        # fname = 'figs2/highway' + str(i) + '.png'
        # pygame.image.save(self.screen, fname)


if __name__ == '__main__':
    v_max = 16.6
    a_max = 1.44
    b_max = 4.61
    T = 1.6
    s0 = 2
    delta = 4
    t0_index = 12

    sim = initialize(v_max, a_max, b_max, T, s0, delta)
    win = Window(sim)
    win.offset = (-145, -100)
    # win.offset = (-145, -95)
    win.zoom = (8, 10)
    win.width = 794
    win.height = 240
    pygame.font.init()
    win.text_font = pygame.font.SysFont('serif', 18)
    minwass = np.load('result/round5_2/minwass.npy')

    # for i in range(1000):
    #     sim.run(1)
    #     win.draw(i)

    for k in range(1, 30):
        sim = initialize(16.66, 2.67, 4.61, 1.86, 2.0, 4.62)
        data = np.loadtxt('observation/round5/obsv' + str(k+t0_index) + '.txt').reshape(-1, 4)
        x0 = data[:, 0:2]
        w0 = np.ones((len(x0),)) / len(x0)
        for i in range(len(data)):
            for j in range(4):
                if data[i][1] == sim.roads[j].start[1]:
                    config = {"path": [j]}
                    theVehicle = Vehicle(config, v_max, a_max, b_max, T, s0)
                    theVehicle.x = data[i, 0]
                    theVehicle.v = data[i, 2]
                    sim.roads[j].vehicles.append(theVehicle)
        win.sim = sim
        win.draw()

        sim = initialize(8.33, 1.44, 4.61, 1.6, 2.0, 4.0)
        data = np.loadtxt('prior/round5/truth' + str(k + t0_index) + '.txt').reshape(-1, 4)
        x1 = data[:, 0:2]
        w1 = np.ones((len(x1),)) / len(x1)
        for i in range(len(data)):
            for j in range(4):
                if data[i][1] == sim.roads[j].start[1]:
                    config = {"path": [j]}
                    theVehicle = Vehicle(config, v_max, a_max, b_max, T, s0)
                    theVehicle.x = data[i, 0]
                    theVehicle.v = data[i, 2]
                    sim.roads[j].vehicles.append(theVehicle)
        # sim.run(10*k)
        win.sim = sim
        # win.draw(k)
        win.draw_vehicles((255, 0, 0), False)  # red
        #
        sim = initialize(16.66, 1.44, 4.61, 1.6, 2.0, 4.0)
        data = np.loadtxt('result/round5_2/' + str(k) + '_' + str(int(minwass[k, 0])) + '_' + str(k) + '.txt').reshape(-1, 2)
        x2 = data[:, 0:2]
        w2 = np.ones((len(x2),)) / len(x2)
        for i in range(len(data)):
            for j in range(4):
                if data[i][1] == sim.roads[j].start[1]:
                    config = {"path": [j]}
                    theVehicle = Vehicle(config, v_max, a_max, b_max, T, s0)
                    theVehicle.x = data[i, 0]
                    sim.roads[j].vehicles.append(theVehicle)
        win.sim = sim
        # sim.run(k * 10)
        win.draw_vehicles((0, 255, 0), False)  # green

        M = np.sqrt(ot.dist(x0, x1))
        emd1 = ot.emd2(w0, w1, M, numItermax=1e7)
        M = np.sqrt(ot.dist(x0, x2))
        emd2 = ot.emd2(w0, w2, M, numItermax=1e7)
        win.draw_status(emd1, emd2)

        fname = 'figs5/highway' + str(k) + '.png'
        pygame.image.save(win.screen, fname)

    # win.draw(1)



