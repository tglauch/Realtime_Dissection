import numpy as np

class Ellipse(object):
    def __init__(self, center_ra, center_dec, settings):
        settings = [float(i) for i in settings]
        self.center_ra = center_ra
        self.center_dec = center_dec
        if len(settings) == 2:
            self.ra_ax = settings[0]
            self.dec_ax = settings[1]
        elif len(settings) == 4:
            arr = np.abs(np.array(settings))
            self.ra_ax = np.sum(arr[[0,1]]) / 2.
            self.center_ra = self.center_ra + (arr[0]-arr[1])/2
            self.dec_ax = np.sum(arr[[2,3]]) / 2.
            self.center_dec = self.center_dec + (arr[2]-arr[3])/2
        else:
            print('No Valid Ellipse')
        if self.dec_ax < self.ra_ax:
            self.rotation = 0
        else:
            self.rotation = 90
        print('The rotation is {}'.format(self.rotation))
        return

    def get_vou_cmd(self, radius):
        return [str(self.center_ra), str(self.center_dec), str(radius), str(np.max([60. * self.ra_ax, 60. *
self.dec_ax])),
                str(np.min([60. * self.ra_ax, 60. * self.dec_ax])), str(self.rotation + 90)]

    def get_max_extension(self):
        return np.max([60. * self.ra_ax, 60. * self.dec_ax ])


class Lightcurve(object):
    def __init__(self, bpath, time_windows):
        self.bpath = bpath
        self.time_windows =  time_windows
        self.time_window_results = []
        return


