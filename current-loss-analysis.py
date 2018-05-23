import openpyxl
import numpy as np
import os
import matplotlib.pyplot as plt

''' Class for storing and analyzing QE Data over wavelength (wl). Finds piecwise
    best fit to estimate current loss.'''
class QEData():

    #Currently loads data from Excel. TODO: accept args from JMP
    def __init__(self, fname):
        wb = openpyxl.load_workbook(fname, read_only = True, data_only = True)
        data = wb.get_sheet_by_name("Subset of Cell QE All")
        self.wl = [r[5].value for r in data.rows]
        self.qe = [r[6].value for r in data.rows]
        self.bandGap = [r[21].value for r in data.rows]

    ''' Returns x,y values for a) slopes between QE,WL data points and
    #   b) slopes between points in a). The idea is to approximate
    #   first and second derivative values in order to determine partition
    #   points *before* doing best fits for each partition'''
    def _derivatives(self):
        # first order derivatives (wl, dqe/dwl)
        slopes = []
        slopes_w = []
        for i, w in enumerate(self.wl[2:96]):
            slopes.append((self.qe[i + 2] - self.qe[i + 1])/(w - self.wl[i + 1]))
            slopes_w.append(w - (w - self.wl[i + 1]) / 2)

        # second order derivatives (wl, d^2qe/dwl^2)
        sec_der = []
        sec_der_w = []
        for i, w in enumerate(slopes_w[1:]):
            sec_der.append((slopes[i + 1] - slopes[i])/(w - slopes_w[i]))
            sec_der_w.append(w - (w - slopes_w[i]) / 2)

        return slopes_w, slopes, sec_der_w, sec_der

    ''' Finds best partition points using hard-coded estimates and derivatives.
        Tries data points for wl values w/in 20nm of estimate for each point.
        On points predicted to be inflection points, uses 'second derivative'
        to pick point. On critical points estimated to join two straight lines,
        uses 'first derivatives'.'''
    def _find_critical_points(self):
        estimates = [410, 520, 580, 950, 1050]
        critical_points = [0, 0, 0, 0, 0]

        # Find points w/in 20nm of estimates
        poss_points = [[], [], [], [], []]
        for i, w in enumerate(self.wl[1:97]):
            for j, e in enumerate(estimates):
                if np.abs(w - e) <= 20:
                    poss_points[j].append(i + 1)

        slopes_w, slopes, sec_der_w, sec_der = self._derivatives()

        # Predicted inflection point (threshold 0.001)
        for i, p in enumerate(poss_points[0]):
            j = get_by_w(self.wl[p], sec_der_w, sec_der)
            critical_points[0] = p
            if sec_der[j] > 0.001:
                break

        #Predicted inflection point (no threshold)
        for i, p in enumerate(poss_points[1]):
            j = get_by_w(self.wl[p], sec_der_w, sec_der)
            critical_points[1] = p
            if sec_der[j] < 0:
                break

        #Predicted inflection point (threshold 0.002)
        for i, p in enumerate(poss_points[2]):
            j = get_by_w(self.wl[p], sec_der_w, sec_der)
            critical_points[2] = p
            if sec_der[j] > -0.002 and sec_der[j] < 0.002:
                break

        #Predicted point joining two lines (threshold 0.08)
        for i, p in enumerate(poss_points[3]):
            j = get_by_w(self.wl[p], slopes_w, slopes)
            critical_points[3] = p
            if slopes[j] < -0.08:
                break

        #Predicted point joining line w/ slope? (threshold 0.18)
        for i, p in enumerate(poss_points[4]):
            j = get_by_w(self.wl[p], slopes_w, slopes)
            critical_points[4] = p
            if slopes[j] < -0.18:
                break

        return critical_points

    ''' Returns piecewise best fit data points'''
    def _best_fits(self):
        critical_points = self._find_critical_points()
        cp0, cp1, cp2, cp3, cp4 = critical_points

        # First curve: parabola
        w = self.wl[1:cp0 + 1]
        qe = self.qe[1:cp0 + 1]

        curve = np.poly1d(np.polyfit(w, qe, 6))
        w_new = np.linspace(w[0], w[cp0 - 1], 50)
        qe_new = curve(w_new)

        # Second curve: parabola
        w = self.wl[cp0:cp1 + 1]
        qe = self.qe[cp0:cp1 + 1]

        curve = np.poly1d(np.polyfit(w, qe, 6))
        w_temp = np.linspace(w[0], w[cp1 - cp0], 50)
        w_new = np.append(w_new, w_temp)
        qe_new = np.append(qe_new, curve(w_temp))

        return w_new, qe_new

    def plot(self):
        #slopes_w, slopes, sec_der_w, sec_der = self._derivatives()
        critical_points_w = [self.wl[cp] for cp in self._find_critical_points()]
        critical_points_qe = [self.qe[cp] for cp in self._find_critical_points()]

        best_fit_x, best_fit_y = self._best_fits()

        #plt.subplot(311)
        QEcurve = plt.plot(self.wl[1:96], self.qe[1:96], '-o',
            critical_points_w, critical_points_qe, 'ro',
            best_fit_x, best_fit_y)
        plt.axis([300, 1400, 0, 100])

        '''plt.subplot(312)
        slopeCurve = plt.plot(slopes_w, slopes, '-o')
        plt.axis([300, 1400, -0.8, 0.8])

        plt.subplot(313)
        secDerCurve = plt.plot(sec_der_w, sec_der, '-o')
        plt.axis([300, 1400, -0.025, 0.025])'''

        plt.grid(True)
        plt.show()

# Returns index for y in y_vals corresponding to wavelength w
def get_by_w(w, w_vals, y_vals):
    for i, y in enumerate(y_vals):
        if w_vals[i] > w:
            return i
    return -1

if __name__ == "__main__":
    QEData = QEData(os.path.join(os.pardir, 'current-loss-analysis/QEData.xlsx'))
    QEData.plot()
