import openpyxl
import numpy as np
import os
import matplotlib.pyplot as plt

class QEData():

    def __init__(self, fname):
        wb = openpyxl.load_workbook(fname, read_only = True, data_only = True)
        data = wb.get_sheet_by_name("Subset of Cell QE All")
        self.wl = [r[5].value for r in data.rows]
        self.qe = [r[6].value for r in data.rows]
        self.bandGap = [r[21].value for r in data.rows]

    def _derivatives(self):
        slopes = []
        slopes_w = []
        for i, w in enumerate(self.wl[2:96]):
            slopes.append((self.qe[i + 2] - self.qe[i + 1])/(w - self.wl[i + 1]))
            slopes_w.append(w - (w - self.wl[i + 1]) / 2)

        sec_der = []
        sec_der_w = []
        for i, w in enumerate(slopes_w[1:]):
            sec_der.append((slopes[i + 1] - slopes[i])/(w - slopes_w[i]))
            sec_der_w.append(w - (w - slopes_w[i]) / 2)

        return slopes_w, slopes, sec_der_w, sec_der

    def _find_critical_points(self):
        estimates = [410, 520, 580, 950, 1050]
        critical_points = [0, 0, 0, 0, 0]

        poss_points = [[], [], [], [], []]
        for i, w in enumerate(self.wl[1:97]):
            for j, e in enumerate(estimates):
                if np.abs(w - e) <= 20:
                    poss_points[j].append(i + 1)

        slopes_w, slopes, sec_der_w, sec_der = self._derivatives()

        for i, p in enumerate(poss_points[0]):
            j = get_by_w(self.wl[p], sec_der_w, sec_der)
            critical_points[0] = p
            if sec_der[j] > 0.001:
                break

        for i, p in enumerate(poss_points[1]):
            j = get_by_w(self.wl[p], sec_der_w, sec_der)
            critical_points[1] = p
            if sec_der[j] < 0:
                break

        for i, p in enumerate(poss_points[2]):
            j = get_by_w(self.wl[p], sec_der_w, sec_der)
            critical_points[2] = p
            if sec_der[j] > -0.002 and sec_der[j] < 0.002:
                break

        for i, p in enumerate(poss_points[3]):
            j = get_by_w(self.wl[p], slopes_w, slopes)
            critical_points[3] = p
            if slopes[j] < -0.08:
                break

        for i, p in enumerate(poss_points[4]):
            j = get_by_w(self.wl[p], slopes_w, slopes)
            critical_points[4] = p
            if slopes[j] < -0.18:
                break

        return critical_points

    def _best_fits(self):
        critical_points = self._find_critical_points()

        cp0 = critical_points[0]
        w = self.wl[1:cp0]
        qe = self.qe[1:cp0]

        curve1 = np.poly1d(np.polyfit(w, qe, 3))
        w_new = np.linspace(w[0], w[cp0 - 2], 50)
        qe_new = curve1(w_new)

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

def get_by_w(w, w_vals, y_vals):
    for i, y in enumerate(y_vals):
        if w_vals[i] > w:
            return i
    return -1

if __name__ == "__main__":
    QEData = QEData(os.path.join(os.pardir, 'current-loss-analysis/QEData.xlsx'))
    QEData.plot()
