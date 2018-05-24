import openpyxl
import numpy as np
import os
import matplotlib.pyplot as plt

''' Class for storing and analyzing QE Data over wavelength (wl). Finds piecwise
    best fit to estimate current loss.'''
class QEData():

    COULOMBS = 6.24150934E18

    #Currently loads data from Excel. TODO: accept args from JMP
    def __init__(self, fname):
        wb = openpyxl.load_workbook(fname, read_only = True, data_only = True)
        data = wb.get_sheet_by_name("Subset of Cell QE All")
        self.wl = [r[5].value for r in data.rows]
        self.qe = [r[6].value for r in data.rows]
        self.device = [r[2].value for r in data.rows]
        bandGaps = [r[21].value for r in data.rows]
        self.bandGap = bandGaps[1]

        AM15G = openpyxl.load_workbook("AM1.5G.xlsx")
        AM15G_Data = AM15G.get_sheet_by_name("Data")
        self.AM_WL = [r[0].value for r in AM15G_Data.rows]
        self.AM_PHOTONS = [r[1].value for r in AM15G_Data.rows]

        self.bad_points = True

    ''' Returns x,y values for a) slopes between QE,WL data points and
    #   b) slopes between points in a). The idea is to approximate
    #   first and second derivative values in order to determine partition
    #   points *before* doing best fits for each partition'''
    def _derivatives(self, device):
        wl = [w for i, w in enumerate(self.wl) if self.device[i] == device]
        # first order derivatives (wl, dqe/dwl)
        slopes = []
        slopes_w = []
        for i, w in enumerate(wl[2:96]):
            slopes.append((self.qe[i + 2] - self.qe[i + 1])/(w - wl[i + 1]))
            slopes_w.append(w - (w - wl[i + 1]) / 2)

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
    def _find_critical_points(self, device):
        estimates = [410, 520, 580, 950, 1050]
        critical_points = [0, 0, 0, 0, 0, 0]

        wl = [w for i, w in enumerate(self.wl) if self.device[i] == device]

        # Find points w/in 20nm of estimates
        poss_points = [[], [], [], [], []]
        for i, w in enumerate(wl):
            for j, e in enumerate(estimates):
                if np.abs(w - e) <= 20:
                    poss_points[j].append(i)

        slopes_w, slopes, sec_der_w, sec_der = self._derivatives(device)

        # Predicted inflection point (threshold 0.001)
        for i, p in enumerate(poss_points[0]):
            j = get_by_w(wl[p] + 1, sec_der_w)
            critical_points[0] = p
            if sec_der[j] > 0.001:
                break

        # Predicted inflection point (no threshold)
        for i, p in enumerate(poss_points[1]):
            j = get_by_w(wl[p] + 1, sec_der_w)
            critical_points[1] = p
            if sec_der[j] < 0:
                break

        # Predicted inflection point (threshold 0.002)
        for i, p in enumerate(poss_points[2]):
            j = get_by_w(wl[p] + 1, sec_der_w)
            critical_points[2] = p
            if sec_der[j] > -0.001 and sec_der[j] < 0.001:
                break

        # Predicted point joining two lines (threshold 0.08)
        for i, p in enumerate(poss_points[3]):
            j = get_by_w(wl[p] + 1, slopes_w)
            critical_points[3] = p
            if slopes[j] < -0.08:
                break

        # Predicted point joining line w/ curve? (threshold 0.18)
        for i, p in enumerate(poss_points[4]):
            j = get_by_w(wl[p] + 1, slopes_w)
            critical_points[4] = p
            if slopes[j] < -0.18:
                break

        #Last critical point is bandGap
        critical_points[5] = get_by_w(1241.52/self.bandGap, wl)

        return critical_points

    ''' Returns piecewise best fit data points'''
    def _best_fits(self, device):
        critical_points = self._find_critical_points(device)
        cp0, cp1, cp2, cp3, cp4, bandGap_i = critical_points

        wl = [w for i, w in enumerate(self.wl) if self.device[i] == device]
        qe = [q for i, q in enumerate(self.qe) if self.device[i] == device]

        w_new = self.AM_WL[1:]

        # First curve: parabola
        w = wl[:cp0]
        q = qe[:cp0]

        curve = np.poly1d(np.polyfit(w, q, 3))
        a, b = get_by_w(w[0], w_new), get_by_w(w[cp0 - 1], w_new)
        w_temp = w_new[a:b]
        qe_new = curve(w_temp)

        # Second curve: parabola
        w = wl[cp0 - 1:cp1]
        q = qe[cp0 - 1:cp1]

        curve = np.poly1d(np.polyfit(w, q, 2))
        a, b = get_by_w(w[0], w_new), get_by_w(w[cp1 - cp0], w_new)
        w_temp = w_new[a:b]
        qe_new = np.append(qe_new, curve(w_temp))

        # Third curve: parabola
        w = wl[cp1 - 1:cp2 - 1]
        q = qe[cp1 - 1:cp2 - 1]

        curve = np.poly1d(np.polyfit(w, q, 2))
        a, b = get_by_w(w[0], w_new), get_by_w(w[cp2 - cp1 - 1], w_new)
        w_temp = w_new[a:b]
        qe_new = np.append(qe_new, curve(w_temp))

        bad_point_a = b

        # Fourth curve: line
        w = wl[cp2 - 1:cp3]
        q = qe[cp2 - 1:cp3]

        curve = np.poly1d(np.polyfit(w, q, 2))
        a, b = get_by_w(w[0], w_new), get_by_w(w[cp3 - cp2], w_new)
        w_temp = w_new[a:b]
        qe_new = np.append(qe_new, curve(w_temp))

        bad_point_b = a

        # Fifth curve: line
        w = wl[cp3 - 1:cp4]
        q = qe[cp3 - 1:cp4]

        curve = np.poly1d(np.polyfit(w, q, 1))
        a, b = get_by_w(w[0], w_new), get_by_w(w[cp4 - cp3], w_new)
        w_temp = w_new[a:b]
        qe_new = np.append(qe_new, curve(w_temp))

        # Sixth curve: parabola. End at bandGap cutoff
        w = wl[cp4 - 1:bandGap_i]
        q = qe[cp4 - 1:bandGap_i]

        curve = np.poly1d(np.polyfit(w, q, 2))
        a, b = get_by_w(w[0], w_new), get_by_w(w[bandGap_i - cp4], w_new)
        w_temp = w_new[a:b]
        qe_new = np.append(qe_new, curve(w_temp))

        # Sixth curve: parabola. bandGap cutoff to end
        end_i = get_by_w(1300, wl)

        w = wl[bandGap_i - 1:end_i + 1]
        q = qe[bandGap_i - 1:end_i + 1]

        curve = np.poly1d(np.polyfit(w, q, 3))
        a = get_by_w(w[0], w_new)
        w_temp = w_new[a:]
        qe_new = np.append(qe_new, curve(w_temp))

        del w_new[bad_point_a:bad_point_b]
        if self.bad_points:
            del self.AM_PHOTONS[bad_point_a + 1:bad_point_b + 1]
            self.bad_points = False
        return w_new, qe_new

    def overallQE(self, device):
        w, q = self._best_fits(device)
        return integral(w, q) / ((1300 - 360) * 100)

    def current_loss(self, device):
        w, q = self._best_fits(device)
        AM15G_convolution = [q[i] * p for i, p in enumerate(self.AM_PHOTONS[1:])]
        plt.plot(w, AM15G_convolution)

        current = integral(w, AM15G_convolution) / (self.COULOMBS * 100)
        return (integral(w, self.AM_PHOTONS[1:]) / self.COULOMBS) - current

    def plot(self, device):
        wl = [w for i, w in enumerate(self.wl) if self.device[i] == device]
        qe = [q for i, q in enumerate(self.qe) if self.device[i] == device]
        if (wl == []):
            raise Exception('No such device ID found in spreadsheet.')

        critical_points_w = [wl[cp - 1] for cp in self._find_critical_points(device)]
        critical_points_qe = [qe[cp - 1] for cp in self._find_critical_points(device)]

        best_fit_x, best_fit_y = self._best_fits(device)

        QEcurve = plt.plot(wl, qe, '-o',
            critical_points_w, critical_points_qe, 'ro',
            best_fit_x, best_fit_y)
        plt.axis([300, 1400, 0, 100])

    def plot_derivatives(self, device):
        wl = [w for i, w in enumerate(self.wl) if self.device[i] == device]
        qe = [q for i, q in enumerate(self.qe) if self.device[i] == device]

        slopes_w, slopes, sec_der_w, sec_der = self._derivatives(device)

        plt.subplot(311)
        QEcurve = plt.plot(wl, qe, '-o')
        plt.axis([300, 1400, 0, 100])

        plt.subplot(312)
        slopeCurve = plt.plot(slopes_w, slopes, '-o')
        plt.axis([300, 1400, -0.8, 0.8])

        plt.subplot(313)
        secDerCurve = plt.plot(sec_der_w, sec_der, '-o')
        plt.axis([300, 1400, -0.025, 0.025])


# Returns index corresponding to wavelength w
def get_by_w(w, w_vals):
    for i, v in enumerate(w_vals):
        if v >= w:
            return i
    return -1

def integral(x, y, a = 360, b = 1300):
    i, j = get_by_w(a, x), get_by_w(b, x)
    sum = 0
    for i, y_val in enumerate(y[i:j]):
        sum += (y_val + y[i + 1]) * (x[i + 1] - x[i]) / 2
    return sum

if __name__ == "__main__":
    QEData = QEData(os.path.join(os.pardir, 'current-loss-analysis/QEData.xlsx'))

    print("Overall QE:", QEData.overallQE(device = 200000))

    plt.figure(1)
    QEData.plot(200000)
    plt.grid(True)

    '''plt.figure(2)
    QEData.plot('2G5')
    plt.grid(True)

    plt.figure(3)
    QEData.plot('2K5')
    plt.grid(True)'''

    plt.figure(2)
    print("Current Loss:", QEData.current_loss(200000), "A/m^2")

    #QEData.plot_derivatives(200000)
    plt.show()
