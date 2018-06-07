import openpyxl
import numpy as np
import os
import matplotlib.pyplot as plt
import datetime

COULOMBS = 6.24150934E18
JMP = True

#JMP Scripts sends data + path for AM15G Excel data
if JMP:
    jmp_wl = eval(input("WL: "))
    jmp_qe = eval(input("QE: "))
    jmp_bandgap = eval(input("Bandgap: "))
    AM15DATAPATH = input("AM1.5G Data Path: ")
else:
    AM15DATAPATH = 'AM1.5G.xlsx'

''' Class for storing and analyzing QE Data over wavelength (wl). Finds piecwise
    best fit to estimate current loss.'''
class QEData():

    # Loads data from Excel and/or JMP
    def __init__(self, fname = ''):
        if JMP:
            self.wl = jmp_wl
            self.qe = jmp_qe
            self.bandGap = jmp_bandgap[0]
        else:
            wb = openpyxl.load_workbook(fname, read_only = True, data_only = True)
            data = wb.get_sheet_by_name('Subset of Cell QE All')
            self.wl = [r[5].value for r in data.rows][1:]
            self.qe = [r[6].value for r in data.rows][1:]
            self.datetime = [r[0].value for r in data.rows][1:]
            bandGaps = [r[21].value for r in data.rows][1:]
            self.bandGap = bandGaps[1]

        AM15G = openpyxl.load_workbook(AM15DATAPATH)
        AM15G_Data = AM15G.get_sheet_by_name('Data')
        self.AM_WL = [r[0].value for r in AM15G_Data.rows][1:]
        self.AM_PHOTONS = [r[1].value for r in AM15G_Data.rows][1:]

        #Fixes points from bad filter by shifting up to make 'continuous'
        self.bad_points = True
        self._fix_points()

    ''' Returns x,y values for a) slopes between QE,WL data points and
    #   b) slopes between points in a). The idea is to approximate
    #   first and second derivative values in order to determine partition
    #   points *before* doing best fits for each partition'''
    def _derivatives(self):
        wl = self.wl
        qe = self.qe
        # first order derivatives (wl, dqe/dwl)
        slopes = []
        slopes_w = []
        for i, w in enumerate(wl[1:]):
            slopes.append((qe[i + 1] - qe[i])/(w - wl[i - 1]))
            slopes_w.append(w - (w - wl[i]) / 2)

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
        wl = self.wl
        qe = self.qe
        estimates = [420, 520, 580, 950, 1050]
        critical_points = [0, 0, 0, 0, 0, 0]

        # Find points w/in 20nm of estimates
        poss_points = [[], [], [], [], []]
        for i, w in enumerate(wl):
            for j, e in enumerate(estimates):
                if np.abs(w - e) <= 20:
                    poss_points[j].append(i)

        slopes_w, slopes, sec_der_w, sec_der = self._derivatives()

        # Predicted inflection point (threshold 0.001)
        for i, p in enumerate(poss_points[0]):
            j = get_by_w(wl[p], sec_der_w) + 1
            critical_points[0] = p + 1
            if sec_der[j - 1] < 0:
                break

        # Predicted inflection point (no threshold)
        for i, p in enumerate(poss_points[1]):
            j = get_by_w(wl[p], sec_der_w) + 1
            critical_points[1] = p
            if sec_der[j - 1] < 0:
                break

        # Predicted inflection point (threshold 0.003)
        for i, p in enumerate(poss_points[2]):
            j = get_by_w(wl[p], sec_der_w) + 1
            critical_points[2] = p
            if sec_der[j] > -0.002 and sec_der[j] < 0.002:
                break

        jump_point = critical_points[2]
        jump = qe[jump_point] - qe[jump_point - 1]

        # Predicted point joining two lines (threshold 0.03)
        for i, p in enumerate(poss_points[3]):
            j = get_by_w(wl[p], slopes_w)
            critical_points[3] = p
            if slopes[j] < -0.03:
                break

        # Predicted point joining line w/ curve? (threshold 0.08)
        for i, p in enumerate(poss_points[4]):
            j = get_by_w(wl[p], slopes_w)
            critical_points[4] = p
            if slopes[j] < -0.08:
                break

        #Last critical point is bandGap
        critical_points[5] = get_by_w(1241.52/self.bandGap, wl)

        return critical_points

    # Fixes points before third critical point by shifting up to make 'continuous'.
    # Accounts for bad filter on machine
    def _fix_points(self):
        if self.bad_points:
            critical_points = self._find_critical_points()
            wl = self.wl

            jp = critical_points[2]
            jump = self.qe[jp] - self.qe[jp - 1]
            jump -= self.qe[jp - 1] - self.qe[jp - 2]
            jump += 2 * self.qe[jp - 2] - self.qe[jp - 3] - self.qe[jp - 1]

            for i, q_val in enumerate(self.qe[0: jp]):
                self.qe[i] = q_val + jump
            self.bad_points = False

    ''' Returns piecewise best fit data points'''
    def _best_fits(self):
        critical_points = self._find_critical_points()
        cp0, cp1, cp2, cp3, cp4, bandGap_i = critical_points

        wl = self.wl
        qe = self.qe

        w_new = self.AM_WL

        # First curve: parabola
        w = wl[:cp0 + 2]
        q = qe[:cp0 + 2]

        curve = np.poly1d(np.polyfit(w, q, 3))
        a, b = get_by_w(w[0], w_new), get_by_w(w[cp0], w_new)
        w_temp = w_new[a:b]
        qe_new = curve(w_temp)

        # Second curve: parabola
        w = wl[cp0 - 1:cp1 + 2]
        q = qe[cp0 - 1:cp1 + 2]

        curve = np.poly1d(np.polyfit(w, q, 4))
        a, b = get_by_w(w[1], w_new), get_by_w(w[cp1 - cp0 + 1], w_new)
        w_temp = w_new[a:b]
        qe_new = np.append(qe_new, curve(w_temp))

        # Third curve: parabola
        w = wl[cp1 - 1:cp2 + 1]
        q = qe[cp1 - 1:cp2 + 1]

        curve = np.poly1d(np.polyfit(w, q, 4))
        a, b = get_by_w(w[1], w_new), get_by_w(w[cp2 - cp1 + 1], w_new)
        w_temp = w_new[a:b]
        qe_new = np.append(qe_new, curve(w_temp))

        # Fourth curve: parabola
        w = wl[cp2 - 1:cp3 + 2]
        q = qe[cp2 - 1:cp3 + 2]

        curve = np.poly1d(np.polyfit(w, q, 4))
        a, b = get_by_w(w[1], w_new), get_by_w(w[cp3 - cp2 + 1], w_new)
        w_temp = w_new[a:b]
        qe_new = np.append(qe_new, curve(w_temp))

        # Fifth curve: line
        w = wl[cp3 - 1:cp4 + 2]
        q = qe[cp3 - 1:cp4 + 2]

        curve = np.poly1d(np.polyfit(w, q, 4))
        a, b = get_by_w(w[1], w_new), get_by_w(w[cp4 - cp3 + 1], w_new)
        w_temp = w_new[a:b]
        qe_new = np.append(qe_new, curve(w_temp))

        # Sixth curve: parabola. End at bandGap cutoff
        w = wl[cp4 - 1:bandGap_i + 2]
        q = qe[cp4 - 1:bandGap_i + 2]

        curve = np.poly1d(np.polyfit(w, q, 4))
        a, b = get_by_w(w[1], w_new), get_by_w(w[bandGap_i - cp4 + 1], w_new)
        w_temp = w_new[a:b]
        qe_new = np.append(qe_new, curve(w_temp))

        # Seventh curve: parabola. bandGap cutoff to end
        end_i = get_by_w(1300, wl)

        w = wl[bandGap_i - 1:end_i + 2]
        q = qe[bandGap_i - 1:end_i + 2]

        curve = np.poly1d(np.polyfit(w, q, 4))
        a = get_by_w(w[1], w_new)
        w_temp = w_new[a:]
        qe_new = np.append(qe_new, curve(w_temp))

        def interp(cp):
            val1 = qe_new[get_by_w(wl[cp], w_new) - 1]
            val2 = qe_new[get_by_w(wl[cp], w_new) + 1]
            qe_new[get_by_w(wl[cp], w_new)] = (val1 + val2) / 2

        interp(cp0)
        interp(cp1)
        interp(cp2)
        interp(cp3)
        interp(cp4)
        interp(bandGap_i)

        return w_new, qe_new.tolist()

    # Takes integral of QE (not used)
    def integralQE(self):
        w, q = self._best_fits()
        return integral(w, q) / ((1300 - 360) * 100)

    # Calculates total current loss by using best fit, convolving w/ Am15G data,
    # taking integral and subtracting from Am1.5G curve
    def current_loss(self):
        w, q = self._best_fits()
        AM15G_convolution = [q[i] / 100 * p for i, p in enumerate(self.AM_PHOTONS)]
        plt.plot(w, AM15G_convolution)

        current = integral(w, AM15G_convolution) / COULOMBS
        return (integral(w, self.AM_PHOTONS / COULOMBS) - current)

    # See above, but also returns best fit data, convolution, current generated,
    # and current data for partitions
    def current_loss_all_data(self):
        critical_points = self._find_critical_points()
        w, q = self._best_fits()
        AM15G_convolution = [q[i] / 100 * p for i, p in enumerate(self.AM_PHOTONS)]
        #plt.plot(w, AM15G_convolution)
        current_partitions_loss = [0, 0, 0, 0, 0, 0, 0]
        current_partitions = [0, 0, 0, 0, 0, 0, 0]
        wl = self.wl

        j = 360
        crit_ps = [wl[cp] for cp in critical_points]
        for i in range(6):
            crit_p = crit_ps[i]
            curr = integral(w, AM15G_convolution, a=j, b=crit_p) / COULOMBS
            curr_loss = (integral(w, self.AM_PHOTONS, a=j, b=crit_p) / COULOMBS ) - curr
            current_partitions[i] = curr
            current_partitions_loss[i] = curr_loss
            j = crit_p

        curr = integral(w, AM15G_convolution, a=j) / COULOMBS
        curr_loss = (integral(w, self.AM_PHOTONS, a = j) / COULOMBS) - curr
        current_partitions[6] = curr
        current_partitions_loss[6] = curr_loss

        current = integral(w, AM15G_convolution) / COULOMBS
        current_loss = (integral(w, self.AM_PHOTONS) / COULOMBS) - current
        return (w, q, AM15G_convolution, current_loss, current, crit_ps,
                current_partitions_loss, current_partitions)

    #Plots critical points, best fit, and original data (w/ fixed points)
    def plot(self):
        wl = self.wl
        qe = self.qe

        best_fit_x, best_fit_y = self._best_fits()

        critical_points_w = [wl[cp] for cp in self._find_critical_points()]
        critical_points_qe = [qe[cp] for cp in self._find_critical_points()]

        QEcurve = plt.plot(wl, qe, '-o',
            critical_points_w, critical_points_qe, 'ro',
            best_fit_x, best_fit_y)
        plt.axis([300, 1400, 0, 100])

        slopes = []
        slopes_w = []
        for i, w in enumerate(best_fit_x):
            slopes.append((best_fit_y[i] - best_fit_y[i - 1])/(w - best_fit_x[i - 1]))
            slopes_w.append(w - (w - best_fit_x[i]) / 2)

        plt.figure(10)
        plt.plot(slopes_w, slopes)

    # Plots qe data, first and decond 'derivatives'
    def plot_derivatives(self):
        wl = self.wl
        qe = self.qe

        slopes_w, slopes, sec_der_w, sec_der = self._derivatives()

        plt.subplot(311)
        QEcurve = plt.plot(wl, qe, '-o')
        plt.axis([300, 1400, 0, 100])
        plt.grid(True)

        plt.subplot(312)
        slopeCurve = plt.plot(slopes_w, slopes, '-o')
        plt.axis([300, 1400, -0.8, 0.8])
        plt.grid(True)

        plt.subplot(313)
        secDerCurve = plt.plot(sec_der_w, sec_der, '-o')
        plt.axis([300, 1400, -0.025, 0.025])
        plt.grid(True)


# Returns index corresponding to wavelength w
def get_by_w(w, w_vals):
    for i, v in enumerate(w_vals):
        if v >= w:
            return i
    return -1

# Performs Riemann sum using mean y-values
def integral(x, y, a = 360, b = 1300):
    i, j = get_by_w(a, x), get_by_w(b, x)
    sum = 0
    for k, y_val in enumerate(y[i:j]):
        sum += (y_val + y[k  + i + 1]) * (x[k + i + 1] - x[k + i]) / 2
    return sum

if __name__ == '__main__':
    # Talks to JMP
    if (JMP):
        QEData = QEData()
        w, q, convolution, current_loss, current, critical_points, \
            current_partitions_loss, current_partitions = QEData.current_loss_all_data()

        print('Wavelengths:', w)
        print('Best Fits:', q)
        print('Convolution:', convolution)
        print('AM15G:', QEData.AM_PHOTONS)
        print('Current Loss:', current_loss, 'A/m^2')
        print('Current Generated:', current, 'A/m^2')
        print('Critical Points:', critical_points)
        print('Current Loss Partitioned:', current_partitions_loss)
        print('Current Generated Partitioned:', current_partitions)

    else:
        #testing stuff
        QEData = QEData(os.path.join(os.pardir, 'current-loss-analysis/QEData.xlsx'))

        print('Integral QE over WL:', QEData.integralQE())

        plt.figure(2)
        QEData.plot()
        plt.grid(True)

        plt.figure(3)
        print('Current Loss:', QEData.current_loss(), 'A/m^2')

        QEData.plot_derivatives()
        plt.show()
