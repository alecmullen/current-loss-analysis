'''# in best fits

    slopes = []
    slopes_w = []
    for i, w in enumerate(w_new):
        slopes.append((qe_new[i] - qe_new[i - 1])/(w - w_new[i - 1]))
        slopes_w.append(w - (w - w_new[i]) / 2)

    def interp2(cp):
        val1 = slopes[get_by_w(wl[cp], w_new) - 1]
        val2 = slopes[get_by_w(wl[cp], w_new) + 1]
        slopes[get_by_w(wl[cp], w_new)] = (val1 + val2) / 2

    interp2(cp0)
    interp2(cp1)
    interp2(cp2)
    interp2(cp3)
    interp2(cp4)
    interp2(bandGap_i)

    c = qe_new[0]
    for i, q in enumerate(qe_new):
        qe_new[i] = integral(slopes_w, slopes, b=slopes_w[i]) + c'''

'''
#Attempt to do continuous best fit by lin. regression. Computer can't handle it
def best_fit(self, datetime):
    wl = [w for i, w in enumerate(self.wl) if self.datetime[i] == datetime]
    qe = [q for i, q in enumerate(self.qe) if self.datetime[i] == datetime]

    def pw_curve(x, b0, c0, d0, x0, y0, c1, d1, x1, y1, c2, d2, x2, y2,
                c3, d3, x3, y3, x4, y4, c4, d4, c5, d5, x5, y5, b6, c6, d6):
        curve0 = lambda x: (y0 - (b0 * x0 + c0 * x0 ** 2 + d0 * x0 ** 3)
                            + b0 * x + c0 * x ** 2 + d0 * x ** 3)
        curve1 = lambda x: (y0 + ((y1 - y0)/(x1 - x0) - c1*(x1 - x0) - d1*(x1 - x0)**2) * (x - x0) +
                            c1 * (x - x0) ** 2 + d1 * (x - x0) ** 3)
        curve2 = lambda x: (y1 + ((y2 - y1)/(x2 - x1) - c2*(x2 - x1) - d2*(x2 - x1)**2) * (x - x1)
                            + c2 * (x - x1) ** 2 + d2 * (x - x1) ** 3)
        curve3 = lambda x: (y2 + ((y3 - y2)/(x3 - x2) - c3*(x3 - x2) - d3*(x3 - x2)**2) * (x - x2)
                            + c3 * (x - x2) ** 2 + d3 * (x - x2) ** 3)
        curve4 = lambda x: (y3 + ((y4 - y3)/(x4 - x3) - c4*(x4 - x3) - d4*(x4 - x3)**2) * (x - x3)
                            + c4 * (x - x3) ** 2 + d4 * (x - x3) ** 3)
        curve5 = lambda x: (y4 + ((y5 - y4)/(x5 - x4) - c5*(x5 - x4) - d5*(x5 - x4)**2) * (x - x4)
                            + c5 * (x - x4) ** 2 + d5 * (x - x4) ** 3)
        curve6 = lambda x: (y5 - (b6 * x5 + c6 * x5 ** 2 + d6 * x5 ** 3)
                            + b6 * x + c6 * x ** 2 + d6 * x ** 3)

        #pdb.set_trace()
        return np.piecewise(x, [x < x0, np.logical_and(x >= x0, x < x1), np.logical_and(x >= x1, x < x2),
                                np.logical_and(x >= x2, x < x3), np.logical_and(x >= x3, x < x4),
                                np.logical_and(x >= x4, x < x5), x >= x5],
                                [curve0, curve1, curve2, curve3, curve4, curve5, curve6])


    scale = 1e1
    scale2 = 1e-2
    p, _ = optimize.curve_fit(pw_curve, np.array(wl)*scale2, scale * np.array(qe), p0 = [1, 1, 1, 420*scale2, 50*scale, 1, 1, 520*scale2,
                                50*scale, 1, 0, 580*scale2, 50*scale, 1, 0, 950*scale2, 50*scale, 0, 0, 1050*scale2, 50*scale, 1, 1, (1241.52/self.bandGap)*scale2,
                                50*scale, 1, 1, 1],
                                bounds = ([-1e90, -1e90, -1e90, 400*scale2, 0, -1e90, -1e90, 500*scale2, 0, -1e90, 0, 560*scale2, 0, -1e90, 0, 930*scale2,
                                             0, 0, 0, 1030*scale2, 0, -1e90, -1e90, (1241.52/self.bandGap - 20)*scale2, 0, -1e90, -1e90, -1e90],
                                         [1e90, 1e90, 1e90, 440*scale2, scale*1e2, 1e90, 1e90, 540*scale2, scale*1e2, 1e90, 0.1, 600*scale2, scale*1e2, 1e90, 0.1, 970*scale2,
                                             scale*1e3, 0.1, 0.1, 1070*scale2, scale*1e2, 1e90, 1e90, (1241.52/self.bandGap + 20)*scale2, scale*1e2, 1e90, 1e90, 1e90]))
    w_new = self.AM_WL[1:]
    print(*p)
    return(w_new, [(1.0/scale) * val for val in pw_curve([w_n*scale2 for w_n in w_new], *p)])'''

'''#handling issues w/ JMP input
from concurrent.futures import ProcessPoolExecutor

def timeout_five(fnc, *args, **kwargs):
        with ProcessPoolExecutor() as p:
            f = p.submit(fnc, *args, **kwargs)
            return f.result(timeout=5)

    def get_from_jmp(s):
        try:
            return timeout_five(input, s)
        except Exception as e:
            raise Exception("Error receiving data sent from JSL script")
'''
