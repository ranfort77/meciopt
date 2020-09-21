#!/usr/bin/env python

import numpy as np


def linesearch_halfstep(phi, phi0):
    ap = None
    for i in range(1, 13):  # 1.0/(2^12) = 0.000244
        ap0 = 1.0/(2.0**i)
        phi_ap0 = phi(ap0)
        if phi_ap0 <= phi0:
            break 
    return ap0, phi_ap0, i
    

def linesearch_armijo(phi, phi0, derphi0):
    # This function is the modified version of the
    # scipy.optimize.linesearch.scalar_search_armijo
    c1 = 1.0e-4
    ap0 = 1.0
    amin = 0.0
    i = 0
    
    phi_ap0 = phi(ap0)  # single point
    i += 1
    if phi_ap0 <= phi0 + c1*ap0*derphi0:
        return ap0, phi_ap0, i

    ap1 = -(derphi0) * ap0**2 / 2.0 / (phi_ap0 - phi0 - derphi0*ap0)
    phi_ap1 = phi(ap1)  # single point
    i += 1
    if phi_ap1 <= phi0 + c1*ap1*derphi0:
        return ap1, phi_ap1, i

    while ap1 > amin:
        factor = ap0**2 * ap1**2 * (ap1 - ap0)
        a = ap0**2 * (phi_ap1 - phi0 - derphi0*ap1) - \
            ap1**2 * (phi_ap0 - phi0 - derphi0*ap0)
        a = a / factor
        b = -ap0**3 * (phi_ap1 - phi0 - derphi0*ap1) + \
             ap1**3 * (phi_ap0 - phi0 - derphi0*ap0)
        b = b / factor

        ap2 = (-b + np.sqrt(abs(b**2 - 3 * a * derphi0))) / (3.0*a)
        phi_ap2 = phi(ap2)  # single point 
        i += 1
        if (phi_ap2 <= phi0 + c1*ap2*derphi0):
            return ap2, phi_ap2, i

        if (ap1 - ap2) > ap1 / 2.0 or (1 - ap2/ap1) < 0.96:
            ap2 = ap1 / 2.0
        ap0 = ap1
        ap1 = ap2
        phi_ap0 = phi_ap1
        phi_ap1 = phi_ap2

    # Failed to find a suitable step length
    return None, phi_ap1, i


if __name__ == '__main__':
    pass 
