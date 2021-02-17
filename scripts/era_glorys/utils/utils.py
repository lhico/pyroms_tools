import numpy as np

def rotate_coords(xm, ym, ang_rot, degrees=False):
    xm_mean = np.mean(xm)
    ym_mean = np.mean(ym)

    xm1 = xm.copy() - xm_mean
    ym1 = ym.copy() - ym_mean

    ang_rot_rad = ang_rot
    if degrees:
        ang_rot_rad = ang_rot/180 *np.pi

    xrot = xm1*np.cos(ang_rot_rad) - ym1*np.sin(ang_rot_rad)
    yrot = xm1*np.sin(ang_rot_rad) + ym1*np.cos(ang_rot_rad)

    xrot += xm_mean
    yrot += ym_mean
    return xrot, yrot
