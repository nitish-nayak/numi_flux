import ROOT
import numpy as np

klp3 = np.array([-102303.94038275,535462.40995564,63064246.35882203])
kon_axis = np.array([0., 0., 1.])
cos_angle = np.dot(klp3, kon_axis)/np.sqrt(np.sum(klp3**2))
tan_phi = klp3[1]/klp3[0]
print(np.arccos(cos_angle))

angles = np.linspace(1.5, 7.5, 7)
z_fn = lambda x : (klp3[2]/cos_angle)*np.cos(x*1.E-3)

for angle in angles:
    z_angle = z_fn(angle)
    r_sqr = np.sum(klp3**2) - z_angle**2
    x_angle = np.sign(klp3[0])*np.sqrt(r_sqr/(1 + tan_phi**2))
    y_angle = x_angle*tan_phi
    angle_n = np.array([x_angle, y_angle, z_angle])
    print('const TVector3 kLP3_%.1fmrad(%s, %s, %s);'% (angle, str(x_angle), str(y_angle), str(z_angle)))


