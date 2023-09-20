from numpy import log,sqrt,exp

        
class gm(object):
    
    def __init__(self):
        self.gmm_type = 'Intensity prediction equation'
        self.gmm_id = 'Allen et al. (2012, hypocentral)'
    
    def get_intensity(r, m,
                      a = 2.085,
                      b = 1.428,#.913,#1.06,
                      c = -1.402,#-1.107,#-0.0010,
                      d = 0.078,#.813,#-3.37,
                      s = 1,
                      m1=-0.209,
                      m2=2.042):

        rm = m1+m2*exp(m-5)

        if r <= 50:
            return a + b*m + c*log(sqrt(r**2 + rm**2))+s

        return a + b*m + c*log(sqrt(r**2 + rm**2))+d*log(r/50)+s

class ipe_allen2012_hyp(gm):
    pass
