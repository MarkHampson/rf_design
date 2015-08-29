import math
import numpy as np
import Touchstone as tch
import matplotlib.pyplot as plt

class Vector:
    def __init__(self, domain=set(), function={}):
        self.d = domain
        self.f = function
    def get_d(self, d):
        """ Vector can be sparse """
        if d in self.f:
            return self.f[d]
        else:
            return 0.0
    def set_d(self, d, f):
        self.f[d] = f

class Sparams:
    def __init__(self, sdict, fmt, z, d):
        """ 
        Parameters:

        Frequency points: d (domain)
        Four vectors: s11, s12, s21, s22 with domain d, encapsulated in a dictionary sdict
        Impedance: z
        Format: fmt (mag/angle, real/imaginary, or db/angle), angles in degrees
        """
        self.s11 = sdict['s11']
        self.s12 = sdict['s12']
        self.s21 = sdict['s21']
        self.s22 = sdict['s22']
        self.fmt = fmt
        self.z = z
        self.d = d

def ri_to_db(Sp):
    """ Convert ri (complex) to (db,angle) s-params """
    r,i = Sp.real, Sp.imag
    dbs = 20.0*math.log10( math.sqrt(r**2 + i**2) )
    rads = math.atan2(i, r)
    degs = rads * (180.0/math.pi)
    db = (dbs,degs)
    return db

def StoABCD(s2p, Zo=50.0):

    abcd = np.zeros( (len(s2p),), dtype=('c16',(2,2)) )
    for m, sp in zip(abcd,s2p):
        s11 = sp[0,0]
        s12 = sp[0,1]
        s21 = sp[1,0]
        s22 = sp[1,1]
        
        
        A = ((1+s11)*(1-s22)+s12*s21)/(2*s21)
        B = Zo*((1+s11)*(1+s22)-s12*s21)/(2*s21)
        C = ((1-s11)*(1-s22)-s12*s21)/(2*Zo*s21)
        D = ((1-s11)*(1+s22)+s12*s21)/(2*s21)

        m[0,0] = A
        m[0,1] = B
        m[1,0] = C
        m[1,1] = D
        
    return abcd

def AtoS(abcd, Zo=50.0):

    s2p = np.zeros( (len(abcd),), dtype=('c16',(2,2)) )
    abcd = abcd.astype(complex) # force type
    
    for m,sp in zip(abcd,s2p):
        A = m[0,0]
        B = m[0,1]
        C = m[1,0]
        D = m[1,1]

        den = A+B/Zo+C*Zo+D
        s11 = (A+B/Zo-C*Zo-D)/den
        s12 = (2*(A*D-B*C))/den
        s21 = 2/den
        s22 = (-A+B/Zo-C*Zo+D)/den
        sp[0,0] = s11
        sp[0,1] = s12
        sp[1,0] = s21
        sp[1,1] = s22

    return s2p

def s2pArray(s2p_ri):
    """S2P dict in RI format """
    s2p_freqs = sorted(s2p_ri.keys())
    dt_s2p = np.dtype( ('c16', (2,2)) )
    s2p = np.zeros( (len(s2p_freqs),), dtype = dt_s2p )

    for freq,sp in zip(s2p_freqs,s2p):
        getS = lambda arg: complex(*s2p_ri[freq][arg])
        sp[0,0] = getS("s11")
        sp[0,1] = getS("s12")
        sp[1,0] = getS("s21")
        sp[1,1] = getS("s22")

    return s2p

# def ri_to_db(s2p_ri):
#     """ Convert ri to db s-params """
#     s2p_db = np.zeros( s2p_ri.shape, s2p_ri.dtype)

    
    
#     return s2p_db
    
    # db = {}
    # for freq in self.sparams:
    #     for Sp in self.sparams[freq]:
    #         r,i = self.sparams[freq][Sp]
    #         dbs = 20.0*math.log10( math.sqrt(r**2 + i**2) )
    #         rads = math.atan2(i, r)
    #         degs = rads * (180.0/math.pi)
    #         if freq not in db:
    #             db[freq] ={}
    #         db[freq][Sp] = (dbs,degs)
    # return db

if __name__=="__main__":

    hfcn = tch.Touchstone("HFCN-3800___Plus25degC.S2P")
    
    s2p_ri = hfcn.ri                 # retreive complex
    s2p_freqs = hfcn.freqs           # retreive freqs
    
    sp_hfcn = s2pArray(s2p_ri)
    abcd_hfcn1 = StoABCD(sp_hfcn)
    abcd_hfcn2 = abcd_hfcn1.copy()

    total_ABCD = np.dot(abcd_hfcn1,abcd_hfcn2)

    total_S = AtoS(total_ABCD)
    
    # total_S_db = ri_to_db(total_S)
    
    

    
    # dt_db_ph = np.dtype( [('db','f8'),('ph','f8')] )
    # dt_s2p = np.dtype( [ ('s11',dt_db_ph),('s12',dt_db_ph),
    #                      ('s21',dt_db_ph),('s22',dt_db_ph)
    #                  ])
    # s2p = np.zeros((len(s2p_freqs), ), dtype=dt_s2p )
    # print(s2p_db[s2p_freqs[0]])
    
    

        

    # print(s2p['s21']['db'])
    # plt.plot( s2p_freqs, s2p['s21']['db'] )
    # plt.show()
