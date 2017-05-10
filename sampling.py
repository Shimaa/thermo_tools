import numpy as np
import matplotlib.pyplot as plt

#
# Define some constants related to butanol species classes
#
ALCOHOL_LIKE_TYPE = 1
ALCOHOL_LIKE_HF_SD = 0.123333333
ALCOHOL_LIKE_S_SD = 0.37

HYDROXY_ALKYL_RADICAL_TYPE = 2
HYDROXY_ALKYL_RADICAL_HF_SD = 0.81
HYDROXY_ALKYL_RADICAL_S_SD = 1.165714286

PEROXY_RADICAL_TYPE = 3
PEROXY_RADICAL_HF_SD = 1.255	
PEROXY_RADICAL_S_SD = 0.9575

HYDROPEROXIDE_TYPE = 4
HYDROPEROXIDE_HF_SD = 2.6225
HYDROPEROXIDE_S_SD = 1.005

ALKOXY_RADICAL_LIKE_PROPOXY_TYPE = 5
ALKOXY_RADICAL_LIKE_PROPOXY_HF_SD = 1.7
ALKOXY_RADICAL_LIKE_PROPOXY_S_SD = 1.135


##############################################
# A funcion to generate random samples that are normally distributed given a mean and SD
##############################################
def doSampling(mu, sigma, numOfSamples): # mean and standard deviation

    # generate Randoms
    s = np.random.normal(mu, sigma, numOfSamples)
    #print (s)

    # verify the mean and the variance
    #assert (abs(mu - np.mean(s)) < 0.01), "Verification Error"
    #assert (abs(sigma - np.std(s, ddof=1)) < 0.01),  "Verification Error"

    '''
    count, bins, ignored = plt.hist(s, 1000, normed=True)
    plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
             np.exp( - (bins - mu)**2 / (2 * sigma**2) ), linewidth=2, color='r')
    plt.show()
    '''

    return s




doSampling (-56.40, 0.123333333, 1000)




