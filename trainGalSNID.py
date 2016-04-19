#!/usr/bin/env python
# D. Jones - 4/19/16
import numpy as np

class txtobj:
    def __init__(self,filename):
        import numpy as np
        fin = open(filename,'r')
        lines = fin.readlines()
        for l in lines:
            if l.startswith('#'):
                l = l.replace('\n','')
                coldefs = l.split()[1:]
                break
        try: coldefs
        except NameError:
            raise RuntimeError('Error : file %s has no header'%filename)
                
        with open(filename) as f:
            reader = [x.split() for x in f if not x.startswith('#')]

        i = 0
        for column in zip(*reader):
            try:
                self.__dict__[coldefs[i]] = np.array(column[:]).astype(float)
            except:
                self.__dict__[coldefs[i]] = np.array(column[:])
            i += 1

class galsnid:
    def __init__(self):
        self.clobber = False
        self.verbose = False

    def add_options(self, parser=None, usage=None, config=None):
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('-v', '--verbose', action="count", dest="verbose",default=0)
        parser.add_option('--debug', default=False, action="store_true",
                          help='debug mode: more output and debug files')
        parser.add_option('--clobber', default=False, action="store_true",
                          help='clobber output image')

        if config:
            parser.add_option('-i','--infile', default=config.get('main','infile'), type="string",
                              help='input file with galaxy data')
            parser.add_option('--galparamfile', default=config.get('main','galparamfile'), type="string",
                              help="""a file that specifies which host galaxy properties to train on
and how many bins for each parameter.  Header should be '# param value'""")
            parser.add_option('-o','--outfile', default=config.get('main','outfile'), type="string",
                              help='output file with SN classifications')
            parser.add_option('--idcol', default=config.get('main','idcol'), type="string",
                              help='name of the ID column in the input file')
            parser.add_option('--spectypecol', default=config.get('main','spectypecol'), type="string",
                              help='name of the ID column in the input file')
            parser.add_option('--snrthresh', default=config.get('main','snrthresh'), type="float",
                              help="""for parameters with uncertainties, require a certain continuum SNR to use them 
in GalSNID (because GalSNID doesn\'t incorporate host galaxy measurement uncertainty.  GalSNID will look for the parameter 
followed by '_csnr' in the input file header to see if parameters have uncertainties""")

        else:
            parser.add_option('-i','--infile', default=None, type="string",
                              help='input file with galaxy data')
            parser.add_option('--galparamfile', default='galparams.txt', type="string",
                              help="""a file that specifies which host galaxy properties to train on
and how many bins for each parameter.  Header should be '# param value'""")
            parser.add_option('-o','--outfile', default=None, type="string",
                              help='output file with SN classifications')
            parser.add_option('--idcol', default='ID', type="string",
                              help='name of the ID column in the input file')
            parser.add_option('--spectypecol', default='type', type="string",
                              help='name of the ID column in the input file')
            parser.add_option('--snrthresh', default=5, type="float",
                              help="""for parameters with uncertainties, require a certain continuum SNR to use them 
in GalSNID (because GalSNID doesn\'t incorporate host galaxy measurement uncertainty.  GalSNID will look for the parameter 
followed by '_csnr' in the input file header to see if parameters have uncertainties""")


        parser.add_option('-c','--configfile', default=None, type="string",
                          help='configuration file with GalSNID options')
        
        return(parser)

    def main(self):
        import os
        import scipy.stats

        sndata = txtobj(self.options.infile)
        gs = txtobj(self.options.galparamfile)

        if not os.path.exists(self.options.outfile) or self.options.clobber:
            fout = open(self.options.outfile,'w')
        else:
            print('file %s exists!  Not clobbering'%self.options.outfile)
            return()

        for p,v in zip(gs.param,gs.value):
            if self.options.snrthresh and sndata.__dict__.has_key('%s_csnr'%p):
                    snr =  sndata.__dict__['%s_csnr'%p]
            else: snr = np.array([99.0]*len(sndata.__dict__[p]))
            if ',' in v:
                v1,v2 = v.split(',')
                v1,v2 = float(v1),float(v2)
                Ia = len(np.where((sndata.__dict__[self.options.spectypecol] == 'Ia') &
                                  (sndata.__dict__[p] > v1) &
                                  (sndata.__dict__[p] < v2) & 
                                  (sndata.__dict__[p] != -99) &
                                  (snr > self.options.snrthresh))[0])
                Ibc = len(np.where((sndata.__dict__[self.options.spectypecol] == 'Ibc') &
                                   (sndata.__dict__[p] > v1) &
                                   (sndata.__dict__[p] < v2) &
                                   (sndata.__dict__[p] != -99) &
                                   (snr > self.options.snrthresh))[0])
                II = len(np.where((sndata.__dict__[self.options.spectypecol] == 'II') &
                                  (sndata.__dict__[p] > v1) &
                                  (sndata.__dict__[p] < v2) &
                                  (sndata.__dict__[p] != -99) &
                                  (snr > self.options.snrthresh))[0])
                totalIa = len(np.where((sndata.__dict__[self.options.spectypecol] == 'Ia') &
                                       (sndata.__dict__[p] != -99) &
                                       (snr > self.options.snrthresh))[0])
                totalIbc = len(np.where((sndata.__dict__[self.options.spectypecol] == 'Ibc') &
                                       (sndata.__dict__[p] != -99) &
                                        (snr > self.options.snrthresh))[0])
                totalII = len(np.where((sndata.__dict__[self.options.spectypecol] == 'II') &
                                       (sndata.__dict__[p] != -99) &
                                       (snr > self.options.snrthresh))[0])
                
            else:
                Ia = len(np.where((sndata.__dict__[self.options.spectypecol] == 'Ia') &
                                  (sndata.__dict__[p] == v) &
                                  (sndata.__dict__[p] != '-99') &
                                  (snr > self.options.snrthresh))[0])
                Ibc = len(np.where((sndata.__dict__[self.options.spectypecol] == 'Ibc') &
                                   (sndata.__dict__[p] == v) &
                                   (sndata.__dict__[p] != '-99') &
                                   (snr > self.options.snrthresh))[0])
                II = len(np.where((sndata.__dict__[self.options.spectypecol] == 'II') &
                                  (sndata.__dict__[p] == v) &
                                  (sndata.__dict__[p] != '-99') &
                                  (snr > self.options.snrthresh))[0])
                totalIa = len(np.where((sndata.__dict__[self.options.spectypecol] == 'Ia') &
                                       (sndata.__dict__[p] != '-99') &
                                       (snr > self.options.snrthresh))[0])
                totalIbc = len(np.where((sndata.__dict__[self.options.spectypecol] == 'Ibc') &
                                       (sndata.__dict__[p] != '-99') &
                                        (snr > self.options.snrthresh))[0])
                totalII = len(np.where((sndata.__dict__[self.options.spectypecol] == 'II') &
                                       (sndata.__dict__[p] != '-99') &
                                       (snr > self.options.snrthresh))[0])
            errIal,errIau = scipy.stats.poisson.interval(0.68,Ia)
            errIbcl,errIbcu = scipy.stats.poisson.interval(0.68,Ibc)
            errIIl,errIIu = scipy.stats.poisson.interval(0.68,II)
            P_Ia_train = (Ia/float(totalIa),errIal/float(totalIa),errIau/float(totalIa))
            P_Ibc_train = (Ibc/float(totalIbc),errIbcl/float(totalIbc),errIbcu/float(totalIbc))
            P_II_train = (II/float(totalII),errIIl/float(totalII),errIIu/float(totalII))

            print >> fout,'%s: %s Ia: %.3f,%.3f,%.3f Ibc: %.3f,%.3f,%.3f II: %.3f,%.3f,%.3f'%(
                p,v,P_Ia_train[0],P_Ia_train[1],P_Ia_train[2],P_Ibc_train[0],P_Ibc_train[1],P_Ibc_train[2],
                P_II_train[0],P_II_train[1],P_II_train[2])

        fout.close()

if __name__ == "__main__":
    usagestring="""Run GalSNID on a SN sample.

USAGE: trainGalSNID.py [options]

"""
    import exceptions
    import os
    import optparse
    import ConfigParser

    gs=galsnid()

    parser = gs.add_options(usage=usagestring)
    options,  args = parser.parse_args()
    if options.configfile:
        config = ConfigParser.ConfigParser()
        config.read(options.configfile)
    else: config=None
    parser = gs.add_options(usage=usagestring,config=config)
    options,  args = parser.parse_args()


    gs.options = options
    gs.verbose = options.verbose
    gs.clobber = options.clobber

    import numpy as np
    import pylab as p

    gs.main()
