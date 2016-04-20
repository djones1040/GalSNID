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
            parser.add_option('-g','--galparams', default=config.get('main','galparams'), type="string",
                              help="""comma-separated list of host galaxy parameters to use 
for classifying; should match header of input file""")
            parser.add_option('-o','--outfile', default=config.get('main','outfile'), type="string",
                              help='output file with SN classifications')
            parser.add_option('-t','--classfile', default=config.get('main','classfile'), type="string",
                              help='file with Bayesian probabilities for each parameter in galparams')
            parser.add_option('-r','--ratespriorfile', default=config.get('main','ratespriorfile'), type="string",
                              help='a file providing the SN rates prior.  Use header "# z frac_Ia frac_Ibc frac_II"')
            parser.add_option('-z','--redshiftcol', default=config.get('main','redshiftcol'), type="string",
                              help='name of the redshift column in the input file')
            parser.add_option('--idcol', default=config.get('main','idcol'), type="string",
                              help='name of the ID column in the input file')
            parser.add_option('--spectypecol', default=config.get('main','spectypecol'), type="string",
                              help='name of the ID column in the input file')
            parser.add_option('--writetypecol', default=config.get('main','writetypecol'), type="int",
                              help='name of the ID column in the input file')
            parser.add_option('--snrthresh', default=map(float,config.get('main','snrthresh'))[0], type="string",
                              help="""for parameters with uncertainties, require a certain continuum SNR to use them 
in GalSNID (because GalSNID doesn\'t incorporate host galaxy measurement uncertainty.  GalSNID will look for the parameter 
followed by '_csnr' in the input file header to see if parameters have uncertainties""")
        else:
            parser.add_option('-i','--infile', default=None, type="string",
                              help='input file with galaxy data')
            parser.add_option('-g','--galparams', default=None, type="string",
                              help="""comma-separated list of host galaxy parameters to use 
for classifying; should match header of input file""")
            parser.add_option('-o','--outfile', default=None, type="string",
                              help='output file with SN classifications')
            parser.add_option('-t','--classfile', default='galsnid_class.txt', type="string",
                              help='file with Bayesian probabilities for each parameter in galparams')
            parser.add_option('-r','--ratespriorfile', default=None, type="string",
                              help='a file providing the SN rates prior.  Use header "# z frac_Ia frac_Ibc frac_II"')
            parser.add_option('-z','--redshiftcol', default='z', type="string",
                              help='name of the redshift column in the input file')
            parser.add_option('--idcol', default='ID', type="string",
                              help='name of the ID column in the input file')
            parser.add_option('--spectypecol', default='ID', type="string",
                              help='name of the ID column in the input file')
            parser.add_option('--writetypecol', default=0, type="int",
                              help='name of the ID column in the input file')
            parser.add_option('--snrthresh', default=5, type="string",
                              help="""for parameters with uncertainties, require a certain continuum SNR to use them 
in GalSNID (because GalSNID doesn\'t incorporate host galaxy measurement uncertainty.  GalSNID will look for the parameter 
followed by '_csnr' in the input file header to see if parameters have uncertainties""")


        parser.add_option('-c','--configfile', default=None, type="string",
                          help='configuration file with GalSNID options')
        
        return(parser)

    def main(self,infile,outfile):
        sndata = txtobj(infile)
        fout = open(outfile,'w')

        if not self.options.writetypecol:
            print >> fout, '# ID z PIa PIa_high PIa_low PIbc PIbc_high PIbc_low PII PII_high PII_low'
        else:
            print >> fout, '# ID z %s PIa PIa_high PIa_low PIbc PIbc_high PIbc_low PII PII_high PII_low'%(self.options.spectypecol)

        pt = PIaTable(trainfile=self.options.classfile)

        for i in range(len(sndata.__dict__[self.options.galparams[0]])):
            PIa,PIbc,PII = 1.0,1.0,1.0
            PIa_low,PIbc_low,PII_low = 1.0,1.0,1.0
            PIa_high,PIbc_high,PII_high = 1.0,1.0,1.0
            for j in range(len(self.options.galparams)):
                # SNR threshold
                if self.options.snrthresh and sndata.__dict__.has_key('%s_csnr'%self.options.galparams[j]):
                    if sndata.__dict__[self.options.galparams[j]][i]/sndata.__dict__['%s_csnr'%self.options.galparams[j]][i] < self.options.snrthresh:
                        if self.options.verbose: print('Measurement SNR for variable %s < threshold of %.1f'%(self.options.galparams[j],self.options.snrthresh))
                        continue

                # default probabilities
                PIa *= pt.getP(self.options.galparams[j],
                               sndata.__dict__[self.options.galparams[j]][i],
                               type='Ia',verbose=self.options.verbose)
                PIbc *= pt.getP(self.options.galparams[j],
                                sndata.__dict__[self.options.galparams[j]][i],
                                type='Ibc',verbose=self.options.verbose)
                PII *= pt.getP(self.options.galparams[j],
                               sndata.__dict__[self.options.galparams[j]][i],
                               type='II',verbose=self.options.verbose)

                # low probabilities (68% CI from training)
                PIa_low *= pt.getP(self.options.galparams[j],
                                   sndata.__dict__[self.options.galparams[j]][i],
                                   type='Ia',low=True,verbose=self.options.verbose)
                PIbc_low *= pt.getP(self.options.galparams[j],
                                    sndata.__dict__[self.options.galparams[j]][i],
                                    type='Ibc',low=True,verbose=self.options.verbose)
                PII_low *= pt.getP(self.options.galparams[j],
                                   sndata.__dict__[self.options.galparams[j]][i],
                                   type='II',low=True,verbose=self.options.verbose)

                # high probabilities (68% CI from training)
                PIa_high *= pt.getP(self.options.galparams[j],
                                    sndata.__dict__[self.options.galparams[j]][i],type='Ia',high=True,verbose=self.options.verbose)
                PIbc_high *= pt.getP(self.options.galparams[j],
                                     sndata.__dict__[self.options.galparams[j]][i],type='Ibc',high=True,verbose=self.options.verbose)
                PII_high *= pt.getP(self.options.galparams[j],
                                    sndata.__dict__[self.options.galparams[j]][i],type='II',high=True,verbose=self.options.verbose)

            # rates prior
            PIa *= self.ratesprior(sndata.__dict__[self.options.redshiftcol][i],'Ia')
            PIa_low *= self.ratesprior(sndata.__dict__[self.options.redshiftcol][i],'Ia')
            PIa_high *= self.ratesprior(sndata.__dict__[self.options.redshiftcol][i],'Ia')
            PIbc *= self.ratesprior(sndata.__dict__[self.options.redshiftcol][i],'Ibc')
            PIbc_low *= self.ratesprior(sndata.__dict__[self.options.redshiftcol][i],'Ibc')
            PIbc_high *= self.ratesprior(sndata.__dict__[self.options.redshiftcol][i],'Ibc')
            PII *= self.ratesprior(sndata.__dict__[self.options.redshiftcol][i],'II')
            PII_low *= self.ratesprior(sndata.__dict__[self.options.redshiftcol][i],'II')
            PII_high *= self.ratesprior(sndata.__dict__[self.options.redshiftcol][i],'II')

            # Ia probs
            try:
                PIa_final = PIa/(PIa + PIbc + PII)
            except:
                import pdb; pdb.set_trace()
            PIa_final_low = PIa_low/(PIa_low + PIbc_low + PII_low)
            PIa_final_high = PIa_high/(PIa_high + PIbc_high + PII_high)

            # Ibc probs
            PIbc_final = PIbc/(PIa + PIbc + PII)
            PIbc_final_low = PIbc_low/(PIa_low + PIbc_low + PII_low)
            PIbc_final_high = PIbc_high/(PIa_high + PIbc_high + PII_high)

            # II probs
            PII_final = PII/(PIa + PIbc + PII)
            PII_final_low = PII_low/(PIa_low + PIbc_low + PII_low)
            PII_final_high = PII_high/(PIa_high + PIbc_high + PII_high)

            if not self.options.writetypecol:
                print >> fout, '%s %.3f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f'%(
                    sndata.__dict__[self.options.idcol][i],
                    sndata.__dict__[self.options.redshiftcol][i],
                    PIa_final,PIa_final_low,PIa_final_high,
                    PIbc_final,PIbc_final_low,PIbc_final_high,
                    PII_final,PII_final_low,PII_final_high)
            else:
                print >> fout, '%s %.3f %s %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f'%(
                    sndata.__dict__[self.options.idcol][i],
                    sndata.__dict__[self.options.redshiftcol][i],
                    sndata.__dict__[self.options.spectypecol][i],
                    PIa_final,PIa_final_low,PIa_final_high,
                    PIbc_final,PIbc_final_low,PIbc_final_high,
                    PII_final,PII_final_low,PII_final_high)


        fout.close()

    def ratesprior(self,z,cls):
        rp = txtobj(self.options.ratespriorfile)
        frac = 'frac_%s'%cls
        return(np.interp(z,rp.__dict__[self.options.redshiftcol],rp.__dict__[frac]))

class PIaTable(object):
    def __init__(self,trainfile=None):
        if trainfile == None:
            raise exceptions.RuntimeError('Error : SN training file not provided!!!')

        fin = open(trainfile,'r')
        for line in fin:
            line = line.replace('\n','')
            key = line.split(':')[0]
            if ',' in line.split()[1]:
                if not self.__dict__.has_key(key):
                    self.__dict__[key.replace('_Ia','')] = {'par_min':float(line.split(' ')[1].split(',')[0]),
                                                            'par_max':float(line.split(' ')[1].split(',')[1]),
                                                            'PIa':float(line.split('Ia: ')[-1].split('CC: ')[0].split(',')[0]),
                                                            'PIa+':float(line.split('Ia: ')[-1].split('CC: ')[0].split(',')[2].split(' ')[0]),
                                                            'PIa-':float(line.split('Ia: ')[-1].split('CC: ')[0].split(',')[1]),
                                                            'PIbc':float(line.split('Ibc: ')[-1].split(',')[0]),
                                                            'PIbc+':float(line.split('Ibc: ')[-1].split(',')[2].split(' ')[0]),
                                                            'PIbc-':float(line.split('Ibc: ')[-1].split(',')[1]),
                                                            'PII':float(line.split('II: ')[-1].split(',')[0]),
                                                            'PII+':float(line.split('II: ')[-1].split(',')[2].split(' ')[0]),
                                                            'PII-':float(line.split('II: ')[-1].split(',')[1])}#,
                                                            #'NSNe':float(line.split(' ')[-1])}
                else:
                    self.__dict__[key] = np.append(self.__dict__[key],{'par_min':float(line.split(' ')[1].split(',')[0]),
                                                                       'par_max':float(line.split(' ')[1].split(',')[1]),
                                                                       'PIa':float(line.split('Ia: ')[-1].split('CC: ')[0].split(',')[0]),
                                                                       'PIa+':float(line.split('Ia: ')[-1].split('CC: ')[0].split(',')[2].split(' ')[0]),
                                                                       'PIa-':float(line.split('Ia: ')[-1].split('CC: ')[0].split(',')[1]),
                                                                       'PIbc':float(line.split('Ibc: ')[-1].split(',')[0]),
                                                                       'PIbc+':float(line.split('Ibc: ')[-1].split(',')[2].split(' ')[0]),
                                                                       'PIbc-':float(line.split('Ibc: ')[-1].split(',')[1]),
                                                                       'PII':float(line.split('II: ')[-1].split(',')[0]),
                                                                       'PII+':float(line.split('II: ')[-1].split(',')[2].split(' ')[0]),
                                                                       'PII-':float(line.split('II: ')[-1].split(',')[1])})#,
                                                                       #'NSNe':float(line.split(' ')[-1])})
            else:
                if not line.startswith('#'):
                    if not self.__dict__.has_key(key):
                        self.__dict__[key] = {'par':filter(None,line.split('%s:'%key)[1].split(' '))[0],
                                              'PIa':float(line.split('Ia: ')[-1].split('CC: ')[0].split(',')[0]),
                                              'PIa+':float(line.split('Ia: ')[-1].split('CC: ')[0].split(',')[2].split(' ')[0]),
                                              'PIa-':float(line.split('Ia: ')[-1].split('CC: ')[0].split(',')[1]),
                                              'PIbc':float(line.split('Ibc: ')[-1].split(',')[0]),
                                              'PIbc+':float(line.split('Ibc: ')[-1].split(',')[2].split(' ')[0]),
                                              'PIbc-':float(line.split('Ibc: ')[-1].split(',')[1]),
                                              'PII':float(line.split('II: ')[-1].split(',')[0]),
                                              'PII+':float(line.split('II: ')[-1].split(',')[2].split(' ')[0]),
                                              'PII-':float(line.split('II: ')[-1].split(',')[1])}#,
                                              #'NSNe':float(line.split(' ')[-1])}
                    else:
                        self.__dict__[key] = np.append(self.__dict__[key],{'par':filter(None,line.split('%s:'%key)[1].split(' '))[0],
                                                                           'PIa':float(line.split('Ia: ')[-1].split('CC: ')[0].split(',')[0]),
                                                                           'PIa+':float(line.split('Ia: ')[-1].split('CC: ')[0].split(',')[2].split(' ')[0]),
                                                                           'PIa-':float(line.split('Ia: ')[-1].split('CC: ')[0].split(',')[1]),
                                                                           'PIbc':float(line.split('Ibc: ')[-1].split(',')[0]),
                                                                           'PIbc+':float(line.split('Ibc: ')[-1].split(',')[2].split(' ')[0]),
                                                                           'PIbc-':float(line.split('Ibc: ')[-1].split(',')[1]),
                                                                           'PII':float(line.split('II: ')[-1].split(',')[0]),
                                                                           'PII+':float(line.split('II: ')[-1].split(',')[2].split(' ')[0]),
                                                                           'PII-':float(line.split('II: ')[-1].split(',')[1])})#,
                                                                           #'NSNe':float(line.split(' ')[-1])})


        fin.close()

    def getP(self,key,val,type='Ia',low=False,high=False,verbose=False):
        for i in range(len(self.__dict__[key])):
            if self.__dict__[key][i].has_key('par'):
                if self.__dict__[key][i]['par'] == val:
                    if not low and not high:
                        return(self.__dict__[key][i]['P%s'%type])
                    elif low:
                        return(self.__dict__[key][i]['P%s-'%type])
                    elif high:
                        return(self.__dict__[key][i]['P%s+'%type])
            else:
                if self.__dict__[key][i]['par_min'] <= val and self.__dict__[key][i]['par_max'] > val:
                    if not low and not high:
                        return(self.__dict__[key][i]['P%s'%type])
                    elif low:
                        return(self.__dict__[key][i]['P%s-'%type])
                    elif high:
                        return(self.__dict__[key][i]['P%s+'%type])

        if verbose: print('Warning : value %s not defined for host parameter %s'%(val,key))
        return(1.0)
        
if __name__ == "__main__":
    usagestring="""Run GalSNID on a SN sample.

USAGE: GalSNID.py [options]

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

    options.galparams = options.galparams.split(',')


    gs.options = options
    gs.verbose = options.verbose
    gs.clobber = options.clobber

    import numpy as np
    import pylab as p

    gs.main(options.infile,options.outfile)
