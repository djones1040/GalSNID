#!/usr/bin/env python
# D. Jones - 4/19/16
import numpy as np

class galsnid:
    def __init__(self):
        self.clobber = False
        self.verbose = False

    def add_options(self, parser=None, usage=None, config=None):
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('-v', '--verbose', action="count", dest="verbose",default=1)
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

        parser.add_option('-c','--configfile', default=None, type="string",
                          help='configuration file with GalSNID options')
        
        return(parser)

    def main(self,infile,outfile):
        sndata = txtobj(infile)

        fout = open(outfile,'w')
        print >> fout, '# ID z PIa PIa_high PIa_low'

        PIa,PIbc,PII = 1.0,1.0,1.0
        PIa_low,PIbc_low,PII_low = 1.0,1.0,1.0
        PIa_high,PIbc_high,PII_high = 1.0,1.0,1.0
        for i in range(len(sndata.__dict__[self.options.galparams[0]])):
            for j in range(len(self.options.galparams)):
                # default probabilities
                for P,type in zip([PIa,PIbc,PII],['Ia','Ibc','II']):
                    P *= self.PIaTable.getP(self.options.galparams[j],sndata.__dict__[self.options.galparams[j]][i],type=type)

                # low probabilities (68% CI from training)
                for P,type in zip([PIa_low,PIbc_low,PII_low],['Ia','Ibc','II']):
                    P *= self.PIaTable.getP(self.options.galparams[j],sndata.__dict__[self.options.galparams[j]][i],low=True)

                # low probabilities (68% CI from training)
                for P,type in zip([PIa_high,PIbc_high,PII_high],['Ia','Ibc','II']):
                    P *= self.PIaTable.getP(self.options.galparams[j],sndata.__dict__[self.options.galparams[j]][i],high=True)

            # rates prior
            for P1,PIb,P2 in zip([PIa,PIa_low,PIa_high],[PIbc,PIbc_low,PIbc_high],[PII,PII_low,PII_high]):
                P1 *= self.ratesprior(sndata.__dict__[self.options.redshift][i],'Ia')
                PIb *= self.ratesprior(sndata.__dict__[self.options.redshift][i],'Ibc')
                P2 *= self.ratesprior(sndata.__dict__[self.options.redshift][i],'II')

            PIa_final = PIa/(PIa + PIbc + PII)
            PIa_final_low = PIa_low/(PIa_low + PIbc_low + PII_low)
            PIa_final_high = PIa_high/(PIa_high + PIbc_high + PII_high)

            print >> fout, '%s %.3f %.4f %.4f %.4f'%(sndata.__dict__[self.options.id][i],
                                                     sndata.__dict__[self.options.redshift][i],
                                                     PIa_final,PIa_final_low,PIa_final_high)

        fout.close()

    def ratesprior(self,z,class):
        rp = txtobj(self.options.ratespriorfile)
        frac = 'frac_%s'%class
        return(np.interp(z,rp.__dict__[self.options.redshift],rp.__dict__[frac]))

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
                                                            'PII-':float(line.split('II: ')[-1].split(',')[1]),
                                                            'NSNe':float(line.split(' ')[-1])}
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
                                                                       'PII-':float(line.split('II: ')[-1].split(',')[1]),
                                                                       'NSNe':float(line.split(' ')[-1])})
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
                                              'PII-':float(line.split('II: ')[-1].split(',')[1]),
                                              'NSNe':float(line.split(' ')[-1])}
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
                                                                           'PII-':float(line.split('II: ')[-1].split(',')[1]),
                                                                           'NSNe':float(line.split(' ')[-1])})


        fin.close()

    def getP(self,key,val,type='Ia',low=False,high=False):
        for i in range(len(self.__dict__[key])):
            if self.__dict__[key][i].has_key('par'):
                if self.__dict__[key][i]['par'] == val:
                    if not low and not high:
                        return(self.__dict__[key][i]['P%s'%type])
                    elif low:
                        return(self.__dict__[key][i]['P%s'%type]-self.__dict__[key][i]['P%s-'%type])
                    elif high:
                        return(self.__dict__[key][i]['P%s'%type]+self.__dict__[key][i]['P%s+'%type])
            else:
                if self.__dict__[key][i]['par_min'] <= val and self.__dict__[key][i]['par_max'] > val:
                    if not low and not high:
                        return(self.__dict__[key][i]['P%s'%type])
                    elif low:
                        return(self.__dict__[key][i]['P%s'%type]-self.__dict__[key][i]['P%s-'%type])
                    elif high:
                        return(self.__dict__[key][i]['P%s'%type]+self.__dict__[key][i]['P%s+'%type])
        print('Warning : value %s not defined for host parameter %s'%(val,key))
        return(1.0)
        
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
        config.read(options.paramfile)
    else: config=None
    parser = gs.add_options(usage=usagestring,config=config)
    options,  args = parser.parse_args()

    options.galparams = options.galparams.split(',')


    gs.options = options
    gs.verbose = options.verbose
    gs.clobber = options.clobber

    import numpy as np
    import pylab as p

    gs.classifyIbc(options.infile,options.outfile)
