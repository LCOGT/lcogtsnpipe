#############################################################################################

def extintion(site):
    kk={}
    print '\n### site = '+site
    kk['lasilla'] =  {'U':0.46,'u':0.46, 'B':0.27, 'g':0.20 ,'V':0.12, 'r':0.09, 'R':0.09, 'i':0.02, 'I':0.02, 'z':0.03,'J':0.07, 'H':0.023,'K':0.045}  
    kk['paranal'] =  {'U':0.43,'u':0.43, 'B':0.22, 'g':0.18 ,'V':0.11, 'r':0.07, 'R':0.07, 'i':0.05, 'I':0.05, 'z':0.06,'J':0.11, 'H':0.06,'K': 0.07}
    kk['roque'] =    {'U':0.46,'u':0.46, 'B':0.22, 'g':0.16 ,'V':0.12, 'r':0.08, 'R':0.08, 'i':0.04, 'I':0.04, 'z':0.06,'J':0.12, 'H':0.06,'K': 0.09}
    kk['calaralto']= {'U':0.47,'u':0.47, 'B':0.35, 'g':0.30 ,'V':0.23, 'r':0.17, 'R':0.17, 'i':0.09, 'I':0.09, 'z':0.10,'J':0.102,'H':0.07,'K': 0.09}  
    kk['siding'] =   {'U':0.63,'u':0.70, 'B':0.32, 'g':0.26 ,'V':0.18, 'r':0.150,'R':0.13, 'i':0.08, 'I':0.07, 'z':0.06,'J':0.0,  'H':0.0, 'K': 0.0}
    kk['mauna'] =    {'U':0.45,'u':0.48, 'B':0.21, 'g':0.16 ,'V':0.12, 'r':0.09, 'R':0.07, 'i':0.04, 'I':0.03, 'z':0.03,'J':0.0,  'H':0.0, 'K': 0.0}
    kk['mcdonald'] = {'U':0.51,'u':0.51, 'B':0.23, 'g':0.20 ,'V':0.15, 'r':0.10, 'R':0.10, 'i':0.05, 'I':0.05, 'z':0.05,'J':0.0,  'H':0.0, 'K': 0.0}
    kk['southafrica'] = {'U':0.51,'u':0.51, 'B':0.23, 'g':0.20 ,'V':0.15, 'r':0.10, 'R':0.10, 'i':0.05, 'I':0.05, 'z':0.05,'J':0.0,  'H':0.0, 'K': 0.0}
    kk['asiago'] =   {'U':0.58,'u':0.58, 'B':0.29, 'g':0.20 ,'V':0.16, 'r':0.12, 'R':0.12, 'i':0.08, 'I':0.08, 'z':0.09,'J':0.16, 'H':0.09,'K': 0.19}
    kk['kait'] =     {'U':0.56,'u':0.56, 'B':0.28, 'g':0.22 ,'V':0.17, 'r':0.13, 'R':0.13, 'i':0.07, 'I':0.07, 'z':0.09,'J':0.0,  'H':0.0, 'K': 0.0}
    kk['bao'] =      {'U':0.76,'u':0.76, 'B':0.41, 'g':0.35 ,'V':0.30, 'r':0.28, 'R':0.28, 'i':0.13, 'I':0.13, 'z':0.15,'J':0.0,  'H':0.0, 'K': 0.0}
    kk['wise'] =     {'U':0.70,'u':0.70, 'B':0.60, 'g':0.55 ,'V':0.52, 'r':0.42, 'R':0.42, 'i':0.27, 'I':0.27, 'z':0.30,'J':0.0,  'H':0.0, 'K': 0.0}
    kk['hct'] =      {'U':0.36,'u':0.0 , 'B':0.21, 'g':0.0  ,'V':0.12, 'r':0.0,  'R':0.09, 'i':0.0,  'I':0.05, 'z':0.0, 'J':0.0,  'H':0.0, 'K': 0.0}
    kk['st'] =       {'U':0.57,'u':0.0 , 'B':0.28 ,'g':0.0  ,'V':0.17, 'r':0.0,  'R':0.11, 'i':0.0,  'I':0.07, 'z':0.0, 'J':0.0,  'H':0.0, 'K': 0.0}
    kk['ctio'] =     {'U':0.50,'u':0.516,'B':0.232,'g':0.203,'V':0.144,'r':0.116,'R':0.083,'i':0.077,'I':0.056,'z':0.04,'J':0.0,  'H':0.0, 'K': 0.0}
    kk['monsec'] =   {'U':0.50,'u':0.516,'B':0.23, 'g':0.23 ,'V':0.15, 'r':0.11, 'R':0.09, 'i':0.04, 'I':0.04, 'z':0.04,'J':0.0,  'H':0.0, 'K': 0.0}
    kk['0'] =        {'U':0.0 ,'u':0.0  ,'B':0.0,  'g':0.0  ,'V':0.0,  'r':0.0 , 'R':0.0,  'i':0.0,  'I':0.0,  'z':0.0, 'J':0.0,  'H':0.0, 'K': 0.0}
    return kk[site]

###############################################################################################

def colfix(instrument,ss='sloan'):
    colorefisso={}
    if ss in ['sloan','landolt']:
        colorefisso['fs01']={'UUB':0.0,'uug':0.0,'BUB':0.0,'BBV':0.0,'VBV':0.0,'VVR':0.0,\
                             'gug':0.0,'ggr':0.105,'RVR':0.0,'RRI':0.0,'rrz':0.0,'zrz':0.0,\
                             'rgr':0.013,'rri':0.029,'iri':0.0874,'iiz':0.0,'IRI':0.0,'ziz':-0.15}
        colorefisso['fs02']={'UUB':0.0,'uug':0.0,'BUB':0.0,'BBV':0.0,'VBV':0.0,'VVR':0.0,'ggi':0.0,'igi':0.0,\
                             'gug':0.0,'ggr':0.105,'RVR':0.0,'RRI':0.0,'rrz':0.0,'zrz':0.0,\
                             'rgr':0.013,'rri':0.029,'iri':0.0874,'iiz':0.0,'IRI':0.0,'ziz':-0.15}
        colorefisso['fs03']={'UUB':0.0,'uug':0.0,'BUB':0.0,'BBV':0.0,'VBV':0.0,'VVR':0.0,'ggi':0.0,'igi':0.0,\
                             'gug':0.0,'ggr':0.105,'RVR':0.0,'RRI':0.0,'rrz':0.0,'zrz':0.0,\
                             'rgr':0.013,'rri':0.029,'iri':0.0874,'iiz':0.0,'IRI':0.0,'ziz':-0.15}
        colorefisso['em03']={'UUB':0.0,'uug':0.0,'BUB':0.0,'BBV':0.0,'VBV':0.0,'VVR':0.0,\
                             'gug':0.0,'ggr':0.105,'RVR':0.0,'RRI':0.0,'rrz':0.0,'zrz':0.0,\
                             'rgr':0.013,'rri':0.029,'iri':0.0874,'iiz':0.0,'IRI':0.0,'ziz':-0.15}
        colorefisso['em01']={'UUB':0.0,'uug':0.0,'BUB':0.0,'BBV':0.0,'VBV':0.0,'VVR':0.0,\
                             'gug':0.0,'ggr':0.105,'RVR':0.0,'RRI':0.0,'rrz':0.0,'zrz':0.0,\
                             'rgr':0.013,'rri':0.029,'iri':0.0874,'iiz':0.0,'IRI':0.0,'ziz':-0.15}
        colorefisso['kb77']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':-0.02,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.034,'rri':0.025,'iri':0.071,'iiz':0.110,'IRI':0.013,'ziz':-0.04}
        colorefisso['kb73']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':-0.02,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.034,'rri':0.025,'iri':0.071,'iiz':0.110,'IRI':0.013,'ziz':-0.04}
        colorefisso['kb78']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':-0.02,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.034,'rri':0.025,'iri':0.071,'iiz':0.110,'IRI':0.013,'ziz':-0.04}
        colorefisso['kb74']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':-0.02,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.034,'rri':0.025,'iri':0.071,'iiz':0.110,'IRI':0.013,'ziz':-0.04}
        colorefisso['kb79']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':-0.02,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.034,'rri':0.025,'iri':0.071,'iiz':0.110,'IRI':0.013,'ziz':-0.04}
        colorefisso['kb70']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':-0.02,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.034,'rri':0.025,'iri':0.071,'iiz':0.110,'IRI':0.013,'ziz':-0.04}
        colorefisso['kb75']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':-0.02,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.034,'rri':0.025,'iri':0.071,'iiz':0.110,'IRI':0.013,'ziz':-0.04}
        colorefisso['kb05']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':-0.02,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.034,'rri':0.025,'iri':0.071,'iiz':0.110,'IRI':0.013,'ziz':-0.04}
        colorefisso['kb76']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':-0.02,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.034,'rri':0.025,'iri':0.071,'iiz':0.110,'IRI':0.013,'ziz':-0.04}
        colorefisso['kb71']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':-0.02,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.034,'rri':0.025,'iri':0.071,'iiz':0.110,'IRI':0.013,'ziz':-0.04}
        colorefisso['fl02']=colorefisso['fl03']=colorefisso['fl04']=colorefisso['kb71']
        colorefisso['fl05']=colorefisso['fl06']=colorefisso['fl07']=colorefisso['kb71']
        colorefisso['fl08']=colorefisso['fl09']=colorefisso['fl10']=colorefisso['kb71']
    elif ss=='apass':
        colorefisso['fs01']={'BBV':-0.0,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['fs02']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['fs03']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['em03']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['em01']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['kb77']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['kb73']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['kb78']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['kb74']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['kb79']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['kb70']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['kb75']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['kb05']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['kb76']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['kb71']={'BBV':-0.13,'VBV':0.0,'VVg':0.0,'gVg':0.0,'ggr':0.0,'rgr':0.0,'rri':0.0,'iri':0.0}
        colorefisso['fl02']=colorefisso['fl03']=colorefisso['fl04']=colorefisso['kb71']
        colorefisso['fl05']=colorefisso['fl06']=colorefisso['fl07']=colorefisso['kb71']
        colorefisso['fl08']=colorefisso['fl09']=colorefisso['fl10']=colorefisso['kb71']
    elif ss=='sloanprime':
        colorefisso['fs01']={'UUB':0.0,'uug':0.0,'BUB':0.0,'BBV':0.0,'VBV':0.0,'VVR':0.0,\
                             'gug':0.0,'ggr':0.0,'RVR':0.0,'RRI':0.0,'rrz':0.0,'zrz':0.0,\
                             'rgr':0.0,'rri':0.0,'iri':0.027,'iiz':0.0,'IRI':0.0,'ziz':0.0}
        colorefisso['fs02']={'UUB':0.0,'uug':0.0,'BUB':0.0,'BBV':0.0,'VBV':0.0,'VVR':0.0,\
                             'gug':0.0,'ggr':0.0,'RVR':0.0,'RRI':0.0,'rrz':0.0,'zrz':0.0,\
                             'rgr':0.0,'rri':0.0,'iri':0.027,'iiz':0.0,'IRI':0.0,'ziz':0.0}
        colorefisso['fs03']={'UUB':0.0,'uug':0.0,'BUB':0.0,'BBV':0.0,'VBV':0.0,'VVR':0.0,\
                             'gug':0.0,'ggr':0.0,'RVR':0.0,'RRI':0.0,'rrz':0.0,'zrz':0.0,\
                             'rgr':0.0,'rri':0.0,'iri':0.027,'iiz':0.0,'IRI':0.0,'ziz':0.0}
        colorefisso['em03']={'UUB':0.0,'uug':0.0,'BUB':0.0,'BBV':0.0,'VBV':0.0,'VVR':0.0,\
                             'gug':0.0,'ggr':0.0,'RVR':0.0,'RRI':0.0,'rrz':0.0,'zrz':0.0,\
                             'rgr':0.0,'rri':0.0,'iri':0.027,'iiz':0.0,'IRI':0.0,'ziz':0.0}
        colorefisso['em01']={'UUB':0.0,'uug':0.0,'BUB':0.0,'BBV':0.0,'VBV':0.0,'VVR':0.0,\
                             'gug':0.0,'ggr':0.0,'RVR':0.0,'RRI':0.0,'rrz':0.0,'zrz':0.0,\
                             'rgr':0.0,'rri':0.0,'iri':0.027,'iiz':0.0,'IRI':0.0,'ziz':0.0}
        colorefisso['kb77']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':0.054,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.003,'rri':-0.007,'iri':0.028,'iiz':0.110,'IRI':0.013,'ziz':-0.16}
        colorefisso['kb73']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':0.054,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.003,'rri':-0.007,'iri':0.028,'iiz':0.110,'IRI':0.013,'ziz':-0.16}
        colorefisso['kb78']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV': 0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':0.054,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.003,'rri':-0.007,'iri':0.028,'iiz':0.110,'IRI':0.013,'ziz':-0.16}
        colorefisso['kb74']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':0.054,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.003,'rri':-0.007,'iri':0.028,'iiz':0.110,'IRI':0.013,'ziz':-0.16}
        colorefisso['kb79']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':0.054,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.003,'rri':-0.007,'iri':0.028,'iiz':0.110,'IRI':0.013,'ziz':-0.16}
        colorefisso['kb70']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':0.054,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.003,'rri':-0.007,'iri':0.028,'iiz':0.110,'IRI':-0.015,'ziz':-0.16}
        colorefisso['kb75']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':0.054,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.003,'rri':-0.007,'iri':0.028,'iiz':0.110,'IRI':0.013,'ziz':-0.16}
        colorefisso['kb05']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':0.054,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.003,'rri':-0.007,'iri':0.028,'iiz':0.110,'IRI':0.013,'ziz':-0.16}
        colorefisso['kb76']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':0.054,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.003,'rri':-0.007,'iri':0.028,'iiz':0.110,'IRI':0.013,'ziz':-0.16}
        colorefisso['kb71']={'UUB':0.059,'uug':0.0,'BUB':-0.095,'BBV':0.06,'VBV':0.03,'VVR':-0.059,\
                             'gug':0.13,'ggr':0.054,'RVR':-0.028,'RRI':-0.033,'rrz':0.0,'zrz':0.0,'ggi':0.0,'igi':0.0,\
                             'rgr':0.003,'rri':-0.007,'iri':0.028,'iiz':0.110,'IRI':0.013,'ziz':-0.16}
        colorefisso['fl02']=colorefisso['fl03']=colorefisso['fl04']=colorefisso['kb71']
        colorefisso['fl05']=colorefisso['fl06']=colorefisso['fl07']=colorefisso['kb71']
        colorefisso['fl08']=colorefisso['fl09']=colorefisso['fl10']=colorefisso['kb71']
    elif ss=='natural':
        colorefisso['em01']={'UUB':0.0,'uug':0.0,'BUB':0.0,'BBV':0.0,'VBV':0.0,'VVR':0.0,\
                             'gug':0.0,'ggr':0.0,'RVR':0.0,'RRI':0.0,'rrz':0.0,'zrz':0.0,\
                             'rgr':0.0,'rri':0.0,'iri':0.0,'iiz':0.0,'IRI':0.0,'ziz':0.0}
        colorefisso['fs01']=colorefisso['fs02']=colorefisso['fs03']=colorefisso['em03']=colorefisso['em01']
        colorefisso['kb77']=colorefisso['kb73']=colorefisso['kb78']=colorefisso['kb74']=colorefisso['em01']
        colorefisso['kb79']=colorefisso['kb70']=colorefisso['kb75']=colorefisso['kb05']=colorefisso['kb71']=colorefisso['kb76']=colorefisso['em01']
        colorefisso['fl02']=colorefisso['fl03']=colorefisso['fl04']=colorefisso['kb71']
        colorefisso['fl05']=colorefisso['fl06']=colorefisso['fl07']=colorefisso['kb71']
        colorefisso['fl08']=colorefisso['fl09']=colorefisso['fl10']=colorefisso['kb71']
    print instrument
    return colorefisso[instrument]

################################################################################

def filterst(telescope):
    _filters={'lsc':{'U':'U','B':'B','V':'V','R':'R','I':'I','u':'up','g':'gp','r':'rp','i':'ip','z':'zs','landolt':'landolt','sloan':'sloan'},
              'elp':{'U':'U','B':'B','V':'V','R':'R','I':'I','u':'up','g':'gp','r':'rp','i':'ip','z':'zs','landolt':'landolt','sloan':'sloan'},
              'ftn':{'U':'U','B':'Bessell-B','V':'Bessell-V','R':'Bessell-R','I':'Bessell-I','u':'up','g':'SDSS-G','r':'SDSS-R','i':'SDSS-I','z':'Pan-Starrs-Z','landolt':'landolt','sloan':'sloan'},
              'fts':{'U':'U','B':'Bessell-B','V':'Bessell-V','R':'Bessell-R','I':'Bessell-I','u':'up','g':'SDSS-G','r':'SDSS-R','i':'SDSS-I','z':'Pan-Starrs-Z','landolt':'landolt','sloan':'sloan'}}
    _filters['1m0-03']=_filters['1m0-04']=_filters['1m0-05']=_filters['1m0-07']=_filters['1m0-08']=_filters['lsc']
    _filters['1m0-09']=_filters['1m0-10']=_filters['1m0-11']=_filters['1m0-12']=_filters['1m0-13']=_filters['2m0-01']=_filters['2m0-02']=_filters['lsc']
    _filters['cpt']=_filters['all']=_filters['ogg']=_filters['coj']=_filters['SDSS']=_filters['extdata']=_filters['lsc']
    return _filters[telescope]

###############################################################################

def filterst1(telescope):
    _filters={'lsc':{'U':'U','B':'B','V':'V','R':'R','I':'I','up':'u','gp':'g','rp':'r','ip':'i','zs':'z','landolt':'landolt','sloan':'sloan'},
              'elp':{'U':'U','B':'B','V':'V','R':'R','I':'I','up':'u','gp':'g','rp':'r','ip':'i','zs':'z','landolt':'landolt','sloan':'sloan'},
              'ftn':{'U':'U','Bessell-B':'B','Bessell-V':'V','Bessell-R':'R','Bessell-I':'I','up':'u','SDSS-G':'g','SDSS-R':'r','SDSS-I':'i','Pan-Starrs-Z':'z','landolt':'landolt','sloan':'sloan'},
              'fts':{'U':'U','Bessell-B':'B','Bessell-V':'V','Bessell-R':'R','Bessell-I':'I','up':'u','SDSS-G':'g','SDSS-R':'r','SDSS-I':'i','Pan-Starrs-Z':'z','landolt':'landolt','sloan':'sloan'}}

    _filters['1m0-03']=_filters['1m0-04']=_filters['1m0-05']=_filters['1m0-07']=_filters['1m0-08']=_filters['lsc']
    _filters['1m0-09']=_filters['1m0-10']=_filters['1m0-11']=_filters['1m0-12']=_filters['1m0-13']=_filters['2m0-01']=_filters['2m0-02']=_filters['lsc']
    _filters['cpt']=_filters['coj']=_filters['ogg']=_filters['lsc']
    _filters['all']={'U':'U','Bessell-B':'B','Bessell-V':'V','Bessell-R':'R','Bessell-I':'I',\
                     'B':'B','V':'V','R':'R','I':'I','up':'u','gp':'g','rp':'r','ip':'i','zs':'z','g':'g','r':'r','i':'i',\
                     'SDSS-G':'g','SDSS-R':'r','SDSS-I':'i','Pan-Starrs-Z':'z','landolt':'landolt','sloan':'sloan'}
    _filters['2m0-02']=_filters['2m0-01']=_filters['extdata']=_filters['SDSS']=_filters['all']
    return _filters[telescope]

#############################################################################

dictcol2={'rsloan':'r','gsloan':'g','isloan':'i','zsloan':'z','B':'B','V':'V','R':'R','I':'I'}

def dictionarycolor():
    dictcol={'rp':{'col0':['gsloan','rsloan'],'col1':['rsloan','isloan'],'mag':'rsloan'},'ip':{'col0':['rsloan'],'col1':['isloan'],'mag':'isloan'},\
             'zp':{'col0':['isloan'],'col1':['zsloan'],'mag':'zsloan'},'gp':{'col0':['gsloan'],'col1':['rsloan'],'mag':'gsloan'},\
             'B':{'col0':['B'],'col1':['V'],'mag':'B'},'V':{'col0':['B'],'col1':['V'],'mag':'V'}}
    return dictcol

############################################################################

def cameraavailable(site=''):
    import numpy as np
    import os,string,re    
    if not site: site=['lsc','elp','coj','cpt','ogg']
    else: site=[site]
    ascifile='http://pydevsba.lco.gtn/configdb/camera_mappings/'
    line='wget '+ascifile+' -O configtel.txt'
    os.system(line)
    data=np.genfromtxt('configtel.txt',str)
    vector=['Site','Observatory','Telescope','Camera',\
                'CameraType', 'Size', 'PScale', 'BinningAvailable', 
            'Std.Binning', 'Overhead', 'Autoguider', 'AutoguiderType', 'Filters']
    aa=zip(*data)
    dict={}
    for i in range(0,len(vector)):
        dict[vector[i]]=np.array(aa[i])
    goodcamera=    [i for i in range(len(dict['CameraType'])) if dict['CameraType'][i] in [\
            '1m0-SciCam-Sinistro','1m0-SciCam-SBIG','2m0-SciCam-Spectral']]

    dict1={}
    for i in dict:
        dict1[i]=dict[i][goodcamera]
    instrument=[]
    for i in site:
        instrument=list(instrument)+list(dict1['Camera'][dict1['Site']==i])
    telescope=[]
    for i in site:
        telescope=list(telescope)+list(dict1['Observatory'][dict1['Site']==i])
    _site=[]
    for i in site:
        _site=list(_site)+list(dict1['Site'][dict1['Site']==i])
    return dict1,instrument,telescope,_site

#######################################################################
