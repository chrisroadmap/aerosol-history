from erf_aprp import driver



#expt = 'amip'
#model = 'E3SM'
#for run in range(1,4):
#    baserun = 'r%di1p1f1' % run
#    pertrun = 'r%di1p1f2' % run
#    driver(baserun, pertrun, model, expt, make_aprp=True, from_pp=False, e3sm=True)
#
#expt = 'histaer'
#model = 'IPSL-CM6A-LR'
#baserun = 'r1i1p1f1'
#pertrun = 'r1i1p1f1'
#driver(baserun, pertrun, model, expt, make_aprp=False, make_od550aer=True)
#
#model = 'GFDL-CM4'
#baserun='r1i1p1f1'
#expt='histaer'
#for pertrun in ['r1i1p1f1','r3i1p1f1']:
#    driver(baserun, pertrun, model, expt, make_aprp=False, make_od550aer=True)
#
#model = 'GISS-E2-1-G'
#baserun='r1i1p1f2'
#for pertrun in ['r1i1p1f2']:
#    for expt in ['histall', 'histghg', 'histnat']:
#    for expt in ['lu']:
#        driver(baserun, pertrun, model, expt, make_aprp=False)
#
model = 'NorESM2-LM'
baserun='r1i1p2f1'
expt='histaer'
for pertrun in ['r1i1p2f1','r2i1p2f1','r3i1p2f1']:
    driver(baserun, pertrun, model, expt, make_aprp=True, make_od550aer=True)
#
#model = 'CanESM5'
#baserun='r1i1p2f1'
#expt='histaer'
#for pertrun in ['r1i1p2f1', 'r2i1p2f1', 'r3i1p2f1']:
#    driver(baserun, pertrun, model, expt, make_aprp=False, make_od550aer=True)
#
#model = 'HadGEM3-GC31-LL'
#baserun='r1i1p1f3'
#expt='histaer'
#for pertrun in ['r1i1p1f3', 'r2i1p1f3', 'r3i1p1f3']:
#     driver(baserun, pertrun, model, expt, from_pp=False, make_aprp=True, make_od550aer=True)



#for pertrun in ['r2i1p1f3', 'r3i1p1f3']:
#    for expt in ['histnat']:
#        driver(baserun, pertrun, model, expt, from_pp=True, make_aprp=False)
#for pertrun in ['r3i1p1f3']:
#    for expt in ['histall']:
#        driver(baserun, pertrun, model, expt, from_pp=True, make_aprp=False)
#
#
#model = 'MIROC6'
#baserun='r1i1p1f1'
#for pertrun in ['r1i1p1f1','r2i1p1f1','r3i1p1f1']:
#    for expt in ['histall', 'histghg', 'histnat']:
#        driver(baserun, pertrun, model, expt, make_aprp=False)
#

#expt='histSST'
#model='NorESM2-LM'
#baserun='r1i1p1f1'
#pertrun='r1i1p1f1'
#driver(baserun, pertrun, model, expt, make_aprp=True, aerchemmip=True, make_od550aer=False, has_rfmip=True)
#
#model = 'UKESM1-0-LL'
#baserun='r1i1p1f2'
#pertrun='r1i1p1f2'
#driver(baserun, pertrun, model, expt, make_aprp=False, aerchemmip=True, make_od550aer=True)
#
#model = 'GFDL-ESM4'
#baserun='r1i1p1f1'
#pertrun='r1i1p1f1'
#driver(baserun, pertrun, model, expt, make_aprp=False, aerchemmip=True, make_od550aer=True)
#
#model = 'GISS-E2-1-G'
#baserun='r1i1p3f1'
#pertrun='r1i1p3f1'
#driver(baserun, pertrun, model, expt, make_aprp=True, aerchemmip=True, make_od550aer=True, has_rfmip=True)
#
#model = 'MIROC6'
#baserun='r1i1p1f1'
#pertrun='r1i1p1f1'
#driver(baserun, pertrun, model, expt, make_aprp=True, aerchemmip=True, make_od550aer=True, has_rfmip=True)
#
#model = 'MRI-ESM2-0'
#baserun = 'r1i1p1f1'
#pertrun = 'r1i1p1f1'
#driver(baserun, pertrun, model, expt, make_aprp=True, aerchemmip=True, make_od550aer=True)

