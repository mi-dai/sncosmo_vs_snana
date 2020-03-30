import copy
import pandas as pd
import sncosmo
import numpy as np
import sys
sys.path.append("/home/mi/Desktop/project-sn-fitter/project/")
from sedfit.utils import get_refmag
from astropy.table import Table

class SNpackageGeneral():
    def __init__():
        pass
    
    def plot_result(lcparams):
        pass

    
class SNCosmo(SNpackageGeneral):
    def __init__():
        pass
    
    def plot_result(lcparams):
        pass


class SNANA(SNpackageGeneral):
    def __init__():
        pass
    
    def plot_result(lcparams):
        pass
    

def sncosmo_fitres_to_df(filenames):
    if isinstance(filenames,str):
        filenames = [filenames]
    dflist = []
    for filename in filenames:
        dflist.append(pd.read_csv(filename))
    return pd.concat(dflist,sort=False)        
        
def snana_fitres_to_df(filenames):
    surveymap = {53:'CfA3',
                 62:'CfA2',
                 5:'CSP'}
    if isinstance(filenames,str):
        filenames = [filenames]
    dflist = []
    for filename in filenames:
        dflist.append(pd.read_csv(filename,comment='#',sep='\s+'))
    df = pd.concat(dflist,sort=False) 
    offset = 0.27
    df['x0_offset'] = df['x0']*np.power(10.,-0.4*(offset))
    survey = [surveymap[x] for x in df['IDSURVEY']]
    df['Survey'] = survey
    return df
    
def get_lcpar_map(source='snana'):
    if source == 'snana':
        return {'x0':'x0_offset',
                'x1':'x1',
                'c':'c',
                't0':'PKMJD',
                'z':'zHEL',
                'mwebv':'MWEBV'}
    elif source == 'sncosmo':
        pardict = {}
        for p in ['x0','x1','c','t0','z','mwebv']:
            pardict[p] = p
        return pardict

def update_sncosmo_model(salt2par,lcpar_map):
    model = sncosmo.Model(source='salt2',
                      effects=[sncosmo.F99Dust()],
                      effect_names=['mw'],
                      effect_frames=['obs'])
    lcpardict = {}
    sncosmo_parlist = ['x0','x1','c','t0','z','mwebv']
    for p in sncosmo_parlist:
        lcpardict[p] = salt2par[lcpar_map[p]]
    model.update(lcpardict)
    return model

def compare_model_sncosmo(data=None,models=[],bands=None,**kwargs):
    if data is None and bands is None:
        raise ValueError("Must provide bands if no data is given")
    fig = sncosmo.plot_lc(data=data,model=models,bands=bands,**kwargs)
    return fig
        
def get_photometry_single_sn(snname,fmeta=None,flc=None):
    meta = pd.read_csv(fmeta)
    lc = pd.read_csv(flc)
    meta_sn = meta.loc[meta.Name_upper==snname.upper()]
    lc_sn = lc.loc[lc.Name_upper==snname.upper()]
    return meta_sn,lc_sn

def convert_lc_for_sncosmo(lc,zp=27.5):
    dflist = []
    for f in lc['Filter'].unique():
        sn_f = lc.set_index('Filter').loc[f].reset_index()
        sn_f['flux'] = np.power(10.,-0.4*(sn_f['Mag']-zp))
        sn_f['flux_err'] = np.absolute(0.921*sn_f['flux']*sn_f['MagErr'])
        s = sn_f['Survey'].unique()[0]
        b = f.split('-')[-1]
        mref,magsys = get_refmag(s,b)
        sn_f['zp'] = zp-mref
        sn_f['zpsys'] = magsys
        dflist.append(sn_f)
    return Table.from_pandas(pd.concat(dflist,sort=False))
    
# Function to find common elements in n arrays 
def commonElements(arr): 
      
    # initialize result with first array as a set 
    result = set(arr[0]) 
  
    # now iterate through list of arrays starting from 
    # second array and take intersection_update() of  
    # each array with result. Every operation will  
    # update value of result with common values in 
    # result set and intersected set 
    for currSet in arr[1:]: 
        result.intersection_update(currSet) 
  
    return list(result)