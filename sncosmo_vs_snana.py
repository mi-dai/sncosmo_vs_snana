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
                 54:'CfA4',
                 56:'Swift',
                 61:'CfA1',
                 62:'CfA2',
                 63:'CfA3',
                 64:'CfA3',
                 65:'CfA4',
                 66:'CfA4',
                 5:'CSP',
                 150:'Foundation',
                 1:'SDSS',
                 15:'PanSTARRS',
                 50:'Hamuy96',
                 51:'LOSS'}
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
        
def get_photometry_single_sn(snname,fmeta=None,flc=None,df_meta=None,df_lc=None):
    meta = df_meta
    lc = df_lc
    if meta is None and fmeta is not None:
        meta = pd.read_csv(fmeta)
    if lc is None and flc is not None:
        lc = pd.read_csv(flc)
    meta_sn = meta.loc[meta.Name_upper==snname.upper()]
    lc_sn = lc.loc[lc.Name_upper==snname.upper()]
    return meta_sn,lc_sn

def convert_lc_for_sncosmo(lc,zp=27.5,select_filts=None):
    dflist = []
    for f in lc['Filter'].unique():
        if isinstance(select_filts,(str,list)) and f.split('-')[1][0] in select_filts:       
            try:
                sn_f = lc.set_index('Filter').loc[f].reset_index()
                sn_f['flux'] = np.power(10.,-0.4*(sn_f['Mag']-zp))
                sn_f['flux_err'] = np.absolute(0.921*sn_f['flux']*sn_f['MagErr'])
                s = sn_f['Survey'].unique()[0]
                b = f.split('-')[-1]
                mref,magsys = get_refmag(s,b)
                sn_f['zp'] = zp-mref
                sn_f['zpsys'] = magsys
                dflist.append(sn_f)
            except:
                print("Error in {} band lc".format(f))
                continue
        else:
            continue
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

def get_mjdrange(group):
    return pd.Series({'first_fitmjd':group['MJD'].min(),'last_fitmjd':group['MJD'].max()})

def get_snana_fitmjd_range(filenames):
    if isinstance(filenames,str):
        filenames = [filenames]
    dflist = []
    for filename in filenames:
        df = pd.read_csv(filename,comment='#',sep='\s+')
        df_mjd = df.loc[df['DATAFLAG']==1].groupby('CID').apply(get_mjdrange)
        df_mjd['Survey'] = filename.split('/')[-1].split('_')[0].strip()
        dflist.append(df_mjd)
    df_res = pd.concat(dflist,sort=False).reset_index() 
    return df_res
    