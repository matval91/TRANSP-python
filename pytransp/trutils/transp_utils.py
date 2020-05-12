import numpy as np

def _fill_dict(ff, keys, varnames):
    """Fills dictionaries
    
    Given the name of the keys and the name of the variables,
    fills the dictionaries

    Parameters:
            keys      (array) : keys of the output dictionaries
        varnames  (array) : name of NETCDF file variables
    Attributes:
        name_dict (dict) : dict associating keys with varnames
        var_dict  (dict) : dict associating actual variables
    Note:

    """
    name_dict = dict.fromkeys(keys)
    var_dict  = dict.fromkeys(keys)

    for ii,kk in enumerate(keys):
        try:
            name_dict[kk] = varnames[ii]
            tmpdata = ff.variables[name_dict[kk]][:]
            var_dict[kk] = tmpdata
        except:
            print('No variable ', kk)
            var_dict[kk] = 0.

    return name_dict, var_dict

def _time_to_ind(time_transp, time, inj_index=[0]):
    """ times to indexes

    Gets the indexes relative to the times requested

    Parameters:
        time (arr) : array with the times   where to plot the lines
    Attributes:
        ind (arr)  : array with the indexes where to plot the lines
    Note:

    """
    if len(inj_index)!=1:
        time_transp = time_transp[inj_index]
    else:
        time_transp = time_transp
    time_transp=np.array(time_transp)
    ind=np.array([])
    try:
        ttt=time[0]
        if ttt[0]!=0:
            for t in ttt:
                ind = np.append(ind, np.argmin(time_transp-t<0.))
            ind=np.array(ind,dtype=int)
        else:
            ind = np.array([0])
            ind = ind.astype(int)
    except:
        ttt=time
        ind=np.append(ind, np.argmin(time_transp-ttt<0.))
            
    return np.array(ind,dtype=int)
