import os

from astropy.io import fits

from . import caldb
from . import env


def make_arf(filename, weights):
    """
    Construct an ARF appropriate for stray light (e.g., the apbgd model)
    
    Parameters
    ----------    
    filename : str
        Path to .pha that you want to generate an ARF for
    
    weights : list
        Four element list giving the weighted value per detector to use.
        
    Returns
    ---------
    
    Path to generated ARF

    """

    
        
    hdr = fits.getheader(filename)
    arf_hdu = fits.open(os.path.join(env.auxildir, 'template.arf'))
    arf = arf_hdu[1].data
    arf_header = arf_hdu[1].header
    # Trim header
    arf_hdu[1].header = arf_header[0:26]
    

    
    # Get DETABS CALDB file for this FPM:
    detabs_file = env._CALDB.getDETABS(hdr['INSTRUME'], 'DET0', hdr['DATE-OBS'])
    beabs_file = env._CALDB.getBEABS('FPM', hdr['DATE-OBS'])
    beabs_file = os.path.join(env._CALDB_PATH, beabs_file)
    detabs_file = os.path.join(env._CALDB_PATH, detabs_file)

    # Load Be absorption file
    print(beabs_file)    
    beabs = fits.getdata(beabs_file, 1)
    
    # Mock up a unit ARF:
    arf['SPECRESP'] = beabs['ATT']
    out_arf = f"{filename.split('.')[0]}_nsb_nodetabs.arf"
    arf_hdu.writeto(out_arf, overwrite=True)

    detabs_hdu = fits.open(detabs_file)
    
    # Loop over each DETABS extension and scale the correctio
    for detind,isc in enumerate(weights):
        print(detind, isc)
        detabs = detabs_hdu[detind+1].data
        if detind == 0:
            arfabs = detabs['DETABS'] * isc
        else:
            arfabs += detabs['DETABS'] * isc
    detabs_hdu.close()
    arf['SPECRESP'] *= arfabs

    out_arf = f"{filename.split('.')[0]}_nsb.arf"

    arf_hdu.writeto(out_arf, overwrite=True)
    arf_hdu.close()
    
    return out_arf
    
