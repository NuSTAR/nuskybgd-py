import astropy.io.fits as pf


class CalDB:

    def __init__(self, path):
        self._CalDBPath = path
        self._index = self._loadCalDBIndex()

    def _loadCalDBIndex(self):
        indexfile = '%s/data/nustar/fpm/caldb.indx' % self._CalDBPath
        _fh = pf.open(indexfile)

    # ------------------------------------
    # Simple look-ups, perform same action
    def _printPath(keyword):
        

    def getPSFEnergy():
        # GRPPSF
        pass

    def getARF():
        # SPECRESP
        pass

    def getVignetting():
        # TVIGNET
        pass

    def getInstrumentMap():
        # INSTRMAP
        pass

    def getPixelLocation():
        # PIXPOS
        pass

    def getAperture():
        # APERTURE
        pass

    # Different
    def getDetectorAbsorption():
        pass

    # ------------------------------------


def getcaldbfile(caldb_path, ):
    pass

