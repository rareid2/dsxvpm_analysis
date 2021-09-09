# EccDipole.py
#
# A tilted, offset dipole model of the Earth's geomagnetic field,
# based on the first eight IGRF coefficients
#
# NOTE:  IGRF coefficients are hard-coded for 2020

# Michael Starks     AFRL/RVBX      10 Feb 2021

from FieldModels.BaseFieldModel import BaseFieldModel


class EccDipole(BaseFieldModel):
    """A tilted, offset dipole model of the Earth's geomagnetic field"""

    # =======================
    # Import required modules
    # =======================
    import numpy as np
    from Tracer.Constants import Constants

    # =====================
    # Class setup activites
    # =====================
    
    # Create the conversion constants 
    # NOTE: Using 2020 IGRF coefficients
    _g10 = -29404.8
    _g11 = -1450.9
    _g20 = -2499.6
    _g21 = 2982.0
    _g22 = 1677.0
    _h11 = 4652.5
    _h21 = -2991.6
    _h22 = -734.6

    # Following Usoskin et al. [2010], except there is an error in the first term of L1
    # It is g11, not g10, as seen in Fraser-Smith [1987] and Olson & Deguen [2012]
    _B0 = np.sqrt(_g10**2+_g11**2+_h11**2)
    _L0 = 2*_g10*_g20+np.sqrt(3.0)*(_g11*_g21+_h11*_h21)
    _L1 = -_g11*_g20+np.sqrt(3.0)*(_g10*_g21+_g11*_g22+_h11*_h22)
    _L2 = -_h11*_g20+np.sqrt(3.0)*(_g10*_h21-_h11*_g22+_g11*_h22)
    _E = (_L0*_g10+_L1*_g11+_L2*_h11)/(4*_B0**2)

    # denominators (for brevity)
    _dm1 = 3*_B0**2
    _dm2 = np.sqrt(_g11**2+_h11**2)

    # offset vector [Mm]
    _d = np.matrix([[(_L1-_g11*_E)/_dm1],[(_L2-_h11*_E)/_dm1],[(_L0-_g10*_E)/_dm1]]) * Constants.EarthRadius

    # transformation matrix
    _a = np.matrix([[_g11*_g10/(_B0*_dm2),_g10*_h11/(_B0*_dm2),-_dm2/_B0],[_h11/_dm2,-_g11/_dm2,0],[-_g11/_B0,-_h11/_B0,-_g10/_B0]])


    # ==================
    # Base Class Methods
    # ==================

    @classmethod
    def GetModelID(cls):
        """Returns a string identifier for the model"""
        return('ecc-dipole')


    @classmethod
    def LocalField(cls,r):
        """Returns the local field intensity [T] and Cartesian direction [Mm] at location geocentric Cartesian location r [Mm]. A row vector is expected."""
        # NOTE: Optimization of Numpy power and transpose functions did not significantly improve performance
        import numpy as np
        from Tracer.Constants import Constants
        from Tracer.Cart2Sph import Cart2Sph
        from Tracer.VectorSph2Cart import VectorSph2Cart

        if (r.ndim < 2):
            r = np.array(r, ndmin=2)
        if (r.shape[1] != 3):
            r = np.reshape(r, (-1,3))
        
        # Convert geocentric to geomagnetic
        rmag = cls.GEOtoMag(r)
                
        # Convert geomagnetic Cartesian to geomagnetic spherical
        rmag_sph = Cart2Sph(rmag)

        # Compute the field magnitude
        a = np.sqrt(1.0 + 3.0*(np.square(np.sin(np.radians(rmag_sph[0,1])))))
        ERbyr3  = np.power(Constants.EarthRadius / rmag_sph[0,0], 3)
        Babs = 0.3012 * ERbyr3 * a / 1.0e4

        # Compute the location in spherical magnetic, convert to Cartesian magnetic, then rotate to geocentric Cartesian
        Bdir_sph = np.array((-2.0 * np.sin(np.radians(rmag_sph[0,1])) / a, np.cos(np.radians(rmag_sph[0,1])) / a, 0.0), ndmin=2)
        Bdir_cart = VectorSph2Cart(rmag, Bdir_sph)
        Bdir = np.transpose(np.matmul(np.transpose(cls._a), np.transpose(Bdir_cart)))

        return Babs,Bdir


    def StoreParameters(self, h5group):
        """Stores model ID and initialization parameters to the provided HDF5 group"""
        import h5py

        h5group.attrs['model-id'] = self.GetModelID()
        return

    def GetParametersDict(self):
        """Returns a dictionary object containing the model ID and initialization parameters"""
        return({"model-id" : self.GetModelID()
               })
    

    # =============
    # Special Methods
    # =============

    @classmethod
    def GEOtoMag(cls,r):
        """Convert a geocentric Cartesian coordinate [Mm] to eccentric dipole magnetic Cartesian coordinate [Mm].  Expects a row vector or an array where each row is a coordinate."""
        # NOTE: Optimization by not using Numpy transpose or multiply resulted in a 20% loss in performance
        import numpy as np

        if (r.ndim < 2):
            r = np.array(r, ndmin=2)
        if (r.shape[1] != 3):
            r = np.reshape(r, (-1,3))

        # Convert geocentric to geomagnetic (keep as a row vector)
        if (r.shape[0] == 1):
            rmag = np.array(np.transpose(np.matmul(cls._a, (np.transpose(r) - cls._d))))
        else:
            rmag = np.empty(r.shape)
            for ind in np.arange(0, r.shape[0]):
                rmag[ind,] = np.array(np.transpose(np.matmul(cls._a, (np.transpose(r[None,ind,:]) - cls._d))))

        return rmag


    @classmethod
    def MagtoGEO(cls, r):
        """Convert a eccentric dipole magnetic Cartesian coordinate [Mm] to a geocentric Cartesian coordinate [Mm].  Expects a row vector or an array where each row is a coordinate."""
        import numpy as np

        if (r.ndim < 2):
            r = np.array(r, ndmin=2)
        if (r.shape[1] != 3):
            r = np.reshape(r, (-1,3))

        # Convert geomagnetic to geocentric (keep as a row vector)
        if (r.shape[0] == 1):
            rgeo = np.array(np.transpose(np.matmul(np.transpose(cls._a),  np.transpose(r)) + cls._d))
        else:
            rgeo= np.empty(r.shape)
            for ind in np.arange(0, r.shape[0]):
                rgeo[ind,:] = np.array(np.transpose(np.matmul(np.transpose(cls._a), np.transpose(r[None,:,ind])) + cls._d))

        return rgeo


    @classmethod
    def ComputeL(cls, r):
        """Compute the L shell of a geocentric Cartesian coordinate [Mm].  Expects a row vector."""
        import numpy as np
        from Tracer.Constants import Constants
        from Tracer.Cart2Sph import Cart2Sph

        if (r.ndim < 2):
            r = np.array(r, ndmin=2)
        if (r.shape[1] != 3):
            r = np.reshape(r, (-1,3))

        rmag = cls.GEOtoMag(r)
        rmag_sph = Cart2Sph(rmag)
        L = rmag_sph[0,0] / Constants.EarthRadius / np.square(np.cos(np.radians(rmag_sph[0,1])))
        # Or use L = rmag_sph[0,0] * np.square(rmag_sph[0,0]) / Constants.EarthRadius / (np.square(r[0,0]) + np.square(r[0,1]))

        return L

    
    @classmethod
    def ComputeRforLMlat(cls, L, mlat):
        """Compute the radial distance [Mm] from the magnetic dipole center of the given L shell at the given magnetic latitude [deg]"""
        import numpy as np
        from Tracer.Constants import Constants

        r = L * Constants.EarthRadius * np.square(np.cos(np.radians(mlat)))

        return r


    @classmethod
    def ComputeRforBMlat(cls, Babs, mlat):
        """Compute the radial distance [Mm] from the magnetic dipole center for a given magnetic field strength [T] and magnetic latitude [deg]"""
        import numpy as np
        from Tracer.Constants import Constants

        a = np.sqrt(1.0 + 3.0 * np.square(np.sin(np.radians(mlat))))
        r = Constants.EarthRadius * np.power(1.0e4 / 0.3012 * Babs / a, -(1.0/3.0))

        return r

