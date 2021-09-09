# -*- coding: utf-8 -*-
#
# OutputReader.py
#
# Provides an object oriented interface to raytracer output data
#

# Michael Starks        AFRL/RVBX           14 March 2021

class OutputReader:
    
    class RayGroup:
        """A class providing a logical interface to a single ray group."""
        
        def __init__(self, f, n):
            """Initialize the ray group object."""
            self.rayGroup = f['{0:d}'.format(n)]
            self.freqGroups = len(self.rayGroup.keys())
            self.frequencies = []
            for freqGroup in self.rayGroup.keys():
                self.frequencies.append(self.rayGroup[freqGroup].attrs['frequency-Hz'])
            try:
                rayComment = self.rayGroup['comment']
            except KeyError:
                self.comment = None
            else:
                self.comment = rayComment
            
        
        def GetNumberOfFrequencies(self):
            """Return the number of frequencies for which this ray group was traced."""
            return self.freqGroups

        
        def GetFrequencies(self):
            """Return the array of frequencies for which this ray group was traced."""
            return self.frequencies


        def GetComment(self):
            """Return the comment associated with the ray group."""
            return self.comment

        
        def GetRay(self, freq):
            """Get a logical interface to a single ray in this group."""
            import numpy as np
            for freqGroup in self.rayGroup.keys():
                this_freq = float(self.rayGroup[freqGroup].attrs['frequency-Hz'])
                req_freq = float(freq)
                if (np.abs(float(self.rayGroup[freqGroup].attrs['frequency-Hz']) - float(freq)) < 1e-6):
                    return OutputReader.Ray(self.rayGroup[freqGroup])
            raise IndexError


    class Ray:
        """A class providing a logical interface to a single ray."""
        
        def __init__(self, freqGroup):
            """Initialize the ray."""
            import numpy as np
            # Read everything directly into Numpy arrays
            self.frequency = float(freqGroup.attrs['frequency-Hz'])
            self.numIters = freqGroup.attrs['number-of-iters']
            self.startTime = freqGroup.attrs['start-time-s']
            self.endTime = freqGroup.attrs['end-time-s']
            self.runTime = self.endTime - self.startTime
            self.termCode = freqGroup.attrs['termination-code']
            if (self.termCode == 'Success'):
                self.termCond = freqGroup.attrs['termination-condition']
            else:
                self.termCond = 'null'
            if (self.numIters > 0):
                self.r = np.empty(freqGroup['r'].shape)
                freqGroup['r'].read_direct(self.r)
                self.n = np.empty(freqGroup['n'].shape)
                freqGroup['n'].read_direct(self.n)
                self.err = np.empty(freqGroup['err'].shape)
                freqGroup['err'].read_direct(self.err)
                self.tg = np.empty(freqGroup['tg'].shape)
                freqGroup['tg'].read_direct(self.tg)
            else:
                self.r = None
                self.n = None
                self.err = None
                self.tg = None


    def __init__(self, inFile):
        """A ray tracer output reader object."""
        import h5py

        self.f = h5py.File(inFile, 'r')
        self.rayGroups = 0
        for key in self.f['ray-data'].keys():
            try:
                _a = int(key)
            except ValueError:
                pass
            else:
                self.rayGroups += 1
        self.inputFile = self.f.attrs['input']
        return

    
    def GetInputFilename(self):
        """Return the input file used to generate the ray tracing."""
        return(self.inputFile)
    
    
    def GetNumberOfRayGroups(self):
        """Return the number of ray groups in the file."""
        return(self.rayGroups)


    def GetRayGroup(self, n):
        """Return a logical interface to a ray group."""
        if (n > self.rayGroups):
            raise IndexError
        return OutputReader.RayGroup(self.f['ray-data'], n)

    
    def GetFieldModel(self):
        """Return the field model used to trace the rays."""
        from FieldModels.BaseFieldModel import BaseFieldModel
        fieldModel = BaseFieldModel.CreateFromHDF5(self.f['models']['field-model'])
        return(fieldModel)

    
    def GetDensityModel(self):
        """Return the density model used to trace the rays."""
        from DensityModels.BaseDensityModel import BaseDensityModel
        densityModel = BaseDensityModel.CreateFromHDF5(self.f['models']['density-model'])
        return(densityModel)

    def GetConstraints(self):
        """Return a list of the constraints objects used in the ray tracing."""
        return(self.constraints)

