# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import numpy as np
import scipy.interpolate

import lsst.afw.cameraGeom as afwCameraGeom

from .hscFilterDict import filterData


class DetectorThroughput(object):
    """
    Simple python class to hold detector throughput as a function of position

    This is intended as a placeholder until something more flexible can
    be implemented.

    This version is for HSC.

    """

    def __init__(self):
        """
        Create a DetectorThroughput instance.
        """

        self._setFilterData()
        self._setQEData()
        self._setMirrorData()
        self._pixelScale = 15e-3

    def getThroughputDetector(self, detector, filterName, lam):
        """
        Get the throughput at the central pixel of a camera detector

        Parameters
        ----------
        detector: lsst.afw.cameraGeom.Detector
        filterName: filter name (short name)
        lam: np.array() with wavelengths (units of nanometers)

        Returns
        -------
        throughput: np.array() with throughput at each lam
        """

        c = detector.getCenter(afwCameraGeom.FOCAL_PLANE)

        return self.getThroughputXY(filterName,
                                    c.getPoint().getX(), c.getPoint().getY(),
                                    lam)

    def getThroughputXY(self, filterName, x, y, lam):
        """
        Get the throughput at an arbitrary x/y on the focal plane

        Parameters
        ----------
        filterName: filter name (short name)
        x: x position (pixels) in focal plane coordinates
        y: y position (pixels) in focal plane coordinates
        lam: np.array() with wavelengths (units of nanometers)

        Returns
        -------
        throughput: np.array() with throughput at each lam
        """

        # Convert the filterData radius in pixels
        filterRadPix = filterData[filterName]['radius'] / self._pixelScale

        # Compute the radius, clipping to bounds
        #  Make sure that we keep within the boundaries by some small epsilon
        radius = np.clip(np.sqrt(x**2. + y**2.), filterRadPix[0] + 1e-7,
                         filterRadPix[-1] - 1e-7)

        # Which input wavelengths are non-zero
        inRange, = np.where((lam > filterData[filterName]['lam'][0]) &
                            (lam < filterData[filterName]['lam'][-1]))

        # set up 2d interpolation (bilinear)
        interpolator = scipy.interpolate.interp2d(filterRadPix,
                                                  filterData[filterName]['lam'],
                                                  filterData[filterName]['T'],
                                                  kind='linear', fill_value=0.0)


        # run interpolator
        # note that we need to flatten the nx1 array
        throughput = interpolator(radius, lam).flatten()

        # now multiply by average CCD QE (unless we had ccd-by-ccd QE)
        # values outside the range are set to 0
        interpolator = scipy.interpolate.interp1d(self.qeLambda, self.qeQe, fill_value=0.0)
        throughput *= interpolator(lam)

        # now multiply by mirror reflectivity
        interpolator = scipy.interpolate.interp1d(self.mirrorLambda,
                                                  self.mirrorReflectivity,
                                                  fill_value=0.0)
        throughput *= interpolator(lam)

        return throughput

    def getThroughputXYOld(self, filterName, x, y, lam):
        # leave this here temporarily

        # convert radius to pixels
        filterRad = self.filterData[filterName]['radius'] / self._pixelScale

        # keep within range.  Not sure what to do about extrapolation
        radius = np.clip(np.sqrt(x**2. + y**2.), filterRad[0], filterRad[-1])

        # and we need the interpolated on10, on50, on80, Tavg, off80, off50, off10
        yvec = np.array([0.0, 10.0, 50.0, 80.0, 0.0,
                         0.0, 80.0, 50.0, 10.0, 0.0])
        xvec = np.zeros_like(yvec)

        # get the on/off x values
        for i, key in enumerate((None, "on10", "on50", "on80",
                                None, None,
                                "off80", "off50", "off10", None)):
            if (key is None):
                continue

            interpolator = self.makeEvenSplineInterpolator(filterRad,
                                                           self.filterData[filterName][key])
            xvec[i] = interpolator(radius)

        # and the Tavg value
        interpolator = self.makeEvenSplineInterpolator(filterRad,
                                                       self.filterData[filterName]['Tavg'])
        yvec[4] = interpolator(radius)
        yvec[5] = interpolator(radius)

        # and linear interpolation to fill out the other points
        xvec[0] = xvec[2] + ((yvec[0] - yvec[2])*(xvec[2] - xvec[1]))/(yvec[2] - yvec[1])
        xvec[4] = xvec[3] + ((yvec[4] - yvec[3])*(xvec[3] - xvec[2]))/(yvec[3] - yvec[2])
        xvec[5] = xvec[7] + ((yvec[5] - yvec[7])*(xvec[7] - xvec[6]))/(yvec[7] - yvec[6])
        xvec[9] = xvec[8] + ((yvec[9] - yvec[8])*(xvec[8] - xvec[7]))/(yvec[8] - yvec[7])

        # now we can start the throughput...
        #  use linear interpolation
        throughput = np.zeros_like(lam, dtype='f8')

        # start with the filter
        interpolator = scipy.interpolate.interp1d(xvec, yvec)
        inRange, = np.where((lam > xvec[0]) & (lam < xvec[-1]))

        # divide by 100 so that 1.0 is proper transmission (not 100%)
        throughput[inRange] = interpolator(lam[inRange])/100.

        # now multiply by average CCD QE (will need to be updated)
        interpolator = scipy.interpolate.interp1d(self.qeLambda, self.qeQe)
        inRange, = np.where((lam > self.qeLambda[0]) & (lam < self.qeLambda[-1]))

        # this will set QE to zero outside the bounds
        qe = np.zeros_like(throughput)
        qe[inRange] = interpolator(lam[inRange])

        throughput *= qe

        # and the mirror reflectivity
        interpolator = scipy.interpolate.interp1d(self.mirrorLambda,
                                                  self.mirrorReflectivity)
        inRange, = np.where((lam > self.mirrorLambda[0]) &
                            (lam < self.mirrorLambda[-1]))

        throughput[inRange] *= interpolator(lam[inRange])

        return throughput

    def makeEvenSplineInterpolator(self, x, y):
        # this was taken from https://github.com/RobertLuptonTheGood/notebooks/blob/master/HSC%20filters.ipynb

        npt = len(x)

        vecs = [x, y]
        for i, v in enumerate(vecs):
            v_o = np.empty(2*npt - 1)
            vecs[i] = v_o

            v_o[0:npt] = v[::-1]
            if i == 0:                      # radius
                v_o[0:npt] *= -1
            v_o[npt:] = v[1:]

        return scipy.interpolate.interp1d(vecs[0], vecs[1], kind='cubic')

    def _setQEData(self):
        # this was taken from https://www.naoj.org/Observing/Instruments/HSC/txt/qe_ccd_HSC.txt
        # wavelength units are angstroms natively ... change to nm
        # note that I added a cushion, extrapolated to zero, outside the range
        # Better than the alternative?
        self.qeLambda = np.array([3400., 3600., 3800., 4000., 4200., 4400., 4600., 4800.,
                                  5000., 5200., 5400., 5600., 5800., 6000., 6200.,
                                  6400., 6600., 6800., 7000., 7200., 7400., 7600.,
                                  7800., 8000., 8200., 8400., 8600., 8800., 9000.,
                                  9200., 9400., 9600., 9800., 10000., 10200., 10400.,
                                  10600., 10800.0])/10.
        self.qeQe = np.array([0.000, 0.285, 0.36, 0.517, 0.647, 0.748, 0.818, 0.864, 0.91,
                              0.936, 0.948, 0.951, 0.958, 0.954, 0.956, 0.95, 0.938,
                              0.929, 0.916, 0.908, 0.899, 0.892, 0.882, 0.863, 0.856,
                              0.829, 0.825, 0.821, 0.78, 0.745, 0.692, 0.609, 0.517,
                              0.378, 0.209, 0.088, 0.047, 0.0])

    def _setMirrorData(self):
        # this was taken from https://www.naoj.org/Observing/Telescope/Parameters/Reflectivity/M1-2010s.txt
        #  wavelength units are nanometers
        self.mirrorLambda = np.array([300., 305., 310., 315., 320., 325., 330., 335.,
                                      340., 345., 350., 355., 360., 365., 370., 375.,
                                      380., 385., 390., 395., 400., 405., 410., 415.,
                                      420., 425., 430., 435., 440., 445., 450., 455.,
                                      460., 465., 470., 475., 480., 485., 490., 495.,
                                      500., 505., 510., 515., 520., 525., 530., 535.,
                                      540., 545., 550., 555., 560., 565., 570., 575.,
                                      580., 585., 590., 595., 600., 605., 610., 615.,
                                      620., 625., 630., 635., 640., 645., 650., 655.,
                                      660., 665., 670., 675., 680., 685., 690., 695.,
                                      700., 705., 710., 715., 720., 725., 730., 735.,
                                      740., 745., 750., 755., 760., 765., 770., 775.,
                                      780., 785., 790., 795., 800., 805., 810., 815.,
                                      820., 825., 830., 835., 840., 845., 850., 855.,
                                      860., 865., 870., 875., 880., 885., 890., 895.,
                                      900., 905., 910., 915., 920., 925., 930., 935.,
                                      940., 945., 950., 955., 960., 965., 970., 975.,
                                      980., 985., 990., 995., 1000., 1005., 1010., 1015.,
                                      1020., 1025., 1030., 1035., 1040., 1045., 1050., 1055.,
                                      1060., 1065., 1070., 1075., 1080., 1085., 1090., 1095.,
                                      1100., 1105., 1110., 1115., 1120., 1125., 1130., 1135.,
                                      1140., 1145., 1150., 1155., 1160., 1165., 1170., 1175.,
                                      1180., 1185., 1190., 1195.])
        self.mirrorReflectivity = np.array([0.894, 0.8976, 0.8992, 0.902, 0.9028, 0.904, 0.9037,
                                            0.903, 0.9034, 0.9048, 0.9067, 0.9081, 0.9086, 0.9077,
                                            0.9077, 0.9085, 0.9079, 0.9089, 0.9085, 0.9099, 0.9104,
                                            0.9101, 0.91, 0.9098, 0.9091, 0.9094, 0.9098, 0.9099,
                                            0.9099, 0.9098, 0.9096, 0.9096, 0.9092, 0.9093, 0.9093,
                                            0.9093, 0.9087, 0.9084, 0.9078, 0.9074, 0.9072, 0.9071,
                                            0.907, 0.907, 0.9067, 0.9062, 0.9061, 0.9058, 0.9054,
                                            0.9054, 0.9047, 0.9043, 0.9039, 0.9034, 0.9031, 0.9026,
                                            0.9024, 0.9018, 0.9014, 0.9006, 0.9003, 0.8997, 0.8993,
                                            0.899, 0.8982, 0.8977, 0.8971, 0.8963, 0.8956, 0.8948,
                                            0.8942, 0.893, 0.8922, 0.8915, 0.8906, 0.8896, 0.8888,
                                            0.8879, 0.8869, 0.8856, 0.8845, 0.8832, 0.882, 0.8807,
                                            0.8795, 0.8782, 0.8771, 0.8765, 0.8746, 0.8734, 0.8719,
                                            0.8698, 0.868, 0.866, 0.8644, 0.8634, 0.8611, 0.8595,
                                            0.8573, 0.8555, 0.8544, 0.8534, 0.8529, 0.8533, 0.853,
                                            0.8524, 0.8525, 0.8542, 0.8575, 0.8628, 0.8675, 0.8732,
                                            0.878, 0.8821, 0.8853, 0.887, 0.8895, 0.8929, 0.8961,
                                            0.8995, 0.9024, 0.9052, 0.9079, 0.9107, 0.9132, 0.9157,
                                            0.9185, 0.921, 0.9236, 0.9254, 0.9271, 0.9285, 0.9299,
                                            0.9315, 0.9329, 0.9343, 0.9355, 0.9369, 0.938, 0.939,
                                            0.94, 0.941, 0.9422, 0.9429, 0.9438, 0.9446, 0.9454,
                                            0.9462, 0.947, 0.9478, 0.9485, 0.9491, 0.9497, 0.9503,
                                            0.9511, 0.9517, 0.9523, 0.9527, 0.9532, 0.9538, 0.9542,
                                            0.9547, 0.9553, 0.9557, 0.9561, 0.9565, 0.9568, 0.9571,
                                            0.9575, 0.958, 0.9584, 0.9587, 0.9591, 0.9594, 0.9596,
                                            0.9599, 0.9603, 0.9605, 0.9606, 0.9608])

    def _setFilterData(self):
        # this was taken from
        # https://github.com/RobertLuptonTheGood/notebooks/blob/master/HSC%20filters.ipynb
        # wavelength units are nanometers

        self.filterData = dict(
            # radius in mm
            g=dict(radius=np.array([0.00, 50.00, 100.00, 150.00, 200.00, 250.00, 270.00]),
                   EW=np.array([143.35, 143.00, 142.90, 142.73, 142.57, 142.33, 142.22]),
                   lambda_bar=np.array([473.07, 473.04, 472.94, 472.78, 472.56, 472.10, 471.99]),
                   peak=np.array([98.18, 97.95, 98.17, 98.09, 98.05, 98.00, 97.95]),
                   on50=np.array([399.65, 399.65, 399.63, 399.56, 399.40, 399.04, 398.96]),
                   off50=np.array([546.79, 546.76, 546.57, 546.34, 546.05, 545.51, 545.38]),
                   on10=np.array([397.02, 397.03, 397.03, 396.97, 396.76, 396.39, 396.30]),
                   off10=np.array([550.84, 550.81, 550.61, 550.30, 549.92, 549.29, 549.08]),
                   on80=np.array([401.02, 401.01, 400.98, 400.92, 400.79, 400.47, 400.38]),
                   off80=np.array([544.43, 544.42, 544.29, 544.15, 543.92, 543.42, 543.30]),
                   Tmin=np.array([95.98, 95.83, 95.78, 95.81, 95.78, 95.60, 95.65]),
                   Tmax=np.array([98.18, 97.95, 98.17, 98.09, 98.05, 98.00, 97.95]),
                   Tavg=np.array([97.14, 96.92, 96.97, 96.96, 96.94, 96.89, 96.85])),
            r=dict(radius=np.array([0.00, 50.00, 100.00, 150.00, 200.00, 250.00, 270.00]),
                   EW=np.array([144.90, 138.90, 135.08, 134.15, 134.40, 135.38, 133.72]),
                   lambda_bar=np.array([620.59, 623.80, 621.14, 622.74, 623.21, 621.77, 623.72]),
                   peak=np.array([96.47, 94.66, 92.39, 91.74, 92.12, 93.19, 91.75]),
                   on50=np.array([542.38, 547.38, 545.37, 547.02, 547.74, 546.63, 548.32]),
                   off50=np.array([697.74, 699.44, 696.49, 697.87, 697.97, 695.92, 697.99]),
                   on10=np.array([537.35, 542.44, 540.43, 542.05, 542.79, 541.69, 543.26]),
                   off10=np.array([702.70, 704.74, 701.74, 703.18, 703.26, 701.16, 703.56]),
                   on80=np.array([546.61, 553.00, 550.31, 551.93, 552.58, 551.37, 553.48]),
                   off80=np.array([694.47, 695.88, 692.84, 694.26, 694.34, 692.39, 694.13]),
                   Tmin=np.array([88.35, 86.91, 84.96, 83.78, 84.31, 85.60, 85.00]),
                   Tmax=np.array([96.25, 94.41, 91.96, 91.74, 92.12, 93.19, 91.75]),
                   Tavg=np.array([93.79, 91.61, 89.39, 89.18, 89.57, 90.67, 89.28])),
            r1=dict(radius=np.array([0.00, 10.00, 20.00, 30.00, 40.00, 50.00, 60.00, 70.00,
                                     80.00, 90.00, 100.00, 110.00, 120.00, 130.00, 140.00, 150.00,
                                     160.00, 170.00, 180.00, 190.00, 200.00, 210.00, 220.00, 230.00,
                                     240.00, 250.00, 260.00, 267.00]),
                    peak=np.array([99.20, 99.68, 99.23, 98.71, 98.92, 98.26, 98.36, 98.23,
                                   97.61, 97.33, 96.51, 96.83, 96.26, 96.18, 96.09, 95.91,
                                   95.87, 95.79, 95.92, 95.73, 96.15, 96.22, 96.67, 96.59,
                                   96.86, 97.06, 96.74, 96.30]),
                    on50=np.array([542.48, 543.30, 546.23, 549.77, 550.49, 549.50, 548.16, 547.44,
                                   547.54, 548.06, 548.63, 549.14, 549.78, 550.27, 550.20, 549.69,
                                   549.17, 548.80, 548.46, 548.25, 548.30, 548.55, 549.08, 549.72,
                                   550.35, 550.78, 551.04, 551.36]),
                    off50=np.array([680.95, 681.41, 683.75, 686.95, 687.74, 686.96, 685.99, 685.54,
                                    685.65, 686.17, 687.10, 688.20, 689.33, 690.24, 690.55, 690.35,
                                    689.99, 689.54, 689.14, 688.91, 688.99, 689.45, 690.22, 691.00,
                                    691.72, 692.09, 692.29, 692.85]),
                    on10=np.array([539.33, 540.09, 541.45, 545.03, 546.90, 545.83, 544.64, 544.22,
                                   544.36, 544.91, 545.41, 546.02, 546.56, 547.15, 547.17, 546.71,
                                   546.26, 545.90, 545.54, 545.36, 545.41, 545.71, 546.20, 546.84,
                                   547.45, 548.00, 548.24, 548.53]),
                    off10=np.array([685.20, 685.55, 688.95, 691.78, 692.32, 691.51, 690.37, 689.71,
                                    689.76, 690.21, 691.04, 692.08, 693.28, 694.18, 694.50, 694.23,
                                    693.82, 693.45, 692.90, 692.65, 692.74, 693.22, 693.97, 694.81,
                                    695.57, 695.97, 696.07, 696.62]),
                    on80=np.array([551.63, 551.87, 552.05, 553.07, 553.51, 552.33, 550.77, 554.03,
                                   549.70, 550.31, 550.84, 551.36, 551.88, 552.45, 557.52, 557.46,
                                   556.96, 556.63, 556.26, 556.23, 556.59, 557.08, 558.02, 558.74,
                                   559.54, 560.09, 560.37, 560.68]),
                    off80=np.array([677.34, 677.87, 679.95, 683.19, 684.12, 683.50, 682.58, 682.22,
                                    682.38, 682.84, 683.61, 684.66, 685.75, 686.75, 687.01, 686.82,
                                    686.43, 686.07, 685.63, 685.38, 685.35, 685.79, 686.32, 687.27,
                                    687.93, 688.17, 688.24, 688.32]),
                    Tmin=np.array([0.93, 0.93, 0.91, 0.93, 0.93, 0.92, 0.93, 0.92,
                                   0.92, 0.92, 0.91, 0.90, 0.90, 0.90, 0.90, 0.90,
                                   0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.91, 0.91,
                                   0.91, 0.91, 0.90, 0.90]),
                    Tmax=np.array([0.99, 1.00, 0.99, 0.99, 0.99, 0.98, 0.98, 0.98,
                                   0.98, 0.97, 0.97, 0.97, 0.96, 0.96, 0.96, 0.96,
                                   0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.97, 0.97,
                                   0.97, 0.97, 0.97, 0.96]),
                    Tavg=np.array([0.96, 0.97, 0.96, 0.96, 0.96, 0.96, 0.96, 0.95,
                                   0.95, 0.95, 0.94, 0.94, 0.94, 0.94, 0.94, 0.93,
                                   0.93, 0.93, 0.93, 0.93, 0.93, 0.94, 0.94, 0.94,
                                   0.94, 0.94, 0.94, 0.94]),
                    EW=np.array([133.29, 133.50, 132.65, 132.12, 131.86, 131.65, 131.74, 131.60,
                                 131.25, 130.73, 130.34, 130.69, 130.90, 131.15, 131.40, 131.33,
                                 131.62, 131.51, 131.38, 131.38, 131.39, 132.07, 132.60, 132.85,
                                 133.07, 133.10, 132.78, 132.35]),
                    lambda_bar=np.array([614.49, 614.87, 616.00, 618.13, 618.82, 617.91, 616.67, 618.12,
                                         616.04, 616.58, 617.23, 618.01, 618.82, 619.60, 622.26, 622.14,
                                         621.69, 621.35, 620.94, 620.81, 620.97, 621.43, 622.17, 623.00,
                                         623.73, 624.13, 624.31, 624.50]),
                    ),
            i=dict(radius=np.array([0.00, 50.00, 100.00, 150.00, 200.00, 250.00, 270.00]),
                   EW=np.array([141.50, 142.01, 143.25, 144.59, 146.72, 148.23, 147.72]),
                   lambda_bar=np.array([771.10, 771.81, 771.80, 770.18, 769.23, 773.33, 774.76]),
                   peak=np.array([95.38, 95.09, 95.00, 94.81, 94.85, 94.61, 93.66]),
                   on50=np.array([696.47, 696.74, 696.01, 693.52, 691.35, 694.64, 696.20]),
                   off50=np.array([844.27, 845.38, 845.84, 845.46, 845.48, 850.22, 852.02]),
                   on10=np.array([688.47, 688.75, 688.03, 685.46, 683.38, 686.58, 688.15]),
                   off10=np.array([857.11, 857.91, 858.34, 857.77, 857.75, 864.25, 865.20]),
                   on80=np.array([701.04, 701.27, 700.50, 697.96, 695.72, 699.06, 700.65]),
                   off80=np.array([839.78, 840.75, 841.29, 840.85, 841.01, 844.89, 845.35]),
                   Tmin=np.array([92.99, 92.59, 92.41, 91.99, 91.93, 92.10, 91.07]),
                   Tmax=np.array([95.33, 94.91, 95.00, 94.56, 94.36, 94.57, 93.66]),
                   Tavg=np.array([94.25, 94.02, 93.93, 93.68, 93.62, 93.68, 92.69])),
            # FIXME: this is a placeholder until I have i2, and is just the middle of the i and flat
            i2=dict(radius=np.array([0.00, 50.00, 100.00, 150.00, 200.00, 250.00, 270.00]),
                    EW=np.array([144.59, 144.59, 144.59, 144.59, 144.59, 144.59, 144.59]),
                    lambda_bar=np.array([770.18, 770.18, 770.18, 770.18, 770.18, 770.18, 770.18]),
                    peak=np.array([94.81, 94.81, 94.81, 94.81, 94.81, 94.81, 94.81]),
                    on50=np.array([693.52, 693.52, 693.52, 693.52, 693.52, 693.52, 693.52]),
                    off50=np.array([845.46, 845.46, 845.46, 845.46, 845.46, 845.46, 845.46]),
                    on10=np.array([685.46, 685.46, 685.46, 685.46, 685.46, 685.46, 685.46]),
                    off10=np.array([857.77, 857.77, 857.77, 857.77, 857.77, 857.77, 857.77]),
                    on80=np.array([697.96, 697.96, 697.96, 697.96, 697.96, 697.96, 697.96]),
                    off80=np.array([840.85, 840.85, 840.85, 840.85, 840.85, 840.85, 840.85]),
                    Tmin=np.array([91.99, 91.99, 91.99, 91.99, 91.99, 91.99, 91.99]),
                    Tmax=np.array([94.56, 94.56, 94.56, 94.56, 94.56, 94.56, 94.56]),
                    Tavg=np.array([93.68, 93.68, 93.68, 93.68, 93.68, 93.68, 93.68])),
            z=dict(radius=np.array([0.00, 50.00, 100.00, 150.00, 200.00, 250.00, 270.00]),
                   EW=np.array([80.01, 80.12, 80.14, 79.37, 78.82, 80.41, 80.26]),
                   lambda_bar=np.array([892.63, 892.33, 892.18, 891.72, 891.99, 891.25, 891.45]),
                   peak=np.array([98.53, 98.82, 99.17, 99.05, 98.59, 99.26, 99.66]),
                   on50=np.array([852.97, 852.75, 852.75, 852.66, 852.77, 852.37, 852.36]),
                   off50=np.array([932.22, 931.92, 931.67, 930.89, 930.99, 930.86, 930.81]),
                   on10=np.array([844.28, 844.11, 844.25, 843.83, 844.40, 843.20, 843.54]),
                   off10=np.array([940.81, 940.58, 940.30, 939.53, 939.57, 939.62, 939.50]),
                   on80=np.array([858.27, 857.95, 858.00, 858.02, 858.10, 857.70, 857.58]),
                   off80=np.array([926.58, 926.29, 926.08, 925.21, 925.36, 925.12, 925.09]),
                   Tmin=np.array([97.07, 96.86, 96.91, 97.21, 97.00, 97.57, 97.78]),
                   Tmax=np.array([98.53, 98.82, 99.17, 99.05, 98.59, 99.26, 99.66]),
                   Tavg=np.array([97.97, 98.07, 98.32, 98.34, 98.03, 98.63, 98.98])),
            y=dict(radius=np.array([0.00, 50.00, 100.00, 150.00, 200.00, 250.00, 270.00]),
                   EW=np.array([141.37, 141.34, 139.06, 141.19, 139.50, 143.70, 144.91]),
                   lambda_bar=np.array([1002.77, 1002.98, 1004.36, 1003.56, 1003.99, 1003.19, 1001.92]),
                   peak=np.array([97.92, 98.21, 97.74, 98.36, 97.59, 98.18, 97.86]),
                   on50=np.array([933.55, 934.12, 935.81, 935.33, 934.85, 932.85, 930.31]),
                   off50=np.array([1072.36, 1072.35, 1072.90, 1072.62, 1072.86, 1073.84, 1073.22]),
                   on10=np.array([923.65, 924.03, 926.13, 925.00, 925.35, 922.85, 920.70]),
                   off10=np.array([1082.85, 1082.86, 1083.39, 1083.36, 1083.38, 1084.44, 1083.83]),
                   on80=np.array([939.26, 939.65, 941.54, 941.23, 941.02, 939.10, 936.69]),
                   off80=np.array([1065.99, 1065.89, 1066.40, 1066.11, 1066.32, 1067.11, 1066.53]),
                   Tmin=np.array([96.25, 96.64, 96.20, 97.17, 96.06, 96.63, 96.31]),
                   Tmax=np.array([97.92, 98.21, 97.74, 98.36, 97.59, 98.18, 97.86]),
                   Tavg=np.array([97.13, 97.42, 97.03, 97.84, 96.85, 97.47, 97.25])),)
