General Remarks:
================

All spectrum_*.txt consist of two columns with the
following contents:
1st column: Energy in keV
2nd column: number of photons per mm^2 at 1m distance with
   an energy around the one given in column 1; these numbers
   can be given as photon density values (per keV) or as
   number of photons per bin (e.g. per 27keV)


keV spectra:
============

spectrum_73keVp_RQA5_xraylib_0.1keV_bins.txt:
The number of photons is given per mAs, per 0.1keV bin.
Parameters used when calling
xrspec(E, kV, target_material, mAs, mdis, degrees, mmpyrex, mmoil, mmlexan, mmAl or mmMo):
xrspec(E, 73, 'W',             1,   1,    12,      2.38,    3.06,  2.66,    1.2 + 22.5)

spectrum_90keVp_xraylib_1keV_bins.txt:
The number of photons is given per mAs, per 1keV bin.
Parameters used when calling
xrspec(E, kV, target_material, mAs, mdis, degrees, mmpyrex, mmoil, mmlexan, mmAl or mmMo):
xrspec(E, 90, 'W',             1,   1,    12,      2.38,    3.06,  2.66,    1.2)

spectrum_120keVp_xraylib_1keV_bins.txt:
The number of photons is given per mAs, per 1keV bin.
Parameters used when calling
xrspec(E, kV, target_material, mAs, mdis, degrees, mmpyrex, mmoil, mmlexan, mmAl or mmMo):
xrspec(E, 120, 'W',             1,   1,    12,      2.38,    3.06,  2.66,    1.2)

spectrum_120keVp_1mmAl_spektr_1keV_bins.txt:
Siewerdsen's Spektr (Med. Phys., 31(11):3057-3067, 2004)
uses the tungsten anode spectra model interpolating
polynomals (tasmip) algorithm.  (See also the spektr_v*
directory in the matlab libraries.)
The number of photons is given per mAs, per 1keV bin.
Parameters used when calling
spektrSpectrum(kVp, mmAl):
spektrSpectrum(120, 1)

spectrum_120keVp_10mmAl_spektr_1keV_bins.txt:
Siewerdsen's Spektr (Med. Phys., 31(11):3057-3067, 2004)
uses the tungsten anode spectra model interpolating
polynomals (tasmip) algorithm.  (See also the spektr_v*
directory in the matlab libraries.)
The number of photons is given per mAs, per 1keV bin.
Parameters used when calling
spektrSpectrum(kVp, mmAl):
spektrSpectrum(120, 10)

Lei_120kVp_10mmAl.txt:
This is a deprecated spectrum originally used by Lei Zhu
which is only kept for reference.  Please refrain from using
it!


MeV spectra:
============

spectrum_2.5MeVp_AlTarget_NoFilter_Varian_27keV_bins.txt:
The number of photons is given per mm^2 (converted from the
original values per cm^2 values), per 27keV bin, and per
second.  (Since the source is pulsed with 360 Hz, the number
of photons per pulse can be obtained by dividing by 360.)
Both data columns received an additional digit compared to
the original Excel file.

spectrum_6MeVp_CuTarget_NoFilter_Varian_62keV_bins.txt:
The number of photons is given per mm^2 (converted from the
original values per cm^2 values), per 62keV bin, and per
second.  (Since the source is pulsed with 360 Hz, the number
of photons per pulse can be obtained by dividing by 360.)
Both data columns received an additional digit compared to
the original Excel file.

spectrum_6MeVp_WCuTarget_ClinacFilter_Varian_62keV_bins.txt:
The number of photons is given per mm^2 (converted from the
original values per cm^2 values), per 62keV bin, and per
second.  (Since the source is pulsed with 360 Hz, the number
of photons per pulse can be obtained by dividing by 360.)
Both data columns received an additional digit compared to
the original Excel file.

Lei_6MVp.txt:
This is a deprecated spectrum originally used by Lei Zhu
which is only kept for reference.  Please refrain from using
it!
