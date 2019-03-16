# Convert text line files into numpy npy files
import numpy as np
import os.path
GHz = 29.9792458     # conversion from cm^-1 to GHz


def convert_nh3_dbs():
    # def readInputFiles(path, im=0, rm=0, vm=0, verbose=False)
    # readInputFiles(path,0,1100,4000,verbose=verbose)
    # global fo, Io, Eo, gammaNH3o, H2HeBroad
    # global fo_rot, Io_rot, Eo_rot, gNH3_rot, gH2_rot, gHe_rot
    # global fo_v2, Io_v2, Eo_v2

    # %% Inversion lines:
    # % fo is frequency in GHz, Io is line intensity in cm^-1/(molecule./cm^2),
    # % Eo is lower state energy in cm^-1, gammaNH3o and H2HeBroad are self and
    # % foreign gas broadening parameters.
    path = 'nh3'
    im = 0
    rm = 1100
    vm = 4000
    filename = os.path.join(path, 'ammonia_inversion.dat')
    fo, Io, Eo, gammaNH3o, H2HeBroad = np.loadtxt(filename, skiprows=1, unpack=True)
    if im > 0:
        fo = fo[:-im]
        Io = Io[:-im]
        Eo = Eo[:-im]
        gammaNH3o = gammaNH3o[:-im]
        H2HeBroad = H2HeBroad[:-im]

    # %% Rotational lines:
    # % fo_rot is frequency in GHz, Io_rot is line intensity in
    # % cm^-1/(molecule./cm^2), Eo_rot is lower state energy in cm^-1, gNH3_rot,
    # % gH2_rot, gHe_rot are broadening parameters for rotational lines.
    filename = os.path.join(path, 'ammonia_rotational.dat')
    fo_rot, Io_rot, Eo_rot, gNH3_rot, gH2_rot, gHe_rot = np.loadtxt(filename, skiprows=1, unpack=True)
    if rm > 0:
        fo_rot = fo_rot[:-rm]
        Io_rot = Io_rot[:-rm]
        Eo_rot = Eo_rot[:-rm]
        gNH3_rot = gNH3_rot[:-rm]
        gH2_rot = gH2_rot[:-rm]
        gHe_rot = gHe_rot[:-rm]

    # %% v2 roto-vibrational lines:
    # % fo_v2 is frequency in GHz, Io_v2 is line intensity in
    # % cm^-1/(molecule./cm^2), Eo_v2 is lower state energy in cm^-1,
    filename = os.path.join(path, 'ammonia_rotovibrational.dat')
    fo_v2, Io_v2, Eo_v2 = np.loadtxt(filename, skiprows=1, unpack=True)
    if vm > 0:
        fo_v2 = fo_v2[:-vm]
        Io_v2 = Io_v2[:-vm]
        Eo_v2 = Eo_v2[:-vm]

    with open('ammonia.npz', 'wb') as ofp:
        np.savez(ofp, fo=fo, Io=Io, Eo=Eo, gammaNH3o=gammaNH3o, H2HeBroad=H2HeBroad, fo_rot=fo_rot, Io_rot=Io_rot, Eo_rot=Eo_rot,
                 gNH3_rot=gNH3_rot, gH2_rot=gH2_rot, gHe_rot=gHe_rot, fo_v2=fo_v2, Io_v2=Io_v2, Eo_v2=Eo_v2)


def convert_fieg(constituent):
    f0 = []
    I0 = []
    E = []
    Gxx = []
    nlin = 0
    if constituent == 'h2s':
        useLinesUpTo = 200
        I0_scale = 1.0
        inputfile = 'h2s/h2s.lin'
        use_these = ['f0', 'I0', 'E', 'Gxx']
    elif constituent == 'ph3':
        useLinesUpTo = 320
        I0_scale = GHz * 1.0E17
        inputfile = 'ph3/ph3jh.lin'
        use_these = ['f0', 'I0', 'E']
    elif constituent == 'nh3_sjs':
        useLinesUpTo = 200
        I0_scale = 1.0
        inputfile = 'nh3/nh3.lin'
        use_these = ['f0', 'I0', 'E', 'Gxx']
    elif constituent == 'co':
        useLinesUpTo = 26
        I0_scale = 1.0
        inputfile = 'co/co.lin'
        use_these = ['f0', 'I0', 'E']
    outputfile = inputfile.split('.')[0] + '.npz'
    with open(inputfile, 'r') as ifp:
        for line in ifp:
            if nlin >= useLinesUpTo:
                break
            nlin += 1
            data = line.split()
            if len(data) == len(use_these):
                if 'f0' in use_these:
                    f0.append(float(data[0]))
                if 'I0' in use_these:
                    I0.append(float(data[1]) / I0_scale)
                if 'E' in use_these:
                    E.append(float(data[2]))
                if 'Gxx' in use_these:
                    Gxx.append(float(data[3]))
            else:
                break
    f0 = np.array(f0)
    I0 = np.array(I0)
    E = np.array(E)
    Gxx = np.array(Gxx)
    # ## Be sure to change this  --------
    #                                    |
    #                                    |
    #                                    V
    with open(outputfile, 'wb') as ofp:
        np.savez(ofp, f0=f0, I0=I0, E=E, G0=Gxx)


def convert_wgt(constituent):
    inputfile = 'ph3/PH3WGT.dat'
    useLinesUpTo = 40
    outputfile = inputfile.split('.')[0] + '.npz'
    WgtI0 = []
    WgtFGB = []
    WgtSB = []
    nwgt = 0
    with open(inputfile, 'r') as ifp:
        for line in ifp:
            if nwgt >= useLinesUpTo:
                break
            nwgt += 1
            data = line.split()
            if len(data) == 4:
                WgtI0.append(float(data[1]))
                WgtFGB.append(float(data[2]))
                WgtSB.append(float(data[3]))
            else:
                break
    pad_ones = np.ones(320 - 40)
    WgtI0 = np.concatenate((np.array(WgtI0), pad_ones))
    WgtFGB = np.concatenate((np.array(WgtFGB), pad_ones))
    WgtSB = np.concatenate((np.array(WgtSB), pad_ones))
    with open(outputfile, 'wb') as ofp:
        np.savez(ofp, WgtI0=WgtI0, WgtFGB=WgtFGB, WgtSB=WgtSB)
