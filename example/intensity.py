from ctypes import Structure
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from os.path import join
from sys import argv
from os  import getcwd
from math import exp


font = {'fontname':'Times New Roman'}
rc('font', **{'family':'Times New Roman', 'serif':['Times New Roman']})
rc('text', usetex=True)
aspect = 0.25
bandwidth = 1
aspect = 0.25

class intensity_plotter(Structure):

    _params = {'filepath': '',
               'plotname': '',
               'outputdir': './',
               'unit': 'a.u.',
               'bandwidth': 2.0, 
               'raman': False,
               'guass': False,
               'ranged': False,
               'freqinf': 0.0,
               'freqsup': 5000.0,
               'band': True,
               'subtract': False,
               'bkgdfile': '',
               'uniintmax': [1.0], # universal max intensity value, set by user
               'margin': 0.0,
               'hide': False,
               'scalef': [1.0],
               'scaleffreq': [1.0],
#               'scalefint': 1.0,
               'ylimf': 1.0,
               'nspectra': 1,
               'datapaths': [],
               'yshift': 0.05,
               'colours': ['b', 'g'],
               'labels': [],
               'legendpos': 'upper right',
               'legendftsize': 'x-large',
               'reverse': False,
              }

    # NOTE: no __init__() method
    x = np.linspace(0, 5000, 100000)
    y = None
    low_id = 0
    up_id  = 0
    modeids     = None
    # wavenumbers as x axis
    wavenumbers = None
    # (IR) intensities as y axis
    intensities = None
    intensities_unit1 = None
    # need to scale the y axis
    # max value of intensities scaled to 1.0
    # for example 0.2 to 1.0 - factor of 5.0
    maxint_value = 0.0 # max intensity value
    mean_int     = False
    count_spectra= 0

    def get_params(self):

        for key, val in self._params.items(): 
            try:
                getattr(self, key)
                print(key, val)
            except:
                setattr(self, key, val)
                print("set", key, val)


    # No longer in use
    def get_argvs(self):

        self.filepath     = argv[1] # file path of IR or Raman file
        self.plotname     = argv[2]
        self.outputdir    = argv[3]
        self.unit         = argv[4]
        self.raman        = bool(int(argv[5])) # if true, plot Raman signatures
        self.gauss   = bool(int(argv[6])) # if true, plot guassians of signatures # default is lorenzian
        self.ranged  = bool(int(argv[7])) # if true, select signatures of frequency between (300, 5000); rule out low frequency signatures
        self.band    = bool(int(argv[8])) # if true, plot vertical band lines
        self.freqinf     = float(argv[9]) 
        self.freqsup      = float(argv[10])     
        self.subtract = bool(int(argv[11]))
        self.bkgdfile = argv[12]
        self.uniintmax = float(argv[13])
        self.margin  = float(argv[14])
        self.hide = bool(int(argv[15]))
        # common scaling factor to replace mean intensity
        # for comparing spectra
        try:
             self.scalef= float(argv[16])
             if self.scalef == 0.0:
                  self.scalef=1.0
        except:
             self.scalef= 1.0
        try:
             self.ylimf = float(argv[17])
             if self.ylimf == 0.0:
                 self.ylimf = 1.0
        except:
             self.ylimf = 1.0

        try:
            self.nspectra    = int(argv[18])
            # if number of spectra is more than 1,
            # then this program will make number of spectra
        except:
            self.nspectra = 1
        print("number of spectra to plot: {}".format(self.nspectra))

        self.yshift = 0.2

        self.datapaths = ['+1/DLF-Raman.dat', '0/DLF-Raman.dat']

        print(argv[:])
        print("a raman plotting: {}".format(self.raman))
        print("plotting gaussian: {}".format(self.gauss))
        print("select frequency range: {}".format(self.ranged))
        print("frquency ranges: ({}, {})".format(self.freqinf, self.freqsup))

        if len(self.filepath)  == 0: 
            self.filepath  = getcwd()

        if len(self.outputdir) == 0: 
            self.outputdir = getcwd()
    
    def read_data_plot_multiple(self):
        # read data file names

        self.count_spectra = 0
        # read multiple sets of data
        figname = self.plotname
        unit='a.u.'
        self.filename = None
        for filename in self.datapaths:
            print("current data file name: {}".format(filename)) 
            self.filename = filename
            self.read_data(filename)
            self.plot_int(figname, unit, background=False)
            self.count_spectra += 1

        plt.legend(loc=self.legendpos, 
                   shadow=False, 
                   fontsize=self.legendftsize)

        plt.savefig(figname+'.pdf', dpi=600, format='pdf')
        
    def read_data(self, filename):
        # read data from file
        self.freqsup = self.freqsup
        self.freqinf = self.freqinf

        with open(filename) as fin:
            self.data = np.loadtxt(filename)
            print("dbg only: check data read in")
            print(self.data.T)
            # wavenumbers in cm-1
            self.modeids     = self.data.T[0]
            self.wavenumbers = self.data.T[1]
            self.wavenumbers *= self.scaleffreq[self.count_spectra]
            # get the lower and upper bound indices in the frequencies list
            #print("dbg only: wave numbers")
            print("### dbg only: wns {}".format(self.wavenumbers))
            if self.freqinf != 0.0:
                for n, freq in enumerate(self.wavenumbers):
                    if self.wavenumbers[n] >= self.freqinf:
                        print("found the lower bound of frquency")
                        self.low_id = n
                        break
                    else:
                        self.low_id = 0
                print("final lower freq id: {}".format(self.low_id))
                        
            if self.freqsup != 0.0:
                for n, freq in enumerate(self.wavenumbers):
                    print(self.wavenumbers[n])
                    if self.wavenumbers[n] >= self.freqsup:
                        print("found the upper bound of frquency")
                        self.up_id = n-1
                        break
                    else:
                        self.up_id = n
                print("final upper freq id: {}".format(self.up_id))
            if self.up_id < self.low_id:
                exit("ERROR: cannot find any signature band within the frequency range: {}".format(self.filename))

            if self.raman: # for Raman
                # a.u.
                self.intensities_unit1    = self.data.T[2]
                self.intensities_unit1_tr = self.data.T[3]
                print("\n wavenumbers (cm-1): \n {} \n".format(self.wavenumbers))
                print("Raman intensities (a.u.): \n {} \n".format(self.intensities_unit1_tr))
                # bohr ang is NWChem's default
                self.intensities_unit2    = self.data.T[4]
                # depolar ratios (polarised) arbitrary unit
                self.intensities_unit3    = self.data.T[5]
                self.intensities_unit3_tr = self.data.T[6]
                # depolar ratios (unpolarised) arbitrary unit
                self.intensities_unit4    = self.data.T[7]
                self.intensities_unit4_tr = self.data.T[8]
            else: # for IR
                # atomic unit
                self.intensities_unit1    = self.data.T[2]
                self.intensities_unit1_tr = self.data.T[3]
                print("\n wavenumbers (cm-1): \n {} \n".format(self.wavenumbers))
                print("IR intensities (a.u.): \n {} \n".format(self.intensities_unit1_tr))
                # debye/angstrom
                self.intensities_unit2    = self.data.T[4]
                self.intensities_unit2_tr = self.data.T[5]
                # KM/mol
                self.intensities_unit3    = self.data.T[6]
                self.intensities_unit3_tr = self.data.T[7]
                # averaged frequency (= a.u.*factor(=3*N atoms/sum int)) arbitrary unit
                self.intensities_unit4    = self.data.T[8]
                self.intensities_unit4_tr = self.data.T[9]

    def plot_int(self, figname, unit, background=False):
        #plt.autoscale(False)
        #xlim = plt.gca().get_xlim()
        #plt.xlim(xlim, )

        if self.raman:
            if unit=='' or unit=='a.u.' or unit=='au':
                #unit = 'a.u.'
                unit = 'arb.u.'
                if self.low_id < self.up_id:
                    self.maxint_value = max(self.intensities_unit1_tr[self.low_id: self.up_id])
                else:
                    self.maxint_value = self.intensities_unit1_tr[self.low_id]

                print("max int value: {}".format(self.maxint_value))
                self.int_data         = self.intensities_unit1_tr

            if unit=='bohr2ang' or unit=='ang': # compare with NWChem only, it reports derivatives only
                self.maxint_value = max(self.intensities_unit2   [self.low_id: self.up_id])
                self.int_data     = self.intensities_unit2
            if unit=='polarised' or unit=='polar':
                self.maxint_value = max(self.intensities_unit3_tr[self.low_id: self.up_id])
                self.int_data     = self.intensities_unit3_tr
            if unit=='unpolarised' or unit=='unpolar':
                self.maxint_value = max(self.intensities_unit4_tr[self.low_id: self.up_id])
                self.int_data     = self.intensities_unit4_tr
        else:
            if unit=='' or unit=='a.u.' or unit=='au':
                #unit = 'a.u.'
                unit = 'arb.u.'
                if self.low_id < self.up_id:
                    self.maxint_value = max(self.intensities_unit1_tr[self.low_id: self.up_id])
                else:
                    self.maxint_value = self.intensities_unit1_tr[self.low_id]

                print("max int value: {}".format(self.maxint_value))
                self.int_data     = self.intensities_unit1_tr
            if unit=='debye/ang' or unit=='debyeang':
                self.maxint_value = max(self.intensities_unit2_tr[self.low_id: self.up_id])
                self.int_data     = self.intensities_unit2_tr
            if unit=='KM/mol' or unit=='kmmol':
                self.maxint_value = max(self.intensities_unit3_tr[self.low_id: self.up_id])
                self.int_data     = self.intensities_unit3_tr
            if unit=='arbitrary' or unit=='arb':
                self.maxint_value = max(self.intensities_unit4_tr[self.low_id: self.up_id])
                self.int_data     = self.intensities_unit4_tr

        if self.mean_int:
             self.int_data = self.int_data/self.int_data.mean()
             self.maxint_value = max(self.int_data[self.low_id: self.up_id])

        if self.uniintmax[self.count_spectra] != 0.0:
            self.maxint_value = self.uniintmax[self.count_spectra]

        if self.scalef[self.count_spectra] != 0.0:
             self.int_data = self.int_data/self.scalef[self.count_spectra]
             self.maxint_value = self.maxint_value/self.scalef[self.count_spectra] 

        if not background:
            plt.ylim(0.0, self.maxint_value*self.ylimf)
            if self.reverse:
                plt.xlim(max(self.x), min(self.x))
            else:
                plt.xlim(min(self.x), max(self.x))

            if self.freqinf !=0 or self.freqsup !=0: 

                range_freqs_wn = self.freqsup - self.freqinf

                if self.reverse:
                    plt.xlim(self.freqsup + self.margin*range_freqs_wn,
                             self.freqinf - self.margin*range_freqs_wn)
                else:
                    plt.xlim(self.freqinf - self.margin*range_freqs_wn,
                             self.freqsup + self.margin*range_freqs_wn)


            # get the mode ids by which the frequencies are within the range
            for i, wn in enumerate(self.wavenumbers):

                if self.ranged:

                    if i > self.up_id or i < self.low_id:
                        continue

                if (self.band):

                    if self.nspectra > 1 and self.count_spectra > 0:
                        shift = self.count_spectra * self.yshift
                    else:
                        shift = 0.0

                    scalefactor = 1.0/(self.maxint_value*self.ylimf)

                    plt.axvline(x=wn,
                                ymin=scalefactor*shift, 
                                ymax=scalefactor*(self.int_data[i]+shift), 
                                linewidth=1.0, 
                                color='red')

            if self.hide:
                #get current axes
                ax = plt.gca()
                #hide y-axis
                #ax.get_yaxis().set_visible(False)
                ax.get_yaxis().set_ticks([])
            plt.xticks(fontsize='x-large')
            plt.yticks(fontsize='x-large')

            csfont = {'fontname':'Comic Sans MS'}
            hfont = {'fontname':'Helvetica'}
            tnrfont = {'fontname': 'Times New Roman'}

            #plt.title('title',**csfont)
            #plt.xlabel('xlabel', **hfont)
            #plt.show()
            #plt.rcParams.update({'font.size': 30})
            plt.xlabel(r"$\tilde{\nu}$ (cm$^{-1}$)",**font, fontsize='x-large')
            plt.ylabel(r"$I$ ({})".format(unit),**font, fontsize='x-large')
            plt.title ("{}".format(self.plotname),**font)#.format(self.sorb))
            plt.xticks(**font)
            plt.yticks(**font)

            #ax = plt.gca()
            #ax.set_box_aspect(aspect)
            #f.set_figheight(height, forward=True)
            #plt.tight_layout()

            #plt.title ("Computational IR intensity of {} NH3-H (CHA)".format(self.sorb), fontsize="large")
            #plt.legend()
            #plt.show()
            ###plt.savefig(figname, dpi=600)
            if (self.gauss):
                self.plot_gaussians(figname)
            else:
                self.plot_lorentzians(figname)
        else: # obtaining the interpolation for backgroud intensities
            if (self.gauss):
                self.y = self.bkgd_gaussians()
            else:
                self.y = self.bkgd_lorentzians()

    def single_lorentzian_func(self, x, amp, wid, cen):
        return(amp*wid**2/((x-cen)**2+wid**2))

    def single_gaussian_func  (self, x, amp, wid, cen):
        y = np.zeros(np.shape(x))
        for n, e_x in enumerate(x):
            y[n] = amp*exp(-(e_x-cen)**2/wid**2)
        return y

    def plot_lorentzians(self, figname):
        y = 0.0
        # loop over all peaks
        for i, wn in enumerate(self.wavenumbers):
            if self.ranged:
                if i > self.up_id or i < self.low_id:
                #if wn > 5000.0 or wn < 300.00:
                    continue
            print("### dbg only: self.int_data[i] = {}".format(self.int_data[i]))
            y += self.single_lorentzian_func(self.x, self.int_data[i], self.bandwidth, wn)

            if self.subtract:
                y = y-self.y

        if self.nspectra > 1 and self.count_spectra != 0: # after the first spectrum, shift up
            y += self.count_spectra * self.yshift

        plt.plot(self.x, y, 
                 self.colours[self.count_spectra], 
                 linewidth=1.0, 
                 label=self.labels[self.count_spectra])

        if self.nspectra == 1:
            plt.savefig(figname, dpi=600)

        #print("\n>>> saving Lorentzian to file...\n")
        #plt.ylim(0.0, self.maxint_value*factor_int)
        #plt.xlim(wn_max, wn_min)
        #print(">>> Lorentzian Y max =", plt.ylim())
        #ax = plt.gca()
        #ax.set_box_aspect(aspect)
        #plt.tight_layout()

        
    def plot_gaussians(self, figname):
        y = 0.0
        # loop over all peaks
        for i, wn in enumerate(self.wavenumbers):
            if self.ranged:
                if i > self.up_id or i < self.low_id:
                #if wn > 5000.0 or wn < 300.00:
                    continue
            y += self.single_gaussian_func(self.x, self.int_data[i], self.bandwidth, wn)
            if self.subtract:
                y = y-self.y

        plt.plot(self.x, y, 'g')
        plt.savefig(figname, dpi=600)

    def bkgd_lorentzians(self): # background Lorenzians
        y = 0.0
        # loop over all peaks
        for i, wn in enumerate(self.wavenumbers):
            if self.ranged:
                if i > self.up_id or i < self.low_id:
                #if wn > 5000.0 or wn < 300.00:
                    continue
            y += self.single_lorentzian_func(self.x, self.int_data[i], self.bandwidth, wn)
        return y
 
    def bkgd_gaussians(self): # background Lorenzians
        y = 0.0
        # loop over all peaks
        for i, wn in enumerate(self.wavenumbers):
            if self.ranged:
                if i > self.up_id or i < self.low_id:
                #if wn > 5000.0 or wn < 300.00:
                    continue
            y += self.single_gaussian_func(self.x, self.int_data[i], self.bandwidth, wn)
        return y

    def read_background_int(self):
        self.read_data(self.bkgdfile)

    def runmain(self):

        self.get_params()

        if self.subtract:
            self.read_data("{}".format(self.bkgdfile))
            self.plot_int('', self.unit, background=True) 

        if self.nspectra > 1:
            self.read_data_plot_multiple()
        else:
            self.read_data("{}".format(self.filepath))#.format(plotter.sorb))
            #plotter.read_data("story/H-NH3/{}/DLF-IR.dat".format(plotter.sorb))
            self.plot_int ("{}.eps".format(join(self.outputdir, self.plotname)), self.unit)#.format(plotter.sorb))

###plotter = intensity_plotter()
###plotter.runmain()
#plotter.plot_int ("H-NH3.{}.IR.png".format(plotter.sorb))
