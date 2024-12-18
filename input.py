from VSCFintensity import intensity_plotter 

plotter = intensity_plotter(filepath='+1/DLF-Raman.dat',
                            plotname='Raman Spectra',
                            outputdir='./',
                            unit='a.u.',
                            bandwidth=4.0,
                            curvewidth=1.0,
                            raman=True,    # plot Raman spectra
                            gauss=False,
                            ranged=True,
                            freqinf=100.0,
                            freqsup=4000.0,
                            band=True,
                            subtract=False, # subtract background
                            bkgdfile='',
                            margin=0.0,     # frequency margin, this is a factor*range gives real margin
                            hide=True,      # hide intensity values
                            ylimf=1.0,     # scaling factor for y limit, max int val * ylimf for ylimf
                            nspectra=3,
#                            datapaths=['_local/vibrational.spectra.vscf.dat', '_normal/vibrational.spectra.vscf.dat', '_harmonic/DLF-IR.dat'],
                            datapaths=['_local/vibrational.spectra.vscf.dat', '_normal/vibrational.spectra.vscf.dat', '_harmonic/DLF-Raman.dat'],
                            yshift=0.30,
#                            uniintmax=[200.0, 200.0, 200.0],   # universal max intensity set by user - force max int val = uniintmax
                            uniintmax=[100000.0, 120000.0, 120000.0], # Raman
#                            scalef=[200.0, 200.0, 200.0], # intensity scaling factor, combining with 'uniintmax' to fit in y axis
                            scalef=[100000.0, 120000.0, 120000.0], # Raman
                            colours=['c', 'g', 'b'], # colour schemes: g:, b--, darkorange, seagreen, see: 'https://matplotlib.org/stable/gallery/color/named_colors.html'
                            labels=['VSCF local modes', 'VSCF normal modes', 'harmonic',  ],
                            legendpos='upper right',
                            vscfspectra=[True, True, False], # if spectra file is by VSCF or not
                            vscfnmodes=45,
)

plotter.runmain()
