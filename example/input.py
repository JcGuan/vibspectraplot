from intensity import intensity_plotter 

plotter = intensity_plotter(filepath='./CYMV-IR.dat',
                            plotname='IR of a single water in ZSM-5',
                            outputdir='./',
                            unit='a.u.',
                            bandwidth=2.0,
                            raman=True,
                            gauss=False,
                            ranged=True,
                            freqinf=1000.0,
                            freqsup=4000.0,
                            band=True,
                            subtract=False, # subtract background
                            bkgdfile='',
                            margin=.01, # frequency margin, this is a factor*range gives real margin
                            hide=True, # hide intensity values
                            ylimf=10.0,       # scaling factor for y limit, max int val * ylimf for ylimf
                            nspectra=1,
                            datapaths=['./CYMV-IR.dat'],
                            yshift=0.8,
                            uniintmax=[1.0],   # universal max intensity set by user - force max int val = uniintmax
                            scalef=[2.0], # intensity scaling factor, combining with 'uniintmax' to fit in y axis
                            scaleffreq=[1.0],
                            colours=['b'], # colour schemes: g:, b--, darkorange, seagreen, see: 'https://matplotlib.org/stable/gallery/color/named_colors.html'
                            labels=['single water'],
                            reverse=True,
                            legendpos='upper left',
                            
                            
                            
)

plotter.runmain()
