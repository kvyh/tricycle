import periodicity as pr
import data

def main():
    allkics = data.select_kics()
    files = pr.Periodicity(allkics[:20])
    #files.autocorrelate_periods(allkics[:80], write = True)
    #files.make_period_files()
    #periods = files.refined_periods(files.pullfile('autocor_all'), files.pullfile('periodogram_all'), correlation = (.95,1.05))
    files.plot_all(filename = 'periodogram_periods_all_interp_0.8')
    #files.lightcurve(allkics[47])
    #files.plot_periodogram(allkics[47], autoccor= True)
    
    
if __name__ == '__main__':
    main()
