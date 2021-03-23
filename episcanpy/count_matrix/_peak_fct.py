# peak loading and normalising matrices
import pandas as pd 

def load_peaks(input_file, sep='\t', header=None):
    """
    Load standard output of macs2 narrowPeak.
    """
    df = pd.read_csv(input_file, sep=sep, header=header)
    dict_peaks = {x.lstrip('chr'):[] for x in set(df[0].tolist())}
    for line in df.values:
        #dict_peaks[line[0]].append([line[1], line[2], "_".join([str(x) for x in line[0:3]])])
        dict_peaks[line[0].lstrip('chr')].append([line[1], line[2]])
    return(dict_peaks)

def norm_peaks(peaks, extension=150):
    """
    takes the middle of the peak and extend it of the number of base pair corresponding
    to the extension parameter(150bp) on both side of the peak.
    """ 
    for chrom, peak_chrom in peaks.items():
        for indiv_peak in peak_chrom:
            average_start = (indiv_peak[0]+indiv_peak[1])//2
            indiv_peak[0] = average_start-extension
            indiv_peak[1] = average_start+extension
            
