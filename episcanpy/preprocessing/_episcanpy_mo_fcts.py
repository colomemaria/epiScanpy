import pandas as pd
import anndata as ad
def load_gtf_file(input_gtf_file,
                  gtf_column_names=['seqname',
                                    'source',
                                    'feature',
                                    'start',
                                    'end',
                                    'score',
                                    'strand',
                                    'attribute',
                                    'other'],
                  comment='#'):
    """
    Load a gtf file. 
    more information about gtf files can be found on the ensembl website
    here: https://www.ensembl.org/info/website/upload/gff.html
    
    Mouse and human gtf files can be found on the gencode website.
    here: https://www.gencodegenes.org
    
    Parameters:
    -----------
    input_gtf_file : str input of the gtf file with path
    
    gtf_column_names : name of the standard columns froming the gtf file
    
    comment: str starting comment lines to skip
    
    Output
    ------
    
    Pandas dataframe containing the gtf file
    """
    gtf_file = pd.read_csv(input_gtf_file, # input gtf file
                       sep='\t', # how to separate the columns 
                       names=gtf_column_names, #column names
                       comment=comment, # skip lines starting with #
                       header=None,
                       index_col=False)
    return(gtf_file)
    
# 1.1 filter the gtf file

def filter_gtf_file(gtf_file, source=None, feature=None, copy=True):
    """
    filter gtf file based on tthe source of the annotations (HAVANA or ENSEMBL),
    and/or the type of feature (transcript, exon, gene)
    if source/feature is not None.
    
    It returns a new filtered pandas.DataFrame.
    
    Parameters
    ----------
    gtf_file : pandas.DataFrame containing the gtf file
    
    source : None or string ('HAVANA', 'ENSEMBL')
    
    feature : None or string ('transcript', 'exon', 'gene')
    
    copy : Boolean. If copy is True, return a new gtf file.
    Otherwise, the function overwrite the input gtf file.
    
    Output
    ------
    
    filtered pandas.DataFrame
    
    """
    filtered_gtf = gtf_file.copy()
    
    # filter the feature
    if feature != None:
        filter_crit =  filtered_gtf['feature']==feature
        filtered_gtf = filtered_gtf[filter_crit]
        filtered_gtf

    # filter the source
    if source != None:
        filter_crit =  filtered_gtf['source']==source
        filtered_gtf = filtered_gtf[filter_crit]
        filtered_gtf
        
    if copy:
        return(filtered_gtf)
    else:
        gtf_file.drop(gtf_file.index[~gtf_file.index.isin(filtered_gtf.index)], axis=0, inplace=True)


def extract_TSS(input_gtf_file):
    """
    extract key information regarding the TSS of a filtered gtf file. 
    """
    tss_dataframe =  {'chromosome': [], 'TSS_pos': [], 'strand' : [], 'annotation': []}
    for line in input_gtf_file.values:
    
        #Add chromosome
        tss_dataframe['chromosome'].append(line[0])
    
        #Add strand info
        tss_dataframe['strand'].append(line[6])
    
        # Add TSS coordinate depending of the strand
        if line[6]=='+':
            tss_dataframe['TSS_pos'].append(line[3])
        else:
            tss_dataframe['TSS_pos'].append(line[4])
        
        # Add gene/transcript name and id 
        tss_dataframe['annotation'].append([x for x in line[8].split(";") if line[2] in x])
    
    return(pd.DataFrame(data=tss_dataframe))


def find_TSS(peak_df, TSS_df):
    """
    extract the distance of peaks to TSS. Given the filtere dataframe of genes (filtered gtf)
    and the adata.var containing peak coordinates as chromosome, start_peak and end_peak
    If the df doesn't contain the right column. I need to adda a couple of command to extract the 
    information from the peak name.
    """    

    # Sort both df by position. 
    peak_infos = peak_df.sort_values('start_peak')
    TSS_chr1 = TSS_df.sort_values('TSS_pos')
    
    index = 0
    column_distance_to_TSS = []
    upstream_downstream = []
    gene_name = []
    peak_number = 0
    for peak in peak_infos.values:
        #print(peak)
        ##########
        if index >= len(TSS_chr1.values):
            TSS = TSS_chr1.values[-1]
            column_distance_to_TSS.append(min([peak[1]-TSS[1],peak[2]-TSS[1]]))
            upstream_downstream.append('+')
            gene_name.append(TSS[-1])
            print('outbound')
        ##########
        for TSS in TSS_chr1.values[index:]:
            peak[1] = int(peak[1])
            peak[2] = int(peak[2])
            if (peak[1] > TSS[1]) and (peak[1] > TSS_chr1.values[-1][1]):
                ############################################################
                #print('here')
                ############################################################
                column_distance_to_TSS.append(min([peak[1]-TSS[1],peak[2]-TSS[1]]))
                upstream_downstream.append('+')
                gene_name.append(TSS[-1])
                break
                
            elif peak[1] > TSS[1]:
                #print('index +1')
                index += 1
                
            elif peak[1] < TSS[1] < peak[2]:
                # peak overlapping the TSS so the distance is zero. 
                column_distance_to_TSS.append(0)
                upstream_downstream.append('+')
                gene_name.append(TSS[-1])
                break
            else:
                
                # these values can be positive or negative (with peak[2]), while peak[1] has to be positive
                distance_index = min([abs(peak[1]-TSS[1]),abs(peak[2]-TSS[1])]) 
                
                # this one has to be positive
                distance_index1 = min([abs(peak[1]-TSS_chr1.values[index-1][1]), 
                                       abs(peak[2]-TSS_chr1.values[index-1][1])])
                
                if index == len(TSS_chr1.values):
                    # these values can be positive or negative (with peak[2]) while peak[1] has to be negative
                    distance_index2 = min([abs(peak[1]-TSS_chr1.values[index+1][1]),
                                          abs(peak[2]-TSS_chr1.values[index+1][1])])
                else:
                    distance_index2=distance_index+1
    
                ##########
                if (distance_index < distance_index1):
                    if (distance_index < distance_index2):
                        # distance_index is the smallest
                        column_distance_to_TSS.append(distance_index)
                        if abs(peak[1]-TSS[1]) < abs(peak[2]-TSS[1]):
                            if peak[1]-TSS[1] < 0:
                                sign = '-'
                            else:
                                sign = '+'
                        else:
                            if peak[2]-TSS[1] < 0:
                                sign = '-'
                            else:
                                sign = '+'
                    else:
                        # distance_index2 is the smallest
                        column_distance_to_TSS.append(distance_index2)
                        if abs(peak[1]-TSS_chr1.values[index+1][1]) < abs(peak[2]-TSS_chr1.values[index+1][1]):
                            if peak[1]-TSS_chr1.values[index+1][1] < 0:
                                sign = '-'
                            else:
                                sign = '+'
                        else:
                            if peak[2]-TSS_chr1.values[index+1][1] < 0:
                                sign = '-'
                            else:
                                sign = '+'
                        index+=1
            
                elif (distance_index1 < distance_index2):
                    # distance_index1 is the smallest
                    column_distance_to_TSS.append(distance_index1)
                    if abs(peak[1]-TSS_chr1.values[index-1][1]) < abs(peak[2]-TSS_chr1.values[index-1][1]):
                        if peak[1]-TSS_chr1.values[index-1][1] < 0:
                            sign = '-'
                        else:
                            sign = '+'
                    else:
                        if peak[2]-TSS_chr1.values[index-1][1] < 0:
                            sign = '-'
                        else:
                            sign = '+'
                    index-=1
                else:
                    # distance_index2 is the smallest
                    column_distance_to_TSS.append(distance_index2)
                    if abs(peak[1]-TSS_chr1.values[index+1][1]) < abs(peak[2]-TSS_chr1.values[index+1][1]):
                        if peak[1]-TSS_chr1.values[index+1][1] < 0:
                            sign = '-'
                        else:
                            sign = '+'
                    else:
                        if peak[2]-TSS_chr1.values[index+1][1] < 0:
                            sign = '-'
                        else:
                            sign = '+'
                    index+=1
                upstream_downstream.append(sign)
                gene_name.append(TSS[-1])
                break
            
            if index >= len(TSS_chr1.values):
                #print('here here')
                index = len(TSS_chr1.values)-1
        
        
    peak_infos['TSS_distance'] = column_distance_to_TSS
    peak_infos['TSS_down_up_stream'] = upstream_downstream
    peak_infos['closest_TSS_name'] = gene_name
    return(peak_infos)

def find_TSS_subset_chromosome(peak_df, TSS_df):
    """
    Run the function find TSS chromosome per chromosome and save it in a new dataframe.
    Right now it returns the dataframe but in the end it should overwrite the anndata.var
    """
    columns = peak_df.columns.tolist()+['TSS_distance', 'TSS_down_up_stream', 'closest_TSS_name']
    peaks_infos_all_chr = pd.DataFrame(columns=columns)
    
    #chromosomes = sorted(list(set(TSS_df['chromosome'])))
    chromosomes = list(set(TSS_df['chromosome']))
    
    for chrom in chromosomes:
        #print(chrom)
        TSS_chr = TSS_df.copy()
        TSS_chr = TSS_chr[TSS_chr['chromosome']==chrom]

        peak_chr = peak_df.copy()
        peak_chr = peak_chr[peak_chr['chromosome']==chrom]
        
        if len(peak_chr.index) == 0:
            print('There is no peaks in ', chrom, '.')
            
        new_chrom_df = find_TSS(peak_chr, TSS_chr)
        peaks_infos_all_chr = peaks_infos_all_chr.append(new_chrom_df, ignore_index=True)
    
    new_index = []
    for row in peaks_infos_all_chr.values:
        new_index.append('_'.join([str(x) for x in row[0:3]]))
    peaks_infos_all_chr.index = new_index
    return(peaks_infos_all_chr)

def tool_distance2TSS(adata,
                      input_gtf_file,
                      source='ENSEMBL',
                      feature='transcript'):
     
    # if TSS infos in the adata.var I need to delete it. 'TSS_distance', 'TSS_down_up_stream'
    if 'TSS_distance' in adata.var.columns.tolist():
        del adata.var['TSS_distance']
    if 'TSS_down_up_stream' in adata.var.columns.tolist():
        del adata.var['TSS_down_up_stream']
    if 'closest_TSS_name' in adata.var.columns.tolist():
        del adata.var['closest_TSS_name']
    # if 'chromosome', 'start_peak', 'end_peak' are missing. I need to create them
    if 'chromosome' not in adata.var.columns.tolist():
        chromosome = []
        start_peak = []
        end_peak = []
        if ":" in adata.var_names[0]:
            for feature_name in adata.var_names:
                feature_name = feature_name.split(':')
                chromosome.append(feature_name[0])
                feature_name = feature_name[1].split('-')
                start_peak.append(int(feature_name[0]))
                end_peak.append(int(feature_name[1]))
        else:
            for feature_name in adata.var_names:
                feature_name = feature_name.split('_')
                chromosome.append(feature_name[0])
                start_peak.append(int(feature_name[1]))
                end_peak.append(int(feature_name[2]))
            
        adata.var['chromosome'] = chromosome
        adata.var['start_peak'] = start_peak
        adata.var['end_peak'] = end_peak
    
    gtf_file = load_gtf_file(input_gtf_file)
    filter_gtf_file(gtf_file, source=source, feature=feature, copy=False)

    TSS_df = extract_TSS(gtf_file)
    del gtf_file
    
    updated_peak_df_all_chr = find_TSS_subset_chromosome(adata.var.copy(), TSS_df)
    #copy_var = pd.merge(adata.var, updated_peak_df_all_chr)
    copy_var = pd.concat([adata.var, updated_peak_df_all_chr], axis=1)
    copy_var.index = adata.var_names.copy()
    adata.var = copy_var.copy()


## Finally ! mean distance per cell at open peaks !
def tool_mean_distance2TSS(adata):
    """
    mean_distance2TSS
    """
    mean_distance_cell = []
    for line in adata.X:
        mean_distance_cell.append(np.mean([adata.var['TSS_distance'][x] for x in line.indices.tolist()]))
    adata.obs['mean_distance_to_TSS'] = mean_distance_cell

