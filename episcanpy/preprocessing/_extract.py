def extract_CG(input_sample, output_variable, head):
    """
    Description function. 
    
    Extract relevant information for cytosines covered in a CG context
    It supposes that it is on the same format at the Methylpy/Ecker data
    
    Parameters
    ----------
    input_sample: file opened that contain the cytosines
    output_variable: the datatype that gather the methylation levels for a
        given type of features
    head: annotation of the file that we are reading (line to skip)
    """
    for cyt_called in input_sample:
        cyt_called = cyt_called.split("\t")
        if (cyt_called != head) and (cyt_called[3][0:2] == 'CG'):
            output_variable[cyt_called[0]].append([int(cyt_called[1]),
                       int(cyt_called[-3]), int(cyt_called[-2])])
    return(output_variable)
    
def extract_CH(sample, output_variable, head):
    """
    Description function. 
    
    Extract relevant information for cytosines covered in a CH context
    It supposes that it is on the same format at the Methylpy/Ecker data

    Parameters
    ----------
    input_sample: file opened that contain the cytosines
    output_variable: the datatype that gather the methylation levels for a
        given type of features
    head: annotation of the file that we are reading (line to skip)
    """
    for cyt_called in sample:
        cyt_called = cyt_called.split("\t")
        if cyt_called != head and cyt_called[3][0:2] != 'CG':
            output_variable[cyt_called[0]].append([int(cyt_called[1]),
                       int(cyt_called[-3]), int(cyt_called[-2])])
    return(output_variable)
