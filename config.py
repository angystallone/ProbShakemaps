import os

filename = 'input_file.txt'

def load_config(file):
    with open(os.path.join(os.getcwd(), "INPUT_FILES", filename), 'r') as file:
        lines = file.readlines()

    # Extract configuration parameters
    tectonicRegionType = lines[0].split(maxsplit=1)[1:][0].rstrip()
    mag_scaling = str(lines[1].split( )[1]) 
    rupture_aratio = int(lines[2].split( )[1])
    ID_Event = str(lines[3].split( )[1])
    if str(lines[4].split( )[1]) == 'None':
        vs30file = None
    else:
        vs30file = str(lines[4].split( )[1])
    CorrelationModel = str(lines[5].split( )[1])
    CrosscorrModel = str(lines[6].split( )[1])
    vs30_clustering = lines[7].split()[1] == 'True'
    truncation_level = float(lines[8].split( )[1])
    
    # Create and return configuration dictionary
    config = {
        'tectonicRegionType': tectonicRegionType,
        'mag_scaling': mag_scaling,
        'rupture_aratio': rupture_aratio,
        'ID_Event': ID_Event,
        'vs30file': vs30file,
        'CorrelationModel': CorrelationModel,
        'CrosscorrModel': CrosscorrModel,
        'vs30_clustering': vs30_clustering,
        'truncation_level': truncation_level
    }
    
    return config
