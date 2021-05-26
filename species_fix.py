def sp_fix():
    '''
    This function changes the short species names in the phtozome details file and changes them
    to the full Genus and species name from the Phytozome_species_index.csv file, and saves the
    output as phytozome_details_mod_all_sp_corrected.txt

    This function takes no arguments
    '''

    import pandas as pd
    from log_func import log

    details_modDF = pd.read_csv('phytozome_details_mod_all.txt', sep='\t')
    speciesDF = pd.read_csv('Phytozome_species_index.csv')
    speciesDF.rename(columns={'Short_species': 'species'}, inplace=True)
    detailsDF = pd.merge(speciesDF, details_modDF, how='inner', on='species')
    del detailsDF['species']
    detailsDF.rename(columns={'Full_species': 'Species'}, inplace=True)
    detailsDF.to_csv('phytozome_details_mod_all_sp_corrected.txt', sep='\t', index=False)

    log(f'phytozome_details_mod_all_sp_corrected.txt was created with {len(detailsDF.index)} rows')
