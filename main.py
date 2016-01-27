def main():
    #RETRIEVING DATA
    handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id="6273291")
    seq_record_fasta = SeqIO.read(handle, "fasta")
    handle.close()
    handler = Entrez.efetch(db="nucleotide", rettype="gbwithparts", retmode="text", id="NC_000919")
    seq_record_genbank = SeqIO.read(handler, "gb") #using "gb" as an alias for "genbank"
    handler.close()
    features = select_features()
    tabular = readtabular()
    sse = select_swissprot_entries()
    sser = swissprot_reviewed()
    print('---Analysis of information---')
    print('|GenBank|')
    print('There are',len(seq_record_genbank),'entries in the file retrieved from GenBank')
    print(len(select_features()),'features are to be analyzed in this work.')
    print(distinguish_feature_type()[0],'entries are CDS and',distinguish_feature_type()[1],'entries are tRNA.')
    print('|SwissProt|')
    print(len(open_swissprot()), 'entries in the file retrieved from SwissProt.')
    print(len(sse),'features from SwissProt are to be analyzed in this work.')
    print('|Protein Data Bank|')    
    print(len(pdb()),'entries from SwissProt exist in the Protein Data Bank.')
    print('|KEGG|')
    print(make_kegg(),'entries from SwissProt exist in KEGG.')
    print()
    print('---Validation of information---')
    print('|Tabular vs GenBank registers|')
    print('There were',validate_tb_gb()[0],'non similiar informations between tabular and genbank registers')
    print('Those non similar informations were:')
    for irre in validate_tb_gb()[1]: print(irre)
    print('|Swiss Prot|')
    print(check_features_in_swiss_prot(),'protein entries from the GenBank file are not featured in the SwissProt register.')
    print(str(round(len(sser)/len(sse)*100,1))+'% of the entries in SwissProt file is reviewed data.')
    interesting_list = ['WP_010881922.1','WP_010881929.1','WP_010881933.1','WP_010881936.1','WP_010881975.1','WP_010881980.1','WP_010881981.1','WP_010881982.1','WP_010881986.1','WP_010882011.1','WP_010882021.1','WP_010882033.1']    
    CDD_results = {'WP_010881922.1':'No conserved domains',\
'WP_010881929.1':'<a href="http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?ascbin=8&maxaln=10&seltype=2&uid=TIGR03546"> CDD links</a>',\
'WP_010881933.1':'<a href="http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?ascbin=8&maxaln=10&seltype=2&uid=pfam04773">CDD Link</a>',\
'WP_010881936.1':'No conserved domains',\
'WP_010881975.1':'<a href="http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?ascbin=8&maxaln=10&seltype=2&uid=COG1394">CDD Link</a>',\
'WP_010881980.1':'<a href="http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?ascbin=8&maxaln=10&seltype=2&uid=COG1269">CDD Link</a>',\
'WP_010881981.1':'<a href="http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?ascbin=8&maxaln=10&seltype=2&uid=COG1527">CDD Link</a>',\
'WP_010881982.1':'No conserved domains',\
'WP_010881986.1':'No conserved domains',\
'WP_010882011.1':'<a href="http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?ascbin=8&maxaln=10&seltype=2&uid=COG1480">CDD Link</a>',\
'WP_010882021.1':'<a href="http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?ascbin=8&maxaln=10&seltype=2&uid=cl20779">CDD Link</a>',\
'WP_010882033.1':'No conserved domains'}
    

if __name__ == '__main__':
    main()
