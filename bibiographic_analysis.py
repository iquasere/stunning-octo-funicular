def protein_literature():
    pp_list=[]
    for ft in features:
        if 'product' in ft.qualifiers.keys():
            if ft.qualifiers['product'][0] != 'hypothetical protein' and 'tRNA' not in ft.qualifiers['product'][0]:
                pp_list.append(ft.qualifiers['product'][0])
    p_literature = dict()
    for product in pp_list:
        handle=Entrez.esearch(db='pubmed',term=product,retmax=2)
        record=Entrez.read(handle)
        idlist=record['IdList']
        Handle=Entrez.efetch(db='pubmed',id=idlist, rettype='medline',retmode='text')
        records=Medline.parse(Handle)
        p_literature['product'] = records
    return p_literature
