def open_swissprot():
    handle = open("uniprot-proteome.txt")
    swiss_record = [record for record in SwissProt.parse(handle)]
    return swiss_record
    
#Saves swissprot register  
def select_swissprot_entries():
    entries = list()
    for ft in features:
        encontrada = False
        if ft.type == 'CDS':
            for entry in open_swissprot():
                for ref in entry.cross_references:
                    if ref[0] == 'RefSeq':
                        if 'protein_id' in ft.qualifiers.keys():           #three entries lack protein_id
                            if ref[1] == ft.qualifiers['protein_id'][0]:
                                entries.append(entry)
    return entries

#Selects reviewed entries
def swissprot_reviewed():
    r_entries = list()
    for entry in sse:
        if entry.data_class == 'Reviewed':
            r_entries.append(entry)
    return r_entries

#Counts number of entries in genbank not found in swissprot
def check_features_in_swiss_prot():
    not_found = 0
    for ft in features:
        if ft.type=='CDS':
            found = False
            for ss in sse:
                if 'protein_id' in ft.qualifiers.keys():
                    found = True
            if found == False:
                not_found +=1
    return not_found


#Search in pdb database entries corresponding to study's proteins
def pdb():
    pdb_ids = list()
    f = open('pdb_database.txt')
    file = f.readlines()
    for ss in sse:
        for access in ss.accessions:
            for line in file:
                if access in line:
                    pdb_id = line.split()[0]
                    for ref in ss.cross_references:
                        if ref[0] == 'RefSeq':
                            pdb_ids.append((access,pdb_id))
    f.close()
    return pdb_ids
    
#Retrieve information from pdb
def create_pdb_files():
    pdbl = PDB.PDBList()
    for match in pdb():
        pdbl.retrieve_pdb_file(match)
