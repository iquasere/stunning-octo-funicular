#counts number of features of CDS and tRNA type
def distinguish_feature_type():
    trna = 0
    protein = 0
    for ft in features:
        if ft.type=='CDS':
            protein += 1
        elif ft.type == 'tRNA':
            trna += 1
    return (protein,trna)

#Selects features with specific locus_tag
def select_features():
    f = seq_record_genbank.features
    features = list()
    for i in range(2345,2951, 5):
        for ft in f:
            if 'locus_tag' in ft.qualifiers.keys():
                if '0'*(5-int(len(str(i))))+str(i) in ft.qualifiers['locus_tag'][0]:
                    features.append(ft)
    return features
    
#Saves tabular register
def readtabular():
    handler = open('tabular.txt') #Opens tabular file
    tabular = list() #Where tabular register is saved
    titles = handler.readline().lstrip('#').rstrip('\n').split('\t') #headlines
    for line in handler.readlines():
        objects = line.rstrip('\n').split('\t')
        for j in range(2345,2951,5):
            if str(j) in objects[7]: #searches for the right locus_tag
                register = dict() #Saves each register as a dictionary
                for k in range(len(objects)):
                    register[titles[k]] = objects[k]
                tabular.append(register)
    return tabular

#compares information in genbank with information in tabular    
def validate_tb_gb():
    irregularities = 0
    irregular_list = list()
    for ft in features:
        if ft.type == 'CDS':
            for t in tabular:
                if t['Locus tag'] == ft.qualifiers['locus_tag'][0]:
                    if 'GeneID:'+t['GeneID'] != str(ft.qualifiers['db_xref'][1]):       #verify geneID
                        irregularities += 1
                    if int(t['Length']) != len(ft.qualifiers['translation'][0]):        #verify length
                        irregularities += 1
                    if ft.qualifiers['product'][0].upper() not in t['Protein name'].upper():            #verify protein name -> problem: some proteins don't have exact same name
                         irregularities += 1
                         irregular_list.append(ft.qualifiers['product'][0]+' =/= '+t['Protein name'])
                    if t['Protein product'] != ft.qualifiers['protein_id'][0]:          #verify protein identity
                        irregularities += 1
                    #t['Replicon Accession'] not identified in features
                    if str(int(t['Start'])-1) != str(int(str(ft.location.start))):      #verify start position -> problem: start in features always 1 less
                        irregularities += 1
                    if t['Stop'] != str(ft.location.end):                               #verify stop position
                        irregularities += 1        
                    if t['Strand'] == '+':                                              #verify strand
                        if ft.location.strand != 1:
                            irregularities += 1
                    else:
                        if ft.location.strand != -1:
                            irregularities += 1
    return irregularities, irregular_list
