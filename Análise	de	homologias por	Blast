
def retrieving_info_blast():   
    for prot in ['WP_010881922.1',
 'WP_010881929.1',
 'WP_010881933.1',
 'WP_010881936.1',
 'WP_010881975.1',
 'WP_010881980.1',
 'WP_010881981.1',
 'WP_010881982.1',
 'WP_010881986.1',
 'WP_010882011.1',
 'WP_010882021.1',
 'WP_010882033.1']:
        result = open('int'+prot+'.xml') 
        blast_record = NCBIXML.read(result)
        for alignment in blast_record.alignments:
            print('Acession:',blast_record.accession)
            print('ID:',blast_record.hit_id)
            print('Description:',blast_record.hit_def)
            print('Alignment:',blast_record)
            for hsp in blast_record.hsps:
                print('e-value:',hsp.expect)
                print('Score:',hsp.score)
                print('Match:',hsp.match)
 
def proteinas_interesse():
    '''lista de proteinas de possivel interesse, ou seja são as proteinas reviewed'''
    ids = []
    for ss in swissprot_reviewed(select_swissprot_entries(select_features(),"uniprot-proteome.txt")):
        for ref in ss.cross_references:
            if ref[0] == 'RefSeq':
                ids.append(ref[1])
    return ids
                        
def poss():
    '''função que cria o fasta para se fazer o blast'''
    interesting_list=proteinas_interesse()
    lista=select_features()
    j=1
    for ft in lista:
        qual=ft.qualifiers
        if 'protein_id' in qual.keys():                        
            if qual['protein_id'][0] in interesting_list:
            #ver se algum das nossas proteinas tem um id dos que estão na lista de interesse
                with open('poss'+str(j)+'.fasta','w') as save_file:
                    qual=ft.qualifiers
                    if 'translation' in qual.keys():
                        save_file.write('>'+qual['protein_id'][0]+'\n')                
                        save_file.write(qual['translation'][0]+'\n')
                        j+=1
                        
def search_interesse():
    '''função do blast das proteinas possiveis de interesse contra humano e e.coli'''
    interesting_list=proteinas_interesse()
    j=1
    for prot in interesting_list:
        open('poss'+str(j)+'.fasta')
        record=SeqIO.read('poss'+str(j)+'.fasta','fasta')
        result_handle = NCBIWWW.qblast('blastp','swissprot',record.format('fasta'),hitlist_size=1, entrez_query='Homo sapiens[organism]'and'Escherichia coli[organism]')
        save_file=open('possib'+str(j)+'.xml','w')
        save_file.write(result_handle.read())
        save_file.close()
        result_handle.close()
        j+=1

def determinar_interesse():
    '''determinação do interesse da proteina, se a identidade do unico hit que foi obtido for menos que 20, é considerado uma proteina de interesse'''
    interesting_list=[]
    list=proteinas_interesse()
    for prot in range(len(list)):
        result = open('possib'+str(prot+1)+'.xml')
        blast_record = NCBIXML.parse(result)
        item=next(blast_record)
        for ali in item.alignments:
            for hsp in ali.hsps:
                if hsp.identities < 20:
                    interessantes.append(prot)
    return interesting_list

def separar_int():
    '''construção de uma lista de proteinas de interesse'''
    num_int=determinar_interesse()
    list_int=proteinas_interesse()
    interesting_list=[]
    for i in num_int:
        a=list_int[i]
        interesting_list.append(a)
    return interesting_list
    
def separar_n_int():
    '''construção de uma lista com as restantes proteinas'''
    interesting_list=separar_int()
    features=select_features()
    n_int=[]
    for ft in features:
        qual=ft.qualifiers
        if 'protein_id' in qual.keys():                        
            if qual['protein_id'][0] not in interesting_list:
                n_int.append(qual['protein_id'][0])
    return n_int

def fasta_files():
    '''função que cria os fastas divisorios entre interessantes e outras'''
    interesting_list=separar_int()
    features=select_features()
    for ft in features:
        qual=ft.qualifiers
        if 'protein_id' in qual.keys():                        
            if qual['protein_id'][0] in interesting_list:
            #ver se algum das nossas proteinas tem um id dos que estão na lista de interesse
                with open('int'+qual['protein_id'][0]+'.fasta','w') as save_file:
                    qual=ft.qualifiers
                    if 'translation' in qual.keys():
                        save_file.write('>'+qual['protein_id'][0]+'\n')                
                        save_file.write(qual['translation'][0]+'\n')
            else:
                with open('ni'+qual['protein_id'][0]+'.fasta','w') as save_file:
                    qual=ft.qualifiers
                    if 'translation' in qual.keys():
                        save_file.write('>'+qual['protein_id'][0]+'\n')                
                        save_file.write(qual['translation'][0]+'\n')

def blast():
    '''blast para todas as proteinas'''
    interesting_list=separar_int()
    features=select_features()    
    for ft in features:
        qual=ft.qualifiers
        if 'protein_id' in qual.keys():                        
            if qual['protein_id'][0] in interesting_list:
                open('int'+qual['protein_id'][0]+'.fasta')
                record=SeqIO.read('int'+qual['protein_id'][0]+'.fasta','fasta')
                result_handle = NCBIWWW.qblast('blastp','swissprot',record.format('fasta'), hitlist_size=10)
                save_file=open('int'+qual['protein_id'][0]+'.xml','w')
                save_file.write(result_handle.read())
                save_file.close()
                result_handle.close()
            else:
                open('ni'+qual['protein_id'][0]+'.fasta')
                record=SeqIO.read('ni'+qual['protein_id'][0]+'.fasta','fasta')
                result_handle = NCBIWWW.qblast('blastp','swissprot',record.format('fasta'), hitlist_size=10)
                save_file=open('ni'+qual['protein_id'][0]+'.xml','w')
                save_file.write(result_handle.read())
                save_file.close()
                result_handle.close()

