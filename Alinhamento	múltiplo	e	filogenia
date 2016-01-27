        
def align_files():
    '''A função vai criar um ficheiro fasta que contenha a nossa sequencia e os seus hits do blast necessário para os alinhamentos multiplos'''
    lista_int=separar_int()  #lista das proteinas de interesse
    lista=select_features()  #lista de features do genbank
    for ft in lista:
        qual=ft.qualifiers
        if 'protein_id' in qual.keys():                        
            if qual['protein_id'][0] in lista_int:
                for i in range(len(lista_int)):
                    with open('align'+str(i+1)+'.fasta','w') as save_file:#vai abrir um ficheiro fasta
                        save_file.write('>'+qual['protein_id'][0]+'\n')                
                        save_file.write(qual['translation'][0]+'\n')#colocar a nossa sequencia
                        result = open('int'+qual['protein_id']+'.xml')
                        blast_record = NCBIXML.read(result)
                        for alignment in blast_record.alignments:
                            for hsp in alignment.hsps:
                                save_file.write('\n'+'>'+alignment.title+'\n')
                                save_file.write(hsp.sbjct+'\n')#e colocara sequencia do hit
