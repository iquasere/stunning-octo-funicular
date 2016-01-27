def tabela():
    with open('nova_tabela.txt','w') as f_nt:
        f_nt.write('<!DOCTYPE html><html><body><table style="width:100%"><tr>')
        cabeçalho = ['NCBI GeneID','G_Acession Number','Locus Tag','Gene name','Strand','UniprotID','Revision degree','NCBI Protein_Acession Number','Molecule name','Sequence length','Cellular location','GeneOntology terms','EC/TC Number','Description']
        linha = str()    
        for cabeça in cabeçalho:
            linha += '<th>'+cabeça+'</th>'
        f_nt.write(linha)
        for ft in features:
            linha = '<tr>'
            if ft.type == 'CDS':                                                                      #exclude gene and tRNA entries
                if 'protein_id' in ft.qualifiers.keys():                                              #three hypotetical proteins lack protein id
                    if ft.qualifiers['protein_id'][0] in interesting_list:
                        linha = '<tr bgcolor="#00e5e6">'
                    for entry in sse:
                        for ref in entry.cross_references:
                            if ref[0] == 'RefSeq':
                                if ref[1] == ft.qualifiers['protein_id'][0]:                          #joins both registers (genbank and unitprot)
                                    found = re.search('GeneID:(.+)', ft.qualifiers['db_xref'][1])
                                    if found != None:
                                        first = found.group(1)
                                    else:
                                        linha += '<td>'+'-'+'</td>'
                                    second = ref[2]
                                    third = ft.qualifiers['locus_tag'][0]
                                    fourth = []
                                    found = re.search('OrderedLocusNames=(TP_[0-9]+);', entry.gene_name)
                                    if found != None:
                                        fourth.append(found.group(1))
                                        found = re.search('Name=(.+?)\s', entry.gene_name)
                                        if found != None:
                                            fourth.append(found.group(1))
                                            found = re.search('Synonyms=(.+?);', entry.gene_name)
                                            
                                    else:
                                        fourthstr = '-'
                                    fourthstr = ''
                                    for name in fourth:
                                        fourthstr  += name + '; '
                                    fifth = ft.location.strand
                                    sixth = entry.accessions[0]
                                    seventh = entry.data_class
                                    eighth = ref[1]
                                    ninth = ft.qualifiers['product'][0]
                                    tenth = entry.sequence_length
                                    if len(entry.comments) > 0:
                                        found = re.search('SUBCELLULAR LOCATION: (.+?) {', entry.comments[0])                                    
                                        if found != None:
                                            eleventh = found.group(1)
                                        else:
                                            eleventh = '-'
                                    else:
                                        eleventh = '-'
                                    twelfth = ''                                    
                                    for ref in entry.cross_references:
                                        if ref[0] == 'GO':
                                            twelfth += ref[1]+' '+ref[2]+' '+ref[3]+'; '
                                    twelfth = twelfth.rstrip('; ')
                                    found = re.search('(EC=.+?) {', entry.description)
                                    if found != None:
                                        thirteenth = found.group(1)
                                    else:
                                        found = re.search('(TC=.+?) {', entry.description)
                                        if found != None:
                                            thirteenth = found.group(1)
                                        else:
                                            thirteenth = '-'
                                    if len(entry.comments) > 0:
                                        found = re.search('FUNCTION: (.+?) {', entry.comments[0])                                    
                                        if found != None:
                                            fourteenth = found.group(1) 
                                        else:
                                            fourteenth = '-'
                                    else:
                                        fourteenth = '-'
                    linha += '<td>'+str(first)+'</td>'                          #Gene ID NCBI
                    linha += '<td>'+second+'</td>'                              #Accession NCBI
                    linha += '<td>'+third+'</td>'                               #locus tag
                    linha += '<td>'+fourthstr.rstrip('; ')+'</td>'              #gene name
                    linha += '<td>'+str(fifth)+'</td>'                          #strand
                    linha += '<td>'+sixth+'</td>'                               #UniProt ID         
                    linha += '<td>'+seventh+'</td>'                             #Revision
                    linha += "<td><a href='blastsite_"+eighth+".html'>"+eighth+"</a></td>"       #Protein NCBI accession
                    linha += '<td>'+ninth+'</td>'                               #Protein name
                    linha += '<td>'+str(tenth)+'</td>'                          #Sequence length
                    linha += '<td>'+eleventh+'</td>'                            #cellular location
                    linha += '<td>'+twelfth+'</td>'                             #gene ontology terms
                    linha += '<td>'+thirteenth+'</td>'                          #EC or TC numbers
                    linha += '<td>'+fourteenth+'</td>'                          #Function
                    f_nt.write(linha)
                else:    
                    found = re.search('GeneID:(.+)', ft.qualifiers['db_xref'][0])
                    if found != None:
                        first = found.group(1)
                        linha += '<td>'+str(first)+'</td>'                          #Gene ID NCBI
                    else:
                        linha += '<td>'+'-'+'</td>'
                    second = '-'
                    linha += '<td>'+second+'</td>'                                  #Accession NCBI
                    third = ft.qualifiers['locus_tag'][0]
                    linha += '<td>'+third+'</td>'                                  #locus tag
                    fourth = '-'
                    linha += '<td>'+str(fourth)+'</td>'                                 #gene name
                    fifth = ft.location.strand
                    linha += '<td>'+str(fifth)+'</td>'                                   #strand
                    sixth = '-'
                    linha += '<td>'+sixth+'</td>'                                   #UniProt ID
                    seventh = '-'
                    linha += '<td>'+seventh+'</td>'                                 #Revision
                    eighth = '-'
                    linha += '<td>'+eighth+'</td>'                                  #Protein NCBI accession
                    ninth = ft.qualifiers['product'][0]
                    linha += '<td>'+ninth+'</td>'                                   #Protein name
                    tenth = int((ft.location.end-ft.location.start-3)/3)
                    linha += '<td>'+str(tenth)+'</td>'                                   #Sequence length
                    eleventh = '-'
                    linha += '<td>'+eleventh+'</td>'                                #cellular location
                    twelfth = '-'
                    linha += '<td>'+twelfth+'</td>'                                #gene ontology terms
                    thirteenth = '-'
                    linha += '<td>'+thirteenth+'</td>'                              #EC or TC numbers
                    fourteenth = ninth
                    linha += '<td>'+fourteenth+'</td>'                              #Function
                    f_nt.write(linha)
            elif ft.type == 'trna':
                found = re.search('GeneID:(.+)', ft.qualifiers['db_xref'][0])            
                if found != None:
                    first = found.group(1)
                    linha += '<td>'+str(first)+'</td>'                          #Gene ID NCBI
                else:
                    linha += '<td>'+'-'+'</td>'
                second = '-'
                linha += '<td>'+second+'</td>'                                  #Accession NCBI
                third = ft.qualifiers['locus_tag'][0]
                linha += '<td>'+third+'</td>'                                  #locus tag
                fourth = '-'
                linha += '<td>'+str(fourth)+'</td>'                                 #gene name
                fifth = ft.location.strand
                linha += '<td>'+str(fifth)+'</td>'                                   #strand
                sixth = '-'
                linha += '<td>'+sixth+'</td>'                                   #UniProt ID
                seventh = '-'
                linha += '<td>'+seventh+'</td>'                                 #Revision
                eighth = '-'
                linha += '<td>'+eighth+'</td>'                                  #Protein NCBI accession
                ninth = ft.qualifiers['product'][0]
                linha += '<td>'+ninth+'</td>'                                   #RNA name
                tenth = int(ft.location.end-ft.location.start)
                linha += '<td>'+str(tenth)+'</td>'                                   #Sequence length
                eleventh = '-'
                linha += '<td>'+eleventh+'</td>'                                #cellular location
                twelfth = '-'
                linha += '<td>'+twelfth+'</td>'                                #gene ontology terms
                thirteenth = '-'
                linha += '<td>'+thirteenth+'</td>'                              #EC or TC numbers
                fourteenth = '-'
                linha += '<td>'+fourteenth+'</td>'                              #Function
                f_nt.write(linha)
            linha += '</tr>'
        f_nt.write('</tr></table></body></html>')            
    f_nt.close()
