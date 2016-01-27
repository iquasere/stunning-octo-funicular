        
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
                                
                        

from Bio import AlignIO
interesting_list = ['WP_010881922.1','WP_010881929.1','WP_010881933.1','WP_010881936.1','WP_010881975.1','WP_010881980.1','WP_010881981.1','WP_010881982.1','WP_010881986.1','WP_010882011.1','WP_010882021.1','WP_010882033.1']

########  0  #########

def create_muscle(): #cria todos os ficheiros fasta obtidos por blast para formato clustal de forma a serem obtidos alinhamentos
    from Bio.Align.Applications._Muscle import MuscleCommandline
    try:
        lista=interesting_list
        muscle_exe = r"C:\workspace_python\Biopython\src\muscle.exe"
        for i in range(len(lista)):
            infile = (r"C:\workspace_python\Biopython\src\align" + str(i+1) + ".fasta")
            outfile = (r"C:\workspace_python\Biopython\src\Malign" + str(i+1) + ".aln")
            cline = MuscleCommandline(muscle_exe, input=infile, out = outfile, clw = True) #organiza de forma a ser um tipico clustalw
            cline()
    except:
        print("MultiAlignment Error")
        
########   1   #########

def readmulti_seq_id_dbxrefs(): #leitura dos 'features' dos alinhamentos multiplos
    try:
        lista=interesting_list
        for i in range(len(lista)):
            align = AlignIO.read(("Malign" + str(i+1) + ".aln"), "clustal")
            print("Alignment " + str(i+1) + " length %i" % align.get_alignment_length()) #tamanho do alinhamentos

            for record in align:
                print("%s - %s" % (record.seq, record.id)) #mostra a seq  + id proteina
                if record.dbxrefs: #caso contenha no ficheiro informacao de outras bases de dados
                    print("/n" "%s %s" % (record.id, record.dbxrefs)+ "/n") #id + referencia de outras bases de dados
    except:
        print("Reading error!")


########   1   #######

def muscle2phy(): #conversao de ficheiros tipo clustal para phylip de forma a obter posteriormente uma arvore de filogenia
    try:
        lista = interesting_list
        for i in range(len(lista)):
            AlignIO.convert(("Malign" + str(i+1) + ".aln"), "clustal", ("Malign" + str(i+1) + ".phy"), "phylip-relaxed")
        print("All MultiAlignments converted!")
    except:
        print("Converting Error!")

########   2   ########

def create_phy(): #criar um ficheiro de formato newick para cada alinhamento multiplo resultante de forma a ser possivel observar as arvores filogeneticas
    from Bio.Phylo.Applications import PhymlCommandline
    try:
        phy_exe = r"C:\workspace_python\Biopython\src\PhyML-3.1\PhyML.exe"
        lista = interesting_list
        for i in range(len(lista)):
            input = (r"C:\workspace_python\Biopython\src\Malign" + str(i+1) + ".phy")
            cmdline = PhymlCommandline(phy_exe, input=input, datatype='aa', model='WAG', alpha='e', bootstrap=20)
            out_log, err_log = cmdline()
        print("Created All Tree files!")
    except:
        print("Error creating tree files!")

###########  3   ###########

def tree(): #obter as arvores filogeneticas em formato ascii
    from Bio import Phylo
    try:
        lista = interesting_list
        for i in range(len(lista)):
            align_tree = Phylo.read(("Malign" + str(i+1) + ".phy_phyml_tree.txt"), "newick")
            print("Happy tree " + str(i+1) +"! " + (interesting_list[i]))
            Phylo.draw_ascii(align_tree)
    except:
        print("Creating Tree error!")

##########   4   #########

def export(): #abre uma janela com a arvore filogenetica resultante e ser possivel guarda-la de forma a exportar para o site
    import pylab
    from Bio import Phylo
    try:
        lista = interesting_list
        for i in range(len(lista)):
            align_tree = Phylo.read(("Malign" + str(i+1) + ".phy_phyml_tree.txt"), "newick")
            Phylo.draw(align_tree)
            pylab.show()
    except:
        print("Vizualizing Tree error!")
