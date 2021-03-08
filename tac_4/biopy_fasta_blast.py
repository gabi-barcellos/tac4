"""
Autora:  Gabriela Barcellos
Date: 05/03/2021
Version: 3
Last update: 07/03/2021
"""

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline

dna=str(input("Digite uma sequencia de DNA: "))
dna=Seq(dna.upper())
print("Essa é a sua sequencia de DNA: {}.".format(dna))
mrna = dna.transcribe()
print("Essa é a sua sequencia de mRNA: {}.".format(mrna))
tradução=mrna.translate()
print("Essa é a sua proteina: {}".format(tradução))

#PARTE FASTA

arquivo_entrada = open(r'C:/Users/gabri/Downloads/sequencias.fasta')
contador = 0
for i in SeqIO.parse(arquivo_entrada, "fasta"):
    contador+=1
    arquivo = open(r'sequencia' + str(contador) + '.fasta', 'w+')
    arquivo.write(str(i.seq))
    arquivo.close()



#PARTE BLAST
sequencia_desconhecida = input("Digite o caminho do arquivo fasta: ")
base_dados = input("Digite o caminho do arquivo multi-fasta de proteínas: ")
meuOutput = r"C:/Users/gabri/PycharmProjects/tac4/tac_4/meuOutput.txt"
blastx_path = r"C:/Arquivos de programas/NCBI/blast-2.11.0+/bin/blastx.exe"

linha_de_comando_blastx = NcbiblastxCommandline(query= sequencia_desconhecida, ## o que quero buscar
                                                subject = base_dados,  ##  onde quero buscar
                                                outfmt=6, ## formato do arquivo de saida
                                                out=meuOutput, ## arquivo de saida
                                                evalue=0.05, ## e-valor
                                                cmd = blastx_path) ## caminho para o executável
print("Executando busca local: ")

stdout, stderr = linha_de_comando_blastx()

print("Fim da busca local...")

# Abrindo o arquivo de output que armazena o resultado
blast_result = open(meuOutput, "r")

##### indíces para os resultados do blast em formato 6 ("outfmt=6")


evalue = 10 # expect value
bitscore = 11 # bit score

## Imprimindo os resultados
results = blast_result.read()
s = results.split('\n')

for linha in s:

    hit = linha.split('\t')
    if len(hit) != 12:
        break

    #hit[bitscore] = list(map(int, hit[bitscore]))
    #print("Score = %s" % hit[bitscore])
    #maximo = hit[bitscore].max()
    print("Score = {}".format(hit[bitscore]))

blast_result.close()