from Bio import SeqIO

def recortar(arquivo_entrada, nome_saida):
    # Abre o arquivo original
    for registro in SeqIO.parse(arquivo_entrada, "fasta"):
        # Recorta os primeiros 265 caracteres
        sequencia_curta = registro.seq[:265]
        
        # Salva o recorte num novo arquivo
        with open(nome_saida, "w") as f:
            f.write(f">{registro.id}_RECORTADO\n{sequencia_curta}\n")
    print(f"Pronto! Arquivo {nome_saida} criado.")

recortar("wuhan.fasta", "wuhan_curto.fasta")
recortar("omicron.fasta", "omicron_curto.fasta")
recortar("alfa.fasta", "alfa_curto.fasta")
recortar("delta.fasta", "delta_curto.fasta")