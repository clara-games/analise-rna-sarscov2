from Bio import SeqIO

def carregar_seq(arquivo):
    try:
        for registro in SeqIO.parse(arquivo, "fasta"):
            return str(registro.seq)
    except FileNotFoundError:
        return None

# Carregando as linhagens evolutivas
wuhan = carregar_seq("wuhan_curto.fasta")
alfa = carregar_seq("alfa_curto.fasta")
delta = carregar_seq("delta_curto.fasta")
omicron = carregar_seq("omicron_curto.fasta")

def rastro_evolutivo(referencia, variante, nome):
    if not referencia or not variante:
        return
    
    mutacoes = 0
    tamanho = min(len(referencia), len(variante))
    
    for i in range(tamanho):
        if referencia[i] != variante[i]:
            mutacoes += 1
            
    print(f"🧬 Linhagem {nome}: {mutacoes} mutações acumuladas em relação à Wuhan.")

print("="*50)
print("RELATÓRIO DE DERIVA GENÉTICA (TRAJETÓRIA EVOLUTIVA)")
print("="*50)

rastro_evolutivo(wuhan, alfa, "ALFA (B.1.1.7)")
rastro_evolutivo(wuhan, delta, "DELTA (B.1.617.2)")
rastro_evolutivo(wuhan, omicron, "OMICRON (BA.1)")
print("="*50) 