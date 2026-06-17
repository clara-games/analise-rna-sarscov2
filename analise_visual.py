import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

# --- 1. DADOS CONSOLIDADOS (Baseado no seu terminal) ---
variantes = ['Wuhan', 'Alfa', 'Delta', 'Omicron']
energias = [-82.10, -70.60, -89.94, -85.30] 
mutacoes = [0, 55, 188, 192] 

# --- 2. CONFIGURAÇÃO VISUAL (LIGHT MODE) ---
plt.style.use('default') # Volta para o padrão fundo branco
# Nova paleta: Azul Suave, Verde Água, Laranja Vibrante, Vermelho Intenso
cores = ['#2E86AB', '#3AB795', '#FF9100', '#FF2D55']

fig = plt.figure(figsize=(16, 10))
grid = fig.add_gridspec(2, 2)

fig.suptitle('RELATÓRIO DE EVOLUÇÃO GENÔMICA: SARS-CoV-2', fontsize=20, fontweight='bold', color='#333333')

# --- GRÁFICO 1: ESTABILIDADE (MFE) ---
ax1 = fig.add_subplot(grid[0, 0])
barras1 = ax1.bar(variantes, energias, color=cores, edgecolor='#555555', linewidth=1)
ax1.set_title('Estabilidade da Estrutura (Energia Livre)', fontsize=14, fontweight='bold')
ax1.set_ylabel('MFE (kcal/mol)')
ax1.grid(axis='y', linestyle='--', alpha=0.7)

for b in barras1:
    yval = b.get_height()
    ax1.text(b.get_x() + b.get_width()/2, yval - 1, f'{yval}', ha='center', va='top', fontweight='bold')

# --- GRÁFICO 2: ACÚMULO DE MUTAÇÕES (ANÁLISE EVOLUTIVA) ---
ax2 = fig.add_subplot(grid[0, 1])
barras2 = ax2.barh(variantes, mutacoes, color=cores, edgecolor='#555555', linewidth=1)
ax2.set_title('Distância Genética (Mutações vs Wuhan)', fontsize=14, fontweight='bold')
ax2.set_xlabel('Quantidade de Trocas de Nucleotídeos')
ax2.grid(axis='x', linestyle='--', alpha=0.7)

for i, v in enumerate(mutacoes):
    ax2.text(v + 3, i, f'{v} mut.', va='center', fontweight='bold')

# --- GRÁFICO 3: COMPOSIÇÃO DO RNA (T -> U CORRIGIDO) ---
ax3 = fig.add_subplot(grid[1, :])
def extrair_rna(caminho):
    for registro in SeqIO.parse(caminho, "fasta"):
        seq_rna = registro.seq.transcribe() # Garante a Uracila (U)
        return {base: seq_rna.count(base) for base in "AUCG"}

arquivos = ["wuhan_curto.fasta", "alfa_curto.fasta", "delta_curto.fasta", "omicron_curto.fasta"]

for i, arq in enumerate(arquivos):
    try:
        dados = extrair_rna(arq)
        ax3.plot(list(dados.keys()), list(dados.values()), marker='o', markersize=10, 
                 linewidth=3, label=variantes[i], color=cores[i])
    except:
        print(f"Erro ao ler {arq}")

ax3.set_title('Perfil Molecular: Transcrição de RNA', fontsize=14, fontweight='bold')
ax3.set_ylabel('Contagem de Bases')
ax3.legend(loc='upper right', frameon=True)
ax3.grid(alpha=0.3)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show() 