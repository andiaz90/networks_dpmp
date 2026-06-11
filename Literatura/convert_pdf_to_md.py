from pathlib import Path
from markitdown import MarkItDown

# Carpeta con los PDFs
carpeta_pdf = Path("/Users/valentinavasquez/Documents/GitHub/networks_dpmp/Literatura/network_labor_mkt_pdf")

# Carpeta donde se guardarán los .md (se crea si no existe)
carpeta_md = carpeta_pdf.parent / "network_labor_mkt_md"
carpeta_md.mkdir(exist_ok=True)

md = MarkItDown()

pdfs = sorted(carpeta_pdf.glob("*.pdf"))
print(f"Se encontraron {len(pdfs)} PDFs\n")

for pdf in pdfs:
    try:
        resultado = md.convert(str(pdf))
        archivo_md = carpeta_md / f"{pdf.stem}.md"
        archivo_md.write_text(resultado.text_content, encoding="utf-8")
        print(f"✓ {pdf.name} -> {archivo_md.name}")
    except Exception as e:
        print(f"✗ Error con {pdf.name}: {e}")

