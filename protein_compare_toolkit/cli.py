import typer
from protein_compare_toolkit.commands import sdi

app = typer.Typer()
app.add_typer(sdi.app, name="sdi")

if __name__ == "__main__":
    app()
