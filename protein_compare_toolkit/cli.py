'''The main CLI module.'''
import typer

from protein_compare_toolkit.commands import logo
from protein_compare_toolkit.commands import stats

app = typer.Typer()

app.add_typer(stats.app, name="stats")
app.add_typer(logo.app)

if __name__ == "__main__":
    app()
